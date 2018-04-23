module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                                  ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec, xtest                       ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22, log_tau          ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2                          ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22                         ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none

    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, T_b, n_b, dx, xstart, xstop, C_r, constant,yp1, ypn, eps, hmin
    logical(lgt) :: use_saha

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo

    ! Spline variables
    yp1 = 1.d30
    ypn = 1.d30

    ! Integration variables
    eps = 1.d-10
    hmin = 0.d0


    ! Allocating arrays
    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(log_tau(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    allocate(xtest(10000))



    !-----------------------------------------------------------------------------
    ! x (rec) grid
    !-----------------------------------------------------------------------------
    dx = (xstop-xstart)/(n-1)
    x_rec(1) = xstart
    do i=1, (n-1)
       x_rec(i+1) = x_rec(i) + dx
    end do
    
    !-----------------------------------------------------------------------------
    ! X_e and n_e at all grid times
    !-----------------------------------------------------------------------------
    use_saha = .true.
    do i = 1, n
       n_b = Omega_b0*rho_c0/(m_H*exp(x_rec(i))**3)

       if (use_saha) then
          !-----------------------------------------------------------------------------
          ! Saha equation for X_e > 0.99
          !-----------------------------------------------------------------------------
          T_b = T_0/exp(x_rec(i))                                                       
          constant = (m_e*k_b*T_b/(2.d0*pi*hbar*hbar))**1.5*exp(-epsilon_0/(k_b*T_b))/n_b ! Constant for quadratic formula
          X_e(i) = (-constant + sqrt(constant*constant + 4.d0*constant))/2.d0             ! Quadratic formula
          if (X_e(i) < saha_limit) use_saha = .false.                                     ! Using Peeble's equation
            
       else
          X_e(i) = X_e(i-1)                                                               ! Integrating forward
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, dx, hmin, dXe_dx, bsstep, output)
       end if
       n_e(i) = X_e(i)*n_b                                                                ! Calculating electron density
    end do
    n_e = log(n_e)                                                                        ! Log for future interpolation benefits

    call spline(x_rec,n_e,yp1,ypn,n_e2)                                                   ! Splined (log of) electron density function

    !-----------------------------------------------------------------------------
    ! Optical depth at all grid points 
    !-----------------------------------------------------------------------------
    tau(n) = 0.d0                                                                         ! Initial tau (today's value)
    log_tau(n) = -18.7d0                                                                  ! Educated guess from the n-1 value of tau
    do i = n-1, 1, -1
        tau(i) = tau(i+1)                                                                 ! Integrating forward
        call odeint(tau(i:i), x_rec(i+1) ,x_rec(i), eps ,dx, hmin, dtau_dx, bsstep, output) 
        log_tau(i) =log(tau(i))                                                           ! Setting log(tau) values for future interpolation
    end do

    !-----------------------------------------------------------------------------
    ! Splined (log of) optical depth
    !-----------------------------------------------------------------------------
    call spline(x_rec, log_tau, yp1, ypn, tau2)

    !-----------------------------------------------------------------------------
    ! Splined second derivative of (log of) optical depth
    !-----------------------------------------------------------------------------
    call spline(x_rec, tau2, yp1, ypn, tau22)


    !-----------------------------------------------------------------------------
    ! Visibility function (g_tilde) at all grid points
    !-----------------------------------------------------------------------------
    do i=1,n      
      g(i) = -get_dtau(x_rec(i))*exp(-tau(i))
    end do

    !-----------------------------------------------------------------------------
    ! Splined visibility function
    !-----------------------------------------------------------------------------
    call spline(x_rec, g, yp1, ypn, g2)

    !-----------------------------------------------------------------------------
    ! Splined second derivative of visibility function
    !-----------------------------------------------------------------------------
    call spline(x_rec,g2,yp1,ypn,g22)

    !-----------------------------------------------
    ! For test of the interpolation
    !-----------------------------------------------
    !xstart = -7.5
    !xstop = 6.0
    !dx = (xstop-xstart)/(9999)
    !do i=1, n
    !    xtest(i) = xstart + (i-1)*dx
    !end do
    

    ! write to file
    open(50, file = "data.dat")
    do i=1,n
       write(50,'(8(E17.8E3))') x_rec(i), X_e(i), get_tau(x_rec(i)), get_dtau(x_rec(i)), get_ddtau(x_rec(i)), get_g(x_rec(i)), &
       get_dg(x_rec(i)), get_ddg(x_rec(i))
    end do
    close(50)

  end subroutine initialize_rec_mod

  subroutine dXe_dx(x, X_e, deriv)
    !-----------------------------------------------
    ! Calculating the derivative of Xe
    !-----------------------------------------------
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x                                                               ! x-value for evaluation
    real(dp), dimension(:), intent(in)  :: X_e                                                             ! Xe-value for evaluation
    real(dp), dimension(:), intent(out) :: deriv                                                           ! Return value
    real(dp) :: lambda_21, lambda_a, Cr, beta, alpha2, n1s, beta2, phi2, H, n_b, T_b    
    T_b = T_0/exp(x)                                                                                       ! Baryon temperature
    n_b = Omega_b0*rho_c0/(m_H*exp(3.d0*x))                                                                ! Baryon number density
    H = get_H(x)                                                                                           ! Hubble constant
    n1s = n_b*(1.d0 - X_e(1))                                                                              ! Number of bound electrons
    phi2 = 0.448d0*log(epsilon_0/(k_b*T_b))                                                                
    lambda_a = H*(3.d0*epsilon_0)**3.d0/((8.d0*pi)**2.d0*n1s*(c*hbar)**3.d0)   ! [s-1]
    
    alpha2 = 64.d0*pi / sqrt(27.d0*pi) * (alpha/m_e)**2.d0 * sqrt(epsilon_0/(k_b*T_b)) * phi2 * hbar**2/c
    beta = alpha2 * (m_e*k_b*T_b / (2.d0*pi*hbar**2.d0))**1.5d0 * exp(-epsilon_0/(k_b*T_b))

    ! combined the betas to avoid infinities    
    beta2 = alpha2 * (m_e*k_b*T_b / (2.d0*pi*hbar**2.d0))**1.5d0 * exp((-epsilon_0)/(4.d0*k_b*T_b))

    Cr = (lambda_21 + lambda_a)/(lambda_21 + lambda_a + beta2) 
    deriv = Cr/H*(beta*(1.d0 - X_e) - n_b*alpha2*X_e**2.d0)
  end subroutine dXe_dx


  function get_n_e(x)
    !-----------------------------------------------
    ! Routine for computing n_e at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    get_n_e = exp(splint(x_rec, n_e, n_e2, x))                                                  
  end function get_n_e

  function get_tau(x)
    !-----------------------------------------------
    ! Routine for computing tau at arbitrary x, using precomputed information
    !-----------------------------------------------
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: get_tau
      get_tau  = exp(splint(x_rec, log_tau, tau2, x))
  end function get_tau

  subroutine dtau_dx(x,tau, deriv)
    !-----------------------------------------------
    ! Calculating the derivative of tau
    !-----------------------------------------------
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: deriv

    deriv = -get_n_e(x)*sigma_T*exp(x)*c/get_H_p(x)
  end subroutine dtau_dx

  function get_dtau(x)
    !-----------------------------------------------
    ! Routine for computing dtau at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    get_dtau = splint_deriv(x_rec, log_tau, tau2, x)
    get_dtau = get_dtau*get_tau(x) 
 end function get_dtau

  function get_ddtau(x)
    !-----------------------------------------------
    ! Routine for computing ddtau at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    get_ddtau = splint(x_rec, tau2, tau22, x)
    get_ddtau = get_ddtau*get_tau(x) + get_dtau(x)**2/get_tau(x)  
  end function get_ddtau

  function get_g(x)
    !-----------------------------------------------
    ! Routine for computing g_tilde at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x)
  end function get_g

  function get_dg(x)
    !-----------------------------------------------
    ! Routine for computing dg_tilde at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec, g, g2, x)
  end function get_dg

  function get_ddg(x)
    !-----------------------------------------------
    ! Routine for computing ddg_tilde at arbitrary x, using precomputed information
    !-----------------------------------------------
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec, g2, g22, x)
  end function get_ddg

end module rec_mod
