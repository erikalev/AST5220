module time_mod
  use healpix_types
  use ode_solver
  use params
  use spline_1D_mod
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta, a_eta       ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2, eta_intp          ! Eta and eta'' at each grid point
  real(dp),    allocatable, dimension(:) :: Omega_m, Omega_b, Omega_r, Omega_lambda, Omega_nu, H
contains

  subroutine initialize_time_mod
    implicit none
    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, eta_init, eps, hmin, yp1
    real(dp)     :: ypn, rho_c
    real(dp), dimension(1)     :: y
    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    
    do i=1,n1                              ! during recombination
       x_t(i) = x_start_rec - (i-1)*(x_start_rec - x_end_rec)/(n1-1)
    end do

    do i=1, n2                             ! after recombination
       x_t(n1+i) = x_end_rec - i*(x_end_rec - x_0)/(n2-1)
    end do

    a_t = exp(x_t)

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))


	! Finding the x_eta-values
    dx = (x_eta2 - x_eta1)/(n_eta - 1)
    x_eta(1) = x_eta1
    do i=2, n_eta
       x_eta(i) = x_eta(i-1) + dx
    end do


    ! Integrating eta-values
    eta(1) = c*a_init/(H_0*sqrt(Omega_r0))                                                  ! Initial eta-value
    dx = abs(1.d-2*(x_eta(2)-x_eta(1)))                                                     ! Step length
    eps = 1.d-10                                                  
    hmin = 0.d0
    do i=2, n_eta
       eta(i) = eta(i-1)										                            ! old value to be updated in next line
       call odeint(eta(i:i), x_eta(i-1), x_eta(i), eps, dx, hmin, derivs, bsstep, output)   ! updating eta(i)-value
    end do
    
    allocate(Omega_r(n_eta))
    allocate(Omega_m(n_eta))
    allocate(Omega_b(n_eta))
    allocate(Omega_lambda(n_eta))
    allocate(Omega_nu(n_eta))
    allocate(H(n_eta))

    ! Calculating Omega-values and H
    do i=1, n_eta
       H(i) = get_H(x_eta(i))
       Omega_r(i) = Omega_r0*exp(x_eta(i))**(-4.d0)*(H_0/H(i))**2
       Omega_nu(i) = Omega_nu0*exp(x_eta(i))**(-4.d0)*(H_0/H(i))**2
       Omega_m(i) = Omega_m0*exp(x_eta(i))**(-3.d0)*(H_0/H(i))**2
       Omega_b(i) = Omega_b0*exp(x_eta(i))**(-3.d0)*(H_0/H(i))**2
       Omega_lambda(i) = 1.d0 - Omega_r(i) - Omega_m(i) - Omega_b(i) - Omega_nu(i)
    end do
    
    yp1 = 1.d30
    ypn = 1.d30
    call spline(x_eta, eta, yp1, ypn, eta2)         ! finding the derivaties of the curve


	! finding the interpolated eta-values
    allocate(eta_intp(5000))
    dx = (x_eta2 - x_eta1)/(n_eta - 1)/5.0			! interpolation step size
	eta_intp(1) = eta(1)						
    do i=2, 5000
       eta_intp(i) = get_eta(x_eta(1) + (i - 1)*dx)
    end do
 

    ! write values to file                                                                                                                                                                                         
    open(54, file = "data.dat")
    do i=1,n_eta
       write(54,'(8(E17.8))') eta(i), x_eta(i), Omega_m(i), Omega_b(i), Omega_r(i), Omega_lambda(i), Omega_nu(i), H(i)
    end do
    close(54)

    open(53, file = "eta_intp.dat")
    do i=1, 5000
       write(53, "(1(E17.8))") eta_intp(i)
    end do
    close(53)
    
  end subroutine initialize_time_mod


  subroutine derivs(x, eta, d_eta_dx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: d_eta_dx
    d_eta_dx = c/get_H_p(x)
  end subroutine derivs


  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output

  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    get_H = H_0*sqrt((Omega_b0 + Omega_m0)*exp(-3*x) + (Omega_r0 + Omega_nu0)*exp(-4*x) + Omega_lambda0)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    get_H_p = exp(x)*get_H(x)

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in)  :: x
    real(dp)              :: get_dH_p
    get_dH_p = exp(x)*(get_H(x) + H_0**2*(-3*(Omega_m0 + Omega_b0)*exp(-3*x) - 4*Omega_r0*exp(-4*x))/(2*get_H(x)))
  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta = splint(x_eta, eta, eta2, x_in)
  end function get_eta

end module time_mod
