module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),     parameter, private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:)   :: ddelta
  real(dp), allocatable, dimension(:,:)   :: ddelta_b
  real(dp), allocatable, dimension(:,:)   :: dv
  real(dp), allocatable, dimension(:,:,:) :: dTheta
  integer(i4b), allocatable, dimension(:) :: h_plot

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current, ck, ckH, ck_current
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, j
    real(dp)     :: dtau, H, dk    

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))

    do i=1, n_k
        ks(i) = k_min + (k_max-k_min)*((i-1.d0)/(n_k-1.d0))**2.d0
    end do

    ! Allocate arrays for perturbation quantities
    allocate(Theta(n_t, 0:lmax_int, n_k))
    allocate(delta(n_t, n_k))
    allocate(delta_b(n_t, n_k))
    allocate(v(n_t, n_k))
    allocate(v_b(n_t, n_k))
    allocate(Phi(n_t, n_k))
    allocate(Psi(n_t, n_k))
    allocate(dPhi(n_t, n_k))
    allocate(dPsi(n_t, n_k))
    allocate(dv_b(n_t, n_k))
    allocate(dTheta(n_t, 0:lmax_int, n_k))
    allocate(ddelta(n_t, n_k))
    allocate(ddelta_b(n_t, n_k))
    allocate(dv(n_t, n_k))
    allocate(h_plot(6))

    !------------------------------------------------------------
    ! Initial conditions for the Boltzmann and Einstein equations
    !------------------------------------------------------------

    Phi(1,:)        = 1.d0
    delta(1,:)      = 1.5d0*Phi(1,:)
    delta_b(1,:)    = delta(1, :)
    Theta(1,0,:)    = 0.5d0*Phi(1, :)
    dPhi(1, :)      = 0.d0
    dPsi(1, :)      = 0.d0
    dv_b(1, :)      = 0.d0
    dTheta(1, :, :) = 0.d0
    dtau = get_dtau(x_t(1))
    H    = get_H_p(x_t(1))
    v(1,:)       =  c*ks/(2.d0*H)*Phi(1, :)
    v_b(1,:)     =  v(1, :)
    Theta(1,1,:) = -c*ks/(6.d0*H)*Phi(1, :)
    Theta(1,2,:) = -20.d0*c*ks/(45.d0*H*dtau)*Theta(1, 1, :)
    do l = 3, lmax_int
       Theta(1,l,:) = -l*c/(2.d0*l + 1.d0)*ks/(H*dtau)*Theta(1, l-1, :)
    end do    
    Psi(1, :) = -Phi(1, :) - 12.d0*H_0**2.d0/(c*ks*exp(x_t(1)))**2.d0*Omega_r0*Theta(1, 2, :)
  end subroutine initialize_perturbation_eqns

  
  subroutine integrate_perturbation_eqns
    implicit none
    integer(i4b) :: i, j, k, l, i_tc, o,p
    real(dp)     :: x1, x2
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2, a, x, H, d_tau, dx, x_test

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, deriv, y_test
    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5

    allocate(y(npar))
    allocate(deriv(npar))
    allocate(y_tight_coupling(7))
    allocate(y_test(npar))

    ! h-values for plotting
    h_plot(1) = 1
    h_plot(2) = 10
    h_plot(3) = 30
    h_plot(4) = 50
    h_plot(5) = 80
    h_plot(6) = 100
 


    ! Propagate each k-mode independently
    do k = 1, n_k
       write(*,*) "k", k

       k_current = ks(k)  ! Store k_current as a global module variable       h1        = 1.d-5
       ck_current = c*k_current

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(1,k)
       y_tight_coupling(2) = delta_b(1,k)
       y_tight_coupling(3) = v(1,k)
       y_tight_coupling(4) = v_b(1,k)
       y_tight_coupling(5) = Phi(1,k)
       y_tight_coupling(6) = Theta(1,0,k)
       y_tight_coupling(7) = Theta(1,1,k)

       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       i = 2
       do while (x_t(i) < x_tc)
           x                = x_t(i)
           a                = a_t(i)
           d_tau            = get_dtau(x)
           H                 = get_H_p(x)
           ckH              = ck_current/H
           call odeint(y_tight_coupling, x_t(i-1), x_t(i), eps, h1, hmin, dydx_tc, bsstep, output)
           
           delta(i, k)      = y_tight_coupling(1) 
           delta_b(i, k)    = y_tight_coupling(2)
           v(i, k)          = y_tight_coupling(3)
           v_b(i, k)        = y_tight_coupling(4)
           Phi(i, k)        = y_tight_coupling(5)
           Theta(i, 0, k)   = y_tight_coupling(6)
           Theta(i, 1, k)   = y_tight_coupling(7)              
           Theta(i, 2, k)   = -20.d0*ckH/(45.d0*d_tau)*Theta(i, 1, k)              

           do l = 3, lmax_int
               Theta(i,l,k) = -l/(2.d0*l + 1.d0)*ckH/d_tau*Theta(i, l-1, k)
           end do    

           Psi(i, k)        = -Phi(i, k) - 12.d0*H_0*H_0/((ck_current*a)*(ck_current*a))*(Omega_r0*Theta(i, 2, k))         

           call dydx_tc(x_t(i), y_tight_coupling, deriv)
           ddelta(i, k)     = deriv(1)
           ddelta_b(i, k)   = deriv(2)
           dv(i, k)         = deriv(3)
           dv_b(i, k)       = deriv(4)
           dPhi(i, k)       = deriv(5)
           dTheta(i, 0, k)  = deriv(6)
           dTheta(i, 1, k)  = deriv(7)

           dPsi(i, k)       = -dPhi(i, k) - 12.d0*H_0**2.d0/(ck_current*a)**2.d0*Omega_r0*(-2.d0*Theta(i, 2, k) + dTheta(i, 2, k))
           i                = i + 1
        end do
        i_tc = i

        ! Task: Set up variables for integration from the end of tight coupling 
        ! until today
        y(1:7) = y_tight_coupling
        y(8)   = Theta(i_tc-1, 2, k) 

        do l = 3, lmax_int
           y(6+l) = Theta(i_tc-1, l, k)
        end do
        x_test = -7.39774538677127
        if (ANY( h_plot==k ) ) then
            write(*,*) "x_tc", x_t(i_tc)
            write(*,*) "--------------------"
        end if
        do i=i_tc, n_t
           call odeint(y, x_t(i-1), x_t(i), eps, h1, hmin, dydx_no_tc, bsstep, output)

           delta(i, k)      = y(1) 
           delta_b(i, k)    = y(2)
           v(i, k)          = y(3)
           v_b(i, k)        = y(4)
           Phi(i, k)        = y(5)

           do l = 0, lmax_int
               Theta(i,l,k) = y(6+l)
           end do           
           Psi(i, k)        = -Phi(i, k) - 12.d0*H_0**2/(ck_current*a_t(i))**2*(Omega_r0*Theta(i, 2, k)) 

           call dydx_no_tc(x, y, deriv)
           dv_b(i, k)       = deriv(4)
           dPhi(i, k)       = deriv(5)
           do l = 0, lmax_int
               dTheta(i,l,k) = deriv(6+l)
           end do
        end do  
    end do
    
    !---------------------------------------------------------
    ! Printing shit to file
    !---------------------------------------------------------

    open(49, file = "delta.dat")  
        write(49,'(6(E17.8E3))') delta(:, 1), delta(:, 10), delta(:, 30), delta(:, 50), delta(:, 80), delta(:, 100) 
    close(49)
    
    open(50, file = "data.dat")  
    do i=1, n_t
        do k=1, 6
            write(50,'(22(E25.8E3))') Phi(i, h_plot(k)), dPhi(i, h_plot(k)), Psi(i, h_plot(k)), dPsi(i, h_plot(k)), &
            delta(i, h_plot(k)), ddelta(i, h_plot(k)), delta_b(i, h_plot(k)), ddelta_b(i, h_plot(k)), v(i, h_plot(k)), &
            dv(i, h_plot(k)), v_b(i, h_plot(k)), dv_b(i, h_plot(k)), Theta(i, 0, h_plot(k)), dTheta(i, 0, h_plot(k)), &
            Theta(i, 1, h_plot(k)), dTheta(i, 1, h_plot(k)), Theta(i, 2, h_plot(k)), dTheta(i, 2, h_plot(k)), &
            Theta(i, 3, h_plot(k)), dTheta(i, 3, h_plot(k)), Theta(i, 4, h_plot(k)), dTheta(i, 4, h_plot(k))
        end do
    end do
    close(50)
    
    open(51, file = "x_data.dat")  
    do i=1, n_t
        write(51,'(1(E17.8E3))') x_t(i)
    end do
    close(51)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)

  subroutine dydx_tc(x, y, deriv)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: deriv
    real(dp)                            :: delta, delta_b, v, v_b, Phi, Theta0, Theta1, Theta2, Psi, dtau, ddtau, Hp, R, a, q, dHp
    real(dp)                            :: dPhi, dTheta0, dTheta1, ddelta, ddelta_b, dv, dv_b
    integer(i4b) :: i

    delta       = y(1)
    delta_b     = y(2)
    v           = y(3)
    v_b         = y(4)
    Phi         = y(5)
    Theta0      = y(6)
    Theta1      = y(7)

    dtau        = get_dtau(x)    
    ddtau       = get_ddtau(x)
    a           = exp(x)
    Hp          = get_H_p(x)
    dHp         = get_dH_p(x)
    ckH         = ck_current/Hp

    Theta2      = -20.d0*ckH/(45.d0*dtau)*Theta1

    Psi         = -Phi - 12.d0*H_0*H_0/(ck_current*a)**2.d0*Omega_r0*Theta2

    R           = 4.d0*Omega_r0/(3.d0*Omega_b0*a)

    dv          = -v - ckH*Psi

    dPhi        = Psi - ckH*ckH/3.d0*Phi + H_0*H_0/(2.d0*Hp*Hp)*(Omega_m0/a*delta + Omega_b0/a*delta_b + &
                  4.d0*Omega_r0/(a*a)*Theta0)

    dTheta0     = -ckH*Theta1 - dPhi


    q = (-((1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau)*(3.d0*Theta1 + v_b) - ckH*Psi + (1.d0 - dHp/Hp)*ckH*(-Theta0 + 2.d0*Theta2)& 
        - ckH*dTheta0)/((1.d0+R)*dtau + dHp/Hp - 1.d0)
    
    ddelta_b    = ckH*v_b - 3.d0*dPhi ! d_deltab

    ddelta      = ckH*v - 3.d0*dPhi ! d_delta

    dv_b        = 1.d0/(1.d0 + R)*(-v_b -ckH*Psi + R*(q + ckH*(-Theta0 + 2.d0*Theta2) - ckH*Psi)) ! d_vb

    dTheta1     = 1.d0/3.d0*(q - dv_b)
    
    deriv(1)    = ddelta
    deriv(2)    = ddelta_b
    deriv(3)    = dv
    deriv(4)    = dv_b
    deriv(5)    = dPhi
    deriv(6)    = dTheta0
    deriv(7)    = dTheta1

  end subroutine dydx_tc

  subroutine dydx_no_tc(x,y, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      integer(i4b) :: l
      real(dp) :: delta,delta_b,ddelta,ddelta_b,v,v_b,dv,dv_b,Phi,Theta0,Theta1,Theta2,Theta3,Psi,dPhi,dTheta0,dTheta1,dTheta2
      real(dp) :: a,H_p,dtau,R

      delta     = y(1)
      delta_b   = y(2)
      v         = y(3)
      v_b       = y(4)
      Phi       = y(5)
      Theta0    = y(6)
      Theta1    = y(7)
      Theta2    = y(8)
      Theta3    = y(9)

      a         = exp(x)
      H_p       = get_H_p(x)
      ckH       = ck_current/H_p
      dtau      = get_dtau(x)

      Psi       = -Phi - 12.d0*H_0**2/(ck_current*a)**2*Omega_r0*Theta2

      R         = 4.d0*Omega_r0/(3.d0*Omega_b0*a)

      dv        = -v - ckH*Psi

      dPhi      = Psi - ckH*ckH*Phi/3.d0 + (H_0*H_0)/(2.d0*H_p*H_p)*(Omega_m0/a*delta + Omega_b0/a*delta_b + &
                  4.d0*Omega_r0/(a*a)*Theta0)

      dTheta0   = -ckH*Theta1 -dPhi

      ddelta_b  = ckH*v_b - 3.d0*dPhi ! d_deltab
      
      ddelta    = ckH*v - 3.d0*dPhi
      
      dv_b      = -v_b - ckH*Psi + dtau*R*(3.d0*Theta1 + v_b) ! d_vb
 
      dTheta1   = ckH*Theta0/3.d0 -2.d0/3.d0*ckH*Theta2 + ckH/3.d0*Psi + dtau*(Theta1 + v_b/3.0)
    
      dTheta2   = 2.d0/5.d0*ckH*Theta1 - 3.d0/5.d0*ckH*Theta3 + dtau*(0.9d0*Theta2)

      dydx(1)   = ddelta
      dydx(2)   = ddelta_b
      dydx(3)   = dv
      dydx(4)   = dv_b
      dydx(5)   = dPhi
      dydx(6)   = dTheta0
      dydx(7)   = dTheta1
      dydx(8)   = dTheta2

      do l=3,lmax_int-1
          dydx(6+l) = l/(2.d0*l + 1.d0)*ckH*y(5 + l) - (l + 1.d0)/(2.d0*l + 1.d0)*ckH*y(7 + l) +dtau*y(6 + l)
      end do

      dydx(6+lmax_int) = ckH*y(6+lmax_int-1) -c*(lmax_int+1.d0)/(H_p*get_eta(x))*y(6+lmax_int) +dtau*y(6+lmax_int)

end subroutine dydx_no_tc

  function get_tight_coupling_time(k)
    implicit none
    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    real(dp)              :: x, d_tau, H
    integer(i4b)          :: i
    do i=1, n_t
        x       = x_t(i)
        d_tau   = get_dtau(x)
        H       = get_H_p(x)
        if ((abs(d_tau) > 10.d0) .and. (abs(c*k/(H*d_tau)) < 0.1d0) .and. (x < x_start_rec)) then
            get_tight_coupling_time = x
        end if
    end do
  end function get_tight_coupling_time

  
end module evolution_mod
