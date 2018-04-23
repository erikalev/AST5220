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
  integer(i4b), parameter          :: n_tc_rec = 1000
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
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
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
    integer(i4b) :: n_tot   
    real(dp)     :: dtau, H, dk, dx  
    real(dp), allocatable, dimension(:):: y1, x_tc_rec
    n_tot = n_tc_rec + n_t    

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do i=1, n_k
        dk = (k_max-k_min)*((i-1)/n_k)**2.d0
        ks(i) = k_min + dk
    end do

    ! Allocate arrays for perturbation quantities
    allocate(Theta(n_tot, 0:lmax_int, n_k))
    allocate(delta(n_tot, n_k))
    allocate(delta_b(n_tot, n_k))
    allocate(v(n_tot, n_k))
    allocate(v_b(n_tot, n_k))
    allocate(Phi(n_tot, n_k))
    allocate(Psi(n_tot, n_k))
    allocate(dPhi(n_tot, n_k))
    allocate(dPsi(n_tot, n_k))
    allocate(dv_b(n_tot, n_k))
    allocate(dTheta(n_tot, 0:lmax_int, n_k))

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(1,:)     = 1.d0
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1, :)
    Theta(1,0,:) = 0.5d0*Phi(1, :)
    dPhi(1, :) = 0.d0
    dPsi(1, :) = 0.d0
    dv_b(1, :) = 0.d0
    dTheta(1, 0, :) = 0.d0
       
    dtau = get_dtau(log(a_init))
    H    = get_H_p(log(a_init))
    do j = 1, n_k
       v(1,j)       =  c*ks(j)/(2.d0*H)*Phi(1, j)
       v_b(1,j)     =  v(1, j)
       Theta(1,1,j) = -c*ks(j)/(6.d0*H)*Phi(1, j)
       Theta(1,2,j) = -20.d0*c*ks(j)/(45.d0*H*dtau)*Theta(1, 0, j)
       do l = 3, lmax_int
           Theta(1,l,j) = -l*c/(2.d0*l + 1.d0)*ks(j)/(H*dtau)*Theta(1, l-1, j)
       end do    
    Psi(1, j) = -Phi(1, j) - 12.d0*H_0**2.d0/(c*ks(j)*a_init)**2.d0*Omega_r0*Theta(1, 2, j)
    end do


    dx = (x_start_rec-x_init)/(n_tc_rec-1)
    allocate(x_tc_rec(n_tc_rec))
    x_tc_rec(1) = x_init 
    x_tc_rec(2) = x_init + dx

    allocate(y1(7))
    y1(1) = delta(1, 1)
    y1(2) = delta_b(1, 1)
    y1(3) = v(1, 1)
    y1(4) = v_b(1, 1)
    y1(5) = Phi(1, 1)
    y1(6) = Theta(1, 0, 1)
    y1(7) = Theta(1, 1, 1)
    call odeint(y1, x_tc_rec(1), x_tc_rec(2), eps, dx, hmin, dydx, bsstep, output)
    do i=1, 7
        write(*,*) y1(i)
    end do
  end subroutine initialize_perturbation_eqns



  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)

  subroutine dydx(x, y, deriv)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: deriv
    real(dp)                            :: delta, delta_b, v, v_b, Phi, Theta0, Theta1, Theta2, Psi, dtau, Hp, R, a
    integer(i4b) :: i

    delta = y(1)
    delta_b = y(2)
    v = y(3)
    v_b = y(4)
    Phi = y(5)
    Theta0 = y(6)
    Theta1 = y(7)
    dtau = get_dtau(x)
    a = exp(x)
    Hp = get_H_p(x)
    Theta2 = -20.d0*c*k_current/(45.d0*Hp*dtau)*Theta1
    Psi = -Phi - 12.d0*H_0**2.d0/(c*k_current*a)**2.d0*Omega_r0*Theta2

    R = 4.d0*Omega_r0/(3.d0*Omega_b0*a)


    deriv(3) = -v - c*k_current/Hp*Psi      ! d_v

    deriv(5) = Psi - (c*k_current)**2.d0/(3.d0*Hp**2.d0)*Phi + H_0**2.d0/(2.d0*Hp**2.d0)*&
    (Omega_m0/a*delta + Omega_b0/a*delta_b + 4.d0*Omega_r0/a**2.d0*Theta0)  !d_phi

    deriv(2) = c*k_current/Hp*v_b - 3.d0*deriv(5) ! d_deltab

    deriv(1) = c*k_current/Hp*v - 3.d0*deriv(5) ! d_delta

    deriv(4) = -v_b -c*k_current/Hp*Psi + dtau*R*(v_b + 3.d0*Theta1) ! d_vb

    deriv(6) = -c*k_current/Hp*Theta1 - deriv(5) !d_Theta0

    deriv(7) = c*k_current/(3.d0*Hp)*Theta0 - 2.d0*c*k_current/(3.d0*Hp)*Theta2 + c*k_current/(3.d0*Hp)*Psi + dtau*&
    (Theta1 + 1.d0/3.d0*v_b)  !d_Theta1

    !do i=1, 7
    !    write(*,*) deriv(i)
    !end do
    !write(*,*) "-----------------"
  end subroutine dydx

  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    real(dp)              :: x, d_tau, H, dx
    integer(i4b)          :: n, i
    dx = (x_start_rec - x_init)/(n_tc_rec-1)

    do i=1, n_tc_rec
        x = log(a_init) + (i-1)*dx
        d_tau = get_dtau(x)
        H = get_H_p(x)
        if ((abs(d_tau) > 10.d0) .and. (abs(c*k/(H*d_tau)) < 0.1d0) .and. (x < x_start_rec)) then
            get_tight_coupling_time = x
        end if
    end do
  end function get_tight_coupling_time

  
end module evolution_mod
