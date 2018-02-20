program harmonic_oscillator
  use healpix_types
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b) :: i, n
  real(dp)     :: dt, x_0, v_0, eps
  real(dp), dimension(2) :: y
  real(dp), allocatable, dimension(:) :: t, v, x, x2, v2

  n = 50                     ! number of points to output functions
  dt = 10.d0 * pi / (n - 1)  ! time step in ode_solver
  x_0 = 0.d0                 ! initial position
  v_0 = 1.d0                 ! initial velocity
  eps = 1.d-8                ! error tolerance in ode_solver function
  
  allocate(t(n))
  allocate(x(n))
  allocate(v(n))
  do i=1,n
     t(i) = (i-1) * dt
  end do
  ! solve ode-system:
  x(1) = x_0
  v(1) = v_0
  y(1) = x(1)
  y(2) = v(1)
  do i=1,n-1
     ! odeint takes arugments ()y, x_start, x_stop, tolerance,  high_t, low_t, function, ode_func, output_name
     call odeint(y, t(i), t(i+1), eps, dt / 10.d0, dt / 10000.d0, derivs, bsstep, output)
     x(i+1) = y(1)
     v(i+1) = y(2)
  end do


  ! output to file
  open(54, file='lowres_data.dat')
  do i=1,n
     write(54, '(3(E17.8))') t(i), x(i), v(i)
  end do
  close(54)
  write(*,*) "Done writing lowres data to file"

  ! spline lowres data
  allocate(x2(n))
  allocate(v2(n))
  ! spline takes arguments (x, y, dy/dx|_1, dy/dx|_n, d^2y/dx^2),
  ! where the last argument is the output of the subroutine.
  ! Setting the two first derivatives to 1e30 corresponds to 
  ! choosing the "natural spline". 
  call spline(t, x, 1d30, 1d30, x2)
  call spline(t, v, 1d30, 1d30, v2)

  ! interpolate to high resolution grid: 
contains

  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)  :: k_over_m
    k_over_m = 1.d0
    dydx(1) = y(2)
    dydx(2) = -k_over_m * y(1)
  end subroutine derivs

  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output
  
end program harmonic_oscillator
