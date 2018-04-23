  subroutine integrate_perturbation_eqns
    implicit none
    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2, a, x, H, d_tau, dx

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, deriv, x_tc_rec

    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(deriv(npar))
    allocate(y_tight_coupling(7))

    allocate(x_tc_rec(n_tc_rec))
    x_tc_rec(1) = x_init 

    ! Propagate each k-mode independently

    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5

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
       !write(*,*) x_tc
       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
        
       dx = (x_start_rec-x_init)/n_tc_rec
       do i=2, n_tc_rec

           x_tc_rec(i) = x_tc_rec(i-1) + dx
           x = x_tc_rec(i)
           a = exp(x)
           d_tau = get_dtau(x)
           H = get_H_p(x)
           if (x < x_tc) then
               write(*,*) "1"
               call odeint(y_tight_coupling, x_tc_rec(i-1), x_tc_rec(i), eps, dx, hmin, dydx, bsstep, output)
               write(*,*) "2"
               delta(i, k)    = y_tight_coupling(1) 
               delta_b(i, k)  = y_tight_coupling(2)
               v(i, k)        = y_tight_coupling(3)
               v_b(i, k)      = y_tight_coupling(4)
               Phi(i, k)      = y_tight_coupling(5)
               Theta(i, 0, k) = y_tight_coupling(6)
               Theta(i, 1, k) = y_tight_coupling(7)
               Theta(i, 2, k) = -20.d0*c*k_current/(45.d0*H*d_tau)*Theta(1, i, k)               
               do l = 3, lmax_int
                   Theta(i,l,k) = -l/(2.d0*l + 1.d0)*c*k_current/(H*d_tau)*Theta(i, l-1, k)
               end do    
               Psi(i, k) = -Phi(i, k) - 12.d0*H_0**2/(c*k_current*a)**2*(Omega_r0*Theta(i, 2, k))         
               call dydx(x, y_tight_coupling, deriv)
               dPhi(i, k)       = deriv(5)
               dv_b(i, k)       = deriv(4)
               dTheta(i, 0, k)  = deriv(6)
               dTheta(i, 1, k)  = deriv(7)
               dTheta(i, 2, k)  = l*c*k_current/(2.d0*l + 1.d0)*Theta(i, 1, k) - (l + 1.d0)*c*k_current/((2.d0*l + 1.d0)*H)&
               *Theta(i, 3, k) + d_tau*(Theta(i, 2, k) - 0.1d0)
            else

               ! Task: Set up variables for integration from the end of tight coupling 
               ! until today
               y(1:7) = y_tight_coupling
               y(8)   = Theta(i, 2, k)
               do l = 3, lmax_int
                  y(6+l) = -l*c/(2.d0*l + 1.d0)*ks(k)/(H*d_tau)*Theta(i, l-1, k)
               end do

            end if 
      
    end do
    end do

