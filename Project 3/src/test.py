           do j = i_tc, n_t
              !Integrate equations from tight coupling to today
              !write(*,*) 'running odeint with j =', j
              call odeint(y, x_t(j-1) ,x_t(j), eps, h1, hmin, dydx_no_tc, bsstep, output)
    
              !Store variables at time step j in global variables
              delta(j,k)   = y(1)
              delta_b(j,k) = y(2)
              v(j,k)       = y(3)
              v_b(j,k)     = y(4)
              Phi(j,k)     = y(5)
              
              do l = 0, lmax_int
                 Theta(j,l,k) = y(6+l)
              end do
              Psi(j,k)     =  - Phi(j,k) - (12.d0*H_0**2.d0)/(ck_current*a_t(j))**2.d0*Omega_r0*Theta(j,2,k)
    
              !Store derivatives that are required for C_l estimation
              call dydx_no_tc(x_t(j),y,deriv)
              dv_b(j,k)     = deriv(4)
              dPhi(j,k)     = deriv(5)
              !if(k==40) then
              !    write(*,*) 'dPhi(',j,k,') =', dPhi(j,k)
              !endif            
    
              do l=0,lmax_int
                  dTheta(j,l,k) = deriv(6+l)
              end do
    
              dPsi(j,k)     = -dPhi(j,k) - 12.d0*H_0**2.d0/(ck_current*a_t(j))**2.d0*&
                               Omega_r0*(-2.d0*Theta(j,2,k)+dTheta(j,2,k)) 
end do
