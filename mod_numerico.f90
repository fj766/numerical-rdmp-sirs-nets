module mod_numerico
   use geraRede
   use types
   implicit none
   	
   contains

        subroutine rk4_grafo(this, t, dt, m1, n1, gam, lam, rho, Iji, Si, Ii, f_Iji, f_Si, f_Ii)
		use types
		use geraRede
		
		class(grafo) :: this
         	real(dp) :: t, dt
        	integer :: m1, n1
        	real(dp) :: gam, lam, rho
        	real(dp) :: Si(n1), Ii(n1), Iji(m1)
        	real(dp) :: Sip(n1), Iip(n1), Ijip(m1)
        	real(dp) :: k1_Iji(m1), k2_Iji(m1), k3_Iji(m1), k4_Iji(m1)
        	
        	real(dp) :: k1_Si(n1), k2_Si(n1), k3_Si(n1), k4_Si(n1)
        	real(dp) :: k1_Ii(n1), k2_Ii(n1), k3_Ii(n1), k4_Ii(n1)
		real(dp) :: argIji(m1), argIi(n1), argSi(n1)  
 		
 		interface
 			function f_Iji(this, t, Iji, Si, m1, n1, rho, lam)
				use types
 				use geraRede
 				
 				class(grafo) :: this
 				real(dp), intent(in) :: t, rho, lam
 				integer :: m1, n1
 				real(dp) :: Iji(m1), Si(n1)
 				real(dp) :: f_Iji(m1)
 			end function
  		
  			function f_Si(this, t, Iji, Si, Ii, m1, n1, gam, lam)
 				use types
 				use geraRede
 				
 				class(grafo) :: this
 				real(dp), intent(in) :: t, gam, lam
 				integer :: m1, n1
 				real(dp) :: Iji(m1), Si(n1), Ii(n1)
 				real(dp) :: f_Si(n1)
 			end function			
  		
  			function f_Ii(this, t, Iji, Si, Ii, m1, n1, rho, lam)
 				use types
 				use geraRede
 				
 				class(grafo) :: this
 				real(dp), intent(in) :: t, rho, lam
 				integer :: m1, n1
 				real(dp) :: Iji(m1), Si(n1), Ii(n1) 
				real(dp) :: f_Ii(n1)
 			end function			
  		end interface	 
  		
  	!####################################################################################
		k1_Iji = dt * f_Iji(this, t, Iji, Si, m1, n1, rho, lam)
		
		k1_Si = dt * f_Si(this, t, Iji, Si, Ii, m1, n1, gam, lam)
		
		k1_Ii = dt * f_Ii(this, t, Iji, Si, Ii, m1, n1, rho, lam)

		!###############################
		argIji = Iji + 1d0/2d0 * k1_Iji
		argSi = Si + 1d0/2d0 * k1_Si
		argIi = Ii + 1d0/2d0 * k1_Ii
		!###############################

		k2_Iji = dt * f_Iji(this, t + (dt/2d0), argIji, argSi, m1, n1, rho, lam)
		
		k2_Si = dt * f_Si(this, t + (dt/2d0), argIji , argSi, argIi, m1, n1, gam, lam)		
		
		k2_Ii = dt * f_Ii(this, t + (dt/2d0), argIji , argSi, argIi, m1, n1, rho, lam)



		!###############################
		argIji = Iji + 1d0/2d0 * k2_Iji
		argSi = Si + 1d0/2d0 * k2_Si
		argIi = Ii + 1d0/2d0 * k2_Ii
		!###############################
		
		k3_Iji = dt * f_Iji(this, t + (dt/2d0), argIji, argSi, m1, n1, rho, lam)

		k3_Si = dt * f_Si(this, t + (dt/2d0), argIji, argSi, argIi, m1, n1, gam, lam)

		k3_Ii = dt * f_Ii(this, t + (dt/2d0), argIji, argSi, argIi, m1, n1, rho, lam)


		!###############################
		argIji = Iji + k3_Iji
		argSi = Si + k3_Si
		argIi = Ii + k3_Ii
		!###############################

		k4_Iji = dt * f_Iji(this, t + dt, argIji, argSi, m1, n1, rho, lam)

		k4_Si = dt * f_Si(this, t + dt, argIji, argSi, argIi, m1, n1, gam, lam)

		k4_Ii = dt * f_Ii(this, t + dt, argIji, argSi, argIi, m1, n1, rho, lam)


		
		Iji = Iji + 1d0/6d0 * (k1_Iji + 2d0 * k2_Iji + 2d0 * k3_Iji + k4_Iji)  		
				
		Si = Si + 1d0/6d0 * (k1_Si + 2d0 * k2_Si + 2d0 * k3_Si + k4_Si)  		
  		
		
		Ii = Ii + 1d0/6d0 * (k1_Ii + 2d0 * k2_Ii + 2d0 * k3_Ii + k4_Ii)  		
  		!####################################################################################
  		! F_Ii
  	  	        		    
        end subroutine		

end module
