module hashimoto
   use geraRede
   use mod_rndgen
   use mod_tools
   use types
   
   implicit none
   real(dp), allocatable :: x(:), y(:)
   integer, allocatable :: link_back2(:)
   real(dp), parameter :: tol = 1.0d-9
   real(dp), allocatable :: fi(:)
   real(dp) :: Y4
! Calcula o autovalor principal da matriz B_hsh
! 
   contains

!=======================================================================
! Subrotinas
!=======================================================================
   
   subroutine aloca_listas_e_matrizHashimoto(this)
      class(grafo) :: this
      integer :: i1, i2, i3, j4
      integer (kind=8) ::  i12, i23

      !########################################################


			if(allocated(x)) deallocate(x)
			allocate(x(this%sumDeg))
			
			if(allocated(y)) deallocate(y)
			allocate(y(this%sumDeg))	   
   
            if(allocated(link_back2)) deallocate(link_back2)
            allocate(link_back2(size(y)))
            
            do i1 = 1, this%nodes
               do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
                  i2 = this%listAdj(i12)
lb_hs:            do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1
                     i3 = this%listAdj(i23)
                     if( i3 == i1 )then
                        link_back2(i12) = i23
                        exit lb_hs
                     endif
                  enddo lb_hs
               enddo
            enddo
            
            
   end subroutine                       
         !#######################################################
                                                                     
         !########################################################

   !====================================================================
   subroutine metodo_potencia_Hashimoto(this, lambda, mu1)
      class(grafo), intent(in) :: this
      real(dp), intent(out) :: lambda
      real(dp), intent(in) :: mu1
	  type(rndgen) :: geni
	  integer :: semente
	  real(dp) :: erro
	  real(dp) :: autovalor, autovalor_0
	  integer :: i1, i2, i3
	  integer :: j1, j2, j3, j4
	  integer(kind = 8) :: j12, j21, j23
      integer (kind=8) :: i12, i21, i23, i32
      real(dp) :: x_norm, y_norm
      real(dp), allocatable :: x_prov(:)      
	  integer :: contador, num_cont
	  real(dp) :: Y4_0
	  real(dp), allocatable :: rhoi_at(:), rhoi_int(:), rhoi_old(:) 
	  real(dp) :: norm_rhoi_at, norm_rhoi_int, norm_rhoi_old
	  real(dp), allocatable :: v1(:)
	  real(dp) :: soma
	  real(dp) :: max_x, max_y, max_x_prov
	  real(dp) :: soma_x, soma_y, soma_xprov
	  real(dp) :: sum_xi, sum_ki_xi, Mu_An
	  semente = 956792214

	  if(allocated(rhoi_at)) deallocate(rhoi_at)
	  allocate(rhoi_at(this%nodes))

	  if(allocated(rhoi_int)) deallocate(rhoi_int)
	  allocate(rhoi_int(this%nodes))

	  if(allocated(rhoi_old)) deallocate(rhoi_old)
	  allocate(rhoi_old(this%nodes))
	  
     !==================================================================	  
	  call geni%init(semente)
     !==================================================================     
      do j1 = 1, this%nodes
         do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
            j2 = this%listAdj(j12)
            j21 = link_back2(j12)
            y(j21) = this%deg(j2)
         enddo
      enddo
      !=================================================================
      y_norm = (sum(y**2.0d0))**0.5d0
      !=================================================================
      if(allocated(x_prov)) deallocate(x_prov)
      allocate(x_prov(size(y)))
      !=================================================================
      ! Apos normalizar o vetor x, a norma ||xp|| = ||x||_{oo}
      ! se torna = 1.0d0.
      !=================================================================
      j1 = 0
      !=================================================================
      num_cont = 2
      Y4_0 = 0.0d0
      autovalor_0 = 0.0d0
      contador = 1
!=======================================================================  
iter: do while(.True.)	
	     !Aqui comeca o algoritmo
         !==============================================================
         ! Y = A X
         !==============================================================
         y = y/y_norm
         x_prov = y
         !==============================================================
         do j3 = 1, num_cont
            x = y
            sit_i1: do i1 = 1, this%nodes
               viz_i2: do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) -1               
                  ! link u-->v
                  !-----------------------------------------------------
                  i21 = link_back2(i12)
                  !-----------------------------------------------------
                  i2 = this%listAdj(i12)
                  soma = 0.0d0
viz_i3:           do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1 ! link v --> w
                     !--------------------------------------------------
                     i3 = this%listAdj(i23)
                     !--------------------------------------------------
                     i32 = link_back2(i23)
                     !--------------------------------------------------
                     !if(i3 == i1) cycle viz_i3                     
                     !soma = soma + x(i23)
                     !--------------------------------------------------
                     if( i32 == i12) cycle viz_i3
                     soma = soma + x(i32)
                     !--------------------------------------------------
                  enddo viz_i3
                  !y(i12) = soma
                  y(i21) = soma                   
	           enddo viz_i2
	        enddo sit_i1      
		 enddo            

         !==============================================================         
         y_norm = (sum(y**2.0d0))**0.5d0
         x_norm = (sum(x**2.0d0))**0.5d0
         !==============================================================
         autovalor = (dot_product(x_prov, y))**(1.0d0/(1.0d0 * num_cont))
         !==============================================================
         ! Essa linha deve ser comentada caso tenha dado errado.
         !==============================================================
         max_x = sum(x)
         max_y = sum(y)
         !==============================================================
         do j1 = 1, this%nodes
            !===========================================================
            soma_x = 0.0d0; soma_y = 0.0d0 !; soma_xprov = 0.0d0
            do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
               !========================================================            
               j21 = link_back2(j12)
               !========================================================
               soma_y = soma_y + y(j21)/max_y
               soma_x = soma_x + x(j21)/max_x
               !========================================================
            enddo
            !===========================================================
            rhoi_at(j1) = soma_y/autovalor
            rhoi_int(j1) = soma_x/autovalor
            !===========================================================
         enddo
   	     !==============================================================
         norm_rhoi_at = (sum(rhoi_at**2.0d0))**0.5d0
   	     norm_rhoi_int = (sum(rhoi_int**2.0d0))**0.5d0
   	     !==============================================================
   	     rhoi_at = rhoi_at/norm_rhoi_at
   	     rhoi_int = rhoi_int/norm_rhoi_int
   	     !==============================================================
         Y4 = sum(((rhoi_at + rhoi_int)/2.0d0)**4.0d0)
         !==============================================================
         erro = max(abs(autovalor - autovalor_0), abs(Y4 - Y4_0))
         autovalor_0 = autovalor
         Y4_0 = Y4
         !==============================================================   
		 if(erro < tol)then
		    sum_xi = sum((rhoi_at + rhoi_int )/2.0d0)
		    sum_ki_xi = dot_product(1.0d0 * this%deg, (rhoi_at + rhoi_int )/2.0d0)
		    Mu_An = sum_ki_xi/sum_xi - 1
		    write(*,*) "Mu_An = ", Mu_An, "e ", "Mu_num = ", autovalor
		    !-----------------------------------------------------------
            if( this%nodes == 10**4 )then
                     
               call calcula_distancias_ao_hub(this)
               
               if(allocated(fi)) deallocate(fi)
               allocate(fi(this%nodes))
               
               fi = (rhoi_at + rhoi_int )/2.0d0
               deallocate(rhoi_at)
               deallocate(rhoi_int)

            endif
		    !-----------------------------------------------------------
		    lambda = mu1/autovalor
		    write(13,*) this%nodes, lambda
			write(*,*) 'Autovalor principal ', autovalor, ' encontrado com erro relativo ', erro, '. '
			
			write(*,*) " "
			write(*,*) 'Foram necessarios ', contador, ' passos.'
			write(*,*) " "
			deallocate(link_back2)
			exit iter
         endif
         
         !######################
         contador = contador + 1
      enddo iter
      close(13)
      close(14)
    
    end subroutine
!########################################################


end module

module mod_SIRS_rDMP_numerico
	use geraRede
	use mod_tools
	use types
	implicit none
	
	real(dp) :: t

	real(dp), allocatable :: Ii(:), Iji(:), Ri(:)
	real(dp), allocatable :: Ii_tild(:), Iji_tild(:), Ri_tild(:)
	real(dp), allocatable :: k1_Ii(:), k1_Iji(:), k1_Ri(:) 
	real(dp), allocatable :: k2_Ii(:), k2_Iji(:), k2_Ri(:)
	real(dp), allocatable :: k3_Ii(:), k3_Iji(:), k3_Ri(:)
	real(dp), allocatable :: k4_Ii(:), k4_Iji(:), k4_Ri(:)
	real(dp), allocatable :: v1(:)
	integer, allocatable :: link_back(:)
	
	contains
	
	
		!#################################################
		!	Modelo SIRS
		!#################################################		

    subroutine aloca_listas_SIRS_rDMP(this)
       class(grafo) :: this
       integer :: i1, i2, i3
       integer(kind=8) :: i12, i13, i21, i23, i31, i32
              
       if(allocated(k1_Ii))then
          deallocate(Ii)
          deallocate(Ri)
          deallocate(Iji)
          
          deallocate(v1)         
          
          deallocate(k1_Ii)
          deallocate(k2_Ii)
          deallocate(k3_Ii)
          deallocate(k4_Ii)
          
          deallocate(k1_Ri)
          deallocate(k2_Ri)
          deallocate(k3_Ri)
          deallocate(k4_Ri)

          deallocate(k1_Iji)
          deallocate(k2_Iji)
          deallocate(k3_Iji)
          deallocate(k4_Iji)  
          deallocate(link_back)       
       endif          
          
       allocate(Ii(this%nodes))
       allocate(k1_Ii(this%nodes))
       allocate(k2_Ii(this%nodes))   
       allocate(k3_Ii(this%nodes))
       allocate(k4_Ii(this%nodes))       

       allocate(Ri(this%nodes))   
       allocate(k1_Ri(this%nodes))
       allocate(k2_Ri(this%nodes))   
       allocate(k3_Ri(this%nodes))
       allocate(k4_Ri(this%nodes))

       allocate(Iji(this%sumDeg))
       allocate(k1_Iji(this%sumDeg))
       allocate(k2_Iji(this%sumDeg))   
       allocate(k3_Iji(this%sumDeg))
       allocate(k4_Iji(this%sumDeg))         
       
       allocate(link_back(this%sumDeg))
       
       do i1 = 1, this%nodes
          do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
             i2 = this%listAdj(i12)
             l_back:do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1
                i3 = this%listAdj(i23)
                if(i3 == i1)then
                   link_back(i12) = i23
                   exit l_back
                endif
             enddo l_back
          enddo
       enddo
       
    end subroutine

    subroutine condicao_inicial_SIRS_rDMP(this, p0)
       class(grafo), intent(in) :: this
       real(dp), intent(in) :: p0
       
       integer :: l1
       integer(kind=8) :: l12
       
       !#################################################
       if( (.not.allocated(Ii)) .and. (.not.allocated(Ri)) .and. (.not.allocated(Iji)) ) stop "Chame a subrotina 'aloca', linha 20"
       !#################################################
       do l1 = 1, this%nodes
          !#################################################
          if(lista_de_clusters(l1) /= i_comp_gigante) cycle
          !#################################################          
          Ii(l1) = p0; Ri(l1) = 0.0d0
          !#################################################
          do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
             Iji(l12) = p0 * ( 1.0_dp - p0 )
          enddo
          !#################################################
       enddo
       t = 0.0_dp
    end subroutine
 
    subroutine rk4_SIRS_rDMP_grafo(this, t, dt, alp, lam, mu) !, f_Si, f_Ii
		    !use geraRede
		
		    class(grafo) :: this
         	real(dp) :: t, dt
        	integer :: n_sites
        	integer (kind=8) :: n_stubs
        	real(dp) :: alp, lam, mu
            
            n_sites = this%nodes
            n_stubs = this%sumDeg
  		!####################################################################################
		
		
		!####################################################################################	
                     
!                    f_Ri(this, t, Ri, Ii, n_sitios, alf, mu)
		k1_Ri = dt * f_Ri(this, t, Ri, Ii, n_sites, alp, mu)
		
		!            f_Ii(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
		k1_Ii = dt * f_Ii(this, t, Ri, Ii, Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
		k1_Iji = dt * f_Iji(this, t, Ri, Ii, Iji, n_sites, n_stubs, mu, lam)
		!###############################################################
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n_sitios, alp, mu)
		k2_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n_sites, alp, mu)
		
!		             f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sites, n_stubs, mu, lam)
       
        !###############################################################
                     !f_Ri(this, t             , Ri                 , Ii                 , n_sitios, alf, mu)
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n_sitios, alp, mu)
		k3_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n_sites, alp, mu)

!                    f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sites, n_stubs, mu, lam)
       !################################################################                   
                    !f_Ri(this, t     , Ri        , Ii        , n_sitios, alf, mu)
!                    f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n_sitios, alp, mu)
		k4_Ri = dt * f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n_sites, alp, mu)
		
!                    f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Ii = dt * f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Iji = dt * f_Iji(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sites, n_stubs, mu, lam)
		!###############################################################
			
		Ri = Ri + (1.0_dp/6.0_dp) * (k1_Ri + 2.0_dp * k2_Ri + 2.0_dp * k3_Ri + k4_Ri)
		Ii = Ii + (1.0_dp/6.0_dp) * (k1_Ii + 2.0_dp * k2_Ii + 2.0_dp * k3_Ii + k4_Ii)
		Iji = Iji + (1.0_dp/6.0_dp) * (k1_Iji + 2.0_dp * k2_Iji + 2.0_dp * k3_Iji + k4_Iji)
  		!###############################################################   
    end subroutine


    subroutine rk4_SIRS_rDMP_grafo_teste(this, t, dt, alp, lam, mu)
        !use geraRede
        class(grafo) :: this
        real(dp) :: t, dt
        integer :: n_sites
        integer (kind=8) :: n_stubs
        real(dp) :: alp, lam, mu
        real(dp), parameter :: fator = 1.0_dp/6.0_dp
        !###############################################################
        n_sites = this%nodes; n_stubs = this%sumDeg
        !###############################################################             
!                    f_Ri(this, t, Ri, Ii, n_sitios, alf, mu)
		k1_Ri = dt * f_Ri(this, t, Ri, Ii, n_sites, alp, mu)
		
		!            f_Ii(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
		k1_Ii = dt * f_Ii(this, t, Ri, Ii, Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
		k1_Iji = dt * f_Iji(this, t, Ri, Ii, Iji, n_sites, n_stubs, mu, lam)
		!###############################################################
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n_sitios, alp, mu)
		k2_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n_sites, alp, mu)
		
!		             f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sites, n_stubs, mu, lam)
       
        !###############################################################
                     !f_Ri(this, t             , Ri                 , Ii                 , n_sitios, alf, mu)
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n_sitios, alp, mu)
		k3_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n_sites, alp, mu)

!                    f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sites, n_stubs, mu, lam)
       !################################################################                   
                    !f_Ri(this, t     , Ri        , Ii        , n_sitios, alf, mu)
!                    f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n_sitios, alp, mu)
		k4_Ri = dt * f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n_sites, alp, mu)
		
!                    f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Ii = dt * f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sites, n_stubs, mu, lam)

!                     f_Iji(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Iji = dt * f_Iji(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sites, n_stubs, mu, lam)
		!###############################################################
		Ri_tild = Ri + fator * (k1_Ri + 2.0_dp * k2_Ri + 2.0_dp * k3_Ri + k4_Ri)
		Ii_tild = Ii + fator * (k1_Ii + 2.0_dp * k2_Ii + 2.0_dp * k3_Ii + k4_Ii)
		Iji_tild = Iji + fator * (k1_Iji + 2.0_dp * k2_Iji + 2.0_dp * k3_Iji + k4_Iji)
  		!###############################################################   
    end subroutine

!#######################################################################
		function f_Iji(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
			!use types
			!use geraRede
			
			class(grafo) :: this
			real(dp), intent(in) :: t, mu, lam
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			real(dp) :: Ri(n_sitios), Ii(n_sitios), Iji(n_arestas)
			real(dp) :: f_Iji(n_arestas)
			real(dp) :: sumIkj
			!###############################
			integer :: l1, l2, l3
            integer(kind=8) :: l12, l21, l13, l31, l23, l32
			!###############################

			! Sitio i
			do l1 = 1, this%nodes
			    if(lista_de_clusters(l1) /= i_comp_gigante) cycle
				do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
					!---------------------------------------------------
                    l2 = this%listAdj(l12)
				    !---------------------------------------------------
				    l21 = link_back(l12)
				    !---------------------------------------------------
					f_Iji(l21) = -1.0_dp * mu * Iji(l21)
				    !---------------------------------------------------
					! Vizinhos k, do sitio j. 
				    !---------------------------------------------------
					sumIkj = 0.0_dp
				    !---------------------------------------------------
					do l23 = this%aux(l2), this%aux(l2) + this%deg(l2) - 1
				        !-----------------------------------------------
					    !l3 = this%listAdj(l23)
					    l32 = link_back(l23)
                        !===============================================
						! O sitio i nao vale como vizinho de j aqui.
                        !===============================================
						if(l32 == l12) cycle
                        !===============================================
						sumIkj = sumIkj + Iji(l32)
						!===============================================
					enddo
					!===================================================
					f_Iji(l21) = f_Iji(l21) + lam * ( 1.0d0 - Ri(l2) - Ii(l2) ) * sumIkj
					!===================================================
				enddo
			enddo
		end function

		!#################################################
		!	F_Si
		!###############################################################
		function f_Ri(this, t, Ri, Ii, n_sitios, alf, mu)
			class(grafo) :: this
			real(dp), intent(in) :: t, alf, mu
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			real(dp) :: Ri(n_sitios), Ii(n_sitios)
			real(dp) :: f_Ri(n_sitios)
			!###########################################################
			integer :: l1, l2
            integer(kind=8) :: l12
			!###########################################################
			! Sitio i
			do l1 = 1, this%nodes
			    !#######################################################
			    if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			    !#######################################################
				f_Ri(l1) = mu * Ii(l1) - alf * Ri(l1)
				!#######################################################
			enddo
		end function
		!###############################################################
		!	F_Ii
		!###############################################################
		function f_Ii(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
			!###########################################################
			class(grafo) :: this
			!###########################################################
			real(dp), intent(in) :: t, mu, lam
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			!###########################################################
			real(dp) :: Ri(n_sitios), Ii(n_sitios), Iji(n_arestas)
			real(dp) :: f_Ii(n_sitios)
			real(dp) :: sumIji
			!###########################################################
			integer :: l1, l2, l3
            integer(kind=8) :: l12, l13, l23, l21
			!###########################################################
			! Sitio i.
			do l1 = 1, this%nodes
				!#################################################
				if(lista_de_clusters(l1) /= i_comp_gigante) cycle
				!#######################################################
				f_Ii(l1) = -1.0_dp * mu * Ii(l1)
				!#######################################################
				sumIji = 0.0_dp
				!#######################################################
				do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
					l21 = link_back(l12)
					sumIji = sumIji +  Iji(l21)
				enddo
				!#######################################################
				f_Ii(l1) = f_Ii(l1) + lam * (1.0d0 - Ii(l1) - Ri(l1) ) * sumIji
				!#######################################################
			enddo
		end function
!#######################################################################
! Define funcoes
end module

program main
   use hashimoto
   use geraRede
   use mod_rndgen
   use mod_tools
   use types
   use mod_SIRS_rDMP_numerico
   
   implicit none   
!#######################################################################   
   type(grafo_PL_UCM) :: rede
!#######################################################################   
   real(dp) :: dlamb
   real(dp) :: lamb0, lambdaf
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   real(dp) :: rho, rho0
   real(dp) :: reco
   real(dp) :: tole
!#######################################################################
   !type(rndgen) :: gen1
   integer :: seed(50)
   integer :: seed1
   type(rndgen) :: ger_inic
!#######################################################################
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
!#######################################################################   
   character(len=500) :: t_vs_Im
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo_rede_real
   character(len=5) :: tam_char
   character(len=5) :: gama_char
   character(len=5) :: indice
!#######################################################################   
   integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9
   integer (kind=8) :: i12
   logical :: T_vs
!#######################################################################
   integer(kind=8) :: sumDeg2
   real(dp) :: t0, tf
   real(dp) :: dt0, dt
   real(dp) :: vec_dt(50)
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   integer :: ntempo
   
   integer :: n_it_rho, n_it_Y4, n_sub, per_conv
   
   real(dp) :: p0
   character(len=300) :: cwd, resultados, tipoCorte
   character(len=1000) :: local
   !##############################################################
   
   real(dp) :: qm, q2m

   integer :: nargus, ind_lamb   
   character(len=3) :: char_ind_lamb
    
   character(len=20) :: buffer
   integer :: divisor

   logical :: teveLeitura
   
   character(len=1000) :: nomeArquivo

   integer :: ind_amostra
   real(dp) ::lamb_C
   real(dp) :: tempoEscrita, tempoEscrita0
   real(dp) :: lamb_Copia
   logical :: existe, existe2
   integer :: st, st2
   logical :: usouCopia
   real(dp) :: lambda_Ultimo_Index
   integer :: int_soCalculaIPR
   logical :: soCalculaIPR
   
   logical :: taAberta
   real(dp) :: Y4_old
   real(dp) :: sumIi2
   logical :: Ii_neg
   logical :: Ii_gt_1
   real(dp) :: LEV_Y4(2)
   real(dp) :: lbd_lido, Y4_lido
   character(len = 100) :: C_lamb
   character(len = 100) :: C_tempo, C_tempoEscrita
   character(len = 100) :: C_I_i, C_I_ji, C_R_i
   integer :: st10
   logical :: existe10
   integer :: gasta_aleatorio
!=======================================================================   
! ----------------Name space dos arquivos de backup---------------------
!=======================================================================
   C_lamb = 'Copia_lambda'
   C_tempo = 'Copia_tempo'
   C_tempoEscrita = 'Copia_tempoEscrita'
   C_I_i = 'Copia_Ii'
   C_I_ji = 'Copia_Iji'
   C_R_i = 'Copia_R_i'      
!=======================================================================
! Esse seed1 estah de acordo com a rotina que calcula o metodo
! estocastico no NSSC, pasta SIRS_Estocastico,
! sob o titulo main_SIRS_Estocastico.f90.
!=======================================================================
   seed1=947361823
!=======================================================================
 
   !tipoCorte ='_Rigido'
   tipoCorte ='_2sqrtN'
   !tipoCorte ='_sqrtN'
   
   resultados = 'Rst_rDMP_Corte'//trim(adjustl(tipoCorte))

   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
       
   local = trim(adjustl(resultados))//'/'

!=======================================================================
   call entradaArgumentos()
!=======================================================================

!#######################################################################
   !====================================================================
   !  Se der ruim no dt, ele eh dividido por dois.
   !====================================================================   
   dt0 = 1.0d-1 
   dt = dt0

   vec_dt(1) = dt0
   do i1 = 2, size(vec_dt)
      vec_dt(i1) = vec_dt(i1 - 1) * 0.5d0
   enddo
   !====================================================================
   ! A principio, t0 = 0, mas, se houver um estado salvo, muda
   ! para o t salvo no arquivo.
   !====================================================================
   t0 = 0.0_dp
   tf = 10000.0_dp
   
   tempoEscrita0 = 1.0d0
   
   ntempo = int( (tf- t0)/dt )
   tole = 1d-7
   per_conv = int( 10.0_dp/dt )

   !====================================================================
   ! Como escolheremos a nossa sequencia de sementes?
   ! Fiz uma cagada nao padronizando direito minhas sementes.
   ! Ou eu arrisco rodar de acordo com os processos estocasticos,
   ! sem nenhuma garantia de estar certo, apenas checando a distribuicao
   ! de graus, clustering, knn e etc. disponivel,
   ! se estiverem identicos aos obtidos na rotina atual,
   ! ou eu rodo tudo de novo! A parte 'boa' eh que pra essa analise,
   ! vou precisar rodar apenas uma amostra por vez, ao inves de ter que
   ! rodar 10 amostras para cada um dos 9 tamanhos, para cada um
   ! dos 3 gamas e para cada um dos 3 valores de alfa, e para cada uma
   ! das 3 abordagens (Estocastico, 2 x Campo medio), totalizando
   ! 10 * 9 * 3 * 3 * 3 = 2430 curvas.
   ! Ao contrario, serao 
   ! 4 tamanhos * 3 gamas * 1 amostra * 3 alfas * 3 abordagens
   ! =  108 curvas. Isso so para as redes PL,
   ! sendo que eu fixei o corte!
   ! Lembrando que eu posso fazer para o grafo estrela, o grafo RRN,
   ! o grafo Estrela dentro do RRN, grafo Roda totalizando
   ! 4 tamanhos * 4 grafos * 1 amostra * 3 alfas * 3 abordagens
   ! = 144.
   ! Totalizando 144 + 108 = 252
   !====================================================================   
   if(.True.)then
      !=================================================================
      ! Vamos usar esse,
      ! que parece ser a versao antiga (mas que pode concordar com algumas
      ! amostras de rede que temos.
      !=================================================================
      call ger_inic%init(seed1)
      i2 = 1
      do i1 = 1, 5000
         if(mod(i1,100) > 0) cycle
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1      
      enddo
   else
      !=================================================================
      ! Ou esse,
      ! que eh da versao mais nova e que concorda com algumas redes.
      !=================================================================   
      i2 = 1
      do i1 = 1, 1000
         if(mod(i1,100) > 0)then
            gasta_aleatorio = ger_inic%int(100000000,999999999)
            cycle
         endif
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1
      enddo
   endif   
   !====================================================================
   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
      
   local = trim(adjustl(trim(adjustl(local))//'gam_'//trim(adjustl(gama_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )   

   local = trim(adjustl(trim(adjustl(local))//'ams_'//trim(adjustl(indice))//'/'))
   
   call system('mkdir -p '//trim(adjustl(local)) )
!##########################################################################################
!				Inicia grafo
!##########################################################################################
  call criaRedeEClassificaClusters(rede, tam_rede, grau_min, grau_max, gama_exp, seed(ind_amostra))
!#######################################################################  
  indice = trim(adjustl(indice))
  indice = trim(adjustl('_'//trim(adjustl(indice)))) 
!#######################################################################  
   write(*,*) "######################Dados da Rede######################"
   write(*,*) ""
   write(*,*) "Tamanho da rede ", rede%nodes, "."
   write(*,*) ""
   write(*,*) "Fracao correspondente aa componente gigante ", 100.0 * comp_gigante/rede%nodes,"%", "."
   write(*,*) ""
   write(*,*) "Expoente da distribuicao de graus da rede ", gama_exp, "."
   write(*,*) ""
   write(*,*) "Grau minimo ", rede%degMin, ".", " Grau maximo ", rede%degMax, "."
   write(*,*) ""
!#######################################################################
   write(alp_char2, '(f9.3)') alp

   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//trim(adjustl(tipoCorte))//'/'))
      
   call system('mkdir -p '//trim(adjustl(local)) )

   call kNN_e_clustering(rede)
!#######################################################################
   write(*,*) "######################Dados temporais######################"
   write(*,*) ""
   write(*,*) "Instante inicial ", t0, "."	
   write(*,*) ""
   write(*,*) "Resolucao ", dt, "."
   write(*,*) ""
   write(*,*) "Instante final ", tf, "."
   write(*,*) ""
   write(*,*) "Quantidade de pontos ", ntempo, "."
   write(*,*) ""
   write(*,*) "Tolerancia considerada no criterio de conv. ", tole, "."
!#######################################################################
   write(*,*) "###############Parametros Epidemicos#####################"
   write(*,*) "Probabilidade de sitio estar infectado ", p0, "."
   write(*,*) ""
   write(*,*) "Taxa de recuperacao ", mu, "."	
   write(*,*) ""
   write(*,*) "Taxa de enfraquecimento imunologico ", alp, "."
   write(*,*) ""
   write(*,*) "Taxa de infeccao inicial ", lamb0, "."
   write(*,*) ""
   write(*,*) "Resolucao da taxa de infeccao ", dlamb, "."
   write(*,*) ""
   write(*,*) "Taxa de infeccao final ", lambdaf, "." 
   write(*,*) ""
   write(*,*) "Numero de pontos de taxa infeccao final ", nlamb, "."
   write(*,*) ""
   write(*,*) "#########################################################"

   call calcula_P_grau(rede)
              
!#######################################################################

!#######################################################################   

!#######################################################################
   
   LEV_Y4 = LEV_e_IPR_NBC_Unc(rede)
   
   nomeArquivo = trim(adjustl(local))//trim(adjustl('N_vs_Y4_rDMP_Teorico'))//'.dat'
   
   open(885, file = trim(adjustl(nomeArquivo)), status = 'unknown')
   
   write(885,*) rede%nodes, LEV_Y4(2)
   close(885)   

   nomeArquivo = trim(adjustl(local))//trim(adjustl('N_vs_1_sobre_Lambda_H_Teorico'))//'.dat'

   open(885, file = trim(adjustl(nomeArquivo)), status = 'unknown')
   write(885,*) rede%nodes, mu/LEV_Y4(1)

   close(885)   
   !stop "Soh calcularah o LEV e IPR Teorico. Venha e edite, se quiser mais de mim."
   !====================================================================
   if( nargus == 10)then
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'_lamb_index_0.dat'
   elseif( nargus == 9)then
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'.dat'
   endif
   !====================================================================
   ! Testo se os arquivos acima existem.
   ! Logo abaixo, faco leituras, caso eles existam.
   ! Em caso negativo, eu chamo as rotinas apropriadas.
   !====================================================================
   inquire( file = trim(adjustl(nomeArquivo)), exist = existe10)
   !====================================================================
   if( .not. existe10 )then
      write(*,*) "Arquivo nao existe"
      call aloca_listas_e_matrizHashimoto(rede)
      call metodo_potencia_Hashimoto(rede, lamb0, mu)
   else
      open(unit = 71, file = trim(adjustl(nomeArquivo)), status = 'old') 
      read(71, *, iostat = st10) lamb0, Y4
      if( st10 /= 0)then
         write(*,*) "Nao houve leitura"
         call aloca_listas_e_matrizHashimoto(rede)
         call metodo_potencia_Hashimoto(rede, lamb0, mu)      
      endif
      close(71)
   endif
   !====================================================================
   if( ( .not. existe10 ) .or. ( st10 /= 0))then
      if( nargus == 10)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'_lamb_index_0.dat'
         open(337, file=trim(adjustl(nomeArquivo)), status = 'unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'_lamb_index_0.dat'
         open(334, file=trim(adjustl(nomeArquivo)), status = 'unknown')       
      elseif( nargus == 9)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'.dat'
         open(337, file=trim(adjustl(nomeArquivo)), access = 'append', status = 'unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'.dat'
         open(334, file=trim(adjustl(nomeArquivo)), access = 'append', status = 'unknown')
      endif
      
      write(*,*) "Ou nao existe ou leu ruim"
      write(334,*) lamb0, 0.0d0
      close(334)
      write(337,*) lamb0, Y4
      close(337)
      if( (nargus == 10) .and. (ind_lamb == 0) ) stop
   endif
   !====================================================================
   call aloca_listas_SIRS_rDMP(rede)
   !====================================================================
   !====================================================================
   ! A principio
   !====================================================================
   usouCopia = .False.

   if( nargus == 9 )then
      !===============================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'.dat'                        

      inquire(unit=334, opened=taAberta)
      if(taAberta) close(334)             
      open(334, file=trim(adjustl(nomeArquivo)), status='unknown')
      !-----------------------------------------------------------------
      lambdaf = 0.0d0
      !-----------------------------------------------------------------      
      lerho:do
         read(334, *, iostat = st) lamb, rho
         if(st /= 0) exit lerho
         if( lamb > lambdaf) lambdaf = lamb
      enddo lerho
      if(lambdaf == 0.0d0)then
         lambdaf = lamb0 + dlamb * 500.0d0
         lamb = lamb0
      else
         lamb = lambdaf + dlamb
         lambdaf = lamb + dlamb * 500.0d0
      endif
      !-----------------------------------------------------------------
      close(334)
      !-----------------------------------------------------------------      
      open(334, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'.dat'
      !-----------------------------------------------------------------
      inquire(unit=335, opened=taAberta)
      if(taAberta) close(335)
                         
      open(335, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
      !===============================================================   
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'.dat'

      inquire(unit=337, opened=taAberta)
      if(taAberta) close(337)
                         
      open(337, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
   elseif( nargus == 10 )then
      inquire(unit = 334, opened = taAberta)
      if(taAberta)close(334)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'       
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo)), status = 'unknown')
      !=================================================================
      inquire(unit = 335, opened = taAberta)
      if(taAberta)close(335)
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat' 
      open(335, file=trim(adjustl(nomeArquivo)), status = 'unknown')
      !=================================================================
      inquire(unit = 337, opened = taAberta)
      if(taAberta)close(337)
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat' 
      open(337, file=trim(adjustl(nomeArquivo)), status = 'unknown')
      !=================================================================
      lambdaf = lamb0 + 1.0d0 * ind_lamb * dlamb
      lamb = lambdaf
   endif
   !====================================================================
   ! Eu testo se os arquivos supracitados existem ou contem
   ! dados apropriados. Em caso negativo, chamo os modulos
   ! apropriados e agora eu os escrevo.
   !====================================================================
   if( rede%nodes == 10**4 )then
      
      if(allocated(fi))then
         call calcula_distancias_ao_hub(rede)
         write(*,*) "Chamou calcula_distancias_ao_hub"
      
         nomeArquivo = trim(adjustl('fi_rDMP_N_10to4_gm_'))//trim(adjustl(gama_char))//'_tam_'//trim(adjustl(tam_char))//'_alp_TODOS_'//'ams_'//trim(adjustl(indice))
         open(999, file = trim(adjustl(nomeArquivo))//'.dat', status = 'unknown')
         do i1 = 1, size(fi)
            write(999,*) rede%deg(i1), lista_distancias(i1), (fi(i1))**2.0d0
         enddo
 
         close(999)
         call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
      
         deallocate(lista_distancias)
      endif
   endif
   !====================================================================
   ! Loop lambda
   !====================================================================
   if(allocated(Ii_tild)) deallocate(Ii_tild)
   allocate(Ii_tild(rede%nodes))

   if(allocated(Ri_tild)) deallocate(Ri_tild)
   allocate(Ri_tild(rede%nodes))

   if(allocated(Iji_tild)) deallocate(Iji_tild)
   allocate(Iji_tild(rede%sumDeg))
   !====================================================================
   ll: do while(lamb <= lambdaf)
           !============================================================
           ! tempo
           !============================================================
           n_it_rho = 0
           n_it_Y4 = 0
           n_sub = 0
           !============================================================
           if( nargus == 10)then              
              inquire(unit = 333, opened = taAberta)
              if(taAberta)close(333)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
              !=========================================================
              open(333, file=trim(adjustl(nomeArquivo)), access = 'append', status='unknown')   
              !=========================================================
              inquire(unit = 336, opened = taAberta)
              if(taAberta)close(336)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
              !=========================================================
              open(336, file=trim(adjustl(nomeArquivo)), access = 'append', status='unknown')   
              !=========================================================
           elseif( nargus == 9)then
              inquire(unit = 333, opened = taAberta)
              if(taAberta)close(333)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_rDMP'))//'_global.dat'
              !============================================================
              open(333, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================

              inquire(unit = 336, opened = taAberta)
              if(taAberta)close(336)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'_global.dat'
              !============================================================
              open(336, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================
           endif           
           !============================================================
           ! Se sair com Ii = 0 da rotina checaBackUp, 
           ! chama subrotina condicao_inicial_SIRS_rDMP.
           !============================================================
           call checaBackUp(rede, trim(adjustl(local)))
           !============================================================
           sumIi2 = sum(Ii)
           if( sumIi2 == 0.0d0)then
              call condicao_inicial_SIRS_rDMP(rede, p0)
              rho0 = p0
              rho = rho0
              reco = 0.01              
              t = 0.0d0
              Y4_old = 1.0_dp/(1.0_dp * comp_gigante)
              Y4 = Y4_old
              tempoEscrita = tempoEscrita0
           else
              !=========================================================
              ! Na subrotina checaBackup, calcula-se rho e Y4,
              ! mas nao rho0 e Y4_old.
              !=========================================================
              rho0 = rho
              Y4_old = Y4
              reco = sum(Ri)/comp_gigante
           endif
           !============================================================
           dt = dt0
           !============================================================
           lt:do while(t <= tf)
              
              call tempoDeEscrever()
              
              write(333, *) t, Y4
              write(336, *) t, rho
              !=========================================================
              write(339, *) t, reco
              !=========================================================
              ! Aqui decidimos se o dt vai ou nao estragar tudo.
              ! Iniciamos Ii_tilda com um valor > 1.0d0
              !=========================================================
              choose_dt:do i1 = 1, size(vec_dt)
                 !======================================================
                 dt = vec_dt(i1)
                 call rk4_SIRS_rDMP_grafo_teste(rede, t, dt, alp, lamb, mu)
                 !======================================================
                 rho = 0.0d0; reco = 0.0d0
                 Ii_neg = .False.
                 Ii_gt_1 = .False.
                 !======================================================
                 conta_Ii:do i2 = 1, size(Ii_tild)   ! Calcula rho
                    if(lista_de_clusters(i2) /= i_comp_gigante) cycle
                    if( (Ii_tild(i2) < 0.0d0) .or. (Ri_tild(i2) < 0.0d0) .or. ( (1.0d0 - Ii_tild(i2) - Ri_tild(i2)) <= 0.0d0 ) )then
                       Ii_neg = .True.
                       exit conta_Ii
                    elseif( (Ii_tild(i2) > 1.0d0) .or. (Ri_tild(i2) > 1.0d0) .or. ( ( 1.0d0 - Ii_tild(i2) - Ri_tild(i2)) > 1.0d0 ))then
                       Ii_gt_1 = .True.
                       exit conta_Ii
                    endif
                    rho = rho + Ii_tild(i2)
                    reco = reco + Ri_tild(i2)
                 enddo conta_Ii
                 !======================================================
                 if(.not. Ii_neg  )then
                    if( .not. Ii_gt_1 )then
                       Ii = Ii_tild
                       Ri = Ri_tild
                       Iji = Iji_tild
                       exit choose_dt
                    endif
                 endif
              enddo choose_dt
              !=========================================================
              !if(i1 > 1) write(*,*) "dt atualizado para = ", dt, "de indice = ", i1
              !=========================================================
              !call rk4(dt, t, rede, rede%nodes, rede%sumDeg, alp, lamb, mu)
              !=========================================================  
              if( ( Ii_neg == .True. ) .or. ( Ii_gt_1 == .True. ) )then
                 stop "Esgotou a precisao da Maquina, sem obter resultados fisicos"
              endif
              !=========================================================
              !call rk4_SIRS_rDMP_grafo(rede, t, dt, alp, lamb, mu)
              !=========================================================
              rho  = rho/( 1.0_dp * comp_gigante )
              reco  = reco/( 1.0_dp * comp_gigante )
              !=========================================================
              sumIi2 = (sum(Ii**2.0_dp))**0.5_dp
         
              Y4 = sum( (Ii/sumIi2)**4.0_dp )
           
              t = t + dt
           
              if( abs(rho - rho0)/rho0 < tole )then
                 n_it_rho = n_it_rho + 1
                 if( n_it_rho >= int(5.0d0/dt) )then
                    write(*,*) "Conv. cm lbd = ", lamb, " e rho = ", rho, "."
                    exit lt
                 endif
              elseif( abs(Y4 - Y4_old)/Y4_old < tole )then
                 n_it_Y4 = n_it_Y4 + 1
                 if( n_it_Y4 >= int(5.0d0/dt) )then
                    write(*,*) "Conv. cm lbd = ", lamb, " e Y4 = ", Y4, "."
                    exit lt
                 endif
              else
                 n_it_Y4 = 0; n_it_rho = 0
              endif
              
              rho0 = rho
              Y4_old = Y4
           enddo lt
           
           close(333)
           close(336)
        
           if( ( (n_it_rho) >= int(5.0d0/dt)) .or. ( (n_it_Y4) >= int(5.0d0/dt)) )then
              write(334, *) lamb, rho           
              write(335, *) lamb, t
              
              write(337, *) lamb, Y4
           endif
           
           lamb = lamb + dlamb

		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_global.dat')
		      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_global.tar.gz')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.tar.gz')
		   endif
        enddo ll
        close(334)
        close(335)
        close(337)
     
   contains
   
      function LEV_e_IPR_NBC_Unc(this)
         class(grafo), intent(in) :: this
         real(dp) :: LEV_e_IPR_NBC_Unc(2)
         integer :: j1, j2
         integer(kind = 8) :: j12, j23
         real(dp) :: numerador1, numerador2, denominador
         real(dp), allocatable :: xi(:)
         real(dp) :: XI_total
         real(dp) :: norm_xi
         !--------------------------------------------------------------
         if(allocated(xi)) deallocate(xi)
         allocate(xi(this%nodes))
         !--------------------------------------------------------------         
         denominador = 0.0d0
         do j1 = 1, this%nodes
            denominador = denominador + 1.0d0 * this%deg(j1) * ( 1.0d0 * this%deg(j1) - 1.0d0)
         enddo
         !--------------------------------------------------------------
         numerador2 = 0.0d0
         do j1 = 1, this%nodes
            !-----------------------------------------------------------
            numerador1 = 0.0d0 
            !-----------------------------------------------------------
            do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
               j2 = this%listAdj(j12)
               numerador1 = numerador1 + 1.0d0 * (this%deg(j2) - 1)
            enddo
            numerador2 = numerador2 + (1.0d0 * (this%deg(j1) - 1)) * numerador1
            !-----------------------------------------------------------
            xi(j1) = numerador1/denominador
            !-----------------------------------------------------------
         enddo
         norm_xi = (sum(xi**2.0d0))**0.5d0
         !--------------------------------------------------------------
         xi = xi/norm_xi
         !--------------------------------------------------------------
         LEV_e_IPR_NBC_Unc(1) = numerador2/denominador
         LEV_e_IPR_NBC_Unc(2) = sum((xi)**4.0d0)
         !--------------------------------------------------------------         
      end function
   !====================================================================
      subroutine criaRedeEClassificaClusters(this, tam_rede1, grau_min1, grau_max1, gama_exp1, seme1)
        type(grafo_PL_UCM) :: this
        integer :: tam_rede1
        integer :: grau_min1
        real(dp) :: grau_max1
        real(dp) :: gama_exp1
        integer :: seme1
!#######################################################################   
        call this%iniciaGrafo(tam_rede1)
!#######################################################################
        call this%inicia(grau_min1, grau_max1, gama_exp1, seme1)
!#######################################################################
        call this%liga(seme1, .False.) 
!#######################################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
!#######################################################################
      end subroutine
!#######################################################################

      subroutine entradaArgumentos()
         nargus = iargc()


         if(nargus == 9)then
            !#############################
            ! Amostra
            !#############################
            write(*,*) "Recebeu 9 arqumentos"
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra

            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf

            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
            
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif               
            
         elseif(nargus == 10)then
            write(*,*) "Recebeu 10 arqumentos"
            !#############################
            ! Amostra
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra

            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf

            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
      
            !#############################
            ! Indice do ponto
            ! Quando fizermos
            ! paralelizacao burra.
            !#############################
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_lamb   
            write(char_ind_lamb, '(I0)') ind_lamb
        
            char_ind_lamb = trim(adjustl(char_ind_lamb))
            
            call getarg(10, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif   
         else
            stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
         endif


         if( trim(adjustl(resultados)) == 'Rst_rDMP_Corte_Rigido' )then
               call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
         elseif(trim(adjustl(resultados)) == 'Rst_rDMP_Corte_sqrtN' )then
               grau_max = (1.0_dp * tam_rede)**(0.5_dp)    
         elseif(trim(adjustl(resultados)) == 'Rst_rDMP_Corte_2sqrtN' )then
               grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
         endif


         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!#######################################################################
!  Indice da semente
!#######################################################################
         write(gama_char,'(f5.2)') gama_exp
!#######################################################################
         write(indice,'(I0)') ind_amostra 
!#######################################################################   
         if(tam_rede == 10)then
            tam_char = '10'
         elseif(tam_rede == 100)then
            tam_char = '100'
         elseif(tam_rede == 1000)then
            tam_char = '1k'
         elseif(tam_rede == 3000)then
            tam_char = '3k'
         elseif(tam_rede == 10000)then
            tam_char = '10k'
         elseif(tam_rede == 30000)then
            tam_char = '30k'
         elseif(tam_rede == 100000)then
            tam_char = '100k'
         elseif(tam_rede == 300000)then
            tam_char = '300k'
         elseif(tam_rede == 1000000)then
            tam_char = '1M'
         elseif(tam_rede == 3000000)then
            tam_char = '3M'
         elseif(tam_rede == 10000000)then
            tam_char = '10M'
         elseif(tam_rede == 30000000)then
            tam_char = '30M'
         elseif(tam_rede == 100000000)then
       	    tam_char = '100M'
         else
            stop 'Escolha um tamanho de rede dentro do catalogo'
         endif         
      end subroutine

!=======================================================================
   subroutine calcula_P_grau(this)      
      class(grafo) :: this
      integer :: i1
!#######################################################################
      if(allocated(P_grau)) deallocate(P_grau)
      allocate(P_grau(this%degMin:this%degMax))
      P_grau = 0_dp
      do i1 = 1, this%nodes
         P_grau(this%deg(i1)) = P_grau(this%deg(i1)) + 1.0_dp
      enddo
      P_grau = P_grau/(1.0_dp * this%nodes)
      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      open(800, file=trim(adjustl(arq_1)), status='unknown')
      do i1 = this%degMin, this%degMax
         if(P_grau(i1) == 0.0_dp) cycle
         write(800,*) i1, P_grau(i1)
      enddo
      close(800)
      deallocate(P_grau)
!#######################################################################      
   end subroutine
!=======================================================================
   subroutine kNN_e_clustering(this)
      class(grafo) :: this

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call calcula_k_nn(this, .True., 800, trim(adjustl(arq_1)))
      close(800)

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call clustering(this,.True., 800, trim(adjustl(arq_1)))
      close(800)
   end subroutine
!=======================================================================   
	subroutine checaBackUp(this, pasta_fonte)
		 !--------------------------------------------------------------
		 class(grafo) :: this
		 character(len = *) :: pasta_fonte
		 character(len = 1000) :: nomeArquivo
		 integer :: i1, i2, i3
		 integer(kind = 8) :: i12, i21, i23, i13, i31, i32
		 logical :: existe
		 integer :: st
		 real(dp), allocatable :: rho_i1(:)
		 real(dp) :: norm_rhoi1
		 character(len = 100) :: C_lamb_f, C_tempo_f, C_tempoEscrita_f, C_R_i_f, C_I_ji_f, C_I_i_f
		 real(dp) :: t_lido
		 real(dp) :: tempoEscrita_Lido
		 !--------------------------------------------------------------
		 pasta_fonte = trim(adjustl(pasta_fonte))
		 !--------------------------------------------------------------
		 if(nargus == 10 )then
			C_lamb_f = trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_tempo_f = trim(adjustl(C_tempo))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_tempoEscrita_f = trim(adjustl(C_tempoEscrita))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_R_i_f = trim(adjustl(C_R_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_I_ji_f = trim(adjustl(C_I_ji))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_I_i_f = trim(adjustl(C_I_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
		 elseif( nargus == 9)then
			C_lamb_f = trim(adjustl(C_lamb))//'_global'
			C_tempo_f = trim(adjustl(C_tempo))//'_global'
			C_tempoEscrita_f = trim(adjustl(C_tempoEscrita))//'_global'
			C_R_i_f = trim(adjustl(C_R_i))//'_global'
			C_I_ji_f = trim(adjustl(C_I_ji))//'_global'
			C_I_i_f = trim(adjustl(C_I_i))//'_global'
		 endif
		 Ii = 0.0d0
		!===========================================================
		!-------------------Leitura de LAMBDA-----------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_lamb_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe ) 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) lbd_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===========================================================
		!-------------------Leitura do TEMPO------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_tempo_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )
        !===============================================================
        ! C_tempo_f nao existe
        !===============================================================

		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		
		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) t_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===========================================================
		!-------------------Leitura do TEMPO DE ESCRITA-------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_tempoEscrita_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )
        !===============================================================
        ! C_tempoEscrita_f nao existe
        !===============================================================
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

		open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
		read(776,*, iostat=st) tempoEscrita_lido
		close(776)
		if( st /= 0 )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif
		!===========================================================
		!-------------------Leitura de I_i--------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_I_i_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe ) 

		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Eu zeraria Ii aqui, mas ele ainda nao foi lido.
		   ! Portanto, ele jah eh zero.
		   !============================================================
		   return
        endif
        !===============================================================
		do i1 = 1, size(Ii)
		   read(776,*, iostat=st) Ii(i1)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !============================================================
		      ! Eh preciso zerar Ii, afim de que meu programa entenda
		      ! que nao foi feita uma leitura apropriada do backup
		      ! e ele chame as condicoes iniciais apropriadas
		      !============================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
		!===========================================================
		!-------------------Leitura de I_ji-------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_I_ji_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
		 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   !============================================================
		   ! Preciso zerar pelo menos o Ii, 
		   ! senao, o programa vai entender que houve leitura correta.
		   !============================================================
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)
        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)
        
        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i12 = 1, size(Iji)
		   read(776,*, iostat=st) Iji(i12)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, pois foi a condicao
		      ! que eu escolhi para declarar que checaBackup
		      ! nao teve sucesso.
		      !=========================================================
			  Ii = 0.0d0
			  return
		   endif
		enddo
		close(776)
		!===========================================================
		!-------------------Leitura de R_i--------------------------
		!===========================================================
		nomeArquivo = trim(adjustl(pasta_fonte))//trim(adjustl(C_R_i_f))
		inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe ) 
		if( .not. existe )then
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
		endif

		call system ('tar -xzvf '//trim(adjustl(nomeArquivo))//'.tar.gz')

		inquire(unit=776, opened=taAberta)
		if(taAberta) close(776)

        !===============================================================
        inquire(file = trim(adjustl(nomeArquivo))//'.dat', exist = existe)

        if(existe)then
  		   open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
        else
		   if(nargus == 9)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		   elseif(nargus == 10)then
		      call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		   endif
		   Ii = 0.0d0
		   return
        endif
        !===============================================================
		do i1 = 1, size(Ri)
		   read(776,*, iostat=st) Ri(i1)
		   if( st /= 0 )then
		      if(nargus == 9)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_global.dat')
		      elseif(nargus == 10)then
		         call system('rm -r '//trim(adjustl(pasta_fonte))//trim(adjustl('Copia'))//'*_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat')
		      endif
		      !=========================================================
		      ! Preciso zerar ao menos Ii, afim de que meu programa
		      ! entenda que checaBackup nao achou backup
		      ! e chame a rotina condicao_inicial.
		      !=========================================================
		      Ii = 0.0d0
		      return
		   endif
		enddo
		close(776)
		!===============================================================
        ! Agora que leu tudo, sera feita a estatistica!
		!===============================================================
		rho = 0.0d0
		do i1 = 1, size(Ii)
		   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		   rho = rho + Ii(i1)
		enddo
		rho = rho/(1.0d0 * comp_gigante)
		!---------------------------------------------------------------    
		!if(allocated(rho_i1)) deallocate(rho_i1)
		!allocate(rho_i1(size(Ii)))
		!---------------------------------------------------------------
		! Calcular o NAV via mensagens que chegam dava um resultado
		! ligeiramente diferente para o Y4. Levemente maior.
		! Isso acontecia porque eu nao estava levando em conta a
		! equacao
		! \rho_i = (lambda/mu) * s_i * \sum_j A_{ij} * I_{ji}.
		! Apos fazer isso, a concordancia ficou muito melhor.
		!---------------------------------------------------------------
		!rho_i1 = 0.0d0
		!---------------------------------------------------------------
		!do i1 = 1, size(Ii)
		!   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		!   do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
		!	  i21 = link_back(i12)
		!	  rho_i1(i1) = rho_i1(i1) + Iji(i21)
		!   enddo
		!   rho_i1(i1) = lbd_lido * (1.0d0 - Ii(i1) - Ri(i1) )* rho_i1(i1)
		!enddo
		!norm_rhoi1 = 0.0d0

		!norm_rhoi1 = (sum(rho_i1**2.d0))**0.5d0
		!rho_i1 = rho_i1/norm_rhoi1
		!Y4 = sum(rho_i1**4.0d0)
		
		!write(*,*) t_lido, Y4

		!===============================================================
		!norm_rhoi1 = (sum(Ii**2.0d0))**0.5d0
		
		!Y4 = sum((Ii/norm_rhoi1)**4.0d0)
		!write(*,*) t_lido, Y4
		!===============================================================
		!deallocate(rho_i1)
		!stop
		!===============================================================
		norm_rhoi1 = (sum(Ii**2.0d0))**0.5d0
		
		Y4 = sum((Ii/norm_rhoi1)**4.0d0)
		!===============================================================
		Y4_old = Y4
		rho0 = rho
		!===============================================================
		t = t_lido
		lamb = lbd_lido
		tempoEscrita = tempoEscrita_Lido
	end subroutine

   subroutine TempoDeEscrever()                        
         if( t >= tempoEscrita)then
            !===========================================================       
            do i1 = 1, size(Ii)
               if( (Ii(i1) > 1.0d0) .or. (Ii(i1) < 0.0d0) ) return
            enddo
            !===========================================================
            do i1 = 1, size(Ri)
               if( ( Ri(i1) > 1.0d0) .or. (Ri(i1) < 0.0d0) )return
            enddo            
            !===========================================================
            do i12 = 1, size(Iji)
               if( ( Iji(i12) > 1.0d0) .or. (Iji(i12) < 0.0d0) )return
            enddo                              
            !===========================================================
            
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_i))//'_global'
            endif
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')            
            !===========================================================            
            do i1 = 1, size(Ii)
               write(777,*) Ii(i1)
            enddo
            !===========================================================
            close(777)                 
            !===========================================================
            
            !===========================================================
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================

            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_R_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_R_i))//'_global'
            endif            
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================            
            do i1 = 1, size(Ri)
               write(777,*) Ri(i1)
            enddo
            !===========================================================
            close(777)
            !===========================================================
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================

            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_ji))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_ji))//'_global'
            endif           
               
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================               
            do i12 = 1, size(Iji)
               write(777,*) Iji(i12)
            enddo
            !===========================================================
            close(777)
            !===========================================================            
            call system('tar -czvf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempo))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempo))//'_global'
            endif           
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) t
            close(777)
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))//'_global'
            endif  
            !===========================================================       
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) lamb
            close(777)
            !===========================================================
            tempoEscrita = tempoEscrita + tempoEscrita0
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempoEscrita))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempoEscrita))//'_global'
            endif           
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) tempoEscrita
            close(777)
            !===========================================================
            !===========================================================
         endif
   end subroutine

	function g_func(k_s, gam1, Nstr)
	  real(dp) :: k_s
	  real(dp) :: gam1
	  integer :: Nstr
	  real(dp) :: g_func
	  real(dp) :: Ap
	  !#################################################################
	  !   Quando g_func for usada pela primeira vez,
	  !   gbuffer = 0.0d0 e kmin = this%degMin,
	  !   k_s recebe o valor que se quer aproximar de kc
	  !   
	  !#################################################################      
	  !g_func = gbuffer
				  
	  Ap = Ap_list(int(k_s))
	  
	  Ap = 1.0_dp/Ap
	  
	  g_func = k_s - Ap**(1.0_dp/gam1) * (1.0_dp * Nstr)**(1.0_dp/gam1)
	  
	end function                          

	subroutine acha_cutoff_rigido(kming, gam_p, N_p)
	  integer :: kming
	  real(dp) :: gam_p
	  integer :: N_p
	  real(dp), parameter :: tol = 5d-5
	  real(dp) :: gminus, gmais
	  real(dp) :: kminus, kmais
	  real(dp) :: k_p
	  real(dp) :: gp
	  real(dp) :: Ap1
	  integer :: kl1, kl2
	  integer, parameter :: N_it = 10**3

	  !#################################################################
	  !   Inicio
	  !#################################################################           
	  kminus = 1.0_dp * kming
	  kmais = 1.5_dp * kming * (1.0_dp * N_p)**(1.0_dp/gam_p)
	  
	  if(allocated(Ap_list)) deallocate(Ap_list)
	  allocate(Ap_list(int(kming):(int(kmais))))
	  
	  gminus = g_func(kminus, gam_p, N_p)
	  
	  if(gminus >= 0.0_dp) stop "Precisamos de um valor de kminus para que gminus < 0"
	   
	  Ap1 = 0.0_dp
	  do kl1 = kming, int(kmais)
		 Ap1 = Ap1 + (1.0_dp * kl1)**(-gam_p)
		 Ap_list(kl1) = Ap1
	  enddo

	  gmais = g_func(kmais, gam_p, N_p)

	  if(gmais <= 0.0_dp) stop "Precisamos de um valor para kmais, tal que gmais > 0"
	  !#################################################################
	  !   Execucao
	  !#################################################################
	  kl1 = 1
	  do while(kl1 <= N_it)
		 k_p = kminus + (kmais - kminus)/2.0_dp
		 gp = g_func(k_p, gam_p, N_p)
		 if((gp == 0.0_dp).or.((kmais-kminus)/2.0_dp <= tol))then
			grau_max = k_p
			write(*,*) "Achou o grau max apos N = ", kl1, " iteracoes."
			exit
		 endif
		 kl1 = kl1 + 1
		 if(gminus * gp > 0.0_dp)then
			kminus = k_p
			gminus = gp
		 else
			kmais = k_p
		 endif
	  enddo      
	  !#################################################################
	  !   Final
	  !#################################################################
	  deallocate(Ap_list)                   
	end subroutine
	 
end program
