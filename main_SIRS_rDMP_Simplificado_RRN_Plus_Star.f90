module hashimoto
   use geraRede
   use mod_tools
   use mod_rndgen
   use types
   
   implicit none
   
   real(dp), allocatable :: x(:), y(:)
   integer, allocatable :: link_back2(:)
   real(dp), allocatable :: fi_rDMP(:)
   real(dp), parameter :: tol = 1.0d-7
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
            y(j12) = this%deg(j2)
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
                  viz_i3: do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1 ! link v --> w
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
		    !-----------------------------------------------------------
		    sum_xi = sum((rhoi_at + rhoi_int )/2.0d0)
		    sum_ki_xi = dot_product(1.0d0 * this%deg, (rhoi_at + rhoi_int )/2.0d0)
		    Mu_An = sum_ki_xi/sum_xi - 1
		    write(*,*) "Mu_An = ", Mu_An, "e ", "Mu_num = ", autovalor
		    !-----------------------------------------------------------
            if( this%nodes == 10**7 )then
               if(allocated(fi_rDMP)) deallocate(fi_rDMP)
               allocate(fi_rDMP(this%nodes))
               fi_rDMP = (rhoi_at + rhoi_int)/2.0d0
            endif
		    !-----------------------------------------------------------
		    lambda = mu1/autovalor
		    write(13,*) this%nodes, lambda
			write(*,*) 'Autovalor principal ', autovalor, ' encontrado com erro relativo ', erro, '. '
			
			write(*,*) " "
			write(*,*) 'Foram necessarios ', contador, ' passos.'
			write(*,*) " "

            deallocate(rhoi_at)
            deallocate(rhoi_int)

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
	
	real(dp), allocatable :: Ii(:), Ri(:), Iji(:)
	real(dp), allocatable :: k1_Ii(:), k1_Ri(:), k1_Iji(:) 
	real(dp), allocatable :: k2_Ii(:), k2_Ri(:), k2_Iji(:)
	real(dp), allocatable :: k3_Ii(:), k3_Ri(:), k3_Iji(:)
	real(dp), allocatable :: k4_Ii(:), k4_Ri(:), k4_Iji(:)
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
l_back:      do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1
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
   type(grafoRRN_Plus_Star) :: rede
!#######################################################################   
   real(dp) :: dlamb
   real(dp) :: lamb0, lambdaf
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   real(dp) :: rho, rho0
   real(dp) :: tole
!#######################################################################
   !type(rndgen) :: gen1
   integer :: seed(10)
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
   real(dp) :: dt
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   integer :: ntempo
   
   integer :: n_it, n_sub, per_conv
   
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

   integer, parameter :: ind_amostra = 1
   real(dp) ::lamb_C
   real(dp) :: tempoEscrita, tempoEscrita0
   real(dp) :: lamb_Copia
   logical :: existe, existe2
   integer :: st, st2
   logical :: usouCopia
   real(dp) :: lambda_Ultimo_Index
   integer :: int_soCalculaIPR
   logical :: soCalculaIPR
   integer :: grauRRN
   logical :: taAberta
   real(dp) :: Y4_old
   real(dp) :: sumIi2
   logical :: Ii_neg
   logical :: Ii_gt_1
   integer :: gasta_aleatorio
   character(len=4) :: char_grau_RRN
   integer :: grauStar
   character(len=4) :: char_grauStar_RRN
   character(len=1000) :: fi_char
   character(len = 100) :: C_lamb
   character(len = 100) :: C_tempo
   character(len = 100) :: C_I_i
   character(len = 100) :: C_I_ji
   character(len = 100) :: C_R_i
!=======================================================================   
   resultados = trim(adjustl('Rst_rDMP_RRN_Plus_Star'))
   call system('mkdir -p '//resultados)
!=======================================================================
! ----------------Name space dos arquivos de backup---------------------
!=======================================================================
   C_lamb = 'Copia_lambda'
   C_tempo = 'Copia_tempo'
   C_I_i = 'Copia_Ii'
   C_I_ji = 'Copia_Iji'
   C_R_i = 'Copia_R_i'
!=======================================================================
   local = trim(adjustl(resultados))//"/"
!=======================================================================   
   call entradaArgumentos()
!=======================================================================
   !====================================================================
   !  Se der ruim no dt, ele eh dividido por dois.
   !====================================================================   
   dt = 1.0d-1 
   !====================================================================
   ! A principio, t0 = 0, mas, se houver um estado salvo, muda
   ! para o t salvo no arquivo.
   !====================================================================
   t0 = 0.0_dp
   tf = 10000.0_dp
   
   tempoEscrita0 = 10.0d0
   
   ntempo = int( (tf- t0)/dt )
   tole = 1d-7
   per_conv = int( 10.0_dp/dt )

   seed1  = 967891968
   
   call ger_inic%init(seed1)

!=======================================================================
   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
      
   call system('mkdir -p '//trim(adjustl(local)) ) 

   write(alp_char2, '(f4.1)') alp   
   alp_char2 = trim(adjustl(alp_char2))
   
   local = trim(adjustl(trim(adjustl(local))//'/alp_'//trim(adjustl(alp_char2))//'/'))
   call system('mkdir -p '//trim(adjustl(local)))
   !====================================================================
   write(*,*) "O valor de alp eh ", alp
   write(*,*) "O valor de mu eh ", mu
   write(*,*) "O valor de lambda0 eh ", lamb0
   !====================================================================
!#######################################################################
!				Inicia grafo
!#######################################################################
   !====================================================================
   call criaRedeEClassificaClusters(rede, tam_rede, grauRRN, grauStar)
   !====================================================================
   write(*,*) ""
   write(*,*) "O tamanho da rede eh", rede%nodes
   write(*,*) ""
   write(*,*)"O grau medio dos sitios da rede eh : ", rede%degMean
   write(*,*) ""
   !====================================================================
   write(*,*) "Gerou a rede"
   write(*,*) ""
   call sub_classifica_clusters(rede,.False., 000, 'sem_arquivo.dat') 
   !====================================================================
   call aloca_listas_SIRS_rDMP(rede)
   !====================================================================
   !====================================================================
   ! A principio
   !====================================================================
   usouCopia = .False.

   call aloca_listas_e_matrizHashimoto(rede)

   call metodo_potencia_Hashimoto(rede, lamb0, mu)

   !--------------------------------------------------------------------
   if( rede%nodes == 10**7 )then
      !-----------------------------------------------------------------
      call calcula_distancias_ao_hub(rede)
      !-----------------------------------------------------------------      
      fi_char = trim(adjustl(local))//trim(adjustl('fi_rDMP_N_10_a_7_RRN_deg_6_Plus_Star_998_alp_TODOS'))
      !-----------------------------------------------------------------
      open(999, file = trim(adjustl(fi_char))//'.dat', status = 'unknown')
      do i1 = 1, rede%nodes
         write(999,*) rede%deg(i1), lista_distancias(i1), (fi_rDMP(i1))**2.0d0
         write(*,*) "i1 = ", i1, "fi = ", fi_rDMP(i1)
      enddo
      deallocate(lista_distancias)
      deallocate(fi_rDMP)
      !-----------------------------------------------------------------
      close(999)
      !-----------------------------------------------------------------
      call system('tar -czf '//trim(adjustl(fi_char))//'.tar.gz '//trim(adjustl(fi_char))//'.dat' )
      call system('rm '//trim(adjustl(fi_char))//'.dat' )
      !-----------------------------------------------------------------      
   endif
   !--------------------------------------------------------------------
   if( nargus == 6 )then
      !===============================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'.dat'                        

      inquire(unit=334, opened=taAberta)
      if(taAberta) close(334)             
      open(334, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
      write(334,*) lamb0, 0.0d0
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'.dat'

      inquire(unit=335, opened=taAberta)
      if(taAberta) close(335)
                         
      open(335, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
      !===============================================================   
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//'.dat'

      inquire(unit=337, opened=taAberta)
      if(taAberta) close(337)
                         
      open(337, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
      write(337,*) lamb0, Y4
      !=================================================================      
      lambdaf = lamb0 + dlamb * 500.0d0
      !=================================================================
   elseif( nargus == 7 )then
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
      lambdaf = lamb + 1.0d0 * ind_lamb * dlamb
      if(ind_lamb == 0) stop
   endif

   !stop "Nao vai rodar mais pontos alem do limiar. Venha e edite."
     !==================================================================
     ! Loop lambda
     !==================================================================
     
ll:     do while(lamb <= lambdaf)          
           !============================================================
           ! tempo
           !============================================================
           
           n_it = 0
           n_sub = 0
           !============================================================
           if( nargus == 7)then              

              inquire(unit = 333, opened = taAberta)
              if(taAberta)close(333)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
              !============================================================
              open(333, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================
           
              inquire(unit = 336, opened = taAberta)
              if(taAberta)close(336)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
              !============================================================
              open(336, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================

           elseif( nargus == 6)then
              inquire(unit = 333, opened = taAberta)
              if(taAberta)close(333)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_rDMP'))//'.dat'              
              !============================================================
              open(333, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================

              inquire(unit = 336, opened = taAberta)
              if(taAberta)close(336)
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'.dat'              
              !============================================================
              open(336, file=trim(adjustl(nomeArquivo)), status='unknown')   
              !============================================================
           endif           
           
           if( .not. usouCopia)then 
              call condicao_inicial_SIRS_rDMP(rede, p0)                 
              
              rho0 = p0
              rho = rho0
              
              Y4_old = 1.0_dp/(1.0_dp * comp_gigante)
              Y4 = Y4_old
              
              tempoEscrita = tempoEscrita0
              !=========================================================
           else   
              usouCopia = .False.
           endif
           
           !============================================================
lt:        do while(t <= tf)
              
              !call tempoDeEscrever()
              
              write(333, *) t, Y4
              write(336, *) t, rho
              
              !=========================================================
              call rk4_SIRS_rDMP_grafo(rede, t, dt, alp, lamb, mu)
              !=========================================================
         
              rho = 0.0_dp
         
              Ii_neg = .False.
              Ii_gt_1 = .False.
              do i1 = 1, size(Ii)   ! Calcula rho
                 if(lista_de_clusters(i1) /= i_comp_gigante) cycle
                 if( Ii(i1) < 0.0_dp)then
                    Ii_neg = .True.
                    exit
                 elseif( Ii(i1) > 1.0_dp )then
                    Ii_gt_1 = .True.
                    exit
                 endif
                 rho = rho + Ii(i1)
              enddo
         
              rho  = rho/( 1.0_dp * comp_gigante )
         
              sumIi2 = (sum(Ii**2.0_dp))**0.5_dp
         
		      Y4 = sum( (Ii/sumIi2)**4.0_dp )
           
              t = t + dt
           
              if( abs(rho - rho0)/rho0 < tol )then
                 n_it = n_it + 1
                 if( n_it >= int(2.0d0/dt) )then
                    write(*,*) "Conv. cm lbd = ", lamb, " e rho = ", rho, "."
                    exit lt
                 endif
              else
                 n_it = 0
              endif
              
              if( ( Ii_neg == .True. ) .or. ( Ii_gt_1 == .True. ) )then                 
                 write(*,*) "Dev. err. num., atlz. dt para ", dt/2.0d0            
                        
                 dt = dt/2.0d0
            
                 if( .not. usouCopia)then 
                    rho0 = p0
         
                    call condicao_inicial_SIRS_rDMP(rede, p0)
         
                    tempoEscrita = tempoEscrita0
                    !=========================================================
                 else   
                    write(*,*) "Atualizou dt e usou uma copia."
                    usouCopia = .False.
                 endif            
            
                 cycle lt
              endif   

              rho0 = rho
              Y4_old = Y4
           enddo lt
           
           
           close(333)
           close(336)
        
           if( n_it >= int(2.0d0/dt) )then
              write(334, *) lamb, rho           
              write(335, *) lamb, t
              
              write(337, *) lamb, Y4
           endif
           
           lamb = lamb + dlamb
           
           call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia*.tar.gz')) )
           call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia*.dat')) )
        enddo ll
        close(334)
        close(335)
        close(337)
     
   contains 
   !====================================================================
      subroutine criaRedeEClassificaClusters(this, tam_rede1, grauRRN1, grauStar1)
        type(grafoRRN_Plus_Star) :: this
        integer :: tam_rede1
        integer :: grauRRN1
        integer :: grauStar1
        
!#######################################################################   
        call this%ligaRRNStar(tam_rede1, grauRRN1, grauStar1)
!#######################################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
!#######################################################################
      end subroutine
!#######################################################################

      subroutine entradaArgumentos()
         nargus = iargc()

         if(nargus == 6)then
            write(*,*) "6 argumentos foram dados de entrada"

            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede
            !#############################
            !  Grau RRN
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauRRN
            write(char_grau_RRN, '(I0)') grauRRN            
            !############################# 
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauStar
            write(char_grauStar_RRN, '(I0)') grauStar
            !#############################
            ! Lambda0
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor

            !#############################
            ! Alfa
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
                              
         elseif(nargus == 7)then
            write(*,*) "Foram recebidos 7 argumentos"
            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede
            !#############################
            !  Grau RRN
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauRRN
            write(char_grau_RRN, '(I0)') grauRRN            
            !############################# 
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauStar
            write(char_grauStar_RRN, '(I0)') grauStar
            !#############################
            ! Lambda0
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor

            !#############################
            ! Alfa
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
            
            !#############################
            ! Ind Lambda
            !#############################
            
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_lamb   
            write(char_ind_lamb, '(I0)') ind_lamb
        
            char_ind_lamb = trim(adjustl(char_ind_lamb))
            
         else
            stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
         endif

         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!#######################################################################
!  Indice da semente
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
		 character(len = 100) :: C_lamb_f, C_tempo_f, C_R_i_f, C_I_ji_f, C_I_i_f
		 real(dp) :: t_lido
		 !--------------------------------------------------------------
		 pasta_fonte = trim(adjustl(pasta_fonte))
		 !--------------------------------------------------------------
		 if(nargus == 10 )then
			C_lamb_f = trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_tempo_f = trim(adjustl(C_tempo))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_R_i_f = trim(adjustl(C_R_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_I_ji_f = trim(adjustl(C_I_ji))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
			C_I_i_f = trim(adjustl(C_I_i))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
		 elseif( nargus == 9)then
			C_lamb_f = trim(adjustl(C_lamb))//'_global'
			C_tempo_f = trim(adjustl(C_tempo))//'_global'
			C_R_i_f = trim(adjustl(C_R_i))//'_global'
			C_I_ji_f = trim(adjustl(C_I_ji))//'_global'
			C_I_i_f = trim(adjustl(C_I_i))
		 endif
		 Ii = 0.0d0; Iji = 0.0d0; Ri = 0.0d0
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
		rho = 0.0d0
		do i1 = 1, size(Ii)
		   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		   rho = rho + Ii(i1)
		enddo
		rho = rho/(1.0d0 * comp_gigante)
		!---------------------------------------------------------------
		Y4 = 0.0d0
		!---------------------------------------------------------------    
		if(allocated(rho_i1)) deallocate(rho_i1)
		allocate(rho_i1(size(Ii)))
		!---------------------------------------------------------------    
		rho_i1 = 0.0d0
		!---------------------------------------------------------------
		do i1 = 1, size(Ii)
		   if( lista_de_clusters(i1) /= i_comp_gigante) cycle
		   do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
			  i21 = link_back(i12)
			  rho_i1(i1) = rho_i1(i1) + Iji(i21)
		   enddo
		enddo
		norm_rhoi1 = (sum(rho_i1**2.d0))**0.5d0
		rho_i1 = rho_i1/norm_rhoi1
		Y4 = sum(rho_i1**4.0d0)
		!===============================================================
		Y4_old = Y4
		rho0 = rho
		!===============================================================
		t = t_lido
		lamb = lbd_lido
        deallocate(rho_i1)
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
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_i))
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
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_R_i))
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
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_I_ji))
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
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_tempo))
            endif           
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) t
            close(777)
            !===========================================================

            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl(C_lamb))
            endif  
            !===========================================================       
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(777,*) lamb
            close(777)
            !===========================================================
            tempoEscrita = tempoEscrita + tempoEscrita0
            !===========================================================
         endif
   end subroutine


end program
