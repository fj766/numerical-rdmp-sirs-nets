module hashimoto
   use geraRede
   use mod_rndgen
   use types
   
   implicit none
   
   real(dp), allocatable :: x(:), y(:)
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
   
   end subroutine                       
         !#######################################################
                                                                     
         !########################################################

   subroutine metodo_potencia_Hashimoto(this, lambda, mu1)
      class(grafo), intent(in) :: this
      real(dp), intent(out) :: lambda
      real(dp), intent(in) :: mu1
	  real(dp) :: xp, yp
	  integer :: ipx, ipy
	  type(rndgen) :: geni
	  integer :: semente
	  real(dp) :: erro
	  real(dp) :: autovalor
	  integer :: i1, i2, i3
      integer (kind=8) :: i12, i23, j4
	  integer :: contador
	  real(dp), allocatable :: v1(:)
	  
	  semente = 956792214	
	  
	  call geni%init(semente)
	  
	  ipx = 1	  
      xp = 0.0d0
      
      do i12 = 1, size(x)
         x(i12) = 1.0_dp * this%deg( this%listAdj(i12) )
         if(abs(x(i12)) > abs(xp))then
            xp = x(i12)
            ipx = i12
         endif
      enddo
      x = x/xp	        

      contador = 0
	   
iter: do while(.True.)	
         !Aqui comeca o algoritmo
         ipy = 1
         yp = 0.0d0
         y = 0.0d0
         
sit_i1:	 do i1 = 1, this%nodes

viz_i2:     do i12 = this%aux(i1), this%aux(i1) + this%deg(i1) -1
                   
               ! link u-->v
               i2 = this%listAdj(i12)
viz_i3:        do i23 = this%aux(i2), this%aux(i2) + this%deg(i2) - 1 ! link v --> w
                      
                  i3 = this%listAdj(i23)
            
                  if(i3 == i1) cycle viz_i3
            
                  y(i12) = y(i12) + x(i23)
			          			          
			      if(abs(y(i12)) > abs(yp))then
			         yp = y(i12)
			         ipy = i12
			      endif	   		
               enddo viz_i3                   
	        enddo viz_i2
	     enddo sit_i1      
			         
         autovalor = y(ipx)
         !#######################
         ! Se yp = 0, para tudo!
         if(yp == 0.0d0)then
            do i12 = 1, size(x)   
			   x(i12) = 2.0d0 * geni%rnd()
			   if(abs(x(i12)) > abs(xp))then
			      xp = x(i12)
			      ipx = i12
	           endif
    	    enddo
		    x = x/xp
		    contador = 0
		    write(*,*) 'Autovalor nulo encontrado'
		    cycle iter			   
         endif
			   
			   
		 erro = abs(maxval(x - y/yp))
		 
		 x = y/yp
		 ipx = ipy
		 xp = 1.0d0
		 	   
		 if(erro < tol)then
		    lambda = mu1/autovalor
		    write(13,*) this%nodes, lambda
			write(*,*) 'Autovalor principal ', autovalor, ' encontrado com erro relativo ', erro, '. '
			
			write(*,*) " "
			write(*,*) 'Foram necessarios ', contador, ' passos.'
			write(*,*) " "
			if(allocated(v1)) deallocate(v1)
			allocate(v1(size(x)))
			
			v1 = x
			
			v1 = v1/(sum(v1**2.0_dp))**0.5_dp
			
			Y4 = sum(v1**4.0_dp)

			deallocate(v1)
			
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
	
	contains
	
	
		!#################################################
		!	Modelo SIRS
		!#################################################		

    subroutine aloca_listas_SIRS_rDMP(this)
       class(grafo) :: this
              
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
        	integer :: n1
        	integer (kind=8) :: n2
        	real(dp) :: alp, lam, mu		 
            
            n1 = this%nodes
            n2 = this%sumDeg
  		!####################################################################################
		
		
		!####################################################################################	
                     
!                    f_Ri(this, t, Ri, Ii, n_sitios, alf, mu)
		k1_Ri = dt * f_Ri(this, t, Ri, Ii, n1, alp, mu)
		
		!            f_Ii(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
		k1_Ii = dt * f_Ii(this, t, Ri, Ii, Iji, n1, n2, mu, lam)

!                     f_Iji(this, t, Ri, Iji, n_sitios, n_arestas, mu, lam)
		k1_Iji = dt * f_Iji(this, t, Ri, Iji, n1, n2, mu, lam)
		
		!####################################################################################
                     
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n_sitios, alp, mu)
		k2_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n1, alp, mu)
		
!		            f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, Iji + 0.5_dp * k1_Iji, n1, n2, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Iji + 0.5_dp * k1_Iji, n_sitios, n_arestas, mu, lam)
		k2_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Iji + 0.5_dp * k1_Iji, n1, n2, mu, lam)
       
        !####################################################################################
                     !f_Ri(this, t             , Ri                 , Ii                 , n_sitios, alf, mu)
!                    f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n_sitios, alp, mu)
		k3_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n1, alp, mu)		

!                    f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, Iji + 0.5_dp * k2_Iji, n1, n2, mu, lam)

!                     f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Iji + 0.5_dp * k2_Iji, n_sitios, n_arestas, mu, lam)
		k3_Iji = dt * f_Iji(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Iji + 0.5_dp * k2_Iji, n1, n2, mu, lam)		
       
       !####################################################################################                     
                    !f_Ri(this, t     , Ri        , Ii        , n_sitios, alf, mu)
!                    f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n_sitios, alp, mu)
		k4_Ri = dt * f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n1, alp, mu)		
		
!                    f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Ii = dt * f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, Iji + k3_Iji, n1, n2, mu, lam)

!                     f_Iji(this, t + dt, Ri + k3_Ri, Iji + k3_Iji, n_sitios, n_arestas, mu, lam)
		k4_Iji = dt * f_Iji(this, t + dt, Ri + k3_Ri, Iji + k3_Iji, n1, n2, mu, lam)		
		
		!####################################################################################
			
		Ri = Ri + (1.0_dp/6.0_dp) * (k1_Ri + 2.0_dp * k2_Ri + 2.0_dp * k3_Ri + k4_Ri)  		
  		
		Ii = Ii + (1.0_dp/6.0_dp) * (k1_Ii + 2.0_dp * k2_Ii + 2.0_dp * k3_Ii + k4_Ii)  		
		
		Iji = Iji + (1.0_dp/6.0_dp) * (k1_Iji + 2.0_dp * k2_Iji + 2.0_dp * k3_Iji + k4_Iji)  
						
  		!####################################################################################

  		        		    
    end subroutine		


		function f_Iji(this, t, Ri, Iji, n_sitios, n_arestas, mu, lam)
			!use types
			!use geraRede
			
			class(grafo) :: this			
			real(dp), intent(in) :: t, mu, lam
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			real(dp) :: Ri(n_sitios), Iji(n_arestas)
			real(dp) :: f_Iji(n_arestas)
			real(dp) :: sumIji
			!###############################
			integer :: l1, l2, l3
                        integer(kind=8) :: l12, l23
			!###############################

			! Sitio i
			do l1 = 1, this%nodes
			    if(lista_de_clusters(l1) /= i_comp_gigante) cycle
				
				do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
					
					f_Iji(l12) = -1.0_dp * mu * Iji(l12)
					
					l2 = this%listAdj(l12)
					
					! Vizinhos k, do sitio j. 
					
					sumIji = 0.0_dp
					
					do l23 = this%aux(l2), this%aux(l2) + this%deg(l2) - 1
						! O sitio i nao vale como vizinho de j aqui.
						l3 = this%listAdj(l23)
						if(l3 == l1) cycle
						sumIji = sumIji + Iji(l23)						
					enddo

					f_Iji(l12) = f_Iji(l12) + lam * (1.0d0 - Ri(l2) - Ii(l2) ) * sumIji
				enddo
			enddo
		end function

		!#################################################
		!	F_Si
		!#################################################		

		function f_Ri(this, t, Ri, Ii, n_sitios, alf, mu)
			
			class(grafo) :: this
			real(dp), intent(in) :: t, alf, mu
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			real(dp) :: Ri(n_sitios), Ii(n_sitios)
			real(dp) :: f_Ri(n_sitios)
			!###############################
			integer :: l1, l2
                        integer(kind=8) :: l12
			!###############################
			
			! Sitio i
			do l1 = 1, this%nodes
			    !#################################################
			    if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			    !#################################################
				f_Ri(l1) = mu * Ii(l1) - alf * Ri(l1)
				!#################################################
			enddo
		end function			
	
		!#################################################
		!	F_Ii
		!#################################################

		function f_Ii(this, t, Ri, Ii, Iji, n_sitios, n_arestas, mu, lam)
			
			class(grafo) :: this
			real(dp), intent(in) :: t, mu, lam
			integer :: n_sitios
			integer (kind=8) :: n_arestas
			real(dp) :: Ri(n_sitios), Ii(n_sitios), Iji(n_arestas)
			real(dp) :: f_Ii(n_sitios)
			real(dp) :: sumIji
			!###############################
			integer :: l1, l2
                        integer (kind=8) :: l12
			!###############################
			
			! Sitio i.
			do l1 = 1, this%nodes
				!#################################################
				if(lista_de_clusters(l1) /= i_comp_gigante) cycle
				!#################################################
				f_Ii(l1) = -1.0_dp * mu * Ii(l1)
				!#################################################
				sumIji = 0.0_dp
				!#################################################
				do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
					sumIji = sumIji +  Iji(l12)
				enddo
				!#################################################				
				f_Ii(l1) = f_Ii(l1) + lam * (1.0d0 - Ii(l1) - Ri(l1) ) * sumIji
			enddo
		end function
!########################################################
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
   character(len=300) :: cwd, resultados
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
!#######################################################################   
   seed1=947361823
!#######################################################################   
   call entradaArgumentos()
!=======================================================================
   
   resultados = 'Rst_rDMP_Corte_2sqrtN'
   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
       
   local = trim(adjustl(resultados))//'/'
!#######################################################################
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
   tole = 1d-5
   per_conv = int( 10.0_dp/dt )
!#######################################################################

!#######################################################################
   call ger_inic%init(seed1)
   i2 = 1
   do i1 = 1, 1000
      if(mod(i1,100) > 0) cycle
      seed(i2)  = ger_inic%int(100000000,999999999)
      write(*,*) i1, seed(i2)
      i2 = i2+1      
   enddo
   
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
   write(alp_char2, '(f9.3)') alp
   
   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//'/'))
      
   call system('mkdir -p '//trim(adjustl(local)) )

   !====================================================================
   call aloca_listas_SIRS_rDMP(rede)
   !====================================================================
   !====================================================================
   ! A principio
   !====================================================================
   usouCopia = .False.
   
   call testaSeTemCopiaEstadoSistema(usouCopia)
   
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
           if( nargus == 10)then              
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
           elseif( nargus == 9)then
              nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_rDMP'))//'.dat'              
           endif           
           
           !============================================================
           open(336, file=trim(adjustl(nomeArquivo)), status='unknown')   
           !============================================================
           
           if( .not. usouCopia)then 
              call condicao_inicial_SIRS_rDMP(rede, p0)                 
              rho0 = p0
              tempoEscrita = tempoEscrita0
              !=========================================================
           else   
              usouCopia = .False.
           endif
           
           !============================================================
lt:        do while(t <= tf)
              
              call tempoDeEscrever()
              
              write(336, *) t, rho
           
              !=========================================================
              call rk4_SIRS_rDMP_grafo(rede, t, dt, alp, lamb, mu)
              !=========================================================
              rho = 0.0d0
           
              do i1 = 1, rede%nodes
                 if(lista_de_clusters(i1) /= i_comp_gigante) cycle
                 rho = rho + Ii(i1)
              enddo
           
              rho = rho/(1.0d0 * comp_gigante)
           
              t = t + dt
           
              if( abs(rho - rho0)/rho0 < tol )then
                 n_it = n_it + 1
                 if( n_it >= int(5.0d0/dt) )then
                    write(*,*) "Convergiu com lambda = ", lamb, " e rho = ", rho, "."
                    exit lt
                 endif
              else
                 n_it = 0
              endif
           
           
!              if( ( abs(rho - rho0)/rho0 > 5.0d0 * tol ) .and. ( rho < 1.0_dp/(10.0_dp * rede%nodes) ) )then
!                 n_sub = n_sub + 1
!                 if ( n_sub >= int(5.0d0/dt) )then
!                    write(*,*) "Abs subcr com lambda = ", lamb, " e rho = ", rho, "."              
!                    exit lt
!                 endif
!              else
!                 n_sub = 0
!              endif
              
              if( rho < 0.0d0 .or. rho > 1.0d0)then                 
                 write(*,*) "Atualizou dt de dt = ", dt, " para dt = ", dt/2.0d0            
            
                 call TestaSeTemCopiaEstadoSistema(usouCopia)
            
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
           enddo lt
           
           
           close(336)
        
           if( n_it >= int(5.0d0/dt) )then
              write(334, *) lamb, rho           
              write(335, *) lamb, t

              if(allocated(v1)) deallocate(v1)
              allocate(v1(size(Ii)))
         
   		      v1 = Ii
			
		      v1 = v1/(sum(v1**2.0_dp))**0.5_dp
			
		      Y4 = sum(v1**4.0_dp)
              
              deallocate(v1)
              
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
      subroutine leAtehOfimUnidimensional(label, nomeArquivo, var1, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1
               !write(*,*) st
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine


   !====================================================================
      subroutine leAtehOfimBidimensional(label, nomeArquivo, var1, var2, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1, var2
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1, var2
               !write(*,*) st
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine

 
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

!#######################################################################
         if(gama_exp < 3.0_dp)then
            grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
         else
            !grau_max = (1.0_dp * tam_rede)**(0.5_dp)
            !grau_max = (1.0_dp * tam_rede)**(1d0/(gama_exp-1.d0))
            call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
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

   recursive subroutine TestaSeTemCopiaEstadoSistema(usouCopia)
      logical, intent(inout) :: usouCopia
      integer :: ind_anterior
      real(dp) :: rho1
      logical ::  existe2
      integer ::  st1, st2
      integer :: ind_arquivo
      logical :: taAberta
      integer :: i10
      !=================================================================
      ! A principio, dizemos que nao usara a Copia.
      ! Apenas com os testes abaixo isso mudara. Ou nao.
      ! UsouCopia eh a saida principal do programa
      !=================================================================

      usouCopia = .False.
      
      !=================================================================
      if(nargus == 10)then

         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))

         inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe ) 
                         
         if(existe)then
            usouCopia = .True.

            inquire(unit=776, opened=taAberta)
            if(taAberta) close(776)
                        
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            read(776,*, iostat=st) lamb0
            close(776)
            
            if( st /= 0)then
               usouCopia = .False.
               call system('rm -r '//trim(adjustl(local))//'Copia*')               
               write(*,*) "Rotina deletou Copia_lambda.dat e se chamou sozinha"
               call TestaSeTemCopiaEstadoSistema(usouCopia)
            endif
         else
            !===============================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))
            
            inquire(file=trim(adjustl(nomeArquivo))//'.dat', exist = existe2)
            teveLeitura = .False.
            !===============================================================                     
            if(existe2)then   
               inquire(unit=334, opened=taAberta)
               if(taAberta) close(334)
                                       
               open(334, file = trim(adjustl(nomeArquivo))//'.dat', status='old')
               read(334, *, iostat=st) lamb0, rho1
               close(334)
               
               if( st == 0 )then			                    
                  if( rho1 == 0.0_dp )then
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))
                     inquire(file=trim(adjustl(nomeArquivo))//'.dat', exist=existe2)
                     
                     st = -1
                     if(existe2)then
                        inquire(unit=444, opened=taAberta)
                        if(taAberta) close(444)                     
                        
                        open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
                        read(444, *, iostat=st) lamb0, Y4
                        if(st == 0) teveLeitura = .True.
                        close(444)
                     endif
                     
                     if( st /= 0 )then        
                        if(soCalculaIPR)then                           
                           call aloca_listas_e_matrizHashimoto(rede)
                           call metodo_potencia_Hashimoto(rede, lamb0, mu)
                           !=====================================================
                           deallocate(x)
                           deallocate(y)               
                           !=====================================================
		   	               open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			               write(444,*) lamb0, Y4
			               close(444)			     
			               stop
			            endif
			         endif			      
			      endif
			   endif
			   close(334)
			      
               !call leAtehOfimBidimensional(334, trim(adjustl(nomeArquivo))//'.dat', lamb0, rho1, teveLeitura)              
               !===============================================================
               !close(334)
               !===============================================================
            endif
            if(.not.teveLeitura)then
               !============================================================
               nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//trim(adjustl('_lamb_index_0'))
               
               !========================================================
               ! Abre arquivo que nao sabemos se ja existia
               ! e quem dira se tinha alguma coisa nele.
               !========================================================
               inquire(unit=334, opened=taAberta)
               if(taAberta) close(334)    
               
               open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')               
               !========================================================
               ! Damos voto de confianca
               !========================================================
               teveLeitura = .True.
               read(334, *, iostat=st2) lamb0, rho1
               if(st2 /=0) teveLeitura = .False. 
               close(334)              
               !========================================================
               if( .not. teveLeitura)then        
                  !=====================================================         
                  call aloca_listas_e_matrizHashimoto(rede)
                  call metodo_potencia_Hashimoto(rede, lamb0, mu)
                                    
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//trim(adjustl('_lamb_index_0'))
                  
                  inquire(unit=444, opened=taAberta)
                  if(taAberta) close(444)    
                  
			      open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			      write(444,*) lamb0, Y4
			      close(444)			      
                  !=====================================================           
                  deallocate(x)
                  deallocate(y)                 
                  !=====================================================
               !========================================================
                  inquire(unit=334, opened=taAberta)
                  
                  if(taAberta) close(334)    
                  
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//trim(adjustl('_lamb_index_0'))               
                  
                  open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')                  
                  
                  write(334,*) lamb0, 0.0d0
                  close(334)
                  !=====================================================            
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'_lamb_index_0' 
                  
                  inquire(unit=335, opened=taAberta)
                  if(taAberta) close(335)

                  open(335, file=trim(adjustl(nomeArquivo))//'.dat', status = 'unknown')
                  write(335,*) lamb0, tf
                  !===========================================================
                  close(335)
               else
                  if(soCalculaIPR)then
                     call aloca_listas_e_matrizHashimoto(rede)
                     call metodo_potencia_Hashimoto(rede, lamb0, mu)
                     !==================================================
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_rDMP'))//trim(adjustl('_lamb_index_0'))
                  
                     inquire(unit=444, opened=taAberta)
                     if(taAberta) close(444)    
                  
			         open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			         write(444,*) lamb0, Y4
   			         close(444)			      
                     !================================================== 

			         deallocate(x,y)
			         stop
			      endif   
               endif
               close(334)                             
            endif   
            lamb0 = lamb0 + (1.0_dp * ind_lamb) * dlamb
         endif
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'       
         !==============================================================
         open(334, file=trim(adjustl(nomeArquivo)), access='append', status = 'unknown')
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat' 
         open(335, file=trim(adjustl(nomeArquivo)), access='append', status = 'unknown')
         !==============================================================       
         !==============================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat' 
         open(337, file=trim(adjustl(nomeArquivo)), access='append', status = 'unknown')
         !==============================================================   


         lambdaf = lamb0        
      elseif(nargus == 9)then
         
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda.dat'))
 
         inquire( file=trim(adjustl(nomeArquivo)), exist=existe ) 

         if(existe)then
            usouCopia = .True.
            inquire(unit=776, opened=taAberta)
            if(taAberta) close(776)             
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            read(776,*, iostat=st) lamb0
            close(776)
            !===========================================================
            if( st1 /= 0)then
               call system('rm -r '//trim(adjustl(local))//'Copia*')
               write(*,*) "Rotina deletou Copia_lambda.dat e se chamou sozinha"
               call TestaSeTemCopiaEstadoSistema(usouCopia)
            endif  
            !=========================================================== 
         else
            
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP_lamb_index_0'))
            inquire( file = trim(adjustl(nomeArquivo))//'.dat', exist = existe2)
            
            !===========================================================
            ! Por padrao, essa variavel sera inicalizada assim
            ! A razao esta nas linhas onde lemos lbd_vs_rho_rDMP.dat.
            !===========================================================
            st2 = -1 
            
            if(existe2)then
               inquire(unit=600, opened=taAberta)
               if(taAberta) close(600)                         
               
               open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               read(600, *, iostat = st2) lambda_Ultimo_Index, rho1
               close(600)
               if( st2 == 0)then            
                  if( rho1 == 0.0_dp )then
                     if(soCalculaIPR)then
                        call aloca_listas_e_matrizHashimoto(rede)
                        call metodo_potencia_Hashimoto(rede, lambda_Ultimo_Index, mu)
                                    
                        nomeArquivo = trim(adjustl(local))//trim(adjustl('N_vs_Y4_rDMP'))
                        inquire(unit=444, opened=taAberta)
                        if(taAberta) close(444)             
			            open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			            write(444,*) lambda_Ultimo_Index, Y4
			            close(444)
			      
                        !=====================================================
                        deallocate(x)
                        deallocate(y)                 
                        !=====================================================
			            stop
			         endif
			      endif               
               
               
               
lams:             do i1 = 1, 10
                     write(buffer, '(I0)') i1
                     buffer = trim(adjustl(buffer))
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(buffer))

                     inquire(unit=600, opened=taAberta)
                     if(taAberta) close(600)             
                  
                     open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
                     read(600, *, iostat = st2) lambda_Ultimo_Index, rho1
                     close(600)
                     if(st2 /= 0)then
                        ind_anterior = i1 - 1
                        write(buffer, '(I0)') ind_anterior
                        buffer = trim(adjustl(buffer))
                        nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'_lamb_index_'//trim(adjustl(buffer))
                        
                        inquire(unit=600, opened=taAberta)
                        if(taAberta) close(600)             
                        
                        open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
                        read(600, *, iostat= st2) lambda_Ultimo_Index, rho1
                        close(600)                        
                        exit lams
                     endif                   
                  enddo lams
               endif
            endif      
            !=================================================================        
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'.dat'
                 
            call leAtehOfimBidimensional(334, nomeArquivo, lamb0, rho, teveLeitura)            
   
            if(.not. teveLeitura)then
               if(st2 /= 0)then              
                  call aloca_listas_e_matrizHashimoto(rede)
                  call metodo_potencia_Hashimoto(rede, lamb0, mu)                                    
                  !=====================================================
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('N_vs_Y4_rDMP'))
            
                  inquire(unit=444, opened=taAberta)
                  if(taAberta) close(444)             

			      open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			      write(444,*) rede%nodes, Y4
			      close(444)
                  !=====================================================           
                  deallocate(x)
                  deallocate(y)
                  !=====================================================        
                  lamb_C = lamb0
   
                  write(*,*) "lambda_C = ", lamb_C
         
                  write(334,*) lamb_C, 0.0d0

                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'.dat'

                  inquire(unit=335, opened=taAberta)
                  if(taAberta) close(335)             
              
                  open(335, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
              
                  write(335,*) lamb_C, tf
                  close(335)
               else
                  lamb0 = lambda_Ultimo_Index  
               endif
            else
               if(lambda_Ultimo_Index > lamb0) lamb0 = lambda_Ultimo_Index                                 
            endif
            lamb0 = lamb0 + dlamb
            !===============================================================
            close(334)      
            !===============================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_rDMP'))//'.dat'                        

            inquire(unit=334, opened=taAberta)
            if(taAberta) close(334)             
            open(334, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
            !=================================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_rDMP'))//'.dat'

            inquire(unit=335, opened=taAberta)
            if(taAberta) close(335)
                         
            open(335, file=trim(adjustl(nomeArquivo)), access='append', status='unknown')
            !===============================================================   
         endif     
      endif
!#######################################################################
   	  lamb = lamb0
!#######################################################################      

      if( ( .not. usouCopia) .or. ( .not. existe ) ) return
      
      if( nargus == 10)then
         !=========================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
                    
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )

         if(existe)then
            usouCopia = .True.
            call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            inquire(unit=777, opened=taAberta)
            if(taAberta) close(777)
                         
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            
            do i1 = 1, rede%nodes
               read(777, *, iostat = st) Ii(i1)
            enddo
            !===========================================================
            close(777)
            !===========================================================                             
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            if(st /= 0)then
               usouCopia = .False.
               return
            endif
            !===========================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
            
            if(existe)then
               !========================================================
               call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

               inquire(unit=778, opened=taAberta)
               if(taAberta) close(778)             
               open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               !========================================================
               
               do i1 = 1, rede%nodes
                  read(778, *, iostat=st) Ri(i1)
               enddo
            endif

            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif
            !======================================================
            close(778)
            !======================================================
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')      
            !========================================================           
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Iji'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )

            if(existe)then
               !========================================================
               call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

               inquire(unit=779, opened=taAberta)
               if(taAberta) close(779)             
               open(779, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               !========================================================
               do i12 = 1, rede%sumDeg
                  read(779, *, iostat=st) Iji(i12)
               enddo
            endif
            
            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif            
            !======================================================
            close(779)
            !======================================================
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')      
            !========================================================           
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))//'.dat'
            inquire( file=trim(adjustl(nomeArquivo)), exist=existe )

            !======================================================
            if(existe)then

               inquire(unit=780, opened=taAberta)
               if(taAberta) close(780)             
               
               open(780, file=trim(adjustl(nomeArquivo)), status='old')                 
               read(780, *, iostat=st) t, dt
            endif
            
            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif                             
            !======================================================
            close(780)
            !======================================================
         endif   
      elseif( nargus == 9)then
         !======================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))
           
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )

         if(existe)then            
            usouCopia = .True.
            
            call system('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            inquire(unit=777, opened=taAberta)
            if(taAberta) close(777)             

            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
                        
            do i1 = 1, rede%nodes
               read(777, *, iostat=st) Ii(i1)
            enddo
            !===========================================================
            close(777)
            !===========================================================         
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            if(st /= 0)then
               usouCopia = .False.
               return
            endif

            !===========================================================                 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))
            inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
         
            if(existe)then
               
               call system('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

               inquire(unit=778, opened=taAberta)
               if(taAberta) close(778)             
               
               open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               
               do i1 = 1, rede%nodes
                  read(778, *, iostat=st) Ri(i1)
               enddo
               close(778)
               call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            endif
            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif            
            !===========================================================                 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Iji'))
            inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
        

            if(existe)then
               
               call system('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

               inquire(unit=779, opened=taAberta)
               if(taAberta) close(779)             

               open(779, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               
               do i12 = 1, rede%sumDeg
                  read(779, *, iostat=st) Iji(i12)
               enddo
               close(779)
               
               call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
               
            endif

            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif
            !===========================================================
         
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t.dat'))
            inquire( file=trim(adjustl(nomeArquivo)), exist=existe )
         
            if(existe)then

               inquire(unit=780, opened=taAberta)
               if(taAberta) close(780)             

               open(780, file=trim(adjustl(nomeArquivo)), status='old')                 
               read(780, *, iostat = st) t, dt
               !======================================================
               close(780)
               !======================================================
            endif

            if( (st /= 0) .or. (.not.existe) )then
               usouCopia = .False.
               return
            endif
         endif           
      endif   
      tempoEscrita = t + tempoEscrita0
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
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))
            endif           
            
            !========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')            
            !========================================================               
            do i1 = 1, size(Ii)
               write(777,*) Ii(i1)
            enddo
            !========================================================
            close(777)                 
            !========================================================
            
            !===========================================================
            call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))
            endif            
            !========================================================
            open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !========================================================               
            do i1 = 1, size(Ri)
               write(778,*) Ri(i1)
            enddo
            !========================================================
            close(778)                 
            !========================================================
                 
            !======================================================
            call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !======================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Iji'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Iji'))
            endif           
               
            !===========================================================
            open(779, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !===========================================================               
            do i12 = 1, size(Iji)
               write(779,*) Iji(i12)
            enddo
            !===========================================================
            close(779)                                         
            !===========================================================
            
            call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !======================================================
                 
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))
            endif           
            !======================================================
            open(780, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(780,*) t, dt
            close(780)                 
            !======================================================
            tempoEscrita = tempoEscrita + tempoEscrita0
            !======================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))
            endif  
            !=========================================================         
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(776,*) lamb
            close(776)
         !=========================================================
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
