program sirsMessagePassing

!####################################################
	use mod_numerico
	use geraRede
	use mod_rndgen
	use types
!####################################################
! Modulos utilizados
	
	implicit none

!####################################################
	real(dp) :: t0, tf, t, dt
!####################################################
! Tempo discretizado e incremento

!########################################################	
	real(dp), allocatable :: Si(:), Ii(:), Ri(:), Iji(:)
	real(dp) :: Si_medio, Ii_medio, Ri_medio, I_jimedio
	real(dp) :: Iji_m, Si_m, Ii_m, Ri_m
	integer ::  n_arestas, n_pts, semente
!########################################################
! Variaveis dinamicas

!########################################################
	real(dp) :: gam, lam, lam0, lamf, dlam, rho
!########################################################
! Taxas de transicao

!########################################################
	type(grafo_PL_UCM) :: rede
!########################################################
! Declara variaveis tipo grafo


!########################################################
	type(rndgen) :: gen
!########################################################
! Gerador de numeros pseudo-aleatorios

!########################################################
	real(dp), parameter :: p0 = 1.d0
	real(dp) :: p_sort
!########################################################
! Probabilidade inicial e numero do sorteio

!########################################################
	integer :: i1, i2, i3, j1, j2, j3, j4
!########################################################
! Variaveis auxiliares	

!########################################################
	character(len=30) :: rowfmt
!########################################################
! Tamanho do espaco requerido para printar cada entrada
! no arquivo

!########################################################
	integer, parameter :: sit_escol = 50
	integer, parameter :: are_escol = 20
	integer, parameter :: n_lam = 5000
	integer :: n_sitios(6)
	real(dp) :: exp_gama(3)
!########################################################
! Como sou estravagante, escolho 100 sitios para plotar
! seus estados

!########################################################
	character(len=50) :: gama_Char, lam_Char, N_char, rede_info, are_escol_Char, formato_padrao
	character(len=1000) :: caminho 
	character(len=2000) :: meu_formato
!########################################################

	!#########################################################################################
	write(are_escol_Char,*) are_escol+1
	are_escol_Char = trim(adjustl(are_escol_Char))
	
	meu_formato = trim(adjustl(""))
	
	formato_padrao = trim(adjustl("F9.7, 1X"))
	
	
	do j1 = 1, are_escol + 1
		meu_formato = trim(adjustl(meu_formato//are_escol_Char//formato_padrao))
	enddo
	
	meu_formato = trim(adjustl(meu_formato))
	!#########################################################################################
	
	caminho = trim(adjustl("/home/jota/S_comp/S_num/SIRS_MP/Rst_vs_t/"))

	
	t0 = 0.0_dp;	tf = 50.0_dp
	
	n_pts = 10000
	
	dt = 1.0_dp * (tf - t0)/n_pts
	
	lam0 = 0.120_dp; lamf = 1.5_dp;
	
	rho = 1.0_dp; gam = 0.5_dp
	
	lam = lam0
	
	dlam = 1.0_dp * (lamf - lam0)/n_lam	
	
	exp_gama(1) = 2.3d0; exp_gama(2) = 2.7d0
	exp_gama(3) = 3.5d0
	
	n_sitios(1) = 1000; n_sitios(2) = 3000;
	n_sitios(3) = 10000; n_sitios(4) = 30000;
	n_sitios(5) = 100000; n_sitios(6) = 300000;
				
!########################################################
! Define parametros


!########################################################
	semente = 995887671


lg: do j2 = 1, 1!size(exp_gama)


ls:  do j3 = 1, 2!size(n_sitios)
 
	call rede%iniciaGrafo(n_sitios(j3))
	
	call rede%inicia(2, sqrt(1.0d0 * n_sitios(j3)), exp_gama(j2), semente)
	
	call rede%liga(semente, .True.)
		
	rede%degMin = minval(rede%deg)		
	rede%degMax = maxval(rede%deg)
	
	n_sitios(j3) = rede%nodes
	
	n_arestas = sum(rede%deg)

	
!	write(*,*) n_arestas
!########################################################
! Gera o grafo que vamos usar	

	
!########################################################
	if(allocated(Iji)) deallocate(Iji)
		allocate(Iji(n_arestas))

	if(allocated(Si)) deallocate(Si)
		allocate(Si(n_sitios(j3)))		

	if(allocated(Ii)) deallocate(Ii)
		allocate(Ii(n_sitios(j3)))
		
	if(allocated(Ri)) deallocate(Ri)
		allocate(Ri(n_sitios(j3)))						
!########################################################
! Aloca as variaveis dinâmicas	
			
!########################################################


!########################################################


!########################################################
! Medias sao inicialmente nulas

!########################################################
ll:    do j4 = 1, n_lam   

	Iji_m = 0.0d0
	Si_m = 0.0d0
	Ii_m = 0.0d0
	Ri_m = 0.0d0


!########################################################
! Abre arquivos importantes

	write(lam_Char, '(F7.5)') lam
	lam_char = trim(adjustl(lam_Char))

	write(N_char, '(I0)') j3
	
	N_char = trim(adjustl(N_char))
	
       write(gama_Char, '(I0)') j2
       gama_char = trim(adjustl(gama_char))


	
	rede_info = trim(adjustl('_l_'//lam_char//'_N_'//N_char//'_g_'//gama_char//'_'))


	open(unit=10, file=trim(adjustl(trim(adjustl(caminho))//'I_ij_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=11, file=trim(adjustl(trim(adjustl(caminho))//'S_i_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=12, file=trim(adjustl(trim(adjustl(caminho))//'I_i_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=13, file=trim(adjustl(trim(adjustl(caminho))//'R_i_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	
	open(unit=14, file=trim(adjustl(trim(adjustl(caminho))//'I_ij_M_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=15, file=trim(adjustl(trim(adjustl(caminho))//'S_i_M_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=16, file=trim(adjustl(trim(adjustl(caminho))//'I_i_M_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')
	open(unit=17, file=trim(adjustl(trim(adjustl(caminho))//'R_i_M_Vs_T'//trim(adjustl(rede_info))//'.dat')), status='unknown')


!########################################################

	t = 0.0d0
	
	Iji = 0.0d0
	Si = 1.0d0
	Ii = 0.0d0
	Ri = 0.0d0
	
!	write(*,*) "O numero de sitios suscetiveis eh ", sum(Si)
		
	!################################################
	! Porcentagem inicial de infectados
	!################################################
			
	call gen%init(semente)
	
!	write(*,*) "O percentual de sitios inicialmente infectados eh ", p0
	
	do i1 = 1, rede%nodes
		p_sort = gen%rnd()

		if(p_sort > p0)cycle

		Ii(i1) = 1.0d0
		Si(i1) = 0.0d0
	 				
	enddo
	
!	write(*,*) "O tamanho da rede eh ", rede%nodes, "o esperado eh ", n_sitios	
	
!	write(*,*) "O numero de sitios infectados eh ", sum(Ii)
!	write(*,*) "O numero de sitios suscetiveis eh ", sum(Si)
!	write(*,*) "O numero de sitios recuperados eh ", sum(Ri)		
	
	do i1 = 1, rede%nodes
		do i2 = rede%aux(i1), rede%aux(i1) + rede%deg(i1) - 1
			j1 = rede%listAdj(i2)
			if(Ii(j1) == 1.0d0) Iji(i2) = 1.0d0
		enddo
	enddo
	
!	write(*,*) "O numero de arestas eh: ", sum(rede%deg)
!	write(*,*) "O numero de arestas infectadas eh: ", sum(Iji)
		
!########################################################
! Setta as condicoes iniciais







lpts:	do i1 = 1, n_pts
		
		!########################################################
		
		Ri = 1.0d0 -1.0d0 * (Si + Ii)
		
		Iji_m = 1.0d0 * sum(Iji)/n_arestas
		Si_m = 1.0d0 * sum(Si)/n_sitios(j3)
		Ii_m = 1.0d0 * sum(Ii)/n_sitios(j3)
		Ri_m = 1.0d0 * sum(Ri)/n_sitios(j3)
		!########################################################
		! Calcula as medias das variaveis dinamicas
		
		!##############################################################

		
		write(10, *) t, (Iji(j1), j1 = 1, are_escol)
		
		write(11, *) t, (Si(j1), j1 = 1, sit_escol)
		write(12, *) t, (Ii(j1), j1 = 1, sit_escol)
		write(13, *) t, (Ri(j1), j1 = 1, sit_escol)
		
		write(14, *) t, Iji_m
		write(15, *) t, Si_m
		write(16, *) t, Ii_m
		write(17, *) t, Ri_m
		
		
		!##############################################################		
		! Escreve nos arquivos.
		
		!##############################################################$$$$$$$$$$$$$$$$
		call rk4_grafo(rede, t, dt, n_arestas, n_sitios(j3), gam, lam, rho, Iji, Si, Ii, f_Iji, f_Si, f_Ii)	
		!##############################################################$$$$$$$$$$$$$$$$
		! Chama o metodo de Runge-Kutta de ordem 4
		
		
		!##############################################################
		t = t + dt
		!##############################################################
		! Incrementa o tempo.
	enddo lpts
	
!########################################################
! Rotina principal


!########################################################
	lam = lam + dlam

	close(10)
	close(11)
	close(12)
	close(13)

	close(14)
	close(15)
	close(16)
	close(17)
	
    enddo ll    
  enddo ls
 enddo lg
 		
!########################################################
! Fecha os arquivos


!########################################################
	stop
!########################################################
! Fecha o programa


!########################################################	
	contains
	
	
		!#################################################
		!	F_Iji
		!#################################################	

		function f_Iji(this, t, Iji, Si, n_arestas, n_sitios, rho, lam)
			use types
			use geraRede
			
			class(grafo) :: this			
			real(dp), intent(in) :: t, rho, lam
			integer :: n_arestas, n_sitios
			real(dp) :: Iji(n_arestas), Si(n_sitios)
			real(dp) :: f_Iji(n_arestas)
			real(dp) :: sumIji
			!###############################
			integer :: i1, i2, i3, j1
			!###############################

			! Sitio i
			do i1 = 1, this%nodes
				! Vizinhos j, do sitio i
				do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
					
					j1 = this%listAdj(i2)
					
					! Vizinhos k, do sitio j. 
					
					sumIji = 0.0d0
					
					do i3 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
						! O sitio i nao vale como vizinho de j aqui.
						if(this%listAdj(i3) /= i1)then
							sumIji = sumIji + Iji(i3)
						endif
					enddo	
					f_Iji(i2) = -1.0d0 * rho * Iji(i2)

					f_Iji(i2) = f_Iji(i2) + lam * Si(j1) * sumIji
				enddo
			enddo
			return									
		end function

		!#################################################
		!	F_Si
		!#################################################		

		function f_Si(this, t, Iji, Si, Ii, n_arestas, n_sitios, gam, lam)
			use types
			use geraRede
			
			class(grafo) :: this
			real(dp), intent(in) :: t, gam, lam
			integer :: n_arestas, n_sitios
			real(dp) :: Iji(n_arestas), Si(n_sitios), Ii(n_sitios)
			real(dp) :: f_Si(n_sitios)
			real(dp) :: sumIji
			!###############################
			integer :: i1, i2, j1
			!###############################
			
			! Sitio i
			do i1 = 1, this%nodes				
				sumIji = 0.0d0
				! Vizinhos do sitio i.
				do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
					sumIji = sumIji + Iji(i2)
				enddo
				f_Si(i1) = gam * (1.0d0 -Si(i1) -Ii(i1)) 				
				f_Si(i1) = f_Si(i1) -lam * Si(i1) * sumIji
			enddo
			return
		end function			
	
		!#################################################
		!	F_Ii
		!#################################################

		function f_Ii(this, t, Iji, Si, Ii, n_arestas, n_sitios, rho, lam)
			use types
			use geraRede
			
			class(grafo) :: this
			real(dp), intent(in) :: t, rho, lam
			integer :: n_arestas, n_sitios
			real(dp) :: Iji(n_arestas), Si(n_sitios), Ii(n_sitios)
			real(dp) :: f_Ii(n)
			real(dp) :: sumIji
			!###############################
			integer :: i1, i2, j1
			!###############################
			
			! Sitio i.
			do i1 = 1, this%nodes
				
				! Vizinhos do sitio i.
				sumIji = 0.0d0
				do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
					sumIji = sumIji +  Iji(i2)
				enddo
				f_Ii(i1) = -1.0d0 * rho * Ii(i1)
				
				f_Ii(i1) = f_Ii(i1) + lam * Si(i1) * sumIji
			enddo
			return
		end function
!########################################################
! Define funcoes

end program
