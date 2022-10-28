program main
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   character(len = 1000) :: cw
   character(len = 1000) :: dir_parent
   character(len = 1000) :: arq_fonte0, arq_fonte1
   integer :: i1, i2, i3
   integer :: io_stat, io_stat2 
   character(len = 2) :: indice
   character(len = 5) :: alfa_char, gama_char
   real(dp), allocatable :: Y4_vs_tam(:), Y4(:)
   real(dp), allocatable :: lbd_vs_tam(:), lbd(:)
   real(dp) :: Y4_mean, Y4_dev, Y4_lido
   real(dp) :: lbd_mean, lbd_dev, lbd_lido
   integer :: n_Y4_lidos(9), n_lbd_lidos(9)
   integer :: tam(9)
   character(len=5) :: tam_char(9)
   integer :: tam0
   logical :: existe
   character(len=10) :: tipoCorte
   integer :: n_ams, n_tam

   !--------------------------------------------------------------------
   alfa_char = '0.100'
   gama_char = '2.30'
   !tipoCorte = '_sqrtN'
   tipoCorte = '_2sqrtN'
   !--------------------------------------------------------------------
   call getcwd(cw)
   dir_parent = trim(adjustl(cw))//trim(adjustl('/Rst_rDMP_Corte'))//trim(adjustl(tipoCorte))
   !--------------------------------------------------------------------
   n_ams = 50
   !--------------------------------------------------------------------
   n_tam = 9
   !--------------------------------------------------------------------   
   write(*,*) trim(adjustl(cw))
   !--------------------------------------------------------------------
   if(allocated(Y4)) deallocate(Y4)
   allocate(Y4(n_ams))
   !--------------------------------------------------------------------
   if(allocated(Y4_vs_tam)) deallocate(Y4_vs_tam)
   allocate(Y4_vs_tam(n_tam))
   !--------------------------------------------------------------------
   if(allocated(lbd)) deallocate(lbd)
   allocate(lbd(n_ams))
   !--------------------------------------------------------------------
   if(allocated(lbd_vs_tam)) deallocate(lbd_vs_tam)
   allocate(lbd_vs_tam(n_tam))

   !--------------------------------------------------------------------
   ! Lista tam nao eh allocatable.
   tam(1) = 10**3; tam(2) = 3*10**3; tam(3) = 10**4; 
   tam(4) = 3*10**4; tam(5) = 1*10**5; tam(6) = 3*10**5;
   tam(7) = 1*10**6; tam(8) = 3*10**6; tam(9) = 10**7;
   !--------------------------------------------------------------------
   ! Lista tam nao eh allocatable.
   tam_char(1) = '1k'; tam_char(2) = '3k'; tam_char(3) = '10k'; 
   tam_char(4) = '30k'; tam_char(5) = '100k'; tam_char(6) = '300k';
   tam_char(7) = '1M'; tam_char(8) = '3M'; tam_char(9) = '10M';   
   !--------------------------------------------------------------------
   n_Y4_lidos = 0
   n_lbd_lidos = 0
   arq_fonte1 = trim(adjustl(cw))//'/N_vs_Y4_DevY4_rDMP_Alfa_TODOS.dat'
   open(777, file = trim(adjustl(arq_fonte1)), access = 'append', status = 'unknown')

   arq_fonte1 = trim(adjustl(cw))//'/N_vs_lbd_DevLbd_rDMP_Alfa_TODOS.dat'
   open(778, file = trim(adjustl(arq_fonte1)), access = 'append', status = 'unknown')   
   !--------------------------------------------------------------------
   do i1 = 1, size(tam)
      !-----------------------------------------------------------------
      arq_fonte0 = trim(adjustl(dir_parent))//'/tam_'//trim(adjustl(tam_char(i1)))//'/gam_'//trim(adjustl(gama_char))
      !-----------------------------------------------------------------
      Y4_mean = 0.0d0
      lbd_mean = 0.0d0
      !-----------------------------------------------------------------
      Y4 = 0.0d0
      lbd = 0.0d0
      lams: do i2 = 1, n_ams
         !--------------------------------------------------------------      
         write(indice, '(I0)') i2
         !--------------------------------------------------------------
         arq_fonte1 = trim(adjustl(arq_fonte0))//'/ams_'//trim(adjustl(indice))//'/alp_'//trim(adjustl(alfa_char))//trim(adjustl(tipoCorte))
         !--------------------------------------------------------------
         inquire(file = trim(adjustl(arq_fonte1))//'/N_vs_Y4_rDMP_Teorico.dat', exist = existe)
         !--------------------------------------------------------------
         if(.not. existe)then
            !-----------------------------------------------------------
            write(*,*) "-----------------------------------------------------------"
            write(*,*) trim(adjustl(arq_fonte1))//'/N_vs_Y4_rDMP_Teorico.dat', 'nao existe'
         else
            open(111, file = trim(adjustl(arq_fonte1))//'/N_vs_Y4_rDMP_Teorico.dat', status = 'old')
               read(111,*, iostat = io_stat) tam0, Y4_lido            
            !-----------------------------------------------------------
            if( io_stat == 0) then
               !--------------------------------------------------------
               Y4(i2) = Y4_lido
               Y4_mean = Y4_mean + Y4(i2)
               n_Y4_lidos(i1) = n_Y4_lidos(i1) + 1
               !--------------------------------------------------------
            endif
            close(111)
         endif
         !--------------------------------------------------------------
         inquire(file = trim(adjustl(arq_fonte1))//'/N_vs_1_sobre_Lambda_H_Teorico.dat', exist = existe)
         !--------------------------------------------------------------

         if( .not. existe)then
            !-----------------------------------------------------------
            write(*,*) "-----------------------------------------------------------"
            write(*,*) trim(adjustl(arq_fonte1))//'/N_vs_1_sobre_Lambda_H_Teorico.dat', 'nao existe'         
            !-----------------------------------------------------------
         else
            open(112, file = trim(adjustl(arq_fonte1))//'/N_vs_1_sobre_Lambda_H_Teorico.dat', status = 'old')
               read(112,*, iostat = io_stat2) tam0, lbd_lido
            
            if(io_stat2 == 0)then
               lbd(i2) = lbd_lido
               lbd_mean = lbd_mean + lbd(i2)
               n_lbd_lidos(i1) = n_lbd_lidos(i1) + 1
            endif
            close(112)
         endif
         !--------------------------------------------------------------         
      enddo lams
      if( n_Y4_lidos(i1) > 1 )then
         !--------------------------------------------------------------
         Y4_mean = Y4_mean/(1.0d0 * n_Y4_lidos(i1))
         Y4_dev = 0.0d0
         !--------------------------------------------------------------
         do i3 = 1, size(Y4)
            if( Y4(i3) > 0.0d0) Y4_dev = Y4_dev + (Y4(i3) - Y4_mean)**2.0d0
         enddo
         !--------------------------------------------------------------
         Y4_dev = ( Y4_dev/( 1.0d0 * n_Y4_lidos(i1) -1.0d0 ) )**0.5d0
         !--------------------------------------------------------------
         write(777,*) tam(i1), Y4_mean, Y4_dev
         !--------------------------------------------------------------
      else
         write(777,*) tam(i1), Y4_mean, 0.0d0
      endif    
      !-----------------------------------------------------------------
      if( n_lbd_lidos(i1) > 1 )then
         !--------------------------------------------------------------
         lbd_mean = lbd_mean/(1.0d0 * n_lbd_lidos(i1))
         lbd_dev = 0.0d0
         !--------------------------------------------------------------
         do i3 = 1, size(lbd)
            if( lbd(i3) > 0.0d0) lbd_dev = lbd_dev + (lbd(i3) - lbd_mean)**2.0d0
         enddo
         !--------------------------------------------------------------
         lbd_dev = ( lbd_dev/( 1.0d0 * n_lbd_lidos(i1) -1.0d0 ) )**0.5d0
         !--------------------------------------------------------------
         write(778,*) tam(i1), lbd_mean, lbd_dev
         !--------------------------------------------------------------
      else
         write(778,*) tam(i1), lbd_mean, 0.0d0
      endif          
   enddo
   !--------------------------------------------------------------------
   close(777)
   close(778)
   stop
   !--------------------------------------------------------------------   
end program
