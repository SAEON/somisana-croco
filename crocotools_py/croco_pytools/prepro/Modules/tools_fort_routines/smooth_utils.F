      subroutine rx_diag_tile(istr,iend,jstr,jend, Lm,Mm, h,msk,
     &                                        iter,ntiles,ifrst_call)

      implicit none

      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h
      real(kind=8)  cff, my_rx,   my_ry
      real(kind=8) :: rx_max=0.D0, ry_max=0.D0

      integer istr,iend,jstr,jend, Lm,Mm, iter,ntiles
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      integer i,j, my_i_rx, my_i_ry, my_j_rx, my_j_ry
      integer :: trd_count=0, ifrst_call,
     &    irx_max=0, jrx_max=0, iry_max=0, jry_max=0

                                          ! This part is purely for
      my_rx=0.D0 ; my_i_rx=0 ; my_j_rx=0  ! diagnostics and causes no
      my_ry=0.D0 ; my_i_ry=0 ; my_j_ry=0  ! effect onto the resultant
                                          ! topography.
      do j=jstr,jend+1
        do i=istr,iend+1
          if (msk(i,j)> 0 .and. msk(i-1,j) > 0) then
            cff=abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
            if (cff > my_rx) then
              my_rx=cff ; my_i_rx=i ; my_j_rx=j
            endif
          endif
          if (msk(i,j) > 0 .and. msk(i,j-1) > 0) then
            cff=abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
            if (cff > my_ry) then
              my_ry=cff ; my_i_ry=i ;  my_j_ry=j
            endif
          endif
        enddo
      enddo
!      write(*,*) 'my_rx:',my_rx,'/ my_ry:',my_ry
C$OMP CRITICAL(rxdiag_cr_rgn)
      if (ifrst_call == 1) then
        ifrst_call=0
        write(*,'(/1x,A,2x,A,1x,A,14x,A,1x,A/1x,2A)')  'iter',
     &          'i,jrx_max', 'rx_max',  'i,jry_max', 'ry_max',
     &                    '---------------------------------',
     &                    '---------------------------------'
      endif
      trd_count=trd_count+1
      if (my_rx > rx_max) then
        rx_max=my_rx
        irx_max=my_i_rx ; jrx_max=my_j_rx
      endif
      if (my_ry > ry_max) then
        ry_max=my_ry
        iry_max=my_i_ry ; jry_max=my_j_ry
      endif
      if (trd_count == ntiles) then
        trd_count=0
        write(*,'(I5,2(2x,I4,I5,1x,A,F10.7))') iter, irx_max,jrx_max,
     &        'rx_max =',rx_max,  iry_max,jry_max, 'ry_max =',ry_max
        rx_max=0.D0
        ry_max=0.D0
      endif
C$OMP END CRITICAL(rxdiag_cr_rgn)
      end



      subroutine neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,A)
      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      real(kind=8) A(0:Lm+1,0:Mm+1)
      if (istr == 1) then
        do j=jstr,jend
          A(istr-1,j)=A(istr,j)
        enddo
      endif
      if (iend == Lm) then
        do j=jstr,jend
          A(iend+1,j)=A(iend,j)
        enddo
      endif
      if (jstr == 1) then
        do i=istr,iend
          A(i,jstr-1)=A(i,jstr)
        enddo
      endif
      if (jend == Mm) then
        do i=istr,iend
          A(i,jend+1)=A(i,jend)
        enddo
      endif
      if (istr == 1 .and. jstr == 1) then
        A(istr-1,jstr-1)=A(istr,jstr)
      endif
      if (istr == 1 .and. jend == Mm) then
        A(istr-1,jend+1)=A(istr,jend)
      endif
      if (iend == Lm .and. jstr == 1) then
        A(iend+1,jstr-1)=A(iend,jstr)
      endif
      if (iend == Lm .and. jend == Mm) then
        A(iend+1,jend+1)=A(iend,jend)
      endif
      end
