#ifndef LAND_SMOOTH
# define TIMING

      subroutine smooth_thread(Lm,Mm, msk, h, Lgh,hmin,hmax,r_max,
     &                                      method,Lgh1)

      implicit none
      integer Lm,Mm
      character(len=64) method
      character(len=64), dimension(4) :: log_list

      real(kind=8),    dimension(0:Lm+1,0:Mm+1) :: h,Lgh,Lgh1
      real(kind=8), allocatable, dimension(:) :: FX,FE, FX1,FE1
      real(kind=8) :: hmin,lcl_hmin, hmax, r_max

      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      integer istr,iend,jstr,jend, i,j, iter,ifrst_call
      integer numthreads,trd, nsub_x,nsub_y, isize,jsize,
     &               tile, my_first,my_last, range, size
      integer iters_cond, iters_lin
      integer(kind=2) :: transform_to_log

C$    integer omp_get_num_threads, omp_get_thread_num
# ifdef TIMING
      integer(kind=4) iclk_start, iclk_end, clk_rate, clk_max
# endif
      log_list = [character(len=64) :: 'LOG_SMOOTHING', 
     &               'LEGACY_LOG_SMOOTH','LOG_SMOOTH_2','LOG_SMOOTH_1']
      transform_to_log=0
      iters_cond=0
      iters_lin=0
      numthreads=1 ; trd=0
      ifrst_call=1
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      if (trd == 0) then
# ifdef TIMING
        call system_clock(iclk_start, clk_rate, clk_max)
# endif
        iters_cond=500 !<-- set number of iterations
        write(*,'(1x,2A,F10.7,3(2x,A,I5))')  'enter smooth_thread ',
     &                 'r_max =', r_max, 'iters_cond =', iters_cond,
     &                                    'iters_lin =', iters_lin
C$   &                                 , 'numthreads =', numthreads
      endif

      call set_tiles(Lm,Mm, nsub_x,nsub_y)

c**   nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      isize=(Lm+nsub_x-1)/nsub_x ;  jsize=(Mm+nsub_y-1)/nsub_y
      size=(isize+4)*(jsize+4)
      allocate(FX(size),FE(size), FX1(size),FE1(size))

      write(*,'(5(2x,A,I4))') 'nsub_x =',nsub_x, 'nsub_y =',nsub_y,
     &               'isize =',isize, 'jsize =',jsize, 'trd =',trd

      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

        if (istr == 1) istr=istr-1 ; if (iend == Lm) iend=iend+1
        if (jstr == 1) jstr=jstr-1 ; if (jend == Mm) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            h(i,j)=min( hmax, max( hmin, h(i,j) ))
            if ( hmin<=0 ) then
                  h(i,j)=h(i,j)-hmin+0.1  ! raise topo when hmin<=0 to avoid log error
                  lcl_hmin=0.1
              else
                  lcl_hmin=hmin
              endif
            
            if ( any( log_list == method ) ) then
              transform_to_log=1
              if (h(i,j) > lcl_hmin) then
                Lgh(i,j)=log(h(i,j)/lcl_hmin)
              else
                Lgh(i,j)=0.D0
              endif
            endif
          enddo
        enddo
      enddo
C$OMP BARRIER

! Conditional smoothing or log-smoothing to reduce r-factors

      do iter=1,iters_cond
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

          call rx_diag_tile(istr,iend,jstr,jend, Lm,Mm, h,msk,
     &                                 iter-1, nsub_x*nsub_y,ifrst_call)
          ifrst_call=0
          if ( method == 'LOG_SMOOTHING') then
              call lsmooth_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     &                           msk, Lgh,Lgh1, FX,FE,FX1,FE1)

          elseif (method == 'LOG_SMOOTH_2') then
            call lsmooth_2_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     &                                msk, Lgh,Lgh1, FX,FE,FE1)

          elseif (method == 'LOG_SMOOTH_1') then
            call lsmooth_1_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     &                                    msk, Lgh,Lgh1, FX,FE)

          elseif (method == 'LEGACY_LOG_SMOOTH') then
            call lsmth_legacy_tile(istr,iend,jstr,jend, Lm,Mm,
     &                  r_max,      Lgh,Lgh1, FX,FE,FX1,FE1)

          else
            call cond_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     &                  r_max, msk, h,  Lgh1, FX,FE,FX1,FE1)
          endif
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)
          if (method == 'LOG_SMOOTHING') then
            call land_lsmth_tile (istr,iend,jstr,jend, Lm,Mm,
     &                  r_max, msk, Lgh1,Lgh, FX,FE,FX1,FE1)

          elseif (method == 'LOG_SMOOTH_2') then
            call land_lsmth_2_tile(istr,iend,jstr,jend, Lm,Mm,
     &                      r_max, msk, Lgh,Lgh1, FX,FE,FE1)

          elseif (method == 'LOG_SMOOTH_1') then
            call land_lsmth_1_tile(istr,iend,jstr,jend, Lm,Mm,
     &                         r_max, msk, Lgh,Lgh1, FX,FE)

          elseif (method == 'LEGACY_LOG_SMOOTH') then
            call lsmth_legacy_tile(istr,iend,jstr,jend, Lm,Mm,
     &                r_max,        Lgh1,Lgh, FX,FE,FX1,FE1)

          else
            call land_cnd_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     &                   r_max, msk, Lgh1,h, FX,FE, FX1,FE1)
          endif

          if (transform_to_log == 1) then 
            if (istr == 1) istr=istr-1 ; if (iend == Lm) iend=iend+1
            if (jstr == 1) jstr=jstr-1 ; if (jend == Mm) jend=jend+1

            do j=jstr,jend                ! this backward conversion is
              do i=istr,iend              ! needed because "rx_diag_tile"
                h(i,j)=lcl_hmin*exp(Lgh(i,j)) ! expects "h" to report "r_max"
              enddo                       ! values achieved at every
            enddo                         ! iteration.
          endif

        enddo   !<-- tile
C$OMP BARRIER
      enddo  !<-- iter

      if ( hmin<=0 ) then
        h=h+hmin-0.1       ! Deepens topo if hmin<0 as we previously raise it   
      endif


! Final linear smoothing

      if (iters_lin > 0 .and. trd == 0) then
        write(*,*) 'Applying linear smoothing.'
      endif

      do iter=1,iters_lin
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                istr,iend,jstr,jend)

          call smooth_tile(istr,iend,jstr,jend, Lm,Mm, msk,
     &                                  h,Lgh1, FX,FE,FE1)
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                istr,iend,jstr,jend)

          call land_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     &                            msk, Lgh1,h, FX,FE,FE1)
        enddo
C$OMP BARRIER

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                istr,iend,jstr,jend)

          call rx_diag_tile(istr,iend,jstr,jend, Lm,Mm,
     &                         h,msk, iter, nsub_x*nsub_y)

        enddo
C$OMP BARRIER
      enddo  !<-- iter
# ifdef TIMING
      if (trd == 0) then
        call system_clock(iclk_end, clk_rate, clk_max)
        if (clk_rate > 0) then
          write(*,'(/1x,2A,F8.2,1x,A,I4,1x,A)')  'Wall Clock time ',
     &   'spent in smoothing', (iclk_end-iclk_start)/dble(clk_rate),
     &                         'sec running', numthreads, 'threads.'
        endif
      endif
# endif
      end
#endif /* LAND_SMOOTH */


