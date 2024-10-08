      subroutine cond_rx0_thread(ncx,ncy, msk, h,hmin, hmax, r_max,wrk)

! Conditions topography to restrict the maximum r-factor via
! deepening-only procedure: scan the grid to find points where
! pair-wise ratio
!                       | h(i+1) - h(i) |
!                  rx= -------------------
!                         h(i+1) + h(i)
!
! exceeds threshold "r_max" and adjust THE SMALLER  value of
! depth toward THE LARGER; iterate until the condition is no
! longer exceeded anywhere.

      implicit none
      integer ncx,ncy
      integer(kind=2) msk(ncx,ncy)
      real(kind=8) h(ncx,ncy), wrk(ncx,ncy)
      integer :: trd_count=0, npts_all=1
      integer :: iters_cond, iters_lin
      real(kind=8) :: hmin, hmax, r_max
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     &                    range, istr,iend,jstr,jend, iter, my_npts

      iters_cond=0
      iters_lin=0
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      
      call set_tiles(ncx,ncy, nsub_x,nsub_y)

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      iter=0
      do while(npts_all > 0)
        iter=iter+1 ; my_npts=0

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     &                           msk, h,wrk, r_max, my_npts)
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                   istr,iend,jstr,jend)
          call cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     &                           msk, wrk,h, r_max, my_npts)
        enddo
C$OMP CRITICAL(cond_rx0_smth)
        if (trd_count == 0) npts_all=0
        trd_count=trd_count+1 ; npts_all=npts_all+my_npts
        if (trd_count == numthreads) then
          trd_count=0
          write(*,'(8x,A,I7,2x,A,I10)') 'iter =', iter,
     &                           'changes =', npts_all
        endif
C$OMP END CRITICAL(cond_rx0_smth)
C$OMP BARRIER
      enddo !<-- while
      end


      subroutine cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     &                       msk, src,targ, r_max, my_npts)

      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, my_npts
      integer(kind=2) msk(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy), r_max
      integer i,j, iw,js,ie,jn
      real(kind=8) ratio,  max_surr
      real(kind=8), parameter :: epsil=1.D-12
      ratio=(1.-r_max)/(1.+r_max)

      do j=jstr,jend
        js=max(j-1,1) ; jn=min(j+1,ncy)
        do i=istr,iend
          iw=max(i-1,1) ; ie=min(i+1,ncx)
          targ(i,j)=src(i,j)      !<-- first copy then adjust

          if (msk(i,j) > 0) then
            max_surr=0.    !<-- initialize
            if (msk(iw,j) > 0) max_surr=max(max_surr, src(iw,j))
            if (msk(i,js) > 0) max_surr=max(max_surr, src(i,js))
            if (msk(ie,j) > 0) max_surr=max(max_surr, src(ie,j))
            if (msk(i,jn) > 0) max_surr=max(max_surr, src(i,jn))

            max_surr=max_surr*ratio
            if (targ(i,j) < max_surr-epsil) then
              targ(i,j)=max_surr
              my_npts=my_npts+1
            endif
          endif  !<-- msk(i,j)
        enddo  !<-- i
      enddo  !<-- j
      end
