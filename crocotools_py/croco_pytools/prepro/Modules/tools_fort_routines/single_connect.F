      subroutine single_connect(Lm,Mm,i0,j0,rmask,outmask)

! Purpose: to enforce single-domain connectedness of unmasked area.
!-------------------------------------------------------------------
! The problem to be addressed is as follows: the landmask is generated
! from topography by declaring that land is where depth is less than
! user-specified minimum depth, this may leave holes in the landmask
! -- lakes which are not connected to the main body of water.
! This program takes a user-specified initial point, indices (i0,j0),
! which is assumed to belong to the main body of water and, starting
! from this point starts filling the adjacent non-masked points as
! "water" (so for each non-masked point on the grid it looks at all the
! immediate neighboring points and if at least one of them is already
! labelled as "water", then the point itself becomes "water").
! Note that at any time during the process there are three types of
! points: 1-masked ("land" according to the initial landmask mask),
! which remain unchanged; 2-unmasked, but non-labelled as "water",
! (initially all unmasked points); and 3-unmasked labelled as "water"
! (initially there is only one such point, i0,j0). The labelling
! continues until the number of "water" points no longer grows, which
! is used as termination signal.  At the end the new landmask is
! defined as the points which are not "water" (i.e., either initially
! masked points, or unmasked, but not labelled as "water").

#define TIMING
#define CLOSE_NARROW_BAYS

      implicit none
      real(kind=8), dimension(Lm+2,Mm+2) :: rmask,outmask
      integer(kind=2), dimension(:,:), allocatable :: imask, mss,mss2
      integer nargs, iargc, ncgrd, i0,j0, Lm,Mm, i,j, ierr, lstr, lgrd
      logical show_changes

      include "netcdf.inc"

cf2py intent(in) Lm,Mm,i0,j0,rmask
cf2py intent(out) outmask
      show_changes=.false.
        if (rmask(i0,j0) > 0.5) then
          allocate(imask(0:Lm+1, 0:Mm+1))
          allocate(mss (-1:Lm+2,-1:Mm+2))
          allocate(mss2(-1:Lm+2,-1:Mm+2))

C$OMP PARALLEL SHARED(Lm,Mm, i0,j0,show_changes, rmask,imask,mss,mss2)
          call sin_con_thread(Lm,Mm, i0,j0, show_changes, rmask,imask,
     &                                                       mss,mss2)
          outmask=rmask
C$OMP END PARALLEL

        else
          write(*,'(/1x,2A,i4,2x,A,i4,1x,A/)')  '### ERROR: selected ',
     &    'point i =', i0, 'j =', j0, 'is on land. Try another point.'
        endif

      end subroutine single_connect

      subroutine sin_con_thread(Lm,Mm, i0,j0,show_changes, rmask,
     &                                           imask, mss,mss2)
      implicit none
      integer :: trd_count,  wtr_pts,  wtr_pts_bak
      integer Lm,Mm, i0,j0
      logical show_changes
      real(kind=8) rmask(0:Lm+1,0:Mm+1)
      integer(kind=2) imask(0:Lm+1,0:Mm+1), mss(-1:Lm+2,-1:Mm+2),
     &                                     mss2(-1:Lm+2,-1:Mm+2)
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last,
     &        range,  tile, istr,iend,jstr,jend, my_sum, iter, i,j
#ifdef CLOSE_NARROW_BAYS
     &                                                        , m
#endif
C$    integer omp_get_num_threads, omp_get_thread_num
#ifdef TIMING
      integer(kind=4) iclk_start, iclk_end, clk_rate, clk_max
      call system_clock(iclk_start, clk_rate, clk_max)
#endif
      trd_count=0
      wtr_pts=0
      wtr_pts_bak=-1
      numthreads=1                     ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      call set_tiles(Lm,Mm, nsub_x,nsub_y)

c**   nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd==0) write(*,'(/1x,2A,2(I4,1x,A))')   'Enforce that all ',
     & 'water points can be reached from point  i =',i0, 'j =',j0,'...'

      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            if (rmask(i,j) > 0.5) then
              imask(i,j)=1
            else
              imask(i,j)=0
            endif
          enddo
        enddo

        if (istr==0) istr=istr-1 ; if (iend==Lm+1) iend=iend+1
        if (jstr==0) jstr=jstr-1 ; if (jend==Mm+1) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            mss(i,j)=0
            mss2(i,j)=0
          enddo
        enddo
        if (istr<=i0 .and. i0<=iend .and.
     &      jstr<=j0 .and. j0<=jend) then
          mss(i0,j0)=1
        endif
      enddo   !<-- tile
C$OMP BARRIER
                                          ! The following while() loop
      iter=0                              ! body consists of two nearly
      do while (wtr_pts /= wtr_pts_bak)   ! identical segments which
        iter=iter+1                       ! differ only by switching
        my_sum=0                          ! arrays "mss" and "mss2".
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          do j=jstr,jend
            do i=istr,iend
              if ( imask(i,j)>0 .and. ( mss(i,j )>0  .or.
     &              mss(i+1,j)>0  .or.  mss(i,j+1)>0   .or.
     &              mss(i-1,j)>0  .or.  mss(i,j-1)>0 )) then
                my_sum=my_sum+1
                mss2(i,j)=1
              endif
            enddo
          enddo
        enddo  !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,50)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     &                     '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

        iter=iter+1
        my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                  istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          do j=jstr,jend
            do i=istr,iend
              if ( imask(i,j)>0 .and. ( mss2(i,j )>0  .or.
     &               mss2(i+1,j)>0 .or. mss2(i,j+1)>0    .or.
     &               mss2(i-1,j)>0 .or. mss2(i,j-1)>0 )) then
                my_sum=my_sum+1
                mss(i,j)=1
              endif
            enddo
          enddo
        enddo  !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,50)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     &                     '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while()


      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        if (show_changes) then
          do j=jstr,jend
            do i=istr,iend
              if (mss(i,j)==0 .and. imask(i,j)>0) imask(i,j)=-1
            enddo
          enddo
        else
          do j=jstr,jend
            do i=istr,iend
              if (mss(i,j)==0) imask(i,j)=0
            enddo
          enddo
        endif
      enddo
C$OMP BARRIER


! The following part is to identify and close unresolved features
! of land mask, such as single-point bays.  For each water point the
! immediately adjacent surroundings is inspected and a weighting value
! is calculated by adding 8 for each water point adjacent on the east,
! west, north, and south sides, while adding 1 for each water point
! adjacent in diagonal direction.  If the resultant value falls below
! a pre-set threshold, the point is set to land.   The procedure is
! repeated iteratively until it is detected that further iterations
! do not result in any progress. In the three examples below "."
! means water, "x" land.
!
!       . . x x x x       . . x x x x x       . . x x x x x x x .
!       . . . . x x       . . . x x x x       . . . x x x x . . .
!       . . . . x x       . . . . o o x       . . . . . . . . . .
!       . . . x x x       . . x x x x x       . . x x x x . . . .
!       . . x x x x       . . x x x x x       . . x x x x x x . .
!
!   all water and         initially water       narrow passage: all
!   land points will      points "o" will       points will be kept
!   be kept as is         be turned to land     as they are


#ifdef CLOSE_NARROW_BAYS
C$OMP CRITICAL(cr_region)
      if (trd_count==0) then
        wtr_pts_bak=wtr_pts+1
      endif

      trd_count=trd_count+1

      if (trd_count==numthreads) then
        trd_count=0
        write(*,'(/1x,A)') 'Closing unresolved narrow bays ...'
      endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

      iter=0
      do while (wtr_pts /= wtr_pts_bak)
        iter=iter+1

        do tile=my_first,my_last
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)
          do j=jstr,jend
            do i=istr,iend
              m=0
              if (imask(i+1,j  )>0) m=m+8
              if (imask(i+1,j+1)>0) m=m+1
              if (imask(i  ,j+1)>0) m=m+8
              if (imask(i-1,j+1)>0) m=m+1
              if (imask(i-1,j  )>0) m=m+8
              if (imask(i-1,j-1)>0) m=m+1
              if (imask(i  ,j-1)>0) m=m+8
              if (imask(i+1,j-1)>0) m=m+1
              mss(i,j)=m
            enddo
          enddo

          if (istr==1) then
            do j=jstr,jend
              m=0
              if (imask(istr-1,j-1)>0) m=m+8
              if (imask(istr  ,j-1)>0) m=m+1
              if (imask(istr  ,j  )>0) m=m+8
              if (imask(istr  ,j+1)>0) m=m+1
              if (imask(istr-1,j+1)>0) m=m+8
              mss(istr-1,j)=m
            enddo
          endif

          if (iend==Lm) then
            do j=jstr,jend
              m=0
              if (imask(iend+1,j+1)>0) m=m+8
              if (imask(iend  ,j+1)>0) m=m+1
              if (imask(iend  ,j  )>0) m=m+8
              if (imask(iend  ,j-1)>0) m=m+1
              if (imask(iend+1,j-1)>0) m=m+8
              mss(iend+1,j)=m
            enddo
          endif

          if (jstr==1) then
            do i=istr,iend
              m=0
              if (imask(i+1,jstr-1)>0) m=m+8
              if (imask(i+1,jstr  )>0) m=m+1
              if (imask(i  ,jstr  )>0) m=m+8
              if (imask(i-1,jstr  )>0) m=m+1
              if (imask(i-1,jstr-1)>0) m=m+8
              mss(i,jstr-1)=m
            enddo
          endif

          if (jend==Mm) then
            do i=istr,iend
              m=0
              if (imask(i-1,jend+1)>0) m=m+8
              if (imask(i-1,jend  )>0) m=m+1
              if (imask(i  ,jend  )>0) m=m+8
              if (imask(i+1,jend  )>0) m=m+1
              if (imask(i+1,jend+1)>0) m=m+8
              mss(i,jend+1)=m
            enddo
          endif

          if (istr==1 .and. jstr==1) then
            m=0
            if (imask(istr  ,jstr-1)>0) m=m+8
            if (imask(istr  ,jstr  )>0) m=m+1
            if (imask(istr-1,jstr  )>0) m=m+8
            mss(istr-1,jstr-1)=m
          endif

          if (istr==1 .and. jend==Mm) then
            m=0
            if (imask(istr-1,jend  )>0) m=m+8
            if (imask(istr  ,jend  )>0) m=m+1
            if (imask(istr  ,jend+1)>0) m=m+8
            mss(istr-1,jend+1)=m
          endif

          if (iend==Lm .and. jstr==1) then
            m=0
            if (imask(iend+1,jstr  )>0) m=m+8
            if (imask(iend  ,jstr  )>0) m=m+1
            if (imask(iend  ,jstr-1)>0) m=m+8
            mss(iend+1,jstr-1)=m
          endif

          if (iend==Lm .and. jend==Mm) then
            m=0
            if (imask(iend  ,jend+1)>0) m=m+8
            if (imask(iend  ,jend  )>0) m=m+1
            if (imask(iend+1,jend  )>0) m=m+8
            mss(iend+1,jend+1)=m
          endif
        enddo  !<-- tile
C$OMP BARRIER

        my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                                 istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          if (show_changes) then
            do j=jstr,jend
              do i=istr,iend
                if (mss(i,j)<12 .and. imask(i,j)>0) imask(i,j)=-2
                if (imask(i,j)>0) my_sum=my_sum+1
              enddo
            enddo
          else
            do j=jstr,jend
              do i=istr,iend
                if (mss(i,j)<12) imask(i,j)=0
                if (imask(i,j)>0) my_sum=my_sum+1
              enddo
            enddo
          endif
        enddo

C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif

        trd_count=trd_count+1

        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,1)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     &                     '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while()
#endif  /* CLOSE_NARROW_BAYS */


      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     &                               istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        if (show_changes) then
          do j=jstr,jend
            do i=istr,iend
              if (imask(i,j) > 0) then
                rmask(i,j)=1.D0            !<-- water stays water
              elseif (imask(i,j)==-1) then
                rmask(i,j)=-1.D0           !<-- water turned into land
              elseif (imask(i,j)==-2) then
                rmask(i,j)=-2.D0           !<-- narrow bays closed
              else
                rmask(i,j)=0.D0
              endif
            enddo
          enddo
        else
          do j=jstr,jend
            do i=istr,iend
              if (imask(i,j) > 0) then
                rmask(i,j)=1.D0
              else
                rmask(i,j)=0.D0
              endif
            enddo
          enddo
        endif
      enddo
C$OMP BARRIER

#ifdef TIMING
      call system_clock(iclk_end, clk_rate, clk_max)
      if (clk_rate > 0) then
C$OMP CRITICAL(cr_region)
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          write(*,'(/ /1x,A,F8.2,1x,A,I4,1x,A/10x,A/)')
     &             'Wall clock time spent in computational part',
     &                      (iclk_end-iclk_start)/dble(clk_rate),
     &                    'sec  running', numthreads, 'threads.',
     &                     '[reading and writing files excluded]'
        endif
C$OMP END CRITICAL(cr_region)
      endif
#endif
      end subroutine sin_con_thread
