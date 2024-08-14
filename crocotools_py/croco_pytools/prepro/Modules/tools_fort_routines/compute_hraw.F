      subroutine compute_hraw (lon,lat, xlon,ylat, htopo, Lm,Mm,
     &                 lon_r,lat_r, pm,pn, h, radius)

! A universal interpolation/coarsening procedure to compute topography
! at model grid.   It takes topography data "htopo" defined at source
! grid (xlon,ylat) and transfers it onto the target grid (lon_r,lat_r)
! using weighted averaging,
!
!          weight(r) = [1-(r/width)^2]^2,    r < width
!
! or, ideally, but computationally costly, using HAN WINDOW,
!
!          weight(r) = 1/2 + 1/2*cos(pi*r/width),    r < width
!
! or its polynomial fit
!
!          weight(r) = 1/2 + 1/2*cos(pi*r/width),    r < width
!
! or its polynomial fit
!
!     weight(r)=[ 1- 1.228..*(r/width)^2 +0.228..*(r/width)^4 ]^2
!
! with coefficient chosen to make weight(width/2)= 1/2 (see below).
!
! where "r" is the distance between target point and data point, and
! "width" is specified averaging width.   The goal here is to achieve
! isotropic response of the algorithm to unresolved spike in the data
! (delta function) in sense that "width", expressed as real distance
! (in meters) should be the same in all directions.  Furthermore,
! "width" controlled by non-dimensional input parameter "radius" as
! multiple of local grid spacing
!
!         width[meters] = radius/sqrt(0.5*(pm^2+pn^2))
!
! Translated into geographical lon-lat coordinates, this means that
! "wdthX" (width along lon, in degrees) and "wdthY" (along lat,
! degrees) are related as wdthY = wdthX*cos(lat), and,
!
!         wdthY[degrees] = width[meters]*180/(pi*Eradius)
!
! where Eradius is radius of the earth. The weighted coarsening is
! activated by CPP-switch WEIGHTED_AVERAGING; otherwise a simple
! bi-linear interpolation is performed.
!
! PARALLELIZATION METHOD: This routine is expected to be called from
! inside an Open MP parallel region and it employs dynamic scheduling
! by processing one j-row of array at a time,
!
!         do while(.true.)
!     C$OMP CRITICAL
!           j=comm_j_index       <-- check out one j-row to process
!           comm_j_index=j+1     <-- increment shared counter to signal
!     C$OMP END CRITICAL             to the other threads that this row
!           if (j > Mm+1) --> DONE   is already taken
!           do i=....
!            .... do useful work
!           enddo
!         enddo
!
! where "comm_j_index" is a shared counter initialized to zero used as
! a signalling variable. This method is preferred to more conventional
! tiling here because the input topographic data array "htopo" may be
! very large, so the algorithm here is designed to keep the threads
! working close to each other in memory space (perhaps synergetically
! benefit from shared shared cache of a multicore CPU by mutually
! reusing data brought into cache by different threads).

! Created and maintained by Alexander Shchepetkin old_galaxy@yahoo.com

#define TIMING
#define WEIGHTED_AVERAGING
c-#define VERBOSE

      implicit none
      integer::  lon,lat,Lm,Mm!, ierr
      real(kind=8), dimension(lon,lat):: htopo
      real(kind=8) ::  radius, xlon(lon), ylat(lat)
      real(kind=8), dimension(0:Lm+1,0:Mm+1):: lon_r,lat_r,pm,pn
      real(kind=8), dimension(0:Lm+1,0:Mm+1):: h
      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

!--> internal variables

      integer, dimension(0:Lm+1) :: ic,jc
      real(kind=8), dimension(0:Lm+1) :: xi,eta
      real(kind=8):: xwest,xeast,ysouth,ynorth, dx_inv, dy_inv, xr,yr
      integer:: numthreads, trd, iwest,ieast,jsouth,jnorth,
     &                            i,j, ir,jr, iter
C$    integer omp_get_num_threads, omp_get_thread_num

#ifdef WEIGHTED_AVERAGING
      logical :: bounded
      integer :: istr,iend,jstr,jend, imin,imax,jmin,jmax, mwx,mwy,
     &         mwx_alloc,mwy_alloc, irc,jrc, itm,jtm, iups,jups
      real(kind=8), dimension(:,:), allocatable :: wgt,wgt_new
      real(kind=8)::  SmScle,wdthX,wdthY, sum,biasX,biasY, xc,yc,
     &                       cff,ctr, cx,cy, dcx,dcy, acx,acy
#endif
#ifdef TIMING
      integer(kind=4):: iclk_start, iclk_end, clk_rate, clk_max
#endif

cf2py intent(in) Lm,Mm,lon,lat,xlon,ylat,htopo,lon_r,lat_r,pm,pn,radius
cf2py intent(out) h 

      include "phys_const.h"
#ifdef WEIGHTED_AVERAGING
      mwx_alloc=0  ; mwy_alloc=0
      dx_inv=dble(lon-1)/(xlon(lon)-xlon(1))
      dy_inv=dble(lat-1)/(ylat(lat)-ylat(1))
#endif
      comm_j1_indx=0
      comm_j2_indx=0
      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()

#ifdef TIMING
      if(trd == 0) call system_clock(iclk_start, clk_rate, clk_max)
#endif
                     ! Check whether the data coverage is sufficient
! Preliminary step:  ! for the entire ROMS  grid. Complain about error
!------------------  ! and signal to quit, if anything goes wrong.

! Find extremes of lon-lat oordinates for topography data set

      xwest=xlon(1)  ;   xeast=xlon(lon)
      ysouth=ylat(1) ;   ynorth=ylat(lat)

! Find extremes of lon-lat coordinates for croco data set
      call roms_grid_geo_bounds(lon_r,lat_r,Lm,Mm,radius,
     &                     lon_min,lon_max,lat_min,lat_max)

C$OMP MASTER
      write(*,'(/1x,2A,2(/4x,A,2F12.6))')  'Geographical extents of ',
     &        'the dataset:',  'Longitude min/max =', xwest,  xeast,
     &                         'Latitude  min/max =', ysouth, ynorth
      write(*,'(/1x,2A,2(/4x,A,2F12.6))')  'Geographical extents of ',
     & 'requested ROMS grid:', 'Longitude min/max =', lon_min,lon_max,
     &                         'Latitude  min/max =', lat_min,lat_max
C$OMP END MASTER

      iter=0
      if (lon_min < xwest) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)')  '### WARNING: ',
     &    'Western extent of model grid is beyond western bound of ',
     &    'dataset.', 'lon_min_data =', xwest, 'lon_min =', lon_min
C$OMP END MASTER
      endif
      if (lon_max > xeast) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)')  '### WARNING: ',
     &    'Eastern extent of model grid is beyond eastern bound of ',
     &    'dataset.', 'lon_max_data =', xeast, 'lon_max =', lon_max
C$OMP END MASTER
      endif
      if (lat_min < ysouth) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)')   '### WARNING: ',
     &   'Southern extent of model grid is beyond southern bound of ',
     &   'dataset.', 'lat_min_data =', ysouth, 'lat_min =', lat_min
C$OMP END MASTER
      endif
      if (lat_max > ynorth) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)')   '### WARNING: ',
     &   'Northern extent of model grid is beyond northern bound of ',
     &   'dataset.',  'lat_max_data =', ynorth,  'lat_max =', lat_max
C$OMP END MASTER
      endif
      if (iter > 0) then  ! Note: Because the error condition checked
C$OMP MASTER              ! here occurs simultaneously on all threads,
!        ierr=iter         ! all exit simultaneously.  The only reason
C$OMP END MASTER          ! for setting "ierr" is to prevent writing
c        return            ! "hraw" into netCDF file in the case of
      endif               ! error detected here.

      iwest=indx_bound2(xlon,lon, lon_min)
      ieast=indx_bound2(xlon,lon, lon_max)
      jsouth=indx_bound2(ylat,lat, lat_min)
      jnorth=indx_bound2(ylat,lat, lat_max)

C$OMP MASTER
      write(*,'(/1x,2A/2(1x,A,I5,I6,1x,A,I6))')   'Bounding indices ',
     &   'for the portion of dataset covering the entire model grid:',
     &                'iwest,ieast =',  iwest,ieast,   'out of', lon,
     &                'jsouth,jnorth =',jsouth,jnorth, 'out of', lat

      if (iwest == 0) write(*,*)    'WARNING: compute_hraw :: ',
     &                          'restricting iwest from 0 to 1'
      if (ieast == lon) write(*,*)  'WARNING: compute_hraw :: ',
     &             'restricting ieast from', lon, ' to', lon-1
      if (jsouth == 0)  write(*,*)  'WARNING: compute_hraw :: ',
     &                         'restricting jsouth from 0 to 1'
      if (jnorth == lat) write(*,*) 'WARNING: compute_hraw :: ',
     &            'restricting jnorth from', lat, ' to', lat-1



#ifdef WEIGHTED_AVERAGING
      write(*,'(/1x,A,F8.5,1x,A)')  'Averaging window width =',
     &                               radius, 'in grid spaces.'
#else
      write(*,'(/1x,A)') 'Bi-linear interpolation.'
#endif
C$OMP END MASTER
C$OMP BARRIER

#define iwest illegal
#define ieast illegal
#define jsouth illegal
#define jnorth illegal


! Start interpolation:   For each point of model grid, "lon_r,lat_r"
!------ --------------   find indices "ic,jc" of the data grid such
! that
!              xlon(ic) <= lon_r < xlon(ic+1)
!              ylat(jc) <= lat_r < ylat(jc+1)
!
! after which compute fractional distances "xi,eta"
!
!                   0 <= xi,eta < 1
!
! which identify position of point "lon_r,lat_r" relatively to the
! rectangle [xlon(ic):xlon(jc+1)]X[ylat(jc):ylat(jc+1)].
! This part is the same for weighted and bi-linear interpolation.
!
! The rest of this code is organized as as single giant j-loop lasting
! all the way to the end.
      do while(.true.)
C$OMP CRITICAL(progress_sync)
        j=comm_j2_indx ; comm_j2_indx=j+1
# ifdef VERBOSE
        write(*,'(8x,2(A,I5,1x),A,I2)') 'processing row j =', j,
     &                             'out of', Mm+1, 'trd =', trd
# else
        if (mod(j,5) == 0) then
          if (mod(j,375) == 0) write(6,*)
          write(6,'(A)',advance='no') '.'
          flush(unit=6)
        endif
# endif
C$OMP END CRITICAL(progress_sync)
        if (j > Mm+1) exit !goto 98     !--> DONE

! Find bounding indices (ic,jc) for each point (i,j) of ROMS grid,
! such that ROMS grid points with coordinates lon_r(i,j),lat_r(i,j)
! falls between ic and ic+1, as well as jc and jc+1 on the data grid
! with coordinates xlon,ylat. Mathematically speaking the "while"
! searches in the longitudinal direction in the code below are not
! necessary because the data grid is always uniform in this direction.
! However, the dataset may store coordinates in single precision which
! may cause round-off-level violations of "xi,eta" bounds triggering
! error messages below.  Latitudinal direction may theoretically use
! non-uniform grid spacing, although in most datasets it does not.

! The code below (starting from here and all the way to the end) is
! adapted to handle the situation that the dataset may not cover the
! entire ROMS grid, and even not to cover the entire rectangular
! portion [xwest,xeast] x [ysouth,ynorth] due to tiled approach for
! storing data in multiple netCDF files in such a way that individual
! tiles use lon,lat coordinates within, but mutial arrangement of the
! tiles is rather arbitrary. For this reason (ic,jc) are set first to
! unrealistic special values, which are then owerwritten with
! meaningful indices only if interpolation is possible.

#ifdef WEIGHTED_AVERAGING
# define margine 2
#else
# define margine 1
#endif
        do i=0,Lm+1
          ic(i)=-1 ; jc(i)=-1 !<-- special values
          if (xwest < lon_r(i,j) .and. lon_r(i,j) < xeast) then
            ir=1+int(dx_inv*(lon_r(i,j)-xwest))
            ir=max(min(ir, lon-margine), margine)
            do while(xlon(ir) > lon_r(i,j) .and. ir > margine)
              ir=ir-1
            enddo
            do while(xlon(ir+1) <= lon_r(i,j) .and. ir < lon-margine)
              ir=ir+1
            enddo
            if ( xlon(ir) <= lon_r(i,j) .and. lon_r(i,j) < xlon(ir+1)
     &                     ) ic(i)=ir
          endif

          if (ysouth < lat_r(i,j) .and. lat_r(i,j) < ynorth) then
            jr=1+int(dy_inv*(lat_r(i,j)-ysouth))
            jr=max(min(jr, lat-margine), margine)
            do while(ylat(jr)> lat_r(i,j) .and. jr > margine)
              jr=jr-1
            enddo
            do while(ylat(jr+1) <= lat_r(i,j) .and. jr < lat-margine)
              jr=jr+1
            enddo
            if ( ylat(jr) <= lat_r(i,j) .and. lat_r(i,j) < ylat(jr+1)
     &                   ) jc(i)=jr
          endif
#undef margine

! Compute compute fractional distances "xi,eta" which locate the
! position of "lon_r,lar_r" within the 2D-cell [ir:ir+1]X[jr:jr+1] on
! the data grid. If everything goes correctly both "xi" and "eta" must
! be bounded as 0 <= xi,eta < 1.

          if (ic(i) > 0 .and. jc(i) > 0) then
            ir=ic(i) ; jr=jc(i)
            xi(i)= (lon_r(i,j)-xlon(ir))/(xlon(ir+1)-xlon(ir))
            eta(i)=(lat_r(i,j)-ylat(jr))/(ylat(jr+1)-ylat(jr))
          else
            xi(i)=0.D0 ; eta(i)=0.D0
          endif
        enddo


#ifdef VERBOSE
        iter=0                     !--> This is just for checking.
        do i=0,Lm+1
          if (xi(i) < 0.D0 .or. xi(i) > 1.D0) then
            write(*,'(1x,2A,3I8,A,F9.5)') '### ERROR: compute_hraw ::',
     &      ' Algorithm failure X: i,j,ic=',i,j,ic(i), ' xi =', xi(i)
            iter=iter+1
          endif
          if (eta(i) < 0.D0 .or. eta(i) > 1.D0) then
            write(*,'(1x,2A,3I8,A,F9.5)') '### ERROR: compute_hraw ::',
     &      ' Algorithm failure Y: i,j,jc=',i,j,jc(i),' eta =',eta(i)
            iter=iter+1
          endif
        enddo
        if (iter > 0) then
C$OMP CRITICAL
!          ierr=ierr+iter
C$OMP END CRITICAL
        endif
#endif



#ifdef WEIGHTED_AVERAGING

! Find bounding indices "istr:iend" and "jstr:jend" enclosing all data
! grid points which lie within the distances "wdthX" (in longitudinal)
! and "wdthY" (in meridional direction) from the point "lon_r,lat_r"
! of the model grid.

        SmScle=radius *180.D0/(pi*Eradius*sqrt(0.5D0))

        do i=0,Lm+1
          bounded=.false.
          if (ic(i) > 0 .and. jc(i) > 0) then
            irc=ic(i) ; xc=lon_r(i,j) ; istr=irc ; iend=irc+1
            jrc=jc(i) ; yc=lat_r(i,j) ; jstr=jrc ; jend=jrc+1

            wdthY=SmScle/sqrt(pm(i,j)**2+pn(i,j)**2)
            wdthX=wdthY/cos(deg2rad*yc)

c*          wdthX=0.1*wdthX                !<-- checking
c*          wdthY=0.1*wdthY                !<-- isotropy


! Expand eastern and western bounding indices until the interval
! xlon(istr):xlon(iend) is wide enough to cover "xc-wdthX:xc+wdthX".
! Similarly for southern and northern bounding indices to make
! ylat(jstr):ylat(jend) cover "yc-wdthY:yc+wdthY".


            do while(xc-wdthX < xlon(istr) .and. istr > 3)
              istr=istr-1
            enddo
            do while(xlon(iend) < xc+wdthX .and. iend < lon-2)
              iend=iend+1
            enddo

            do while(yc-wdthY < ylat(jstr) .and. jstr > 3)
              jstr=jstr-1
            enddo
            do while(ylat(jend) < yc+wdthY .and. jend < lat-2)
              jend=jend+1
            enddo

            if (xlon(istr) < xc-wdthX .and. xc+wdthX < xlon(iend) .and.
     &          ylat(jstr) < yc-wdthY .and. yc+wdthY < ylat(jend)) then
              bounded=.true.
            endif
          endif  !<--  ic(i) > 0 .and. jc(i) > 0



          if (bounded) then
            mwx=max(iend-irc+1,irc-istr+1)
            mwy=max(jend-jrc+1,jrc-jstr+1)
            if (mwx > mwx_alloc .or. mwy > mwy_alloc) then
              mwx_alloc=mwx+1            ! add 1 merely to reduce
              mwy_alloc=mwy+1            ! probability of the need
              if (allocated(wgt)) then   ! to allocate arrays of a
                deallocate(wgt)          ! larger size later
                deallocate(wgt_new)
              endif
              allocate(wgt(-mwx_alloc:mwx_alloc,-mwy_alloc:mwy_alloc),
     &             wgt_new(-mwx_alloc:mwx_alloc,-mwy_alloc:mwy_alloc))
# ifdef VERBOSE
C$OMP CRITICAL (wgt_alloc_sync)
              write(*,'(32x,A,2(1x,A,I4),2x,A,I2)') 'allocated window',
     &       'array:', 2*mwx_alloc+1, 'x', 2*mwy_alloc+1, 'trd =',trd
C$OMP END CRITICAL (wgt_alloc_sync)
# endif
            endif

! Construct interpolation/weighting coefficients "wgt":
!---------- ------------- --------- ------------ -------
! The intended shape of weight is a cosine function,
!
!            f(x) = 1/2 + (1/2)*cos(pi*x) ,  |x|<1
!
! Below it is approximated by f(x)=(1- 1.228..*x^2 +0.228..*x^4)^2
! where 1.22876383367175=5-8*sqrt(2)/3 is chosen to make f(1/2)=1/2
! like it is for cosine function above.  The difference between the
! two functions within the entire interval |x|<1 is between -1.1e-3
! to +5e-4.

! Note that in the case when width of averaging window is narrower
! than the distance between adjacent data points, there is a chance
! that all weights computed by the algorithm below are all zeros. In
! anticipation of this possibility weights are set initially to delta
! function at the closest discrete point, and are then only if the
! window is wide enough for averaging to occur (hence the assignment
! of "weight" occurs conditionally in the 2D ir-jr loop below).
! Subsequently, in the case when "weight" remain delta-function, the
! exact coefficients for bi-linear interpolation are recovered by the
! iterative iterative centering procedure.

            do jr=jstr-1,jend+1
              do ir=istr-1,iend+1
                wgt(ir-irc,jr-jrc)=0.D0 ; wgt_new(ir-irc,jr-jrc)=0.D0
              enddo
            enddo
            if (xi(i) < 0.5D0) then
              ir=0
            else
              ir=1
            endif
            if (eta(i) < 0.5D0) then
              jr=0
            else
              jr=1
            endif
            wgt(ir,jr)=1.D0    !<-- delta function

            dcx=1.D0/wdthX  ;  dcy=1.D0/wdthY
            sum=0.D0 ; biasX=0.D0 ;  biasY=0.D0


            do jr=jstr,jend
              do ir=istr,iend
                xr=xlon(ir)  ;  yr=ylat(jr)
                cff=(dcx*(xr-xc))**2 + (dcy*(yr-yc))**2
                if (cff < 1.D0) then
                  wgt(ir-irc,jr-jrc)=( 1.D0-cff*( 1.22876383367175D0
     &                                  -0.22876383367175D0*cff ))**2
                endif
                sum  =sum   + wgt(ir-irc,jr-jrc)
                biasX=biasX + wgt(ir-irc,jr-jrc)*xr
                biasY=biasY + wgt(ir-irc,jr-jrc)*yr
              enddo
            enddo

            if (sum > 0.D0) then
              cx=dx_inv*(xc-biasX/sum)      ! initial estimate of
              cy=dy_inv*(yc-biasY/sum)      ! off-centering errors

              if (abs(cx) > 1.D0 .or. abs(cy) > 1.D0) then
                 write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     &              'failure, cx =', cx, ' cy =',cy, ' istr,iend =',
     &                        istr,iend, ' jstr,jend =', jstr,jend
                 stop
              endif
            else
              write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     &                                      'failure, line 483.'
              stop
            endif

! Weight-centering iterations:  Although the initial setting of "wgt"
!----------------- -----------  is centered already close to (xi,eta)
! defined as the fractional distance between points (ir:ir+1,jr:jr+1),
! this centering may not be accurate enough if too few points are
! participating in averaging. The iterative procedure below compares
! the actual location of the center of gravity with the desired one,
! (xi,eta), and modifies the weights using a first-order upstream
! advection scheme. In essence, it computes the fractional shift
! (cx,cy) which moves the center of gravity  exactly to (xi,eta).

            do iter=1,6
              if (cx > 0.D0) then
                iups=+1 ; acx=cx  ; imin=istr ; imax=iend+1
              else
                iups=-1 ; acx=-cx ; imin=istr-1 ; imax=iend
              endif
              if (cy > 0.D0) then
                jups=+1 ; acy=cy  ; jmin=jstr ; jmax=jend+1
              else
                jups=-1 ; acy=-cy ; jmin=jstr-1 ; jmax=jend
              endif

              ctr=(1.D0-acx)*(1.D0-acy) ; dcx=acx*(1.D0-acy)
              dcy=(1.D0-acx)*      acy  ; cff=acx  *    acy

              sum=0.D0 ; biasX=0.D0 ; biasY=0.D0

              do jr=jmin,jmax
                do ir=imin,imax
                  itm=ir-irc ; jtm=jr-jrc

                  wgt_new(itm,jtm)=ctr*wgt(itm     ,jtm     )
     &                            +dcx*wgt(itm-iups,jtm     )
     &                            +dcy*wgt(itm     ,jtm-jups)
     &                            +cff*wgt(itm-iups,jtm-jups)

                  sum  = sum  + wgt_new(itm,jtm)
                  biasX=biasX + wgt_new(itm,jtm)*xlon(ir)
                  biasY=biasY + wgt_new(itm,jtm)*ylat(jr)
                enddo
              enddo

              if (sum > 0.D0) then
                dcx=dx_inv*(xc-biasX/sum)      ! normalized residual
                dcy=dy_inv*(yc-biasY/sum)      ! off-centering errors

c**             write(*,'(I2,4F19.16)') iter,cx,dcx, cy,dcy

                cx=cx+dcx  ;  cy=cy+dcy       !<-- adjust cx,cy

                if (abs(cx) > 1.D0 .or. abs(cy) > 1.D0) then
                   write(*,*) '### ERROR: compute_hraw :: cx =',cx,
     &                                                  ' cy =',cy
                   stop
                endif
              else
                write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     &                                      'failure at line 542.'
                stop
              endif                            ! After this adjustment
            enddo !<-- iter                    ! "wgt_new" is centered
                                               ! exactly at (xi,eta);
            sum=0.D0 ; biasX=0.D0 ; biasY=1.D0 ! use it to interpolate
            do jr=jmin,jmax                    ! htopo--> hraw.
              do ir=imin,imax
                if (htopo(ir,jr) < 1.E+32) then
                  sum  = sum +wgt_new(ir-irc,jr-jrc)
                  biasX=biasX+wgt_new(ir-irc,jr-jrc)*dble(htopo(ir,jr))
                else
                  biasY=-1.D0
                endif
              enddo
            enddo
            if (biasY < 0.D0) then
              h(i,j)=1.D+33
            elseif (sum > 0.D0) then
              h(i,j)=biasX/sum
            else
              write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     &                                    'failure at line 562.'
              stop
            endif

          else   ! <-- .not.bounded
            h(i,j)=0.D0
          endif
        enddo   ! <-- i

# else

! Alternative to WEIGHTED_AVERAGING: Bi-linear interpolation using
! four nearest data points.   Note that in the case when the target
! grid is coarser than the source grid this approach results in
! sub-sampling and ultimately aliasing error.  This part of code is
! kept for algorithm verification purposes only and not to be used
! in practice.

        do i=0,Lm+1
          if (ic(i) > 0 .and. jc(i) > 0) then
            ir=ic(i) ; jr=jc(i)

            h(i,j)=(1.D0-xi(i))*(1.D0-eta(i))*dble(htopo(ir  ,jr  ))
     &                  +xi(i) *(1.D0-eta(i))*dble(htopo(ir+1,jr  ))
     &            +(1.D0-xi(i))*      eta(i) *dble(htopo(ir  ,jr+1))
     &                  +xi(i) *      eta(i) *dble(htopo(ir+1,jr+1))
          else
            h(i,j)=0.D0
          endif
        enddo
#endif

        do i=0,Lm+1                     
            h(i,j)=-h(i,j)  ! Convert into CROCO standart

                            ! but do not set h<0 to 0
                            ! as it can serve when doing wet_and_dry
                                         ! Convert topographic data
!          if (h(i,j) < 0.D0) then       ! into ROMS standard: ETOPO
!            h(i,j)=-h(i,j)              ! and SRTM conventions imply
!          else                          ! that "htopo" is elevation,
!            h(i,j)=0.D0                 ! hence it is positive above
!          endif                         ! and negative below sea
        enddo                           ! level.  ROMS convention
      enddo   !<-- j, while(.true.)     ! anticipates positive depth.
!  98  continue
C$OMP BARRIER
#ifdef TIMING
      if (trd == 0) then
        call system_clock(iclk_end, clk_rate, clk_max)
        if (clk_rate>0) then
          write(*,'(/ /1x,2A,F8.2,1x,A,I4,1x,A)')   'Wall Clock time ',
     &   'spent to compute hraw', (iclk_end-iclk_start)/dble(clk_rate),
     &                           'sec running', numthreads, 'threads.'
        endif
      endif
#endif
      
      contains
        integer function indx_bound2(x,n, x0)
        integer n, i      ! bounded by
        real(kind=8) x(n), x0         !         x(i) <= x0 < x(i+1)

        if (x0 < x(1)) then
          i=0                         ! if x0 is outside the full range
        elseif (x0 > x(n)) then       ! of x(1) ... x(n), then return
          i=n                         ! i=0 or i=n.
        else
          i=int( ( x(n)-x0 +n*(x0-x(1)) )/(x(n)-x(1)) )
          if (x(i+1)<x0) then
            do while (x(i+1) < x0)    ! This algorithm computes "i" as
              i=i+1                   ! linear interpolation between x(1)
            enddo                     ! and x(n) which should yield the
          elseif (x(i) > x0) then     ! correct value for "i" right a way
            do while (x(i) > x0)      ! because array elements x(i) are
              i=i-1                   ! equidistantly spaced.  The while
            enddo                     ! loops are here merely to address
          endif                       ! possible roundoff errors.
 
          if (x(i+1)-x0 < 0 .or. x0-x(i) < 0) then
            write(*,'(1x,A,5F12.6)') '### ERROR: indx_bound :: ',
     &                   x(i), x0, x(i+1), x0-x(i), x(i+1)-x0
            stop
          endif
        endif
        indx_bound2=i
        end function indx_bound2
     
      end subroutine compute_hraw


      subroutine roms_grid_geo_bounds(lon_r,lat_r, Lm,Mm,radius,
     &           lon_min,lon_max,lat_min,lat_max)

! Determine geographical limits of data needed to generate ROMS grid
! topography. Algorithmically this is done by going along the perimeter
! of the grid and slightly extrapolating nx-ny coordinates outside to
! accommodate for the data points needed for averaging to compute the
! topography for the points on the perimeter of ROMS grid. The distance
! for extrapolation depends on ROMS grid spacing and the smoothing
! "radius" (with an extra 50% margin for safety of the logic).
!
!                   input:  lon_r,lat_r, Lm,Mm, radius,
!                  output:  lon_min,lon_max,lat_min,lat_max

!      use comm_vars_hraw
      implicit none
      integer Lm,Mm, i,j
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: lon_r,lat_r
      real(kind=8) radius,  cff, lon_ext, lat_ext

      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

cf2py intent(in) lon_r,lat_r, Lm,Mm,radius
cf2py intent(out) lon_min,lon_max,lat_min,lat_max
      comm_j1_indx=0
      comm_j2_indx=0
      cff=1.5D0*radius  !<-- safety margin

      lon_min=lon_r(0,0)+cff*(lon_r(0,0)-lon_r(1,1))
      lon_max=lon_min
      lat_min=lat_r(0,0)+cff*(lat_r(0,0)-lat_r(1,1))
      lat_max=lat_min

      do i=1,Lm
        lon_ext=lon_r(i,0)+cff*(lon_r(i,0)-lon_r(i,1))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(i,0)+cff*(lat_r(i,0)-lat_r(i,1))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(Lm+1,0)+cff*(lon_r(Lm+1,0)-lon_r(Lm,1))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(Lm+1,0)+cff*(lat_r(Lm+1,0)-lat_r(Lm,1))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do j=1,Mm
        lon_ext=lon_r(Lm+1,j)+cff*(lon_r(Lm+1,j)-lon_r(Lm,j))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(Lm+1,j)+cff*(lat_r(Lm+1,j)-lat_r(Lm,j))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(Lm+1,Mm+1)+cff*(lon_r(Lm+1,Mm+1)-lon_r(Lm,Mm))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(Lm+1,Mm+1)+cff*(lat_r(Lm+1,Mm+1)-lat_r(Lm,Mm))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do i=Lm,1,-1
        lon_ext=lon_r(i,Mm+1)+cff*(lon_r(i,Mm+1)-lon_r(i,Mm))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(i,Mm+1)+cff*(lat_r(i,Mm+1)-lat_r(i,Mm))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(0,Mm+1)+cff*(lon_r(0,Mm+1)-lon_r(1,Mm))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(0,Mm+1)+cff*(lat_r(0,Mm+1)-lat_r(1,Mm))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do j=Mm,1,-1
        lon_ext=lon_r(0,j)+cff*(lon_r(0,j)-lon_r(1,j))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(0,j)+cff*(lat_r(0,j)-lat_r(1,j))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      end subroutine roms_grid_geo_bounds
