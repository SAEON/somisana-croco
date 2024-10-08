      subroutine r2r_match_topo (OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH, 
     &           ic,jc,iip,jjp,wdth,xc,yc,hchd,xi,xp,yp,hp,hwg)


!xp lon_rho parent
!yp lat_rho parent
!hp topo parent
!
!xc lon_rho child
!yc lat_rho child
!hwg child topo              
!xi mask_rho child
!wght weight child
!hc new htopo


! Modifies child-grid bottom topography in the vicinity of forced
! open boundaries to exactly match topography interpolated from the
! parent grid right at the boundary, while at the same time, time
! making smooth transition into the interior of child-grid domain
! where the topography is expected to be less smooth because of
! finer resolution.

! The method is essentially generating a 2D weighting function,
! wgt=wgt(i,j), 0 <= wgt <= 1, hence
!
!             h_matched = (1-wgt) * h + wgt * h_parent
!
! where wgt=1 at open boundaries, and wgt=0 in the interior beyond
! merging zone, with smooth transition within - sounds quite simple,
! however is should be noted that merging should be avoided in the
! areas closed to the edge of the child domain, but blocked by the
! land mask. This in its turn leads to a rather complicated logical
! "etching" to determine the exact coverage of merging zone.
! Interpolation of parent-grid topography is by bi-cubic spline
! using Hermite polynomials.

! Created/maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

c--#define LEGACY_METHOD
#define SMOOTH_JUNCTIONS
      implicit none
      character(len=4) merging_flags
      character(len=64) prntgrd, chldgrd, str
      integer :: ic,jc,iip,jjp
      real(kind=8), dimension(ic+2,jc+2) :: xc,yc,xi,eta
     &                                     ,wgt,hwg,hchd
      real(kind=8), dimension(iip+2,jjp+2) :: hp,xp,yp
      
      real(kind=8), allocatable, dimension(:,:) :: hc,srX,srY,sXY,sYX
      real(kind=4), allocatable, dimension(:,:) :: mwgt
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp
      integer(kind=2), allocatable, dimension(:,:) ::mask,mgz,ms1,ms2
      integer nargs, wdth, net_alloc_mem, nx,ny, ncx,ncy, ncsrc,nctarg,
     &        ncid, ierr, ipvar,jpvar, xivar,etavar, hpvar,hvar, dhvar,
     &        wgtvar,mwgvar, r2dgrd(2), i,j,k, iter, lprnt,lchld
      real(kind=8) xmin,xmax, xcmin,xcmax, cff
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH, rename_hvar

       integer rad,rad2, ict,jct

      include "netcdf.inc"

Cf2py intent(in) OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH,ic,jc,iip,jjp,wdth,xc,yc,hchd,xi,xp,yp,hp
Cf2py intent(out) hwg

      net_alloc_mem=0 !<-- allocated memory expressed in 4-byte numbers
      merging_flags='    '; rename_hvar=.false.

! Check whether all the arguments are specified correctly and both
! netCDF files can be opened (at least in read-only mode); diagnose
! errors and quite if an errors occurs; write help page and quit if
! the program is called with insufficient number of arguments.

        k=0
        if (OBC_WEST) then                    ! prepare signature
          k=k+1 ; merging_flags(k:k)='W'      ! string to be saved
        endif                                 ! as an attribute in
        if (OBC_EAST) then                    ! netCDF file
         k=k+1 ; merging_flags(k:k)='E'
        endif
        if (OBC_SOUTH) then
          k=k+1 ; merging_flags(k:k)='S'
        endif
        if (OBC_NORTH) then
          k=k+1 ; merging_flags(k:k)='N'
        endif

! Allocate arrays and read coordinates first for the target grid,
! then for the source.  This leads to a more optimal memory use as
! xp,yp arrays for the source grid can be deallocated after computing
! ip,jp-indices and fractional offsets xi,eta, the xp,yp are allocated
! at the end to be deallocated to free memory for arrays associated
! with spline interpolation.

      ncx=ic+2 ;  ncy=jc+2

      if (ierr == 0) allocate( mwgt(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( ip(ncx,ncy),jp(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( mask(ncx,ncy), mgz(ncx,ncy),
     &                       ms1(ncx,ncy),ms2(ncx,ncy), stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +ncx*ncy*(6*2 + 1*1 + 2*1)
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif

      do j=1,ncy
        do i=1,ncx
          if (xi(i,j) > 0.5D0) then   ! temporarily use array "xi"
            mask(i,j)=1               ! to read mask from the file;
          else                        ! thereafter this program uses
            mask(i,j)=0               ! only integer(kind=2) version
          endif                       ! of mask
        enddo
      enddo

      xcmin=xc(1,1) ; xcmax=xc(1,1)
      do j=1,ncy
        do i=1,ncx 
          if (xc(i,j) < xcmin) then
            xcmin=xc(i,j)
          elseif (xc(i,j) > xcmax) then
            xcmax=xc(i,j)
          endif
        enddo
      enddo 

! The same for parent grid

      nx=iip+2 ; ny=jjp+2

      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+(3*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
      stop
      endif

! Interpolate parent-grid topography onto child grid...

      xmin=xp(1,1) ; xmax=xp(1,1)
      do j=1,ny
        do i=1,nx
          if (xp(i,j) < xmin) then
            xmin=xp(i,j)
          elseif (xp(i,j) > xmax) then
            xmax=xp(i,j)
          endif
        enddo
      enddo
      write(*,*) 'Parent grid xmin,xmax =', xmin,xmax
      if (xmin > xcmax) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin-360.D0,xmax-360.D0
      elseif (xmax < xcmin) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin+360.D0, xmax+360.D0
      endif

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
        call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                              ip,jp, xi,eta)
C$OMP END PARALLEL

      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)

      net_alloc_mem=net_alloc_mem-4*nx*ny

      allocate(srX(nx,ny),srY(nx,ny), sXY(nx,ny),sYX(nx,ny))
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (4*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)

      deallocate(sYX) ; net_alloc_mem=net_alloc_mem-2*nx*ny
      allocate(hc(ncx,ncy), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (1*2)*ncx*ncy
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif

C$OMP PARALLEL SHARED(nx,ny, hp, srX,srY,sXY, ncx,ncy,
C$OMP&                               ip,jp, xi,eta, hc)
      call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY,
     &                      ncx,ncy, ip,jp, xi,eta, hc)
C$OMP END PARALLEL


! Start forming merging zone and weights: the two methods below differ
!------ --------------- -------- -------- only by handling land mask
! and produce identical results in the case of no land.  Either way,
! the outcome is "mgz" assuming values 0 or 1 to define merging zone
! (unmasked water points only) and "wgt" initialized as constant-slope
! function wgt=wdth at each open boundary water points from which it
! descends linearly into the interior of the domain by decreasing by
! 1 for each row of grid points until it vanishes to 0.  In the case
! of non-LEGACY_METHOD "wgt" also decreases by 1 every step when it
! goes around the corner of the coastline -- the minimal path from
! the costline to the given point is counted as the number of steps
! connecting unmasked water points by C-grid masking rules rather than
! along a straight line to the nearest boundary point, so the resultant
! weights and the extent of merging zone expected to be somewhat
! smaller in areas behind land.

#ifdef LEGACY_METHOD

! Initialize constant-slope function "wgt" and the initial version of
! merging zone "mgz" where wgt > 0 in a straightforward analytical way,
! then apply land mask which potentially may break up "mgz" into
! several unconnected areas some of which cannot be reached by water
! points from the open boundary, then exclude such points by applying
! algorithm similar to enforcing connectivity of unmasked water area
! by single connect.

      do j=1,ncy
        do i=1,ncx
          k=0
          if (OBC_WEST)  k=max(k, wdth-i+1)
          if (OBC_EAST)  k=max(k, wdth-ncx+i)
          if (OBC_SOUTH) k=max(k, wdth-j+1)
          if (OBC_NORTH) k=max(k, wdth-ncy+j)

          if (k > 0 .and. mask(i,j) > 0) then
            mgz(i,j)=1 ; wgt(i,j)=dble(k)
          else
            mgz(i,j)=0 ; wgt(i,j)=0.D0
          endif
        enddo
      enddo

C$OMP PARALLEL SHARED(ncx,ncy, mgz,ms1,ms2)
      call mrg_zone_cont_thread(ncx,ncy, mgz,ms1,ms2)
C$OMP END PARALLEL

#else
      do j=1,ncy                     ! Initialize "mgz" as a single
        do i=1,ncx                   ! row of perimeter points along
          mgz(i,j)=0                 ! the unmasked parts of open
        enddo                        ! boundaries, then "etch" inward,
      enddo                          ! while obeying the connectivity
      if (OBC_WEST) then             ! rules of land-water masking:
        do j=1,ncy                   ! the by the construction "mgz"
          mgz(1,j)=mask(1,j)         ! has the property that each of
        enddo                        ! its points can be reached from
      endif                          ! the open boundary by water.
      if (OBC_SOUTH) then
        do i=1,ncx
          mgz(i,1)=mask(i,1)
        enddo
      endif
      if (OBC_EAST) then
        do j=1,ncy
          mgz(ncx,j)=mask(ncx,j)
        enddo
      endif
      if (OBC_NORTH) then
        do i=1,ncx
          mgz(i,ncy)=mask(i,ncy)
        enddo
      endif

C$OMP PARALLEL SHARED(ncx,ncy, mask,ms1,ms2, wdth)
      call etch_mgz_weights_thread(ncx,ncy, mask,mgz,ms2, wdth)
C$OMP END PARALLEL

      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=dble(mgz(i,j)) ; mgz(i,j)=1
          else
            wgt(i,j)=0.D0
          endif
        enddo
      enddo
#endif /* LEGACY_METHOD */

! The following part makes round junctions between two adjacent open
! boundaries by extending merging zone inward. This is useful to avoid
! steep gradient near inner corner of merging zone when weighting
! finction "wgt" is subjected to Laplacian smoothing.  Note that all
! loops below are non-reversible as they process points in row-by-row
! manner (both in i and j directions) with checking that the previous
! row is set to 1.

#ifdef SMOOTH_JUNCTIONS
      rad=(wdth+1)/3 ; rad2=rad**2
      if (OBC_NORTH .and. OBC_WEST) then
        ict=wdth+rad+1 ; jct=ncy-wdth-rad  !<-- center of the circle
        do j=jct+rad,jct,-1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i-1,j) > 0 .and. mgz(i,j+1) > 0
     &                         .and. mgz(i-1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_WEST .and. OBC_SOUTH) then
        ict=wdth+rad+1 ; jct=wdth+rad+1
        do j=jct-rad,jct,+1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i-1,j) > 0 .and. mgz(i,j-1) > 0
     &                         .and. mgz(i-1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_SOUTH .and. OBC_EAST) then
        ict=ncx-wdth-rad ; jct=wdth+rad+1  !<-- center of the circle
        do j=jct-rad,jct,+1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i+1,j) > 0 .and. mgz(i,j-1) > 0
     &                         .and. mgz(i+1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_EAST .and. OBC_NORTH) then
        ict=ncx-wdth-rad ; jct=ncy-wdth-rad
        do j=jct+rad,jct,-1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     &                   .and. mgz(i+1,j) > 0 .and. mgz(i,j+1) > 0
     &                         .and. mgz(i+1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
#endif

! Etch the merging zone into the land mask, in such a way that
! any possibility of reaching water behind mask is excluded.

C$OMP PARALLEL SHARED(ncx,ncy, mask, mgz,ms1,ms2)
      call etch_mgz_into_land_thread(ncx,ncy, mask, mgz,ms1,ms2)
C$OMP END PARALLEL

! At this stage "ms1 > 0" defines the area to which the merging zone
! should be allowed to grow. This area consists of all water points
! of the merging zone which can be reached by water starting from
! the boundary row of non-masked points, and the adjacent land area
! to which to merging zone can be expanded to ensure smoothness of
! the merged topography, but without interfering with water interior.
! The next stage is to etch the actual merging function into the
! land-masked area.

      do j=1,ncy                            ! Set "mgz" to be merging
        do i=1,ncx                          ! zone mask (positive "ms1"
          if (ms1(i,j) > 0) then            ! area), which will remain
            if (mask(i,j) > 0) then         ! unchanged from now on,
              ms1(i,j)=int(wgt(i,j)+0.5D0)  ! then initialize both
            else                            ! "ms1" and "ms2" to the
              ms1(i,j)=0                    ! constant-slope merging
            endif                           ! function (which later
            mgz(i,j)=1                      ! be used as the argument
          else                              ! for the actual merging
            ms1(i,j)=0                      ! function), using integer
            mgz(i,j)=0                      ! numbers from the range
          endif                             ! from 0 to "wdth".
          ms2(i,j)=ms1(i,j)
        enddo
      enddo

C$OMP PARALLEL SHARED(wdth, ncx,ncy, mgz,ms1,ms2, OBC_WEST,OBC_EAST,
C$OMP&                                           OBC_SOUTH,OBC_NORTH)
      call etch_weights_into_land_thread(wdth, ncx,ncy, mgz,ms1,ms2,
     &                         OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
C$OMP END PARALLEL

! At this moment if there would be no land "ms1" is constant-slope
! function: the same maximum value at the outermost row of points and
! descending linearly to zero within the band "wdth" into the domain -
! so smoothing iterations below (subject to Dirichlet B.C) would not
! actually cause any change to it with the exception of (i) in the
! vicinity of corners between adjacent open boundaries (resulting in
! sharp bend along the 45-degree row points where two side slopes
! come together), and (ii) irregular shapes caused by "etching" of
! "ms1" into the land.

      do j=1,ncy
        do i=1,ncx
          wgt(i,j)=dble(ms1(i,j)) ; hwg(i,j)=wgt(i,j)
        enddo
      enddo
      iter=4*wdth
C$OMP PARALLEL SHARED(iter, ncx,ncy, mgz,wgt,hwg)
      call smooth_wgt_thread(iter, ncx,ncy, mgz,wgt,hwg)
C$OMP END PARALLEL

! Convert the originally constant-slope --> etched --> smoothed "wgt"
! into merging function, then perform the actual merging of topography.
! Note that the incoming "wgt" is within the range of 0 <= wgt <= wdth
! nominally changing by 1 from one point to the next (nominally because
! masking-etching-smoothing may modify this). A slightly different
! rescaling,
!               -1/(2*(wdth-1)) <= wgt <= (wdth-1/2)(wdth-1)
!
! subject to limiting 0 <= wgt <= 1 is used below with the rationale
! to take out discretely jagged inner border of merging zone (obtained
! by etching), and to make cos() function a bit closer to 1 at the
! second row of points near the open boundary.
      hwg=hchd

      write(*,*) 'Merging topography...'
      cff=1.D0/dble(wdth-1)
      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=min(1.D0, cff*(wgt(i,j)-0.5D0))
            if (wgt(i,j) < 0.D0) then
              wgt(i,j)=0.D0 ; mgz(i,j)=0
            endif
            wgt(i,j)=0.5D0-0.5D0*cos(3.14159265358979323D0*wgt(i,j))
            hwg(i,j)=hwg(i,j)*(1.D0-wgt(i,j))+hc(i,j)*wgt(i,j)
          endif

          if (mgz(i,j) > 0) then
            if (mask(i,j) > 0) then
              mwgt(i,j)=1.
            else                          ! Array "mwgt" is merely
              mwgt(i,j)=0.8               ! for 4-color illustration
            endif                         ! of the layout of merging
          else                            ! zone relatively to land
            wgt(i,j)=0.D0                 ! mask:
            if (mask(i,j) > 0) then       !
              mwgt(i,j)=0.                !    mwgt=1 or 0.8 merge
            else                          !    mwgt=0 or 0.2 intact
              mwgt(i,j)=0.2               !
            endif                         !    mwgt=0  or  1 water
          endif                           !    mwgt= 0.8 or 0.2 land
        enddo
      enddo
      end subroutine
