      subroutine srtopo(Lm,Mm,srcdir,lon_r,lat_r,pm,pn,radius,hraw)

! Purpose: Read "raw" topography from SRTM30 dataset -- a directory
!          containing 33 tiles in netCDF format -- and interpolate or
!          coarsen it onto ROMS model grid. The outcome is written as
!          netCDF variable "hraw".

! NOTES: (1) the content of this file is just a driver: the actual
!            interpolation/coarsening routine is in "compute_hraw.F".

!        (2) nothing is done here about setting mask, which needs to
!            be processed separately via "rmask.F". This is  to avoid
!            accidental overwriting of land mask field "rmask" which
!            may be subject to elaborate hand editing.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

!      use roms_grid_vars
!      use comm_vars_hraw

      implicit none
      real(kind=8) radius
      integer nx_alc,ny_alc, nx,ny,  nx_lon,ny_lat, nxy_pts,
     &        nargs, Lm,Mm, ncid, ndims, varid, tile, id, size,
     &        i,j, ierr, lsrc, lgrd, ldim,llat,llon
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: lon_r,lat_r,pm,pn
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: hraw
      real(kind=8), allocatable :: x(:),y(:), xlon(:),ylat(:)
      real(kind=8), allocatable :: htopo(:,:)
      integer(kind=2), allocatable :: bffr(:)
!! roms_grid_vars
      logical curv_grid
!! comm_vars_hraw
      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

      integer, parameter :: ntiles=33
      character(len=7), dimension(ntiles) :: srtile = (/
     &    'w180s60','w120s60','w060s60','w000s60','e060s60','e120s60',
     &    'w180s10','w140s10','w100s10','w060s10','w020s10',
     &                        'e020s10','e060s10','e100s10','e140s10',
     &    'w180n40','w140n40','w100n40','w060n40','w020n40',
     &                        'e020n40','e060n40','e100n40','e140n40',
     &    'w180n90','w140n90','w100n90','w060n90','w020n90',
     &                       'e020n90', 'e060n90','e100n90','e140n90'/)

      character(len=7) curr_srtile
      real(kind=8), dimension(ntiles) :: xW,xE, yS,yN
      real(kind=8)  xmin,xmax,ymin,ymax
      integer, dimension(ntiles) :: tndx, istr,iend,jstr,jend
      integer, dimension(ntiles) :: ewshft, iwest,ieast,jsouth
      integer, dimension(ntiles) :: jnorth
      integer :: ntls, ishft,jshft,imin,imax,jmin,jmax

      character(len=160) :: srcdir
      character(len=16) :: str, dname, lon_name, lat_name
      real(kind=8), parameter :: spv=-99999.D0

      include "netcdf.inc"

cf2py intent(in) Lm,Mm,srcdir,lon_r,lat_r,pm,pn,radius
cf2py intent(out) hraw


      comm_j1_indx=0
      comm_j2_indx=0

      call lenstr(srcdir,lsrc)
      if (srcdir(lsrc:lsrc) /= '/') then
        lsrc=lsrc+1
        srcdir(lsrc:lsrc)='/'
      endif


! Reset everything:       ! "nx,ny_alc" are the actual allocated sizes
! -----------------       ! of tile-sized buffer arrays [they will grow
      nx_alc=0 ; ny_alc=0 ! if the tiles are of different sizes]; 
      ntls=0              ! number of active tiles;
      do tile=1,ntiles
        tndx(tile)=0
        iwest(tile)=0  ; ieast(tile)=0
        jsouth(tile)=0 ; jnorth(tile)=0
      enddo

! Open and read ROMS grid file. Find geographical limits of the grid.

      call roms_grid_geo_bounds(lon_r,lat_r, Lm,Mm, radius,
     &           lon_min,lon_max,lat_min,lat_max)

      write(*,'(1x,2A,2F16.8/22x,A,2F16.8 )')  'roms grid extremes: ',
     &   'longitude:', lon_min, lon_max, 'latitude:', lat_min,lat_max

! Open and scan topography data files.

      do tile=1,ntiles
        ierr=nf_open(srcdir(1:lsrc)/ /srtile(tile)/ /'.nc',
     &                                   nf_nowrite, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_ndims (ncid, ndims)
          if (ierr == nf_noerr) then
            nx=0 ; ny=0
            do id=1,ndims
              dname='       '
              ierr=nf_inq_dim (ncid, id, dname, size)
              if (ierr == nf_noerr) then
                call lenstr(dname,ldim)
                if( (ldim==1 .and. dname(1:ldim)=='x')  .or.
     &              (ldim==3 .and. dname(1:ldim)=='nx')  .or.
     &              (ldim==9 .and. dname(1:ldim)=='longitude') ) then
                  lon_name=dname(1:ldim)
                  nx=size
                elseif( (ldim==1 .and.dname(1:ldim)=='y') .or.
     &                  (ldim==3 .and.dname(1:ldim)=='ny') .or.
     &                  (ldim==8 .and.dname(1:ldim)=='latitude') ) then
                  lat_name=dname(1:ldim)
                  ny=size
                endif
              else
                write(*,*) '### ERROR: dimension id =', id, '?'
              endif
            enddo
          else
            write(*,*) '### ERROR: nf_inq_ndims?'
          endif
!#ifdef VERBOSE
!          write(*,'(1x,A,3(1x,A),I6,2x,2(1x,A),I6)') srtile(tile),
!     &    ': dimensions:',  lon_name(1:lenstr(lon_name)), '=', nx,
!     &                      lat_name(1:lenstr(lat_name)), '=', ny
!#endif

          if (nx > 0 .and. ny > 0) then
            if (nx > nx_alc) then              ! NOTE: Since the data
              if (allocated(x)) deallocate(x)  ! files are expected to
              allocate(x(nx))                  ! be what is called
              nx_alc=nx                        ! "CF-compliant" its
            endif                              ! coordinate variables
                                               ! should have the same
            if (ny > ny_alc) then              ! names as their
              if (allocated(y)) deallocate(y)  ! corresponding
              allocate(y(ny))                  ! dimensions.
              ny_alc=ny
            endif

            ierr=nf_inq_varid (ncid, lon_name, varid)
            if (ierr == nf_noerr) then
              ierr=nf_get_var_double (ncid, varid, x)
              if (ierr == nf_noerr) then
                ierr=nf_inq_varid (ncid, lat_name, varid)
                if (ierr == nf_noerr) then
                  ierr=nf_get_var_double (ncid, varid, y)
                  if (ierr == nf_noerr) then

! Note that fiction "indx_bound" returns 0 on "ny" if the test value
! "lon_mix/max" is outside the range of "y" for this tile, so the first
! if-condition indicates that that there is an overlap in y-coordinate
! (latitude of data) between ROMS grid and the tile.

                    jmin=indx_bound3(y,ny, lat_min)
                    jmax=indx_bound3(y,ny, lat_max)
                    if (jmin<ny .and. jmax>0) then
                      if (jmin == 0) jmin=1 
                      if (jmin > 1 ) jmin=jmin-1
                      if (jmax < ny) jmax=jmax+1

! Longitude coordinate is defined with 360 degree periodicity, so if
! the tile dies not fit right a way, try to shift it east of west by
! 360 degrees and then check again whether it has overlap. 

                      ishft=0 
                      imin=indx_bound3(x,nx, lon_min)
                      imax=indx_bound3(x,nx, lon_max)
                      if (imin==nx) then
                        ishft=+1
                        imin=indx_bound3(x,nx, lon_min-360.D0)
                        imax=indx_bound3(x,nx, lon_max-360.D0)
                      elseif (imax==0) then
                        ishft=-1
                        imin=indx_bound3(x,nx, lon_min+360.D0)
                        imax=indx_bound3(x,nx, lon_max+360.D0)
                      endif
                      if (imin<nx .and. imax>0) then 
                        if (imin==0) imin=1
                        if (imin >1) imin=imin-1
                        if (imax<nx) imax=imax+1

! Once it passed all the logical check above, add the tile to the list
! of tiles to be read, and record all its attributes.

                        ntls=ntls+1
                        tndx(ntls)=tile   ; ewshft(ntls)=ishft
                        iwest(ntls)=imin  ; jsouth(ntls)=jmin
                        ieast(ntls)=imax  ; jnorth(ntls)=jmax

                        if (ishft > 0) then
                          xW(ntls)=x(imin) +360.D0
                          xE(ntls)=x(imax) +360.D0
                        elseif (ishft < 0) then
                          xW(ntls)=x(imin) -360.D0
                          xE(ntls)=x(imax) -360.D0
                        else
                          xW(ntls)=x(imin)
                          xE(ntls)=x(imax)
                        endif
                        yS(ntls)=y(jsouth(ntls))
                        yN(ntls)=y(jnorth(ntls))

                      endif    !<-- imin<nx .and. imax>0
                    endif    !<--   jmin<ny .and. jmax>0

                  else
                    call lenstr(lat_name,llat)
                    write(*,'(/1x,6A/12x,A/)')    '### ERROR: Cannot ',
     &                'read variable ''', lat_name(1:llat),
     &                            ''' from ''', srtile(tile), '.nc''.',
     &                                               nf_strerror(ierr)
                  endif
                else
                  call lenstr(lat_name,llat)
                  write(*,'(/1x,6A/12x,A/)')  '### ERROR: Cannot get ',
     &                                     'netCDF variable ID for ''',
     &                     lat_name(1:llat), ''' from  ''',
     &                      srtile(tile), '.nc''.',  nf_strerror(ierr)
                endif
              else
                call lenstr(lon_name,llon)
                write(*,'(/1x,6A/12x,A/)')   '### ERROR: Cannot read ',
     &                     'variable ''', lon_name(1:llon),
     &         ''' from ''', srtile(tile), '.nc''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(/1x,6A/12x,A/)')      '### ERROR: Cannot get ',
     &       'netCDF variable ID for ''', lon_name(1:llon),
     &        ''' from ''',  srtile(tile), '.nc''.', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/)')        '### ERROR: Cannot determine ',
     &                   'dimension sizes for topography data file ''',
     &                                          srtile(tile), '.nc''.'
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot open netCDF ',
     &             'file ''', srtile(tile), '.nc'' in read-only mode.',
     &                                               nf_strerror(ierr)
        endif
        ierr=nf_close(ncid)
      enddo
      if (allocated(x)) deallocate(x)
      if (allocated(y)) deallocate(y)

#define i illegal

! Summarize the results of tile scan.


      write(*,'(/ /1x,2A,I3,3x,A/1x,2A/)')      'SUMMARY: number of ',
     &  'selected tiles, ntls =', ntls,  'Their file names, bounding',
     &  'indices of actually used portions, and geographical limits ',
     &  'of used portions:'

      write(*,'(2x,A,1x,A,3x,A,1x,A,1x,A,1x,A,2x,A,11x,A,10x,A,9x,A)')
     &  '#', 'file', 'iwst','iest','jsth','jnrth', 'xW','xE','yS','yN'
      write(*,'(1x,2A)') '---------------------------------------',
     &                   '---------------------------------------'
      do tile=1,ntls
        write(*,'(I3,1x,A,4I5,2F13.7,2F11.7)') tile, srtile(tndx(tile)),
     &               iwest(tile),ieast(tile),jsouth(tile),jnorth(tile),
     &             xW(tile),xE(tile),yS(tile),yN(tile)
      enddo

      xmin=xW(1) ; xmax=xE(1)
      ymin=YS(1) ; ymax=yN(1)
      do tile=2,ntls
        if (xW(tile)<xmin) xmin=xW(tile)
        if (xE(tile)>xmax) xmax=xE(tile)
        if (yS(tile)<ymin) ymin=yS(tile)
        if (yN(tile)>ymax) ymax=yN(tile)
      enddo

      write(*,'(/1x,2A,/1x,A,2F13.8,2x,A,2F12.8/)')  'Geographical ',
     &              'extremes of the used portions among all tiles:', 
     &               'xmin,max =',xmin,xmax, 'ymin,max =',ymin,ymax


! The following code segment determines bounding indices of each tile
! (actually used portion thereof) within the index coordinates of would
! be global data array.  This is needed for correct placement of data
! from each individual tile.
! The mapping is as follows:
!                              iwest:ieast  --> istr:iend
!                             jsouth:jnorth --> jstr:jend

! The algorithm is essentially to convert geographical limits of each
! used tile portion into integer numbers implicitly assuming that the
! grid spacing is globally uniform among all the tiles. [Note that 0.5
! is added inside each int() to counter roundoff error - mathematically
! the result of division (ieast-iwest)/(xE-xW) yields the inverse of
! grid spacing while (xW)-xmin) should be an integer number of grid
! spaces, so the whole expression inside int() without 0.5 should be
! an integer number. However  minute roundoff error may cause it to
! be 1 less than it should.]   
   
! The dimensions of global data arrays are defined as the maximum
! iend and jend among all the tiles.

      nxy_pts=0 ; nx_lon=0 ; ny_lat=0

      do tile=1,ntls
        istr(tile)=1+int( 0.5D0 + dble( ieast(tile)-iwest(tile) )
     &                      *(xW(tile)-xmin)/(xE(tile)-xW(tile)) )
        iend(tile)=istr(tile) + ieast(tile)-iwest(tile)
        if (iend(tile)>nx_lon) nx_lon=iend(tile)

        jstr(tile)=1+int( 0.5D0 + dble( jnorth(tile)-jsouth(tile) )
     &                       *(yS(tile)-ymin)/(yN(tile)-yS(tile)) )
        jend(tile)=jstr(tile) +jnorth(tile)-jsouth(tile)
        if (jend(tile)>ny_lat) ny_lat=jend(tile)

        nxy_pts=nxy_pts + (iend(tile)-istr(tile)+1)
     &                   *(jend(tile)-jstr(tile)+1)
      enddo

      write(*,'(/2(1x,A,I7)/15x,A,I12/42x,A,I12/66x,A/47x,A,I12/)')
     &   'Required dimensions of array to hold the entire used data:',
     &                                            nx_lon, 'x', ny_lat,
     &   'Aggregate number of points in all selected tiles:', nxy_pts,
     &                        'product of dimensions:', nx_lon*ny_lat,
     &   '-----------',    'their difference:', nxy_pts-nx_lon*ny_lat


      write(*,'(/1x,2A/1x,A/)') 'Starting and ending indices for ',
     &                   'the used portions of each tile as defined',
     &         'within the logical coordinates of global data array:'

      write(*,'(10x,A,1x,A,2x,A/4x,2A)') '#', 'file', 'ew_shft',
     &                   '------------------------------------',
     &                   '------------------------------------'
      do tile=1,ntls
        write(*,'(8x,I3,1x,A,1x,I3,1x,A,2I7,2x,A,2I7)') tile,
     &     srtile(tndx(tile)), ewshft(tile),
     &                      'istr,iend =', istr(tile),iend(tile),
     &                      'jstr,jend =', jstr(tile),jend(tile)
      enddo

      if (nxy_pts < nx_lon*ny_lat) then
        write(*,'(/1x,2A/)')   '### ERROR: Available data tiles do ',
     &                             'not cover the entire ROMS grid.'
        stop
      elseif (nxy_pts > nx_lon*ny_lat) then
        write(*,'(/1x,2A/)') '### ERROR: Possible overlapping tiles ',
     &                         'or tile selection algorithm failure.'
        stop
      endif
#undef i

      write(*,'(/1x,A/)') 'Reading topographic data..'

! Note that the tile-sized coordinate arrays x,y were deallocated above
! and are allocated again just below, but after the arrays covering the
! entire grid area. This is to avoid memory fragmentation by keeping 
! the tile-sized x,y, and bffr at the end of allocated memory, are they 
! will be deallocated as soon as reading of tiled data is complete, and
! a new array to hold interpolated topography will be allocated.

      allocate(xlon(nx_lon))  ;  allocate(ylat(ny_lat))
      allocate(htopo(nx_lon,ny_lat))

      allocate(x(nx_alc))     ;  allocate(y(ny_alc))
      allocate(bffr(nx_alc*ny_alc))

      do i=1,nx_lon     ! Initialize coordinate arrays with a special
        xlon(i)=spv     ! value. This is needed to check consistency of
      enddo             ! coordinates stored in different tiles.  Each
      do j=1,ny_lat     ! coordinate value will be recorded only once,
        ylat(j)=spv     ! and thereafter checked that the value from a
      enddo             ! different tile matches the already recorded.

      do tile=1,ntls
        curr_srtile=srtile(tndx(tile))
        ierr=nf_open(srcdir(1:lsrc)/ /curr_srtile/ /'.nc',
     &                                  nf_nowrite, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_dimid (ncid, lon_name, id)
          if (ierr == nf_noerr) then
            ierr=nf_inq_dimlen (ncid, id, nx)
            if (ierr == nf_noerr) then
              ierr=nf_inq_dimid (ncid, lat_name, id)
              if (ierr == nf_noerr) then
                ierr=nf_inq_dimlen (ncid, id, ny)
                if (ierr == nf_noerr) then
                  call lenstr(lon_name,llon)
                  call lenstr(lat_name,llat)
                  write(*,'(I3,4(1x,A),I6,2x,2(1x,A),I6)',advance='no')
     &                           tile,   curr_srtile,   'dimensions:',
     &                           lon_name(1:llon), '=', nx,
     &                           lat_name(1:llat), '=', ny
                else
                  call lenstr(lat_name,llat)
                  write(*,'(1x,2A,I3,1x,5A/12x,A)')      '### ERROR: ',
     &                     'Cannot determine length of dimension ', id,
     &                        'named ''', lat_name(1:llat),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
                endif
              else
                call lenstr(lat_name,llat)
                write(*,'(1x,6A/12x,A)')      '### ERROR: Cannot get ',
     &             'dimension ID for ''', lat_name(1:llat),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(1x,2A,I3,1x,5A/12x,A)')   '### ERROR: Cannot ',
     &                'determine length of dimension ', id, 'named ''',
     &                      lon_name(1:llon), ''' from ''',
     &                           curr_srtile, '''.', nf_strerror(ierr)
            endif
          else
            call lenstr(lon_name,llon)
            write(*,'(1x,6A/12x,A)')          '### ERROR: Cannot get ',
     &             'dimension ID for ''', lon_name(1:llon),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
          endif

          if (ierr == nf_noerr) then
            ierr=nf_inq_varid (ncid, lon_name, varid)
            if (ierr == nf_noerr) then
              ierr=nf_get_var_double (ncid, varid, x)
              if (ierr == nf_noerr) then
                ierr=nf_inq_varid (ncid, lat_name, varid)
                if (ierr == nf_noerr) then
                  ierr=nf_get_var_double (ncid, varid, y)
                  if (ierr /= nf_noerr) then
                    call lenstr(lat_name,llat)
                    write(*,'(1x,2A,I3,1x,5A/12x,A)')    '### ERROR: ',
     &                       'Cannot read coordinate variable ', varid,
     &                        'named ''', lat_name(1:llat),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)

#ifdef VERBOSE
                  else
                    write(*,'(1x,A)',advance='no') 
     &                  'read coordinate variables'
#endif
                  endif
                else
                  call lenstr(lat_name,llat)
                  write(*,'(1x,6A/12x,A)')    '### ERROR: Cannot get ',
     &              'variable ID for ''', lat_name(1:llat),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
                endif
              else
                call lenstr(lon_name,llon)
                write(*,'(1x,2A,I3,1x,5A/12x,A)') '### ERROR: Cannot ',
     &                  'read coordinate variable ', varid, 'named ''',
     &                      lon_name(1:llon), ''' from ''',
     &                           curr_srtile, '''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(1x,6A/12x,A)')  '### ERROR: Cannot determine ',
     &              'variable ID for ''', lon_name(1:llon),
     &             ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
            endif
          endif

          if (ierr == nf_noerr) then
            ierr=nf_inq_varid (ncid, 'topo', varid)
            if (ierr == nf_noerr) then
              write(*,'(12x,A)') 'found variable ''topo''.'
            else
              ierr=nf_inq_varid (ncid, 'elevation', varid)
              if (ierr == nf_noerr) then
                write(*,'(12x,A)') 'found variable ''elevation''.'
              else
                ierr=nf_inq_varid (ncid, 'z', varid)
#ifdef VERBOSE
                if (ierr == nf_noerr) write(*,'(1x,A)',advance='no')
     &                                           'found topography'
#endif
              endif
            endif
          endif

          if (ierr == nf_noerr) then
            ierr=nf_get_var_int2 (ncid, varid, bffr)
            if (ierr == nf_noerr) then
              write(*,'(1x,A)')  'retrieved topographic data'
            else
              write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot read ',
     &                  'topographic data from ''', curr_srtile, '''.',
     &                                              nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot determine ',
     &                 'netCDF variable ID for topography field in ''',
     &                          curr_srtile, '''.',  nf_strerror(ierr)
          endif

          if (ierr == nf_noerr) then
            ishft=istr(tile)-iwest(tile)
            if (ewshft(tile)>0) then
              do i=iwest(tile),ieast(tile)
                x(i)=x(i) +360.D0
              enddo
            elseif (ewshft(tile)<0) then
              do i=iwest(tile),ieast(tile)
                x(i)=x(i) -360.D0
              enddo
            endif
            ierr=0
            do i=iwest(tile),ieast(tile)
              if (xlon(i+ishft) > spv) then
                if (xlon(i+ishft) /= x(i)) ierr=ierr+1
              else
                xlon(i+ishft)=x(i)
              endif
            enddo
            if (ierr>0) write(*,'(/1x,2A/)') '### ERROR: Conflicting ',
     &                      'longitude coordinate data between tiles.'

            jshft=jstr(tile)-jsouth(tile)
            ierr=0
            do j=jsouth(tile),jnorth(tile)
              if (ylat(j+jshft) > spv) then
                if (ylat(j+jshft) /= y(j)) ierr=ierr+1
              else
                ylat(j+jshft)=y(j)
              endif
            enddo
            if (ierr>0) write(*,'(/1x,2A/)') '### ERROR: Conflicting ',
     &                       'latitude coordinate data between tiles.'


            do j=jsouth(tile),jnorth(tile)
              do i=iwest(tile),ieast(tile)
                htopo(i+ishft,j+jshft)=bffr(i+nx*(j-1))
              enddo
            enddo

#ifdef VISUALIZE_TILING
            do j=jsouth(tile),jnorth(tile)
              htopo(iwest(tile)+ishft,j+jshft)=1.0E+6
              htopo(ieast(tile)+ishft,j+jshft)=1.0E+6
            enddo
            do i=iwest(tile),ieast(tile)
              htopo(i+ishft,jsouth(tile)+jshft)=1.0E+6
              htopo(i+ishft,jnorth(tile)+jshft)=1.0E+6
            enddo
#endif

          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot open netCDF ',
     &             'file ''', srtile(tile), '.nc'' in read-only mode.',
     &                                               nf_strerror(ierr)
        endif
        ierr=nf_close(ncid)
      enddo

      deallocate(bffr) ; deallocate(y) ; deallocate(x)

! 3. Once both the input topography data and the target roms grid
! file are successfully opened and the relevant data is read, do the
! actual interpolation:

      if (ierr == nf_noerr) then
!        allocate (hraw(0:Lm+1,0:Mm+1))  ; ierr=0

C$OMP PARALLEL SHARED(nx_lon,ny_lat, xlon,ylat, htopo, Lm,Mm,
C$OMP&                        lon_r,lat_r, pm,pn, hraw, radius, ierr)
        call compute_hraw(nx_lon,ny_lat,  xlon,ylat,  htopo,
     &                 Lm,Mm, lon_r,lat_r, pm,pn, hraw, radius)
C$OMP END PARALLEL
!        if (ierr==0) call write_hraw(grid, srcdir,radius, Lm,Mm,hraw)
      endif

      contains
        integer function indx_bound3(x,n, x0)
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
        indx_bound3=i
        end function indx_bound3
  

      end subroutine srtopo
