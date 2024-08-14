! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
!                        | h(i+1) - h(i) |
!                  rx = -------------------  <  rx_max
!                          h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine.     Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
!                    Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition 
!
!      | Lgh(i+1,j) - Lgh(i+1,j) |  <  log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth".  Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
!              lsmooth hmin hmax r_max file.nc
! where
!        #1  hmin   desired minimum depth limit, meters, positive;
!        #2  hmax   --/ /-- maximum        --/ /--
!        #3  r_max  desired maximum r-factor, nondimensional;
!        #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.

# ifdef LOG_SMOOTHING
      subroutine lsmooth (Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
#elif defined LEGACY_LOG_SMOOTH
      subroutine lsmooth_legacy(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# elif defined LOG_SMOOTH_2
      subroutine lsmooth2(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# elif defined LOG_SMOOTH_1
      subroutine lsmooth1(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
#elif defined NON_DECREASING_RX_COND
      subroutine cond_rx0_topo(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# else
      subroutine smooth (Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# endif

      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2) 
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0
      
!      call lenstr(grid,lgrd)

!      ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )
      
!      call get_var_by_name_double(ncid, 'mask_rho', wrk)
      
      do j=1,ncy
        do i=1,ncx                      ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then    ! to read mask from the file;
            msk(i,j)=1                  ! from now on this code will
          else                          ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j)                      !<-- make another copy
        enddo
      enddo

#if defined LOG_SMOOTHING
      method='LOG_SMOOTHING'
#elif defined LEGACY_LOG_SMOOTH
      method='LEGACY_LOG_SMOOTH'
# elif defined LOG_SMOOTH_2
      method='LOG_SMOOTH_2'
# elif defined LOG_SMOOTH_1
      method='LOG_SMOOTH_1'
#elif defined NON_DECREASING_RX_COND
      method='RX_COND'
#else
      method='SMOOTH'
#endif

# ifndef NON_DECREASING_RX_COND
C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL
# else
C$OMP PARALLEL SHARED(ncx,ncy, msk, h,wrk)
      call cond_rx0_thread(ncx,ncy, msk, h,hmin,hmax,r_max,wrk)
C$OMP END PARALLEL
# endif

!      call write_topo(ncid, grid, ncx,ncy, h)
      end

