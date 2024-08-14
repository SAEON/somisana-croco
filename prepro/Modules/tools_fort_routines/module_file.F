      module smooth_topo_params
       character(len=64) method
       real(kind=8) :: hmin, hmax, r_max
       integer :: iters_cond=0, iters_lin=0
      end

      module rx_comm_vars
      real(kind=8) :: rx_max=0.D0, ry_max=0.D0
      integer :: trd_count=0, ifrst_call=1,
     &    irx_max=0, jrx_max=0, iry_max=0, jry_max=0
      end module rx_comm_vars


! The following module is designed to be completely initialized by
! "read_roms_grid", which includes both allocation arrays with proper
! dimensions matching the actual grid file and filling them with data.
! Note that "angle" is no longer part of the module because all what
! is needed in most cases when grid file is read is cos and sin of
! angle to rotate vector componets, but not angle itself.

      module roms_grid_vars
       real(kind=8), allocatable, dimension(:,:) :: lon_r,lat_r,
     &            lon_p,lat_p, pm,pn, f, rmask, csA,snA, orterr
       logical curv_grid
      end module roms_grid_vars

! The following module contains parameter for grid generation and is
! shared by only by "read_grid_params.F" (where they are initialized),
! "grid.F", and "read_roms_grid" (later in this file, used only as
! input to leave signature in netCDF file).

      module roms_grid_params
       real(kind=8) lon_min,lon_max, lat_min,lat_max,
     &    size_x,size_y, cent_lat, psi0,theta0,alpha
       integer nx,ny, flip_xy
      end module

! This module is to monitor the total amount of data read and written.

      module mod_io_size_acct
       integer(kind=8), save :: read_clk=0,  sz_read_acc=0,
     &                          write_clk=0, sz_write_acc=0
      end module mod_io_size_acct



!!!! Module r2r_vert_int !!!!

      module comp_uvbar_vars
      integer, save :: alloc_zw_size=0
C$OMP THREADPRIVATE(alloc_zw_size)
      real(kind=8), allocatable, dimension(:) :: D,z_w,Hz
C$OMP THREADPRIVATE(D,z_w,Hz)
      end module comp_uvbar_vars

      module r2r_vertint_vars
       integer, save :: trd_count=0
       real(kind=8), save :: vert_int_error=0.
       integer, save :: alloc_zc_size=0, allc_zpr_size=0
C$OMP THREADPRIVATE(allc_zpr_size, alloc_zc_size)
       real(kind=8), allocatable, dimension(:) :: zc,zp_r,drv
C$OMP THREADPRIVATE(zc,zp_r,drv)
      end module r2r_vertint_vars

