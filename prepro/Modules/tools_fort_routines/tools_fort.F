!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROMS ROUTINES
!!
!! copied from actual ROMS scripts
!!
!! compile with:
!! "cpp R_tools_fort.F R_tools_fort.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort R_tools_fort.f" for python use
!!
!! print R_tools_fort.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# include "sigma_to_z_intr.F"
# include "sigma_to_z_intr_bot.F"
# include "sigma_to_z_intr_bot_2d.F"
# include "sigma_to_z_intr_bounded.F"
# include "sigma_to_z_intr_sfc.F"
# include "single_connect.F"
!!! diverse !!!
# include "roms_read_write.F"
# include "srtopo.F"
!!! compute hraw !!!
# include "compute_hraw.F"

!!! smooth !!!
# include "topo_smooth_subs.F"
# include "tools_topo.F"
# include "smooth.F"

!!! lsmooth !!!
# define LOG_SMOOTHING
# include "smooth.F"
# undef LOG_SMOOTHING

!!! legacy_lsmooth !!!
# define LEGACY_LOG_SMOOTH
# include "smooth.F"
# undef LEGACY_LOG_SMOOTH

!!! lsmooth2 !!!
# define LOG_SMOOTH_2
# include "smooth.F"
# undef LOG_SMOOTH_2

!!! lsmooth1 !!!
# define LOG_SMOOTH_1
# include "smooth.F"
# undef LOG_SMOOTH_1

!!! non_decreasing_rx_cond !!!
# define NON_DECREASING_RX_COND
# include "cond_rx0_topo.F"
# include "smooth.F"
# undef NON_DECREASING_RX_COND

# include "r2r_match_topo.F"


# include "r2r_init.F"
# include "r2r_bry.F"
# include "etch_into_land.F"
# include "set_depth.F"
# include "def_roms_file.F"
# include "def_bry_file.F"
# include "r2r_subs.F"
# include "r2r_vert_interp.F"
# include "r2r_rotate.F"
# include "r2r_bry_rotate.F"

!!! tools used for smoothing !!!
# include "smooth_utils.F"
# include "tiling.F"
# include "lenstr.F"
# include "mrg_zone_subs.F"
# include "r2r_interp_init.F"
# include "r2r_bry_interp.F"
!!!!!!!!!!

