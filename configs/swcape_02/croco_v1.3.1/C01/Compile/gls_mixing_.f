      SUBROUTINE gls_mixing (tile)
      IMPLICIT NONE
      INTEGER*4         :: tile, trd
      INTEGER*4         :: omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=128,  MMm0=446,  N=30)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=2,  NP_ETA=8,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=5)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,7,0:NPP-1)
      common/private_scratch/ A2d,A3d
      real u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hz_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Huon(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hvom(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Wi(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real rho(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=omp_get_thread_num()
      call gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                    A3d(1, 1,trd), A3d(1, 2,trd), A3d(1, 3,trd),
     &                    A3d(1, 4,trd), A3d(1, 5,trd), A2d(1, 2,trd),
     &                    A2d(1, 3,trd), A2d(1, 4,trd), A2d(1, 5,trd),
     &                    A2d(1, 6,trd), A2d(1, 7,trd), A2d(1, 8,trd))
      RETURN
      END
      SUBROUTINE gls_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                          tke_new, gls_new, shear2, vort2, work,
     &                                diss, ustar_sfc_sq,ustar_bot_sq,
     &                                                 DC, FC, CF, RH)
      IMPLICIT NONE
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=128,  MMm0=446,  N=30)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=2,  NP_ETA=8,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=5)
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real diff4_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff4(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff4_sponge/diff4_sponge
      common /mixing_diff4/diff4
      real diff3d_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real diff3d_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /mixing_diff3d_u/diff3d_u
      common /mixing_diff3d_v/diff3d_v
      real dRdx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real dRde(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real idRz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_dRdx/dRdx
      common /mixing_dRde/dRde
      common /mixing_idRz/idRz
      real Rslope_max,Gslope_max
      parameter (Gslope_max=1.D0, Rslope_max=0.05D0)
      integer*4 ismooth
      real csmooth
      common /mixing_csmooth/ csmooth
      common /mixing_ismooth/ ismooth
      real Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real Akv_old(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt_old(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_Akvold/Akv_old /mixing_Aktold/Akt_old
      real bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real trb(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2,NGLS)
      common /gls_trb/trb
      real Lscale(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /gls_lsc/Lscale
      real Eps_gls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /gls_eps/Eps_gls
      integer*4 kbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /gls_kbl/ kbl
      real hbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E  )
      common /gls_hbl/ hbl
      real cm0
      common /gls_cm0/ cm0
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      integer*4 nwrtsurf
      integer*4 ntssurf_avg, nwrtsurf_avg
      logical ldefhis
      logical got_tini(NT)
      logical ldefsurf
      logical ldefsurf_avg
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefsurf, nwrtsurf
     &                      , ldefsurf_avg
     &                      , nwrtsurf_avg
     &                      , ntssurf_avg
     &                      , ldefhis
      real Akv_bak
      common /scalars_akv/ Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akt_bak
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      INTEGER*4         ::   Istr, Iend, Jstr, Jend
      INTEGER*4         ::   i,       j,    k, tind, kref
      INTEGER*4         ::   imin, imax, jmin, jmax
      INTEGER*4         ::   ig,    ig1,  ig2
      REAL            ::  tke_new     (Istr-2:Iend+2,Jstr-2:Jend+2,0:N  
     &                                                                 )
      REAL            ::  gls_new     (Istr-2:Iend+2,Jstr-2:Jend+2,0:N  
     &                                                                 )
      REAL            ::  work        (Istr-2:Iend+2,Jstr-2:Jend+2,0:N  
     &                                                                 )
      REAL            ::  shear2      
     &                               (Istr-2:Iend+2,Jstr-2:Jend+2,0:N-1)
      REAL            ::  vort2       
     &                               (Istr-2:Iend+2,Jstr-2:Jend+2,0:N-1)
      REAL            ::  diss        (Istr-2:Iend+2,1:N-1)
      REAL            ::  ustar_sfc_sq(Istr-2:Iend+2,Jstr-2:Jend+2      
     &                                                                 )
      REAL            ::  ustar_bot_sq(Istr-2:Iend+2,Jstr-2:Jend+2      
     &                                                                 )
      REAL            ::  DC          (Istr-2:Iend+2,0:N  )
      REAL            ::  FC          (Istr-2:Iend+2,0:N  )
      REAL            ::  CF          (Istr-2:Iend+2,1:N-1)
      REAL            ::  RH          (Istr-2:Iend+2,1:N-1)
      REAL            ::  cff , cff1 , cff2, cff3m, cff3p
      REAL            ::  invk, invG, Bprod, Sprod, epsilon
      REAL            ::  alpha_n, alpha_m, c_mu, c_mu_prim
      REAL            ::  alpha_n_min, alpha_m_max, cm0inv2, gls
      REAL            ::  flux_top, flux_bot, lgthsc, L_lim, du,dv,dw
      REAL            ::  trb_new, trb_sfc, trb_bot, z0_s, z0_b, gls_min
      REAL            ::  HUon_w, HVom_w, trb_min(2), Denom
      REAL            ::  su_r,sv_r
      REAL, PARAMETER ::  eps_min =  1.0D-12
      REAL, PARAMETER ::  tke_min =  1.0D-10
      REAL, PARAMETER ::  eps     =  1.0D-14
      REAL, PARAMETER ::  nuws    =  0.1D-05
      REAL, PARAMETER ::  nuwm    =  1.0D-05
      REAL, PARAMETER ::  galp    =  0.53D0
      REAL, PARAMETER ::  chk     =  1400.D0/g
      REAL, PARAMETER ::  Zosmin  =  1.D-2
      REAL, PARAMETER ::  Zobmin  =  1.D-4
      real            :: rp,    rm,    rn
      real            :: beta1, beta2, beta3m, beta3p
      real            :: OneOverSig(2)
      parameter( rp    = 3.0D0 , rm    = 1.5D0 , rn     = -1.0D0        
     &                                                                 )
      parameter( beta1 = 1.44D0, beta2 = 1.92D0, beta3m = -0.4D0, 
     &                                                   beta3p = 1.0D0)
      parameter( OneOverSig = (/ 1.0D0, 0.7692D0 /) )
      REAL, PARAMETER :: e1 =  3.0D0 + 1.D0*rp / rn
      REAL, PARAMETER :: e2 =  1.5D0 + 1.D0*rm / rn
      REAL, PARAMETER :: e3 = -1.0D0 / rn
      REAL, PARAMETER :: smth_a = 1.D0/12.D0
      REAL, PARAMETER :: smth_b = 3.D0/16.D0
       REAL ::  c1   ,c2    ,c3    ,c4    ,c5    , c6
       REAL :: cb1   ,cb2   ,cb3   ,cb4   ,cb5   ,cbb
       REAL :: a1    ,a2    ,a3    ,a5    ,nn
       REAL :: ab1   ,ab2   ,ab3   ,ab5   ,nb
       REAL :: sf_d0 ,sf_d1 ,sf_d2 ,sf_d3 ,sf_d4 , sf_d5
       REAL :: sf_n0 ,sf_n1 ,sf_n2
       REAL :: sf_nb0,sf_nb1,sf_nb2
       REAL :: lim_am0,lim_am1,lim_am2,lim_am3,lim_am4,lim_am5,lim_am6
       PARAMETER(c1=5.D0   ,
     &           c2=0.8D0  ,
     &           c3=1.968D0,
     &           c4=1.136D0,
     &           c5=0.D0   ,
     &           c6=0.4D0   )
       PARAMETER(cb1=5.95D0  ,
     &           cb2=0.6D0   ,
     &           cb3=1.D0    ,
     &           cb4=0.D0    ,
     &           cb5=0.3333D0,
     &           cbb=0.72D0   )
       PARAMETER(  a1 = 0.66666666667D0 - 0.5D0*c2 )
       PARAMETER(  a2 = 1.D0            - 0.5D0*c3 )
       PARAMETER(  a3 = 1.D0            - 0.5D0*c4 )
       PARAMETER(  a5 = 0.5D0           - 0.5D0*c6 )
       PARAMETER( ab1 = 1.D0 - cb2               )
       PARAMETER( ab2 = 1.D0 - cb3               )
       PARAMETER( ab3 = 2.D0*(1.D0-cb4)            )
       PARAMETER( ab5 = 2.D0*cbb*(1.D0-cb5)        )
       PARAMETER( nn  = 0.5D0*c1                 )
       PARAMETER( nb  = cb1                    )
       PARAMETER( sf_d0 = 36.0D0*nn*nn*nn*nb*nb                         
     &                                                                 )
       PARAMETER( sf_d1 = 84.0D0*a5*ab3*nn*nn*nb+36.0D0*ab5*nn*nn*nn*nb 
     &                                                                 )
       PARAMETER( sf_d2 = 9.0D0*(ab2*ab2-ab1*ab1)*nn*nn*nn
     &                  - 12.0D0*(a2*a2-3.D0*a3*a3)*nn*nb*nb)
       PARAMETER( sf_d3 = 12.0D0*a5*ab3*(a2*ab1-3.0D0*a3*ab2)* nn
     &                    + 12.0D0*a5*ab3*(    a3*a3-a2*a2)* nb
     &                    + 12.0D0*   ab5*(3.0D0*a3*a3-a2*a2)*nn*nb     
     &                                                                 )
       PARAMETER( sf_d4 = 48.0D0*a5*a5*ab3*ab3*nn + 
     &                                         36.0D0*a5*ab3*ab5*nn*nn )
       PARAMETER( sf_d5 = 3.0D0*(a2*a2-3.0D0*a3*a3)
     &                       *(ab1*ab1-ab2*ab2)*nn    )
       PARAMETER( sf_n0  = 36.0D0*a1*nn*nn*nb*nb )
       PARAMETER( sf_n1  = - 12.0D0*a5*ab3*(ab1+ab2)*nn*nn
     &                    + 8.0D0*a5*ab3*(6.0D0*a1-a2-3.0D0*a3)*nn*nb
     &                    + 36.0D0*a1*ab5*nn*nn*Nb )
       PARAMETER( sf_n2  = 9.0D0*a1*(ab2*ab2-ab1*ab1)*nn*nn )
       PARAMETER( sf_nb0 = 12.0D0*ab3*nn*nn*nn*nb  )
       PARAMETER( sf_nb1 = 12.0D0*a5*ab3*ab3*nn*nn )
       PARAMETER( sf_nb2 = 9.0D0*a1*ab3*(ab1-ab2)*nn*nn + ( 
     &                                            6.0D0*a1*(a2-3.0D0*a3)
     &                               - 4.0D0*(a2*a2-3.0D0*a3*a3) )*ab3 
     &                                                        * nn * nb)
       PARAMETER( lim_am0 = sf_d0*sf_n0               )
       PARAMETER( lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0 )
       PARAMETER( lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0 )
       PARAMETER( lim_am3 = sf_d4*sf_n1               )
       PARAMETER( lim_am4 = sf_d2*sf_n0               )
       PARAMETER( lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0   )
       PARAMETER( lim_am6 = sf_d3*sf_n1               )
      real u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hz_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Huon(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hvom(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Wi(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real rho(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real svstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      real    sms_origin_date_in_sec
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_origin_date_in_sec
      common /smsdat3/ sms_cycle, sms_scale
      common /smsdat4/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat5/ lsvsgrd,sms_tid, susid, svsid
      real wspd_cfb(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /smsdat_wspd_cfb/ wspd_cfb
      integer*4 lwgrd, wid
      common /smsdat5/ lwgrd, wid
      real bustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real bvstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real bvstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real shflx_rsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rsw/shflx_rsw
      real shflx_rlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rlw/shflx_rlw
      real shflx_lat(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_lat/shflx_lat
      real shflx_sen(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_sen/shflx_sen
      real stflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      REAL(kind=8) :: stf_origin_date_in_sec
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ stf_origin_date_in_sec
      common /stfdat3/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat4/  stf_tid, stf_id
      real btflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_btflx/btflx
      real tair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rhum(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real prate(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real radlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real radsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real wspd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real uwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /bulk_tair/ tair
      common /bulk_rhum/ rhum
      common /bulk_prate/ prate
      common /bulk_radlw/ radlw
      common /bulk_radsw/ radsw
      common /bulk_wspd/ wspd
      common /bulk_uwnd/ uwnd
      common /bulk_vwnd/ vwnd
      real tairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real rhumg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real prateg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real radlwg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real radswg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real uwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real vwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real wspdg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulkdat_tairg/tairg
      common /bulkdat_rhumg/rhumg
      common /bulkdat_prateg/prateg
      common /bulkdat_radlwg/radlwg
      common /bulkdat_radswg/radswg
      common /bulk_uwndg/uwndg
      common /bulk_vwndg/vwndg
      common /bulkdat_wspdg/wspdg
      real    tairp(2),rhump(2),pratep(2),radlwp(2),radswp(2)
      real    uwndp(2),vwndp(2)
      real    bulk_time(2), bulk_cycle
      integer*4 tair_id,rhum_id,prate_id,radlw_id,radsw_id
      integer*4 ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      REAL(kind=8) :: blk_origin_date_in_sec
      integer*4 uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      integer*4 itbulk,bulk_ncycle,bulk_rec,bulk_tid
      integer*4 bulkunused
      common /bulkdat1_for/ tair_id,rhum_id,prate_id,radlw_id,radsw_id
      common /bulkdat1_grd
     &                 / ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      common /bulkdat1_tim/ itbulk, bulk_ncycle, bulk_rec, bulk_tid
      common /bulkdat1_uns/ bulkunused
      common /bulkdat1_wnd/ uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      common /bulkdat2_for/ tairp,rhump,pratep,radlwp,radswp
      common /bulkdat2_tim/ bulk_time, bulk_cycle, 
     &                                            blk_origin_date_in_sec
      common /bulkdat2_wnd/ uwndp,vwndp
      real srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_srflx/srflx
      real srflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      REAL(kind=8) :: srf_origin_date_in_sec
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ srf_origin_date_in_sec
      common /srfdat3/ itsrf,srf_ncycle,srf_rec,lsrfgrd,srf_tid,srf_id
      real h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_angler/angler
      real latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv
      real pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real zob(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /Z0B_VAR/zob
      real zeta_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real ubar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vbar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_zeta/zeta_avg
     &       /avg_ubar/ubar_avg
     &       /avg_vbar/vbar_avg
      real bostr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_bostr/bostr_avg
      real bustr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_bustr/bustr_avg
      real bvstr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_bvstr/bvstr_avg
      real wstr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_wstr/wstr_avg
      real sustr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_sustr/sustr_avg
      real svstr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_svstr/svstr_avg
      real srflx_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_srflx/srflx_avg
      real u_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real v_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real t_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real rho_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real bvf_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real omega_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real w_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &       /avg_rho/rho_avg /avg_omega/omega_avg
     &       /avg_bvf/bvf_avg
     &       /avg_w/w_avg
      real stflx_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /avg_stflx/stflx_avg
      real btflx_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /avg_btflx/btflx_avg
      real hbl_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_hbl/hbl_avg
      real tke_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real gls_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Lscale_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /avg_tke/tke_avg
      common /avg_gls/gls_avg
      common /avg_Lscale/Lscale_avg
      real shflx_rsw_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real shflx_rlw_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real shflx_lat_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real shflx_sen_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_shflx_rsw/shflx_rsw_avg
      common /avg_shflx_rlw/shflx_rlw_avg
      common /avg_shflx_lat/shflx_lat_avg
      common /avg_shflx_sen/shflx_sen_avg
      real diff3d_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /avg_diff3d/diff3d_avg
      real Akv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /avg_Akv/Akv_avg /avg_Akt/Akt_avg
      if (.not.WEST_INTER) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (.not.EAST_INTER) then
        imax=Iend
      else
        imax=Iend+1
      endif
      if (.not.SOUTH_INTER) then
        jmin=Jstr
      else
        jmin=Jstr-1
      endif
      if (.not.NORTH_INTER) then
        jmax=Jend
      else
        jmax=Jend+1
      endif
      cm0     =  ( (a2*a2 - 3.0D0*a3*a3 + 3.0D0*a1*nn)/(3.0D0*nn*nn) 
     &                                                         )**0.25D0
      cm0inv2 = 1.D0/cm0**2
      alpha_n_min = 0.5D0*( - ( sf_d1 + sf_nb0 )
     &             + sqrt(  ( sf_d1 + sf_nb0 )**2
     &            - 4.D0 * sf_d0 *( sf_d4 + sf_nb1 ) ) )
     &                            / ( sf_d4 + sf_nb1 )
      cff     = (cm0**3 )*(tke_min**1.5D0) / eps_min
      gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )
      trb_min(itke) = tke_min
      trb_min(igls) = gls_min
      ig=itke
      DO k=1,N-1
        DO j=jmin,jmax
          DO i=imin,imax+1
            HUon_w = 0.5D0*(Huon(i,j,k)+Huon(i,j,k+1)) * umask(i,j)
            ustar_sfc_sq(i,j)=
     &            trb(i-1,j,k,nstp,ig)*max(HUon_w,0.D0)
     &           +trb(i  ,j,k,nstp,ig)*min(HUon_w,0.D0)
          ENDDO
        ENDDO
        DO j=jmin,jmax+1
          DO i=imin,imax
            HVom_w = 0.5D0*(Hvom(i,j,k)+Hvom(i,j,k+1)) * vmask(i,j)
            ustar_bot_sq(i,j)=
     &            trb(i,j-1,k,nstp,ig)*max(HVom_w,0.D0)
     &           +trb(i,j  ,k,nstp,ig)*min(HVom_w,0.D0)
          ENDDO
        ENDDO
        DO j=jmin,jmax
          DO i=imin,imax
            cff    = 2.D0 / ( Hz(i,j,k)+ Hz(i,j,k+1) )
            tke_new(i,j,k)= trb(i,j,k,nstp,ig)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                       *(  ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                         + ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
          ENDDO
        ENDDO
      ENDDO
      ig=igls
      DO k=1,N-1
        DO j=jmin,jmax
          DO i=imin,imax+1
            HUon_w = 0.5D0*(Huon(i,j,k)+Huon(i,j,k+1)) * umask(i,j)
            ustar_sfc_sq(i,j)=
     &            trb(i-1,j,k,nstp,ig)*max(HUon_w,0.D0)
     &           +trb(i  ,j,k,nstp,ig)*min(HUon_w,0.D0)
          ENDDO
        ENDDO
        DO j=jmin,jmax+1
          DO i=imin,imax
            HVom_w = 0.5D0*(Hvom(i,j,k)+Hvom(i,j,k+1)) * vmask(i,j)
            ustar_bot_sq(i,j)=
     &            trb(i,j-1,k,nstp,ig)*max(HVom_w,0.D0)
     &           +trb(i,j  ,k,nstp,ig)*min(HVom_w,0.D0)
          ENDDO
        ENDDO
        DO j=jmin,jmax
          DO i=imin,imax
            cff    = 2.D0 / ( Hz(i,j,k)+ Hz(i,j,k+1) )
            gls_new(i,j,k)= trb(i,j,k,nstp,ig)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                       *(  ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                         + ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
          ENDDO
        ENDDO
      ENDDO
      ig=itke
      DO j=jmin,jmax
        DO k=1,N
          DO i=imin,imax
            cff=0.5D0*(We(i,j,k)+We(i,j,k-1))
            FC(i,k)=trb(i,j,k-1,nstp,ig)*max(cff,0.D0)
     &             +trb(i,j,k  ,nstp,ig)*min(cff,0.D0)
          ENDDO
        ENDDO
        DO k=1,N-1
          DO i=imin,imax
            cff=2.D0/(Hz(i,j,k)+Hz(i,j,k+1))
            tke_new(i,j,k)= tke_new(i,j,k)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                        *(FC(i,k+1)-FC(i,k))
          ENDDO
        ENDDO
      ENDDO
      ig=igls
      DO j=jmin,jmax
        DO k=1,N
          DO i=imin,imax
            cff=0.5D0*(We(i,j,k)+We(i,j,k-1))
            FC(i,k)=trb(i,j,k-1,nstp,ig)*max(cff,0.D0)
     &             +trb(i,j,k  ,nstp,ig)*min(cff,0.D0)
          ENDDO
        ENDDO
        DO k=1,N-1
          DO i=imin,imax
            cff=2.D0/(Hz(i,j,k)+Hz(i,j,k+1))
            gls_new(i,j,k)= gls_new(i,j,k)
     &                     -dt*cff*pm(i,j)*pn(i,j)
     &                        *(FC(i,k+1)-FC(i,k))
          ENDDO
        ENDDO
      ENDDO
      tind  = nrhs
      DO k=1,N-1
         DO j=jmin,jmax
            DO i=imin,imax
               cff = 1.D0 / ( Hz( i, j, k ) + Hz( i, j, k+1 ) )
               du  = cff*( u(i, j, k+1,tind)+u(i+1, j, k+1,tind)
     &                    -u(i, j, k  ,tind)-u(i+1, j, k  ,tind))
               dv  = cff*( v(i, j, k+1,tind)+v(i, j+1, k+1,tind)
     &                    -v(i, j, k  ,tind)-v(i, j+1, k  ,tind))
               shear2(i,j,k) = du*du + dv*dv
           ENDDO
         ENDDO
      ENDDO
      DO j=jmin,jmax
         DO i=imin,imax
            su_r=0.5D0*(sustr(i,j)+sustr(i+1,j))
            sv_r=0.5D0*(svstr(i,j)+svstr(i,j+1))
            ustar_sfc_sq( i, j ) = sqrt(su_r**2+sv_r**2)
            ustar_bot_sq( i, j ) =
     &                        sqrt( (0.5D0*(bustr(i,j)+bustr(i+1,j)))**2
     &                             
     &                           +(0.5D0*(bvstr(i,j)+bvstr(i,j+1)))**2 )
         ENDDO
      ENDDO
      DO j=jmin,jmax
         DO i=imin,imax
            DO k=1,N-1
               cff       = (cm0**e1) * ( trb( i,j,k,nstp,itke )**e2 )
     &                               * ( trb( i,j,k,nstp,igls )**e3 )
               diss(i,k) = MAX( cff , eps_min )
            ENDDO
         ENDDO
         DO ig = 1,ngls
            cff=-0.5D0*dt
            DO k=2,N-1
               DO i=imin,imax
                  FC(i,k)= cff* OneOverSig(ig)*
     &                   ( Akv_old(i,j,k)+Akv_old(i,j,k-1) ) / Hz(i,j,k)
               ENDDO
            ENDDO
            DO i=imin,imax
               FC(i,1)=0.D0
               FC(i,N)=0.D0
            ENDDO
            DO k=1,N-1
               DO i=imin,imax
                  ig1   = (igls-ig)
                  ig2   = (ig-itke)
                  invk  =     1.D0 / trb( i,j,k,nstp,itke )
                  gls   =          trb( i,j,k,nstp,igls )
                  invG  =  ig1+ig2*(1.D0/gls)
                  cff1  =  ig1+ig2*beta1   * invk*gls
                  cff2  = (ig1+ig2*beta2 ) * invk
                  cff3m =  ig1+ig2*beta3m  * invk*gls
                  cff3p =  ig1+ig2*beta3p  * invk*gls
                  Sprod =  cff1*Akv_old(i,j,k) * shear2(i,j,k)
                  Bprod = -Akt_old(i,j,k)*( cff3m*MAX(bvf(i,j,k),0.D0)
     &                                   +  cff3p*MIN(bvf(i,j,k),0.D0) )
                  cff   =       0.5D0*(Hz(i,j,k)+Hz(i,j,k+1))
                  IF( ig == itke ) THEN
                    trb_new=tke_new(i,j,k)
                  ELSE
                    trb_new=gls_new(i,j,k)
                  ENDIF
                  IF( (Bprod + Sprod) .gt. 0.D0) THEN
                     RH(i,k) = cff*( trb_new + dt*(Bprod+Sprod) )
                     DC(i,k) = cff*(1.D0+dt*cff2*diss(i,k))
     &                                                -FC(i,k)-FC(i,k+1)
                  ELSE
                     RH(i,k) = cff*( trb_new + dt*Sprod  )
                     DC(i,k) = cff*(1.D0+dt*(cff2*diss(i,k)
     &                              -invG*Bprod)) - FC(i,k) - FC(i,k+1)
                  ENDIF
               ENDDO
            ENDDO
            IF( ig == itke ) THEN
               DO i=imin,imax
                  trb_sfc = MAX( tke_min, cm0inv2*ustar_sfc_sq(i,j) )
                  flux_top = 0.D0
                  trb_bot = MAX( tke_min, cm0inv2*ustar_bot_sq(i,j) )
                  flux_bot = 0.D0
                  RH(i,1  ) = RH(i,  1) + dt*flux_bot
                  RH(i,N-1) = RH(i,N-1) + dt*flux_top
                  tke_new(i,j,N) = trb_sfc
                  tke_new(i,j,0) = trb_bot
               ENDDO
            ELSE
               DO i=imin,imax
                  z0_s=MAX( Zosmin , chk*ustar_sfc_sq(i,j) )
                  cff = 0.5D0*(tke_new(i,j,N-1)+tke_new(i,j,N))
                  lgthsc = vonKar*(0.5D0*Hz(i,j,N)+z0_s)
                  trb_sfc = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                              *(cff**rm))
                  flux_top = -rn*cm0**(rp+1.D0)
     &                             *vonKar*OneOverSig(igls)
     &                             *(cff**(rm+0.5D0))*(lgthsc**rn)
                  z0_b = MAX( Zob(i,j) , Zobmin )
                  cff = 0.5D0*(tke_new(i,j,1)+tke_new(i,j,0))
                  lgthsc = vonKar*(0.5D0*Hz(i,j,1)+z0_b)
                  trb_bot = MAX(gls_min,(cm0**rp)*(lgthsc**rn)
     &                                   *(tke_new(i,j,0)**rm))
                  flux_bot =-rn*cm0**(rp+1.D0)
     &                          *vonKar*OneOverSig(igls)
     &                          *(cff**(rm+0.5D0))*(lgthsc**rn)
                  RH(i,  1) = RH(i,  1) + dt*flux_bot
                  RH(i,N-1) = RH(i,N-1) + dt*flux_top
                  gls_new(i,j,N) = trb_sfc
                  gls_new(i,j,0) = trb_bot
               ENDDO
            ENDIF
            DO i=imin,imax
               cff       =  1.D0/DC(i,N-1)
               CF(i,N-1) = cff*FC(i,N-1)
               RH(i,N-1) = cff*RH(i,N-1)
            ENDDO
            DO k=N-2,1,-1
               DO i=imin,imax
                  cff     =   1.D0/(DC(i,k)-CF(i,k+1)*FC(i,k+1))
                  CF(i,k) = cff*FC(i,k)
                  RH(i,k) = cff*( RH(i,k)-FC(i,k+1)*RH(i,k+1))
               ENDDO
            ENDDO
            IF( ig == itke ) THEN
              DO i=imin,imax
                tke_new(i,j,1) = MAX( RH(i,1), trb_min(ig) )
              ENDDO
              DO k=2,N-1
                DO i=imin,imax
                  RH(i,k) = RH(i,k)-CF(i,k)*RH(i,k-1)
                  tke_new(i,j,k) = MAX( RH(i,k), trb_min(ig) )
                ENDDO
              ENDDO
            ELSE
              DO i=imin,imax
                gls_new(i,j,1) = MAX( RH(i,1), trb_min(ig) )
              ENDDO
              DO k=2,N-1
                DO i=imin,imax
                  RH(i,k) = RH(i,k)-CF(i,k)*RH(i,k-1)
                  gls_new( i,j,k) = MAX( RH(i,k), trb_min(ig) )
                ENDDO
              ENDDO
            ENDIF
         ENDDO
      ENDDO
      ig = igls
      DO k=0,N
      if (.not.WEST_INTER) then
        do j=jmin,jmax
          gls_new(Istr-1,j,k)=gls_new(Istr,j,k)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=jmin,jmax
          gls_new(Iend+1,j,k)=gls_new(Iend,j,k)
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=imin,imax
          gls_new(i,Jstr-1,k)=gls_new(i,Jstr,k)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=imin,imax
          gls_new(i,Jend+1,k)=gls_new(i,Jend,k)
        enddo
      endif
      if (.not.WEST_INTER.and..not.SOUTH_INTER) then
        gls_new(Istr-1,Jstr-1,k)=gls_new(Istr,Jstr,k)
      endif
      if (.not.WEST_INTER.and..not.NORTH_INTER) then
        gls_new(Istr-1,Jend+1,k)=gls_new(Istr,Jend,k)
      endif
      if (.not.EAST_INTER.and..not.SOUTH_INTER) then
        gls_new(Iend+1,Jstr-1,k)=gls_new(Iend,Jstr,k)
      endif
      if (.not.EAST_INTER.and..not.NORTH_INTER) then
        gls_new(Iend+1,Jend+1,k)=gls_new(Iend,Jend,k)
      endif
         DO j=jstr-1,jend+1
            DO i=istr,iend+1
               ustar_sfc_sq (i,j  )=( gls_new(i  ,j,k)
     &                   -  gls_new(i-1,j,k) )
     &                             *umask(i,j)
            ENDDO
         ENDDO
         DO j=jstr,jend+1
            DO i=istr-1,iend+1
               work(i,j,0)=( gls_new(i,j  ,k)
     &                    - gls_new(i,j-1,k) )
     &                             *vmask(i,j)
            ENDDO
            DO i=istr,iend
              ustar_bot_sq(i,j)=work(i,j,0)
     &                + smth_a*( ustar_sfc_sq(i+1,j)+ustar_sfc_sq(i  
     &                                                             ,j-1)
     &                          -ustar_sfc_sq(i  
     &                                        ,j)-ustar_sfc_sq(i+1,j-1))
            ENDDO
         ENDDO
         DO j=jstr,jend
            DO i=istr,iend+1
              ustar_sfc_sq(i,j)=ustar_sfc_sq(i,j  )
     &                + smth_a*( work(i,j+1,0)+work(i-1,j  ,0)
     &                          -work(i,j  ,0)-work(i-1,j+1,0))
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=gls_new(i,j,k)
     &                        + smth_b*( 
     &                             ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                                  
     &                          +ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
     &                                               *rmask(i,j)
            ENDDO
         ENDDO
      ENDDO
        ig = itke
      DO k=0,N
      if (.not.WEST_INTER) then
        do j=jmin,jmax
          tke_new(Istr-1,j,k)=tke_new(Istr,j,k)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=jmin,jmax
          tke_new(Iend+1,j,k)=tke_new(Iend,j,k)
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=imin,imax
          tke_new(i,Jstr-1,k)=tke_new(i,Jstr,k)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=imin,imax
          tke_new(i,Jend+1,k)=tke_new(i,Jend,k)
        enddo
      endif
      if (.not.WEST_INTER.and..not.SOUTH_INTER) then
        tke_new(Istr-1,Jstr-1,k)=tke_new(Istr,Jstr,k)
      endif
      if (.not.WEST_INTER.and..not.NORTH_INTER) then
        tke_new(Istr-1,Jend+1,k)=tke_new(Istr,Jend,k)
      endif
      if (.not.EAST_INTER.and..not.SOUTH_INTER) then
        tke_new(Iend+1,Jstr-1,k)=tke_new(Iend,Jstr,k)
      endif
      if (.not.EAST_INTER.and..not.NORTH_INTER) then
        tke_new(Iend+1,Jend+1,k)=tke_new(Iend,Jend,k)
      endif
         DO j=jstr-1,jend+1
            DO i=istr,iend+1
               ustar_sfc_sq (i,j  )=( tke_new(i  ,j,k)
     &                   -  tke_new(i-1,j,k) )
     &                             *umask(i,j)
            ENDDO
         ENDDO
         DO j=jstr,jend+1
            DO i=istr-1,iend+1
               work(i,j,0)=( tke_new(i,j  ,k)
     &                    - tke_new(i,j-1,k) )
     &                             *vmask(i,j)
            ENDDO
            DO i=istr,iend
              ustar_bot_sq(i,j)=work(i,j,0)
     &                + smth_a*( ustar_sfc_sq(i+1,j)+ustar_sfc_sq(i  
     &                                                             ,j-1)
     &                          -ustar_sfc_sq(i  
     &                                        ,j)-ustar_sfc_sq(i+1,j-1))
            ENDDO
         ENDDO
         DO j=jstr,jend
            DO i=istr,iend+1
              ustar_sfc_sq(i,j)=ustar_sfc_sq(i,j  )
     &                + smth_a*( work(i,j+1,0)+work(i-1,j  ,0)
     &                          -work(i,j  ,0)-work(i-1,j+1,0))
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=tke_new(i,j,k)
     &                        + smth_b*( 
     &                             ustar_sfc_sq(i+1,j)-ustar_sfc_sq(i,j)
     &                                  
     &                          +ustar_bot_sq(i,j+1)-ustar_bot_sq(i,j) )
     &                                               *rmask(i,j)
            ENDDO
         ENDDO
      ENDDO
      DO j=jstr,jend
         DO k=1,N-1
            DO i=istr,iend
               L_lim = galp * sqrt( 2.D0* trb(i,j,k,nnew,itke)) /
     &                            ( sqrt(max(eps, bvf(i,j,k)))  )
               cff = (cm0**rp)*(L_lim**rn)*(trb(i,j,k,nnew,itke)**rm)
               trb( i,j,k,nnew,igls ) = MAX( trb( i,j,k,nnew,igls ),cff)
               epsilon = (cm0**e1) * ( trb( i,j,k,nnew,itke )**e2 )
     &                             * ( trb( i,j,k,nnew,igls )**e3 )
               epsilon = MAX(epsilon, eps_min)
               cff     = ( trb(i,j,k,nnew,itke)/epsilon )**2
               alpha_m     = cff*  shear2(i,j,k)
               alpha_n     = cff*     bvf(i,j,k)
               alpha_n     = MIN(  MAX( 0.73D0*alpha_n_min , alpha_n ) 
     &                                                        , 1.0D10 )
               alpha_m_max = ( lim_am0 + lim_am1 * alpha_n
     &                                 + lim_am2 * alpha_n**2
     &                                 + lim_am3 * alpha_n**3) /
     &                         ( lim_am4 + lim_am5 * alpha_n
     &                                   + lim_am6 * alpha_n**2 )
               alpha_m = MIN(alpha_m , alpha_m_max)
               Denom = sf_d0  + sf_d1*alpha_n +  sf_d2*alpha_m
     &                        + sf_d3*alpha_n*alpha_m
     &                        + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
               c_mu      = (sf_n0  + sf_n1*alpha_n  + sf_n2*alpha_m)
     &                                                        /Denom
               c_mu_prim = (sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)
     &                                                        /Denom
               cff = trb( i,j,k,nnew,itke )**2 / epsilon
               Akv(i,j,k) = MAX( cff*c_mu , nuwm ) * rmask(i,j)
               Akt(i,j,k,itemp)= MAX( cff*c_mu_prim , nuws )
     &                                             * rmask(i,j)
               Akt(i,j,k,isalt)= MAX( cff*c_mu_prim , nuws )
     &                                             * rmask(i,j)
               Lscale( i, j , k ) =  cm0 * cm0 * cm0 * cff
     &                              / sqrt( trb( i,j,k,nnew,itke ) )
     &                                             * rmask(i,j)
               Eps_gls(i,j,k) = epsilon * rmask(i,j)
            ENDDO
         ENDDO
         DO i=istr,iend
           Akv(i,j,N) = MAX( 1.5D0*Akv(i,j,N-1)
     &                      -0.5D0*Akv(i,j,N-2), nuwm)
           Akv(i,j,0) = MAX( 1.5D0*Akv(i,j,  1)
     &                      -0.5D0*Akv(i,j,  2), nuwm)
           Akt(i,j,N,itemp) = MAX(  1.5D0*Akt(i,j,N-1,itemp)
     &                             -0.5D0*Akt(i,j,N-2,itemp), nuws )
           Akt(i,j,0,itemp) = MAX(  1.5D0*Akt(i,j,  1,itemp)
     &                             -0.5D0*Akt(i,j,  2,itemp), nuws )
           Akt(i,j,N,isalt)=Akt(i,j,N,itemp)
           Akt(i,j,0,isalt)=Akt(i,j,0,itemp)
         ENDDO
      ENDDO
      if (.not.WEST_INTER) then
        do j=jstr,jend
          do k=0,N
            trb(istr-1,j,k,nnew,itke)=trb(istr,j,k,nnew,itke)
            trb(istr-1,j,k,nnew,igls)=trb(istr,j,k,nnew,igls)
            Akv(istr-1,j,k      )=Akv(istr,j,k      )
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
          enddo
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=jstr,jend
          do k=0,N
            trb(iend+1,j,k,nnew,itke)=trb(iend,j,k,nnew,itke)
            trb(iend+1,j,k,nnew,igls)=trb(iend,j,k,nnew,igls)
            Akv(iend+1,j,k      )=Akv(iend,j,k      )
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
          enddo
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=istr,iend
          do k=0,N
            trb(i,jstr-1,k,nnew,itke)=trb(i,jstr,k,nnew,itke)
            trb(i,jstr-1,k,nnew,igls)=trb(i,jstr,k,nnew,igls)
            Akv(i,jstr-1,k      )=Akv(i,jstr,k      )
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
          enddo
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=istr,iend
          do k=0,N
            trb(i,jend+1,k,nnew,itke)=trb(i,jend,k,nnew,itke)
            trb(i,jend+1,k,nnew,igls)=trb(i,jend,k,nnew,igls)
            Akv(i,jend+1,k)=Akv(i,jend,k)
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
          enddo
        enddo
      endif
      if (.not.WEST_INTER .and. .not.SOUTH_INTER) then
        do k=0,N
          trb(istr-1,jstr-1,k,nnew,itke)=trb(istr,jstr,k,nnew,itke)
          trb(istr-1,jstr-1,k,nnew,igls)=trb(istr,jstr,k,nnew,igls)
          Akv(istr-1,jstr-1,k      )=Akv(istr,jstr,k      )
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
        enddo
      endif
      if (.not.WEST_INTER .and. .not.NORTH_INTER) then
        do k=0,N
          trb(istr-1,jend+1,k,nnew,itke)=trb(istr,jend,k,nnew,itke)
          trb(istr-1,jend+1,k,nnew,igls)=trb(istr,jend,k,nnew,igls)
          Akv(istr-1,jend+1,k      )=Akv(istr,jend,k      )
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
        enddo
      endif
      if (.not.EAST_INTER .and. .not.SOUTH_INTER) then
        do k=0,N
          trb(iend+1,jstr-1,k,nnew,itke)=trb(iend,jstr,k,nnew,itke)
          trb(iend+1,jstr-1,k,nnew,igls)=trb(iend,jstr,k,nnew,igls)
          Akv(iend+1,jstr-1,k      )=Akv(iend,jstr,k      )
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
        enddo
      endif
      if (.not.EAST_INTER .and. .not.NORTH_INTER) then
        do k=0,N
          trb(iend+1,jend+1,k,nnew,itke)=trb(iend,jend,k,nnew,itke)
          trb(iend+1,jend+1,k,nnew,igls)=trb(iend,jend,k,nnew,igls)
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
        enddo
      endif
      do j=Jstr,Jend
        do i=Istr,Iend
          kbl(i,j)=1
          hbl(i,j)=z_w(i,j,N)-z_r(i,j,1)
          hbl(i,j)=hbl(i,j)*rmask(i,j)
        enddo
      enddo
      kref=max(1,N-3)
      do k=1,kref
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=rho1(i,j,k)-rho1(i,j,kref)
            if (cff.gt.0.01D0) then
              kbl(i,j)=k
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,k)
              hbl(i,j)=hbl(i,j)*rmask(i,j)
            endif
          enddo
        enddo
      enddo
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(-1,-1,0,nnew,itke))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        trb(-1,-1,0,nnew,igls))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(-1,-1,0,itemp))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(-1,-1,0,isalt))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Lscale(-1,-1,0))
      return
      end
