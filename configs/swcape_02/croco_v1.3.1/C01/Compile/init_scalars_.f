      subroutine init_scalars (ierr)
      implicit none
      integer*4 ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3
      integer*4 omp_get_num_threads
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4 indxS
      parameter (indxS=indxV+ntrc_temp+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxsurft,indxsurfs,indxsurfz,indxsurfu,
     &        indxsurfv
      parameter (indxsurft=indxV+ntrc_temp+ntrc_salt
     &                           +ntrc_pas+ntrc_bio+ntrc_sed
     &                           +ntrc_diats+ntrc_diauv+ntrc_diavrt
     &                           +ntrc_diaek+ntrc_diapv+ntrc_diaeddy+40
     &                                                                0,
     &           indxsurfs=indxsurft+1,
     &           indxsurfz=indxsurfs+1,
     &           indxsurfu=indxsurfz+1,
     &           indxsurfv=indxsurfu+1)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkv+ntrc_temp+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxTke
      parameter (indxTke=indxAkv+ntrc_temp+7)
      integer*4 indxGls
      parameter (indxGls=indxAkv+ntrc_temp+8)
      integer*4 indxLsc
      parameter (indxLsc=indxAkv+ntrc_temp+9)
      integer*4 indxAkk
      parameter (indxAkk=indxAkv+ntrc_temp+10)
      integer*4 indxAkp
      parameter (indxAkp=indxAkv+ntrc_temp+11)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWSPD,indxTAIR,indxRHUM,indxRADLW,indxRADSW,
     &        indxPRATE,indxUWND,indxVWND,indxPATM
      parameter (indxWSPD=indxSST+3,  indxTAIR=indxSST+4,
     &           indxRHUM=indxSST+5,  indxRADLW=indxSST+6,
     &           indxRADSW=indxSST+7, indxPRATE=indxSST+8,
     &           indxUWND=indxSST+9,  indxVWND=indxSST+10,
     &           indxPATM=indxSST+11)
      integer*4 indxShflx_rlw,indxShflx_lat,indxShflx_sen
      parameter (indxShflx_rlw=indxSST+12,
     &           indxShflx_lat=indxSST+13, indxShflx_sen=indxSST+14)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstAkv,rstAkt
      integer*4 rstAks
      integer*4 rstTke,rstGls
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
     &      , hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
      integer*4 hisT(NT)
      integer*4 ncidsurf, nrecsurf, nrpfsurf
     &      , surfTime, surfTime2, surfTstep
     &      , surf_surft(2), surf_surfs(2),  surf_surfz(2)
     &      , surf_surfu(2), surf_surfv(2)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
      integer*4 avgT(NT)
      integer*4 avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
       integer*4 ncidsurf_avg, nrecsurf_avg, nrpfsurf_avg
     &      , surfTime_avg, surfTime2_avg, surfTstep_avg
     &      , surf_surft_avg(2), surf_surfs_avg(2), surf_surfz_avg(2)
     &      , surf_surfu_avg(2), surf_surfv_avg(2)
      logical wrthis(800+NT)
     &      , wrtavg(800+NT)
     &      , wrtsurf(3)
     &      , wrtsurf_avg(3)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstAkv,rstAkt
     &      , rstAks
     &      , rstTke,rstGls
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
     &      , hisShflx_rlw
     &      , hisShflx_lat, hisShflx_sen
     &      , ncidsurf, nrecsurf, nrpfsurf
     &      , surfTime, surfTime2, surfTstep
     &      , surf_surft, surf_surfs,  surf_surfz
     &      , surf_surfu, surf_surfv
     &      , ncidsurf_avg, nrecsurf_avg, nrpfsurf_avg
     &      , surfTime_avg, surfTime2_avg, surfTstep_avg
     &      , surf_surft_avg, surf_surfs_avg,  surf_surfz_avg
     &      , surf_surfu_avg, surf_surfv_avg
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
     &      , avgShflx_rlw
     &      , avgShflx_lat, avgShflx_sen
     &      , wrthis
     &      , wrtavg
     &      , wrtsurf
     &      , wrtsurf_avg
      character*80 date_str, title
      character*80 origin_date, start_date_run, xios_origin_date
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8) :: origin_date_in_sec, xios_origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,  surfname
     &                                ,  surfname_avg
     &                                ,   clmname
      character*75  vname(20, 800)
      common /cncscrum/   date_str,   title
     &         ,   origin_date, start_date_run
     &         ,   xios_origin_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   xios_origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,  surfname
     &                                ,  surfname_avg
     &                                ,   clmname
     &                                ,   vname
      real bry_time(2)
      REAL(kind=8) :: bry_origin_date_in_sec
      common /bry_indices_array/ bry_time,
     &        bry_origin_date_in_sec
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
      real ssh(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /climat_ssh/ssh
      real Znudgcof(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /climat_Znudgcof/Znudgcof
      real sshg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /climat_sshg/sshg
      real    ssh_time(2)
      real    ssh_cycle
      integer*4 itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id
      REAL(kind=8) :: ssh_origin_date_in_sec
      common /climat_zdat1/ ssh_time, ssh_origin_date_in_sec
      common /climat_zdat2/ ssh_cycle
      common /climat_zdat3/
     &        itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id
      real tclm(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      common /climat_tclm/tclm
      real Tnudgcof(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      common /climat_Tnudgcof/Tnudgcof
      real tclima(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,2,NT)
      common /climat_tclima/tclima
      real tclm_time(2,NT)
      real tclm_cycle(NT)
      integer*4 ittclm(NT), tclm_ncycle(NT), tclm_rec(NT),
     &        tclm_tid(NT), tclm_id(NT)
      logical got_tclm(NT)
      REAL(kind=8) :: tclm_origin_date_in_sec
      common /climat_tdat/  tclm_time,       tclm_cycle,
     &        ittclm,       tclm_ncycle,     tclm_rec,
     &                      tclm_tid,        tclm_id,
     &                                       got_tclm,
     &                        tclm_origin_date_in_sec
      real ubclm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vbclm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /climat_ubclm/ubclm /climat_vbclm/vbclm
      real uclm(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real vclm(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /climat_uclm/uclm /climat_vclm/vclm
      real M2nudgcof(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /climat_M2nudgcof/M2nudgcof
      real ubclima(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real vbclima(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /climat_ubclima/ubclima /climat_vbclima/vbclima
      real M3nudgcof(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /climat_M3nudgcof/M3nudgcof
      real uclima(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,2)
      real vclima(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,2)
      common /climat_uclima/uclima /climat_vclima/vclima
      real     uclm_time(2)
      real     uclm_cycle
      integer*4 ituclm, uclm_ncycle, uclm_rec, uclm_tid,
     &        ubclm_id, vbclm_id, uclm_id, vclm_id
      REAL(kind=8) :: uclm_origin_date_in_sec
      common /climat_udat1/  uclm_time, uclm_origin_date_in_sec
      common /climat_udat2/  uclm_cycle
      common /climat_udat3/
     &             ituclm,   uclm_ncycle, uclm_rec,
     &             uclm_tid, ubclm_id,    vbclm_id,
     &             uclm_id,  vclm_id
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
      if (mynode.eq.0) write(stdout,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
      if (numthreads.gt.NPP) then
        if (mynode.eq.0) write(stdout,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds NPP =', NPP
        ierr=ierr+1
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
        if (mynode.eq.0) write(stdout,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1
      endif
      time=0.D0
      tdays=0.D0
      PREDICTOR_2D_STEP=.FALSE.
      iic=0
      kstp=1
      krhs=1
      knew=1
      ntstart=1
      nstp=1
      nrhs=1
      nnew=1
      nfast=1
      synchro_flag=.true.
      first_time=0
      may_day_flag=0
      do j=0,NPP
        do i=0,31
          proc(i,j)=0
          CPU_time(i,j)=0.D0
        enddo
      enddo
      trd_count=0
      nrecrst=0
      nrechis=0
      nrecavg=0
      tile_count=0
      bc_count=0
      avgke=0.D0
      avgpe=0.D0
      avgkp=0.D0
      volume=0.D0
      hmin=+1.D+20
      hmax=-1.D+20
      grdmin=+1.D+20
      grdmax=-1.D+20
      Cu_min=+1.D+20
      Cu_max=-1.D+20
      lonmin=+1.D+20
      lonmax=-1.D+20
      latmin=+1.D+20
      latmax=-1.D+20
      bc_crss=0.D0
      rx0=-1.D+20
      rx1=-1.D+20
      ncidrst=-1
      ncidhis=-1
      ncidavg=-1
      ncidfrc=-1
      ncidbulk=-1
      ncidclm=-1
      ncidqbar=-1
      ncidbtf=-1
      ncidsurf=-1
      nrecsurf=0
      ncidsurf_avg=-1
      nrecsurf_avg=0
       bry_id=-1
      call get_date (date_str)
      csmooth=dt/(86400.D0)
      ismooth=int(1.D0/csmooth)
      vname(:,:)='  '
      vname(1,indxTime)='scrum_time'
      vname(2,indxTime)='time since initialization'
      vname(3,indxTime)='second'
      vname(4,indxTime)='time, scalar, series'
      vname(5,indxTime)='time'
      vname(6,indxTime)='    '
      vname(7,indxTime)='T'
      vname(8,indxTime)='.NO.'
      vname(9,indxTime)=' '
      vname(1,indxTime2)='time'
      vname(2,indxTime2)='time since initialization'
      vname(3,indxTime2)='second'
      vname(4,indxTime2)='time, scalar, series'
      vname(5,indxTime2)='time'
      vname(6,indxTime2)='  '
      vname(7,indxTime2)='T '
      vname(8,indxTime2)='.NO.'
      vname(9,indxTime2)=' '
      vname(1,indxZ)='zeta                                        '
      vname(2,indxZ)='free-surface                                '
      vname(3,indxZ)='meter                                       '
      vname(4,indxZ)='free-surface, scalar, series                '
      vname(5,indxZ)='sea_surface_height                          '
      vname(6,indxZ)='lat_rho lon_rho                             '
      vname(7,indxZ)='                                            '
      vname(1,indxUb)='ubar                                       '
      vname(2,indxUb)='vertically integrated u-momentum component '
      vname(3,indxUb)='meter second-1                             '
      vname(4,indxUb)='ubar-velocity, scalar, series              '
      vname(5,indxUb)='barotropic_sea_water_x_'
     &                                  // 'velocity_at_u_location'
      vname(6,indxUb)='lat_u lon_u                                '
      vname(7,indxUb)='                                           '
      vname(1,indxVb)='vbar                                       '
      vname(2,indxVb)='vertically integrated v-momentum component '
      vname(3,indxVb)='meter second-1                             '
      vname(4,indxVb)='vbar-velocity, scalar, series              '
      vname(5,indxVb)='barotropic_sea_water_y_'
     &                                  // 'velocity_at_v_location'
      vname(6,indxVb)='lat_v lon_v                                '
      vname(7,indxVb)='                                           '
      vname(11,indxVb)=' '
      vname(1,indxBostr)='bostr                                   '
      vname(2,indxBostr)='Kinematic bottom stress                 '
      vname(3,indxBostr)='N/m2                                    '
      vname(4,indxBostr)='                                        '
      vname(5,indxBostr)='                                        '
      vname(6,indxBostr)='lat_rho lon_rho                         '
      vname(7,indxBostr)='                                        '
      vname(1,indxBustr)='bustr                                   '
      vname(2,indxBustr)='Kinematic u bottom stress component     '
      vname(3,indxBustr)='N/m2                                    '
      vname(4,indxBustr)='                                        '
      vname(5,indxBustr)='                                        '
      vname(6,indxBustr)='lat_u lon_u                             '
      vname(7,indxBustr)='                                        '
      vname(1,indxBvstr)='bvstr                                   '
      vname(2,indxBvstr)='Kinematic v bottom stress component     '
      vname(3,indxBvstr)='N/m2                                    '
      vname(4,indxBvstr)='                                        '
      vname(5,indxBvstr)='                                        '
      vname(6,indxBvstr)='lat_v lon_v                             '
      vname(7,indxBvstr)='                                        '
      vname(1,indxWstr)='wstr                                     '
      vname(2,indxWstr)='Kinematic wind stress                    '
      vname(3,indxWstr)='N/m2                                     '
      vname(4,indxWstr)='                                         '
      vname(5,indxWstr)='magnitude_of_surface_downward_stress     '
      vname(6,indxWstr)='lat_rho lon_rho                          '
      vname(7,indxWstr)='                                         '
      vname(1,indxUWstr)='sustr                                   '
      vname(2,indxUWstr)='Kinematic u wind stress component       '
      vname(3,indxUWstr)='N/m2                                    '
      vname(4,indxUWstr)='                                        '
      vname(5,indxUWstr)='surface_downward_eastward_stress        '
      vname(6,indxUWstr)='lat_u lon_u                             '
      vname(7,indxUWstr)='                                        '
      vname(8,indxUWstr)='                                        '
      vname(1,indxVWstr)='svstr                                   '
      vname(2,indxVWstr)='Kinematic v wind stress component       '
      vname(3,indxVWstr)='N/m2                                    '
      vname(4,indxVWstr)='                                        '
      vname(5,indxVWstr)='surface_downward_northward_stress       '
      vname(6,indxVWstr)='lat_v lon_v                             '
      vname(7,indxVWstr)='                                        '
      vname(1,indxU)='u                                           '
      vname(2,indxU)='u-momentum component                        '
      vname(3,indxU)='meter second-1                              '
      vname(4,indxU)='u-velocity, scalar, series                  '
      vname(5,indxU)='sea_water_x_velocity_at_u_location          '
      vname(6,indxU)='lat_u lon_u                                 '
      vname(7,indxU)='                                            '
      vname(1,indxV)='v                                           '
      vname(2,indxV)='v-momentum component                        '
      vname(3,indxV)='meter second-1                              '
      vname(4,indxV)='v-velocity, scalar, series                  '
      vname(5,indxV)='sea_water_y_velocity_at_v_location          '
      vname(6,indxV)='lat_v lon_v                                 '
      vname(7,indxV)='                                            '
      vname(1,indxT)='temp                                        '
      vname(2,indxT)='potential temperature                       '
      vname(3,indxT)='Celsius                                     '
      vname(4,indxT)='temperature, scalar, series                 '
      vname(5,indxT)='sea_water_potential_temperature             '
      vname(6,indxT)='lat_rho lon_rho                             '
      vname(7,indxT)='                                            '
      vname(1,indxS)='salt                                        '
      vname(2,indxS)='salinity                                    '
      vname(3,indxS)='PSU                                         '
      vname(4,indxS)='salinity, scalar, series                    '
      vname(5,indxS)='sea_water_salinity                          '
      vname(6,indxS)='lat_rho lon_rho                             '
      vname(7,indxS)='                                            '
      vname(1,indxShflx)='shflux                                  '
      vname(2,indxShflx)='surface net heat flux                   '
      vname(3,indxShflx)='Watts meter-2                           '
      vname(4,indxShflx)='surface heat flux, scalar, series       '
      vname(6,indxShflx)='lat_rho lon_rho                         '
      vname(7,indxShflx)='                                        '
      vname(1,indxSwflx)='swflux                                  '
      vname(2,indxSwflx)='surface freshwater flux (E-P)           '
      vname(3,indxSwflx)='centimeter day-1                        '
      vname(4,indxSwflx)='surface freshwater flux, scalar, series '
      vname(5,indxSwflx)='                                        '
      vname(6,indxSwflx)='lat_rho lon_rho                         '
      vname(7,indxSwflx)='                                        '
      vname(8,indxSwflx)='                                        '
      vname(1,indxShflx_rsw)='radsw                               '
      vname(2,indxShflx_rsw)='Short-wave surface radiation        '
      vname(3,indxShflx_rsw)='Watts meter-2                       '
      vname(4,indxShflx_rsw)='                                    '
      vname(5,indxShflx_rsw)='                                    '
      vname(6,indxShflx_rsw)='lat_rho lon_rho                     '
      vname(7,indxShflx_rsw)='                                    '
      vname(1,indxShflx_rlw)='shflx_rlw                           '
      vname(2,indxShflx_rlw)='Long-wave surface radiation         '
      vname(3,indxShflx_rlw)='Watts meter-2                       '
      vname(5,indxShflx_rlw)='                                    '
      vname(6,indxShflx_rlw)='lat_rho lon_rho                     '
      vname(7,indxShflx_rlw)='                                    '
      vname(1,indxShflx_lat)='shflx_lat                           '
      vname(2,indxShflx_lat)='Latent surface heat flux            '
      vname(3,indxShflx_lat)='Watts meter-2                       '
      vname(4,indxShflx_sen)='                                    '
      vname(5,indxShflx_lat)='                                    '
      vname(6,indxShflx_lat)='lat_rho lon_rho                     '
      vname(7,indxShflx_lat)='                                    '
      vname(1,indxShflx_sen)='shflx_sen                           '
      vname(2,indxShflx_sen)='Sensible surface heat flux          '
      vname(3,indxShflx_sen)='Watts meter-2                       '
      vname(4,indxShflx_sen)='                                    '
      vname(4,indxShflx_sen)='                                    '
      vname(6,indxShflx_sen)='lat_rho lon_rho                     '
      vname(7,indxShflx_sen)='                                    '
      vname(7,indxShflx_sen)='  '
      vname(1,indxO)='omega                                       '
      vname(2,indxO)='S-coordinate vertical momentum component    '
      vname(3,indxO)='meter second-1                              '
      vname(4,indxO)='omega, scalar, series                       '
      vname(5,indxO)='                                            '
      vname(6,indxO)='lat_rho lon_rho                             '
      vname(7,indxO)='                                            '
      vname(1,indxW)='w                                           '
      vname(2,indxW)='vertical momentum component                 '
      vname(3,indxW)='meter second-1                              '
      vname(4,indxW)='w-velocity, scalar, series                  '
      vname(5,indxW)='upward_sea_water_velocity                   '
      vname(6,indxW)='lat_rho lon_rho                             '
      vname(1,indxR)='rho                                         '
      vname(2,indxR)='density anomaly                             '
      vname(3,indxR)='kilogram meter-3                            '
      vname(4,indxR)='density, scalar, series                     '
      vname(5,indxR)='sea_water_sigma_t                           '
      vname(7,indxR)='                                            '
      vname(1,indxbvf)='bvf                                       '
      vname(2,indxbvf)='Brunt Vaisala Frequency                   '
      vname(3,indxbvf)='second-1                                  '
      vname(4,indxbvf)='bvf, scalar, series                       '
      vname(5,indxbvf)='brunt_vaisala_frequency                   '
      vname(7,indxbvf)='                                          '
      vname(1,indxDiff)='diff3d                                   '
      vname(2,indxDiff)='horizontal diffusivity coefficient       '
      vname(3,indxDiff)='meter4 second-1                          '
      vname(5,indxDiff)='ocean_tracer_xy_biharmonic_diffusivity   '
      vname(4,indxDiff)='diff3d, scalar, series                   '
      vname(6,indxDiff)='lat_rho lon_rho                          '
      vname(7,indxDiff)='                                         '
      vname(1,indxAkv)='AKv                                       '
      vname(2,indxAkv)='vertical viscosity coefficient            '
      vname(3,indxAkv)='meter2 second-1                           '
      vname(4,indxAkv)='AKv, scalar, series                       '
      vname(5,indxAkv)='ocean_vertical_momentum_diffusivity_'
     &                                          // 'at_w_location '
      vname(6,indxAkv)='lat_rho lon_rho                           '
      vname(7,indxAkv)='                                          '
      vname(1,indxAkt)='AKt                                       '
      vname(2,indxAkt)='temperature vertical diffusion coefficient'
      vname(3,indxAkt)='meter2 second-1                           '
      vname(4,indxAkt)='AKt, scalar, series                       '
      vname(5,indxAkt)='ocean_vertical_heat_diffusivity_'
     &                                         //  'at_w_location '
      vname(6,indxAkt)='lat_rho lon_rho                           '
      vname(7,indxAkt)='                                          '
      vname(1,indxAks)='AKs                                       '
      vname(2,indxAks)='salinity vertical diffusion coefficient   '
      vname(3,indxAks)='meter2 second-1                           '
      vname(4,indxAks)='AKs, scalar, series                       '
      vname(5,indxAks)='ocean_vertical_salt_diffusivity_'
     &               / /  'at_w_location                          '
      vname(6,indxAks)='lat_rho lon_rho                           '
      vname(7,indxAks)='                                          '
      vname(1,indxHbl)='hbl                                       '
      vname(2,indxHbl)='depth of planetary boundary layer         '
      vname(3,indxHbl)='meter                                     '
      vname(4,indxHbl)='hbl, scalar, series                       '
      vname(5,indxHbl)='ocean_mixed_layer_thickness_defined_'
     &                                       // 'by_mixing_scheme '
      vname(6,indxHbl)='lat_rho lon_rho                           '
      vname(7,indxHbl)='                                          '
      vname(1,indxTke)='tke                                       '
      vname(2,indxTke)='turbulent kinetic energy                  '
      vname(3,indxTke)='meter2 second-2                           '
      vname(4,indxTke)='tke, scalar, series                       '
      vname(5,indxTke)='ocean_turbulent_kinetic_energy'
     &                                         // '_at_w_location '
      vname(6,indxTke)='lat_rho lon_rho                           '
      vname(7,indxTke)='                                          '
      vname(1,indxGls)='gls                                       '
      vname(2,indxGls)='turbulent generic length scale            '
      vname(3,indxGls)='meter3 second-2                           '
      vname(4,indxGls)='gls, scalar, series                       '
      vname(5,indxGls)='ocean_generic_length_scale'
     &                                         // '_at_w_location '
      vname(6,indxGls)='lat_rho lon_rho                           '
      vname(7,indxGls)='                                          '
      vname(1,indxLsc)='Lscale                                    '
      vname(2,indxLsc)='vertical mixing turbulent length scale    '
      vname(3,indxLsc)='meter                                     '
      vname(4,indxLsc)='Lscale, scalar, series                    '
      vname(5,indxLsc)='ocean_turbulent_length_scale'
     &                                         // '_at_w_location '
      vname(6,indxLsc)='lat_rho lon_rho                           '
      vname(7,indxLsc)='                                          '
      vname(1,indxSSH)='SSH                                       '
      vname(2,indxSSH)='sea surface height                        '
      vname(3,indxSSH)='meter                                     '
      vname(4,indxSSH)='SSH, scalar, series                       '
      vname(5,indxSSH)='sea_surface_height_above_sea_level        '
      vname(6,indxSSH)='lat_rho lon_rho                           '
      vname(7,indxSSH)='                                          '
      vname(1,indxSUSTR)='sustr                                   '
      vname(2,indxSUSTR)='surface u-momentum stress               '
      vname(3,indxSUSTR)='Newton meter-2                          '
      vname(4,indxSUSTR)='surface u-mom. stress, scalar, series   '
      vname(5,indxSUSTR)='surface_downward_x_stress               '
      vname(6,indxSUSTR)='lat_u lon_u                             '
      vname(7,indxSUSTR)='                                        '
      vname(1,indxSVSTR)='svstr                                   '
      vname(2,indxSVSTR)='surface v-momentum stress               '
      vname(3,indxSVSTR)='Newton meter-2                          '
      vname(4,indxSVSTR)='surface v-mom. stress, scalar, series   '
      vname(5,indxSVSTR)='surface_downward_y_stress               '
      vname(6,indxSVSTR)='lat_v lon_v                             '
      vname(7,indxSVSTR)='                                        '
      vname(1,indxWSPD)='wspd                                     '
      vname(2,indxWSPD)='surface wind speed 10 m                  '
      vname(3,indxWSPD)='meter second-1                           '
      vname(4,indxWSPD)='surface wind speed, scalar, series       '
      vname(5,indxWSPD)='wind_speed                               '
      vname(6,indxWSPD)='lat_rho lon_rho                          '
      vname(7,indxWSPD)='                                         '
      vname(1,indxTAIR)='tair                                     '
      vname(2,indxTAIR)='surface air temperature 2m               '
      vname(3,indxTAIR)='Celsius                                  '
      vname(4,indxTAIR)='surface air temperature, scalar, series  '
      vname(5,indxTAIR)='air_temperature_at_2m                    '
      vname(6,indxTAIR)='lat_rho lon_rho                          '
      vname(7,indxTAIR)='                                         '
      vname(1,indxRHUM)='rhum                                     '
      vname(2,indxRHUM)='surface air relative humidity 2m         '
      vname(3,indxRHUM)='fraction                                 '
      vname(4,indxRHUM)='surface relative humidity, scalar, series'
      vname(5,indxRHUM)='relative_humidity_at_2m                  '
      vname(6,indxRHUM)='lat_rho lon_rho                          '
      vname(7,indxRHUM)='                                         '
      vname(1,indxRADLW)='radlw_in                                '
      vname(2,indxRADLW)='downward longwave radiation             '
      vname(3,indxRADLW)='Watts meter-2                           '
      vname(4,indxRADLW)='downward longwave, scalar, series       '
      vname(5,indxRADLW)='surface_net_downward_longwave_flux      '
      vname(6,indxRADLW)='                                        '
      vname(7,indxRADLW)='                                        '
      vname(1,indxPRATE)='prate                                   '
      vname(2,indxPRATE)='surface precipitation rate              '
      vname(3,indxPRATE)='Kg meter-2 second-1                     '
      vname(4,indxPRATE)='precipitation rate, scalar, series      '
      vname(5,indxPRATE)='                                        '
      vname(6,indxPRATE)='lat_rho lon_rho                         '
      vname(7,indxPRATE)='                                        '
      vname(1,indxUWND)='uwnd                                     '
      vname(2,indxUWND)='surface u-wind speed 10 m                '
      vname(3,indxUWND)='meter second-1                           '
      vname(4,indxUWND)='surface wind speed, scalar, series       '
      vname(5,indxUWND)='x_wind                                   '
      vname(6,indxUWND)='lat_u lon_u                              '
      vname(7,indxUWND)='                                         '
      vname(1,indxVWND)='vwnd                                     '
      vname(2,indxVWND)='surface v-wind speed 10 m                '
      vname(3,indxVWND)='meter second-1                           '
      vname(4,indxVWND)='surface wind speed, scalar, series       '
      vname(5,indxVWND)='y_wind                                   '
      vname(6,indxVWND)='lat_v lon_v                              '
      vname(7,indxVWND)='                                         '
       write(vname1,*) 'temp'
       vname(1,indxsurft)=vname1
       vname(2,indxsurft)='SST                                 '
       write(vname3,*) 'Celsius'
       vname(3,indxsurft)=vname3
       write(vname1,*) 'salt'
       vname(1,indxsurfs)=vname1
       vname(2,indxsurfs)='SSS                                 '
       write(vname3,*) 'PSU'
       vname(3,indxsurfs)=vname3
       write(vname1,*) 'zeta'
       vname(1,indxsurfz)=vname1
       vname(2,indxsurfz)='SSH                                 '
       write(vname3,*) 'm'
       vname(3,indxsurfz)=vname3
       write(vname1,*) 'u'
       vname(1,indxsurfu)=vname1
       vname(2,indxsurfu)='surface u                           '
       write(vname3,*) 'meter second-1                         '
       vname(3,indxsurfu)=vname3
       write(vname1,*) 'v'
       vname(1,indxsurfv)=vname1
       vname(2,indxsurfv)='surface v                           '
       write(vname3,*) trim(unitt),'meter second-1             '
       vname(3,indxsurfv)=vname3
      do itrc=1,NT
         got_tclm(itrc)=.false.
      enddo
      return
      end
