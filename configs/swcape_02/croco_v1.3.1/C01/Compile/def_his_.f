      subroutine def_his (ncid, total_rec, ierr)
      implicit none
      logical create_new_file
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , b3dgrd(4)
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc,itrv
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
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_64bit_data
      integer*4 nf_cdf5
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_format_64bit_offset
      integer*4 nf_format_64bit_data
      integer*4 nf_format_cdf5
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_int64
      external        nf_put_att_int64
      integer*4         nf_get_att_int64
      external        nf_get_att_int64
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_compact
      parameter (nf_compact = 2)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_szip
      external nf_def_var_szip
      integer*4 nf_inq_var_szip
      external nf_inq_var_szip
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_def_var_filter
      external nf_def_var_filter
      integer*4 nf_inq_var_filter
      external nf_inq_var_filter
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer*4 nf_set_log_level
      external nf_set_log_level
      character*70 text
      ierr=0
      lstr=lenstr(hisname)
      if (nrpfhis.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfhis))
        call insert_time_index (hisname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create(hisname(1:lstr),nf_64bit_offset, ncid)
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) write(stdout,'(/3(1x,A)/)')
     &           'ERROR in def_his/avg:',
     &           'Cannot create netCDF file:', hisname(1:lstr)
          goto 99
        endif
        if (nrpfhis.eq.0) total_rec=0
        call put_global_atts (ncid, ierr)
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    6,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        b3dgrd(1)=r2dgrd(1)
        b3dgrd(2)=r2dgrd(2)
        b3dgrd(4)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        if (total_rec.le.1) then
          call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd
     &                                  ,r3dgrd, w3dgrd)
        endif
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 hisTstep)
        ierr=nf_put_att_text (ncid, hisTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, hisTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, hisTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, hisTime, indxTime, 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, hisTime2)
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, hisTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, hisTime2, indxTime2, 5,
     &       NF_REAL, ierr)
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                              NF_REAL, 3, r2dgrd, hisZ)
          lvar=lenstr(vname(2,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'long_name', lvar,
     &                                  vname(2,indxZ)(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
          lvar=lenstr(vname(4,indxZ))
          ierr=nf_put_att_text (ncid, hisZ, 'field',     lvar,
     &                                  vname(4,indxZ)(1:lvar))
          call nf_add_attribute(ncid, hisZ, indxZ, 5, NF_REAL,
     &         ierr)
        endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, hisUb)
          lvar=lenstr(vname(2,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'long_name', lvar,
     &                                  vname(2,indxUb)(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
          lvar=lenstr(vname(4,indxUb))
          ierr=nf_put_att_text (ncid, hisUb, 'field',    lvar,
     &                                  vname(4,indxUb)(1:lvar))
          call nf_add_attribute(ncid, hisUb, indxUb, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                              NF_REAL, 3, v2dgrd, hisVb)
          lvar=lenstr(vname(2,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'long_name', lvar,
     &                                  vname(2,indxVb)(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
          lvar=lenstr(vname(4,indxVb))
          ierr=nf_put_att_text (ncid, hisVb, 'field',     lvar,
     &                                  vname(4,indxVb)(1:lvar))
          call nf_add_attribute(ncid, hisVb, indxVb, 5, NF_REAL,
     &         ierr)
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar),
     &                             NF_REAL, 4, u3dgrd, hisU)
          lvar=lenstr(vname(2,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'long_name', lvar,
     &                                  vname(2,indxU)(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
          lvar=lenstr(vname(4,indxU))
          ierr=nf_put_att_text (ncid, hisU, 'field',     lvar,
     &                                  vname(4,indxU)(1:lvar))
          call nf_add_attribute(ncid, hisU, indxU, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar),
     &                             NF_REAL, 4, v3dgrd, hisV)
          lvar=lenstr(vname(2,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'long_name', lvar,
     &                                  vname(2,indxV)(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
          lvar=lenstr(vname(4,indxV))
          ierr=nf_put_att_text (ncid, hisV, 'field',     lvar,
     &                                  vname(4,indxV)(1:lvar))
          call nf_add_attribute(ncid, hisV, indxV, 5, NF_REAL, ierr)
        endif
        do itrc=1,NT
          if (wrthis(indxV+itrc)) then
            lvar=lenstr(vname(1,indxV+itrc))
            ierr=nf_def_var (ncid, vname(1,indxV+itrc)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisT(itrc))
            lvar=lenstr(vname(2,indxV+itrc))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'long_name',
     &                         lvar, vname(2,indxV+itrc)(1:lvar))
            lvar=lenstr(vname(3,indxV+itrc))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'units', lvar,
     &                               vname(3,indxV+itrc)(1:lvar))
            lvar=lenstr(vname(4,indxV+itrc))
            ierr=nf_put_att_text (ncid, hisT(itrc), 'field', lvar,
     &                               vname(4,indxV+itrc)(1:lvar))
            call nf_add_attribute(ncid,hisT(itrc),indxV+itrc,5,
     &           NF_REAL, ierr)
          endif
        enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncid, vname(1,indxR)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisR)
          lvar=lenstr(vname(2,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'long_name', lvar,
     &                                  vname(2,indxR)(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
          lvar=lenstr(vname(4,indxR))
          ierr=nf_put_att_text (ncid, hisR, 'field',     lvar,
     &                                  vname(4,indxR)(1:lvar))
          call nf_add_attribute(ncid, hisR, indxR, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_def_var (ncid, vname(1,indxbvf)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisbvf)
          lvar=lenstr(vname(2,indxbvf))
          ierr=nf_put_att_text (ncid, hisbvf, 'long_name', lvar,
     &                                  vname(2,indxbvf)(1:lvar))
          lvar=lenstr(vname(3,indxbvf))
          ierr=nf_put_att_text (ncid, hisbvf, 'units',     lvar,
     &                                  vname(3,indxbvf)(1:lvar))
          lvar=lenstr(vname(4,indxbvf))
          ierr=nf_put_att_text (ncid, hisR, 'field',     lvar,
     &                                  vname(4,indxbvf)(1:lvar))
          call nf_add_attribute(ncid, hisbvf, indxbvf, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncid, vname(1,indxO)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisO)
          lvar=lenstr(vname(2,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'long_name', lvar,
     &                                  vname(2,indxO)(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
          lvar=lenstr(vname(4,indxO))
          ierr=nf_put_att_text (ncid, hisO, 'field',     lvar,
     &                                  vname(4,indxO)(1:lvar))
          call nf_add_attribute(ncid, hisO, indxO, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar),
     &                               NF_REAL, 4, r3dgrd, hisW)
          lvar=lenstr(vname(2,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'long_name', lvar,
     &                                  vname(2,indxW)(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
          lvar=lenstr(vname(4,indxW))
          ierr=nf_put_att_text (ncid, hisW, 'field',     lvar,
     &                                  vname(4,indxW)(1:lvar))
          call nf_add_attribute(ncid, hisW, indxW, 5, NF_REAL, ierr)
        endif
        if (wrthis(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_def_var (ncid, vname(1,indxBostr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisBostr)
          lvar=lenstr(vname(2,indxBostr))
          ierr=nf_put_att_text (ncid, hisBostr, 'long_name', lvar,
     &                                 vname(2,indxBostr)(1:lvar))
          lvar=lenstr(vname(3,indxBostr))
          ierr=nf_put_att_text (ncid, hisBostr, 'units',     lvar,
     &                                 vname(3,indxBostr)(1:lvar))
          call nf_add_attribute(ncid, hisBostr,indxBostr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxBustr)) then
          lvar=lenstr(vname(1,indxBustr))
          ierr=nf_def_var (ncid, vname(1,indxBustr)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, hisBustr)
          lvar=lenstr(vname(2,indxBustr))
          ierr=nf_put_att_text (ncid, hisBustr, 'long_name', lvar,
     &                                 vname(2,indxBustr)(1:lvar))
          lvar=lenstr(vname(3,indxBustr))
          ierr=nf_put_att_text (ncid, hisBustr, 'units',     lvar,
     &                                 vname(3,indxBustr)(1:lvar))
          call nf_add_attribute(ncid, hisBustr,indxBustr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxBvstr)) then
          lvar=lenstr(vname(1,indxBvstr))
          ierr=nf_def_var (ncid, vname(1,indxBvstr)(1:lvar),
     &                             NF_REAL, 3, v2dgrd, hisBvstr)
          lvar=lenstr(vname(2,indxBvstr))
          ierr=nf_put_att_text (ncid, hisBvstr, 'long_name', lvar,
     &                                 vname(2,indxBvstr)(1:lvar))
          lvar=lenstr(vname(3,indxBvstr))
          ierr=nf_put_att_text (ncid, hisBvstr, 'units',     lvar,
     &                                 vname(3,indxBvstr)(1:lvar))
          call nf_add_attribute(ncid, hisBvstr,indxBvstr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_def_var (ncid, vname(1,indxWstr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisWstr)
          lvar=lenstr(vname(2,indxWstr))
          ierr=nf_put_att_text (ncid, hisWstr, 'long_name', lvar,
     &                                 vname(2,indxWstr)(1:lvar))
          lvar=lenstr(vname(3,indxWstr))
          ierr=nf_put_att_text (ncid, hisWstr, 'units',     lvar,
     &                                 vname(3,indxWstr)(1:lvar))
          call nf_add_attribute(ncid, hisWstr, indxWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_def_var (ncid, vname(1,indxUWstr)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, hisUWstr)
          lvar=lenstr(vname(2,indxUWstr))
          ierr=nf_put_att_text (ncid, hisUWstr, 'long_name', lvar,
     &                                 vname(2,indxUWstr)(1:lvar))
          lvar=lenstr(vname(3,indxUWstr))
          ierr=nf_put_att_text (ncid, hisUWstr, 'units',     lvar,
     &                                 vname(3,indxUWstr)(1:lvar))
          call nf_add_attribute(ncid, hisUWstr, indxUWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_def_var (ncid, vname(1,indxVWstr)(1:lvar),
     &                             NF_REAL, 3, v2dgrd, hisVWstr)
          lvar=lenstr(vname(2,indxVWstr))
          ierr=nf_put_att_text (ncid, hisVWstr, 'long_name', lvar,
     &                                 vname(2,indxVWstr)(1:lvar))
          lvar=lenstr(vname(3,indxVWstr))
          ierr=nf_put_att_text (ncid, hisVWstr, 'units',     lvar,
     &                                 vname(3,indxVWstr)(1:lvar))
          call nf_add_attribute(ncid, hisVWstr, indxVWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_def_var (ncid, vname(1,indxDiff)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, hisDiff)
          lvar=lenstr(vname(2,indxDiff))
          ierr=nf_put_att_text (ncid, hisDiff, 'long_name', lvar,
     &                                  vname(2,indxDiff)(1:lvar))
          lvar=lenstr(vname(4,indxDiff))
          ierr=nf_put_att_text (ncid, hisDiff, 'field',     lvar,
     &                                  vname(4,indxDiff)(1:lvar))
          call nf_add_attribute(ncid, hisDiff, indxDiff, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, hisAkv)
          lvar=lenstr(vname(2,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'long_name', lvar,
     &                                  vname(2,indxAkv)(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
          lvar=lenstr(vname(4,indxAkv))
          ierr=nf_put_att_text (ncid, hisAkv, 'field',     lvar,
     &                                  vname(4,indxAkv)(1:lvar))
          call nf_add_attribute(ncid, hisAkv, indxAkv, 5, NF_REAL,
     &                                                          ierr)
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAkt)
          lvar=lenstr(vname(2,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'long_name', lvar,
     &                                  vname(2,indxAkt)(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
          lvar=lenstr(vname(4,indxAkt))
          ierr=nf_put_att_text (ncid, hisAkt, 'field',     lvar,
     &                                  vname(4,indxAkt)(1:lvar))
          call nf_add_attribute(ncid, hisAkt, indxAkt, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisAks)
          lvar=lenstr(vname(2,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'long_name', lvar,
     &                                  vname(2,indxAks)(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
          lvar=lenstr(vname(4,indxAks))
          ierr=nf_put_att_text (ncid, hisAks, 'field',     lvar,
     &                                  vname(4,indxAks)(1:lvar))
          call nf_add_attribute(ncid, hisAks, indxAks, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrthis(indxHbl)) then
          lvar=lenstr(vname(1,indxHbl))
          ierr=nf_def_var (ncid, vname(1,indxHbl)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisHbl)
          lvar=lenstr(vname(2,indxHbl))
          ierr=nf_put_att_text (ncid, hisHbl, 'long_name', lvar,
     &                                  vname(2,indxHbl)(1:lvar))
          lvar=lenstr(vname(3,indxHbl))
          ierr=nf_put_att_text (ncid, hisHbl, 'units',     lvar,
     &                                  vname(3,indxHbl)(1:lvar))
          lvar=lenstr(vname(4,indxHbl))
          ierr=nf_put_att_text (ncid, hisHbl, 'field',     lvar,
     &                                  vname(4,indxHbl)(1:lvar))
          call nf_add_attribute(ncid, hisHbl, indxHbl, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrthis(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_def_var (ncid, vname(1,indxTke)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisTke)
          lvar=lenstr(vname(2,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'long_name', lvar,
     &                                  vname(2,indxTke)(1:lvar))
          lvar=lenstr(vname(3,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'units',     lvar,
     &                                  vname(3,indxTke)(1:lvar))
          lvar=lenstr(vname(4,indxTke))
          ierr=nf_put_att_text (ncid, hisTke, 'field',     lvar,
     &                                  vname(4,indxTke)(1:lvar))
          call nf_add_attribute(ncid, hisTke, indxTke, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrthis(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_def_var (ncid, vname(1,indxGls)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisGls)
          lvar=lenstr(vname(2,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'long_name', lvar,
     &                                  vname(2,indxGls)(1:lvar))
          lvar=lenstr(vname(3,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'units',     lvar,
     &                                  vname(3,indxGls)(1:lvar))
          lvar=lenstr(vname(4,indxGls))
          ierr=nf_put_att_text (ncid, hisGls, 'field',     lvar,
     &                                  vname(4,indxGls)(1:lvar))
          call nf_add_attribute(ncid, hisGls, indxGls, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrthis(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_def_var (ncid, vname(1,indxLsc)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, hisLsc)
          lvar=lenstr(vname(2,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'long_name', lvar,
     &                                  vname(2,indxLsc)(1:lvar))
          lvar=lenstr(vname(3,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'units',     lvar,
     &                                  vname(3,indxLsc)(1:lvar))
          lvar=lenstr(vname(4,indxLsc))
          ierr=nf_put_att_text (ncid, hisLsc, 'field',     lvar,
     &                                  vname(4,indxLsc)(1:lvar))
          call nf_add_attribute(ncid, hisLsc, indxLsc, 5, NF_REAL,
     &                                                             ierr)
        endif
        if (wrthis(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_def_var (ncid, vname(1,indxShflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisShflx)
          lvar=lenstr(vname(2,indxShflx))
          ierr=nf_put_att_text (ncid, hisShflx, 'long_name', lvar,
     &                                 vname(2,indxShflx)(1:lvar))
          lvar=lenstr(vname(3,indxShflx))
          ierr=nf_put_att_text (ncid, hisShflx, 'units',     lvar,
     &                                 vname(3,indxShflx)(1:lvar))
          call nf_add_attribute(ncid, hisShflx, indxShflx, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrthis(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_def_var (ncid, vname(1,indxSwflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, hisSwflx)
          lvar=lenstr(vname(2,indxSwflx))
          ierr=nf_put_att_text (ncid, hisSwflx, 'long_name', lvar,
     &                                 vname(2,indxSwflx)(1:lvar))
          lvar=lenstr(vname(3,indxSwflx))
          ierr=nf_put_att_text (ncid, hisSwflx, 'units',     lvar,
     &                                 vname(3,indxSwflx)(1:lvar))
          call nf_add_attribute(ncid, hisSwflx, indxSwflx, 5,
     &                          NF_REAL, ierr)
        endif
      if (wrthis(indxShflx_rsw)) then
        lvar=lenstr(vname(1,indxShflx_rsw))
        ierr=nf_def_var (ncid, vname(1,indxShflx_rsw)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, hisShflx_rsw)
          lvar=lenstr(vname(2,indxShflx_rsw))
          ierr=nf_put_att_text (ncid, hisShflx_rsw, 'long_name', lvar,
     &                                 vname(2,indxShflx_rsw)(1:lvar))
          lvar=lenstr(vname(3,indxShflx_rsw))
          ierr=nf_put_att_text (ncid, hisShflx_rsw, 'units',     lvar,
     &                                 vname(3,indxShflx_rsw)(1:lvar))
        call nf_add_attribute(ncid, hisShflx_rsw, indxShflx_rsw,5,
     &                                                   NF_REAL, ierr)
      endif
        if (wrthis(indxShflx_rlw)) then
          lvar=lenstr(vname(1,indxShflx_rlw))
          ierr=nf_def_var (ncid, vname(1,indxShflx_rlw)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, hisShflx_rlw)
          lvar=lenstr(vname(2,indxShflx_rlw))
          ierr=nf_put_att_text (ncid, hisShflx_rlw, 'long_name', lvar,
     &                               vname(2,indxShflx_rlw)(1:lvar))
          lvar=lenstr(vname(3,indxShflx_rlw))
          ierr=nf_put_att_text (ncid, hisShflx_rlw, 'units',     lvar,
     &                               vname(3,indxShflx_rlw)(1:lvar))
          call nf_add_attribute(ncid, hisShflx_rlw, indxShflx_rlw,5,
     &         NF_REAL, ierr)
        endif
      if (wrthis(indxShflx_lat)) then
        lvar=lenstr(vname(1,indxShflx_lat))
        ierr=nf_def_var (ncid, vname(1,indxShflx_lat)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, hisShflx_lat)
          lvar=lenstr(vname(2,indxShflx_lat))
          ierr=nf_put_att_text (ncid, hisShflx_lat, 'long_name', lvar,
     &                               vname(2,indxShflx_lat)(1:lvar))
          lvar=lenstr(vname(3,indxShflx_lat))
          ierr=nf_put_att_text (ncid, hisShflx_lat, 'units',     lvar,
     &                                 vname(3,indxShflx_lat)(1:lvar))
          call nf_add_attribute(ncid, hisShflx_lat, indxShflx_lat, 5,
     &         NF_REAL, ierr)
        endif
      if (wrthis(indxShflx_sen)) then
        lvar=lenstr(vname(1,indxShflx_sen))
        ierr=nf_def_var (ncid, vname(1,indxShflx_sen)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, hisShflx_sen)
          lvar=lenstr(vname(2,indxShflx_sen))
          ierr=nf_put_att_text (ncid, hisShflx_sen, 'long_name', lvar,
     &                                 vname(2,indxShflx_sen)(1:lvar))
          lvar=lenstr(vname(3,indxShflx_sen))
          ierr=nf_put_att_text (ncid, hisShflx_sen, 'units',     lvar,
     &                                 vname(3,indxShflx_sen)(1:lvar))
        call nf_add_attribute(ncid, hisShflx_sen, indxShflx_sen, 5,
     &       NF_REAL, ierr)
      endif
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &                'DEF_HIS/AVG - Created ',
     &                'new netCDF file ''', hisname(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (hisname(1:lstr), nf_write, ncid)
        if (ierr .eq. nf_noerr) then
          ierr=checkdims (ncid, hisname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfhis.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfhis))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_HIS/AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  hisname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfhis.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            if (mynode.eq.0) write(stdout,'(/1x,4A,2x,A,I4/)')
     &                  'DEF_HIS/AVG ERROR: ',
     &                  'Cannot open file ''', hisname(1:lstr), '''.'
     &                   ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', hisTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', hisname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),hisTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), hisname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),hisTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), hisname(1:lstr)
          goto 99
        endif
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), hisZ)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxZ)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), hisUb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxUb)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), hisVb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxVb)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_inq_varid (ncid,vname(1,indxBostr)(1:lvar),
     &                                                   hisBostr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxBostr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxBustr)) then
          lvar=lenstr(vname(1,indxBustr))
          ierr=nf_inq_varid (ncid,vname(1,indxBustr)(1:lvar),
     &                                                   hisBustr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBustr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxBvstr)) then
          lvar=lenstr(vname(1,indxBvstr))
          ierr=nf_inq_varid (ncid,vname(1,indxBvstr)(1:lvar),
     &                                                   hisBvstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBvstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxWstr)(1:lvar),
     &                                                   hisWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxUWstr)(1:lvar),
     &                                                   hisUWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxUWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxVWstr)(1:lvar),
     &                                                   hisVWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxVWstr)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), hisU)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxU)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), hisV)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxV)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        do itrc=1,NT
          if (wrthis(indxV+itrc)) then
            lvar=lenstr(vname(1,indxV+itrc))
            ierr=nf_inq_varid (ncid, vname(1,indxV+itrc)(1:lvar),
     &                                                 hisT(itrc))
            if (ierr .ne. nf_noerr) then
              write(stdout,1) vname(1,indxV+itrc)(1:lvar),
     &                                       hisname(1:lstr)
              goto 99
            endif
          endif
        enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncid, vname(1,indxR)(1:lvar), hisR)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxR)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_inq_varid (ncid, vname(1,indxbvf)(1:lvar), hisbvf)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxbvf)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncid, vname(1,indxO)(1:lvar), hisO)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxO)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncid, vname(1,indxW)(1:lvar), hisW)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxW)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_inq_varid (ncid, vname(1,indxDiff)(1:lvar), hisDiff)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxDiff)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), hisAkv)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkv)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), hisAkt)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkt)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), hisAks)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAks)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxHbl)) then
          lvar=lenstr(vname(1,indxHbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbl)(1:lvar), hisHbl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxHbl)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_inq_varid (ncid,vname(1,indxTke)(1:lvar), hisTke)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTke)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_inq_varid (ncid,vname(1,indxGls)(1:lvar), hisGls)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxGls)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_inq_varid (ncid,vname(1,indxLsc)(1:lvar), hisLsc)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxLsc)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx)(1:lvar),
     &                                                   hisShflx)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_inq_varid (ncid,vname(1,indxSwflx)(1:lvar),
     &                                                   hisSwflx)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxSwflx)(1:lvar), hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx_rsw)) then
          lvar=lenstr(vname(1,indxShflx_rsw))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_rsw)(1:lvar),
     &                                                hisShflx_rsw)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), 
     &                                                   hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx_rlw)) then
          lvar=lenstr(vname(1,indxShflx_rlw))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_rlw)(1:lvar),
     &                                               hisShflx_rlw)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_rlw)(1:lvar), 
     &                                                   hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx_lat)) then
          lvar=lenstr(vname(1,indxShflx_lat))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_lat)(1:lvar),
     &                                               hisShflx_lat)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_lat)(1:lvar), 
     &                                                   hisname(1:lstr)
            goto 99
          endif
        endif
        if (wrthis(indxShflx_sen)) then
          lvar=lenstr(vname(1,indxShflx_sen))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_sen)(1:lvar),
     &                                               hisShflx_sen)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_sen)(1:lvar), 
     &                                                   hisname(1:lstr)
            goto 99
          endif
        endif
      if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_HIS/AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      else
        ierr=nf_open (hisname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_HIS/AVG ERROR: ',
     &                'Cannot open file ''', hisname(1:lstr), '''.'
     &                 ,' mynode =', mynode
          goto 99
        endif
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_HIS/AVG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_HIS/AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
        if (total_rec.le.1) call wrt_grid (ncid, hisname, lstr)
  99  return
      end
      subroutine def_avg (ncid, total_rec, ierr)
      implicit none
      logical create_new_file
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , b3dgrd(4)
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc,itrv
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
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_64bit_data
      integer*4 nf_cdf5
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_format_64bit_offset
      integer*4 nf_format_64bit_data
      integer*4 nf_format_cdf5
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_int64
      external        nf_put_att_int64
      integer*4         nf_get_att_int64
      external        nf_get_att_int64
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_compact
      parameter (nf_compact = 2)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_szip
      external nf_def_var_szip
      integer*4 nf_inq_var_szip
      external nf_inq_var_szip
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_def_var_filter
      external nf_def_var_filter
      integer*4 nf_inq_var_filter
      external nf_inq_var_filter
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer*4 nf_set_log_level
      external nf_set_log_level
      character*70 text
      ierr=0
      lstr=lenstr(avgname)
      if (nrpfavg.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfavg))
        call insert_time_index (avgname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create(avgname(1:lstr),nf_64bit_offset, ncid)
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) write(stdout,'(/3(1x,A)/)')
     &           'ERROR in def_his/avg:',
     &           'Cannot create netCDF file:', avgname(1:lstr)
          goto 99
        endif
        if (nrpfavg.eq.0) total_rec=0
        call put_global_atts (ncid, ierr)
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    6,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        b3dgrd(1)=r2dgrd(1)
        b3dgrd(2)=r2dgrd(2)
        b3dgrd(4)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        if (total_rec.le.1) then
          call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd
     &                                  ,r3dgrd, w3dgrd)
        endif
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 avgTstep)
        ierr=nf_put_att_text (ncid, avgTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, avgTime)
        text='averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, avgTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, avgTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, avgTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, avgTime, indxTime, 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, avgTime2)
        text='averaged '/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, avgTime2, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, avgTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, avgTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, avgTime2, indxTime2, 5,
     &       NF_REAL, ierr)
        if (wrtavg(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                              NF_REAL, 3, r2dgrd, avgZ)
          text='averaged '/ /vname(2,indxZ)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgZ, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, avgZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
          lvar=lenstr(vname(4,indxZ))
          ierr=nf_put_att_text (ncid, avgZ, 'field',     lvar,
     &                                  vname(4,indxZ)(1:lvar))
          call nf_add_attribute(ncid, avgZ, indxZ, 5, NF_REAL,
     &         ierr)
        endif
        if (wrtavg(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, avgUb)
          text='averaged '/ /vname(2,indxUb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgUb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncid, avgUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
          lvar=lenstr(vname(4,indxUb))
          ierr=nf_put_att_text (ncid, avgUb, 'field',    lvar,
     &                                  vname(4,indxUb)(1:lvar))
          call nf_add_attribute(ncid, avgUb, indxUb, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                              NF_REAL, 3, v2dgrd, avgVb)
          text='averaged '/ /vname(2,indxVb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgVb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncid, avgVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
          lvar=lenstr(vname(4,indxVb))
          ierr=nf_put_att_text (ncid, avgVb, 'field',     lvar,
     &                                  vname(4,indxVb)(1:lvar))
          call nf_add_attribute(ncid, avgVb, indxVb, 5, NF_REAL,
     &         ierr)
        endif
        if (wrtavg(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar),
     &                             NF_REAL, 4, u3dgrd, avgU)
          text='averaged '/ /vname(2,indxU)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgU, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, avgU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
          lvar=lenstr(vname(4,indxU))
          ierr=nf_put_att_text (ncid, avgU, 'field',     lvar,
     &                                  vname(4,indxU)(1:lvar))
          call nf_add_attribute(ncid, avgU, indxU, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar),
     &                             NF_REAL, 4, v3dgrd, avgV)
          text='averaged '/ /vname(2,indxV)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgV, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, avgV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
          lvar=lenstr(vname(4,indxV))
          ierr=nf_put_att_text (ncid, avgV, 'field',     lvar,
     &                                  vname(4,indxV)(1:lvar))
          call nf_add_attribute(ncid, avgV, indxV, 5, NF_REAL, ierr)
        endif
        do itrc=1,NT
          if (wrtavg(indxV+itrc)) then
            lvar=lenstr(vname(1,indxV+itrc))
            ierr=nf_def_var (ncid, vname(1,indxV+itrc)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, avgT(itrc))
            text='averaged '/ /vname(2,indxV+itrc)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncid, avgT(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
            lvar=lenstr(vname(3,indxV+itrc))
            ierr=nf_put_att_text (ncid, avgT(itrc), 'units', lvar,
     &                               vname(3,indxV+itrc)(1:lvar))
            lvar=lenstr(vname(4,indxV+itrc))
            ierr=nf_put_att_text (ncid, avgT(itrc), 'field', lvar,
     &                               vname(4,indxV+itrc)(1:lvar))
            call nf_add_attribute(ncid,avgT(itrc),indxV+itrc,5,
     &           NF_REAL, ierr)
          endif
        enddo
        if (wrtavg(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncid, vname(1,indxR)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, avgR)
          text='averaged '/ /vname(2,indxR)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgR, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncid, avgR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
          lvar=lenstr(vname(4,indxR))
          ierr=nf_put_att_text (ncid, avgR, 'field',     lvar,
     &                                  vname(4,indxR)(1:lvar))
          call nf_add_attribute(ncid, avgR, indxR, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_def_var (ncid, vname(1,indxbvf)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, avgbvf)
          text='averaged '/ /vname(2,indxbvf)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgbvf, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxbvf))
          ierr=nf_put_att_text (ncid, avgbvf, 'units',     lvar,
     &                                  vname(3,indxbvf)(1:lvar))
          lvar=lenstr(vname(4,indxbvf))
          ierr=nf_put_att_text (ncid, avgR, 'field',     lvar,
     &                                  vname(4,indxbvf)(1:lvar))
          call nf_add_attribute(ncid, avgbvf, indxbvf, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncid, vname(1,indxO)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, avgO)
          text='averaged '/ /vname(2,indxO)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgO, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncid, avgO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
          lvar=lenstr(vname(4,indxO))
          ierr=nf_put_att_text (ncid, avgO, 'field',     lvar,
     &                                  vname(4,indxO)(1:lvar))
          call nf_add_attribute(ncid, avgO, indxO, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar),
     &                               NF_REAL, 4, r3dgrd, avgW)
          text='averaged '/ /vname(2,indxW)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgW, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, avgW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
          lvar=lenstr(vname(4,indxW))
          ierr=nf_put_att_text (ncid, avgW, 'field',     lvar,
     &                                  vname(4,indxW)(1:lvar))
          call nf_add_attribute(ncid, avgW, indxW, 5, NF_REAL, ierr)
        endif
        if (wrtavg(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_def_var (ncid, vname(1,indxBostr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, avgBostr)
          text='averaged '/ /vname(2,indxBostr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgBostr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxBostr))
          ierr=nf_put_att_text (ncid, avgBostr, 'units',     lvar,
     &                                 vname(3,indxBostr)(1:lvar))
          call nf_add_attribute(ncid, avgBostr,indxBostr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxBustr)) then
          lvar=lenstr(vname(1,indxBustr))
          ierr=nf_def_var (ncid, vname(1,indxBustr)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, avgBustr)
          text='averaged '/ /vname(2,indxBustr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgBustr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxBustr))
          ierr=nf_put_att_text (ncid, avgBustr, 'units',     lvar,
     &                                 vname(3,indxBustr)(1:lvar))
          call nf_add_attribute(ncid, avgBustr,indxBustr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxBvstr)) then
          lvar=lenstr(vname(1,indxBvstr))
          ierr=nf_def_var (ncid, vname(1,indxBvstr)(1:lvar),
     &                             NF_REAL, 3, v2dgrd, avgBvstr)
          text='averaged '/ /vname(2,indxBvstr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgBvstr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxBvstr))
          ierr=nf_put_att_text (ncid, avgBvstr, 'units',     lvar,
     &                                 vname(3,indxBvstr)(1:lvar))
          call nf_add_attribute(ncid, avgBvstr,indxBvstr,5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_def_var (ncid, vname(1,indxWstr)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, avgWstr)
          text='averaged '/ /vname(2,indxWstr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgWstr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxWstr))
          ierr=nf_put_att_text (ncid, avgWstr, 'units',     lvar,
     &                                 vname(3,indxWstr)(1:lvar))
          call nf_add_attribute(ncid, avgWstr, indxWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_def_var (ncid, vname(1,indxUWstr)(1:lvar),
     &                             NF_REAL, 3, u2dgrd, avgUWstr)
          text='averaged '/ /vname(2,indxUWstr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgUWstr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxUWstr))
          ierr=nf_put_att_text (ncid, avgUWstr, 'units',     lvar,
     &                                 vname(3,indxUWstr)(1:lvar))
          call nf_add_attribute(ncid, avgUWstr, indxUWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_def_var (ncid, vname(1,indxVWstr)(1:lvar),
     &                             NF_REAL, 3, v2dgrd, avgVWstr)
          text='averaged '/ /vname(2,indxVWstr)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgVWstr, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxVWstr))
          ierr=nf_put_att_text (ncid, avgVWstr, 'units',     lvar,
     &                                 vname(3,indxVWstr)(1:lvar))
          call nf_add_attribute(ncid, avgVWstr, indxVWstr, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_def_var (ncid, vname(1,indxDiff)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, avgDiff)
          text='averaged '/ /vname(2,indxDiff)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgDiff, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(4,indxDiff))
          ierr=nf_put_att_text (ncid, avgDiff, 'field',     lvar,
     &                                  vname(4,indxDiff)(1:lvar))
          call nf_add_attribute(ncid, avgDiff, indxDiff, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar),
     &                             NF_REAL, 4, w3dgrd, avgAkv)
          text='averaged '/ /vname(2,indxAkv)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgAkv, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, avgAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
          lvar=lenstr(vname(4,indxAkv))
          ierr=nf_put_att_text (ncid, avgAkv, 'field',     lvar,
     &                                  vname(4,indxAkv)(1:lvar))
          call nf_add_attribute(ncid, avgAkv, indxAkv, 5, NF_REAL,
     &                                                          ierr)
        endif
        if (wrtavg(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, avgAkt)
          text='averaged '/ /vname(2,indxAkt)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgAkt, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, avgAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
          lvar=lenstr(vname(4,indxAkt))
          ierr=nf_put_att_text (ncid, avgAkt, 'field',     lvar,
     &                                  vname(4,indxAkt)(1:lvar))
          call nf_add_attribute(ncid, avgAkt, indxAkt, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrtavg(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, avgAks)
          text='averaged '/ /vname(2,indxAks)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgAks, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, avgAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
          lvar=lenstr(vname(4,indxAks))
          ierr=nf_put_att_text (ncid, avgAks, 'field',     lvar,
     &                                  vname(4,indxAks)(1:lvar))
          call nf_add_attribute(ncid, avgAks, indxAks, 5, NF_REAL,
     &                                                           ierr)
        endif
        if (wrtavg(indxHbl)) then
          lvar=lenstr(vname(1,indxHbl))
          ierr=nf_def_var (ncid, vname(1,indxHbl)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, avgHbl)
          text='averaged '/ /vname(2,indxHbl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgHbl, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbl))
          ierr=nf_put_att_text (ncid, avgHbl, 'units',     lvar,
     &                                  vname(3,indxHbl)(1:lvar))
          lvar=lenstr(vname(4,indxHbl))
          ierr=nf_put_att_text (ncid, avgHbl, 'field',     lvar,
     &                                  vname(4,indxHbl)(1:lvar))
          call nf_add_attribute(ncid, avgHbl, indxHbl, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrtavg(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_def_var (ncid, vname(1,indxTke)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, avgTke)
          text='averaged '/ /vname(2,indxTke)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgTke, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxTke))
          ierr=nf_put_att_text (ncid, avgTke, 'units',     lvar,
     &                                  vname(3,indxTke)(1:lvar))
          lvar=lenstr(vname(4,indxTke))
          ierr=nf_put_att_text (ncid, avgTke, 'field',     lvar,
     &                                  vname(4,indxTke)(1:lvar))
          call nf_add_attribute(ncid, avgTke, indxTke, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrtavg(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_def_var (ncid, vname(1,indxGls)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, avgGls)
          text='averaged '/ /vname(2,indxGls)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgGls, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxGls))
          ierr=nf_put_att_text (ncid, avgGls, 'units',     lvar,
     &                                  vname(3,indxGls)(1:lvar))
          lvar=lenstr(vname(4,indxGls))
          ierr=nf_put_att_text (ncid, avgGls, 'field',     lvar,
     &                                  vname(4,indxGls)(1:lvar))
          call nf_add_attribute(ncid, avgGls, indxGls, 5, NF_REAL,
     &                                                            ierr)
        endif
        if (wrtavg(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_def_var (ncid, vname(1,indxLsc)(1:lvar),
     &                              NF_REAL, 4, w3dgrd, avgLsc)
          text='averaged '/ /vname(2,indxLsc)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgLsc, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxLsc))
          ierr=nf_put_att_text (ncid, avgLsc, 'units',     lvar,
     &                                  vname(3,indxLsc)(1:lvar))
          lvar=lenstr(vname(4,indxLsc))
          ierr=nf_put_att_text (ncid, avgLsc, 'field',     lvar,
     &                                  vname(4,indxLsc)(1:lvar))
          call nf_add_attribute(ncid, avgLsc, indxLsc, 5, NF_REAL,
     &                                                             ierr)
        endif
        if (wrtavg(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_def_var (ncid, vname(1,indxShflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, avgShflx)
          text='averaged '/ /vname(2,indxShflx)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgShflx, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxShflx))
          ierr=nf_put_att_text (ncid, avgShflx, 'units',     lvar,
     &                                 vname(3,indxShflx)(1:lvar))
          call nf_add_attribute(ncid, avgShflx, indxShflx, 5,
     &                          NF_REAL, ierr)
        endif
        if (wrtavg(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_def_var (ncid, vname(1,indxSwflx)(1:lvar),
     &                             NF_REAL, 3, r2dgrd, avgSwflx)
          text='averaged '/ /vname(2,indxSwflx)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgSwflx, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxSwflx))
          ierr=nf_put_att_text (ncid, avgSwflx, 'units',     lvar,
     &                                 vname(3,indxSwflx)(1:lvar))
          call nf_add_attribute(ncid, avgSwflx, indxSwflx, 5,
     &                          NF_REAL, ierr)
        endif
      if (wrtavg(indxShflx_rsw)) then
        lvar=lenstr(vname(1,indxShflx_rsw))
        ierr=nf_def_var (ncid, vname(1,indxShflx_rsw)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, avgShflx_rsw)
          text='averaged '/ /vname(2,indxShflx_rsw)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgShflx_rsw, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxShflx_rsw))
          ierr=nf_put_att_text (ncid, avgShflx_rsw, 'units',     lvar,
     &                                 vname(3,indxShflx_rsw)(1:lvar))
        call nf_add_attribute(ncid, avgShflx_rsw, indxShflx_rsw,5,
     &                                                   NF_REAL, ierr)
      endif
        if (wrtavg(indxShflx_rlw)) then
          lvar=lenstr(vname(1,indxShflx_rlw))
          ierr=nf_def_var (ncid, vname(1,indxShflx_rlw)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, avgShflx_rlw)
          text='averaged '/ /vname(2,indxShflx_rlw)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgShflx_rlw, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxShflx_rlw))
          ierr=nf_put_att_text (ncid, avgShflx_rlw, 'units',     lvar,
     &                               vname(3,indxShflx_rlw)(1:lvar))
          call nf_add_attribute(ncid, avgShflx_rlw, indxShflx_rlw,5,
     &         NF_REAL, ierr)
        endif
      if (wrtavg(indxShflx_lat)) then
        lvar=lenstr(vname(1,indxShflx_lat))
        ierr=nf_def_var (ncid, vname(1,indxShflx_lat)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, avgShflx_lat)
          text='averaged '/ /vname(2,indxShflx_lat)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgShflx_lat, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxShflx_lat))
          ierr=nf_put_att_text (ncid, avgShflx_lat, 'units',     lvar,
     &                                 vname(3,indxShflx_lat)(1:lvar))
          call nf_add_attribute(ncid, avgShflx_lat, indxShflx_lat, 5,
     &         NF_REAL, ierr)
        endif
      if (wrtavg(indxShflx_sen)) then
        lvar=lenstr(vname(1,indxShflx_sen))
        ierr=nf_def_var (ncid, vname(1,indxShflx_sen)(1:lvar),
     &                           NF_REAL, 3, r2dgrd, avgShflx_sen)
          text='averaged '/ /vname(2,indxShflx_sen)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, avgShflx_sen, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxShflx_sen))
          ierr=nf_put_att_text (ncid, avgShflx_sen, 'units',     lvar,
     &                                 vname(3,indxShflx_sen)(1:lvar))
        call nf_add_attribute(ncid, avgShflx_sen, indxShflx_sen, 5,
     &       NF_REAL, ierr)
      endif
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &                'DEF_HIS/AVG - Created ',
     &                'new netCDF file ''', avgname(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (avgname(1:lstr), nf_write, ncid)
        if (ierr .eq. nf_noerr) then
          ierr=checkdims (ncid, avgname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfavg.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfavg))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_HIS/AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  avgname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfavg.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            if (mynode.eq.0) write(stdout,'(/1x,4A,2x,A,I4/)')
     &                  'DEF_HIS/AVG ERROR: ',
     &                  'Cannot open file ''', avgname(1:lstr), '''.'
     &                   ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', avgTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', avgname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),avgTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), avgname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),avgTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), avgname(1:lstr)
          goto 99
        endif
        if (wrtavg(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), avgZ)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxZ)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), avgUb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxUb)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), avgVb)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxVb)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxBostr)) then
          lvar=lenstr(vname(1,indxBostr))
          ierr=nf_inq_varid (ncid,vname(1,indxBostr)(1:lvar),
     &                                                   avgBostr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxBostr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxBustr)) then
          lvar=lenstr(vname(1,indxBustr))
          ierr=nf_inq_varid (ncid,vname(1,indxBustr)(1:lvar),
     &                                                   avgBustr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBustr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxBvstr)) then
          lvar=lenstr(vname(1,indxBvstr))
          ierr=nf_inq_varid (ncid,vname(1,indxBvstr)(1:lvar),
     &                                                   avgBvstr)
          if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxBvstr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxWstr)) then
          lvar=lenstr(vname(1,indxWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxWstr)(1:lvar),
     &                                                   avgWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxWstr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxUWstr)) then
          lvar=lenstr(vname(1,indxUWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxUWstr)(1:lvar),
     &                                                   avgUWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxUWstr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxVWstr)) then
          lvar=lenstr(vname(1,indxVWstr))
          ierr=nf_inq_varid (ncid,vname(1,indxVWstr)(1:lvar),
     &                                                   avgVWstr)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxVWstr)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), avgU)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxU)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), avgV)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxV)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        do itrc=1,NT
          if (wrtavg(indxV+itrc)) then
            lvar=lenstr(vname(1,indxV+itrc))
            ierr=nf_inq_varid (ncid, vname(1,indxV+itrc)(1:lvar),
     &                                                 avgT(itrc))
            if (ierr .ne. nf_noerr) then
              write(stdout,1) vname(1,indxV+itrc)(1:lvar),
     &                                       avgname(1:lstr)
              goto 99
            endif
          endif
        enddo
        if (wrtavg(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncid, vname(1,indxR)(1:lvar), avgR)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxR)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxbvf)) then
          lvar=lenstr(vname(1,indxbvf))
          ierr=nf_inq_varid (ncid, vname(1,indxbvf)(1:lvar), avgbvf)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxbvf)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncid, vname(1,indxO)(1:lvar), avgO)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxO)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncid, vname(1,indxW)(1:lvar), avgW)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxW)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxDiff)) then
          lvar=lenstr(vname(1,indxDiff))
          ierr=nf_inq_varid (ncid, vname(1,indxDiff)(1:lvar), avgDiff)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxDiff)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), avgAkv)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkv)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), avgAkt)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAkt)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), avgAks)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxAks)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxHbl)) then
          lvar=lenstr(vname(1,indxHbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbl)(1:lvar), avgHbl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxHbl)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxTke)) then
          lvar=lenstr(vname(1,indxTke))
          ierr=nf_inq_varid (ncid,vname(1,indxTke)(1:lvar), avgTke)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTke)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxGls)) then
          lvar=lenstr(vname(1,indxGls))
          ierr=nf_inq_varid (ncid,vname(1,indxGls)(1:lvar), avgGls)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxGls)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxLsc)) then
          lvar=lenstr(vname(1,indxLsc))
          ierr=nf_inq_varid (ncid,vname(1,indxLsc)(1:lvar), avgLsc)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxLsc)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxShflx)) then
          lvar=lenstr(vname(1,indxShflx))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx)(1:lvar),
     &                                                   avgShflx)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxSwflx)) then
          lvar=lenstr(vname(1,indxSwflx))
          ierr=nf_inq_varid (ncid,vname(1,indxSwflx)(1:lvar),
     &                                                   avgSwflx)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxSwflx)(1:lvar), avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxShflx_rsw)) then
          lvar=lenstr(vname(1,indxShflx_rsw))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_rsw)(1:lvar),
     &                                                avgShflx_rsw)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), 
     &                                                   avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxShflx_rlw)) then
          lvar=lenstr(vname(1,indxShflx_rlw))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_rlw)(1:lvar),
     &                                               avgShflx_rlw)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_rlw)(1:lvar), 
     &                                                   avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxShflx_lat)) then
          lvar=lenstr(vname(1,indxShflx_lat))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_lat)(1:lvar),
     &                                               avgShflx_lat)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_lat)(1:lvar), 
     &                                                   avgname(1:lstr)
            goto 99
          endif
        endif
        if (wrtavg(indxShflx_sen)) then
          lvar=lenstr(vname(1,indxShflx_sen))
          ierr=nf_inq_varid (ncid,vname(1,indxShflx_sen)(1:lvar),
     &                                               avgShflx_sen)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxShflx_sen)(1:lvar), 
     &                                                   avgname(1:lstr)
            goto 99
          endif
        endif
      if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_HIS/AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      else
        ierr=nf_open (avgname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          if (mynode.eq.0) write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_HIS/AVG ERROR: ',
     &                'Cannot open file ''', avgname(1:lstr), '''.'
     &                 ,' mynode =', mynode
          goto 99
        endif
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_HIS/AVG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_HIS/AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
        if (total_rec.le.1) call wrt_grid (ncid, avgname, lstr)
  99  return
      end
