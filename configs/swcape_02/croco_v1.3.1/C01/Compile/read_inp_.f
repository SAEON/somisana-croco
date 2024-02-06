      subroutine read_inp (ierr)
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
      include 'mpif.h'
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      integer*4 kwsize, testunit, input
      parameter (kwsize=32, testunit=40, input=15)
      character end_signal*3, keyword*32, fname*180
      parameter (end_signal='end')
      integer*4 ierr, iargc, is,ie, kwlen, lstr, lenstr
     &                                       , itrc
      logical dumboolean
      fname='croco.in'
      if (mynode.eq.0 .and. iargc().GT.0) call getarg(1,fname)
      call MPI_Bcast(fname,64,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
      wrthis(indxTime)=.false.
      wrtavg(indxTime)=.false.
      wrtsurf(3)=.false.
      wrtsurf_avg(3)=.false.
      ierr=0
      call setup_kwds (ierr)
      open (input,file=fname,status='old',form='formatted',err=97)
   1  keyword='                                '
      read(input,'(A)',err=1,end=99) keyword
      if (ichar(keyword(1:1)).eq.33) goto 1
      is=1
   2  if (is.eq.kwsize) then
        goto 1
      elseif (keyword(is:is).eq.' ') then
        is=is+1
        goto 2
      endif
      ie=is
   3  if (keyword(ie:ie).eq.':') then
        keyword(ie:ie)=' '
        goto 4
      elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
        ie=ie+1
        goto 3
      endif
      goto 1
   4  kwlen=ie-is
      if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
      if (keyword(1:kwlen).eq.'title') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) title
        lstr=lenstr(title)
        if (mynode.eq.0) write(stdout,'(/1x,A)') title(1:lstr)
      elseif (keyword(1:kwlen).eq.'time_stepping') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntimes,dt,ndtfast, ninfo
        if (mynode.eq.0) write(stdout,
     & '(I10,2x,A,1x,A /F10.2,2x,A,2(/I10,2x,A,1x,A)/F10.4,2x,A)')
     &    ntimes,  'ntimes   Total number of timesteps for',
     &                                            '3D equations.',
     &    dt,      'dt       Timestep [sec] for 3D equations',
     &    ndtfast, 'ndtfast  Number of 2D timesteps within each',
     &                                                 '3D step.',
     &    ninfo,   'ninfo    Number of timesteps between',
     &                                     'runtime diagnostics.'
        dtfast=dt/float(ndtfast)
        if (NWEIGHT.lt.(2*ndtfast-1)) then
          write(stdout,'(a,i3)')
     &    ' Error - Number of 2D timesteps (2*ndtfast-1): ',
     &    2*ndtfast-1
          write(stdout,'(a,i3)')
     &    '           exceeds barotopic weight dimension: ',NWEIGHT
          goto 95
        endif
      elseif (keyword(1:kwlen).eq.'S-coord') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) theta_s, theta_b, Tcline
        if (mynode.eq.0) write(stdout,
     &                        '(3(1pe10.3,2x,A,1x,A/),32x,A)')
     &    theta_s, 'theta_s  S-coordinate surface control',
     &                                               'parameter.',
     &    theta_b, 'theta_b  S-coordinate bottom control',
     &                                               'parameter.',
     &    Tcline,  'Tcline   S-coordinate surface/bottom layer',
     &  'width used in', 'vertical coordinate stretching, meters.'
      elseif (keyword(1:kwlen).eq.'initial') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrrec
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          ininame=fname(1:lstr)
          if (mynode.eq.0) write(stdout,'(1x,A,2x,A,4x,A,I3)')
     &     'Initial State File:', ininame(1:lstr), 'Record:',nrrec
      elseif (keyword(1:kwlen).eq.'grid') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open(testunit,file=fname(1:lstr), status='old', err=97)
        close(testunit)
        grdname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,'(10x,A,2x,A)')
     &                   'Grid File:', grdname(1:lstr)
        elseif (keyword(1:kwlen).eq.'bulk_forcing') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          bulkname=fname(1:lstr)
          if (mynode.eq.0) write(stdout,'(2x,A,2x,A)')
     &               '   Bulk Data File:', bulkname(1:lstr)
      elseif (keyword(1:kwlen).eq.'climatology') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        clmname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,'(3x,A,2x,A)')
     &            'Climatology File:', clmname(1:lstr)
      elseif (keyword(1:kwlen).eq.'restart') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrst, nrpfrst
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
        rstname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,4x,A,I6,4x,A,I4)')
     &             'Restart File:', rstname(1:lstr),
     &             'nrst =', nrst, 'rec/file: ', nrpfrst
      elseif (keyword(1:kwlen).eq.'history') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefhis, nwrt, nrpfhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        hisname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,2x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &       'History File:', hisname(1:lstr),  'Create new:',
     &       ldefhis, 'nwrt =', nwrt, 'rec/file =', nrpfhis
      elseif (keyword(1:kwlen).eq.'averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntsavg, navg, nrpfavg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        avgname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &         '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      ntsavg, 'ntsavg      Starting timestep for the',
     &         'accumulation of output', 'time-averaged data.',
     &      navg,   'navg        Number of timesteps between',
     &     'writing of time-averaged','data into averages file.',
     &     'Averages File:', avgname(1:lstr),
     &     'rec/file =', nrpfavg
      elseif (keyword(1:kwlen).eq.'surf') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefsurf, nwrtsurf, nrpfsurf
        if (nwrtsurf.eq.0) nwrtsurf = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        surfname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &    '(5x,A,2x,A/,32x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &    'Surface outputs File:', surfname(1:lstr),
     &    'Create new:', ldefsurf, 'nwrt =', nwrtsurf,
     &    'rec/file =', nrpfsurf
      elseif (keyword(1:kwlen).eq.'surf_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefsurf_avg, ntssurf_avg,
     &                      nwrtsurf_avg,  nrpfsurf_avg
        if (nwrtsurf_avg.eq.0) nwrtsurf_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        surfname_avg=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &    '(5x,A,2x,A/,32x,A,1x,L1,2x,A,I4,2x,A,I3,/32x,A,I10)')
     &    'Surface outputs AVG File:',surfname_avg(1:lstr),
     &    'Create new:', ldefsurf_avg,'nwrt =',nwrtsurf_avg,
     &    'rec/file =',nrpfsurf_avg,
     &    'Starting timestep = ',ntssurf_avg
      elseif (keyword(1:kwlen).eq.'primary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                       ,  wrthis(indxVb)
     &                    ,  wrthis(indxU),  wrthis(indxV)
     &                    , (wrthis(itrc), itrc=indxV+1,indxV+NT)
        if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
     &                        .or. wrthis(indxU) .or. wrthis(indxV)
     &     ) wrthis(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &    , wrthis(indxZ),  'write zeta ', 'free-surface.'
     &    , wrthis(indxUb), 'write UBAR ', '2D U-momentum component.'
     &    , wrthis(indxVb), 'write VBAR ', '2D V-momentum component.'
     &    , wrthis(indxU),  'write U    ', '3D U-momentum component.'
     &    , wrthis(indxV),  'write V    ', '3D V-momentum component.'
        do itrc=1,NT
          if (wrthis(indxV+itrc)) wrthis(indxTime)=.true.
          if (mynode.eq.0) write(stdout, '(6x,L1,2x,A,I2,A,I2,A)')
     &                     wrthis(indxV+itrc), 'write T(', itrc,
     &                              ')  Tracer of index', itrc,'.'
        enddo
      elseif (keyword(1:kwlen).eq.'auxiliary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                                             wrthis(indxR)
     &                                          ,  wrthis(indxO)
     &                                          ,  wrthis(indxW)
     &                                          ,  wrthis(indxAkv)
     &                                          ,  wrthis(indxAkt)
     &                                          ,  wrthis(indxAks)
     &                                          ,  wrthis(indxbvf)
     &                                          ,  dumboolean
     &                                          ,  wrthis(indxDiff)
     &                                          ,  wrthis(indxHbl)
     &                                          ,  dumboolean
     &                                          ,  wrthis(indxBostr)
     &                                          ,  wrthis(indxBustr)
     &                                          ,  wrthis(indxBvstr)
     &                                          ,  wrthis(indxWstr)
     &                                          ,  wrthis(indxUWstr)
     &                                          ,  wrthis(indxVWstr)
     &                                          ,  wrthis(indxShflx)
     &                                          ,  wrthis(indxSwflx)
     &                                          ,  wrthis(indxShflx_rsw)
     &                                          ,  wrthis(indxShflx_rlw)
     &                                          ,  wrthis(indxShflx_lat)
     &                                          ,  wrthis(indxShflx_sen)
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
        if ( wrthis(indxR)
     &                                        .or. wrthis(indxO)
     &                                        .or. wrthis(indxW)
     &                                        .or. wrthis(indxAkv)
     &                                        .or. wrthis(indxAkt)
     &                                        .or. wrthis(indxAks)
     &                                        .or. wrthis(indxbvf)
     &                                        .or. wrthis(indxDiff)
     &                                        .or. wrthis(indxHbl)
     &                                        .or. wrthis(indxBostr)
     &                                        .or. wrthis(indxBustr)
     &                                        .or. wrthis(indxBvstr)
     &                                        .or. wrthis(indxWstr)
     &                                        .or. wrthis(indxUWstr)
     &                                        .or. wrthis(indxVWstr)
     &                                        .or. wrthis(indxShflx)
     &                                        .or. wrthis(indxSwflx)
     &                                        .or. wrthis(indxShflx_rsw)
     &                                        .or. wrthis(indxShflx_rlw)
     &                                        .or. wrthis(indxShflx_lat)
     &                                        .or. wrthis(indxShflx_sen)
     &     ) wrthis(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrthis(indxR),    'write RHO  ', 'Density anomaly.'
     &  , wrthis(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrthis(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrthis(indxAkv),  'write Akv  ', 'Vertical viscosity.'
     &  , wrthis(indxAkt),  'write Akt  ',
     &                      'Vertical diffusivity for temperature.'
     &  , wrthis(indxAks),  'write Aks  ',
     &                      'Vertical diffusivity for salinity.'
     &  , wrthis(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
     &  , wrthis(indxDiff),  'write Visc3d', 'Horizontal diffusivity.'
     &  , wrthis(indxHbl),  'write Hbl  ',
     &                      'Depth of model boundary layer.'
     &  , wrthis(indxShflx_rlw), 'write shflx_rlw [W/m2]',
     &                                 'Long Wave heat flux.'
     &  , wrthis(indxShflx_lat), 'write shflx_lat [W/m2]',
     &                                 'Latent heat flux.'
     &  , wrthis(indxShflx_sen), 'write shflx_sen [W/m2]',
     &                                 'Sensible heat flux'
     &  , wrthis(indxBostr), 'write Bostr', 'Bottom Stress.'
     &  , wrthis(indxBustr), 'write Bustr', 'U-Bottom Stress.'
     &  , wrthis(indxBvstr), 'write Bvstr', 'V-Bottom Stress.'
     &  , wrthis(indxWstr),  'write Wstress', 'Wind Stress.'
     &  , wrthis(indxUWstr), 'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrthis(indxVWstr), 'write V-Wstress comp.', 'V-Wind Stress.'
     &  , wrthis(indxShflx), 'write Shflx [W/m2]',
     &                       'Surface net heat flux'
     &  , wrthis(indxSwflx), 'write Swflx [cm/day]',
     &                       'Surface freshwater flux (E-P)'
     &  , wrthis(indxShflx_rsw),'write Shflx_rsw [W/m2]',
     &                          'Short-wave surface radiation'
      elseif (keyword(1:kwlen).eq.'gls_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxTke),  wrthis(indxGls)
     &                    ,  wrthis(indxLsc)
        if (wrthis(indxTke) .or. wrthis(indxGls) .or. wrthis(indxLsc)
     &     ) wrthis(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &   , wrthis(indxTke), 'write TKE ', 'turbulent kinetic energy.  '
     &   , wrthis(indxGls), 'write GLS ', 'generic length scale.'
     &   , wrthis(indxLsc), 'write Lscale ',
     &                                  'vertical mixing length scale.'
      elseif (keyword(1:kwlen).eq.'primary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxZ),  wrtavg(indxUb)
     &                                    ,  wrtavg(indxVb)
     &                    ,  wrtavg(indxU),  wrtavg(indxV)
     &                    , (wrtavg(itrc), itrc=indxV+1,indxV+NT)
        if ( wrtavg(indxZ) .or. wrtavg(indxUb) .or. wrtavg(indxVb)
     &                     .or. wrtavg(indxU)  .or. wrtavg(indxV)
     &     ) wrtavg(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &  'Fields to be saved in averages file: (T/F)'
     &  , wrtavg(indxZ),  'write zeta ', 'free-surface.'
     &  , wrtavg(indxUb), 'write UBAR ', '2D U-momentum component.'
     &  , wrtavg(indxVb), 'write VBAR ', '2D V-momentum component.'
     &  , wrtavg(indxU),  'write U    ', '3D U-momentum component.'
     &  , wrtavg(indxV),  'write V    ', '3D V-momentum component.'
         do itrc=1,NT
          if (wrtavg(indxV+itrc)) wrtavg(indxTime)=.true.
          if (mynode.eq.0) write(stdout,
     &                     '(6x,L1,2x,A,I2,A,2x,A,I2,A)')
     &                      wrtavg(indxV+itrc), 'write T(',
     &                      itrc,')', 'Tracer of index', itrc,'.'
        enddo
      elseif (keyword(1:kwlen).eq.'auxiliary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxR), wrtavg(indxO)
     &        ,  wrtavg(indxW),  wrtavg(indxAkv)
     &                                          ,  wrtavg(indxAkt)
     &                                          ,  wrtavg(indxAks)
     &                                          ,  wrtavg(indxbvf)
     &                                          ,  dumboolean
     &                                          ,  wrtavg(indxDiff)
     &                                          ,  wrtavg(indxHbl)
     &                                          ,  dumboolean
     &                                          ,  wrtavg(indxBostr)
     &                                          ,  wrtavg(indxBustr)
     &                                          ,  wrtavg(indxBvstr)
     &                                          ,  wrtavg(indxWstr)
     &                                          ,  wrtavg(indxUWstr)
     &                                          ,  wrtavg(indxVWstr)
     &                                          ,  wrtavg(indxShflx)
     &                                          ,  wrtavg(indxSwflx)
     &                                          ,  wrtavg(indxShflx_rsw)
     &                                          ,  wrtavg(indxShflx_rlw)
     &                                          ,  wrtavg(indxShflx_lat)
     &                                          ,  wrtavg(indxShflx_sen)
     &                                          , dumboolean
     &                                          , dumboolean
     &                                          ,  dumboolean
        if ( wrtavg(indxR) .or. wrtavg(indxO) .or. wrtavg(indxW)
     &                   .or. wrtavg(indxAkv)
     &                                        .or. wrtavg(indxAkt)
     &                                        .or. wrtavg(indxAks)
     &                                        .or. wrtavg(indxbvf)
     &                                        .or. wrtavg(indxDiff)
     &                                        .or. wrtavg(indxHbl)
     &                                        .or. wrtavg(indxBostr)
     &                                        .or. wrtavg(indxBustr)
     &                                        .or. wrtavg(indxBvstr)
     &                                        .or. wrtavg(indxWstr)
     &                                        .or. wrtavg(indxUWstr)
     &                                        .or. wrtavg(indxVWstr)
     &                                        .or. wrtavg(indxShflx)
     &                                        .or. wrtavg(indxSwflx)
     &                                        .or. wrtavg(indxShflx_rsw)
     &                                        .or. wrtavg(indxShflx_rlw)
     &                                        .or. wrtavg(indxShflx_lat)
     &                                        .or. wrtavg(indxShflx_sen)
     &     ) wrtavg(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrtavg(indxR),    'write RHO  ', 'Density anomaly'
     &  , wrtavg(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrtavg(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrtavg(indxAkv),  'write Akv  ', 'Vertical viscosity'
     &  , wrtavg(indxAkt),  'write Akt  ',
     &                      'Vertical diffusivity for temperature.'
     &  , wrtavg(indxAks),  'write Aks  ',
     &                         'Vertical diffusivity for salinity.'
     &  , wrtavg(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
     &  , wrtavg(indxDiff),'write diff3d', 'Horizontal diffusivity'
     &  , wrtavg(indxHbl),  'write Hbl  ',
     &                          'Depth of model boundary layer'
     &  , wrtavg(indxShflx_rlw), 'write shflx_rlw [W/m2]',
     &                                 'Long Wave heat flux.'
     &  , wrtavg(indxShflx_lat), 'write shflx_lat[W/m2] ',
     &                                 'Latente heat flux.'
     &  , wrtavg(indxShflx_sen), 'write shflx_sen [W/m2]',
     &                                 'Sensible heat flux.'
     &  , wrtavg(indxBostr),'write Bostr', 'Bottom Stress.'
     &  , wrtavg(indxBustr),'write Bustr', 'U-Bottom Stress.'
     &  , wrtavg(indxBvstr),'write Bvstr', 'V-Bottom Stress.'
     &  , wrtavg(indxWstr), 'write Wstr', 'Wind Stress.'
     &  , wrtavg(indxUWstr),'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrtavg(indxVWstr),'write V-Wstress comp.', 'V-Wind Stress.'
     &  , wrtavg(indxShflx),'write Shflx [W/m2]',
     &                      'Surface net heat flux.'
     &  , wrtavg(indxSwflx),'write Swflx [cm/day]',
     &                      'Surface freshwater flux (E-P)'
     &  , wrtavg(indxShflx_rsw),'write Shflx_rsw [W/m2]',
     &                      'Short-wave surface radiation.'
      elseif (keyword(1:kwlen).eq.'gls_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxTke),  wrtavg(indxGls)
     &                    ,  wrtavg(indxLsc)
        if ( wrtavg(indxAkk) .or. wrtavg(indxAkp) .or. wrtavg(indxTke)
     &                       .or. wrtavg(indxGls) .or. wrtavg(indxLsc)
     &     ) wrtavg(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in average file: (T/F)'
     &   , wrtavg(indxTke), 'write TKE ', 'turbulent kinetic energy.  '
     &   , wrtavg(indxGls), 'write GLS ', 'generic length scale.'
     &   , wrtavg(indxLsc), 'write Lscale ',
     &                                  'vertical mixing length scale.'
      elseif (keyword(1:kwlen).eq.'surf_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  (wrtsurf(itrc), itrc=1,1)
         do itrc=1,1
           if (wrtsurf(itrc)) wrtsurf(3)=.true.
           if (mynode.eq.0) write(stdout, '(6x,L1,2x,A,2x,A,I2)')
     &         wrtsurf(itrc),
     &         'write surface outputs ',
     &         ' Momentum of index', itrc
        enddo
      elseif (keyword(1:kwlen).eq.'surf_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtsurf_avg(itrc),itrc=1,1)
        do itrc=1,1
          if (wrtsurf_avg(itrc)) wrtsurf_avg(3)=.true.
          if (mynode.eq.0) write(stdout, '(6x,L1,2x,A,2x,A,I2)')
     &        wrtsurf_avg(itrc),
     &        'write averaged surface outputs ',
     &        ' Momentum of index', itrc
       enddo
      elseif (keyword(1:kwlen).eq.'rho0') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rho0
        if (mynode.eq.0) write(stdout,'(F10.4,2x,A,1x,A)')
     &        rho0, 'rho0     Boussinesq approximation',
     &                           'mean density, kg/m3.'
      elseif (keyword(1:kwlen).eq.'lateral_visc') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) visc2, visc4
        if (mynode.eq.0) write(stdout,9) visc2
   9    format(1pe10.3,2x,'visc2    Horizontal Laplacian ',
     &       'mixing coefficient [m2/s]',/,32x,'for momentum.')
      elseif (keyword(1:kwlen).eq.'bottom_drag') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rdrg, rdrg2, Zobt, Cdb_min, Cdb_max
        if (mynode.eq.0) write(stdout,'(5(1pe10.3,2x,A/))')
     &     rdrg, 'rdrg     Linear bottom drag coefficient (m/si).',
     &    rdrg2, 'rdrg2    Quadratic bottom drag coefficient.',
     &     Zobt, 'Zobt     Bottom roughness for logarithmic law (m).',
     &  Cdb_min, 'Cdb_min  Minimum bottom drag coefficient.',
     &  Cdb_max, 'Cdb_max  Maximum bottom drag coefficient.'
      elseif (keyword(1:kwlen).eq.'gamma2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) gamma2
        if (mynode.eq.0) write(stdout,'(f10.2,2x,A,1x,A)')
     &     gamma2, 'gamma2   Slipperiness parameter:',
     &                     'free-slip +1, or no-slip -1.'
      elseif (keyword(1:kwlen).eq.'tracer_diff2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (tnu2(itrc),itrc=1,NT)
        do itrc=1,NT
          if (mynode.eq.0) write(stdout,7) tnu2(itrc), itrc, itrc
   7      format(1pe10.3,'  tnu2(',i2,')  Horizontal Laplacian '
     &     ,'mixing coefficient (m2/s)',/,32x,'for tracer ',i2,'.')
        enddo
      elseif (keyword(1:kwlen).eq.'tracer_diff4') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (tnu4(itrc),itrc=1,NT)
        do itrc=1,NT
          if (mynode.eq.0) write(stdout,8) tnu4(itrc), itrc, itrc
   8      format(1pe10.3,'  tnu4(',i2,')  Horizontal biharmonic'
     &    ,' mixing coefficient [m4/s]',/,32x,'for tracer ',i2,'.')
        enddo
      elseif (keyword(1:kwlen).eq.'sponge') then
         call cancel_kwd (keyword(1:kwlen), ierr)	
	     if (mynode.eq.0) write(stdout,'(/,1x,A,/,25x,A/)')
     &   'SPONGE_GRID is defined: x_sponge parameter in sponge/nudging',
     &   'layer is set generically in set_nudgcof.F routine'
      elseif (keyword(1:kwlen).eq.'nudg_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) tauT_in,tauT_out,tauM_in,tauM_out
          tauT_in =1.D0/(tauT_in *86400.D0)
          tauT_out=1.D0/(tauT_out*86400.D0)
          tauM_in =1.D0/(tauM_in *86400.D0)
          tauM_out=1.D0/(tauM_out*86400.D0)
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauT_in,'tauT_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &       tauT_out,'tauT_out Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauM_in,'tauM_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A/)')
     &       tauM_out,'tauM_out Nudging coefficients [sec^-1]'
      else
        if (mynode.eq.0) write(stdout,'(/3(1x,A)/)')
     &                  'WARNING: Unrecognized keyword:',
     &                   keyword(1:kwlen),' --> DISREGARDED.'
      endif
      if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
  95  write(stdout,'(/1x,4A/)') 'READ_INP ERROR while reading block',
     &                    ' with keyword ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  write(stdout,'(/1x,4A/)') 'READ_INP ERROR: Cannot find input ',
     &                                'file ''', fname(1:lstr), '''.'
      ierr=ierr+1
  99  close (input)
      if (ierr.eq.0) then
        call check_kwds (ierr)
        call check_srcs
        call check_switches1 (ierr)
        call check_switches2 (ierr)
      endif
      if (ierr.ne.0) then
        write(stdout,'(/1x,2A,I3,1x,A/)') 'READ_INP ERROR: ',
     & 'A total of', ierr, 'configuration errors discovered.'
        return
      endif
      return
      end
