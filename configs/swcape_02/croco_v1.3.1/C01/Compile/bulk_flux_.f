      subroutine bulk_flux (tile)
      implicit none
      integer*4 tile, trd, omp_get_thread_num
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
      call bulk_flux_tile (Istr,Iend,Jstr,Jend,
     &                               A2d(1,1,trd),A2d(1,2,trd))
      return
      end
      subroutine bulk_flux_tile (Istr,Iend,Jstr,Jend, aer,cer)
      implicit none
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
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      INTEGER*4              :: i, j, m, iter
      INTEGER*4              :: Istr, Iend, Jstr, Jend
      INTEGER*4              :: imin, imax, jmin, jmax
      INTEGER*4, PARAMETER   :: IterFl = 3
      INTEGER*4, PARAMETER :: mb = 2
      REAL           :: cff  , cff2
      REAL           :: rho0i, cpi  , RH   , Hlv     , wspd0
      REAL           :: TseaC, TseaK, Qsea , spec_hum, Tvstar
      REAL           :: TairC, TairK, Q    , rhoAir  , Bf
      REAL           :: delW(2) , delT , delQ , ddelW(2) , ddelT , ddelQ
      REAL           :: Cd   , Ch   , Ce   , qsat    , iexns , iexna
      REAL           :: Wstar(2), Tstar, Qstar, Wstar0  , Tstar0, Qstar0
      REAL           :: psi_u(2), ZoLu(2) , psi_t, ZoLt    , patm
      REAL           :: hfsen, hflat, hflw , upvel   , evap
      REAL           :: wspd0_cfb
      REAL, PARAMETER ::  blk_ZW =  10.0D0
      REAL, PARAMETER ::  blk_ZT =   2.0D0
      REAL, PARAMETER ::  blk_ZToZW = blk_ZT/blk_ZW
      REAL            ::  logus10,logts10
      REAL,PARAMETER :: Log10oLogZw = LOG(  10.0D0*10000.D0)/
     &                                LOG(blk_ZW*10000.D0)
      REAL           :: bulk_psiu_coare, bulk_psit_coare, air_visc
      REAL           :: iZoW , iZoT  , CC, Ch10, Ribcu , charn
      REAL           :: iZo10, iZoT10, Ri, Rr  , VisAir
      REAL           ::  aer(Istr-2:Iend+2,Jstr-2:Jend+2)
      REAL           ::  cer(Istr-2:Iend+2,Jstr-2:Jend+2)
      REAL           ::  stau(Istr-2:Iend+2,Jstr-2:Jend+2)
      if (.not.WEST_INTER) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (.not.EAST_INTER) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
      if (.not.SOUTH_INTER) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (.not.NORTH_INTER) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
      rho0i=1.0D0/rho0
      cpi=1.0D0/cp
      iexns  = (psurf*ip00)**(-rdocpd)
      DO j=jmin,jmax
        DO i=imin,imax
          wspd0     = wspd(i,j)
          wspd0     = MAX ( wspd0 , 0.1D0 * MIN(10.D0, blk_ZW) )
          wspd0_cfb = wspd_cfb(i,j)
          wspd0_cfb = MAX ( wspd0_cfb , 0.1D0 * MIN(10.D0, blk_ZW) )
          TairC = tair(i,j)
          TairK = TairC + CtoK
          TseaC = t(i,j,N,nrhs,itemp)
          TseaK = TseaC + CtoK
          RH     = rhum(i,j)
          Q      = spec_hum (RH,psurf,TairC)
          CALL exner_patm_from_tairabs(iexna,patm,Q,TairK,blk_ZT,psurf)
          rhoAir = patm*(1.D0+Q) / ( blk_Rgas*TairK*(1.D0+MvoMa*Q) )
          Qsea   = qsat(TseaK,psurf,0.98D0)
          delW(1) = SQRT(wspd0*wspd0+0.25D0)
          delW(2) = SQRT(wspd0_cfb*wspd0_cfb+0.25D0)
          delQ = Q-Qsea
          cff  = CtoK*(iexna-iexns)
          delT = TairC*iexna - TseaC*iexns + cff
          Wstar(1:mb)  = 0.035D0*delW(1:mb)*Log10oLogZw
          VisAir = air_visc ( TairC )
          charn  = 0.011D0
          Ch10   = 0.00115D0
          Ribcu  = -    blk_ZW / ( blk_Zabl * 0.004D0 * blk_beta**3 )
          DO m=1,mb
            iZo10  = g*Wstar(m) /
     &               (charn*Wstar(m)*Wstar(m)*Wstar(m)+0.11D0*g*VisAir)
            iZoT10 = 0.1D0 * exp(vonKar*vonKar /
     &                     ( Ch10*LOG( 10.0D0*iZo10 ) ) )
            CC     = LOG( blk_ZW*iZo10 )*LOG( blk_ZW*iZo10 )/
     &               LOG( blk_ZT*iZoT10 )
            Ri     =  g * blk_ZW * ( delT+cpvir*TairK*delQ )/
     &                             ( TairK*delW(m)*delW(m) )
            IF ( Ri < 0.0D0 ) THEN
              ZoLu(m)=CC*Ri/(1.0D0+Ri/Ribcu)
            ELSE
              ZoLu(m)=CC*Ri/(1.0D0+3.0D0*Ri/CC)
            ENDIF
            psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
            logus10   = LOG(blk_ZW* iZo10)
            Wstar(m)  = delW(m)*vonKar/(logus10-psi_u(m))
          ENDDO
          ZoLt      = ZoLu(mb)*blk_ZToZW
          psi_t     = bulk_psit_coare ( ZoLt    )
          logts10   = LOG(blk_ZT*iZoT10)
          cff       =      vonKar/(logts10-psi_t)
          Tstar     = delT*cff
          Qstar     = delQ*cff
          charn = 0.011D0
          IF     ( delW(1) > 18.0D0 ) then
            charn = 0.018D0
          ELSEIF ( delW(1) > 10.0D0 ) then
            charn = 0.011D0+0.125D0*(0.018D0-0.011D0)*(delW(1)-10.D0)
          ENDIF
          DO Iter=1,IterFl
            DO m=1,mb
              iZoW    = g*Wstar(m) / ( charn*Wstar(m)*Wstar(m)*Wstar(m)
     &                                 +0.11D0*g*VisAir )
              Rr      = Wstar(m)/(iZow*VisAir)
              iZoT    = MAX(8695.65D0,18181.8D0*(Rr**0.6D0))
              ZoLu(m)    = vonKar*g*blk_ZW*
     &             (Tstar*(1.0D0+cpvir*Q)+cpvir*TairK*Qstar)/
     &             (TairK*Wstar(m)*Wstar(m)*(1.0D0+cpvir*Q)+eps)
              psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
              logus10 = LOG(blk_ZW*iZoW)
              Wstar(m)   = delW(m)*vonKar/(logus10-psi_u(m))
            ENDDO
            ZoLt      = ZoLu(mb)*blk_ZToZW
            psi_t     = bulk_psit_coare ( ZoLt )
            logts10   = LOG(blk_ZT*iZoT)
            cff       = vonKar/(logts10-psi_t)
            Tstar     = delT*cff
            Qstar     = delQ*cff
              Bf=-g/TairK*Wstar(1)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0D0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2D0
              endif
              delW(1)  = SQRT(wspd0*wspd0+cff*cff)
              Bf=-g/TairK*Wstar(2)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0D0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2D0
              endif
              delW(2)  = SQRT(wspd0_cfb*wspd0_cfb+cff*cff)
            ddelW = delW
          ENDDO
          hfsen = - blk_Cpa*rhoAir*Wstar(mb)*Tstar
          Hlv   = (2.5008D0 - 0.0023719D0*TseaC)*1.0D+6
          hflat = - Hlv*rhoAir*Wstar(mb)*Qstar
          hflw   = radlw(i,j)
     &           - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK
          upvel=-1.61D0*Wstar(mb)*Qstar-(1.0D0+1.61D0*Q)*Wstar(mb)
     &                                           *     Tstar/TairK
          hflat=hflat+rhoAir*Hlv*upvel*Q
          hflat=-hflat*rho0i*cpi
          hfsen=-hfsen*rho0i*cpi
          stflx(i,j,itemp)=srflx(i,j)+hflw+hflat+hfsen
          evap=-cp*hflat/Hlv
          stflx(i,j,isalt)=(evap-prate(i,j))*t(i,j,N,nrhs,isalt)
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
          aer(i,j)  = rhoAir*delW(1)
          Cd        = (Wstar(1)/ddelW(1))**2
          cer(i,j)  = Cd
          stau(i,j) = cfb_slope * delW(1) + cfb_offset
          shflx_rsw(i,j)=srflx(i,j)
          shflx_lat(i,j)=hflat
          shflx_sen(i,j)=hfsen
          shflx_rlw(i,j)=hflw
        ENDDO
      ENDDO
      do j=jmin,jmax
        do i=imin+1,imax
          cff =0.5D0*(cer(i-1,j)+cer(i,j))
          cff2=0.5D0*(aer(i-1,j)+aer(i,j))
          sustr(i,j)= ( cff*cff2*uwnd(i,j)
     &     + 0.5D0*(stau(i,j)+stau(i-1,j))*u(i,j,N,nrhs)
     &                             )*rho0i
          sustr(i,j)=sustr(i,j)*umask(i,j)
        enddo
      enddo
      do j=jmin+1,jmax
        do i=imin,imax
          cff =0.5D0*(cer(i,j-1)+cer(i,j))
          cff2=0.5D0*(aer(i,j-1)+aer(i,j))
          svstr(i,j)=( cff*cff2*vwnd(i,j)
     &   + 0.5D0*(stau(i,j)+stau(i,j-1))*v(i,j,N,nrhs)
     &                                        )*rho0i
          svstr(i,j)=svstr(i,j)*vmask(i,j)
        enddo
      enddo
      return
      end
      FUNCTION spec_hum (RH,psfc,TairC)
      IMPLICIT NONE
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      REAL       ::  RH   , spec_hum , cff
      REAL       ::  psfc, TairC
      cff=(1.0007D0+3.46D-6*0.01D0*psfc)*6.1121D0*
     &        exp(17.502D0*TairC/(240.97D0+TairC))
      IF (RH.lt.2.0D0) then
         cff=cff*RH
         spec_hum=MvoMa*(cff/(psfc*0.01D0-0.378D0*cff))
      ELSE
         spec_hum=0.001D0*RH
      ENDIF
      END FUNCTION spec_hum
      SUBROUTINE exner_patm_from_tairabs (iexn,pair,q,tairabs,z,psfc)
      IMPLICIT NONE
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      REAL,INTENT(  out)   :: iexn,pair
      REAL,INTENT(in   )   :: q, tairabs, z, psfc
      REAL                 :: xm,q_sat,qsat
      REAL, PARAMETER      ::  g    = 9.80665D0
      INTEGER*4              :: iter
      INTEGER*4, PARAMETER   :: Niter = 3
      pair = psfc
      DO Iter = 1, Niter
        q_sat = qsat(tairabs, pair, 1.D0)
        xm    =  mm_dryair + (q/q_sat) * ( mm_water - mm_dryair )
        pair  = psfc * EXP( -g * xm * z / ( r_gas * tairabs ) )
      ENDDO
      iexn =  (pair*ip00)**(-rdocpd)
      return
      END SUBROUTINE exner_patm_from_tairabs
      FUNCTION qsat (TairK, patm, coeff)
      IMPLICIT NONE
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      REAL                 ::  qsat
      REAL                 ::  TairK, patm, coeff
      REAL                 ::  psat
      REAL, PARAMETER      ::  alpw    = 60.2227554D0
      REAL, PARAMETER      ::  betaw   = 6822.40088D0
      REAL, PARAMETER      ::  gamw    = 5.13926744D0
      REAL, PARAMETER      ::  alpi    = 32.62117980819471D0
      REAL, PARAMETER      ::  betai   = 6295.421338904806D0
      REAL, PARAMETER      ::  gami    = 0.5631331575423155D0
      IF (TairK .LE. CtoK) then
        psat = EXP( alpi - betai/TairK - gami*LOG(TairK) )
      ELSE
        psat = EXP( alpw - betaw/TairK - gamw*LOG(TairK) )
      ENDIF
      psat = coeff * psat
      qsat = (MvoMa*psat)/(patm+(MvoMa-1.0D0)*psat)
      return
      END FUNCTION qsat
      FUNCTION air_visc(TairC)
      REAL                 :: air_visc,cff
      REAL, PARAMETER      :: c0 = 1.326D-5
      REAL, PARAMETER      :: c1 = 6.542D-3
      REAL, PARAMETER      :: c2 = 8.301D-6
      REAL, PARAMETER      :: c3 = 4.84D-9
      cff      = TairC*TairC
      air_visc = c0*(1.D0+c1*TairC+c2*cff-c3*cff*TairC)
      return
      END FUNCTION air_visc
      FUNCTION bulk_psiu_coare (ZoL)
      IMPLICIT NONE
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      REAL                 ::  bulk_psiu_coare
      REAL                 ::  ZoL
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
      IF (ZoL <= 0.0D0) then
        chik = (1.0D0-15.0D0*ZoL)**0.25D0
        psik = 2.0D0*LOG(0.5D0*(1.0D0+chik))+LOG(0.5D0*(1.0D0+chik**2))
     &                                  -2.0D0*ATAN(chik)+pis2
        chic = (1.0D0-10.15D0*ZoL)**r3
        psic  = 1.5D0*LOG(r3*(chic**2+chic+1.0D0))
     &          - sqr3*ATAN((2.0D0*chic+1.0D0)/sqr3)+2.D0*pis2osqr3
        bulk_psiu_coare=psic+(psik-psic)/(1.0D0+ZoL**2)
      ELSE
        chic=-MIN(50.0D0,0.35D0*ZoL)
        bulk_psiu_coare=-((1.0D0+ZoL)+0.6667D0*(ZoL-14.28D0)*EXP(chic)+
     &                                                          8.525D0)
      ENDIF
      return
      END FUNCTION bulk_psiu_coare
      FUNCTION bulk_psit_coare (ZoL)
      IMPLICIT NONE
      REAL, PARAMETER ::  CtoK      =  273.16D0
      REAL, PARAMETER ::  blk_Rgas  =  287.0596736665907D0
      REAL, PARAMETER ::  blk_Rvap  =  461.5249933083879D0
      REAL, PARAMETER ::  blk_Cpa   = 1004.708857833067D0
      REAL, PARAMETER ::  ip00      =    1.D-5
      REAL, PARAMETER ::  p00       =    1.D+5
      REAL, PARAMETER ::  rdocpd    = blk_Rgas/blk_Cpa
      REAL, PARAMETER ::  cpdord    = blk_Cpa/blk_Rgas
      REAL, PARAMETER ::  r_gas     = 8.314510D0
      REAL, PARAMETER ::  mm_dryair = 28.9644D-3
      REAL, PARAMETER ::  mm_water  = 18.0153D-3
      REAL, PARAMETER ::  cpvir     = blk_Rvap/blk_Rgas - 1.D0
      REAL, PARAMETER ::  MvoMa     = mm_water/mm_dryair
      REAL, PARAMETER ::  eps       = 1.D-8
      REAL, PARAMETER ::  r3        = 1.D0/3.D0
      REAL, PARAMETER ::  pis2      = 2.D0*ATAN(1.D0)
      REAL, PARAMETER ::  sqr3      = SQRT(3.D0)
      REAL, PARAMETER ::  pis2osqr3 = pis2/sqr3
      REAL, PARAMETER ::  LMOmin    = -200.D0
      REAL, PARAMETER ::  LMOmax    = 0.25D0
      REAL, PARAMETER ::  blk_beta  = 1.2D0
      REAL, PARAMETER ::  blk_Zabl  = 600.D0
      REAL, PARAMETER ::  dWstar0   = 1.D-6
      REAL, PARAMETER ::  dTstar0   = 1.D-6
      REAL, PARAMETER ::  dQstar0   = 1.D-9
      REAL, PARAMETER ::  emiss_lw  =  0.985D0
      REAL, PARAMETER ::  SigmaSB   =  5.6697D-8
      REAL            :: psurf
      PARAMETER (psurf=100000.0D0)
      real swparam
      parameter (swparam=0.3D0)
      real cfb_slope, cfb_offset
      parameter (cfb_slope=-0.0029D0)
      parameter (cfb_offset=0.008D0)
      REAL                 ::  bulk_psit_coare
      REAL                 ::  ZoL
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
      IF (ZoL < 0.0D0) THEN
        chik  = (1.0D0-15.0D0*ZoL)**0.25D0
        psik  = 2.0D0*LOG(0.5D0*(1.0D0+chik**2))
        chic  = (1.0D0-34.15D0*ZoL)**r3
        psic  = 1.5D0*LOG((chic**2+chic+1.0D0)*r3)
     &               -sqr3*ATAN((2.0D0*chic+1.0D0)/sqr3)
     &               +2.D0*pis2osqr3
        bulk_psit_coare = psic+(psik-psic)/(1.0D0+ZoL**2)
      ELSE
        chic=-MIN(50.0D0,0.35D0*ZoL)
        bulk_psit_coare = -((1.0D0+2.0D0*ZoL/3.0D0)**1.5D0+
     &            0.6667D0*(ZoL-14.28D0)*EXP(chic)+8.525D0)
      ENDIF
      return
      END FUNCTION bulk_psit_coare
