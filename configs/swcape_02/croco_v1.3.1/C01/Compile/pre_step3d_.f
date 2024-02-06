      subroutine pre_step3d (tile)
      implicit none
      integer*4 tile,  trd,omp_get_thread_num
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
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
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
      call pre_step3d_tile (Istr,Iend,Jstr,Jend,
     &                A3d(1,1,trd), A3d(1,2,trd), A3d(1,3,trd),
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd)
     &                                          , A3d(1,4,trd)
     &                                                        )
      return
      end
      subroutine pre_step3d_tile (Istr,Iend,Jstr,Jend, ru,rv,rw,
     &                                                 FC,CF,DC,
     &                                               FX,FE,WORK
     &                                                 ,Hz_half
     &                                                         )
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
      integer*4 Istr,Iend,Jstr,Jend, itrc, i,j,k, indx
     &       ,imin,imax,jmin,jmax,nadv,iAkt
      real   ru(Istr-2:Iend+2,Jstr-2:Jend+2,N),    cff,
     &       rv(Istr-2:Iend+2,Jstr-2:Jend+2,N),    cff1,
     &       rw(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &       FC(Istr-2:Iend+2,0:N),  cff2,
     &       CF(Istr-2:Iend+2,0:N),
     &       DC(Istr-2:Iend+2,0:N),  gamma,
     &       FX(Istr-2:Iend+2,Jstr-2:Jend+2),      epsil,
     &       FE(Istr-2:Iend+2,Jstr-2:Jend+2),        cdt,
     &     WORK(Istr-2:Iend+2,Jstr-2:Jend+2)
      real Hz_half(Istr-2:Iend+2,Jstr-2:Jend+2,N)
      parameter (gamma=1.D0/6.D0, epsil=1.D-16)
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
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real rho(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real rhoA(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rhoS(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rvfrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rufrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real rvfrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DU_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DV_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DU_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DV_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
      real zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
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
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (.not.WEST_INTER) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (.not.EAST_INTER) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (.not.SOUTH_INTER) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (.not.NORTH_INTER) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      indx=3-nstp
      nadv=  nstp
      if (iic.eq.ntstart) then
        cff=0.5D0*dt
        cff1=1.D0
        cff2=0.D0
      else
        cff=(1.D0-gamma)*dt
        cff1=0.5D0+gamma
        cff2=0.5D0-gamma
      endif
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            Hz_half(i,j,k)=cff1*Hz(i,j,k)+cff2*Hz_bak(i,j,k)
     &        -cff*pm(i,j)*pn(i,j)*( Huon(i+1,j,k)-Huon(i,j,k)
     &                              +Hvom(i,j+1,k)-Hvom(i,j,k)
     &                                  +We(i,j,k)-We(i,j,k-1)
     &                                                       )
          enddo
        enddo
      enddo
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
      do itrc=1,NT
        do k=1,N
          if (WEST_INTER) then
            imin=Istr-1
          else
            imin=max(Istr-1,1)
          endif
          if (EAST_INTER) then
            imax=Iend+2
          else
            imax=min(Iend+2,Lmmpi+1)
          endif
          if (SOUTH_INTER) then
            jmin=Jstr-1
          else
            jmin=max(Jstr-1,1)
          endif
          if (NORTH_INTER) then
            jmax=Jend+2
          else
            jmax=min(Jend+2,Mmmpi+1)
          endif
          do j=Jstr,Jend
            do i=imin,imax
              FX(i,j)=(t(i,j,k,nadv,itrc)-t(i-1,j,k,nadv,itrc))
     &                                               *umask(i,j)
            enddo
          enddo
          if (.not.WEST_INTER) then
            do j=Jstr,Jend
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (.not.EAST_INTER) then
            do j=Jstr,Jend
              FX(Lmmpi+2,j)=FX(Lmmpi+1,j)
            enddo
          endif
          do j=Jstr,Jend
            do i=Istr-1,Iend+1
              WORK(i,j)=0.5D0*(FX(i+1,j)+FX(i,j))
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FX(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                     -0.333333333333D0*(WORK(i,j)-WORK(i-1,j))
     &                                                )*Huon(i,j,k)
            enddo
          enddo
          do j=jmin,jmax
            do i=Istr,Iend
              FE(i,j)=(t(i,j,k,nadv,itrc)-t(i,j-1,k,nadv,itrc))
     &                                               *vmask(i,j)
            enddo
          enddo
          if (.not.SOUTH_INTER) then
            do i=Istr,Iend
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (.not.NORTH_INTER) then
            do i=Istr,Iend
              FE(i,Mmmpi+2)=FE(i,Mmmpi+1)
            enddo
          endif
          do j=Jstr-1,Jend+1
            do i=Istr,Iend
              WORK(i,j)=0.5D0*(FE(i,j+1)+FE(i,j))
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              FE(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                     -0.333333333333D0*(WORK(i,j)-WORK(i,j-1))
     &                                               )*Hvom(i,j,k)
            enddo
          enddo
          if (iic.eq.ntstart) then
            cff=0.5D0*dt
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                      -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                            +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          else
            cff=(1.D0-gamma)*dt
            cff1=0.5D0+gamma
            cff2=0.5D0-gamma
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=cff1*Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                            +cff2*Hz_bak(i,j,k)*t(i,j,k,indx,itrc)
     &                         -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                               +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          endif
        enddo
      enddo
      if (iic.eq.ntstart) then
            cdt=0.5D0*dt
      else
            cdt=(1.D0-gamma)*dt
      endif
      do j=Jstr,Jend
        do i=Istr,Iend
           do k=1,N
             DC(i,k)=1.D0/Hz_half(i,j,k)
           enddo
           DC(i,0)=cdt*pn(i,j)*pm(i,j)
        enddo
        do itrc=1,NT
          do i=Istr,Iend
            FC(i,0)=0.D0
            CF(i,0)=0.D0
          enddo
          do k=1,N-1,+1
            do i=Istr,Iend
              cff    = 1.D0
     &                    /(2.D0*Hz(i,j,k+1)+Hz(i,j,k)*(2.D0-FC(i,k-1)))
              FC(i,k)= cff*Hz(i,j,k+1)
              CF(i,k)= cff*( 6.D0*( t(i,j,k+1,nadv,itrc)
     &                           -t(i,j,k  ,nadv,itrc) 
     &                                             )-Hz(i,j,k)*CF(i,k-1)
     &                                                                  
     &                                                                 )
            enddo
          enddo
          do i=Istr,Iend
            CF(i,N)=0.D0
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            enddo
          enddo
          cff=1.D0/3.D0
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=We(i,j,k)*( t(i,j,k,nadv,itrc)+cff*Hz(i,j,k)
     &                                  *(CF(i,k)+0.5D0*CF(i,k-1))
     &                                                            )
            enddo
          enddo
          do i=Istr,Iend
            FC(i,N)=0.D0
            FC(i,0)=0.D0
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo
          do k=1,N
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=DC(i,k)*( t(i,j,k,nnew,itrc)
     &               -DC(i,0)*(FC(i,k)-FC(i,k-1)))
            enddo
          enddo
        enddo
        do i=IstrU,Iend
          DC(i,0)=pm_u(i,j)*pn_u(i,j)
        enddo
        if (iic.eq.ntstart) then
          do k=1,N
            do i=IstrU,Iend
              cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i-1,j,k))
              u(i,j,k,nnew)=(u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                           +cdt*DC(i,0)*ru(i,j,k)
     &                      )*cff
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                         Hz(i-1,j,k))
            enddo
          enddo
        else
          cff1=0.5D0+gamma
          cff2=0.5D0-gamma
          do k=1,N
            do i=IstrU,Iend
              cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i-1,j,k))
              u(i,j,k,nnew)=( cff1*u(i,j,k,nstp)*0.5D0*(Hz(i  ,j,k)+
     &                                                Hz(i-1,j,k))
     &                       +cff2*u(i,j,k,indx)*0.5D0*(Hz_bak(i  ,j,k)+
     &                                                Hz_bak(i-1,j,k))
     &                       +cdt*DC(i,0)*ru(i,j,k)
     &                                                     )*cff
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                         Hz(i-1,j,k))
            enddo
          enddo
         endif
        if (j.ge.JstrV) then
          do i=Istr,Iend
            DC(i,0)=pm_v(i,j)*pn_v(i,j)
          enddo
          if (iic.eq.ntstart) then
            do k=1,N
              do i=Istr,Iend
                cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i,j-1,k))
                v(i,j,k,nnew)=(v(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+Hz(i,j-1,
     &                                                               k))
     &                                           +cdt*DC(i,0)*rv(i,j,k)
     &                           )*cff
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5D0*(Hz(i,j  ,k)+
     &                                           Hz(i,j-1,k))
              enddo
            enddo
          else
            cff1=0.5D0+gamma
            cff2=0.5D0-gamma
            do k=1,N
              do i=Istr,Iend
                cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i,j-1,k))
                v(i,j,k,nnew)=( cff1*v(i,j,k,nstp)*0.5D0*(Hz(i,j  ,k)+
     &                                                  Hz(i,j-1,k))
     &                         +cff2*v(i,j,k,indx)*0.5D0*(Hz_bak(i,j  
     &                                                              ,k)+
     &                                                  Hz_bak(i,j-1,k))
     &                         +cdt*DC(i,0)*rv(i,j,k)
     &                                                           )*cff
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                           Hz(i,j-1,k))
              enddo
            enddo
          endif
        endif
      enddo
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(-1,-1,1,nnew,itrc))
      enddo
      call u3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(-1,-1,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(-1,-1,1,nnew))
      do k=0,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Akv_old(i,j,k)=Akv(i,j,k)
            Akt_old(i,j,k)=Akt(i,j,k,itemp)
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,knew)=Zt_avg1(i,j)
        enddo
      enddo
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(-1,-1,knew))
      return
      end
