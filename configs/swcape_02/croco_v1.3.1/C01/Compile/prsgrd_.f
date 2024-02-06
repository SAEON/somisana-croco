      subroutine prsgrd (tile)
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
      call prsgrd_tile (Istr,Iend,Jstr,Jend,
     &                  A3d(1,1,trd), A3d(1,2,trd), A3d(1,3,trd),
     &                  A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                  A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd))
      return
      end
      subroutine prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                                       dR,dZ, FC,dZx,rx,dRx)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, imin,imax,jmin,jmax
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
      real ru(Istr-2:Iend+2,Jstr-2:Jend+2,N),   OneFifth,
     &     rv(Istr-2:Iend+2,Jstr-2:Jend+2,N),   OneTwelfth,
     &      P(Istr-2:Iend+2,Jstr-2:Jend+2,N),   epsil,dpth,cff1,cff2,
     &     dR(Istr-2:Iend+2,0:N), cff, GRho,
     &     dZ(Istr-2:Iend+2,0:N), cfr, HalfGRho,
     &     FC(Istr-2:Iend+2,Jstr-2:Jend+2),
     &    dZx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &     rx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &    dRx(Istr-2:Iend+2,Jstr-2:Jend+2)
      parameter (OneFifth=0.2D0, OneTwelfth=1.D0/12.D0, epsil=0.D0)
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
      GRho=g/rho0
      HalfGRho=0.5D0*GRho
      do j=JstrV-1,Jend
        do k=1,N-1
          do i=IstrU-1,Iend
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
            dpth=z_w(i,j,N)-0.5D0*(z_r(i,j,k+1)+z_r(i,j,k))
            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)
     &              +(qp1(i,j,k+1)-qp1(i,j,k))
     &                     *dpth*(1.D0-qp2*dpth)
          enddo
        enddo
        do i=IstrU-1,Iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1
          do i=IstrU-1,Iend
            cff=2.D0*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))
            cfr=2.D0*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.D0
            endif
            dpth=z_w(i,j,N)-z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1.D0-2.D0*qp2*dpth)
          enddo
        enddo
        do i=IstrU-1,Iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     &            +0.5D0*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,
     &                                                               N))
     &              /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=IstrU-1,Iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                    *(z_r(i,j,k+1)-z_r(i,j,k))
     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                             -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )
     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                             -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                            ))
          enddo
        enddo
      enddo
      if (.not.WEST_INTER) then
        imin=IstrU
      else
        imin=IstrU-1
      endif
      if (.not.EAST_INTER) then
        imax=Iend
      else
        imax=Iend+1
      endif
      do k=N,1,-1
        do j=Jstr,Jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
     &                              *umask(i,j)
            dpth=0.5D0*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))
            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)
     &                +(qp1(i,j,k)-qp1(i-1,j,k))
     &                     *dpth*(1.D0-qp2*dpth) )
     &                              *umask(i,j)
          enddo
        enddo
        if (.not.WEST_INTER) then
          do j=Jstr,Jend
            FC(imin-1,j)=FC(imin,j)
            rx(imin-1,j)=rx(imin,j)
          enddo
        endif
        if (.not.EAST_INTER) then
          do j=Jstr,Jend
            FC(imax+1,j)=FC(imax,j)
            rx(imax+1,j)=rx(imax,j)
          enddo
        endif
        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cff=2.D0*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.D0
            endif
            cfr=2.D0*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.D0
            endif
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.D0-2.D0*qp2*(z_w(i,j,N)-z_r(i,j,k)))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*(
     &                                             P(i-1,j,k)-P(i,j,k)
     &                           -HalfGRho*( (rho(i,j,k)+rho(i-1,j,k))
     &                                      *(z_r(i,j,k)-z_r(i-1,j,k))
     &     -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                             -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )
     &                -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                             -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                             ))
     &    )
          enddo
        enddo
        if (.not.SOUTH_INTER) then
          jmin=JstrV
        else
          jmin=JstrV-1
        endif
        if (.not.NORTH_INTER) then
          jmax=Jend
        else
          jmax=Jend+1
        endif
        do j=jmin,jmax
          do i=Istr,Iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
     &                              *vmask(i,j)
            dpth=0.5D0*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))
            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)
     &                +(qp1(i,j,k)-qp1(i,j-1,k))
     &                     *dpth*(1.D0-qp2*dpth) )
     &                              *vmask(i,j)
          enddo
        enddo
        if (.not.SOUTH_INTER) then
          do i=Istr,Iend
            FC(i,jmin-1)=FC(i,jmin)
            rx(i,jmin-1)=rx(i,jmin)
          enddo
        endif
        if (.not.NORTH_INTER) then
          do i=Istr,Iend
            FC(i,jmax+1)=FC(i,jmax)
            rx(i,jmax+1)=rx(i,jmax)
          enddo
        endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
            cff=2.D0*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.D0
            endif
            cfr=2.D0*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.D0
            endif
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.D0-2.D0*qp2*(z_w(i,j,N)-z_r(i,j,k)))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*(
     &                                             P(i,j-1,k)-P(i,j,k)
     &                           -HalfGRho*( (rho(i,j,k)+rho(i,j-1,k))
     &                                      *(z_r(i,j,k)-z_r(i,j-1,k))
     &     -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                             -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )
     &                -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                             -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                             ))
     &                                                               )
          enddo
        enddo
      enddo
      return
      end
