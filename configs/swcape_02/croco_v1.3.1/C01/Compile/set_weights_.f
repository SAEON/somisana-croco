      subroutine set_weights
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
      integer*4 i,j, iter
      real p,q,r, scale
      real*8 sum, shft, cff
      nfast=0
      do i=0,2*ndtfast
        weight(1,i)=0.D0
        weight(2,i)=0.D0
      enddo
      p=2.0D0
      q=4.0D0
      r=0.25D0
      scale=(p+1.D0)*(p+q+1.D0)/((p+2.D0)*(p+q+2.D0)*float(ndtfast))
      do iter=1,16
        nfast=0
        do i=1,2*ndtfast
          cff=scale*float(i)
          weight(1,i)=cff**p - cff**(p+q) - r*cff
          if (i.gt.ndtfast .and. weight(1,i).lt.0.D0) then
            weight(1,i)=0.D0
          else
            nfast=i
          endif
        enddo
        sum=0.D0
        shft=0.D0
        do i=1,nfast
          sum=sum+weight(1,i)
          shft=shft+weight(1,i)*float(i)
        enddo
        scale=scale * shft/(sum*float(ndtfast))
      enddo
      do iter=1,ndtfast
        sum=0.D0
        shft=0.D0
        do i=1,nfast
          sum=sum+weight(1,i)
          shft=shft+float(i)*weight(1,i)
        enddo
        shft=shft/sum
        cff=float(ndtfast)-shft
        if (cff .gt. 1.D0) then
          nfast=nfast+1
          do i=nfast,2,-1
            weight(1,i)=weight(1,i-1)
          enddo
          weight(1,1)=0.D0
        elseif (cff .gt. 0.D0) then
          sum=1.D0-cff
          do i=nfast,2,-1
            weight(1,i)=sum*weight(1,i)+cff*weight(1,i-1)
          enddo
          weight(1,1)=sum*weight(1,1)
        elseif (cff .lt. -1.D0) then
          nfast=nfast-1
          do i=1,nfast,+1
            weight(1,i)=weight(1,i+1)
          enddo
          weight(1,nfast+1)=0.D0
        elseif (cff .lt. 0.D0) then
          sum=1.D0+cff
          do i=1,nfast-1,+1
            weight(1,i)=sum*weight(1,i)-cff*weight(1,i+1)
          enddo
          weight(1,nfast)=sum*weight(1,nfast)
        endif
      enddo
      do j=1,nfast
        cff=weight(1,j)
        do i=1,j
          weight(2,i)=weight(2,i)+cff
        enddo
      enddo
      sum=0.D0
      cff=0.D0
      do i=1,nfast
        sum=sum+weight(1,i)
        cff=cff+weight(2,i)
      enddo
      sum=1.D0/sum
      cff=1.D0/cff
      do i=1,nfast
        weight(1,i)=sum*weight(1,i)
        weight(2,i)=cff*weight(2,i)
      enddo
      if (mynode.eq.0) write(stdout,'(/1x,A,I3,4x,A,I3/)')
     &   'Time splitting: ndtfast =', ndtfast, 'nfast =', nfast
      return
      end
