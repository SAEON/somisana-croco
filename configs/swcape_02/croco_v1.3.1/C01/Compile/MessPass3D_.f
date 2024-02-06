      subroutine MessPass3D_tile (Istr,Iend,Jstr,Jend, A, nmax)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=2)
      integer*4 nmax
      real A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmax)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr
      integer*4 iter, mdii, mdjj
      integer*4 sub_X,size_X, sub_E,size_E, size_Z
      parameter (size_Z=Npts*Npts*(N+1),
     &     sub_X=(Lm+NSUB_X-1)/NSUB_X, size_X=(N+1)*Npts*(sub_X+2*Npts),
     &     sub_E=(Mm+NSUB_E-1)/NSUB_E, size_E=(N+1)*Npts*(sub_E+2*Npts))
      real buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &     buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &    jbuf_sndW(size_E),                      jbuf_sndE(size_E),
     &    jbuf_revW(size_E),                      jbuf_revE(size_E),
     &     buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &     buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (ii.eq.0 .and. Istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
      ishft=imax-imin+1
      if (jj.eq.0 .and. Jstr.eq.1) then
        jmin=Jstr-1
      else
        jmin=Jstr
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts*nmax
      isize=Npts*ishft*nmax
      jsize=Npts*jshft*nmax
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if ((WEST_INTER2)) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndW(k+nmax*(j-jmin+(ipts-1)*jshft))=A(ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER2) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndE(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                                       A(Lmmpi-Npts+ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER2) then
            ibuf_snds = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndS(k+nmax*(i-imin+(jpts-1)*ishft))=A(i,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER2) then
            ibuf_sndn = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndN(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                                         A(i,Mmmpi-Npts+jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (CORNER_SW) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd1(k+nmax*(ipts-1+Npts*(jpts-1)))=A(ipts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (CORNER_NE) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd2(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                            A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (CORNER_SE) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd3(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(Lmmpi+ipts-Npts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                8, MPI_COMM_WORLD, req(7),ierr)
            call MPI_Send  (buf_snd3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                7, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (CORNER_NW) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd4(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(ipts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                7, MPI_COMM_WORLD, req(8),ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                8, MPI_COMM_WORLD,        ierr)
          endif
        endif
      enddo
      if (WEST_INTER2) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
             A(ipts-Npts,j,k)=jbuf_revW(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (WEST_INTER .and. .not. WEST_INTER2) then
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(ipts-Npts,j,k)=A(ipts,j,k)
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER2) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=
     &                         jbuf_revE(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER .and. .not. EAST_INTER2) then
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=A(Lmmpi+ipts-Npts,j,k)
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER2) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=ibuf_revS(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. .not. SOUTH_INTER2) then
       do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=A(i,jpts,k)
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER2) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=
     &                         ibuf_revN(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER .and. .not. NORTH_INTER2) then
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=A(i,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
      endif
      if (CORNER_SW) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           buf_rev1(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SW  .and.
     &   SOUTH_INTER .and.  WEST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           A(ipts,jpts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_NE) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &                           buf_rev2(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NE  .and.
     &   NORTH_INTER .and.  EAST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &              A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_SE) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,jpts-Npts,k)=
     &                           buf_rev3(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SE .and.
     &   SOUTH_INTER .and.  EAST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
               A(Lmmpi+ipts,jpts-Npts,k)=
     &                           A(Lmmpi+ipts-Npts,jpts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_NW) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=
     &                           buf_rev4(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NW  .and.
     &   NORTH_INTER .and.  WEST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=A(ipts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif
      return
      end
      subroutine MessPass3D_3pts_tile (Istr,Iend,Jstr,Jend, A, nmax)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=3)
      integer*4 nmax
      real A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmax)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr
      integer*4 iter, mdii, mdjj
      integer*4 sub_X,size_X, sub_E,size_E, size_Z
      parameter (size_Z=Npts*Npts*(N+1),
     &     sub_X=(Lm+NSUB_X-1)/NSUB_X, size_X=(N+1)*Npts*(sub_X+2*Npts),
     &     sub_E=(Mm+NSUB_E-1)/NSUB_E, size_E=(N+1)*Npts*(sub_E+2*Npts))
      real buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &     buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &    jbuf_sndW(size_E),                      jbuf_sndE(size_E),
     &    jbuf_revW(size_E),                      jbuf_revE(size_E),
     &     buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &     buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (ii.eq.0 .and. Istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
      ishft=imax-imin+1
      if (jj.eq.0 .and. Jstr.eq.1) then
        jmin=Jstr-1
      else
        jmin=Jstr
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts*nmax
      isize=Npts*ishft*nmax
      jsize=Npts*jshft*nmax
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if ((WEST_INTER2)) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndW(k+nmax*(j-jmin+(ipts-1)*jshft))=A(ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER2) then
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndE(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                                       A(Lmmpi-Npts+ipts,j,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER2) then
            ibuf_snds = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndS(k+nmax*(i-imin+(jpts-1)*ishft))=A(i,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER2) then
            ibuf_sndn = 0.D0
            do k=1,nmax
              do i=imin,imax
                do jpts=1,Npts
                  ibuf_sndN(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                                         A(i,Mmmpi-Npts+jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (CORNER_SW) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd1(k+nmax*(ipts-1+Npts*(jpts-1)))=A(ipts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (CORNER_NE) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd2(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                            A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2, ksize, MPI_DOUBLE_PRECISION, p_NE,
     &                                  6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (CORNER_SE) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd3(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(Lmmpi+ipts-Npts,jpts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                8, MPI_COMM_WORLD, req(7),ierr)
            call MPI_Send  (buf_snd3, ksize, MPI_DOUBLE_PRECISION, p_SE,
     &                                7, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (CORNER_NW) then
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd4(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(ipts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                7, MPI_COMM_WORLD, req(8),ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                                8, MPI_COMM_WORLD,        ierr)
          endif
        endif
      enddo
      if (WEST_INTER2) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
             A(ipts-Npts,j,k)=jbuf_revW(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (WEST_INTER .and. .not. WEST_INTER2) then
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(ipts-Npts,j,k)=A(ipts,j,k)
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER2) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=
     &                         jbuf_revE(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
      endif
      if (EAST_INTER .and. .not. EAST_INTER2) then
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=A(Lmmpi+ipts-Npts,j,k)
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER2) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=ibuf_revS(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. .not. SOUTH_INTER2) then
       do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=A(i,jpts,k)
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER2) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=
     &                         ibuf_revN(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
      endif
      if (NORTH_INTER .and. .not. NORTH_INTER2) then
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=A(i,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
      endif
      if (CORNER_SW) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           buf_rev1(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SW  .and.
     &   SOUTH_INTER .and.  WEST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           A(ipts,jpts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_NE) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &                           buf_rev2(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NE  .and.
     &   NORTH_INTER .and.  EAST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &              A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_SE) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,jpts-Npts,k)=
     &                           buf_rev3(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SE .and.
     &   SOUTH_INTER .and.  EAST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
               A(Lmmpi+ipts,jpts-Npts,k)=
     &                           A(Lmmpi+ipts-Npts,jpts,k)
            enddo
          enddo
        enddo
       endif
      if (CORNER_NW) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=
     &                           buf_rev4(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NW  .and.
     &   NORTH_INTER .and.  WEST_INTER ) then
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=A(ipts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif
      return
      end
