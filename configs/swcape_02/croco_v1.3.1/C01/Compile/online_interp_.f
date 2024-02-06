      SUBROUTINE linterp2d (LBx, UBx, LBy, UBy,
     &                      Xinp, Yinp, Finp,
     &                      Istr, Iend, Jstr, Jend,
     &                      Xout, Yout, Fout)
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
      integer*4, intent(in) :: LBx, UBx, LBy, UBy
      integer*4, intent(in) :: Istr, Iend, Jstr, Jend
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
      real(kind=8), intent(in) :: Xout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8), intent(in) :: Yout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8), intent(out) :: Fout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 i, i1, i2, j, j1, j2, ii, jj
      real(kind=8) cff, x, x1, x2, y, y1, y2
      DO j=Jstr,Jend
        DO i=Istr,Iend
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     &           (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2'
           goto 100
 10        continue
           DO jj=LBy,(UBy-1)
             if ((Yinp(jj).le.Yout(i,j)).and.
     &           (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
 20        continue
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
            x1=Xinp(i1)
            x2=Xinp(i2)
            y1=Yinp(j1)
            y2=Yinp(j2)
            x=Xout(i,j)
            y=Yout(i,j)
            cff= Finp(i1,j1)*(x2-x )*(y2-y )
     &          +Finp(i2,j1)*(x -x1)*(y2-y )
     &          +Finp(i1,j2)*(x2-x )*(y -y1)
     &          +Finp(i2,j2)*(x -x1)*(y -y1)
            Fout(i,j)=cff/((x2-x1)*(y2-y1))
          END IF
        END DO
      END DO
      RETURN
 100  continue
      print*, 'error in linterp2d'
      END SUBROUTINE linterp2d
      SUBROUTINE cinterp2d (LBx, UBx, LBy, UBy,
     &                      Xinp, Yinp, Finp,
     &                      Istr, Iend, Jstr, Jend,
     &                      Xout, Yout, Fout)
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
      integer*4, intent(in) :: LBx, UBx, LBy, UBy
      integer*4, intent(in) :: Istr, Iend, Jstr, Jend
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
      real(kind=8), intent(in) :: Xout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8), intent(in) :: Yout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8), intent(out) :: Fout(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 i, ic, iter, i1, i2, j, jc, j1, j2, ii, jj
      real(kind=8) :: a11, a12, a21, a22
      real(kind=8) :: e11, e12, e21, e22
      real(kind=8) :: cff, d1, d2, dfc, dx, dy, eta, xi, xy, yx
      real(kind=8) :: f0, fx, fxx, fxxx, fxxy, fxy, fxyy, fy, fyy, fyyy
      real(kind=8), parameter :: C01 = 1.0D0/48.0D0
      real(kind=8), parameter :: C02 = 1.0D0/32.0D0
      real(kind=8), parameter :: C03 = 0.0625D0
      real(kind=8), parameter :: C04 = 1.0D0/6.0D0
      real(kind=8), parameter :: C05 = 0.25D0
      real(kind=8), parameter :: C06 = 0.5D0
      real(kind=8), parameter :: C07 = 0.3125D0
      real(kind=8), parameter :: C08 = 0.625D0
      real(kind=8), parameter :: C09 = 1.5D0
      real(kind=8), parameter :: C10 = 13.0D0/24.0D0
      real(kind=8), parameter :: LIMTR = 3.0D0
      real(kind=8), parameter :: spv = 0.0D0
      real(kind=8), dimension(-1:2,-1:2) :: dfx, dfy, ff
      DO j=Jstr,Jend
        DO i=Istr,Iend
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     &           (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2',
     &           Istr,Iend,Jstr,Jend,i,j,Xout(i,j),Xout(i-1,j)
           goto 100
 10        continue
           DO jj=LBy,UBy-1
             if ((Yinp(jj).le.Yout(i,j)).and.
     &           (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
 20        continue
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
            xy=Xinp(i2)-Xinp(i1)-Xinp(i2)+Xinp(i1)
            yx=Yinp(j2)-Yinp(j2)-Yinp(j1)+Yinp(j1)
            dx=Xout(i,j)-0.25D0*(Xinp(i2)+Xinp(i1)+
     &                         Xinp(i2)+Xinp(i1))
            dy=Yout(i,j)-0.25D0*(Yinp(j2)+Yinp(j2)+
     &                         Yinp(j1)+Yinp(j1))
            e11=0.5D0*(Xinp(i2)-Xinp(i1)+Xinp(i2)-Xinp(i1))
            e12=0.5D0*(Xinp(i2)+Xinp(i1)-Xinp(i2)-Xinp(i1))
            e21=0.5D0*(Yinp(j2)-Yinp(j2)+Yinp(j1)-Yinp(j1))
            e22=0.5D0*(Yinp(j2)+Yinp(j2)-Yinp(j1)-Yinp(j1))
            cff=1.0D0/(e11*e22-e12*e21)
            xi=cff*(e22*dx-e12*dy)
            eta=cff*(e11*dy-e21*dx)
            DO iter=1,4
              d1=dx-e11*xi-e12*eta-xy*xi*eta
              d2=dy-e21*xi-e22*eta-yx*xi*eta
              a11=e11+xy*eta
              a12=e12+xy*xi
              a21=e21+yx*eta
              a22=e22+yx*xi
              cff=1.0D0/(a11*a22-a12*a21)
              xi =xi +cff*(a22*d1-a12*d2)
              eta=eta+cff*(a11*d2-a21*d1)
            END DO
            DO jc=-1,2
              DO ic=-1,2
                ff(ic,jc)=Finp(MAX(1,MIN(UBx,i1+ic)),
     &                         MAX(1,MIN(UBy,j1+jc)))
              END DO
            END DO
            f0=C07*(ff(1,1)+ff(1,0)+ff(0,1)+ff(0,0))-
     &         C02*(ff(2,0)+ff(2,1)+ff(1,2)+ff(0,2)+
     &              ff(-1,1)+ff(-1,0)+ff(0,-1)+ff(1,-1))
            fx=C08*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))-
     &         C01*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     &         C03*(ff(1,2)-ff(0,2)+ff(1,-1)-ff(0,-1))
            fy=C08*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))-
     &         C01*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     &         C03*(ff(2,1)-ff(2,0)+ff(-1,1)-ff(-1,0))
            fxy=ff(1,1)-ff(1,0)-ff(0,1)+ff(0,0)
            fxx=C05*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)+
     &               ff(2,0)-ff(1,0)-ff(0,0)+ff(-1,0))
            fyy=C05*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)+
     &               ff(0,2)-ff(0,1)-ff(0,0)+ff(0,-1))
            fxxx=C06*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     &           C09*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))
            fyyy=C06*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     &           C09*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))
            fxxy=C06*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)-
     &                ff(2,0)+ff(1,0)+ff(0,0)-ff(-1,0))
            fxyy=C06*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)-
     &                ff(0,2)+ff(0,1)+ff(0,0)-ff(0,-1))
            Fout(i,j)=f0+
     &                fx*xi+
     &                fy*eta+
     &                C06*fxx*xi*xi+
     &                fxy*xi*eta+
     &                C06*fyy*eta*eta+
     &                C04*fxxx*xi*xi*xi+
     &                C06*fxxy*xi*xi*eta+
     &                C04*fyyy*eta*eta*eta+
     &                C06*fxyy*xi*eta*eta
          END IF
        END DO
      END DO
      RETURN
 100  continue
      print*, 'error in cinterp2d'
      END SUBROUTINE cinterp2d
