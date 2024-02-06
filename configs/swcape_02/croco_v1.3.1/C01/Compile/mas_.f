      module parameters
      integer*4,parameter :: rsh=4,rlg=8
      real,parameter    :: pi=3.141592654D0
      real(kind=rsh),parameter    :: valmanq=999.D0
      end module parameters
      subroutine ma7(flag,seci,hn,haut,k)
      use parameters
      implicit none
      logical :: flag
      integer*4                      :: idk,l,k
      real               :: d1,haut,fk,d,hh,d2
      real,dimension(130):: hn
      real               :: seci
      fk=seci/60.0D0/11.25D0+2.0D0
      k=fk+0.5D0
      d=fk-k
      if(k >= 130)then
       idk=k-129
       k=k-idk
       d=d+idk
      endif
      if (flag) write(6,*)"ma7 seci",seci
      if (flag) write(6,*)"ma7 k",k
      if (flag) write(6,*)"ma7 hn",hn
      l=k
      if(l == 130)l=129
      hh=hn(l)+hn(l)
      d1=hn(l+1)-hn(l-1)
      d2=hn(l+1)+hn(l-1)-hh
      haut=(hh+d*(d*d2+d1))/800.0D0
      if (flag) write(6,*)"ma7 haut",haut,l,hh,hn(l+1),hn(l-1),d1,d2
      end subroutine ma7
      subroutine ma1(flag,my,ii,jj,t0,
     &               jour,mois,ia,hn,rr,gg,nno,fr1,fr2,nam)
      use parameters
      implicit none
      logical :: flag
      integer*4 :: l,nj,i,julien2,nd0,n,m,jour,mois,ia,kg,ng,j,k
      integer*4 :: ll,kk,ii,jj,my
      integer*4,dimension(8)              :: nam
      integer*4,dimension(30)             :: noa
      integer*4,dimension(30,8)           :: nno
      real                    :: t0,fq,co,so
      double precision        :: fnjd
      real,dimension(30)      :: f,q,q0,r1,r2,v0
      real,dimension(31)      :: g,r
      real,dimension(130)     :: hn
      real,dimension(362)     :: c
      real,dimension(91)      :: cs=(/
     &  1.000000D0,0.999848D0,0.999391D0,0.998630D0,0.997564D0,
     &                                            0.996195D0,0.994522D0,
     &  0.992546D0,0.990268D0,0.987688D0,0.984808D0,0.981627D0,
     &                                            0.978148D0,0.974370D0,
     &  0.970296D0,0.965926D0,0.961262D0,0.956305D0,0.951057D0,
     &                                            0.945519D0,0.939693D0,
     &  0.933580D0,0.927184D0,0.920505D0,0.913545D0,0.906308D0,
     &                                            0.898794D0,0.891007D0,
     &  0.882948D0,0.874620D0,0.866025D0,0.857167D0,0.848048D0,
     &                                            0.838671D0,0.829038D0,
     &  0.819152D0,0.809017D0,0.798636D0,0.788011D0,0.777146D0,
     &                                            0.766044D0,0.754710D0,
     &  0.743145D0,0.731354D0,0.719340D0,0.707107D0,0.694658D0,
     &                                            0.681998D0,0.669131D0,
     &  0.656059D0,0.642788D0,0.629320D0,0.615661D0,0.601815D0,
     &                                            0.587785D0,0.573576D0,
     &  0.559193D0,0.544639D0,0.529919D0,0.515038D0,0.500000D0,
     &                                            0.484810D0,0.469472D0,
     &  0.453991D0,0.438371D0,0.422618D0,0.406737D0,0.390731D0,
     &                                            0.374607D0,0.358368D0,
     &  0.342020D0,0.325568D0,0.309017D0,0.292372D0,0.275637D0,
     &                                            0.258819D0,0.241922D0,
     &  0.224951D0,0.207912D0,0.190809D0,0.173648D0,0.156434D0,
     &                                            0.139173D0,0.121869D0,
     &  0.104528D0,0.087156D0,0.069756D0,0.052336D0,0.034899D0,
     &                                           0.017452D0,0.000000D0/)
       real,dimension(30,8)    :: rr,gg,fr1,fr2
       real,dimension(11,3)    :: x,y
      equivalence (c(1),cs(1))
      l=1
      x(:,:)=0.0D0
      y(:,:)=0.0D0
      c(271)=0.0D0
      do i=1,90
       c(182-i)=-c(i)
       c(362-i)=c(i)
       c(180+i)=-c(i)
      end do
      if (flag) write(6,*)"ma1 c",c(:)
      nj=julien2(ia,jour,mois)-2415021
      fnjd=real(nj+l,8)-DBLE(1.4921875D0)
      nj=mod(nj,7)+1
      if(nj <= 0)nj=nj+7
      nd0=-1
      n=0
      m=1
      if (flag) write(6,*)"ma1 nj",nj
      if (flag) write(6,*)"ma1 fnjd",fnjd
      do kg=1,8
       n=nam(kg)
       if(n /= 0) then
        do i=1,n
         r(i)=rr(i,kg)
         g(i)=gg(i,kg)
         noa(i)=nno(i,kg)
         r1(i)=fr1(i,kg)
         r2(i)=fr2(i,kg)
         if (flag)
     &  write(6,*)"ma1 before ma3",i,r(i),g(i),r1(i),r2(i),noa(i)
        end do
        call ma3(f,v0,q,n,r,g,fnjd,ng,t0,r1,r2,noa)
        do i=1,n
         r(i)=r(i)*f(i)
         q0(i)=v0(i)-g(i)+360.0D0
         q(i)=mod(q(i),360.0D0)
         q0(i)=mod(q0(i),360.0D0)
         if(q0(i) < 0.0D0)q0(i)=q0(i)+360.0D0
         if (flag) write(6,*)"ma1 after ma3",i,f(i),q0(i),v0(i),q(i),
     &                      r(i),g(i),ng,r1(i),r2(i),noa(i)
        end do
        do i=l,3
         do j=1,n
          fq=q0(j)
          k=fq+1
          fq=fq-k+1
          co=((c(k+1)-c(k))*fq+c(k))*r(j)
          if (flag) write(6,*)"x,y",j,k,fq,co
          k=k-90
          if(k <= 0) k=k+360
          so=((c(k+1)-c(k))*fq+c(k))*r(j)
          x(ng,i)=x(ng,i)+co
          y(ng,i)=y(ng,i)+so
          q0(j)=q0(j)+q(j)
          if(q0(j) > 360.0D0)q0(j)=q0(j)-360.0D0
          if (flag) write(6,*)"x,y2",j,so,x(ng,i),y(ng,i),q0
         end do
        end do
       end if
      end do
      if (flag) write(6,*)'ma1 before ma4',x,y
      call ma4(hn,l,x,y)
      if (flag) write(6,*)'ma1 after ma4',hn
      end subroutine ma1
      subroutine ma4(hn,l,hm1,hm2)
      use parameters
      implicit none
      integer*4                        :: k,l,i
      real                 :: fnm
      real,dimension(130)  :: hn
      real,dimension(11,3) :: hm1,hm2
      real,dimension(128,3):: h
      complex,dimension(128)         :: x
      do  k=l,3
       h(:,k)=0.0D0
       x(:)=0.0D0
       x(1)=hm1(1,k)
       fnm=x(1)*2.0D0
       hm2(1,k)=0.0D0
       hm1(1,k)=0.0D0
       do i=2,11
        x(i)=cmplx(hm1(i,k),hm2(i,k))
        hm1(i,k)=0.0D0
        hm2(i,k)=0.0D0
        x(130-i)=conjg(x(i))
       end do
       call fft(x)
       h(:,k)=x(:)
      end do
      call ma5(h,hn)
      end subroutine ma4
      subroutine ma5(h,x)
      use parameters
      implicit none
      integer*4                        :: k
      real                 :: b,hh
      real,dimension(128,3):: h
      real,dimension(130)  :: x
      b=-0.5078125D0
      x(1)=x(129)
      x(2)=x(130)
      do k=1,128
       hh=h(k,2)+h(k,2)
       b=b+0.0078125D0
       x(k+2)=((h(k,1)+h(k,3)-hh)*b+h(k,3)-h(k,1))*b+hh
      end do
      if(x(1) == -valmanq)then
       x(2)=3.0D0*(x(3)-x(4))+x(5)
       x(1)=6.0D0*x(3)-8.0D0*x(4)+3.0D0*x(5)
      endif
      end subroutine ma5
      subroutine fft(a)
      use parameters
      implicit none
      integer*4                :: m,nv2,nm1,n,i,j,ip,l,k,le,le1
      real         :: fij
      complex,dimension(128) :: a
      complex                :: u,w,t
      m=7
      nv2=64
      nm1=127
      n=128
      j=1
      do  i=1,nm1
      if(i < j) then
       t=a(j)
       a(j)=a(i)
       a(i)=t
      end if
      k=nv2
    6 continue
      if(k < j) then
       j=j-k
       k=k/2
       goto 6
      end if
      j=j+k
      end do
      do l=1,m
      le=2**l
      le1=le/2
      u=cmplx(1.0D0,0.0D0)
      fij=pi/le1
      w=cexp(cmplx(0.0D0,fij))
      do j=1,le1
       do i=j,n,le
        ip=i+le1
        t=a(ip)*u
        a(ip)=a(i)-t
        a(i)=a(i)+t
       end do
       u=u*w
      end do
      end do
      end subroutine fft
      subroutine ma3(f,v0,q,n,r,g,fnjd,ng,t0,r1,r2,noa)
      use parameters
       implicit none
       integer*4                     :: i,ng,n
       integer*4,dimension(6)        :: nd
       integer*4,dimension(30)       :: noa
       real              :: fa1,fb1,fa2,fb2,a1,a2,asfov,t0
       double precision            :: fnjd
       real,dimension(30):: q,v0,f,r,g,r1,r2
       double precision            :: dq
       complex                     :: fc
      do i=1,n
        nd(6)=noa(i)
        nd(5)=noa(i)/10
        nd(4)=noa(i)/100
        nd(3)=noa(i)/1000
        nd(2)=noa(i)/10000
        nd(1)=noa(i)/100000
        nd(6)=nd(6)-nd(5)*10
        nd(5)=nd(5)-nd(4)*10
        nd(4)=nd(4)-nd(3)*10
        nd(3)=nd(3)-nd(2)*10
        nd(2)=nd(2)-nd(1)*10
        ng=nd(1)+1
        if(nd(1) == 0) r(i)=r(i)*2.0D0
        call masfo(nd,dq)
        q(i)=dq
        g(i)=g(i)+q(i)*t0/24.0D0
        g(i)=mod(g(i),360.0D0)
        if(g(i) < 0.0D0) g(i)=g(i)+360.0D0
        v0(i)=asfov(nd)+mod(REAL(dq*fnjd,8),DBLE(360.D0))
        fc=1.0D0
        if( r1(i) /= 0.0D0 .or. r2(i) /= 0.0D0 ) then
          fb2=9.242202D-4
          fb1=-fb2
          fa2=1.760045D0
          fa1=4.523139D0
          if(noa(i) == 275555)then
           fb1=2*fb2
           fa1=2*fa2
          endif
          a1=fa1+fb1*fnjd
          a2=fa2+fb2*fnjd
          fc=1.0D0+r1(i)*cmplx(cos(a1),sin(a1))+
     &    r2(i)*cmplx(cos(a2),sin(a2))
        end if
        v0(i)=v0(i)+atan2(aimag(fc),real(fc))*57.29578D0
        f(i)=cabs(fc)
        if(noa(i) == 355555.or.noa(i) == 382555 .or.noa(i) == 164555)
     &          v0(i)=v0(i)-90.0D0
      end do
      end subroutine ma3
      subroutine masfo(nd,asfo)
      use parameters
      implicit none
      integer*4,dimension(6) :: nd
      double precision       :: asfo
      asfo =(DBLE(360.D0)-DBLE(12.19074939D0))*(nd(1))
     &     +DBLE(13.17639673D0)*(nd(2)-5)+DBLE(0.98564734D0)*(nd(3)-5)
     &     +DBLE(0.11140408D0)*(nd(4)-5)+DBLE(0.05295392D0)*(nd(5)-5)+
     &      DBLE(0.00004707D0)*(nd(6)-5)
      end subroutine masfo
      function asfov(nd)
      use parameters
      implicit none
      integer*4,dimension(6) :: nd
      real       :: fov,asfov
      fov=280.1895D0*(nd(1)+nd(3)-5)+
     &  277.0248D0*(nd(2)-nd(1)-5)+334.3853D0*(nd(4)-5)
     &  +100.8432D0*(nd(5)-5)+281.2209D0*(nd(6)-5)+(mod(nd(1),2))*90.0D0
      asfov=(mod(fov,360.0D0))
      end function asfov
      function julien2(ia,jou,moi)
      use parameters
      implicit none
      integer*4                 :: 
     &                           i,julien2,ia,jou,moi,ibs,jour,mois,iy,m
      integer*4,dimension(12,2),PARAMETER :: jo=
     &     reshape((/0,31,59,90,120,151,181,212,243,273,304,334,
     &     0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
      real          :: a,b
      jour=jou
      mois=moi
      if(mois==1) then
      ibs=1
      if( ia >= 1582 .and.(ia /= 1582 .or. jou >= 277) ) then
       ibs=mod(ia,4)+2
       if(ibs /= 2)ibs=1
      end if
       do i=1,12
        mois=12-i+1
        jour=jou-jo(mois,ibs)
        if(jour > 0) goto 20
       end do
   20 continue
      end if
      a=ia+mois/100.0D0+jour/10000.0D0
      b=0.0D0
      iy=ia
      m=mois
      if(mois <= 2)then
        iy=iy-1
        m=m+12
      endif
      if(a >= 1582.1015D0)then
       b=2-int(iy/100)+int(iy/100)/4
      endif
      julien2=int(365.25D0*iy)+int(30.6001D0*(m+1))+jour+1720995+b
      end function julien2
