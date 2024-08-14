!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot_2d (Lm,N, nz, z_r, z_w, rmask, var,
     &                    z_lev, var_zlv, below,imin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev. 
!


      implicit none

      integer Lm,Mm,N, nz,  imin,imax, kmin, i,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,imin:Lm+1), 
     &       z_r(N,0:Lm+1),   rmask(0:Lm+1),
     &       z_w(0:N,0:Lm+1), z_lev(nz,imin:Lm+1),
     &       FillValue, var_zlv(nz,imin:Lm+1),
     &       zz(0:N+1,0:Lm+1), dpth, below
#ifdef CUBIC
     &     , dz(kmin-1:N,0:Lm+1), FC(kmin-1:N,0:Lm+1), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, below, imin,kmin, FillValue
Cf2py intent(out) var_zlv



        if (kmin.eq.1) then
          if (imin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,i)
              zz(N+1,i)=z_w(N,i)
            enddo
          elseif (imin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,i)+z_r(k,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,i-1)+z_w(0,i))
              zz(N+1,i)=0.5D0*(z_w(N,i-1)+z_w(N,i))
            enddo
          endif
        else
          if (imin.eq.0 ) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,i)
              enddo
            enddo
          elseif (imin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,i)+z_w(k,i-1))
              enddo
            enddo
          endif
        endif
#ifdef CUBIC
        do k=kmin,N-1
          do i=imin,imax
            dz(k,i)=zz(k+1,i)-zz(k,i)
            FC(k,i)=var(k+1,i)-var(k,i)
          enddo
        enddo
        do i=imin,imax
          dz(kmin-1,i)=dz(kmin,i)
          FC(kmin-1,i)=FC(kmin,i)

          dz(N,i)=dz(N-1,i)
          FC(N,i)=FC(N-1,i)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=imin,imax
            cff=FC(k,i)*FC(k-1,i)
            if (cff.gt.0.D0) then
              FC(k,i)=cff*(dz(k,i)+dz(k-1,i))/( (FC(k,i)+FC(k-1,i))
     &                                          *dz(k,i)*dz(k-1,i) )
            else
              FC(k,i)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
       

          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(N,i)-zz(0,i) 
              if (rmask(i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(m,i)-zz(N,i)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,i)).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)     
              if (rmask(i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(m,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(z_lev(m,i)-zz(N,i)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-below-z_lev(m,i)).gt.0.) then
                km(i)=-3     !<-- below bottom 
              elseif (dpth*(zz(0,i)-z_lev(m,i)).gt.0.) then
                km(i)=-2     !<-- below bottom but close 
              elseif (dpth*(zz(1,i)-z_lev(m,i)).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,i))*(z_lev(m,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,i)=FillValue             !<-- masked out
            elseif (km(i).eq.-2) then
#ifdef INTERP_BELOW
              var_zlv(m,i)=var(kmin,i)  !    
#ifdef CUBIC
     &           -FC(kmin,i)*(zz(kmin,i)-z_lev(m,i))
#else
     &          -(zz(kmin,i)-z_lev(m,i))*(var(kmin+1,i)-var(kmin,i))
     &                                  /(zz(kmin+1,i)-zz(kmin,i))
#endif
#else
              var_zlv(m,i)=FillValue      !<-- below bottom
#endif
            elseif (km(i).eq.N+2) then
#ifdef INTERP_ABOVE
              var_zlv(m,i)=var(N,i)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(N,i)*(z_lev(m,i)-zz(N,i))
#else
     &              +(z_lev(m,i)-zz(N,i))*(var(N,i)-var(N-1,i))
     &                                   /(zz(N,i)-zz(N-1,i))
#endif
#else
              var_zlv(m,i)=-FillValue     !<-- above surface
#endif
            elseif (km(i).eq.N) then
              var_zlv(m,i)=var(N,i)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(N,i)*(z_lev(m,i)-zz(N,i))
#else
     &              +(z_lev(m,i)-zz(N,i))*(var(N,i)-var(N-1,i))
     &                                   /(zz(N,i)-zz(N-1,i))
#endif
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(m,i)=var(kmin,i)  !     but above bottom  
#ifdef CUBIC
     &           -FC(kmin,i)*(zz(kmin,i)-z_lev(m,i))
#else
     &          -(zz(kmin,i)-z_lev(m,i))*(var(kmin+1,i)-var(kmin,i))
     &                                  /(zz(kmin+1,i)-zz(kmin,i))
#endif
            else
              k=km(i)
              !write(*,*) k,km
#ifdef CUBIC
              cff=1.D0/(zz(k+1,i)-zz(k,i))
              p=z_lev(m,i)-zz(k,i)
              q=zz(k+1,i)-z_lev(m,i)

              var_zlv(m,i)=cff*( q*var(k,i) + p*var(k+1,i)
     &             -cff*p*q*( cff*(q-p)*(var(k+1,i)-var(k,i))
     &                                  +p*FC(k+1,i) -q*FC(k,i) )
     &                                                          )
#else
              var_zlv(m,i)=( var(k,i)*(zz(k+1,i)-z_lev(m,i))
     &                        +var(k+1,i)*(z_lev(m,i)-zz(k,i))
     &                               )/(zz(i,k+1)-zz(i,k))
#endif


            !write(*,*) 'bof',i,k,zz(i,k), zz(i,k+1), z_lev(i,m), m



            endif
          enddo
        enddo ! <-- m


      return
      end
 
