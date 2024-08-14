!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bounded (Lm,Mm,N, nz, z_r, z_w, 
     &                    rmask, var,z_lev, var_zlv, 
     &                    imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev. 
!


      implicit none

      integer Lm,Mm,N, nz,  imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,jmin:Mm+1,imin:Lm+1), 
     &       z_r(N,0:Mm+1,0:Lm+1),   rmask(0:Mm+1,0:Lm+1),
     &       z_w(0:N,0:Mm+1,0:Lm+1), z_lev(nz,jmin:Mm+1,imin:Lm+1),
     &       FillValue, var_zlv(nz,jmin:Mm+1,imin:Lm+1),
     &       zz(0:N+1,0:Lm+1), dpth
#ifdef CUBICNO
     &     , dz(kmin-1:N,0:Lm+1), FC(kmin-1:N,0:Lm+1), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max(  trd   *chunk_size  -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv


#ifdef VERBOSE
      write(*,'(2x,2A,3I3,2x,A,3I4)') 'sigma_to_z_intr: imin,jmin,',
     &   'kmin =', imin,jmin,kmin, 'trd,jstr,jend =', trd,jstr,jend
      write(*,*) Lm,Mm,N
#endif


      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,j,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,j,i)
              zz(N+1,i)=z_w(N,j,i)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,j,i)+z_r(k,j,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i-1)+z_w(0,j,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i-1)+z_w(N,j,i))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5*(z_r(k,j,i)+z_r(k,j-1,i))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i)+z_w(0,j-1,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i)+z_w(N,j-1,i))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_r(k,j,i)+z_r(k,j,i-1)
     &                      +z_r(k,j-1,i)+z_r(k,j-1,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.25D0*( z_w(0,j,i)+z_w(0,j,i-1)
     &                     +z_w(0,j-1,i)+z_w(0,j-1,i-1))

              zz(N+1,i)=0.25D0*( z_w(N,j,i)+z_w(N,j,i-1)
     &                      +z_w(N,j-1,i)+z_w(N,j-1,i-1))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,j,i)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,j,i)+z_w(k,j,i-1))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5*(z_w(k,j,i)+z_w(k,j-1,i))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_w(k,j,i)+z_w(k,j,i-1)
     &                      +z_w(k,j-1,i)+z_w(k,j-1,i-1))
              enddo
            enddo
          endif
        endif
#ifdef CUBICNO
        do k=kmin,N-1
          do i=imin,imax
            dz(k,i)=zz(k+1,i)-zz(k,i)
            FC(k,i)=var(k+1,j,i)-var(k,j,i)
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
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
       

          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(N,i)-zz(0,i) 
              if (rmask(j,i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)     
              if (rmask(j,i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2     !<-- below bottom
              elseif (dpth*(zz(1,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,j,i))*(z_lev(m,j,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,j,i)=FillValue             !<-- masked out
            elseif (km(i).eq.-2) then
              var_zlv(m,j,i)=var(1,j,i)      !<-- below bottom
            elseif (km(i).eq.N+2) then
              var_zlv(m,j,i)=var(N,j,i)     !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(m,j,i)=var(N,j,i)     !-> R-point, above z_r(N)
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(m,j,i)=var(1,j,i)      !<-- below bottom
            else
              k=km(i)
              !write(*,*) k,km
#ifdef CUBICNO
              cff=1.D0/(zz(k+1,i)-zz(k,i))
              p=z_lev(m,j,i)-zz(k,i)
              q=zz(k+1,i)-z_lev(m,j,i)

              var_zlv(m,j,i)=cff*( q*var(k,j,i) + p*var(k+1,j,i)
     &             -cff*p*q*( cff*(q-p)*(var(k+1,j,i)-var(k,j,i))
     &                                  +p*FC(k+1,i) -q*FC(k,i) )
     &                                                          )
#else
              var_zlv(m,j,i)=( var(k,j,i)*(zz(k+1,i)-z_lev(m,j,i))
     &                        +var(k+1,j,i)*(z_lev(m,j,i)-zz(k,i))
     &                               )/(zz(k+1,i)-zz(k,i))
#endif


            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m



#ifdef VERBOSE
              if ((zz(k+1,i)-z_lev(m,j,i))*(z_lev(m,j,i)-zz(k,i))
     & .lt.0.D0) then

                write(*,*) zz(k,i), zz(k+1,i), z_lev(m,j,i), m
                write(*,'(/1x,2A/)') '### ERROR: sigma_to_z_intr ',
     &                               ':: Search algoritm failure.'
                stop 
              endif
#endif
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
 
