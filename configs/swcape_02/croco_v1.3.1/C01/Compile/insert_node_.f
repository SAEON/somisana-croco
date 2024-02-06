      subroutine insert_node (name, lstr, node, nnodes, ierr)
      implicit none
      character*(*) name, sffx*16
      integer*4 lstr, ierr, i,j,k, lsffx, digits, power, ndots, idot(3)
     &                                        , node,  nnodes
      parameter (digits=5)
      logical leading_dots
      ndots=0
      leading_dots=.true.
      do i=1,lstr
        if (name(i:i).eq.'.') then
         if (.not.leading_dots) then
            if (ndots.lt.3) then
              ndots=ndots+1
              idot(ndots)=i
            else
              write(*,'(/1x,4A/)') 'INSERT_NODE/INDEX ERROR: too ',
     &             'many dots in file name ''', name(1:lstr), '''.'
              ierr=ierr+1
              return
            endif
          endif
        else
          leading_dots=.false.
        endif
      enddo
      if(name(lstr-3:lstr-1) .eq. 'nc.') ndots=ndots-1
      lsffx=0
      if (ndots.gt.0) then
        i=idot(ndots)+1
   1    k=ichar(name(i:i))-48
        if ((k.lt.0 .or. k.gt.9) .and.
     &    name(i:i).ne.'*' .and. name(i:i).ne.'?') then
          lsffx=lstr-idot(ndots)+1
        elseif (i.lt.lstr) then
          i=i+1
          goto 1
        endif
      endif
      do j=1,ndots-1
        i=idot(j)+1
   2    k=ichar(name(i:i))-48
        if (k.lt.0 .or. k.gt.9) then
          if (name(i:i).ne.'*' .and. name(i:i).ne. '?') then
            write(*,'(/1x,2A/20x,3A/)') 'INSERT_NODE/INDEX ERROR: ',
     &                    'a non-digital character found in index ',
     &                    'segment of name ''',  name(1:lstr), '''.'
            ierr=ierr+1
          endif
        elseif (i.lt.idot(j+1)-1) then
          i=i+1
          goto 2
        endif
      enddo
      if (ierr.ne.0) return
      if (ndots.eq.0) then
        i=lstr+1
        j=lstr+1
        name(i:i)='.'
      else
        if (ndots.eq.1) then
          i=idot(1)
        elseif (ndots.eq.2) then
          if (idot(2)-idot(1).le.digits) then
            i=idot(1)
          else
            i=idot(2)
          endif
        else
          i=idot(ndots-1)
        endif
        if (lsffx.gt.0) then
          j=idot(ndots)
        else
          j=lstr+1
        endif
      endif
      lsffx=lstr+1-j
      if (lsffx.gt.0) sffx(1:lsffx)=name(j:lstr)
      k=node
      power=10
   3  if (nnodes.gt.power) then
        power=10*power
        goto 3
      endif
      if (power .ge. 10**digits) then
        write(*,'(/1x,2A/6x,2A/6x,A/)')  'INSERT_NODE/INDEX ERROR: ',
     &   'Possible ambiguity between MPI-node segment',    'length ',
     &   'and time index segment length. To fix: increase parameter',
     &   '''digits'' in file "insert_node.F" and recompile.'
        ierr=ierr+1
        return
      endif
   4  power=power/10
       i=i+1
       j=k/power
       name(i:i)=char(48+j)
       k=k-j*power
       if (power.gt.1) goto 4
      if (lsffx.gt.0) name(i+1:i+lsffx)=sffx(1:lsffx)
      lstr=i+lsffx
      return
      end
      subroutine   insert_time_index (name, lstr, indx, ierr)
      implicit none
      character*(*) name, sffx*16
      integer*4 lstr, ierr, i,j,k, lsffx, digits, power, ndots, idot(3)
     &                                        , indx
      parameter (digits=5)
      logical leading_dots
      ndots=0
      leading_dots=.true.
      do i=1,lstr
        if (name(i:i).eq.'.') then
         if (.not.leading_dots) then
            if (ndots.lt.3) then
              ndots=ndots+1
              idot(ndots)=i
            else
              write(*,'(/1x,4A/)') 'INSERT_NODE/INDEX ERROR: too ',
     &             'many dots in file name ''', name(1:lstr), '''.'
              ierr=ierr+1
              return
            endif
          endif
        else
          leading_dots=.false.
        endif
      enddo
      if(name(lstr-3:lstr-1) .eq. 'nc.') ndots=ndots-1
      lsffx=0
      if (ndots.gt.0) then
        i=idot(ndots)+1
   1    k=ichar(name(i:i))-48
        if ((k.lt.0 .or. k.gt.9) .and.
     &    name(i:i).ne.'*' .and. name(i:i).ne.'?') then
          lsffx=lstr-idot(ndots)+1
        elseif (i.lt.lstr) then
          i=i+1
          goto 1
        endif
      endif
      do j=1,ndots-1
        i=idot(j)+1
   2    k=ichar(name(i:i))-48
        if (k.lt.0 .or. k.gt.9) then
          if (name(i:i).ne.'*' .and. name(i:i).ne. '?') then
            write(*,'(/1x,2A/20x,3A/)') 'INSERT_NODE/INDEX ERROR: ',
     &                    'a non-digital character found in index ',
     &                    'segment of name ''',  name(1:lstr), '''.'
            ierr=ierr+1
          endif
        elseif (i.lt.idot(j+1)-1) then
          i=i+1
          goto 2
        endif
      enddo
      if (ierr.ne.0) return
      if (ndots.eq.0) then
        i=lstr+1
        j=lstr+1
        name(i:i)='.'
      else
        i=idot(1)
        if (ndots.eq.1) then
          if (lsffx.gt.0 .or. lstr-idot(1).lt.digits) then
            j=idot(1)
          else
            j=lstr+1
          endif
        elseif (ndots.eq.2 .and. idot(2)-idot(1).le.digits) then
          j=idot(1)
        else
          j=idot(2)
        endif
      endif
      lsffx=lstr+1-j
      if (lsffx.gt.0) sffx(1:lsffx)=name(j:lstr)
      k=indx
      power=10**digits
   4  power=power/10
       i=i+1
       j=k/power
       name(i:i)=char(48+j)
       k=k-j*power
       if (power.gt.1) goto 4
      if (lsffx.gt.0) name(i+1:i+lsffx)=sffx(1:lsffx)
      lstr=i+lsffx
      return
      end
      subroutine extract_time_index (name, lstr, indx, ierr)
      implicit none
      character*(*) name, sffx*16
      integer*4 lstr, ierr, i,j,k, lsffx, digits, power, ndots, idot(3)
     &                                        , indx
      parameter (digits=5)
      logical leading_dots
      ndots=0
      leading_dots=.true.
      do i=1,lstr
        if (name(i:i).eq.'.') then
         if (.not.leading_dots) then
            if (ndots.lt.3) then
              ndots=ndots+1
              idot(ndots)=i
            else
              write(*,'(/1x,4A/)') 'INSERT_NODE/INDEX ERROR: too ',
     &             'many dots in file name ''', name(1:lstr), '''.'
              ierr=ierr+1
              return
            endif
          endif
        else
          leading_dots=.false.
        endif
      enddo
      if(name(lstr-3:lstr-1) .eq. 'nc.') ndots=ndots-1
      lsffx=0
      if (ndots.gt.0) then
        i=idot(ndots)+1
   1    k=ichar(name(i:i))-48
        if ((k.lt.0 .or. k.gt.9) .and.
     &    name(i:i).ne.'*' .and. name(i:i).ne.'?') then
          lsffx=lstr-idot(ndots)+1
        elseif (i.lt.lstr) then
          i=i+1
          goto 1
        endif
      endif
      do j=1,ndots-1
        i=idot(j)+1
   2    k=ichar(name(i:i))-48
        if (k.lt.0 .or. k.gt.9) then
          if (name(i:i).ne.'*' .and. name(i:i).ne. '?') then
            write(*,'(/1x,2A/20x,3A/)') 'INSERT_NODE/INDEX ERROR: ',
     &                    'a non-digital character found in index ',
     &                    'segment of name ''',  name(1:lstr), '''.'
            ierr=ierr+1
          endif
        elseif (i.lt.idot(j+1)-1) then
          i=i+1
          goto 2
        endif
      enddo
      if (ierr.ne.0) return
      if (ndots.eq.1 .and. lsffx.eq.0) then
        i=idot(1)+1
        j=lstr
      elseif (ndots.gt.1) then
        i=idot(1)+1
        j=idot(2)-1
      else
        i=0
        j=0
      endif
      indx=0
      if (j-i+1.ge.digits) then
        do k=i,j
          indx=10*indx + ichar(name(k:k))-48
        enddo
      endif
      return
      end
