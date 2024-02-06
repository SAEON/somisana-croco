      subroutine cancel_kwd (keyword,ierr)
      implicit none
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      character*(*) keyword
      integer*4 ierr, is,i,ie, lenkw,lenstr
      lenkw=lenstr(keyword)
      is=0
   1   is=is+1
       if (Coptions(is:is).eq.' ' .and. is.lt.max_opt_size) goto 1
      ie=is
   2   ie=ie+1
       if (Coptions(ie:ie).ne.' ' .and. ie.lt.max_opt_size) goto 2
      if (lenkw.eq.ie-is .and. Coptions(is:ie-1).eq.keyword) then
        do i=is,ie-1
          Coptions(i:i)=' '
        enddo
        return
      elseif (is.lt.max_opt_size) then
        is=ie
        goto 1
      endif
      write(*,'(2(/1x,A,1x,A,1x,A)/)') 'CANCEL_KW ERROR:',
     &        'Can not cancel keyword:',  keyword(1:lenkw),
     &        'Check SCRUM/ROMS input script for possible',
     &        'duplicated keyword.'
      ierr=ierr+1
      return
      end
      subroutine check_kwds (ierr)
      implicit none
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      integer*4 ierr, is,ie
      is=0
   1   is=is+1
       if (is .gt. max_opt_size) return
       if (Coptions(is:is) .eq. ' ') goto 1
      ie=is
   2   ie=ie+1
       if (Coptions(ie:ie).ne.' ' .and. ie.lt.max_opt_size) goto 2
      ierr=ierr+1
      write(*,'(/2(1x,A)/)') 'ERROR: keyword not found:',
     &                                  Coptions(is:ie-1)
      is=ie
      goto 1
      end
