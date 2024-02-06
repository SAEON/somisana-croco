      program srcscheck
      implicit none
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      integer*4 input,iout, max_string
      parameter (input=11, iout=12, max_string=128)
      character tab*1, backslash*1, string*128
      integer*4 lstr, i, is,ie,indx,iocheck
      logical end_of_file, sources
      tab=char(9)
      backslash=char(92)
      sources=.false.
      end_of_file=.false.
      write(*,'(/1x,A,1x,A/)') 'This is SRCSCHECK:',
     &          'Creating new version of check_srcs.F.'
      open (unit=iout,file='check_srcs.F',form='formatted')
      write(iout,'(A/, /6x,A/, /3(A,1x,A/),A,29x,A/)')
     &     '#include "cppdefs.h"',    'subroutine check_srcs',
     &     '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                              'CODE, DO NOT EDIT! !!!!!!',
     &     '!!!!!! This file needs to be updated only if',
     &                               'the new files     !!!!!!',
     &     '!!!!!! were introduced into or deleted from',
     &                              'the list of source !!!!!!',
     &     '!!!!!! codes SRCS in the Makefile.',       '!!!!!!'
      write(iout,'(6x,A/6x,A/A/6x,A/8x,A/6x,A)')
     &     'implicit none', 'integer i', '#include "strings.h"',
     &     'do i=1,max_opt_size',  'srcs(i:i)='' ''',   'enddo'
      indx=1
      open(unit=input,file='Makefile')
   1  string=' '
       read(input,'(A)',iostat=iocheck,end=2) string
        goto 3
   2    end_of_file=.true.
   3    lstr=max_string
   4    if (string(lstr:lstr).eq.' ') then
          lstr=lstr-1
          if (lstr.gt.0) goto 4
        endif
        if (lstr.eq.0 .and. .not.end_of_file) goto 1
        do i=1,lstr
          if (string(i:i).eq.tab .or.
     &        string(i:i).eq.'=') string(i:i)=' '
        enddo
        is=1
   5    if (string(is:is).eq.' ') then
          is=is+1
          if (is.lt.lstr) goto 5
        endif
        if (string(is:is+5).eq.'SRCS90') then
          sources=.true.
          string(is:is+5)='    '
        endif
        if (string(is:is+3).eq.'SRCS') then
          sources=.true.
          string(is:is+3)='    '
        endif
        if (sources) then
          if (string(lstr:lstr).eq.backslash) then
              lstr=lstr-1
            if (lstr .gt. 0 ) then
   8          if (string(lstr:lstr).eq.' ') then
                lstr=lstr-1
                if (lstr.gt.0) goto 8
              endif
            endif
            if (lstr.eq.0 .and. .not.end_of_file) goto 1
          else
            sources=.false.
          endif
   6      if (string(is:is).eq.' ') then
            is=is+1
            if (is.lt.lstr) goto 6
          endif
          ie=is
   7      if (string(ie+1:ie+1).ne.' ') then
            ie=ie+1
            if (ie.lt.lstr) goto 7
          endif
          if (indx.eq.1) then
            write(iout,'(6x,A,I4,A1,I4,3A)') 'srcs(', indx, ':',
     &                   indx+ie-is, ')=''', string(is:ie), ''''
            indx=indx+ie-is
          else
            write(iout,'(6x,A,I4,A1,I4,3A)') 'srcs(', indx, ':',
     &                indx+ie-is+1, ')='' ', string(is:ie), ''''
            indx=indx+ie-is+1
          endif
          if (indx.gt.max_opt_size) then
            write(iout,'(A/2(A,1x,A/),A)') '*****', '*****   ERROR:',
     &           'parameter max_opt_size in file "strings.h" is not',
     &           '*****          sufficient to accomodate SCRS list',
     &           'from Makefile.', '*****'
            goto 99
          endif
          indx=indx+1
          is=ie+1
          if (is.le.lstr) goto 6
        endif
       if (.not.end_of_file) goto 1
      close(unit=input)
  99  write(iout,'(6x,A/6x,A)') 'return', 'end'
      close(iout)
      stop
      end
