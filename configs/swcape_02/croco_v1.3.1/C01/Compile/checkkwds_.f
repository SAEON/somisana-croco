      program checkkwds
      implicit none
      integer*4 input,iout, maxstring, nwords
      parameter (input=11, iout=12, maxstring=80, nwords=16)
      character string*80, backslash*1
      logical cont_switch, end_of_file
      integer*4 is,i,ie, kwlen, lstring, iocheck, count
      backslash=char(92)
      write(*,'(/1x,A,1x,A/)') 'This is CHECKKWDS: Creating',
     &                        'new version of "setup_kwds.F".'
      open(iout,file='setup_kwds.F', form='formatted')
      write(iout,'(A/,/6x,A/3(/A,1x,A)/A,8x,A/A,1x,A/A,39x,A)')
     &  '#include "cppdefs.h"',    'subroutine setup_kwds (ierr)',
     &  '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                                   'CODE, DO NOT EDIT! !!!!!!',
     &  '!!!!!! This file needs to be updated only if',
     &                                    'new keywords were !!!!!!',
     &  '!!!!!! introduced into "read_inp.F".',
     &                            'To create or refresh this !!!!!!',
     &  '!!!!!! file use compile and execute "checkkwds.F" as an',
     &  '!!!!!!', '!!!!!! independent program, or use',
     &                          'commands "make checkkwds"   !!!!!!',
     &  '!!!!!! or "make depend".',                         '!!!!!!'
      write(iout,'(2(/6x,A), 5(/A) /6x,A /8x,A, 2(/6x,A))')
     &  'implicit none', 'integer ierr, is,ie', '#include "param.h"',
     &  '#include "strings.h"', '#ifdef MPI','# include "scalars.h"',
     &  '#endif',   'do is=1,max_opt_size',  'Coptions(is:is)='' ''',
     &  'enddo',    'is=1'
      open(input,file='read_inp.F',status='old',form='formatted')
      count=0
      cont_switch=.false.
      end_of_file=.false.
   1   count=count+1
        do i=1,maxstring
         string(i:i)=' '
        enddo
        read(input,'(A)',iostat=iocheck,end=2) string
        goto 3
   2    end_of_file=.true.
   3    lstring=maxstring
   4    if (string(lstring:lstring).eq.' ') then
          lstring=lstring-1
          if (lstring.gt.0) goto 4
        endif
        if (lstring.eq.0 .and. end_of_file) goto 11
        if (lstring.eq.0       .or. string(1:1).eq.'!' .or.
     &      string(1:1).eq.'C' .or. string(1:1).eq.'c') goto 1
        if (string(1:1).eq.'#') then
         is=2
   5      if (string(is:is).ne.'i' .and. is.lt.lstring-6) then
            is=is+1
            goto 5
          elseif (string(is:is+6).eq.'include') then
            goto 1
          endif
        endif
        if (string(1:1).eq.'#' .or. cont_switch) then
          write(iout,'(A)') string(1:lstring)
          if (string(lstring:lstring).eq.backslash) then
            cont_switch=.true.
          else
            cont_switch=.false.
          endif
          goto 1
        endif
        is=7
   6     if (string(is:is).ne.'k' .and. is.lt.lstring-6) then
           is=is+1
           goto 6
         elseif (string(is:is+6).eq.'keyword') then
           is=is+10
   7       if (string(is:is).ne.'k' .and. is.lt.lstring-4) then
             is=is+1
             goto 7
           elseif (string(is:is+4).eq.'kwlen') then
             is=is+6
   8         if (string(is:is).ne.'.' .and. is.lt.lstring-4) then
               is=is+1
               goto 8
             elseif (string(is:is+3).eq. '.eq.') then
               is=is+4
   9           i=ichar(string(is:is))
               if (i.ne.39 .and. is.lt.lstring) then
                 is=is+1
                 goto 9
               elseif (i.eq.39) then
                 ie=is+1
  10             i=ichar(string(ie:ie))
                 if (i.ne.39 .and. ie.lt.lstring) then
                   ie=ie+1
                   goto 10
                 elseif (i.eq.39) then
                   kwlen=ie-is-1
                   write(iout,'(6x,A,I2/6x,A/6x,A,A/6x,A/6x,A)')
     &                      'ie=is +', kwlen,
     &                      'if (ie.ge.max_opt_size) goto 99',
     &                      'Coptions(is:ie)=',     string(is:ie),
     &                      'Coptions(ie+1:ie+1)='' ''', 'is=ie+2'
                 endif
               endif
             endif
           endif
         endif
       goto 1
  11  close(input)
      write(iout,'(6x,A)') 'return'
      write(iout,'(2x,A/5x,A1,2x,A,1x,A/5x,A1,2x,A,1x,A)')
     &   '99  MPI_master_only write(stdout,''(/1x,A,A/14x,A)'')',
     &   '&', '''SETUP_KWDS ERROR: Unsufficient size of string',
     &     'Coptions'',', '&', '''in file "strings.h".'',',
     &                   '''Increase the size it and recompile.'''
      write(iout,'(6x,A,2(/6x,A))') 'ierr=ierr+1', 'return', 'end'
      close(iout)
      stop
      end
