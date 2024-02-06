      program cppcheck
      implicit none
      integer*4 input,iout, maxstring, lstring
      parameter (input=11, iout=12, maxstring=80)
      character*80 string
      integer*4 nwords                    , nswitches,nexample
     &                                  , max_switches
      parameter (nwords=16              , max_switches=1024)
      integer*4 istart(nwords), is        , size(max_switches)
     &       ,  iend(nwords), ie,ln     , line(max_switches)
      logical macro(nwords)             , example
      character*32                        switch(max_switches)
      integer*4 count, iocheck, i,k,n
      logical end_of_file, comment, word, command, new
      write(*,'(/1x,A,1x,A/)') 'This is CPPCHECK: Creating',
     &                   'new version of check_switches1.F.'
      example=.true.
      nexample=0
      nswitches=0
      do i=1,max_switches
        size(i)=0
        switch(i)='                                '
      enddo
      open(input,file='mergcpp.txt',status='old',form='formatted')
      count=0
      end_of_file=.false.
   1   count=count+1
        do i=1,maxstring
         string(i:i)=' '
        enddo
        read(input,'(A)',iostat=iocheck,end=2) string
        if (string(1:1).ne.'#') goto 1
        goto 3
   2    end_of_file=.true.
   3    lstring=maxstring+1
   4     lstring=lstring-1
         if ((string(lstring:lstring).eq.' ').AND.(lstring.GT.1)) goto 4
        n=0
        comment=.false.
        i=1
   5     i=i+1
          if (.not.comment .and. string(i:i+1).eq.'/*') then
            comment=.true.
            n=n+1
            istart(n)=i
          elseif (comment  .and. string(i:i+1).eq.'*/') then
            comment=.false.
            iend(n)=i+1
          endif
         if (i+1.lt.lstring) goto 5
        if (comment) then
          lstring=istart(n)-1
          n=n-1
        endif
        do k=1,n
          do i=istart(k),iend(k)
            string(i:i)=' '
          enddo
        enddo
        do i=1,lstring
          k=ichar(string(i:i))
          if (k.lt.48 .or. (k.gt.57 .and. k.lt.65) .or. (k.gt.90
     &    .and. k.lt.95) .or. k.eq.96 .or. k.gt.122) string(i:i)=' '
        enddo
        n=0
        word=.false.
        i=1
   6     i=i+1
          if (string(i:i).ne.' ' .and. .not.word) then
            word=.true.
            n=n+1
            istart(n)=i
          elseif (string(i:i).eq.' ' .and.  word) then
            word=.false.
            iend(n)=i-1
          endif
         if (i.lt.lstring) goto 6
        if (word) iend(n)=i
        command=.false.
        do k=1,n
          macro(k)=.false.
          is=istart(k)
          ie=iend(k)
          ln=ie-is+1
          if (ln.eq.6 .and. string(is:ie).eq.'define') then
            command=.true.
          elseif (ln.eq.5 .and. string(is:ie).eq.'undef') then
            command=.true.
          elseif (ln.eq.2 .and. string(is:ie).eq.'if') then
            command=.true.
            example=.false.
          elseif (ln.eq.5 .and. string(is:ie).eq.'ifdef') then
            command=.true.
            example=.false.
          elseif (ln.eq.7 .and. string(is:ie).eq.'defined') then
            command=.true.
            example=.false.
          elseif (ln.eq.4 .and. string(is:ie).eq.'elif') then
            command=.true.
            example=.false.
          elseif (ln.eq.4 .and. string(is:ie).eq.'else') then
          elseif (ln.eq.5 .and. string(is:ie).eq.'endif') then
          elseif (ln.eq.7 .and. string(is:ie).eq.'include') then
          elseif (command) then
            command=.false.
            macro(k)=.true.
          endif
        enddo
        do k=1,n
          if (macro(k)) then
            is=istart(k)
            ie=iend(k)
            ln=ie-is+1
            new=.true.
            do i=1,nswitches
              if (ln.eq.size(i)) then
                if (string(is:ie).eq.switch(i)(1:ln)) new=.false.
              endif
            enddo
            if (new) then
              nswitches=nswitches+1
              size(nswitches)=ln
              switch(nswitches)(1:ln)=string(is:ie)
              line(nswitches)=count
              if (example) nexample=nexample+1
            endif
          endif
        enddo
      if (.not.end_of_file) goto 1
      close(unit=input)
      open (unit=iout,file='check_switches1.F',form='formatted')
      write(iout,'(A/)')  '#include "cppdefs.h"'
      write(iout,'(/6x,A/)') 'subroutine check_switches1 (ierr)'
      write(iout,'(4(A,1x,A/),A,14x,A/A,1x,A,3x,A/A,14x,A/A,22x,A)')
     &  '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                                   'CODE, DO NOT EDIT! !!!!!!',
     &  '!!!!!! This file needs to be updated only if',
     &                                    'new CPP-switches  !!!!!!',
     &  '!!!!!! were introduced into "cppdefs.h".',
     &                                ' NO ACTION IS NEEDED  !!!!!!',
     &  '!!!!!! if changes in "cppdefs.h" are limited',
     &                                    'to activation or  !!!!!!',
     &  '!!!!!! deactivation of previously known switches.','!!!!!!',
     &  '!!!!!! To refresh this file compile and execute',
     &                                      '"cppcheck.F"', '!!!!!!',
     &  '!!!!!! as an independent program, or use commands','!!!!!!',
     &  '!!!!!! "make checkdefs" or "make depend".',        '!!!!!!'
      write(iout,'(A,20x,I3,1x,A/A,23x,I3,1x,A)')
     &  '!!!!!! Number of Configuration Choices:',nexample, '!!!!!!',
     &  '!!!!!! Total number of CPP-switches:', nswitches,  '!!!!!!'
      write(iout,'(2(/6x,A), 5(/A) /6x,A /5x,A1,6x,A, 5(/6x,A))')
     &  'implicit none',        'integer ierr, is,ie, iexample',
     &  '#include "param.h"',   '#include "strings.h"',
     &  '#ifdef MPI',           '# include "scalars.h"',    '#endif',
     &  'MPI_master_only write(stdout,''(/1x,A/)'')',       '&',
     &  '''Activated C-preprocessing Options:''',
     &  'do is=1,max_opt_size', '  Coptions(is:is)='' ''',  'enddo',
     &                                       'iexample=0',  'is=1'
      do i=1,nswitches
        ln=size(i)
        write(iout,'(A,1x,A)') '#ifdef', switch(i)(1:ln)
        if (i.le.nexample) write(iout,'(6x,A)') 'iexample=iexample+1'
        write(iout,'(6x,A,1x,A1,A,A1)')
     &         'MPI_master_only write(stdout,''(10x,A)'')',
     &                       '''', switch(i)(1:ln), ''''
        write(iout,'(6x,A7,I2/6x,A/6x,A,A,A1/6x,A/6x,A/A)')
     &       'ie=is +', ln-1, 'if (ie.ge.max_opt_size) goto 99',
     &       'Coptions(is:ie)=''', switch(i)(1:ln), '''',
     &       'Coptions(ie+1:ie+1)='' ''', 'is=ie+2', '#endif'
      enddo
      write(iout,'(6x,A/6x,A/8x,A/5x,A1,1x,A/8x,A/6x,A)')
     &     'MPI_master_only write(stdout,''(/)'')',
     &     'if (iexample.eq.0) then',
     &     'MPI_master_only write(stdout,''(1x,A)'')', '&',
     &   '''ERROR in "cppdefs.h": no configuration is specified.''',
     &     'ierr=ierr+1',  'elseif (iexample.gt.1) then'
      write(iout,'(8x,A/5x,A1,1x,A/8x,A/6x,A/6x,A)')
     &     'MPI_master_only write(stdout,''(1x,A)'')', '&',
     &   '''ERROR: more than one configuration in "cppdefs.h".''',
     &     'ierr=ierr+1', 'endif', 'return'
      write(iout,'(2x,A/5x,A1,2x,A,1x,A/5x,A1,2x,A,1x,A)')
     &   '99  MPI_master_only write(stdout,''(/1x,A,A/14x,A)'')',
     &   '&', '''CHECKDEFS -- ERROR: Unsufficient size of string',
     &     'Coptions'',', '&', '''in file "strings.h".'',',
     &                   '''Increase the size it and recompile.'''
      write(iout,'(6x,A,2(/6x,A))') 'ierr=ierr+1', 'return', 'end'
      close(unit=iout)
      stop
      end
