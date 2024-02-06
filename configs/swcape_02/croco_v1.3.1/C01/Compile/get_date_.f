      subroutine get_date (date_str)
      implicit none
      character*(*) date_str
      integer*4 year,hour, minute,sec, half, iday,imon, dstat,tstat,
     &        nday, lmonth(12), lday(31), len1, len2, len3, lenstr
      character*3  ampm(0:1)
      character*9 day(0:6), month(12)
      data lmonth/7,8,5,5,3,4,4,6,9,7,8,8/ ampm/' AM',' PM'/
     &    lday /9*1,22*2/ day /'Sunday   ', 'Monday   ', 'Tuesday  ',
     &            'Wednesday', 'Thursday ', 'Friday   ', 'Saturday '/
     &      month/'January  ', 'February ', 'March    ', 'April    ',
     &            'May      ', 'June     ', 'July     ', 'August   ',
     &            'September', 'October  ', 'November ', 'December '/
      character*11 ctime*11, today*18, fmt*20, wkday*44
      hour=0
      minute=0
      sec=0
      nday=1
      dstat=1
      tstat=1
      wkday=' '
      today=' '
      ctime=' '
      if (tstat.eq.0) then
        half=hour/12
        hour=hour-half*12
        if (hour.eq.0) hour=12
        if (half.eq.2) half=0
      endif
      if (dstat.eq.0) then
        write(fmt,10) lmonth(imon), lday(nday)
  10    format('(a',i1,',1x,i',i1,',1h,,1x,i4)')
        write(today,fmt) month(imon),nday,year
        wkday=day(iday)
      endif
      if(tstat.eq.0) then
        write(ctime,20) hour, minute, sec, ampm(half)
  20    format(i2,':',i2.2,':',i2.2,a3)
      endif
      len1=lenstr(wkday)
      len2=lenstr(today)
      len3=lenstr(ctime)
      date_str=wkday(1:len1)
      if (len2.gt.0) then
        len1=lenstr(date_str)
        date_str=date_str(1:len1)/ /' - '/ /today(1:len2)
      endif
      if (len3.gt.0) then
        len1=lenstr(date_str)
        date_str=date_str(1:len1)/ /' - '/ /ctime(1:len3)
      endif
      return
      end
