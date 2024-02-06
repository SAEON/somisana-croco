      integer*4 function lenstr (string)
      implicit none
      integer*4 is, ie
      character*(*) string
      ie=len(string)
   1   if (string(ie:ie).eq.' ') then
         ie=ie-1
         if (ie.gt.0) goto 1
       endif
      is=0
   2   is=is+1
       if (string(is:is).eq.' ' .and. is.lt.ie) goto 2
      lenstr=ie-is+1
      if (is.gt.1) string=string(is:ie)
      return
      end
