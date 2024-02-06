










! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!







                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      

                      
                      
                      
                      
                      
                      

                      


                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                     
                      
!    RVTK test (Restartability or Parallel reproducibility)                








                      
                      
                      
                      
                      


! $Id: set_global_definitions.h 1616 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!


 










 



 
  

























































































 
 

 
      



! $Id: set_global_definitions.h 1618 2014-12-18 14:39:51Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!





















!









!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh



 


 


MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys

   !!======================================================================
END MODULE p4zlys
