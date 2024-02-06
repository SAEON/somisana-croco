










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



 


 


MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_pisces             ! Empty routine
   END SUBROUTINE trc_ini_pisces

   !!======================================================================
END MODULE trcini_pisces
