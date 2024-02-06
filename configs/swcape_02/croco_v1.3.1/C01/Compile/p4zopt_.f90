










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



 


 


MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :  1.0  !  2004     (O. Aumont) Original code
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!            3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!            3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt

   !!======================================================================
END MODULE p4zopt
