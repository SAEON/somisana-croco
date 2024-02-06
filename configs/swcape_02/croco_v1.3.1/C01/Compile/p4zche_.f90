










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



 


 


MODULE p4zche
   !!======================================================================
   !!                         ***  MODULE p4zche  ***
   !! TOP :   PISCES Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!             3.6  !  2016-03  (O. Aumont) Change chemistry to MOCSY standards
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_che                   ! Empty routine
      WRITE(*,*) 'p4z_che: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_che


   !!======================================================================
END MODULE p4zche
