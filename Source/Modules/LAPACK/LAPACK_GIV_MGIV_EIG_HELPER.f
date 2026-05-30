! ##################################################################################################################################

      MODULE LAPACK_GIV_MGIV_EIG_HELPER

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE_HELPER ( RANGE, INFO )

      USE PARAMS, ONLY                :  SUPINFO

      character range

      integer info

      Write(err,9902)
      if (supinfo == 'N') then
         Write(f06,9902)
      endif

      if      ((info == 1) .or. (info == 3) .and. (range /= 'I')) then
         Write(err,99021)
         Write(f06,99021)
      else if ((info == 2) .or. (info == 3) .and. (range == 'I')) then
         Write(err,99022)
         Write(f06,99022)
      else if (( info == 4) .and. (range == 'I')) then
         Write(err,803)
         Write(f06,803)
         fatal_err = fatal_err + 1
         call outa_here ( 'Y' )
      endif

 9902 format(' *INFORMATION: SOME OR ALL OF THE EIGENVALUES FAILED TO CO
     &NVERGE OR WERE NOT COMPUTED IN LAPACK SUBROUTINE DSTEBZ:')

99021 format(15x,'BISECTION FAILED TO CONVERGE FOR SOME EIGENVALUES; THE
     &SE EIGENVALUES ARE FLAGGED BY A NEGATIVE BLOCK NUMBER.',/,15X,
     &'THE EFFECT IS THAT THE EIGENVALUES MAY NOT BE AS ACCURATE AS THE
     &ABSOLUTE AND RELATIVE TOLERANCES.',/,15X,
     &'THIS IS GENERALLY CAUSED BY UNEXPECTEDLY INACCURATE ARITHMETIC.'
     &,/)

99022 format(15x,'NOT ALL OF THE EIGENVALUES IN THE RANGE REQUESTED WERE
     & FOUND:',/,15X,
     &'CAUSE: NON-MONOTONIC ARITHMETIC, CAUSING THE STURM SEQUENCE TO BE
     & NON-MONOTONIC.',/,15X,
     &'CURE : RECALCULATE, REQUESTING ALL EIGENVALUES',/)

  803 format(' *ERROR   803: PROGRAMMING ERROR IN SUBROUTINE DSTEBZ.'
     &,/,15X,'NO EIGENVALUES WERE COMPUTED BY LAPACK SUBROUTINE DSTEBZ.
     &THE GERSHGORIN INTERVAL INITIALLY USED WAS TOO SMALL.',/,15X,
     &'PROBABLE CAUSE: YOUR MACHINE HAS SLOPPY FLOATING-POINT ARITHMETIC
     &',/,15X,'CURE          : INCREASE THE PARAMETER "FUDGE" IN LAPACK
     &SUBROUTINE DSTEBZ, RECOMPILE, AND TRY AGAIN',/)

      END SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE_HELPER

      END MODULE LAPACK_GIV_MGIV_EIG_HELPER
