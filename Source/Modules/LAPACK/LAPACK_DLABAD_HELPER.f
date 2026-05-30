! ##################################################################################################################################

      MODULE LAPACK_DLABAD_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DLABAD_HELPER( SMALL, LARGE )

      REAL(DOUBLE)   LARGE, SMALL

      INTRINSIC          LOG10, SQRT

      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF

      RETURN

      END SUBROUTINE DLABAD_HELPER

      END MODULE LAPACK_DLABAD_HELPER
