! ##################################################################################################################################

      MODULE LAPACK_DLAPY2_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      DOUBLE PRECISION FUNCTION DLAPY2_HELPER( X, Y )

      REAL(DOUBLE)   X, Y
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D0 )
      REAL(DOUBLE)   W, XABS, YABS, Z

      INTRINSIC          ABS, MAX, MIN, SQRT

      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2_HELPER = W
      ELSE
         DLAPY2_HELPER = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN

      END FUNCTION DLAPY2_HELPER

      END MODULE LAPACK_DLAPY2_HELPER
