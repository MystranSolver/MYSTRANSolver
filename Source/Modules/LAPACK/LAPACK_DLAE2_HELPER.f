! ##################################################################################################################################

      MODULE LAPACK_DLAE2_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DLAE2_HELPER( A, B, C, RT1, RT2 )

      REAL(DOUBLE)   A, B, C, RT1, RT2
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D0 )
      REAL(DOUBLE)   TWO
      PARAMETER          ( TWO = 2.0D0 )
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      REAL(DOUBLE)   HALF
      PARAMETER          ( HALF = 0.5D0 )
      REAL(DOUBLE)   AB, ACMN, ACMX, ADF, DF, RT, SM, TB

      INTRINSIC          ABS, SQRT

      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN

      END SUBROUTINE DLAE2_HELPER

      END MODULE LAPACK_DLAE2_HELPER
