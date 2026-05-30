! ##################################################################################################################################

      MODULE LAPACK_DLAEV2_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DLAEV2_HELPER( A, B, C, RT1, RT2, CS1, SN1 )

      REAL(DOUBLE)   A, B, C, CS1, RT1, RT2, SN1
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D0 )
      REAL(DOUBLE)   TWO
      PARAMETER          ( TWO = 2.0D0 )
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      REAL(DOUBLE)   HALF
      PARAMETER          ( HALF = 0.5D0 )
      INTEGER            SGN1, SGN2
      REAL(DOUBLE)   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN

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
         SGN1 = -1
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF

      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN

      END SUBROUTINE DLAEV2_HELPER

      END MODULE LAPACK_DLAEV2_HELPER
