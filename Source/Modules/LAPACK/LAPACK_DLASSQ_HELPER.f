! ##################################################################################################################################

      MODULE LAPACK_DLASSQ_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DLASSQ_HELPER( N, X, INCX, SCALE, SUMSQ )

      INTEGER            INCX, N
      REAL(DOUBLE)   SCALE, SUMSQ
      REAL(DOUBLE)   X( * )
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            IX
      REAL(DOUBLE)   ABSXI

      INTRINSIC          ABS

      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN

      END SUBROUTINE DLASSQ_HELPER

      END MODULE LAPACK_DLASSQ_HELPER
