! ##################################################################################################################################

MODULE LAPACK_DLARF_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARF_HELPER( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )

CHARACTER          SIDE
INTEGER            INCV, LDC, M, N
REAL(DOUBLE)       TAU
REAL(DOUBLE)       C( LDC, * ), V( * ), WORK( * )
REAL(DOUBLE)       ONE, ZERO
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
LOGICAL            LSAME
EXTERNAL           LSAME
EXTERNAL           DGEMV, DGER

IF( LSAME( SIDE, 'L' ) ) THEN
   IF( TAU.NE.ZERO ) THEN
      CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                  WORK, 1 )
      CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
   END IF
ELSE
   IF( TAU.NE.ZERO ) THEN
      CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                  WORK, 1 )
      CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
   END IF
END IF

END SUBROUTINE DLARF_HELPER

END MODULE LAPACK_DLARF_HELPER
