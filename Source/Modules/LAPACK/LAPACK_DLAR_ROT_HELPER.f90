! ##################################################################################################################################

MODULE LAPACK_DLAR_ROT_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARGV_HELPER( N, X, INCX, Y, INCY, C, INCC )

INTEGER            INCC, INCX, INCY, N
REAL(DOUBLE)       C( * ), X( * ), Y( * )
REAL(DOUBLE)       ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
INTEGER            I, IC, IX, IY
REAL(DOUBLE)       F, G, T, TT

INTRINSIC          ABS, SQRT

IX = 1
IY = 1
IC = 1
DO 10 I = 1, N
   F = X( IX )
   G = Y( IY )
   IF( G.EQ.ZERO ) THEN
      C( IC ) = ONE
   ELSE IF( F.EQ.ZERO ) THEN
      C( IC ) = ZERO
      Y( IY ) = ONE
      X( IX ) = G
   ELSE IF( ABS( F ).GT.ABS( G ) ) THEN
      T = G / F
      TT = SQRT( ONE+T*T )
      C( IC ) = ONE / TT
      Y( IY ) = T*C( IC )
      X( IX ) = F*TT
   ELSE
      T = F / G
      TT = SQRT( ONE+T*T )
      Y( IY ) = ONE / TT
      C( IC ) = T*Y( IY )
      X( IX ) = G*TT
   END IF
   IC = IC + INCC
   IY = IY + INCY
   IX = IX + INCX
10 CONTINUE

END SUBROUTINE DLARGV_HELPER

! ##################################################################################################################################

SUBROUTINE DLARTV_HELPER( N, X, INCX, Y, INCY, C, S, INCC )

INTEGER            INCC, INCX, INCY, N
REAL(DOUBLE)       C( * ), S( * ), X( * ), Y( * )
INTEGER            I, IC, IX, IY
REAL(DOUBLE)       XI, YI

IX = 1
IY = 1
IC = 1
DO 10 I = 1, N
   XI = X( IX )
   YI = Y( IY )
   X( IX ) = C( IC )*XI + S( IC )*YI
   Y( IY ) = C( IC )*YI - S( IC )*XI
   IX = IX + INCX
   IY = IY + INCY
   IC = IC + INCC
10 CONTINUE

END SUBROUTINE DLARTV_HELPER

END MODULE LAPACK_DLAR_ROT_HELPER
