! ##################################################################################################################################

      MODULE LAPACK_BLAS_AUX
! --- lapack_surgery begin --- !

! This is  the set of LAPACK auxiliary routines called by other LAPACK subroutines

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                  :  ERR, F06, SC1
      USE SCONTR, ONLY                  :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                  :  HOUR, MINUTE, SEC,
     &                                     SFRAC, TSEC
      USE PARAMS, ONLY                  :  NOCOUNTS
      USE LAPACK_DLACON_HELPER
      USE LAPACK_DLAGTS_HELPER
      USE LAPACK_DLAN_HELPER
      USE LAPACK_DLABAD_HELPER
      USE LAPACK_DLACPY_HELPER
      USE LAPACK_DLAE2_HELPER
      USE LAPACK_DLAEV2_HELPER
      USE LAPACK_DLARF_HELPER
      USE LAPACK_DLARFB_HELPER
      USE LAPACK_DLARFG_HELPER
      USE LAPACK_DLARFT_HELPER
      USE LAPACK_DLAR_ROT_HELPER
      USE LAPACK_DLARTG_HELPER
      USE LAPACK_DLAPY2_HELPER
      USE LAPACK_DLASCL_HELPER
      USE LAPACK_DLAS_MISC_HELPER
      USE LAPACK_DLASRT_HELPER
      USE LAPACK_DLASSQ_HELPER
      USE LAPACK_DISNAN_HELPER

      USE OUTA_HERE_Interface

      character(1*byte), parameter      :: cr13_lba = char(13)

      CONTAINS

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  BLAS routines

! ##################################################################################################################################
! 001 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return

      END FUNCTION DASUM

! ##################################################################################################################################
! 002 LAPACK_BLAS_AUX

      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return

      END SUBROUTINE DAXPY

! ##################################################################################################################################
! 003 LAPACK_BLAS_AUX

      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return

      END SUBROUTINE DCOPY

! ##################################################################################################################################
! 004 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return

      END FUNCTION DDOT

! ##################################################################################################################################
! 005 LAPACK_BLAS_AUX

      SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA, BETA
      INTEGER            INCX, INCY, KL, KU, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  KL     - INTEGER.
*           On entry, KL specifies the number of sub-diagonals of the
*           matrix A. KL must satisfy  0 .le. KL.
*           Unchanged on exit.
*
*  KU     - INTEGER.
*           On entry, KU specifies the number of super-diagonals of the
*           matrix A. KU must satisfy  0 .le. KU.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry, the leading ( kl + ku + 1 ) by n part of the
*           array A must contain the matrix of coefficients, supplied
*           column by column, with the leading diagonal of the matrix in
*           row ( ku + 1 ) of the array, the first super-diagonal
*           starting at position 2 in row ku, the first sub-diagonal
*           starting at position 1 in row ( ku + 2 ), and so on.
*           Elements in the array A that do not correspond to elements
*           in the band matrix (such as the top left ku by ku triangle)
*           are not referenced.
*           The following program segment will transfer a band matrix
*           from conventional full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    K = KU + 1 - J
*                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
*                       A( K + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( kl + ku + 1 ).
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL(DOUBLE).
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL(DOUBLE) array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY,
     $                   LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( KL.LT.0 )THEN
         INFO = 4
      ELSE IF( KU.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( KL + KU + 1 ) )THEN
         INFO = 8
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 10
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGBMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the band part of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KUP1 = KU + 1
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  K    = KUP1 - J
                  DO 50, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( I ) = Y( I ) + TEMP*A( K + I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  K    = KUP1 - J
                  DO 70, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( IY ) = Y( IY ) + TEMP*A( K + I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               IF( J.GT.KU )
     $            KY = KY + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               K    = KUP1 - J
               DO 90, I = MAX( 1, J - KU ), MIN( M, J + KL )
                  TEMP = TEMP + A( K + I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               K    = KUP1 - J
               DO 110, I = MAX( 1, J - KU ), MIN( M, J + KL )
                  TEMP = TEMP + A( K + I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
               IF( J.GT.KU )
     $            KX = KX + INCX
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGBMV .
*
      END SUBROUTINE DGBMV

! ##################################################################################################################################
! 008 LAPACK_BLAS_AUX

      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END SUBROUTINE DGER

! ##################################################################################################################################
! 009 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )

*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      REAL(DOUBLE)                  X( * )
*     ..
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      REAL(DOUBLE)      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      REAL(DOUBLE)      ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END FUNCTION DNRM2

! ##################################################################################################################################
! 010 LAPACK_BLAS_AUX

      SUBROUTINE  DROT (N,DX,INCX,DY,INCY,C,S)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dy(*),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return

      END SUBROUTINE DROT

! ##################################################################################################################################
! 011 LAPACK_BLAS_AUX

      SUBROUTINE DSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY, col_num )

*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA, BETA
      INTEGER            INCX, INCY, K, LDA, N, col_num
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSBMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric band matrix, with k super-diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the band matrix A is being supplied as
*           follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  being supplied.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  being supplied.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry, K specifies the number of super-diagonals of the
*           matrix A. K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the symmetric matrix, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer the upper
*           triangular part of a symmetric band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the symmetric matrix, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer the lower
*           triangular part of a symmetric band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL(DOUBLE).
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  Y      - REAL(DOUBLE) array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.

! col_num - Integer
!           The col num corresponding to the vector y being multiplied
!           (for use when there are multiple calls to dsbmv to mult
!           a matrix times each col of another matrix).
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( K.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSBMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN

*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of the array A
*     are accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when upper triangle of A is stored.
*
         KPLUS1 = K + 1
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF (NOCOUNTS .NE. 'Y') THEN
                  write(sc1,12345,advance='no') j,n,col_num,cr13_lba
               ENDIF
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               L     = KPLUS1 - J
               DO 50, I = MAX( 1, J - K ), J - 1
                  Y( I ) = Y( I ) + TEMP1*A( L + I, J )
                  TEMP2  = TEMP2  + A( L + I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               L     = KPLUS1 - J
               DO 70, I = MAX( 1, J - K ), J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( L + I, J )
                  TEMP2   = TEMP2   + A( L + I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               IF( J.GT.K )THEN
                  KX = KX + INCX
                  KY = KY + INCY
               END IF
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when lower triangle of A is stored.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               IF (NOCOUNTS .NE. 'Y') THEN
                  write(sc1,12345,advance='no') j,n,col_num,cr13_lba
               ENDIF
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( 1, J )
               L      = 1            - J
               DO 90, I = J + 1, MIN( N, J + K )
                  Y( I ) = Y( I ) + TEMP1*A( L + I, J )
                  TEMP2  = TEMP2  + A( L + I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( 1, J )
               L       = 1             - J
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, MIN( N, J + K )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( L + I, J )
                  TEMP2   = TEMP2   + A( L + I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
12345 format(7X,'mult row ',i8,' of ',i8,' times col ',i8, a)

*     End of DSBMV .
*
      END SUBROUTINE DSBMV

! ##################################################################################################################################
! 014 LAPACK_BLAS_AUX

      SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )

*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL(DOUBLE).
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END SUBROUTINE DSYMV

! ##################################################################################################################################
! 015 LAPACK_BLAS_AUX

      SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )

*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in upper triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYR  .
*
      END SUBROUTINE DSYR

! ##################################################################################################################################
! 016 LAPACK_BLAS_AUX

      SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )

*     .. Scalar Arguments ..
      REAL(DOUBLE)   ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYR2 .
*
      END SUBROUTINE DSYR2

! ##################################################################################################################################
! 017 LAPACK_BLAS_AUX

      SUBROUTINE DSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )

*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL(DOUBLE)   ALPHA, BETA
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR2K  performs one of the symmetric rank 2k operations
*
*     C := alpha*A*B' + alpha*B*A' + beta*C,
*
*  or
*
*     C := alpha*A'*B + alpha*B'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*  matrices in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
*                                        beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrices  A and B.  K must be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL(DOUBLE) array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL(DOUBLE).
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - REAL(DOUBLE) array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL(DOUBLE)   TEMP1, TEMP2
*     .. Parameters ..
      REAL(DOUBLE)   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2K', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  C := alpha*A*B' + alpha*B*A' + C.
*
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*A'*B + alpha*B'*A + C.
*
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYR2K.
*
      END SUBROUTINE DSYR2K

! ##################################################################################################################################
! 018 LAPACK_BLAS_AUX

      SUBROUTINE DSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )

*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      REAL(DOUBLE)   ALPHA, BETA
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYRK  performs one of the symmetric rank k operations
*
*     C := alpha*A*A' + beta*C,
*
*  or
*
*     C := alpha*A'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*  in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrix  A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL(DOUBLE).
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - REAL(DOUBLE) array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL(DOUBLE)   TEMP
*     .. Parameters ..
      REAL(DOUBLE)   ONE ,         ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYRK ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  C := alpha*A*A' + beta*C.
*
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*A'*A + beta*C.
*
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYRK .
*
      END SUBROUTINE DSYRK

! ##################################################################################################################################
! 019 LAPACK_BLAS_AUX

      SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX,
     &                   dtbsv_msg )     ! my addition

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DTBSV'
*     .. Scalar Arguments ..
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
      character*1        dtbsv_msg
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTBSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U' or 'u', K specifies the number of
*           super-diagonals of the matrix A.
*           On entry with UPLO = 'L' or 'l', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer an upper
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer a lower
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that when DIAG = 'U' or 'u' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
!  dtbsv_msg (input) CHARACTER
!           = 'Y', have subr DTBSV print Fwd, Back pass messages
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*


! **********************************************************************************************************************************
**     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTBSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  if (dtbsv_msg == 'Y' .AND. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,22345,advance='no') j,cr13_lba
                  endif
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  if (dtbsv_msg == 'Y' .AND. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,22345,advance='no') j,cr13_lba
                  endif
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  if (dtbsv_msg == 'Y' .and. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,12345,advance='no') j,n,cr13_lba
                  endif
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  if (dtbsv_msg == 'Y' .and. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,12345,advance='no') j,n,cr13_lba
                  endif
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A')*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  if (dtbsv_msg == 'Y' .and. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,12345,advance='no') j,n,cr13_lba
                  endif
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  if (dtbsv_msg == 'Y' .and. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,12345,advance='no') j,n,cr13_lba
                  endif
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )
     $               KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  if (dtbsv_msg == 'Y' .AND. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,22345,advance='no') j,cr13_lba
                  endif
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  if (dtbsv_msg == 'Y' .AND. (NOCOUNTS .NE. 'Y')) then
                     write(sc1,22345,advance='no') j,cr13_lba
                  endif
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )
     $               KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
12345 format(7X,'Forward  pass, row ',i8,' to ',i8, a)

22345 format(7X,'Backward pass, row ',i8,' to        1 ', a)

*     End of DTBSV .
*
! **********************************************************************************************************************************
 9000 continue            ! My lines

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE DTBSV

! ##################################################################################################################################
! 020 LAPACK_BLAS_AUX

      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )

*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL(DOUBLE)   ALPHA
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL(DOUBLE).
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL(DOUBLE) array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL(DOUBLE)   TEMP
*     .. Parameters ..
      REAL(DOUBLE)   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*A'*B.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END SUBROUTINE DTRMM

! ##################################################################################################################################
! 021 LAPACK_BLAS_AUX

      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )

*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL(DOUBLE) array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - REAL(DOUBLE) array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL(DOUBLE)   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END SUBROUTINE DTRMV

! ##################################################################################################################################
! 023 LAPACK_BLAS_AUX

      INTEGER FUNCTION IDAMAX(N,DX,INCX)

c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(DOUBLE) dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return

      END FUNCTION IDAMAX

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Auxiliary routines

! ##################################################################################################################################
! 024 LAPACK_BLAS_AUX

      SUBROUTINE DLABAD( SMALL, LARGE )
      REAL(DOUBLE)   LARGE, SMALL

      CALL DLABAD_HELPER( SMALL, LARGE )

      END SUBROUTINE DLABAD

! ##################################################################################################################################
! 025 LAPACK_BLAS_AUX

      SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE, itmax )  ! My itmax

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DLACON'

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      integer            itmax
      REAL(DOUBLE)   EST
      INTEGER            ISGN( * )
      REAL(DOUBLE)   V( * ), X( * )
      CALL DLACON_HELPER( N, V, X, ISGN, EST, KASE, itmax )
      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE DLACON

! ##################################################################################################################################
! 026 LAPACK_BLAS_AUX

      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
      REAL(DOUBLE)   A( LDA, * ), B( LDB, * )
      CALL DLACPY_HELPER( UPLO, M, N, A, LDA, B, LDB )
      RETURN
*
*     End of DLACPY
*
      END SUBROUTINE DLACPY

! ##################################################################################################################################
! 027 LAPACK_BLAS_AUX

      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(DOUBLE)   A, B, C, RT1, RT2
      CALL DLAE2_HELPER( A, B, C, RT1, RT2 )
      RETURN

      END SUBROUTINE DLAE2

! ##################################################################################################################################
! 028 LAPACK_BLAS_AUX

      SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,
     $                   RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,
     $                   NAB, WORK, IWORK, INFO )
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN,
     $                   NITMAX
      REAL(DOUBLE)   ABSTOL, PIVMIN, RELTOL
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      REAL(DOUBLE)   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),
     $                   WORK( * )

      CALL DLAEBZ_HELPER( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,
     $                    RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,
     $                    NAB, WORK, IWORK, INFO )

      RETURN

      END SUBROUTINE DLAEBZ

! **********************************************************************************************************************************

      SUBROUTINE DLAEBZ_HELPER( IJOB, NITMAX, N, MMAX, MINP, NBMIN,
     $                          ABSTOL, RELTOL, PIVMIN, D, E, E2, NVAL,
     $                          AB, C, MOUT, NAB, WORK, IWORK, INFO )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
      REAL(DOUBLE)   ABSTOL, PIVMIN, RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      REAL(DOUBLE)   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAEBZ contains the iteration loops which compute and use the
*  function N(w), which is the count of eigenvalues of a symmetric
*  tridiagonal matrix T less than or equal to its argument  w.  It
*  performs a choice of two types of loops:
*
*  IJOB=1, followed by
*  IJOB=2: It takes as input a list of intervals and returns a list of
*          sufficiently small intervals whose union contains the same
*          eigenvalues as the union of the original intervals.
*          The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
*          The output interval (AB(j,1),AB(j,2)] will contain
*          eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
*
*  IJOB=3: It performs a binary search in each input interval
*          (AB(j,1),AB(j,2)] for a point  w(j)  such that
*          N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
*          the search.  If such a w(j) is found, then on output
*          AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
*          (AB(j,1),AB(j,2)] will be a small interval containing the
*          point where N(w) jumps through NVAL(j), unless that point
*          lies outside the initial interval.
*
*  Note that the intervals are in all cases half-open intervals,
*  i.e., of the form  (a,b] , which includes  b  but not  a .
*
*  To avoid underflow, the matrix should be scaled so that its largest
*  element is no greater than  overflow**(1/2) * underflow**(1/4)
*  in absolute value.  To assure the most accurate computation
*  of small eigenvalues, the matrix should be scaled to be
*  not much smaller than that, either.
*
*  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*  Matrix", Report CS41, Computer Science Dept., Stanford
*  University, July 21, 1966
*
*  Note: the arguments are, in general, *not* checked for unreasonable
*  values.
*
*  Arguments
*  =========
*
*  IJOB    (input) INTEGER
*          Specifies what is to be done:
*          = 1:  Compute NAB for the initial intervals.
*          = 2:  Perform bisection iteration to find eigenvalues of T.
*          = 3:  Perform bisection iteration to invert N(w), i.e.,
*                to find a point which has a specified number of
*                eigenvalues of T to its left.
*          Other values will cause DLAEBZ to return with INFO=-1.
*
*  NITMAX  (input) INTEGER
*          The maximum number of "levels" of bisection to be
*          performed, i.e., an interval of width W will not be made
*          smaller than 2^(-NITMAX) * W.  If not all intervals
*          have converged after NITMAX iterations, then INFO is set
*          to the number of non-converged intervals.
*
*  N       (input) INTEGER
*          The dimension n of the tridiagonal matrix T.  It must be at
*          least 1.
*
*  MMAX    (input) INTEGER
*          The maximum number of intervals.  If more than MMAX intervals
*          are generated, then DLAEBZ will quit with INFO=MMAX+1.
*
*  MINP    (input) INTEGER
*          The initial number of intervals.  It may not be greater than
*          MMAX.
*
*  NBMIN   (input) INTEGER
*          The smallest number of intervals that should be processed
*          using a vector loop.  If zero, then only the scalar loop
*          will be used.
*
*  ABSTOL  (input) REAL(DOUBLE)
*          The minimum (absolute) width of an interval.  When an
*          interval is narrower than ABSTOL, or than RELTOL times the
*          larger (in magnitude) endpoint, then it is considered to be
*          sufficiently small, i.e., converged.  This must be at least
*          zero.
*
*  RELTOL  (input) REAL(DOUBLE)
*          The minimum relative width of an interval.  When an interval
*          is narrower than ABSTOL, or than RELTOL times the larger (in
*          magnitude) endpoint, then it is considered to be
*          sufficiently small, i.e., converged.  Note: this should
*          always be at least radix*machine epsilon.
*
*  PIVMIN  (input) REAL(DOUBLE)
*          The minimum absolute value of a "pivot" in the Sturm
*          sequence loop.  This *must* be at least  max |e(j)**2| *
*          safe_min  and at least safe_min, where safe_min is at least
*          the smallest number that can divide one without overflow.
*
*  D       (input) REAL(DOUBLE) array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
*  E       (input) REAL(DOUBLE) array, dimension (N)
*          The offdiagonal elements of the tridiagonal matrix T in
*          positions 1 through N-1.  E(N) is arbitrary.
*
*  E2      (input) REAL(DOUBLE) array, dimension (N)
*          The squares of the offdiagonal elements of the tridiagonal
*          matrix T.  E2(N) is ignored.
*
*  NVAL    (input/output) INTEGER array, dimension (MINP)
*          If IJOB=1 or 2, not referenced.
*          If IJOB=3, the desired values of N(w).  The elements of NVAL
*          will be reordered to correspond with the intervals in AB.
*          Thus, NVAL(j) on output will not, in general be the same as
*          NVAL(j) on input, but it will correspond with the interval
*          (AB(j,1),AB(j,2)] on output.
*
*  AB      (input/output) REAL(DOUBLE) array, dimension (MMAX,2)
*          The endpoints of the intervals.  AB(j,1) is  a(j), the left
*          endpoint of the j-th interval, and AB(j,2) is b(j), the
*          right endpoint of the j-th interval.  The input intervals
*          will, in general, be modified, split, and reordered by the
*          calculation.
*
*  C       (input/output) REAL(DOUBLE) array, dimension (MMAX)
*          If IJOB=1, ignored.
*          If IJOB=2, workspace.
*          If IJOB=3, then on input C(j) should be initialized to the
*          first search point in the binary search.
*
*  MOUT    (output) INTEGER
*          If IJOB=1, the number of eigenvalues in the intervals.
*          If IJOB=2 or 3, the number of intervals output.
*          If IJOB=3, MOUT will equal MINP.
*
*  NAB     (input/output) INTEGER array, dimension (MMAX,2)
*          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).
*          If IJOB=2, then on input, NAB(i,j) should be set.  It must
*             satisfy the condition:
*             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),
*             which means that in interval i only eigenvalues
*             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,
*             NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with
*             IJOB=1.
*             On output, NAB(i,j) will contain
*             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of
*             the input interval that the output interval
*             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the
*             the input values of NAB(k,1) and NAB(k,2).
*          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),
*             unless N(w) > NVAL(i) for all search points  w , in which
*             case NAB(i,1) will not be modified, i.e., the output
*             value will be the same as the input value (modulo
*             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)
*             for all search points  w , in which case NAB(i,2) will
*             not be modified.  Normally, NAB should be set to some
*             distinctive value(s) before DLAEBZ is called.
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (MMAX)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (MMAX)
*          Workspace.
*
*  INFO    (output) INTEGER
*          = 0:       All intervals converged.
*          = 1--MMAX: The last INFO intervals did not converge.
*          = MMAX+1:  More than MMAX intervals were generated.
*
*  Further Details
*  ===============
*
*      This routine is intended to be called only by other LAPACK
*  routines, thus the interface is less user-friendly.  It is intended
*  for two purposes:
*
*  (a) finding eigenvalues.  In this case, DLAEBZ should have one or
*      more initial intervals set up in AB, and DLAEBZ should be called
*      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.
*      Intervals with no eigenvalues would usually be thrown out at
*      this point.  Also, if not all the eigenvalues in an interval i
*      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.
*      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest
*      eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX
*      no smaller than the value of MOUT returned by the call with
*      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1
*      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the
*      tolerance specified by ABSTOL and RELTOL.
*
*  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).
*      In this case, start with a Gershgorin interval  (a,b).  Set up
*      AB to contain 2 search intervals, both initially (a,b).  One
*      NVAL element should contain  f-1  and the other should contain  l
*      , while C should contain a and b, resp.  NAB(i,1) should be -1
*      and NAB(i,2) should be N+1, to flag an error if the desired
*      interval does not lie in (a,b).  DLAEBZ is then called with
*      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --
*      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while
*      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r
*      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and
*      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and
*      w(l-r)=...=w(l+k) are handled similarly.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
*     ..
*     .. Local Scalars ..
      INTEGER            ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL,
     $                   KLNEW
      REAL(DOUBLE)   TMP1, TMP2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Check for Errors
*
      INFO = 0
      IF( IJOB.LT.1 .OR. IJOB.GT.3 ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Initialize NAB
*
      IF( IJOB.EQ.1 ) THEN
*
*        Compute the number of eigenvalues in the initial intervals.
*
         MOUT = 0
*DIR$ NOVECTOR
         DO 30 JI = 1, MINP
            DO 20 JP = 1, 2
               TMP1 = D( 1 ) - AB( JI, JP )
               IF( ABS( TMP1 ).LT.PIVMIN )
     $            TMP1 = -PIVMIN
               NAB( JI, JP ) = 0
               IF( TMP1.LE.ZERO )
     $            NAB( JI, JP ) = 1
*
               DO 10 J = 2, N
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP )
                  IF( ABS( TMP1 ).LT.PIVMIN )
     $               TMP1 = -PIVMIN
                  IF( TMP1.LE.ZERO )
     $               NAB( JI, JP ) = NAB( JI, JP ) + 1
   10          CONTINUE
   20       CONTINUE
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 )
   30    CONTINUE
         RETURN
      END IF
*
*     Initialize for loop
*
*     KF and KL have the following meaning:
*        Intervals 1,...,KF-1 have converged.
*        Intervals KF,...,KL  still need to be refined.
*
      KF = 1
      KL = MINP
*
*     If IJOB=2, initialize C.
*     If IJOB=3, use the user-supplied starting point.
*
      IF( IJOB.EQ.2 ) THEN
         DO 40 JI = 1, MINP
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
   40    CONTINUE
      END IF
*
*     Iteration loop
*
      DO 130 JIT = 1, NITMAX
*
*        Loop over intervals
*
         IF( KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0 ) THEN
*
*           Begin of Parallel Version of the loop
*
            DO 60 JI = KF, KL
*
*              Compute N(c), the number of eigenvalues less than c
*
               WORK( JI ) = D( 1 ) - C( JI )
               IWORK( JI ) = 0
               IF( WORK( JI ).LE.PIVMIN ) THEN
                  IWORK( JI ) = 1
                  WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
               END IF
*
               DO 50 J = 2, N
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI )
                  IF( WORK( JI ).LE.PIVMIN ) THEN
                     IWORK( JI ) = IWORK( JI ) + 1
                     WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
                  END IF
   50          CONTINUE
   60       CONTINUE
*
            IF( IJOB.LE.2 ) THEN
*
*              IJOB=2: Choose all intervals containing eigenvalues.
*
               KLNEW = KL
               DO 70 JI = KF, KL
*
*                 Insure that N(w) is monotone
*
                  IWORK( JI ) = MIN( NAB( JI, 2 ),
     $                          MAX( NAB( JI, 1 ), IWORK( JI ) ) )
*
*                 Update the Queue -- add intervals if both halves
*                 contain eigenvalues.
*
                  IF( IWORK( JI ).EQ.NAB( JI, 2 ) ) THEN
*
*                    No eigenvalue in the upper interval:
*                    just use the lower interval.
*
                     AB( JI, 2 ) = C( JI )
*
                  ELSE IF( IWORK( JI ).EQ.NAB( JI, 1 ) ) THEN
*
*                    No eigenvalue in the lower interval:
*                    just use the upper interval.
*
                     AB( JI, 1 ) = C( JI )
                  ELSE
                     KLNEW = KLNEW + 1
                     IF( KLNEW.LE.MMAX ) THEN
*
*                       Eigenvalue in both intervals -- add upper to
*                       queue.
*
                        AB( KLNEW, 2 ) = AB( JI, 2 )
                        NAB( KLNEW, 2 ) = NAB( JI, 2 )
                        AB( KLNEW, 1 ) = C( JI )
                        NAB( KLNEW, 1 ) = IWORK( JI )
                        AB( JI, 2 ) = C( JI )
                        NAB( JI, 2 ) = IWORK( JI )
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
   70          CONTINUE
               IF( INFO.NE.0 )
     $            RETURN
               KL = KLNEW
            ELSE
*
*              IJOB=3: Binary search.  Keep only the interval containing
*                      w   s.t. N(w) = NVAL
*
               DO 80 JI = KF, KL
                  IF( IWORK( JI ).LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = C( JI )
                     NAB( JI, 1 ) = IWORK( JI )
                  END IF
                  IF( IWORK( JI ).GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = C( JI )
                     NAB( JI, 2 ) = IWORK( JI )
                  END IF
   80          CONTINUE
            END IF
*
         ELSE
*
*           End of Parallel Version of the loop
*
*           Begin of Serial Version of the loop
*
            KLNEW = KL
            DO 100 JI = KF, KL
*
*              Compute N(w), the number of eigenvalues less than w
*
               TMP1 = C( JI )
               TMP2 = D( 1 ) - TMP1
               ITMP1 = 0
               IF( TMP2.LE.PIVMIN ) THEN
                  ITMP1 = 1
                  TMP2 = MIN( TMP2, -PIVMIN )
               END IF
*
*              A series of compiler directives to defeat vectorization
*              for the next loop
*
*$PL$ CMCHAR=' '
CDIR$          NEXTSCALAR
C$DIR          SCALAR
CDIR$          NEXT SCALAR
CVD$L          NOVECTOR
CDEC$          NOVECTOR
CVD$           NOVECTOR
*VDIR          NOVECTOR
*VOCL          LOOP,SCALAR
CIBM           PREFER SCALAR
*$PL$ CMCHAR='*'
*
               DO 90 J = 2, N
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1
                  IF( TMP2.LE.PIVMIN ) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN( TMP2, -PIVMIN )
                  END IF
   90          CONTINUE
*
               IF( IJOB.LE.2 ) THEN
*
*                 IJOB=2: Choose all intervals containing eigenvalues.
*
*                 Insure that N(w) is monotone
*
                  ITMP1 = MIN( NAB( JI, 2 ),
     $                    MAX( NAB( JI, 1 ), ITMP1 ) )
*
*                 Update the Queue -- add intervals if both halves
*                 contain eigenvalues.
*
                  IF( ITMP1.EQ.NAB( JI, 2 ) ) THEN
*
*                    No eigenvalue in the upper interval:
*                    just use the lower interval.
*
                     AB( JI, 2 ) = TMP1
*
                  ELSE IF( ITMP1.EQ.NAB( JI, 1 ) ) THEN
*
*                    No eigenvalue in the lower interval:
*                    just use the upper interval.
*
                     AB( JI, 1 ) = TMP1
                  ELSE IF( KLNEW.LT.MMAX ) THEN
*
*                    Eigenvalue in both intervals -- add upper to queue.
*
                     KLNEW = KLNEW + 1
                     AB( KLNEW, 2 ) = AB( JI, 2 )
                     NAB( KLNEW, 2 ) = NAB( JI, 2 )
                     AB( KLNEW, 1 ) = TMP1
                     NAB( KLNEW, 1 ) = ITMP1
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
*
*                 IJOB=3: Binary search.  Keep only the interval
*                         containing  w  s.t. N(w) = NVAL
*
                  IF( ITMP1.LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = TMP1
                     NAB( JI, 1 ) = ITMP1
                  END IF
                  IF( ITMP1.GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  END IF
               END IF
  100       CONTINUE
            KL = KLNEW
*
*           End of Serial Version of the loop
*
         END IF
*
*        Check for convergence
*
         KFNEW = KF
         DO 110 JI = KF, KL
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) )
            TMP2 = MAX( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) )
            IF( TMP1.LT.MAX( ABSTOL, PIVMIN, RELTOL*TMP2 ) .OR.
     $          NAB( JI, 1 ).GE.NAB( JI, 2 ) ) THEN
*
*              Converged -- Swap with position KFNEW,
*                           then increment KFNEW
*
               IF( JI.GT.KFNEW ) THEN
                  TMP1 = AB( JI, 1 )
                  TMP2 = AB( JI, 2 )
                  ITMP1 = NAB( JI, 1 )
                  ITMP2 = NAB( JI, 2 )
                  AB( JI, 1 ) = AB( KFNEW, 1 )
                  AB( JI, 2 ) = AB( KFNEW, 2 )
                  NAB( JI, 1 ) = NAB( KFNEW, 1 )
                  NAB( JI, 2 ) = NAB( KFNEW, 2 )
                  AB( KFNEW, 1 ) = TMP1
                  AB( KFNEW, 2 ) = TMP2
                  NAB( KFNEW, 1 ) = ITMP1
                  NAB( KFNEW, 2 ) = ITMP2
                  IF( IJOB.EQ.3 ) THEN
                     ITMP1 = NVAL( JI )
                     NVAL( JI ) = NVAL( KFNEW )
                     NVAL( KFNEW ) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  110    CONTINUE
         KF = KFNEW
*
*        Choose Midpoints
*
         DO 120 JI = KF, KL
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
  120    CONTINUE
*
*        If no more intervals to refine, quit.
*
         IF( KF.GT.KL )
     $      GO TO 140
  130 CONTINUE
*
*     Converged
*
  140 CONTINUE
      INFO = MAX( KL+1-KF, 0 )
      MOUT = KL
*
      RETURN
*
*     End of DLAEBZ_HELPER
*
      END SUBROUTINE DLAEBZ_HELPER

! ##################################################################################################################################
! 029 LAPACK_BLAS_AUX

      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(DOUBLE)   A, B, C, CS1, RT1, RT2, SN1
      CALL DLAEV2_HELPER( A, B, C, RT1, RT2, CS1, SN1 )
      RETURN
*
*     End of DLAEV2
*
      END SUBROUTINE DLAEV2

! ##################################################################################################################################
! 030 LAPACK_BLAS_AUX

      SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, JOB, N
      REAL(DOUBLE)   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      REAL(DOUBLE)   A( * ), B( * ), C( * ), D( * ), Y( * )
*     ..
*
      CALL DLAGTS_HELPER( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
      END SUBROUTINE DLAGTS

! ##################################################################################################################################
! 032 LAPACK_BLAS_AUX
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T

      CALL DLAMC1_HELPER( BETA, T, RND, IEEE1 )

      RETURN

      END SUBROUTINE DLAMC1

! **********************************************************************************************************************************

      SUBROUTINE DLAMC1_HELPER( BETA, T, RND, IEEE1 )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL(DOUBLE)   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1_HELPER
*
      END SUBROUTINE DLAMC1_HELPER
*
! ##################################################################################################################################
! 033 LAPACK_BLAS_AUX
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL(DOUBLE)   EPS, RMAX, RMIN

      CALL DLAMC2_HELPER( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

      RETURN

      END SUBROUTINE DLAMC2

! **********************************************************************************************************************************

      SUBROUTINE DLAMC2_HELPER( BETA, T, RND, EPS, EMIN, RMIN, EMAX,
     $                          RMAX )

*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL(DOUBLE)   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) REAL(DOUBLE)
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) REAL(DOUBLE)
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) REAL(DOUBLE)
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      REAL(DOUBLE)   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2_HELPER
*
      END SUBROUTINE DLAMC2_HELPER
*
! ##################################################################################################################################
! 034 LAPACK_BLAS_AUX
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(DOUBLE)   A, B
      DLAMC3 = DLAMC3_HELPER( A, B )
*
      RETURN
*
*     End of DLAMC3
*
      END FUNCTION DLAMC3
*
! ##################################################################################################################################
! 035 LAPACK_BLAS_AUX
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
      INTEGER            BASE, EMIN
      REAL(DOUBLE)   START

      CALL DLAMC4_HELPER( EMIN, START, BASE )

      RETURN

      END SUBROUTINE DLAMC4

! **********************************************************************************************************************************

      SUBROUTINE DLAMC4_HELPER( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL(DOUBLE)   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) REAL(DOUBLE)
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      REAL(DOUBLE)   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4_HELPER
*
      END SUBROUTINE DLAMC4_HELPER
*
! ##################################################################################################################################
! 036 LAPACK_BLAS_AUX
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL(DOUBLE)   RMAX

      CALL DLAMC5_HELPER( BETA, P, EMIN, IEEE, EMAX, RMAX )

      RETURN

      END SUBROUTINE DLAMC5

! **********************************************************************************************************************************

      SUBROUTINE DLAMC5_HELPER( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL(DOUBLE)   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) REAL(DOUBLE)
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL(DOUBLE)   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5_HELPER
*
      END SUBROUTINE DLAMC5_HELPER

! ##################################################################################################################################
! 037 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DLANSB( NORM, UPLO, N, K, AB, LDAB,
     $                 WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
      REAL(DOUBLE)   AB( LDAB, * ), WORK( * )
      DLANSB = DLANSB_HELPER( NORM, UPLO, N, K, AB, LDAB, WORK )
      RETURN
*
*     End of DLANSB

      END FUNCTION DLANSB

! ##################################################################################################################################
! 039 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
      REAL(DOUBLE)   A( LDA, * ), WORK( * )
      DLANSY = DLANSY_HELPER( NORM, UPLO, N, A, LDA, WORK )
      RETURN
*
*     End of DLANSY
*
      END FUNCTION DLANSY

! ##################################################################################################################################
! 040 LAPACK_BLAS_AUX

      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL(DOUBLE)   X, Y
      DLAPY2 = DLAPY2_HELPER( X, Y )
      RETURN

      END FUNCTION DLAPY2

! ##################################################################################################################################
! 041 LAPACK_BLAS_AUX

      SUBROUTINE DLAQSB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED
     &                 , thresh, small, large )            ! My line
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            KD, LDAB, N
      REAL(DOUBLE)   AMAX, SCOND
      REAL(DOUBLE)   thresh, small, large              ! My line
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   AB( LDAB, * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAQSB equilibrates a symmetric band matrix A using the scaling
*  factors in the vector S.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) REAL(DOUBLE) array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U'*U or A = L*L' of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*          The scale factors for A.                               --------------------
*
*  SCOND   (input) REAL(DOUBLE)
*          Ratio of the smallest S(i) to the largest S(i).
*
*  AMAX    (input) REAL(DOUBLE)
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies whether or not equilibration was done.
*          = 'N':  No equilibration.
*          = 'Y':  Equilibration was done, i.e., A has been replaced by
*                  diag(S) * A * diag(S).

!  thresh  (input) REAL(DOUBLE)
!          threshold value (see comment below)    ! My line

!  small   (output) REAL(DOUBLE)
!          threshold value (see comment below)    ! My line

!  large   (output) REAL(DOUBLE)
!          threshold value (see comment below)    ! My line
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if scaling should be done
*  based on the ratio of the scaling factors.  If SCOND < THRESH,
*  scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if scaling should
*  be done based on the absolute size of the largest matrix element.
*  If AMAX > LARGE or AMAX < SMALL, scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL(DOUBLE)   CJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL(DOUBLE)       DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
      ELSE
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A is stored in band format.
*
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = MAX( 1, J-KD ), J
                  AB( KD+1+I-J, J ) = CJ*S( I )*AB( KD+1+I-J, J )
   10          CONTINUE
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, MIN( N, J+KD )
                  AB( 1+I-J, J ) = CJ*S( I )*AB( 1+I-J, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of DLAQSB
*
      END SUBROUTINE DLAQSB

! ##################################################################################################################################
! 042 LAPACK_BLAS_AUX

      SUBROUTINE DLAR2V( N, X, Y, Z, INCX, C, S, INCC )
      INTEGER            INCC, INCX, N
      REAL(DOUBLE)       C( * ), S( * ), X( * ), Y( * ), Z( * )

      CALL DLAR2V_HELPER( N, X, Y, Z, INCX, C, S, INCC )
      RETURN

      END SUBROUTINE DLAR2V

      SUBROUTINE DLAR2V_HELPER( N, X, Y, Z, INCX, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   C( * ), S( * ), X( * ), Y( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAR2V applies a vector of real plane rotations from both sides to
*  a sequence of 2-by-2 real symmetric matrices, defined by the elements
*  of the vectors x, y and z. For i = 1,2,...,n
*
*     ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
*     ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be applied.
*
*  X       (input/output) REAL(DOUBLE) array,
*                         dimension (1+(N-1)*INCX)
*          The vector x.
*
*  Y       (input/output) REAL(DOUBLE) array,
*                         dimension (1+(N-1)*INCX)
*          The vector y.
*
*  Z       (input/output) REAL(DOUBLE) array,
*                         dimension (1+(N-1)*INCX)
*          The vector z.
*
*  INCX    (input) INTEGER
*          The increment between elements of X, Y and Z. INCX > 0.
*
*  C       (input) REAL(DOUBLE) array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  S       (input) REAL(DOUBLE) array, dimension (1+(N-1)*INCC)
*          The sines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C and S. INCC > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IX
      REAL(DOUBLE)   CI, SI, T1, T2, T3, T4, T5, T6, XI, YI, ZI
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IX )
         ZI = Z( IX )
         CI = C( IC )
         SI = S( IC )
         T1 = SI*ZI
         T2 = CI*ZI
         T3 = T2 - SI*XI
         T4 = T2 + SI*YI
         T5 = CI*XI + T1
         T6 = CI*YI - T1
         X( IX ) = CI*T5 + SI*T4
         Y( IX ) = CI*T6 - SI*T3
         Z( IX ) = CI*T4 - SI*T5
         IX = IX + INCX
         IC = IC + INCC
   10 CONTINUE
*
*     End of DLAR2V
*
      RETURN

      END SUBROUTINE DLAR2V_HELPER

! ##################################################################################################################################
! 043 LAPACK_BLAS_AUX

      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL(DOUBLE)   TAU
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   C( LDC, * ), V( * ), WORK( * )
*     ..
*
      CALL DLARF_HELPER( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      END SUBROUTINE DLARF

! ##################################################################################################################################
! 044 LAPACK_BLAS_AUX

      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
! B 02/07/04 //////////////////////////////////////////////////////////B
      REAL(DOUBLE)   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
! Double the first dimension

!     REAL(DOUBLE)   C( LDC, * ), T( LDT, * ), V( LDV, * ),
!    $                   WORK( 2*LDWORK, * )
! E ////////////////////////////////////////////////////////////////////E
*     ..
*
*  Purpose
*  =======
*
*  DLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) REAL(DOUBLE) array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) REAL(DOUBLE) array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) REAL(DOUBLE) array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Executable Statements ..
*
      CALL DLARFB_HELPER( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                    T, LDT, C, LDC, WORK, LDWORK )
      RETURN
      END SUBROUTINE DLARFB

! ##################################################################################################################################
! 045 LAPACK_BLAS_AUX

      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL(DOUBLE)   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   X( * )
*     ..
*
      CALL DLARFG_HELPER( N, ALPHA, X, INCX, TAU )
      END SUBROUTINE DLARFG

! ##################################################################################################################################
! 046 LAPACK_BLAS_AUX

      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
      CALL DLARFT_HELPER( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
      END SUBROUTINE DLARFT

! ##################################################################################################################################
! 047 LAPACK_BLAS_AUX

      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   C( * ), X( * ), Y( * )
*     ..
*
      CALL DLARGV_HELPER( N, X, INCX, Y, INCY, C, INCC )
      END SUBROUTINE DLARGV

! ##################################################################################################################################
! 048 LAPACK_BLAS_AUX

      SUBROUTINE DLARNV( IDIST, ISEED, N, X )
      INTEGER            IDIST, N
      INTEGER            ISEED( 4 )
      REAL(DOUBLE)   X( * )

      CALL DLARNV_HELPER( IDIST, ISEED, N, X )

      RETURN

      END SUBROUTINE DLARNV

! **********************************************************************************************************************************

      SUBROUTINE DLARNV_HELPER( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL(DOUBLE)   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARNV returns a vector of n random real numbers from a uniform or
*  normal distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  uniform (0,1)
*          = 2:  uniform (-1,1)
*          = 3:  normal (0,1)
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated.
*
*  X       (output) REAL(DOUBLE) array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine DLARUV to generate random
*  real numbers from a uniform (0,1) distribution, in batches of up to
*  128 using vectorisable code. The Box-Muller method is used to
*  transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      REAL(DOUBLE)   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
*     ..
*     .. Local Arrays ..
      REAL(DOUBLE)   U( LV )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
*     ..
*     .. External Subroutines ..
*     ..
*     .. Executable Statements ..
*
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
*
*        Call DLARUV to generate IL2 numbers from a uniform (0,1)
*        distribution (IL2 <= LV)
*
         CALL DLARUV( ISEED, IL2, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           Copy generated numbers
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           Convert generated numbers to uniform (-1,1) distribution
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           Convert generated numbers to normal (0,1) distribution
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
*
*     End of DLARNV_HELPER
*
      END SUBROUTINE DLARNV_HELPER

! ##################################################################################################################################
! 049 LAPACK_BLAS_AUX

      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      REAL(DOUBLE)   CS, F, G, R, SN
*     ..
*
      CALL DLARTG_HELPER( F, G, CS, SN, R )
      END SUBROUTINE DLARTG

! ##################################################################################################################################
! 050 LAPACK_BLAS_AUX

      SUBROUTINE DLARTV( N, X, INCX, Y, INCY, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   C( * ), S( * ), X( * ), Y( * )
*     ..
*
      CALL DLARTV_HELPER( N, X, INCX, Y, INCY, C, S, INCC )
      END SUBROUTINE DLARTV

! ##################################################################################################################################
! 051 LAPACK_BLAS_AUX

      SUBROUTINE DLARUV( ISEED, N, X )
      INTEGER            N
      INTEGER            ISEED( 4 )
      REAL(DOUBLE)   X( * )

      CALL DLARUV_HELPER( ISEED, N, X )

      RETURN

      END SUBROUTINE DLARUV

! **********************************************************************************************************************************

      SUBROUTINE DLARUV_HELPER( ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL(DOUBLE)   X( N )
*     ..
*
*  Purpose
*  =======
*
*  DLARUV returns a vector of n random real numbers from a uniform (0,1)
*  distribution (n <= 128).
*
*  This is an auxiliary routine called by DLARNV and ZLARNV.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated. N <= 128.
*
*  X       (output) REAL(DOUBLE) array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE
      PARAMETER          ( ONE = 1.0D0 )
      INTEGER            LV, IPW2
      REAL(DOUBLE)   R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
*     ..
*     .. Local Arrays ..
      INTEGER            MM( LV, 4 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, MOD
*     ..
*     .. Data statements ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508,
     $                   2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754,
     $                   1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766,
     $                   2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572,
     $                   305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893,
     $                   3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307,
     $                   1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297,
     $                   3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966,
     $                   2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758,
     $                   3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598,
     $                   1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406,
     $                   1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922,
     $                   3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038,
     $                   2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934,
     $                   1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091,
     $                   541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451,
     $                   2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580,
     $                   949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958,
     $                   2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055,
     $                   1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507,
     $                   4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078,
     $                   2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273,
     $                   3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17,
     $                   3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854,
     $                   3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916,
     $                   3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971,
     $                   409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889,
     $                   2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831,
     $                   1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621,
     $                   3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541,
     $                   1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893,
     $                   2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736,
     $                   1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992,
     $                   3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787,
     $                   3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125,
     $                   77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364,
     $                   3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460,
     $                   2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257,
     $                   1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574,
     $                   3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912,
     $                   225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216,
     $                   85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248,
     $                   3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401,
     $                   3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124,
     $                   3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762,
     $                   1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149,
     $                   2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245,
     $                   413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166,
     $                   65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466,
     $                   1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018,
     $                   697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399,
     $                   3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190,
     $                   3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879,
     $                   1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153,
     $                   3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320,
     $                   2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18,
     $                   929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712,
     $                   533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159,
     $                   2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318,
     $                   4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091,
     $                   721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443,
     $                   2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510,
     $                   2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449,
     $                   2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956,
     $                   2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201,
     $                   245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137,
     $                   1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399,
     $                   1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321,
     $                   3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271,
     $                   997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667,
     $                   1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703,
     $                   2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629,
     $                   1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365,
     $                   981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431,
     $                   2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113,
     $                   941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922,
     $                   2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554,
     $                   197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184,
     $                   2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099,
     $                   285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228,
     $                   1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012,
     $                   2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921,
     $                   3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452,
     $                   909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901,
     $                   2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572,
     $                   421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309,
     $                   4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171,
     $                   2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817,
     $                   2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039,
     $                   1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696,
     $                   1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256,
     $                   1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715,
     $                   81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077,
     $                   1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019,
     $                   2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497,
     $                   2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101,
     $                   129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717,
     $                   1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51,
     $                   249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981,
     $                   3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978,
     $                   2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813,
     $                   3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881,
     $                   2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76,
     $                   2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846,
     $                   3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694,
     $                   1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682,
     $                   345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124,
     $                   2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660,
     $                   1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997,
     $                   3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479,
     $                   2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141,
     $                   157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886,
     $                   2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514,
     $                   3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301,
     $                   1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604,
     $                   2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888,
     $                   2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836,
     $                   3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990,
     $                   361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058,
     $                   2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692,
     $                   3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194,
     $                   2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20,
     $                   3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285,
     $                   3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046,
     $                   3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107,
     $                   517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508,
     $                   3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525,
     $                   2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801,
     $                   1537 /
*     ..
*     .. Executable Statements ..
*
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
*
      DO 10 I = 1, MIN( N, LV )
*
*        Multiply the seed by i-th power of the multiplier modulo 2**48
*
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) +
     $         I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
*
*        Convert 48-bit integer to a real number in the interval (0,1)
*
         X( I ) = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $            DBLE( IT4 ) ) ) )
   10 CONTINUE
*
*     Return final value of seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
*
*     End of DLARUV_HELPER
*
      END SUBROUTINE DLARUV_HELPER

! ##################################################################################################################################
! 052 LAPACK_BLAS_AUX

      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL(DOUBLE)   CFROM, CTO
      REAL(DOUBLE)   A( LDA, * )
      CALL DLASCL_HELPER( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
      RETURN
*
*     End of DLASCL
*
      END SUBROUTINE DLASCL

! ##################################################################################################################################
! 053 LAPACK_BLAS_AUX

      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      REAL(DOUBLE)   ALPHA, BETA
      REAL(DOUBLE)   A( LDA, * )
      CALL DLASET_HELPER( UPLO, M, N, ALPHA, BETA, A, LDA )
*
      RETURN
*
*     End of DLASET
*
      END SUBROUTINE DLASET

! ##################################################################################################################################
! 054 LAPACK_BLAS_AUX

      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
      REAL(DOUBLE)   A( LDA, * ), C( * ), S( * )

      CALL DLASR_HELPER( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )

      RETURN

      END SUBROUTINE DLASR

! **********************************************************************************************************************************

      SUBROUTINE DLASR_HELPER( SIDE, PIVOT, DIRECT, M, N, C, S, A,
     $                         LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), C( * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n real matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  This version vectorises across rows of the array A when SIDE = 'L'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) REAL(DOUBLE) arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL(DOUBLE)   CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLASR_HELPER
*
      END SUBROUTINE DLASR_HELPER

! ##################################################################################################################################
! 055 LAPACK_BLAS_AUX

      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
      REAL(DOUBLE)   D( * )
      CALL DLASRT_HELPER( ID, N, D, INFO )
      RETURN
*
*     End of DLASRT
*
      END SUBROUTINE DLASRT

! ##################################################################################################################################
! 056 LAPACK_BLAS_AUX

      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL(DOUBLE)   SCALE, SUMSQ
      REAL(DOUBLE)   X( * )
      CALL DLASSQ_HELPER( N, X, INCX, SCALE, SUMSQ )
      RETURN
*
*     End of DLASSQ
*
      END SUBROUTINE DLASSQ

! ##################################################################################################################################
! 057 LAPACK_BLAS_AUX

      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)   A( LDA, * )
      CALL DLASWP_HELPER( N, A, LDA, K1, K2, IPIV, INCX )
*
      RETURN
*
*     End of DLASWP
*
      END SUBROUTINE DLASWP

! ##################################################################################################################################
! 058 LAPACK_BLAS_AUX

      SUBROUTINE DLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,
     $                   SCALE, CNORM, INFO, iter_num, dtbsv_msg )
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      character*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, N
      integer            iter_num
      REAL(DOUBLE)   SCALE
      REAL(DOUBLE)   AB( LDAB, * ), CNORM( * ), X( * )

      CALL DLATBS_HELPER( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB,
     $                    X, SCALE, CNORM, INFO, iter_num, dtbsv_msg )

      RETURN

      END SUBROUTINE DLATBS

! **********************************************************************************************************************************

      SUBROUTINE DLATBS_HELPER( UPLO, TRANS, DIAG, NORMIN, N, KD, AB,
     $                          LDAB, X, SCALE, CNORM, INFO, iter_num,
     $                          dtbsv_msg )
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DLATBS_HELPER'
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      character*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, N
      integer            iter_num
      REAL(DOUBLE)   SCALE
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   AB( LDAB, * ), CNORM( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLATBS solves one of the triangular systems
*
*     A *x = s*b  or  A'*x = s*b
*
*  with scaling to prevent overflow, where A is an upper or lower
*  triangular band matrix.  Here A' denotes the transpose of A, x and b
*  are n-element vectors, and s is a scaling factor, usually less than
*  or equal to 1, chosen so that the components of x will be less than
*  the overflow threshold.  If the unscaled problem will not cause
*  overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A
*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
*  non-trivial solution to A*x = 0 is returned.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  Solve A * x = s*b  (No transpose)
*          = 'T':  Solve A'* x = s*b  (Transpose)
*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  NORMIN  (input) CHARACTER*1
*          Specifies whether CNORM has been set or not.
*          = 'Y':  CNORM contains the column norms on entry
*          = 'N':  CNORM is not set on entry.  On exit, the norms will
*                  be computed and stored in CNORM.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of subdiagonals or superdiagonals in the
*          triangular matrix A.  KD >= 0.
*
*  AB      (input) REAL(DOUBLE) array, dimension (LDAB,N)
*          The upper or lower triangular band matrix A, stored in the
*          first KD+1 rows of the array. The j-th column of A is stored
*          in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  X       (input/output) REAL(DOUBLE) array, dimension (N)
*          On entry, the right hand side b of the triangular system.
*          On exit, X is overwritten by the solution vector x.
*
*  SCALE   (output) REAL(DOUBLE)
*          The scaling factor s for the triangular system
*             A * x = s*b  or  A'* x = s*b.
*          If SCALE = 0, the matrix A is singular or badly scaled, and
*          the vector x is an exact or approximate solution to A*x = 0.
*
*  CNORM   (input or output) REAL(DOUBLE) array, dimension (N)
*
*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*          must be greater than or equal to the 1-norm.
*
*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of A.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
!  iter_num (input) INTEGER
!           This subr may be called iteratively (by DPBCON when
!           estimating RCOND (cond number of AB)

!  dtbsv_msg (input) CHARACTER
!           = 'Y', have subr DTBSV print Fwd, Back pass messages

*  Further Details
*  ======= =======
*
*  A rough bound on x is computed; if that is less than overflow, DTBSV
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*  if A is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / A(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*       end
*
*  Define bounds on the components of x after j iterations of the loop:
*     M(j) = bound on x[1:j]
*     G(j) = bound on x[j+1:n]
*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*
*  Then for iteration j+1 we have
*     M(j+1) <= G(j) / | A(j+1,j+1) |
*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*
*  where CNORM(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of A, not counting the diagonal.  Hence
*
*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*                                   1<=i< j
*
*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the
*  reciprocal of the largest M(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  The bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  If
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*
*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
*  algorithm for A upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*       end
*
*  We simultaneously compute two bounds
*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       M(j) = bound on x(i), 1<=i<=j
*
*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*  Then the bound on x(j) is
*
*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*
*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*                      1<=i<=j
*
*  and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater
*  than max(underflow, 1/overflow).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND
      REAL(DOUBLE)   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS,
     $                   TMAX, TSCAL, USCAL, XBND, XJ, XMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME

      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH

*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*


! **********************************************************************************************************************************
*      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( KD.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATBS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            DO 10 J = 1, N
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = DASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N
               JLEN = MIN( KD, N-J )
               IF( JLEN.GT.0 ) THEN
                  CNORM( J ) = DASUM( JLEN, AB( 2, J ), 1 )
               ELSE
                  CNORM( J ) = ZERO
               END IF
   20       CONTINUE
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM.
*
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine DTBSV can be used.
*
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              M(j) = G(j-1) / abs(A(j,j))
*
               TJJ = ABS( AB( MAIND, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A' * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
               TJJ = ABS( AB( MAIND, J ) )
               IF( XJ.GT.TJJ )
     $            XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL DTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1
     &             , dtbsv_msg )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 110 J = JFIRST, JLAST, JINC

               write(sc1,12345,advance='no') iter_num,j,jlast,cr13_lba
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = AB( MAIND, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE )
     $               GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
*                                             x(j)* A(max(1,j-kd):j-1,j)
*
                     JLEN = MIN( KD, J-1 )
                     CALL DAXPY( JLEN, -X( J )*TSCAL,
     $                           AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE IF( J.LT.N ) THEN
*
*                 Compute the update
*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
*                                          x(j) * A(j+1:min(j+kd,n),j)
*
                  JLEN = MIN( KD, N-J )
                  IF( JLEN.GT.0 )
     $               CALL DAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1,
     $                           X( J+1 ), 1 )
                  I = J + IDAMAX( N-J, X( J+1 ), 1 )
                  XMAX = ABS( X( I ) )
               END IF
  110       CONTINUE
*
         ELSE
*
*           Solve A' * x = b
*
            DO 160 J = JFIRST, JLAST, JINC

               write(sc1,12345,advance='no') iter_num,j,jlast,cr13_lba
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call DDOT to perform the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     SUMJ = DDOT( JLEN, AB( KD+1-JLEN, J ), 1,
     $                      X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     IF( JLEN.GT.0 )
     $                  SUMJ = DDOT( JLEN, AB( 2, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     DO 120 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( KD+I-JLEN, J )*USCAL )*
     $                         X( J-JLEN-1+I )
  120                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 130 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
  130                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.TSCAL ) THEN
*
*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 150
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to A'*x = 0.
*
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - sumj if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
*
*     End of DLATBS_HELPER
*
12345 format(5X,'Iteration number ',i4,' : J = ',i8,' to ',i8, a)


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE DLATBS_HELPER

! ##################################################################################################################################
! 059 LAPACK_BLAS_AUX

      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
      REAL(DOUBLE)   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )

      CALL DLATRD_HELPER( UPLO, N, NB, A, LDA, E, TAU, W, LDW )

      RETURN

      END SUBROUTINE DLATRD

! **********************************************************************************************************************************

      SUBROUTINE DLATRD_HELPER( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  DLATRD reduces NB rows and columns of a real symmetric matrix A to
*  symmetric tridiagonal form by an orthogonal similarity
*  transformation Q' * A * Q, and returns the matrices V and W which are
*  needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by DSYTRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements above the diagonal
*            with the array TAU, represent the orthogonal matrix Q as a
*            product of elementary reflectors;
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            with the array TAU, represent the  orthogonal matrix Q as a
*            product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= (1,N).
*
*  E       (output) REAL(DOUBLE) array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) REAL(DOUBLE) array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) REAL(DOUBLE) array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*  and tau in TAU(i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*  and tau in TAU(i).
*
*  The elements of the vectors v together form the n-by-nb matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a symmetric rank-2k update of the form:
*  A := A - V*W' - W*V'.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  a   a   a   v4  v5 )              (  d                  )
*    (      a   a   v4  v5 )              (  1   d              )
*    (          a   1   v5 )              (  v1  1   a          )
*    (              d   1  )              (  v1  v2  a   a      )
*    (                  d  )              (  v1  v2  a   a   a  )
*
*  where d denotes a diagonal element of the reduced matrix, a denotes
*  an element of the original matrix that is unchanged, and vi denotes
*  an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      REAL(DOUBLE)   ALPHA
*     ..
*     .. External Subroutines ..
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DLATRD_HELPER
*
      END SUBROUTINE DLATRD_HELPER

! ##################################################################################################################################
! 060 LAPACK_BLAS_AUX

      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
      INTEGER            INFO, K, LDA, M, N
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )

      CALL DORG2L_HELPER( M, N, K, A, LDA, TAU, WORK, INFO )

      RETURN

      END SUBROUTINE DORG2L

! **********************************************************************************************************************************

      SUBROUTINE DORG2L_HELPER( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORG2L generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the last n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by DGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQLF in the last k columns of its array
*          argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL(DOUBLE) array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQLF.
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2L_HELPER
*
      END SUBROUTINE DORG2L_HELPER

! ##################################################################################################################################
! 061 LAPACK_BLAS_AUX

      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
      INTEGER            INFO, K, LDA, M, N
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )

      CALL DORG2R_HELPER( M, N, K, A, LDA, TAU, WORK, INFO )

      RETURN

      END SUBROUTINE DORG2R

! **********************************************************************************************************************************

      SUBROUTINE DORG2R_HELPER( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORG2R generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL(DOUBLE) array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2R_HELPER
*
      END SUBROUTINE DORG2R_HELPER

! ##################################################################################################################################
! 062 LAPACK_BLAS_AUX

      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )

      CALL DORGQL_HELPER( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

      RETURN

      END SUBROUTINE DORGQL

! **********************************************************************************************************************************

      SUBROUTINE DORGQL_HELPER( M, N, K, A, LDA, TAU, WORK, LWORK,
     $                          INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGQL generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the last N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by DGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQLF in the last k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL(DOUBLE) array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQLF.
*
*  WORK    (workspace/output) REAL(DOUBLE) array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQL_HELPER
*
      END SUBROUTINE DORGQL_HELPER

! ##################################################################################################################################
! 063 LAPACK_BLAS_AUX

      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )

      CALL DORGQR_HELPER( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

      RETURN

      END SUBROUTINE DORGQR

! **********************************************************************************************************************************

      SUBROUTINE DORGQR_HELPER( M, N, K, A, LDA, TAU, WORK, LWORK,
     $                          INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL(DOUBLE) array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace/output) REAL(DOUBLE) array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR_HELPER
*
      END SUBROUTINE DORGQR_HELPER

! ##################################################################################################################################
! 064 LAPACK_BLAS_AUX

      SUBROUTINE DRSCL( N, SA, SX, INCX )
      INTEGER            INCX, N
      REAL(DOUBLE)   SA
      REAL(DOUBLE)   SX( * )

      CALL DRSCL_HELPER( N, SA, SX, INCX )

      RETURN

      END SUBROUTINE DRSCL

! **********************************************************************************************************************************

      SUBROUTINE DRSCL_HELPER( N, SA, SX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL(DOUBLE)   SA
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   SX( * )
*     ..
*
*  Purpose
*  =======
*
*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.
*  This is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of components of the vector x.
*
*  SA      (input) REAL(DOUBLE)
*          The scalar a which is used to divide each component of x.
*          SA must be >= 0, or the subroutine will divide by zero.
*
*  SX      (input/output) REAL(DOUBLE) array, dimension
*                         (1+(N-1)*abs(INCX))
*          The n-element vector x.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector SX.
*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      REAL(DOUBLE)   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Functions ..
      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply X by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector X by MUL
*
      CALL DSCAL( N, MUL, SX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DRSCL_HELPER
*
      END SUBROUTINE DRSCL_HELPER

! ##################################################################################################################################
! 065 LAPACK_BLAS_AUX

      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      REAL(DOUBLE)   A( LDA, * ), D( * ), E( * ), TAU( * )

      CALL DSYTD2_HELPER( UPLO, N, A, LDA, D, E, TAU, INFO )

      RETURN

      END SUBROUTINE DSYTD2

! **********************************************************************************************************************************

      SUBROUTINE DSYTD2_HELPER( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   A( LDA, * ), D( * ), E( * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
*  form T by an orthogonal similarity transformation: Q' * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) REAL(DOUBLE) array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) REAL(DOUBLE) array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) REAL(DOUBLE) array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      REAL(DOUBLE)   ALPHA, TAUI
*     ..
*     .. External Subroutines ..
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of DSYTD2_HELPER
*
      END SUBROUTINE DSYTD2_HELPER

! #################################################################################################################################
! 069 LAPACK_BLAS_AUX

      LOGICAL FUNCTION DISNAN(DIN)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN      (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
      DISNAN = DISNAN_HELPER(DIN)
      RETURN
      END FUNCTION DISNAN

! #################################################################################################################################
! 070 LAPACK_BLAS_AUX

      LOGICAL FUNCTION DLAISNAN(DIN1,DIN2)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN1,DIN2
*     ..
*
*  Purpose
*  =======
*
*  This routine is not for general use.  It exists solely to avoid
*  over-optimization in DISNAN.
*
*  DLAISNAN checks for NaNs by comparing its two arguments for
*  inequality.  NaN is the only floating-point value where NaN != NaN
*  returns .TRUE.  To check for NaNs, pass the same variable as both
*  arguments.
*
*  Strictly speaking, Fortran does not allow aliasing of function
*  arguments. So a compiler must assume that the two arguments are
*  not the same variable, and the test will not be optimized away.
*  Interprocedural or whole-program optimization may delete this
*  test.  The ISNAN functions will be replaced by the correct
*  Fortran 03 intrinsic once the intrinsic is widely available.
*
*  Arguments
*  =========
*
*  DIN1     (input) DOUBLE PRECISION
*  DIN2     (input) DOUBLE PRECISION
*          Two numbers to compare for inequality.
*
*  =====================================================================
*
      DLAISNAN = DLAISNAN_HELPER(DIN1,DIN2)
      RETURN
      END FUNCTION DLAISNAN

! --- lapack_surgery end --- !
      END MODULE LAPACK_BLAS_AUX
