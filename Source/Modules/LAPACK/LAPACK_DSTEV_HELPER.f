! ##################################################################################################################################

      MODULE LAPACK_DSTEV_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE

      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )

      CHARACTER          JOBZ
      INTEGER            INFO, LDZ, N
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      LOGICAL            WANTZ
      INTEGER            IMAX, ISCALE
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM,
     $                   TNRM
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
      EXTERNAL           DSCAL, DSTEQR, DSTERF, XERBLA
      INTRINSIC          SQRT

      WANTZ = LSAME( JOBZ, 'V' )

      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -6
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEV ', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 )
     $   RETURN

      IF( N.EQ.1 ) THEN
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )

      ISCALE = 0
      TNRM = DLANST( 'M', N, D, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, SIGMA, D, 1 )
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
      END IF

      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, D, E, INFO )
      ELSE
         CALL DSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
      END IF

      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, D, 1 )
      END IF

      RETURN
      END SUBROUTINE DSTEV
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DSTEV_HELPER
