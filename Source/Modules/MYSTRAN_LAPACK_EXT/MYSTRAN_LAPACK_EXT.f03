! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

! ##################################################################################################################################
!
! MYSTRAN_LAPACK_EXT
!
! Salvaged MYSTRAN-specific LAPACK extensions. These routines have
! mathematical deviations from the corresponding upstream LAPACK
! reference implementations and must therefore stay compiled into
! MYSTRAN regardless of which BLAS+LAPACK provider is being linked.
!
! Every other routine that used to live under Source/Modules/LAPACK has
! been removed; consumers now call standard LAPACK directly (resolved
! at link time against the chosen provider).
!
! Members:
!   DPTTRF_MYSTRAN              Modified DPTTRF: tests D(I) == 0 instead
!                               of D(I) <= 0, so the L*D*L^T factorisation
!                               is produced even with negative diagonal
!                               entries. Required by EST_NUMBER_OF_EIGENS
!                               (LINK4) which counts negative diagonals
!                               to bound the eigenvalue spectrum.
!
!   DSBGVX_GIV_MGIV             Renamed DSBGVX. Adds three output
!                               arguments (mlam / eig_num / mvec) and
!                               accepts a `method` selector to drive
!                               the GIV / MGIV eigensolvers. Calls
!                               DSTEBZ_MYSTRAN below.
!
!   DSTEBZ_MYSTRAN              Renamed DSTEBZ. Produces two extra
!                               output integers (lowest_mode_num,
!                               highest_mode_num) needed by
!                               DSBGVX_GIV_MGIV to populate eig_num.
!                               No other algorithmic deviation.
!
!   EIGENVALUE_CONVERGENCE_FAILURE
!                               MYSTRAN-specific helper used by
!                               DSBGVX_GIV_MGIV to surface DSTEBZ
!                               convergence failures through MYSTRAN's
!                               error-reporting infrastructure.
!
! ##################################################################################################################################

      MODULE MYSTRAN_LAPACK_EXT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, SOL_NAME

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  LINK_MESSAGE

      IMPLICIT NONE

      CONTAINS

! ##################################################################################################################################
! DPTTRF_MYSTRAN

      SUBROUTINE DPTTRF_MYSTRAN( N, D, E, INFO )

!  This is a MYSTRAN-specific modification of LAPACK subroutine DPTTRF.
!  The original singularity test
!         IF( D( I ).LE.ZERO ) THEN
!  has been replaced with
!         IF( D( I ).EQ.ZERO ) THEN
!  so that the L*D*L^T decomposition is produced even when diagonal
!  elements are negative. EST_NUMBER_OF_EIGENS in LINK4 needs the
!  number of negative diagonals after factorisation to estimate the
!  number of eigenvalues below a given shift.

!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      REAL(DOUBLE)   D( * ), E( * )
!     ..
!     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I4
      REAL(DOUBLE)   EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DPTTRF', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 ) RETURN

      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE

      DO 20 I = I4 + 1, N - 4, 4
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI

         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 1
            GO TO 30
         END IF
         EI = E( I+1 )
         E( I+1 ) = EI / D( I+1 )
         D( I+2 ) = D( I+2 ) - E( I+1 )*EI

         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 2
            GO TO 30
         END IF
         EI = E( I+2 )
         E( I+2 ) = EI / D( I+2 )
         D( I+3 ) = D( I+3 ) - E( I+2 )*EI

         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 3
            GO TO 30
         END IF
         EI = E( I+3 )
         E( I+3 ) = EI / D( I+3 )
         D( I+4 ) = D( I+4 ) - E( I+3 )*EI
   20 CONTINUE

      IF( D( N ).LE.ZERO ) INFO = N

   30 CONTINUE
      RETURN

      END SUBROUTINE DPTTRF_MYSTRAN

! ##################################################################################################################################
! DSBGVX_GIV_MGIV
!
! MYSTRAN-specific driver: thin wrapper around the same algorithm as
! upstream LAPACK DSBGVX, but augmented with three extra outputs
! (mlam / eig_num / mvec) and a `method` selector. The body still
! drives the same chain of LAPACK helpers (DPBSTF, DSBGST, DSBTRD,
! DSTERF, DSTEQR, DSTEIN, DCOPY, DGEMV, DSWAP, DLACPY) plus the
! MYSTRAN-specific DSTEBZ_MYSTRAN below.

      SUBROUTINE DSBGVX_GIV_MGIV ( JOBZ, RANGE, UPLO, N, KA, KB, AB,    &
     &                             LDAB, BB, LDBB, Q, LDQ, VL, VU,      &
     &                             IL, IU, ABSTOL, mlam, W, Z, LDZ,     &
     &                             WORK, IWORK, IFAIL, INFO,            &
     &                             method, eig_num, mvec )

!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      character(LEN=8)   method
      INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ,    &
     &                   mlam, n, eig_num(n), mvec
      REAL(DOUBLE)   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL(DOUBLE)   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),         &
     &                   W( * ), WORK( * ), Z( LDZ, * )
!     ..
!     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, UPPER, VALEIG, WANTZ
      CHARACTER          ORDER, VECT
      INTEGER            I, IINFO, INDD, INDE, INDEE, INDIBL, INDISP,   &
     &                   INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT
      integer            lowest_mode_num, highest_mode_num
      REAL(DOUBLE)   TMP1
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, DPBSTF, DSBGST, DSBTRD, DSTERF, DSTEQR,&
     &                   DLACPY, DCOPY, DSTEIN, DGEMV, DSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..

! Initialize eig_num
      do i=1,n
         eig_num(i) = 0
      enddo

!     Test the input parameters.
      WANTZ  = LSAME( JOBZ,  'V' )
      UPPER  = LSAME( UPLO,  'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KA.LT.0 ) THEN
         INFO = -5
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -8
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 ) THEN
         INFO = -12
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -14
      ELSE IF( INDEIG .AND. IL.LT.0 ) THEN
         INFO = -15
      ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) ) ) THEN
         INFO = -16
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -21
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBGVX_GIV_MGIV', -INFO )
      END IF

!     Form a split Cholesky factorization of B.
      if (method(1:3) == 'GIV') then
         if (sol_name(1:8) == 'BUCKLING') then
            CALL LINK_MESSAGE(                                          &
     &       '  CHOLESKY FACTORIZATION OF DIFFER STIFF MATRIX')
         else
            CALL LINK_MESSAGE(                                          &
     &       '  CHOLESKY FACTORIZATION OF MASS MATRIX')
         endif
      else if (method(1:4) == 'MGIV') then
         CALL LINK_MESSAGE('  CHOLESKY FACTORIZATION OF STIFF MATRIX')
      endif
      CALL DPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF

!     Transform problem to standard eigenvalue problem.
      CALL LINK_MESSAGE('  TRANSFORM TO STANDARD EIGENVALUE PROBLEM')
      CALL DSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ,   &
     &             WORK, IINFO )

!     Reduce symmetric band matrix to tridiagonal form.
      CALL LINK_MESSAGE('  REDUCE SYMM BAND MATRIX TO TRIDIAG FORM')
      INDD = 1
      INDE = INDD + N
      INDWRK = INDE + N
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL DSBTRD( VECT, UPLO, N, KA, AB, LDAB, WORK( INDD ),           &
     &             WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO )

!     If all eigenvalues are desired and ABSTOL <= 0, try DSTERF/DSTEQR
!     first; fall back to DSTEBZ_MYSTRAN+DSTEIN if that fails.
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.ge.N ) ) .AND.  &
     &    ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
         IF( .NOT.WANTZ ) THEN
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,             &
     &                   WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            END IF
         END IF

         IF( INFO.EQ.0 ) THEN
            mlam = N
            mvec = N
            GO TO 30
         END IF
         INFO = 0
      END IF

!     Otherwise, call DSTEBZ_MYSTRAN and, if eigenvectors are desired,
!     call DSTEIN.
      if (info > 0) then
         Write(err,9901) info
         Write(f06,9901) info
      endif
 9901 format(' *INFORMATION: LAPACK SUBR DSTERF OR DSTEQR HAS FAILED TO &
     &FIND ALL OF THE EIGENVALUES IN A TOTAL OF 30*NDOFA ITERATIONS'    &
     &,/,14X,' A TOTAL OF ',I8,' SUB-DIAGONAL ELEMENTS OF THE TRIDIAGONA&
     &L MATRIX E HAVE NOT CONVERGED TO ZERO.'                           &
     &,/,14X,' LAPACK WILL ATTEMPT TO USE SUBR DSTEBZ TO FIND THE EIGENV&
     &ALUES.',/)

      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ_MYSTRAN( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL,     &
     &             WORK( INDD ), WORK( INDE ), mlam, NSPLIT, W,         &
     &             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),    &
     &             IWORK( INDIWO ), INFO,                               &
     &             lowest_mode_num, highest_mode_num )

      if (info > 0) then
         call eigenvalue_convergence_failure ( range, info )
         if ((info == 1) .or. (info == 3) .and. (range /= 'I')) then
            Write(f06,9903)
            do i=1,mlam
               if (iwork(indibl-1+i) < 0) then
                  Write(f06,9904) i,w(i)
               endif
            enddo
            Write(f06,*)
         endif
      endif
 9903 format(15x,'THE EIGENVALUES OF QUESTIONABLE VALUE ARE',/,         &
     &       15X,'   INDEX     EIGENVALUE')
 9904 FORMAT(15X,I8,1ES15.6)

      do i=1,mlam
         if (method(1:3) == 'GIV') then
            eig_num(i) =  lowest_mode_num + (i - 1)
         else
            eig_num(i) = (n + 1) - (lowest_mode_num + (i - 1))
         endif
      enddo

      IF( WANTZ ) THEN
         mvec = mlam
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), mvec, W,           &
     &                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,         &
     &                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )

!        Apply transformation matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
         DO 20 J = 1, mvec
            CALL DCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL DGEMV( 'N', N, N, ONE, Q, LDQ, WORK, 1, ZERO,          &
     &                  Z( 1, J ), 1 )
   20    CONTINUE
      END IF

   30 CONTINUE

!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
      IF( WANTZ ) THEN
         DO 50 J = 1, mvec - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, mvec
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   40       CONTINUE

            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   50    CONTINUE
      END IF

      RETURN

      END SUBROUTINE DSBGVX_GIV_MGIV

! ##################################################################################################################################
! EIGENVALUE_CONVERGENCE_FAILURE

      SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE ( RANGE, INFO )

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

 9902 format(' *INFORMATION: SOME OR ALL OF THE EIGENVALUES FAILED TO CO&
     &NVERGE OR WERE NOT COMPUTED IN LAPACK SUBROUTINE DSTEBZ:')

99021 format(15x,'BISECTION FAILED TO CONVERGE FOR SOME EIGENVALUES; THE&
     &SE EIGENVALUES ARE FLAGGED BY A NEGATIVE BLOCK NUMBER.',/,15X,    &
     &'THE EFFECT IS THAT THE EIGENVALUES MAY NOT BE AS ACCURATE AS THE &
     &ABSOLUTE AND RELATIVE TOLERANCES.',/,15X,                         &
     &'THIS IS GENERALLY CAUSED BY UNEXPECTEDLY INACCURATE ARITHMETIC.' &
     &,/)

99022 format(15x,'NOT ALL OF THE EIGENVALUES IN THE RANGE REQUESTED WERE&
     & FOUND:',/,15X,                                                   &
     &'CAUSE: NON-MONOTONIC ARITHMETIC, CAUSING THE STURM SEQUENCE TO BE&
     & NON-MONOTONIC.',/,15X,                                           &
     &'CURE : RECALCULATE, REQUESTING ALL EIGENVALUES',/)

  803 format(' *ERROR   803: PROGRAMMING ERROR IN SUBROUTINE DSTEBZ.'   &
     &,/,15X,'NO EIGENVALUES WERE COMPUTED BY LAPACK SUBROUTINE DSTEBZ. &
     &THE GERSHGORIN INTERVAL INITIALLY USED WAS TOO SMALL.',/,15X,     &
     &'PROBABLE CAUSE: YOUR MACHINE HAS SLOPPY FLOATING-POINT ARITHMETIC&
     &',/,15X,'CURE          : INCREASE THE PARAMETER "FUDGE" IN LAPACK &
     &SUBROUTINE DSTEBZ, RECOMPILE, AND TRY AGAIN',/)

      END SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE

! ##################################################################################################################################
! DSTEBZ_MYSTRAN
!
! Renamed copy of LAPACK DSTEBZ. The only deviation from upstream is
! that two extra outputs (lowest_mode_num, highest_mode_num) are
! returned at the end of the argument list. The numerical algorithm is
! unchanged.

      SUBROUTINE DSTEBZ_MYSTRAN( RANGE, ORDER, N, VL, VU, IL, IU,       &
     &                   ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT,    &
     &                   WORK, IWORK, INFO,                             &
     &                   lowest_mode_num, highest_mode_num )

!     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      integer            lowest_mode_num, highest_mode_num
      REAL(DOUBLE)   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      REAL(DOUBLE)   D( * ), E( * ), W( * ), WORK( * )
!     ..
!     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,      &
     &                   HALF = 1.0D0 / TWO )
      REAL(DOUBLE)   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO,   &
     &                   IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX,     &
     &                   ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL,   &
     &                   NWU
      REAL(DOUBLE)   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN,       &
     &                   TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
!     ..
!     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME

      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH

      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAEBZ, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..

      INFO = 0

!     Decode RANGE
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF

!     Decode ORDER
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF

!     Check for Errors
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( IORDER.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.2 ) THEN
         IF( VL.GE.VU )                                                 &
     &      INFO = -5
      ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )   &
     &          THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )  &
     &          THEN
         INFO = -7
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEBZ', -INFO )
         RETURN
      END IF

      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.

      M = 0
      IF( N.EQ.0 ) RETURN

      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N )                     &
     &   IRANGE = 1

      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 ) NB = 0

!     Special Case when N=1
      IF( N.EQ.1 ) THEN
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         IF( IRANGE.EQ.2 .AND. ( VL.GE.D( 1 ) .OR. VU.LT.D( 1 ) ) ) THEN
            M = 0
         ELSE
            W( 1 ) = D( 1 )
            IBLOCK( 1 ) = 1
            M = 1
         END IF
         highest_mode_num = 1
         lowest_mode_num  = 1
         RETURN
      END IF

!     Compute Splitting Points
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE

!DIR$ NOVECTOR
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         IF( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            ISPLIT( NSPLIT ) = J - 1
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN

      IF( IRANGE.EQ.3 ) THEN
         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO

         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
   20    CONTINUE

         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN

         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /           &
     &           LOG( TWO ) ) + 2
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF

         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU

         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,&
     &                WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, &
     &                IWORK, W, IBLOCK, IINFO )

         IF( IWORK( 6 ).EQ.IU ) THEN
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         END IF

         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),                      &
     &           ABS( D( N ) )+ABS( E( N-1 ) ) )

         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+          &
     &              ABS( E( J ) ) )
   30    CONTINUE

         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF

         IF( IRANGE.EQ.2 ) THEN
            WL = VL
            WU = VU
         ELSE
            WL = ZERO
            WU = ZERO
         END IF
      END IF

      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0

      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF

         IF( IN.EQ.1 ) THEN
            IF( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN )             &
     &         NWL = NWL + 1
            IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN )             &
     &         NWU = NWU + 1
            IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE.&
     &          D( IBEGIN )-PIVMIN ) ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            END IF
         ELSE
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO

            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE

            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN

            IF( ABSTOL.LE.ZERO ) THEN
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            ELSE
               ATOLI = ABSTOL
            END IF

            IF( IRANGE.GT.1 ) THEN
               IF( GU.LT.WL ) THEN
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               END IF
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU )                                           &
     &            GO TO 70
            END IF

            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,     &
     &                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),      &
     &                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,     &
     &                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )

            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )

            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) /        &
     &              LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
     &                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),      &
     &                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,   &
     &                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )

            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )

               IF( J.GT.IOUT-IINFO ) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF,                       &
     &                 IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE

            M = M + IM
         END IF
   70 CONTINUE

      IF( IRANGE.EQ.3 ) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU

         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
            DO 80 JE = 1, M
               IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
                  IDISCL = IDISCL - 1
               ELSE IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
   80       CONTINUE
            M = IM
         END IF
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.                        &
     &                   ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.                        &
     &                   ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
  110             CONTINUE
                  IBLOCK( IW ) = 0
  120          CONTINUE
            END IF
            IM = 0
            DO 130 JE = 1, M
               IF( IBLOCK( JE ).NE.0 ) THEN
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
  130       CONTINUE
            M = IM
         END IF
         IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
            TOOFEW = .TRUE.
         END IF
      END IF

      IF( IORDER.EQ.1 .AND. NSPLIT.GT.1 ) THEN
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               IF( W( J ).LT.TMP1 ) THEN
                  IE = J
                  TMP1 = W( J )
               END IF
  140       CONTINUE

            IF( IE.NE.0 ) THEN
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            END IF
  150    CONTINUE
      END IF

      highest_mode_num = nwu
      lowest_mode_num = highest_mode_num - m + 1

      INFO = 0
      IF( NCNVRG ) INFO = INFO + 1
      IF( TOOFEW )  INFO = INFO + 2

      RETURN

      END SUBROUTINE DSTEBZ_MYSTRAN

      END MODULE MYSTRAN_LAPACK_EXT
