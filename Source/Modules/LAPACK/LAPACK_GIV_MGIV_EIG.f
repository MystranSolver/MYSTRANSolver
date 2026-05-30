! ##################################################################################################################################

      MODULE LAPACK_GIV_MGIV_EIG
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LINKNO
      USE TIMDAT, ONLY                :  TSEC
      USE LAPACK_BLAS_AUX
      USE LAPACK_MISCEL                                    ! This contains DSTEQR, used in this module
      USE LAPACK_GIV_MGIV_EIG_HELPER
      USE LAPACK_DPBSTF_HELPER

      USE OURTIM_Interface
      USE OUTA_HERE_Interface
      USE LINK_MESSAGE_Interface

      character(1*byte), parameter      :: cr13_lge = char(13)

! This is a set of LAPACK routines for solving for all, or some, eigenvalues and, possibly, some eigenvectors of:

!                             Ax = (Lambda)Bx             (1)

! where A and B are real and symmetric.

! This module contains the major LAPACK subroutines described below (and contained in this module):

!           DSBGVX_GIV_MGIV : Main driver for solving (1) for eigenvalues and eigenvectors, This driver calls:

!           DPBSTF: to form split Cholesky factorization of B.

!           DSBGST: to transform problem from from Ax = (Lambda)Bx to std eigen problem Cy = (Lambda)y. This subr calls:

!           DSBTRD: to reduce the symmetric band matrix C to a symmetric tridiagonal form by orthogonal similarity transforms.

! In the case where all eigenvalues and, possibly all eigenvecs are requested, the driver calls 1 or all of the following subr's:

!           DSTERF: to compute all eigenvalues of the tridiagonal matrix.
!                   NOTE: DSTERF is not in this module, it is in module LAPACK_MISCEL, since it is used in several modules
! or
!           DSTEQR: to compute all eigenvalues and all eigenvectors of the tridiagonal matrix (if eigenvectors are requested).
!                   NOTE: DSTEQR is not in this module, it is in module LAPACK_MISCEL, since it is used in several modules

! or, in the case where some eigenvalues and, possibly some eigenvectors are requested, the driver calls one or both of the subr's:

!           DSTEBZ: to compute some eigenvalues of the tridiagonal matrix

! then, if vectors are requested:

!           DSTEIN: to compute the eigenvectors corresponding to the eigenvalues computed by subr DSTEBZ. This subr calls:

!                   DLAGTF: to factorize the matrix (T - Lambds*I) where T is the tridiagonal matrix

! In addition, other LAPACK (BLAS) procedures are called by the above and are contained in module LAPACK_BLAS_AUX_1

      CONTAINS

! ##################################################################################################################################
! 001 LAPACK_FIV_MGIV_EIG

      SUBROUTINE DSBGVX_GIV_MGIV ( JOBZ, RANGE, UPLO, N, KA, KB, AB,
     &                             LDAB, BB, LDBB, Q, LDQ, VL, VU,
     &                             IL, IU, ABSTOL, mlam, W, Z, LDZ,
     &                             WORK, IWORK, IFAIL, INFO,
     &                             method, eig_num, mvec )
*
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE PARAMS, ONLY                :  SUPINFO, NOCOUNTS
      USE SCONTR, ONLY                :  SOL_NAME

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DSBGVX_GIV_MGIV'

*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      character*8        method
      INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ,
     $                   mlam, n, eig_num(n), mvec
      REAL(DOUBLE)   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL(DOUBLE)   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),
     $                   W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSBGVX computes selected eigenvalues, and optionally, eigenvectors
*  of a real generalized symmetric-definite banded eigenproblem, of
*  the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
*  and banded, and B is also positive definite.  Eigenvalues and
*  eigenvectors can be selected by specifying either all eigenvalues,
*  a range of values or a range of indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  KA      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
*
*  KB      (input) INTEGER
*          The number of superdiagonals of the matrix B if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
*
*  AB      (input/output) REAL(DOUBLE) array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first ka+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
*
*          On exit, the contents of AB are destroyed.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KA+1.
*
*  BB      (input/output) REAL(DOUBLE) array, dimension (LDBB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix B, stored in the first kb+1 rows of the array.  The
*          j-th column of B is stored in the j-th column of the array BB
*          as follows:
*          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
*
*          On exit, the factor S from the split Cholesky factorization
*          B = S**T*S, as returned by DPBSTF.
*
*  LDBB    (input) INTEGER
*          The leading dimension of the array BB.  LDBB >= KB+1.
*
*  Q       (output) REAL(DOUBLE) array, dimension (LDQ, N)
*          If JOBZ = 'V', the n-by-n matrix used in the reduction of
*          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,
*          and consequently C to tridiagonal form.
*          If JOBZ = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  If JOBZ = 'N',
*          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
*
*  VL      (input) REAL(DOUBLE)
*  VU      (input) REAL(DOUBLE)
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          0 <= IL <= IU; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'V'.
*
*  ABSTOL  (input) REAL(DOUBLE)
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*DLAMCH('S').
*
*  mlam    (output) INTEGER
*          The total number of eigenvalues found.  0 <= mlam <= N.
*          If RANGE = 'I', mlam = IU-IL+1.
*
*  W       (output) REAL(DOUBLE) array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) REAL(DOUBLE) array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
*          eigenvectors, with the i-th column of Z holding the
*          eigenvector associated with W(i).  The eigenvectors are
*          normalized so Z**T*B*Z = I.
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) REAL(DOUBLE) array, dimension (7N)
*
*  IWORK   (workspace/output) INTEGER array, dimension (5N)
*
*  IFAIL   (input) INTEGER array, dimension (mlam)
*          If JOBZ = 'V', then if INFO = 0, the first mlam elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvalues that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0 : successful exit
*          < 0 : if INFO = -i, the i-th argument had an illegal value
*          <= N: if INFO = i, then i eigenvectors failed to converge.
*                  Their indices are stored in IFAIL.
*          > N : DPBSTF returned an error code; i.e.,
*                if INFO = N + i, for 1 <= i <= N, then the leading
*                minor of order i of B is not positive definite.
*                The factorization of B could not be completed and
*                no eigenvalues or eigenvectors were computed.

!  method  (input) character*8 = 'GIV' or 'MGIV'

!  mvec    (output) INTEGER, number of eigenvectors found

!  eig_num (output) Integer array (dim N) of mode numbers

*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, UPPER, VALEIG, WANTZ
      CHARACTER          ORDER, VECT
      character*6        sub_nam
      INTEGER            I, IINFO, INDD, INDE, INDEE, INDIBL, INDISP,
     $                   INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT
      integer            lowest_mode_num, highest_mode_num, pm, qm, m1
      REAL(DOUBLE)   TMP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..

! Initialize eig_num

      do i=1,n
         eig_num(i) = 0
      enddo
*
*     Test the input parameters.
*
      WANTZ  = LSAME( JOBZ,  'V' )
      UPPER  = LSAME( UPLO,  'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
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
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBGVX_GIV_MGIV', -INFO )
      END IF
*
*     Form a split Cholesky factorization of B.
*
      if (method(1:3) == 'GIV') then
         if (sol_name(1:8) == 'BUCKLING') then
            CALL LINK_MESSAGE(
     $       '  CHOLESKY FACTORIZATION OF DIFFER STIFF MATRIX')
         else
            CALL LINK_MESSAGE(
     $       '  CHOLESKY FACTORIZATION OF MASS MATRIX')
         endif
      else if (method(1:4) == 'MGIV') then
         CALL LINK_MESSAGE('  CHOLESKY FACTORIZATION OF STIFF MATRIX')
      endif
      CALL DPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem.
*
      CALL LINK_MESSAGE('  TRANSFORM TO STANDARD EIGENVALUE PROBLEM')
      CALL DSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ,
     $             WORK, IINFO )
*
*     Reduce symmetric band matrix to tridiagonal form.
*
      CALL LINK_MESSAGE('  REDUCE SYMM BAND MATRIX TO TRIDIAG FORM')
      INDD = 1
      INDE = INDD + N
      INDWRK = INDE + N
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL DSBTRD( VECT, UPLO, N, KA, AB, LDAB, WORK( INDD ),
     $             WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO )
*
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then call:
!        DSTERF (eigenvalues only) or
!        DSTEQR (eigenvalues and eigenvectors).
!     If this fails for some eigenvalue, then try DSTEBZ.
*
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.ge.N ) ) .AND.
     $    ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
         IF( .NOT.WANTZ ) THEN
            sub_nam = 'DSTERF'
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
            sub_nam = 'DSTEQR'
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,
     $                   WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            END IF
         END IF

         do i=1,n
            eig_num = i
         enddo

         IF( INFO.EQ.0 ) THEN
            mlam = N
            mvec = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call DSTEBZ and, if eigenvectors are desired,
*     call DSTEIN.
! /////////////////////////////////////////////////////////////////////B
* If we got here, then 1 of 2 things happened:
!  (a) all eigenvalues were requested and ABSTOL <= 0, and subrs DSTERF
!      (or DSTEQR), above, failed to find all of them (INFO > 0).
!      If INFO > 0, LAPACK will attempt to use subrs DSTEBZ and DSTEIN
!      to find the eigenvalues/vectors
!  (b) all eigenvalues were not requested and subrs DSTEBZ and DSTEIN
!      will be called to find them
! At any rate, print message if INFO > 0

      if (info > 0) then
         Write(err,9901) info
         if (supinfo == 'N') then
            Write(f06,9901) info
         endif
      endif
 9901 format(' *INFORMATION: LAPACK SUBR DSTERF OR DSTEQR HAS FAILED TO
     &FIND ALL OF THE EIGENVALUES IN A TOTAL OF 30*NDOFA ITERATIONS'
     &,/,14X,' A TOTAL OF ',I8,' SUB-DIAGONAL ELEMENTS OF THE TRIDIAGONA
     &L MATRIX E HAVE NOT CONVERGED TO ZERO.'
     &,/,14X,' LAPACK WILL ATTEMPT TO USE SUBR DSTEBZ TO FIND THE EIGENV
     &ALUES.',/)
! /////////////////////////////////////////////////////////////////////E
!
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      sub_nam = 'DSTEBZ'
      CALL DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL,
     $             WORK( INDD ), WORK( INDE ), mlam, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), INFO,
     &             lowest_mode_num, highest_mode_num )

! /////////////////////////////////////////////////////////////////////B
      if (info > 0) then                 ! Call routine to print message
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
 9903 format(15x,'THE EIGENVALUES OF QUESTIONABLE VALUE ARE',/,
     &       15X,'   INDEX     EIGENVALUE')
 9904 FORMAT(15X,I8,1ES15.6)
! /////////////////////////////////////////////////////////////////////E

      do i=1,mlam
         if (method(1:3) == 'GIV') then
            eig_num(i) =  lowest_mode_num + (i - 1)
         else
            eig_num(i) = (n + 1) - (lowest_mode_num + (i - 1))
         endif
      enddo
*
      IF( WANTZ ) THEN
         sub_nam = 'DSTEIN'
         mvec = mlam
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), mvec, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
*
*        Apply transformation matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSTEIN.
*
         DO 20 J = 1, mvec
            CALL DCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL DGEMV( 'N', N, N, ONE, Q, LDQ, WORK, 1, ZERO,
     $                  Z( 1, J ), 1 )
   20    CONTINUE
      END IF

   30 CONTINUE

*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
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
*
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
*
*     End of DSBGVX
*
      END SUBROUTINE DSBGVX_GIV_MGIV

! ##################################################################################################################################
! 002 LAPACK_GIV_MGIV_EIG

      SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE ( RANGE, INFO )
      character range
      integer info

      CALL EIGENVALUE_CONVERGENCE_FAILURE_HELPER ( RANGE, INFO )

      END SUBROUTINE EIGENVALUE_CONVERGENCE_FAILURE
! ##################################################################################################################################
! 003 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DPBSTF( UPLO, N, KD, AB, LDAB, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       AB( LDAB, * )

      CALL DPBSTF_HELPER( UPLO, N, KD, AB, LDAB, INFO )

      END SUBROUTINE DPBSTF

! ##################################################################################################################################
! 004 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,
     $                   LDX, WORK, INFO )

      CHARACTER          UPLO, VECT
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N
      REAL(DOUBLE)       AB( LDAB, * ), BB( LDBB, * ), WORK( * ),
     $                   X( LDX, * )

      CALL DSBGST_HELPER( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB,
     $                    X, LDX, WORK, INFO )

      END SUBROUTINE DSBGST

      SUBROUTINE DSBGST_HELPER( VECT, UPLO, N, KA, KB, AB, LDAB, BB,
     $                          LDBB, X, LDX, WORK, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DSBGST'
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   AB( LDAB, * ), BB( LDBB, * ), WORK( * ),
     $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DSBGST reduces a real symmetric-definite banded generalized
*  eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,
*  such that C has the same bandwidth as A.
*
*  B must have been previously factorized as S**T*S by DPBSTF, using a
*  split Cholesky factorization. A is overwritten by C = X**T*A*X, where
*  X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the
*  bandwidth of A.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form the transformation matrix X;
*          = 'V':  form X.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  KA      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
*
*  KB      (input) INTEGER
*          The number of superdiagonals of the matrix B if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.
*
*  AB      (input/output) REAL(DOUBLE) array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first ka+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
*
*          On exit, the transformed matrix X**T*A*X, stored in the same
*          format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KA+1.
*
*  BB      (input) REAL(DOUBLE) array, dimension (LDBB,N)
*          The banded factor S from the split Cholesky factorization of
*          B, as returned by DPBSTF, stored in the first KB+1 rows of
*          the array.
*
*  LDBB    (input) INTEGER
*          The leading dimension of the array BB.  LDBB >= KB+1.
*
*  X       (output) REAL(DOUBLE) array, dimension (LDX,N)
*          If VECT = 'V', the n-by-n matrix X.
*          If VECT = 'N', the array X is not referenced.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.
*          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPDATE, UPPER, WANTX
      INTEGER            I, I0, I1, I2, INCA, J, J1, J1T, J2, J2T, K,
     $                   KA1, KB1, KBT, L, M, NR, NRT, NX
      integer            loopno, phase, phase_index
      REAL(DOUBLE)   BII, RA, RA1, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..

*
*     Test the input parameters
*
      WANTX = LSAME( VECT, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      KA1 = KA + 1
      KB1 = KB + 1
      INFO = 0
      IF( .NOT.WANTX .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KA.LT.0 ) THEN
         INFO = -4
      ELSE IF( KB.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -7
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -9
      ELSE IF( LDX.LT.1 .OR. WANTX .AND. LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBGST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      INCA = LDAB*KA1
*
*     Initialize X to the unit matrix, if needed
*
      IF( WANTX )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, X, LDX )
*
*     Set M to the splitting point m. It must be the same value as is
*     used in DPBSTF. The chosen value allows the arrays WORK and RWORK
*     to be of dimension (N).
*
      M = ( N+KB ) / 2
*
*     The routine works in two phases, corresponding to the two halves
*     of the split Cholesky factorization of B as S**T*S where
*
*     S = ( U    )
*         ( M  L )
*
*     with U upper triangular of order m, and L lower triangular of
*     order n-m. S has the same bandwidth as B.
*
*     S is treated as a product of elementary matrices:
*
*     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)
*
*     where S(i) is determined by the i-th row of S.
*
*     In phase 1, the index i takes the values n, n-1, ... , m+1;
*     in phase 2, it takes the values 1, 2, ... , m.
*
*     For each value of i, the current matrix A is updated by forming
*     inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside
*     the band of A. The bulge is then pushed down toward the bottom of
*     A in phase 1, and up toward the top of A in phase 2, by applying
*     plane rotations.
*
*     There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
*     of them are linearly independent, so annihilating a bulge requires
*     only 2*kb-1 plane rotations. The rotations are divided into a 1st
*     set of kb-1 rotations, and a 2nd set of kb rotations.
*
*     Wherever possible, rotations are generated and applied in vector
*     operations of length NR between the indices J1 and J2 (sometimes
*     replaced by modified values NRT, J1T or J2T).
*
*     The cosines and sines of the rotations are stored in the array
*     WORK. The cosines of the 1st set of rotations are stored in
*     elements n+2:n+m-kb-1 and the sines of the 1st set in elements
*     2:m-kb-1; the cosines of the 2nd set are stored in elements
*     n+m-kb+1:2*n and the sines of the second set in elements m-kb+1:n.
*
*     The bulges are not formed explicitly; nonzero elements outside the
*     band are created only when they are required for generating new
*     rotations; they are stored in the array WORK, in positions where
*     they are later overwritten by the sines of the rotations which
*     annihilate them.
*
*     **************************** Phase 1 *****************************
*
*     The logical structure of this phase is:
*
*     UPDATE = .TRUE.
*     DO I = N, M + 1, -1
*        use S(i) to update A and create a new bulge
*        apply rotations to push all bulges KA positions downward
*     END DO
*     UPDATE = .FALSE.
*     DO I = M + KA + 1, N - 1
*        apply rotations to push all bulges KA positions downward
*     END DO
*
*     To avoid duplicating code, the two loops are merged.
*
      UPDATE = .TRUE.
      I = N + 1
      write(sc1,*)
      phase = 1
      phase_index = 0
   10 CONTINUE
      IF( UPDATE ) THEN
         I = I - 1
         phase_index = phase_index + 1
         write(sc1,12345,advance='no') phase,phase_index,i-m,cr13_lge
         KBT = MIN( KB, I-1 )
         I0 = I - 1
         I1 = MIN( N, I+KA )
         I2 = I - KBT + KA1
         IF( I.LT.M+1 ) THEN
            UPDATE = .FALSE.
            I = I + 1
            I0 = M
            IF( KA.EQ.0 )
     $         GO TO 480
            GO TO 10
         END IF
      ELSE
         I = I + KA
         IF( I.GT.N-1 )
     $      GO TO 480
      END IF
*
      IF( UPPER ) THEN
*
*        Transform A, working with the upper triangle
*
         IF( UPDATE ) THEN
*
*           Form  inv(S(i))**T * A * inv(S(i))
*
            BII = BB( KB1, I )
            DO 20 J = I, I1
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII
   20       CONTINUE
            DO 30 J = MAX( 1, I-KA ), I
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII
   30       CONTINUE
            DO 60 K = I - KBT, I - 1
               DO 40 J = I - KBT, K
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) -
     $                               BB( J-I+KB1, I )*AB( K-I+KA1, I ) -
     $                               BB( K-I+KB1, I )*AB( J-I+KA1, I ) +
     $                               AB( KA1, I )*BB( J-I+KB1, I )*
     $                               BB( K-I+KB1, I )
   40          CONTINUE
               DO 50 J = MAX( 1, I-KA ), I - KBT - 1
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) -
     $                               BB( K-I+KB1, I )*AB( J-I+KA1, I )
   50          CONTINUE
   60       CONTINUE
            DO 80 J = I, I1
               DO 70 K = MAX( J-KA, I-KBT ), I - 1
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) -
     $                               BB( K-I+KB1, I )*AB( I-J+KA1, J )
   70          CONTINUE
   80       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by inv(S(i))
*
               CALL DSCAL( N-M, ONE / BII, X( M+1, I ), 1 )
               IF( KBT.GT.0 )
     $            CALL DGER( N-M, KBT, -ONE, X( M+1, I ), 1,
     $                       BB( KB1-KBT, I ), 1, X( M+1, I-KBT ), LDX )
            END IF
*
*           store a(i,i1) in RA1 for use in next loop over K
*
            RA1 = AB( I-I1+KA1, I1 )
         END IF
*
*        Generate and apply vectors of rotations to chase all the
*        existing bulges KA positions down toward the bottom of the
*        band
*
!        write(sc1,*)
         DO 130 K = 1, KB - 1
            loopno = 130
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            IF( UPDATE ) THEN
*
*              Determine the rotations which would annihilate the bulge
*              which has in theory just been created
*
               IF( I-K+KA.LT.N .AND. I-K.GT.1 ) THEN
*
*                 generate rotation to annihilate a(i,i-k+ka+1)
*
                  CALL DLARTG( AB( K+1, I-K+KA ), RA1,
     $                         WORK( N+I-K+KA-M ), WORK( I-K+KA-M ),
     $                         RA )
*
*                 create nonzero element a(i-k,i-k+ka+1) outside the
*                 band and store it in WORK(i-k)
*
                  T = -BB( KB1-K, I )*RA1
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T -
     $                          WORK( I-K+KA-M )*AB( 1, I-K+KA )
                  AB( 1, I-K+KA ) = WORK( I-K+KA-M )*T +
     $                              WORK( N+I-K+KA-M )*AB( 1, I-K+KA )
                  RA1 = RA
               END IF
            END IF
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            IF( UPDATE ) THEN
               J2T = MAX( J2, I+2*KA-K+1 )
            ELSE
               J2T = J2
            END IF
            NRT = ( N-J2T+KA ) / KA1
            DO 90 J = J2T, J1, KA1
*
*              create nonzero element a(j-ka,j+1) outside the band
*              and store it in WORK(j-m)
*
               WORK( J-M ) = WORK( J-M )*AB( 1, J+1 )
               AB( 1, J+1 ) = WORK( N+J-M )*AB( 1, J+1 )
   90       CONTINUE
*
*           generate rotations in 1st set to annihilate elements which
*           have been created outside the band
*
            IF( NRT.GT.0 )
     $         CALL DLARGV( NRT, AB( 1, J2T ), INCA, WORK( J2T-M ), KA1,
     $                      WORK( N+J2T-M ), KA1 )
            IF( NR.GT.0 ) THEN
*
*              apply rotations in 1st set from the right
*
               DO 100 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J2 ), INCA,
     $                         AB( KA-L, J2+1 ), INCA, WORK( N+J2-M ),
     $                         WORK( J2-M ), KA1 )
  100          CONTINUE
*
*              apply rotations in 1st set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( KA1, J2 ), AB( KA1, J2+1 ),
     $                      AB( KA, J2+1 ), INCA, WORK( N+J2-M ),
     $                      WORK( J2-M ), KA1 )
*
            END IF
*
*           start applying rotations in 1st set from the left
*
            DO 110 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA,
     $                         AB( L+1, J2+KA1-L ), INCA,
     $                         WORK( N+J2-M ), WORK( J2-M ), KA1 )
  110       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 1st set
*
               DO 120 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1,
     $                       WORK( N+J-M ), WORK( J-M ) )
  120          CONTINUE
            END IF
  130    CONTINUE
*
         IF( UPDATE ) THEN
            IF( I2.LE.N .AND. KBT.GT.0 ) THEN
*
*              create nonzero element a(i-kbt,i-kbt+ka+1) outside the
*              band and store it in WORK(i-kbt)
*
               WORK( I-KBT ) = -BB( KB1-KBT, I )*RA1
            END IF
         END IF
*
!        write(sc1,*)
         DO 170 K = KB, 1, -1
            loopno = 170
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            IF( UPDATE ) THEN
               J2 = I - K - 1 + MAX( 2, K-I0+1 )*KA1
            ELSE
               J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            END IF
*
*           finish applying rotations in 2nd set from the left
*
            DO 140 L = KB - K, 1, -1
               NRT = ( N-J2+KA+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J2-L+1 ), INCA,
     $                         AB( L+1, J2-L+1 ), INCA, WORK( N+J2-KA ),
     $                         WORK( J2-KA ), KA1 )
  140       CONTINUE
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            DO 150 J = J1, J2, -KA1
               WORK( J ) = WORK( J-KA )
               WORK( N+J ) = WORK( N+J-KA )
  150       CONTINUE
            DO 160 J = J2, J1, KA1
*
*              create nonzero element a(j-ka,j+1) outside the band
*              and store it in WORK(j)
*
               WORK( J ) = WORK( J )*AB( 1, J+1 )
               AB( 1, J+1 ) = WORK( N+J )*AB( 1, J+1 )
  160       CONTINUE
            IF( UPDATE ) THEN
               IF( I-K.LT.N-KA .AND. K.LE.KBT )
     $            WORK( I-K+KA ) = WORK( I-K )
            END IF
  170    CONTINUE
*
!        write(sc1,*)
         DO 210 K = KB, 1, -1
            loopno = 210
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            IF( NR.GT.0 ) THEN
*
*              generate rotations in 2nd set to annihilate elements
*              which have been created outside the band
*
               CALL DLARGV( NR, AB( 1, J2 ), INCA, WORK( J2 ), KA1,
     $                      WORK( N+J2 ), KA1 )
*
*              apply rotations in 2nd set from the right
*
               DO 180 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J2 ), INCA,
     $                         AB( KA-L, J2+1 ), INCA, WORK( N+J2 ),
     $                         WORK( J2 ), KA1 )
  180          CONTINUE
*
*              apply rotations in 2nd set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( KA1, J2 ), AB( KA1, J2+1 ),
     $                      AB( KA, J2+1 ), INCA, WORK( N+J2 ),
     $                      WORK( J2 ), KA1 )
*
            END IF
*
*           start applying rotations in 2nd set from the left
*
            DO 190 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA,
     $                         AB( L+1, J2+KA1-L ), INCA, WORK( N+J2 ),
     $                         WORK( J2 ), KA1 )
  190       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 2nd set
*
               DO 200 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1,
     $                       WORK( N+J ), WORK( J ) )
  200          CONTINUE
            END IF
  210    CONTINUE
*
!        write(sc1,*)
         DO 230 K = 1, KB - 1
            loopno = 230
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
*
*           finish applying rotations in 1st set from the left
*
            DO 220 L = KB - K, 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J2+KA1-L ), INCA,
     $                         AB( L+1, J2+KA1-L ), INCA,
     $                         WORK( N+J2-M ), WORK( J2-M ), KA1 )
  220       CONTINUE
  230    CONTINUE
*
         IF( KB.GT.1 ) THEN
            DO 240 J = N - 1, I - KB + 2*KA + 1, -1
               WORK( N+J-M ) = WORK( N+J-KA-M )
               WORK( J-M ) = WORK( J-KA-M )
  240       CONTINUE
         END IF
*
      ELSE
*
*        Transform A, working with the lower triangle
*
         IF( UPDATE ) THEN
*
*           Form  inv(S(i))**T * A * inv(S(i))
*
            BII = BB( 1, I )
            DO 250 J = I, I1
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII
  250       CONTINUE
            DO 260 J = MAX( 1, I-KA ), I
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII
  260       CONTINUE
            DO 290 K = I - KBT, I - 1
               DO 270 J = I - KBT, K
                  AB( K-J+1, J ) = AB( K-J+1, J ) -
     $                             BB( I-J+1, J )*AB( I-K+1, K ) -
     $                             BB( I-K+1, K )*AB( I-J+1, J ) +
     $                             AB( 1, I )*BB( I-J+1, J )*
     $                             BB( I-K+1, K )
  270          CONTINUE
               DO 280 J = MAX( 1, I-KA ), I - KBT - 1
                  AB( K-J+1, J ) = AB( K-J+1, J ) -
     $                             BB( I-K+1, K )*AB( I-J+1, J )
  280          CONTINUE
  290       CONTINUE
            DO 310 J = I, I1
               DO 300 K = MAX( J-KA, I-KBT ), I - 1
                  AB( J-K+1, K ) = AB( J-K+1, K ) -
     $                             BB( I-K+1, K )*AB( J-I+1, I )
  300          CONTINUE
  310       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by inv(S(i))
*
               CALL DSCAL( N-M, ONE / BII, X( M+1, I ), 1 )
               IF( KBT.GT.0 )
     $            CALL DGER( N-M, KBT, -ONE, X( M+1, I ), 1,
     $                       BB( KBT+1, I-KBT ), LDBB-1,
     $                       X( M+1, I-KBT ), LDX )
            END IF
*
*           store a(i1,i) in RA1 for use in next loop over K
*
            RA1 = AB( I1-I+1, I )
         END IF
*
*        Generate and apply vectors of rotations to chase all the
*        existing bulges KA positions down toward the bottom of the
*        band
*
         DO 360 K = 1, KB - 1
            IF( UPDATE ) THEN
*
*              Determine the rotations which would annihilate the bulge
*              which has in theory just been created
*
               IF( I-K+KA.LT.N .AND. I-K.GT.1 ) THEN
*
*                 generate rotation to annihilate a(i-k+ka+1,i)
*
                  CALL DLARTG( AB( KA1-K, I ), RA1, WORK( N+I-K+KA-M ),
     $                         WORK( I-K+KA-M ), RA )
*
*                 create nonzero element a(i-k+ka+1,i-k) outside the
*                 band and store it in WORK(i-k)
*
                  T = -BB( K+1, I-K )*RA1
                  WORK( I-K ) = WORK( N+I-K+KA-M )*T -
     $                          WORK( I-K+KA-M )*AB( KA1, I-K )
                  AB( KA1, I-K ) = WORK( I-K+KA-M )*T +
     $                             WORK( N+I-K+KA-M )*AB( KA1, I-K )
                  RA1 = RA
               END IF
            END IF
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            IF( UPDATE ) THEN
               J2T = MAX( J2, I+2*KA-K+1 )
            ELSE
               J2T = J2
            END IF
            NRT = ( N-J2T+KA ) / KA1
            DO 320 J = J2T, J1, KA1
*
*              create nonzero element a(j+1,j-ka) outside the band
*              and store it in WORK(j-m)
*
               WORK( J-M ) = WORK( J-M )*AB( KA1, J-KA+1 )
               AB( KA1, J-KA+1 ) = WORK( N+J-M )*AB( KA1, J-KA+1 )
  320       CONTINUE
*
*           generate rotations in 1st set to annihilate elements which
*           have been created outside the band
*
            IF( NRT.GT.0 )
     $         CALL DLARGV( NRT, AB( KA1, J2T-KA ), INCA, WORK( J2T-M ),
     $                      KA1, WORK( N+J2T-M ), KA1 )
            IF( NR.GT.0 ) THEN
*
*              apply rotations in 1st set from the left
*
               DO 330 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J2-L ), INCA,
     $                         AB( L+2, J2-L ), INCA, WORK( N+J2-M ),
     $                         WORK( J2-M ), KA1 )
  330          CONTINUE
*
*              apply rotations in 1st set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ),
     $                      INCA, WORK( N+J2-M ), WORK( J2-M ), KA1 )
*
            END IF
*
*           start applying rotations in 1st set from the right
*
            DO 340 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA,
     $                         AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ),
     $                         WORK( J2-M ), KA1 )
  340       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 1st set
*
               DO 350 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1,
     $                       WORK( N+J-M ), WORK( J-M ) )
  350          CONTINUE
            END IF
  360    CONTINUE
*
         IF( UPDATE ) THEN
            IF( I2.LE.N .AND. KBT.GT.0 ) THEN
*
*              create nonzero element a(i-kbt+ka+1,i-kbt) outside the
*              band and store it in WORK(i-kbt)
*
               WORK( I-KBT ) = -BB( KBT+1, I-KBT )*RA1
            END IF
         END IF
*
         DO 400 K = KB, 1, -1
            IF( UPDATE ) THEN
               J2 = I - K - 1 + MAX( 2, K-I0+1 )*KA1
            ELSE
               J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            END IF
*
*           finish applying rotations in 2nd set from the right
*
            DO 370 L = KB - K, 1, -1
               NRT = ( N-J2+KA+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J2-KA ), INCA,
     $                         AB( KA1-L, J2-KA+1 ), INCA,
     $                         WORK( N+J2-KA ), WORK( J2-KA ), KA1 )
  370       CONTINUE
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            DO 380 J = J1, J2, -KA1
               WORK( J ) = WORK( J-KA )
               WORK( N+J ) = WORK( N+J-KA )
  380       CONTINUE
            DO 390 J = J2, J1, KA1
*
*              create nonzero element a(j+1,j-ka) outside the band
*              and store it in WORK(j)
*
               WORK( J ) = WORK( J )*AB( KA1, J-KA+1 )
               AB( KA1, J-KA+1 ) = WORK( N+J )*AB( KA1, J-KA+1 )
  390       CONTINUE
            IF( UPDATE ) THEN
               IF( I-K.LT.N-KA .AND. K.LE.KBT )
     $            WORK( I-K+KA ) = WORK( I-K )
            END IF
  400    CONTINUE
*
         DO 440 K = KB, 1, -1
            J2 = I - K - 1 + MAX( 1, K-I0+1 )*KA1
            NR = ( N-J2+KA ) / KA1
            J1 = J2 + ( NR-1 )*KA1
            IF( NR.GT.0 ) THEN
*
*              generate rotations in 2nd set to annihilate elements
*              which have been created outside the band
*
               CALL DLARGV( NR, AB( KA1, J2-KA ), INCA, WORK( J2 ), KA1,
     $                      WORK( N+J2 ), KA1 )
*
*              apply rotations in 2nd set from the left
*
               DO 410 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J2-L ), INCA,
     $                         AB( L+2, J2-L ), INCA, WORK( N+J2 ),
     $                         WORK( J2 ), KA1 )
  410          CONTINUE
*
*              apply rotations in 2nd set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( 1, J2 ), AB( 1, J2+1 ), AB( 2, J2 ),
     $                      INCA, WORK( N+J2 ), WORK( J2 ), KA1 )
*
            END IF
*
*           start applying rotations in 2nd set from the right
*
            DO 420 L = KA - 1, KB - K + 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA,
     $                         AB( KA1-L, J2+1 ), INCA, WORK( N+J2 ),
     $                         WORK( J2 ), KA1 )
  420       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 2nd set
*
               DO 430 J = J2, J1, KA1
                  CALL DROT( N-M, X( M+1, J ), 1, X( M+1, J+1 ), 1,
     $                       WORK( N+J ), WORK( J ) )
  430          CONTINUE
            END IF
  440    CONTINUE
*
         DO 460 K = 1, KB - 1
            J2 = I - K - 1 + MAX( 1, K-I0+2 )*KA1
*
*           finish applying rotations in 1st set from the right
*
            DO 450 L = KB - K, 1, -1
               NRT = ( N-J2+L ) / KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J2 ), INCA,
     $                         AB( KA1-L, J2+1 ), INCA, WORK( N+J2-M ),
     $                         WORK( J2-M ), KA1 )
  450       CONTINUE
  460    CONTINUE
*
         IF( KB.GT.1 ) THEN
            DO 470 J = N - 1, I - KB + 2*KA + 1, -1
               WORK( N+J-M ) = WORK( N+J-KA-M )
               WORK( J-M ) = WORK( J-KA-M )
  470       CONTINUE
         END IF
*
      END IF
*
      GO TO 10
*
  480 CONTINUE
*
*     **************************** Phase 2 *****************************
*
*     The logical structure of this phase is:
*
*     UPDATE = .TRUE.
*     DO I = 1, M
*        use S(i) to update A and create a new bulge
*        apply rotations to push all bulges KA positions upward
*     END DO
*     UPDATE = .FALSE.
*     DO I = M - KA - 1, 2, -1
*        apply rotations to push all bulges KA positions upward
*     END DO
*
*     To avoid duplicating code, the two loops are merged.
*
      UPDATE = .TRUE.
      I = 0
      write(sc1,*)
      phase = 1
  490 CONTINUE
      IF( UPDATE ) THEN
         write(sc1,12345,advance='no') phase,i,m,cr13_lge
         I = I + 1
         KBT = MIN( KB, M-I )
         I0 = I + 1
         I1 = MAX( 1, I-KA )
         I2 = I + KBT - KA1
         IF( I.GT.M ) THEN
            UPDATE = .FALSE.
            I = I - 1
            I0 = M + 1
            IF( KA.EQ.0 )
     $         RETURN
            GO TO 490
         END IF
      ELSE
         I = I - KA
         IF( I.LT.2 )
     $      RETURN
      END IF
*
      IF( I.LT.M-KBT ) THEN
         NX = M
      ELSE
         NX = N
      END IF
*
      IF( UPPER ) THEN
*
*        Transform A, working with the upper triangle
*
         IF( UPDATE ) THEN
*
*           Form  inv(S(i))**T * A * inv(S(i))
*
            BII = BB( KB1, I )
            DO 500 J = I1, I
               AB( J-I+KA1, I ) = AB( J-I+KA1, I ) / BII
  500       CONTINUE
            DO 510 J = I, MIN( N, I+KA )
               AB( I-J+KA1, J ) = AB( I-J+KA1, J ) / BII
  510       CONTINUE
            DO 540 K = I + 1, I + KBT
               DO 520 J = K, I + KBT
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) -
     $                               BB( I-J+KB1, J )*AB( I-K+KA1, K ) -
     $                               BB( I-K+KB1, K )*AB( I-J+KA1, J ) +
     $                               AB( KA1, I )*BB( I-J+KB1, J )*
     $                               BB( I-K+KB1, K )
  520          CONTINUE
               DO 530 J = I + KBT + 1, MIN( N, I+KA )
                  AB( K-J+KA1, J ) = AB( K-J+KA1, J ) -
     $                               BB( I-K+KB1, K )*AB( I-J+KA1, J )
  530          CONTINUE
  540       CONTINUE
            DO 560 J = I1, I
               DO 550 K = I + 1, MIN( J+KA, I+KBT )
                  AB( J-K+KA1, K ) = AB( J-K+KA1, K ) -
     $                               BB( I-K+KB1, K )*AB( J-I+KA1, I )
  550          CONTINUE
  560       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by inv(S(i))
*
               CALL DSCAL( NX, ONE / BII, X( 1, I ), 1 )
               IF( KBT.GT.0 )
     $            CALL DGER( NX, KBT, -ONE, X( 1, I ), 1, BB( KB, I+1 ),
     $                       LDBB-1, X( 1, I+1 ), LDX )
            END IF
*
*           store a(i1,i) in RA1 for use in next loop over K
*
            RA1 = AB( I1-I+KA1, I )
         END IF
*
*        Generate and apply vectors of rotations to chase all the
*        existing bulges KA positions up toward the top of the band
*
!        write(sc1,*)
         DO 610 K = 1, KB - 1
            loopno = 610
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            IF( UPDATE ) THEN
*
*              Determine the rotations which would annihilate the bulge
*              which has in theory just been created
*
               IF( I+K-KA1.GT.0 .AND. I+K.LT.M ) THEN
*
*                 generate rotation to annihilate a(i+k-ka-1,i)
*
                  CALL DLARTG( AB( K+1, I ), RA1, WORK( N+I+K-KA ),
     $                         WORK( I+K-KA ), RA )
*
*                 create nonzero element a(i+k-ka-1,i+k) outside the
*                 band and store it in WORK(m-kb+i+k)
*
                  T = -BB( KB1-K, I+K )*RA1
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T -
     $                               WORK( I+K-KA )*AB( 1, I+K )
                  AB( 1, I+K ) = WORK( I+K-KA )*T +
     $                           WORK( N+I+K-KA )*AB( 1, I+K )
                  RA1 = RA
               END IF
            END IF
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            IF( UPDATE ) THEN
               J2T = MIN( J2, I-2*KA+K-1 )
            ELSE
               J2T = J2
            END IF
            NRT = ( J2T+KA-1 ) / KA1
            DO 570 J = J1, J2T, KA1
*
*              create nonzero element a(j-1,j+ka) outside the band
*              and store it in WORK(j)
*
               WORK( J ) = WORK( J )*AB( 1, J+KA-1 )
               AB( 1, J+KA-1 ) = WORK( N+J )*AB( 1, J+KA-1 )
  570       CONTINUE
*
*           generate rotations in 1st set to annihilate elements which
*           have been created outside the band
*
            IF( NRT.GT.0 )
     $         CALL DLARGV( NRT, AB( 1, J1+KA ), INCA, WORK( J1 ), KA1,
     $                      WORK( N+J1 ), KA1 )
            IF( NR.GT.0 ) THEN
*
*              apply rotations in 1st set from the left
*
               DO 580 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J1+L ), INCA,
     $                         AB( KA-L, J1+L ), INCA, WORK( N+J1 ),
     $                         WORK( J1 ), KA1 )
  580          CONTINUE
*
*              apply rotations in 1st set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( KA1, J1 ), AB( KA1, J1-1 ),
     $                      AB( KA, J1 ), INCA, WORK( N+J1 ),
     $                      WORK( J1 ), KA1 )
*
            END IF
*
*           start applying rotations in 1st set from the right
*
            DO 590 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J1T ), INCA,
     $                         AB( L+1, J1T-1 ), INCA, WORK( N+J1T ),
     $                         WORK( J1T ), KA1 )
  590       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 1st set
*
               DO 600 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1,
     $                       WORK( N+J ), WORK( J ) )
  600          CONTINUE
            END IF
  610    CONTINUE
*
         IF( UPDATE ) THEN
            IF( I2.GT.0 .AND. KBT.GT.0 ) THEN
*
*              create nonzero element a(i+kbt-ka-1,i+kbt) outside the
*              band and store it in WORK(m-kb+i+kbt)
*
               WORK( M-KB+I+KBT ) = -BB( KB1-KBT, I+KBT )*RA1
            END IF
         END IF
*
!        write(sc1,*)
         DO 650 K = KB, 1, -1
            loopno = 650
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            IF( UPDATE ) THEN
               J2 = I + K + 1 - MAX( 2, K+I0-M )*KA1
            ELSE
               J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            END IF
*
*           finish applying rotations in 2nd set from the right
*
            DO 620 L = KB - K, 1, -1
               NRT = ( J2+KA+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J1T+KA ), INCA,
     $                         AB( L+1, J1T+KA-1 ), INCA,
     $                         WORK( N+M-KB+J1T+KA ),
     $                         WORK( M-KB+J1T+KA ), KA1 )
  620       CONTINUE
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            DO 630 J = J1, J2, KA1
               WORK( M-KB+J ) = WORK( M-KB+J+KA )
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA )
  630       CONTINUE
            DO 640 J = J1, J2, KA1
*
*              create nonzero element a(j-1,j+ka) outside the band
*              and store it in WORK(m-kb+j)
*
               WORK( M-KB+J ) = WORK( M-KB+J )*AB( 1, J+KA-1 )
               AB( 1, J+KA-1 ) = WORK( N+M-KB+J )*AB( 1, J+KA-1 )
  640       CONTINUE
            IF( UPDATE ) THEN
               IF( I+K.GT.KA1 .AND. K.LE.KBT )
     $            WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K )
            END IF
  650    CONTINUE
*
!        write(sc1,*)
         DO 690 K = KB, 1, -1
            loopno = 690
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            IF( NR.GT.0 ) THEN
*
*              generate rotations in 2nd set to annihilate elements
*              which have been created outside the band
*
               CALL DLARGV( NR, AB( 1, J1+KA ), INCA, WORK( M-KB+J1 ),
     $                      KA1, WORK( N+M-KB+J1 ), KA1 )
*
*              apply rotations in 2nd set from the left
*
               DO 660 L = 1, KA - 1
                  CALL DLARTV( NR, AB( KA1-L, J1+L ), INCA,
     $                         AB( KA-L, J1+L ), INCA,
     $                         WORK( N+M-KB+J1 ), WORK( M-KB+J1 ), KA1 )
  660          CONTINUE
*
*              apply rotations in 2nd set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( KA1, J1 ), AB( KA1, J1-1 ),
     $                      AB( KA, J1 ), INCA, WORK( N+M-KB+J1 ),
     $                      WORK( M-KB+J1 ), KA1 )
*
            END IF
*
*           start applying rotations in 2nd set from the right
*
            DO 670 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J1T ), INCA,
     $                         AB( L+1, J1T-1 ), INCA,
     $                         WORK( N+M-KB+J1T ), WORK( M-KB+J1T ),
     $                         KA1 )
  670       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 2nd set
*
               DO 680 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1,
     $                       WORK( N+M-KB+J ), WORK( M-KB+J ) )
  680          CONTINUE
            END IF
  690    CONTINUE
*
!        write(sc1,*)
         DO 710 K = 1, KB - 1
            loopno = 710
!           write(sc1,22345,advance='no') loopno,k,kb,cr13_lge
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
*
*           finish applying rotations in 1st set from the right
*
            DO 700 L = KB - K, 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( L, J1T ), INCA,
     $                         AB( L+1, J1T-1 ), INCA, WORK( N+J1T ),
     $                         WORK( J1T ), KA1 )
  700       CONTINUE
  710    CONTINUE
*
         IF( KB.GT.1 ) THEN
            DO 720 J = 2, MIN( I+KB, M ) - 2*KA - 1
               WORK( N+J ) = WORK( N+J+KA )
               WORK( J ) = WORK( J+KA )
  720       CONTINUE
         END IF
*
      ELSE
*
*        Transform A, working with the lower triangle
*
         IF( UPDATE ) THEN
*
*           Form  inv(S(i))**T * A * inv(S(i))
*
            BII = BB( 1, I )
            DO 730 J = I1, I
               AB( I-J+1, J ) = AB( I-J+1, J ) / BII
  730       CONTINUE
            DO 740 J = I, MIN( N, I+KA )
               AB( J-I+1, I ) = AB( J-I+1, I ) / BII
  740       CONTINUE
            DO 770 K = I + 1, I + KBT
               DO 750 J = K, I + KBT
                  AB( J-K+1, K ) = AB( J-K+1, K ) -
     $                             BB( J-I+1, I )*AB( K-I+1, I ) -
     $                             BB( K-I+1, I )*AB( J-I+1, I ) +
     $                             AB( 1, I )*BB( J-I+1, I )*
     $                             BB( K-I+1, I )
  750          CONTINUE
               DO 760 J = I + KBT + 1, MIN( N, I+KA )
                  AB( J-K+1, K ) = AB( J-K+1, K ) -
     $                             BB( K-I+1, I )*AB( J-I+1, I )
  760          CONTINUE
  770       CONTINUE
            DO 790 J = I1, I
               DO 780 K = I + 1, MIN( J+KA, I+KBT )
                  AB( K-J+1, J ) = AB( K-J+1, J ) -
     $                             BB( K-I+1, I )*AB( I-J+1, J )
  780          CONTINUE
  790       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by inv(S(i))
*
               CALL DSCAL( NX, ONE / BII, X( 1, I ), 1 )
               IF( KBT.GT.0 )
     $            CALL DGER( NX, KBT, -ONE, X( 1, I ), 1, BB( 2, I ), 1,
     $                       X( 1, I+1 ), LDX )
            END IF
*
*           store a(i,i1) in RA1 for use in next loop over K
*
            RA1 = AB( I-I1+1, I1 )
         END IF
*
*        Generate and apply vectors of rotations to chase all the
*        existing bulges KA positions up toward the top of the band
*
         DO 840 K = 1, KB - 1
            IF( UPDATE ) THEN
*
*              Determine the rotations which would annihilate the bulge
*              which has in theory just been created
*
               IF( I+K-KA1.GT.0 .AND. I+K.LT.M ) THEN
*
*                 generate rotation to annihilate a(i,i+k-ka-1)
*
                  CALL DLARTG( AB( KA1-K, I+K-KA ), RA1,
     $                         WORK( N+I+K-KA ), WORK( I+K-KA ), RA )
*
*                 create nonzero element a(i+k,i+k-ka-1) outside the
*                 band and store it in WORK(m-kb+i+k)
*
                  T = -BB( K+1, I )*RA1
                  WORK( M-KB+I+K ) = WORK( N+I+K-KA )*T -
     $                               WORK( I+K-KA )*AB( KA1, I+K-KA )
                  AB( KA1, I+K-KA ) = WORK( I+K-KA )*T +
     $                                WORK( N+I+K-KA )*AB( KA1, I+K-KA )
                  RA1 = RA
               END IF
            END IF
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            IF( UPDATE ) THEN
               J2T = MIN( J2, I-2*KA+K-1 )
            ELSE
               J2T = J2
            END IF
            NRT = ( J2T+KA-1 ) / KA1
            DO 800 J = J1, J2T, KA1
*
*              create nonzero element a(j+ka,j-1) outside the band
*              and store it in WORK(j)
*
               WORK( J ) = WORK( J )*AB( KA1, J-1 )
               AB( KA1, J-1 ) = WORK( N+J )*AB( KA1, J-1 )
  800       CONTINUE
*
*           generate rotations in 1st set to annihilate elements which
*           have been created outside the band
*
            IF( NRT.GT.0 )
     $         CALL DLARGV( NRT, AB( KA1, J1 ), INCA, WORK( J1 ), KA1,
     $                      WORK( N+J1 ), KA1 )
            IF( NR.GT.0 ) THEN
*
*              apply rotations in 1st set from the right
*
               DO 810 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ),
     $                         INCA, WORK( N+J1 ), WORK( J1 ), KA1 )
  810          CONTINUE
*
*              apply rotations in 1st set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( 1, J1 ), AB( 1, J1-1 ),
     $                      AB( 2, J1-1 ), INCA, WORK( N+J1 ),
     $                      WORK( J1 ), KA1 )
*
            END IF
*
*           start applying rotations in 1st set from the left
*
            DO 820 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA,
     $                         AB( KA1-L, J1T-KA1+L ), INCA,
     $                         WORK( N+J1T ), WORK( J1T ), KA1 )
  820       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 1st set
*
               DO 830 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1,
     $                       WORK( N+J ), WORK( J ) )
  830          CONTINUE
            END IF
  840    CONTINUE
*
         IF( UPDATE ) THEN
            IF( I2.GT.0 .AND. KBT.GT.0 ) THEN
*
*              create nonzero element a(i+kbt,i+kbt-ka-1) outside the
*              band and store it in WORK(m-kb+i+kbt)
*
               WORK( M-KB+I+KBT ) = -BB( KBT+1, I )*RA1
            END IF
         END IF
*
         DO 880 K = KB, 1, -1
            IF( UPDATE ) THEN
               J2 = I + K + 1 - MAX( 2, K+I0-M )*KA1
            ELSE
               J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            END IF
*
*           finish applying rotations in 2nd set from the left
*
            DO 850 L = KB - K, 1, -1
               NRT = ( J2+KA+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J1T+L-1 ), INCA,
     $                         AB( KA1-L, J1T+L-1 ), INCA,
     $                         WORK( N+M-KB+J1T+KA ),
     $                         WORK( M-KB+J1T+KA ), KA1 )
  850       CONTINUE
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            DO 860 J = J1, J2, KA1
               WORK( M-KB+J ) = WORK( M-KB+J+KA )
               WORK( N+M-KB+J ) = WORK( N+M-KB+J+KA )
  860       CONTINUE
            DO 870 J = J1, J2, KA1
*
*              create nonzero element a(j+ka,j-1) outside the band
*              and store it in WORK(m-kb+j)
*
               WORK( M-KB+J ) = WORK( M-KB+J )*AB( KA1, J-1 )
               AB( KA1, J-1 ) = WORK( N+M-KB+J )*AB( KA1, J-1 )
  870       CONTINUE
            IF( UPDATE ) THEN
               IF( I+K.GT.KA1 .AND. K.LE.KBT )
     $            WORK( M-KB+I+K-KA ) = WORK( M-KB+I+K )
            END IF
  880    CONTINUE
*
         DO 920 K = KB, 1, -1
            J2 = I + K + 1 - MAX( 1, K+I0-M )*KA1
            NR = ( J2+KA-1 ) / KA1
            J1 = J2 - ( NR-1 )*KA1
            IF( NR.GT.0 ) THEN
*
*              generate rotations in 2nd set to annihilate elements
*              which have been created outside the band
*
               CALL DLARGV( NR, AB( KA1, J1 ), INCA, WORK( M-KB+J1 ),
     $                      KA1, WORK( N+M-KB+J1 ), KA1 )
*
*              apply rotations in 2nd set from the right
*
               DO 890 L = 1, KA - 1
                  CALL DLARTV( NR, AB( L+1, J1 ), INCA, AB( L+2, J1-1 ),
     $                         INCA, WORK( N+M-KB+J1 ), WORK( M-KB+J1 ),
     $                         KA1 )
  890          CONTINUE
*
*              apply rotations in 2nd set from both sides to diagonal
*              blocks
*
               CALL DLAR2V( NR, AB( 1, J1 ), AB( 1, J1-1 ),
     $                      AB( 2, J1-1 ), INCA, WORK( N+M-KB+J1 ),
     $                      WORK( M-KB+J1 ), KA1 )
*
            END IF
*
*           start applying rotations in 2nd set from the left
*
            DO 900 L = KA - 1, KB - K + 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA,
     $                         AB( KA1-L, J1T-KA1+L ), INCA,
     $                         WORK( N+M-KB+J1T ), WORK( M-KB+J1T ),
     $                         KA1 )
  900       CONTINUE
*
            IF( WANTX ) THEN
*
*              post-multiply X by product of rotations in 2nd set
*
               DO 910 J = J1, J2, KA1
                  CALL DROT( NX, X( 1, J ), 1, X( 1, J-1 ), 1,
     $                       WORK( N+M-KB+J ), WORK( M-KB+J ) )
  910          CONTINUE
            END IF
  920    CONTINUE
*
         DO 940 K = 1, KB - 1
            J2 = I + K + 1 - MAX( 1, K+I0-M+1 )*KA1
*
*           finish applying rotations in 1st set from the left
*
            DO 930 L = KB - K, 1, -1
               NRT = ( J2+L-1 ) / KA1
               J1T = J2 - ( NRT-1 )*KA1
               IF( NRT.GT.0 )
     $            CALL DLARTV( NRT, AB( KA1-L+1, J1T-KA1+L ), INCA,
     $                         AB( KA1-L, J1T-KA1+L ), INCA,
     $                         WORK( N+J1T ), WORK( J1T ), KA1 )
  930       CONTINUE
  940    CONTINUE
*
         IF( KB.GT.1 ) THEN
            DO 950 J = 2, MIN( I+KB, M ) - 2*KA - 1
               WORK( N+J ) = WORK( N+J+KA )
               WORK( J ) = WORK( J+KA )
  950       CONTINUE
         END IF
*
      END IF
*
      GO TO 490


      RETURN
*
*     End of DSBGST

! **********************************************************************************************************************************
22345 format(5X,'Loop ',i8,': K = ',i8,' of ',i8, a)

12345 format(5X,'Phase ',i1,': Updating from index ',i8,' to ',i8, a)

! **********************************************************************************************************************************

      END SUBROUTINE DSBGST_HELPER

! ##################################################################################################################################
! 005 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
     $                   WORK, INFO )

      CHARACTER          UPLO, VECT
      INTEGER            INFO, KD, LDAB, LDQ, N
      REAL(DOUBLE)       AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),
     $                   WORK( * )

      CALL DSBTRD_HELPER( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
     $                    WORK, INFO )

      END SUBROUTINE DSBTRD

      SUBROUTINE DSBTRD_HELPER( VECT, UPLO, N, KD, AB, LDAB, D, E, Q,
     $                          LDQ, WORK, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DSBTRD'
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KD, LDAB, LDQ, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSBTRD reduces a real symmetric band matrix A to symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form Q;
*          = 'V':  form Q;
*          = 'U':  update a matrix X, by forming X*Q.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) REAL(DOUBLE) array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the diagonal elements of AB are overwritten by the
*          diagonal elements of the tridiagonal matrix T; if KD > 0, the
*          elements on the first superdiagonal (if UPLO = 'U') or the
*          first subdiagonal (if UPLO = 'L') are overwritten by the
*          off-diagonal elements of T; the rest of AB is overwritten by
*          values generated during the reduction.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  D       (output) REAL(DOUBLE) array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
*  E       (output) REAL(DOUBLE) array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
*  Q       (input/output) REAL(DOUBLE) array, dimension (LDQ,N)
*          On entry, if VECT = 'U', then Q must contain an N-by-N
*          matrix X; if VECT = 'N' or 'V', then Q need not be set.
*
*          On exit:
*          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;
*          if VECT = 'U', Q contains the product X*Q;
*          if VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  Modified by Linda Kaufman, Bell Labs.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            INITQ, UPPER, WANTQ
      INTEGER            I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J,
     $                   J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1,
     $                   KDM1, KDN, L, LAST, LEND, NQ, NR, NRT
      REAL(DOUBLE)   TEMP
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..

*
*     Test the input parameters
*
      INITQ = LSAME( VECT, 'V' )
      WANTQ = INITQ .OR. LSAME( VECT, 'U' )
      UPPER = LSAME( UPLO, 'U' )
      KD1 = KD + 1
      KDM1 = KD - 1
      INCX = LDAB - 1
      IQEND = 1
*
      INFO = 0
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD1 ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSBTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize Q to the unit matrix, if needed
*
      IF( INITQ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
*
*     Wherever possible, plane rotations are generated and applied in
*     vector operations of length NR over the index set J1:J2:KD1.
*
*     The cosines and sines of the plane rotations are stored in the
*     arrays D and WORK.
*
      INCA = KD1*LDAB
      KDN = MIN( N-1, KD )
      IF( UPPER ) THEN
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to tridiagonal form, working with upper triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            write(sc1,*)
            DO 90 I = 1, N - 2
               IF (NOCOUNTS .NE. 'Y') THEN
                  write(sc1,12345,advance='no') i,n-2,cr13_lge
               ENDIF
*
*              Reduce i-th row of matrix to tridiagonal form
*
               DO 80 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL DLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ),
     $                            KD1, D( J1 ), KD1 )
*
*                    apply rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    DLARTV or DROT is used
*
                     IF( NR.GE.2*KD-1 ) THEN
                        DO 10 L = 1, KD - 1
                           CALL DLARTV( NR, AB( L+1, J1-1 ), INCA,
     $                                  AB( L, J1 ), INCA, D( J1 ),
     $                                  WORK( J1 ), KD1 )
   10                   CONTINUE
*
                     ELSE
                        JEND = J1 + ( NR-1 )*KD1
                        DO 20 JINC = J1, JEND, KD1
                           CALL DROT( KDM1, AB( 2, JINC-1 ), 1,
     $                                AB( 1, JINC ), 1, D( JINC ),
     $                                WORK( JINC ) )
   20                   CONTINUE
                     END IF
                  END IF
*
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i,i+k-1)
*                       within the band
*
                        CALL DLARTG( AB( KD-K+3, I+K-2 ),
     $                               AB( KD-K+2, I+K-1 ), D( I+K-1 ),
     $                               WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP
*
*                       apply rotation from the right
*
                        CALL DROT( K-3, AB( KD-K+4, I+K-2 ), 1,
     $                             AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL DLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ),
     $                            AB( KD, J1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the left
*
                  IF( NR.GT.0 ) THEN
                     IF( 2*KD-1.LT.NR ) THEN
*
*                    Dependent on the the number of diagonals either
*                    DLARTV or DROT is used
*
                        DO 30 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL DLARTV( NRT, AB( KD-L, J1+L ), INCA,
     $                                     AB( KD-L+1, J1+L ), INCA,
     $                                     D( J1 ), WORK( J1 ), KD1 )
   30                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 40 JIN = J1, J1END, KD1
                              CALL DROT( KD-1, AB( KD-1, JIN+1 ), INCX,
     $                                   AB( KD, JIN+1 ), INCX,
     $                                   D( JIN ), WORK( JIN ) )
   40                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL DROT( LEND, AB( KD-1, LAST+1 ), INCX,
     $                                AB( KD, LAST+1 ), INCX, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 50 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL DROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), WORK( J ) )
   50                   CONTINUE
                     ELSE
*
                        DO 60 J = J1, J2, KD1
                           CALL DROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), WORK( J ) )
   60                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 70 J = J1, J2, KD1
*
*                    create nonzero element a(j-1,j+kd) outside the band
*                    and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           copy off-diagonal elements to E
*
            DO 100 I = 1, N - 1
               E( I ) = AB( KD, I+1 )
  100       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 110 I = 1, N - 1
               E( I ) = ZERO
  110       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 120 I = 1, N
            D( I ) = AB( KD1, I )
  120    CONTINUE
*
      ELSE
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to tridiagonal form, working with lower triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            DO 210 I = 1, N - 2
*
*              Reduce i-th column of matrix to tridiagonal form
*
               DO 200 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL DLARGV( NR, AB( KD1, J1-KD1 ), INCA,
     $                            WORK( J1 ), KD1, D( J1 ), KD1 )
*
*                    apply plane rotations from one side
*
*
*                    Dependent on the the number of diagonals either
*                    DLARTV or DROT is used
*
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 130 L = 1, KD - 1
                           CALL DLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA,
     $                                  AB( KD1-L+1, J1-KD1+L ), INCA,
     $                                  D( J1 ), WORK( J1 ), KD1 )
  130                   CONTINUE
                     ELSE
                        JEND = J1 + KD1*( NR-1 )
                        DO 140 JINC = J1, JEND, KD1
                           CALL DROT( KDM1, AB( KD, JINC-KD ), INCX,
     $                                AB( KD1, JINC-KD ), INCX,
     $                                D( JINC ), WORK( JINC ) )
  140                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i+k-1,i)
*                       within the band
*
                        CALL DLARTG( AB( K-1, I ), AB( K, I ),
     $                               D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( K-1, I ) = TEMP
*
*                       apply rotation from the left
*
                        CALL DROT( K-3, AB( K-2, I+1 ), LDAB-1,
     $                             AB( K-1, I+1 ), LDAB-1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL DLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ),
     $                            AB( 2, J1-1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    DLARTV or DROT is used
*
                  IF( NR.GT.0 ) THEN
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 150 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL DLARTV( NRT, AB( L+2, J1-1 ), INCA,
     $                                     AB( L+1, J1 ), INCA, D( J1 ),
     $                                     WORK( J1 ), KD1 )
  150                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 160 J1INC = J1, J1END, KD1
                              CALL DROT( KDM1, AB( 3, J1INC-1 ), 1,
     $                                   AB( 2, J1INC ), 1, D( J1INC ),
     $                                   WORK( J1INC ) )
  160                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL DROT( LEND, AB( 3, LAST-1 ), 1,
     $                                AB( 2, LAST ), 1, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
*
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 170 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL DROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), WORK( J ) )
  170                   CONTINUE
                     ELSE
*
                        DO 180 J = J1, J2, KD1
                           CALL DROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), WORK( J ) )
  180                   CONTINUE
                     END IF
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 190 J = J1, J2, KD1
*
*                    create nonzero element a(j+kd,j-1) outside the
*                    band and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = D( J )*AB( KD1, J )
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           copy off-diagonal elements to E
*
            DO 220 I = 1, N - 1
               E( I ) = AB( 2, I )
  220       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 230 I = 1, N - 1
               E( I ) = ZERO
  230       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 240 I = 1, N
            D( I ) = AB( 1, I )
  240    CONTINUE
      END IF
*

      RETURN
*
*     End of DSBTRD
*
! **********************************************************************************************************************************
12345 format(5X,'Loop ',i8,' of ',i8, a)

! **********************************************************************************************************************************

      END SUBROUTINE DSBTRD_HELPER

! ##################################################################################################################################
! 006 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,
     $                   M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,
     $                   INFO, lowest_mode_num, highest_mode_num  )

! /////////////////////////////////////////////////////////////////////B
      use pentium_ii_kind

      use outa_here_interface

      character(len=len(blnk_sub_nam)):: subr_name = 'DSTEBZ'
! /////////////////////////////////////////////////////////////////////E
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      integer            lowest_mode_num, highest_mode_num
      REAL(DOUBLE)   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      REAL(DOUBLE)   D( * ), E( * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
*  matrix T.  The user may ask for all eigenvalues, all eigenvalues
*  in the half-open interval (VL, VU], or the IL-th through IU-th
*  eigenvalues.
*
*  To avoid overflow, the matrix must be scaled so that its
*  largest element is no greater than overflow**(1/2) *
*  underflow**(1/4) in absolute value, and for greatest
*  accuracy, it should not be much smaller than that.
*
*  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*  Matrix", Report CS41, Computer Science Dept., Stanford
*  University, July 21, 1966.
*
*  Arguments
*  =========
*
*  RANGE   (input) CHARACTER
*          = 'A': ("All")   all eigenvalues will be found.
*          = 'V': ("Value") all eigenvalues in the half-open interval
*                           (VL, VU] will be found.
*          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*                           entire matrix) will be found.
*
*  ORDER   (input) CHARACTER
*          = 'B': ("By Block") the eigenvalues will be grouped by
*                              split-off block (see IBLOCK, ISPLIT) and
*                              ordered from smallest to largest within
*                              the block.
*          = 'E': ("Entire matrix")
*                              the eigenvalues for the entire matrix
*                              will be ordered from smallest to
*                              largest.
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  VL      (input) REAL(DOUBLE)
*  VU      (input) REAL(DOUBLE)
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues.  Eigenvalues less than or equal
*          to VL, or greater than VU, will not be returned.  VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL(DOUBLE)
*          The absolute tolerance for the eigenvalues.  An eigenvalue
*          (or cluster) is considered to be located if it has been
*          determined to lie in an interval whose width is ABSTOL or
*          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
*          will be used, where |T| means the 1-norm of T.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*
*  D       (input) REAL(DOUBLE) array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (input) REAL(DOUBLE) array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*
*  M       (output) INTEGER
*          The actual number of eigenvalues found. 0 <= M <= N.
*          (See also the description of INFO=2,3.)
*
*  NSPLIT  (output) INTEGER
*          The number of diagonal blocks in the matrix T.
*          1 <= NSPLIT <= N.
*
*  W       (output) REAL(DOUBLE) array, dimension (N)
*          On exit, the first M elements of W will contain the
*          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
*          workspace.)
*
*  IBLOCK  (output) INTEGER array, dimension (N)
*          At each row/column j where E(j) is zero or small, the
*          matrix T is considered to split into a block diagonal
*          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
*          block (from 1 to the number of blocks) the eigenvalue W(i)
*          belongs.  (DSTEBZ may use the remaining N-M elements as
*          workspace.)
*
*  ISPLIT  (output) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*          (Only the first NSPLIT elements will actually be used, but
*          since the user cannot know a priori what value NSPLIT will
*          have, N words must be reserved for ISPLIT.)
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  some or all of the eigenvalues failed to converge or
*                were not computed:
*                =1 or 3: Bisection failed to converge for some
*                        eigenvalues; these eigenvalues are flagged by a
*                        negative block number.  The effect is that the
*                        eigenvalues may not be as accurate as the
*                        absolute and relative tolerances.  This is
*                        generally caused by unexpectedly inaccurate
*                        arithmetic.
*                =2 or 3: RANGE='I' only: Not all of the eigenvalues
*                        IL:IU were found.
*                        Effect: M < IU+1-IL
*                        Cause:  non-monotonic arithmetic, causing the
*                                Sturm sequence to be non-monotonic.
*                        Cure:   recalculate, using RANGE='A', and pick
*                                out eigenvalues IL:IU.  In some cases,
*                                increasing the PARAMETER "FUDGE" may
*                                make things work.
*                = 4:    RANGE='I', and the Gershgorin interval
*                        initially used was too small.  No eigenvalues
*                        were computed.
*                        Probable cause: your machine has sloppy
*                                        floating-point arithmetic.
*                        Cure: Increase the PARAMETER "FUDGE",
*                              recompile, and try again.
*
*  Internal Parameters
*  ===================
*
*  RELFAC  REAL(DOUBLE), default = 2.0e0
*          The relative tolerance.  An interval (a,b] lies within
*          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
*          where "ulp" is the machine precision (distance from 1 to
*          the next larger floating point number.)
*
*  FUDGE   REAL(DOUBLE), default = 2
*          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
*          a value of 1 should work, but on machines with sloppy
*          arithmetic, this needs to be larger.  The default for
*          publicly released versions should be large enough to handle
*          the worst machine around.  Note that this has no effect
*          on accuracy of the solution.
*
! lowest_mode_num:  (output) INTEGER mode number of lowest  mode calculated

! highest_mode_num: (output) INTEGER mode number of highest mode calculated

*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
      REAL(DOUBLE)   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO,
     $                   IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX,
     $                   ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL,
     $                   NWU
      REAL(DOUBLE)   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN,
     $                   TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
*     ..
*     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
*     ..
*     .. External Functions ..

      LOGICAL            LSAME
      EXTERNAL           LSAME

      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH

      INTEGER            ILAENV
      EXTERNAL           ILAENV

*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..

*
      INFO = 0
*
*     Decode RANGE
*
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
*
*     Decode ORDER
*
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
*
*     Check for Errors
*
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( IORDER.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.2 ) THEN
         IF( VL.GE.VU )
     $      INFO = -5
      ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $          THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )
     $          THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEBZ', -INFO )
         RETURN
      END IF
*
*     Initialize error flags
*
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
*     Simplifications:
*
      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N )
     $   IRANGE = 1
*
*     Get machine constants
*     NB is the minimum vector length for vector bisection, or 0
*     if only scalar is to be done.
*
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 )
     $   NB = 0
*
*     Special Case when N=1
*
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
*
*     Compute Splitting Points
*
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE
*
*DIR$ NOVECTOR
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
*
*     Compute Interval and ATOLI
*
      IF( IRANGE.EQ.3 ) THEN
*
*        RANGE='I': Compute the interval containing eigenvalues
*                   IL through IU.
*
*        Compute Gershgorin interval for entire (split) matrix
*        and use it as the initial interval
*
         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO
*
         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
   20    CONTINUE
*
         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
*
*        Compute Iteration parameters
*
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
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
*
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
     $                WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
     $                IWORK, W, IBLOCK, IINFO )
! /////////////////////////////////////////////////////////////////////B
         if (iinfo > 0) then
            if ((iinfo >=1) .and. (iinfo <=2)) then
               Write(err,801) subr_name,iinfo
               Write(f06,801) subr_name,iinfo
               fatal_err = fatal_err + 1
            else
               Write(err,802) subr_name,iinfo
               Write(f06,802) subr_name,iinfo
               fatal_err = fatal_err + 1
            endif
         endif
! /////////////////////////////////////////////////////////////////////E
*
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
*
         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
*
*        RANGE='A' or 'V' -- Set ATOLI
*
         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $           ABS( D( N ) )+ABS( E( N-1 ) ) )
*
         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+
     $              ABS( E( J ) ) )
   30    CONTINUE
*
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
         IF( IRANGE.EQ.2 ) THEN
            WL = VL
            WU = VU
         ELSE
            WL = ZERO
            WU = ZERO
         END IF
      END IF
*
*     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
*     NWL accumulates the number of eigenvalues .le. WL,
*     NWU accumulates the number of eigenvalues .le. WU
*
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
*
      write(sc1,*)
      DO 70 JB = 1, NSPLIT
         write(sc1,12345,advance='no') jb,nsplit,cr13_lge
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
*
         IF( IN.EQ.1 ) THEN
*
*           Special Case -- IN=1
*
            IF( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN )
     $         NWL = NWL + 1
            IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN )
     $         NWU = NWU + 1
            IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE.
     $          D( IBEGIN )-PIVMIN ) ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            END IF
         ELSE
*
*           General Case -- IN > 1
*
*           Compute Gershgorin Interval
*           and use it as the initial interval
*
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO
*
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE
*
            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
*
*           Compute ATOLI for the current submatrix
*
            IF( ABSTOL.LE.ZERO ) THEN
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            ELSE
               ATOLI = ABSTOL
            END IF
*
            IF( IRANGE.GT.1 ) THEN
               IF( GU.LT.WL ) THEN
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               END IF
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU )
     $            GO TO 70
            END IF
*
*           Set Up Initial Interval
*
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
*
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )
*
*           Compute Eigenvalues
*
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) /
     $              LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )

! /////////////////////////////////////////////////////////////////////B
            if (iinfo > 0) then
               if ((iinfo >=1) .and. (iinfo <=in)) then
                  Write(err,801) subr_name,iinfo
                  Write(f06,801) subr_name,iinfo
                  fatal_err = fatal_err + 1
               else
                  Write(err,802) subr_name,iinfo
                  Write(f06,802) subr_name,iinfo
                  fatal_err = fatal_err + 1
               endif
            endif
! /////////////////////////////////////////////////////////////////////E
*
*           Copy Eigenvalues Into W and IBLOCK
*           Use -JB for block number for unconverged eigenvalues.
*
            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
*
*              Flag non-convergence.
*
               IF( J.GT.IOUT-IINFO ) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF,
     $                 IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
*
            M = M + IM
         END IF
   70 CONTINUE
*
*     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
*     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
*
      IF( IRANGE.EQ.3 ) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
*
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
*
*           Code to deal with effects of bad arithmetic:
*           Some low eigenvalues to be discarded are not in (WL,WLU],
*           or high eigenvalues to be discarded are not in (WUL,WU]
*           so just kill off the smallest IDISCL/largest IDISCU
*           eigenvalues, by simply finding the smallest/largest
*           eigenvalue(s).
*
*           (If N(w) is monotone non-decreasing, this should never
*               happen.)
*
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
*
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
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
*
*     If ORDER='B', do nothing -- the eigenvalues are already sorted
*        by block.
*     If ORDER='E', sort the eigenvalues from smallest to largest
*
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
*
            IF( IE.NE.0 ) THEN
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            END IF
  150    CONTINUE
      END IF
*
      highest_mode_num = nwu
      lowest_mode_num = highest_mode_num - m + 1

      INFO = 0
      IF( NCNVRG )
     $   INFO = INFO + 1
      IF( TOOFEW )
     $   INFO = INFO + 2


      RETURN
*
*     End of DSTEBZ
*
! **********************************************************************************************************************************
12345 format(5X,'Block ',i8,' of ',i8, a)

  801 format(' *ERROR   801: ERROR ATTEMPTING TO SOLVE FOR EIGENVALUES',
     &'. IN SUBROUTINE DLAEBZ, (CALLED FROM SUBR ',A,')'
     &   ,/,14X,' THE LAST ',I8,' INTERVALS DID NOT CONVERGE')

  802 format(' *ERROR   802: ERROR ATTEMPTING TO SOLVE FOR EIGENVALUES',
     &'. IN SUBROUTINE DLAEBZ, (CALLED FROM SUBR ',A,')'
     &   ,/,14X,' TOO MANY INTERVALS WERE GENERATED. MAX IS ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE DSTEBZ

! ##################################################################################################################################
! 007 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
     $                   IWORK, IFAIL, INFO )
*
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DSTEIN'

*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
     $                   IWORK( * )
      REAL(DOUBLE)   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEIN computes the eigenvectors of a real symmetric tridiagonal
*  matrix T corresponding to specified eigenvalues, using inverse
*  iteration.
*
*  The maximum number of iterations allowed for each eigenvector is
*  specified by an internal parameter MAXITS (currently set to 5).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input) REAL(DOUBLE) array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (input) REAL(DOUBLE) array, dimension (N)
*          The (n-1) subdiagonal elements of the tridiagonal matrix
*          T, in elements 1 to N-1.  E(N) need not be set.
*
*  M       (input) INTEGER
*          The number of eigenvectors to be found.  0 <= M <= N.
*
*  W       (input) REAL(DOUBLE) array, dimension (N)
*          The first M elements of W contain the eigenvalues for
*          which eigenvectors are to be computed.  The eigenvalues
*          should be grouped by split-off block and ordered from
*          smallest to largest within the block.  ( The output array
*          W from DSTEBZ with ORDER = 'B' is expected here. )
*
*  IBLOCK  (input) INTEGER array, dimension (N)
*          The submatrix indices associated with the corresponding
*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
*          the first submatrix from the top, =2 if W(i) belongs to
*          the second submatrix, etc.  ( The output array IBLOCK
*          from DSTEBZ is expected here. )
*
*  ISPLIT  (input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to
*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*          through ISPLIT( 2 ), etc.
*          ( The output array ISPLIT from DSTEBZ is expected here. )
*
*  Z       (output) REAL(DOUBLE) array, dimension (LDZ, M)
*          The computed eigenvectors.  The eigenvector associated
*          with the eigenvalue W(i) is stored in the i-th column of
*          Z.  Any vector which fails to converge is set to its current
*          iterate after MAXITS iterations.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= max(1,N).
*
*  WORK    (workspace) REAL(DOUBLE) array, dimension (5*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  IFAIL   (output) INTEGER array, dimension (M)
*          On normal exit, all elements of IFAIL are zero.
*          If one or more eigenvectors fail to converge after
*          MAXITS iterations, then their indices are stored in
*          array IFAIL.
*
*  INFO    (output) INTEGER
*          = 0: successful exit.
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, then i eigenvectors failed to converge
*               in MAXITS iterations.  Their indices are stored in
*               array IFAIL.
*
*  Internal Parameters
*  ===================
*
*  MAXITS  INTEGER, default = 5
*          The maximum number of iterations performed.
*
*  EXTRA   INTEGER, default = 2
*          The number of iterations performed after norm growth
*          criterion is satisfied, should be at least 1.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1,
     $                   ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1,
     $                   INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1,
     $                   JBLK, JMAX, NBLK, NRMCHK
      REAL(DOUBLE)   DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL,
     $                   SCL, SEP, TOL, XJ, XJM, ZTR
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. External Functions ..

      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*

*     Test the input parameters.
*
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) )
     $           THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEIN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      EPS = DLAMCH( 'Precision' )
*
*     Initialize seed for random number generator DLARNV.
*
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
*
*     Initialize pointers.
*
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
*
*     Compute eigenvectors of matrix blocks.
*
      J1 = 1
      DO 160 NBLK = 1, IBLOCK( M )
*
*        Find starting and ending indices of block nblk.
*
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 )
     $      GO TO 60
         GPIND = B1
*
*        Compute reorthogonalization criterion and stopping criterion.
*
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+
     $               ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
*
         DTPCRT = SQRT( ODM1 / BLKSIZ )
*
*        Loop through eigenvalues of block nblk.
*
   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 160
            END IF
            JBLK = JBLK + 1
            XJ = W( J )
*
*           Skip all the work if the block size is one.
*
            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            END IF
*
*           If eigenvalues j and j-1 are too close, add a relatively
*           small perturbation.
*
            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL )
     $            XJ = XJM + PERTOL
            END IF
*
            ITS = 0
            NRMCHK = 0
*
*           Get random starting vector.
*
            CALL DLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
*
*           Copy the matrix T so it won't be destroyed in factorization.
*
            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
*
*           Compute LU factors with partial pivoting  ( PT = LU )
*
            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK,
     $                   IINFO )
*
*           Update iteration count.
*
   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS )
     $         GO TO 100
*
*           Normalize and scale the righthand side vector Pb.
*
            SCL = BLKSIZ*ONENRM*MAX( EPS,
     $            ABS( WORK( INDRV4+BLKSIZ ) ) ) /
     $            DASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
*
*           Solve the system LU = Pb.
*
            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK,
     $                   WORK( INDRV1+1 ), TOL, IINFO )
*
*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
*           close enough.
*
            IF( JBLK.EQ.1 )
     $         GO TO 90
            IF( ABS( XJ-XJM ).GT.ORTOL )
     $         GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ),
     $                  1 )
                  CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1,
     $                        WORK( INDRV1+1 ), 1 )
   80          CONTINUE
            END IF
*
*           Check the infinity norm of the iterate.
*
   90       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
*
*           Continue for additional iterations after norm reaches
*           stopping criterion.
*
            IF( NRM.LT.DTPCRT )
     $         GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 )
     $         GO TO 70
*
            GO TO 110
*
*           If stopping criterion was not satisfied, update info and
*           store eigenvector number in array ifail.
*
  100       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
*
*           Accept iterate as jth eigenvector.
*
  110       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO )
     $         SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  120       CONTINUE
            DO 130 I = 1, N
               Z( I, J ) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z( B1+I-1, J ) = WORK( INDRV1+I )
  140       CONTINUE
*
*           Save the shift to check eigenvalue spacing at next
*           iteration.
*
            XJM = XJ
*
  150    CONTINUE
  160 CONTINUE
*
      RETURN
*
*     End of DSTEIN
*
      END SUBROUTINE DSTEIN

! ##################################################################################################################################
! 008 LAPACK_GIV_MGIV_EIG

      SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DLAGTF'
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
      REAL(DOUBLE)   LAMBDA, TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      REAL(DOUBLE)   A( * ), B( * ), C( * ), D( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n
*  tridiagonal matrix and lambda is a scalar, as
*
*     T - lambda*I = PLU,
*
*  where P is a permutation matrix, L is a unit lower tridiagonal matrix
*  with at most one non-zero sub-diagonal elements per column and U is
*  an upper triangular matrix with at most two non-zero super-diagonal
*  elements per column.
*
*  The factorization is obtained by Gaussian elimination with partial
*  pivoting and implicit row scaling.
*
*  The parameter LAMBDA is included in the routine so that DLAGTF may
*  be used, in conjunction with DLAGTS, to obtain eigenvectors of T by
*  inverse iteration.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix T.
*
*  A       (input/output) REAL(DOUBLE) array, dimension (N)
*          On entry, A must contain the diagonal elements of T.
*
*          On exit, A is overwritten by the n diagonal elements of the
*          upper triangular matrix U of the factorization of T.
*
*  LAMBDA  (input) REAL(DOUBLE)
*          On entry, the scalar lambda.
*
*  B       (input/output) REAL(DOUBLE) array, dimension (N-1)
*          On entry, B must contain the (n-1) super-diagonal elements of
*          T.
*
*          On exit, B is overwritten by the (n-1) super-diagonal
*          elements of the matrix U of the factorization of T.
*
*  C       (input/output) REAL(DOUBLE) array, dimension (N-1)
*          On entry, C must contain the (n-1) sub-diagonal elements of
*          T.
*
*          On exit, C is overwritten by the (n-1) sub-diagonal elements
*          of the matrix L of the factorization of T.
*
*  TOL     (input) REAL(DOUBLE)
*          On entry, a relative tolerance used to indicate whether or
*          not the matrix (T - lambda*I) is nearly singular. TOL should
*          normally be chose as approximately the largest relative error
*          in the elements of T. For example, if the elements of T are
*          correct to about 4 significant figures, then TOL should be
*          set to about 5*10**(-4). If TOL is supplied as less than eps,
*          where eps is the relative machine precision, then the value
*          eps is used in place of TOL.
*
*  D       (output) REAL(DOUBLE) array, dimension (N-2)
*          On exit, D is overwritten by the (n-2) second super-diagonal
*          elements of the matrix U of the factorization of T.
*
*  IN      (output) INTEGER array, dimension (N)
*          On exit, IN contains details of the permutation matrix P. If
*          an interchange occurred at the kth step of the elimination,
*          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)
*          returns the smallest positive integer j such that
*
*             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,
*
*          where norm( A(j) ) denotes the sum of the absolute values of
*          the jth row of the matrix A. If no such j exists then IN(n)
*          is returned as zero. If IN(n) is returned as positive, then a
*          diagonal element of U is small, indicating that
*          (T - lambda*I) is singular or nearly singular,
*
*  INFO    (output) INTEGER
*          = 0   : successful exit
*          .lt. 0: if INFO = -k, the kth argument had an illegal value
*
* =====================================================================
*
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            K
      REAL(DOUBLE)   EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
*     ..
*     .. Executable Statements ..

*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLAGTF', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      A( 1 ) = A( 1 ) - LAMBDA
      IN( N ) = 0
      IF( N.EQ.1 ) THEN
         IF( A( 1 ).EQ.ZERO )
     $      IN( 1 ) = 1
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
*
      TL = MAX( TOL, EPS )
      SCALE1 = ABS( A( 1 ) ) + ABS( B( 1 ) )
      DO 10 K = 1, N - 1
         A( K+1 ) = A( K+1 ) - LAMBDA
         SCALE2 = ABS( C( K ) ) + ABS( A( K+1 ) )
         IF( K.LT.( N-1 ) )
     $      SCALE2 = SCALE2 + ABS( B( K+1 ) )
         IF( A( K ).EQ.ZERO ) THEN
            PIV1 = ZERO
         ELSE
            PIV1 = ABS( A( K ) ) / SCALE1
         END IF
         IF( C( K ).EQ.ZERO ) THEN
            IN( K ) = 0
            PIV2 = ZERO
            SCALE1 = SCALE2
            IF( K.LT.( N-1 ) )
     $         D( K ) = ZERO
         ELSE
            PIV2 = ABS( C( K ) ) / SCALE2
            IF( PIV2.LE.PIV1 ) THEN
               IN( K ) = 0
               SCALE1 = SCALE2
               C( K ) = C( K ) / A( K )
               A( K+1 ) = A( K+1 ) - C( K )*B( K )
               IF( K.LT.( N-1 ) )
     $            D( K ) = ZERO
            ELSE
               IN( K ) = 1
               MULT = A( K ) / C( K )
               A( K ) = C( K )
               TEMP = A( K+1 )
               A( K+1 ) = B( K ) - MULT*TEMP
               IF( K.LT.( N-1 ) ) THEN
                  D( K ) = B( K+1 )
                  B( K+1 ) = -MULT*D( K )
               END IF
               B( K ) = TEMP
               C( K ) = MULT
            END IF
         END IF
         IF( ( MAX( PIV1, PIV2 ).LE.TL ) .AND. ( IN( N ).EQ.0 ) )
     $      IN( N ) = K
   10 CONTINUE
      IF( ( ABS( A( N ) ).LE.SCALE1*TL ) .AND. ( IN( N ).EQ.0 ) )
     $   IN( N ) = N
*
      RETURN

*
*     End of DLAGTF
*
      END SUBROUTINE DLAGTF

! --- lapack_surgery end --- !
      END MODULE LAPACK_GIV_MGIV_EIG
