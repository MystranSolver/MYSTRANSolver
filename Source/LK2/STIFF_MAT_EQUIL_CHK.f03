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

   MODULE REDUCTION_CHECKS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: STIFF_MAT_EQUIL_CHK

   CONTAINS

      SUBROUTINE STIFF_MAT_EQUIL_CHK ( OUTPUT, X_SET, SYM_KIN, NROWS, NTERM_KIN, I_KIN, J_KIN, KIN, KIN_DIAG, KIN_MAX_DIAG, RBMAT )

! Performs a stiffness matrix equilibrium check on input matrix KIN by calculating KIN*RBMAT = PRB where RBMAT is a rigid body
! displacement matrix (NROWS x 6). Each column of RBMAT represents a rigid body displacement for one of the 6 components of motion:
! T1, T2, T3, R1, R2, R3. In order for the check to make any sense, KIN must be not have any rigid body modes restrained (as for
! example would be the case if it were grounded - e.g. a cantilevered beam has rigid body modes restrained)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, SC1

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSPOINT, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES, ONLY            :  TDOFI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND
      USE PARAMS, ONLY                :  EPSIL, EQCHK_NORM, SUPWARN, SUPINFO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE LAPACK_ADAPTERS, ONLY       :  BANDGEN_LAPACK_DPB, BANDSIZ
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_LAPACK_MAT, DEALLOCATE_LAPACK_MAT
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_MATRIX_BY_ROWS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'STIFF_MAT_EQUIL_CHK'
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_KIN             ! Whether KIN is stored as symmetric (only nonzero terms on and above
!                                                              the diag) or as nonsym (all nonzero terms). Passed as arg in a CALL
      CHARACTER(LEN=*), INTENT(IN)    :: X_SET               ! The displ set for the equil check (e,g, 'G ', etc
      CHARACTER( 2*BYTE)              :: DIR(6)      =         (/'T1','T2','T3','R1','R2','R3'/)
      CHARACTER( 3*BYTE)              :: MAT_NAME            ! Character name for output purposes
      CHARACTER(20*BYTE)              :: MAT_DESCR           ! Character name for output purposes
      CHARACTER(LEN=LEN(EQCHK_NORM))  :: NORM                ! Character name for output purposes

      INTEGER(LONG), INTENT(IN)       :: NROWS               ! Number of rows in KIN
      INTEGER(LONG), INTENT(IN)       :: NTERM_KIN           ! Number of nonzero terms in KIN
      INTEGER(LONG), INTENT(IN)       :: I_KIN(NROWS+1)      ! Row start indices for KIN
      INTEGER(LONG), INTENT(IN)       :: J_KIN(NTERM_KIN)    ! Col numbers of terms in KIN
      INTEGER(LONG), INTENT(IN)       :: OUTPUT              ! =1, output PRB, =2 output RB_STRN_ENRGY, =3 output both
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: KIN_SDIA            ! No. of superdiags in KIN upper triangle
      INTEGER(LONG)                   :: ROW,COL             ! Row/col where max term in RB_STRAIN_ENERGY exists
      INTEGER(LONG)                   :: SA_SET_COL          ! Col no. in array TDOF where the  SA-set is (from subr TDOF_COL_NUM)


      CHARACTER( 1*BYTE)              :: FLAG(NROWS)         ! Character to designate whether PRBN was normalized to diag KIN or not

      REAL(DOUBLE), INTENT(IN)        :: KIN(NTERM_KIN)      ! Nonzero terms in KIN
      REAL(DOUBLE), INTENT(IN)        :: KIN_DIAG(NROWS)     ! Diagonal of KIN
      REAL(DOUBLE), INTENT(IN)        :: KIN_MAX_DIAG        ! Max diag term from KIN
      REAL(DOUBLE), INTENT(IN)        :: RBMAT(NROWS,6)      ! Rigid body displacement matrix (6 rigid body modes)
      REAL(DOUBLE)                    :: EPS1                ! Small number to compare real zero to
      REAL(DOUBLE)                    :: INV_KIN_DIAGI       ! 1.0/KIN_DIAG(I) or 1.0/KIN_MAX_DIAG
      REAL(DOUBLE)                    :: MAX_ABS             ! Max abs value from RB_STRAIN_ENERGY
      REAL(DOUBLE)                    :: PRB(NROWS,6)        ! KIN*RBMAT (forces resulting from equil check on KIN)
      REAL(DOUBLE)                    :: PRB_COL(NROWS)      ! One column of PRB
      REAL(DOUBLE)                    :: PRBN(NROWS,6)       ! KIN*RBMAT w/ i-th row normalized on i-th diag term of KIN
      REAL(DOUBLE)                    :: RBMAT_COL(NROWS)    ! One row of RBMAT
      REAL(DOUBLE)                    :: RB_STRN_ENRGY(6,6)  ! Strain energy for the RBGLOBAL displs

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! If NSPOINT > 0 MYSTRAN code, as of 05/05/07, will not accomodate equil check on models with any SPOINT's

      IF (NSPOINT > 0) THEN
         WARN_ERR = WARN_ERR + 1
         IF (DEBUG(23) == 0) THEN
            WRITE(ERR,1001) NSPOINT, X_SET
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1001) NSPOINT, X_SET
            ENDIF
            RETURN
         ELSE
            WRITE(ERR,1002) NSPOINT, X_SET
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1002) NSPOINT, X_SET
            ENDIF
            CONTINUE
         ENDIF
      ENDIF

      EPS1 = EPSIL(1)

! Make sure that the max diagonal term in KIN is not zero. If so, we cannot normalize to get PRBN

      IF (DABS(KIN_MAX_DIAG) > EPS1) THEN
         NORM = EQCHK_NORM
      ELSE
         NORM = 'N'
         IF (EQCHK_NORM == 'Y') THEN
            WRITE(ERR,101) X_SET, KIN_MAX_DIAG
            IF (SUPINFO == 'N') THEN
               WRITE(F06,101) X_SET, KIN_MAX_DIAG
            ENDIF
         ENDIF
      ENDIF

      IF      (X_SET == 'G ') THEN
         MAT_NAME = 'KGG'
         MAT_DESCR = 'STIFFNESS MATRIX KGG'
      ELSE IF (X_SET == 'N ') THEN
         MAT_NAME = 'KNN'
         MAT_DESCR = 'STIFFNESS MATRIX KNN'
      ELSE IF (X_SET == 'F ') THEN
         MAT_NAME = 'KFF'
         MAT_DESCR = 'STIFFNESS MATRIX KFF'
      ELSE IF (X_SET == 'A ') THEN
         MAT_NAME = 'KAA'
         MAT_DESCR = 'STIFFNESS MATRIX KAA'
      ELSE IF (X_SET == 'L ') THEN
         MAT_NAME = 'KLL'
         MAT_DESCR = 'STIFFNESS MATRIX KLL'
      ELSE
         WRITE(ERR,970) SUBR_NAME, X_SET
         WRITE(F06,970) SUBR_NAME, X_SET
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Write RBMAT out, if requested

      IF (DEBUG(22) == 1) THEN
         CALL WRITE_RB_MATS ( 'DISPL' )
      ENDIF

! Calculate grid forces due to rigid body motion

      IF (OUTPUT > 0) THEN

         IF (DEBUG(21) == 0) THEN                          ! Use MATMULT_SFF to mult KIN*RBMAT
            CALL MATMULT_SFF ( MAT_NAME, NROWS, NROWS, NTERM_KIN, SYM_KIN, I_KIN, J_KIN, KIN, 'RBMAT', NROWS, 6, RBMAT, 'Y', 'PRB',&
                               ONE, PRB )

         ELSE                                              ! Use DSBMV to mult KIN*RBMAT

            WRITE(SC1,2191) MAT_NAME
            CALL BANDSIZ ( NROWS, NTERM_KIN, I_KIN, J_KIN, KIN_SDIA )

            WRITE(SC1,2192) MAT_NAME
            CALL ALLOCATE_LAPACK_MAT ( 'ABAND', KIN_SDIA+1, NROWS, SUBR_NAME )

            WRITE(SC1,2193) MAT_NAME
            CALL BANDGEN_LAPACK_DPB ( MAT_NAME, NROWS, KIN_SDIA, NTERM_KIN, I_KIN, J_KIN, KIN, ABAND, SUBR_NAME )

            IF ((DEBUG(34) == 1) .OR. (DEBUG(34) == 3)) THEN
               CALL WRITE_MATRIX_BY_ROWS ( MAT_DESCR, ABAND, KIN_SDIA+1, NROWS, F06 )
            ENDIF
            WRITE(SC1,2194)
            WRITE(SC1, * )
            DO J=1,6
               DO I=1,NROWS
                  RBMAT_COL(I) = RBMAT(I,J)
                  PRB_COL(I)   = ZERO
               ENDDO
               CALL DSBMV ( 'U', NROWS, KIN_SDIA, 1.0D0, ABAND, KIN_SDIA+1, RBMAT_COL, 1, 0.0D0, PRB_COL, 1 )
               DO I=1,NROWS
                  PRB(I,J) = PRB_COL(I)
               ENDDO
            ENDDO

            CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )

         ENDIF

         CALL TDOF_COL_NUM ( 'SA',  SA_SET_COL )
         IF (NORM == 'Y') THEN
            DO I=1,NROWS
               IF ((DABS(KIN_DIAG(I)) > EPS1) .AND. (TDOFI(I,SA_SET_COL) == 0)) THEN
                  INV_KIN_DIAGI = ONE/KIN_DIAG(I)
                  FLAG(I) = ' '
               ELSE
                  INV_KIN_DIAGI = ONE/KIN_MAX_DIAG
                  FLAG(I) = '*'
               ENDIF
               DO J=1,6
                  PRBN(I,J) = PRB(I,J)*INV_KIN_DIAGI
               ENDDO
            ENDDO
         ENDIF
      ENDIF

! Write grid forces if requested

      IF ((OUTPUT == 1) .OR. (OUTPUT == 3)) THEN
         CALL WRITE_RB_MATS ( 'FORCE' )
      ENDIF

! Output strain energy if requested

      IF ((OUTPUT == 2) .OR. (OUTPUT == 3)) THEN

         CALL MATMULT_FFF_T ( RBMAT, PRB, NROWS, 6, 6, RB_STRN_ENRGY )

         MAX_ABS = ZERO
         ROW = 0
         COL = 0
         DO I=1,6
            DO J=1,6
               IF(DABS(RB_STRN_ENRGY(I,J)) > MAX_ABS) THEN
                  MAX_ABS = DABS(RB_STRN_ENRGY(I,J))
                  ROW = I
                  COL = J
               ENDIF
            ENDDO
         ENDDO

         IF (MAX_ABS == ZERO) THEN
            ROW = 1
            COL = 1
            MAX_ABS = DABS(RB_STRN_ENRGY(1,1))
         ENDIF
            WRITE(F06,2107) X_SET, MAX_ABS, ROW, COL

         DO I=1,6
            WRITE(F06,2108) DIR(I), (RB_STRN_ENRGY(I,J),J=1,6)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT('*INFORMATION: CANNOT NORMALIZE THE EQUILIBRIUM CHECK FOR THE ',A,'-SET. THE MAX DIAGONAL OF THAT STIFFNESS MATRIX IS'&
                   ,1ES15.6,' AND IS TOO SMALL'                                                                                    &
                   ,/,14X,' THE UNNORMALIZED EQUILIBRIUM CHECK WILL BE OUTPUT INSTEAD')

  970 FORMAT(' *ERROR   970: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' REQUEST FOR EQUILIBRIUM CHECK ON THE ',A,' SET CANNOT BE DONE. THE INPUT SET MUST BE G, N, F OR A')

 1001 FORMAT(' *WARNING    : DUE TO THE PRESENCE OF ',I8,' SPOINT''s, MYSTRAN CANNOT PERFORM THE REQUESTED EQUIL CHECK ON THE ',A, &
                           ' SET STIFF MATRIX.'                                                                                    &
                    ,/,14X,' USER CAN OVERRIDE THIS WITH B.D ENTRY DEBUG 23 > 0')

 1002 FORMAT(' *WARNING    : EVEN THOUGH THERE ARE ',I8,' SPOINT''s IN THE MODEL, THE USER HAS REQUESTED A STIFF MATRIX EQUIL',    &
                           ' CHECK FOR THE ',A,' SET.'                                                                             &
                    ,/,14X,' THE RESULTS MAY BE UNPREDICTABLE DUE TO THE FACT THAT SPOINT FORCES ARE ALWAYS REPORTED AS COMP T1',/)

 2107 FORMAT(33X,' STRAIN ENERGY DUE TO RBGLOGAL ',A,'-SET RIGID BODY DISPLACEMENTS',/,                                            &
             42X,'(max absolute value = ',1ES9.2,' at row',I2,', col',I2,')'                                                       &
                  ,/,26X,'T1             T2             T3             R1             R2             R3',/)

 2108 FORMAT(15X,A,1X,6(1ES15.6))

 2191 FORMAT(8X,'CALCULATE BANDWIDTH OF ',A,' MATRIX')

 2192 FORMAT(8X,'ALLOCATE ARRAYS FOR LAPACK BAND FORM OF ',A)

 2193 FORMAT(8X,'PUT ',A,' MATRIX IN LAPACK BAND FORM')

 2194 FORMAT(8X,'CALC FORCES DUE TO R.B. DISPLS FOR EACH OF THE 6 COLS OF RBGLOBAL')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE WRITE_RB_MATS ( WHAT )

! Write rigid body matrices (RB displ, forces due to RB disp)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  NDOFG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOFI
      USE PARAMS, ONLY                :  EQCHK_REF_GRID, EQCHK_TINY
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      IMPLICIT NONE

      CHARACTER(LEN=* )               :: WHAT              ! What to write to F06

      INTEGER(LONG)                   :: OLD_GRID          ! Grid ID (for output deliniation)
      INTEGER(LONG)                   :: X_SET_COL         ! Col no. in array TDOFI where the X_SET is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: X_SET_COMP        ! X_SET component number
      INTEGER(LONG)                   :: X_SET_DOF         ! X_SET DOF number
      INTEGER(LONG)                   :: X_SET_GRID        ! X_SET grid number

      REAL(DOUBLE)                    :: AFORCE(6)         ! Absolute value of PRB or PRBN forces in one row
      REAL(DOUBLE)                    :: MAX_AFORCE(6)     ! largest value of AFORCE(J) for J=1,6

      INTRINSIC                       :: DABS

! **********************************************************************************************************************************
      CALL TDOF_COL_NUM ( X_SET, X_SET_COL )

      IF (WHAT == 'DISPL') THEN

         IF (EQCHK_REF_GRID == 0) THEN
            WRITE(F06,2901) X_SET
         ELSE
            WRITE(F06,2902) X_SET, EQCHK_REF_GRID
         ENDIF
         OLD_GRID = 0
         DO I=1,NDOFG
            X_SET_DOF  = TDOFI(I,X_SET_COL)
            IF (X_SET_DOF > 0) THEN
               X_SET_GRID = TDOFI(I,1)
               X_SET_COMP = TDOFI(I,2)
               IF (DEBUG(88) == 0) THEN
                  IF (OLD_GRID /= X_SET_GRID) THEN
                     WRITE(F06,*)
                  ENDIF
               ENDIF
               OLD_GRID = X_SET_GRID
               WRITE(F06,2904) X_SET_GRID, X_SET_COMP, (RBMAT(X_SET_DOF,J),J=1,6)
            ENDIF
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == 'FORCE') THEN

         IF (EQCHK_NORM == 'N') THEN
            WRITE(F06,2911) X_SET, MAT_NAME, X_SET
         ELSE
            WRITE(F06,2912) X_SET, MAT_NAME, X_SET
         ENDIF
         OLD_GRID = 0
         DO I=1,6
!xx         MAX_AFORCE(I) = ZERO
            MAX_AFORCE(I) = -MACH_LARGE_NUM
         ENDDO
         DO I=1,NDOFG
            X_SET_DOF  = TDOFI(I,X_SET_COL)
            IF (X_SET_DOF > 0) THEN
               X_SET_GRID = TDOFI(I,1)
               X_SET_COMP = TDOFI(I,2)
               IF (EQCHK_NORM == 'Y') THEN
                  DO J=1,6
                     AFORCE(J) = DABS(PRBN(X_SET_DOF,J))
                     IF (AFORCE(J) > MAX_AFORCE(J)) THEN
                        MAX_AFORCE(J) = AFORCE(J)
                     ENDIF
                  ENDDO
                  IF ((AFORCE(1) > EQCHK_TINY) .OR. (AFORCE(2) > EQCHK_TINY) .OR. (AFORCE(3) > EQCHK_TINY) .OR.                    &
                      (AFORCE(4) > EQCHK_TINY) .OR. (AFORCE(5) > EQCHK_TINY) .OR. (AFORCE(6) > EQCHK_TINY)) THEN
                      IF (DEBUG(88) == 0) THEN
                        IF (OLD_GRID /= X_SET_GRID) THEN
                           WRITE(F06,*)
                        ENDIF
                     ENDIF
                     OLD_GRID = X_SET_GRID
                     WRITE(F06,2903) X_SET_GRID, X_SET_COMP, FLAG(X_SET_DOF), (PRBN(X_SET_DOF,J),J=1,6)
                  ENDIF
               ELSE
                  DO J=1,6
                     AFORCE(J) = DABS(PRB(X_SET_DOF,J))
                     IF (AFORCE(J) > MAX_AFORCE(J)) THEN
                        MAX_AFORCE(J) = AFORCE(J)
                     ENDIF
                  ENDDO
                  IF ((AFORCE(1) > EQCHK_TINY) .OR. (AFORCE(2) > EQCHK_TINY) .OR. (AFORCE(3) > EQCHK_TINY) .OR.                    &
                      (AFORCE(4) > EQCHK_TINY) .OR. (AFORCE(5) > EQCHK_TINY) .OR. (AFORCE(6) > EQCHK_TINY)) THEN
                      IF (DEBUG(88) == 0) THEN
                        IF (OLD_GRID /= X_SET_GRID) THEN
                           WRITE(F06,*)
                        ENDIF
                     ENDIF
                     OLD_GRID = X_SET_GRID
                     WRITE(F06,2904) X_SET_GRID, X_SET_COMP, (PRB(X_SET_DOF,J) ,J=1,6)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,2905)
         WRITE(F06,2906) (MAX_AFORCE(J),J=1,6)

         IF (EQCHK_NORM == 'Y') THEN
            WRITE(F06,2907) KIN_MAX_DIAG
         ENDIF

      ENDIF

! **********************************************************************************************************************************
 2901 FORMAT(38X,A,'SET RBGLOBAL RIGID BODY DISPLACEMENTS'                                                                         &
                  ,/,46X,'RELATIVE TO BASIC ORIGIN'                                                                                &
                  ,/,26X,'T1             T2             T3             R1             R2             R3')

 2902 FORMAT(38X,A,'SET RBGLOBAL RIGID BODY DISPLACEMENTS'                                                                         &
                  ,/,46X,'RELATIVE TO GRID ',I8                                                                                    &
                  ,/,26X,'T1             T2             T3             R1             R2             R3')

 2911 FORMAT(6X,A,'SET EQUILIBRIUM CHECK: ',A,'*RBGLOBAL = ',A,'SET GRID POINT FORCES DUE TO RBGLOBAL RIGID BODY DISPLACEMENTS'    &
                  ,/,46X,'UNNORMALIZED GRID FORCES',/,46X,'(should be all zero''s)'                                                &
                  ,/,26X,'T1             T2             T3             R1             R2             R3')

 2912 FORMAT(6X,A,'SET EQUILIBRIUM CHECK: ',A,'*RBGLOBAL = ',A,'SET GRID POINT FORCES DUE TO RBGLOBAL RIGID BODY DISPLACEMENTS'    &
                  ,/,36X,'GRID FORCES NORMALIZED TO DIAGONAL STIFFNESS',/,46X,'(should be all zero''s)'                            &
                  ,/,26X,'T1             T2             T3             R1             R2             R3')

 2903 FORMAT(1X,2I8,A1,6(1ES15.6))

 2904 FORMAT(1X,2I8,1X,6(1ES15.6))

 2905 FORMAT(20X,'-------------  -------------  -------------  -------------  -------------  -------------')

 2906 FORMAT(7X,'ABS MAX',4X,6(1ES15.6),/)

 2907 FORMAT(7X,'* stiffness for row is zreo, so row was normalized to the max diagonal term = ',1ES13.6,//)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_RB_MATS

      END SUBROUTINE STIFF_MAT_EQUIL_CHK

   END MODULE REDUCTION_CHECKS
