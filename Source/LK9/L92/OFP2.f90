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

      SUBROUTINE OFP2 ( JVEC, WHAT, SC_OUT_REQ, ZERO_GEN_STIFF, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )

      ! Processes SPC and MPC force output requests for 1 subcase.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, OT4

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, GROUT_SPCF_BIT, GROUT_MPCF_BIT, GROUT_GPFO_BIT, IBIT, INT_SC_NUM,&
                                         MELGP, MOGEL, NGRID, NDOFF, NDOFG, NDOFM, NDOFN, NDOFS, NDOFSA, NTERM_GMN,                &
                                         NTERM_HMN, NTERM_KFS, NTERM_KFSD, NTERM_LMN, NTERM_MFS, NTERM_QS, SOL_NAME

      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START, TDOFI
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, GEN_MASS, MEFFMASS, MPFACTOR_N6
      USE MODEL_STUF, ONLY            :  ANY_SPCF_OUTPUT, ANY_MPCF_OUTPUT, GRID, GRID_ID, GROUT, MEFFMASS_CALC, MPFACTOR_CALC
      USE PARAMS, ONLY                :  AUTOSPC_SPCF, EPSIL, MEFMCORD, OTMSKIP, PRTNEU

      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SPARSE_MATRICES, ONLY       :  I_GMN  , J_GMN  , GMN    , I_GMNt  , J_GMNt , GMNt   , I_HMN, J_HMN, HMN,                 &
                                         I_KSF  , J_KSF  , KSF    , I_KSFD  , J_KSFD , KSFD   ,                                    &
                                         I_LMN  , J_LMN  , LMN    , I_MSF   , J_MSF  , MSF    ,                                    &
                                         SYM_GMN, SYM_HMN, SYM_KFS, SYM_KFSD, SYM_MFS, SYM_LMN

      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, MAXREQ, OGEL
      USE COL_VECS, ONLY              :  UF_COL, UG_COL, UN_COL, PHIXG_COL, PHIXN_COL, PM_COL, PS_COL,                             &
                                         QGm_COL, QGs_COL, QM_COL, QN_COL, QS_COL, QSYS_COL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE OUTPUT4_MATRICES, ONLY      :  OTM_MPCF, OTM_SPCF, TXT_MPCF, TXT_SPCF
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  MPCF_OUT, SPCF_OUT

      USE OFP2_USE_IFs

      IMPLICIT NONE

      LOGICAL                         :: WRITE_F06, WRITE_OP2, WRITE_PCH   ! flag
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP2'
      CHARACTER(LEN=*) , INTENT(IN)   :: WHAT              ! Indicator of whether to process output requests for SPC or MPC forces
      CHARACTER(LEN=*) , INTENT(IN)   :: ZERO_GEN_STIFF    ! Indicator of whether there are zero gen stiffs (can't calc MEFFMASS)
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRD_PRT_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: SPCF_ALL_SAME_CID ! Indicator of whether all grids, for the output set, have the same
!                                                            global coord sys for SPC force output
      CHARACTER( 1*BYTE)              :: SPCF_MEFM_MPF     ! Indicator of whether to calc MEFM, MPF via SPC force method
      CHARACTER( 1*BYTE)              :: MPCF_ALL_SAME_CID ! Indicator of whether all grids, for the output set, have the same
!                                                            global coord sys for MPC force output
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER( 1*BYTE)              :: PRINT_TOTALS      ! This will be set to 'Y' if DEBUG(92) > 0 so OLOAD, SPCF, MPCF force
!                                                            totals will be printed even if ALL_SAME_CID = 'N'
      CHARACTER(1*BYTE)               :: WRITE_OGEL(NGRID) ! 'Y'/'N' as to whether to write OGEL for a grid (used to avoid writing
!                                                            constr forces in subr WRITE_GRD_PRT_OUTPUTS for grids w no constr force

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITG               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: SC_OUT_REQ        ! If > 0, then req1uests for WHAT are to be output
      INTEGER(LONG), INTENT(INOUT)    :: OT4_GROW          ! Row number in OT4 file for grid related OTM descriptors
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            !
      LOGICAL, INTENT(INOUT)          :: NEW_RESULT        ! is this the first result of a table
      INTEGER(LONG)                   :: AGRID             ! An actual grid number
      INTEGER(LONG)                   :: G_SET_COL         ! Col number in TDOF where the G-set DOF's exist
      INTEGER(LONG)                   :: M_SET_COL         ! Col number in TDOF where the M-set DOF's exist
      INTEGER(LONG)                   :: N_SET_COL         ! Col number in TDOF where the N-set DOF's exist
      INTEGER(LONG)                   :: F_SET_COL         ! Col number in TDOF where the F-set DOF's exist
      INTEGER(LONG)                   :: S_SET_COL         ! Col number in TDOF where the S-set DOF's exist
      INTEGER(LONG)                   :: SA_SET_COL        ! Col number in TDOF where the SA-set DOF's exist
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IB1               ! If > 0, there are SPC or MPC force output requests
      INTEGER(LONG)                   :: IB2               ! If > 0, there are GP force balance requests so calc SPC or MPC forces
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IROW_FILE         ! Row number in text file
      INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: NREQ              ! Number of user requested outputs of displ/force
      INTEGER(LONG)                   :: NUM               ! Count of the number of rows added to array OGEL
      INTEGER(LONG)                   :: NUM_COMPS         ! Either 6 or 1 depending on whether grid is a physical grid or a SPOINT
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: SDOF              ! S-set DOF number
      INTEGER(LONG)                   :: SADOF             ! SA-set DOF number

      INTEGER(LONG)                   :: TDOF_ROW          ! Row no. in array TDOF to find GDOF DOF number

      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: MPF               ! Intermediate variable
      REAL(DOUBLE)                    :: QSK_COL(NDOFS)    ! KSF*UF
      REAL(DOUBLE)                    :: QSM_COL(NDOFS)    ! -EIGEN_VAL*MFS*UF
      REAL(DOUBLE)                    :: QMK_COL(NDOFM)    ! HMN*UF
      REAL(DOUBLE)                    :: QMM_COL(NDOFM)    ! -EIGEN_VAL*LMN*UF
      REAL(DOUBLE)                    :: QGs_MEFM(NDOFG)   ! QGs_COL transformed from global to coord system MEMFCORD
      REAL(DOUBLE)                    :: QGs_MEFM_SUM(6)   ! QGs_COL transformed from global to coord system MEMFCORD
      REAL(DOUBLE)                    :: QSA_MAX_ABS(6)    ! Max abs value of any QS force on an AUTOSPC'd DOF
      REAL(DOUBLE)                    :: QSA_SUM(6)        ! Sum of all QS forces on AUTOSPC'd DOF's
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND
      WRITE(ERR,9000) "OFP2 - SPC and MPC force"
 9000 FORMAT(' *DEBUG:    RUNNING=', A)
 9003 FORMAT(' *DEBUG:    ITABLE BAD=', i4)


      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
      DO I=1,MAXREQ
         DO J=1,MOGEL
            OGEL(I,J) = ZERO
         ENDDO
      ENDDO

      ! Initialize WRITE_OGEL
      DO I=1,NGRID
         WRITE_OGEL(I) = 'Y'
      ENDDO

! ---------------------------------------------------------------------------------------------------------------------------------
      ! Process SPC force requests
      NEW_RESULT = .TRUE.
      IF (WHAT == 'SPCF') THEN
      WRITE(ERR,9000) "OFP2 - SPC"
      WRITE(ERR,9003) ITABLE

         SPCF_MEFM_MPF = 'N'
         IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
            IF (SOL_NAME(1:12) /= 'GEN CB MODEL') THEN
               SPCF_MEFM_MPF = 'Y'
            ENDIF
         ENDIF

         IROW_FILE = 0
         IROW_MAT  = 0
         OT4_DESCRIPTOR = 'SPC force'

         CALL ALLOCATE_COL_VEC ('QS_COL', NDOFS, SUBR_NAME)! Initialize the array for SPC forces for this solution vector
         DO I=1,NDOFS
            QS_COL(I) = ZERO
         ENDDO

         CALL ALLOCATE_COL_VEC ('UF_COL', NDOFF, SUBR_NAME)! Get UF_COL to calc (KFS - EIGEN_VAL*MFS)*UF
         IF ((NTERM_KFS > 0) .OR. (NTERM_MFS > 0)) THEN
            CALL TDOF_COL_NUM ( 'F ', F_SET_COL )
            K = 0
            DO I=1,NDOFG
               IF (TDOFI(I,F_SET_COL) > 0) THEN
                  K = K + 1
                  UF_COL(K) = UG_COL(I)
               ENDIF
            ENDDO
         ENDIF

         QSK_COL = ZERO
         QSM_COL = ZERO

         IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            IF (NTERM_KFSD > 0) THEN                          ! Calc QSK = KSFD*UF
               CALL MATMULT_SFF ( 'KSFD', NDOFS, NDOFF, NTERM_KFSD, SYM_KFSD, I_KSFD, J_KSFD, KSFD, 'UF', NDOFF, 1, UF_COL, 'Y',   &
                                  'QSK', ONE, QSK_COL )
            ENDIF
         ELSE
            IF (NTERM_KFS  > 0) THEN                          ! Calc QSK = KSF*UF
               CALL MATMULT_SFF ( 'KSF ', NDOFS, NDOFF, NTERM_KFS , SYM_KFS , I_KSF , J_KSF , KSF , 'UF', NDOFF, 1, UF_COL, 'Y',   &
                                  'QSK', ONE, QSK_COL )
            ENDIF
         ENDIF
                                                           ! Don't do the following for CB soln (vecs are CB vecs, NOT eigen vecs,
!                                                            so no inertia effect is to be added for SPC forces due to CB vecs)
         IF (SOL_NAME(1:5) == 'MODES') THEN
            IF (NTERM_MFS > 0) THEN                        ! Calc QSM = -EIGEN_VAL*MSF*UF
               CALL MATMULT_SFF ( 'MSF', NDOFS, NDOFF, NTERM_MFS, SYM_MFS, I_MSF, J_MSF, MSF, 'UF', NDOFF, 1, UF_COL, 'Y',         &
                                  'QSM', -EIGEN_VAL(JVEC), QSM_COL )
            ENDIF
         ENDIF
         CALL DEALLOCATE_COL_VEC ( 'UF_COL' )

         DO I=1,NDOFS                                      ! Add (QSYS - PS) to QS to get  final QS
            IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               QS_COL(I) = QSK_COL(I)
            ELSE
               QS_COL(I) = QSK_COL(I) + QSM_COL(I) + QSYS_COL(I) - PS_COL(I)
            ENDIF
         ENDDO

         SPCF_ALL_SAME_CID = 'Y'                           ! Check if all grids, for which there will be output, have same coord sys
         DO I=1,NGRID-1                                    ! If not, then we won't write SPC output totals
            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_SPCF_BIT))
            IF (IB1 > 0) THEN
               IF (GRID(I+1,3) /= GRID(I,3)) THEN
                  SPCF_ALL_SAME_CID = 'N'
                  EXIT
               ENDIF
           ENDIF
         ENDDO

         NREQ = 0
         DO I=1,NGRID
            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_SPCF_BIT))
            IF (IB1 > 0) THEN
               NREQ = NREQ + 1
            ENDIF
         ENDDO

         NUM  = 0                                          ! Put QS into 2-d output array OGEL for this subcase (NDOFS x 6).
         DO I=1,NGRID

            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_SPCF_BIT))
            IB2 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_GPFO_BIT))
            IF ((IB1 > 0) .OR. (IB2 > 0) .OR. (SPCF_MEFM_MPF == 'Y')) THEN
               NUM = NUM + 1
               GID_OUT_ARRAY(NUM,1)       = GRID(I,1)
               GID_OUT_ARRAY(NUM,2)       = GRID(I,3)
               GID_OUT_ARRAY(NUM,MELGP+1) = GRID(I,5)

!xx            CALL CALC_TDOF_ROW_NUM ( GRID_ID(I), ROW_NUM_START, 'N' )
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(I), IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )

!xx            WRITE_OGEL(NUM) = 'N'                       ! Set WRITE_OGEL to 'Y' for all grids that have a component in S-set
!xx            DO J=1,NUM_COMPS
!xx               CALL TDOF_COL_NUM ( 'S ', S_SET_COL )
!xx               TDOF_ROW = ROW_NUM_START + J - 1
!xx               SDOF = TDOF(TDOF_ROW,S_SET_COL)
!xx               IF (SDOF > 0) THEN
!xx                  WRITE_OGEL(NUM) = 'Y'
!xx                  EXIT doj1
!xx               ENDIF
!xx            ENDDO doj1

               DO J=1,NUM_COMPS                            ! Calc SPC forces for all grids in requested output set (not only ones
                  CALL TDOF_COL_NUM ( 'S ', S_SET_COL )    ! that have a component in S-set)
                  CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
                  TDOF_ROW = ROW_NUM_START + J - 1
                  SDOF = TDOF(TDOF_ROW,S_SET_COL)
                  GDOF = TDOF(TDOF_ROW,G_SET_COL)
                  IF (SDOF > 0) THEN
                     OGEL(NUM,J)   = QS_COL(SDOF)
                     QGs_COL(GDOF) = QS_COL(SDOF)
                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (IB1 > 0)) THEN
                        IROW_FILE = IROW_FILE + 1
                        IROW_MAT  = IROW_MAT + 1
                        OTM_SPCF(IROW_MAT,JVEC) = QS_COL(SDOF)
                        IF (JVEC == 1) THEN
                           WRITE(TXT_SPCF(IROW_FILE), 9191) IROW_MAT, OT4_DESCRIPTOR, GRID(I,1), J
                        ENDIF
                     ENDIF
                  ELSE
                     OGEL(NUM,J) = ZERO
                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (IB1 > 0)) THEN
                        IROW_FILE = IROW_FILE + 1
                        IROW_MAT  = IROW_MAT + 1
                        OTM_SPCF(IROW_MAT,JVEC) = ZERO
                        IF (JVEC == 1) THEN
                           WRITE(TXT_SPCF(IROW_FILE), 9191) IROW_MAT, OT4_DESCRIPTOR, GRID(I,1), J
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               WRITE_OGEL(NUM) = 'Y'

               IF ((NUM == NREQ) .AND. (SC_OUT_REQ > 0)) THEN

                  WRITE_F06 = (SPCF_OUT(1:1) == 'Y')
                  WRITE_OP2 = (SPCF_OUT(2:2) == 'Y')
                  WRITE_PCH = (SPCF_OUT(3:3) == 'Y')
                  IF (WRITE_OP2) THEN
                     CALL WRITE_GRD_OP2_OUTPUTS ( JVEC, NUM, WHAT, ITABLE, NEW_RESULT )
                     NEW_RESULT = .FALSE.
                  ENDIF

                  IF (WRITE_PCH) THEN
                     CALL WRITE_GRD_PCH_OUTPUTS ( JVEC, NUM, WHAT )
                  ENDIF

                  IF (WRITE_F06) THEN
                     CALL CHK_OGEL_ZEROS ( NUM )
                     CALL WRITE_GRD_PRT_OUTPUTS ( JVEC, NUM, WHAT, IHDR, SPCF_ALL_SAME_CID, WRITE_OGEL )
                  ENDIF

                  EXIT
               ENDIF

               IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (IROW_FILE >= 1)) THEN
                  ! Write OTMSKIP blank separator lines
                  DO J=1,OTMSKIP
                     IROW_FILE = IROW_FILE + 1
                     WRITE(TXT_SPCF(IROW_FILE), 9199)
                  ENDDO
               ENDIF

            ENDIF

         ENDDO

! Calc modal effective masses and/or modal participation factors for current eigenector, if requested (but only if not CB soln)
! Requests for MPF, MEFFMASS for CB are handled elsewhere

         IF (SPCF_MEFM_MPF == 'Y') THEN                    ! NOTE: this would not be true if CB soln due to its value set above

            IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
               IF (ZERO_GEN_STIFF == 'N') THEN
                  CALL CONVERT_VEC_COORD_SYS ( 'Eigenvector', QGs_COL, QGs_MEFM, MEFMCORD )
                  DO J=1,6
                     QGs_MEFM_SUM(J) = ZERO
                  ENDDO
                  K = 0
                  DO I=1,NGRID
                     CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
                     DO J=1,NUM_COMPS
                        K = K + 1
                        QGs_MEFM_SUM(J) = QGs_MEFM_SUM(J) + QGs_MEFM(K)
                     ENDDO
                  ENDDO
                  DEN = EIGEN_VAL(JVEC)*GEN_MASS(JVEC)
                  DO J=1,6
                     MPF = QGs_MEFM_SUM(J)/DEN
                     IF (MPFACTOR_CALC == 'Y') THEN
                        MPFACTOR_N6(JVEC,J) = MPF
                     ENDIF
                     IF (MEFFMASS_CALC == 'Y') THEN
                        MEFFMASS(JVEC,J) = GEN_MASS(JVEC)*MPF*MPF
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

         ENDIF

! Calc largest and sum of QS forces for AUTOSPC DOF's. Write SPC forces on the SA DOF's if PARAM AUTOSPC_SPCF = 'Y'

         IF (NDOFSA > 0) THEN

            WRITE(F06,*)
            IF (AUTOSPC_SPCF == 'Y') THEN
               WRITE(F06,1101)
            ENDIF

            DO K=1,6
               QSA_MAX_ABS = ZERO
               QSA_SUM     = ZERO
            ENDDO
            DO I=1,NGRID
!xx            CALL CALC_TDOF_ROW_NUM ( GRID_ID(I), ROW_NUM_START, 'N' )
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(I), IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  CALL TDOF_COL_NUM ( 'S' , S_SET_COL )
                  CALL TDOF_COL_NUM ( 'SA', SA_SET_COL )
                  TDOF_ROW = ROW_NUM_START + J - 1
                  SDOF  = TDOF(TDOF_ROW,S_SET_COL)
                  SADOF = TDOF(TDOF_ROW,SA_SET_COL)
                  IF (SADOF > 0) THEN
                     AGRID = GRID_ID(I)
                     IF (AUTOSPC_SPCF == 'Y') THEN
                        WRITE(F06,1105) AGRID, J, QS_COL(SDOF)
                     ENDIF
                     IF (DABS(QS_COL(SDOF)) > QSA_MAX_ABS(J)) THEN
                        QSA_MAX_ABS(J) = DABS(QS_COL(SDOF))
                     ENDIF
                     QSA_SUM(J) = QSA_SUM(J) + QS_COL(SDOF)
                  ENDIF
               ENDDO
            ENDDO

            IF (AUTOSPC_SPCF == 'Y') THEN
               WRITE(F06,*)
               WRITE(F06,*)
            ENDIF

            IF (DEBUG(92) == 0) THEN
               PRINT_TOTALS = SPCF_ALL_SAME_CID
            ELSE
               PRINT_TOTALS = 'Y'
            ENDIF

            WRITE(F06,1102)
            IF (AUTOSPC_SPCF == 'N') THEN
               WRITE(F06,1103)
            ENDIF
            WRITE(F06,1104)

            WRITE(F06,9121) (QSA_MAX_ABS(J),J=1,6)

            IF (PRINT_TOTALS == 'Y') THEN
               WRITE(F06,9123) (QSA_SUM(J),J=1,6)
            ELSE
               WRITE(F06,9133)
            ENDIF

         ENDIF

         IF (WRITE_NEU .AND. (ANY_SPCF_OUTPUT > 0)) THEN
            CALL WRITE_FEMAP_GRID_VECS ( QGs_COL, FEMAP_SET_ID, 'SPCF' )
         ENDIF

         CALL DEALLOCATE_COL_VEC ( 'QS_COL' )

! ---------------------------------------------------------------------------------------------------------------------------------
!     Process MPC force requests

      ELSE IF (WHAT == 'MPCF') THEN
         WRITE(ERR,9000) "OFP2 - MPC"
         !IF (.NOT. NEW_RESULT) THEN
         !IF (NEW_RESULT .EQV. .FALSE.) THEN
         !IF (NEW_RESULT .EQ. .FALSE.) THEN    ! bad
         IF (ITABLE .NE. -1) THEN   ! not equal
           ! cleanup op2
           NEW_RESULT = .TRUE.
           ITABLE = -1
         ENDIF
         WRITE(ERR,9003) ITABLE

         ! Initialize the array for MPC forces for this solution vector
         IROW_FILE = 0
         IROW_MAT  = 0
         OT4_DESCRIPTOR = 'MPC force'

         CALL ALLOCATE_COL_VEC ( 'QM_COL', NDOFM, SUBR_NAME )
         DO I=1,NDOFM
            QM_COL(I) = ZERO
         ENDDO

         ! Partition NG-set vectors from G-set
         IF ((NTERM_HMN > 0) .OR. (NTERM_LMN > 0)) THEN
            CALL ALLOCATE_COL_VEC ( 'UN_COL', NDOFN, SUBR_NAME )
            CALL TDOF_COL_NUM ( 'N ', N_SET_COL )
            K = 0
            DO I=1,NDOFG
               IF (TDOFI(I,N_SET_COL) > 0) THEN
                  K = K + 1
                  UN_COL(K) = UG_COL(I)
               ENDIF
            ENDDO
         ENDIF

         DO I=1,NDOFM
            QMK_COL(I) = ZERO
            QMM_COL(I) = ZERO
         ENDDO

         ! Get 1st portion of QM: HMN*UN
         IF (NTERM_HMN > 0) THEN
            CALL MATMULT_SFF ( 'HMN', NDOFM, NDOFN, NTERM_HMN, SYM_HMN, I_HMN, J_HMN, HMN, 'UN', NDOFN, 1, UN_COL, 'Y',            &
                                  'QMK', ONE, QMK_COL )
         ENDIF

         ! Don't do the following for CB soln (vecs are CB vecs, NOT eigen vecs,
         ! so no inertia effect is to be added for MPC forces due to CB vecs)
         ! Get 1st portion of QM: LMN*UN
         IF (NTERM_LMN > 0) THEN
            IF (SOL_NAME(1:5) == 'MODES') THEN
               CALL MATMULT_SFF ( 'LMN', NDOFM, NDOFN, NTERM_LMN, SYM_LMN, I_LMN, J_LMN, LMN, 'UN', NDOFN, 1, UN_COL, 'Y',         &
                                  'QMM', -EIGEN_VAL(JVEC), QMM_COL )
            ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               CALL ALLOCATE_COL_VEC ( 'PHIXN_COL', NDOFN, SUBR_NAME )
               CALL TDOF_COL_NUM ( 'N ', N_SET_COL )
               K = 0
               DO I=1,NDOFG
                  IF (TDOFI(I,N_SET_COL) > 0) THEN
                     K = K + 1
                     PHIXN_COL(K) = PHIXG_COL(I)
                  ENDIF
               ENDDO
               CALL MATMULT_SFF ( 'LMN', NDOFM, NDOFN, NTERM_LMN, SYM_LMN, I_LMN, J_LMN, LMN, 'PHIXN', NDOFN, 1, PHIXN_COL, 'Y',   &
                                  'QMM', ONE, QMM_COL )
               CALL DEALLOCATE_COL_VEC ( 'PHIXN_COL' )
            ENDIF
         ENDIF

         CALL DEALLOCATE_COL_VEC ( 'UN_COL' )

         ! Add (HMN*UN - PM) to get QM
         DO I=1,NDOFM
            QM_COL(I) = QMK_COL(I) + QMM_COL(I) - PM_COL(I)
         ENDDO

         CALL DEALLOCATE_COL_VEC ( 'PM_COL' )

                                                           ! Calculate QN_COL = -GMNt*QM_COL
         CALL ALLOCATE_COL_VEC ( 'QN_COL', NDOFN, SUBR_NAME )
         CALL MATMULT_SFF ( 'GMNt', NDOFN, NDOFM, NTERM_GMN, SYM_GMN, I_GMNt, J_GMNt, GMNt, 'QM', NDOFM, 1, QM_COL, 'Y',           &
                            'QN', -ONE, QN_COL )

         CALL TDOF_COL_NUM ( 'G ', G_SET_COL )             ! Merge QN and QM to get QGm
         CALL TDOF_COL_NUM ( 'N ', N_SET_COL )
         CALL TDOF_COL_NUM ( 'M ', M_SET_COL )
         CALL MERGE_COL_VECS ( N_SET_COL, NDOFN, QN_COL, M_SET_COL, NDOFM, QM_COL, G_SET_COL, NDOFG, QGm_COL )
         CALL DEALLOCATE_COL_VEC ( 'QN_COL' )

         MPCF_ALL_SAME_CID = 'Y'                           ! Check if all grids, for which there will be output, have same coord sys
         DO I=1,NGRID-1                                    ! If not, then we won't write MPC output totals
            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_MPCF_BIT))
            IF (IB1 > 0) THEN
               IF (GRID(I+1,3) /= GRID(I,3)) THEN
                  MPCF_ALL_SAME_CID = 'N'
                  EXIT
               ENDIF
           ENDIF
         ENDDO

         NREQ = 0
         DO I=1,NGRID
            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_MPCF_BIT))
            IF (IB1 > 0) THEN
               NREQ = NREQ + 1
            ENDIF
         ENDDO

         ! Put QGm into 2-d output array OGEL for this subcase (NDOFS x 6).
         NUM = 0
         DO I=1,NGRID
            IB1 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_MPCF_BIT))
            IB2 = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_MPCF_BIT))
            IF ((IB1 > 0) .OR. (IB2 > 0)) THEN
               NUM = NUM + 1
               GID_OUT_ARRAY(NUM,1)       = GRID(I,1)
               GID_OUT_ARRAY(NUM,2)       = GRID(I,3)
               GID_OUT_ARRAY(NUM,MELGP+1) = GRID(I,5)
!xx            CALL CALC_TDOF_ROW_NUM ( GRID_ID(I), ROW_NUM_START, 'N' )
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(I), IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
                  TDOF_ROW = ROW_NUM_START + J - 1
                  GDOF = TDOF(TDOF_ROW,G_SET_COL)
!                 IF (GDOF > 0) THEN                       ! 08/01/05: Remove test on GDOF. It is always > 0
                     OGEL(NUM,J) = QGm_COL(GDOF)
                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                        IROW_FILE = IROW_FILE + 1
                        IROW_MAT  = IROW_MAT + 1
                        OTM_MPCF(IROW_MAT,JVEC) = QGm_COL(GDOF)
                        IF (JVEC == 1) THEN
                           WRITE(TXT_MPCF(IROW_FILE), 9191) IROW_MAT, OT4_DESCRIPTOR, GRID(I,1), J
                        ENDIF
                     ENDIF
!                 ELSE
!                    OGEL(NUM,J) = ZERO
!                 ENDIF
               ENDDO
               WRITE_OGEL(NUM) = 'N'
               DO J=1,NUM_COMPS
                  IF (OGEL(NUM,J) /= ZERO) THEN
                     WRITE_OGEL(NUM) = 'Y'
                  ENDIF
               ENDDO

               IF ((NUM == NREQ) .AND. (SC_OUT_REQ > 0)) THEN

                  WRITE_F06 = (MPCF_OUT(1:1) == 'Y')
                  WRITE_OP2 = (MPCF_OUT(2:2) == 'Y')
                  WRITE_PCH = (MPCF_OUT(3:3) == 'Y')
                  IF (WRITE_OP2) THEN
                     CALL WRITE_GRD_OP2_OUTPUTS ( JVEC, NUM, WHAT, ITABLE, NEW_RESULT )
                     NEW_RESULT = .FALSE.
                  ENDIF

                  IF (WRITE_PCH) THEN
                     CALL WRITE_GRD_PCH_OUTPUTS ( JVEC, NUM, WHAT )
                  ENDIF

                  IF (WRITE_F06) THEN
                     CALL CHK_OGEL_ZEROS ( NUM )
                     CALL WRITE_GRD_PRT_OUTPUTS ( JVEC, NUM, WHAT, IHDR, MPCF_ALL_SAME_CID, WRITE_OGEL )
                  ENDIF

                  EXIT

               ENDIF

               IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (IROW_MAT >= 1)) THEN
                  DO J=1,OTMSKIP                           ! Write OTMSKIP blank separator lines
                     IROW_FILE = IROW_FILE + 1
                     WRITE(TXT_MPCF(IROW_FILE), 9199)
                  ENDDO
               ENDIF

            ENDIF

         ENDDO

         IF (WRITE_NEU .AND. (ANY_MPCF_OUTPUT > 0)) THEN
            CALL WRITE_FEMAP_GRID_VECS ( QGm_COL, FEMAP_SET_ID, 'MPCF' )
         ENDIF

         CALL DEALLOCATE_COL_VEC ( 'QM_COL' )

! ---------------------------------------------------------------------------------------------------------------------------------
! Coding error: illegal input WHAT

      ELSE

         WRITE(ERR,9100) WHAT
         WRITE(ERR,9100) WHAT
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(55X,'SPC Forces on AUTOSPC''d DOF''s',/,48X,'(in global coordinate system at each grid)'                              &
          ,/,58X,'GRID   COMP      VALUE',//)

 1102 FORMAT(52X,'SPC Force Summary on AUTOSPC''d DOF''s',/,49X,'(in global coordinate system at each grid)')

 1103 FORMAT(8X,'For detailed information on the SPC forces for the AUTOSPC DOF''s, set field 7 (AUTOSPC_SPCF) on the ',           &
                'PARAM AUTOSPC entry to Y',/)

 1104 FORMAT(30X,'T1             T2             T3             R1             R2             R3',/)

 1105 FORMAT(54X,I8,5X,I2,1ES15.6)

 9100 FORMAT(' *ERROR  9100: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ILLEGAL INPUT FOR VARIABLE "WHAT" = ',A)

 9121 FORMAT(1X,'ABS AUTOSPC FORCES  :',6(ES15.6))

 9123 FORMAT(1X,'AUTOSPC FORCE TOTALS:',6(ES15.6),/,5X,'(for output set)')

 9133 FORMAT(10X,'             -------------- -------------- -------------- -------------- -------------- --------------',/,       &
             1X,'AUTOSPC FORCE TOTALS: not printed since all grids do not have the same global coordinate system')

 9191 FORMAT(I8,1X,A,I8,I8)

 9199 FORMAT(' ')


! **********************************************************************************************************************************

      END SUBROUTINE OFP2
