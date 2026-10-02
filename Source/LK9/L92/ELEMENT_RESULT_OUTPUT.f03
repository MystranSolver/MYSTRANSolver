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

   MODULE ELEMENT_RESULT_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OFP3_ELFE_1D, OFP3_ELFE_2D, OFP3_ELFN, OFP3_STRE_NO_PCOMP, OFP3_STRN_NO_PCOMP

   CONTAINS

      SUBROUTINE OFP3_ELFE_1D ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element engr force output requests for 1D (ELAS, BUSH, ROD, BAR) elements for one subcase. Results go into array OGEL
! for later output in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_ELFE_BIT, FATAL_ERR, IBIT, INT_SC_NUM, MBUG, MOGEL,&
                                         NELE, NCBAR, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCROD, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  ANY_ELFE_OUTPUT, BUSH_CID, BUSH_VVEC, EDAT, ELAS_COMP, ELEM_LEN_12, ELEM_LEN_AB, EPNT,    &
                                         ETYPE, EID, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS, OFFDIS_GA_GB, OFFDIS_L,            &
                                         PE_GA_GB, PEL, PLY_NUM, STRESS, TE, TE_GA_GB, TYPE, XEL
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFE, TXT_ELFE

      USE EMG_MOD, ONLY               :  EMG
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: CALC_ELEM_NODE_FORCES
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE ELEMENT_STRESS_RECOVERY, ONLY: CALC_ELEM_STRESSES
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE ELEMENT_OUTPUT_WRITERS, ONLY:  WRITE_ELEM_ENGR_FORCE
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_ELFO_VECS

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_FORCE_OUTPUT, ONLY      :  SET_OEF_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_ELFE_1D'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER(20*BYTE)              :: FORCE_ITEM(8)     ! Char description of element engineering forces
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_ELFE        ! If > 0, there are ELFORCE(ENGR) requests for some elems
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: IERROR       = 0  ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NDUM              ! An arg passed to CALC_ELEM_STRESSES
      INTEGER(LONG)                   :: NUM_ELEM          ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_ELEM since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file


      REAL(DOUBLE)                    :: DUM0(6,12)        ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM1(6)           ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM21(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM22(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM31(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM32(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: EEF(6)            ! Element engineering force for BUSH
      REAL(DOUBLE)                    :: DX,DY,DZ          ! Offset dist1
      REAL(DOUBLE)                    :: FORCES(12)        ! Forces at the grid points
      REAL(DOUBLE)                    :: LENGTH

      ! OP2 parameters
      INTEGER(LONG)                   :: ITABLE            ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME        ! the op2 table name

      REAL(DOUBLE)                    :: TET(3,3)          ! Transpose of TE
      REAL(DOUBLE)                    :: TET_GA_GB(3,3)    ! Transpose of TE_GA_GB
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND

! **********************************************************************************************************************************
!     Initialize
      TABLE_NAME = "OEF ERR "
      ITABLE = 0


      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element engineering force requests for BAR, BUSH, ELAS, ROD. Use subr CALC_ELEM_NODE_FORCES and then convert the node
! forces to engineering forces (see equations below after subr CALC_ELEM_NODE_FORCES is called)

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'Y'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'Y'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      FORCE_ITEM(1) = 'M1a: Mom Plane1 EndA'
      FORCE_ITEM(2) = 'M1b: Mom Plane2 EndA'
      FORCE_ITEM(3) = 'M2a: Mom Plane1 EndB'
      FORCE_ITEM(4) = 'M2b: Mom Plane2 EndB'
      FORCE_ITEM(5) = 'V1 : Shear Plane1   '
      FORCE_ITEM(6) = 'V2 : Shear Plane2   '
      FORCE_ITEM(7) = 'FX : Axial force    '
      FORCE_ITEM(8) = 'T  : Torque         '

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the no. requests/elem.
         NELREQ(I) = 0
      ENDDO

      DO I=1,METYPE
         DO J=1,NELE
            IF ((ETYPE(J)(1:3) == 'BAR') .OR. (ETYPE(J)(1:4) == 'BUSH') .OR. (ETYPE(J)(1:4) == 'ELAS') .OR.                        &
                (ETYPE(J)(1:3) == 'ROD'))THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                  IF (ELOUT_ELFE > 0) THEN
                     NELREQ(I) = NELREQ(I) + 1
                   ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element engineering force, ELFO'
reqs2:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs2
         NUM_ELEM  = 0
         NUM_OGEL = 0

elems_2: DO J = 1,NELE
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF ((ETYPE(J)(1:3) == 'BAR') .OR. (ETYPE(J)(1:4) == 'BUSH') .OR. (ETYPE(J)(1:4) == 'ELAS') .OR.                        &
                (ETYPE(J)(1:3) == 'ROD'))THEN

               IF (ETYPE(J) == ELMTYP(I)) THEN
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                  IF (ELOUT_ELFE > 0) THEN

                     PLY_NUM = 0                           ! 'N' in call to EMG means do not write to BUG file
                     CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                     IF (NUM_EMG_FATAL_ERRS > 0) THEN
                        IERROR = IERROR + 1
                        CYCLE elems_2
                     ENDIF
                     CALL ELMDIS

                     CALL CALC_ELEM_NODE_FORCES            ! Use NODE to get engr forces (SE matrices don't have torque)

                     NUM_OGEL = NUM_OGEL + 1
                     IF (NUM_OGEL > MAXREQ) THEN
                        WRITE(ERR,9200) SUBR_NAME, MAXREQ
                        WRITE(F06,9200) SUBR_NAME, MAXREQ
                        FATAL_ERR = FATAL_ERR + 1
                        CALL OUTA_HERE ( 'Y' )             ! Coding error (dim of array OGEL too small), so quit
                     ENDIF

!                    ---------------------------------------------------------------------------------------------------------------
                     IF (ETYPE(J)(1:4) == 'ELAS') THEN     ! Set engr forces based on the node force values
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        NDUM = 0
                        CALL CALC_ELEM_STRESSES ( 1, NDUM, 0, 'N', 'N' )
                        OGEL(NUM_OGEL,1) = STRESS(1)       ! ELAS engr force is stored in the stress array
                     ! end elas
!                    ---------------------------------------------------------------------------------------------------------------
                     ELSE IF (ETYPE(J)(1:4) == 'BUSH') THEN
                        DO K=1,3                           ! Calculate element forces in GA-GB axes (x along line from GA to GB)
                          DUM21(K) = ZERO
                          DUM22(K) = ZERO
                          DUM31(K) = ZERO
                          DUM32(K) = ZERO
                        ENDDO
                        IF (ELEM_LEN_12 > .0001D0) THEN    ! Element has a GA-GB axis so start with PE_GA_GB
                           DX = ABS(OFFDIS_GA_GB(2,1))
                           DY =    (OFFDIS_GA_GB(2,2))
                           DZ =    (OFFDIS_GA_GB(2,3))
                           DUM21(1) =  PE_GA_GB(7)
                           DUM21(2) =  PE_GA_GB(8)
                           DUM21(3) =  PE_GA_GB(9)
                           DUM31(1) =  PE_GA_GB(8)*DZ - PE_GA_GB(9)*DY + PE_GA_GB(10)
                           DUM31(2) = -PE_GA_GB(7)*DZ - PE_GA_GB(9)*DX + PE_GA_GB(11)
                           DUM31(3) =  PE_GA_GB(7)*DY + PE_GA_GB(8)*DX + PE_GA_GB(12)
                           DO K=1,3
                              EEF(K)   = DUM21(K)
                              EEF(K+3) = DUM31(K)
                           ENDDO
                                                           ! There is a local elem coord system (via CID or v-vec)
                        IF ((BUSH_CID >= 0) .OR. (BUSH_VVEC /= 0)) THEN

                              DO K=1,3                     ! Transform elem forces from GA-GB axes to basic
                                 DO L=1,3
                                    TET_GA_GB(K,L) = TE_GA_GB(L,K)
                                 ENDDO
                              ENDDO

                              DO K=1,3
                              DUM22(K) = ZERO
                              DUM32(K) = ZERO
                              ENDDO

                              CALL MATMULT_FFF ( TET_GA_GB, DUM21, 3, 3, 1, DUM22 )
                              CALL MATMULT_FFF ( TET_GA_GB, DUM31, 3, 3, 1, DUM32 )
                              DO K=1,3
                                 EEF(K)   = DUM22(K)
                                 EEF(K+3) = DUM32(K)
                              ENDDO
                           ENDIF
                                                           ! Transform elem forces from basic to local
                           IF ((BUSH_CID > 0) .OR. (BUSH_VVEC /= 0)) THEN

                              DO K=1,3
                                 DO L=1,3
                                    TET(K,L) = TE(L,K)
                                 ENDDO
                              ENDDO

                              CALL MATMULT_FFF ( TE, DUM22, 3, 3, 1, DUM21 )
                              CALL MATMULT_FFF ( TE, DUM32, 3, 3, 1, DUM31 )
                              DO K=1,3
                                 EEF(K)   = DUM21(K)
                                 EEF(K+3) = DUM31(K)
                              ENDDO

                           ENDIF

                        ELSE                               ! Element has GA, GB coincident so element loads are in PEL

                           DO K=1,6
                              EEF(K) = PEL(K+6)
                           ENDDO

                        ENDIF

                     OGEL(NUM_OGEL,1) = EEF(1)             ! Now set OGEL output equal to the correct EEF's
                     OGEL(NUM_OGEL,2) = EEF(2)
                     OGEL(NUM_OGEL,3) = EEF(3)
                     OGEL(NUM_OGEL,4) = EEF(4)
                     OGEL(NUM_OGEL,5) = EEF(5)
                     OGEL(NUM_OGEL,6) = EEF(6)

                     ! end bush
!                    ---------------------------------------------------------------------------------------------------------------
                     ELSE IF (ETYPE(J)(1:3) == 'ROD') THEN
                        OGEL(NUM_OGEL,7) = -PEL(1)         ! Fx  (axial force for ROD)
                        OGEL(NUM_OGEL,8) = -PEL(4)         ! T   (torque for ROD)
                     !end rod
!                    ---------------------------------------------------------------------------------------------------------------
                     ELSE IF (ETYPE(J)(1:3) == 'BAR') THEN
                        LENGTH = ELEM_LEN_AB
                        OGEL(NUM_OGEL,1) = -PEL(6)                 ! M1a (bending moment, plane 1, end a for BAR)
                        OGEL(NUM_OGEL,2) =  PEL(5)                 ! M2a (bending moment, plane 2, end a for BAR)
                        OGEL(NUM_OGEL,3) = -PEL(6) + PEL(2)*LENGTH ! M1b (bending moment, plane 1, end b for BAR)
                        OGEL(NUM_OGEL,4) =  PEL(5) + PEL(3)*LENGTH ! M2b (bending moment, plane 2, end b for BAR)
                        OGEL(NUM_OGEL,5) = -PEL(2)                 ! V1  (plane 1 shear for BAR)
                        OGEL(NUM_OGEL,6) = -PEL(3)                 ! V2  (plane 2 shear for BAR)
                        OGEL(NUM_OGEL,7) = -PEL(1)                 ! Fx  (axial force for BAR)
                        OGEL(NUM_OGEL,8) = -PEL(4)                 ! T   (torque for BAR)
                     ENDIF !end bar

                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                        IF (ETYPE(J)(1:4) == 'ELAS') THEN
                           DO K=1,1                        ! ELAS has only 1 engr force
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF (ETYPE(J)(1:4) == 'BUSH') THEN
                           DO K=1,6
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF (ETYPE(J)(1:3) == 'ROD') THEN
                           DO K=7,8
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF (ETYPE(J)(1:3) == 'BAR') THEN
                           DO K=1,8
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF (ETYPE(J)(1:4) == 'BEAM') THEN
                           FATAL_ERR = FATAL_ERR + 1
                           WRITE(ERR,963) SUBR_NAME, ETYPE(J)
                           WRITE(ERR,963) SUBR_NAME, ETYPE(J)

                        ENDIF

                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                     ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_ELFE(OT4_EROW), 9199)
                        ENDDO
                     ENDIF

                     NUM_ELEM = NUM_ELEM + 1
                     EID_OUT_ARRAY(NUM_ELEM,1) = EID
                     IF (NUM_ELEM == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100                    FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                        WRITE(ERR,100) "F1",ETYPE(J),TABLE_NAME,ITABLE
                        CALL SET_OEF_TABLE_NAME(ETYPE(J), TABLE_NAME, ITABLE)
                        WRITE(ERR,100) "F2",ETYPE(J),TABLE_NAME,ITABLE
                        CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_ELEM, IHDR, 1, ITABLE )
                        EXIT
                     ENDIF

                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_2

      ENDDO reqs2
 10   FORMAT("*DEBUG:      OEF_END 1D:    TABLE_NAME",A)
      WRITE(ERR,10) TABLE_NAME
      IF ((TABLE_NAME .NE. "OEF ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF

      IF (WRITE_NEU .AND. (ANY_ELFE_OUTPUT > 0)) THEN

! bar    ---------------------------------------------------------------------------------------------------------------------------
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBAR, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out BAR engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:3) == 'BAR') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               LENGTH = ELEM_LEN_AB
               CALL ELMDIS
               CALL CALC_ELEM_NODE_FORCES
               FEMAP_EL_VECS(NUM_FROWS,1) = -PEL(6)                 ! M1a (bending moment, plane 1, end a for BAR)
               FEMAP_EL_VECS(NUM_FROWS,2) = -PEL(6) + PEL(2)*LENGTH ! M1b (bending moment, plane 1, end b for BAR)
               FEMAP_EL_VECS(NUM_FROWS,3) =  PEL(5)                 ! M2a (bending moment, plane 2, end a for BAR)
               FEMAP_EL_VECS(NUM_FROWS,4) =  PEL(5) + PEL(3)*LENGTH ! M2b (bending moment, plane 2, end b for BAR)
               FEMAP_EL_VECS(NUM_FROWS,5) = -PEL(2)                 ! V1  (plane 1 shear for BAR)
               FEMAP_EL_VECS(NUM_FROWS,6) = -PEL(3)                 ! V2  (plane 2 shear for BAR)
               FEMAP_EL_VECS(NUM_FROWS,7) = -PEL(1)                 ! Fx  (axial force for BAR or ROD)
               FEMAP_EL_VECS(NUM_FROWS,8) = -PEL(4)                 ! T   (torque for BAR or ROD)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'BAR     ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! bush   ---------------------------------------------------------------------------------------------------------------------------
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBUSH, 6, SUBR_NAME )
         DO J=1,NELE                                       ! Write out BUSH engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:4) == 'BUSH') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL CALC_ELEM_NODE_FORCES

               IF (ETYPE(J)(1:4) == 'BUSH') THEN

                  DO K=1,3                           ! Calculate element forces in GA-GB axes (x along line from GA to GB)
                    DUM21(K) = ZERO
                    DUM22(K) = ZERO
                    DUM31(K) = ZERO
                    DUM32(K) = ZERO
                  ENDDO

                  IF (ELEM_LEN_12 > .0001D0) THEN    ! Element has a GA-GB axis so start with PE_GA_GB

                     DX = ABS(OFFDIS_GA_GB(2,1))
                     DY =    (OFFDIS_GA_GB(2,2))
                     DZ =    (OFFDIS_GA_GB(2,3))

                     DUM21(1) =  PE_GA_GB(7)
                     DUM21(2) =  PE_GA_GB(8)
                     DUM21(3) =  PE_GA_GB(9)
                     DUM31(1) =  PE_GA_GB(8)*DZ - PE_GA_GB(9)*DY + PE_GA_GB(10)
                     DUM31(2) = -PE_GA_GB(7)*DZ - PE_GA_GB(9)*DX + PE_GA_GB(11)
                     DUM31(3) =  PE_GA_GB(7)*DY + PE_GA_GB(8)*DX + PE_GA_GB(12)

                     DO K=1,3
                        EEF(K)   = DUM21(K)
                        EEF(K+3) = DUM31(K)
                     ENDDO
                                                     ! There is a local elem coord system (via CID or v-vec)
                     IF ((BUSH_CID >= 0) .OR. (BUSH_VVEC /= 0)) THEN

                        DO K=1,3                     ! Transform elem forces from GA-GB axes to basic
                           DO L=1,3
                              TET_GA_GB(K,L) = TE_GA_GB(L,K)
                           ENDDO
                        ENDDO

                        DO K=1,3
                        DUM22(K) = ZERO
                        DUM32(K) = ZERO
                        ENDDO

                        CALL MATMULT_FFF ( TET_GA_GB, DUM21, 3, 3, 1, DUM22 )
                        CALL MATMULT_FFF ( TET_GA_GB, DUM31, 3, 3, 1, DUM32 )
                        DO K=1,3
                           EEF(K)   = DUM22(K)
                           EEF(K+3) = DUM32(K)
                        ENDDO

                     ENDIF
                                                     ! Transform elem forces from basic to local
                     IF ((BUSH_CID > 0) .OR. (BUSH_VVEC /= 0)) THEN

                        DO K=1,3
                           DO L=1,3
                              TET(K,L) = TE(L,K)
                           ENDDO
                        ENDDO

                        CALL MATMULT_FFF ( TE, DUM22, 3, 3, 1, DUM21 )
                        CALL MATMULT_FFF ( TE, DUM32, 3, 3, 1, DUM31 )
                        DO K=1,3
                           EEF(K)   = DUM21(K)
                           EEF(K+3) = DUM31(K)
                        ENDDO

                     ENDIF

                  ELSE                               ! Element has GA, GB coincident so element loads are in PEL

                     DO K=1,6
                        EEF(K) = PEL(K+6)
                     ENDDO

                  ENDIF

               ENDIF

               FEMAP_EL_VECS(NUM_FROWS,1) = EEF(1)
               FEMAP_EL_VECS(NUM_FROWS,2) = EEF(2)
               FEMAP_EL_VECS(NUM_FROWS,3) = EEF(3)
               FEMAP_EL_VECS(NUM_FROWS,4) = EEF(4)
               FEMAP_EL_VECS(NUM_FROWS,5) = EEF(5)
               FEMAP_EL_VECS(NUM_FROWS,6) = EEF(6)

            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'BUSH    ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! rod    ---------------------------------------------------------------------------------------------------------------------------
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCROD, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out ROD engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:3) == 'ROD') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL CALC_ELEM_NODE_FORCES
               FEMAP_EL_VECS(NUM_FROWS,7) = -PEL(1)                 ! Fx  (axial force for BAR or ROD)
               FEMAP_EL_VECS(NUM_FROWS,8) = -PEL(4)                 ! T   (torque for BAR or ROD)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'ROD     ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! elas1  ---------------------------------------------------------------------------------------------------------------------------
         NDUM = 0
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS1, 2, SUBR_NAME )
         DO J=1,NELE                                       ! Write out ELAS1 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:5) == 'ELAS1') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               CALL CALC_ELEM_STRESSES ( NCELAS1, NDUM, NUM_FROWS, 'N', 'Y' )
               FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'ELAS1   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! elas2  ---------------------------------------------------------------------------------------------------------------------------
         NDUM = 0
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS2, 2, SUBR_NAME )
         DO J=1,NELE                                       ! Write out ELAS2 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:5) == 'ELAS2') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               CALL CALC_ELEM_STRESSES ( NCELAS2, NDUM, NUM_FROWS, 'N', 'Y' )
               FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'ELAS2   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! elas3  ---------------------------------------------------------------------------------------------------------------------------
         NDUM = 0
         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS3, 2, SUBR_NAME )
         DO J=1,NELE                                       ! Write out ELAS3 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:5) == 'ELAS3') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               CALL CALC_ELEM_STRESSES ( NCELAS3, NDUM, NUM_FROWS, 'N', 'Y' )
               FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'ELAS3   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

! elas4  ---------------------------------------------------------------------------------------------------------------------------
         NDUM = 0
         NUM_FROWS= 0                                      ! 'N' in call to EMG means do not write to BUG file
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS4, 2, SUBR_NAME )
         DO J=1,NELE                                       ! Write out ELAS4 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:5) == 'ELAS4') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               CALL CALC_ELEM_STRESSES ( NCELAS4, NDUM, NUM_FROWS, 'N', 'Y' )
               FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'ELAS4   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT ENGINEERING FORCE'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
  963 FORMAT(' *ERROR   946: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NO CODE FOR 1D ELEMENT TYPE "',A,'"')

 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9199 FORMAT(' ')

 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE OFP3_ELFE_1D


      SUBROUTINE OFP3_ELFE_2D ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element engr force output requests for 2D (TRIA3, QUAD4, SHEAR) elements for one subcase. Results go into array OGEL
! for later output in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_ELFE_BIT, FATAL_ERR, IBIT, INT_SC_NUM, MBUG, MOGEL,                   &
                                         WARN_ERR, NELE, NCQUAD4, NCQUAD4K, NCSHEAR, NCTRIA3, NCTRIA3K, SOL_NAME, MAX_STRESS_POINTS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      use model_stuf, only            :  pcomp_props
      USE MODEL_STUF, ONLY            :  ANY_ELFE_OUTPUT, EDAT, EPNT, ETYPE, FCONV, EID, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,&
                                         PLY_NUM, TYPE, STRESS, SHELL_STR_ANGLE, NUM_SEi, ELGP, AGRID
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  FORC_LOC
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFE, TXT_ELFE

      USE VECTOR_GEOMETRY, ONLY       :  PLANE_COORD_TRANS_21
      USE RESULT_COORDINATES, ONLY    :  TRANSFORM_SHELL_STR
      USE COMPOSITE_SHELL_PREPARATION, ONLY:  IS_ELEM_PCOMP_PROPS
      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: SHELL_ENGR_FORCE_OGEL
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: POLYNOM_FIT_STRE_STRN
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE ELEMENT_OUTPUT_WRITERS, ONLY:  WRITE_ELEM_ENGR_FORCE
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_ELFO_VECS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_FORCE_OUTPUT, ONLY      :  SET_OEF_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_ELFE_2D'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER(20*BYTE)              :: FORCE_ITEM(8)     ! Char description of element engineering forces
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
     !INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the op2 subtable number, should be -3, -5, ...
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_ELFE        ! If > 0, there are ELFORCE(ENGR) requests for some elems
      INTEGER(LONG)                   :: I,J,K,M           ! DO loop indices
      INTEGER(LONG)                   :: IERROR       = 0  ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elem points processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff force points for one element
      integer(long)                   :: num_pcomp_elems   ! number of elements that are composites (used to prevent output of engr
!                                                            forces for PCOMP elems until I fix that output)

                                                           ! Stress index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRESS_OUT_ERR_INDEX(MAX_STRESS_POINTS)

                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRESS_OUT_PCT_ERR(MAX_STRESS_POINTS)

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRESS_OUT_PCT_ERR

      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.
                                                           ! Array of values from array STRESS for all stress points
      REAL(DOUBLE)                    :: STRESS_RAW(9,MAX_STRESS_POINTS)
                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRESS_OUT(9,MAX_STRESS_POINTS)


      ! OP2 parameters
      INTEGER(LONG)                   :: ITABLE            ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME        ! the op2 table name

      LOGICAL                         :: WRITE_NEU
      INTRINSIC IAND

! **********************************************************************************************************************************
!     Initialize
      TABLE_NAME = "OEF ERR "
      ITABLE = 0

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element engineering force requests for plate and USERIN elements.
! For the MIN3,4 elements, the stiffness matrix has to be generated to get FCONV(3). Therefore, set OPT(4) to 'N' initially, and if
! the element being processed is a MIN3 or MIN4, reset OPT(4) to 'Y' in the calculation loop prior to EMG call.

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      FORCE_ITEM(1) = 'Nxx: Normal x Force '
      FORCE_ITEM(2) = 'Nyy: Normal y Force '
      FORCE_ITEM(3) = 'Nxy: Shear xy Force '
      FORCE_ITEM(4) = 'Mxx: Moment x Plane '
      FORCE_ITEM(5) = 'Myy: Moment y Plane '
      FORCE_ITEM(6) = 'Mxy: Twist Mom xy   '
      FORCE_ITEM(7) = 'Qx : Transv Shear x '
      FORCE_ITEM(8) = 'Qy : Transv Shear y '

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the number of requests per element
         NELREQ(I) = 0
      ENDDO

      num_pcomp_elems = 0                                  ! Remove lower case code when I fix engr force output for PCOMP's
      DO I=1,METYPE
         DO J=1,NELE
            IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'QUAD8') .OR.                     &
               (ETYPE(J)(1:5) == 'SHEAR') .OR. (ETYPE(J)(1:6) == 'USERIN')) THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  call is_elem_pcomp_props ( j )
                  if (pcomp_props == 'N') then
                     IF ((FORC_LOC == 'CORNER  ') .OR.                                                                             &
                         (ETYPE(J)(1:5) == 'QUAD8')) THEN
                        NUM_PTS(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS(I) = 1
                     ENDIF
                     ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NELREQ(I) = NELREQ(I) + NUM_PTS(I)
                     ENDIF
                  else
                     num_pcomp_elems = num_pcomp_elems + 1
                  endif
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!xx   DO I=1,METYPE                                        ! Engr force requests for ELAS elems not honored. Give message
!xx      IF ((ELMTYP(I)(1:4) == 'ELAS') .AND. (NELREQ(I) > 0)) THEN
!xx         WARN_ERR = WARN_ERR + 1
!xx         WRITE(ERR,9204) NELREQ(I), ELMTYP(I)
!xx         IF (SUPWARN == 'N') THEN
!xx            WRITE(F06,9204) NELREQ(I), ELMTYP(I)
!xx         ENDIF
!xx      ENDIF
!xx   ENDDO

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element engineering force'
reqs3:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs3
         NUM_OGEL_ROWS = 0
         NUM_OGEL = 0

elems_3: DO J = 1,NELE
            call is_elem_pcomp_props ( j )                 ! Remove lower case code when I fix engr force output for PCOMP's
            if (pcomp_props == 'N') then
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'QUAD8') .OR.                  &
                  (ETYPE(J)(1:5) == 'SHEAR') .OR. (ETYPE(J)(1:6) == 'USERIN')) THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'QUAD8') .OR.         &
                           (ETYPE(J)(1:5) == 'SHEAR')) THEN
                           OPT(4) = 'Y'
                        ENDIF
                        DO K=0,MBUG-1
                           WRT_BUG(K) = 0
                        ENDDO
                        PLY_NUM = 0                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE elems_3
                        ENDIF
                        OPT(4) = 'N'
                        CALL ELMDIS


                        DO M=1,NUM_PTS(I)                  ! Gauss point stress
                           CALL ELEM_STRE_STRN_ARRAYS ( M )
                           STRESS_RAW(:,M) = STRESS(:)
                        ENDDO

                        STRESS_OUT(:,1) = STRESS_RAW(:,1)  ! Set STRAIN_OUT for NUM_PTS(I) = 1

                        IF ((FORC_LOC == 'CORNER  ') .OR.                                                                          &
                            (ETYPE(J)(1:5) == 'QUAD8')) THEN

                           IF (TYPE(1:5) == 'QUAD4') THEN

                                                           ! Extrapolate stress to corners
                              CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,              &
                                    STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                           ELSEIF (ETYPE(J)(1:5) == 'QUAD8') THEN

                                                           ! Extrapolate stress to corners
                              CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,              &
                                STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transfrom stress to element coordinate system
                              DO M=2,NUM_PTS(I)
                                 CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                                 CALL TRANSFORM_SHELL_STR( TEL, STRESS_OUT(:,M), ONE)
                              ENDDO

                                                           ! Center stress is the average of corner stress in element coordinates.
                                                           ! In MSC, the center force and moment resultants are the average but
                                                           ! this is equivalent.
                              STRESS_OUT(:,1) = (STRESS_OUT(:,2) + STRESS_OUT(:,3) + STRESS_OUT(:,4) + STRESS_OUT(:,5)) / FOUR

                           ENDIF

                        ENDIF

                        DO M=1,NUM_PTS(I)                  ! Calculate forces and moments from stresses
                           STRESS(:) = STRESS_OUT(:,M)
                           CALL SHELL_ENGR_FORCE_OGEL ( NUM_OGEL )

                           NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                           EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
                           GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                           DO K=1,ELGP
                              GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                           ENDDO

                        ENDDO


                        IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                           DO K=1,8
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                           DO K=1,OTMSKIP                     ! Write OTMSKIP blank separator lines
                              OT4_EROW = OT4_EROW + 1
                              WRITE(TXT_ELFE(OT4_EROW), 9199)
                           ENDDO
                        ENDIF

                        IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                           CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100                       FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                           WRITE(ERR,100) "F3",ETYPE(J),TABLE_NAME,ITABLE
                           CALL SET_OEF_TABLE_NAME(ETYPE(J), TABLE_NAME, ITABLE)
                           WRITE(ERR,100) "F4",ETYPE(J),TABLE_NAME,ITABLE
                           CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_PTS(I), ITABLE )
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO elems_3
      ENDDO reqs3

 10   FORMAT("*DEBUG:      OEF_END 2D:    TABLE_NAME",A)
      WRITE(ERR,10) TABLE_NAME
      IF ((TABLE_NAME .NE. "OEF ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF

      IF (WRITE_NEU .AND. (ANY_ELFE_OUTPUT > 0)) THEN

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out TRIA3K engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'TRIA3K  ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out TRIA3 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'TRIA3   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out QUAD4K engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'QUAD4K  ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out QUAD4 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'QUAD4 ') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'QUAD4   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out SHEAR engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'SHEAR') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'SHEAR   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT ENGINEERING FORCE'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

 9204 FORMAT(' *WARNING    : ELEMENT ENGINEERING FORCE OUTPUT REQUESTS FOR ',I8,1X,A,' NOT ALLOWED'                                &
                    ,/,14X,' REQUEST STRESS OUTPUT IN CASE CONTROL FOR THESE ELEMENTS')


! **********************************************************************************************************************************

      END SUBROUTINE OFP3_ELFE_2D


      SUBROUTINE OFP3_ELFN ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element node force output requests for one subcase, all element types

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_FIJ, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_ELFN_BIT, ELDT_BUG_U_P_BIT, ELDT_F25_U_P_BIT, FATAL_ERR,NELE, IBIT,   &
                                         INT_SC_NUM, MBUG, MOGEL, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  ELFORCEN, OTMSKIP
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, ETYPE, AGRID, EID, ELDT, ELGP, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,     &
                                         PEB, PEG, PEL, PLY_NUM, TYPE, SCNUM, BGRID
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, EID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFN, TXT_ELFN

      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS, TRANSFORM_NODE_FORCES
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: CALC_ELEM_NODE_FORCES
      USE ELEMENT_DIAGNOSTICS, ONLY   :  ELMOUT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE ELEMENT_OUTPUT_WRITERS, ONLY:  WRITE_ELEM_NODE_FORCE
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_FIJFIL

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_ELFN'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: DUM_BUG(0:MBUG-1) ! Values from WRT_BUG sent to subr ELMOUT in a particular call
      INTEGER(LONG)                   :: ELOUT_ELFN        ! If > 0, there are ELFORCE(NODE) requests for some elems
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: IERROR      = 0   ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: I2                ! A calculated index into an array
      INTEGER(LONG)                   :: NBUG(METYPE)      ! Count of the no. of requests for ELDATA print requests of UEL, PEL
      INTEGER(LONG)                   :: NDISK(METYPE)     ! Count of the no. of requests for ELDATA disk file requests of UEL, PEL
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_COMPS         ! Either 6 or 1 depending on whether grid is a physical grid or a SPOINT
      INTEGER(LONG)                   :: NUM_ELEM          ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_ELEM since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file


      INTRINSIC IAND



! **********************************************************************************************************************************
! Process element node force requests for all elements

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'Y'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'Y'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing no. requests/elem.
         NELREQ(I) = 0
         NBUG(I)   = 0
         NDISK(I)  = 0
      ENDDO

      DO I=1,METYPE
         DO J=1,NELE
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_ELFN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFN_BIT))
               WRT_BUG(6) = IAND(ELDT(J)            ,IBIT(ELDT_BUG_U_P_BIT))
               WRT_FIJ(5) = IAND(ELDT(J)            ,IBIT(ELDT_F25_U_P_BIT))
               IF (ELOUT_ELFN > 0) THEN
                  NELREQ(I) = NELREQ(I) + 1
               ENDIF
               IF (WRT_BUG(6) > 0) THEN
                  NBUG(I)   = NBUG(I)   + 1
               ENDIF
               IF (WRT_FIJ(5) > 0) THEN
                  NDISK(I)  = NDISK(I)  + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element nodal force'

      DO I=0,MBUG-1
         DUM_BUG(I) = 0
      ENDDO

reqs1:DO I=1,METYPE
         IF ((NELREQ(I) + NBUG(I) + NDISK(I)) == 0) CYCLE reqs1
         NUM_ELEM = 0
         NUM_OGEL = 0

elems_1: DO J = 1,NELE
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_ELFN = IAND(ELOUT(J,INT_SC_NUM), IBIT(ELOUT_ELFN_BIT))
               DUM_BUG(6) = IAND(ELDT(J)            , IBIT(ELDT_BUG_U_P_BIT))
               WRT_FIJ(5) = IAND(ELDT(J), IBIT(ELDT_F25_U_P_BIT))

               IF ((ELOUT_ELFN > 0) .OR. (DUM_BUG(6) > 0) .OR. (WRT_FIJ(5) > 0)) THEN
                  PLY_NUM = 0                              ! 'Y' in call to EMG means write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'Y' )
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE elems_1
                  ENDIF
                  CALL ELMDIS
                  CALL CALC_ELEM_NODE_FORCES
                  IF      (ELFORCEN == 'BASIC') THEN       ! Transform to basic or global coords depending on ELFORCEN
                     CALL TRANSFORM_NODE_FORCES ( 'B' )
                  ELSE IF (ELFORCEN == 'GLOBAL') THEN
                     CALL TRANSFORM_NODE_FORCES ( 'G' )
                  ENDIF
                  IF (DUM_BUG(6) > 0) THEN                 ! Output grid point displs and loads in local elem or basic coords
                     CALL ELMOUT ( J, DUM_BUG, SCNUM(INT_SC_NUM), OPT )
                  ENDIF

                  IF (ELOUT_ELFN > 0) THEN
                     NUM_ELEM = NUM_ELEM + 1
                     EID_OUT_ARRAY(NUM_ELEM,1) = EID
                     DO K=1,ELGP
                        GID_OUT_ARRAY(NUM_ELEM,K) = AGRID(K)
                     ENDDO

                     I2 = 0
                     DO K=1,ELGP
                        NUM_OGEL = NUM_OGEL + 1
                        IF (NUM_OGEL > MAXREQ) THEN
                           WRITE(ERR,9200) SUBR_NAME, MAXREQ
                           WRITE(F06,9200) SUBR_NAME, MAXREQ
                           FATAL_ERR = FATAL_ERR + 1
                           CALL OUTA_HERE ( 'Y' )          ! Coding error (dim of array OGEL too small), so quit
                        ENDIF
                        DO L=1,6
                           OGEL(NUM_OGEL,L) = ZERO
                        ENDDO
                        CALL GET_GRID_NUM_COMPS ( BGRID(K), NUM_COMPS, SUBR_NAME )
                        DO L=1,NUM_COMPS
                           I2 = I2 + 1
                           IF      (ELFORCEN == 'LOCAL') THEN
                              OGEL(NUM_OGEL,L) = PEL(I2)
                           ELSE IF (ELFORCEN == 'BASIC') THEN
                              OGEL(NUM_OGEL,L) = PEB(I2)
                           ELSE IF (ELFORCEN == 'GLOBAL') THEN
                              OGEL(NUM_OGEL,L) = PEG(I2)
                           ENDIF
                           IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFN(OT4_EROW,JVEC) = OGEL(NUM_OGEL,L)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFN(OT4_EROW), 9191) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, AGRID(K), L
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO

                     IF (NUM_ELEM == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )
                        CALL WRITE_ELEM_NODE_FORCE ( JVEC, ELGP, NUM_ELEM, IHDR )
                     ENDIF
                  ENDIF

                  IF (WRT_FIJ(5) > 0) THEN
                     CALL WRITE_FIJFIL ( 5, JVEC )
                  ENDIF

                  IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                     DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                        OT4_EROW = OT4_EROW + 1
                        WRITE(TXT_ELFN(OT4_EROW), 9199)
                     ENDDO
                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_1

      ENDDO reqs1

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT NODE FORCE'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9191 FORMAT(I8,1X,A,A8,I8,I8,I8)

 9199 FORMAT(' ')

 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE OFP3_ELFN


      SUBROUTINE OFP3_STRE_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element stress output requests for non PCOMP elements for one subcase. Also write Output Transformation Matrices (OTM's)
! for stresses for Craig-Bampton models)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRE_BIT, FATAL_ERR, IBIT, INT_SC_NUM,                                &
                                         MAX_STRESS_POINTS, MBUG, MOGEL,                                                           &
                                         NELE, NCBAR, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCHEXA8, NCHEXA20, NCPENTA6,     &
                                         NCPENTA15,NCTETRA4, NCTETRA10, NCQUAD4, NCQUAD4K, NCROD, NCSHEAR, NCTRIA3, NCTRIA3K,      &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  AGRID, ANY_STRE_OUTPUT, EDAT, EPNT, ETYPE, EID, ELGP, ELMTYP, ELOUT,                      &
                                         METYPE, NUM_SEi, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PLY_NUM, STRESS, TYPE, SHELL_STR_ANGLE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, STRE_OPT
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRE, TXT_STRE

      USE VECTOR_GEOMETRY, ONLY       :  PLANE_COORD_TRANS_21
      USE RESULT_COORDINATES, ONLY    :  TRANSFORM_SHELL_STR
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE COMPOSITE_SHELL_PREPARATION, ONLY:  IS_ELEM_PCOMP_PROPS
      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: POLYNOM_FIT_STRE_STRN
      USE ELEMENT_STRESS_RECOVERY, ONLY: CALC_ELEM_STRESSES
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE ELEMENT_OUTPUT_WRITERS, ONLY:  WRITE_ELEM_STRESSES
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_STRE_VECS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_STRESS_OUTPUT, ONLY     :  SET_OES_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRE_NO_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRESS_ITEM(20)   ! Char description of element stresses

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRE        ! If > 0, there are STRESS   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Value initialized to zero and used in call to CALC_ELEM_STRESSES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)

      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff stress points for one element (3rd dim in arrays SEi, STEi)

                                                           ! Stress index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRESS_OUT_ERR_INDEX(MAX_STRESS_POINTS)



                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRESS_OUT_PCT_ERR(MAX_STRESS_POINTS)

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRESS_OUT_PCT_ERR

                                                           ! Array of values from array STRESS for all stress points
      REAL(DOUBLE)                    :: STRESS_RAW(9,MAX_STRESS_POINTS)

                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRESS_OUT(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.

      ! OP2 stuff
      CHARACTER(8*BYTE)               :: TABLE_NAME   ! name of the op2 table name
      INTEGER(LONG)                   :: ITABLE       ! the subtable
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND
      ITABLE = 0
      TABLE_NAME = "OES ERR "

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element stress output (STRESS) requests for all elems except composite shells

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff


! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the no. requests/elem.
         NELREQ(I) = 0
      ENDDO

      DO I=1,METYPE
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  IF ((STRE_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRE_LOC == 'GAUSS   ') .OR.                                                                                &
                      (ETYPE(J)(1:4) == 'HEXA') .OR.                                                                               &
                      (ETYPE(J)(1:5) == 'PENTA') .OR.                                                                              &
                      (ETYPE(J)(1:5) == 'TETRA') .OR.                                                                              &
                      (ETYPE(J)(1:5) == 'QUAD8')) THEN
                     NUM_PTS(I) = NUM_SEi(I)
                  ELSE
                     NUM_PTS(I) = 1
                  ENDIF
                  ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
                  IF (ELOUT_STRE > 0) THEN
                     NELREQ(I) = NELREQ(I) + NUM_PTS(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      OGEL = ZERO

! 101  FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE_INT=",I8,"; TABLE_NAME=",A)
!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element stress'
reqs5:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs5
         NUM_OGEL_ROWS = 0
         NUM_OGEL      = 0

elems_5: DO J = 1,NELE

            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
               IF (ELOUT_STRE > 0) THEN
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE elems_5
                  ENDIF
                  CALL ELMDIS

                  DO M=1,NUM_PTS(I)
                     CALL ELEM_STRE_STRN_ARRAYS ( M )
                     STRESS_RAW(:,M) = STRESS(:)
                  ENDDO

                  STRESS_OUT(:,1) = STRESS(:)            ! Set STRESS_OUT for NUM_PTS(I) = 1

                  IF ((STRE_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRE_LOC == 'GAUSS   ') .OR.                                                                                &
                      (TYPE(1:4) == 'HEXA') .OR.                                                                                   &
                      (TYPE(1:5) == 'PENTA') .OR.                                                                                  &
                      (TYPE(1:5) == 'TETRA') .OR.                                                                                  &
                      (TYPE(1:5) == 'QUAD8')) THEN

                     IF (TYPE(1:5) == 'QUAD4') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,                    &
                                                     STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                     ELSE IF (TYPE(1:5) == 'QUAD8') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,                    &
                                                     STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transform stress from the cartesian local coordinate system to
                                                           ! the element coordinate system
                        DO M=2,NUM_PTS(I)
                           CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                           CALL TRANSFORM_SHELL_STR( TEL, STRESS_OUT(:,M), ONE)
                        ENDDO
                                                           ! Center stress is the average of corner stress in element coordinates.
                                                           ! This is how MSC does it.
                        STRESS_OUT(:,1) = (STRESS_OUT(:,2) + STRESS_OUT(:,3) + STRESS_OUT(:,4) + STRESS_OUT(:,5)) / FOUR

                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR.                                                                           &
                              (TYPE(1:5) == 'PENTA') .OR.                                                                          &
                              (TYPE(1:5) == 'TETRA')) THEN
! Stresses are directly evaluated at the corner grid points. If they are going to be evaluated at Gauss points
! then extrapolated to grid points, that should be done here, in POLYNOM_FIT_STRE_STRN, or in an equivalent subroutine.
                        STRESS_OUT(:,:) = STRESS_RAW(:,:)

                     ENDIF

                  ENDIF

do_stress_pts:    DO M=1,NUM_PTS(I)

                     DO K=1,9
                        STRESS(K) = STRESS_OUT(K,M)
                     ENDDO

                     CALL CALC_ELEM_STRESSES ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                                                           ! If CB soln, write rows of OGEL, from CALC_ELEM_STRESSES, to OTM_STRE
                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                        CALL GET_STRESS_ITEM_DATA

                        IF ((TYPE == 'BAR     ') .OR. (TYPE == 'TRIA3   ') .OR. (TYPE == 'QUAD4   ') .OR. (TYPE == 'SHEAR   ')) THEN
                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRE(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
                                       IF (M == 1) THEN
                                          IF (TYPE(1:5) == 'QUAD4') THEN
                                             WRITE(TXT_STRE(OT4_EROW),9190) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                   &
                                                                           STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ELSE
                                             WRITE(TXT_STRE(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                  &
                                                                             STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ENDIF
                                       ELSE
                                          WRITE(TXT_STRE(OT4_EROW),9191) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, AGRID(M-1),          &
                                                                         STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                       ENDIF
                                    ELSE
                                       WRITE(TXT_STRE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                        &
                                                                       STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ELSE
                           DO K=1,NUM_OTM_ENTRIES
                              OT4_EROW = OT4_EROW + 1
                              OTM_STRE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_STRE(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, STRESS_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_STRE(OT4_EROW), 9199)
                        ENDDO
                     ENDIF

                     NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                     EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
                     GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                     IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
                        IF (TYPE(1:5) == 'QUAD4') THEN
                           POLY_FIT_ERR(NUM_OGEL_ROWS)       = STRESS_OUT_PCT_ERR(M)
                           POLY_FIT_ERR_INDEX(NUM_OGEL_ROWS) = STRESS_OUT_ERR_INDEX(M)
                        ENDIF
                     ENDIF
                     DO K=1,ELGP
                        GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                     ENDDO

                  ENDDO do_stress_pts

                  IF (ETYPE(J)(1:5) /='USER1') THEN
                     IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )
 100                    FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                        WRITE(ERR,100) "A",TYPE,TABLE_NAME,ITABLE
                        CALL SET_OES_TABLE_NAME(TYPE, TABLE_NAME, ITABLE)
                        WRITE(ERR,100) "B",TYPE,TABLE_NAME,ITABLE
                        CALL WRITE_ELEM_STRESSES ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_PTS(I), ITABLE )
                        EXIT
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_5

      ENDDO reqs5

      IF ((TABLE_NAME .NE. "OES ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF
!===========================
      IF (WRITE_NEU .AND. (ANY_STRE_OUTPUT > 0)) THEN

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out BUSH stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBUSH, 6, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BUSH') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCBUSH, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'BUSH   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS1 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS1, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS1') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS1, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS1   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS2 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS2, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS2') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS2, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS2   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS3 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS3, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS3') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS4 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS4, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS4') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS4   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ROD stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCROD, 4, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:3) == 'ROD') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCROD, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ROD     ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out BAR stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBAR, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:3) == 'BAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCBAR, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'BAR     ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3K stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTRIA3K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TRIA3K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TRIA3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4K stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCQUAD4K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'QUAD4K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCQUAD4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'QUAD4   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA8 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA8, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA8 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCHEXA8, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'HEXA8   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA20 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA20, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA20') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCHEXA20, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'HEXA20  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA6 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA6, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA6 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCPENTA6, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PENTA6  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA15 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA15, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA15') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCPENTA15, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PENTA15 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA4 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA4, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTETRA4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TETRA4  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA10 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA10, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA10') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTETRA10, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TETRA10 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out SHEAR stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCSHEAR, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'SHEAR   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT STRESS'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9190 FORMAT(I8,1X,A,A8,I8,4X,'CENTER',9X,A20)

 9191 FORMAT(I8,1X,A,A8,I8,4X,'GRID',I8,3X,A20)

 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9193 FORMAT(I8,1X,A,A8,I8,19X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_STRESS_ITEM_DATA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II               ! DO loop index

! **********************************************************************************************************************************
      DO II=1,18
         STRESS_ITEM(II)(1:) = ' '
      ENDDO

      IF       (TYPE(1:5) == 'ELAS1') THEN
         NUM_OTM_ENTRIES = 1
         STRESS_ITEM( 1) = 'Spring elem stress  '

      ELSE IF  (TYPE(1:3) == 'ROD'  ) THEN
         NUM_OTM_ENTRIES = 4
         STRESS_ITEM( 1) = 'Axial Stress        '
         STRESS_ITEM( 2) = 'MS - Axial          '
         STRESS_ITEM( 3) = 'Torsional Stress    '
         STRESS_ITEM( 4) = 'MS - Torsion        '

      ELSE IF  (TYPE(1:3) == 'BAR'  ) THEN
         NUM_OTM_ENTRIES = 9
         STRESS_ITEM( 1) = 'SA1: Stress Pt1 EndA'  ;  STRESS_ITEM(10) = 'SB1: Stress Pt1 EndB'
         STRESS_ITEM( 2) = 'SA2: Stress Pt2 EndA'  ;  STRESS_ITEM(11) = 'SB2: Stress Pt2 EndB'
         STRESS_ITEM( 3) = 'SA3: Stress Pt3 EndA'  ;  STRESS_ITEM(12) = 'SB3: Stress Pt3 EndB'
         STRESS_ITEM( 4) = 'SA4: Stress Pt4 EndA'  ;  STRESS_ITEM(13) = 'SB4: Stress Pt4 EndB'
         STRESS_ITEM( 5) = 'Axial Stress        '  ;  STRESS_ITEM(14) = 'Axial stress        '
         STRESS_ITEM( 6) = 'SA-Max              '  ;  STRESS_ITEM(15) = 'SB-Max              '
         STRESS_ITEM( 7) = 'SA-Min              '  ;  STRESS_ITEM(16) = 'SB-Min              '
         STRESS_ITEM( 8) = 'MS-Tension          '  ;  STRESS_ITEM(17) = 'MS-Compression      '
         STRESS_ITEM( 9) = 'Torsional Stress    '  ;  STRESS_ITEM(18) = 'MS-Torsion          '

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN
         NUM_OTM_ENTRIES = 10
         STRESS_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRESS_ITEM(11) = 'Fibre Dist      +Z1 '
         STRESS_ITEM( 2) = 'Normal X Stress -Z1 '  ;  STRESS_ITEM(12) = 'Normal X Stress +Z1 '
         STRESS_ITEM( 3) = 'Normal Y Stress -Z1 '  ;  STRESS_ITEM(13) = 'Normal Y Stress +Z1 '
         STRESS_ITEM( 4) = 'Shear XY Stress -Z1 '  ;  STRESS_ITEM(14) = 'Shear XY Stress +Z1 '
         STRESS_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRESS_ITEM(15) = 'Princ Angle  at +Z1 '
         STRESS_ITEM( 6) = 'Major Stress at -Z1 '  ;  STRESS_ITEM(16) = 'Major Stress at +Z1 '
         STRESS_ITEM( 7) = 'Minor Stress at -Z1 '  ;  STRESS_ITEM(17) = 'Minor Stress at +Z1 '
         IF      (STRE_OPT(1:8) == 'VONMISES') THEN
            STRESS_ITEM( 8) = 'von Mises XY at -Z1 '  ;  STRESS_ITEM(18) = 'von Mises XY at +Z1 '
         ELSE IF (STRE_OPT(1:4) == 'MAXS'    ) THEN
            STRESS_ITEM( 8) = 'Max Shear XY at -Z1 '  ;  STRESS_ITEM(18) = 'Max Shear XY at +Z1 '
         ELSE
            STRESS_ITEM( 8) = '**** undefined **** '  ;  STRESS_ITEM(18) = '**** undefined **** '
         ENDIF
         STRESS_ITEM( 9) = 'Shear XZ Stress avg '  ;  STRESS_ITEM(19) = 'Shear XZ Stress avg '
         STRESS_ITEM(10) = 'Shear YZ Stress avg '  ;  STRESS_ITEM(20) = 'Shear XZ Stress avg '

      ELSE IF  (TYPE(1:6) == 'USERIN') THEN
         NUM_OTM_ENTRIES = 9
         STRESS_ITEM( 1) = '*** Not Defined ****'
         STRESS_ITEM( 2) = '*** Not Defined ****'
         STRESS_ITEM( 3) = '*** Not Defined ****'
         STRESS_ITEM( 4) = '*** Not Defined ****'
         STRESS_ITEM( 5) = '*** Not Defined ****'
         STRESS_ITEM( 6) = '*** Not Defined ****'
         STRESS_ITEM( 7) = '*** Not Defined ****'
         STRESS_ITEM( 8) = '*** Not Defined ****'
         STRESS_ITEM( 9) = '*** Not Defined ****'

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         NUM_OTM_ENTRIES = 8
         STRESS_ITEM( 1) = 'Normal x Stress     '
         STRESS_ITEM( 2) = 'Normal y Stress     '
         STRESS_ITEM( 3) = 'Normal z Stress     '
         STRESS_ITEM( 4) = 'Shear xy Stress     '
         STRESS_ITEM( 5) = 'Shear yz Stress     '
         STRESS_ITEM( 6) = 'Shear zx Stress     '
         STRESS_ITEM( 7) = 'Oct Direct Stress   '
         STRESS_ITEM( 8) = 'Oct Shear Stress    '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRESS_ITEM_DATA

!====================================================================================================

      END SUBROUTINE OFP3_STRE_NO_PCOMP


      SUBROUTINE OFP3_STRN_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element strain output requests for non PCOMP elements for one subcase. Also write Output Transformation Matrices (OTM's)
! for strains for Craig-Bampton models)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRN_BIT, FATAL_ERR, IBIT, INT_SC_NUM,                                &
                                         MAX_STRESS_POINTS, MBUG, MOGEL,                                                           &
                                         NELE, NCBAR, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCHEXA8, NCHEXA20, NCPENTA6,     &
                                         NCPENTA15,NCTETRA4, NCTETRA10, NCQUAD4, NCQUAD4K, NCROD, NCSHEAR, NCTRIA3, NCTRIA3K,      &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  AGRID, ANY_STRN_OUTPUT, EDAT, EPNT, ETYPE, EID, ELGP, ELMTYP, ELOUT,                      &
                                         METYPE, NUM_SEi, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PLY_NUM, STRAIN, TYPE, SHELL_STR_ANGLE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRN_OPT
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRN, TXT_STRN

      USE VECTOR_GEOMETRY, ONLY       :  PLANE_COORD_TRANS_21
      USE RESULT_COORDINATES, ONLY    :  TRANSFORM_SHELL_STR
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE COMPOSITE_SHELL_PREPARATION, ONLY:  IS_ELEM_PCOMP_PROPS
      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: POLYNOM_FIT_STRE_STRN
      USE ELEMENT_STRAIN_RECOVERY, ONLY: CALC_ELEM_STRAINS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE ELEMENT_OUTPUT_WRITERS, ONLY:  WRITE_ELEM_STRAINS
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_STRN_VECS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_STRESS_OUTPUT, ONLY     :  SET_OST_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRN_NO_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRAIN_ITEM(22)   ! Char description of element strains

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRN        ! If > 0, there are STRAIN   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Dummy valye needed in call to CALC_ELEM_ENFR_FORCES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff strain points for one element (3rd dim in arrays SEi, STEi)

                                                           ! Strain index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRAIN_OUT_ERR_INDEX(MAX_STRESS_POINTS)



                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRAIN_OUT_PCT_ERR(MAX_STRESS_POINTS)

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRAIN_OUT_PCT_ERR

                                                           ! Array of values from array STRAIN for all stress points
      REAL(DOUBLE)                    :: STRAIN_RAW(9,MAX_STRESS_POINTS)

                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRAIN_OUT(9,MAX_STRESS_POINTS)
      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.

      ! OP2 stuff
      CHARACTER(8*BYTE)               :: TABLE_NAME   ! name of the op2 table name
      INTEGER(LONG)                   :: ITABLE       ! the subtable
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND
      ITABLE = 0
      TABLE_NAME = "OES ERR "

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element strain output (STRAIN) requests for all elems except composite shells

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the no. requests/elem.
         NELREQ(I) = 0
      ENDDO

      DO I=1,METYPE                                        ! Only count requests for elem types that can have strain output
         IF((ELMTYP(I)(1:5) == 'TRIA3') .OR. (ELMTYP(I)(1:5) == 'QUAD4') .OR. (ELMTYP(I)(1:5) == 'SHEAR') .OR.                     &
            (ELMTYP(I)(1:4) == 'HEXA' ) .OR. (ELMTYP(I)(1:5) == 'PENTA') .OR. (ELMTYP(I)(1:5) == 'TETRA') .OR.                     &
            (ELMTYP(I)(1:4) == 'BUSH' ) .OR. (ELMTYP(I)(1:5) == 'QUAD8')) THEN
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'N') THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                  IF ((STRN_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRN_LOC == 'GAUSS   ') .OR.                                                                                &
                      (ETYPE(J)(1:4) == 'HEXA') .OR.                                                                               &
                      (ETYPE(J)(1:5) == 'PENTA') .OR.                                                                              &
                      (ETYPE(J)(1:5) == 'TETRA') .OR.                                                                              &
                      (ETYPE(J)(1:5) == 'QUAD8')) THEN
                        NUM_PTS(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS(I) = 1
                     ENDIF
                     ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
                     IF (ELOUT_STRN > 0) THEN
                        NELREQ(I) = NELREQ(I) + NUM_PTS(I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element strain'
reqs7:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs7
         NUM_OGEL_ROWS = 0
         NUM_OGEL = 0
elems_7: DO J = 1,NELE

            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
               IF (ELOUT_STRN > 0) THEN
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )! Calc SEi matrices
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE elems_7
                  ENDIF
                  CALL ELMDIS

                  DO M=1,NUM_PTS(I)
                     CALL ELEM_STRE_STRN_ARRAYS ( M )
                     DO K=1,9
                        STRAIN_RAW(K,M) = STRAIN(K)
                     ENDDO
                  ENDDO

                  STRAIN_OUT(:,1) = STRAIN(:)              ! Set STRAIN_OUT for NUM_PTS(I) = 1

                  IF ((STRN_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRN_LOC == 'GAUSS   ') .OR.                                                                                &
                      (TYPE(1:4) == 'HEXA') .OR.                                                                                   &
                      (TYPE(1:5) == 'PENTA') .OR.                                                                                  &
                      (TYPE(1:5) == 'TETRA') .OR.                                                                                  &
                      (TYPE(1:5) == 'QUAD8')) THEN

                     IF (TYPE(1:5) == 'QUAD4') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW, 9, NUM_PTS(I), STRAIN_OUT, STRAIN_OUT_PCT_ERR,                    &
                                                     STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX )

                     ELSE IF (TYPE(1:5) == 'QUAD8') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW, 9, NUM_PTS(I), STRAIN_OUT, STRAIN_OUT_PCT_ERR,                    &
                                                     STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transform strain from the cartesian local coordinate system to
                                                           ! the element coordinate system
                        DO M=1,NUM_PTS(I)
                           CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                           CALL TRANSFORM_SHELL_STR( TEL, STRAIN_OUT(:,M), TWO)
                        ENDDO

                                                           ! Center strain is the average of corner strain in element coordinates.
                                                           ! This is how MSC does it.
                        STRAIN_OUT(:,1) = (STRAIN_OUT(:,2) + STRAIN_OUT(:,3) + STRAIN_OUT(:,4) + STRAIN_OUT(:,5)) / FOUR

                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR.                                                                           &
                              (TYPE(1:5) == 'PENTA') .OR.                                                                          &
                              (TYPE(1:5) == 'TETRA')) THEN
! Strains are directly evaluated at the corner grid points. If they are going to be evaluated at Gauss points
! then extrapolated to grid points, that should be done here, in POLYNOM_FIT_STRE_STRN, or in an equivalent subroutine.
                        STRAIN_OUT(:,:) = STRAIN_RAW(:,:)

                     ENDIF

                  ENDIF

do_strain_pts:    DO M=1,NUM_PTS(I)

                     DO K=1,9
                        STRAIN(K) = STRAIN_OUT(K,M)
                     ENDDO
                     CALL CALC_ELEM_STRAINS ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                                                              ! If CB soln, write rows of OGEL, from CALC_ELEM_STRAINS, to OTM_STRN
                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                        CALL GET_STRAIN_ITEM_DATA

                        IF ((TYPE == 'TRIA3   ') .OR. (TYPE == 'QUAD4   ') .OR. (TYPE == 'SHEAR   ')) THEN
                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRN(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
                                       IF (M == 1) THEN
                                          IF (TYPE(1:5) == 'QUAD4') THEN
                                             WRITE(TXT_STRN(OT4_EROW),9190) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                   &
                                                                           STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ELSE
                                             WRITE(TXT_STRN(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                  &
                                                                             STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ENDIF
                                       ELSE
                                          WRITE(TXT_STRN(OT4_EROW),9191) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, AGRID(M-1),          &
                                                                         STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                       ENDIF
                                    ELSE
                                       WRITE(TXT_STRN(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                        &
                                                                       STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ELSE
                           DO K=1,NUM_OTM_ENTRIES
                              OT4_EROW = OT4_EROW + 1
                              OTM_STRN(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_STRN(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, STRAIN_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_STRN(OT4_EROW), 9199)
                        ENDDO
                     ENDIF

                     NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                     EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
                     GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                     IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
                        IF (TYPE(1:5) == 'QUAD4') THEN
                           POLY_FIT_ERR(NUM_OGEL_ROWS)       = STRAIN_OUT_PCT_ERR(M)
                           POLY_FIT_ERR_INDEX(NUM_OGEL_ROWS) = STRAIN_OUT_ERR_INDEX(M)
                        ENDIF
                     ENDIF
                     DO K=1,ELGP
                        GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                     ENDDO

                  ENDDO do_strain_pts

                  IF (ETYPE(J)(1:5) /='USER1') THEN
                     IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )
 100                    FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                        WRITE(ERR,100) "A",TYPE,TABLE_NAME,ITABLE
                        CALL SET_OST_TABLE_NAME(TYPE, TABLE_NAME, ITABLE)
                        WRITE(ERR,100) "B",TYPE,TABLE_NAME,ITABLE
                        CALL WRITE_ELEM_STRAINS ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_PTS(I), ITABLE )
                        EXIT
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_7

      ENDDO reqs7

      IF ((TABLE_NAME .NE. "OES ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF
!===========================
      IF (WRITE_NEU .AND. (ANY_STRN_OUTPUT > 0)) THEN

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out BUSH strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBUSH, 6, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BUSH') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCBUSH, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'BUSH  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3K strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTRIA3K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TRIA3K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TRIA3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4K strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCQUAD4K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'QUAD4K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCQUAD4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'QUAD4   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA8 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA8, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA8 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCHEXA8, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'HEXA8   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA20 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA20, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA20') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCHEXA20, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'HEXA20  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA6 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA6, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA6 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCPENTA6, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'PENTA6  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA15 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA15, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA15') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCPENTA15, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'PENTA15 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA4 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA4, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTETRA4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TETRA4  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA10 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA10, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA10') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTETRA10, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TETRA10 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0                                      ! Write out SHEAR strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCSHEAR , NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'SHEAR   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT STRAIN'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9190 FORMAT(I8,1X,A,A8,I8,4X,'CENTER',9X,A20)

 9191 FORMAT(I8,1X,A,A8,I8,4X,'GRID',I8,3X,A20)

 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9193 FORMAT(I8,1X,A,A8,I8,19X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_STRAIN_ITEM_DATA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II               ! DO loop index

! **********************************************************************************************************************************
      DO II=1,18
         STRAIN_ITEM(II)(1:) = ' '
      ENDDO

      IF  (TYPE(1:3) == 'ROD'  ) THEN
         NUM_OTM_ENTRIES = 4
         STRAIN_ITEM( 1) = 'Axial Strain        '
         STRAIN_ITEM( 2) = 'MS - Axial          '
         STRAIN_ITEM( 3) = 'Torsional Strain    '
         STRAIN_ITEM( 4) = 'MS - Torsion        '

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN
         NUM_OTM_ENTRIES = 10
         STRAIN_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRAIN_ITEM(11) = 'Fibre Dist      +Z1 '
         STRAIN_ITEM( 2) = 'Normal X Strain -Z1 '  ;  STRAIN_ITEM(12) = 'Normal X Strain +Z1 '
         STRAIN_ITEM( 3) = 'Normal Y Strain -Z1 '  ;  STRAIN_ITEM(13) = 'Normal Y Strain +Z1 '
         STRAIN_ITEM( 4) = 'Shear XY Strain -Z1 '  ;  STRAIN_ITEM(14) = 'Shear XY Strain +Z1 '
         STRAIN_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRAIN_ITEM(15) = 'Princ Angle  at +Z1 '
         STRAIN_ITEM( 6) = 'Major Strain at -Z1 '  ;  STRAIN_ITEM(16) = 'Major Strain at +Z1 '
         STRAIN_ITEM( 7) = 'Minor Strain at -Z1 '  ;  STRAIN_ITEM(17) = 'Minor Strain at +Z1 '
         IF      (STRN_OPT(1:8) == 'VONMISES') THEN
            STRAIN_ITEM( 8) = 'von Mises XY at -Z1 '  ;  STRAIN_ITEM(18) = 'von Mises XY at +Z1 '
         ELSE IF (STRN_OPT(1:4) == 'MAXS'    ) THEN
            STRAIN_ITEM( 8) = 'Max Shear XY at -Z1 '  ;  STRAIN_ITEM(18) = 'Max Shear XY at +Z1 '
         ELSE
            STRAIN_ITEM( 8) = '**** undefined **** '  ;  STRAIN_ITEM(18) = '**** undefined **** '
         ENDIF
         STRAIN_ITEM( 9) = 'Transv Shear XZ avg '  ;  STRAIN_ITEM(19) = 'Transv Shear XZ avg '
         STRAIN_ITEM(10) = 'Transv Shear YZ avg '  ;  STRAIN_ITEM(20) = 'Transv Shear YZ avg '

      ELSE IF  (TYPE(1:6) == 'USERIN') THEN
         NUM_OTM_ENTRIES = 9
         STRAIN_ITEM( 1) = '*** Not Defined ****'
         STRAIN_ITEM( 2) = '*** Not Defined ****'
         STRAIN_ITEM( 3) = '*** Not Defined ****'
         STRAIN_ITEM( 4) = '*** Not Defined ****'
         STRAIN_ITEM( 5) = '*** Not Defined ****'
         STRAIN_ITEM( 6) = '*** Not Defined ****'
         STRAIN_ITEM( 7) = '*** Not Defined ****'
         STRAIN_ITEM( 8) = '*** Not Defined ****'
         STRAIN_ITEM( 9) = '*** Not Defined ****'

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         NUM_OTM_ENTRIES = 8
         STRAIN_ITEM( 1) = 'Normal x Strain     '
         STRAIN_ITEM( 2) = 'Normal y Strain     '
         STRAIN_ITEM( 3) = 'Normal z Strain     '
         STRAIN_ITEM( 4) = 'Shear xy Strain     '
         STRAIN_ITEM( 5) = 'Shear yz Strain     '
         STRAIN_ITEM( 6) = 'Shear zx Strain     '
         STRAIN_ITEM( 7) = 'Oct Direct Strain   '
         STRAIN_ITEM( 8) = 'Oct Shear Strain    '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRAIN_ITEM_DATA

      END SUBROUTINE OFP3_STRN_NO_PCOMP

   END MODULE ELEMENT_RESULT_OUTPUT
