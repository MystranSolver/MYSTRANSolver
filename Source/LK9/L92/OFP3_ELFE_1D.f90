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
                                         ETYPE, EID, ELMTYP, ELOUT, FCONV, METYPE, NUM_EMG_FATAL_ERRS, OFFDIS_GA_GB, OFFDIS_L,     &
                                         PE_GA_GB, PEL, PLY_NUM, STRESS, TE, TE_GA_GB, TYPE, XEL
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFE, TXT_ELFE
!--- cbeam_add begin --- ! 
      USE SCONTR, ONLY                :  NCBEAM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  AGRID, CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_XL, PBEAM_NSTATIONS, PRESS
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT,GID_OUT_ARRAY
!--- cbeam_add end --- ! 

      USE OFP3_ELFE_1D_USE_IFs

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
! --- cbeam_add begin --- !
      INTEGER(LONG)                   :: NREQ_ELEM         ! Actual number of engineering-force rows for current element
      INTEGER(LONG)                   :: NSTA             ! Active number of CBEAM stations for the current element
      INTEGER(LONG)                   :: ISTA             ! CBEAM station loop index
      INTEGER(LONG)                   :: WRITE_NUM_PTS    ! NUM_PTS argument for WRITE_ELEM_ENGR_FORCE
      INTEGER(LONG)                   :: OGEL_ENDA_ROW    ! OGEL row for the x/L = 0.0 CBEAM endpoint
      INTEGER(LONG)                   :: OGEL_ENDB_ROW    ! OGEL row for the x/L = 1.0 CBEAM endpoint
!                                                            (this can be > NUM_ELEM since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file
! --- cbeam_add end --- !

      REAL(DOUBLE)                    :: DUM0(6,12)        ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM1(6)           ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM21(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM22(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM31(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM32(3)          ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: EEF(6) = ZERO     ! Element engineering force for BUSH/1D force work arrays
      REAL(DOUBLE)                    :: DX,DY,DZ          ! Offset dist1
      REAL(DOUBLE)                    :: FORCES(12)        ! Forces at the grid points
      REAL(DOUBLE)                    :: LENGTH
! --- cbeam_add begin --- !
      REAL(DOUBLE)                    :: DT
      REAL(DOUBLE)                    :: GRAD
      REAL(DOUBLE)                    :: M1_CORR
      REAL(DOUBLE)                    :: M2_CORR
      REAL(DOUBLE)                    :: P1_LOC
      REAL(DOUBLE)                    :: P2_LOC
      REAL(DOUBLE)                    :: T_CUT
      REAL(DOUBLE)                    :: V1_CORR
      REAL(DOUBLE)                    :: V2_CORR
      REAL(DOUBLE)                    :: X1L
      REAL(DOUBLE)                    :: X2L
      REAL(DOUBLE)                    :: XI_REL
      REAL(DOUBLE)                    :: XI_SPAN
      REAL(DOUBLE)                    :: XI_STA
      REAL(DOUBLE)                    :: XI_TOL
! --- cbeam_add end --- !
      ! OP2 parameters
      INTEGER(LONG)                   :: ITABLE            ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME        ! the op2 table name

      REAL(DOUBLE)                    :: TET(3,3)          ! Transpose of TE
      REAL(DOUBLE)                    :: TET_GA_GB(3,3)    ! Transpose of TE_GA_GB
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND

! **********************************************************************************************************************************
! --- op2_upgraded begin --- !
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
      OPT(5) = 'Y'                                         ! OPT(5) is for calc of PPE 
! --- cbeam_add begin --- ! >> change to 'Y' from 'N' PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff
! --- cbeam_add end --- !
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
! added CBEAM >> (ETYPE(J)(1:4) == 'BEAM') 

            IF ((ETYPE(J)(1:3) == 'BAR') .OR. (ETYPE(J)(1:4) == 'BEAM') .OR. (ETYPE(J)(1:4) == 'BUSH') .OR.                       &
                (ETYPE(J)(1:4) == 'ELAS') .OR. (ETYPE(J)(1:3) == 'ROD'))THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                  IF (ELOUT_ELFE > 0) THEN

! --- cbeam_add begin --- !
                     NREQ_ELEM = 1
                     IF (ETYPE(J)(1:4) == 'BEAM') THEN
                        NREQ_ELEM = PBEAM_NSTATIONS(EDAT(EPNT(J)+1))
                        IF (NREQ_ELEM <= 0) NREQ_ELEM = 5
                     ENDIF
                     NELREQ(I) = NELREQ(I) + NREQ_ELEM 
! --- cbeam_end begin --- ! 
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
! --- cbeam_add begin --- !
         WRITE_NUM_PTS = 1
! --- cbeam_add end --- !
elems_2: DO J = 1,NELE
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
! added CBEAM 
            IF ((ETYPE(J)(1:3) == 'BAR') .OR. (ETYPE(J)(1:4) == 'BEAM') .OR. (ETYPE(J)(1:4) == 'BUSH') .OR.                       &
                (ETYPE(J)(1:4) == 'ELAS') .OR. (ETYPE(J)(1:3) == 'ROD'))THEN

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
! --- cbeam_add begin --- ! 
                     NSTA = 1
                     IF (ETYPE(J)(1:4) == 'BEAM') THEN
                        NSTA = CBEAM_ACTIVE_NSTATIONS
                        IF (NSTA <= 0) NSTA = 5
                        IF (DEBUG(233) > 0) THEN
                           WRITE(ERR,'("DBG_BEAM_PRE EID=",I8," NSTA=",I8," NUM_OGEL=",I8," NUM_ELEM=",I8)') EID, NSTA, NUM_OGEL, NUM_ELEM
                        ENDIF
                     ENDIF
! --- cbeam_add end --- !
!                    ---------------------------------------------------------------------------------------------------------------
                     IF (ETYPE(J)(1:4) == 'ELAS') THEN     ! Set engr forces based on the node force values	
                        NUM_OGEL = NUM_OGEL + 1
                        IF (NUM_OGEL > MAXREQ) THEN
                           WRITE(ERR,9200) SUBR_NAME, MAXREQ
                           WRITE(F06,9200) SUBR_NAME, MAXREQ
                           FATAL_ERR = FATAL_ERR + 1
                           CALL OUTA_HERE ( 'Y' )             ! Coding error (dim of array OGEL too small), so quit
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        NDUM = 0
                        CALL CALC_ELEM_STRESSES ( 1, NDUM, 0, 'N', 'N' )
                        IF (FCONV(1) > ZERO) THEN          ! ELAS engr force is stress/FCONV
                           OGEL(NUM_OGEL,1) = STRESS(1)/FCONV(1)
                        ELSE
                           OGEL(NUM_OGEL,1) = ZERO
                        ENDIF
                     ! end elas
!                    ---------------------------------------------------------------------------------------------------------------
                     ELSE IF (ETYPE(J)(1:4) == 'BUSH') THEN
! added OGEL for not crashing	
                        NUM_OGEL = NUM_OGEL + 1
                        IF (NUM_OGEL > MAXREQ) THEN
                           WRITE(ERR,9200) SUBR_NAME, MAXREQ
                           WRITE(F06,9200) SUBR_NAME, MAXREQ
                           FATAL_ERR = FATAL_ERR + 1
                           CALL OUTA_HERE ( 'Y' )
                        ENDIF
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
! added OGEL for not crashing	
                        NUM_OGEL = NUM_OGEL + 1
                        IF (NUM_OGEL > MAXREQ) THEN
                           WRITE(ERR,9200) SUBR_NAME, MAXREQ
                           WRITE(F06,9200) SUBR_NAME, MAXREQ
                           FATAL_ERR = FATAL_ERR + 1
                           CALL OUTA_HERE ( 'Y' )
                        ENDIF
                        OGEL(NUM_OGEL,7) = -PEL(1)         ! Fx  (axial force for ROD)
                        OGEL(NUM_OGEL,8) = -PEL(4)         ! T   (torque for ROD)
                     !end rod
!                    ---------------------------------------------------------------------------------------------------------------
                     ELSE IF (ETYPE(J)(1:3) == 'BAR') THEN
! added OGEL for not crashing	
                        NUM_OGEL = NUM_OGEL + 1
                        IF (NUM_OGEL > MAXREQ) THEN
                           WRITE(ERR,9200) SUBR_NAME, MAXREQ
                           WRITE(F06,9200) SUBR_NAME, MAXREQ
                           FATAL_ERR = FATAL_ERR + 1
                           CALL OUTA_HERE ( 'Y' )
                        ENDIF
                        LENGTH = ELEM_LEN_AB
                        OGEL(NUM_OGEL,1) = -PEL(6)                 ! M1a (bending moment, plane 1, end a for BAR)
                        OGEL(NUM_OGEL,2) =  PEL(5)                 ! M2a (bending moment, plane 2, end a for BAR)
                        OGEL(NUM_OGEL,3) = -PEL(6) + PEL(2)*LENGTH ! M1b (bending moment, plane 1, end b for BAR)
                        OGEL(NUM_OGEL,4) =  PEL(5) + PEL(3)*LENGTH ! M2b (bending moment, plane 2, end b for BAR)
                        OGEL(NUM_OGEL,5) = -PEL(2)                 ! V1  (plane 1 shear for BAR)
                        OGEL(NUM_OGEL,6) = -PEL(3)                 ! V2  (plane 2 shear for BAR)
                        OGEL(NUM_OGEL,7) = -PEL(1)                 ! Fx  (axial force for BAR)
                        OGEL(NUM_OGEL,8) = -PEL(4)                 ! T   (torque for BAR)
! --- cbeam_add begin --- !
! --- cbeam_stations begin--- !
                     ELSE IF (ETYPE(J)(1:4) == 'BEAM') THEN
                        LENGTH = ELEM_LEN_AB
                        DO ISTA=1,NSTA
                           NUM_OGEL = NUM_OGEL + 1
                           IF (NUM_OGEL > MAXREQ) THEN
                              WRITE(ERR,9200) SUBR_NAME, MAXREQ
                              WRITE(F06,9200) SUBR_NAME, MAXREQ
                              FATAL_ERR = FATAL_ERR + 1
                              CALL OUTA_HERE ( 'Y' )
                           ENDIF
                           XI_STA = CBEAM_ACTIVE_XL(ISTA)
                           IF (NSTA == 1) XI_STA = ZERO
                           M1_CORR = ZERO
                           M2_CORR = ZERO
                           V1_CORR = ZERO
                           V2_CORR = ZERO

                           P1_LOC = PRESS(1,INT_SC_NUM)
                           P2_LOC = PRESS(2,INT_SC_NUM)
                           X1L    = PRESS(3,INT_SC_NUM)
                           X2L    = PRESS(4,INT_SC_NUM)
                           IF (PRESS(25,INT_SC_NUM) > HALF) THEN
                              P1_LOC = TE(2,2)*P1_LOC
                              P2_LOC = TE(2,2)*P2_LOC
                           ENDIF
                           IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                              IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                                 V1_CORR = V1_CORR - P1_LOC
                                 M1_CORR = M1_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                              ELSE
                                 T_CUT = MIN(XI_STA, X2L)
                                 IF (T_CUT > X1L) THEN
                                    DT     = T_CUT - X1L
                                    XI_REL = XI_STA - X1L
                                    XI_SPAN= X2L - X1L
                                    GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                                    V1_CORR = V1_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                                    M1_CORR = M1_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT)                        &
                                              + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                                 ENDIF
                              ENDIF
                           ENDIF

                           P1_LOC = PRESS(5,INT_SC_NUM)
                           P2_LOC = PRESS(6,INT_SC_NUM)
                           X1L    = PRESS(7,INT_SC_NUM)
                           X2L    = PRESS(8,INT_SC_NUM)
                           IF (PRESS(26,INT_SC_NUM) > HALF) THEN
                              P1_LOC = TE(2,3)*P1_LOC
                              P2_LOC = TE(2,3)*P2_LOC
                           ENDIF
                           IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                              IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                                 V1_CORR = V1_CORR - P1_LOC
                                 M1_CORR = M1_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                              ELSE
                                 T_CUT = MIN(XI_STA, X2L)
                                 IF (T_CUT > X1L) THEN
                                    DT     = T_CUT - X1L
                                    XI_REL = XI_STA - X1L
                                    XI_SPAN= X2L - X1L
                                    GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                                    V1_CORR = V1_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                                    M1_CORR = M1_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT)                        &
                                              + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                                 ENDIF
                              ENDIF
                           ENDIF

                           P1_LOC = PRESS(5,INT_SC_NUM)
                           P2_LOC = PRESS(6,INT_SC_NUM)
                           X1L    = PRESS(7,INT_SC_NUM)
                           X2L    = PRESS(8,INT_SC_NUM)
                           IF (PRESS(26,INT_SC_NUM) > HALF) THEN
                              P1_LOC = TE(3,3)*P1_LOC
                              P2_LOC = TE(3,3)*P2_LOC
                           ENDIF
                           IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                              IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                                 V2_CORR = V2_CORR - P1_LOC
                                 M2_CORR = M2_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                              ELSE
                                 T_CUT = MIN(XI_STA, X2L)
                                 IF (T_CUT > X1L) THEN
                                    DT     = T_CUT - X1L
                                    XI_REL = XI_STA - X1L
                                    XI_SPAN= X2L - X1L
                                    GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                                    V2_CORR = V2_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                                    M2_CORR = M2_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT)                        &
                                              + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                                 ENDIF
                              ENDIF
                           ENDIF

                           P1_LOC = PRESS(1,INT_SC_NUM)
                           P2_LOC = PRESS(2,INT_SC_NUM)
                           X1L    = PRESS(3,INT_SC_NUM)
                           X2L    = PRESS(4,INT_SC_NUM)
                           IF (PRESS(25,INT_SC_NUM) > HALF) THEN
                              P1_LOC = TE(3,2)*P1_LOC
                              P2_LOC = TE(3,2)*P2_LOC
                           ENDIF
                           IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                              IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                                 V2_CORR = V2_CORR - P1_LOC
                                 M2_CORR = M2_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                              ELSE
                                 T_CUT = MIN(XI_STA, X2L)
                                 IF (T_CUT > X1L) THEN
                                    DT     = T_CUT - X1L
                                    XI_REL = XI_STA - X1L
                                    XI_SPAN= X2L - X1L
                                    GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                                    V2_CORR = V2_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                                    M2_CORR = M2_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT)                        &
                                              + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                                 ENDIF
                              ENDIF
                           ENDIF

                            OGEL(NUM_OGEL,1) = -PEL(6) + PEL(2)*LENGTH*XI_STA + M1_CORR
                            OGEL(NUM_OGEL,2) =  PEL(5) + PEL(3)*LENGTH*XI_STA + M2_CORR
                            OGEL(NUM_OGEL,3) = OGEL(NUM_OGEL,1)
                            OGEL(NUM_OGEL,4) = OGEL(NUM_OGEL,2)
                            OGEL(NUM_OGEL,5) = -PEL(2) + V1_CORR
                            OGEL(NUM_OGEL,6) = -PEL(3) + V2_CORR
                           OGEL(NUM_OGEL,7) = -PEL(1)
                           OGEL(NUM_OGEL,8) = -PEL(4)
                           CBEAM_XL_OUT(NUM_OGEL) = XI_STA
                           EID_OUT_ARRAY(NUM_OGEL,1) = EID
                           GID_OUT_ARRAY(NUM_OGEL,1) = 0
                           GID_OUT_ARRAY(NUM_OGEL,2) = AGRID(1)
                           GID_OUT_ARRAY(NUM_OGEL,3) = AGRID(2)
                           IF ((DEBUG(233) > 0) .AND. ((EID == 14) .OR. (EID == 15))) THEN
                              WRITE(F06,9233) EID, ISTA, XI_STA, PEL(2), PEL(3), M1_CORR, M2_CORR, V1_CORR, V2_CORR,                &
                                               OGEL(NUM_OGEL,1), OGEL(NUM_OGEL,2), OGEL(NUM_OGEL,5), OGEL(NUM_OGEL,6)
                           ENDIF
                        ENDDO
! --- cbeam_add end --- !
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
! --- cbeam_add begin --- !
                        IF (ETYPE(J)(1:4) == 'BEAM') THEN
                           DO K=1,8
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL-NSTA+1,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF
! --- cbeam_add end --- !
                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                     ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_ELFE(OT4_EROW), 9199)
                        ENDDO
                     ENDIF
! --- cbeam_add begin --- !
                     IF (ETYPE(J)(1:4) == 'BEAM') THEN
                        NUM_ELEM = NUM_ELEM + NSTA
                        WRITE_NUM_PTS = NSTA
                        IF (DEBUG(233) > 0) THEN
                           WRITE(ERR,'("DBG_BEAM_POST EID=",I8," NSTA=",I8," NUM_OGEL=",I8," NUM_ELEM=",I8)') EID, NSTA, NUM_OGEL, NUM_ELEM
                        ENDIF
                     ELSE
                        NUM_ELEM = NUM_ELEM + 1
                        EID_OUT_ARRAY(NUM_ELEM,1) = EID

                       IF (NUM_ELEM == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )

                        WRITE(ERR,100) "F1",ETYPE(J),TABLE_NAME,ITABLE
                        CALL SET_OEF_TABLE_NAME(ETYPE(J), TABLE_NAME, ITABLE)
                        WRITE(ERR,100) "F2",ETYPE(J),TABLE_NAME,ITABLE
                        CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_ELEM, IHDR, 1, ITABLE )
                        EXIT
                       ENDIF
                     ENDIF

                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_2
! added OGEL fix  
         IF (NUM_ELEM > 0) THEN
            CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100        FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
            TYPE = ELMTYP(I)
            IF (DEBUG(200) > 0) THEN
               WRITE(ERR,100) "F1",TYPE,TABLE_NAME,ITABLE
            ENDIF
            CALL SET_OEF_TABLE_NAME(TYPE, TABLE_NAME, ITABLE)
            IF (DEBUG(200) > 0) THEN
               WRITE(ERR,100) "F2",TYPE,TABLE_NAME,ITABLE
            ENDIF
            CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_ELEM, IHDR, WRITE_NUM_PTS, ITABLE )
         ENDIF

      ENDDO reqs2
 10   FORMAT("*DEBUG:      OEF_END 1D:    TABLE_NAME",A)
      IF (DEBUG(200) > 0) THEN
         WRITE(ERR,10) TABLE_NAME
      ENDIF
      IF ((TABLE_NAME .NE. "OEF ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF
! --- op2_upgraded end --- !

      IF (WRITE_NEU .AND. (ANY_ELFE_OUTPUT > 0)) THEN

! beam
! --- cbeam_add begin ---!  
         NUM_FROWS= 0
! --- neu_upgrade begin --- !
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBEAM, 14, SUBR_NAME )
         DO J=1,NELE
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:4) == 'BEAM') THEN
               NUM_FROWS= NUM_FROWS+ 1
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL CALC_ELEM_NODE_FORCES
               LENGTH = ELEM_LEN_AB
               XI_STA = ZERO
               M1_CORR = ZERO
               M2_CORR = ZERO
               V1_CORR = ZERO
               V2_CORR = ZERO

               P1_LOC = PRESS(1,INT_SC_NUM)
               P2_LOC = PRESS(2,INT_SC_NUM)
               X1L    = PRESS(3,INT_SC_NUM)
               X2L    = PRESS(4,INT_SC_NUM)
               IF (PRESS(25,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(2,2)*P1_LOC
                  P2_LOC = TE(2,2)*P2_LOC
               ENDIF

               P1_LOC = PRESS(5,INT_SC_NUM)
               P2_LOC = PRESS(6,INT_SC_NUM)
               X1L    = PRESS(7,INT_SC_NUM)
               X2L    = PRESS(8,INT_SC_NUM)
               IF (PRESS(26,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(2,3)*P1_LOC
                  P2_LOC = TE(2,3)*P2_LOC
               ENDIF

               XI_STA = 1.0D0
               M1_CORR = ZERO
               M2_CORR = ZERO
               V1_CORR = ZERO
               V2_CORR = ZERO

               P1_LOC = PRESS(1,INT_SC_NUM)
               P2_LOC = PRESS(2,INT_SC_NUM)
               X1L    = PRESS(3,INT_SC_NUM)
               X2L    = PRESS(4,INT_SC_NUM)
               IF (PRESS(25,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(2,2)*P1_LOC
                  P2_LOC = TE(2,2)*P2_LOC
               ENDIF
               IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                  IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                     V1_CORR = V1_CORR - P1_LOC
                     M1_CORR = M1_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                  ELSE
                     T_CUT = MIN(XI_STA, X2L)
                     IF (T_CUT > X1L) THEN
                        DT     = T_CUT - X1L
                        XI_REL = XI_STA - X1L
                        XI_SPAN= X2L - X1L
                        GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                        V1_CORR = V1_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                        M1_CORR = M1_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT) + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                     ENDIF
                  ENDIF
               ENDIF

               P1_LOC = PRESS(5,INT_SC_NUM)
               P2_LOC = PRESS(6,INT_SC_NUM)
               X1L    = PRESS(7,INT_SC_NUM)
               X2L    = PRESS(8,INT_SC_NUM)
               IF (PRESS(26,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(2,3)*P1_LOC
                  P2_LOC = TE(2,3)*P2_LOC
               ENDIF
               IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                  IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                     V1_CORR = V1_CORR - P1_LOC
                     M1_CORR = M1_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                  ELSE
                     T_CUT = MIN(XI_STA, X2L)
                     IF (T_CUT > X1L) THEN
                        DT     = T_CUT - X1L
                        XI_REL = XI_STA - X1L
                        XI_SPAN= X2L - X1L
                        GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                        V1_CORR = V1_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                        M1_CORR = M1_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT) + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                     ENDIF
                  ENDIF
               ENDIF

               P1_LOC = PRESS(5,INT_SC_NUM)
               P2_LOC = PRESS(6,INT_SC_NUM)
               X1L    = PRESS(7,INT_SC_NUM)
               X2L    = PRESS(8,INT_SC_NUM)
               IF (PRESS(26,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(3,3)*P1_LOC
                  P2_LOC = TE(3,3)*P2_LOC
               ENDIF
               IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                  IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                     V2_CORR = V2_CORR - P1_LOC
                     M2_CORR = M2_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                  ELSE
                     T_CUT = MIN(XI_STA, X2L)
                     IF (T_CUT > X1L) THEN
                        DT     = T_CUT - X1L
                        XI_REL = XI_STA - X1L
                        XI_SPAN= X2L - X1L
                        GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                        V2_CORR = V2_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                        M2_CORR = M2_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT) + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                     ENDIF
                  ENDIF
               ENDIF

               P1_LOC = PRESS(1,INT_SC_NUM)
               P2_LOC = PRESS(2,INT_SC_NUM)
               X1L    = PRESS(3,INT_SC_NUM)
               X2L    = PRESS(4,INT_SC_NUM)
               IF (PRESS(25,INT_SC_NUM) > HALF) THEN
                  P1_LOC = TE(3,2)*P1_LOC
                  P2_LOC = TE(3,2)*P2_LOC
               ENDIF
               IF ((X1L >= ZERO) .AND. (XI_STA > X1L)) THEN
                  IF (DABS(X2L - X1L) <= 1.0D-12) THEN
                     V2_CORR = V2_CORR - P1_LOC
                     M2_CORR = M2_CORR + P1_LOC*LENGTH*(XI_STA - X1L)
                  ELSE
                     T_CUT = MIN(XI_STA, X2L)
                     IF (T_CUT > X1L) THEN
                        DT     = T_CUT - X1L
                        XI_REL = XI_STA - X1L
                        XI_SPAN= X2L - X1L
                        GRAD   = (P2_LOC - P1_LOC)/XI_SPAN
                        V2_CORR = V2_CORR - LENGTH*(P1_LOC*DT + 0.5D0*GRAD*DT*DT)
                        M2_CORR = M2_CORR + LENGTH*LENGTH*(P1_LOC*(XI_REL*DT - 0.5D0*DT*DT) + GRAD*(0.5D0*XI_REL*DT*DT - DT*DT*DT/3.0D0))
                     ENDIF
                  ENDIF
               ENDIF

               FEMAP_EL_VECS(NUM_FROWS,1)  = -PEL(6)
               FEMAP_EL_VECS(NUM_FROWS,2)  =  PEL(5)
               FEMAP_EL_VECS(NUM_FROWS,3)  = -PEL(6) + PEL(2)*LENGTH + M1_CORR
               FEMAP_EL_VECS(NUM_FROWS,4)  =  PEL(5) + PEL(3)*LENGTH + M2_CORR
               FEMAP_EL_VECS(NUM_FROWS,5)  = -PEL(2)
               FEMAP_EL_VECS(NUM_FROWS,6)  = -PEL(3)
               FEMAP_EL_VECS(NUM_FROWS,7)  = -PEL(2) + V1_CORR
               FEMAP_EL_VECS(NUM_FROWS,8)  = -PEL(3) + V2_CORR
               FEMAP_EL_VECS(NUM_FROWS,9)  = -PEL(1)
               FEMAP_EL_VECS(NUM_FROWS,10) = -PEL(1)
               FEMAP_EL_VECS(NUM_FROWS,11) = -PEL(4)
               FEMAP_EL_VECS(NUM_FROWS,12) = -PEL(4)
               FEMAP_EL_VECS(NUM_FROWS,13) = ZERO
               FEMAP_EL_VECS(NUM_FROWS,14) = ZERO

! --- CBEAM_standard begin --- !
               OGEL_ENDA_ROW = 0
               OGEL_ENDB_ROW = 0
               XI_TOL        = 1.0D-10
               DO K=1,MIN(NUM_OGEL,MAXREQ)
                  IF (EID_OUT_ARRAY(K,1) == EID) THEN
                     IF (DABS(CBEAM_XL_OUT(K) - ZERO ) <= XI_TOL) OGEL_ENDA_ROW = K
                     IF (DABS(CBEAM_XL_OUT(K) - 1.0D0) <= XI_TOL) OGEL_ENDB_ROW = K
                  ENDIF
               ENDDO
               IF (OGEL_ENDA_ROW > 0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,1)  = OGEL(OGEL_ENDA_ROW,1)
                  FEMAP_EL_VECS(NUM_FROWS,2)  = OGEL(OGEL_ENDA_ROW,2)
                  FEMAP_EL_VECS(NUM_FROWS,5)  = OGEL(OGEL_ENDA_ROW,5)
                  FEMAP_EL_VECS(NUM_FROWS,6)  = OGEL(OGEL_ENDA_ROW,6)
                  FEMAP_EL_VECS(NUM_FROWS,9)  = OGEL(OGEL_ENDA_ROW,7)
                  FEMAP_EL_VECS(NUM_FROWS,11) = OGEL(OGEL_ENDA_ROW,8)
               ENDIF
               IF (OGEL_ENDB_ROW > 0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,3)  = OGEL(OGEL_ENDB_ROW,1)
                  FEMAP_EL_VECS(NUM_FROWS,4)  = OGEL(OGEL_ENDB_ROW,2)
                  FEMAP_EL_VECS(NUM_FROWS,7)  = OGEL(OGEL_ENDB_ROW,5)
                  FEMAP_EL_VECS(NUM_FROWS,8)  = OGEL(OGEL_ENDB_ROW,6)
                  FEMAP_EL_VECS(NUM_FROWS,10) = OGEL(OGEL_ENDB_ROW,7)
                  FEMAP_EL_VECS(NUM_FROWS,12) = OGEL(OGEL_ENDB_ROW,8)
               ENDIF
! --- CBEAM_standard end --- !
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'BEAM    ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA
! --- neu_upgrade end --- !
! --- cbeam_add begin ---! 

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

! For ELAS we need to calculate elem engr forces from the stresses since there is no "local" elem coord system

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
               IF (FCONV(1) > 0.D0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)/FCONV(1)
               ELSE

               ENDIF
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
               IF (FCONV(1) > 0.D0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)/FCONV(1)
               ELSE

               ENDIF
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
               IF (FCONV(1) > 0.D0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)/FCONV(1)
               ELSE

               ENDIF
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
               IF (FCONV(1) > 0.D0) THEN
                  FEMAP_EL_VECS(NUM_FROWS,1) = STRESS(1)/FCONV(1)
               ELSE

               ENDIF
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

 9233 FORMAT(' DEBUG233 CBEAM ELFO EID=',I8,' ISTA=',I4,' XI=',1ES14.6,' PEL23=',2(1X,1ES14.6),                                   &
             ' CORR(M1,M2,V1,V2)=',4(1X,1ES14.6),' OUT(M1,M2,V1,V2)=',4(1X,1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE OFP3_ELFE_1D



