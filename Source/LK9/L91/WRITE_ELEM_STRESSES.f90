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

      SUBROUTINE WRITE_ELEM_STRESSES ( JSUB, NUM, IHEADER, NUM_PTS, ITABLE )

      ! Writes blocks of element stresses for one subcase and one element type for elements that do not have PCOMP properties, including
      ! all 1-D, 2-D, 3-D elements.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, BARTOR, INT_SC_NUM, MAX_NUM_STR, NDOFR, NUM_CB_DOFS,             &
                                         NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  STR_CID
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE MODEL_STUF, ONLY            :  ELEM_ONAME, ELMTYP, LABEL, SCNUM, STITLE, TITLE, TYPE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, STRE_OPT, STRE_OUT

      USE WRITE_ELEM_STRESSES_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_ELEM_STRESSES'
      CHARACTER(LEN=*), INTENT(IN)    :: IHEADER           ! Indicator of whether to write an output header

                                                           ! Array of different notes to write regarding poly fit errors
      CHARACTER( 50*BYTE)             :: ERR_INDEX_NOTE(MAX_NUM_STR)
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format
      CHARACTER(LEN=LEN(ELEM_ONAME))  :: ONAME             ! Element name to write out in F06 file
      CHARACTER( 1*BYTE)              :: WRITE_NOTES = 'N' ! Indicator of whether to write any WRT_ERR_INDEX_NOTE(i)

                                                           ! Indicators of whether to write note on indices of POLY_FIT_ERR
      CHARACTER( 1*BYTE)              :: WRT_ERR_INDEX_NOTE(MAX_NUM_STR)

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS           ! Num diff stress points for one element (3rd dim in arrays SEi, STEi)
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J,L             ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: IBEG, IEND, IELEM, ISTA, NSTA_ELEM, NROW_ELEM
      INTEGER(LONG)                   :: NCOLS             ! Num of cols to write out


      REAL(DOUBLE)                    :: ABS_ANS(11)       ! Max ABS for all element output
      REAL(DOUBLE)                    :: MAX_ANS(11)       ! Max for all element output
      REAL(DOUBLE)                    :: MIN_ANS(11)       ! Min for all element output
      REAL(DOUBLE)                    :: TINT, XI_STD, XI0, XI1
      REAL(DOUBLE)                    :: XI_RAW(11), SXC_RAW(11), SXD_RAW(11), SXE_RAW(11), SXF_RAW(11), SMAX_RAW(11),            &
                                         SMIN_RAW(11), MST_RAW(11), MSC_RAW(11)
      REAL(DOUBLE), ALLOCATABLE       :: BEAM_XI(:,:), BEAM_SXC(:,:), BEAM_SXD(:,:), BEAM_SXE(:,:), BEAM_SXF(:,:),               &
                                         BEAM_SMAX(:,:), BEAM_SMIN(:,:), BEAM_MST(:,:), BEAM_MSC(:,:)
      INTEGER(LONG), ALLOCATABLE      :: BEAM_EID(:), BEAM_GRID(:,:)

      ! op2 info
      CHARACTER( 8*BYTE)              :: TABLE_NAME             ! the name of the op2 table
      INTEGER(LONG)                   :: NNODES                 ! number of nodes for the element

      ! table -3 info
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      LOGICAL                         :: WRITE_F06, WRITE_OP2   ! flag
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: STITLEI                ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: STRESS_CODE            ! flag for type of stress; see GET_STRESS_CODE

!     op2 specific flags
      INTEGER(LONG)                   :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)                   :: NUM_WIDE         ! the number of "words" for an element
      INTEGER(LONG)                   :: NVALUES          ! the number of "words" for all the elments
      INTEGER(LONG)                   :: NTOTAL           ! the number of bytes for all NVALUES
      INTEGER(LONG)                   :: ISUBCASE         ! the subcase ID
      INTEGER(LONG)                   :: NELEMENTS
      INTEGER(LONG)                   :: ISUBCASE_INDEX   ! the index into SCNUM
      INTEGER(LONG)                   :: CID              ! coordinate system
      CHARACTER(4*BYTE)               :: CEN_WORD         ! the word "CEN/" (we need to cast the length)
      CHARACTER(139*BYTE)             :: CLINE_BUF        ! Pre-assembled CENTER line for solid stresses (matches FORMAT 1303)
      CHARACTER(139*BYTE)             :: GLINE_BUF        ! Pre-assembled GRD    line for solid stresses (matches FORMAT 1306)



! **********************************************************************************************************************************
      ! Initialize
      DEVICE_CODE = 1  ! PLOT
      STRESS_CODE = 0
 1    FORMAT("WRITE OES F06/OP2; ITABLE=",I8," (should be -4, -6, ...)")
      WRITE(ERR,1) ITABLE
      FILL(1:) = ' '

      DO I=1,MAX_NUM_STR
         WRT_ERR_INDEX_NOTE(I) = 'N'
      ENDDO

      ERR_INDEX_NOTE(1) = ' (1): Polynomial fit error is in X  normal  stress'
      ERR_INDEX_NOTE(2) = ' (2): Polynomial fit error is in Y  normal  stress'
      ERR_INDEX_NOTE(3) = ' (3): Polynomial fit error is in XY shear   stress'
      ERR_INDEX_NOTE(4) = ' (4): Polynomial fit error is in X  bending stress'
      ERR_INDEX_NOTE(5) = ' (5): Polynomial fit error is in Y  bending stress'
      ERR_INDEX_NOTE(6) = ' (6): Polynomial fit error is in XY twist   stress'
      ERR_INDEX_NOTE(7) = ' (7): Polynomial fit error is in XZ shear   stress '
      ERR_INDEX_NOTE(8) = ' (8): Polynomial fit error is in YZ shear   stress'
      ERR_INDEX_NOTE(9) = ''

      FILL(1:) = ' '

      ! Get element output name
      ONAME(1:) = ' '
      CALL GET_ELEM_ONAME ( ONAME )

      ! Write output headers if this is not the first use of this subr.
      ANALYSIS_CODE = -1
      FIELD_5_INT_FLAG = .TRUE.
      FIELD5_INT_MODE = 0
      !FIELD5_FLOAT_TIME_FREQ = 0.0
      FIELD6_EIGENVALUE = 0.0
      WRITE_F06 = (STRE_OUT(1:1) == 'Y')
      INQUIRE ( UNIT=OP2, OPENED=WRITE_OP2 )

      IF (IHEADER == 'Y') THEN
         IF (WRITE_F06) THEN
             WRITE(F06,*)
             WRITE(F06,*)
         ENDIF
         ! -- F06 header: OUTPUT FOR SUBCASE, EIGENVECTOR or CRAIG-BAMPTON DOF
         ISUBCASE_INDEX = 0
         IF    (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JSUB
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = 1  ! temp
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)
         ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
            ISUBCASE_INDEX = 1  ! statics
            ANALYSIS_CODE = 10
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
            ISUBCASE_INDEX = 1  ! statics
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            ISUBCASE_INDEX = 2  ! modes
            ANALYSIS_CODE = 7
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ISUBCASE_INDEX = 1  ! modes
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:8) == 'MFREQ') THEN
            ISUBCASE_INDEX = INT_SC_NUM
            ANALYSIS_CODE = 5
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            ISUBCASE_INDEX = 1  ! modes
            IF ((JSUB <= NDOFR) .OR. (JSUB >= NDOFR+NVEC)) THEN
               IF (JSUB <= NDOFR) THEN
                  BDY_DOF_NUM = JSUB
               ELSE
                  BDY_DOF_NUM = JSUB-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

            IF (WRITE_F06) THEN
               IF (JSUB <= NDOFR) THEN
                   WRITE(F06,103) JSUB, NUM_CB_DOFS, 'acceleration', BDY_GRID, BDY_COMP
               ELSE IF ((JSUB > NDOFR) .AND. (JSUB <= NDOFR+NVEC)) THEN
                   WRITE(F06,104) JSUB, NUM_CB_DOFS, JSUB-NDOFR
               ELSE
                   WRITE(F06,103) JSUB, NUM_CB_DOFS, 'displacement', BDY_GRID, BDY_COMP
               ENDIF
            ENDIF  ! write f06

         ENDIF
         ISUBCASE = SCNUM(ISUBCASE_INDEX)

         ! -- F06 header for TITLE, SUBTITLE, LABEL (but only to F06)
         TITLEI = TITLE(INT_SC_NUM)
         STITLEI = STITLE(INT_SC_NUM)
         LABELI = LABEL(INT_SC_NUM)

         IF (WRITE_F06) THEN
            IF (TITLE(INT_SC_NUM)(1:)  /= ' ') THEN
               WRITE(F06,201) TITLE(INT_SC_NUM)
            ENDIF

            IF (STITLE(INT_SC_NUM)(1:) /= ' ') THEN
               WRITE(F06,201) STITLE(INT_SC_NUM)
            ENDIF

            IF (LABEL(INT_SC_NUM)(1:)  /= ' ') THEN
               WRITE(F06,201) LABEL(INT_SC_NUM)
            ENDIF
            WRITE(F06,*)

           ! -- F06 1st 2 header lines for stress output description
            IF     ((TYPE(1:3) == 'BAR') .OR. (TYPE(1:4) == 'BEAM')) THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF (TYPE(1:4) == 'BUSH') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 11)
               ENDIF
               WRITE(F06,401) FILL(1: 40), ONAME

            ELSE IF (TYPE(1:4) == 'ELAS') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 11)
               ENDIF
               WRITE(F06,401) FILL(1: 40), ONAME

            ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                     IF(STR_CID == -2) THEN
                        WRITE(F06,312) FILL(1: 20)
                     ELSE
                        WRITE(F06,302) FILL(1: 15)
                     ENDIF
                  ELSE
                     IF(STR_CID == -2) THEN
                        WRITE(F06,311) FILL(1: 32)
                     ELSE
                        WRITE(F06,301) FILL(1: 27)
                     ENDIF
                  ENDIF
                  WRITE(F06,401) FILL(1: 55), ONAME
               ELSE
                  IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                     IF(STR_CID == -2) THEN
                        WRITE(F06,312) FILL(1: 27)
                     ELSE
                        WRITE(F06,302) FILL(1: 22)
                     ENDIF
                  ELSE
                     IF(STR_CID == -2) THEN
                        WRITE(F06,311) FILL(1: 38)
                     ELSE
                        WRITE(F06,301) FILL(1: 33)
                     ENDIF
                  ENDIF
                  WRITE(F06,401) FILL(1: 61), ONAME
               ENDIF

            ELSE IF (((TYPE(1:5) == 'QUAD4') ) .OR. (TYPE(1:5) == 'QUAD8')) THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 42)
               ENDIF
               WRITE(F06,401) FILL(1: 71), ONAME

            ELSE IF (TYPE(1:3) == 'ROD') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF (TYPE(1:5) == 'SHEAR') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF (TYPE(1:5) == 'TRIA3') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 36)
               ENDIF
               WRITE(F06,401) FILL(1: 65), ONAME
            ENDIF

            ! -- F06 header lines describing stress columns
            IF      (TYPE == 'BAR     ') THEN
               IF (BARTOR == 'Y') THEN
                  WRITE(F06,1101) FILL(1:1), FILL(1:1)
               ELSE
                  WRITE(F06,1102) FILL(1:1), FILL(1:1)
               ENDIF

            ELSE IF (TYPE == 'BEAM    ') THEN
               IF (BARTOR == 'Y') THEN
                  WRITE(F06,1104) FILL(1:1), FILL(1:1)
               ELSE
                  WRITE(F06,1105) FILL(1:1), FILL(1:1)
               ENDIF

            ELSE IF (TYPE(1:4) == 'ELAS') THEN
               WRITE(F06,1201) FILL(1:1), FILL(1:1)

            ELSE IF((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1301) FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1302) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF (((TYPE(1:5) == 'QUAD4') ) .OR. (TYPE(1:5) == 'QUAD8')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1401) FILL(1: 1), FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1402) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF  (TYPE == 'ROD     ') THEN
               WRITE(F06,1501) FILL(1: 1), FILL(1: 1)

            ELSE IF (TYPE(1:5) == 'SHEAR') THEN
               WRITE(F06,1601) FILL(1: 1), FILL(1: 1)

            ELSE IF (TYPE(1:5) == 'TRIA3') THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1701) FILL(1: 1), FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1702) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF  (TYPE == 'BUSH    ') THEN
               WRITE(F06,1801) FILL(1: 1), FILL(1: 1)

            ELSE IF  (TYPE == 'USERIN  ') THEN
               WRITE(F06,1901) FILL(1: 1), FILL(1: 1)

            ENDIF


         ENDIF  ! write f06






      ENDIF

      ! Write the element stress output
! --- cbeam_stations begin --- !
      IF      (TYPE == 'BAR     ') THEN
         CALL WRITE_BAR(NUM, FILL(1:1), ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                        FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_F06)

      ELSE IF (TYPE == 'BEAM    ') THEN
! --- CBEAM_standard begin --- !
         IF (WRITE_OP2) THEN
            NELEMENTS = 0
            I = 1
            DO WHILE (I <= NUM)
               NELEMENTS = NELEMENTS + 1
               J = EID_OUT_ARRAY(I,1)
               DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == J))
                  I = I + 1
               ENDDO
            ENDDO

            ALLOCATE ( BEAM_EID(NELEMENTS), BEAM_GRID(NELEMENTS,11), BEAM_XI(NELEMENTS,11), BEAM_SXC(NELEMENTS,11),               &
                       BEAM_SXD(NELEMENTS,11), BEAM_SXE(NELEMENTS,11), BEAM_SXF(NELEMENTS,11), BEAM_SMAX(NELEMENTS,11),            &
                       BEAM_SMIN(NELEMENTS,11), BEAM_MST(NELEMENTS,11), BEAM_MSC(NELEMENTS,11) )

            BEAM_GRID(:,:) = 0
            BEAM_XI(:,:)   = 0.0D0
            BEAM_SXC(:,:)  = 0.0D0
            BEAM_SXD(:,:)  = 0.0D0
            BEAM_SXE(:,:)  = 0.0D0
            BEAM_SXF(:,:)  = 0.0D0
            BEAM_SMAX(:,:) = 0.0D0
            BEAM_SMIN(:,:) = 0.0D0
            BEAM_MST(:,:)  = 0.0D0
            BEAM_MSC(:,:)  = 0.0D0

            I = 1
            IELEM = 0
            DO WHILE (I <= NUM)
               IELEM = IELEM + 1
               BEAM_EID(IELEM) = EID_OUT_ARRAY(I,1)
               IBEG = I
               DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == BEAM_EID(IELEM)))
                  I = I + 1
               ENDDO
               IEND = I - 1
               NROW_ELEM = IEND - IBEG + 1
               NSTA_ELEM = NROW_ELEM
               IF (NSTA_ELEM > 11) NSTA_ELEM = 11

! --- CBEAM_standard begin --- !
               ! Populate station grid ids for OP2 CBEAM stress records.
               ! Keep interior stations at zero and identify only the beam
               ! endpoint stations with the element end grids.
               BEAM_GRID(IELEM,:) = 0
               BEAM_GRID(IELEM,1)  = GID_OUT_ARRAY(IBEG,2)
               BEAM_GRID(IELEM,11) = GID_OUT_ARRAY(IBEG,3)
! --- CBEAM_standard end --- !

               DO ISTA=1,NSTA_ELEM
! --- CBEAM_standard begin --- !
                  ! For beam stress, EID_OUT_ARRAY/CBEAM_XL_OUT are stored once per
                  ! station, but OGEL stores two rows per station (top/bottom style
                  ! beam section stress output). Convert the station index into the
                  ! corresponding OGEL row pair explicitly.
                  K = 2*(IBEG + ISTA - 2) + 1
                  XI_RAW  (ISTA) = CBEAM_XL_OUT(IBEG + ISTA - 1)
                  SXC_RAW (ISTA) = OGEL(K    ,1)
                  SXD_RAW (ISTA) = OGEL(K    ,2)
                  SXE_RAW (ISTA) = OGEL(K    ,3)
                  SXF_RAW (ISTA) = OGEL(K    ,4)
                  SMAX_RAW(ISTA) = OGEL(K    ,6)
                  SMIN_RAW(ISTA) = OGEL(K    ,7)
                  MST_RAW (ISTA) = OGEL(K    ,8)
                  MSC_RAW (ISTA) = OGEL(K + 1,8)
                  IF (MST_RAW(ISTA) <= -0.999D0) MST_RAW(ISTA) = 0.0D0
                  IF (MSC_RAW(ISTA) <= -0.999D0) MSC_RAW(ISTA) = 0.0D0
! --- CBEAM_standard end --- !
               ENDDO

               DO ISTA=1,11
                  XI_STD = DBLE(ISTA - 1)/10.0D0
                  BEAM_XI(IELEM,ISTA) = XI_STD
                  IF (NSTA_ELEM <= 1) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (1)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (1)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (1)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (1)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(1)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(1)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (1)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (1)
                  ELSE IF (XI_STD <= XI_RAW(1)) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (1)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (1)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (1)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (1)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(1)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(1)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (1)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (1)
                  ELSE IF (XI_STD >= XI_RAW(NSTA_ELEM)) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (NSTA_ELEM)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (NSTA_ELEM)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (NSTA_ELEM)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (NSTA_ELEM)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(NSTA_ELEM)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(NSTA_ELEM)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (NSTA_ELEM)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (NSTA_ELEM)
                  ELSE
                     DO K=1,NSTA_ELEM-1
                        XI0 = XI_RAW(K)
                        XI1 = XI_RAW(K+1)
                        IF ((XI_STD >= XI0) .AND. (XI_STD <= XI1)) THEN
                           TINT = (XI_STD - XI0)/(XI1 - XI0)
                           BEAM_SXC (IELEM,ISTA) = (1.0D0 - TINT)*SXC_RAW (K) + TINT*SXC_RAW (K+1)
                           BEAM_SXD (IELEM,ISTA) = (1.0D0 - TINT)*SXD_RAW (K) + TINT*SXD_RAW (K+1)
                           BEAM_SXE (IELEM,ISTA) = (1.0D0 - TINT)*SXE_RAW (K) + TINT*SXE_RAW (K+1)
                           BEAM_SXF (IELEM,ISTA) = (1.0D0 - TINT)*SXF_RAW (K) + TINT*SXF_RAW (K+1)
                           BEAM_SMAX(IELEM,ISTA) = (1.0D0 - TINT)*SMAX_RAW(K) + TINT*SMAX_RAW(K+1)
                           BEAM_SMIN(IELEM,ISTA) = (1.0D0 - TINT)*SMIN_RAW(K) + TINT*SMIN_RAW(K+1)
                           BEAM_MST (IELEM,ISTA) = (1.0D0 - TINT)*MST_RAW (K) + TINT*MST_RAW (K+1)
                           BEAM_MSC (IELEM,ISTA) = (1.0D0 - TINT)*MSC_RAW (K) + TINT*MSC_RAW (K+1)
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO

            ELEMENT_TYPE = 2
            NUM_WIDE = 111
            NVALUES = NELEMENTS*NUM_WIDE
            CALL GET_STRESS_CODE( STRESS_CODE, 1, 0, 0 )
            CALL WRITE_OES3_STATIC_AC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
            WRITE(OP2) NVALUES
            WRITE(OP2) (BEAM_EID(IELEM)*10+DEVICE_CODE,                                                                           &
                        (BEAM_GRID(IELEM,ISTA), REAL(BEAM_XI(IELEM,ISTA),4), REAL(BEAM_SXC(IELEM,ISTA),4),                      &
                          REAL(BEAM_SXD(IELEM,ISTA),4), REAL(BEAM_SXE(IELEM,ISTA),4), REAL(BEAM_SXF(IELEM,ISTA),4),            &
                          REAL(BEAM_SMAX(IELEM,ISTA),4), REAL(BEAM_SMIN(IELEM,ISTA),4),                                          &
                          REAL(BEAM_MST(IELEM,ISTA),4), REAL(BEAM_MSC(IELEM,ISTA),4), ISTA=1,11),                               &
                        IELEM=1,NELEMENTS)
            DEALLOCATE ( BEAM_EID, BEAM_GRID, BEAM_XI, BEAM_SXC, BEAM_SXD, BEAM_SXE, BEAM_SXF, BEAM_SMAX, BEAM_SMIN,              &
                         BEAM_MST, BEAM_MSC )
         ENDIF
! --- CBEAM_standard end --- !
         CALL WRITE_CBEAM_STRESS(NUM, WRITE_F06)
! --- cbeam_stations end --- !

      ELSE IF (TYPE(1:4) == 'ELAS') THEN
         IF (WRITE_OP2) THEN
             CALL GET_SPRING_OP2_ELEMENT_TYPE(ELEMENT_TYPE)

             NUM_WIDE = 2 ! eid, spring_stress
             NVALUES = NUM_WIDE * NUM

             DEVICE_CODE = 1   ! PLOT

             !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
             CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         0)
             CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

             WRITE(OP2) NVALUES
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), I=1,NUM)
         ENDIF   ! end of op2

         IF(WRITE_F06) WRITE(F06,1103) (FILL(1:1), EID_OUT_ARRAY(I,1), OGEL(I,1),I=1,NUM)

      ELSE IF((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:4) == 'PYRA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         !       12345
         ! 39  : CTETRA
         ! 67  : CHEXA
         ! 68  : CPENTA
         ! 255 : CPYRAM
         IF (TYPE(1:4) == "HEXA") THEN
             ELEMENT_TYPE = 67
             NNODES = 9
         ELSE IF (TYPE(1:5) == "TETRA") THEN
             ELEMENT_TYPE = 39
             NNODES = 5
         ELSE IF (TYPE(1:5) == "PENTA") THEN
             ELEMENT_TYPE = 68
             NNODES = 7
         ELSE IF (TYPE(1:4) == "PYRA") THEN
             ELEMENT_TYPE = 255
             NNODES = 6
         ENDIF
         NUM_WIDE = 4 + 21*NNODES
         NELEMENTS = NUM / NUM_PTS
         NVALUES = NUM_WIDE * NELEMENTS

         IF (WRITE_OP2) THEN
           !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
           CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         0)
           CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                  TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
           WRITE(OP2) NVALUES
           CEN_WORD = "CEN/"

          ! See the CHEXA, CPENTA, or CTETRA entry for the definition of the element coordinate systems.
          ! The material coordinate system (CORDM) may be the basic system (0 or blank), any defined system
          ! (Integer > 0), or the standard internal coordinate system of the element designated as:
          ! -1: element coordinate system (-1)
          ! -2: element system based on eigenvalue techniques to insure non bias in the element formulation.

          ! TODO hardcoded
           CID = -1

          ! setting:
          !  - CTETRA: [element_device, cid, 'CEN/', 4]
          !  - CPYRAM: [element_device, cid, 'CEN/', 5]
          !  - CPENTA: [element_device, cid, 'CEN/', 6]
          !  - CHEXA:  [element_device, cid, 'CEN/', 8]

           !                 1             2             3            4            5               6             7
           !  Element    Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx      von Mises
           !     ID

           WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, CID, CEN_WORD, NNODES-1,                                       &
                       (GID_OUT_ARRAY(I,J),                                                                              &
                        REAL(OGEL(I+J-1,1),4), REAL(OGEL(I+J-1,4),4), REAL(OGEL(I+J-1,9), 4), 0., 0., 0.,              &
                        REAL(OGEL(I+J-1,12),4), REAL(OGEL(I+J-1,7),4),                                                    &
                        REAL(OGEL(I+J-1,2),4), REAL(OGEL(I+J-1,5),4), REAL(OGEL(I+J-1,10),4), 0., 0., 0.,              &
                        REAL(OGEL(I+J-1,3),4), REAL(OGEL(I+J-1,6),4), REAL(OGEL(I+J-1,11),4), 0., 0., 0.,              &
                        J=1,NNODES), I=1,NUM,NUM_PTS)
         ENDIF  ! end of op2

         IF (STRE_OPT == 'VONMISES') THEN
            NCOLS = 7
         ELSE
            NCOLS = 8
         ENDIF

         IF (WRITE_F06) THEN
            K = 0
            DO I=1,NUM,NUM_PTS
               K = K + 1
               ! Center
               WRITE(F06,1303) EID_OUT_ARRAY(I,1),(OGEL(K,J),J=1,NCOLS)
               ! Corner
               DO L=1,NUM_PTS-1
                  K = K + 1
                  WRITE(F06,1306) FILL(1: 0), GID_OUT_ARRAY(I,L+1),(OGEL(K,J),J=1,NCOLS)
               ENDDO
            ENDDO
         ENDIF

         CALL GET_MAX_MIN_ABS_STR ( NUM, NCOLS, 'N', MAX_ANS, MIN_ANS, ABS_ANS )

         IF (WRITE_F06) THEN
            IF (STRE_OPT == 'VONMISES') THEN
               WRITE(F06,1304) (MAX_ANS(J),J=1,7), (MIN_ANS(J),J=1,7), (ABS_ANS(J),J=1,7)
            ELSE
               WRITE(F06,1305) (MAX_ANS(J),J=1,8), (MIN_ANS(J),J=1,8), (ABS_ANS(J),J=1,8)
            ENDIF
         ENDIF

      ELSE IF (((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) .OR. (TYPE(1:5) == 'QUAD8')) THEN
         !CALL WRITE_OES_CQUAD4 ( NUM, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI )

         !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
         CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)

         IF (WRITE_OP2) THEN
           IF ((STRE_LOC == 'CENTER  ') .AND. (TYPE(1:5) /= 'QUAD8')) THEN
              ! CQUAD4-33
               !(eid_device,
              ! fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              ! fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out; n=17
              NUM_WIDE = 17
              ELEMENT_TYPE = 33
              NVALUES = NUM_WIDE * NUM
              CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                     TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
              !NUM_PTS = 1
              ! just a copy of the CTRIA3 code
              ! op2 version of the upper & lower layers all in one call, but without the transverse shear
              WRITE(OP2) NVALUES
              WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8), I=1,NUM)
           ELSE
              ! CQUAD4-144
               ELEMENT_TYPE = 144
              NUM_WIDE = 87 ! 2 + 17 * (4+1)  ! 4 nodes + 1 centroid

              ! TODO: probably wrong...divide NUM by NUM_PTS?
              NELEMENTS = NUM / NUM_PTS
              NVALUES = NUM_WIDE * NELEMENTS
              ! NUM=  10 NUM_PTS=   5
              !(eid_device, "CEN/", 4, # "CEN/4"
              ! fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              ! fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = n = 17+2
              !
              ! (grid,
              !  fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              !  fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)*4 = n = 17*4
              CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                     TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
              WRITE(OP2) NVALUES
              ! see the CQUAD4-33 stress/strain (the IF part of this IF-ELSE block)
              ! writing before trying to understand this...
              !
              ! basically a one-liner version of the F06 writing
              ! we broke out the L=1,NUM_PTS-1 loop to 4 lines (the GID_OUT_ARRAY lines)
              ! to avoid an additional hard to write loop
              WRITE(OP2) (EID_OUT_ARRAY(5*I+1,1)*10+DEVICE_CODE, "CEN/", 4,                                           &
                                                  (REAL(OGEL(10*I+1,J),4), J=1,8), (REAL(OGEL(10*I+2,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,2), (REAL(OGEL(10*I+3,J),4), J=1,8), (REAL(OGEL(10*I+4,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,3), (REAL(OGEL(10*I+5,J),4), J=1,8), (REAL(OGEL(10*I+6,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,4), (REAL(OGEL(10*I+7,J),4), J=1,8), (REAL(OGEL(10*I+8,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,5), (REAL(OGEL(10*I+9,J),4), J=1,8), (REAL(OGEL(10*(I+1),J),4), J=1,8), &
                          I=0,NELEMENTS-1)
           ENDIF
         ENDIF  ! end of op2

         K = 0
         DO I=1,NUM,NUM_PTS
            K = K + 1
            IF (WRITE_F06) WRITE(F06,*)
            IF (WRITE_F06) WRITE(F06,1403) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(K,J),J=1,10)
            K = K + 1
            IF (WRITE_F06) WRITE(F06,1404) FILL(1: 0), (OGEL(K,J),J=1,8)

            DO L=1,NUM_PTS-1
               K = K + 1
               IF (WRITE_F06) WRITE(F06,*)
               IF (DABS(POLY_FIT_ERR(I+L)) >= 0.01D0) THEN
                  IF (WRITE_F06) THEN
                     WRITE(F06,1405) FILL(1: 0), GID_OUT_ARRAY(I,L+1),(OGEL(K,J),J=1,10), POLY_FIT_ERR(I+L),          &
                                     POLY_FIT_ERR_INDEX(I+L)
                  ENDIF
                  WRT_ERR_INDEX_NOTE(POLY_FIT_ERR_INDEX(I+L)) = 'Y'
               ELSE
                  IF (WRITE_F06) WRITE(F06,1406) FILL(1: 0), GID_OUT_ARRAY(I,L+1),(OGEL(K,J),J=1,10), POLY_FIT_ERR(I+L)
               ENDIF

               K = K + 1
               IF (WRITE_F06) WRITE(F06,1407) FILL(1: 0), (OGEL(K,J),J=1,8)

            ENDDO
         ENDDO

         CALL GET_MAX_MIN_ABS_STR ( NUM, 10, 'Y', MAX_ANS, MIN_ANS, ABS_ANS )

             ! Get max POLY_FIT_ERR
         MAX_ANS(11) = ZERO
         K = 0
         DO I=1,NUM
            K = K + 1
            IF (POLY_FIT_ERR(I) > MAX_ANS(11)) THEN
               MAX_ANS(11) = POLY_FIT_ERR(I)
            ENDIF
            K = K + 1
         ENDDO

         MIN_ANS(11) = MAX_ANS(11)

             ! Get min POLY_FIT_ERR
         K = 0
         DO I=1,NUM
            K = K + 1
            IF (POLY_FIT_ERR(I) < MIN_ANS(11)) THEN
               MIN_ANS(11) = POLY_FIT_ERR(I)
            ENDIF
            K = K + 1
         ENDDO
             ! Get abs POLY_FIT_ERR
         ABS_ANS(11) = MAX( DABS(MAX_ANS(11)), DABS(MIN_ANS(11)) )

         IF (WRITE_F06) THEN
            IF ((STRE_LOC == 'CORNER  ') .OR. (TYPE(1:5) == 'QUAD8')) THEN
               WRITE(F06,1408) FILL(1: 0),                                                                                      &
                               FILL(1: 0),       MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),  &
                                                 MAX_ANS(10),MAX_ANS(11),                                                       &
                               FILL(1: 0),       MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),  &
                                                 MIN_ANS(10),MIN_ANS(11),                                                       &
                               FILL(1: 0),       ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),  &
                                                 ABS_ANS(10),ABS_ANS(11), FILL(1: 0)
            ELSE
               WRITE(F06,1408) FILL(1: 0),                                                                                      &
                               FILL(1: 0),       MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),  &
                                                 MAX_ANS(10),MAX_ANS(11),                                                       &
                               FILL(1: 0),       MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),  &
                                                 MIN_ANS(10),MIN_ANS(11),                                                       &
                               FILL(1: 0),       ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),  &
                                                 ABS_ANS(10),ABS_ANS(11), FILL(1: 0)
            ENDIF
         ENDIF

         WRITE_NOTES = 'N'
         DO I=1,MAX_NUM_STR
            IF (WRT_ERR_INDEX_NOTE(I) == 'Y') THEN
               WRITE_NOTES = 'Y'
            ENDIF
         ENDDO

         IF ((WRITE_NOTES == 'Y') .AND. (WRITE_F06)) THEN
            WRITE(F06,1498)
            DO I=1,MAX_NUM_STR
               IF (WRT_ERR_INDEX_NOTE(I) == 'Y') THEN
                  WRITE(F06,1499) ERR_INDEX_NOTE(I)
               ENDIF
            ENDDO
         ENDIF

      ELSE IF (TYPE == 'ROD     ') THEN
         CALL WRITE_ROD (ISUBCASE, NUM, FILL(1:1), ITABLE, TITLEI, STITLEI, LABELI,  &
                         FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_F06, WRITE_OP2 )

      ELSE IF (TYPE(1:5) == 'SHEAR') THEN
         CALL WRITE_OES_CSHEAR(NUM, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                               FIELD5_INT_MODE, FIELD6_EIGENVALUE,                   &
                               WRITE_F06, WRITE_OP2)

      ELSE IF (TYPE(1:5) == 'TRIA3') THEN
         CALL WRITE_OES_CTRIA3(NUM, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                               FIELD5_INT_MODE, FIELD6_EIGENVALUE,                   &
                               WRITE_F06, WRITE_OP2)

      ELSE IF (TYPE == 'BUSH    ') THEN
         IF (WRITE_OP2) THEN
             ELEMENT_TYPE = 102 ! CBUSH
             NUM_WIDE = 7       ! eid, tx, ty, tz, rx, ry, rz
             STRESS_CODE = 1    ! dunno
             NVALUES = NUM * NUM_WIDE

             CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
             WRITE(OP2) NVALUES
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE,(REAL(OGEL(I,J),4),J=1,6), I=1,NUM)
         ENDIF

         IF (WRITE_F06) THEN
            DO I=1,NUM
               WRITE(F06,1802) EID_OUT_ARRAY(I,1), (OGEL(I,J),J=1,6)
            ENDDO
         ENDIF

      ELSE IF (TYPE == 'USERIN  ') THEN
         IF (WRITE_F06) THEN
            DO I=1,NUM
               WRITE(F06,1902) EID_OUT_ARRAY(I,1), (OGEL(I,J),J=1,6)
            ENDDO
         ENDIF

      ELSE
         WRITE(ERR,9300) SUBR_NAME,TYPE
         WRITE(F06,9300) SUBR_NAME,TYPE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (elem type not valid) , so quit
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

  102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

  103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

  104 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

  201 FORMAT(1X,A)

  301 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M')

  302 FORMAT(A,'C B   E L E M E N T   S T R E S S E S   O T M   I N   L O C A L   E L E M E N T   C O O R D I N A T E',            &
  '   S Y S T E M')

  311 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   M A T E R I A L   C O O R D I N A T E   S Y S T E M')

  312 FORMAT(A,'C B   E L E M E N T   S T R E S S E S   O T M   I N   M A T E R I A L   C O O R D I N A T E   S Y S T E M')

  401 FORMAT(A,'F O R   E L E M E N T   T Y P E   ',A11)



! BAR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1101 FORMAT(                                                                                                                      &
          1X,A,'Element      SA1           SA2           SA3           SA4           Axial        SA-Max        SA-Min      M.S.-T'&
          ,'     Torsional'                                                                                                        &
       ,/,1X,A,'   ID        SB1           SB2           SB3           SB4          Stress        SB-Max        SB-Min      M.S.-C'&
          ,'   Stress/Margin')

 1102 FORMAT(  &
         1X,A,'Element      SA1           SA2           SA3           SA4          Axial         SA-Max        SA-Min      M.S.-T' &
      ,/,1X,A,'   ID        SB1           SB2           SB3           SB4          Stress        SB-Max        SB-Min      M.S.-C')

 1104 FORMAT(                                                                                                                      &
         1X,A,'                       S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )'                   &
      ,/,1X,A,'GRID   ELEMENT-ID       SXC           SXD           SXE           SXF         S-MAX          S-MIN         M.S.-T'      &
         ,'      M.S.-C'                                                                                                            &
      ,/,10X,'x/L')

 1105 FORMAT(                                                                                                                      &
         1X,A,'                       S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )'                   &
      ,/,1X,A,'GRID   ELEMENT-ID       SXC           SXD           SXE           SXF         S-MAX          S-MIN         M.S.-T'      &
         ,'      M.S.-C'                                                                                                            &
      ,/,10X,'x/L')

! ELAS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1201 FORMAT(1X,A,'Element     Stress     Element     Stress     Element     Stress     Element     Stress     Element     Stress' &
          ,/,1X,A,'   ID                     ID                     ID                     ID                     ID')

 1103 FORMAT(5(A,I8,1ES14.6))

! 3D Elems >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 1301 FORMAT(1X,A,'  Elem  Location            Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx    ', &
             '   von Mises'                                                                                                        &
          ,/,1X,A,'   ID')

 1302 FORMAT(1X,A,'  Elem  Location            Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx    ', &
             '      Octahedral Stress'                                                                                             &
          ,/,1X,A,'   ID',109X,'Direct        Shear')

 1303 FORMAT(1X,I8,2X,'CENTER  ',8X,8(1ES14.6))

 1304 FORMAT(28X,'------------- ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             16X,'MAX* :     ',7(ES14.6),/,                                                                                        &
             16X,'MIN* :     ',7(ES14.6),//,                                                                                       &
             16X,'ABS* :     ',7(ES14.6),/                                                                                         &
             16X,'* for output set')

 1305 FORMAT(27X,' ------------- ------------- ------------- ------------- ------------- ------------- -------------',             &
                 ' -------------',/,                                                                                               &
             16X,'MAX* :     ',8(ES14.6),/,                                                                                        &
             16X,'MIN* :     ',8(ES14.6),//,                                                                                       &
             16X,'ABS* :     ',8(ES14.6),/                                                                                         &
             16X,'* for output set')

 1306 FORMAT(1X,A,10X,'GRD',I8,5X,8(1ES14.6))

! QUAD4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1401 FORMAT(1X,A,'  Elem  Location       Fibre       Stresses In Element Coord System     Principal Stresses (Zero Shear)',       &
  '               Transverse   Transverse   % Poly',/,1X,A,                                                                        &
  '   ID                 Distance    Normal-X     Normal-Y     Shear-XY     Angle     Major        Minor      von Mises',          &
  '    Shear-XZ     Shear-YZ    Fit Err',/,1X,A,119X,'(max through thickness)')

 1402 FORMAT(1X,A,'Elem  Location         Fibre       Stresses In Element Coord System     Principal Stresses (Zero Shear)',       &
  '       Max     Transverse   Transverse   % Poly',/,1X,A,                                                                        &
  ' ID                   Distance    Normal-X     Normal-Y     Shear-XY     Angle     Major        Minor      Shear-XY',           &
  '     Shear-XZ     Shear-YZ    Fit Err',/,1X,A,119X,'(max through thickness)')

 1403 FORMAT(1X,A,I8,2X,'CENTER  ',3X,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5))

 1404 FORMAT(1X,A,21X,1ES11.3,3(1ES13.5),0PF8.2,3(1ES13.5))

 1405 FORMAT(1X,A,10X,'GRD',I8,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5),E9.1,'(',I1,')')

 1406 FORMAT(1X,A,10X,'GRD',I8,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5),E9.1)

 1407 FORMAT(1X,A,21X,1ES11.3,3(1ES13.5),0PF8.2,3(1ES13.5))

 1408 FORMAT(1X,A,32X,' ------------ ------------ ------------         ------------ ------------ ------------ ------------',       &
                 ' ------------ --------',/,                                                                                       &
             1X,A,'MAX* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,/,                                                                     &
             1X,A,'MIN* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,//,                                                                    &
             1X,A,'ABS* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,/,                                                                     &
             1X,A,'*for output set')

 1498 FORMAT(' NOTE: Explanation of errors in the polynomial fit to extrapolate element corner point stresses from values at the', &
                   ' Gauss points:')

 1499 FORMAT(6X,A)

! ROD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1501 FORMAT(  &
          1X,A,'Element     Axial       Safety     Torsional     Safety    Element    Axial       Safety     Torsional     Safety' &
       ,/,1X,A,'   ID       Stress      Margin       Stress      Margin       ID      Stress      Margin       Stress      Margin')

! SHEAR ----------------------------------------------------------------------------------------------------------------------------
 1601 FORMAT(1X,A,'Element              S t r e s s e s                           Element              S t r e s s e s'   &
          ,/,1X,A,'   ID      Normal-X      Normal-Y      Shear-XY                   ID      Normal-X      Normal-Y'               &
                 ,'      Shear-XY')


! TRIA3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1701 FORMAT(1X,A,'Element    Location      Fibre        Stresses In Element Coord System       Principal Stresses (Zero Shear)',  &
                '               Transverse   Transverse'                                                                           &
          ,/,1X,A,'   ID                   Distance     Normal-X     Normal-Y      Shear-XY     Angle     Major        Minor'      &
          ,'      von Mises     Shear-XZ     Shear-YZ'                                                                             &
          ,/,1X,A,123X,'(max through thickness)')

 1702 FORMAT(1X,A,'Element    Location      Fibre        Stresses In Element Coord System       Principal Stresses (Zero Shear)',  &
  '      Max      Transverse   Transverse'                                                                                         &
          ,/,1X,A,'   ID                   Distance     Normal-X     Normal-Y      Shear-XY     Angle     Major        Minor',     &
          '      Shear-XY     Shear-XZ     Shear-YZ',/,1X,123X,'(max through thickness)')

 1703 FORMAT(1X,I8,4X,'Anywhere',2X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1704 FORMAT(13X,'in elem',3X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1705 FORMAT(37X,'------------ ------------ ------------          ------------ ------------ ------------ ------------',            &
                 ' ------------',/,                                                                                                &
             1X,'MAX* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'MIN* : ',28x,3(ES13.5),9X,5(ES13.5),//,                                                                           &
             1X,'ABS* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'*for output set')

! BUSH >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1801 FORMAT(20X,A,'Element   Stress-1      Stress-2      Stress-3      Stress-4      Stress-5      Stress-6'                      &
          ,/,20X,A,'   ID')

 1802 FORMAT(19X,I8,8(1ES14.6))

! USERIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1901 FORMAT(20X,A,'Element   Stress-1      Stress-2      Stress-3      Stress-4      Stress-5      Stress-6'                      &
          ,/,20X,A,'   ID')

 1902 FORMAT(19X,I8,8(1ES14.6))

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 9300 FORMAT(' *ERROR  9300: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NO OUTPUT FORMAT AVAILABLE FOR ELEMENT TYPE = ',A)
! **********************************************************************************************************************************
      END SUBROUTINE WRITE_ELEM_STRESSES
!==============================================================================

      SUBROUTINE WRITE_OES_CSHEAR(NUM, FILL, ISUBCASE, ITABLE, TITLE, SUBTITLE, LABEL, &
                                  FIELD5_INT_MODE, FIELD6_EIGENVALUE,                  &
                                  WRITE_F06, WRITE_OP2)
!     TODO: calculate margin
!
      USE PENTIUM_II_KIND, ONLY     :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY :  F06, OP2, ERR
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, OGEL
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_Value, IEEE_QUIET_NAN
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32
      IMPLICIT NONE
      !
      INTEGER(LONG), INTENT(IN)       :: NUM               ! the number of elements
      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      LOGICAL, INTENT(IN)             :: WRITE_F06, WRITE_OP2
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format

      INTEGER(LONG), INTENT(INOUT) :: ITABLE          ! the current subtable number
      !LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      !REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6

      INTEGER(LONG)               :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)               :: ANALYSIS_CODE = 1 ! static fallback for OP2 table-3 header
      INTEGER(LONG)               :: NUM_WIDE = 4     ! the number of "words" for an element
      INTEGER(LONG)               :: NVALUES          ! the number of "words" for all the elments
      INTEGER(LONG)               :: NTOTAL           ! the number of bytes for all NVALUES
      INTEGER(LONG)               :: ELEMENT_TYPE = 4 ! the OP2 flag for the element
      INTEGER(LONG)               :: STRESS_CODE      ! the OP2 flag for the stress
      REAL(DOUBLE)                :: ABS_ANS(3)       ! Max ABS for output
      REAL(DOUBLE)                :: MAX_ANS(3)       ! Max for output
      REAL(DOUBLE)                :: MIN_ANS(3)       ! Min for output
      INTEGER(LONG)               :: I, J             ! DO loop indices
      REAL(REAL32)  :: NAN
      NAN = IEEE_VALUE(NAN, IEEE_QUIET_NAN)

      IF (WRITE_OP2) THEN
          DEVICE_CODE = 1   ! PLOT
          NVALUES = NUM * NUM_WIDE
          NTOTAL = NVALUES * 4

          ! the below call is a placeholder to prevent a memory bug
          CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
          ! eid, max_shear, avg_shear, margin
          CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                 TITLE, SUBTITLE, LABEL, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

 100      FORMAT("*DEBUG: WRITE_CSHEAR    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
 101      FORMAT("*DEBUG: WRITE_CSHEAR    ITABLE=",I8," (should be -5, -7,...)")
          NVALUES = NUM * NUM_WIDE
          NTOTAL = NVALUES * 4
          WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL
          WRITE(OP2) NVALUES

          ! Nastran OP2 requires this write call be a one liner...so it's a little weird...
          ! translating:
          !    DO I=1,NUM
          !        WRITE(OP2) EID_OUT_ARRAY(I,1)*10+DEVICE_CODE  ! Nastran is weird and requires scaling the ELEMENT_ID
          !
          !        convert from float64 (double precision) to float32 (single precision)
          !        RE1 = REAL(OGEL(I,1), 4)
          !        RE2 = REAL(OGEL(I,2), 4)
          !        RE3 = REAL(OGEL(I,3), 4)
          !
          !        write the max_shear, avg_shear,
          !        WRITE(OP2) RE1, RE2, RE3
          !    ENDDO
          !
          ! write the CSHEAR stress/strain data
          !Normal-X      Normal-Y      Shear-XY -> max_shear, avg_shear, margin
          WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,3), 4), REAL(OGEL(I,3), 4), &
                                                     NAN, I=1,NUM)
          WRITE(ERR,100) ITABLE
      ENDIF  ! write op2

      DO I=1,NUM,2
         IF (I+1 <= NUM) THEN
            WRITE(F06,1603) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(I,J),J=1,3), EID_OUT_ARRAY(I+1,1),(OGEL(I+1,J),J=1,3)
         ELSE
            WRITE(F06,1603) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(I,J),J=1,3)
         ENDIF
      ENDDO

      CALL GET_MAX_MIN_ABS_STR ( NUM, 3, 'N', MAX_ANS, MIN_ANS, ABS_ANS )

      WRITE(F06,1604) FILL(1: 0), FILL(1: 0), MAX_ANS(1),MAX_ANS(2),MAX_ANS(3),                                                 &
                      FILL(1: 0),             MIN_ANS(1),MIN_ANS(2),MIN_ANS(3),                                                 &
                      FILL(1: 0),             ABS_ANS(1),ABS_ANS(2),ABS_ANS(3)


 1603 FORMAT(1X,A,I8,3(1ES14.6),13X,I8,3(1ES14.6))
 1604 FORMAT(1X,A,'         ------------- ------------- ------------- ',20X,' ------------- ------------- ------------- ',/,       &
             1X,A,'MAX* : ',1X,3(ES14.6),/,                                                                                        &
             1X,A,'MIN* : ',1X,3(ES14.6),//,                                                                                       &
             1X,A,'ABS* : ',1X,3(ES14.6),/,                                                                                        &
             1X,A,'*for output set')
      END SUBROUTINE WRITE_OES_CSHEAR

!==============================================================================
      SUBROUTINE WRITE_OES_CTRIA3 ( NUM, FILL, ISUBCASE, ITABLE, TITLE, SUBTITLE, LABEL, &
                                    FIELD5_INT_MODE, FIELD6_EIGENVALUE ,                 &
                                    WRITE_F06, WRITE_OP2)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE GET_MAX_MIN_ABS_STR_Interface
      IMPLICIT NONE
      !
      INTEGER(LONG), INTENT(IN)       :: NUM               ! the number of elements
      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      LOGICAL, INTENT(IN)             :: WRITE_F06, WRITE_OP2  ! flags
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format

      INTEGER(LONG), INTENT(INOUT) :: ITABLE           ! the current subtable number
      !LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      !REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6

      INTEGER(LONG)               :: DEVICE_CODE = 1   ! PLOT, PRINT, PUNCH flag; set as PLOT
      INTEGER(LONG)               :: ANALYSIS_CODE = 1 ! static fallback for OP2 table-3 header
      INTEGER(LONG)               :: NUM_WIDE          ! the number of "words" for an element
      INTEGER(LONG)               :: NVALUES           ! the number of "words" for all the elments
      INTEGER(LONG)               :: NTOTAL            ! the number of bytes for all NVALUES
      INTEGER(LONG)               :: ELEMENT_TYPE      ! the OP2 flag for the element
      INTEGER(LONG)               :: STRESS_CODE = 1   ! the OP2 flag for the stress; preallocate
      REAL(DOUBLE)                :: ABS_ANS(11)       ! Max ABS for output
      REAL(DOUBLE)                :: MAX_ANS(11)       ! Max for output
      REAL(DOUBLE)                :: MIN_ANS(11)       ! Min for output
      INTEGER(LONG)               :: I, J, K, L, NELEMENTS, NUM_PTS_TRI           ! DO loop indices

      ! [eid, fiber_dist/curvature, oxx, oyy, txy, angle, omax, omin, ovm/max_shear,   ! upper
      !       fiber_dist/curvature, oxx, oyy, txy, angle, omax, omin, ovm/max_shear,   ! lower
      !]
      K = 0

!100  FORMAT("*DEBUG: WRITE_CTRIA3    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
!101  FORMAT("*DEBUG: WRITE_CTRIA3    ITABLE=",I8," (should be -5, -7,...)")
      IF (STRE_LOC == 'CENTER  ') THEN
         NUM_WIDE = 17
         ELEMENT_TYPE = 74
         NVALUES = NUM * NUM_WIDE
      ELSE
         NUM_PTS_TRI = 1
         NUM_WIDE = 70
         ELEMENT_TYPE = 70
         NELEMENTS = NUM
         NVALUES = NELEMENTS * NUM_WIDE
      ENDIF
      NTOTAL = NVALUES * 4
!      WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL

      IF (WRITE_OP2) THEN
          !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
          CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
          CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                 TITLE, SUBTITLE, LABEL, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
          WRITE(OP2) NVALUES

!1702     FORMAT(1X,A,'Element    Location      Fibre        Stresses In Element Coord System       Principal Stresses (Zero Shear)', &
!  '          Max      Transverse   Transverse'                                                                                       &
!              ,/,1X,A,'   ID                   Distance     Normal-X     Normal-Y      Shear-XY     Angle     Major        Minor',   &
!              '      Shear-XY     Shear-XZ     Shear-YZ',/,1X,123X,'(max through thickness)')

          IF (STRE_LOC == 'CENTER  ') THEN
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(2*I-1,J),4), J=1,8), &
                        (REAL(OGEL(2*I,J),4), J=1,8), I=1,NUM)
          ELSE
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, "CEN/", 3,                                                      &
                         (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8),                                      &
                         GID_OUT_ARRAY(I,2), (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8),                 &
                         GID_OUT_ARRAY(I,3), (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8),                 &
                         GID_OUT_ARRAY(I,4), (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8),                 &
                         I=1,NELEMENTS)
          ENDIF
      ENDIF  ! write op2

 1703 FORMAT(1X,I8,4X,'Anywhere',2X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1704 FORMAT(13X,'in elem',3X,4(1ES13.5),0PF9.3,5(1ES13.5))
 1706 FORMAT(1X,A,I8,4X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1705 FORMAT(37X,'------------ ------------ ------------          ------------ ------------ ------------ ------------',            &
                 ' ------------',/,                                                                                                &
             1X,'MAX* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'MIN* : ',28x,3(ES13.5),9X,5(ES13.5),//,                                                                           &
             1X,'ABS* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'*for output set')

      IF (STRE_LOC == 'CENTER  ') THEN
         DO I=1,NUM
            K = K + 1
            WRITE(F06,*)
            WRITE(F06,1703) EID_OUT_ARRAY(I,1),(OGEL(K,J),J=1,10)
            K = K + 1
            WRITE(F06,1704) (OGEL(K,J),J=1,8)
         ENDDO
      ELSE
         DO I=1,NUM
            K = 2*I - 1
            WRITE(F06,*)
            WRITE(F06,1703) EID_OUT_ARRAY(I,1),(OGEL(K,J),J=1,10)
            K = K + 1
            WRITE(F06,1704) (OGEL(K,J),J=1,8)
            DO L=1,3
               WRITE(F06,*)
               WRITE(F06,1706) FILL(1: 0), GID_OUT_ARRAY(I,L+1),(OGEL(2*I-1,J),J=1,10)
               WRITE(F06,1704) (OGEL(K,J),J=1,8)
            ENDDO
         ENDDO
      ENDIF

      CALL GET_MAX_MIN_ABS_STR ( NUM, 10, 'Y', MAX_ANS, MIN_ANS, ABS_ANS )

      WRITE(F06,1705) MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),MAX_ANS(10),                 &
                      MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),MIN_ANS(10),                 &
                      ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),ABS_ANS(10)

      END SUBROUTINE WRITE_OES_CTRIA3

!==============================================================================
      SUBROUTINE GET_SPRING_OP2_ELEMENT_TYPE(ELEMENT_TYPE)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY      :  ERR
      USE MODEL_STUF, ONLY  :  TYPE
      IMPLICIT NONE
      INTEGER(LONG)    :: ELEMENT_TYPE ! the OP2 flag for the element
      !               12345
      ! 11 : CELAS1 - ELAS1
      ! 12 : CELAS2
      ! 13 : CELAS3
      ! 14 : CELAS4
      IF (TYPE(5:5) == "1") THEN
          ELEMENT_TYPE = 11
      ELSE IF (TYPE(5:5) == "2") THEN
          ELEMENT_TYPE = 12
      ELSE IF (TYPE(5:5) == "3") THEN
          ELEMENT_TYPE = 13
      ELSE IF (TYPE(5:5) == "4") THEN
          ELEMENT_TYPE = 14
      ELSE
 42       FORMAT("TYPE(4:4)=",A," TYPE(5:5)=",A," TYPE(6:6)=",A)
          WRITE(ERR,42) TYPE(4:4),TYPE(5:5),TYPE(6:6)
          ELEMENT_TYPE = -1
      ENDIF
      END SUBROUTINE GET_SPRING_OP2_ELEMENT_TYPE
