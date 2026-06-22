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

      SUBROUTINE WRITE_ELEM_ENGR_FORCE ( JSUB, NUM, IHDR, NUM_PTS, ITABLE )

      ! Writes blocks of element engineering force output for one element type, one
      ! subcase. Elements that can have engineering force output are the ones
      ! enumerated below fin the IF(TYPE == ???)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, OP2
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, INT_SC_NUM, NDOFR, NUM_CB_DOFS, NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL
      USE MODEL_STUF, ONLY            :  ELEM_ONAME, LABEL, SCNUM, STITLE, TITLE, TYPE
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  FORC_LOC, FORC_OUT
      USE WRITE_ELEM_ENGR_FORCE_USE_IFs



      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_ELEM_ENGR_FORCE'
      CHARACTER(LEN=*), INTENT(IN)    :: IHDR              ! Indicator of whether to write an output header
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format
      CHARACTER(LEN=LEN(ELEM_ONAME))  :: ONAME             ! Element name to write out in F06 file

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS           ! Num diff stress points for one element
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J,J1,K,L        ! DO loop indices or counters
      INTEGER(LONG)                   :: IBEG, IEND, IELEM, ISTA, NSTA_ELEM, NELEMENTS
      INTEGER(LONG)                   :: GRID_ID
      INTEGER(LONG)                   :: NUM_TERMS         ! Number of terms to write out for shell elems

      LOGICAL                         :: WRITE_F06, WRITE_OP2   ! flag

      REAL(DOUBLE)                    :: ABS_ANS(8)       ! Max ABS for all element output
      REAL(DOUBLE)                    :: MAX_ANS(8)       ! Max for all element output
      REAL(DOUBLE)                    :: MIN_ANS(8)       ! Min for all element output
      REAL(DOUBLE)                    :: STA_XL
      REAL(DOUBLE)                    :: TINT, XI_STD, XI0, XI1
      REAL(DOUBLE)                    :: XI_RAW(11), BM1_RAW(11), BM2_RAW(11), V1_RAW(11), V2_RAW(11), AX_RAW(11), TRQ_RAW(11)
      REAL(DOUBLE), ALLOCATABLE       :: BEAM_XI(:,:), BEAM_BM1(:,:), BEAM_BM2(:,:), BEAM_V1(:,:), BEAM_V2(:,:), BEAM_AX(:,:),   &
                                         BEAM_TRQ(:,:), BEAM_WTRQ(:,:)
      INTEGER(LONG), ALLOCATABLE      :: BEAM_EID(:), BEAM_GRID(:,:)
            

      ! op2 info
      CHARACTER( 8*BYTE)              :: TABLE_NAME             ! the name of the op2 table

      ! table -3 info
      INTEGER(LONG)                   :: ISUBCASE_INDEX         ! the index into SCNUM
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: STITLEI                ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: FIELD5_INT_MODE
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE

!     op2 specific flags
      INTEGER(LONG)                   :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)                   :: NUM_WIDE     ! the number of "words" for an element
      INTEGER(LONG)                   :: NVALUES      ! the number of "words" for all the elments
      INTEGER(LONG)                   :: NTOTAL       ! the number of bytes for all NVALUES
      INTEGER(LONG)                   :: ISUBCASE     ! the subcase ID
      
      ! initialize
      ANALYSIS_CODE = -1



! **********************************************************************************************************************************
      ! initialize
      DEVICE_CODE = 1  ! PLOT
      ANALYSIS_CODE = -1
      FILL(1:) = ' '

      ! Get element output name
      ONAME(1:) = ' '
      CALL GET_ELEM_ONAME ( ONAME )

      ! Write output headers.
      ANALYSIS_CODE = -1
      FIELD5_INT_MODE = 0
      FIELD6_EIGENVALUE = 0.0

      WRITE_F06 = (FORC_OUT(1:1) == 'Y')
      INQUIRE ( UNIT=OP2, OPENED=WRITE_OP2 )


headr:IF (IHDR == 'Y') THEN

         !--- Subcase num, TITLE, SUBT, LABEL:
         !IF(WRITE_F06) WRITE(F06,*)
         !IF(WRITE_F06) WRITE(F06,*)
         CALL WRITE_SUBCASE_EIGENVEC_HEADER(JSUB, WRITE_F06)
         ISUBCASE_INDEX = 0
         IF    (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JSUB ! statics
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF(WRITE_F06) WRITE(F06,101) SCNUM(JSUB)
         ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 10
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF(WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF(WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            ISUBCASE_INDEX = 2
            ANALYSIS_CODE = 7
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF(WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF(WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:8) == 'MFREQ') THEN
            ISUBCASE_INDEX = INT_SC_NUM
            ANALYSIS_CODE = 5
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF(WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            ! Write info on what CB DOF the output is for
            ISUBCASE_INDEX = 1
            IF ((JSUB <= NDOFR) .OR. (JSUB >= NDOFR+NVEC)) THEN
               IF (JSUB <= NDOFR) THEN
                  BDY_DOF_NUM = JSUB
               ELSE
                  BDY_DOF_NUM = JSUB-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

            IF(WRITE_F06) THEN
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

         TITLEI = TITLE(INT_SC_NUM)
         STITLEI = STITLE(INT_SC_NUM)
         LABELI = LABEL(INT_SC_NUM)

         IF(WRITE_F06) THEN
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

             !--- 1st 2 lines of element specific headers - general info on what type of output:
             IF      (TYPE(1:3) == 'BAR') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:33)
                ELSE
                   WRITE(F06,301) FILL(1:39)
                ENDIF
                WRITE(F06,401) FILL(1:45), ONAME

             ELSE IF (TYPE(1:4) == 'BEAM') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:33)
                ELSE
                   WRITE(F06,301) FILL(1:39)
                ENDIF
                WRITE(F06,401) FILL(1:45), ONAME

             ELSE IF (TYPE(1:4) == 'BUSH') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:19)
                ELSE
                   WRITE(F06,301) FILL(1:24)
                ENDIF
                WRITE(F06,401) FILL(1:29), ONAME

             ELSE IF (TYPE(1:4) == 'ELAS') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:27)
                ELSE
                   WRITE(F06,301) FILL(1:33)
                ENDIF
                WRITE(F06,401) FILL(1:37), ONAME

             ELSE IF (TYPE(1:3) == 'ROD') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:27)
                ELSE
                   WRITE(F06,301) FILL(1:33)
                ENDIF
                WRITE(F06,401) FILL(1:37), ONAME

             ELSE IF (TYPE(1:5) == 'SHEAR') THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:22)
                ELSE
                   WRITE(F06,301) FILL(1:29)
                ENDIF
                WRITE(F06,401) FILL(1:32), ONAME

             ELSE IF((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8')) THEN
                IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                   WRITE(F06,302) FILL(1:33)
                ELSE
                   WRITE(F06,301) FILL(1:39)
                ENDIF
                WRITE(F06,401) FILL(1:43), ONAME
             ENDIF

             !--- Header lines describing columns of output for an element type:
             IF      (TYPE(1:3) == 'BAR'  ) THEN
                WRITE(F06,1101) FILL(1: 0), FILL(1: 0)

             ELSE IF (TYPE(1:4) == 'BEAM') THEN

             ELSE IF (TYPE(1:4) == 'ELAS') THEN
                WRITE(F06,1201) FILL(1: 0), FILL(1: 0)

             ELSE IF (TYPE(1:3) == 'ROD') THEN
                WRITE(F06,1301) FILL(1: 0), FILL(1: 0)

             ELSE IF (TYPE(1:5) == 'SHEAR') THEN
                WRITE(F06,1401) FILL(1: 0), FILL(1: 0)

             ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8')) THEN
                WRITE(F06,1501) FILL(1: 0), FILL(1: 0), FILL(1: 0)

             ELSE IF (TYPE(1:4) == 'BUSH') THEN
                WRITE(F06,1601) FILL(1: 0), FILL(1: 0)
             ENDIF

         ENDIF ! write f06

      ENDIF headr

      ! Write element force output
      IF      (TYPE == 'BAR     ') THEN

         CALL GET_MAX_MIN_ABS ( 1, 8 )

         ! (1) PRINT, (2) PLOT, (3) PUNCH, (4) NEU, (5) CSV
         IF (WRITE_OP2)  THEN  ! op2/plot
           ELEMENT_TYPE = 34
           NUM_WIDE = 9  ! eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
           NVALUES = NUM_WIDE * NUM
           CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                  TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
           WRITE(OP2) NVALUES
           WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(I,J), 4), J=1,8), I=1,NUM)
         ENDIF

         IF (WRITE_F06)  THEN  ! f06/print
           DO I=1,NUM
              WRITE(F06,1102) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(I,J),J=1,8)
           ENDDO
           !CALL GET_MAX_MIN_ABS ( 1, 8 )
           WRITE(F06,1103) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,8), FILL(1: 0), (MIN_ANS(J),J=1,8), FILL(1: 0),                 &
                                                   (ABS_ANS(J),J=1,8), FILL(1: 0)
         ENDIF
!        IF (FORC_OUT(3:3) == 'Y')  CALL WRITE_GRD_PCH_OUTPUTS(JVEC, NUM, WHAT)  ! pch/punch
!        IF (FORC_OUT(4:4) == 'Y')  CALL WRITE_GRD_NEU_OUTPUTS(JVEC, NUM, WHAT)  ! NEU
!        IF (FORC_OUT(5:5) == 'Y')  CALL WRITE_GRD_CSV_OUTPUTS(JVEC, NUM, WHAT)  ! CSV

      ELSE IF (TYPE == 'BEAM    ') THEN

         CALL GET_MAX_MIN_ABS ( 1, 8 )

! --- CBEAM_standard begin --- !
         IF (WRITE_OP2)  THEN
            NELEMENTS = 0
            I = 1
            DO WHILE (I <= NUM)
               NELEMENTS = NELEMENTS + 1
               J = EID_OUT_ARRAY(I,1)
               DO WHILE (I <= NUM)
                  IF (EID_OUT_ARRAY(I,1) /= J) EXIT
                  I = I + 1
               ENDDO
            ENDDO

            ALLOCATE ( BEAM_EID(NELEMENTS), BEAM_GRID(NELEMENTS,11), BEAM_XI(NELEMENTS,11), BEAM_BM1(NELEMENTS,11),              &
                       BEAM_BM2(NELEMENTS,11), BEAM_V1(NELEMENTS,11), BEAM_V2(NELEMENTS,11), BEAM_AX(NELEMENTS,11),               &
                       BEAM_TRQ(NELEMENTS,11), BEAM_WTRQ(NELEMENTS,11) )

            BEAM_GRID(:,:) = 0
            BEAM_XI(:,:)   = 0.0D0
            BEAM_BM1(:,:)  = 0.0D0
            BEAM_BM2(:,:)  = 0.0D0
            BEAM_V1(:,:)   = 0.0D0
            BEAM_V2(:,:)   = 0.0D0
            BEAM_AX(:,:)   = 0.0D0
            BEAM_TRQ(:,:)  = 0.0D0
            BEAM_WTRQ(:,:) = 0.0D0

            I = 1
            IELEM = 0
            DO WHILE (I <= NUM)
               IELEM = IELEM + 1
               BEAM_EID(IELEM) = EID_OUT_ARRAY(I,1)
               IBEG = I
               DO WHILE (I <= NUM)
                  IF (EID_OUT_ARRAY(I,1) /= BEAM_EID(IELEM)) EXIT
                  I = I + 1
               ENDDO
               IEND = I - 1
               NSTA_ELEM = IEND - IBEG + 1
               IF (NSTA_ELEM > 11) NSTA_ELEM = 11

! --- CBEAM_standard begin --- !
               ! Populate station grid ids for OP2 CBEAM records.
               ! Use the end grids at x/L=0 and x/L=1 and leave interior
               ! stations as zero, which matches classic beam-station
               ! conventions better than an all-zero grid list.
               BEAM_GRID(IELEM,:) = 0
               BEAM_GRID(IELEM,1)  = GID_OUT_ARRAY(IBEG,2)
               BEAM_GRID(IELEM,11) = GID_OUT_ARRAY(IBEG,3)
! --- CBEAM_standard end --- !

               DO ISTA=1,NSTA_ELEM
                  XI_RAW (ISTA) = CBEAM_XL_OUT(IBEG + ISTA - 1)
                  BM1_RAW(ISTA) = OGEL(IBEG + ISTA - 1,1)
                  BM2_RAW(ISTA) = OGEL(IBEG + ISTA - 1,2)
                  V1_RAW (ISTA) = OGEL(IBEG + ISTA - 1,5)
                  V2_RAW (ISTA) = OGEL(IBEG + ISTA - 1,6)
                  AX_RAW (ISTA) = OGEL(IBEG + ISTA - 1,7)
                  TRQ_RAW(ISTA) = OGEL(IBEG + ISTA - 1,8)
               ENDDO

               DO ISTA=1,11
                  XI_STD = DBLE(ISTA - 1)/10.0D0
                  BEAM_XI(IELEM,ISTA) = XI_STD
                  IF (NSTA_ELEM <= 1) THEN
                     BEAM_BM1(IELEM,ISTA) = BM1_RAW(1)
                     BEAM_BM2(IELEM,ISTA) = BM2_RAW(1)
                     BEAM_V1 (IELEM,ISTA) = V1_RAW(1)
                     BEAM_V2 (IELEM,ISTA) = V2_RAW(1)
                     BEAM_AX (IELEM,ISTA) = AX_RAW(1)
                     BEAM_TRQ(IELEM,ISTA) = TRQ_RAW(1)
                  ELSE IF (XI_STD <= XI_RAW(1)) THEN
                     BEAM_BM1(IELEM,ISTA) = BM1_RAW(1)
                     BEAM_BM2(IELEM,ISTA) = BM2_RAW(1)
                     BEAM_V1 (IELEM,ISTA) = V1_RAW(1)
                     BEAM_V2 (IELEM,ISTA) = V2_RAW(1)
                     BEAM_AX (IELEM,ISTA) = AX_RAW(1)
                     BEAM_TRQ(IELEM,ISTA) = TRQ_RAW(1)
                  ELSE IF (XI_STD >= XI_RAW(NSTA_ELEM)) THEN
                     BEAM_BM1(IELEM,ISTA) = BM1_RAW(NSTA_ELEM)
                     BEAM_BM2(IELEM,ISTA) = BM2_RAW(NSTA_ELEM)
                     BEAM_V1 (IELEM,ISTA) = V1_RAW(NSTA_ELEM)
                     BEAM_V2 (IELEM,ISTA) = V2_RAW(NSTA_ELEM)
                     BEAM_AX (IELEM,ISTA) = AX_RAW(NSTA_ELEM)
                     BEAM_TRQ(IELEM,ISTA) = TRQ_RAW(NSTA_ELEM)
                  ELSE
                     DO K=1,NSTA_ELEM-1
                        XI0 = XI_RAW(K)
                        XI1 = XI_RAW(K+1)
                        IF ((XI_STD >= XI0) .AND. (XI_STD <= XI1)) THEN
                           TINT = (XI_STD - XI0)/(XI1 - XI0)
                           BEAM_BM1(IELEM,ISTA) = (1.0D0 - TINT)*BM1_RAW(K) + TINT*BM1_RAW(K+1)
                           BEAM_BM2(IELEM,ISTA) = (1.0D0 - TINT)*BM2_RAW(K) + TINT*BM2_RAW(K+1)
                           BEAM_V1 (IELEM,ISTA) = (1.0D0 - TINT)*V1_RAW (K) + TINT*V1_RAW (K+1)
                           BEAM_V2 (IELEM,ISTA) = (1.0D0 - TINT)*V2_RAW (K) + TINT*V2_RAW (K+1)
                           BEAM_AX (IELEM,ISTA) = (1.0D0 - TINT)*AX_RAW (K) + TINT*AX_RAW (K+1)
                           BEAM_TRQ(IELEM,ISTA) = (1.0D0 - TINT)*TRQ_RAW(K) + TINT*TRQ_RAW(K+1)
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO

            ELEMENT_TYPE = 2
            NUM_WIDE = 100
            NVALUES = NELEMENTS*NUM_WIDE
            CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                   TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
            WRITE(OP2) NVALUES
            WRITE(OP2) (BEAM_EID(IELEM)*10+DEVICE_CODE,                                                                             &
                         (BEAM_GRID(IELEM,ISTA), REAL(BEAM_XI(IELEM,ISTA),4), REAL(BEAM_BM1(IELEM,ISTA),4),                        &
                          REAL(BEAM_BM2(IELEM,ISTA),4), REAL(BEAM_V1(IELEM,ISTA),4), REAL(BEAM_V2(IELEM,ISTA),4),                 &
                          REAL(BEAM_AX(IELEM,ISTA),4), REAL(BEAM_TRQ(IELEM,ISTA),4), REAL(BEAM_WTRQ(IELEM,ISTA),4), ISTA=1,11),  &
                         IELEM=1,NELEMENTS)
            DEALLOCATE ( BEAM_EID, BEAM_GRID, BEAM_XI, BEAM_BM1, BEAM_BM2, BEAM_V1, BEAM_V2, BEAM_AX, BEAM_TRQ, BEAM_WTRQ )
         ENDIF
! --- CBEAM_standard end --- !

         IF (WRITE_F06)  THEN
            WRITE(F06,'(A,/,A,/,A)') '                          F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )', &
                                      '         ELEMENT-ID        - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL       TOTAL        WARPING', &
                                      '    GRID   STAT X/L       PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE       TORQUE       TORQUE'
            I = 1
            DO WHILE (I <= NUM)
               IBEG = I
               IELEM = EID_OUT_ARRAY(I,1)
               DO WHILE (I <= NUM)
                  IF (EID_OUT_ARRAY(I,1) /= IELEM) EXIT
                  I = I + 1
               ENDDO
               IEND = I - 1
               NSTA_ELEM = IEND - IBEG + 1

               WRITE(F06,*)
               WRITE(F06,1114) 0, IELEM
               DO ISTA=1,NSTA_ELEM
                  GRID_ID = 0
                  IF (ISTA == 1) GRID_ID = GID_OUT_ARRAY(IBEG,2)
                  IF (ISTA == NSTA_ELEM) GRID_ID = GID_OUT_ARRAY(IBEG,3)
                  STA_XL = CBEAM_XL_OUT(IBEG + ISTA - 1)
                  WRITE(F06,1115) FILL(1: 0), GRID_ID, STA_XL,                                                    &
                                  OGEL(IBEG + ISTA - 1,1), OGEL(IBEG + ISTA - 1,2),                               &
                                  OGEL(IBEG + ISTA - 1,5), OGEL(IBEG + ISTA - 1,6),                               &
                                  OGEL(IBEG + ISTA - 1,7), OGEL(IBEG + ISTA - 1,8), ZERO
               ENDDO
            ENDDO
            WRITE(F06,1116) FILL(1: 0), MAX_ANS(1), MAX_ANS(2), MAX_ANS(5), MAX_ANS(6), MAX_ANS(7), MAX_ANS(8), ZERO
            WRITE(F06,1117) FILL(1: 0), MIN_ANS(1), MIN_ANS(2), MIN_ANS(5), MIN_ANS(6), MIN_ANS(7), MIN_ANS(8), ZERO
            WRITE(F06,1118) FILL(1: 0), ABS_ANS(1), ABS_ANS(2), ABS_ANS(5), ABS_ANS(6), ABS_ANS(7), ABS_ANS(8), ZERO
         ENDIF

      ELSE IF (TYPE(1:4) == 'ELAS') THEN
           ! Engr force for ELAS was put into OGEL(I,1)

!          CALL GET_SPRING_OP2_ELEMENT_TYPE(ELEMENT_TYPE)
!          NUM_WIDE = 2 ! eid, spring_force
!          NVALUES = NUM_WIDE * NUM
!          CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
!                                 TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!          WRITE(OP2) NVALUES
!          WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), I=1,NUM)
!
!          ! TODO: what's going on with this loop having the 1,NUM,5??? and the J=J1,J1+4???
! =======
!
!xx      WRITE(F06,1202) FILL(1: 0), (EID_OUT_ARRAY(I,1),OGEL(I,1),I=1,NUM)

         IF (WRITE_F06)  THEN  ! f06/print
           J1 = 1
           DO I=1,NUM,5
              IF (J1+4 <= NUM) THEN
                 WRITE(F06,1202) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,1), J=J1,J1+4)
                 J1 = J1 + 5
              ELSE
                 WRITE(F06,1202) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,1), J=J1,NUM)
              ENDIF
           ENDDO
           CALL GET_MAX_MIN_ABS ( 1, 1 )
           WRITE(F06,1203) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,1), FILL(1: 0), (MIN_ANS(J),J=1,1), FILL(1: 0),                 &
                                                   (ABS_ANS(J),J=1,1), FILL(1: 0)
         ENDIF

      ELSE IF (TYPE == 'ROD     ') THEN
         IF (WRITE_OP2)  THEN  ! op2/plot
           !CALL WRITE_OEF_ROD ( ISUBCASE, NUM, FILL(1:1), FILL(1:16), ITABLE, TITLEI, STITLEI, LABELI )
           ELEMENT_TYPE = 1
           NUM_WIDE = 3  ! eid, axial, torsion
           NVALUES = NUM_WIDE * NUM
           CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                  TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

           ! TODO: why does fields 7/8 write out the axial and torsion?
           WRITE(OP2) NVALUES
           WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,7), 4), REAL(OGEL(I,8), 4), I=1,NUM)
         ENDIF

         IF (WRITE_F06)  THEN  ! f06/print
           J1 = 1
           DO I=1,NUM,3
              IF (J1+2 <= NUM) THEN
                 WRITE(F06,1302) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,7), OGEL(J,8), J=J1,J1+2)
                 J1 = J1 + 3
              ELSE
                 WRITE(F06,1302) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,7), OGEL(J,8), J=J1,NUM)
              ENDIF
           ENDDO
           CALL GET_MAX_MIN_ABS ( 7, 8 )
           WRITE(F06,1303) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=7,8), FILL(1: 0), (MIN_ANS(J),J=7,8), FILL(1: 0),  &
                                                   (ABS_ANS(J),J=7,8), FILL(1: 0)
         ENDIF

      ELSE IF (TYPE == 'SHEAR   ') THEN
         IF (WRITE_OP2)  THEN  ! op2/plot
           !CALL WRITE_SHEAR_OEF()
           ELEMENT_TYPE = 4  ! CSHEAR
           ! eid,[
           !  force41, force21, force12, force32, force23, force43,
           !  force34, force14,
           !  kick_force1, shear12, kick_force2, shear23,
           !  kick_force3, shear34, kick_force4, shear41,
           !
           NUM_WIDE = 17
           !CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
           !                       TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
           NVALUES = NUM * NUM_WIDE
           !WRITE(OP2) NVALUES
           ! write the CSHEAR force data
           !WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,3), 4), REAL(OGEL(I,3), 4), &
           !                                               NAN, I=1,NUM)
         ENDIF

         IF (WRITE_F06)  THEN  ! f06/print
           J1 = 1
           DO I=1,NUM,2
              IF      (J1+1 <= NUM) THEN
                 WRITE(F06,1402) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,1), OGEL(J,2), OGEL(J,3), J=J1,J1+1)
                 J1 = J1 + 2
              ELSE
                 WRITE(F06,1402) FILL(1: 0), (EID_OUT_ARRAY(J,1), OGEL(J,1), OGEL(J,2), OGEL(J,3), J=J1,NUM)
              ENDIF
           ENDDO
           CALL GET_MAX_MIN_ABS ( 1, 3 )
           WRITE(F06,1403) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,3), FILL(1: 0), (MIN_ANS(J),J=1,3), FILL(1: 0),                 &
                                                   (ABS_ANS(J),J=1,3), FILL(1: 0)
         ENDIF

      ELSE IF ((TYPE == 'TRIA3K  ') .OR. (TYPE == 'QUAD4K  ')) THEN
         IF (WRITE_F06) THEN
             DO I=1,NUM
                WRITE(F06,1512) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(I,J),J=1,6)
             ENDDO
             CALL GET_MAX_MIN_ABS ( 1, 8 )
             WRITE(F06,1513) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,6), FILL(1: 0), (MIN_ANS(J),J=1,6), FILL(1: 0),  &
                                                     (ABS_ANS(J),J=1,6), FILL(1: 0)
         ENDIF
         NUM_TERMS = 6

      ELSE IF ((TYPE == 'TRIA3   ') .OR. ((TYPE == 'QUAD4   ') ) .OR. (TYPE == 'QUAD8   ')) THEN
        IF (WRITE_OP2)  THEN
          IF ((FORC_LOC == 'CENTER  ') .AND. (TYPE /= 'QUAD8   ')) THEN
             IF (TYPE == 'TRIA3   ') THEN
                 ELEMENT_TYPE = 74
             ELSE
                 ELEMENT_TYPE = 33
             ENDIF
             NUM_WIDE = 9
             NVALUES = NUM * NUM_WIDE
             CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
             WRITE(OP2) NVALUES
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(I,J),4),J=1,8), I=1,NUM)
          ELSE
             IF (TYPE == 'TRIA3   ') THEN
                ELEMENT_TYPE = 70
                NUM_WIDE = 38
                NELEMENTS = NUM
                NVALUES = NUM_WIDE * NELEMENTS
                CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                       TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
                WRITE(OP2) NVALUES
                WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, "CEN/",                                                    &
                            0, (REAL(OGEL(I,J),4),J=1,8),                                                                  &
                            GID_OUT_ARRAY(I,2), (REAL(OGEL(I,J),4),J=1,8),                                                 &
                            GID_OUT_ARRAY(I,3), (REAL(OGEL(I,J),4),J=1,8),                                                 &
                            GID_OUT_ARRAY(I,4), (REAL(OGEL(I,J),4),J=1,8),                                                 &
                            I=1,NELEMENTS)
             ELSE
                NELEMENTS = NUM / NUM_PTS
                IF (TYPE == 'QUAD8   ') THEN
                   ELEMENT_TYPE = 64
                ELSE
                   ELEMENT_TYPE = 144
                ENDIF
                NUM_WIDE = 47
                NVALUES = NUM_WIDE * NELEMENTS
                CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                       TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
                WRITE(OP2) NVALUES
                IF (TYPE == 'QUAD8   ') THEN
                   WRITE(OP2) (EID_OUT_ARRAY(5*I+1,1)*10+DEVICE_CODE, "CEN/",                                           &
                               0,                     (REAL(OGEL(5*I+1,J),4),J=1,8),                                      &
                               GID_OUT_ARRAY(5*I+1,2),(REAL(OGEL(5*I+2,J),4),J=1,8),                                      &
                               GID_OUT_ARRAY(5*I+1,3),(REAL(OGEL(5*I+3,J),4),J=1,8),                                      &
                               GID_OUT_ARRAY(5*I+1,4),(REAL(OGEL(5*I+4,J),4),J=1,8),                                      &
                               GID_OUT_ARRAY(5*I+1,5),(REAL(OGEL(5*I+5,J),4),J=1,8),                                      &
                               I=0,NELEMENTS-1)
                ELSE
                   WRITE(OP2) (EID_OUT_ARRAY(4*I+1,1)*10+DEVICE_CODE, "CEN/",                                           &
                               0, (REAL((OGEL(4*I+1,J)+OGEL(4*I+2,J)+OGEL(4*I+3,J)+OGEL(4*I+4,J))/4.0D0,4),J=1,8),     &
                               GID_OUT_ARRAY(4*I+1,2), (REAL(OGEL(4*I+1,J),4),J=1,8),                                     &
                               GID_OUT_ARRAY(4*I+1,3), (REAL(OGEL(4*I+2,J),4),J=1,8),                                     &
                               GID_OUT_ARRAY(4*I+1,4), (REAL(OGEL(4*I+3,J),4),J=1,8),                                     &
                               GID_OUT_ARRAY(4*I+1,5), (REAL(OGEL(4*I+4,J),4),J=1,8),                                     &
                               I=0,NELEMENTS-1)
                ENDIF
             ENDIF
          ENDIF
        ENDIF

        IF (WRITE_F06)  THEN  ! f06
          K = 0
          IF (TYPE == 'TRIA3   ') THEN
             DO I=1,NUM
                WRITE(F06,1524) FILL(1: 0), EID_OUT_ARRAY(I,1), '        ', (OGEL(I,J),J=1,8)
                DO L=2,4
                   WRITE(F06,1525) FILL(1: 0), GID_OUT_ARRAY(I,L),(OGEL(I,J),J=1,8)
                ENDDO
             ENDDO
          ELSE
             DO I=1,NUM,NUM_PTS
                K = K + 1
                                                           ! Center forces
                IF(TYPE == 'QUAD8   ') THEN
                  WRITE(F06,1524) FILL(1: 0), EID_OUT_ARRAY(I,1), 'CENTER  ', (OGEL(K,J),J=1,8)
                ELSE
                  WRITE(F06,1524) FILL(1: 0), EID_OUT_ARRAY(I,1), '        ', (OGEL(K,J),J=1,8)
                ENDIF

                DO L=2,NUM_PTS                                ! Corner forces
                  K = K + 1
                  WRITE(F06,1525) FILL(1: 0), GID_OUT_ARRAY(I,L),(OGEL(K,J),J=1,8)
                ENDDO
             ENDDO
          ENDIF
          CALL GET_MAX_MIN_ABS ( 1, 8 )
          WRITE(F06,1523) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,8), FILL(1: 0), (MIN_ANS(J),J=1,8), FILL(1: 0),  &
                                                  (ABS_ANS(J),J=1,8), FILL(1: 0)
        ENDIF
        NUM_TERMS = 8

      ELSE IF (TYPE(1:4) == 'BUSH') THEN
         ! Engr force for BUSH was put into OGEL(I,1-6)
         IF (WRITE_OP2)  THEN  ! op2/plot
           ELEMENT_TYPE = 102 ! CBUSH
           NUM_WIDE = 7       ! eid, tx, ty, tz, rx, ry, rz
           NVALUES = NUM * NUM_WIDE
           CALL WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, &
                                  TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
           WRITE(OP2) NVALUES
           WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE,(REAL(OGEL(I,J),4),J=1,6), I=1,NUM)
         ENDIF

         IF (WRITE_F06)  THEN  ! f06/print
           DO I=1,NUM
              WRITE(F06,1602) FILL(1: 0), EID_OUT_ARRAY(I,1),(OGEL(I,J),J=1,6)
           ENDDO
           CALL GET_MAX_MIN_ABS ( 1, 6 )
           WRITE(F06,1603) FILL(1: 0), FILL(1: 0), (MAX_ANS(J),J=1,6), FILL(1: 0), (MIN_ANS(J),J=1,6), FILL(1: 0),  &
                                                   (ABS_ANS(J),J=1,6), FILL(1: 0)
         ENDIF

      ENDIF


      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

  102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

  103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

  104 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

  201 FORMAT(1X,A)

  301 FORMAT(16X,A,'E L E M E N T   E N G I N E E R I N G   F O R C E S')

  302 FORMAT(16X,A,'C B   E L E M E N T   E N G I N E E R I N G   F O R C E   O T M')

  401 FORMAT(16X,A,'F O R   E L E M E N T   T Y P E   ',A11)

! BAR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1101 FORMAT(16X,A,' Element       Bend-Moment End A           Bend-Moment End B              - Shear -              Axial'        &
          ,'         Torque'  &
          ,/,16X,A,'    ID       Plane 1       Plane 2       Plane 1       Plane 2      Plane 1       Plane 2        Force')

 1102 FORMAT(16X,A,I8,8(1ES14.6))

! --- cbeam_stations begin --- !
 1114 FORMAT(I1,8X,I8)

 1115 FORMAT(1X,A,I8,2X,F7.3,1X,7(1ES14.6))

 1116 FORMAT(11X,A,'MAX* :  ',7(ES14.6))

 1117 FORMAT(11X,A,'MIN* :  ',7(ES14.6))

 1118 FORMAT(11X,A,'ABS* :  ',7(ES14.6))
! --- cbeam_stations end --- !

 1103 FORMAT(1X,A,'         ------------- ------------- ------------- ------------- ------------- ------------- -------------',    &
                        ' -------------',/,                                                                                        &
             16X,A,'MAX* :  ',8(ES14.6),/,                                                                                         &
             16X,A,'MIN* :  ',8(ES14.6),//,                                                                                        &
             16X,A,'ABS* :  ',8(ES14.6),/,                                                                                         &
             16X,A,'*for output set')

! ELAS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1201 FORMAT(16X,A,' Element     Force      Element     Force      Element     Force      Element     Force      Element     Force'&
          ,/,16X,A,'    ID                     ID                     ID                     ID                     ID')

 1202 FORMAT(16X,A,5(I8,1ES14.6))

 1203 FORMAT(16X,A,'         -------------',/,                                                                                     &
             16X,A,'MAX* :  ',1(ES14.6),/,                                                                                         &
             16X,A,'MIN* :  ',1(ES14.6),//,                                                                                        &
             16X,A,'ABS* :  ',1(ES14.6),/,                                                                                         &
             16X,A,'*for output set')

! ROD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1301 FORMAT(16X,A,' Element     Axial        Torque      Element     Axial        Torque      Element     Axial        Torque'    &
          ,/,16X,A,'    ID       Force                       ID       Force                       ID       Force')

 1302 FORMAT(16X,A,3(I8,1ES14.6,1ES14.6))

 1303 FORMAT(16X,A,'         ------------- -------------',/,                                                                       &
             16X,A,'MAX* :  ',2(1ES14.6),/,                                                                                        &
             16X,A,'MIN* :  ',2(1ES14.6),//,                                                                                       &
             16X,A,'ABS* :  ',2(1ES14.6),/,                                                                                        &
             16X,A,'*for output set')

! SHEAR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1401 FORMAT(16X,A,' Element        N o r m a l   F o r c e s           Element        N o r m a l   F o r c e s          '        &
          ,/,16X,A,'    ID        Nxx           Nyy           Nxy          ID        Nxx           Nyy           Nxy')

 1402 FORMAT(1X,A,2(I8,3(1ES14.6),1X))

 1403 FORMAT(16X,A,'         ------------- ------------- -------------',/,                                                         &
             16X,A,'MAX* :  ',3ES14.6,/,                                                                                           &
             16X,A,'MIN* :  ',3ES14.6,//,                                                                                          &
             16X,A,'ABS* :  ',3ES14.6,/,                                                                                           &
             16X,A,'*for output set')

! SHELL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1501 FORMAT(1X,A,' Element Location               N o r m a l   F o r c e s                       M o m e n t s'                  &
          ,19X,'T r a n s v e r s e',/1X,A,'    ID', 104X,'S h e a r   F o r c e s'                                                &
          ,/,16X,A,'              Nxx           Nyy           Nxy           Mxx           Myy           Mxy            Qx         '&
          ,'  Qy')

!            WRITE(F06,1501) FILL(1: 0), FILL(1: 0), FILL(1: 0)
 1512 FORMAT(1X,A,I8,15X,6(1ES14.6))

 1513 FORMAT(16X,A,'          ------------- ------------- ------------- ------------- ------------- -------------',/,              &
             16X,A,'MAX* :  ',6(ES14.6),/,                                                                                         &
             16X,A,'MIN* :  ',6(ES14.6),//,                                                                                        &
             16X,A,'ABS* :  ',6(ES14.6),/,                                                                                         &
             16X,A,'*for output set')

 1523 FORMAT(16X,A,'          ------------- ------------- ------------- ------------- ------------- ------------- -------------',  &
                        ' -------------',/,                                                                                        &
             16X,A,'MAX* :  ',8(ES14.6),/,                                                                                         &
             16X,A,'MIN* :  ',8(ES14.6),//,                                                                                        &
             16X,A,'ABS* :  ',8(ES14.6),/,                                                                                         &
             16X,A,'*for output set')

 1524 FORMAT(1X,A,I8,2X,A,5X,8(1ES14.6))

 1525 FORMAT(1X,A,10X,'GRD',I8,2X,8(1ES14.6))


! BUSH >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1601 FORMAT(16X,A,' Element      Force         Force         Force        Moment        Moment        Moment'                     &
          ,/,16X,A,'    ID         XE            YE            ZE            XE            YE            ZE')

 1602 FORMAT(16X,A,I8,6(1ES14.6))

 1603 FORMAT(16X,A,'          ------------- ------------- ------------- ------------- ------------- ------------- ',/,             &
             16X,A,'MAX* :  ',6(ES14.6),/,                                                                                         &
             16X,A,'MIN* :  ',6(ES14.6),//,                                                                                        &
             16X,A,'ABS* :  ',6(ES14.6),/,                                                                                         &
             16X,A,'*for output set')

! **********************************************************************************************************************************

      CONTAINS


! ##################################################################################################################################

      SUBROUTINE GET_MAX_MIN_ABS ( BEG_COL, END_COL )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: BEG_COL           ! Col number in OGEL where to beg for averaging to get max, min, abs
      INTEGER(LONG), INTENT(IN)       :: END_COL           ! Col number in OGEL where to end for averaging to get max, min, abs
      INTEGER(LONG)                   :: II,JJ             ! DO loop indices or counters

! **********************************************************************************************************************************
      ! Get MAX, MIN, ABS values
      DO JJ=BEG_COL,END_COL
         MAX_ANS(JJ) = -MACH_LARGE_NUM
      ENDDO

      DO II=1,NUM
         DO JJ=BEG_COL,END_COL
            IF (OGEL(II,JJ) > MAX_ANS(JJ)) THEN
               MAX_ANS(JJ) = OGEL(II,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO JJ=BEG_COL,END_COL
         MIN_ANS(JJ) = MAX_ANS(JJ)
      ENDDO

      DO II=1,NUM
         DO JJ=BEG_COL,END_COL
            IF (OGEL(II,JJ) < MIN_ANS(JJ)) THEN
               MIN_ANS(JJ) = OGEL(II,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO II=BEG_COL,END_COL
         ABS_ANS(II) = MAX( DABS(MAX_ANS(II)), DABS(MIN_ANS(II)) )
      ENDDO

      END SUBROUTINE GET_MAX_MIN_ABS

      END SUBROUTINE WRITE_ELEM_ENGR_FORCE

