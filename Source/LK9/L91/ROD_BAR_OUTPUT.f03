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

   MODULE ROD_BAR_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_ROD, WRITE_BAR, ROD_MARGIN, BAR_MARGIN

   CONTAINS

      SUBROUTINE WRITE_ROD ( ISUBCASE, NUM, FILL_F06, ITABLE,           &
                             TITLE, SUBTITLE, LABEL,                    &
                             FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_F06, WRITE_OP2)

! Routine for writing output to text files F06 for ROD element stresses. Up to 2 elements written per line of output.
! Data is first written to character variables and then that character variable is output the F06.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE CONSTANTS_1, ONLY           :  ZERO
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MSPRNT, OGEL
      USE OP2_STRESS_OUTPUT, ONLY     :  WRITE_OES3_STATIC
      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      CHARACTER(LEN=*), INTENT(IN)    :: FILL_F06          ! Padding for output format
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      INTEGER(LONG), INTENT(IN)       :: FIELD5_INT_MODE
      REAL(DOUBLE),  INTENT(IN)       :: FIELD6_EIGENVALUE
      LOGICAL, INTENT(IN)             :: WRITE_OP2         ! write the op2
      LOGICAL, INTENT(IN)             :: WRITE_F06         ! write the f06

      CHARACTER(  1*BYTE)             :: MSFLAG            ! If margin is negative, MSFLAG is an *
      CHARACTER(118*BYTE)             :: RLINE_F06         ! Result of concatenating char. variables below to make a line of
!                                                            stress output for 1 or 2 CROD's

      CHARACTER(  8*BYTE)             :: REID1             ! Internal file: element ID of 1st CROD
      CHARACTER( 14*BYTE)             :: RSTR11            ! Internal file: axial stress in 1st CROD
      CHARACTER( 10*BYTE)             :: RMS11             ! Internal file: M.S. for axial stress in 1st CROD
      CHARACTER(  1*BYTE)             :: RMSF11            ! Internal file: MSFLAG for axial M.S. in 1st CROD
      CHARACTER( 14*BYTE)             :: RSTR12            ! Internal file: torsional stress in 1st CROD
      CHARACTER( 10*BYTE)             :: RMS12             ! Internal file: M.S. for torsional stress in 1st CROD
      CHARACTER(  1*BYTE)             :: RMSF12            ! Internal file: MSFLAG for torsional M.S. in 1st CROD
      CHARACTER( 14*BYTE)             :: RMS31             ! Internal file: M.S. for axial stress in 1st CROD
      CHARACTER( 14*BYTE)             :: RMS32             ! Internal file: M.S. for axial stress in 1st CROD
      CHARACTER( 14*BYTE)             :: RMS41             ! Internal file: M.S. for axial stress in 1st CROD
      CHARACTER( 14*BYTE)             :: RMS42             ! Internal file: M.S. for axial stress in 1st CROD

      CHARACTER(  8*BYTE)             :: REID2             ! Internal file: element ID of 2nd CROD
      CHARACTER( 14*BYTE)             :: RSTR21            ! Internal file: axial stress in 2nd CROD
      CHARACTER( 10*BYTE)             :: RMS21             ! Internal file: M.S. for axial stress in 2nd CROD
      CHARACTER(  1*BYTE)             :: RMSF21            ! Internal file: MSFLAG for axial M.S. in 2nd CROD
      CHARACTER( 14*BYTE)             :: RSTR22            ! Internal file: torsional stress in 2nd CROD
      CHARACTER( 10*BYTE)             :: RMS22             ! Internal file: M.S. for torsional stress in 2nd CROD
      CHARACTER(  1*BYTE)             :: RMSF22            ! Internal file: MSFLAG for torsional M.S. in 2nd CROD

      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: ABS_ANS(4)        ! Max ABS for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MAX_ANS(4)        ! Max for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MIN_ANS(4)        ! Min for all grids output for each of the 6 disp components

      INTEGER(LONG)                   :: ELEM_TYPE         ! should be 1=CROD, 3=CTUBE, 10=CONROD


! **********************************************************************************************************************************
      ELEM_TYPE = 1  ! CROD
      !ELEM_TYPE = 3  ! CTUBE
      !ELEM_TYPE = 10  ! CONROD
      CALL OUTPUT2_WRITE_OES_ROD(ISUBCASE, ELEM_TYPE, NUM, ITABLE, TITLE, SUBTITLE, LABEL, &
                                 FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_OP2)
      DO I=1,NUM,2

         RLINE_F06(1:)  = ' '

         REID1(1:)  = ' '
         RSTR11(1:) = ' '
         RMS11(1:)  = ' '
         RMSF11(1:) = ' '
         RSTR12(1:) = ' '
         RMS12(1:)  = ' '
         RMSF12(1:) = ' '

         REID2(1:)  = ' '
         RSTR21(1:) = ' '
         RMS21(1:)  = ' '
         RMSF21(1:) = ' '
         RSTR22(1:) = ' '
         RMS22(1:)  = ' '
         RMSF22(1:) = ' '

         RMS31(1:)  = ' '
         RMS32(1:)  = ' '
         RMS41(1:)  = ' '
         RMS42(1:)  = ' '

         WRITE(REID1,2201) EID_OUT_ARRAY(I,1)

! Write axial stress output to a temporary internal file for one element

         WRITE(RSTR11,2202) OGEL(I,1)
         MSFLAG = ' '
         IF (MSPRNT(I,1) /= '0') THEN
            IF (OGEL(I,2) < ZERO) THEN
               MSFLAG = '*'
            ENDIF
            WRITE(RMS11,2203) OGEL(I,2)
            WRITE(RMS31,2213) OGEL(I,2)
            WRITE(RMSF11,2204) MSFLAG
         ENDIF

! Write torsional stress output to a temporary internal file for one element

         WRITE(RSTR12,2202) OGEL(I,3)
         MSFLAG = ' '
         IF (MSPRNT(I,2) /= '0') THEN
            IF (OGEL(I,4) < ZERO) THEN
               MSFLAG = '*'
            ENDIF
            WRITE(RMS12,2203) OGEL(I,4)
            WRITE(RMS32,2213) OGEL(I,4)
            WRITE(RMSF12,2204) MSFLAG
         ENDIF

! Write axial stress output to a temporary internal file for another element

         IF ((I+1) <= NUM) THEN

            WRITE(REID2,2201) EID_OUT_ARRAY(I+1,1)

            WRITE(RSTR21,2202) OGEL(I+1,1)
            MSFLAG = ' '
            IF (MSPRNT(I+1,1) /= '0') THEN
               IF (OGEL(I+1,2) < ZERO) THEN
                  MSFLAG = '*'
               ENDIF
               WRITE(RMS21,2203) OGEL(I+1,2)
               WRITE(RMS41,2213) OGEL(I+1,2)
               WRITE(RMSF21,2204) MSFLAG
            ENDIF

! Write torsional stress output to a temporary internal file for another element

            WRITE(RSTR22,2202) OGEL(I+1,3)
            MSFLAG = ' '
            IF (MSPRNT(I+1,2) /= '0') THEN
               IF (OGEL(I+1,4) < ZERO) THEN
                  MSFLAG = '*'
               ENDIF
               WRITE(RMS22,2203) OGEL(I+1,4)
               WRITE(RMS42,2213) OGEL(I+1,4)
               WRITE(RMSF22,2204) MSFLAG
            ENDIF

         ENDIF

         IF (WRITE_F06) THEN
                                                           ! Write a line of output, consisting of stress output for one
                                                           ! or two elements, to the output file
            RLINE_F06 = REID1//RSTR11//RMS11//RMSF11//RSTR12//RMS12//RMSF12//REID2//RSTR21//RMS21//RMSF21//RSTR22//RMS22//RMSF22
            WRITE(F06,2205) FILL_F06, RLINE_F06
         ENDIF

      ENDDO

      IF (WRITE_F06) THEN
         CALL GET_MAX_MIN_ABS ( 1, 4 )
         WRITE(F06,9104) (MAX_ANS(J),J=1,4),(MIN_ANS(J),J=1,4),(ABS_ANS(J),J=1,4)
      ENDIF


      RETURN

! **********************************************************************************************************************************
 2201 FORMAT(I8)

 2202 FORMAT(ES14.6)

 2203 FORMAT(ES10.2)

 2213 FORMAT(4X,ES10.2)

 2204 FORMAT(A,A)

 2205 FORMAT(A,A)

9104 FORMAT( 1X,'         ------------- ---------  ------------- ---------',/,                                                     &
             1X,'MAX* :  ',ES14.6,ES10.2,1X,ES14.6,ES10.2/,                                                                        &
             1X,'MIN* :  ',ES14.6,ES10.2,1X,ES14.6,ES10.2//,                                                                       &
             1X,'ABS* :  ',ES14.6,ES10.2,1X,ES14.6,ES10.2,/,                                                                       &
             1X,'*for output set')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_MAX_MIN_ABS ( BEG_COL, END_COL )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: BEG_COL           ! Col number in OGEL where to beg for averaging to get max, min, abs
      INTEGER(LONG), INTENT(IN)       :: END_COL           ! Col number in OGEL where to end for averaging to get max, min, abs
      INTEGER(LONG)                   :: II,JJ,KK          ! DO loop indices or counters

! **********************************************************************************************************************************
! Get MAX, MIN, ABS values

      DO JJ=1,END_COL-BEG_COL+1
         MAX_ANS(JJ) = -MACH_LARGE_NUM
      ENDDO

      DO II=1,NUM
         KK = 0
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(II,JJ) > MAX_ANS(JJ)) THEN
               MAX_ANS(KK) = OGEL(II,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO JJ=1,END_COL-BEG_COL+1
         MIN_ANS(JJ) = MAX_ANS(JJ)
      ENDDO

      DO II=1,NUM
         KK = 0
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(II,JJ) < MIN_ANS(JJ)) THEN
               MIN_ANS(KK) = OGEL(II,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO II=1,END_COL-BEG_COL+1
         ABS_ANS(II) = MAX( DABS(MAX_ANS(II)), DABS(MIN_ANS(II)) )
      ENDDO

      END SUBROUTINE GET_MAX_MIN_ABS

      END SUBROUTINE WRITE_ROD
!==================================================================================================
      SUBROUTINE OUTPUT2_WRITE_OES_ROD(ISUBCASE, ELEM_TYPE, NUM, ITABLE, TITLE, SUBTITLE, LABEL, &
                                       FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_OP2)
!     writes the CROD/CTUBE/CONROD stress/strain results.
!     Data is first written to character variables and then that character variable is output the F06.
!
!     Parameters
!     ==========
!     ELEM_TYPE : int
!       flag for the element type
!       - 1 : CROD
!       - 3 : CTUBE
!       - 10 : CONROD
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  OP2, ERR
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, OGEL
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_Value, IEEE_QUIET_NAN
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32
      USE OP2_STRESS_OUTPUT, ONLY     :  WRITE_OES3_STATIC
      INTEGER(LONG) :: ANALYSIS_CODE        ! static, time, frequency, modal, etc. flag
      INTEGER(LONG), INTENT(IN) :: ISUBCASE  ! subcase id
      INTEGER(LONG), INTENT(IN) :: NUM       ! the number of elements in OGEL to write
      INTEGER(LONG), INTENT(IN) :: ELEM_TYPE
      LOGICAL,       INTENT(IN) :: WRITE_OP2               ! is this a PLOT result -> OP2
      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL
      INTEGER(LONG), INTENT(IN) :: FIELD5_INT_MODE
      REAL(DOUBLE),  INTENT(IN) :: FIELD6_EIGENVALUE

      INTEGER(LONG) :: ITABLE       ! the current subtable number
      INTEGER(LONG) :: NUM_WIDE     ! the number of "words" for a single element
      INTEGER(LONG) :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG) :: STRESS_CODE  ! flag for what the output result means (e.g., fiber_distance/curvature, max_shear, von_mises)
      INTEGER(LONG) :: NVALUES      ! the number of "words" for the OP2 data
      INTEGER(LONG) :: NTOTAL       ! the number of bytes corresponding to nvalues
      INTEGER(LONG) :: I
      REAL(REAL32)  :: NAN
      LOGICAL       :: IS_PRINT     ! is this a PRINT result -> F06

      NAN = IEEE_VALUE(NAN, IEEE_QUIET_NAN)
      ! TODO: assuming PLOT
      DEVICE_CODE = 1

! **********************************************************************************************************************************
!******OP2
      ! we're already setup for subtable -3
!==================================================================================================

      IF (WRITE_OP2) THEN
        ! eid, axisl_stress, axial_margin, torsional stress, torsional_margin
        NUM_WIDE = 5

        ! dunno???
        STRESS_CODE = 1
        CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE, &
                               TITLE, SUBTITLE, LABEL,                                          &
                               FIELD5_INT_MODE, FIELD6_EIGENVALUE)
        ! ITABLE = -4, -6, ...
        !NWORDS = NUM * NUM_WIDE
        !NTOTAL = NBYTES_PER_WORD * NWORDS

 100    FORMAT("*DEBUG:    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
        NVALUES = NUM * NUM_WIDE
        NTOTAL = NVALUES * 4
        WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL
        WRITE(OP2) NVALUES

 101    FORMAT ("*OP2:     EID_DEVICE=",I8," STRESS=",ES10.2,ES10.2)
        DO I=1,NUM
          WRITE(ERR,101) EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), REAL(OGEL(I,2), 4)
        ENDDO
      ! Nastran OP2 requires this write call be a one liner...so it's a little weird...
      ! translating:
      !    DO I=1,NUM
      !        WRITE(OP2) EID_OUT_ARRAY(I,1)*10+DEVICE_CODE  ! Nastran is weird and requires scaling the ELEMENT_ID
      !
      !        convert from float64 (double precision) to float32 (single precision)
      !        RE1 = REAL(OGEL(I,1), 4)
      !        RE2 = REAL(OGEL(I,2), 4)
      !        RE3 = REAL(OGEL(I,3), 4)
      !        RE4 = REAL(OGEL(I,4), 4)
      !
      !        write the axisl_stress, axial_margin, torsional stress, torsional_margin
      !        WRITE(OP2) RE1, RE2, RE3, RE4
      !    ENDDO
      !
!      TODO: figure out what's going on with the shear stress
!      TODO: calculate axial/torision margin
        ! write the rod stress/strain data
        !WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), REAL(OGEL(I,2), 4), &
        !                                               REAL(OGEL(I,3), 4), REAL(OGEL(I,4), 4), I=1,NUM)
        WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), NAN, &
                                                       REAL(OGEL(I,2), 4), NAN, I=1,NUM)
 102    FORMAT("*DEBUG: WRITE_CROD      ITABLE=",I8," (should be -5, -7,...)")
        WRITE(ERR,102) ITABLE
      ENDIF

!      IF (IS_PRINT) THEN
!        F06
!      ENDIF
      END SUBROUTINE OUTPUT2_WRITE_OES_ROD


      SUBROUTINE WRITE_BAR (NUM, FILL_F06, ISUBCASE, ITABLE,            &
                            TITLE, SUBTITLE, LABEL,                     &
                            FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_F06 )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2
      USE SCONTR, ONLY                :  BARTOR, BLNK_SUB_NAM, MOGEL
      USE CONSTANTS_1, ONLY           :  ZERO
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MAXREQ, MSPRNT, OGEL

      USE RESULT_FORMATTING, ONLY     :  WRT_REAL_TO_CHAR_VAR

      USE OP2_STRESS_OUTPUT, ONLY     :  WRITE_OES3_STATIC
      USE OP2_STRESS_OUTPUT, ONLY     :  GET_STRESS_CODE, WRITE_OES3_STATIC
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_BAR'

      CHARACTER(LEN=*), INTENT(IN)    :: FILL_F06          ! Padding for output format
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      INTEGER(LONG), INTENT(IN)       :: FIELD5_INT_MODE
      REAL(DOUBLE),  INTENT(IN)       :: FIELD6_EIGENVALUE
      LOGICAL,       INTENT(IN)       :: WRITE_F06

      CHARACTER(133*BYTE)             :: BLINE1A           ! Result of concatenating char. variables BOUT1, BMS1, BMSF1, BTOR to
!                                                            make the 1st line of stress output for a CBAR with torsional stress
      CHARACTER(133*BYTE)             :: BLINE1B           ! Result of concatenating char. variables BOUT1, BMS1, BMSF1, BMS2, BMSF2
!                                                            to make the 1st line of stress output for a CBAR w/o torsional stress
      CHARACTER(133*BYTE)             :: BLINE2A           ! Result of concatenating char. variables BOUT2, BMS2, BMSF2
!                                                            to make the 2nd line of stress output for a CBAR
      CHARACTER(133*BYTE)             :: BLINE2B           ! Result of concatenating char. variables BOUT2, BMS2, BMSF2
!                                                            to make the 2nd line of stress output for a CBAR
      CHARACTER(107*BYTE)             :: BOUT1             ! Internal file: element ID and elem stresses for end A of element
      CHARACTER(107*BYTE)             :: BOUT2             ! Internal file: element ID and elem stresses for end B of element
      CHARACTER( 10*BYTE)             :: BMS1              ! Internal file: margin of safety in tension for element
      CHARACTER( 10*BYTE)             :: BMS2              ! Internal file: margin of safety in tension for element
      CHARACTER( 14*BYTE)             :: BMS3              ! Internal file: margin of safety in tension for element
      CHARACTER(  1*BYTE)             :: BMSF1             ! Internal file: margin of safety flag (MSFLAG) for tension M.S.
      CHARACTER(  1*BYTE)             :: BMSF2             ! Internal file: margin of safety flag (MSFLAG) for compr M.S.
      CHARACTER(  1*BYTE)             :: BMSF3             ! Internal file: margin of safety flag (MSFLAG) for torque M.S.
      CHARACTER( 14*BYTE)             :: BTOR              ! Internal file: torsional stress

      CHARACTER(14*BYTE)              :: ABS_ANS_CHAR(16)  ! Character variable that contains the 6 grid abs  outputs
      CHARACTER(14*BYTE)              :: MAX_ANS_CHAR(16)  ! Character variable that contains the 6 grid max  outputs
      CHARACTER(14*BYTE)              :: MIN_ANS_CHAR(16)  ! Character variable that contains the 6 grid min  outputs
      CHARACTER(  1*BYTE)             :: MSFLAG            ! If margin is negative, MSFLAG is an *
      CHARACTER(14*BYTE)              :: OGEL_CHAR(MOGEL)  ! Char representation of 1 row of OGEL outputs

      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! The subcase ID
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE)                    :: ABS_ANS(16)       ! Max ABS for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MAX_ANS(16)       ! Max for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MIN_ANS(16)       ! Min for all grids output for each of the 6 disp components
      ! op2
      INTEGER(LONG)                   :: ELEMENT_TYPE = 34 ! CBAR-34; constant
      INTEGER(LONG)                   :: NUM_WIDE = 16     ! number of fields for the element; constant
      INTEGER(LONG)                   :: DEVICE_CODE = 1   ! PLOT; constant
      INTEGER(LONG)                   :: STRESS_CODE
      INTEGER(LONG)                   :: NVALUES           ! number of values in the op2 block
      STRESS_CODE = 1



! **********************************************************************************************************************************
      ! Routine for writing output to text file F06 for BAR element stresses. Up to 2 elements written per line of output.
      ! Data is first written to character variables and then that character variable is output the F06.
      ! op2_headers = ['s1a', 's2a', 's3a', 's4a', 'axial', 'smaxa', 'smina', 'MS_tension',
      !                's1b', 's2b', 's3b', 's4b',          'smaxb', 'sminb', 'MS_compression']

      NVALUES = NUM_WIDE * NUM

      !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
      CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         0)
      CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                             TITLE, SUBTITLE, LABEL, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

      WRITE(OP2) NVALUES
      !K = 2*I+1
      !
      ! TODO: this is likely not right, but it should be the right shape...
      WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE,   &
                  (REAL(OGEL(2*I-1,J), 4), J=1,8),     & ! s1a-smina; K = 2*I+1
                  (REAL(OGEL(2*I,  J), 4), J=1,7),     & ! s1b-sminb; K = 2*I
                  I=1,NUM)
! ******************************
      K = 0
      DO I=1,NUM
         BLINE1A(1:) = ' '
         BLINE1B(1:) = ' '
         BLINE2A(1:) = ' '
         BLINE2B(1:) = ' '

         BOUT1(1:)   = ' '
         BMS1(1:)    = ' '
         BMSF1       = ' '

         BOUT2(1:)   = ' '
         BMS2(1:)    = ' '
         BMSF2       = ' '

         BTOR(1:)    = ' '

         ! Write first line of output for one element to a temporary internal file
         K = K + 1
         CALL WRT_REAL_TO_CHAR_VAR ( OGEL, MAXREQ, MOGEL, K, OGEL_CHAR )
         WRITE(BOUT1,9011) EID_OUT_ARRAY(I,1), (OGEL_CHAR(J),J=1,7)

         MSFLAG = ' '
         IF (MSPRNT(K,1) /= '0') THEN
            IF (OGEL(K,8) < ZERO) THEN
               MSFLAG = '*'
            ENDIF
            WRITE(BMS1, 9021) OGEL(K,8)
            WRITE(BMSF1,9031) MSFLAG
         ELSE
            WRITE(BMS1, 9022)
            WRITE(BMSF1,9032)
         ENDIF
         IF (BARTOR == 'Y') THEN
            WRITE(BTOR, 9041) OGEL(K,9)
         ENDIF

        ! Write second line of output for one element to a temporary internal file
         K = K + 1
         CALL WRT_REAL_TO_CHAR_VAR ( OGEL, MAXREQ, MOGEL, K, OGEL_CHAR )
         WRITE(BOUT2,9012) (OGEL_CHAR(J),J=1,4), (OGEL_CHAR(J),J=6,7)

         MSFLAG = ' '
         IF (MSPRNT(K,2) /= '0') THEN
            IF (OGEL(K,8) < ZERO) THEN
               MSFLAG = '*'
            ENDIF
            WRITE(BMS2, 9021) OGEL(K,8)
            WRITE(BMSF2,9031) MSFLAG
         ELSE
            WRITE(BMS2, 9022)
            WRITE(BMSF2,9032)
         ENDIF

         MSFLAG = ' '
         IF (MSPRNT(K,3) /= '0') THEN
            IF (OGEL(K,9) < ZERO) THEN
               MSFLAG = '*'
            ENDIF
            WRITE(BMS3, 9023) OGEL(K,9)
            WRITE(BMSF3,9031) MSFLAG
         ELSE
            WRITE(BMS3, 9024)
            WRITE(BMSF3,9032)
         ENDIF

         ! Write the two lines of stress output for one element to F06
         IF (WRITE_F06) THEN
            WRITE(F06,*)
            IF (BARTOR == 'Y') THEN
               BLINE1A = BOUT1//BMS1//BMSF1//BTOR
               BLINE2A = BOUT2//BMS2//BMSF2//BMS3//BMSF3
               WRITE(F06,9031) BLINE1A
               WRITE(F06,9031) BLINE2A
            ELSE
               BLINE1B = BOUT1//BMS1//BMSF1
               BLINE2B = BOUT2//BMS2//BMSF2
               WRITE(F06,9031) BLINE1B
               WRITE(F06,9031) BLINE2B
            ENDIF
         ENDIF

      ENDDO

      IF (WRITE_F06) THEN
         CALL GET_MAX_MIN_ABS ( 1, 8 )
         WRITE(F06,9108) (MAX_ANS_CHAR(J),J=1,7), MAX_ANS(8), (MAX_ANS_CHAR(J),J=9,15), MAX_ANS(16),                               &
                         (MIN_ANS_CHAR(J),J=1,7), MIN_ANS(8), (MIN_ANS_CHAR(J),J=9,15), MIN_ANS(16),                               &
                         (ABS_ANS_CHAR(J),J=1,7), ABS_ANS(8), (ABS_ANS_CHAR(J),J=9,15), ABS_ANS(16)
      ENDIF


      RETURN

! **********************************************************************************************************************************
 9011 FORMAT(1X,I8,7A)

 9012 FORMAT(1X,8X,4A,14X,2A)

 9021 FORMAT(ES10.2)

 9022 FORMAT('          ')

 9023 FORMAT(4X,ES10.2)

 9024 FORMAT('              ')

 9031 FORMAT(A)

 9032 FORMAT(' ')

 9041 FORMAT(ES14.6)

 9108 FORMAT( 1X,'         ------------- ------------- ------------- ------------- ------------- ------------- -------------',     &
                        ' ---------',/,                                                                                            &
             1X,'MAX* :  ',7A,1ES10.2,/,                                                                                           &
             1X,'MAX* :  ',7A,1ES10.2,//,                                                                                          &
             1X,'MIN* :  ',7A,1ES10.2,/,                                                                                           &
             1X,'MIN* :  ',7A,1ES10.2,//,                                                                                          &
             1X,'ABS* :  ',7A,1ES10.2,/,                                                                                           &
             1X,'ABS* :  ',7A,1ES10.2,/,                                                                                           &
             1X,'*for output set')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_MAX_MIN_ABS ( BEG_COL, END_COL )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: BEG_COL           ! Col number in OGEL where to beg for averaging to get max, min, abs
      INTEGER(LONG), INTENT(IN)       :: END_COL           ! Col number in OGEL where to end for averaging to get max, min, abs
      INTEGER(LONG)                   :: II,JJ,KK,LL       ! DO loop indices or counters

! **********************************************************************************************************************************
! Get MAX, MIN, ABS values

      DO JJ=1,2*(END_COL-BEG_COL+1)
         MAX_ANS(JJ) = -MACH_LARGE_NUM
      ENDDO

      LL = 0
      DO II=1,NUM
         KK = 0
         LL = LL + 1
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(LL,JJ) > MAX_ANS(KK)) THEN
               MAX_ANS(KK) = OGEL(LL,JJ)
            ENDIF
         ENDDO
         LL = LL + 1
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(2*II,JJ) > MAX_ANS(KK)) THEN
               MAX_ANS(KK) = OGEL(LL,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO JJ=1,2*(END_COL-BEG_COL+1)
         MIN_ANS(JJ) = MAX_ANS(JJ)
      ENDDO

      LL = 0
      DO II=1,NUM
         KK = 0
         LL = LL + 1
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(LL,JJ) < MIN_ANS(KK)) THEN
               MIN_ANS(KK) = OGEL(LL,JJ)
            ENDIF
         ENDDO
         LL = LL + 1
         DO JJ=BEG_COL,END_COL
            KK = KK + 1
            IF (OGEL(LL,JJ) < MIN_ANS(KK)) THEN
               MIN_ANS(KK) = OGEL(LL,JJ)
            ENDIF
         ENDDO
      ENDDO

      DO II=1,2*(END_COL-BEG_COL+1)
         ABS_ANS(II) = MAX( DABS(MAX_ANS(II)), DABS(MIN_ANS(II)) )
      ENDDO

      DO II=1,16

         IF (ABS_ANS(II) == 0.0) THEN
            WRITE(ABS_ANS_CHAR(II),'(A)') '  0.0         '
         ELSE
            WRITE(ABS_ANS_CHAR(II),'(1ES14.6)') ABS_ANS(II)
         ENDIF

         IF (MAX_ANS(II) == 0.0) THEN
            WRITE(MAX_ANS_CHAR(II),'(A)') '  0.0         '
         ELSE
            WRITE(MAX_ANS_CHAR(II),'(1ES14.6)') MAX_ANS(II)
         ENDIF

         IF (MIN_ANS(II) == 0.0) THEN
            WRITE(MIN_ANS_CHAR(II),'(A)') '  0.0         '
         ELSE
            WRITE(MIN_ANS_CHAR(II),'(1ES14.6)') MIN_ANS(II)
         ENDIF

      ENDDO

      END SUBROUTINE GET_MAX_MIN_ABS

      END SUBROUTINE WRITE_BAR


      SUBROUTINE ROD_MARGIN (ICOL, S1, S2, MS1, MS2, MSP1, MSP2 )

! Calculates margins of safety for ROD element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, ONEPM6, ONEPP10
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  ULT_STRE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ROD_MARGIN'
      CHARACTER(LEN=*), INTENT(OUT)   :: MSP1,MSP2         ! If '1',  print margins in F06 file. If '0', do not print.

      INTEGER(LONG), INTENT(IN)       :: ICOL              ! Column no. from ULT_STRE to get max allow. stresses


      REAL(DOUBLE), INTENT(OUT)       :: MS1               ! Calculated margin of safety
      REAL(DOUBLE), INTENT(OUT)       :: MS2               ! Calculated margin of safety
      REAL(DOUBLE), INTENT(IN)        :: S1                ! An input stress for which margins are calculated
      REAL(DOUBLE), INTENT(IN)        :: S2                ! An input stress for which margins are calculated
      REAL(DOUBLE)                    :: EPS5              ! Small value
      REAL(DOUBLE)                    :: S1M,S2M           ! DABS(S1), etc
      REAL(DOUBLE)                    :: SCP               !-DABS(SC)
      REAL(DOUBLE)                    :: ST,SC,SS          ! Material allowables in tension, compr., shear
      REAL(DOUBLE)                    :: STM,SCM,SSM       ! DABS(ST), etc

      INTRINSIC                       :: DABS, DMIN1



! **********************************************************************************************************************************
      EPS5 = EPSIL(5)

! Initialize outputs

      MSP1 = ' '
      MSP2 = ' '
      MS1  = ZERO
      MS2  = ZERO

! Calc outputs

      ST   =  ULT_STRE(1,ICOL)
      SC   =  ULT_STRE(2,ICOL)
      SS   =  ULT_STRE(3,ICOL)
      STM  =  DABS(ST)
      SCM  =  DABS(SC)
      SSM  =  DABS(SS)
      SCP  = -SCM

! Calculate margins for ROD1 element. S1 is axial stress, S2 is torsional stress (S3, S4 not used)
! MS1 is the margin for axial stress and MS2 is the margin for torsion.

      S1M = DABS(S1)

      IF (S1 >= ZERO) THEN                                 ! ROD margin for positive axial stress

         IF (S1M < EPS5) THEN
            MS1  = ONEPP10
            MSP1 = '0'
         ELSE
            MS1  = (STM/S1) - ONE
            MSP1 = '1'
         ENDIF
         IF (STM < EPS5) THEN
            MSP1 = '0'
         ENDIF

      ELSE                                                 ! ROD margin for negative axial stress

         IF (S1M < EPS5) THEN
            MS1  = ONEPP10
            MSP1 = '0'
         ELSE
            MS1  = (SCP/S1) - ONE
            MSP1 = '1'
         ENDIF
         IF (SCM < EPS5) THEN
            MSP1 = '0'
         ENDIF

      ENDIF

      S2M = DABS(S2)                                       ! ROD torsional margin (positive or negative torsional stress)

      IF (S2M < EPS5) THEN
         MS2  = ONEPP10
         MSP2 = '0'
      ELSE
         MS2  = (SSM/S2M) - ONE
         MSP2 = '1'
      ENDIF
      IF (SSM < EPS5) THEN
         MSP2 = '0'
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ROD_MARGIN


      SUBROUTINE BAR_MARGIN ( ICOL, S1, S2, S3, S4, S5, MS1, MS2, MS3, MSP1, MSP2, MSP3 )

! Calculates margins of safety for BAR element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, ONEPM6, ONEPM15, ONEPP10
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  ULT_STRE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BAR_MARGIN'
      CHARACTER(LEN=*), INTENT(OUT)   :: MSP1              ! If '1',  print margins in F06 file. If '0', do not print.
      CHARACTER(LEN=*), INTENT(OUT)   :: MSP2              ! If '1',  print margins in F06 file. If '0', do not print.
      CHARACTER(LEN=*), INTENT(OUT)   :: MSP3              ! If '1',  print margins in F06 file. If '0', do not print.

      INTEGER(LONG), INTENT(IN)       :: ICOL              ! Column no. from ULT_STRE to get max allow. stresses


      REAL(DOUBLE), INTENT(OUT)       :: MS1               ! Calculated margin of safety
      REAL(DOUBLE), INTENT(OUT)       :: MS2               ! Calculated margin of safety
      REAL(DOUBLE), INTENT(OUT)       :: MS3               ! Calculated margin of safety
      REAL(DOUBLE), INTENT(IN)        :: S1                ! An input stress for which margins are calculated
      REAL(DOUBLE), INTENT(IN)        :: S2                ! An input stress for which margins are calculated
      REAL(DOUBLE), INTENT(IN)        :: S3                ! An input stress for which margins are calculated
      REAL(DOUBLE), INTENT(IN)        :: S4                ! An input stress for which margins are calculated
      REAL(DOUBLE), INTENT(IN)        :: S5                ! An input stress for which margins are calculated
      REAL(DOUBLE)                    :: EPS5              ! Small value
      REAL(DOUBLE)                    :: EPS6              ! Small value
      REAL(DOUBLE)                    :: S1U               ! = S1 unless S1 is zero and then it is (+/-)EPS6
      REAL(DOUBLE)                    :: S2U               ! = S2 unless S2 is zero and then it is (+/-)EPS6
      REAL(DOUBLE)                    :: S3U               ! = S3 unless S3 is zero and then it is (+/-)EPS6
      REAL(DOUBLE)                    :: S4U               ! = S4 unless S4 is zero and then it is (+/-)EPS6
      REAL(DOUBLE)                    :: S1M               ! DABS(S1)
      REAL(DOUBLE)                    :: S2M               ! DABS(S2)
      REAL(DOUBLE)                    :: S3M               ! DABS(S3)
      REAL(DOUBLE)                    :: S4M               ! DABS(S4)
      REAL(DOUBLE)                    :: S5M               ! DABS(S5)
      REAL(DOUBLE)                    :: SCP               !-DABS(SC)
      REAL(DOUBLE)                    :: ST,SC,SS          ! Material allowables in tension, compr., shear
      REAL(DOUBLE)                    :: STM,SCM,SSM       ! DABS(ST), etc

      INTRINSIC                       :: DABS, DMIN1



! **********************************************************************************************************************************
      EPS5 = EPSIL(5)
      EPS6 = EPSIL(6)

! Initialize outputs

      MSP1 = ' '
      MSP2 = ' '
      MSP3 = ' '
      MS1  = ZERO
      MS2  = ZERO
      MS3  = ZERO

! Calc outputs

      ST   =  ULT_STRE(1,ICOL)
      SC   =  ULT_STRE(2,ICOL)
      SS   =  ULT_STRE(3,ICOL)
      STM  =  DABS(ST)
      SCM  =  DABS(SC)
      SSM  =  DABS(SS)
      SCP  = -SCM

! Calculate margins for BAR1 element. S1 is max stress at end A, S2 is max stress at end B. S3 is min stress at
! end A, S4 is min stress at end B.

      S1M = DABS(S1)
      S2M = DABS(S2)
      S3M = DABS(S3)
      S4M = DABS(S4)

      S1U = S1
      IF (S1M <= EPS5) THEN
         IF (S1 >= ZERO) THEN
            S1U =  EPS6
         ELSE
            S1U = -EPS6
         ENDIF
      ENDIF

      S2U = S2
      IF (S2M <= EPS5) THEN
         IF (S2 >= ZERO) THEN
            S2U =  EPS6
         ELSE
            S2U = -EPS6
         ENDIF
      ENDIF

      S3U = S3
      IF (S3M <= EPS5) THEN
         IF (S3 >= ZERO) THEN
            S3U =  EPS6
         ELSE
            S3U = -EPS6
         ENDIF
      ENDIF

      S4U = S4
      IF (S4M <= EPS5) THEN
         IF (S4 >= ZERO) THEN
            S4U =  EPS6
         ELSE
            S4U = -EPS6
         ENDIF
      ENDIF

      IF ((S1M < EPS5) .AND. (S2M < EPS5)) THEN            ! Margin for tension
         MS1  = ONEPP10   !(i.e. 1.E+10)
         MSP1 = '0'
      ELSE IF ((S1 > ZERO) .AND. (S2 <= ZERO)) THEN
         MS1  = STM/S1U - ONE
         MSP1 = '1'
      ELSE IF ((S1 <= ZERO) .AND. (S2 > ZERO)) THEN
         MS1  = STM/S2U - ONE
         MSP1 = '1'
      ELSE
         MS1  = DMIN1(STM/S1U,STM/S2U) - ONE
         IF ((S1 > ZERO) .AND. (S2 > ZERO)) THEN
            MSP1 = '1'
         ELSE
            MSP1 = '0'
         ENDIF
      ENDIF

      IF (STM < EPS5) THEN
         MSP1 = '0'
      ENDIF

      IF ((S3M < EPS5) .AND. (S4M < EPS5)) THEN            ! Margin for compression
         MS2  = ONEPP10   !(i.e. 1.E+10)
         MSP2 = '0'
      ELSE IF ((S3 >= ZERO) .AND. (S4 < ZERO)) THEN
         MS2  = SCP/S4U - ONE
         MSP2 = '1'
      ELSE IF ((S3 < ZERO) .AND. (S4 >= ZERO)) THEN
         MS2 = SCP/S3U - ONE
         MSP2 = '1'
      ELSE
         MS2 = DMIN1(SCP/S3U,SCP/S4U) - ONE
         IF ((S3 < ZERO) .AND. (S4 < ZERO)) THEN
            MSP2 = '1'
         ELSE
            MSP2 = '0'
         ENDIF
      ENDIF

      IF (SCM < EPS5) THEN
         MSP2 = '0'
      ENDIF

      S5M = DABS(S5)                                       ! ROD1 torsional margin (positive or negative torsional stress)

      IF (S5M < EPS5) THEN
         MS3  = ONEPP10   !(i.e. 1.E+10)
         MSP3 = '0'
      ELSE
         MS3  = (SSM/S5M) - ONE
         MSP3 = '1'
      ENDIF
      IF (SSM < EPS5) THEN
         MSP3 = '0'
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BAR_MARGIN

   END MODULE ROD_BAR_OUTPUT
