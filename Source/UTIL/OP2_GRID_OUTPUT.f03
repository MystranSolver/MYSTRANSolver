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

   MODULE OP2_GRID_OUTPUT

   USE OP2_GEOMETRY_OUTPUT, ONLY :  WRITE_ITABLE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_OUG3_STATIC, WRITE_OUG3_EIGN

   CONTAINS
! ##################################################################################################################################
      SUBROUTINE WRITE_OUG3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NEW_RESULT, &
                                   TITLE, SUBTITLE, LABEL)
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN)        :: ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE
      CHARACTER(LEN=128), INTENT(IN)   :: TITLE            ! Solution title
      CHARACTER(LEN=128), INTENT(IN)   :: SUBTITLE         ! Subcase subtitle
      CHARACTER(LEN=128), INTENT(IN)   :: LABEL            ! Subcase label
      LOGICAL, INTENT(INOUT)           :: NEW_RESULT
      INTEGER(LONG)                    :: FORMAT_CODE
      INTEGER(LONG)                    :: NUM_WIDE
      INTEGER(LONG)                    :: LSDVMN  ! field 5
      REAL(DOUBLE)                     :: FIELD6

!      1 = OP2 is the output file
!      DEVICE_CODE = 1

!     real
      FORMAT_CODE = 1
      NUM_WIDE = 8

      LSDVMN = 1
      FIELD6 = 0.0
      CALL WRITE_OUG3(ITABLE, ISUBCASE, FORMAT_CODE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NUM_WIDE, NEW_RESULT, &
                      TITLE, SUBTITLE, LABEL, LSDVMN, FIELD6)
      END SUBROUTINE WRITE_OUG3_STATIC


! ##################################################################################################################################
    SUBROUTINE WRITE_OUG3_EIGN(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NEW_RESULT, &
                               TITLE, SUBTITLE, LABEL, MODE, EIGENVALUE)
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN)        :: ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE
      INTEGER(LONG), INTENT(INOUT)     :: TABLE_CODE
      CHARACTER(LEN=128), INTENT(IN)   :: TITLE            ! Solution title
      CHARACTER(LEN=128), INTENT(IN)   :: SUBTITLE         ! Subcase subtitle
      CHARACTER(LEN=128), INTENT(IN)   :: LABEL            ! Subcase label
      LOGICAL, INTENT(INOUT)           :: NEW_RESULT
      INTEGER(LONG)                    :: FORMAT_CODE
      INTEGER(LONG)                    :: NUM_WIDE
      INTEGER(LONG)                    :: NINT = 4

      INTEGER(LONG), INTENT(IN)        :: MODE       ! field 5
      REAL(DOUBLE), INTENT(IN)         :: EIGENVALUE ! field 6

!      1 = OP2 is the output file
!      DEVICE_CODE = 1

      ! eigenvector
      IF (TABLE_CODE == 1) THEN
          ! displacement -> eigenvector
          ! don't update spc/mpc forces
          TABLE_CODE = 7
      ENDIF

      ! real
      FORMAT_CODE = 1
      NUM_WIDE = 8

      CALL WRITE_OUG3(ITABLE, ISUBCASE, FORMAT_CODE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NUM_WIDE, NEW_RESULT, &
                      TITLE, SUBTITLE, LABEL, MODE, EIGENVALUE)
      END SUBROUTINE WRITE_OUG3_EIGN

! ##################################################################################################################################
      SUBROUTINE WRITE_OUG3(ITABLE, ISUBCASE, FORMAT_CODE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NUM_WIDE, NEW_RESULT, &
                            TITLE, SUBTITLE, LABEL, FIELD5, FIELD6)
!      Parameters
!      ==========
!      isubcase : int
!         the subcase id
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY           :  ERR, OP2
      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN)        :: ITABLE, ISUBCASE, FORMAT_CODE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NUM_WIDE
      CHARACTER(LEN=128), INTENT(IN)   :: TITLE            ! Solution title
      CHARACTER(LEN=128), INTENT(IN)   :: SUBTITLE         ! Subcase subtitle
      CHARACTER(LEN=128), INTENT(IN)   :: LABEL            ! Subcase label
      INTEGER(LONG), INTENT(IN)        :: FIELD5
      REAL(DOUBLE), INTENT(IN)         :: FIELD6
      REAL(DOUBLE)                     :: FIELD7
      LOGICAL, INTENT(INOUT)           :: NEW_RESULT
      CHARACTER(LEN=128) :: TITLE2                 ! the subcase TITLE
      CHARACTER(LEN=128) :: SUBTITLE2              ! the subcase SUBTITLE
      CHARACTER(LEN=128) :: LABEL2                 ! the subcase LABEL

      INTEGER(LONG) :: THERMAL, SORT_CODE, &
        RANDOM_CODE, ACOUSTIC_FLAG, OCODE, &
        APPROACH_CODE, TCODE
      INTEGER(LONG)                    :: NINT = 4

      TITLE2 = TITLE(1:100)
      SUBTITLE2 = SUBTITLE(1:67)
      LABEL2 = LABEL(1:100)

      WRITE(*,*) NEW_RESULT, itable
      IF(NEW_RESULT .AND. (ITABLE .NE. -3)) THEN
!        header = [
!          4, 146, 4,
!        ]
        WRITE(OP2) 146
        NEW_RESULT = .FALSE.
      ELSE
!        header = [
!          4, itable, 4,
!          4, 1, 4,
!          4, 0, 4,
!          4, 146, 4,
!        ]
        CALL WRITE_ITABLE(ITABLE)
        !WRITE(OP2) ITABLE
        !WRITE(OP2) 1
        !WRITE(OP2) 0
        WRITE(OP2) 146
      ENDIF
      ! displacement only
      !if TABLE_NAME .EQ. "OUG1"
      !TABLE_CODE = 1
      THERMAL = 0

!       static
      !ANALYSIS_CODE = 1

!      sort1
!      should this be 1???
      SORT_CODE = 0

!      other minor assumptions
      RANDOM_CODE = 0
      ACOUSTIC_FLAG = 0
      OCODE = 0

!      table_code = tCode % 1000
!      sort_code = tCode // 1000
!      approach, tcode, int3, isubcase
      TCODE = SORT_CODE * 1000 + TABLE_CODE

!      if isinstance(acoustic_flag, float_types):
!          ftable3 = set_table3_field(ftable3, 12, b'f') # field 11

      IF(ANALYSIS_CODE == 1) THEN
!        what actually is this??
!        field5 = self.lsdvmns[itime]
!       statics
!        FIELD5 = 1
!        FIELD6 = 0
        FIELD7 = 0.0
      ELSE IF(ANALYSIS_CODE == 2) THEN
!       modes - real eigenvalues
        !FIELD5 = MODE             ! int
        !FIELD6 = EIGN             ! float

!       abs_freqs = np.sqrt(np.abs(self.eigns)) / (2 * np.pi)
!       column_names.append('Freq')
!       column_values.append(abs_freqs)
!       column_names.append('Eigenvalue')
!       column_values.append(times)
!       column_names.append('Radians')
!       column_values.append(abs_freqs * 2 * np.pi)  ! radians
!       ABS_FREQ = SQRT(ABS(FIELD6)) / (2 * PI)
!       FIELD7 = ABS_FREQ * 2 * PI
        FIELD7 = SQRT(ABS(FIELD6))  ! SQRT(ABS(EIGN))

        !FIELD7 = MODE_CYCLE_FREQ  ! float; radians
      ELSE IF(ANALYSIS_CODE == 7) THEN
!       pre-buckling????
!        op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
        !FIELD5 = 1

        ! do nothing statement
        APPROACH_CODE = 0

        ! TODO gotta figure out the appropriate value
        FIELD7 = 0.0
      ELSE
        WRITE(ERR,100) ANALYSIS_CODE
        STOP
!        raise NotImplementedError(self.analysis_code)
      ENDIF

      APPROACH_CODE = ANALYSIS_CODE*10 + DEVICE_CODE
!      table3 = [
!          approach_code, table_code, 0, isubcase, field5,
!          field6, field7, random_code, format_code, num_wide,
!          oCode, acoustic_flag, 0, 0, 0,
!          0, 0, 0, 0, 0,
!          0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
!          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
!          0, 0, 0, 0,
!          title, subtitle, label,
!      ]
!      WRITE(OP2)
      WRITE(OP2) APPROACH_CODE, TABLE_CODE, 0, ISUBCASE, FIELD5, &
          REAL(FIELD6, 4), REAL(FIELD7, 4),                      &
          RANDOM_CODE, FORMAT_CODE, NUM_WIDE, &
          OCODE, ACOUSTIC_FLAG, 0, 0, 0, &
          0, 0, 0, 0, 0, &
          0, 0, THERMAL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
          0, 0, 0, 0, &
          TITLE2, SUBTITLE2, LABEL2
!      assert table3[22] == thermal
!      print *, "finished OUGV1 table=3"

      CALL WRITE_ITABLE(ITABLE)

 100  FORMAT(" analysis_code=", i4, " is not supported...stopping")
      END SUBROUTINE WRITE_OUG3

   END MODULE OP2_GRID_OUTPUT
