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

   MODULE OP2_FORCE_OUTPUT

   USE OP2_GEOMETRY_OUTPUT, ONLY :  END_OP2_TABLE, WRITE_ITABLE, WRITE_TABLE_HEADER

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SET_OEF_TABLE_NAME, WRITE_OEF3_STATIC

   CONTAINS

!==================================================================================================
      SUBROUTINE SET_OEF_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE)
      ! initializes the OEF table name
      ! updates ITABLE and TABLE_NAME
      !
      ! TODO: mak sure the TABLE_NAME_NEW is correct...
      !
      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      USE IOUNT1, ONLY                 :  ERR

      CHARACTER(8*BYTE), INTENT(IN)    :: ETYPE          ! name of element type
      CHARACTER(8*BYTE), INTENT(INOUT) :: TABLE_NAME     ! name of the op2 table name
      CHARACTER(8*BYTE)                :: TABLE_NAME_NEW ! name of the op2 table name
      INTEGER(LONG)                    :: ITABLE         ! the subtable
      LOGICAL                          :: RETURN_FLAG    ! return from the subroutine early

 1    FORMAT("*DEBUG:      OUTPUT2_WRITE_FORCE:  ", A)
 2    FORMAT("*DEBUG:      OUTPUT2_WRITE_FORCE:  ", A, "; TABLE_NAME= ", A,";ITABLE=", I8)
 3    FORMAT("*DEBUG:      OUTPUT2_WRITE_FORCE:  ", A, " ",A)
      RETURN_FLAG = .TRUE.
      IF      ((ETYPE == 'BAR     ') .OR. (ETYPE == 'BEAM    ')) THEN
        TABLE_NAME_NEW= "OEF1X   "
        RETURN_FLAG = .FALSE.
      ELSE IF ((ETYPE == 'ELAS1   ') .OR. (ETYPE == 'ELAS2   ') .OR. (ETYPE == 'ELAS3   ') .OR. (ETYPE == 'ELAS4   ') .OR.         &
               (ETYPE == 'BUSH    ') .OR. (ETYPE == 'ROD     ') .OR.                                                               &
               (ETYPE == 'TRIA3   ') .OR. (ETYPE == 'QUAD4   ') .OR. (ETYPE == 'QUAD8   ') .OR. (ETYPE == 'SHEAR   ') .OR.         &
               (ETYPE == 'HEXA8   ') .OR. (ETYPE == 'PENTA6  ') .OR. (ETYPE == 'TETRA4  ') .OR.                                    &
               (ETYPE == 'HEXA20  ') .OR. (ETYPE == 'PENTA15 ') .OR. (ETYPE == 'TETRA10 ')) THEN
        TABLE_NAME_NEW= "OEF1X   "
        WRITE(ERR,3) "OEF1X found",ETYPE
        RETURN_FLAG = .FALSE.
      ELSE
        WRITE(ERR,3) "ERROR STATE",ETYPE
        ! we're now in an error state
        ! also let's close the old table
        TABLE_NAME_NEW= "OEF ERR "
        IF (ITABLE < -1) THEN
          WRITE(ERR,2) "closing force table", TABLE_NAME,ITABLE
          CALL END_OP2_TABLE(ITABLE)   ! close the previous
        ENDIF
!        WRITE(ERR,2) "invalidated tableA",TABLE_NAME,ITABLE
        ITABLE = 0
      ENDIF

      IF (RETURN_FLAG) THEN
        !WRITE(ERR,1) "early return"
        TABLE_NAME = TABLE_NAME_NEW
      ELSE
        !WRITE(ERR,1) "will write table header"
        IF (TABLE_NAME /= TABLE_NAME_NEW) THEN
        ! first let's find out if we need to close off the previous table
          IF (ITABLE < -1) THEN
            WRITE(ERR,2) "closing force table",TABLE_NAME,ITABLE
            CALL END_OP2_TABLE(ITABLE)   ! close the previous
          ENDIF
          ! we're now at the beginning
          TABLE_NAME = TABLE_NAME_NEW
          ITABLE = -1
          WRITE(ERR,2) "will create force table",TABLE_NAME,ITABLE
        ENDIF

        ! if we started/restarted, we need to write the TABLE_NAME
        IF (ITABLE == -1) THEN
          WRITE(ERR,2) "creating force table",TABLE_NAME,ITABLE
          CALL WRITE_TABLE_HEADER(TABLE_NAME)
          ITABLE = -3
        ENDIF
      ENDIF

100   END SUBROUTINE SET_OEF_TABLE_NAME


! ##################################################################################################################################
      SUBROUTINE WRITE_OEF3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEM_TYPE, NUM_WIDE, &
                                   TITLE, LABEL, SUBTITLE,                                            &
                                   FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!      Parameters
!      ==========
!      op2 : file int
!         the file pointer
!      itable : int
!         the subtable number
!          - starts as -3, -5, ...
!          - ends as -4, -6, ...
!      isubcase : int
!         the subcase id (from the bulk data deck)
!      etype : int
!         element type
!      num_wide : int
!         the number of fields in each element
!
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT) :: ITABLE               ! an OP2 subtable counter
      INTEGER(LONG), INTENT(IN) :: ISUBCASE                ! the subcase ID
      INTEGER(LONG), INTENT(IN) :: DEVICE_CODE             ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG), INTENT(IN) :: ELEM_TYPE               ! the element type ID
      INTEGER(LONG), INTENT(IN) :: NUM_WIDE                ! the number of words per element
      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL
      INTEGER(LONG), INTENT(IN)      :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)       :: FIELD6_EIGENVALUE
      INTEGER(LONG)                  :: FORMAT_CODE, ANALYSIS_CODE

!      static is real
      FORMAT_CODE = 1
      CALL WRITE_OEF3(ITABLE, ANALYSIS_CODE, ISUBCASE, DEVICE_CODE, FORMAT_CODE, ELEM_TYPE, NUM_WIDE,       &
                      TITLE, LABEL, SUBTITLE,                                                               &
                      FIELD5_INT_MODE, FIELD6_EIGENVALUE)
      END SUBROUTINE WRITE_OEF3_STATIC

! ##################################################################################################################################
      SUBROUTINE WRITE_OEF3(ITABLE, ANALYSIS_CODE, ISUBCASE, DEVICE_CODE, FORMAT_CODE, ELEM_TYPE, NUM_WIDE, &
                            TITLE, LABEL, SUBTITLE,                                                         &
                            FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!      Parameters
!      ==========
!      analysis_code
!        the solution type flag
!      approach_code
!        ???
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR,OP2
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT) :: ITABLE               ! an OP2 subtable counter
      INTEGER(LONG), INTENT(IN) :: ANALYSIS_CODE           ! static, modal, time, freq, etc. flag
      INTEGER(LONG), INTENT(IN) :: ISUBCASE                ! the subcase ID
      INTEGER(LONG), INTENT(IN) :: DEVICE_CODE             ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG), INTENT(IN) :: ELEM_TYPE               ! the element type ID
      INTEGER(LONG), INTENT(IN) :: FORMAT_CODE
      INTEGER(LONG), INTENT(IN) :: NUM_WIDE                ! the number of words per element
      !INTEGER(LONG), INTENT(IN) :: STRESS_CODE             ! flag for von_mises/max_shear/octehedral and fiber_distance/strain_curvature

      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL
      CHARACTER(LEN=128) :: TITLE2                 ! the subcase TITLE
      CHARACTER(LEN=128) :: SUBTITLE2              ! the subcase SUBTITLE
      CHARACTER(LEN=128) :: LABEL2                 ! the subcase LABEL

      INTEGER(LONG), INTENT(IN) :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)  :: FIELD6_EIGENVALUE
      REAL(DOUBLE)              :: FIELD7
      INTEGER(LONG)             :: APPROACH_CODE
      INTEGER(LONG)             :: TABLE_CODE, LOAD_SET, THERMAL, ACOUSTIC_FLAG
      TITLE2 = TITLE(1:100)
      SUBTITLE2 = SUBTITLE(1:67)
      LABEL2 = LABEL(1:100)

      CALL WRITE_ITABLE(ITABLE)  ! write the -3, -5, ... subtable header
 1    FORMAT("WRITE_OEF3: ITABLE_START=",I8)
      WRITE(ERR,1) ITABLE

      IF ((ANALYSIS_CODE == 1) .OR. (ANALYSIS_CODE == 10)) THEN
        ! statics
        FIELD7 = 0.0
      ELSE
         ! frequency in radians
         FIELD7 = SQRT(ABS(FIELD6_EIGENVALUE))
      ENDIF

      WRITE(OP2) 146
      ! force only
      TABLE_CODE = 4

      ! ???
      LOAD_SET = 1

      ! we're not doing acoustic
      ACOUSTIC_FLAG = 0

      ! not always 0 for stress, but for now
      THERMAL = 0

      APPROACH_CODE = ANALYSIS_CODE * 10 + DEVICE_CODE
2     FORMAT(" APPROACH_CODE=",I4," TABLE_CODE=",I4," ELEM_TYPE=",I4," ISUBCASE=",I4)
      ! 584 bytes
      WRITE(ERR,2) APPROACH_CODE, TABLE_CODE, ELEM_TYPE, ISUBCASE
      WRITE(OP2) APPROACH_CODE, TABLE_CODE, ELEM_TYPE, ISUBCASE, FIELD5_INT_MODE, &
            REAL(FIELD6_EIGENVALUE, 4), REAL(FIELD7, 4),                          &
            LOAD_SET, FORMAT_CODE, NUM_WIDE, &
            0, ACOUSTIC_FLAG, 0, 0, 0, &
            0, 0, 0, 0, 0, &
            0, 0, THERMAL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, &
            TITLE2, SUBTITLE2, LABEL2

      ITABLE = ITABLE - 1        ! flip it to -4, -6, ... so we don't have to do this later
 3    FORMAT("WRITE_OEF3: ITABLE_END=",I8)
      WRITE(ERR,3) ITABLE
      CALL WRITE_ITABLE(ITABLE)
      ITABLE = ITABLE - 1
      END SUBROUTINE WRITE_OEF3



! ##################################################################################################################################

   END MODULE OP2_FORCE_OUTPUT
