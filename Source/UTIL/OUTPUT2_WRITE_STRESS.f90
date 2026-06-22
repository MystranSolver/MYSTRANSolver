
!==================================================================================================
      SUBROUTINE SET_OST_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE)
      ! strain table writer
      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      CHARACTER(8*BYTE), INTENT(IN)    :: ETYPE          ! name of element type
      CHARACTER(8*BYTE), INTENT(INOUT) :: TABLE_NAME     ! name of the op2 table name
      INTEGER(LONG),INTENT(INOUT)      :: ITABLE         ! the subtable
      CALL SET_OES_OST_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE, "OSTR1X  ", "OSTR1X  ")
      END SUBROUTINE SET_OST_TABLE_NAME
!==================================================================================================

      SUBROUTINE SET_OES_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE)
      ! stress table writer
      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      CHARACTER(8*BYTE), INTENT(IN)    :: ETYPE          ! name of element type
      CHARACTER(8*BYTE), INTENT(INOUT) :: TABLE_NAME     ! name of the op2 table name
      INTEGER(LONG),INTENT(INOUT)      :: ITABLE         ! the subtable
      CALL SET_OES_OST_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE, "OES1X   ", "OES1X1  ")
      END SUBROUTINE SET_OES_TABLE_NAME
!==================================================================================================

      SUBROUTINE SET_OES_OST_TABLE_NAME (ETYPE, TABLE_NAME, ITABLE, TABLE_NAME_BAR, TABLE_NAME_SHELL_SOLID)
      ! initializes the OES/OST table name
      ! updates ITABLE and TABLE_NAME
      !
      ! TODO: mak sure the TABLE_NAME_NEW is correct...
      !
      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      USE IOUNT1, ONLY                 :  ERR

      CHARACTER(8*BYTE), INTENT(IN)    :: ETYPE          ! name of element type
      CHARACTER(8*BYTE), INTENT(INOUT) :: TABLE_NAME     ! name of the op2 table name
      CHARACTER(8*BYTE), INTENT(IN)    :: TABLE_NAME_BAR ! name of table for bars
      CHARACTER(8*BYTE), INTENT(IN)    :: TABLE_NAME_SHELL_SOLID ! name of table for shells and solids
      INTEGER(LONG),INTENT(INOUT)      :: ITABLE         ! the subtable
      CHARACTER(8*BYTE)                :: TABLE_NAME_NEW ! name of the op2 table name
      LOGICAL                          :: RETURN_FLAG    ! return from the subroutine early

 1    FORMAT("*DEBUG:      OUTPUT2_WRITE_STRESS: ", A)
 2    FORMAT("*DEBUG:      OUTPUT2_WRITE_STRESS: ", A, "; TABLE_NAME= ", A,";ITABLE=", I8)
 3    FORMAT("*DEBUG:      OUTPUT2_WRITE_STRESS: ", A, " ",A)
      RETURN_FLAG = .TRUE.
      IF      ((ETYPE == 'BAR     ') .OR. (ETYPE == 'BEAM    ')) THEN
        TABLE_NAME_NEW= TABLE_NAME_BAR !"OES1X   "
        RETURN_FLAG = .FALSE.

       ELSE IF ((ETYPE == 'ELAS1   ') .OR. (ETYPE == 'ELAS2   ') .OR. (ETYPE == 'ELAS3   ') .OR. (ETYPE == 'ELAS4   ') .OR.        &
                (ETYPE == 'BUSH    ') .OR. (ETYPE == 'ROD     ') .OR.                                                              &
                (ETYPE == 'TRIA3   ') .OR. (ETYPE == 'QUAD4   ') .OR. (ETYPE == 'SHEAR   ') .OR.                                   &
                (ETYPE == 'HEXA8   ') .OR. (ETYPE == 'PENTA6  ') .OR. (ETYPE == 'TETRA4  ') .OR.                                   &
                (ETYPE == 'HEXA20  ') .OR. (ETYPE == 'PENTA15 ') .OR. (ETYPE == 'TETRA10 ')) THEN

        TABLE_NAME_NEW= TABLE_NAME_SHELL_SOLID !"OES1X1  "
        WRITE(ERR,3) "OES1X1 found",ETYPE
        RETURN_FLAG = .FALSE.
      ELSE
        WRITE(ERR,3) "ERROR STATE",ETYPE
        ! we're now in an error state
        ! also let's close the old table
        TABLE_NAME_NEW= "OES ERR "
        IF (ITABLE < -1) THEN
          WRITE(ERR,2) "closing stress table", TABLE_NAME,ITABLE
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
            WRITE(ERR,2) "closing stress table",TABLE_NAME,ITABLE
            CALL END_OP2_TABLE(ITABLE)   ! close the previous
          ENDIF
          ! we're now at the beginning
          TABLE_NAME = TABLE_NAME_NEW
          ITABLE = -1
          WRITE(ERR,2) "will create stress table",TABLE_NAME,ITABLE
        ENDIF

        ! if we started/restarted, we need to write the TABLE_NAME
        IF (ITABLE == -1) THEN
          WRITE(ERR,2) "creating stress table",TABLE_NAME,ITABLE	
          CALL WRITE_TABLE_HEADER(TABLE_NAME)
          ITABLE = -3
        ENDIF
      ENDIF

100   END SUBROUTINE SET_OES_OST_TABLE_NAME


! ##################################################################################################################################
      SUBROUTINE SET_OESC_TABLE_NAME(TABLE_NAME, ITABLE)
      ! initializes the OES1C table name
      ! updates ITABLE and TABLE_NAME
      !
      ! TODO: mak sure the TABLE_NAME_NEW is correct...
      !
      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      !USE IOUNT1, ONLY                 :  ERR

      !CHARACTER(8*BYTE), INTENT(IN)    :: ETYPE          ! name of element type
      CHARACTER(8*BYTE), INTENT(INOUT) :: TABLE_NAME     ! name of the op2 table name
      !CHARACTER(8*BYTE)                :: TABLE_NAME_NEW ! name of the op2 table name
      INTEGER(LONG)                    :: ITABLE         ! the subtable
      !LOGICAL                          :: RETURN_FLAG    ! return from the subroutine early
      IF (ITABLE == 0) THEN
          CALL WRITE_TABLE_HEADER(TABLE_NAME)
          ITABLE = -3
      ENDIF
      END SUBROUTINE SET_OESC_TABLE_NAME

! ##################################################################################################################################
      SUBROUTINE WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE, &
                                   TITLE, LABEL, SUBTITLE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!
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
!      stress_code : int
!         0 : stress???
!         1 : strain???
!
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT) :: ITABLE               ! an OP2 subtable counter
      INTEGER(LONG), INTENT(IN) :: ISUBCASE                ! the subcase ID
      INTEGER(LONG), INTENT(IN) :: DEVICE_CODE             ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG), INTENT(IN) :: ELEM_TYPE               ! the element type ID
      INTEGER(LONG), INTENT(IN) :: NUM_WIDE                ! the number of words per element
      INTEGER(LONG), INTENT(IN) :: STRESS_CODE             ! flag for von_mises/max_shear/octehedral and fiber_distance/strain_curvature
      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL
      INTEGER(LONG), INTENT(IN)      :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)       :: FIELD6_EIGENVALUE
      INTEGER(LONG)                  :: FORMAT_CODE, ANALYSIS_CODE
!      we assumed static
      ANALYSIS_CODE = 1

      ! LSDVMN
      !FIELD5 = 0

!      static is real
      FORMAT_CODE = 1
      CALL WRITE_OES3(ITABLE, ANALYSIS_CODE, ISUBCASE, DEVICE_CODE, FORMAT_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE,       &
                      TITLE, LABEL, SUBTITLE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
      END SUBROUTINE WRITE_OES3_STATIC

! ##################################################################################################################################
      SUBROUTINE WRITE_OES3_STATIC_AC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE, &
                                      TITLE, LABEL, SUBTITLE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!     Removed in DEV..
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT) :: ITABLE
      INTEGER(LONG), INTENT(IN) :: ISUBCASE
      INTEGER(LONG), INTENT(IN) :: DEVICE_CODE
      INTEGER(LONG), INTENT(IN) :: ANALYSIS_CODE
      INTEGER(LONG), INTENT(IN) :: ELEM_TYPE
      INTEGER(LONG), INTENT(IN) :: NUM_WIDE
      INTEGER(LONG), INTENT(IN) :: STRESS_CODE
      CHARACTER(LEN=128), INTENT(IN) :: TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL
      INTEGER(LONG), INTENT(IN)      :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)       :: FIELD6_EIGENVALUE
      INTEGER(LONG)                  :: FORMAT_CODE

      FORMAT_CODE = 1
      CALL WRITE_OES3(ITABLE, ANALYSIS_CODE, ISUBCASE, DEVICE_CODE, FORMAT_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE, &
                      TITLE, LABEL, SUBTITLE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
      END SUBROUTINE WRITE_OES3_STATIC_AC

! ##################################################################################################################################
      SUBROUTINE WRITE_OES3(ITABLE, ANALYSIS_CODE, ISUBCASE, DEVICE_CODE, FORMAT_CODE, ELEM_TYPE, NUM_WIDE, STRESS_CODE, &
                            TITLE, LABEL, SUBTITLE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
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
      INTEGER(LONG), INTENT(IN) :: STRESS_CODE             ! flag for von_mises/max_shear/octehedral and fiber_distance/strain_curvature

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
 1    FORMAT("WRITE_OES3: ITABLE_START=",I8)
      WRITE(ERR,1) ITABLE	

      IF ((ANALYSIS_CODE == 1) .OR. (ANALYSIS_CODE == 10)) THEN
        ! statics
        FIELD7 = 0.0
      ELSE
         ! frequency in radians
         FIELD7 = SQRT(ABS(FIELD6_EIGENVALUE))
      ENDIF

      WRITE(OP2) 146
      ! stress/strain only
      TABLE_CODE = 5

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
            STRESS_CODE, ACOUSTIC_FLAG, 0, 0, 0, &
            0, 0, 0, 0, 0, &
            0, 0, THERMAL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, &
            TITLE2, SUBTITLE2, LABEL2

      ITABLE = ITABLE - 1        ! flip it to -4, -6, ... so we don't have to do this later
 3    FORMAT("WRITE_OES3: ITABLE_END=",I8)
      WRITE(ERR,3) ITABLE
      CALL WRITE_ITABLE(ITABLE)
      ITABLE = ITABLE - 1
      END SUBROUTINE WRITE_OES3



! ##################################################################################################################################
      SUBROUTINE GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
      ! bit description
      ! 1   is_von_mises (vs. is_max_shear)
      ! 2   is_strain (vs is_stress)
      ! 3   is_strain_curvature (vs. is_fiber_distance)
      ! 4   is_strain (vs is_stress)
      ! 5   material coordinate system flag
      ! stress_bits[1] = 0 -> is_max_shear=True       isVonMises=False
      ! stress_bits[1] = 1 -> is_max_shear=False      isVonMises=True
      ! stress_bits[2] = 0 -> is_stress=True        is_strain=False
      ! stress_bits[3] = 0 -> isFiberCurvature=True isFiberDistance=False
      ! stress_bits[4] = 0 -> duplicate of Bit[1] (stress/strain)
      ! stress_bits[5] = 0 -> material coordinate system flag

      ! s_code =  0 -> stress_bits = [0,0,0,0,0]
      ! s_code =  1 -> stress_bits = [0,0,0,0,1]
      ! s_code =  2 -> stress_bits = [0,0,0,1,0]
      ! s_code =  3 -> stress_bits = [0,0,0,1,1]
      ! etc.
      ! s_code = 32 -> stress_bits = [1,1,1,1,1]
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT)          :: STRESS_CODE
      INTEGER(LONG), INTENT(IN)             :: IS_VON_MISES
      INTEGER(LONG), INTENT(IN)             :: IS_STRAIN
      INTEGER(LONG), INTENT(IN)             :: IS_FIBER_DISTANCE
      INTEGER(LONG)                         :: MATERIAL_FLAG
      MATERIAL_FLAG = 0
      STRESS_CODE = IS_VON_MISES + 2*IS_STRAIN + 4*IS_FIBER_DISTANCE + 8*IS_STRAIN + 16*MATERIAL_FLAG
      END SUBROUTINE GET_STRESS_CODE

! ##################################################################################################################################
