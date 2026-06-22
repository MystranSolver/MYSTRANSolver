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

      SUBROUTINE WRITE_FEMAP_GRID_VECS ( GRID_VEC, FEMAP_SET_ID, WHAT )

! Writes grid related vectors to FEMAP neutral file (displ, applied load, SPC and MPC forces)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, NEU
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD, NDOFG, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  CORD, GRID, GRID_ID, INV_GRID_SEQ

      USE WRITE_FEMAP_GRID_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_GRID_VECS'
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT              ! Indicator if GRID_VEC is DISP, OLOA, SPCF or MPCF
      CHARACTER(LEN= 5*BYTE)          :: TITLE1(4,2)       ! Titles for vectors written to NEU
      CHARACTER(LEN=20*BYTE)          :: TITLE2(2)         ! Titles for vectors written to NEU

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! FEMAP set ID to write out
      INTEGER(LONG)                   :: ACID_G            ! Actual coordinate system ID for a grid
      INTEGER(LONG)                   :: GRID_MAX          ! Grid ID where vector is max
      INTEGER(LONG)                   :: GRID_MIN          ! Grid ID where vector is min
      INTEGER(LONG), ALLOCATABLE      :: GRID_NUMS(:)      ! Grid ID's in global order
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICID              ! Internal coord sys no. corresponding to an actual coord sys no.
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID for a grid in array GRID_NUMS
      INTEGER(LONG), ALLOCATABLE      :: IARRAY(:)         ! Original GRID_NUMS array
      INTEGER(LONG)                   :: IDOFG             ! A G-set DOF number
      INTEGER(LONG)                   :: ID(20)            ! Vector ID's for FEMAP output
      INTEGER(LONG)                   :: J                 ! Counter
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: VEC_ID_OFFSET     ! Offset in determining output vector ID
      INTEGER(LONG)                   :: VEC_ID            ! Vector ID for FEMAP output


      REAL(DOUBLE) , INTENT(IN)       :: GRID_VEC(NDOFG)   ! G-set Vector to process
      REAL(DOUBLE)                    :: DIS(3)            ! Array of 3 translation components
      REAL(DOUBLE)                    :: ROT(3)            ! Array of 3 rotation components
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE), ALLOCATABLE       :: TOTR_VEC(:)       ! RSS of 6 rotation    components in  GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: TOTT_VEC(:)       ! RSS of 6 translation components in  GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: T1_VEC(:)         ! T1 translation component from GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: T2_VEC(:)         ! T2 translation component from GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: T3_VEC(:)         ! T3 translation component from GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: R1_VEC(:)         ! R1 rotation    component from GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: R2_VEC(:)         ! R2 rotation    component from GRID_VEC
      REAL(DOUBLE), ALLOCATABLE       :: R3_VEC(:)         ! R3 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: T0G(3,3)           ! Matrix to transform offsets from global to basic  coords
      REAL(DOUBLE)                    :: VEC_ABS           ! Abs value in vector
      REAL(DOUBLE)                    :: VEC_MAX           ! Max value in vector
      REAL(DOUBLE)                    :: VEC_MIN           ! Min value in vector



! **********************************************************************************************************************************
      TITLE1(1,1) = 'RSS'
      TITLE1(2,1) = 'T1'
      TITLE1(3,1) = 'T2'
      TITLE1(4,1) = 'T3'
      TITLE1(1,2) = 'RSS'
      TITLE1(2,2) = 'R1'
      TITLE1(3,2) = 'R2'
      TITLE1(4,2) = 'R3'

      IF      (WHAT == 'DISP') THEN
! --- neu_upgrade begin --- !
         VEC_ID_OFFSET = 0
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Translation'
         TITLE2(2) = ' Rotation'
! --- neu_upgrade end --- !
      ELSE IF (WHAT == 'OLOA') THEN
         VEC_ID_OFFSET = 20000
          TITLE1(1,1) = 'Total'
          TITLE1(1,2) = 'Total'
          TITLE2(1) = ' Applied Force'
          TITLE2(2) = ' Applied Moment'
      ELSE IF (WHAT == 'SPCF') THEN
         VEC_ID_OFFSET = 30000
! --- neu_upgrade begin --- !
          TITLE1(1,1) = 'Total'
          TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Constraint Force'
         TITLE2(2) = ' Constraint Moment'
! --- neu_upgrade end --- !
      ELSE IF (WHAT == 'MPCF') THEN
         VEC_ID_OFFSET = 40000
         TITLE2(1) = ' MPC force'
         TITLE2(2) = ' MPC moment'
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,939) SUBR_NAME, WHAT
         WRITE(F06,939) SUBR_NAME, WHAT
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      ALLOCATE ( GRID_NUMS(NGRID), IARRAY(NGRID), TOTT_VEC(NGRID), TOTR_VEC(NGRID), T1_VEC(NGRID), T2_VEC(NGRID), T3_VEC(NGRID), &
                 R1_VEC(NGRID), R2_VEC(NGRID), R3_VEC(NGRID) )

      IDOFG = 0
      DO I=1,NGRID

         T1_VEC(I)    = ZERO
         T2_VEC(I)    = ZERO
         T3_VEC(I)    = ZERO
         R1_VEC(I)    = ZERO
         R2_VEC(I)    = ZERO
         R3_VEC(I)    = ZERO

         GRID_NUMS(I) = GRID_ID(INV_GRID_SEQ(I))
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS == 6) THEN                          ! Grid point, 6 DOF
            IDOFG = IDOFG + 1   ;  T1_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1   ;  T2_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1   ;  T3_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1   ;  R1_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1   ;  R2_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1   ;  R3_VEC(I)    = GRID_VEC(IDOFG)
         ELSE                                              ! Scalar point, only 1 DOF
            IDOFG = IDOFG + 1   ;  T1_VEC(I)    = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 0   ;  T2_VEC(I)    = ZERO
            IDOFG = IDOFG + 0   ;  T3_VEC(I)    = ZERO
            IDOFG = IDOFG + 0   ;  R1_VEC(I)    = ZERO
            IDOFG = IDOFG + 0   ;  R2_VEC(I)    = ZERO
            IDOFG = IDOFG + 0   ;  R3_VEC(I)    = ZERO
         ENDIF

         TOTR_VEC(I)  = DSQRT( R1_VEC(I)*R1_VEC(I) + R2_VEC(I)*R2_VEC(I) + R3_VEC(I)*R3_VEC(I) )
         TOTT_VEC(I)  = DSQRT( T1_VEC(I)*T1_VEC(I) + T2_VEC(I)*T2_VEC(I) + T3_VEC(I)*T3_VEC(I) )

         IF (WHAT == 'DISP') THEN                          ! 10/01/14: Code to transform displs from global to basic
            DIS(1) = T1_VEC(I)
            DIS(2) = T2_VEC(I)
            DIS(3) = T3_VEC(I)
            ROT(1) = R1_VEC(I)
            ROT(2) = R2_VEC(I)
            ROT(3) = R3_VEC(I)

            IGRID  = INV_GRID_SEQ(I)                       ! Get transformation matrix for displs from global to basic (TG0)
            ACID_G = GRID(IGRID,3)                         ! Get global coord sys for this grid
            IF (ACID_G /= 0) THEN                          ! Global is not basic so need to transform offset from basic to global
               ICID = 0
               DO J=1,NCORD
                  IF (ACID_G == CORD(J,2)) THEN
                     ICID = J
                     EXIT
                  ENDIF
               ENDDO
               CALL GEN_T0L ( IGRID, ICID, THETAD, PHID, T0G )
               T1_VEC(I) = T0G(1,1)*DIS(1) + T0G(1,2)*DIS(2) + T0G(1,3)*DIS(3)
               T2_VEC(I) = T0G(2,1)*DIS(1) + T0G(2,2)*DIS(2) + T0G(2,3)*DIS(3)
               T3_VEC(I) = T0G(3,1)*DIS(1) + T0G(3,2)*DIS(2) + T0G(3,3)*DIS(3)
               R1_VEC(I) = T0G(1,1)*ROT(1) + T0G(1,2)*ROT(2) + T0G(1,3)*ROT(3)
               R2_VEC(I) = T0G(2,1)*ROT(1) + T0G(2,2)*ROT(2) + T0G(2,3)*ROT(3)
               R3_VEC(I) = T0G(3,1)*ROT(1) + T0G(3,2)*ROT(2) + T0G(3,3)*ROT(3)
            ELSE                                           ! Global was basic so no transformation of coords needed
               T1_VEC(I) = DIS(1)
               T2_VEC(I) = DIS(2)
               T3_VEC(I) = DIS(3)
               R1_VEC(I) = ROT(1)
               R2_VEC(I) = ROT(2)
               R3_VEC(I) = ROT(3)
            ENDIF

         ENDIF

      ENDDO

! Write total translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 1
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(1,1), TITLE2(1)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, TOTT_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = VEC_ID + 1
      ID(2) = VEC_ID + 2
      ID(3) = VEC_ID + 3
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, TOTT_VEC', NGRID, IARRAY, TOTT_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), TOTT_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write T1 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 2
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(2,1), TITLE2(1)
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, T1_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = VEC_ID
      ID(2) = 0
      ID(3) = 0
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, T1_VEC', NGRID,  IARRAY, T1_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), T1_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write T2 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 3
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(3,1), TITLE2(1)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, T2_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = 0
      ID(2) = VEC_ID
      ID(3) = 0
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, T2_VEC', NGRID,  IARRAY, T2_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), T2_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write T3 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 4
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(4,1), TITLE2(1)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, T3_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(2) = 0
      ID(1) = 0
      ID(3) = VEC_ID
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, T3_VEC', NGRID,  IARRAY, T3_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), T3_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write total rotation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 5
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(1,2), TITLE2(2)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, TOTR_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = VEC_ID + 1
      ID(2) = VEC_ID + 2
      ID(3) = VEC_ID + 3
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, TOTR_VEC', NGRID,  IARRAY, TOTR_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), TOTR_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write R1 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 6
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(2,2), TITLE2(2)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, R1_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = VEC_ID
      ID(2) = 0
      ID(3) = 0
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, R1_VEC', NGRID,  IARRAY, R1_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), R1_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write R2 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 7
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(3,2), TITLE2(2)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, R2_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(1) = 0
      ID(2) = VEC_ID
      ID(3) = 0
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, R2_VEC', NGRID,  IARRAY, R2_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), R2_VEC(I)
      ENDDO
      WRITE(NEU,1007)

! Write R3 translation vector output to FEMAP neutral file

      VEC_ID = VEC_ID_OFFSET + 8
      WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
      WRITE(NEU,1002) TITLE1(4,2), TITLE2(2)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, R3_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
      ID(2) = 0
      ID(1) = 0
      ID(3) = VEC_ID
      DO I=4,20
         ID(I) = 0
      ENDDO
      WRITE(NEU,1004) (ID(I),I= 1,10)
      WRITE(NEU,1004) (ID(I),I=11,20)
      WRITE(NEU,1005) GRID_MIN, GRID_MAX
      DO I=1,NGRID
         IARRAY(I) = GRID_NUMS(I)
      ENDDO
      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, R3_VEC', NGRID,  IARRAY, R3_VEC )
      DO I=1,NGRID
         WRITE(NEU,1006) IARRAY(I), R3_VEC(I)
      ENDDO
      WRITE(NEU,1007)



      IF (ALLOCATED(GRID_NUMS)) DEALLOCATE(GRID_NUMS)
      IF (ALLOCATED(IARRAY   )) DEALLOCATE(IARRAY)
      IF (ALLOCATED(TOTT_VEC )) DEALLOCATE(TOTT_VEC)
      IF (ALLOCATED(TOTR_VEC )) DEALLOCATE(TOTR_VEC)
      IF (ALLOCATED(T1_VEC   )) DEALLOCATE(T1_VEC)
      IF (ALLOCATED(T2_VEC   )) DEALLOCATE(T2_VEC)
      IF (ALLOCATED(T3_VEC   )) DEALLOCATE(T3_VEC)
      IF (ALLOCATED(R1_VEC   )) DEALLOCATE(R1_VEC)
      IF (ALLOCATED(R2_VEC   )) DEALLOCATE(R2_VEC)
      IF (ALLOCATED(R3_VEC   )) DEALLOCATE(R3_VEC)

      RETURN

! **********************************************************************************************************************************
  939 FORMAT(' *ERROR   939: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' WRONG VALUE = ',A,' FOR ARGUMENT "WHAT"')

 1001 FORMAT(2(I8,','),'       1,')

 1002 FORMAT(2A)

 1003 FORMAT(3(1ES17.6,','))

 1004 FORMAT(10(I8,','))

 1005 FORMAT(2(I8,','),'       1,       7,',/,'       1,       1,       1')

 1006 FORMAT(I8,',',1ES17.6,',')

 1007 FORMAT('      -1,     0.          ,')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_FEMAP_GRID_VECS
