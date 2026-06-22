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
 
      SUBROUTINE WRITE_FEMAP_ELFO_VECS ( ELEM_TYP, NUM_FEMAP_ROWS, FEMAP_SET_ID )
 
! Writes element engineering forces to FEMAP neutral file for ROD, BAR,TRIA3, QUAD4, SHEAR   

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, NEU
      USE PARAMS, ONLY                :  SUPWARN
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NGRID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
 
      USE WRITE_FEMAP_ELFO_VECS_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_ELFO_VECS'
      CHARACTER( 1*BYTE)              :: CALC_WARN              ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: CENT_TOTAL             ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: COMP_DIR               ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: ENT_TYPE = '8'         ! FEMAP value for record 6
      CHARACTER( 1*BYTE)              :: OUT_TYPE = '3'         ! FEMAP value for record 6
      CHARACTER(LEN=*), INTENT(IN)    :: ELEM_TYP               ! Element type
      CHARACTER(LEN=LEN(ELEM_TYP))    :: ELEM_NAME              ! ELEM_TYP with trailing blanks stripped
      CHARACTER(25*BYTE)              :: TITLE_E(12)            ! Titles for vectors written to NEU

      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS         ! Number of rows of FEMAP data to write
      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID           ! FEMAP set ID to write out
      INTEGER(LONG)                   :: ELEM_MAX               ! Grid ID where vector is max
      INTEGER(LONG)                   :: ELEM_MIN               ! Grid ID where vector is min

                                                                ! Col from FEMAP_EL_NUMS (elem ID's)
      INTEGER(LONG), ALLOCATABLE      :: ELEM_NUMS(:)

      INTEGER(LONG)                   :: ELEM_NAME_LEN          ! Length of ELEM_TYP without trailing blanks
      INTEGER(LONG)                   :: I,J                    ! DO loop indices
      INTEGER(LONG)                   :: ID(20)                 ! Vector ID's for FEMAP output
      INTEGER(LONG)                   :: VEC_ID_OFFSET          ! Offset in determining output vector ID
      INTEGER(LONG)                   :: VEC_ID                 ! Vector ID for FEMAP output

                                                                ! Columns from FEMAP_EL_VECS
      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VECS(:,:)

                                                                ! One column from FEMAP_EL_VECS
      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VEC(:)

      REAL(DOUBLE)                    :: VEC_ABS                ! Abs value in vector
      REAL(DOUBLE)                    :: VEC_MAX                ! Max value in vector
      REAL(DOUBLE)                    :: VEC_MIN                ! Min value in vector
 

! **********************************************************************************************************************************
      ELEM_NAME_LEN = LEN(ELEM_TYP)
      ELEM_NAME(1:) = ELEM_TYP(1:)
      DO I=LEN(ELEM_TYP),1,-1
         IF (ELEM_TYP(I:I) == ' ') THEN
            CYCLE
         ELSE
            ELEM_NAME_LEN = I
            EXIT
         ENDIF
      ENDDO

      ALLOCATE ( ELEM_NUMS(NUM_FEMAP_ROWS) )
      ALLOCATE ( ELEM_VECS(NUM_FEMAP_ROWS,12) )
      ALLOCATE ( ELEM_VEC(NUM_FEMAP_ROWS) )

      IF      (ELEM_TYP == 'ROD     ') THEN
         VEC_ID_OFFSET = 50100
      ELSE IF (ELEM_TYP == 'BAR     ') THEN
         VEC_ID_OFFSET = 50200
      ELSE IF (ELEM_TYP == 'TRIA3K  ') THEN
         VEC_ID_OFFSET = 50300
      ELSE IF (ELEM_TYP == 'TRIA3   ') THEN
         VEC_ID_OFFSET = 50400
      ELSE IF (ELEM_TYP == 'QUAD4K  ') THEN
         VEC_ID_OFFSET = 50500
      ELSE IF (ELEM_TYP == 'QUAD4   ') THEN
         VEC_ID_OFFSET = 50600
      ELSE IF (ELEM_TYP == 'SHEAR   ') THEN
         VEC_ID_OFFSET = 50700
      ELSE IF (ELEM_TYP == 'ELAS1   ') THEN
         VEC_ID_OFFSET = 50800
      ELSE IF (ELEM_TYP == 'ELAS2   ') THEN
         VEC_ID_OFFSET = 50900
      ELSE IF (ELEM_TYP == 'ELAS3   ') THEN
         VEC_ID_OFFSET = 51000
      ELSE IF (ELEM_TYP == 'ELAS4   ') THEN
         VEC_ID_OFFSET = 51100
      ELSE IF (ELEM_TYP == 'BUSH    ') THEN
         VEC_ID_OFFSET = 51200
!--- cbeam add begin --- !
      ELSE IF (ELEM_TYP == 'BEAM    ') THEN
         VEC_ID_OFFSET = 51300
!--- cbeam add end --- !
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,943) TRIM(ELEM_TYP), 'ELEM FORCE', TRIM(SUBR_NAME)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,943) TRIM(ELEM_TYP), 'ELEM FORCE', TRIM(SUBR_NAME)
         ENDIF
      ENDIF

! Process BAR and ROD elements
!--- cbeam add begin --- !
      IF (ELEM_TYP == 'BEAM    ') THEN

         TITLE_E( 1) = 'EndA Plane1 Moment'
         TITLE_E( 2) = 'EndA Plane2 Moment'
         TITLE_E( 3) = 'EndB Plane1 Moment'
         TITLE_E( 4) = 'EndB Plane2 Moment'
         TITLE_E( 5) = 'EndA Pl1 Shear Force'
         TITLE_E( 6) = 'EndA Pl2 Shear Force'
         TITLE_E( 7) = 'EndB Pl1 Shear Force'
         TITLE_E( 8) = 'EndB Pl2 Shear Force'
         TITLE_E( 9) = 'EndA Axial Force'
         TITLE_E(10) = 'EndB Axial Force'
         TITLE_E(11) = 'EndA Torque'
         TITLE_E(12) = 'EndB Torque'

         DO J=1,12
            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            DO I=1,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            CALC_WARN  = '0'
            COMP_DIR   = '3'
            CENT_TOTAL = '1'
            WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)
         ENDDO
!--- cbeam add end --- !
      ELSE IF ((ELEM_TYP == 'BAR     ') .OR. (ELEM_TYP == 'ROD     ')) THEN

         TITLE_E( 1) = 'EndA Plane1 Moment'
         TITLE_E( 2) = 'EndB Plane1 Moment'
         TITLE_E( 3) = 'EndA Plane2 Moment'
         TITLE_E( 4) = 'EndB Plane2 Moment'

         TITLE_E( 5) = 'EndA Pl1 Shear Force'
         TITLE_E( 6) = 'EndB Pl1 Shear Force'
         TITLE_E( 7) = 'EndA Pl2 Shear Force'
         TITLE_E( 8) = 'EndB Pl2 Shear Force'

         TITLE_E( 9) = 'EndA Axial Force'
         TITLE_E(10) = 'EndB Axial Force'

         TITLE_E(11) = 'EndA Torque'
         TITLE_E(12) = 'EndB Torque'

         DO I=1,NUM_FEMAP_ROWS

            ELEM_VECS(I, 1) = FEMAP_EL_VECS(I,1)         ! M1a
            ELEM_VECS(I, 2) = FEMAP_EL_VECS(I,2)         ! M1b
            ELEM_VECS(I, 3) = FEMAP_EL_VECS(I,3)         ! M2a
            ELEM_VECS(I, 4) = FEMAP_EL_VECS(I,4)         ! M2b

            ELEM_VECS(I, 5) = FEMAP_EL_VECS(I,5)         ! V1a
            ELEM_VECS(I, 6) = FEMAP_EL_VECS(I,5)         ! V1b
            ELEM_VECS(I, 7) = FEMAP_EL_VECS(I,6)         ! V2a
            ELEM_VECS(I, 8) = FEMAP_EL_VECS(I,6)         ! V2b

            ELEM_VECS(I, 9) = FEMAP_EL_VECS(I,7)         ! Fa
            ELEM_VECS(I,10) = FEMAP_EL_VECS(I,7)         ! Fb

            ELEM_VECS(I,11) = FEMAP_EL_VECS(I,8)         ! Ta
            ELEM_VECS(I,12) = FEMAP_EL_VECS(I,8)         ! Tb

         ENDDO

         IF (ELEM_TYP == 'ROD     ') THEN

            DO J=9,12,2
               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
               DO I=1,NUM_FEMAP_ROWS
                  ELEM_VEC(I)  = ELEM_VECS(I,J)
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO
               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               ID(1) = VEC_ID
               ID(2) = VEC_ID + 1
               DO I=3,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               CALC_WARN  = '0'
               COMP_DIR   = '0'
               CENT_TOTAL = '1'
               WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)

               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID+1
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J+1)
               DO I=1,NUM_FEMAP_ROWS
                  ELEM_VEC(I)  = ELEM_VECS(I,J+1)
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO
               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               ID(1) = VEC_ID
               ID(2) = VEC_ID + 1
               DO I=3,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               CALC_WARN  = '0'
               COMP_DIR   = '0'
               CENT_TOTAL = '1'
               WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)
            ENDDO

         ELSE IF (ELEM_TYP == 'BAR     ') THEN

            DO J=1,12,2

               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
               DO I=1,NUM_FEMAP_ROWS
                  ELEM_VEC(I)  = ELEM_VECS(I,J)
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO
               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               ID(1) = VEC_ID
               ID(2) = VEC_ID + 1
               DO I=3,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               CALC_WARN  = '0'
               COMP_DIR   = '3'
               CENT_TOTAL = '1'
               WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)

               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID+1
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J+1)
               DO I=1,NUM_FEMAP_ROWS
                  ELEM_VEC(I)  = ELEM_VECS(I,J+1)
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO
               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               ID(1) = VEC_ID
               ID(2) = VEC_ID + 1
               DO I=3,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               CALC_WARN  = '0'
               COMP_DIR   = '3'
               CENT_TOTAL = '1'
               WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)

            ENDDO

         ENDIF

! Process plate elements

      ELSE IF ((ELEM_TYP == 'TRIA3K  ') .OR. (ELEM_TYP == 'TRIA3   ') .OR.                                                         &
               (ELEM_TYP == 'QUAD4K  ') .OR. (ELEM_TYP == 'QUAD4   ') .OR. (ELEM_TYP == 'SHEAR   ')) THEN

         TITLE_E( 1) = 'X  Membrane Force'
         TITLE_E( 2) = 'Y  Membrane Force'
         TITLE_E( 3) = 'XY Membrane Force'
         TITLE_E( 4) = 'X  Moment'
         TITLE_E( 5) = 'Y  Moment'
         TITLE_E( 6) = 'XY Moment'
         TITLE_E( 7) = 'X  TransShear Force'
         TITLE_E( 8) = 'Y  TransShear Force'

         DO J=1,8
            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            DO I=1,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            CALC_WARN  = '0'
            COMP_DIR   = '0'
            CENT_TOTAL = '1'
            WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)
         ENDDO

! Process ELAS elements

      ELSE IF ((ELEM_TYP == 'ELAS1   ') .OR. (ELEM_TYP == 'ELAS2   ') .OR.                                                         &
               (ELEM_TYP == 'ELAS3   ') .OR. (ELEM_TYP == 'ELAS4   ')) THEN

         TITLE_E( 1) = 'Spring Force'

         VEC_ID = VEC_ID_OFFSET + 1
         WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
         WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(1)
         DO I=1,NUM_FEMAP_ROWS
            ELEM_VEC(I)  = FEMAP_EL_VECS(I,1)
            ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
         ENDDO
         CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
         WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
         DO I=1,20
            ID(I) = 0
         ENDDO
         WRITE(NEU,1004) (ID(I),I= 1,10)
         WRITE(NEU,1004) (ID(I),I=11,20)
         WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
         CALC_WARN  = '0'
         COMP_DIR   = '0'
         CENT_TOTAL = '1'
         WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
         DO I=1,NUM_FEMAP_ROWS
            WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
         ENDDO
         WRITE(NEU,1008)

! Process BUSH elements

      ELSE IF (ELEM_TYP == 'BUSH    ') THEN

         TITLE_E( 1) = 'Force XE'
         TITLE_E( 2) = 'Force YE'
         TITLE_E( 3) = 'Force ZE'
         TITLE_E( 4) = 'Moment XE'
         TITLE_E( 5) = 'Moment YE'
         TITLE_E( 6) = 'Moment ZE'

         DO J=1,6
            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            DO I=1,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            CALC_WARN  = '0'
            COMP_DIR   = '0'
            CENT_TOTAL = '1'
            WRITE(NEU,1006) CALC_WARN, COMP_DIR, CENT_TOTAL
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)
         ENDDO

      ENDIF


      IF (ALLOCATED(ELEM_NUMS)) THEN
         DEALLOCATE ( ELEM_NUMS )
      ENDIF
      IF (ALLOCATED(ELEM_VECS)) THEN
         DEALLOCATE ( ELEM_VECS )
      ENDIF
      IF (ALLOCATED(ELEM_VEC)) THEN
         DEALLOCATE ( ELEM_VEC )
      ENDIF

      RETURN

! **********************************************************************************************************************************
  943 FORMAT(' *WARNING    : ELEMENT TYPE = "',A,'" FOR FEMAP ',A,' OUTPUT IN SUBROUTINE ',A,' HAS NOT BEEN PROGRAMMED')

 1001 FORMAT(2(I8,','),'       1,')

 1002 FORMAT(A,1X,A)

 1003 FORMAT(3(1ES17.6,','))

 1004 FORMAT(10(I8,','))

 1005 FORMAT(2(I8,','),2(7X,A,','))

 1006 FORMAT(3(7X,A,','))

 1007 FORMAT(I8,',',1ES17.6,',')

 1008 FORMAT('      -1,     0.          ,')

! **********************************************************************************************************************************
 
      END SUBROUTINE WRITE_FEMAP_ELFO_VECS

