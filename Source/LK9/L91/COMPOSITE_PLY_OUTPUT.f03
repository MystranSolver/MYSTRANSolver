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

   MODULE COMPOSITE_PLY_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_PLY_STRAINS, WRITE_PLY_STRESSES

   CONTAINS

      SUBROUTINE WRITE_PLY_STRAINS ( JSUB, NUM, IHDR, ETYPE, ITABLE  )

! Writes blocks of element ply strains for one subcase one element type for elements with PCOMP properties.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, OP2
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, BARTOR, INT_SC_NUM, LPCOMP_PLIES, NDOFR, NUM_CB_DOFS,            &
                                         NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, FTNAME, OGEL
      USE MODEL_STUF, ONLY            :  ANY_FAILURE_THEORY, ELEM_ONAME, LABEL, PCOMP, SCNUM, STITLE, TITLE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_OPT
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELEM_ONAME
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP
      USE MODAL_OUTPUT_WRITERS, ONLY  :  WRITE_SUBCASE_EIGENVEC_HEADER

      USE OP2_STRESS_OUTPUT, ONLY     :  GET_STRESS_CODE, WRITE_OES3_STATIC
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_PLY_STRAINS'
      CHARACTER(LEN=*), INTENT(IN)    :: IHDR              ! Indicator of whether to write an output header
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format
      CHARACTER(LEN=LEN(ELEM_ONAME))  :: ONAME             ! Element name to write out in F06 file

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      CHARACTER(8*BYTE), INTENT(IN)   :: ETYPE             ! the name of the element
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the op2 subtable number
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: ABS_ANS(10)       ! Max ABS for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MAX_ANS(10)       ! Max for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MIN_ANS(10)       ! Min for all grids output for each of the 6 disp components

      ! table -3 info
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      INTEGER(LONG)                   :: STRESS_CODE            ! flag for op2
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: STITLEI                ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: FIELD5_INT_MODE
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE

      ! op2 specific flags
      INTEGER(LONG)                   :: ISUBCASE_INDEX         ! the index into SCNUM
      INTEGER(LONG)                   :: DEVICE_CODE            ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)                   :: NUM_WIDE               ! the number of "words" for an element
      INTEGER(LONG)                   :: NVALUES                ! the number of "words" for all the elments
      INTEGER(LONG)                   :: NTOTAL                 ! the number of bytes for all NVALUES
      INTEGER(LONG)                   :: ISUBCASE               ! the subcase ID

      INTRINSIC                       :: MAX, MIN, DABS

! **********************************************************************************************************************************
      ! initialize
      ANALYSIS_CODE = -1
      ELEMENT_TYPE = -1
      DEVICE_CODE = 1



! **********************************************************************************************************************************
      FILL(1:) = ' '

      DO J=1,10
         ABS_ANS(J) = ZERO
         MAX_ANS(J) = -MACH_LARGE_NUM
         MIN_ANS(J) = ZERO
      ENDDO

! Get element output name

      ONAME(1:) = ' '
      CALL GET_ELEM_ONAME ( ONAME )

      ! Write output headers if this is not the first use of this subr.
      ANALYSIS_CODE = -1
      FIELD5_INT_MODE = 0
      FIELD6_EIGENVALUE = 0.0

      IF (IHDR == 'Y') THEN

         ! -- F06 header: OUTPUT FOR SUBCASE, EIGENVECTOR or CRAIG-BAMPTON DOF
         CALL WRITE_SUBCASE_EIGENVEC_HEADER(JSUB, .TRUE.)
         ISUBCASE_INDEX = 0
         IF      (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JSUB
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
         ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 10
            FIELD5_INT_MODE = SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            ISUBCASE_INDEX = 2
            ANALYSIS_CODE = 7
            FIELD5_INT_MODE = JSUB

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN   ! Write info on what CB DOF the output is for
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB

            IF ((JSUB <= NDOFR) .OR. (JSUB >= NDOFR+NVEC)) THEN
               IF (JSUB <= NDOFR) THEN
                  BDY_DOF_NUM = JSUB
               ELSE
                  BDY_DOF_NUM = JSUB-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

            IF       (JSUB <= NDOFR) THEN
               WRITE(F06,103) JSUB, NUM_CB_DOFS, 'acceleration', BDY_GRID, BDY_COMP
            ELSE IF ((JSUB > NDOFR) .AND. (JSUB <= NDOFR+NVEC)) THEN
               WRITE(F06,104) JSUB, NUM_CB_DOFS, JSUB-NDOFR
            ELSE
               WRITE(F06,103) JSUB, NUM_CB_DOFS, 'displacement', BDY_GRID, BDY_COMP
            ENDIF
         ENDIF
         ISUBCASE = SCNUM(ISUBCASE_INDEX)

         TITLEI = TITLE(INT_SC_NUM)
         STITLEI = STITLE(INT_SC_NUM)
         LABELI = LABEL(INT_SC_NUM)

         ! -- F06 1st 2 header lines for strain output description
         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            WRITE(F06,302) FILL(1: 0)
         ELSE
            WRITE(F06,301)  FILL(1: 0)
         ENDIF
         IF (ANY_FAILURE_THEORY == 'N') THEN
            IF (STRN_OPT == 'VONMISES') THEN
               WRITE(F06,1401) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
            ELSE
               WRITE(F06,1402) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
            ENDIF
         ELSE
            IF (STRN_OPT == 'VONMISES') THEN
               WRITE(F06,1403) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
            ELSE
               WRITE(F06,1404) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
            ENDIF
         ENDIF

      ENDIF

! Write the element strain output
      ! 95:  cquad4
      ! 96:  cquad8
      ! 97:  ctria3
      ! 98:  ctria6
      ! 232: cquadr
      ! 233: ctriar
      IF      (ETYPE == "QUAD4   ") THEN
          ELEMENT_TYPE = 95
      ELSE IF (ETYPE == "TRIA3   ") THEN
          ELEMENT_TYPE = 97
      ENDIF

      ! op2 output
      IF (STRN_OPT == 'VONMISES') THEN
         ! von Mises
         !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
         CALL GET_STRESS_CODE( STRESS_CODE, 1,            1,         1)
      ELSE
         ! max shear
         CALL GET_STRESS_CODE( STRESS_CODE, 0,            1,         1)
      ENDIF

      ! elementi = [eid_device, layer,
      !             o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
      NUM_WIDE = 11
      NVALUES = NUM * NUM_WIDE
      CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                             TITLEI, STITLEI, LABELI, &
                             FIELD5_INT_MODE, FIELD6_EIGENVALUE)

      WRITE(OP2) NVALUES
      ! TODO: dropping the failure theory
      WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,2), (REAL(OGEL(I,J),4),J=1,9), I=1,NUM)

      ! F06 output
      DO I=1,NUM

         IF (EID_OUT_ARRAY(I,2) == 1) THEN
            WRITE(F06,*)
            IF (ANY_FAILURE_THEORY == 'Y') THEN
               IF (FTNAME(I) /= 'None') THEN
                  WRITE(F06,1405) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,10), FTNAME(I)
               ELSE
                  WRITE(F06,1406) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9), FTNAME(I)
               ENDIF
            ELSE
               WRITE(F06,1406) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
            ENDIF
         ELSE
            IF (ANY_FAILURE_THEORY == 'Y') THEN
               IF (FTNAME(I) /= 'None') THEN
                  WRITE(F06,1407) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,10)
               ELSE
                  WRITE(F06,1408) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
               ENDIF
            ELSE
               WRITE(F06,1408) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
            ENDIF
         ENDIF

      ENDDO

      DO I=1,NUM

         DO J=1,10
            IF (OGEL(I,J) > MAX_ANS(J)) THEN
               MAX_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO

      ENDDO

      DO J=1,10
         MIN_ANS(J) = MAX_ANS(J)
      ENDDO

      DO I=1,NUM

         DO J=1,10
            IF (OGEL(I,J) < MIN_ANS(J)) THEN
               MIN_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO

      ENDDO

      DO I=1,10                                            ! Get absolute max strains
         ABS_ANS(I) = MAX( DABS(MAX_ANS(I)), DABS(MIN_ANS(I)) )
      ENDDO

      WRITE(F06,1409) FILL(1: 0), FILL(1 :0), (MAX_ANS(I),I=1,10),                                                                 &
                      FILL(1: 0)            , (MIN_ANS(I),I=1,10),                                                                 &
                      FILL(1: 0)            , (ABS_ANS(I),I=1,10), FILL(1: 0)




      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

  102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

  103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

  104 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

  201 FORMAT(1X,A)

  301 FORMAT(45X,A,' S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S')

  302 FORMAT(9X,A,'C B   E L E M E N T   S T R A I N   O T M   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M')

 1401 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Strains in fiber and matrix directions     Inter-laminar strains      Principal strains  (zero shear)',&
           '      von',/,                                                                                                          &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Mises')

 1402 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Strains in fiber and matrix directions     Inter-laminar strains      Principal strains  (zero shear)',&
           '      Max',/,                                                                                                          &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Shear')

 1403 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Strains in fiber and matrix directions     Inter-laminar strains      Principal strains  (zero shear)',&
           '      von          Failure',/,                                                                                         &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Mises      Index   Theory')

 1404 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Strains in fiber and matrix directions     Inter-laminar strains      Principal strains  (zero shear)',&
           '      Max          Failure',/,                                                                                         &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Shear      Index   Theory')

!5001 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),12X,A)
 1406 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),12X,A)

!5002 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5))
 1408 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5))

!5003 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X,A)
 1405 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X,A)

!5004 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X)
 1407 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X)

 1409 FORMAT(14X,A,'  ------------ ------------ ------------    ------------  ------------  ------- ------------ ------------',    &
                 ' ------------ ---------',/,                                                                                      &
             1X,A,'MAX* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,/,                                              &
             1X,A,'MIN* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,//,                                             &
             1X,A,'ABS* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,/,                                              &
             1X,A,'*for output set')



! **********************************************************************************************************************************

      END SUBROUTINE WRITE_PLY_STRAINS


      SUBROUTINE WRITE_PLY_STRESSES ( JSUB, NUM, IHDR, ETYPE, ITABLE )

      ! Writes blocks of element ply stresses for one subcase one
      ! element type for elements with PCOMP properties.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, OP2
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, BARTOR, INT_SC_NUM, LPCOMP_PLIES, NDOFR, NUM_CB_DOFS,            &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, FTNAME, OGEL
      USE MODEL_STUF, ONLY            :  ANY_FAILURE_THEORY, ELEM_ONAME, LABEL, PCOMP, SCNUM, STITLE, TITLE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_OPT
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELEM_ONAME
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP
      USE MODAL_OUTPUT_WRITERS, ONLY  :  WRITE_SUBCASE_EIGENVEC_HEADER

      USE OP2_STRESS_OUTPUT, ONLY     :  GET_STRESS_CODE, WRITE_OES3_STATIC
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_PLY_STRESSES'
      CHARACTER(LEN=*), INTENT(IN)    :: IHDR              ! Indicator of whether to write an output header
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format
      CHARACTER(LEN=LEN(ELEM_ONAME))  :: ONAME             ! Element name to write out in F06 file

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      CHARACTER(8*BYTE), INTENT(IN)   :: ETYPE             ! the name of the element
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the op2 subtable number
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: ABS_ANS(10)       ! Max ABS for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MAX_ANS(10)       ! Max for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MIN_ANS(10)       ! Min for all grids output for each of the 6 disp components

      ! table -3 info
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      INTEGER(LONG)                   :: STRESS_CODE            ! flag for op2
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: STITLEI                ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: FIELD5_INT_MODE
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE

      ! op2 specific flags
      INTEGER(LONG)                   :: ISUBCASE_INDEX         ! the index into SCNUM
      INTEGER(LONG)                   :: DEVICE_CODE            ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)                   :: NUM_WIDE               ! the number of "words" for an element
      INTEGER(LONG)                   :: NVALUES                ! the number of "words" for all the elments
      INTEGER(LONG)                   :: NTOTAL                 ! the number of bytes for all NVALUES
      INTEGER(LONG)                   :: ISUBCASE               ! the subcase ID
      LOGICAL                         :: DEBUG_OP2              ! flag

      INTRINSIC                       :: MAX, MIN, DABS

! **********************************************************************************************************************************
      ! initialize
      ANALYSIS_CODE = -1
      ELEMENT_TYPE = -1
      DEVICE_CODE = 1
      DEBUG_OP2 = .FALSE.



! **********************************************************************************************************************************
      FILL(1:) = ' '

      DO J=1,10
         ABS_ANS(J) = ZERO
         MAX_ANS(J) = -MACH_LARGE_NUM
         MIN_ANS(J) = ZERO
      ENDDO

! Get element output name

      ONAME(1:) = ' '
 90   FORMAT("*DEBUG:   WRITE_PLY_STRESSES:  ELEM_ONAME=",A," ONAME=",A," NUM=",I8)
      CALL GET_ELEM_ONAME ( ONAME )
      IF (DEBUG_OP2)  WRITE(ERR,90) ETYPE, ONAME, NUM

      ! Write output headers if this is not the first use of this subr.
      ANALYSIS_CODE = -1
      FIELD5_INT_MODE = 0
      FIELD6_EIGENVALUE = 0.0

      IF (IHDR == 'Y') THEN

         ! -- F06 header: OUTPUT FOR SUBCASE, EIGENVECTOR or CRAIG-BAMPTON DOF
         CALL WRITE_SUBCASE_EIGENVEC_HEADER(JSUB, .TRUE.)
         ISUBCASE_INDEX = 0
         IF      (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JSUB
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
         ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 10
            FIELD5_INT_MODE = SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            ISUBCASE_INDEX = 2
            ANALYSIS_CODE = 7
            FIELD5_INT_MODE = JSUB

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN   ! Write info on what CB DOF the output is for
            ISUBCASE_INDEX = 1
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB

            IF ((JSUB <= NDOFR) .OR. (JSUB >= NDOFR+NUM_CB_DOFS)) THEN
               IF (JSUB <= NDOFR) THEN
                  BDY_DOF_NUM = JSUB
               ELSE
                  BDY_DOF_NUM = JSUB-(NDOFR+NUM_CB_DOFS)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

            IF       (JSUB <= NDOFR) THEN
               WRITE(F06,103) JSUB, NUM_CB_DOFS, 'acceleration', BDY_GRID, BDY_COMP
            ELSE IF ((JSUB > NDOFR) .AND. (JSUB <= NDOFR+NUM_CB_DOFS)) THEN
               WRITE(F06,104) JSUB, NUM_CB_DOFS, JSUB-NDOFR
            ELSE
               WRITE(F06,103) JSUB, NUM_CB_DOFS, 'displacement', BDY_GRID, BDY_COMP
            ENDIF
         ENDIF
         ISUBCASE = SCNUM(ISUBCASE_INDEX)

         TITLEI  = TITLE(INT_SC_NUM)
         STITLEI = STITLE(INT_SC_NUM)
         LABELI  = LABEL(INT_SC_NUM)

         ! -- F06 1st 2 header lines for stress output description
         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            WRITE(F06,302) FILL(1: 0)
         ELSE
            WRITE(F06,301)  FILL(1: 0)
         ENDIF
         WRITE(F06,*) 'ANY_FAILURE_THEORY=',ANY_FAILURE_THEORY
         WRITE(F06,*) 'STRE_OPT=',STRE_OPT
         IF (ANY_FAILURE_THEORY == 'N') THEN
            IF (STRE_OPT == 'VONMISES') THEN
               WRITE(F06,*) 'A'
               WRITE(F06,1401) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
               WRITE(F06,1499) FILL(1: 0), FILL(1: 0), FILL(1: 0), FILL(1: 0)
            ELSE
               WRITE(F06,*) 'B'
               WRITE(F06,1402) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
               WRITE(F06,1499) FILL(1: 0), FILL(1: 0), FILL(1: 0), FILL(1: 0)
            ENDIF
         ELSE
            IF (STRE_OPT == 'VONMISES') THEN
               WRITE(F06,*) 'C'
               WRITE(F06,1403) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
               WRITE(F06,1499) FILL(1: 0), FILL(1: 0), FILL(1: 0), FILL(1: 0)
            ELSE
               WRITE(F06,*) 'D'
               WRITE(F06,1404) FILL(1: 0), ONAME, FILL(1: 0), FILL(1: 0)
               WRITE(F06,1499) FILL(1: 0), FILL(1: 0), FILL(1: 0), FILL(1: 0)
            ENDIF
         ENDIF

      ENDIF

! Write the element stress output
      ! 95:  cquad4
      ! 96:  cquad8
      ! 97:  ctria3
      ! 98:  ctria6
      ! 232: cquadr
      ! 233: ctriar
      IF      (ETYPE == "QUAD4   ") THEN
          ELEMENT_TYPE = 95
      ELSE IF (ETYPE == "TRIA3   ") THEN
          ELEMENT_TYPE = 97
      ENDIF

      ! op2 output
      IF (STRE_OPT == 'VONMISES') THEN
         ! von Mises
         !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
         CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
      ELSE
         ! max shear
         CALL GET_STRESS_CODE( STRESS_CODE, 0,            0,         1)
      ENDIF

      ! elementi = [eid_device, layer,
      !             o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
      NUM_WIDE = 11
      NVALUES = NUM * NUM_WIDE
      CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                             TITLEI, STITLEI, LABELI, &
                             FIELD5_INT_MODE, FIELD6_EIGENVALUE)

      WRITE(OP2) NVALUES
      ! TODO: dropping the failure theory
      WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,2), (REAL(OGEL(I,J),4),J=1,9), I=1,NUM)

      ! F06 output
      DO I=1,NUM

         IF (EID_OUT_ARRAY(I,2) == 1) THEN
            WRITE(F06,*)
            IF (ANY_FAILURE_THEORY == 'Y') THEN
               IF (FTNAME(I) /= 'None') THEN
                  WRITE(F06,1405) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,10), FTNAME(I)
               ELSE
                  WRITE(F06,1406) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9), FTNAME(I)
               ENDIF
            ELSE
               WRITE(F06,1406) FILL(1: 0), EID_OUT_ARRAY(I,1), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
            ENDIF
         ELSE
            IF (ANY_FAILURE_THEORY == 'Y') THEN
               IF (FTNAME(I) /= 'None') THEN
                  WRITE(F06,1407) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,10)
               ELSE
                  WRITE(F06,1408) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
               ENDIF
            ELSE
               WRITE(F06,1408) FILL(1: 0), EID_OUT_ARRAY(I,2), (OGEL(I,J),J=1,9)
            ENDIF
         ENDIF

      ENDDO

      DO I=1,NUM

         DO J=1,10
            IF (OGEL(I,J) > MAX_ANS(J)) THEN
               MAX_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO

      ENDDO

      DO J=1,10
         MIN_ANS(J) = MAX_ANS(J)
      ENDDO

      DO I=1,NUM

         DO J=1,10
            IF (OGEL(I,J) < MIN_ANS(J)) THEN
               MIN_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO

      ENDDO

      DO I=1,10                                          ! Get absolute max stresses
         ABS_ANS(I) = MAX( DABS(MAX_ANS(I)), DABS(MIN_ANS(I)) )
      ENDDO

      WRITE(F06,1409) FILL(1: 0), FILL(1 :0), (MAX_ANS(I),I=1,10),                                                                 &
                      FILL(1: 0)            , (MIN_ANS(I),I=1,10),                                                                 &
                      FILL(1: 0)            , (ABS_ANS(I),I=1,10), FILL(1: 0)




      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

  102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

  103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

  104 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

  201 FORMAT(1X,A)

  301 FORMAT(45X,A,'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S')

  302 FORMAT(9X,A,'C B   E L E M E N T   S T R E S S   O T M   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M')

 1401 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Stresses in fiber and matrix directions    Inter-laminar stresses     Principal stresses (zero shear)',&
           '      von',/,                                                                                                          &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Mises')

 1402 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Stresses in fiber and matrix directions    Inter-laminar stresses     Principal stresses (zero shear)',&
           '      Max',/,                                                                                                          &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Shear')

 1403 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Stresses in fiber and matrix directions    Inter-laminar stresses     Principal stresses (zero shear)',&
           '      von          Failure',/,                                                                                         &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Mises      Index   Theory')

 1404 FORMAT(53X,A,'F O R   E L E M E N T   T Y P E   ',A11,/,                                                                     &
      1X,A,' Element   Ply  Stresses in fiber and matrix directions    Inter-laminar stresses     Principal stresses (zero shear)',&
           '      Max          Failure',/,                                                                                         &
      1X,A,'    ID     Num     Normal-1     Normal-2     Shear-12      Shear-13      Shear-23      Angle     Major        Minor',  &
           '       Shear      Index   Theory')

 1499 FORMAT(15X,                                                                                                                  &
 '--------------------------------------------------------------------------------------------------------------------' ,/,14X,A,  &
'| N O T E:  TRANSVERSE SHEAR STRESSES IN ALL PLIES ARE BASED ON THE TRANSVERSE SHEAR STRAINS AND PARAMETER PCMPTSTM. |',/,14X,A,  &
'|                     THIS DOES NOT TAKE INTO CONSIDERATION ANY THEORETICAL (e.g. PARABOLIC) DISTRIBUTION            |',/,14X,A,  &
'|                                       THAT MAY EXIST FOR THE TRANSVERSE SHEAR STRESSES                             |',/,15X,A,  &
 '--------------------------------------------------------------------------------------------------------------------')

!5001 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),12X,A)
 1406 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),12X,A)

!5002 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5))
 1408 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5))

!5003 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X,A)
 1405 FORMAT(1X,A,I8,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X,A)

!5004 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X)
 1407 FORMAT(1X,A,8X,I6,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,2X)

 1409 FORMAT(14X,A,'  ------------ ------------ ------------    ------------  ------------  ------- ------------ ------------',    &
                 ' ------------ ---------',/,                                                                                      &
             1X,A,'MAX* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,/,                                              &
             1X,A,'MIN* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,//,                                             &
             1X,A,'ABS* :  ',6X,3(1ES13.5),2X,2(1ES14.5),0PF9.3,3(1ES13.5),1ES10.2,/,                                              &
             1X,A,'*for output set')


! **********************************************************************************************************************************

      END SUBROUTINE WRITE_PLY_STRESSES

   END MODULE COMPOSITE_PLY_OUTPUT
