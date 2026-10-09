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

   MODULE GRID_OUTPUT_WRITERS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_GRD_PRT_OUTPUTS, WRITE_GRD_OP2_OUTPUTS, WRITE_FEMAP_GRID_VECS

   CONTAINS

      SUBROUTINE WRITE_GRD_PRT_OUTPUTS ( JVEC, NUM, WHAT, IHDR, ALL_SAME_CID, WRITE_OGEL )

! Writes printed output for grid point related quantities (accels, displacements, eigenvectors, applied loads and SPC, MPC forces)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, INT_SC_NUM, MELGP, MOGEL, NDOFR, NVEC, NUM_CB_DOFS,              &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, MAXREQ, OGEL
      USE MODEL_STUF, ONLY            :  SCNUM, SUBLOD
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP
      USE MODAL_OUTPUT_WRITERS, ONLY :  WRITE_SUBCASE_EIGENVEC_HEADER
      USE TEXT_FIELD_UTILS, ONLY      :  FMT_ES14_6, FMT_I8_RJ

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_GRD_PRT_OUTPUTS'
      CHARACTER(LEN=*) , INTENT(IN)   :: IHDR              ! Indicator of whether to write an output header in this use of this subr
      CHARACTER(LEN=*) , INTENT(IN)   :: WHAT              ! Indicator whether to process displ or force output requests
      CHARACTER(1*BYTE), INTENT(IN)   :: ALL_SAME_CID      ! Indicator of whether all grids, for the output set, have the same
!                                                            global coord sys
      CHARACTER(14*BYTE)              :: OGEL_CHAR(MOGEL)  ! Char representation of 1 row of OGEL outputs

      CHARACTER( 1*BYTE)              :: PRINT_TOTALS      ! This will be set to 'Y' if DEBUG(92) > 0 so OLOAD, SPCF, MPCF force
!                                                            totals will be printed even if ALL_SAME_CID = 'N'
      CHARACTER(14*BYTE)              :: ABS_ANS_CHAR(6)   ! Character variable that contains the 6 grid abs  outputs
      CHARACTER(14*BYTE)              :: MAX_ANS_CHAR(6)   ! Character variable that contains the 6 grid max  outputs
      CHARACTER(14*BYTE)              :: MIN_ANS_CHAR(6)   ! Character variable that contains the 6 grid min  outputs
      CHARACTER(14*BYTE)              :: TOTALS_CHAR(6)    ! Character variable that contains the 6 grid tot  outputs
      CHARACTER(108*BYTE)             :: LINE_BUF          ! Pre-assembled per-grid output line (matches FORMAT 9902 layout)

      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Sol'n vector num. Can be internal subcase number or eigenvector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: LINES_WRITTEN     ! Number of lines written for the grids


      REAL(DOUBLE)                    :: ABS_ANS(6)        ! Max Abs for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MAX_ANS(6)        ! Max for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: MIN_ANS(6)        ! Min for all grids output for each of the 6 disp components
      REAL(DOUBLE)                    :: TOTALS(6)         ! Totals of each of the 6 components output

      CHARACTER(1*BYTE), INTENT(IN)   :: WRITE_OGEL(NUM)   ! 'Y'/'N' as to whether to write OGEL for a grid (used to avoid writing
!                                                            constr forces in subr this subr for grids that have no constr force
      INTRINSIC                       :: MAX, MIN, DABS



! **********************************************************************************************************************************
!  Make sure that WHAT is a valid value

      IF ((WHAT == 'ACCE') .OR. (WHAT == 'DISP') .OR. (WHAT == 'OLOAD') .OR. (WHAT == 'SPCF') .OR. (WHAT == 'MPCF')) THEN
         CONTINUE
      ELSE
         WRITE(ERR,9100) WHAT
         WRITE(F06,9100) WHAT
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Write output headers.

      IF (IHDR == 'Y') THEN

         CALL WRITE_SUBCASE_EIGENVEC_HEADER(JVEC, .TRUE.)
         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN   ! Write info on what CB DOF the output is for

            IF ((JVEC <= NDOFR) .OR. (JVEC >= NDOFR+NVEC)) THEN
               IF (JVEC <= NDOFR) THEN
                  BDY_DOF_NUM = JVEC
               ELSE
                  BDY_DOF_NUM = JVEC-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

               IF       (JVEC <= NDOFR) THEN
                  WRITE(F06,9013) JVEC, NUM_CB_DOFS, 'acceleration', BDY_GRID, BDY_COMP
               ELSE IF ((JVEC > NDOFR) .AND. (JVEC <= NDOFR+NVEC)) THEN
                  WRITE(F06,9015) JVEC, NUM_CB_DOFS, JVEC-NDOFR
               ELSE
                  WRITE(F06,9013) JVEC, NUM_CB_DOFS, 'displacement', BDY_GRID, BDY_COMP
               ENDIF
         ENDIF

         IF (WHAT == 'ACCE') THEN
            IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               WRITE(F06,9314)
            ENDIF

         ELSE IF (WHAT == 'DISP') THEN
            IF    ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
               WRITE(F06,9322)
            ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
               WRITE(F06,9322)
            ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               WRITE(F06,9323)
            ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
               WRITE(F06,9323)
            ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               WRITE(F06,9324)
            ENDIF

         ELSE IF (WHAT == 'OLOAD') THEN
            WRITE(F06,9331)
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
               WRITE(F06,9332)
            ENDIF

         ELSE IF (WHAT == 'SPCF') THEN
            IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               WRITE(F06,9342)
            ELSE
               WRITE(F06,9341)
            ENDIF

         ELSE IF (WHAT == 'MPCF') THEN
            IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               WRITE(F06,9352)
            ELSE
               WRITE(F06,9351)
            ENDIF
         ENDIF
         WRITE(F06,9501)


      ENDIF

! Get MAX, MIN, ABS values

      DO J=1,6
         MAX_ANS(J) = -MACH_LARGE_NUM
      ENDDO

      DO I=1,NUM
         DO J=1,6
            IF (OGEL(I,J) > MAX_ANS(J)) THEN
               MAX_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO
      ENDDO

      DO J=1,6
         MIN_ANS(J) = MAX_ANS(J)
      ENDDO

      DO I=1,NUM
         DO J=1,6
            IF (OGEL(I,J) < MIN_ANS(J)) THEN
               MIN_ANS(J) = OGEL(I,J)
            ENDIF
         ENDDO
      ENDDO

      DO I=1,6
         ABS_ANS(I) = MAX( DABS(MAX_ANS(I)), DABS(MIN_ANS(I)) )
      ENDDO

      DO I=1,6
         IF (ABS(ABS_ANS(I)) == ZERO) THEN
            ABS_ANS_CHAR(I) = '  0.0         '
         ELSE
            CALL FMT_ES14_6 ( ABS_ANS(I), ABS_ANS_CHAR(I) )
         ENDIF
         IF (ABS(MAX_ANS(I)) == ZERO) THEN
            MAX_ANS_CHAR(I) = '  0.0         '
         ELSE
            CALL FMT_ES14_6 ( MAX_ANS(I), MAX_ANS_CHAR(I) )
         ENDIF
         IF (ABS(MIN_ANS(I)) == ZERO) THEN
            MIN_ANS_CHAR(I) = '  0.0         '
         ELSE
            CALL FMT_ES14_6 ( MIN_ANS(I), MIN_ANS_CHAR(I) )
         ENDIF
      ENDDO

! Write accels, displ's, applied forces or SPC forces (also calc TOTALS for forces if that is being output)
! TOTALS(J) is summation of G.P. values of applied forces, SPC forces, or MFC forces, for each of the J=1,6 components.

      TOTALS = ZERO

! Pre-fill the fixed-whitespace positions of LINE_BUF that match FORMAT 9902 = (6X,2(1X,I8),6A).
! Variable fields (GID, COORD, 6x14-char values) are overwritten per iteration in the loop below.
      LINE_BUF = ' '

      LINES_WRITTEN = 0
      DO I=1,NUM

         IF ((WHAT == 'OLOAD') .OR. (WHAT == 'SPCF') .OR. (WHAT == 'MPCF')) THEN
            DO J=1,6
               TOTALS(J) = TOTALS(J) + OGEL(I,J)
               IF (ABS(TOTALS(J)) == ZERO) THEN
                  TOTALS_CHAR(J) = '  0.0         '
               ELSE
                  CALL FMT_ES14_6 ( TOTALS(J), TOTALS_CHAR(J) )
               ENDIF
            ENDDO
         ENDIF

         IF (WRITE_OGEL(I) == 'Y') THEN

!           Assemble the per-grid output line directly into LINE_BUF, then emit with a single A-format WRITE.
!           Layout (matches FORMAT 9902): 6X | 1X | I8 | 1X | I8 | 6 * A14  ==> 108 chars total.
            CALL FMT_I8_RJ ( GID_OUT_ARRAY(I,1), LINE_BUF( 8:15) )
            CALL FMT_I8_RJ ( GID_OUT_ARRAY(I,2), LINE_BUF(17:24) )
            DO J=1,6
               IF (ABS(OGEL(I,J)) == ZERO) THEN
                  LINE_BUF(25 + (J-1)*14 : 24 + J*14) = '  0.0         '
               ELSE
                  CALL FMT_ES14_6 ( OGEL(I,J), LINE_BUF(25 + (J-1)*14 : 24 + J*14) )
               ENDIF
            ENDDO
            WRITE(F06,'(A)') LINE_BUF

            IF (GID_OUT_ARRAY(I,MELGP+1) > 0) THEN
               DO J=1,GID_OUT_ARRAY(I,MELGP+1)
                  WRITE(F06,*)
               ENDDO
            ENDIF

            LINES_WRITTEN = LINES_WRITTEN + 1

         ENDIF

      ENDDO

      IF (LINES_WRITTEN > 2) THEN
         WRITE(F06,9601) (MAX_ANS_CHAR(J),J=1,6), (MIN_ANS_CHAR(J),J=1,6), (ABS_ANS_CHAR(J),J=1,6)
      ENDIF

      IF (DEBUG(92) == 0) THEN
         PRINT_TOTALS = ALL_SAME_CID
      ELSE
         PRINT_TOTALS = 'Y'
      ENDIF

      IF (LINES_WRITTEN > 1) THEN
         IF (PRINT_TOTALS == 'Y') THEN
            IF (WHAT == 'OLOAD') THEN
               WRITE(F06,9701) (TOTALS_CHAR(J),J=1,6)
            ELSE IF (WHAT == 'SPCF' ) THEN
               WRITE(F06,9702) (TOTALS_CHAR(J),J=1,6)
            ELSE IF (WHAT == 'MPCF' ) THEN
               WRITE(F06,9703) (TOTALS_CHAR(J),J=1,6)
            ENDIF
         ELSE
            IF (WHAT == 'OLOAD') THEN
               WRITE(F06,9711)
            ELSE IF (WHAT == 'SPCF' ) THEN
               WRITE(F06,9712)
            ELSE IF (WHAT == 'MPCF' ) THEN
               WRITE(F06,9713)
            ENDIF
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9009 FORMAT(1X,A)

 9011 FORMAT(' OUTPUT FOR SUBCASE ',I8)

 9012 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

 9013 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

 9015 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

 9100 FORMAT(' *ERROR  9100: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ILLEGAL INPUT FOR VARIABLE "WHAT" = ',A)

 9314 FORMAT(1X,'                                                C B   A C C E L E R A T I O N   O T M',/,                         &
             1X,'                                             (in global coordinate system at each grid)')

 9322 FORMAT(1X,'                                                      D I S P L A C E M E N T S',/,                               &
             1X,'                                             (in global coordinate system at each grid)')

 9323 FORMAT(1X,'                                                        E I G E N V E C T O R',/,                                 &
             1X,'                                             (in global coordinate system at each grid)')

 9324 FORMAT(1X,'                                                C B   D I S P L A C E M E N T   O T M',/,                         &
             1X,'                                             (in global coordinate system at each grid)')

 9331 FORMAT(1X,'                                                    A P P L I E D    F O R C E S',/,                              &
             1X,'                                             (in global coordinate system at each grid)')

 9332 FORMAT(1X,'                                                (including equivalent thermal loads)')

 9341 FORMAT(1X,'                                                         S P C   F O R C E S',/,                                  &
             1X,'                                             (in global coordinate system at each grid)')

 9342 FORMAT(1X,'                                                   C B   S P C   F O R C E   O T M',/,                            &
             1X,'                                             (in global coordinate system at each grid)')

 9351 FORMAT(1X,'                                                         M P C   F O R C E S',/,                                  &
             1X,'                                             (in global coordinate system at each grid)')

 9352 FORMAT(1X,'                                                   C B   M P C   F O R C E   O T M',/,                            &
             1X,'                                             (in global coordinate system at each grid)')

 9501 FORMAT(11X,'GRID     COORD      T1            T2            T3            R1            R2            R3',/,                 &
             11X,'          SYS')

 9601 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             16X,'MAX* :  ',6A14,/,                                                                                                &
             16X,'MIN* :  ',6A14,//,                                                                                               &
             16X,'ABS* :  ',6A14,/,                                                                                                &
             16X,'*for output set')

 9701 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'APPLIED FORCE TOTALS:  ',6A14,/,3X,'(for output set)')

 9702 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'    SPC FORCE TOTALS:  ',6A14,/,5X,'(for output set)')

 9703 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'    MPC FORCE TOTALS:  ',6A14,/,5X,'(for output set)')

 9711 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'APPLIED FORCE TOTALS: not printed since all grids do not have the same global coordinate system')

 9712 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'    SPC FORCE TOTALS: not printed since all grids do not have the same global coordinate system')

 9713 FORMAT(11X,'              ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             1X,'    MPC FORCE TOTALS: not printed since all grids do not have the same global coordinate system')

 9902 FORMAT(6X,2(1X,I8),6A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_GRD_PRT_OUTPUTS


      SUBROUTINE WRITE_GRD_OP2_OUTPUTS ( JSUB, NUM, WHAT, ITABLE, NEW_RESULT )
!      Writes "plot" output for grid point related quantities:
!        - accels
!        - displacements
!        - eigenvectors
!        - applied loads
!        - SPC / MPC forces
!        - velocity????

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, INT_SC_NUM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, OGEL
      USE MODEL_STUF, ONLY            :  GRID, LABEL, SCNUM, SUBLOD, STITLE, TITLE
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL

      USE FILE_LIFECYCLE, ONLY        :  OUTA_HERE
      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE, WRITE_ITABLE, WRITE_TABLE_HEADER
      USE OP2_GRID_OUTPUT, ONLY       :  WRITE_OUG3_EIGN, WRITE_OUG3_STATIC
      IMPLICIT NONE
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)) :: SUBR_NAME = 'WRITE_GRD_OP2_OUTPUTS'
      CHARACTER(LEN=*), INTENT(IN)     :: WHAT               ! Indicator whether to process displ or
                                                            ! force output requests
!     CHARACTER(LEN=1)                 :: G_OR_S            ! 'G' if a grid point or 'S' if a scalar point
      CHARACTER(LEN=8)                 :: TABLE_NAME        ! Name of the op2 table that we're writing
      CHARACTER(LEN=128)                :: TITLEI            ! Solution title
      CHARACTER(LEN=128)                :: STITLEI           ! Subcase subtitle
      CHARACTER(LEN=128)                :: LABELI            ! Subcase label

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the OP2 subtable
      LOGICAL,       INTENT(INOUT)    :: NEW_RESULT        ! Is this a new result?
      INTEGER(LONG)                   :: I,J               ! DO loop indices
!      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = WRITE_GRD_OP2_OUTPUTS_BEGEND
      INTEGER(LONG)                   :: ISUBCASE          ! the current subcase ID
      INTEGER(LONG)                   :: TABLE_CODE        ! flag for the type of table
      INTEGER(LONG)                   :: ANALYSIS_CODE     ! flag for the solution type
      INTEGER(LONG), DIMENSION(NUM)   :: G_OR_S            ! flag for the type of point
      INTEGER(LONG)                   :: DEVICE_CODE       ! flag for PLOT,PRINT,PUNCH
      INTEGER(LONG)                   :: MODE              ! mode number for an eigenvector solution
      REAL(DOUBLE)                    :: EIGENVALUE        ! the eigenvalue for an eigenvector solution
      INTEGER(LONG)                   :: THERMAL_FLAG      ! flag for a temperature result
      INTEGER(LONG)                   :: NTOTAL            ! the number of total bytes for all the "words"
      INTEGER(LONG)                   :: NUM_WIDE          ! the width in bytes of a result
      INTEGER(LONG)                   :: NVALUES           ! the width in "words" of a result
      INTEGER(LONG)                   :: ISUBCASE_INDEX    ! the index into SCNUM

! **********************************************************************************************************************************
      ! TODO: assuming PLOT
      DEVICE_CODE = 1


! **********************************************************************************************************************************
      ! Make sure that WHAT is a valid value
      IF ((WHAT == 'ACCE') .OR. (WHAT == 'DISP') .OR. (WHAT == 'OLOAD') .OR. &
          (WHAT == 'SPCF') .OR. (WHAT == 'MPCF')) THEN
         CONTINUE
      ELSE
         WRITE(ERR,9100) WHAT
         WRITE(F06,9100) WHAT
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

!     Write output headers.
      THERMAL_FLAG = 0 ! 1 for heat transfer, 0 otherwise

      ! TODO: is an eigenvector classified as displacement?
      ! TODO: where is velocity???
      ! TODO: can we return two values from a subroutine, so we don't have these
      !       hanging out?
      TABLE_NAME = 'OUG ERR '
      TABLE_CODE = -1 ! error
      CALL GET_TABLE_NAME_OUG(WHAT, TABLE_NAME, TABLE_CODE)
      IF (ITABLE .EQ. -1) THEN
          CALL WRITE_TABLE_HEADER(TABLE_NAME)
          ITABLE = -3
      ENDIF
      CALL WRITE_ITABLE(ITABLE)
      ITABLE = ITABLE - 1


      EIGENVALUE = 0.0
      MODE = 0
      ISUBCASE_INDEX = 0
      CALL GET_ANALYSIS_CODE_FIELD5_FIELD6(JSUB, ANALYSIS_CODE, MODE, EIGENVALUE, ISUBCASE_INDEX)

      TITLEI = TITLE(INT_SC_NUM)
      STITLEI = STITLE(INT_SC_NUM)
      LABELI = LABEL(INT_SC_NUM)

      ISUBCASE = SCNUM(ISUBCASE_INDEX)
      IF ((ANALYSIS_CODE == 1) .OR. (ANALYSIS_CODE == 10)) THEN
          ! static
          CALL WRITE_OUG3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NEW_RESULT, &
                                 TITLEI, STITLEI, LABELI)
      ELSE
          CALL WRITE_OUG3_EIGN(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, TABLE_CODE, NEW_RESULT, &
                               TITLEI, STITLEI, LABELI, MODE, EIGENVALUE)
      ENDIF

      ITABLE = ITABLE - 1
      ! Write accels, displ's, applied forces or SPC forces (also calc TOTALS for forces if that is being output)
      ! TOTALS(J) is summation of G.P. values of applied forces, SPC forces, or MFC forces, for each of the J=1,6 components.

      ! fill the G_OR_S array
      CALL GET_G_OR_S ( NUM, G_OR_S )

      ! write the real "displacment" data
      NUM_WIDE = 8
 100  FORMAT("*DEBUG:    NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
!      NGRID = NUM - 3
      NVALUES = NUM * NUM_WIDE
      NTOTAL = NVALUES * 4
      WRITE(ERR,100) NUM,NVALUES,NTOTAL
      WRITE(OP2) NVALUES
      ! Nastran OP2 requires this write call be a one liner...so it's a little weird...
      ! translating:
      !    DO I=1,NUM
      !        WRITE(OP2) GID_OUT_ARRAY(I,1)*10+DEVICE_CODE  ! Nastran is weird and requires scaling the NODE_ID
      !        WRITE(OP2) G_OR_S(I)                          ! GRID, SPOINT flag
      !
      !        write the TX, TY, TZ, RX, RY, RZ
      !        DO J=1,6
      !            FLOAT_VAL = REAL(OGEL(I,J), 4)   ! convert from float64 (double precision) to float32 (single precision)
      !            WRITE(OP2) FLOAT_VAL
      !        ENDDO
      !    ENDDO
      !
      WRITE(OP2) (GID_OUT_ARRAY(I,1)*10+DEVICE_CODE, G_OR_S(I), (REAL(OGEL(I,J),4), J=1,6), I=1,NUM)
      CALL END_OP2_TABLE(ITABLE)


      RETURN

! **********************************************************************************************************************************

! 9006 FORMAT('$',A,52X,I8)

 9100 FORMAT(' *ERROR  9100: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ILLEGAL INPUT FOR VARIABLE "WHAT" = ',A)


! **********************************************************************************************************************************

      END SUBROUTINE WRITE_GRD_OP2_OUTPUTS

!==============================================================================
      SUBROUTINE GET_TABLE_NAME_OUG ( WHAT, TABLE_NAME, TABLE_CODE )
      USE PENTIUM_II_KIND, ONLY     :  BYTE, LONG
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)  :: WHAT   ! Indicator whether to process displ or
                                              ! force output requests
      CHARACTER(LEN=8)  :: TABLE_NAME         ! Name of the op2 table that we're writing
      INTEGER(LONG)     :: TABLE_CODE         ! flag for the type of table

      IF (WHAT == 'DISP') THEN
        TABLE_NAME = 'OUGV1   '
        TABLE_CODE = 1
      ELSE IF (WHAT == 'VELO') THEN
        TABLE_NAME = 'OUGV1   '
        TABLE_CODE = 10
      ELSE IF (WHAT == 'ACCE') THEN
        TABLE_NAME = 'OUGV1   '
        TABLE_CODE = 11

      ELSE IF (WHAT == 'OLOAD') THEN
        TABLE_NAME = 'OPG1    '  ! TODO: should this be OPGV1?
        TABLE_CODE = 2

      ELSE IF (WHAT == 'SPCF') THEN
        TABLE_NAME = 'OQGV1   '
        TABLE_CODE = 3
      ELSE IF (WHAT == 'MPCF') THEN
        TABLE_NAME = 'OQMG1   '   ! TODO: should this be OQGV1/39?
        TABLE_CODE = 39

      ELSE
        TABLE_NAME = 'OUG ERR '
        TABLE_CODE = -1 ! error
      ENDIF
      END SUBROUTINE GET_TABLE_NAME_OUG

!==============================================================================
      SUBROUTINE GET_G_OR_S ( NUM, G_OR_S )
      USE MODEL_STUF, ONLY       :  GRID
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG
      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN)  :: NUM  ! The number of rows of OGEL to write out

      INTEGER(LONG)                  :: I       ! DO loop index
      INTEGER(LONG), DIMENSION(NUM)  :: G_OR_S  ! flag for the type of point

      ! putting this G/S calc into an array
      DO I=1,NUM
         ! type
         ! 0 - H / SECTOR/HARMONIC/RING POINT
         ! 1 - G / GRID
         ! 2 - S / SPOINT
         ! 3 - E / EXTRA POINT
         ! 4 - M / MODAL POINT
         ! 7 - L / RIGID POINT (e.g. RBE3)
         IF (GRID(I,6) == 1) THEN
            G_OR_S(I) = 2
         ELSE IF (GRID(I,6) == 6) THEN
            G_OR_S(I) = 1
         ELSE
            G_OR_S(I) = -1 ! error
         ENDIF
      ENDDO
      END SUBROUTINE GET_G_OR_S

!==============================================================================
      SUBROUTINE GET_ANALYSIS_CODE_FIELD5_FIELD6(JSUB, ANALYSIS_CODE, MODE, EIGENVALUE, ISUBCASE_INDEX)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR
      USE SCONTR, ONLY                :  SOL_NAME
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP

      CHARACTER(LEN=128)                :: LABELI            ! Subcase label

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: ANALYSIS_CODE     ! flag for the solution type
      INTEGER(LONG), INTENT(INOUT)    :: MODE              ! mode number for an eigenvector solution
      REAL(DOUBLE), INTENT(INOUT)     :: EIGENVALUE        ! the eigenvalue for an eigenvector solution
      INTEGER(LONG), INTENT(INOUT)    :: ISUBCASE_INDEX    ! the index into SCNUM

      IF (SOL_NAME(1:7) == 'STATICS') THEN
        ISUBCASE_INDEX = JSUB
        ANALYSIS_CODE = 1  ! statics
      ELSE IF((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
        ISUBCASE_INDEX = 1
        ANALYSIS_CODE = 2 ! eigenvectors
        EIGENVALUE = EIGEN_VAL(JSUB)
        MODE = JSUB
      ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
        ISUBCASE_INDEX = 1
        ANALYSIS_CODE = 1 ! statics

      ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
        ISUBCASE_INDEX = 2
        ANALYSIS_CODE = 7 ! pre-buckling
        EIGENVALUE = EIGEN_VAL(JSUB)
        MODE = JSUB
!      ELSE IF ???
!        ANALYSIS_CODE = 5 ! frequency
!      ELSE IF ???
!        ANALYSIS_CODE = 6 ! transient
!      ELSE IF ???
!        ANALYSIS_CODE = 9 ! complex eigenvectors
      ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
        ISUBCASE_INDEX = 1
        ANALYSIS_CODE = 10 ! nonlinear statics
      ELSE
        ANALYSIS_CODE = -1 ! error
 99     FORMAT("*ERROR: ANALYSIS_CODE=-1; SOL_NAME =",A)
        WRITE(ERR,99) SOL_NAME
      ENDIF
      END SUBROUTINE GET_ANALYSIS_CODE_FIELD5_FIELD6


      SUBROUTINE WRITE_FEMAP_GRID_VECS ( GRID_VEC, FEMAP_SET_ID, WHAT )

! Writes grid related vectors to FEMAP neutral file (displ, applied load, SPC and MPC forces)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, NEU
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD, NDOFG, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  CORD, GRID, GRID_ID, INV_GRID_SEQ

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE VECTOR_METRICS, ONLY        :  GET_VEC_MIN_MAX_ABS
      USE SORTING, ONLY               :  SORT_INT1_REAL1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_GRID_VECS'
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT              ! Indicator if GRID_VEC is DISP, OLOA, SPCF or MPCF
      CHARACTER(LEN= 3*BYTE)          :: TITLE1(4,2)       ! Titles for vectors written to NEU
      CHARACTER(LEN=20*BYTE)          :: TITLE2(2)         ! Titles for vectors written to NEU

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! FEMAP set ID to write out
      INTEGER(LONG)                   :: ACID_G            ! Actual coordinate system ID for a grid
      INTEGER(LONG)                   :: GRID_MAX          ! Grid ID where vector is max
      INTEGER(LONG)                   :: GRID_MIN          ! Grid ID where vector is min
      INTEGER(LONG)                   :: GRID_NUMS(NGRID)  ! Grid ID's in global order
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICID              ! Internal coord sys no. corresponding to an actual coord sys no.
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID for a grid in array GRID_NUMS
      INTEGER(LONG)                   :: IARRAY(NGRID)     ! Original GRID_NUMS array
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
      REAL(DOUBLE)                    :: TOTR_VEC(NGRID)   ! RSS of 6 rotation    components in  GRID_VEC
      REAL(DOUBLE)                    :: TOTT_VEC(NGRID)   ! RSS of 6 translation components in  GRID_VEC
      REAL(DOUBLE)                    :: T1_VEC(NGRID)     ! T1 translation component from GRID_VEC
      REAL(DOUBLE)                    :: T2_VEC(NGRID)     ! T2 translation component from GRID_VEC
      REAL(DOUBLE)                    :: T3_VEC(NGRID)     ! T3 translation component from GRID_VEC
      REAL(DOUBLE)                    :: R1_VEC(NGRID)     ! R1 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: R2_VEC(NGRID)     ! R2 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: R3_VEC(NGRID)     ! R3 rotation    component from GRID_VEC
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
         VEC_ID_OFFSET = 10000
         TITLE2(1) = ' translation'
         TITLE2(2) = ' rotation'
      ELSE IF (WHAT == 'OLOA') THEN
         VEC_ID_OFFSET = 20000
         TITLE2(1) = ' applied force'
         TITLE2(2) = ' applied moment'
      ELSE IF (WHAT == 'SPCF') THEN
         VEC_ID_OFFSET = 30000
         TITLE2(1) = ' SPC force'
         TITLE2(2) = ' SPC moment'
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

   END MODULE GRID_OUTPUT_WRITERS
