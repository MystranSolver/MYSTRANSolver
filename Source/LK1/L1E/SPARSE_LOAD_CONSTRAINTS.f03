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

   MODULE SPARSE_LOAD_CONSTRAINTS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SPARSE_PG, SPARSE_RMG

   CONTAINS

      SUBROUTINE SPARSE_PG

! Convert full array SYS_LOAD of all loads for all subcases to sparse array PG and write data to file LINK1E. The data that is
! written (for nonzero loads) is: G-set DOF number, internal subcase number, non-zero load value

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1E, L1E_MSG, L1ESTAT, LINK1E, SC1, WRT_ERR
      USE SCONTR, ONLY                :  FATAL_ERR, NDOFG, NSUB, NTERM_PG, BLNK_SUB_NAM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL, PRTFOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP, NL_NUM_LOAD_STEPS
      USE MODEL_STUF, ONLY            :  SYS_LOAD
      USE SPARSE_MATRICES, ONLY       :  I_PG, J_PG, PG

      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_PG'

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: KTERM_PG          ! Count of the number of terms written to file L1E for PG loads
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Open L1E to write G-set loads to.

      OUNT(1) = ERR
      OUNT(2) = F06
      CALL FILE_OPEN ( L1E, LINK1E, OUNT, 'REPLACE', L1E_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! Count the nonzero's in SYS_LOAD

      DO I=1,NDOFG
         DO J=1,NSUB
            IF (DABS(SYS_LOAD(I,J)) > EPS1) THEN
               NTERM_PG = NTERM_PG + 1
            ENDIF
         ENDDO
      ENDDO

! Allocate PG sparse arrays

      CALL ALLOCATE_SPARSE_MAT ( 'PG', NDOFG, NTERM_PG, SUBR_NAME )

! Formulate I_PG

      I_PG(1) = 1
      DO I=1,NDOFG
         I_PG(I+1) = I_PG(I)
         DO J=1,NSUB
            IF (DABS(SYS_LOAD(I,J)) > EPS1) THEN
               I_PG(I+1) = I_PG(I+1) + 1
!xx?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
!xx         ELSE
!xx            I_PG(I+1) = I_PG(I)
!xx 11/06/08: THESE LINES SHOULD BE COMMENTED OUT SINCE, IF SYS_LOAD(I,J) IS ZERO, WE DON'T WANT THAT TERM IN PG. ALLOWING THE ABOVE
!xx           TO BE IN THE CODE MADE A SERIOUS ERROR IN RUN TRUSS-2-CRODS THAT I DIDN'T DISCOVER UNTIL THIS DATE
!xx?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
            ENDIF
         ENDDO
      ENDDO

! Write the systems loads matrix to L1E in sparse form

      WRITE(L1E) NTERM_PG
      KTERM_PG = 0
      WRITE(SC1, * )
      CALL COUNTER_INIT('      Write L1E G-set DOF', NDOFG)
      DO I=1,NDOFG                                         ! Inner loop must be over subcases so that READ_MATRIX_1 will be reading
         DO J=1,NSUB                                       ! one row of sparse PG after another
            IF (DABS(SYS_LOAD(I,J)) > EPS1) THEN
               KTERM_PG = KTERM_PG + 1
               J_PG(KTERM_PG) = J
               IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
                  IF (NL_NUM_LOAD_STEPS > 0) THEN
                     SYS_LOAD(I,J) = SYS_LOAD(I,J)*LOAD_ISTEP/NL_NUM_LOAD_STEPS
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,965) SUBR_NAME, NL_NUM_LOAD_STEPS
                     WRITE(F06,965) SUBR_NAME, NL_NUM_LOAD_STEPS
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ENDIF
               PG(KTERM_PG) = SYS_LOAD(I,J)
               WRITE(L1E) I, J, SYS_LOAD(I,J)
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

      IF (KTERM_PG /= NTERM_PG) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1627) KTERM_PG, NTERM_PG
         WRITE(F06,1627) KTERM_PG, NTERM_PG
      ENDIF


      IF (NTERM_PG > 0) THEN
         CALL FILE_CLOSE ( L1E, LINK1E, 'KEEP' )
      ELSE
         CALL FILE_CLOSE ( L1E, LINK1E, L1ESTAT )
      ENDIF

      IF (PRTFOR(1) == 1) THEN                             ! Print PG if requested
         IF (NTERM_PG > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'G-SET LOADS, MATRIX PG', 'G ', 'SUBCASE', NTERM_PG, NDOFG, I_PG, J_PG, PG )
         ENDIF
      ENDIF

! **********************************************************************************************************************************
  965 FORMAT(' *ERROR   965: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NONLINEAR PARAMETER NL_NUM_LOAD_STEPS MUST BE > 0 BUT VALUE WAS ',I8)

 1627 FORMAT(' *ERROR  1614: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF G-SET LOAD MATRIX RECORDS WRITTEN TO FILE:'                                             &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' WAS KTERM_PG = ',I12,'. IT SHOULD HAVE BEEN NTERM_PG = ',I12)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_PG


      SUBROUTINE SPARSE_RMG

! Reads RMG constraint terms from file LINK1J. Zero terms are stripped and rows are sorted in numerical order. The final sparse RMG
! constraint matrix is written to file LINK1J in format: i, j, RMG(i,j)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1J, LINK1J, L1J_MSG
      USE SCONTR, ONLY                :  NDOFM, NTERM_RMG, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE SPARSE_MATRICES, ONLY       :  I_RMG, J_RMG, RMG

      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN, READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE SORTING, ONLY               :  SORT_INT2_REAL1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_RMG'

      INTEGER(LONG)                   :: I,K               ! DO loop indices or counters

      INTEGER(LONG)                   :: I2_RMG(NTERM_RMG) ! Row numbers of all terms in matrix RMG. Coming into this subr there is
!                                                            a conservative est of NTERM_RMG. The actual value will be determined
!                                                            by a count, herein, of the actual number of terms written to L1J

      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: IRMG              ! Row number for RMG
      INTEGER(LONG)                   :: IRMG_OLD          ! Row number for RMG
      INTEGER(LONG)                   :: JRMG              ! Col number for RMG
      INTEGER(LONG)                   :: KTERM_RMG         ! Count of number of terms in RMG
      INTEGER(LONG)                   :: NTERM_ROW_I       ! Number of nonzero terms in row I of RMG
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: RRMG              ! Real value for RMG

      INTRINSIC DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Read RMG constraint matrix. The matrix is not in DOF order

      IF (NDOFM > 0) THEN

         CALL FILE_OPEN ( L1J, LINK1J, OUNT, 'OLD', L1J_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

         NTERM_RMG = 0                                     ! First, calc NTERM_RMG
         REC_NO = 0
nterm:   DO
            IRMG = 0
            JRMG = 0
            RRMG = ZERO
            READ(L1J,IOSTAT=IOCHK) IRMG, JRMG, RRMG
            REC_NO = REC_NO + 1
            IF      (IOCHK == 0) THEN
               IF (DABS(RRMG) > EPS1) THEN
                  NTERM_RMG = NTERM_RMG + 1
                  CYCLE nterm
               ENDIF
            ELSE IF (IOCHK > 0) THEN
               CALL READERR ( IOCHK, LINK1J, L1J_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading RMG file, so quit
            ELSE
               EXIT nterm
            ENDIF
         ENDDO nterm
         CALL FILE_CLOSE ( L1J, LINK1J, 'KEEP' )
                                                           ! Allocate memory for RMG sparse arrays
         CALL ALLOCATE_SPARSE_MAT ( 'RMG', NDOFM, NTERM_RMG, SUBR_NAME )

         CALL FILE_OPEN ( L1J, LINK1J, OUNT, 'OLD', L1J_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

         KTERM_RMG = 0
         REC_NO    = 0
read_l1j:DO                                                ! Now calc sparse arrays for RMG
            IRMG      = 0
            JRMG      = 0
            RRMG      = ZERO
            READ(L1J,IOSTAT=IOCHK) IRMG, JRMG, RRMG
            REC_NO = REC_NO + 1
            IF      (IOCHK == 0) THEN
               IF (DABS(RRMG) > EPS1) THEN
                  KTERM_RMG         = KTERM_RMG + 1
                  IF (KTERM_RMG > NTERM_RMG) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM_RMG, 'I2_RMG, J_RMG, RMG' )
                  I2_RMG(KTERM_RMG) = IRMG
                   J_RMG(KTERM_RMG) = JRMG
                     RMG(KTERM_RMG) = RRMG
                  CYCLE read_l1j
               ENDIF
            ELSE IF (IOCHK > 0) THEN
               CALL READERR ( IOCHK, LINK1J, L1J_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                              ! Error reading RMG file, so quit
            ELSE
               EXIT read_l1j
            ENDIF
         ENDDO read_l1j

         I_RMG(1) = 1
         DO I=2,NTERM_RMG
            NTERM_ROW_I = 1
            IF (I2_RMG(I) /= I2_RMG(I)) THEN
               I_RMG(I) = I_RMG(I) + NTERM_ROW_I
            ENDIF
         ENDDO

         CALL FILE_CLOSE ( L1J, LINK1J, 'DELETE' )

! Sort RMG so that the rows (M-set DOF's) are in numerically increasing order

         CALL SORT_INT2_REAL1 ( SUBR_NAME, 'I2_RMG, J_RMG, RMG', NTERM_RMG, I2_RMG, J_RMG, RMG )

! Rewrite RMG matrix in the format i, j, RMG(i,j) to L1J

         CALL FILE_OPEN ( L1J, LINK1J, OUNT, 'REPLACE', L1J_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         WRITE(L1J) NTERM_RMG
         DO I=1,NTERM_RMG
            WRITE(L1J) I2_RMG(I),J_RMG(I),RMG(I)
         ENDDO
         CALL FILE_CLOSE ( L1J, LINK1J, 'KEEP' )

      ENDIF

      IRMG_OLD    = 0
      I_RMG(1) = 1
      DO K = 1,NTERM_RMG
         IRMG = I2_RMG(K)
         DO I=IRMG_OLD+1,IRMG
            I_RMG(I+1) = I_RMG(I)
         ENDDO
         IRMG_OLD = IRMG
         I_RMG(IRMG+1) = I_RMG(IRMG+1) + 1
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1616 FORMAT(' *ERROR  1616: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF TERMS IN THE RMG MATRIX   = ',I12,' BUT SHOULD BE NTERM_RMG = ',I12,' IN FILE:'         &
                    ,/,15X,A)



! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_RMG

   END MODULE SPARSE_LOAD_CONSTRAINTS
