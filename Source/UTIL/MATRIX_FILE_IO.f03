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

   MODULE MATRIX_FILE_IO

   USE FILE_LIFECYCLE, ONLY :  FILE_CLOSE, FILE_OPEN, READERR
   USE PROGRESS_COUNTERS, ONLY :  COUNTER_INIT, COUNTER_PROGRESS
   USE DOF_ARRAY_INDEXING, ONLY :  ARRAY_SIZE_ERROR_1
   USE DIAGNOSTICS_MEMORY_REPORTING, ONLY :  GET_GRID_AND_COMP

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: READ_MATRIX_1, READ_MATRIX_2, WRITE_MATRIX_1, WRITE_SPARSE_CRS

   CONTAINS

      SUBROUTINE READ_MATRIX_1 ( FILNAM, UNT, OPND, CLOSE_IT, CLOSE_STAT, MESSAG, NAME, NTERM, READ_NTERM, NROWS  &
                               , I_MATOUT, J_MATOUT, MATOUT )

! Reads matrix data from an unformatted file into a sparse format described below. The format of the data in the file must be:

! If READ_NTERM = 'Y':

!    Record 1  : NTERM       = No. nonzero terms in matrix (also noumber of records following this one)
!    Record 1+i: i, j, value = row no., col no., nonzero value for matrix MATOUT

! If READ_NTERM = 'N':

!    Record i:   i, j, value = row no., col no., nonzero value for matrix MATOUT

! The matrix must be read in one row at a time and the rows MUST be in numerical order

! The sparse matrix format output from this subroutine is the matrix described in compressed row storage format:

!             I_MATOUT(1 to NROWS+1) : i-th value is index in MATOUT where matrix row i begins
!             J_MATOUT(1 to NTERM)   : k-th value is the matrix col no. of the k-th term in array MATOUT
!               MATOUT(1 to NTERM)   : k-th value is the k-th nonzero value in the matrix


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE FILE_LIFECYCLE, ONLY: READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_MATRIX_1'
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_IT          ! ='Y'/'N' whether to close UNT or note
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER(LEN=*), INTENT(IN)    :: READ_NTERM        ! If 'Y', read NTERM from file before reading matrix
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Matrix name
      CHARACTER(LEN=*), INTENT(IN)    :: OPND              ! If 'Y', then do not open UNT, If 'N', open it

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MATOUT
      INTEGER(LONG), INTENT(IN)       :: NTERM             ! Number of matrix terms that should be in FILNAM
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number of FILNAM
      INTEGER(LONG), INTENT(OUT)      :: I_MATOUT(NROWS+1) ! Row numbers for terms in matrix MATOUT
      INTEGER(LONG), INTENT(OUT)      :: J_MATOUT(NTERM)   ! Col numbers for terms in matrix MATOUT
      INTEGER(LONG)                   :: READ_ERR    = 0   ! Error count
      INTEGER(LONG)                   :: I,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: IROW              ! Integer row value read from FILNAM
      INTEGER(LONG)                   :: IROW_OLD          ! Previous value of IROW
      INTEGER(LONG)                   :: JCOL              ! Col number for MATOUT
      INTEGER(LONG)                   :: KTERM             ! Count of number of nonzero terms read from FILNAM
      INTEGER(LONG)                   :: NUM_TERMS         ! Head rec read from files that denotes how many records in FILNAM
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading FILNAM


      REAL(DOUBLE) , INTENT(OUT)      :: MATOUT(NTERM)     ! Real values for matrix MATOUT
      REAL(DOUBLE)                    :: RVAL              ! Real values read from FILNAM

      CHARACTER(LEN=7+LEN(NAME)+LEN(": read row")) :: COUNTER_TEMPLATE

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Quick return if there are no terms in the matrix

      IF (NTERM == 0) RETURN

! Initialize outputs

      DO I=1,NROWS+1
         I_MATOUT(I) = 0
      ENDDO

      DO I=1,NTERM
         J_MATOUT(I) = 0
           MATOUT(I) = ZERO
      ENDDO

      OUNT(1) = ERR
      OUNT(2) = F06

      IF (OPND == 'N') THEN
         CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
      ENDIF

! Should we read NTERM from file before reading matrix?

      IF (READ_NTERM == 'Y') THEN
         READ(UNT,IOSTAT=IOCHK) NUM_TERMS
         IF (IOCHK /= 0) THEN
            REC_NO = 1
            CALL READERR ( IOCHK, FILNAM, MESSAG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                                 ! Can't read NUM_TERMS from file, so quit
         ENDIF
         IF (NUM_TERMS /= NTERM) THEN
            WRITE(ERR, 924) SUBR_NAME, NAME, NUM_TERMS, NTERM, FILNAM
            WRITE(F06, 924) SUBR_NAME, NAME, NUM_TERMS, NTERM, FILNAM
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                                 ! Coding error (wrong no. terms in matrix), so quit
         ENDIF
      ENDIF

! If we got here, there was no problem with NTERM, so read matrix

      KTERM       = 0
      IROW_OLD    = 0
      I_MATOUT(1) = 1

!xx   WRITE(SC1, * )
      WRITE(COUNTER_TEMPLATE, 12345) NAME
      !CALL COUNTER_INIT(COUNTER_TEMPLATE, NROWS)
k_do1:DO K = 1,NTERM
         READ(UNT,IOSTAT=IOCHK) IROW,JCOL,RVAL
         IF (IOCHK /= 0) THEN
            IF (READ_NTERM == 'Y') THEN
               REC_NO = K + 1
            ELSE
               REC_NO = K
            ENDIF
            CALL READERR ( IOCHK, FILNAM, MESSAG, REC_NO, OUNT )
            READ_ERR = READ_ERR + 1                        ! Error reading IROW, JCOL, RVAL record from unit UNT
            CYCLE k_do1
         ELSE
            IF (IROW > IROW_OLD) THEN
               !CALL COUNTER_PROGRESS(IROW)
               DO I=IROW_OLD+1,IROW
                  I_MATOUT(I+1) = I_MATOUT(I)
               ENDDO
               IROW_OLD = IROW
            ELSE IF (IROW < IROW_OLD) THEN
               WRITE(ERR,926) SUBR_NAME, NAME, IROW_OLD, IROW, FILNAM
               WRITE(F06,926) SUBR_NAME, NAME, IROW_OLD, IROW, FILNAM
               READ_ERR = READ_ERR + 1                     ! Coding error (matrix stored incorrectly), so quit
               FATAL_ERR = FATAL_ERR + 1
               CYCLE k_do1
            ENDIF
         ENDIF
         I_MATOUT(IROW+1) = I_MATOUT(IROW+1) + 1
         KTERM            = KTERM + 1
         J_MATOUT(KTERM)  = JCOL
           MATOUT(KTERM)  = RVAL
      ENDDO k_do1
      WRITE(SC1,*) CR13

      IF (READ_ERR /= 0) THEN
         WRITE(ERR,9996) SUBR_NAME,READ_ERR
         WRITE(F06,9996) SUBR_NAME,READ_ERR
         CALL OUTA_HERE ( 'Y' )                            ! Quit due to above errors in k_do1 loop
      ENDIF

      IF (IROW < NROWS) THEN                               ! Fill out remainder of I_MATOUT, if needed
         DO I=IROW+1,NROWS
            I_MATOUT(I+1) = I_MATOUT(I)
         ENDDO
      ENDIF

      IF (CLOSE_IT == 'Y') THEN
         CALL FILE_CLOSE ( UNT, FILNAM, CLOSE_STAT )
      ENDIF

! Check sensibility of I_MATOUT

      IF (DEBUG(101) >= 2) THEN
         CALL CHECK_SPARSE_CRS_I ( NAME, SUBR_NAME, NROWS, NTERM, I_MATOUT, 101 )
      ENDIF

! Write I_A if requested

      IF (DEBUG(101) >= 1) THEN
         WRITE(F06,101) NAME, NAME
         WRITE(F06,102) (I_MATOUT(I),I=1,NROWS+1)
         WRITE(F06,*)
         WRITE(F06,103) NAME, NAME, NAME, NROWS
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(38X,'SPARSE ROW INDICATOR ARRAY I_',A,'(I) FOR MATRIX ',A,':',/)

  102 FORMAT(10I12)

  103 FORMAT(' NOTE: THE DIFFERENCE   I_',A,'(I+1) - I_',A,'(I)   SHOULD EQUAL THE NUMBER OF NONZERO TERMS IN ROW I OF ',A         &
          ,/,'       THE MATRIX HAS ',I8,' ROWS AND IS LISTED ABOVE WITH ELEMENTS 1-10 ON THE FIRST LINE, ELEMENTS 11-20 ON THE',  &
                   ' 2ND LINE, ETC.'/)

  924 FORMAT(' *ERROR   924: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF TERMS IN MATRIX ',A,' = ',I12,' BUT THE NUMBER OF TERMS SHOULD BE = ',I12,' IN FILE:'   &
                    ,/,15X,A)

  926 FORMAT(' *ERROR   926: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX ',A,' IS STORED INCORRECTLY. THE MATRIX MUST BE STORED BY ROWS IN NUMERICAL ORDER '            &
                    ,/,14X,' HOWEVER, ROW ',I12,' IS STORED BEFORE ROW ',I12,' IN FILE:'                                           &
                    ,/,15X,A)

 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

12345 FORMAT("       ", A, ': read row')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CHECK_SPARSE_CRS_I ( MAT_A_NAME, CALLING_SUBR, NROWS_A, NTERM_A, I_A, DEBUG_NUM )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=*)                :: CALLING_SUBR
      CHARACTER(LEN=*)                :: MAT_A_NAME
      CHARACTER( 1*BYTE)              :: QUIT = 'N'

      INTEGER(LONG) ,INTENT(IN)       :: DEBUG_NUM
      INTEGER(LONG) ,INTENT(IN)       :: I_A(NROWS+1)
      INTEGER(LONG) ,INTENT(IN)       :: NROWS_A
      INTEGER(LONG) ,INTENT(IN)       :: NTERM_A
      INTEGER(LONG)                   :: KTERM_A
      INTEGER(LONG)                   :: NTERMS_A_ROW_I    !
      INTEGER(LONG)                   :: NUM_ROW_ERRS

! **********************************************************************************************************************************
      NUM_ROW_ERRS = 0
      KTERM_A      = 0
      DO I=1,NROWS_A
         NTERMS_A_ROW_I = I_A(I+1) - I_A(I)
         IF (NTERMS_A_ROW_I < 0) THEN                      ! Error. This indicates row I has < 0 number of terms in it.
            NUM_ROW_ERRS = NUM_ROW_ERRS + 1
         ELSE
            KTERM_A = KTERM_A + NTERMS_A_ROW_I
         ENDIF
      ENDDO

      IF (NUM_ROW_ERRS > 0) THEN
         WRITE(ERR,913) SUBR_NAME, MAT_A_NAME, MAT_A_NAME, NUM_ROW_ERRS, DEBUG_NUM, MAT_A_NAME, CALLING_SUBR
         WRITE(F06,913) SUBR_NAME, MAT_A_NAME, MAT_A_NAME, NUM_ROW_ERRS, DEBUG_NUM, MAT_A_NAME, CALLING_SUBR
         IF (DEBUG(DEBUG_NUM) >= 3) THEN
            QUIT = 'Y'
         ENDIF
      ENDIF

      IF (KTERM_A /= NTERM_A) THEN
         WRITE(ERR,928) SUBR_NAME, MAT_A_NAME, MAT_A_NAME, NROWS_A, MAT_A_NAME, MAT_A_NAME, KTERM_A, NTERM_A, DEBUG_NUM,           &
                                   MAT_A_NAME, CALLING_SUBR
         WRITE(F06,928) SUBR_NAME, MAT_A_NAME, MAT_A_NAME, NROWS_A, MAT_A_NAME, MAT_A_NAME, KTERM_A, NTERM_A, DEBUG_NUM,           &
                                   MAT_A_NAME, CALLING_SUBR
         IF (DEBUG(DEBUG_NUM) >= 3) THEN
            QUIT = 'Y'
         ENDIF
      ENDIF

      IF (QUIT == 'Y') THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      RETURN

! **********************************************************************************************************************************
  913 FORMAT(' *ERROR   913: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX ',A,' HAS A SPARSE I_',A,' THAT IS INCORRECT.'                                                 &
                    ,/,14X,' THERE ARE AT LEAST ',I8,' ROWS WHERE THERE ARE < 0 NUMBER OF TERMS INDICATED.'                        &
                    ,/,14X,' IN ADDITION, THERE SHOULD BE NO ZERO TERMS IN THE MATRIX'              &
                    ,/,14X,' USE BULK DATA ENTRY DEBUG ',I3,' WITH VALUE >=2 TO GET A LISTING OF I_',A                             &
                    ,/,14X,' THIS SUBR WAS CALLED BY SUBR ',A,/)

  928 FORMAT(' *ERROR   928: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX ',A,' HAS A SPARSE I_',A,' THAT IS INCORRECT. THE SUM OVER I FROM 1 TO THE NUMBER OF ROWS = ', &
                             I8,' OF THE MATRIX TERMS:'                                                                            &
                    ,/,14X,' I_',A,'(I+1) - I_',A,'(I) = ',I8,' SHOULD EQUAL THE NUMBER OF NONZERO TERMS IN THE MATRIX WHICH = ',I8&
                    ,/,14X,' USE BULK DATA ENTRY DEBUG ',I3,' WITH VALUE >=3 TO GET A LISTING OF I_',A                             &
                    ,/,14X,' THIS SUBR WAS CALLED BY SUBR ',A,/)

! **********************************************************************************************************************************

      END SUBROUTINE CHECK_SPARSE_CRS_I

      END SUBROUTINE READ_MATRIX_1


      SUBROUTINE READ_MATRIX_2 ( FILNAM, UNT, OPND, CLOSE_IT, CLOSE_STAT, MESSAG, NAME, NROWS, NTERMS, READ_NTERM  &
                               , I2_MATOUT, J_MATOUT, MATOUT )

! Reads matrix data from an unformatted file into a sparse format described below The format of the data in the file must be:

! If READ_NTERM = 'Y':

!    Record 1  : NTERMS = number nonzero terms in matrix (also number records following this one)
!    Record 1+i: row number, col number, nonzero value for matrix MATOUT

! If READ_NTERM = 'N':

!    Record i:   row number, col number, nonzero value for matrix MATOUT

! The matrix must be read in one row at a time and the rows MUST be in numerical order

! The output from this subroutine is the matrix described in row, col, nonzero value format:

!            I2_MATOUT(1 - NTERMS) : k-th value is the matrix row number of the k-th term in array MATOUT
!             J_MATOUT(1 - NTERMS) : k-th value is the matrix col number of the k-th term in array MATOUT
!               MATOUT(1 - NTERMS) : k-th value is the k-th nonzero value in the matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE FILE_LIFECYCLE, ONLY: READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_MATRIX_2'
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_IT          ! ='Y'/'N' whether to close UNT or note
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER(LEN=*), INTENT(IN)    :: READ_NTERM        ! If 'Y', read NTERM from file before reading matrix
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Matrix name
      CHARACTER(LEN=*), INTENT(IN)    :: OPND              ! If 'Y', then do not open UNT, If 'N', open it

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in the matrix
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of matrix terms that should be in FILNAM
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number of FILNAM
      INTEGER(LONG), INTENT(OUT)      :: I2_MATOUT(NTERMS) ! Row numbers for terms in matrix MATOUT
      INTEGER(LONG), INTENT(OUT)      :: J_MATOUT(NTERMS)  ! Col numbers for terms in matrix MATOUT
      INTEGER(LONG)                   :: IERROR    = 0     ! Error count
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: K                 ! DO loop index
      INTEGER(LONG)                   :: NUM_TERMS         ! Head rec read from files that denotes how many records in FILNAM
      INTEGER(LONG)                   :: OLD_ROW_NUM       ! A variable used to tell when a new row of MATOUT is being read
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading FILNAM


      REAL(DOUBLE) , INTENT(OUT)      :: MATOUT(NTERMS)    ! Real values for matrix MATOUT

      CHARACTER(LEN=7+LEN(NAME)+LEN(": read row")) :: COUNTER_TEMPLATE

      INTRINSIC DABS



! **********************************************************************************************************************************
! Initialize outputs

      DO K=1,NTERMS
         I2_MATOUT(K) = 0
          J_MATOUT(K) = 0
            MATOUT(K) = ZERO
      ENDDO

! Calc outputs

      OUNT(1) = ERR
      OUNT(2) = F06

      IF (OPND == 'N') THEN
         CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
      ENDIF

! Should we read NTERMS from file before reading matrix?

      IF (READ_NTERM == 'Y') THEN
         READ(UNT,IOSTAT=IOCHK) NUM_TERMS
         IF (IOCHK /= 0) THEN
            REC_NO = 1
            CALL READERR ( IOCHK, FILNAM, MESSAG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                                 ! Can't read NUM_TERMS fro file, so quit
         ENDIF
         IF (NUM_TERMS /= NTERMS) THEN
            WRITE(ERR, 924) SUBR_NAME, NAME,NUM_TERMS, NTERMS, FILNAM
            WRITE(F06, 924) SUBR_NAME, NAME,NUM_TERMS, NTERMS, FILNAM
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                                 ! Coding error (wrong number terms in matrix), so quit
         ENDIF
      ENDIF

! If we got here, there was no problem with NTERMS, so read matrix

!xx   WRITE(SC1, * )
      OLD_ROW_NUM = 0
      IERROR  = 0
      !WRITE(COUNTER_TEMPLATE, 12345) NAME
      !CALL COUNTER_INIT(COUNTER_TEMPLATE, NROWS)
      DO K = 1,NTERMS
         READ(UNT,IOSTAT=IOCHK) I2_MATOUT(K), J_MATOUT(K), MATOUT(K)
         IF (IOCHK /= 0) THEN
            IF (READ_NTERM == 'Y') THEN
               REC_NO = K + 1
            ELSE
               REC_NO = K
            ENDIF
            CALL READERR ( IOCHK, FILNAM, MESSAG, REC_NO, OUNT )
            IERROR = IERROR + 1
         ENDIF
         IF (I2_MATOUT(K) > OLD_ROW_NUM) THEN
            !CALL COUNTER_PROGRESS(I2_MATOUT(K))
         ENDIF
      ENDDO
      WRITE(SC1,*) CR13

      IF (IERROR /= 0) THEN
         WRITE(ERR,9996) SUBR_NAME,IERROR
         WRITE(F06,9996) SUBR_NAME,IERROR
         CALL OUTA_HERE ( 'Y' )                                    ! Quit due to above errors reading matrix
      ENDIF

      IF (CLOSE_IT == 'Y') THEN
         CALL FILE_CLOSE ( UNT, FILNAM, CLOSE_STAT )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  924 FORMAT(' *ERROR   924: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF TERMS IN MATRIX ',A,' = ',I12,' BUT THE NUMBER OF TERMS SHOULD BE = ',I12,' IN FILE:'   &
                    ,/,15X,A)

 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

12345 FORMAT("       ",A,': read row')

! **********************************************************************************************************************************

      END SUBROUTINE READ_MATRIX_2


      SUBROUTINE WRITE_MATRIX_1 ( FILNAM, UNT, CLOSE_IT, CLOSE_STAT, MESSAG, NAME, NTERM, NROWS, I_MATIN, J_MATIN, MATIN )

! Writes sparse (compressed row storage format) matrix data to an unformatted file in the format:

!    Record 1  : NTERM       = No. nonzero terms in matrix (also no. records following this one)
!    Record 1+i: i, j, value = row no., col no., nonzero value for matrix MATOUT


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE

      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MATRIX_1'
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_IT          ! ='Y'/'N' whether to close UNT or not
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Matrix name
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERM             ! Number of matrix terms that should be in FILNAM
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number of FILNAM
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)  ! Row numbers for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM)    ! Col numbers for terms in matrix MATIN
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NTERM_ROW_I       ! Number of terms in row I of MATIN
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM)      ! Real values for matrix MATIN

      CHARACTER(LEN=LEN(NAME)+7+LEN(": writing row")) :: COUNTER_TEMPLATE

      INTRINSIC DABS



! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'REPLACE', MESSAG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! Write sparse (compressed row storage) matrix to file in i, j, val format:

      WRITE(UNT) NTERM
      K = 0
!xx   WRITE(SC1, * )
      WRITE(COUNTER_TEMPLATE, 12345) NAME
      CALL COUNTER_INIT(COUNTER_TEMPLATE, NROWS)
      DO I=1,NROWS
         NTERM_ROW_I = I_MATIN(I+1) - I_MATIN(I)
         DO J=1,NTERM_ROW_I
            K = K + 1
            IF (K > NTERM) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM, NAME)
            WRITE(UNT) I,J_MATIN(K),MATIN(K)
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

      IF (CLOSE_IT == 'Y') THEN
         CALL FILE_CLOSE ( UNT, FILNAM, CLOSE_STAT )
      ENDIF



      RETURN

! **********************************************************************************************************************************
12345 FORMAT("       ",A,': writing row')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MATRIX_1


      SUBROUTINE WRITE_SPARSE_CRS ( MAT_NAME, ROW_SET, COL_SET, NTERM_A, NROWS_A, I_AXX, J_AXX, AXX )

! Writes a matrix that is in sparse CRS format to the F06 output file based on user request via Bulk Data PARAM PRTijk entries

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SPARSTOR, TINY

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_SPARSE_CRS'
      CHARACTER(LEN=*), INTENT(IN)    :: COL_SET           ! Set designator for cols of matrix
      CHARACTER(LEN=*), INTENT(IN)    :: ROW_SET           ! Set designator for rows of matrix
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Input matrix descriptor
      CHARACTER(132*BYTE)             :: LINE_OUT          ! Line to print out (to describe matrix) that is centered

      INTEGER(LONG), INTENT(IN)       :: NTERM_A           ! No. of terms in sparse matrix
      INTEGER(LONG), INTENT(IN)       :: NROWS_A           ! No. of rows  in sparse matrix
      INTEGER(LONG), INTENT(IN)       :: I_AXX(NROWS_A+1)  ! Array of starting indices for the 1-st term in rows of AXX
      INTEGER(LONG), INTENT(IN)       :: J_AXX(NTERM_A)    ! Array of col no's for terms in matrix AXX
      INTEGER(LONG)                   :: COL_COMP    = 0   ! Component number returned from subr GET_GRID_AND_COMP
      INTEGER(LONG)                   :: COL_GRID    = 0   ! Grid number returned from subr GET_GRID_AND_COMP
      INTEGER(LONG)                   :: I,J               ! DO loop index
      INTEGER(LONG)                   :: INDEX             ! Index into character array LINE_OUT
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: MAT_NAME_LEN      ! Length of char array MAT_NAME. On input, it is the length as defined
!                                                            in the calling subr. In this subr, MAT_NAME is striped of trailing
!                                                            blanks to get only the actual message. On exit MAT_NAME_LEN is the
!                                                            length of the finite message in MAT_NAME (i.e. without trailing blanks)
      INTEGER(LONG)                   :: NTERM_ROW_I       ! Number of terms in row I of AXX
      INTEGER(LONG)                   :: NULL_ROWS_A       ! Number of null rows in input matrix
      INTEGER(LONG)                   :: ROW_COMP    = 0   ! Component number returned from subr GET_GRID_AND_COMP
      INTEGER(LONG)                   :: ROW_GRID    = 0   ! Grid number returned from subr GET_GRID_AND_COMP


      REAL(DOUBLE) , INTENT(IN)       :: AXX(NTERM_A)      ! Array of terms in matrix AXX



! **********************************************************************************************************************************
! Strip out trailing blanks from MAT_NAME and put remainder centered in array LINE_OUT

      MAT_NAME_LEN = LEN(MAT_NAME)                         ! This is the length of MAT_NAME as input (includes trailing blanks)

      DO  I=MAT_NAME_LEN,1,-1                              ! Calc length of description in MAT_NAME, excluding trailing blanks
         IF (MAT_NAME(I:I) == ' ') THEN
            CYCLE
         ELSE
            MAT_NAME_LEN = I
            EXIT
         ENDIF
      ENDDO

      LINE_OUT(1:) = ' '                                   ! Center MAT_NAME (w/0 trailing blanks) in LINE_OUT
      INDEX = (LEN(LINE_OUT) - MAT_NAME_LEN)/2 -12         ! Subt 12 to allow for SPARSTOR msg
      LINE_OUT(INDEX:) = MAT_NAME(1:MAT_NAME_LEN) // ' (with PARAM SPARSTOR = ' // SPARSTOR // ')'

! Determine how many null rows there are in the input matrix (used in output description of matrix)

      NULL_ROWS_A = 0
      DO I=1,NROWS_A
         NTERM_ROW_I = I_AXX(I+1) - I_AXX(I)
         IF (NTERM_ROW_I == 0) THEN
            NULL_ROWS_A = NULL_ROWS_A + 1
         ENDIF
      ENDDO

! Write the matrix out

      IF (COL_SET == 'SUBCASE') THEN                       ! Input matrix has as many cols as there are subcases

         IF (TINY == ZERO) THEN
            WRITE(F06,101) LINE_OUT, NROWS_A, NTERM_A, NULL_ROWS_A
         ELSE
            WRITE(F06,102) LINE_OUT, TINY, NROWS_A, NTERM_A, NULL_ROWS_A
         ENDIF
         K = 0
         DO I=1,NROWS_A

            NTERM_ROW_I = I_AXX(I+1) - I_AXX(I)
            IF (NTERM_ROW_I == 0) CYCLE
            ROW_GRID = 0
            ROW_COMP = 0
            IF (ROW_SET(1:) /= ' ') THEN
               CALL GET_GRID_AND_COMP ( ROW_SET, I, ROW_GRID, ROW_COMP )
            ENDIF

            IF ((ROW_GRID > 0) .AND. (ROW_COMP > 0)) THEN

               DO J=1,NTERM_ROW_I
                  K = K+1
                  IF (DABS(AXX(K)) > TINY) THEN
                     WRITE(F06,11) K, ROW_GRID, ROW_COMP, J_AXX(K), AXX(K)
                  ENDIF
               ENDDO

            ELSE

               DO J=1,NTERM_ROW_I
                  IF (DABS(AXX(K)) > TINY) THEN
                     WRITE(F06,12) K, I, J_AXX(K),AXX(K)
                  ENDIF
               ENDDO

            ENDIF

            WRITE(F06,*)

         ENDDO

         WRITE(F06,*)

      ELSE                                                 ! Input matrix has column numbers output as integers, not grid-comp

         IF (TINY == ZERO) THEN
            WRITE(F06,201) LINE_OUT, NROWS_A, NTERM_A, NULL_ROWS_A
         ELSE
            WRITE(F06,202) LINE_OUT, TINY, NROWS_A, NTERM_A, NULL_ROWS_A
         ENDIF
         K = 0
         DO I=1,NROWS_A

            NTERM_ROW_I = I_AXX(I+1) - I_AXX(I)
            IF (NTERM_ROW_I == 0) CYCLE
            ROW_GRID = 0
            ROW_COMP = 0
            IF (ROW_SET(1:) /= ' ') THEN
               CALL GET_GRID_AND_COMP ( ROW_SET, I, ROW_GRID, ROW_COMP )
            ENDIF

            IF ((ROW_GRID > 0) .AND. (ROW_COMP > 0)) THEN

               DO J=1,NTERM_ROW_I
                  K = K+1
                  COL_GRID = 0
                  COL_COMP = 0
                  IF (COL_SET(1:) /= ' ') THEN
                     CALL GET_GRID_AND_COMP ( COL_SET, J_AXX(K), COL_GRID, COL_COMP )
                  ENDIF
                  IF (( COL_GRID > 0) .AND. (COL_COMP > 0)) THEN
                     IF (DABS(AXX(K)) > TINY) THEN
                        WRITE(F06,21) K, ROW_GRID, ROW_COMP, COL_GRID, COL_COMP, AXX(K)
                     ENDIF
                  ELSE
                     IF (DABS(AXX(K)) > TINY) THEN
                        WRITE(F06,22) K, ROW_GRID, ROW_COMP, J_AXX(K), AXX(K)
                     ENDIF
                  ENDIF
               ENDDO

            ELSE

               DO J=1,NTERM_ROW_I
                  K = K+1
                  COL_GRID = 0
                  COL_COMP = 0
                  IF (COL_SET(1:) /= ' ') THEN
                     CALL GET_GRID_AND_COMP ( COL_SET, J_AXX(K), COL_GRID, COL_COMP )
                  ENDIF
                  IF (( COL_GRID > 0) .AND. (COL_COMP > 0)) THEN
                     IF (DABS(AXX(K)) > TINY) THEN
                        WRITE(F06,23) K, I, COL_GRID, COL_COMP, AXX(K)
                     ENDIF
                  ELSE
                     IF (DABS(AXX(K)) > TINY) THEN
                        WRITE(F06,24) K, I, J_AXX(K), AXX(K)
                     ENDIF
                  ENDIF
               ENDDO

            ENDIF

            WRITE(F06,*)

         ENDDO

         WRITE(F06,*)

      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(/,'                                    NONZERO TERMS OF SPARSE CRS (Compressed Row Storage) MATRIX:'                  &
            ,/,A                                                                                                                   &
            ,/,'                         (Matrix has ',I8,' rows and ',I8,' nonzeros. ',I8,' null rows not written)'               &
            ,/,'                            Term           Row                    Col                     Value'                   &
            ,/,'                                  (Grid-Comp or row no)  (Internal subcase no)',/)

  102 FORMAT(/,'                                    NONZERO TERMS OF SPARSE CRS (Compressed Row Storage) MATRIX:'                  &
            ,/,A                                                                                                                   &
            ,/,'                                    (Terms with abs value <= PARAM TINY =',1ES10.3,' not output)'                  &
            ,/,'                         (Matrix has ',I8,' rows and ',I8,' nonzeros. ',I8,' null rows not written)'               &
            ,/,'                            Term           Row                    Col                     Value'                   &
            ,/,'                                  (Grid-Comp or row no)  (Internal subcase no)',/)

  201 FORMAT(/,'                                    NONZERO TERMS OF SPARSE CRS (Compressed Row Storage) MATRIX:'                  &
            ,/,A                                                                                                                   &
            ,/,'                         (Matrix has ',I8,' rows and ',I8,' nonzeros. ',I8,' null rows not written)'               &
            ,/,'                            Term           Row                    Col                     Value'                   &
            ,/,'                                  (Grid-Comp or row no)  (Grid-Comp or col no)',/)

  202 FORMAT(/,'                                    NONZERO TERMS OF SPARSE CRS (Compressed Row Storage) MATRIX:'                  &
            ,/,A                                                                                                                   &
            ,/,'                                    (Terms with abs value <= PARAM TINY =',1ES10.3,' not output)'                  &
            ,/,'                         (Matrix has ',I8,' rows and ',I8,' nonzeros. ',I8,' null rows not written)'               &
            ,/,'                            Term           Row                    Col                     Value'                   &
            ,/,'                                  (Grid-Comp or row no)  (Grid-Comp or col no)',/)

   11 FORMAT(20X,I12,5X,I8,'-',I1,13X,I8,14X,1ES21.14,'     11')

   12 FORMAT(20X,I12,5X,I8,15X,I8,14X,1ES21.14,'     12')

   21 FORMAT(20X,I12,5X,I8,'-',I1,13X,I8,'-',I1,12X,1ES21.14,'     21')

   22 FORMAT(20X,I12,5X,I8,'-',I1,13X,I8,14X,1ES21.14,'     22')

   23 FORMAT(20X,I12,5X,I8,15X,I8,'-',I1,12X,1ES14.6,'     23')

   24 FORMAT(20X,I12,5X,I8,15X,I8,14X,1ES21.14,'     24')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_SPARSE_CRS

   END MODULE MATRIX_FILE_IO
