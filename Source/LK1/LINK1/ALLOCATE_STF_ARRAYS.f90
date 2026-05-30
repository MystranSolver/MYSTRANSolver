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

      SUBROUTINE ALLOCATE_STF_ARRAYS ( NAME, CALLING_SUBR )

!  Allocate some arrays for use in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, ONEPP6
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LINKNO, LTERM_KGG, LTERM_KGGD, NDOFG, SOL_NAME,      &
                                         TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  YEAR, MONTH, DAY, HOUR, MINUTE, SEC, SFRAC, STIME, TSEC
      USE PARAMS, ONLY                :  MEMAFAC, MXALLOCA, SUPINFO, WINAMEM
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE STF_ARRAYS, ONLY            :  STF, STFCOL, STFKEY, STFPNT, STF3

      USE ALLOCATE_STF_ARRAYS_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_STF_ARRAYS'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name (used for output error message)
      CHARACTER( 1*BYTE)              :: ALLOC_SUCCESS     ! 'Y' if an allocation attempt was successful

      INTEGER(LONG)                   :: ALLOC_ATTEMPT_NUM ! The number of times an attempt has been made to allocate an array
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: LTERM             ! Count of number of estimated terms in KGG or KGGD
      INTEGER(LONG)                   :: NROWS             ! Number of rows  for matrix NAME
      INTEGER(LONG)                   :: NTERMS            ! Number of terms for matrix NAME


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of memory allocated for each array to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: MB_ALLOC_THIS_TIME! Megabytes of memory allocated for all arrays in this call
      REAL(DOUBLE)                    :: MB_NEEDED         ! Megabytes of memory needed for allocation
      REAL(DOUBLE)                    :: RDOUBLE           ! Real value of DOUBLE
      REAL(DOUBLE)                    :: RLONG             ! Real value of LONG

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
! Set LTERM, which will be the size allocated to the G-set stiffness matrix, to the appropriate value

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         LTERM = LTERM_KGGD
      ELSE
         LTERM = LTERM_KGG
      ENDIF

! Allocate arrays

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      RDOUBLE = REAL(DOUBLE)
      RLONG   = REAL(LONG)

      MB_ALLOCATED       = ZERO
      MB_ALLOC_THIS_TIME = ZERO
      JERR               = 0

      IF      (NAME == 'STFKEY') THEN

         NROWS = NDOFG
         IF (ALLOCATED(STFKEY)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STFKEY(NDOFG),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NDOFG)/ONEPP6
            MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               WRITE(SC1,12345,ADVANCE='NO') NAME, NDOFG, ' rows', CR13
               DO I=1,NDOFG
                  STFKEY(I) = 0
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STFCOL') THEN

         NROWS = LTERM
         IF (ALLOCATED(STFCOL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STFCOL(LTERM),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(LTERM)/ONEPP6
            MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               WRITE(SC1,12345,ADVANCE='NO') NAME, LTERM, ' terms', CR13
               DO I=1,LTERM
                  STFCOL(I) = 0
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STFPNT') THEN

         NROWS = LTERM
         IF (ALLOCATED(STFPNT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STFPNT(LTERM),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(LTERM)/ONEPP6
            MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               WRITE(SC1,12345,ADVANCE='NO') NAME, LTERM, ' terms', CR13
               DO I=1,LTERM
                  STFPNT(I) = 0
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STF   ') THEN

         NROWS = LTERM
         IF (ALLOCATED(STF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STF(LTERM),STAT=IERR)
            MB_ALLOCATED = RDOUBLE*REAL(LTERM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
               WRITE(SC1,12345,ADVANCE='NO') NAME, LTERM, ' terms', CR13
               DO I=1,LTERM
                  STF(I) = ZERO
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STF3  ') THEN

         NTERMS = LTERM
         MB_NEEDED = RDOUBLE*REAL(NTERMS)/ONEPP6 + TWO*RLONG*REAL(NTERMS)/ONEPP6
! --- BANDED_optimizisation -begin-- !
         IF ((WINAMEM > ZERO) .AND. (MB_NEEDED >= WINAMEM)) THEN
! --- BANDED_optimizisation -end-- !
            NTERMS = MEMAFAC*(WINAMEM/MB_NEEDED)*NTERMS
         ENDIF
         ALLOC_ATTEMPT_NUM = 1
         ALLOC_SUCCESS     = 'N'
         IF (ALLOCATED(STF3)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STF3(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               ALLOC_SUCCESS = 'Y'
               MB_ALLOCATED = RDOUBLE*REAL(NTERMS)/ONEPP6 + TWO*RLONG*REAL(NTERMS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
               WRITE(SC1,12345,ADVANCE='NO') NAME, NTERMS, ' terms', CR13
               WRITE(SC1,*) CR13
            ELSE
i_do:          DO
                  NTERMS = MEMAFAC*NTERMS
                  ALLOCATE (STF3(NTERMS),STAT=IERR)
                  IF (ALLOC_ATTEMPT_NUM <= MXALLOCA) THEN
                     MB_ALLOCATED = RDOUBLE*REAL(NTERMS)/ONEPP6 + TWO*RLONG*REAL(NTERMS)/ONEPP6
                     ALLOC_ATTEMPT_NUM = ALLOC_ATTEMPT_NUM + 1
                     IF (IERR == 0) THEN
                        ALLOC_SUCCESS = 'Y'
                        WRITE(SC1,32345,ADVANCE='NO') ALLOC_ATTEMPT_NUM, MB_ALLOCATED, NAME,' was successful', CR13
                        WRITE(SC1,*) CR13
                        CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
                        MB_ALLOC_THIS_TIME = MB_ALLOC_THIS_TIME + MB_ALLOCATED
                        EXIT i_do
                     ELSE
                        WRITE(SC1,32345,ADVANCE='NO') ALLOC_ATTEMPT_NUM, MB_ALLOCATED, NAME,' failed        ', CR13
                        WRITE(SC1,*) CR13
                        CYCLE i_do
                     ENDIF
                  ELSE
                     ALLOC_SUCCESS = 'N'
                     EXIT i_do
                  ENDIF
               ENDDO i_do
               WRITE(SC1,*) CR13
            ENDIF
         ENDIF

         IF (ALLOC_SUCCESS == 'Y') THEN
            DO I=1,NTERMS
               STF3(I)%Col_1 = 0
               STF3(I)%Col_2 = 0
               STF3(I)%Col_3 = ZERO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
            IF (ALLOC_ATTEMPT_NUM >= MXALLOCA) THEN
               WRITE(ERR,999) ALLOC_ATTEMPT_NUM, 'STF3', 'MXALLOCA'
               WRITE(F06,999) ALLOC_ATTEMPT_NUM, 'STF3', 'MXALLOCA'
            ENDIF
         ENDIF


      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

      WRITE(SC1,22345,ADVANCE='NO') MB_ALLOC_THIS_TIME, NAME, '                    ', CR13
      WRITE(SC1,*) CR13
      WRITE(ERR, 1702) MB_ALLOC_THIS_TIME, NAME
      IF (SUPINFO == 'N') THEN
         WRITE(F06, 1702) MB_ALLOC_THIS_TIME, NAME
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

  999 FORMAT('               THERE WERE ',I3,' ATTEMPTS TO ALLOCATE MEMORY TO ARRAY ',A                                            &
                    ,/,14X,' THE MAX ALLOWABLE ATTEMPTS CAN BE INCREASED VIA BULK DATA PARAM ',A)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

 1702 FORMAT('               ALLOCATED   ',1ES9.2,' MB MEMORY TO   ARRAY ',A)

12345 FORMAT(5X,'Initializing ',A,' with ',I12, A, A)

22345 FORMAT(5X,'Allocated   ',1ES9.2,' MB of mem to   array: ',A, A, A)

32345 FORMAT(5X,'Attempt ', I3,' to alloc ', F9.3,' MB memory to array ',A, A, A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_STF_ARRAYS
