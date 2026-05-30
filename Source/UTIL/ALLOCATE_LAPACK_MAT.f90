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

      SUBROUTINE ALLOCATE_LAPACK_MAT ( NAME, NROWS, NCOLS, CALLING_SUBR )

! Allocate matrices used in LAPACK band form

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, ONEPP6
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  SUPINFO, WINAMEM
      USE ARPACK_MATRICES_1 , ONLY    :  IWORK, RFAC, RESID, SELECT, VBAS, WORKD, WORKL
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, BBAND, LAPACK_S, RES

      USE ALLOCATE_LAPACK_MAT_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_LAPACK_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Name of matrix to be allocated
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in array to be allocated
      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols in array to be allocated
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: RNCOLS            ! Real value of NCOLS
      REAL(DOUBLE)                    :: RNROWS            ! Real value of NROWS

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      RNROWS = REAL(NROWS)
      RNCOLS = REAL(NCOLS)

      MB_ALLOCATED = ZERO
      JERR = 0

! LAPACK matrices

      IF (NAME == 'ABAND') THEN

         IF (ALLOCATED(ABAND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            WRITE(SC1,9199) NAME, MB_ALLOCATED
! --- BANDED_optimizisation -begin-- !
            IF (WINAMEM > ZERO) THEN
               WRITE(ERR,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
               IF (SUPINFO == 'N') WRITE(F06,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
            ELSE
               WRITE(ERR,942) NAME, NROWS, NCOLS, MB_ALLOCATED
               IF (SUPINFO == 'N') WRITE(F06,942) NAME, NROWS, NCOLS, MB_ALLOCATED
            ENDIF
            IF ((WINAMEM > ZERO) .AND. (MB_ALLOCATED > WINAMEM)) THEN
               WRITE(ERR,940) MB_ALLOCATED, NAME, WINAMEM
               WRITE(F06,940) MB_ALLOCATED, NAME, WINAMEM
               CALL OUTA_HERE ( 'Y' )
            ENDIF
! --- BANDED_optimizisation -end-- !
            ALLOCATE (ABAND(NROWS,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
         !xx   WRITE(SC1, * )                              ! Advance 1 line for screen messages
               CALL COUNTER_INIT("Init ABAND row", NROWS)
               DO I=1,NROWS
                  DO J=1,NCOLS
                     ABAND(I,J) = ZERO
                  ENDDO
                  CALL COUNTER_PROGRESS(I)
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'BBAND') THEN

         IF (ALLOCATED(BBAND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            WRITE(SC1,9199) NAME, MB_ALLOCATED
! --- BANDED_optimizisation -begin-- !
            IF (WINAMEM > ZERO) THEN
               WRITE(ERR,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
               IF (SUPINFO == 'N') WRITE(F06,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
            ELSE
               WRITE(ERR,942) NAME, NROWS, NCOLS, MB_ALLOCATED
               IF (SUPINFO == 'N') WRITE(F06,942) NAME, NROWS, NCOLS, MB_ALLOCATED
            ENDIF
            IF ((WINAMEM > ZERO) .AND. (MB_ALLOCATED > WINAMEM)) THEN
               WRITE(ERR,940) MB_ALLOCATED, NAME, WINAMEM
               WRITE(F06,940) MB_ALLOCATED, NAME, WINAMEM
               CALL OUTA_HERE ( 'Y' )
            ENDIF
! --- BANDED_optimizisation -end-- !
            ALLOCATE (BBAND(NROWS,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
         !xx   WRITE(SC1, * )                              ! Advance 1 line for screen messages
               CALL COUNTER_INIT("Init BBAND row", NROWS)
               DO I=1,NROWS
                  DO J=1,NCOLS
                     BBAND(I,J) = ZERO
                  ENDDO
                  CALL COUNTER_PROGRESS(I)
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RES') THEN

         IF (ALLOCATED(RES)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS/ONEPP6
            ALLOCATE (RES(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  RES(I) = ONE
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'LAPACK_S') THEN

         IF (ALLOCATED(LAPACK_S)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS/ONEPP6
            ALLOCATE (LAPACK_S(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  LAPACK_S(I) = ONE
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ARPACK matrices

      ELSE IF (NAME == 'IWORK') THEN

         IF (ALLOCATED(IWORK)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(LONG*NROWS)/ONEPP6
            ALLOCATE (IWORK(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  IWORK(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'SELECT') THEN

         IF (ALLOCATED(SELECT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(BYTE*NROWS)/ONEPP6
            ALLOCATE (SELECT(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  SELECT(I) = .FALSE.
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RESID') THEN

         IF (ALLOCATED(RESID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS/100000.
            ALLOCATE (RESID(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  RESID(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RFAC') THEN

         IF (ALLOCATED(RFAC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            WRITE(SC1,9199) NAME, MB_ALLOCATED
! --- BANDED_optimizisation -begin-- !
            IF (WINAMEM > ZERO) THEN
               WRITE(ERR,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
               IF (SUPINFO == 'N') WRITE(F06,941) NAME, NROWS, NCOLS, MB_ALLOCATED, WINAMEM
            ELSE
               WRITE(ERR,942) NAME, NROWS, NCOLS, MB_ALLOCATED
               IF (SUPINFO == 'N') WRITE(F06,942) NAME, NROWS, NCOLS, MB_ALLOCATED
            ENDIF
            IF ((WINAMEM > ZERO) .AND. (MB_ALLOCATED > WINAMEM)) THEN
               WRITE(ERR,940) MB_ALLOCATED, NAME, WINAMEM
               WRITE(F06,940) MB_ALLOCATED, NAME, WINAMEM
               CALL OUTA_HERE ( 'Y' )
            ENDIF
! --- BANDED_optimizisation -end-- !
            ALLOCATE (RFAC(NROWS,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
         !xx   WRITE(SC1, * )                              ! Advance 1 line for screen messages
               CALL COUNTER_INIT("Init RFAC row", NROWS)
               DO I=1,NROWS
                  DO J=1,NCOLS
                     RFAC(I,J) = ZERO
                  ENDDO
                  CALL COUNTER_PROGRESS(I)
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'VBAS') THEN

         IF (ALLOCATED(VBAS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            WRITE(SC1,9199) NAME, MB_ALLOCATED
            ALLOCATE (VBAS(NROWS,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
         !xx   WRITE(SC1, * )                              ! Advance 1 line for screen messages
               CALL COUNTER_INIT("Init VBAS row", NROWS)
               DO I=1,NROWS
                  DO J=1,NCOLS
                     VBAS(I,J) = ZERO
                  ENDDO
                  CALL COUNTER_PROGRESS(I)
               ENDDO
               WRITE(SC1,*) CR13
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'WORKD') THEN

         IF (ALLOCATED(WORKD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = (REAL(DOUBLE))*(REAL(NROWS))/100000.
            ALLOCATE (WORKD(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  WORKD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'WORKL') THEN

         IF (ALLOCATED(WORKL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            MB_ALLOCATED = (REAL(DOUBLE))*(REAL(NROWS))/100000.
            ALLOCATE (WORKL(NROWS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  WORKL(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  940 FORMAT(' *ERROR   940: ATTEMPT TO ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' EXCEEDS PARAM WINAMEM LIMIT OF ',          &
                             F10.3,' MB')

! --- BANDED_optimizisation -begin-- !
  941 FORMAT(' *INFORMATION: BANDED ALLOCATION REQUEST FOR ',A,': ROWS = ',I12,', COLS = ',I12,', ESTIMATED MB = ',F13.3,          &
                    /,14X,' PARAM WINAMEM LIMIT IS ACTIVE AT ',F13.3,' MB; EXCEEDANCE IS FATAL FOR BANDED ALLOCATION')

  942 FORMAT(' *INFORMATION: BANDED ALLOCATION REQUEST FOR ',A,': ROWS = ',I12,', COLS = ',I12,', ESTIMATED MB = ',F13.3,          &
                    /,14X,' PARAM WINAMEM LIMIT IS DISABLED; ALLOCATE/STAT WILL REPORT MEMORY FAILURE')
! --- BANDED_optimizisation -end-- !

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

 9199 FORMAT('      Memory needed for ',A,' = ',F13.6,' MB')

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_LAPACK_MAT
