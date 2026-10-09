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

   MODULE LINK1_WORKSPACE_LIFECYCLE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ALLOCATE_EMS_ARRAYS, ALLOCATE_L1_MGG, ALLOCATE_STF_ARRAYS, ALLOCATE_TEMPLATE, DEALLOCATE_EMS_ARRAYS, DEALLOCATE_L1_MGG, DEALLOCATE_STF_ARRAYS, DEALLOCATE_TEMPLATE

   CONTAINS

      SUBROUTINE ALLOCATE_EMS_ARRAYS ( CALLING_SUBR )

!  Allocate some arrays for use in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, ONEPP6
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LINKNO, LTERM_MGGE, NDOFG, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  YEAR, MONTH, DAY, HOUR, MINUTE, SEC, SFRAC, STIME, TSEC
      USE EMS_ARRAYS, ONLY            :  EMS, EMSCOL, EMSKEY, EMSPNT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_EMS_ARRAYS'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER( 6*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Number of rows for matrix NAME


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: RDOUBLE           ! Real value of DOUBLE
      REAL(DOUBLE)                    :: RLONG             ! Real value of LONG

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      RDOUBLE = REAL(DOUBLE)
      RLONG   = REAL(LONG)

      MB_ALLOCATED = ZERO
      JERR = 0

! Allocate array EMSKEY

      NAME = 'EMSKEY'
      NROWS = NDOFG
      IF (ALLOCATED(EMSKEY)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (EMSKEY(NDOFG),STAT=IERR)
         MB_ALLOCATED = RLONG*REAL(NDOFG)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            WRITE(SC1,22345,ADVANCE='NO') NAME, NDOFG, CR13
            DO I=1,NDOFG
!!             WRITE(SC1,12345,ADVANCE='NO') NAME, I, NDOFG, CR13
               EMSKEY(I) = 0
            ENDDO
            WRITE(SC1,*) CR13
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate array EMSCOL

      NAME = 'EMSCOL'
      NROWS = LTERM_MGGE
      IF (ALLOCATED(EMSCOL)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (EMSCOL(LTERM_MGGE),STAT=IERR)
         MB_ALLOCATED = RLONG*REAL(LTERM_MGGE)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            WRITE(SC1,22345,ADVANCE='NO') NAME, LTERM_MGGE, CR13
            DO I=1,LTERM_MGGE
!!             WRITE(SC1,12345,ADVANCE='NO') NAME, I, LTERM_MGGE, CR13
               EMSCOL(I) = 0
            ENDDO
            WRITE(SC1,*) CR13
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate array EMSPNT

      NAME = 'EMSPNT'
      NROWS = LTERM_MGGE
      IF (ALLOCATED(EMSPNT)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (EMSPNT(LTERM_MGGE),STAT=IERR)
         MB_ALLOCATED = RLONG*REAL(LTERM_MGGE)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            WRITE(SC1,22345,ADVANCE='NO') NAME, LTERM_MGGE, CR13
            DO I=1,LTERM_MGGE
!!             WRITE(SC1,12345,ADVANCE='NO') NAME, I, LTERM_MGGE, CR13
               EMSPNT(I) = 0
            ENDDO
            WRITE(SC1,*) CR13
        ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate array EMS

      NAME = 'EMS   '
      NROWS = LTERM_MGGE
      IF (ALLOCATED(EMS)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (EMS(LTERM_MGGE),STAT=IERR)
         MB_ALLOCATED = RDOUBLE*REAL(LTERM_MGGE)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            WRITE(SC1,22345,ADVANCE='NO') NAME, LTERM_MGGE, CR13
            DO I=1,LTERM_MGGE
!!             WRITE(SC1,12345,ADVANCE='NO') NAME, I, LTERM_MGGE, CR13
               EMS(I) = ZERO
            ENDDO
            WRITE(SC1,*) CR13
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

      WRITE(SC1,*) CR13

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

12345 FORMAT(5X,'array ',A6,' row ',I8,' of ',I8, A, A)

22345 FORMAT(5X,'Initializing ',A,' with ',I12,' rows', A, A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_EMS_ARRAYS



      SUBROUTINE ALLOCATE_L1_MGG ( NAME, CALLING_SUBR )

!  Allocate some arrays for use in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NTERM_MGG, NTERM_MGGC, NTERM_MGGE, NTERM_MGGS,            &
                                         TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_MGG, I2_MGG, J_MGG, MGG, I_MGGC, J_MGGC, MGGC, I_MGGE, J_MGGE, MGGE, I_MGGS, J_MGGS, MGGS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_L1_MGG'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Name of matrix to be allocated
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(6*BYTE)               :: NAMEO             ! Array name (used for output error message)

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Number of rows in array


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: RDOUBLE           ! Real value of DOUBLE
      REAL(DOUBLE)                    :: RLONG             ! Real value of LONG

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      RDOUBLE = REAL(DOUBLE)
      RLONG   = REAL(LONG)

      MB_ALLOCATED = ZERO
      JERR = 0

      IF (NAME == 'I2_MGG') THEN                           ! Allocate arrays for MGG

         NAMEO = 'I2_MGG'
         NROWS = NTERM_MGG
         IF (ALLOCATED(I2_MGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I2_MGG(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  I2_MGG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MGGC') THEN                         ! Allocate arrays for MGGC

         NAMEO = 'I_MGGC'
         NROWS = NDOFG + 1
         IF (ALLOCATED(I_MGGC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MGGC(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  I_MGGC(I) = 1
               ENDDO
            ELSE
               NAMEO = 'I_MGGC'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'J_MGGC'
         NROWS = NTERM_MGGC
         IF (ALLOCATED(J_MGGC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MGGC(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  J_MGGC(I) = 0
               ENDDO
            ELSE
               NAMEO = 'J_MGGC'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'MGGC'
         NROWS = NTERM_MGGC
         IF (ALLOCATED(MGGC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MGGC(NROWS),STAT=IERR)
            MB_ALLOCATED = RDOUBLE*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  MGGC(I) = ZERO
               ENDDO
            ELSE
               NAMEO = 'MGGC'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MGGE') THEN                         ! Allocate arrays for MGGE

         NAMEO = 'I_MGGE'
         NROWS = NDOFG+1
         IF (ALLOCATED(I_MGGE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MGGE(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  I_MGGE(I) = 1
               ENDDO
            ELSE
               NAMEO = 'I_MGGE'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'J_MGGE'
         NROWS = NTERM_MGGE
         IF (ALLOCATED(J_MGGE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MGGE(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  J_MGGE(I) = 0
               ENDDO
            ELSE
               NAMEO = 'J_MGGE'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'MGGE'
         NROWS = NTERM_MGGE
         IF (ALLOCATED(MGGE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MGGE(NROWS),STAT=IERR)
            MB_ALLOCATED = RDOUBLE*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  MGGE(I) = ZERO
               ENDDO
            ELSE
               NAMEO = 'MGGE'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MGGS') THEN                         ! Allocate arrays for MGGS

         NAMEO = 'I_MGGS'
         NROWS = NDOFG+1
         IF (ALLOCATED(I_MGGS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MGGS(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  I_MGGS(I) = 1
               ENDDO
            ELSE
               NAMEO = 'I_MGGS'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'J_MGGS'
         NROWS = NTERM_MGGS
         IF (ALLOCATED(J_MGGS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MGGS(NROWS),STAT=IERR)
            MB_ALLOCATED = RLONG*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  J_MGGS(I) = 0
               ENDDO
            ELSE
               NAMEO = 'J_MGGS'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAMEO = 'MGGS'
         NROWS = NTERM_MGGS
         IF (ALLOCATED(MGGS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MGGS(NROWS),STAT=IERR)
            MB_ALLOCATED = RDOUBLE*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS
                  MGGS(I) = ZERO
               ENDDO
            ELSE
               NAMEO = 'MGGS'
               WRITE(ERR,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAMEO,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1

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

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_L1_MGG


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

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

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
         IF (MB_NEEDED >= WINAMEM) THEN                 ! Reduce request for memory to
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



      SUBROUTINE ALLOCATE_TEMPLATE ( CALLING_SUBR )

!  Allocate some arrays for use in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE STF_TEMPLATE_ARRAYS, ONLY   :  CROW, TEMPLATE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_TEMPLATE'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Nunber of rows in array NAME being allocated


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      MB_ALLOCATED = ZERO
      NROWS = NDOFG
      JERR = 0

! Allocate array TEMPLATE

      NAME = 'TEMPLATE                  '
      IF (ALLOCATED(TEMPLATE)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (TEMPLATE(NDOFG,NDOFG),STAT=IERR)
         MB_ALLOCATED = REAL(BYTE)*REAL(NDOFG)*REAL(NDOFG)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,NDOFG
               DO J=1,NDOFG
                  TEMPLATE(I,J) = .FALSE.
               ENDDO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate array CROW

      NAME = 'CROW                  '
      IF (ALLOCATED(CROW)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (CROW(NDOFG),STAT=IERR)
         MB_ALLOCATED = REAL(BYTE)*REAL(LEN(CROW))*REAL(NDOFG)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,NDOFG
               CROW(I) = ' '
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)
! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_TEMPLATE



      SUBROUTINE DEALLOCATE_EMS_ARRAYS

!  Deallocate some arrays used in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE EMS_ARRAYS, ONLY            :  EMSCOL, EMSKEY, EMSPNT, EMS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_EMS_ARRAYS'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! Deallocate EMSCOL

      IF (ALLOCATED(EMSCOL)) THEN
         DEALLOCATE (EMSCOL,STAT=IERR)
         NAME = 'EMSCOL'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Deallocate EMSPNT

      IF (ALLOCATED(EMSPNT)) THEN
         DEALLOCATE (EMSPNT,STAT=IERR)
         NAME = 'EMSPNT'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Deallocate EMSKEY

      IF (ALLOCATED(EMSKEY)) THEN
         DEALLOCATE (EMSKEY,STAT=IERR)
         NAME = 'EMSKEY'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Deallocate EMS

      IF (ALLOCATED(EMS)) THEN
         DEALLOCATE (EMS,STAT=IERR)
         NAME = 'EMS'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_EMS_ARRAYS


      SUBROUTINE DEALLOCATE_L1_MGG ( NAME_IN )

!  Deallocate some arrays used in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SPARSE_MATRICES, ONLY       :  I_MGG, I2_MGG, J_MGG, MGG, I_MGGC, J_MGGC, MGGC, I_MGGE, J_MGGE, MGGE, I_MGGS, J_MGGS, MGGS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_L1_MGG'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Name of matrix to be allocated
      CHARACTER(6*BYTE)               :: NAME              ! Name of matrix to be allocated

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF (NAME_IN == 'I2_MGG') THEN

         IF (ALLOCATED(I2_MGG)) THEN
            DEALLOCATE (I2_MGG,STAT=IERR)
            NAME = 'I2_MGG'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MGGC') THEN

         IF (ALLOCATED(I_MGGC)) THEN
            DEALLOCATE (I_MGGC,STAT=IERR)
            NAME = 'I_MGGC'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(J_MGGC)) THEN
            DEALLOCATE (J_MGGC,STAT=IERR)
            NAME = 'J_MGGC'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(MGGC)) THEN
            DEALLOCATE (MGGC,STAT=IERR)
            NAME = 'MGGC'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MGGE') THEN

         IF (ALLOCATED(I_MGGE)) THEN
            DEALLOCATE (I_MGGE,STAT=IERR)
            NAME = 'I_MGGE'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(J_MGGE)) THEN
            DEALLOCATE (J_MGGE,STAT=IERR)
            NAME = 'J_MGGE'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(MGGE)) THEN
            DEALLOCATE (MGGE,STAT=IERR)
            NAME = 'MGGE'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(F06,992) NAME, SUBR_NAME
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MGGS') THEN

         IF (ALLOCATED(I_MGGS)) THEN
            DEALLOCATE (I_MGGS,STAT=IERR)
            NAME = 'I_MGGS'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(J_MGGS)) THEN
            DEALLOCATE (J_MGGS,STAT=IERR)
            NAME = 'J_MGGS'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (ALLOCATED(MGGS)) THEN
            DEALLOCATE (MGGS,STAT=IERR)
            NAME = 'MGGS'
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME, SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(F06,992) NAME, SUBR_NAME
            ENDIF
         ENDIF

      ELSE

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED' ,NAME_IN
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_L1_MGG


      SUBROUTINE DEALLOCATE_STF_ARRAYS ( NAME )

!  Deallocate some arrays used in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE STF_ARRAYS, ONLY            :  STFCOL, STFKEY, STFPNT, STF, STF3

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_STF_ARRAYS'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF      (NAME == 'STFKEY') THEN

         IF (ALLOCATED(STFKEY)) THEN
            DEALLOCATE (STFKEY,STAT=IERR)
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STFCOL') THEN

         IF (ALLOCATED(STFCOL)) THEN
            DEALLOCATE (STFCOL,STAT=IERR)
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STFPNT') THEN

         IF (ALLOCATED(STFPNT)) THEN
            DEALLOCATE (STFPNT,STAT=IERR)
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STF   ') THEN

         IF (ALLOCATED(STF)) THEN
            DEALLOCATE (STF,STAT=IERR)
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'STF3') THEN

         IF (ALLOCATED(STF3)) THEN
            DEALLOCATE (STF3,STAT=IERR)
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ELSE

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      WRITE(SC1,12345,ADVANCE='NO') CUR_MB_ALLOCATED, NAME, CR13
      WRITE(F06,1702) CUR_MB_ALLOCATED, NAME



      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

 1702 FORMAT('               DEALLOCATED ',1ES9.2,' MB MEMORY FROM ARRAY ',A)

12345 FORMAT(5X,'Deallocated ',1ES9.2,' MB of mem from array: ',A, A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_STF_ARRAYS


      SUBROUTINE DEALLOCATE_TEMPLATE

!  Deallocate some arrays used in LINK1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE STF_TEMPLATE_ARRAYS, ONLY   :  CROW, TEMPLATE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_TEMPLATE'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! Deallocate TEMPLATE

      IF (ALLOCATED(TEMPLATE)) THEN
         DEALLOCATE (TEMPLATE,STAT=IERR)
         NAME = 'TEMPLATE'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Deallocate CROW

      IF (ALLOCATED(CROW)) THEN
         DEALLOCATE (CROW,STAT=IERR)
         NAME = 'CROW'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_TEMPLATE

   END MODULE LINK1_WORKSPACE_LIFECYCLE
