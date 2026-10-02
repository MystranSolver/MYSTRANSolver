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

   MODULE EIGEN_MATRIX_LIFECYCLE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT

   CONTAINS

      SUBROUTINE ALLOCATE_EIGEN1_MAT ( NAME, NROWS, NCOLS, CALLING_SUBR )

! Allocate arrays for use in eigenvalue extraction

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, TOT_MB_MEM_ALLOC
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL, EIGEN_VEC, GEN_MASS, MODE_NUM, MEFFMASS, MPFACTOR_N6, MPFACTOR_NR

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_EIGEN1_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name of the matrix to be allocated in sparse format
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows to allocate to matrix NAME
      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols to allocate to matrix NAME
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

      IF (NAME == 'EIGEN_VAL') THEN

         IF (ALLOCATED(EIGEN_VAL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (EIGEN_VAL(NROWS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  EIGEN_VAL(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'EIGEN_VEC') THEN

         IF (ALLOCATED(EIGEN_VEC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (EIGEN_VEC(NROWS,NCOLS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  DO J=1,NCOLS
                     EIGEN_VEC(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'GEN_MASS') THEN

         IF (ALLOCATED(GEN_MASS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GEN_MASS(NROWS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  GEN_MASS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MODE_NUM') THEN

         IF (ALLOCATED(MODE_NUM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MODE_NUM(NROWS),STAT=IERR)
            MB_ALLOCATED = REAL(LONG)*REAL(NROWS)/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  MODE_NUM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MEFFMASS') THEN

         IF (ALLOCATED(MEFFMASS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MEFFMASS(NROWS,NCOLS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  DO J=1,NCOLS
                     MEFFMASS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MPFACTOR_N6') THEN

         IF (ALLOCATED(MPFACTOR_N6)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPFACTOR_N6(NROWS,NCOLS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  DO J=1,NCOLS
                     MPFACTOR_N6(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MPFACTOR_NR') THEN

         IF (ALLOCATED(MPFACTOR_NR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPFACTOR_NR(NROWS,NCOLS),STAT=IERR)
            MB_ALLOCATED = REAL(DOUBLE)*RNROWS*RNCOLS/ONEPP6
            IF (IERR == 0) THEN
               DO I=1,NROWS
                  DO J=1,NCOLS
                     MPFACTOR_NR(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

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

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_EIGEN1_MAT


      SUBROUTINE DEALLOCATE_EIGEN1_MAT ( NAME )

! Deallocate arrays in used in eigenvalue extraction

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL, EIGEN_VEC, GEN_MASS, MODE_NUM, MEFFMASS, MPFACTOR_N6, MPFACTOR_NR

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_EIGEN1_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name of the matrix to be allocated in sparse format

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF (NAME == 'EIGEN_VAL') THEN                        ! Deallocate arrays for EIGEN_VAL

        IF (ALLOCATED(EIGEN_VAL)) THEN
           DEALLOCATE (EIGEN_VAL,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'EIGEN_VEC') THEN                   ! Deallocate arrays for EIGEN_VEC

        IF (ALLOCATED(EIGEN_VEC)) THEN
           DEALLOCATE (EIGEN_VEC,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'GEN_MASS') THEN                    ! Deallocate arrays for GEN_MASS

        IF (ALLOCATED(GEN_MASS)) THEN
           DEALLOCATE (GEN_MASS,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'MODE_NUM') THEN                    ! Deallocate arrays for MODE_NUM

        IF (ALLOCATED(MODE_NUM)) THEN
           DEALLOCATE (MODE_NUM,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'MEFFMASS') THEN                    ! Deallocate arrays for MEFFMASS

        IF (ALLOCATED(MEFFMASS)) THEN
           DEALLOCATE (MEFFMASS,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'MPFACTOR_N6') THEN                 ! Deallocate arrays for MPFACTOR

        IF (ALLOCATED(MPFACTOR_N6)) THEN
           DEALLOCATE (MPFACTOR_N6,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE IF (NAME == 'MPFACTOR_NR') THEN                 ! Deallocate arrays for MPFACTOR

        IF (ALLOCATED(MPFACTOR_NR)) THEN
           DEALLOCATE (MPFACTOR_NR,STAT=IERR)
           IF (IERR /= 0) THEN
              WRITE(ERR,992) NAME,SUBR_NAME
              WRITE(F06,992) NAME,SUBR_NAME
              JERR = JERR + 1
           ENDIF
        ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_EIGEN1_MAT

   END MODULE EIGEN_MATRIX_LIFECYCLE
