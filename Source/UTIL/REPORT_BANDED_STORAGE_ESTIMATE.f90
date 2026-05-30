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

      SUBROUTINE REPORT_BANDED_STORAGE_ESTIMATE ( MAT_NAME, NROWS, NTERMS, I_MAT, J_MAT, BAND_WIDTH, SUPINFO_IN )

! --- BANDED_optimizisation -begin-- !
! Reports Stage 1 storage estimates for banded solver planning. This routine does not change solver dispatch.
! --- BANDED_optimizisation -end-- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: SUPINFO_IN

      INTEGER(LONG), INTENT(IN)       :: BAND_WIDTH
      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NTERMS
      INTEGER(LONG), INTENT(IN)       :: I_MAT(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MAT(NTERMS)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: PROFILE_HEIGHT
      INTEGER(LONG)                   :: PROFILE_MAX_HEIGHT
      INTEGER(LONG), ALLOCATABLE      :: PROFILE_FIRST_ROW(:)

      REAL(DOUBLE)                    :: BAND_MB_EST
      REAL(DOUBLE)                    :: BAND_TERMS_EST
      REAL(DOUBLE)                    :: CSR_MB_EST
      REAL(DOUBLE)                    :: PROFILE_MB_EST
      REAL(DOUBLE)                    :: PROFILE_TERMS_EST

! **********************************************************************************************************************************
      BAND_TERMS_EST     = REAL(BAND_WIDTH,DOUBLE)*REAL(NROWS,DOUBLE)
      BAND_MB_EST        = REAL(DOUBLE,DOUBLE)*BAND_TERMS_EST/ONEPP6
      CSR_MB_EST         = (REAL(LONG,DOUBLE)*REAL(NROWS+1,DOUBLE) + REAL(LONG+DOUBLE,DOUBLE)*REAL(NTERMS,DOUBLE))/ONEPP6
      PROFILE_TERMS_EST  = ZERO
      PROFILE_MB_EST     = ZERO
      PROFILE_MAX_HEIGHT = 0

      ALLOCATE (PROFILE_FIRST_ROW(NROWS),STAT=IERR)
      IF (IERR == 0) THEN
         DO I=1,NROWS
            PROFILE_FIRST_ROW(I) = I
         ENDDO

         K = 0
         DO I=1,NROWS
            DO J=1,I_MAT(I+1)-I_MAT(I)
               K = K + 1
               IF ((K <= NTERMS) .AND. (J_MAT(K) >= 1) .AND. (J_MAT(K) <= NROWS)) THEN
                  IF (I < PROFILE_FIRST_ROW(J_MAT(K))) THEN
                     PROFILE_FIRST_ROW(J_MAT(K)) = I
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         DO I=1,NROWS
            PROFILE_HEIGHT = I - PROFILE_FIRST_ROW(I) + 1
            IF (PROFILE_HEIGHT < 1) PROFILE_HEIGHT = 1
            PROFILE_TERMS_EST = PROFILE_TERMS_EST + REAL(PROFILE_HEIGHT,DOUBLE)
            IF (PROFILE_HEIGHT > PROFILE_MAX_HEIGHT) PROFILE_MAX_HEIGHT = PROFILE_HEIGHT
         ENDDO

         PROFILE_MB_EST = REAL(DOUBLE,DOUBLE)*PROFILE_TERMS_EST/ONEPP6
         DEALLOCATE (PROFILE_FIRST_ROW)
      ENDIF

      WRITE(ERR,3095) MAT_NAME, NROWS, NTERMS, BAND_WIDTH, BAND_TERMS_EST, BAND_MB_EST, CSR_MB_EST,                              &
                      PROFILE_TERMS_EST, PROFILE_MB_EST, PROFILE_MAX_HEIGHT
      IF (SUPINFO_IN == 'N') THEN
         WRITE(F06,3095) MAT_NAME, NROWS, NTERMS, BAND_WIDTH, BAND_TERMS_EST, BAND_MB_EST, CSR_MB_EST,                           &
                         PROFILE_TERMS_EST, PROFILE_MB_EST, PROFILE_MAX_HEIGHT
      ENDIF

      RETURN

! **********************************************************************************************************************************
 3095 FORMAT(' *INFORMATION: STORAGE ESTIMATE FOR MATRIX ',A11,': N = ',I12,', STORED NNZ = ',I12,                               &
                    /,14X,' BAND WIDTH = ',I12,', COMPACT BAND TERMS = ',1ES13.6,', COMPACT BAND MB = ',F13.3,                  &
                    /,14X,' CSR MEMORY ESTIMATE MB = ',F13.3,', SKYLINE/PROFILE TERMS = ',1ES13.6,', SKYLINE/PROFILE MB = ',F13.3,&
                    /,14X,' MAX SKYLINE/PROFILE HEIGHT = ',I12,' (ESTIMATOR ONLY; NO SOLVER DISPATCH CHANGE)',/)

! **********************************************************************************************************************************

      END SUBROUTINE REPORT_BANDED_STORAGE_ESTIMATE
