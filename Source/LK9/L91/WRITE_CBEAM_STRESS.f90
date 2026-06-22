! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE WRITE_CBEAM_STRESS (NUM, WRITE_F06)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL

      USE WRITE_CBEAM_STRESS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_CBEAM_STRESS'

      INTEGER(LONG), INTENT(IN)       :: NUM
      LOGICAL,       INTENT(IN)       :: WRITE_F06
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IBEG
      INTEGER(LONG)                   :: IEND
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: NSTA
      INTEGER(LONG)                   :: ELEMENT_ID
      INTEGER(LONG)                   :: GRID_ID

! --- cbeam_stations begin --- !
      IF (.NOT. WRITE_F06) THEN
         RETURN
      ENDIF

      I = 1
      DO WHILE (I <= NUM)
         IBEG = I
         ELEMENT_ID = EID_OUT_ARRAY(I,1)
         DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == ELEMENT_ID))
            I = I + 1
         ENDDO
         IEND = I - 1
         NSTA = IEND - IBEG + 1

         IF (NSTA >= 1) THEN
            WRITE(F06,*)
            WRITE(F06,9001) 0, ELEMENT_ID
            DO ISTA = 1, NSTA
               IF (ISTA == 1) THEN
                   GRID_ID = GID_OUT_ARRAY(IBEG,2)
               ELSE IF (ISTA == NSTA) THEN
                  GRID_ID = GID_OUT_ARRAY(IBEG,3)
               ELSE
                  GRID_ID = 0
               ENDIF

               ! Beam stress/strain results are stored in OGEL as a pair of
               ! rows per output station. The first row contains the SXC..M.S.-T
               ! values and the second paired row carries the compressive margin.
               K = 2*(IBEG + ISTA - 2) + 1

               WRITE(F06,9002) GRID_ID, CBEAM_XL_OUT(IBEG + ISTA - 1),                                         &
                               OGEL(K,1), OGEL(K,2), OGEL(K,3), OGEL(K,4),                                       &
                               OGEL(K,6), OGEL(K,7), OGEL(K,8), OGEL(K + 1,8)
            ENDDO
         ENDIF
      ENDDO
! --- cbeam_stations end --- !

      RETURN

 9001 FORMAT(I1,8X,I8)
 9002 FORMAT(1X,I8,2X,F7.3,1X,8(1ES14.6))

      END SUBROUTINE WRITE_CBEAM_STRESS
