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
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
! EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
! THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE WRITE_MATRIX_MARKET_VECTOR ( VEC_NAME, NUM, VEC, ISUB )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, INFILE, LEN_INPUT_FNAME
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MATRIX_MARKET_VECTOR'
      CHARACTER(LEN=*), INTENT(IN)    :: VEC_NAME
      INTEGER(LONG), INTENT(IN)       :: NUM
      INTEGER(LONG), INTENT(IN)       :: ISUB
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: LBASE
      INTEGER(LONG)                   :: LNAME
      INTEGER(LONG)                   :: UNT
      REAL(DOUBLE), INTENT(IN)        :: VEC(NUM)
      CHARACTER(256*BYTE)             :: FILNAM
      CHARACTER(32*BYTE)              :: VEC_TAG

      UNT = 92
      FILNAM = ' '
      VEC_TAG = ' '

      LBASE = LEN_INPUT_FNAME - 1
      IF (LBASE < 1) THEN
         LBASE = LEN_TRIM(INFILE)
      ENDIF

      LNAME = MIN(LEN_TRIM(VEC_NAME), LEN(VEC_TAG))
      VEC_TAG(1:LNAME) = VEC_NAME(1:LNAME)

      DO I=1,LNAME
         IF (VEC_TAG(I:I) == ' ') VEC_TAG(I:I) = '_'
      ENDDO

      FILNAM(1:LBASE) = INFILE(1:LBASE)
      FILNAM(LBASE+1:LBASE+1) = '_'
      FILNAM(LBASE+2:LBASE+1+LNAME) = VEC_TAG(1:LNAME)
      WRITE(FILNAM(LBASE+2+LNAME:LBASE+18+LNAME),'("_subcase_",I4.4,".mtx")') ISUB

      OPEN(UNIT=UNT,FILE=FILNAM,STATUS='REPLACE',FORM='FORMATTED',ACTION='WRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1001) SUBR_NAME, VEC_NAME, FILNAM(1:LEN_TRIM(FILNAM))
         WRITE(F06,1001) SUBR_NAME, VEC_NAME, FILNAM(1:LEN_TRIM(FILNAM))
         RETURN
      ENDIF

      WRITE(UNT,'(A)') '%%MatrixMarket matrix array real general'
      WRITE(UNT,'(A)') '% MYSTRAN vector export as N x 1 dense array'
      WRITE(UNT,'(I0,1X,I0)') NUM, 1

      DO I=1,NUM
         WRITE(UNT,'(ES24.16E3)') VEC(I)
      ENDDO

      CLOSE(UNT)

      WRITE(ERR,1002) VEC_NAME, ISUB, FILNAM(1:LEN_TRIM(FILNAM))
      WRITE(F06,1002) VEC_NAME, ISUB, FILNAM(1:LEN_TRIM(FILNAM))

 1001 FORMAT(' *WARNING    : ',A,' COULD NOT OPEN MATRIX MARKET FILE FOR ',A,' : ',A)
 1002 FORMAT(' *INFORMATION: MATRIX MARKET EXPORT WRITTEN FOR ',A,' SUBCASE ',I0,' TO FILE ',A)

      END SUBROUTINE WRITE_MATRIX_MARKET_VECTOR
