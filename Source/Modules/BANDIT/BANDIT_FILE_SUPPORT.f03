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

   MODULE BANDIT_FILE_SUPPORT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BANDIT_FILES

   CONTAINS

      SUBROUTINE BANDIT_FILES ( IOU6, IOU7, IOU8, IOU9, IOU11, IOU12, IOU13, IOU14, IOU15, IOU16, IOU17 )

! Opens all output files and closes and deletes them so that no confusion about files if MYSTRAN aborts

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, F06, SC1

      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST        ! T/F depending on whether a file exists

      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_OUT        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F07        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F08        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F09        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F11        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F14        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F15        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F16        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F17        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F12        ! Bandit file name
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BANDIT_F13        ! Bandit file name

      CHARACTER(132*BYTE)             :: OUT_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F07_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F08_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F09_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F11_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F14_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F15_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F16_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F17_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F12_MSG           ! Description of Bandit file
      CHARACTER(132*BYTE)             :: F13_MSG           ! Description of Bandit file

      INTEGER(LONG), INTENT(IN)       :: IOU6              ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU7              ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU8              ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU9              ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU11             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU12             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU13             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU14             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU15             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU16             ! Bandit file unit number
      INTEGER(LONG), INTENT(IN)       :: IOU17             ! Bandit file unit number

      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR

! **********************************************************************************************************************************
! Default units for writing errors the screen (until LINK1A is read) and set filename length

      OUNT(1) = SC1
      OUNT(2) = SC1

! Formatted files.

      BANDIT_OUT(1:) = ' '         ;   OUT_MSG(1:) = ' '
      BANDIT_F07(1:) = ' '         ;   F07_MSG(1:) = ' '
      BANDIT_F08(1:) = ' '         ;   F08_MSG(1:) = ' '
      BANDIT_F09(1:) = ' '         ;   F09_MSG(1:) = ' '
      BANDIT_F11(1:) = ' '         ;   F11_MSG(1:) = ' '
      BANDIT_F14(1:) = ' '         ;   F14_MSG(1:) = ' '
      BANDIT_F15(1:) = ' '         ;   F15_MSG(1:) = ' '
      BANDIT_F16(1:) = ' '         ;   F16_MSG(1:) = ' '
      BANDIT_F17(1:) = ' '         ;   F17_MSG(1:) = ' '

      BANDIT_OUT = 'bandit.out'    ;   OUT_MSG = 'Bandit printed output file'
      BANDIT_F07 = 'bandit.f07'    ;   F07_MSG = 'SEQGP card images'
      BANDIT_F08 = 'bandit.f08'    ;   F08_MSG = 'MYSTRAN input file with SEQGP card images'
      BANDIT_F09 = 'bandit.f09'    ;   F09_MSG = 'References in NASNUM, SPRING, and FINISH'
      BANDIT_F11 = 'bandit.f11'    ;   F11_MSG = 'References in DOLLAR, NASNUM, and NEWIN'
      BANDIT_F14 = 'bandit.f14'    ;   F14_MSG = 'Output new element list for frontal solution'
      BANDIT_F15 = 'bandit.f15'    ;   F15_MSG = 'Component and grid list (CM must be run to get this list)'
      BANDIT_F16 = 'bandit.f16'    ;   F16_MSG = 'Connection table output ($TABLE YES)'
      BANDIT_F17 = 'bandit.f17'    ;   F17_MSG = 'Some run-time messages to benefit interactive running'

      INQUIRE ( FILE=BANDIT_OUT, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU6 , BANDIT_OUT, OUNT,'REPLACE', OUT_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU6 , BANDIT_OUT,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F07, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU7 , BANDIT_F07, OUNT,'REPLACE', F07_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU7 , BANDIT_F07,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F08, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU8 , BANDIT_F08, OUNT,'REPLACE', F08_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU8 , BANDIT_F08,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F09, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU9 , BANDIT_F09, OUNT,'REPLACE', F09_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU9 , BANDIT_F09,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F11, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU11, BANDIT_F11, OUNT,'REPLACE', F11_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU11, BANDIT_F11,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F14, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU14, BANDIT_F14, OUNT,'REPLACE', F14_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU14, BANDIT_F14,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F15, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU15, BANDIT_F15, OUNT,'REPLACE', F15_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU15, BANDIT_F15,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F16, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU16, BANDIT_F16, OUNT,'REPLACE', F16_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU16, BANDIT_F16,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F17, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU17, BANDIT_F17, OUNT,'REPLACE', F17_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU17, BANDIT_F17,'DELETE')
      ENDIF


! Unformatted files

      BANDIT_F12(1:) = ' '         ;   F12_MSG(1:) = ' '
      BANDIT_F13(1:) = ' '         ;   F13_MSG(1:) = ' '

      BANDIT_F12 = 'bandit.f12'    ;   F12_MSG = 'References in MPC, RESET, RIGID, and TIGER'
      BANDIT_F13 = 'bandit.f13'    ;   F13_MSG = 'Store elems for generating elem ordering for frontal solver'

      INQUIRE ( FILE=BANDIT_F12, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU12, BANDIT_F12, OUNT,'REPLACE', F12_MSG,'NEITHER','UNFORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU12, BANDIT_F12,'DELETE')
      ENDIF

      INQUIRE ( FILE=BANDIT_F13, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         CALL FILE_OPEN ( IOU13, BANDIT_F13, OUNT,'REPLACE', F13_MSG,'NEITHER','UNFORMATTED','READWRITE','REWIND','N','N')
         CALL FILE_CLOSE( IOU13, BANDIT_F13,'DELETE')
      ENDIF

! **********************************************************************************************************************************
  150 FORMAT(/,' >> MYSTRAN BEGIN  : ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,' The input file is ',A,/)

  170 FORMAT(/,' >> MYSTRAN RESTART: ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,' The input file is ',A,/)

! **********************************************************************************************************************************

      END SUBROUTINE BANDIT_FILES

   END MODULE BANDIT_FILE_SUPPORT
