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

   MODULE INI_FILE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: READ_INI

   CONTAINS

      SUBROUTINE READ_INI ( INI_EXIST )

! Processes MYSTRAN.INI file, which contains default values for things such as the default drive, directory that the input file
! (and all output) will go. Also can change where some output files will go (to screen or printer rather than disk). This later
! feature is useful for debugging.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, DEFDIR, INIFIL, SC1, MOU4, WRT_ERR

      USE IOUNT1, ONLY                :  BUG,     ERR,     F06,     IN0,     IN1,     INI,     L1A,     NEU,                       &
                                         SEQ,     SPC,                                                                             &
                                         F21,     F22,     F23,     F24,     F25,                                                  &
                                         L1B,     L1C,     L1D,     L1E,     L1F,     L1G,     L1H,     L1I,     L1J,     L1K,     &
                                         L1L,     L1M,     L1N,     L1O,     L1P,     L1Q,     L1R,     L1S,     L1T,     L1U,     &
                                         L1V,     L1W,     L1X,     L1Y,     L1Z,                                                  &
                                         L2A,     L2B,     L2C,     L2D,     L2E,     L2F,     L2G,     L2H,     L2I,     L2J,     &
                                         L2K,     L2L,     L2M,     L2N,     L2O,     L2P,     L2Q,     L2R,     L2S,     L2T,     &
                                         L3A,     L4A,     L4B,     L4C,     L4D,     L5A,     L5B,     OP2,     OU4

      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, BUGSTAT, ERRSTAT, F06STAT, IN0STAT, IN1STAT,                            &
                                         L1ASTAT, NEUSTAT, SEQSTAT, SPCSTAT,                                                       &
                                         F21STAT, F22STAT, F23STAT, F24STAT, F25STAT,                                              &
                                         L1BSTAT, L1CSTAT, L1DSTAT, L1ESTAT, L1FSTAT, L1GSTAT, L1HSTAT, L1ISTAT, L1JSTAT, L1KSTAT, &
                                         L1LSTAT, L1MSTAT, L1NSTAT, L1OSTAT, L1PSTAT, L1QSTAT, L1RSTAT, L1SSTAT, L1TSTAT, L1USTAT, &
                                         L1VSTAT, L1WSTAT, L1XSTAT, L1YSTAT, L1ZSTAT,                                              &
                                         L2ASTAT, L2BSTAT, L2CSTAT, L2DSTAT, L2ESTAT, L2FSTAT, L2GSTAT, L2HSTAT, L2ISTAT, L2JSTAT, &
                                         L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, L2OSTAT, L2PSTAT, L2QSTAT, L2RSTAT, L2SSTAT, L2TSTAT, &
                                         L3ASTAT, L4ASTAT, L4BSTAT, L4CSTAT, L4DSTAT, L5ASTAT, L5BSTAT, OP2STAT, OU4STAT

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, IERRFL, INI_ENTRY_LEN, JF, LINKNO_START, PRINTENV,                    &
                                         PROG_NAME

      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD_08
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, SET_FILE_CLOSE_STAT, WRITE_FILNAM

      IMPLICIT NONE

      LOGICAL                         :: LEXIST

      CHARACTER( 1*BYTE), INTENT(OUT) :: INI_EXIST         ! 'Y' if file MYSTRAN.INI exists or 'N' otherwise
      CHARACTER( 8*BYTE)              :: ALL_CLOSE_STAT    ! Status to use when closing MYSTRAN files when no longer needed.
      CHARACTER( 8*BYTE)              :: CHAR8             ! Character field read from line in INI file
      CHARACTER(LEN=INI_ENTRY_LEN)    :: CARD              ! A card image from file MYSTRAN.INI
      CHARACTER( 8*BYTE)              :: DUMSTAT           ! Dummy status used for file SC1 (screen) that cannot have a close status
      CHARACTER( 1*BYTE)              :: FLD_ERR_MSG(10)   ! 'Y'/'N' designator of whether error msg has been written for a field
      CHARACTER( 8*BYTE)              :: JCARD_08(10)      ! The 10 fields of 8 characters making up CARD
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: MYSTRAN_DIR       ! Directory where MYSTRAN executable (and INI file) exist
      CHARACTER(LEN=LEN(ECHO))        :: PERM_ECHO         ! Set equal to ECHO, initially, so we can use ECHO = NONE for this subr
!                                                            At end, ECHO is reset by ECHO = PERM_ECHO
      CHARACTER( 1*BYTE)              :: RESPONSE          ! 'Y'/'N' response from user
      CHARACTER( 1*BYTE)              :: WRT_CARD          ! 'Y', 'N' indicator if line from INI file was written to show errors
      CHARACTER( 1*BYTE)              :: WRT_HDR           ! 'Y', 'N' indicator if header msg written when there are errors

      INTEGER(LONG)                   :: DUMUNIT           ! Dummy unit number
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IBEGIN            ! Counter used to find where DEFDIR begins on CARD
      INTEGER(LONG)                   :: IERR              ! Error count when reading MYSTRAN.INI file
      INTEGER(LONG)                   :: INIFIL_NAME_LEN   ! Length of INI file name (incl path)
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from OPEN/READ
      INTEGER(LONG)                   :: LINE_NUMBER       ! Line number in the INI file
      INTEGER(LONG)                   :: MYSTRAN_DIR_LEN   ! Length of MYSTRAN_DIR (not including trailing blanks)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to

! **********************************************************************************************************************************
! Default units for writing errors the screen (until LINK1A is read)

      OUNT(1) = SC1
      OUNT(2) = SC1

! Set default close status for files

      ALL_CLOSE_STAT = 'DELETE  '

! Initialize DEFDIR and INIFIL

      DEFDIR(1:FILE_NAM_MAXLEN) = ' '
      INIFIL(1:FILE_NAM_MAXLEN) = ' '

! Use ECHO = NONE when reading MYSTRAN.INI so subr CRDERR will print card errors. ECHO will be reset below.

      PERM_ECHO = ECHO
      ECHO = 'NONE    '

! Initialize array that will say whether an error message has been written for a card field

      DO I=1,10
         FLD_ERR_MSG = 'N'
      ENDDO

! Get directory where executable and INI file exist and check its (non-blank) length:

      CALL GET_MYSTRAN_DIR ( MYSTRAN_DIR, MYSTRAN_DIR_LEN )

! Get the INI file name including complete path

      CALL GET_INI_FILNAM ( MYSTRAN_DIR, MYSTRAN_DIR_LEN, INIFIL_NAME_LEN )

! Look for file the INI. If it exists, open and get default values.

      IERR = 0
      INI_EXIST = 'Y'
      INQUIRE (FILE=INIFIL,EXIST=LEXIST)
      IF (.NOT.LEXIST) THEN
         INI_EXIST = 'N'
         CALL WRITE_INIFIL_MSG ( '1', INIFIL_NAME_LEN )
      ELSE
         OPEN (INI,FILE=INIFIL,STATUS='OLD',IOSTAT=IOCHK)
         IF (IOCHK /= 0) THEN
            CALL OPNERR ( IOCHK, INIFIL, OUNT )
            DO
               WRITE(SC1,* ) ' Cannot open MYSTRAN.INI file. Continue? (Y/N)'
               WRITE(SC1,* ) ' If Y, Then default values will be used'
               WRITE(SC1,* )
               READ ( *, * ) RESPONSE
               IF ((RESPONSE == 'Y') .OR. (RESPONSE == 'y')) THEN
                  EXIT
               ELSE IF ((RESPONSE == 'N') .OR. (RESPONSE == 'n')) THEN
                  WRITE(SC1,*) ' Processing terminated based on user input'
                  CALL OUTA_HERE ( 'Y' )
               ELSE
                  CYCLE
               ENDIF
            ENDDO
         ENDIF

! Process data in MYSTRAN.INI

         WRT_HDR     = 'N'
         LINE_NUMBER = 0
         DO

            WRT_CARD    = 'N'
            READ(INI,101,IOSTAT=IOCHK) CARD
            LINE_NUMBER = LINE_NUMBER + 1

! Give message and CYCLE if error occurs. EXIT if encounter EOF/EOR. Otherwise read CARD.

            IF (IOCHK < 0) THEN                            ! EOF/EOR (no more lines to read in the file) so exit loop
               EXIT
            ELSE IF (IOCHK > 0) THEN
               IF (WRT_HDR  == 'N') THEN
                  CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
                  WRT_HDR  = 'Y'
               ENDIF
               IF (WRT_CARD == 'N') THEN
                  WRITE(SC1,'(1X,A79)') CARD(1:79)
                  WRT_CARD = 'Y'
               ENDIF
               WRITE(SC1,1501) LINE_NUMBER
         !xx   WRITE(SC1, * )
               IERR = IERR + 1
               CYCLE
            ELSE IF ((CARD(1:1) /= '$') .AND. (CARD(1:) /= ' ')) THEN

               CALL MKJCARD_08 ( CARD, JCARD_08 )

               IF (CARD(1:8) ==  'DEF DIR ') THEN
                  IBEGIN = 9
                  DO I=9,FILE_NAM_MAXLEN+9
                     IF (CARD(I:I) == ' ') THEN            ! Checks for leading blanks (starting in col 9) before DEFDIR chars begin
                        IBEGIN = IBEGIN + 1
                        CYCLE
                     ELSE
                        EXIT                               ! Once we get to col where DEFDIR begins, exit loop
                     ENDIF                                 !                                 ------
                  ENDDO
                  DEFDIR(1:) = CARD(IBEGIN:)

               ELSE IF (CARD(1:8) == 'LINKNOST') THEN
                  CALL I4FLD0 ( JCARD_08(2), JF(2), LINKNO_START, WRT_HDR, WRT_CARD, FLD_ERR_MSG )
                  CALL CARD_FLDS_NOT_BLANK0 ( JCARD_08, 0,3,4,5,6,7,8,9, WRT_HDR, WRT_CARD )
                  CALL CRDERR0 ( CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:8) == 'PRT ENV ') THEN
                  CALL C8FLD0 ( JCARD_08(2), JF(2), CHAR8)
                  PRINTENV = CHAR8(1:1)
                  CALL CARD_FLDS_NOT_BLANK0 ( JCARD_08, 0,3,4,5,6,7,8,9, WRT_HDR, WRT_CARD )
                  CALL CRDERR0 ( CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:8) == 'ALLFILES') THEN
                  DUMUNIT = 0
                  CALL READ_INI_LINE ( 'ALLFILES', DUMUNIT, ALL_CLOSE_STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )
                                                           ! Check for all MYSTRAN files except INI, SCR

               ELSE IF (CARD(1:3) == 'SC1') THEN           ! 01
                  DUMSTAT(1:) = ' '
                  CALL READ_INI_LINE ( 'SC1', SC1, DUMSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'BUG') THEN           ! 03 (NOTE: BUGSTAT not allowed to be changed here)
                  CALL READ_INI_LINE ( 'BUG', BUG, BUGSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'ERR') THEN           ! 04
                  CALL READ_INI_LINE ( 'ERR', ERR, ERRSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F06') THEN           ! 06 (NOTE: F06STAT not allowed to be changed here)
                  DUMSTAT(1:) = ' '
                  CALL READ_INI_LINE ( 'F06', F06, DUMSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'IN0') THEN           ! 07
                  CALL READ_INI_LINE ( 'IN0', IN0, IN0STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'IN1') THEN           ! 08 (NOTE: IN1STAT not allowed to be changed here)
                  DUMSTAT(1:) = ' '
                  CALL READ_INI_LINE ( 'IN1', IN1, DUMSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1A') THEN           ! 09
                  CALL READ_INI_LINE ( 'L1A', L1A, L1ASTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'NEU') THEN           ! 10 (NOTE: NEUSTAT not allowed to be changed here)
                  DUMSTAT(1:) = ' '
                  CALL READ_INI_LINE ( 'NEU', NEU, DUMSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'SEQ') THEN           ! 11
                  CALL READ_INI_LINE ( 'SEQ', SEQ, SEQSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'SPC') THEN           ! 12 (NOTE: SPCSTAT not allowed to be changed here)
                  DUMSTAT(1:) = ' '
                  CALL READ_INI_LINE ( 'SPC', SPC, DUMSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F21') THEN           ! 13
                  CALL READ_INI_LINE ( 'F21', F21, F21STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F22') THEN           ! 14
                  CALL READ_INI_LINE ( 'F22', F22, F22STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F23') THEN           ! 15
                  CALL READ_INI_LINE ( 'F23', F23, F23STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F24') THEN           ! 16
                  CALL READ_INI_LINE ( 'F24', F24, F24STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'F25') THEN           ! 17
                  CALL READ_INI_LINE ( 'F25', F25, F25STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1B') THEN           ! 18
                  CALL READ_INI_LINE ( 'L1B', L1B, L1BSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1C') THEN           ! 19
                  CALL READ_INI_LINE ( 'L1C', L1C, L1CSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1D') THEN           ! 20
                  CALL READ_INI_LINE ( 'L1D', L1D, L1DSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1E') THEN           ! 21
                  CALL READ_INI_LINE ( 'L1E', L1E, L1ESTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1F') THEN           ! 22
                  CALL READ_INI_LINE ( 'L1F', L1F, L1FSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1G') THEN           ! 23
                  CALL READ_INI_LINE ( 'L1G', L1G, L1GSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1H') THEN           ! 24
                  CALL READ_INI_LINE ( 'L1H', L1H, L1HSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1I') THEN           ! 25
                  CALL READ_INI_LINE ( 'L1I', L1I, L1ISTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1J') THEN           ! 26
                  CALL READ_INI_LINE ( 'L1J', L1J, L1JSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1K') THEN           ! 27
                  CALL READ_INI_LINE ( 'L1K', L1K, L1KSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1L') THEN           ! 28
                  CALL READ_INI_LINE ( 'L1L', L1L, L1LSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1M') THEN           ! 29
                  CALL READ_INI_LINE ( 'L1M', L1M, L1MSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1N') THEN           ! 30
                  CALL READ_INI_LINE ( 'L1N', L1N, L1NSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1O') THEN           ! 31
                  CALL READ_INI_LINE ( 'L1O', L1O, L1OSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1P') THEN           ! 32
                  CALL READ_INI_LINE ( 'L1P', L1P, L1PSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1Q') THEN           ! 33
                  CALL READ_INI_LINE ( 'L1Q', L1Q, L1QSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1R') THEN           ! 34
                  CALL READ_INI_LINE ( 'L1R', L1R, L1RSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1S') THEN           ! 35
                  CALL READ_INI_LINE ( 'L1S', L1S, L1SSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1T') THEN           ! 36
                  CALL READ_INI_LINE ( 'L1T', L1T, L1TSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1U') THEN           ! 37
                  CALL READ_INI_LINE ( 'L1U', L1U, L1USTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1V') THEN           ! 38
                  CALL READ_INI_LINE ( 'L1V', L1V, L1VSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1W') THEN           ! 39
                  CALL READ_INI_LINE ( 'L1W', L1W, L1WSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1X') THEN           ! 40
                  CALL READ_INI_LINE ( 'L1X', L1X, L1XSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1Y') THEN           ! 41
                  CALL READ_INI_LINE ( 'L1Y', L1Y, L1YSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L1Z') THEN           ! 42
                  CALL READ_INI_LINE ( 'L1Z', L1Z, L1ZSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2A') THEN           ! 43
                  CALL READ_INI_LINE ( 'L2A', L2A, L2ASTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2B') THEN           ! 44
                  CALL READ_INI_LINE ( 'L2B', L2B, L2BSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2C') THEN           ! 45
                  CALL READ_INI_LINE ( 'L2C', L2C, L2CSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2D') THEN           ! 46
                  CALL READ_INI_LINE ( 'L2D', L2D, L2DSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2E') THEN           ! 47
                  CALL READ_INI_LINE ( 'L2E', L2E, L2ESTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2F') THEN           ! 48
                  CALL READ_INI_LINE ( 'L2F', L2F, L2FSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2G') THEN           ! 49
                  CALL READ_INI_LINE ( 'L2G', L2G, L2GSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2H') THEN           ! 50
                  CALL READ_INI_LINE ( 'L2H', L2H, L2HSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2I') THEN           ! 51
                  CALL READ_INI_LINE ( 'L2I', L2I, L2ISTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2J') THEN           ! 52
                  CALL READ_INI_LINE ( 'L2J', L2J, L2JSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2K') THEN           ! 53
                  CALL READ_INI_LINE ( 'L2K', L2K, L2KSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2L') THEN           ! 54
                  CALL READ_INI_LINE ( 'L2L', L2L, L2LSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2M') THEN           ! 55
                  CALL READ_INI_LINE ( 'L2M', L2M, L2MSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2O') THEN           ! 56
                  CALL READ_INI_LINE ( 'L2O', L2O, L2OSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2P') THEN           ! 57
                  CALL READ_INI_LINE ( 'L2P', L2P, L2PSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2Q') THEN           ! 58
                  CALL READ_INI_LINE ( 'L2Q', L2Q, L2QSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L2N') THEN           ! 59
                  CALL READ_INI_LINE ( 'L2N', L2N, L2NSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L3A') THEN           ! 60
                  CALL READ_INI_LINE ( 'L3A', L3A, L3ASTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L4A') THEN           ! 60
                  CALL READ_INI_LINE ( 'L4A', L4A, L4ASTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L4B') THEN           ! 60
                  CALL READ_INI_LINE ( 'L4B', L4B, L4BSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L4C') THEN           ! 60
                  CALL READ_INI_LINE ( 'L4C', L4C, L4CSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L4D') THEN           ! 60
                  CALL READ_INI_LINE ( 'L4D', L4D, L4DSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L5A') THEN           ! 61
                  CALL READ_INI_LINE ( 'L5A', L5A, L5ASTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'L5B') THEN           ! 62
                  CALL READ_INI_LINE ( 'L5B', L5B, L5BSTAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE IF (CARD(1:3) == 'OP2') THEN           ! 63
                  CALL READ_INI_LINE ( 'OP2', OP2, OP2STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

               ELSE
                  IF (WRT_HDR == 'N') THEN
                     CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
                     WRT_HDR = 'Y'
                  ENDIF
                  WRITE(SC1,'(1X,A79)') CARD(1:79)
                  WRITE(SC1,1502)
            !xx   WRITE(SC1, * )
                  CYCLE

               ENDIF

               DO I=1,10
                  IF (IERRFL(I) == 'Y') THEN
                     IERR = IERR + 1
                  ENDIF
               ENDDO
               DO I=1,10
                   IERRFL(I)      = 'N'
                   FLD_ERR_MSG(I) = 'N'
               ENDDO

            ENDIF

         ENDDO

         CALL FILE_CLOSE ( INI, INIFIL, 'KEEP' )

      ENDIF

! Reset ECHO

      ECHO = PERM_ECHO

! Reset close status of files, if ALL_CLOSE_STAT is KEEP

      IF (ALL_CLOSE_STAT == 'KEEP    ') THEN
         CALL SET_FILE_CLOSE_STAT ( ALL_CLOSE_STAT )
      ENDIF

! Check IERR and ask user if we should continue

      IF (IERR > 0) THEN
         DO
            WRITE(SC1,1511)
            CALL WRITE_FILNAM ( INIFIL, SC1, 1 )
            WRITE(SC1,1512)
      !xx   WRITE(SC1, * )
            READ(*,*) RESPONSE
            IF      ((RESPONSE == 'Y') .OR. (RESPONSE == 'y')) THEN
               EXIT
            ELSE IF ((RESPONSE == 'N') .OR. (RESPONSE == 'n')) THEN
               WRITE(SC1,1513)
               CALL OUTA_HERE ( 'N' )
            ELSE
               CYCLE
            ENDIF
         ENDDO
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A80)

  201 FORMAT(1X,A79)

 1501 FORMAT(' Line number ',I3,' cannot be read.')

 1502 FORMAT(' Field 1 of above line was not recognized. Line is ignored.')

 1511 FORMAT(' Errors reading:')

 1512 FORMAT(' Continue? (Y/N) (If Y, then default values will be used)')

 1513 FORMAT(' Processing terminated based on user input')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

!***********************************************************************************************************************************

      SUBROUTINE READ_INI_LINE ( FNAME, UNIT, CLOSE_STAT, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

! Reads data in fields 2 and 3 from one line of INI file for lines that contain the unit number and close status change for a file

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                  :  MOT4, MOU4, INI, OT4, OU4, SCR
      USE SCONTR, ONLY                  :  PROG_NAME

      IMPLICIT NONE

      CHARACTER( 3*BYTE), INTENT(IN)    :: FNAME           ! A 3 char designator for the file name (e.g. L1A for file LINK1A)
      CHARACTER( 8*BYTE), INTENT(INOUT) :: CLOSE_STAT      ! Close status assigned to file (or input value if read wrong)
      CHARACTER( 1*BYTE), INTENT(INOUT) :: FLD_ERR_MSG(10) ! 'Y'/'N' designator of whether error msg has been written for a field
      CHARACTER( 1*BYTE), INTENT(INOUT) :: WRT_CARD        ! 'Y', 'N' indicator if line from INI file was written to show errors
      CHARACTER( 1*BYTE), INTENT(INOUT) :: WRT_HDR         ! 'Y', 'N' indicator if header msg written when there are errors
      CHARACTER( 8*BYTE)                :: CHAR8           ! Close status read from field 3 of INI line

      INTEGER(LONG)     , INTENT(INOUT) :: UNIT            ! File unit number assigned to file (or input value if read wrong)
      INTEGER(LONG)                     :: II              ! DO loop index
      INTEGER(LONG)                     :: IERR     = 0    ! Local error count
      INTEGER(LONG)                     :: INT4            ! File unit number read from field 2 of INI line

!***********************************************************************************************************************************
      CALL LEFT_ADJ_BDFLD0 (JCARD_08(2))
      IF (JCARD_08(2) /= '        ') THEN                  ! Read file unit number
         CALL I4FLD0 ( JCARD_08(2), JF(2), INT4, WRT_HDR, WRT_CARD, FLD_ERR_MSG )
         IF (IERRFL(2) == 'N') THEN
            IF ((INT4 /= INI) .AND. (INT4 < SCR(1)) .AND. (INT4 > 0)) THEN
               DO II=1,MOU4                                ! make sure unit number is not an OU4 unit number
                  IF (INT4 /= OU4(II)) THEN
                     CYCLE
                  ELSE
                     IERR = IERR + 1
                     EXIT
                  ENDIF
               ENDDO
               DO II=1,MOT4                                ! make sure unit number is not an OT4 unit number
                  IF (INT4 /= OT4(II)) THEN
                     CYCLE
                  ELSE
                     IERR = IERR + 1
                     EXIT
                  ENDIF
               ENDDO
               IF (IERR == 0) THEN
                  UNIT = INT4
               ELSE
                  IF (WRT_HDR  == 'N') THEN
                     CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
                     WRT_HDR  = 'Y'
                  ENDIF
                  IF (WRT_CARD == 'N') THEN
                     WRITE(SC1,'(1X,A79)') CARD(1:79)
                     WRT_CARD = 'Y'
                  ENDIF
                  WRITE(SC1,1510) INT4, FNAME
                  FLD_ERR_MSG(2) = 'Y'
               ENDIF
            ELSE
               IF (WRT_HDR  == 'N') THEN
                  CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
                  WRT_HDR  = 'Y'
               ENDIF
               IF (WRT_CARD == 'N') THEN
                  WRITE(SC1,'(1X,A79)') CARD(1:79)
                  WRT_CARD = 'Y'
               ENDIF
               WRITE(SC1,1510) INT4, FNAME
               FLD_ERR_MSG(2) = 'Y'
            ENDIF
         ENDIF
      ENDIF

      CALL LEFT_ADJ_BDFLD0 (JCARD_08(3))
      IF (JCARD_08(3) /= '        ') THEN                  ! Read file close status
         CALL C8FLD0 ( JCARD_08(3), JF(3), CHAR8 )
         CALL LEFT_ADJ_BDFLD0 ( CHAR8 )
         IF ((CHAR8 == 'KEEP    ') .OR. (CHAR8 == 'DELETE  ')) THEN
            CLOSE_STAT = CHAR8
         ELSE
            IERRFL(3) = 'Y'
            IF (WRT_HDR  == 'N') THEN
               CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
               WRT_HDR  = 'Y'
            ENDIF
            IF (WRT_CARD == 'N') THEN
               WRITE(SC1,'(1X,A79)') CARD(1:79)
               WRT_CARD = 'Y'
            ENDIF
            WRITE(SC1,1511) CHAR8, FNAME
            FLD_ERR_MSG(3) = 'Y'
         ENDIF
      ENDIF

      CALL BD_IMBEDDED_BLANK0  ( JCARD_08, 2, 3, 0, 0, 0, 0, 0, 0, WRT_HDR, WRT_CARD, FLD_ERR_MSG )
      CALL CARD_FLDS_NOT_BLANK0( JCARD_08, 0, 0, 4, 5, 6, 7, 8, 9, WRT_HDR, WRT_CARD )
      CALL CRDERR0 ( CARD, FLD_ERR_MSG )

      IF (WRT_HDR == 'Y') THEN
   !xx   WRITE(SC1, * )
      ENDIF

!***********************************************************************************************************************************
 1510 FORMAT(' Cannot use ',I3,' for file unit ',A)

 1511 FORMAT(' Cannot use ',A,' for close status for file unit ',A)

!***********************************************************************************************************************************

      END SUBROUTINE READ_INI_LINE

! ##################################################################################################################################

      SUBROUTINE I4FLD0 ( JCARDI_08, IFLD, I4INP, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

! Reads 8 column field of INTEGER*4 data

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                 :  IERRFL

      IMPLICIT NONE

      CHARACTER( 8*BYTE), INTENT(IN)   :: JCARDI_08        ! The field of 8 characters to read
      CHARACTER( 1*BYTE), INTENT(INOUT):: FLD_ERR_MSG(10)  ! 'Y'/'N' designator of whether error msg has been written for a field
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_CARD         ! 'Y', 'N' indicator if line from INI file was written to show errors
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_HDR          ! 'Y', 'N' indicator if header msg written when there are errors
      CHARACTER( 1*BYTE)               :: DEC_PT           ! 'Y'/'N' indicator of whether a decimal point was founr in JCARDI_08

      INTEGER(LONG), INTENT(IN)        :: IFLD             ! Field (2 - 9) of a Bulk Data card to read
      INTEGER(LONG), INTENT(OUT)       :: I4INP            ! The 4 byte integer value read
      INTEGER(LONG)                    :: I                ! DO loop index
      INTEGER(LONG)                    :: IOCHK            ! IOSTAT error value from READ

! **********************************************************************************************************************************
! Initialize outputs

      I4INP = 0

      READ(JCARDI_08,'(I8)',IOSTAT=IOCHK) I4INP

      IF      (IOCHK < 0) THEN                             ! EOF/EOR during read

         IERRFL(IFLD) = 'Y'

      ELSE IF (IOCHK == 0) THEN                            ! READ was OK

         IERRFL(IFLD) = 'N'

      ELSE IF (IOCHK > 0) THEN                             ! Error during READ

         IERRFL(IFLD) = 'Y'

      ENDIF

! Scan to make sure there was not a decimal point.

      IF (JCARDI_08 /= '        ') THEN
         DEC_PT = 'N'
         DO I=1,8
            IF (JCARDI_08(I:I) == '.') THEN
               DEC_PT = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (DEC_PT == 'Y') THEN
            IF (WRT_HDR  == 'N') THEN
               CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
               WRT_HDR  = 'Y'
            ENDIF
            IF (WRT_CARD == 'N') THEN
               WRITE(SC1,'(1X,A79)') CARD(1:79)
               WRT_CARD = 'Y'
            ENDIF
            IERRFL(IFLD) = 'Y'
            WRITE(SC1,1500) IFLD
            FLD_ERR_MSG = 'Y'
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1500 FORMAT(' A decimal pt was found in what is supposed to be an integer number in field ',I2)

! **********************************************************************************************************************************

      END SUBROUTINE I4FLD0

! ##################################################################################################################################

      SUBROUTINE BD_IMBEDDED_BLANK0 ( JCARD_08, CF2, CF3, CF4, CF5, CF6, CF7, CF8, CF9, WRT_HDR, WRT_CARD, FLD_ERR_MSG )

! Prepares message when some fields of a Bulk data card imbedded blanks when they should not

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                 :  PROG_NAME

      IMPLICIT NONE

      CHARACTER( 8*BYTE), INTENT(IN)   :: JCARD_08(10)     ! The 10 fields of 8 characters making up CARD
      CHARACTER( 1*BYTE), INTENT(INOUT):: FLD_ERR_MSG(10)  ! 'Y'/'N' designator of whether error msg has been written for a field
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_CARD         ! 'Y'/'N' indicator if line from INI file was written to show errors
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_HDR          ! 'Y'/'N' indicator if header msg written when there are errors
      CHARACTER( 1*BYTE)               :: ERRORS           ! 'Y'/'N' indicator of errors
      CHARACTER( 1*BYTE)               :: FOUND_DATA       ! 'Y'/'N' indicator of data found in a Bulk Data card field
      CHARACTER( 1*BYTE)               :: IMB_BLANK(2:9)   ! 'Y'/'N' indicator of whether fields 2-9 have imbedded blanks

      INTEGER(LONG), INTENT(IN)        :: CF2              ! = 2 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF3              ! = 3 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF4              ! = 4 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF5              ! = 5 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF6              ! = 6 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF7              ! = 7 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF8              ! = 8 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)        :: CF9              ! = 9 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG)                    :: CHK_FLD(2:9)     ! Array containing CF2 through CF9
      INTEGER(LONG)                    :: I,J              ! DO loop indices
      INTEGER(LONG)                    :: NUMBER(2:9)      ! Number of imbedded blanks found in a Bulk Data card field

! **********************************************************************************************************************************

! Load CF2 through CF9 into array CHK_FLD

      CHK_FLD(2) = CF2
      CHK_FLD(3) = CF3
      CHK_FLD(4) = CF4
      CHK_FLD(5) = CF5
      CHK_FLD(6) = CF6
      CHK_FLD(7) = CF7
      CHK_FLD(8) = CF8
      CHK_FLD(9) = CF9

! Initialize IMB_BLANK array

      DO I=2,9
         IMB_BLANK(I) = 'N'
      ENDDO

! Check fields for any imbedded blanks and set error if so

      ERRORS = 'N'
      DO I=2,9
         NUMBER(I) = 0
         IF      (CHK_FLD(I) == 0) THEN                    ! We don't want to check field I
            CYCLE
         ELSE IF (CHK_FLD(I) == I) THEN                    ! We do want to check field I
            FOUND_DATA   = 'N'
            IMB_BLANK(I) = 'N'
            DO J=8,1,-1
               IF(JCARD_08(I)(J:J) /= ' ') THEN
                  FOUND_DATA = 'Y'
               ELSE
                  IF(FOUND_DATA == 'Y') THEN
                     IMB_BLANK(I) = 'Y'
                     NUMBER(I) = NUMBER(I) + 1
                     ERRORS = 'Y'
                  ELSE
                     CYCLE
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (ERRORS == 'Y') THEN
         IF (WRT_HDR  == 'N') THEN
            CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
            WRT_HDR  = 'Y'
         ENDIF
         IF (WRT_CARD == 'N') THEN
            WRITE(SC1,'(1X,A79)') CARD(1:79)
            WRT_CARD = 'Y'
         ENDIF
      ENDIF

! Write error if fields checked have imbedded blanks

      DO I=2,9
         IF ((IMB_BLANK(I) == 'Y') .AND. (FLD_ERR_MSG(I) == 'N')) THEN
            IERRFL(I) = 'Y'
            WRITE(SC1,1722) NUMBER(I),I
            FLD_ERR_MSG(I) = 'Y'
            CYCLE
         ENDIF
      ENDDO

      RETURN

! **********************************************************************************************************************************
 1722 FORMAT(' There were ',I2,' imbedded blanks (not allowed) found in field ',I2)

! **********************************************************************************************************************************

      END SUBROUTINE BD_IMBEDDED_BLANK0

! ##################################################################################################################################

      SUBROUTINE CARD_FLDS_NOT_BLANK0 ( JCARD_08, FLD2, FLD3, FLD4, FLD5, FLD6, FLD7, FLD8, FLD9, WRT_HDR, WRT_CARD )

! Prepares message when some fields of a Bulk data card that should be blank, aren't

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                 :  PROG_NAME

      IMPLICIT NONE

      CHARACTER( 8*BYTE), INTENT(IN)   :: JCARD_08(10)      ! The 10 fields of 8 characters making up CARD
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_CARD          ! 'Y', 'N' indicator if line from INI file was written to show errors
      CHARACTER( 1*BYTE), INTENT(INOUT):: WRT_HDR           ! 'Y', 'N' indicator if header msg written when there are errors
      CHARACTER( 1*BYTE)               :: COMMENT           ! 'Y' or 'N' depending on whether non-blank fields are a comment
      CHARACTER( 8*BYTE)               :: MSSG8             ! Message with all fields that are not blank that should be blank
      CHARACTER( 1*BYTE)               :: MSSG1             ! Message that has the field number in it

      INTEGER(LONG), INTENT(IN)        :: FLD2              ! Refers to field 2 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD3              ! Refers to field 3 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD4              ! Refers to field 4 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD5              ! Refers to field 5 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD6              ! Refers to field 6 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD7              ! Refers to field 7 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD8              ! Refers to field 8 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)        :: FLD9              ! Refers to field 9 of a Bulk Data card. If /= 0, then check this field
      INTEGER(LONG)                    :: ALL_FLDS(2:9)     ! Array of the FLDi (2 through 9)
      INTEGER(LONG)                    :: I,J               ! Do loop indices

! **********************************************************************************************************************************
! Set ALL_FLDS

      ALL_FLDS(2) = FLD2
      ALL_FLDS(3) = FLD3
      ALL_FLDS(4) = FLD4
      ALL_FLDS(5) = FLD5
      ALL_FLDS(6) = FLD6
      ALL_FLDS(7) = FLD7
      ALL_FLDS(8) = FLD8
      ALL_FLDS(9) = FLD9

! Scan for 1st non-blank field that should not have data and see if the 1st non-blank character is '$'.
! If it is, then the non-blank fields have a comment, so return

      COMMENT = 'N'
i_do: DO I=2,9
         IF (ALL_FLDS(I) /= I) THEN                        ! CYCLE until we find 1st field to check
            CYCLE
         ELSE                                              ! Found a field that should be blank (unless it begins a comment)
j_do:       DO J=1,8                                       ! CYCLE through 8 chars of this field to see if 1st non-blank char is '$'
               IF (JCARD_08(I)(J:J) == ' ') THEN
                  CYCLE j_do
               ELSE
                  IF (JCARD_08(I)(J:J) == '$') THEN           ! Found a '$' (with blanks preceeding it) so we found a comment
                     COMMENT = 'Y'
                     EXIT i_do                             ! Set COMMENT = 'Y' and EXIT outer loop
                  ENDIF
               ENDIF
            ENDDO j_do
         ENDIF
      ENDDO i_do

! Write message if fields not blank that should be.

      IF (COMMENT == 'N') THEN

         MSSG8 = '        '
         IF ((JCARD_08(2) /= '        ') .AND. (FLD2 == 2)) THEN
            MSSG1 = '2'
            MSSG8 =                  MSSG1(1:1) // MSSG8(2:8)
         ENDIF
         IF ((JCARD_08(3) /= '        ') .AND. (FLD3 == 3)) THEN
            MSSG1 = '3'
            MSSG8 = MSSG8(1:1) // MSSG1(1:1) // MSSG8(3:8)
         ENDIF
         IF ((JCARD_08(4) /= '        ') .AND. (FLD4 == 4)) THEN
            MSSG1 = '4'
            MSSG8 = MSSG8(1:2) // MSSG1(1:1) // MSSG8(4:8)
         ENDIF
         IF ((JCARD_08(5) /= '        ') .AND. (FLD5 == 5)) THEN
            MSSG1 = '5'
            MSSG8 = MSSG8(1:3) // MSSG1(1:1) // MSSG8(5:8)
         ENDIF
         IF ((JCARD_08(6) /= '        ') .AND. (FLD6 == 6)) THEN
            MSSG1 = '6'
            MSSG8 = MSSG8(1:4) // MSSG1(1:1) // MSSG8(6:8)
         ENDIF
         IF ((JCARD_08(7) /= '        ') .AND. (FLD7 == 7)) THEN
            MSSG1 = '7'
            MSSG8 = MSSG8(1:5) // MSSG1(1:1) // MSSG8(7:8)
         ENDIF
         IF ((JCARD_08(8) /= '        ') .AND. (FLD8 == 8)) THEN
            MSSG1 = '8'
            MSSG8 = MSSG8(1:6) // MSSG1(1:1) // MSSG8(8:8)
         ENDIF
         IF ((JCARD_08(9) /= '        ') .AND. (FLD9 == 9)) THEN
            MSSG1 = '9'
            MSSG8 = MSSG8(1:7) // MSSG1(1:1)
         ENDIF

         IF (MSSG8 /= '        ') THEN
            IF (WRT_HDR  == 'N') THEN
               CALL WRITE_INIFIL_MSG ( '2', INIFIL_NAME_LEN )
               WRT_HDR  = 'Y'
            ENDIF
            IF (WRT_CARD == 'N') THEN
               WRITE(SC1,'(1X,A79)') CARD(1:79)
               WRT_CARD = 'Y'
            ENDIF
            WRITE(SC1,1726) MSSG8
         ENDIF

      ENDIF

      RETURN

! **********************************************************************************************************************************
 1726 FORMAT(' Field(s) ',A8,' should be blank and are ignored')

! **********************************************************************************************************************************

      END SUBROUTINE CARD_FLDS_NOT_BLANK0

! ##################################################################################################################################

      SUBROUTINE CRDERR0 ( CARD, FLD_ERR_MSG )

! Prints Bulk Data card errors and warnings

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  IERRFL

      IMPLICIT NONE

      CHARACTER(LEN=*)  , INTENT(IN)  :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE), INTENT(IN)  :: FLD_ERR_MSG(10)   ! 'Y'/'N' designator of whether error msg has been written for a field
      CHARACTER( 1*BYTE)              :: CPN               ! = 'Y' if IERRFL is 'Y' for any Bulk Data card field

      INTEGER(LONG)                   :: I                 ! DO loop index

! **********************************************************************************************************************************
      CPN = 'N'
      DO I=1,10
         IF (IERRFL(I) == 'Y') THEN
            CPN = 'Y'
         ENDIF
      ENDDO

      IF (CPN == 'Y') THEN
         WRITE(SC1,1501) CARD(1:79)
         DO I=1,10
            IF ((IERRFL(I) == 'Y') .AND. (FLD_ERR_MSG(I) == 'N')) THEN
               WRITE(SC1,1502) I
            ENDIF
         ENDDO
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1501 FORMAT(1X,A)

 1502 FORMAT(' Format error in field',I3)

! **********************************************************************************************************************************

      END SUBROUTINE CRDERR0

! ##################################################################################################################################

      SUBROUTINE LEFT_ADJ_BDFLD0 ( CHR8_FLD )

! Shifts an 8 character string so that it is left adjusted

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(8*BYTE), INTENT(INOUT):: CHR8_FLD          ! Char field to left adjust and return
      CHARACTER(8*BYTE)               :: TCHR8_FLD         ! Temporary char field

      INTEGER(LONG)                   :: I                 ! DO loop index

! **********************************************************************************************************************************
      IF (CHR8_FLD(1:1) == ' ') THEN                       ! We need to shift:

         TCHR8_FLD(1:8) = CHR8_FLD(1:8)                    ! Set temporary field to CHR8_FLD

         DO I = 2,8                                        ! Perform shift
            IF (CHR8_FLD(I:I) /= ' ') THEN
               TCHR8_FLD(1:) = CHR8_FLD(I:)
               EXIT
            ENDIF
         ENDDO

         CHR8_FLD(1:8) = TCHR8_FLD(1:8)                    ! Reset CHR*_FLD and return

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE LEFT_ADJ_BDFLD0

! ##################################################################################################################################

      SUBROUTINE C8FLD0 ( JCARDI_08, IFLD, C8INP )

! Reads a field of CHARACTER data that can be 1 to 8 chars in length

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR

      IMPLICIT NONE

      CHARACTER( 8*BYTE), INTENT(IN)  :: JCARDI_08         ! The field of 8 characters to read
      CHARACTER(8*BYTE) , INTENT(OUT) :: C8INP             ! The character variable to read

      INTEGER(LONG), INTENT(IN)       :: IFLD              ! Field (2 - 9) of a Bulk Data card to read
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ

! **********************************************************************************************************************************
      READ(JCARDI_08,'(A8)',IOSTAT=IOCHK) C8INP

      IF      (IOCHK < 0) THEN                             ! EOF/EOR during read

         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1

      ELSE IF (IOCHK == 0) THEN                            ! READ was OK

         IERRFL(IFLD) = 'N'

      ELSE IF (IOCHK > 0) THEN                             ! Error during READ

         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1

      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE C8FLD0

! ##################################################################################################################################

      SUBROUTINE WRITE_INIFIL_MSG ( WHICH_MSG, INIFIL_NAME_LEN )

! Writes message to screen about INIFIL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE iount1, ONLY                :  FILE_NAM_MAXLEN

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: NUM_MSGS    = 2     ! Number of messages to write

      CHARACTER( 1*BYTE), INTENT(IN)  :: WHICH_MSG           ! Which message to write to screen
      CHARACTER(79*BYTE)              :: MSG(NUM_MSGS)       ! Message to write

      INTEGER(LONG)                   :: INIFIL_NAME_LEN     ! Length of INI file name (incl path)

! **********************************************************************************************************************************
      DO I=1,NUM_MSGS
         MSG(I)(1:) = ' '
      ENDDO

      WRITE(SC1,1001) TRIM(INIFIL)

      IF      (WHICH_MSG == '1') THEN
         MSG(1) = 'does not exist. Internal defaults will be used instead.'
         MSG(2) = 'If you want to use the INI file, see the MYSTRAN Installation and Run Manual'
         DO I=1,NUM_MSGS
            WRITE(SC1,2000) MSG(I)
         ENDDO
      ELSE IF (WHICH_MSG == '2') THEN
         MSG(1) = 'has errors (file line with errors is displayed followed by the error messages):'
         WRITE(SC1,2000) MSG(1)
      ENDIF

! **********************************************************************************************************************************
 1001 FORMAT(/,/' Optional Initialization file: ',/,1X,A)

 2000 FORMAT(1X,A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_INIFIL_MSG

      END SUBROUTINE READ_INI


! ##################################################################################################################################

      SUBROUTINE GET_INI_FILNAM  ( MYSTRAN_DIR, MYSTRAN_DIR_LEN, INIFIL_NAME_LEN )

! Gets name (incl path) of the MYSTRAN.INI initialization file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, INIFIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, PROG_NAME
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(PROG_NAME)+4) :: FILNAM            ! File name for MYSTRAN initialization file (not including path)

                                                           ! Directory where MYSTRAN executable (and INI file) exist
      CHARACTER(FILE_NAM_MAXLEN*BYTE), INTENT(IN):: MYSTRAN_DIR
      CHARACTER :: OS_DIR_SEP = '/'                        ! Directory separator

      INTEGER(LONG), INTENT(IN)       :: MYSTRAN_DIR_LEN   ! Length of MYSTRAN_DIR (not including trailing blanks)
      INTEGER(LONG), INTENT(OUT)      :: INIFIL_NAME_LEN   ! Length of INI file name (incl path)

! **********************************************************************************************************************************
      FILNAM( 1:LEN(PROG_NAME)) = PROG_NAME
      FILNAM(LEN(PROG_NAME)+1:LEN(PROG_NAME)+4) = '.INI'

      IF (MYSTRAN_DIR_LEN > 0) THEN
         INIFIL(1:MYSTRAN_DIR_LEN)  = MYSTRAN_DIR
         INIFIL(MYSTRAN_DIR_LEN+1:MYSTRAN_DIR_LEN+2)  = OS_DIR_SEP
         INIFIL(MYSTRAN_DIR_LEN+2:) = FILNAM
         INIFIL_NAME_LEN = MYSTRAN_DIR_LEN +1 + LEN(PROG_NAME) + 4
      ELSE
         INIFIL = FILNAM
         INIFIL_NAME_LEN = LEN(PROG_NAME) + 4
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_INI_FILNAM


      SUBROUTINE GET_MYSTRAN_DIR ( MYSTRAN_DIR, MYSTRAN_DIR_LEN )

! Gets the environment variable MYSTRAN_DIR that tells Windows where the MYSTRAN executable is located. The user must have set this
! environment variable on their computer

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN

      IMPLICIT NONE

      CHARACTER(FILE_NAM_MAXLEN*BYTE), INTENT(OUT) :: MYSTRAN_DIR       ! Directory where executable (and INI file) exist

      INTEGER(LONG), INTENT(OUT)                   :: MYSTRAN_DIR_LEN   ! Length of MYSTRAN_DIR (not including trailing blanks)
      INTEGER(LONG)                                :: I                 ! DO loop index

      INTRINSIC                                    :: GETENV

! **********************************************************************************************************************************
      CALL GETENV ( 'MYSTRAN_directory', MYSTRAN_DIR )
      MYSTRAN_DIR_LEN = FILE_NAM_MAXLEN
      DO I=FILE_NAM_MAXLEN,1,-1
         IF (MYSTRAN_DIR(I:I) /= ' ') THEN
            EXIT
         ELSE
            MYSTRAN_DIR_LEN = MYSTRAN_DIR_LEN - 1
            CYCLE
         ENDIF
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE GET_MYSTRAN_DIR



   END MODULE INI_FILE
