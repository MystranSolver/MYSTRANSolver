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

   MODULE STARTUP_FILES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MYSTRAN_FILES, PROCESS_INCLUDE_FILES, READ_CL, READ_INPUT_FILE_NAME

   CONTAINS

      SUBROUTINE MYSTRAN_FILES ( START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC)

! Sets all MYSTRAN file names. Opens all files and closes and deletes them so that no confusion about files if MYSTRAN aborts
! Reopen BUG, ERR, F06

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, MOT4, MOU4, WRT_BUG, WRT_ERR, LEN_INPUT_FNAME,                  &
                                         LEN_RESTART_FNAME, RESTART_FILNAM

      USE IOUNT1, ONLY                :  OU4_EXT, OT4_EXT

      USE IOUNT1, ONLY                :  BUG,     EIN,     ENF,     ERR,     F06,     IN0,     SC1,                                &
                                         SEQ,     SPC,                                                                             &
                                         L1A,     L1B,     L1C,     L1D,     L1E,     L1F,     L1G,     L1H,     L1I,     L1J,     &
                                         L1K,     L1L,     L1M,     L1N,     L1O,     L1P,     L1Q,     L1R,     L1S,     L1T,     &
                                         L1U,     L1V,     L1W,     L1X,     L1Y,     L1Z,                                         &
                                         L2A,     L2B,     L2C,     L2D,     L2E,     L2F,     L2G,     L2H,     L2I,     L2J,     &
                                         L2K,     L2L,     L2M,     L2N,     L2O,     L2P,     L2Q,     L2R,     L2S,     L2T,     &
                                         L3A,     L4A,     L4B,     L4C,     L4D,     L5A,     L5B,                                &
                                         NEU,     F21,     F22,     F23,     F24,     F25,     OP2,     OT4,     OU4

      USE IOUNT1, ONLY                :  BUGFIL,  EINFIL,  ENFFIL,  ERRFIL,  F06FIL,  IN0FIL,  INFILE,                             &
                                         OT4FIL,  SEQFIL,  SPCFIL,                                                                 &
                                         LINK1A,  LINK1B,  LINK1C,  LINK1D,  LINK1E,  LINK1F,  LINK1G,  LINK1H,  LINK1I,  LINK1J,  &
                                         LINK1K,  LINK1L,  LINK1M,  LINK1N,  LINK1O,  LINK1P,  LINK1Q,  LINK1R,  LINK1S,  LINK1T,  &
                                         LINK1U,  LINK1V,  LINK1W,  LINK1X,  LINK1Y,  LINK1Z,                                      &
                                         LINK2A,  LINK2B,  LINK2C,  LINK2D,  LINK2E,  LINK2F,  LINK2G,  LINK2H,  LINK2I,  LINK2J,  &
                                         LINK2K,  LINK2L,  LINK2M,  LINK2N,  LINK2O,  LINK2P,  LINK2Q,  LINK2R,  LINK2S,  LINK2T,  &
                                         LINK3A,  LINK4A,  LINK4B,  LINK4C,  LINK4D,  LINK5A,  LINK5B,                             &
                                         NEUFIL,  F21FIL,  F22FIL,  F23FIL,  F24FIL,  F25FIL,  OP2FIL,  OT4FIL,  OU4FIL

      USE IOUNT1, ONLY                :  BUG_MSG, EIN_MSG, ENF_MSG, ERR_MSG, F06_MSG, IN0_MSG, OT4_MSG,                            &
                                         SEQ_MSG, L1A_MSG, L1B_MSG, L1C_MSG, L1D_MSG, L1E_MSG, L1F_MSG, L1G_MSG, L1H_MSG, L1I_MSG, &
                                         L1J_MSG, L1K_MSG, L1L_MSG, L1M_MSG, L1N_MSG, L1O_MSG, L1P_MSG, L1Q_MSG, L1R_MSG, L1S_MSG, &
                                         L1T_MSG, L1U_MSG, L1V_MSG, L1W_MSG, L1X_MSG, L1Y_MSG, L1Z_MSG,                            &
                                         L2A_MSG, L2B_MSG, L2C_MSG, L2D_MSG, L2E_MSG, L2F_MSG, L2G_MSG, L2H_MSG, L2I_MSG, L2J_MSG, &
                                         L2K_MSG, L2L_MSG, L2M_MSG, L2N_MSG, L2O_MSG, L2P_MSG, L2Q_MSG, L2R_MSG, L2S_MSG, L2T_MSG, &
                                         L3A_MSG, L4A_MSG, L4B_MSG, L4C_MSG, L4D_MSG, L5A_MSG, L5B_MSG,                            &
                                         NEU_MSG, F21_MSG, F22_MSG, F23_MSG, F24_MSG, F25_MSG, OP2_MSG, OT4_MSG, OU4_MSG, SPC_MSG

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, RESTART
      USE TIMDAT, ONLY                :  TSEC, stime

      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST        ! T/F depending on whether a file exists

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME     = 'MYSTRAN_FILES'
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: FILNAM            ! Name of the input file or restart file


      INTEGER(LONG), INTENT(IN)       :: START_HOUR        ! The hour     when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_MINUTE      ! The minute   when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_SEC         ! The second   when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_SFRAC       ! The sec frac when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_YEAR        ! The year     when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_MONTH       ! The month    when MYSTRAN started.
      INTEGER(LONG), INTENT(IN)       :: START_DAY         ! The day      when MYSTRAN started.
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I1                ! Filename (less extension) length
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR

! **********************************************************************************************************************************
! Default units for writing errors the screen (until LINK1A is read) and set filename length

      OUNT(1) = SC1
      OUNT(2) = SC1

      IF (RESTART == 'Y') THEN
         I1     = LEN_RESTART_FNAME
         FILNAM = RESTART_FILNAM
      ELSE
         I1     = LEN_INPUT_FNAME
         FILNAM = INFILE
      ENDIF

! Form file names then open them and close and delete them. In this way, we can get rid of all of these files before
! we begin and start anew.

! Formatted files. Note: for ERR, F06, BUG reopen them after deleting any old version and write STIME



! **********************************************************************************************************************************
      EINFIL(1:I1)  = FILNAM(1:I1)
      EINFIL(I1+1:) = 'EIN'

      ENFFIL(1:I1)  = FILNAM(1:I1)
      ENFFIL(I1+1:) = 'ENF'

      F06FIL(1:I1)  = FILNAM(1:I1)
      F06FIL(I1+1:) = 'F06'
      IF (F06 /= SC1) THEN
         INQUIRE ( FILE=F06FIL, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN
               CALL FILE_OPEN ( F06, F06FIL, OUNT,'OLD    ', F06_MSG,'NEITHER'    ,'FORMATTED','READWRITE','APPEND','N','N')
               WRITE(F06,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
               CALL FILE_CLOSE ( F06, F06FIL,'KEEP')
               CALL FILE_OPEN ( F06, F06FIL, OUNT,'OLD    ', F06_MSG,'NEITHER'    ,'FORMATTED','READWRITE','APPEND','N','N')
            ELSE
               CALL FILE_OPEN ( F06, F06FIL, OUNT,'UNKNOWN', F06_MSG,'NEITHER'    ,'FORMATTED','READWRITE','REWIND','N','N')
               CALL FILE_CLOSE ( F06, F06FIL,'DELETE')
               CALL FILE_OPEN ( F06, F06FIL, OUNT,'NEW'    , F06_MSG,'WRITE_STIME','FORMATTED','WRITE'    ,'REWIND','Y','Y')
               WRITE(F06,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ELSE
            CALL FILE_OPEN ( F06, F06FIL, OUNT,'NEW', F06_MSG,'WRITE_STIME','FORMATTED','WRITE','REWIND','Y','Y')
            IF(RESTART == 'Y') THEN
               WRITE(F06,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ELSE
               WRITE(F06,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ENDIF
      ENDIF

      IN0FIL(1:I1)  = FILNAM(1:I1)
      IN0FIL(I1+1:) = 'IN0'
      IF (IN0 /= SC1) THEN
         INQUIRE ( FILE=IN0FIL, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN
               CALL FILE_OPEN ( IN0, IN0FIL, OUNT,'OLD    ', IN0_MSG,'NEITHER','FORMATTED','WRITE','APPEND','N','N')
               CALL FILE_CLOSE ( IN0, IN0FIL,'KEEP')
               CALL FILE_OPEN ( IN0, IN0FIL, OUNT,'OLD    ', IN0_MSG,'NEITHER','FORMATTED','WRITE','APPEND','N','N')
            ELSE
               CALL FILE_OPEN ( IN0, IN0FIL, OUNT,'REPLACE', IN0_MSG,'NEITHER','FORMATTED','WRITE','REWIND','N','N')
               CALL FILE_CLOSE ( IN0, IN0FIL,'DELETE')
               CALL FILE_OPEN ( IN0, IN0FIL, OUNT,'NEW', IN0_MSG,'NEITHER','FORMATTED','WRITE','REWIND','N','N')
            ENDIF
         ELSE
            CALL FILE_OPEN ( IN0, IN0FIL, OUNT,'NEW', IN0_MSG,'NEITHER','FORMATTED','WRITE','REWIND','N','N')
         ENDIF
      ENDIF

      BUGFIL(1:I1)  = FILNAM(1:I1)
      BUGFIL(I1+1:) = 'BUG'
      IF (BUG /= SC1) THEN
         INQUIRE ( FILE=BUGFIL, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN
               CALL FILE_OPEN ( BUG, BUGFIL, OUNT,'OLD    ', BUG_MSG,'NEITHER','FORMATTED','READWRITE','APPEND','N','N')
               WRITE(BUG,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
               CALL FILE_CLOSE ( BUG, BUGFIL,'KEEP')
               CALL FILE_OPEN ( BUG, BUGFIL, OUNT,'OLD    ', BUG_MSG,'NEITHER','FORMATTED','READWRITE','APPEND','N','N')
            ELSE
               CALL FILE_OPEN ( BUG, BUGFIL, OUNT,'REPLACE', BUG_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
               CALL FILE_CLOSE ( BUG, BUGFIL,'DELETE')
               CALL FILE_OPEN ( BUG, BUGFIL, OUNT,'NEW', BUG_MSG,'WRITE_STIME','FORMATTED','WRITE','REWIND','Y','Y')
               WRITE(BUG,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ELSE
            CALL FILE_OPEN ( BUG, BUGFIL, OUNT,'NEW', BUG_MSG,'WRITE_STIME','FORMATTED','WRITE','REWIND','Y','Y')
            IF(RESTART == 'Y') THEN
               WRITE(BUG,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ELSE
               WRITE(BUG,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ENDIF
      ENDIF

      ERRFIL(1:I1)  = FILNAM(1:I1)
      ERRFIL(I1+1:) = 'ERR'
      IF (ERR /= SC1) THEN
         INQUIRE ( FILE=ERRFIL, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN
               CALL FILE_OPEN ( ERR, ERRFIL, OUNT,'OLD    ', ERR_MSG,'NEITHER','FORMATTED','READWRITE','APPEND','N','N')
               WRITE(ERR,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
               CALL FILE_CLOSE ( ERR, ERRFIL,'KEEP')
               CALL FILE_OPEN ( ERR, ERRFIL, OUNT,'OLD    ', ERR_MSG,'NEITHER','FORMATTED','READWRITE','APPEND','N','N')
            ELSE
               CALL FILE_OPEN ( ERR, ERRFIL, OUNT,'REPLACE', ERR_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
               CALL FILE_CLOSE ( ERR, ERRFIL,'DELETE')
               CALL FILE_OPEN ( ERR, ERRFIL, OUNT,'NEW', ERR_MSG,'WRITE_STIME','FORMATTED','WRITE','REWIND','Y','Y')
               WRITE(ERR,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ELSE
            CALL FILE_OPEN ( ERR, ERRFIL, OUNT,'NEW', ERR_MSG,'WRITE_STIME','FORMATTED','WRITE','REWIND','Y','Y')
            IF(RESTART == 'Y') THEN
               WRITE(ERR,170) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ELSE
               WRITE(ERR,150) START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC, INFILE
            ENDIF
         ENDIF
      ENDIF

      LINK1A(1:I1)  = FILNAM(1:I1)
      LINK1A(I1+1:) = 'L1A'
      INQUIRE ( FILE=LINK1A, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1A, LINK1A, OUNT,'OLD    ', L1A_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1A, LINK1A,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1A, LINK1A, OUNT,'REPLACE', L1A_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1A, LINK1A,'DELETE')
         ENDIF
      ENDIF

      NEUFIL(1:I1)  = FILNAM(1:I1)
      NEUFIL(I1+1:) = 'NEU'
      INQUIRE ( FILE=NEUFIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( NEU, NEUFIL, OUNT,'OLD    ', NEU_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( NEU, NEUFIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( NEU, NEUFIL, OUNT,'REPLACE', NEU_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( NEU, NEUFIL,'DELETE')
         ENDIF
      ENDIF

      SEQFIL(1:I1)  = FILNAM(1:I1)
      SEQFIL(I1+1:) = 'SEQ'
      INQUIRE ( FILE=SEQFIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( SEQ, SEQFIL, OUNT,'OLD    ', SEQ_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( SEQ, SEQFIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( SEQ, SEQFIL, OUNT,'REPLACE', SEQ_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( SEQ, SEQFIL,'DELETE')
         ENDIF
      ENDIF

      SPCFIL(1:I1)  = FILNAM(1:I1)
      SPCFIL(I1+1:) = 'SPC'
      INQUIRE ( FILE=SPCFIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( SPC, SPCFIL, OUNT,'OLD    ', SPC_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( SPC, SPCFIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( SPC, SPCFIL, OUNT,'REPLACE', SPC_MSG,'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
            CALL FILE_CLOSE ( SPC, SPCFIL,'DELETE')
         ENDIF
      ENDIF

      DO I=1,MOT4
         OT4FIL(I)(1:I1)  = FILNAM(1:I1)
         OT4FIL(I)(I1+1:) = OT4_EXT(I)
         INQUIRE ( FILE=OT4FIL(I), EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN
               CALL FILE_OPEN (OT4(I), OT4FIL(I), OUNT,'OLD    ', OT4_MSG(I),'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
               CALL FILE_CLOSE (OT4(I), OT4FIL(I),'KEEP')
            ELSE
               CALL FILE_OPEN (OT4(I), OT4FIL(I), OUNT,'REPLACE', OT4_MSG(I),'NEITHER','FORMATTED','READWRITE','REWIND','N','N')
               CALL FILE_CLOSE (OT4(I), OT4FIL(I),'DELETE')
            ENDIF
         ENDIF
      ENDDO

! Unformatted files

      OP2FIL(1:I1)  = FILNAM(1:I1)
      OP2FIL(I1+1:) = 'OP2'
      INQUIRE ( FILE=OP2FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( OP2, OP2FIL, OUNT,'OLD    ', OP2_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( OP2, OP2FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( OP2, OP2FIL, OUNT,'REPLACE', OP2_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( OP2, OP2FIL,'DELETE')
         ENDIF
      ENDIF

      F21FIL(1:I1)  = FILNAM(1:I1)
      F21FIL(I1+1:) = 'F21'
      INQUIRE ( FILE=F21FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( F21, F21FIL, OUNT,'OLD    ', F21_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F21, F21FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( F21, F21FIL, OUNT,'REPLACE', F21_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F21, F21FIL,'DELETE')
         ENDIF
      ENDIF

      F22FIL(1:I1)  = FILNAM(1:I1)
      F22FIL(I1+1:) = 'F22'
      INQUIRE ( FILE=F22FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( F22, F22FIL, OUNT,'OLD    ', F22_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F22, F22FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( F22, F22FIL, OUNT,'REPLACE', F22_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F22, F22FIL,'DELETE')
         ENDIF
      ENDIF

      F23FIL(1:I1)  = FILNAM(1:I1)
      F23FIL(I1+1:) = 'F23'
      INQUIRE ( FILE=F23FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( F23, F23FIL, OUNT,'OLD    ', F23_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F23, F23FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( F23, F23FIL, OUNT,'REPLACE', F23_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F23, F23FIL,'DELETE')
         ENDIF
      ENDIF

      F24FIL(1:I1)  = FILNAM(1:I1)
      F24FIL(I1+1:) = 'F24'
      INQUIRE ( FILE=F24FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( F24, F24FIL, OUNT,'OLD    ', F24_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F24, F24FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( F24, F24FIL, OUNT,'REPLACE', F24_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F24, F24FIL,'DELETE')
         ENDIF
      ENDIF

      F25FIL(1:I1)  = FILNAM(1:I1)
      F25FIL(I1+1:) = 'F25'
      INQUIRE ( FILE=F25FIL, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( F25, F25FIL, OUNT,'OLD    ', F25_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F25, F25FIL,'KEEP')
         ELSE
            CALL FILE_OPEN ( F25, F25FIL, OUNT,'REPLACE', F25_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( F25, F25FIL,'DELETE')
         ENDIF
      ENDIF

      LINK1B(1:I1)  = FILNAM(1:I1)
      LINK1B(I1+1:) = 'L1B'
      INQUIRE ( FILE=LINK1B, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1B, LINK1B, OUNT,'OLD    ', L1B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1B, LINK1B,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1B, LINK1B, OUNT,'REPLACE', L1B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1B, LINK1B,'DELETE')
         ENDIF
      ENDIF

      LINK1C(1:I1)  = FILNAM(1:I1)
      LINK1C(I1+1:) = 'L1C'
      INQUIRE ( FILE=LINK1C, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1C, LINK1C, OUNT,'OLD    ', L1C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1C, LINK1C,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1C, LINK1C, OUNT,'REPLACE', L1C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1C, LINK1C,'DELETE')
         ENDIF
      ENDIF

      LINK1D(1:I1)  = FILNAM(1:I1)
      LINK1D(I1+1:) = 'L1D'
      INQUIRE ( FILE=LINK1D, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1D, LINK1D, OUNT,'OLD    ', L1D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1D, LINK1D,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1D, LINK1D, OUNT,'REPLACE', L1D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1D, LINK1D,'DELETE')
         ENDIF
      ENDIF

      LINK1E(1:I1)  = FILNAM(1:I1)
      LINK1E(I1+1:) = 'L1E'
      INQUIRE ( FILE=LINK1E, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1E, LINK1E, OUNT,'OLD    ', L1E_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1E, LINK1E,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1E, LINK1E, OUNT,'REPLACE', L1E_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1E, LINK1E,'DELETE')
         ENDIF
      ENDIF

      LINK1F(1:I1)  = FILNAM(1:I1)
      LINK1F(I1+1:) = 'L1F'
      INQUIRE ( FILE=LINK1F, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1F, LINK1F, OUNT,'OLD    ', L1F_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1F, LINK1F,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1F, LINK1F, OUNT,'REPLACE', L1F_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1F, LINK1F,'DELETE')
         ENDIF
      ENDIF

      LINK1G(1:I1)  = FILNAM(1:I1)
      LINK1G(I1+1:) = 'L1G'
      INQUIRE ( FILE=LINK1G, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1G, LINK1G, OUNT,'OLD    ', L1G_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1G, LINK1G,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1G, LINK1G, OUNT,'REPLACE', L1G_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1G, LINK1G,'DELETE')
         ENDIF
      ENDIF

      LINK1H(1:I1)  = FILNAM(1:I1)
      LINK1H(I1+1:) = 'L1H'
      INQUIRE ( FILE=LINK1H, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1H, LINK1H, OUNT,'OLD    ', L1H_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1H, LINK1H,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1H, LINK1H, OUNT,'REPLACE', L1H_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1H, LINK1H,'DELETE')
         ENDIF
      ENDIF

      LINK1I(1:I1)  = FILNAM(1:I1)
      LINK1I(I1+1:) = 'L1I'
      INQUIRE ( FILE=LINK1I, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1I, LINK1I, OUNT,'OLD    ', L1I_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1I, LINK1I,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1I, LINK1I, OUNT,'REPLACE', L1I_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1I, LINK1I,'DELETE')
         ENDIF
      ENDIF

      LINK1J(1:I1)  = FILNAM(1:I1)
      LINK1J(I1+1:) = 'L1J'
      INQUIRE ( FILE=LINK1J, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1J, LINK1J, OUNT,'OLD    ', L1J_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1J, LINK1J,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1J, LINK1J, OUNT,'REPLACE', L1J_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1J, LINK1J,'DELETE')
         ENDIF
      ENDIF

      LINK1K(1:I1)  = FILNAM(1:I1)
      LINK1K(I1+1:) = 'L1K'
      INQUIRE ( FILE=LINK1K, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1K, LINK1K, OUNT,'OLD    ', L1K_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1K, LINK1K,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1K, LINK1K, OUNT,'REPLACE', L1K_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1K, LINK1K,'DELETE')
         ENDIF
      ENDIF

      LINK1L(1:I1)  = FILNAM(1:I1)
      LINK1L(I1+1:) = 'L1L'
      INQUIRE ( FILE=LINK1L, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1L, LINK1L, OUNT,'OLD    ', L1L_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1L, LINK1L,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1L, LINK1L, OUNT,'REPLACE', L1L_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1L, LINK1L,'DELETE')
         ENDIF
      ENDIF

      LINK1M(1:I1)  = FILNAM(1:I1)
      LINK1M(I1+1:) = 'L1M'
      INQUIRE ( FILE=LINK1M, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1M, LINK1M, OUNT,'OLD    ', L1M_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1M, LINK1M,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1M, LINK1M, OUNT,'REPLACE', L1M_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1M, LINK1M,'DELETE')
         ENDIF
      ENDIF

      LINK1N(1:I1)  = FILNAM(1:I1)
      LINK1N(I1+1:) = 'L1N'
      INQUIRE ( FILE=LINK1N, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1N, LINK1N, OUNT,'OLD    ', L1N_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1N, LINK1N,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1N, LINK1N, OUNT,'REPLACE', L1N_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1N, LINK1N,'DELETE')
         ENDIF
      ENDIF

      LINK1O(1:I1)  = FILNAM(1:I1)
      LINK1O(I1+1:) = 'L1O'
      INQUIRE ( FILE=LINK1O, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1O, LINK1O, OUNT,'OLD    ', L1O_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1O, LINK1O,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1O, LINK1O, OUNT,'REPLACE', L1O_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1O, LINK1O,'DELETE')
         ENDIF
      ENDIF

      LINK1P(1:I1)  = FILNAM(1:I1)
      LINK1P(I1+1:) = 'L1P'
      INQUIRE ( FILE=LINK1P, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1P, LINK1P, OUNT,'OLD    ', L1P_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1P, LINK1P,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1P, LINK1P, OUNT,'REPLACE', L1P_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1P, LINK1P,'DELETE')
         ENDIF
      ENDIF

      LINK1Q(1:I1)  = FILNAM(1:I1)
      LINK1Q(I1+1:) = 'L1Q'
      INQUIRE ( FILE=LINK1Q, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1Q, LINK1Q, OUNT,'OLD    ', L1Q_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Q, LINK1Q,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1Q, LINK1Q, OUNT,'REPLACE', L1Q_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Q, LINK1Q,'DELETE')
         ENDIF
      ENDIF

      LINK1R(1:I1)  = FILNAM(1:I1)
      LINK1R(I1+1:) = 'L1R'
      INQUIRE ( FILE=LINK1R, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1R, LINK1R, OUNT,'OLD    ', L1R_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1R, LINK1R,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1R, LINK1R, OUNT,'REPLACE', L1R_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1R, LINK1R,'DELETE')
         ENDIF
      ENDIF

      LINK1S(1:I1)  = FILNAM(1:I1)
      LINK1S(I1+1:) = 'L1S'
      INQUIRE ( FILE=LINK1S, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1S, LINK1S, OUNT,'OLD    ', L1S_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1S, LINK1S,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1S, LINK1S, OUNT,'REPLACE', L1S_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1S, LINK1S,'DELETE')
         ENDIF
      ENDIF

      LINK1T(1:I1)  = FILNAM(1:I1)
      LINK1T(I1+1:) = 'L1T'
      INQUIRE ( FILE=LINK1T, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1T, LINK1T, OUNT,'OLD    ', L1T_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1T, LINK1T,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1T, LINK1T, OUNT,'REPLACE', L1T_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1T, LINK1T,'DELETE')
         ENDIF
      ENDIF

      LINK1U(1:I1)  = FILNAM(1:I1)
      LINK1U(I1+1:) = 'L1U'
      INQUIRE ( FILE=LINK1U, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1U, LINK1U, OUNT,'OLD    ', L1U_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1U, LINK1U,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1U, LINK1U, OUNT,'REPLACE', L1U_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1U, LINK1U,'DELETE')
         ENDIF
      ENDIF

      LINK1V(1:I1)  = FILNAM(1:I1)
      LINK1V(I1+1:) = 'L1V'
      INQUIRE ( FILE=LINK1V, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1V, LINK1V, OUNT,'OLD    ', L1V_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1V, LINK1V,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1V, LINK1V, OUNT,'REPLACE', L1V_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1V, LINK1V,'DELETE')
         ENDIF
      ENDIF

      LINK1W(1:I1)  = FILNAM(1:I1)
      LINK1W(I1+1:) = 'L1W'
      INQUIRE ( FILE=LINK1W, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1W, LINK1W, OUNT,'OLD    ', L1W_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1W, LINK1W,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1W, LINK1W, OUNT,'REPLACE', L1W_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1W, LINK1W,'DELETE')
         ENDIF
      ENDIF

      LINK1X(1:I1)  = FILNAM(1:I1)
      LINK1X(I1+1:) = 'L1X'
      INQUIRE ( FILE=LINK1X, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1X, LINK1X, OUNT,'OLD    ', L1X_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1X, LINK1X,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1X, LINK1X, OUNT,'REPLACE', L1X_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1X, LINK1X,'DELETE')
         ENDIF
      ENDIF

      LINK1Y(1:I1)  = FILNAM(1:I1)
      LINK1Y(I1+1:) = 'L1Y'
      INQUIRE ( FILE=LINK1Y, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1Y, LINK1Y, OUNT,'OLD    ', L1Y_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Y, LINK1Y,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1Y, LINK1Y, OUNT,'REPLACE', L1Y_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Y, LINK1Y,'DELETE')
         ENDIF
      ENDIF

      LINK1Z(1:I1)  = FILNAM(1:I1)
      LINK1Z(I1+1:) = 'L1Z'
      INQUIRE ( FILE=LINK1Z, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L1Z, LINK1Z, OUNT,'OLD    ', L1Z_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Z, LINK1Z,'KEEP')
         ELSE
            CALL FILE_OPEN ( L1Z, LINK1Z, OUNT,'REPLACE', L1Z_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L1Z, LINK1Z,'DELETE')
         ENDIF
      ENDIF

      LINK2A(1:I1)  = FILNAM(1:I1)
      LINK2A(I1+1:) = 'L2A'
      INQUIRE ( FILE=LINK2A, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2A, LINK2A, OUNT,'OLD    ', L2A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2A, LINK2A,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2A, LINK2A, OUNT,'REPLACE', L2A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2A, LINK2A,'DELETE')
         ENDIF
      ENDIF

      LINK2B(1:I1)  = FILNAM(1:I1)
      LINK2B(I1+1:) = 'L2B'
      INQUIRE ( FILE=LINK2B, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2B, LINK2B, OUNT,'OLD    ', L2B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2B, LINK2B,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2B, LINK2B, OUNT,'REPLACE', L2B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2B, LINK2B,'DELETE')
         ENDIF
      ENDIF

      LINK2C(1:I1)  = FILNAM(1:I1)
      LINK2C(I1+1:) = 'L2C'
      INQUIRE ( FILE=LINK2C, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2C, LINK2C, OUNT,'OLD    ', L2C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2C, LINK2C,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2C, LINK2C, OUNT,'REPLACE', L2C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2C, LINK2C,'DELETE')
         ENDIF
      ENDIF

      LINK2D(1:I1)  = FILNAM(1:I1)
      LINK2D(I1+1:) = 'L2D'
      INQUIRE ( FILE=LINK2D, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2D, LINK2D, OUNT,'OLD    ', L2D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2D, LINK2D,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2D, LINK2D, OUNT,'REPLACE', L2D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2D, LINK2D,'DELETE')
         ENDIF
      ENDIF

      LINK2E(1:I1)  = FILNAM(1:I1)
      LINK2E(I1+1:) = 'L2E'
      INQUIRE ( FILE=LINK2E, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2E, LINK2E, OUNT,'OLD    ', L2E_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2E, LINK2E,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2E, LINK2E, OUNT,'REPLACE', L2E_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2E, LINK2E,'DELETE')
         ENDIF
      ENDIF

      LINK2F(1:I1)  = FILNAM(1:I1)
      LINK2F(I1+1:) = 'L2F'
      INQUIRE ( FILE=LINK2F, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2F, LINK2F, OUNT,'OLD    ', L2F_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2F, LINK2F,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2F, LINK2F, OUNT,'REPLACE', L2F_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2F, LINK2F,'DELETE')
         ENDIF
      ENDIF

      LINK2G(1:I1)  = FILNAM(1:I1)
      LINK2G(I1+1:) = 'L2G'
      INQUIRE ( FILE=LINK2G, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2G, LINK2G, OUNT,'OLD    ', L2G_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2G, LINK2G,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2G, LINK2G, OUNT,'REPLACE', L2G_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2G, LINK2G,'DELETE')
         ENDIF
      ENDIF

      LINK2H(1:I1)  = FILNAM(1:I1)
      LINK2H(I1+1:) = 'L2H'
      INQUIRE ( FILE=LINK2H, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2H, LINK2H, OUNT,'OLD    ', L2H_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2H, LINK2H,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2H, LINK2H, OUNT,'REPLACE', L2H_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2H, LINK2H,'DELETE')
         ENDIF
      ENDIF

      LINK2I(1:I1)  = FILNAM(1:I1)
      LINK2I(I1+1:) = 'L2I'
      INQUIRE ( FILE=LINK2I, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2I, LINK2I, OUNT,'OLD    ', L2I_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2I, LINK2I,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2I, LINK2I, OUNT,'REPLACE', L2I_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2I, LINK2I,'DELETE')
         ENDIF
      ENDIF

      LINK2J(1:I1)  = FILNAM(1:I1)
      LINK2J(I1+1:) = 'L2J'
      INQUIRE ( FILE=LINK2J, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2J, LINK2J, OUNT,'OLD    ', L2J_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2J, LINK2J,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2J, LINK2J, OUNT,'REPLACE', L2J_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2J, LINK2J,'DELETE')
         ENDIF
      ENDIF

      LINK2K(1:I1)  = FILNAM(1:I1)
      LINK2K(I1+1:) = 'L2K'
      INQUIRE ( FILE=LINK2K, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2K, LINK2K, OUNT,'OLD    ', L2K_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2K, LINK2K,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2K, LINK2K, OUNT,'REPLACE', L2K_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2K, LINK2K,'DELETE')
         ENDIF
      ENDIF

      LINK2L(1:I1)  = FILNAM(1:I1)
      LINK2L(I1+1:) = 'L2L'
      INQUIRE ( FILE=LINK2L, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2L, LINK2L, OUNT,'OLD    ', L2L_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2L, LINK2L,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2L, LINK2L, OUNT,'REPLACE', L2L_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2L, LINK2L,'DELETE')
         ENDIF
      ENDIF

      LINK2M(1:I1)  = FILNAM(1:I1)
      LINK2M(I1+1:) = 'L2M'
      INQUIRE ( FILE=LINK2M, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2M, LINK2M, OUNT,'OLD    ', L2M_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2M, LINK2M,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2M, LINK2M, OUNT,'REPLACE', L2M_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2M, LINK2M,'DELETE')
         ENDIF
      ENDIF

      LINK2N(1:I1)  = FILNAM(1:I1)
      LINK2N(I1+1:) = 'L2N'
      INQUIRE ( FILE=LINK2N, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2N, LINK2N, OUNT,'OLD    ', L2N_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2N, LINK2N,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2N, LINK2N, OUNT,'REPLACE', L2N_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2N, LINK2N,'DELETE')
         ENDIF
      ENDIF

      LINK2O(1:I1)  = FILNAM(1:I1)
      LINK2O(I1+1:) = 'L2O'
      INQUIRE ( FILE=LINK2O, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2O, LINK2O, OUNT,'OLD    ', L2O_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2O, LINK2O,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2O, LINK2O, OUNT,'REPLACE', L2O_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2O, LINK2O,'DELETE')
         ENDIF
      ENDIF

      LINK2P(1:I1)  = FILNAM(1:I1)
      LINK2P(I1+1:) = 'L2P'
      INQUIRE ( FILE=LINK2P, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2P, LINK2P, OUNT,'OLD    ', L2P_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2P, LINK2P,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2P, LINK2P, OUNT,'REPLACE', L2P_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2P, LINK2P,'DELETE')
         ENDIF
      ENDIF

      LINK2Q(1:I1)  = FILNAM(1:I1)
      LINK2Q(I1+1:) = 'L2Q'
      INQUIRE ( FILE=LINK2Q, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2Q, LINK2Q, OUNT,'OLD    ', L2Q_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2Q, LINK2Q,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2Q, LINK2Q, OUNT,'REPLACE', L2Q_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2Q, LINK2Q,'DELETE')
         ENDIF
      ENDIF

      LINK2R(1:I1)  = FILNAM(1:I1)
      LINK2R(I1+1:) = 'L2R'
      INQUIRE ( FILE=LINK2R, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2R, LINK2R, OUNT,'OLD    ', L2R_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2R, LINK2R,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2R, LINK2R, OUNT,'REPLACE', L2R_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2R, LINK2R,'DELETE')
         ENDIF
      ENDIF

      LINK2S(1:I1)  = FILNAM(1:I1)
      LINK2S(I1+1:) = 'L2S'
      INQUIRE ( FILE=LINK2S, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2S, LINK2S, OUNT,'OLD    ', L2S_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2S, LINK2S,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2S, LINK2S, OUNT,'REPLACE', L2S_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2S, LINK2S,'DELETE')
         ENDIF
      ENDIF

      LINK2T(1:I1)  = FILNAM(1:I1)
      LINK2T(I1+1:) = 'L2T'
      INQUIRE ( FILE=LINK2T, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L2T, LINK2T, OUNT,'OLD    ', L2T_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2T, LINK2T,'KEEP')
         ELSE
            CALL FILE_OPEN ( L2T, LINK2T, OUNT,'REPLACE', L2T_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L2T, LINK2T,'DELETE')
         ENDIF
      ENDIF

      LINK3A(1:I1)  = FILNAM(1:I1)
      LINK3A(I1+1:) = 'L3A'
      INQUIRE ( FILE=LINK3A, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L3A, LINK3A, OUNT,'OLD    ', L3A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L3A, LINK3A,'KEEP')
         ELSE
            CALL FILE_OPEN ( L3A, LINK3A, OUNT,'REPLACE', L3A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L3A, LINK3A,'DELETE')
         ENDIF
      ENDIF

      LINK4A(1:I1)  = FILNAM(1:I1)
      LINK4A(I1+1:) = 'L4A'
      INQUIRE ( FILE=LINK4A, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L4A, LINK4A, OUNT,'OLD    ', L4A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4A, LINK4A,'KEEP')
         ELSE
            CALL FILE_OPEN ( L4A, LINK4A, OUNT,'REPLACE', L4A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4A, LINK4A,'DELETE')
         ENDIF
      ENDIF

      LINK4B(1:I1)  = FILNAM(1:I1)
      LINK4B(I1+1:) = 'L4B'
      INQUIRE ( FILE=LINK4B, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L4B, LINK4B, OUNT,'OLD    ', L4B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4B, LINK4B,'KEEP')
         ELSE
            CALL FILE_OPEN ( L4B, LINK4B, OUNT,'REPLACE', L4B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4B, LINK4B,'DELETE')
         ENDIF
      ENDIF

      LINK4C(1:I1)  = FILNAM(1:I1)
      LINK4C(I1+1:) = 'L4C'
      INQUIRE ( FILE=LINK4C, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L4C, LINK4C, OUNT,'OLD    ', L4C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4C, LINK4C,'KEEP')
         ELSE
            CALL FILE_OPEN ( L4C, LINK4C, OUNT,'REPLACE', L4C_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4C, LINK4C,'DELETE')
         ENDIF
      ENDIF

      LINK4D(1:I1)  = FILNAM(1:I1)
      LINK4D(I1+1:) = 'L4D'
      INQUIRE ( FILE=LINK4D, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L4D, LINK4D, OUNT,'OLD    ', L4D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4D, LINK4D,'KEEP')
         ELSE
            CALL FILE_OPEN ( L4D, LINK4D, OUNT,'REPLACE', L4D_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L4D, LINK4D,'DELETE')
         ENDIF
      ENDIF

      LINK5A(1:I1)  = FILNAM(1:I1)
      LINK5A(I1+1:) = 'L5A'
      INQUIRE ( FILE=LINK5A, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L5A, LINK5A, OUNT,'OLD    ', L5A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L5A, LINK5A,'KEEP')
         ELSE
            CALL FILE_OPEN ( L5A, LINK5A, OUNT,'REPLACE', L5A_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L5A, LINK5A,'DELETE')
         ENDIF
      ENDIF

      LINK5B(1:I1)  = FILNAM(1:I1)
      LINK5B(I1+1:) = 'L5B'
      INQUIRE ( FILE=LINK5B, EXIST=FILE_EXIST )
      IF (FILE_EXIST) THEN
         IF (RESTART == 'Y') THEN
            CALL FILE_OPEN ( L5B, LINK5B, OUNT,'OLD    ', L5B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L5B, LINK5B,'KEEP')
         ELSE
            CALL FILE_OPEN ( L5B, LINK5B, OUNT,'REPLACE', L5B_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
            CALL FILE_CLOSE ( L5B, LINK5B,'DELETE')
         ENDIF
      ENDIF

      DO I=1,MOU4
         OU4FIL(I)(1:I1)  = FILNAM(1:I1)
         OU4FIL(I)(I1+1:) = OU4_EXT(I)
         INQUIRE ( FILE=OU4FIL(I), EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            IF (RESTART == 'Y') THEN                       ! Keep these until after reading Exec Cont and finding out which ones
               CALL FILE_OPEN ( OU4(I), OU4FIL(I), OUNT,'OLD    ', OU4_MSG(I),'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
               CALL FILE_CLOSE ( OU4(I), OU4FIL(I),'KEEP') ! will be needed. Then close those and 'KEEP' them
            ELSE
               CALL FILE_OPEN ( OU4(I), OU4FIL(I), OUNT,'REPLACE', OU4_MSG(I),'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
               CALL FILE_CLOSE ( OU4(I), OU4FIL(I),'DELETE')
            ENDIF
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
  150 FORMAT(/,' >> MYSTRAN BEGIN  : ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,' The input file is ',A,/)

  170 FORMAT(/,' >> MYSTRAN RESTART: ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,' The input file is ',A,/)

! **********************************************************************************************************************************

      END SUBROUTINE MYSTRAN_FILES


      SUBROUTINE PROCESS_INCLUDE_FILES ( NUM_INCL_FILES )

! PROCESS_INCLUDE_FILES reads in the complete Bulk Data file, checks for INCLUDE files and calls 2 subrs to process them.
! At completion, a new data file is created (IN0FIL on unit IN0) which has the original Bulk Data file plus the contents of all of
! the INCLUDE files. This file then becomes the input file by re-opening it as INFILE (the normal input file)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN0, IN1, INC, INFILE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, READ_INCLUDE_FILNAM, RW_INCLUDE_FILES

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PROCESS_INCLUDE_FILES'
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD              ! Exec Control deck card
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1

      INTEGER(LONG), INTENT(OUT)      :: NUM_INCL_FILES    ! Number of INCLUDE files in the Bulk Data file
      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator.
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1




! **********************************************************************************************************************************
! Initialize

      NUM_INCL_FILES = 0

main: DO

         READ(IN1,101,IOSTAT=IOCHK) CARD
                                                           ! Exit when the complete Bulk Data file has been read
         IF (IOCHK < 0) THEN
            EXIT main
         ENDIF

         IF (IOCHK > 0) THEN                               ! Check if error occurs during read
            WRITE(ERR,1010)
            WRITE(F06,1010)
            WRITE(F06,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE main
         ELSE
            WRITE(IN0,101,IOSTAT=IOCHK) CARD
         ENDIF

         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )

         IERR = 0
         IF (CARD1(1:7) == 'INCLUDE' )  THEN
            NUM_INCL_FILES = NUM_INCL_FILES + 1
            CALL READ_INCLUDE_FILNAM ( CARD1, IERR )
            IF (IERR == 0) THEN
               CALL RW_INCLUDE_FILES ( INC, IN0 )
            ENDIF
         ENDIF

      ENDDO main



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ENTRY FROM BULK DATA FILE WHILE SEARCHING FOR INCLUDE FILES . ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

! **********************************************************************************************************************************

      END SUBROUTINE PROCESS_INCLUDE_FILES


! ##################################################################################################################################

      SUBROUTINE READ_CL ( FILNAM, NC_FILNAM )

! Gets command line string, which is the name of a file (FILNAM), and counts the number (NC_FILNAM) of characters
! (leading blanks ignored) in the name.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(OUT)   :: FILNAM            ! File name on command line

      INTEGER(LONG), INTENT(OUT)      :: NC_FILNAM         ! Length, in chars, of FILNAM (with leading blanks removed)
      INTEGER(LONG)                   :: STATUS            ! Status from GETARG. If /= -1, it is the length of the argument

      INTRINSIC                       :: GET_COMMAND_ARGUMENT

! **********************************************************************************************************************************
! Initialize outputs

      NC_FILNAM = 0
      FILNAM(1:FILE_NAM_MAXLEN) = ' '

! Get command line string which should contain the file name of the input data deck.

      CALL GET_COMMAND_ARGUMENT ( 1, FILNAM, NC_FILNAM, STATUS )
! if retrieval was successful, set STATUS to equal length, as per comments
      IF ( STATUS /= -1 ) THEN
        STATUS = NC_FILNAM
      END IF

      RETURN

      END SUBROUTINE READ_CL


      SUBROUTINE READ_INPUT_FILE_NAME ( INI_EXIST )

! Gets the command line when MYSTRAN is invoked. This should contain the filname of the input file to be run. The default for the
! extension of rhe input file is DAT. If the default extension is being used, it does not have to be specified on the command line.
! If any extension other than the default is used it must be supplied with the filename on the command line. Examples of the command
! line are:

! If the input file is filename.dat, the user can invoke MYSTRAN from the directory in which that input file exists with:

!    MYSTRAN filenam

! If the input file is filename.bdf, then that complete name must be supplied

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, DEFDIR, INFILE,                                        &
                                         LEN_INPUT_FNAME, SC1
      USE SCONTR, ONLY                :  PROG_NAME

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  WRITE_FILNAM

      IMPLICIT NONE

      LOGICAL                         :: LEXIST

      CHARACTER( 1*BYTE), INTENT(IN)  :: INI_EXIST         ! 'Y' if file MYSTRAN.INI exists or 'N' otherwise
      CHARACTER( 1*BYTE)              :: CEXT              ! = 'Y' if there is an extension following a decimal point in FILNAM
      CHARACTER(LEN=LEN(INFILE))      :: FILNAM            ! File name
      CHARACTER(LEN=LEN(INFILE))      :: DUMFIL            ! File name
      CHARACTER( 1*BYTE)              :: POINT             ! = 'Y' if we find a decimal point in INFILE (FILNAM)

      INTEGER(LONG)                   :: LEXT              ! Length (chars) of input file extension
      INTEGER(LONG)                   :: NC_TOT            ! Total length (chars) of input file name including directory
      INTEGER(LONG)                   :: NC_DIR            ! Total length (chars) of input file name's directory
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NC_FILNAM         ! Number of chars in FILNAM read on command line by subr READ_CL

      INTRINSIC                       :: INDEX

! **********************************************************************************************************************************
! If directory for files was input in MYSTRAN.INI, load it into the front end of INFILE. NC_TOT is a count on the
! total number of characters in INFILE (including the dir from MYSTRAN.INI as well as FILNAM from subr READ_CL).

      NC_TOT     = 0
      NC_DIR     = 0
      INFILE(1:) = ' '

      IF (INI_EXIST == 'Y') THEN                           ! Get DEFDIR length. Leading blanks were stripped of in subr READ_INI

         NC_DIR = LEN_TRIM ( DEFDIR )

         IF (NC_DIR > 0) THEN
            IF (NC_DIR+1 > FILE_NAM_MAXLEN) THEN           ! make sure that DEFDIR part of INFILE not too long up to this point
               INFILE(1:FILE_NAM_MAXLEN) = DEFDIR(1:FILE_NAM_MAXLEN)
               WRITE(SC1,1002) FILE_NAM_MAXLEN
               WRITE(SC1,'(A)') TRIM(INFILE)
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            INFILE(1:NC_DIR) = DEFDIR(1:NC_DIR)            ! Set INFILE to DEFDIR
            IF (INFILE(NC_DIR:NC_DIR) /= '\') THEN         ! Add a '\' if INFILE (which is DEFDIR at this point) does not have one
               NC_DIR = NC_DIR + 1
               INFILE(NC_DIR:NC_DIR) = '\'
            ENDIF
         ENDIF
         NC_TOT = NC_DIR                                   ! This is length of DEFDIR including '\' (if needed)

      ENDIF

      NC_FILNAM = 0
      FILNAM(1:) = ' '
      CALL READ_CL ( FILNAM, NC_FILNAM )
      IF (NC_FILNAM == 0) THEN
         WRITE(SC1,1006)
         WRITE(SC1,*)
         READ(*,'(A)') FILNAM
         DO I=FILE_NAM_MAXLEN,1,-1
            IF (FILNAM(I:I) == ' ') THEN
               CYCLE
            ELSE
               NC_FILNAM = I
               EXIT
            ENDIF
         ENDDO
      ENDIF

outer:DO                                                   ! Loop which sets filename, check if it exists and cycles back if not
         IF (NC_TOT + NC_FILNAM > FILE_NAM_MAXLEN) THEN
            WRITE(SC1,1002) FILE_NAM_MAXLEN
            WRITE(SC1,'(A)') TRIM(INFILE)
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         INFILE(NC_TOT+1:) = FILNAM(1:)
         NC_TOT = NC_TOT + NC_FILNAM

         CEXT  = 'N'                                       ! If a file extension was not included, add '.DAT'
         POINT = 'N'
inner_1: DO I=NC_TOT,1,-1
            IF (INFILE(I:I) == ' ') THEN
               CYCLE inner_1
            ELSE IF (INFILE(I:I) == '.') THEN              ! '.' indicates there is an extension if non white space after it
               POINT = 'Y'
               IF (INFILE(I+1:I+1) == ' ') THEN            ! Check if whote space or not after '.'
                  CEXT  = 'N'                              ! All white space after '.', so no file ext. is in INFILE at this point
               ELSE
                  CEXT = 'Y'                               ! Non-white space after '.', so file extension is in INFILE
                  EXIT inner_1
               ENDIF
            ENDIF
         ENDDO inner_1

         IF (CEXT == 'N') THEN                             ! If there was no file name extension, add default extension 'DAT'
            IF ((NC_TOT+4) > FILE_NAM_MAXLEN) THEN
               WRITE(SC1,1002) FILE_NAM_MAXLEN
               CALL WRITE_FILNAM ( INFILE, SC1, 1 )
               CALL OUTA_HERE ( 'Y' )
            ELSE
               IF (POINT == 'Y') THEN
                  INFILE(NC_TOT+1:) = 'DAT'       ! Add extension
                  NC_TOT = NC_TOT + 3
               ELSE IF (POINT == 'N') THEN
                  INFILE(NC_TOT+1:) = '.' // 'DAT'
                  NC_TOT = NC_TOT + 4
               ENDIF
            ENDIF
         ENDIF

         DUMFIL(1:) = ' '
inner_2: DO
            INQUIRE (FILE=INFILE,EXIST=LEXIST)             ! Check whether INFILE exists
            IF (.NOT.LEXIST) THEN                          ! INFILE seems to not exist, but maybe ext. is different so check
               IF (CEXT == 'N') THEN

                  DUMFIL(1:NC_TOT-3) = INFILE(1:NC_TOT-3)

                  DUMFIL(NC_TOT-2:NC_TOT) = 'dat'
                  INQUIRE (FILE=DUMFIL,EXIST=LEXIST)
                  IF (LEXIST) THEN
                     INFILE(1:) = DUMFIL(1:)
                     EXIT outer
                  ENDIF

                  DUMFIL(NC_TOT-2:NC_TOT) = 'BDF'
                  INQUIRE (FILE=DUMFIL,EXIST=LEXIST)
                  IF (LEXIST) THEN
                     INFILE(1:) = DUMFIL(1:)
                     EXIT outer
                  ENDIF

                  DUMFIL(NC_TOT-2:NC_TOT) = 'bdf'
                  INQUIRE (FILE=DUMFIL,EXIST=LEXIST)
                  IF (LEXIST) THEN
                     INFILE(1:) = DUMFIL(1:)
                     EXIT outer
                  ENDIF

               ENDIF
               WRITE(SC1,1004)
               CALL WRITE_FILNAM ( INFILE, SC1, 1 )
               WRITE(SC1,1005)
               READ(*,'(A)') INFILE                        ! Get user input of complete file name FILNAM
               CYCLE inner_2
            ELSE
               EXIT outer
            ENDIF
         ENDDO inner_2

      ENDDO outer

! Count length of extension of INFILE (after '.'). File name must have an extension length of at least 1 character.

      LEXT = 0
      DO I=NC_TOT,1,-1
         IF (INFILE(I:I) == '.') THEN
            EXIT
         ELSE
            LEXT = LEXT + 1
         ENDIF
      ENDDO
      LEN_INPUT_FNAME = NC_TOT - LEXT

      RETURN

! **********************************************************************************************************************************
 1002 FORMAT(' Input data file name was too long. max length is ',I8,' File name input was:')

 1004 FORMAT(' Input data file:')

 1005 FORMAT(' does not exist. Specify complete file name with drive and extension',/,                                             &
             ' or hit Ctrl-Break to start over.',/)

 1006 FORMAT(' Input the input data file name:')

! **********************************************************************************************************************************

      END SUBROUTINE READ_INPUT_FILE_NAME

   END MODULE STARTUP_FILES
