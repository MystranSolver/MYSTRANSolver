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

   MODULE FILE_LIFECYCLE

   USE OP2_GEOMETRY_OUTPUT, ONLY :  END_OP2_TABLES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, SET_FILE_CLOSE_STAT, CLOSE_LIJFILES, CLOSE_OUTFILES, WRITE_FILNAM, FILERR, READERR, STMERR, OPNERR, OUTA_HERE, WRITE_L1A

   CONTAINS

      SUBROUTINE FILE_CLOSE ( UNIT, FILNAM, CLOSE_STAT )

! Closes files and writes message if the close fails

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  SC1

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=*)   , INTENT(IN) :: FILNAM            ! File name

      CHARACTER(LEN=*)   , INTENT(IN) :: CLOSE_STAT        ! Status for close

      INTEGER(LONG), INTENT(IN)       :: UNIT              ! File unit number
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when closing a file




! **********************************************************************************************************************************
      CLOSE ( UNIT,STATUS=CLOSE_STAT,IOSTAT=IOCHK )

      IF (IOCHK /= 0) THEN
         WRITE(SC1,903) IOCHK
         CALL WRITE_FILNAM ( FILNAM, SC1, 15 )
         WRITE(SC1,9232)
         STOP
      ENDIF


      RETURN

! **********************************************************************************************************************************
  903 FORMAT(' *ERROR   903: ERROR ENCOUNTERED WITH IOSTAT = ',I8,' CLOSING FILE:')

 9232 FORMAT('               IT MAY BE OPEN IN ANOTHER PROGRAM. IF SO, CLOSE IT & START AGAIN.')

! **********************************************************************************************************************************

      END SUBROUTINE FILE_CLOSE


      SUBROUTINE FILE_INQUIRE ( MESSAGE )

! Inquires about whether files are opened. Writes results to file F06

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN,  MOT4,    MOU4,    OU4_EXT, OT4_EXT

      USE IOUNT1, ONLY                :  BUG,     EIN,     ENF,     ERR,     F06,     IN0,     IN1,     NEU,                       &
                                         SEQ,     SC1,     SPC,                                                                    &
                                         F21,     F22,     F23,     F24,     F25,                                                  &
                                         L1A,     L1B,     L1C,     L1D,     L1E,     L1F,     L1G,     L1H,     L1I,     L1J,     &
                                         L1K,     L1L,     L1M,     L1N,     L1O,     L1P,     L1Q,     L1R,     L1S,     L1T,     &
                                         L1U,     L1V,     L1W,     L1X,     L1Y,     L1Z,                                         &
                                         L2A,     L2B,     L2C,     L2D,     L2E,     L2F,     L2G,     L2H,     L2I,     L2J,     &
                                         L2K,     L2L,     L2M,     L2N,     L2O,     L2P,     L2Q,     L2R,     L2S,     L2T,     &
                                         L3A,     L4A,     L4B,     L4C,     L4D,     L5A,     L5B,     OP2,     OT4,     OU4,     &
                                         MAX_FIL

      USE IOUNT1, ONLY                :  BUGFIL,  EINFIL,  ENFFIL,  ERRFIL,  F06FIL,  IN0FIL,  INFILE,  NEUFIL,                    &
                                         SEQFIL,  SPCFIL,                                                                          &
                                         F21FIL,  F22FIL,  F23FIL,  F24FIL,  F25FIL,                                               &
                                         LINK1A,  LINK1B,  LINK1C,  LINK1D,  LINK1E,  LINK1F,  LINK1G,  LINK1H,  LINK1I,  LINK1J,  &
                                         LINK1K,  LINK1L,  LINK1M,  LINK1N,  LINK1O,  LINK1P,  LINK1Q,  LINK1R,  LINK1S,  LINK1T,  &
                                         LINK1U,  LINK1V,  LINK1W,  LINK1X,  LINK1Y,  LINK1Z,                                      &
                                         LINK2A,  LINK2B,  LINK2C,  LINK2D,  LINK2E,  LINK2F,  LINK2G,  LINK2H,  LINK2I,  LINK2J,  &
                                         LINK2K,  LINK2L,  LINK2M,  LINK2N,  LINK2O,  LINK2P,  LINK2Q,  LINK2R,  LINK2S,  LINK2T,  &
                                         LINK3A,  LINK4A,  LINK4B,  LINK4C,  LINK4D,  LINK5A,  LINK5B,  OP2FIL,  OT4FIL,  OU4FIL

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      LOGICAL                         :: LEXIST            ! True if file exists
      LOGICAL                         :: LOPND             ! True if file is opened

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FILE_INQUIRE'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAGE           ! Message written when this subr is called
      CHARACTER( 3*BYTE)              :: FIL(100)          ! Descriptor of a MYSTRAN file
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: FILNAM(100)       ! Filename
      CHARACTER(14*BYTE)              :: ANSE              ! Message set based on value of LEXIST
      CHARACTER(10*BYTE)              :: ANSO              ! Message set based on value of LOPND

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: UNT(100)          ! Unit number of a MYSTRAN file




! **********************************************************************************************************************************
      FIL(  1) = 'SC1'   ;   UNT(  1) =  SC1               ! SC1 - don't need to do INQUIRE on it
      FIL(  3) = 'BUG'   ;   UNT(  3) =  BUG   ;   FILNAM(  3) = BUGFIL
      FIL(  4) = 'EIN'   ;   UNT(  4) =  EIN   ;   FILNAM(  4) = EINFIL
      FIL(  5) = 'ENF'   ;   UNT(  5) =  ENF   ;   FILNAM(  5) = ENFFIL
      FIL(  6) = 'ERR'   ;   UNT(  6) =  ERR   ;   FILNAM(  6) = ERRFIL
      FIL(  8) = 'F06'   ;   UNT(  8) =  F06   ;   FILNAM(  8) = F06FIL
      FIL(  9) = 'L1A'   ;   UNT(  9) =  L1A   ;   FILNAM(  9) = LINK1A
      FIL( 10) = 'IN0'   ;   UNT( 10) =  IN0   ;   FILNAM( 10) = IN0FIL
      FIL( 11) = 'NEU'   ;   UNT( 11) =  NEU   ;   FILNAM( 11) = NEUFIL
      FIL( 13) = 'SEQ'   ;   UNT( 13) =  SEQ   ;   FILNAM( 13) = SEQFIL
      FIL( 14) = 'SPC'   ;   UNT( 14) =  SPC   ;   FILNAM( 14) = SPCFIL
      FIL( 15) = 'F21'   ;   UNT( 15) =  F21   ;   FILNAM( 15) = F21FIL
      FIL( 16) = 'F22'   ;   UNT( 16) =  F22   ;   FILNAM( 16) = F22FIL
      FIL( 17) = 'F23'   ;   UNT( 17) =  F23   ;   FILNAM( 17) = F23FIL
      FIL( 18) = 'F24'   ;   UNT( 18) =  F24   ;   FILNAM( 18) = F24FIL
      FIL( 19) = 'F25'   ;   UNT( 19) =  F25   ;   FILNAM( 19) = F25FIL
      FIL( 20) = 'L1B'   ;   UNT( 20) =  L1B   ;   FILNAM( 20) = LINK1B
      FIL( 21) = 'L1C'   ;   UNT( 21) =  L1C   ;   FILNAM( 21) = LINK1C
      FIL( 22) = 'L1D'   ;   UNT( 22) =  L1D   ;   FILNAM( 22) = LINK1D
      FIL( 23) = 'L1E'   ;   UNT( 23) =  L1E   ;   FILNAM( 23) = LINK1E
      FIL( 24) = 'L1F'   ;   UNT( 24) =  L1F   ;   FILNAM( 24) = LINK1F
      FIL( 25) = 'L1G'   ;   UNT( 25) =  L1G   ;   FILNAM( 25) = LINK1G
      FIL( 26) = 'L1H'   ;   UNT( 26) =  L1H   ;   FILNAM( 26) = LINK1H
      FIL( 27) = 'L1I'   ;   UNT( 27) =  L1I   ;   FILNAM( 27) = LINK1I
      FIL( 28) = 'L1J'   ;   UNT( 28) =  L1J   ;   FILNAM( 28) = LINK1J
      FIL( 29) = 'L1K'   ;   UNT( 29) =  L1K   ;   FILNAM( 29) = LINK1K
      FIL( 30) = 'L1L'   ;   UNT( 30) =  L1L   ;   FILNAM( 30) = LINK1L
      FIL( 31) = 'L1M'   ;   UNT( 31) =  L1M   ;   FILNAM( 31) = LINK1M
      FIL( 32) = 'L1N'   ;   UNT( 32) =  L1N   ;   FILNAM( 32) = LINK1N
      FIL( 33) = 'L1O'   ;   UNT( 33) =  L1O   ;   FILNAM( 33) = LINK1O
      FIL( 34) = 'L1P'   ;   UNT( 34) =  L1P   ;   FILNAM( 34) = LINK1P
      FIL( 35) = 'L1Q'   ;   UNT( 35) =  L1Q   ;   FILNAM( 35) = LINK1Q
      FIL( 36) = 'L1R'   ;   UNT( 36) =  L1R   ;   FILNAM( 36) = LINK1R
      FIL( 37) = 'L1S'   ;   UNT( 37) =  L1S   ;   FILNAM( 37) = LINK1S
      FIL( 38) = 'L1T'   ;   UNT( 38) =  L1T   ;   FILNAM( 38) = LINK1T
      FIL( 39) = 'L1U'   ;   UNT( 39) =  L1U   ;   FILNAM( 39) = LINK1U
      FIL( 40) = 'L1V'   ;   UNT( 40) =  L1V   ;   FILNAM( 40) = LINK1V
      FIL( 41) = 'L1W'   ;   UNT( 41) =  L1W   ;   FILNAM( 41) = LINK1W
      FIL( 42) = 'L1X'   ;   UNT( 42) =  L1X   ;   FILNAM( 42) = LINK1X
      FIL( 43) = 'L1Y'   ;   UNT( 43) =  L1Y   ;   FILNAM( 43) = LINK1Y
      FIL( 44) = 'L1Z'   ;   UNT( 44) =  L1Z   ;   FILNAM( 44) = LINK1Z
      FIL( 45) = 'L2A'   ;   UNT( 45) =  L2A   ;   FILNAM( 45) = LINK2A
      FIL( 46) = 'L2B'   ;   UNT( 46) =  L2B   ;   FILNAM( 46) = LINK2B
      FIL( 47) = 'L2C'   ;   UNT( 47) =  L2C   ;   FILNAM( 47) = LINK2C
      FIL( 48) = 'L2D'   ;   UNT( 48) =  L2D   ;   FILNAM( 48) = LINK2D
      FIL( 49) = 'L2E'   ;   UNT( 49) =  L2E   ;   FILNAM( 49) = LINK2E
      FIL( 50) = 'L2F'   ;   UNT( 50) =  L2F   ;   FILNAM( 50) = LINK2F
      FIL( 51) = 'L2G'   ;   UNT( 51) =  L2G   ;   FILNAM( 51) = LINK2G
      FIL( 52) = 'L2H'   ;   UNT( 52) =  L2H   ;   FILNAM( 52) = LINK2H
      FIL( 53) = 'L2I'   ;   UNT( 53) =  L2I   ;   FILNAM( 53) = LINK2I
      FIL( 54) = 'L2J'   ;   UNT( 54) =  L2J   ;   FILNAM( 54) = LINK2J
      FIL( 55) = 'L2K'   ;   UNT( 55) =  L2K   ;   FILNAM( 55) = LINK2K
      FIL( 56) = 'L2L'   ;   UNT( 56) =  L2L   ;   FILNAM( 56) = LINK2L
      FIL( 57) = 'L2M'   ;   UNT( 57) =  L2M   ;   FILNAM( 57) = LINK2M
      FIL( 58) = 'L2N'   ;   UNT( 58) =  L2N   ;   FILNAM( 58) = LINK2N
      FIL( 59) = 'L2O'   ;   UNT( 59) =  L2O   ;   FILNAM( 59) = LINK2O
      FIL( 60) = 'L2P'   ;   UNT( 60) =  L2P   ;   FILNAM( 60) = LINK2P
      FIL( 61) = 'L2Q'   ;   UNT( 61) =  L2Q   ;   FILNAM( 61) = LINK2Q
      FIL( 62) = 'L2R'   ;   UNT( 62) =  L2R   ;   FILNAM( 62) = LINK2R
      FIL( 63) = 'L2S'   ;   UNT( 63) =  L2S   ;   FILNAM( 63) = LINK2S
      FIL( 64) = 'L2T'   ;   UNT( 64) =  L2T   ;   FILNAM( 64) = LINK2T
      FIL( 65) = 'L3A'   ;   UNT( 65) =  L3A   ;   FILNAM( 65) = LINK3A
      FIL( 66) = 'L4A'   ;   UNT( 66) =  L4A   ;   FILNAM( 66) = LINK4A
      FIL( 67) = 'L4B'   ;   UNT( 67) =  L4B   ;   FILNAM( 67) = LINK4B
      FIL( 68) = 'L4C'   ;   UNT( 68) =  L4C   ;   FILNAM( 68) = LINK4C
      FIL( 69) = 'L4D'   ;   UNT( 69) =  L4D   ;   FILNAM( 69) = LINK4D
      FIL( 70) = 'L5A'   ;   UNT( 70) =  L5A   ;   FILNAM( 70) = LINK5A
      FIL( 71) = 'L5B'   ;   UNT( 71) =  L5B   ;   FILNAM( 71) = LINK5B
      FIL( 72) = 'OP2'   ;   UNT( 72) =  OP2   ;   FILNAM( 71) = OP2FIL

      IF ( 72 > MAX_FIL) THEN
         WRITE(ERR,944) SUBR_NAME, MAX_FIL
         WRITE(F06,944) SUBR_NAME, MAX_FIL
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO I=1,MOU4
         FIL(MAX_FIL+I)      = OU4_EXT(I);   UNT(MAX_FIL+I)      =  OU4(I)   ;   FILNAM(MAX_FIL+I)      = OU4FIL(I)
      ENDDO
      DO I=1,MOT4
         FIL(MAX_FIL+MOU4+I) = OT4_EXT(I);   UNT(MAX_FIL+MOU4+I) =  OT4(I)   ;   FILNAM(MAX_FIL+MOU4+I) = OT4FIL(I)
      ENDDO

      WRITE(ERR,1)
      WRITE(ERR,2) MESSAGE

      WRITE(F06,1)
      WRITE(F06,2) MESSAGE

      DO I=2,MAX_FIL+MOU4+MOT4                                  ! Start at 2 since SC1 was 1 and we do not do INQUIRE on it
         IF (UNT(I) /= 0) THEN
           INQUIRE(FILE=FILNAM(I),EXIST=LEXIST,OPENED=LOPND)
           IF (LEXIST) THEN
               ANSE = 'exists'
            ELSE
               ANSE = 'does not exist'
            ENDIF
            IF (LOPND) THEN
               ANSO = '    opened'
            ELSE
               ANSO = 'not opened'
            ENDIF
            WRITE(ERR,2999) FIL(I), UNT(I), ANSE, ANSO, FILNAM(I)
            WRITE(F06,2999) FIL(I), UNT(I), ANSE, ANSO, FILNAM(I)
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
    1 FORMAT(' *******************************************************************************************************************')

    2 FORMAT(' Messages from subr FILE_INQUIRE on whether files are opened ',A)

  944 FORMAT(' *ERROR   944: PROGRAMMING ERROR IN SUBROUTINE ',A,/,                                                                &
                        14X,'ATTEMPT TO EXCEED MAX_FIL = ',I4,' NUMBER OF FILES IN WRITING FILE STATUS')

 2999 FORMAT(1X,A3,' on unit ',I4,1X,A14,' and is ',A10,' with file name = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE FILE_INQUIRE








































































      SUBROUTINE FILE_OPEN (UNIT, FILNAM, OUNT, STATUS, MESSAG, RW_STIME, FORMAT, ACTION, POSITION, WRITE_L1A, WRITE_VER)

      ! Opens formatted files that have STIME for read or write. If open for read, check STIME. If open for write, write STIME
      ! If file needs to be opened for READWRITE, this subr needs to be called twice:
      !     (1) called 1st time for READ (to open file and to read and check STIME)
      !     (2) called 2nd time (after closing by calling subr) to position at end for subsequent writing after returning
      !
      ! Parameters
      ! ----------
      ! UNIT : int
      !    the file number to open
      ! FILNAM : str
      !    the filename to open
      ! OUNT : tuple[int, int]
      !    I think this is Bill's black magic for handling fortran files...
      ! STATUS : str8
      !    not checked....
      !    'NEW': creates a new file
      !    'OLD': opens an existing file
      !    'REPLACE' : If the file doesn't exist, it creates it. If it does exist, then it replaces it.
      !    'UNKNOWN' : ???
      ! MESSAG : str
      !    used for error messages
      ! RW_STIME : ???
      !    'READ_STIME'  : read the stime?
      !    'WRITE_STIME' : write the stime?
      !    'NEITHER'     : skip this
      ! FORMAT : str
      !    'FORMATTED'   : for ASCII files
      !    'UNFORMATTED' : for binary files that have records
      !    'BINARY'      : binary files, but you get to handle everything; no ACCESS or RECL
      !    'SYSTEM'      : ???; no ACCESS or RECL
      !    Opening a file with FORM='BINARY' has roughly the same effect as FORM='UNFORMATTED', except
      !    that no record lengths are embedded in the file. Without this data, there is no way to tell
      !    where one record begins, or ends. Thus, it is impossible to BACKSPACE a FORM='BINARY' file,
      !    because there is no way of telling where to backspace to. A READ on a 'BINARY' file will read
      !    as much data as needed to fill the variables on the input list.
      ! ACTION : str8
      !    'READ' : open the file for reading
      !    'WRITE': open the file for writing
      !    'READWRITE': ???
      ! POSITION : str8
      !    'ASIS'   : ???
      !    'APPEND' : open the file in append mode
      !    'REWIND' : start from the beginning of the file
      ! RW_STIME : ???
      !    ???
      ! WRITE_L1A : str1
      !    Y/N : write to L1A???
      ! WRITE_VER : str1
      !    Y/N : write to VER???
      !
      ! Unused
      ! -------
      ! ACCESS : str
      !    'SEQUENTIAL' : ???; FILE_OPEN always uses this
      !    'DIRECT'     : ???
      !    'STREAM'     : ???
      ! DIRECT="NO"
      ! RECL : int / None
      !    this is not modifyable per block; record_length
      !    int : the size of the record length; if you go over the size, it'll 0-pad
      !    None: automatic
      !
      ! Examples
      ! --------
      ! FILE_OPEN ( OP2, OP2FIL, OUNT,'OLD    ', OP2_MSG,'NEITHER','UNFORMATTED','WRITE','REWIND','N','N')
      !
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06, IN1, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, PROG_NAME
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE DEBUG_PARAMETERS
      USE MYSTRAN_Version, ONLY       :  MYSTRAN_VER_NUM, MYSTRAN_VER_MONTH, MYSTRAN_VER_DAY, MYSTRAN_VER_YEAR, MYSTRAN_AUTHOR,  &
                                         MYSTRAN_COMMENT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FILE_OPEN'
      CHARACTER(LEN=*), INTENT(IN)    :: ACTION            ! File description
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name
      CHARACTER(LEN=*), INTENT(IN)    :: FORMAT            ! File format
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description
      CHARACTER(LEN=*), INTENT(IN)    :: POSITION          ! File error message
      CHARACTER(LEN=*), INTENT(IN)    :: STATUS            ! File status indicator (NEW, OLD, REPLACE)
      CHARACTER(LEN=*), INTENT(IN)    :: RW_STIME          ! Indicator of whether to read or write STIME
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_L1A         ! 'Y'/'N' Arg passed to subr OUTA_HERE
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_VER         ! 'Y'/'N' Arg to tell whether to write MYSTRAN version info
      CHARACTER( 9*BYTE)              :: NAM_ACT
      CHARACTER(11*BYTE)              :: NAM_FOR
      CHARACTER( 6*BYTE)              :: NAM_POS
      CHARACTER( 7*BYTE)              :: NAM_STA
      CHARACTER(11*BYTE)              :: NAM_RWS

      INTEGER(LONG), INTENT(IN)       :: UNIT              ! Unit number file is attached to
      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to. Input to subr FILE_OPEN
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Error count
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from file




! **********************************************************************************************************************************
      ! Check inputs for sensibility (coding errors if wrong)
      IERR = 0

      ! Check FORM
      IF ((FORMAT /= 'FORMATTED') .AND. (FORMAT /= 'UNFORMATTED') .AND. (FORMAT /= 'BINARY'))THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'FORMAT', 'FORMATTED, UNFORMATTED or BINARY', FORMAT
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'FORMAT', 'FORMATTED, UNFORMATTED or BINARY', FORMAT
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      ! Check ACTION
      IF ((ACTION /= 'READ') .AND. (ACTION /= 'WRITE') .AND. (ACTION /= 'READWRITE')) THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'ACTION', 'READ or WRITE', ACTION
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'ACTION', 'READ or WRITE', ACTION
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      ! Check POSITION
      IF ((POSITION /= 'ASIS') .AND. (POSITION /= 'APPEND') .AND. (POSITION /= 'REWIND')) THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'POSITION', 'ASIS or APPEND or REWIND', POSITION
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'POSITION', 'ASIS or APPEND or REWIND', POSITION
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      ! Check STATUS
!     IF ((STATUS /= 'NEW') .AND.(STATUS /= 'OLD') .AND.(STATUS /= 'REPLACE')) THEN
!        IERR = IERR + 1
!        DO I=1,2
!           IF (OUNT(I) > 0) THEN
!              WRITE(OUNT(I),909) SUBR_NAME, 'STATUS', 'NEW or OLD or REPLACE', STATUS
!           ELSE
!              WRITE(SC1,909) SUBR_NAME, 'STATUS', 'NEW or OLD or REPLACE', STATUS
!           ENDIF
!           IF (OUNT(2) == OUNT(1)) EXIT
!        ENDDO
!     ENDIF

      ! Check WRITE_L1A
      IF ((WRITE_L1A /= 'Y') .AND. (WRITE_L1A /= 'N'))THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'WRITE_L1A', 'Y or N', WRITE_L1A
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'WRITE_L1A', 'Y or N', WRITE_L1A
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      ! Check WRITE_STIME
      IF ((RW_STIME /= 'READ_STIME') .AND. (RW_STIME /= 'WRITE_STIME') .AND. (RW_STIME /= 'NEITHER'))THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'RW_STIME', 'READ_STIME, WRITE_STIME or NEITHER', RW_STIME
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'RW_STIME', 'READ, WRITE or NEITHER', RW_STIME
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      IF ((WRITE_VER /= 'Y') .AND. (WRITE_VER /= 'N'))THEN
         IERR = IERR + 1
         DO I=1,2
            IF (OUNT(I) > 0) THEN
               WRITE(OUNT(I),909) SUBR_NAME, 'WRITE_VER', 'Y or N', WRITE_VER
            ELSE
               WRITE(SC1,909) SUBR_NAME, 'WRITE_VER', 'Y or N', WRITE_VER
            ENDIF
            IF (OUNT(2) == OUNT(1)) EXIT
         ENDDO
      ENDIF

      IF (IERR > 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         IF (WRITE_L1A == 'Y') THEN
            CALL OUTA_HERE ( 'Y' )
         ELSE
            CALL OUTA_HERE ( 'N' )
         ENDIF
      ENDIF

      ! No problems with inputs, so proceed to open file
      IERR = 0
      OPEN ( UNIT, FILE=FILNAM, STATUS=STATUS, FORM=FORMAT, ACCESS='SEQUENTIAL', ACTION=ACTION, POSITION=POSITION, IOSTAT=IOCHK )

      IF (DEBUG(52) > 0) THEN
         CALL DEB_FILE_OPEN
      ENDIF

      IF (IOCHK == 0) THEN
         IF      (RW_STIME == 'READ_STIME') THEN           ! Read and check STIME
            IF (FORMAT == 'FORMATTED') THEN
               READ(UNIT,'(1X,I11)',IOSTAT=IOCHK) XTIME
            ELSE
               READ(UNIT,IOSTAT=IOCHK) XTIME
            ENDIF
            IF (IOCHK /= 0) THEN
               REC_NO = 1
               CALL READERR ( IOCHK, FILNAM, MESSAG, REC_NO, OUNT )
               IERR = IERR + 1
            ELSE
               IF (XTIME /= STIME) THEN
                  CALL STMERR ( XTIME, FILNAM, OUNT )
                  IERR = IERR +1
               ENDIF
            ENDIF
         ELSE IF (RW_STIME == 'WRITE_STIME') THEN          ! Write STIME
            IF (FORMAT == 'FORMATTED') THEN
               WRITE(UNIT,'(1X,I11)') STIME
            ELSE
               WRITE(UNIT) STIME
            ENDIF
            IF (WRITE_VER == 'Y') THEN
               WRITE(UNIT,117) PROG_NAME, MYSTRAN_VER_NUM, MYSTRAN_VER_MONTH, MYSTRAN_VER_DAY, MYSTRAN_VER_YEAR,                &
                               MYSTRAN_AUTHOR, MYSTRAN_COMMENT
            ENDIF
         ENDIF
      ELSE
         CALL OPNERR ( IOCHK, FILNAM, OUNT )
         IERR = IERR +1
      ENDIF

      IF (IERR > 0) THEN
         CALL FILERR ( OUNT )
         CALL OUTA_HERE ( 'Y' )
      ENDIF


      RETURN

! **********************************************************************************************************************************
  909 FORMAT(' *ERROR   909: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ILLEGAL INPUT FOR ARGUMENT ',A,'. MUST BE ',A,' BUT IS "',A,'". ERROR OCCURRED OPENING FILE:')

  117 FORMAT(/,1X,A,' Version',5(1X,A),/,A)

  118 FORMAT(/,1X,A,' Version',5(1X,A),/)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEB_FILE_OPEN

      IMPLICIT NONE

! **********************************************************************************************************************************
      WRITE(F06,98720)

      WRITE(F06,'(A,I8)') ' In file OPEN, UNIT     = ', UNIT
      WRITE(F06,'(A,A) ') ' In file OPEN, FILE     = ', FILNAM(1:40)
      WRITE(F06,'(A,A) ') ' In file OPEN, STATUS   = ', STATUS
      WRITE(F06,'(A,A) ') ' In file OPEN, FORM     = ', FORMAT
      WRITE(F06,'(A,A) ') ' In file OPEN, ACTION   = ', ACTION
      WRITE(F06,'(A,A) ') ' In file OPEN, POSITION = ', POSITION
      WRITE(F06,'(A,I8)') ' In file OPEN, IOCHK    = ', IOCHK
      WRITE(F06,*)
      WRITE(F06,*) ' Open ', filnam(1:32), 'RW_STIME, STIME = "', rw_stime, '"   ', stime

      WRITE(F06,*)

      WRITE(F06,98799)

      WRITE(F06,*)

! **********************************************************************************************************************************
98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::::START DEBUG(52) OUTPUT FROM SUBROUTINE FILE_OPEN:::::::::::::::::::::::::',&
             ':::::::::::::::::',/)

98799 FORMAT(' ::::::::::::::::::::::::::::::::::::::::::END DEBUG(52) OUTPUT FROM SUBROUTINE FILE_OPEN::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEB_FILE_OPEN


      END SUBROUTINE FILE_OPEN


      SUBROUTINE SET_FILE_CLOSE_STAT ( CLOSE_STAT )

! This subr gets called by subrs LINK0 and READ_INI to set a common close status for all files.
! LINK0 does this if the run is check pointed
! READ_INI does this if the MYSTRAN.INI file sets a common close status (like KEEP) for ALLFILES

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      USE IOUNT1, ONLY                :  WRT_ERR, ERRSTAT, SEQSTAT, SPCSTAT, L1ASTAT,                                              &
                                         L1BSTAT, L1CSTAT, L1DSTAT, L1ESTAT, L1FSTAT, L1GSTAT, L1HSTAT, L1ISTAT, L1TSTAT, L1JSTAT, &
                                         L1KSTAT, L1LSTAT, L1MSTAT, L1NSTAT, L1OSTAT, L1PSTAT, L1QSTAT, L1RSTAT, L1SSTAT, L1USTAT, &
                                         L1VSTAT, L1WSTAT, L1XSTAT, L1YSTAT, L1ZSTAT,                                              &
                                         L2ASTAT, L2BSTAT, L2CSTAT, L2DSTAT, L2ESTAT, L2FSTAT, L2GSTAT, L2HSTAT, L2ISTAT, L2JSTAT, &
                                         L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, L2OSTAT, L2PSTAT, L2QSTAT, L2RSTAT, L2SSTAT, L2TSTAT, &
                                         L3ASTAT, L4ASTAT, L4BSTAT, L4CSTAT, L4DSTAT, L5ASTAT, L5BSTAT
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! Close status of files

! **********************************************************************************************************************************

      L1ASTAT = CLOSE_STAT
      L1BSTAT = CLOSE_STAT
      L1CSTAT = CLOSE_STAT
      L1DSTAT = CLOSE_STAT
      L1ESTAT = CLOSE_STAT
      L1FSTAT = CLOSE_STAT
      L1GSTAT = CLOSE_STAT
      L1HSTAT = CLOSE_STAT
      L1ISTAT = CLOSE_STAT
      L1JSTAT = CLOSE_STAT
      L1KSTAT = CLOSE_STAT
      L1LSTAT = CLOSE_STAT
      L1MSTAT = CLOSE_STAT
      L1NSTAT = CLOSE_STAT
      L1OSTAT = CLOSE_STAT
      L1PSTAT = CLOSE_STAT
      L1QSTAT = CLOSE_STAT
      L1RSTAT = CLOSE_STAT
      L1SSTAT = CLOSE_STAT
      L1TSTAT = CLOSE_STAT
      L1USTAT = CLOSE_STAT
      L1VSTAT = CLOSE_STAT
      L1WSTAT = CLOSE_STAT
      L1XSTAT = CLOSE_STAT
      L1YSTAT = CLOSE_STAT
      L1ZSTAT = CLOSE_STAT
      L2ASTAT = CLOSE_STAT
      L2BSTAT = CLOSE_STAT
      L2CSTAT = CLOSE_STAT
      L2DSTAT = CLOSE_STAT
      L2ESTAT = CLOSE_STAT
      L2FSTAT = CLOSE_STAT
      L2GSTAT = CLOSE_STAT
      L2HSTAT = CLOSE_STAT
      L2ISTAT = CLOSE_STAT
      L2JSTAT = CLOSE_STAT
      L2KSTAT = CLOSE_STAT
      L2LSTAT = CLOSE_STAT
      L2MSTAT = CLOSE_STAT
      L2NSTAT = CLOSE_STAT
      L2OSTAT = CLOSE_STAT
      L2PSTAT = CLOSE_STAT
      L2QSTAT = CLOSE_STAT
      L2RSTAT = CLOSE_STAT
      L2SSTAT = CLOSE_STAT
      L2TSTAT = CLOSE_STAT
      L3ASTAT = CLOSE_STAT
      L3ASTAT = CLOSE_STAT
      L4ASTAT = CLOSE_STAT
      L4BSTAT = CLOSE_STAT
      L4CSTAT = CLOSE_STAT
      L4DSTAT = CLOSE_STAT
      L5ASTAT = CLOSE_STAT
      L5BSTAT = CLOSE_STAT

! **********************************************************************************************************************************

      END SUBROUTINE SET_FILE_CLOSE_STAT


      SUBROUTINE CLOSE_LIJFILES ( CLOSE_STAT )

! Closes Lij unformatted files

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, RESTART
      USE IOUNT1, ONLY                :  MOU4,    WRT_ERR, ERR, F06,                                                      &
                                         L1B,     L1C,     L1D,     L1E,     L1F,     L1G,     L1H,     L1I,     L1J,     L1K,     &
                                         L1L,     L1M,     L1N,     L1O,     L1P,     L1Q,     L1R,     L1S,     L1T,     L1U,     &
                                         L1V,     L1W,     L1X,     L1Y,     L1Z,                                                  &
                                         L2A,     L2B,     L2C,     L2D,     L2E,     L2F,     L2G,     L2H,     L2I,     L2J,     &
                                         L2K,     L2L,     L2M,     L2N,     L2O,     L2P,     L2Q,     L2R,     L2S,     L2T,     &
                                         L3A,     L4A,     L4B,     L4C,     L4D,     L5A,     L5B,     OU4

      USE IOUNT1, ONLY                :  WRT_ERR,                                                                &
                                         LINK1B,  LINK1C,  LINK1D,  LINK1E,  LINK1F,  LINK1G,  LINK1H,  LINK1I,  LINK1J,  LINK1K,  &
                                         LINK1L,  LINK1M,  LINK1N,  LINK1O,  LINK1P,  LINK1Q,  LINK1R,  LINK1S,  LINK1T,  LINK1U,  &
                                         LINK1V,  LINK1W,  LINK1X,  LINK1Y,  LINK1Z,                                               &
                                         LINK2A,  LINK2B,  LINK2C,  LINK2D,  LINK2E,  LINK2F,  LINK2G,  LINK2H,  LINK2I,  LINK2J,  &
                                         LINK2K,  LINK2L,  LINK2M,  LINK2N,  LINK2O,  LINK2P,  LINK2Q,  LINK2R,  LINK2S,  LINK2T,  &
                                         LINK3A,  LINK4A,  LINK4B,  LINK4C,  LINK4D,  LINK5A,  LINK5B,  OU4FIL

      USE IOUNT1, ONLY                :  WRT_ERR,                                                                &
                                         L1BSTAT, L1CSTAT, L1DSTAT, L1ESTAT, L1FSTAT, L1GSTAT, L1HSTAT, L1ISTAT, L1JSTAT, L1KSTAT, &
                                         L1LSTAT, L1MSTAT, L1NSTAT, L1OSTAT, L1PSTAT, L1QSTAT, L1RSTAT, L1SSTAT, L1TSTAT, L1USTAT, &
                                         L1VSTAT, L1WSTAT, L1XSTAT, L1YSTAT, L1ZSTAT,                                              &
                                         L2ASTAT, L2BSTAT, L2CSTAT, L2DSTAT, L2ESTAT, L2FSTAT, L2GSTAT, L2HSTAT, L2ISTAT, L2JSTAT, &
                                         L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, L2OSTAT, L2PSTAT, L2QSTAT, L2RSTAT, L2SSTAT, L2TSTAT, &
                                         L3ASTAT, L4ASTAT, L4BSTAT, L4CSTAT, L4DSTAT, L5ASTAT, L5BSTAT, OU4STAT

      USE TIMDAT, ONLY                :  STIME, TSEC

      IMPLICIT NONE

      LOGICAL                         :: LEXIST            ! Result from INQUIRE regarding whether a file exists
      LOGICAL                         :: LOPND             ! Result from INQUIRE regarding whether a file is open

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CLOSE_LIJFILES'
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! Indicator of what to do with file when it is closed

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN

! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

      IF ((CLOSE_STAT /= 'DELETE') .AND. (CLOSE_STAT /= 'FILE_STAT')) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,954) SUBR_NAME, CLOSE_STAT
         WRITE(F06,954) SUBR_NAME, CLOSE_STAT
         STOP
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1B, LINK1B, L1BSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1B, LINK1B, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1C, LINK1C, L1CSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1C, LINK1C, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1D, LINK1D, L1DSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1D, LINK1D, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1E, LINK1E, L1ESTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1E, LINK1E, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1F, LINK1F, L1FSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1F, LINK1F, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1G, LINK1G, L1GSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1G, LINK1G, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1H, LINK1H, L1HSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1H, LINK1H, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1I, LINK1I, L1ISTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1I, LINK1I, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1J, LINK1J, L1JSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1J, LINK1J, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1K, LINK1K, L1KSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1K, LINK1K, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1L, LINK1L, L1LSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1L, LINK1L, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1M, LINK1M, L1MSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1M, LINK1M, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1N, LINK1N, L1NSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1N, LINK1N, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1O, LINK1O, L1OSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1O, LINK1O, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1P, LINK1P, L1PSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1P, LINK1P, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1Q, LINK1Q, L1QSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1Q, LINK1Q, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1R, LINK1R, L1RSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1R, LINK1R, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1S, LINK1S, L1SSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1S, LINK1S, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1T, LINK1T, L1TSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1T, LINK1T, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1U, LINK1U, L1USTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1U, LINK1U, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1V, LINK1V, L1VSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1V, LINK1V, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1W, LINK1W, L1WSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1W, LINK1W, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1X, LINK1X, L1XSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1X, LINK1X, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L1Y, LINK1Y, L1YSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L1Y, LINK1Y, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         IF (RESTART == 'N') THEN
            CALL CLOSE_THIS_FILE ( L1Z, LINK1Z, L1ZSTAT )
         ELSE
            CALL CLOSE_THIS_FILE ( L1Z, LINK1Z, 'KEEP' )
         ENDIF
      ELSE
         IF (RESTART == 'N') THEN
            CALL CLOSE_THIS_FILE ( L1Z, LINK1Z, CLOSE_STAT )
         ELSE
            CALL CLOSE_THIS_FILE ( L1Z, LINK1Z, 'KEEP' )
         ENDIF
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2A, LINK2A, L2ASTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2A, LINK2A, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2B, LINK2B, L2BSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2B, LINK2B, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2C, LINK2C, L2CSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2C, LINK2C, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2D, LINK2D, L2DSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2D, LINK2D, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2E, LINK2E, L2ESTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2E, LINK2E, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2F, LINK2F, L2FSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2F, LINK2F, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2G, LINK2G, L2GSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2G, LINK2G, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2H, LINK2H, L2HSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2H, LINK2H, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2I, LINK2I, L2ISTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2I, LINK2I, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2J, LINK2J, L2JSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2J, LINK2J, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2K, LINK2K, L2KSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2K, LINK2K, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2L, LINK2L, L2LSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2L, LINK2L, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2M, LINK2M, L2MSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2M, LINK2M, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2N, LINK2N, L2NSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2N, LINK2N, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2O, LINK2O, L2OSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2O, LINK2O, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2P, LINK2P, L2PSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2P, LINK2P, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2Q, LINK2Q, L2QSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2Q, LINK2Q, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2R, LINK2R, L2RSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2R, LINK2R, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2S, LINK2S, L2SSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2S, LINK2S, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L2T, LINK2T, L2TSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L2T, LINK2T, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L5A, LINK5A, L5ASTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L5A, LINK5A, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L3A, LINK3A, L3ASTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L3A, LINK3A, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L4A, LINK4A, L4ASTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L4A, LINK4A, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L4B, LINK4B, L4BSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L4B, LINK4B, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L4C, LINK4C, L4CSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L4C, LINK4C, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L4D, LINK4D, L4DSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L4D, LINK4D, CLOSE_STAT )
      ENDIF

      IF (CLOSE_STAT == 'FILE_STAT') THEN
         CALL CLOSE_THIS_FILE ( L5B, LINK5B, L5BSTAT )
      ELSE
         CALL CLOSE_THIS_FILE ( L5B, LINK5B, CLOSE_STAT )
      ENDIF

      DO I=1,MOU4
         INQUIRE(FILE=OU4FIL(I),EXIST=LEXIST,OPENED=LOPND)
         IF (LEXIST) THEN
            IF (LOPND) THEN
               CALL FILE_CLOSE ( OU4(I), OU4FIL(I), OU4STAT(I) )
            ELSE
               OPEN (OU4(I),FILE=OU4FIL(I),STATUS='OLD',IOSTAT=IOCHK)
               IF (IOCHK /= 0) THEN
                  CALL OPNERR ( IOCHK, OU4FIL(I), OUNT )
                  CALL OUTA_HERE ( 'Y' )
               ELSE
                  CALL FILE_CLOSE ( OU4(I), OU4FIL(I), OU4STAT(I) )
               ENDIF
            ENDIF
         ENDIF
      ENDDO

! **********************************************************************************************************************************
  954 FORMAT(' *ERROR   954: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INTENT(IN) ARG CLOSE_STAT MUST BE "DELETE" OR "FILE_STAT" BUT IS "',A,'"')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CLOSE_THIS_FILE ( UNT, FILNAM, STATUS )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM
      CHARACTER(LEN=*), INTENT(IN)    :: STATUS

      INTEGER(LONG)   , INTENT(IN)    :: UNT

! **********************************************************************************************************************************

      INQUIRE(FILE=FILNAM,EXIST=LEXIST,OPENED=LOPND)
      IF (LEXIST) THEN
         IF (LOPND) THEN
            CALL FILE_CLOSE ( UNT, FILNAM, STATUS )
         ELSE
            OPEN (UNT,FILE=FILNAM,STATUS='OLD',IOSTAT=IOCHK)
            IF (IOCHK /= 0) THEN
               CALL OPNERR ( IOCHK, FILNAM, OUNT )
               CALL OUTA_HERE ( 'Y' )
            ELSE
               CALL FILE_CLOSE ( UNT, FILNAM, STATUS )
            ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE CLOSE_THIS_FILE

      END SUBROUTINE CLOSE_LIJFILES


      SUBROUTINE CLOSE_OUTFILES ( BUG_CLOSE_STAT, ERR_CLOSE_STAT, OP2_CLOSE_STAT )

! Closes BUGFIL, ERRFIL, F06FIL

      USE PENTIUM_II_KIND, ONLY       :  BYTE
      USE IOUNT1, ONLY                :  BUG   , ERR   , F06   , OP2   , SC1, BUGFIL, ERRFIL, F06FIL, OP2FIL

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: BUG_CLOSE_STAT    ! Input value for close status for BUG
      CHARACTER(LEN=*), INTENT(IN)    :: ERR_CLOSE_STAT    ! Input value for close status for ERR
      CHARACTER(LEN=*), INTENT(IN)    :: OP2_CLOSE_STAT    ! Input value for close status for OP2

! **********************************************************************************************************************************
      ! close standard output files first
      IF (F06 /= SC1) THEN
         CALL FILE_CLOSE ( F06, F06FIL, 'KEEP' )
      ENDIF

      IF (OP2 /= SC1) THEN
         CALL END_OP2_TABLES()
         CALL FILE_CLOSE ( OP2, OP2FIL, OP2_CLOSE_STAT )
      ENDIF

      ! close error/log files last
      IF (BUG /= SC1) THEN
         CALL FILE_CLOSE ( BUG, BUGFIL, BUG_CLOSE_STAT )
      ENDIF

      IF (ERR /= SC1) THEN
         CALL FILE_CLOSE ( ERR, ERRFIL, ERR_CLOSE_STAT )
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE CLOSE_OUTFILES


      SUBROUTINE OPEN_OUTFILES

! Opens BUGFIL, ERRFIL, F06FIL and, after checking STIME, closes the file so it can be reopened with APPEND.
! This subr is intended for opening these files

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG    , ERR    , F06    , SC1, BUGOUT, FILE_NAM_MAXLEN,                                  &
                                         BUGFIL , ERRFIL , F06FIL ,                                                                &
                                         BUG_MSG, ERR_MSG, F06_MSG
      USE TIMDAT, ONLY                :  STIME, TSEC

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: QUIT              ! If 'Y' quit

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR(4)           ! Error indicators
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to


! **********************************************************************************************************************************
! Default units for writing errors the screen (until LINK1A is read)

      OUNT(1) = SC1
      OUNT(2) = SC1

      QUIT = 'N'

      DO I=1,4
         IERR(I) = 0
      ENDDO

! Open BUG, ERR, F06 files, check STIME, and position file at end

      CALL OPENIT ( BUGFIL, BUG, 'BUG', BUG_MSG, IERR(2) )

      CALL OPENIT ( ERRFIL, ERR, 'ERR', ERR_MSG, IERR(3) )

      CALL OPENIT ( F06FIL, F06, 'F06', F06_MSG, IERR(4) )

! If there were any errors based on opening above files, quit.

      DO I=1,4
         IF (IERR(I) > 0) THEN
            CALL FILERR ( OUNT )
            QUIT = 'Y'
         ENDIF
      ENDDO

      IF (QUIT == 'Y') THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE OPENIT ( FILNAM, UNIT, UNIT_NAME, FILE_MSG, IERR0 )

! Opens formatted files

      USE IOUNT1, ONLY                        :  FILE_NAM_MAXLEN
      USE SCONTR, ONLY                        :  LINKNO, FATAL_ERR
      USE TIMDAT, ONLY                        :  MONTH,DAY,YEAR,HOUR,MINUTE,SEC,SFRAC

      IMPLICIT NONE

      LOGICAL                                 :: FILE_EXIST! Result from INQUIRE regarding whether a file is open
      LOGICAL                                 :: FILE_OPND ! Result from INQUIRE regarding whether a file is open

      CHARACTER(FILE_NAM_MAXLEN*BYTE), INTENT(IN)                                                                                  &
                                              :: FILNAM    ! File name (needs to be same length as other file names in MYSTRAN

      CHARACTER(LEN=*)                        :: FILE_MSG  ! Char message describing file
      CHARACTER(LEN=*)                        :: UNIT_NAME ! Char name of unit being opened (e.g. 'BUG')

      INTEGER(LONG)                           :: IERR0     ! Error indicator
      INTEGER(LONG), INTENT(IN)               :: UNIT      ! Unit number for the file to be opened

! **********************************************************************************************************************************
! FILNAM is opened, if not already opened, and positioned at the file end so we can append data.
! First, though, we need to make sure that
!   1) there is no problem opening the file and
!   2) that there is no error reading the time stamp and
!   3) the time stamp is the correct value.
! If this is the case, we close the file and reopen it with STATUS = 'APPEND'

      IERR0 = 0

      INQUIRE (FILE=FILNAM,EXIST=FILE_EXIST)                ! Make sure file exists (will be true if file OR unit number exist)
      IF (FILE_EXIST) THEN
         INQUIRE (FILE=FILNAM,OPENED=FILE_OPND)             ! If it is opened we assume it is already positioned at the end
         IF (.NOT.FILE_OPND) THEN
            IF (UNIT /= SC1) THEN                           ! Open file, check STIME, close file and then reopen as APPEND
               CALL FILE_OPEN ( UNIT, FILNAM, OUNT, 'OLD', FILE_MSG, 'READ_STIME', 'FORMATTED', 'READ'     , 'REWIND','Y','N')
               CALL FILE_CLOSE ( UNIT, FILNAM, 'KEEP' )
               CALL FILE_OPEN ( UNIT, FILNAM, OUNT, 'OLD', FILE_MSG, 'NEITHER'   , 'FORMATTED', 'READWRITE', 'APPEND','Y','N')
            ENDIF
         ENDIF
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         IERR0     = IERR0 + 1
         WRITE(ERR,939) UNIT_NAME,UNIT
         CALL WRITE_FILNAM ( FILNAM, ERR, 15 )
         WRITE(F06,939) UNIT_NAME,UNIT
         CALL WRITE_FILNAM ( FILNAM, F06, 15 )
      ENDIF

! **********************************************************************************************************************************
  939 FORMAT(' *ERROR   939: THE FOLLOWING OUTPUT FILE (CONNECTED TO FILE UNIT ',A,' = ',I3,') DOES NOT EXIST:')

! **********************************************************************************************************************************

      END SUBROUTINE OPENIT

      END SUBROUTINE OPEN_OUTFILES


      SUBROUTINE WRITE_FILNAM ( FILNAM, UNT, BLEN )

! Writes file name, FILNAM, to unit UNT with a format that depends on the actual length, in characters, of FILNAM (used
! to avoid using 4 lines of screen to write a short file name when the allowable length of FILNAM is 256 bytes)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, SC1

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name
      CHARACTER(15*BYTE)              :: LEADING_BLANKS    !

      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number of file FILNAM to write to
      INTEGER(LONG), INTENT(IN)       :: BLEN              ! Length, in char's, of blanks to print before filnam
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: FLEN              ! Length, in char's, of FILNAM (non-blanks)

! **********************************************************************************************************************************
! Find actual length of file name

      LEADING_BLANKS(1:) = ' '

! Print filename

      IF (UNT == SC1) THEN
         WRITE(SC1,1001) LEADING_BLANKS(1:BLEN), TRIM(FILNAM)
      ELSE
         WRITE(UNT,1001) LEADING_BLANKS(1:BLEN), TRIM(FILNAM)
      ENDIF

! **********************************************************************************************************************************
 1001 FORMAT(A,A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_FILNAM


      SUBROUTINE FILERR ( OUNT )

! Writes message about file open errors and errors reading data

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
      DO I=1,2
         WRITE(OUNT(I),900)
         IF (OUNT(2) == OUNT(1)) EXIT
      ENDDO

! **********************************************************************************************************************************
  900 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE FILE OPEN/READ ERRORS.')



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE FILERR


      SUBROUTINE READERR (IOCHK, FILNAM, MESSAG, REC_NO, OUNT )

! Writes message about errors encountered when reading files

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  SC1
      USE SCONTR, ONLY                :  FATAL_ERR

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description. Used for error messaging
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name

      INTEGER(LONG), INTENT(IN)       :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to
      INTEGER(LONG), INTENT(IN)       :: REC_NO            ! Indicator of record number when error encountered reading file
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IEND              ! End col for MESSAG




! **********************************************************************************************************************************
! First, truncate trailing blanks in MESSAG

      DO I=LEN(MESSAG),1,-1
         IF (MESSAG(I:I) == ' ') THEN
            CYCLE
         ELSE
            IEND = I
            EXIT
         ENDIF
      ENDDO

! IOCHK < 0 is due to EOF/EOR during read. IOCHK > 0 is due to error during read.

      DO I=1,2

         IF (OUNT(I) /= SC1) THEN

            IF (IOCHK < 0) THEN
               IF (REC_NO > 0) THEN
                  WRITE(OUNT(I),904) IOCHK, REC_NO, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15 )
               ELSE
                  WRITE(OUNT(I),905) IOCHK, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ENDIF
            ELSE IF (IOCHK > 0) THEN
               IF (REC_NO > 0) THEN
                  WRITE(OUNT(I),906) IOCHK, REC_NO, MESSAG(1:IEND)
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ELSE
                  WRITE(OUNT(I),907) IOCHK, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ENDIF
            ENDIF

         ELSE

            IF (IOCHK < 0) THEN
               IF (REC_NO > 0) THEN
                  WRITE(SC1,904) IOCHK, REC_NO, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ELSE
                  WRITE(SC1,905) IOCHK, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ENDIF
            ELSE IF (IOCHK > 0) THEN
               IF (REC_NO > 0) THEN
                  WRITE(SC1,906) IOCHK, REC_NO, MESSAG(1:IEND)
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ELSE
                  WRITE(SC1,907) IOCHK, MESSAG
                  CALL WRITE_FILNAM(FILNAM,OUNT(I), 15)
               ENDIF
            ENDIF

         ENDIF

         IF (OUNT(2) == OUNT(1)) THEN
            EXIT
         ENDIF

      ENDDO

      FATAL_ERR = FATAL_ERR + 1



      RETURN

! **********************************************************************************************************************************
  904 FORMAT(' *ERROR   904: EOF/EOR WITH IOSTAT = ',I8,' ENCOUNTERED READING FROM THE FOLLOWING FILE ON RECORD NUMBER ',I8,       &
                           ' FOR DATA NAMED "',A,'"')

  905 FORMAT(' *ERROR   905: EOF/EOR WITH IOSTAT = ',I8,' ENCOUNTERED READING FROM THE FOLLOWING FILE FOR DATA NAMED "',A,'"')

  906 FORMAT(' *ERROR   906: ERR WITH IOSTAT = ',I8,' ENCOUNTERED READING FROM THE FOLLOWING FILE ON RECORD NUMBER ',I8,           &
                           ' FOR DATA NAMED "',A,'"')

  907 FORMAT(' *ERROR   907: ERR WITH IOSTAT = ',I8,' ENCOUNTERED READING FROM THE FOLLOWING FILE FOR DATA NAMED "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE READERR


      SUBROUTINE STMERR ( XTIME, FILNAM, OUNT )

! Prints error messages when the wrong time stamp, STIME, is read as the first record in a file that has been opened

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN
      USE SCONTR, ONLY                :  FATAL_ERR
      USE TIMDAT, ONLY                :  STIME

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to
      INTEGER(LONG), INTENT(IN)       :: XTIME             ! Time stamp read from file LINK1A
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IEND              ! Index




! **********************************************************************************************************************************
      DO I=FILE_NAM_MAXLEN,1,-1
         IF (FILNAM(I:I) /= ' ') THEN
            IEND = I
            EXIT
         ENDIF
      ENDDO

      DO I=1,2
         WRITE(OUNT(I),908) XTIME, STIME, FILNAM(1:IEND)
         IF (OUNT(2) == OUNT(1)) EXIT
      ENDDO

      FATAL_ERR = FATAL_ERR + 1



      RETURN

! **********************************************************************************************************************************
  908 FORMAT(' *ERROR   908: WRONG SYSTEM TIME STAMP = ',I12,'. SHOULD BE = ',I12                                                  &
                    ,/,14X,' READ FROM FILE:'                                                                                      &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' IT COULD BE THAT THE FILE DOES NOT EXIST')

! **********************************************************************************************************************************
       END SUBROUTINE STMERR


      SUBROUTINE OPNERR ( IOCHK, FILNAM, OUNT )

! Prints error messages when IOSTAT is not zero on a file OPEN.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  FATAL_ERR, RESTART

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST        ! True if FILNAM exists
      LOGICAL                         :: FILE_OPENED       ! True if FILNAM is open

      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG), INTENT(IN)       :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to




! **********************************************************************************************************************************
! IOCHK < 0 is due to EOF/EOR during open. IOCHK > 0 is due to error during open.

      INQUIRE (FILE=FILNAM,OPENED=FILE_OPENED)
      INQUIRE (FILE=FILNAM, EXIST=FILE_EXIST)
      DO I=1,2

         WRITE(OUNT(I), * )

         IF (IOCHK < 0) THEN

            WRITE(OUNT(I),902) IOCHK
            CALL WRITE_FILNAM ( FILNAM, OUNT(I), 15 )

         ELSE

            IF (.NOT.FILE_EXIST) THEN

               WRITE(OUNT(I),903) IOCHK
               CALL WRITE_FILNAM ( FILNAM, OUNT(I), 15 )
               IF (RESTART == 'Y') THEN
                  WRITE(OUNT(I),9223)
               ELSE
                  WRITE(OUNT(I),9222)
               ENDIF

            ELSE IF (FILE_OPENED) THEN

               WRITE(OUNT(I),903) IOCHK
               CALL WRITE_FILNAM ( FILNAM, OUNT(I), 15 )
               WRITE(OUNT(I),9232)

            ELSE

               WRITE(OUNT(I),903) IOCHK
               CALL WRITE_FILNAM ( FILNAM, OUNT(I), 15 )
               WRITE(OUNT(I),9242)

            ENDIF

         ENDIF

         WRITE(OUNT(I), * )

         IF (OUNT(2) == OUNT(1)) EXIT

      ENDDO

      FATAL_ERR = FATAL_ERR + 1

      IF (FILE_OPENED) THEN
         CALL OUTA_HERE ( 'N' )
      ENDIF


      RETURN

! **********************************************************************************************************************************
  902 FORMAT(' *ERROR   902: EOF/EOR ENCOUNTERED WITH IOSTAT = ',I8,' OPENING FILE:')

  903 FORMAT(' *ERROR   903: ERROR ENCOUNTERED WITH IOSTAT = ',I8,' OPENING FILE:')

 9222 FORMAT('               THE FILE DOES NOT EXIST. THIS IS A PROGRAMMING ERROR.',/)

 9223 FORMAT('               THE FILE COULD NOT BE FOUND AND IS REQUIRED IN A RESTART. MAKE SURE THE ORIGINAL RUN PRODUCED AND',   &
                           ' SAVED THIS FILE',/)

 9232 FORMAT('               IT MAY BE OPEN IN ANOTHER PROGRAM. IF SO, CLOSE IT & START AGAIN.')

 9242 FORMAT('               THE FILE EXISTS AND IS NOT OPENED.  THIS IS A PROGRAMMING ERROR.')

! **********************************************************************************************************************************

      END SUBROUTINE OPNERR


      SUBROUTINE OUTA_HERE ( WRITE_TO_L1A )

! Routine called when MYSTRAN has to stop due to some error. Writes the termination message and stops.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  BUGOUT, F06, F06FIL, SC1,                                                 &
                                         BUGSTAT, BUGSTAT_OLD, ERRSTAT, ERRSTAT_OLD,                          &
                                         OP2STAT, L1A, LINK1A, L1ASTAT

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LINKNO, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE


      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OUTA_HERE'
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_TO_L1A      ! Y/N indicator of whether to call subr WRITE_L1A
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT status while deleting LINK1A
      LOGICAL                         :: L1A_OPEN          ! Whether LINK1A is currently open
      LOGICAL                         :: L1A_EXIST         ! Whether LINK1A exists on disk
      LOGICAL                         :: F06_OPEN          ! Whether the F06 file is open



! **********************************************************************************************************************************


! Write termination messages and close output files

      IF (WRITE_TO_L1A == 'Y') THEN

         IF (FATAL_ERR > 0) THEN                           ! Check fatal error flag and write message
            WRITE(SC1,9992) FATAL_ERR
            INQUIRE ( UNIT=F06, OPENED=F06_OPEN )          ! Every stop for fatal errors ends the F06 with the same line
            IF (F06_OPEN) WRITE(F06,9991) FATAL_ERR
         ENDIF

         IF (WARN_ERR > 0) THEN                            ! Check warning flag and write message
            IF (SUPWARN == 'Y') THEN
               WRITE(SC1,9994) WARN_ERR
            ELSE
               WRITE(SC1,9993) WARN_ERR
            ENDIF
         ENDIF

         WRITE(SC1,9999)
         CALL WRITE_FILNAM ( F06FIL, SC1, 1 )


! Set close status for output files

         IF (BUGOUT == 'Y') THEN
            BUGSTAT = 'KEEP'
         ELSE
            IF (BUGSTAT_OLD == 'KEEP    ') THEN
               BUGSTAT = 'KEEP'
            ELSE
               BUGSTAT = 'DELETE'
            ENDIF
         ENDIF

         IF ((FATAL_ERR > 0) .OR. (WARN_ERR > 0)) THEN
            ERRSTAT = 'KEEP'
         ELSE
            IF (ERRSTAT_OLD == 'KEEP    ') THEN
               ERRSTAT = 'KEEP'
            ELSE
               ERRSTAT = 'DELETE'
            ENDIF
         ENDIF


      ELSE

         WRITE(SC1,9999)
         CALL WRITE_FILNAM ( F06FIL, SC1, 1 )

      ENDIF

      INQUIRE ( UNIT=L1A, OPENED=L1A_OPEN )
      INQUIRE ( FILE=LINK1A, EXIST=L1A_EXIST )
      IF (L1A_OPEN) THEN
         CLOSE ( L1A, STATUS='DELETE', IOSTAT=IOCHK )
      ELSE IF (L1A_EXIST) THEN
         OPEN ( L1A, FILE=LINK1A, STATUS='OLD', IOSTAT=IOCHK )
         IF (IOCHK == 0) THEN
            CLOSE ( L1A, STATUS='DELETE', IOSTAT=IOCHK )
         ENDIF
      ENDIF

      CALL CLOSE_OUTFILES ( BUGSTAT, ERRSTAT, OP2STAT )

      STOP

! **********************************************************************************************************************************
 9991 FORMAT(/,' *** PROCESSING TERMINATED: ',I0,' FATAL ERROR(S), LISTED ABOVE')

 9992 FORMAT(' CHECK F06 OUTPUT FILE FOR ',I8,' FATAL MESSAGE(S)')

 9993 FORMAT(' CHECK F06 OUTPUT FILE FOR ',I8,' WARNING MESSAGE(S)')

 9994 FORMAT(' CHECK ERR OUTPUT FILE FOR ',I8,' WARNING MESSAGE(S)')

 9999 FORMAT(' Processing terminated. Check for error messages in output file:')

! **********************************************************************************************************************************

      END SUBROUTINE OUTA_HERE


      SUBROUTINE WRITE_L1A ( CLOSE_STAT, CALL_OUTA_HERE )

! Writes data to file LINK1A at the end of each LINK. This is read by all LINK's after LINK1, as they begin. This text file contains
! the names of files opened for a run, the "counter" info (e.g. NGRID, number of grids, etc), solution number, PARAM's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  MOT4,    MOU4,    WRT_ERR

      USE IOUNT1, ONLY                :  BUG,     EIN,     ENF,     ERR,     F06,     IN0,     IN1,     INI,                       &
                                         L1A,     NEU,     OT4,     SEQ,     SPC,     SC1,                                &
                                         F21,     F22,     F23,     F24,     F25,                                                  &
                                         L1B,     L1C,     L1D,     L1E,     L1F,     L1G,     L1H,     L1I,     L1J,     L1K,     &
                                         L1L,     L1M,     L1N,     L1O,     L1P,     L1Q,     L1R,     L1S,     L1T,     L1U,     &
                                         L1V,     L1W,     L1X,     L1Y,     L1Z,                                                  &
                                         L2A,     L2B,     L2C,     L2D,     L2E,     L2F,     L2G,     L2H,     L2I,     L2J,     &
                                         L2K,     L2L,     L2M,     L2N,     L2O,     L2P,     L2Q,     L2R,     L2S,     L2T,     &
                                         L3A,     L4A,     L4B,     L4C,     L4D,     L5A,     L5B,     OP2,     OU4

      USE IOUNT1, ONLY                :  BUGSTAT, EINSTAT, ENFSTAT, ERRSTAT, F06STAT, IN0STAT, IN1STAT, INISTAT,                   &
                                         L1ASTAT, NEUSTAT, OT4STAT, SEQSTAT, SPCSTAT,                                              &
                                         F21STAT, F22STAT, F23STAT, F24STAT, F25STAT,                                              &
                                         L1BSTAT, L1CSTAT, L1DSTAT, L1ESTAT, L1FSTAT, L1GSTAT, L1HSTAT, L1ISTAT, L1JSTAT, L1KSTAT, &
                                         L1LSTAT, L1MSTAT, L1NSTAT, L1OSTAT, L1PSTAT, L1QSTAT, L1RSTAT, L1SSTAT, L1TSTAT, L1USTAT, &
                                         L1VSTAT, L1WSTAT, L1XSTAT, L1YSTAT, L1ZSTAT,                                              &
                                         L2ASTAT, L2BSTAT, L2CSTAT, L2DSTAT, L2ESTAT, L2FSTAT, L2GSTAT, L2HSTAT, L2ISTAT, L2JSTAT, &
                                         L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, L2OSTAT, L2PSTAT, L2QSTAT, L2RSTAT, L2SSTAT, L2TSTAT, &
                                         L3ASTAT, L4ASTAT, L4BSTAT, L4CSTAT, L4DSTAT, L5ASTAT, L5BSTAT, OP2STAT, OU4STAT

      USE IOUNT1, ONLY                :  BUGFIL,  EINFIL,  ENFFIL,  ERRFIL,  F06FIL,  IN0FIL,  INIFIL,  LINK1A,                    &
                                         NEUFIL,  OT4FIL,  SEQFIL,  SPCFIL,  F21FIL,  F22FIL,  F23FIL,  F24FIL,  F25FIL,           &
                                         LINK1A,  LINK1B,  LINK1C,  LINK1D,  LINK1E,  LINK1F,  LINK1G,  LINK1H,  LINK1I,  LINK1J,  &
                                         LINK1K,  LINK1L,  LINK1M,  LINK1N,  LINK1O,  LINK1P,  LINK1Q,  LINK1R,  LINK1S,  LINK1T,  &
                                         LINK1U,  LINK1V,  LINK1W,  LINK1X,  LINK1Y,  LINK1Z,                                      &
                                         LINK2A,  LINK2B,  LINK2C,  LINK2D,  LINK2E,  LINK2F,  LINK2G,  LINK2H,  LINK2I,  LINK2J,  &
                                         LINK2K,  LINK2L,  LINK2M,  LINK2N,  LINK2O,  LINK2P,  LINK2Q,  LINK2R,  LINK2S,  LINK2T,  &
                                         LINK3A,  LINK4A,  LINK4B,  LINK4C,  LINK4D,  LINK5A,  LINK5B,  OP2FIL,  OU4FIL

      USE IOUNT1, ONLY                :  BUG_MSG, EIN_MSG, ENF_MSG, ERR_MSG, F06_MSG, IN0_MSG, IN1_MSG, INI_MSG,                   &
                                         L1A_MSG, NEU_MSG, OT4_MSG, SEQ_MSG, SPC_MSG,                                              &
                                         F21_MSG, F22_MSG, F23_MSG, F24_MSG, F25_MSG,                                              &
                                         L1B_MSG, L1C_MSG, L1D_MSG, L1E_MSG, L1F_MSG, L1G_MSG, L1H_MSG, L1I_MSG, L1J_MSG, L1K_MSG, &
                                         L1L_MSG, L1M_MSG, L1N_MSG, L1O_MSG, L1P_MSG, L1Q_MSG, L1R_MSG, L1S_MSG, L1T_MSG, L1U_MSG, &
                                         L1V_MSG, L1W_MSG, L1X_MSG, L1Y_MSG, L1Z_MSG,                                              &
                                         L2A_MSG, L2B_MSG, L2C_MSG, L2D_MSG, L2E_MSG, L2F_MSG, L2G_MSG, L2H_MSG, L2I_MSG, L2J_MSG, &
                                         L2K_MSG, L2L_MSG, L2M_MSG, L2N_MSG, L2O_MSG, L2P_MSG, L2Q_MSG, L2R_MSG, L2S_MSG, L2T_MSG, &
                                         L3A_MSG, L4A_MSG, L4B_MSG, L4C_MSG, L4D_MSG, L5A_MSG, L5B_MSG, OP2_MSG, OU4_MSG
      USE SCONTR
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, ELFORCEN, HEXAXIS, IORQ1B, IORQ1M, IORQ1S, IORQ2B, IORQ2T,&
                                         MATSPARS, MIN4TRED, QUAD4TYP, QUADAXIS, SPARSTOR

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_L1A'
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! STATUS when closing file LINK1A
      CHARACTER(LEN=*), INTENT(IN)    :: CALL_OUTA_HERE    ! 'Y'/'N' indicator of whether to call OUTA_HERE (this should be 'Y'
!                                                             except when this subr is called by OUTA_HERE

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to




! **********************************************************************************************************************************
! Units for writing open errors

      OUNT(1) = ERR
      OUNT(2) = F06

! Open L1A file and write STIME

      OPEN (L1A,FILE=LINK1A,STATUS='REPLACE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR (IOCHK, LINK1A, OUNT )
         CALL FILERR ( OUNT )
         IF (CALL_OUTA_HERE == 'Y') THEN
            CALL OUTA_HERE ( 'N' )
         ENDIF
      ENDIF
      WRITE(L1A,110) STIME

! Write current LINK number

      WRITE(L1A,110) LINKNO

! Write solution name

      WRITE(L1A,120) SOL_NAME

! Write I/0 unit numbers, close status and names

      WRITE(L1A,140) SC1                                   !   1

      WRITE(L1A,151) BUG,BUGSTAT,BUG_MSG,BUGFIL            !   3
      WRITE(L1A,151) EIN,EINSTAT,EIN_MSG,EINFIL            !   4
      WRITE(L1A,151) ENF,ENFSTAT,ENF_MSG,ENFFIL            !   5
      WRITE(L1A,151) ERR,ERRSTAT,ERR_MSG,ERRFIL            !   6
      WRITE(L1A,151) F06,F06STAT,F06_MSG,F06FIL            !   8
      WRITE(L1A,151) IN0,IN0STAT,IN0_MSG,IN0FIL            !   9
      WRITE(L1A,151) L1A,L1ASTAT,L1A_MSG,LINK1A            !  10
      WRITE(L1A,151) NEU,NEUSTAT,NEU_MSG,NEUFIL            !  11

      WRITE(L1A,151) SEQ,SEQSTAT,SEQ_MSG,SEQFIL            !  13
      WRITE(L1A,151) SPC,SPCSTAT,SPC_MSG,SPCFIL            !  14

      WRITE(L1A,151) F21,F21STAT,F21_MSG,F21FIL            ! 15
      WRITE(L1A,151) F22,F22STAT,F22_MSG,F22FIL            ! 16
      WRITE(L1A,151) F23,F23STAT,F23_MSG,F23FIL            ! 17
      WRITE(L1A,151) F24,F24STAT,F24_MSG,F24FIL            ! 18
      WRITE(L1A,151) F25,F25STAT,F25_MSG,F25FIL            ! 19
      WRITE(L1A,151) L1B,L1BSTAT,L1B_MSG,LINK1B            ! 20
      WRITE(L1A,151) L1C,L1CSTAT,L1C_MSG,LINK1C            ! 21
      WRITE(L1A,151) L1D,L1DSTAT,L1D_MSG,LINK1D            ! 22
      WRITE(L1A,151) L1E,L1ESTAT,L1E_MSG,LINK1E            ! 23
      WRITE(L1A,151) L1F,L1FSTAT,L1F_MSG,LINK1F            ! 24
      WRITE(L1A,151) L1G,L1GSTAT,L1G_MSG,LINK1G            ! 25
      WRITE(L1A,151) L1H,L1HSTAT,L1H_MSG,LINK1H            ! 26
      WRITE(L1A,151) L1I,L1ISTAT,L1I_MSG,LINK1I            ! 27
      WRITE(L1A,151) L1J,L1JSTAT,L1J_MSG,LINK1J            ! 28
      WRITE(L1A,151) L1K,L1KSTAT,L1K_MSG,LINK1K            ! 29
      WRITE(L1A,151) L1L,L1LSTAT,L1L_MSG,LINK1L            ! 30
      WRITE(L1A,151) L1M,L1MSTAT,L1M_MSG,LINK1M            ! 31
      WRITE(L1A,151) L1N,L1NSTAT,L1N_MSG,LINK1N            ! 32
      WRITE(L1A,151) L1O,L1OSTAT,L1O_MSG,LINK1O            ! 33
      WRITE(L1A,151) L1P,L1PSTAT,L1P_MSG,LINK1P            ! 34
      WRITE(L1A,151) L1Q,L1QSTAT,L1Q_MSG,LINK1Q            ! 35
      WRITE(L1A,151) L1R,L1RSTAT,L1R_MSG,LINK1R            ! 36
      WRITE(L1A,151) L1S,L1SSTAT,L1S_MSG,LINK1S            ! 37
      WRITE(L1A,151) L1T,L1TSTAT,L1T_MSG,LINK1T            ! 38
      WRITE(L1A,151) L1U,L1USTAT,L1U_MSG,LINK1U            ! 39
      WRITE(L1A,151) L1V,L1VSTAT,L1V_MSG,LINK1V            ! 40
      WRITE(L1A,151) L1W,L1WSTAT,L1W_MSG,LINK1W            ! 41
      WRITE(L1A,151) L1W,L1WSTAT,L1X_MSG,LINK1X            ! 42
      WRITE(L1A,151) L1Y,L1YSTAT,L1Y_MSG,LINK1Y            ! 43
      WRITE(L1A,151) L1Z,L1ZSTAT,L1Z_MSG,LINK1Z            ! 44
      WRITE(L1A,151) L2A,L2ASTAT,L2A_MSG,LINK2A            ! 45
      WRITE(L1A,151) L2B,L2BSTAT,L2B_MSG,LINK2B            ! 46
      WRITE(L1A,151) L2C,L2CSTAT,L2C_MSG,LINK2C            ! 47
      WRITE(L1A,151) L2D,L2DSTAT,L2D_MSG,LINK2D            ! 48
      WRITE(L1A,151) L2E,L2ESTAT,L2E_MSG,LINK2E            ! 49
      WRITE(L1A,151) L2F,L2FSTAT,L2F_MSG,LINK2F            ! 50
      WRITE(L1A,151) L2G,L2GSTAT,L2G_MSG,LINK2G            ! 51
      WRITE(L1A,151) L2H,L2HSTAT,L2H_MSG,LINK2H            ! 52
      WRITE(L1A,151) L2I,L2ISTAT,L2I_MSG,LINK2I            ! 53
      WRITE(L1A,151) L2J,L2JSTAT,L2J_MSG,LINK2J            ! 54
      WRITE(L1A,151) L2K,L2KSTAT,L2K_MSG,LINK2K            ! 55
      WRITE(L1A,151) L2L,L2LSTAT,L2L_MSG,LINK2L            ! 56
      WRITE(L1A,151) L2M,L2MSTAT,L2M_MSG,LINK2M            ! 57
      WRITE(L1A,151) L2N,L2NSTAT,L2N_MSG,LINK2N            ! 58
      WRITE(L1A,151) L2O,L2OSTAT,L2O_MSG,LINK2O            ! 59
      WRITE(L1A,151) L2P,L2PSTAT,L2P_MSG,LINK2P            ! 60
      WRITE(L1A,151) L2Q,L2QSTAT,L2Q_MSG,LINK2Q            ! 61
      WRITE(L1A,151) L2R,L2QSTAT,L2R_MSG,LINK2R            ! 62
      WRITE(L1A,151) L2S,L2QSTAT,L2S_MSG,LINK2S            ! 63
      WRITE(L1A,151) L2T,L2QSTAT,L2T_MSG,LINK2T            ! 64
      WRITE(L1A,151) L3A,L3ASTAT,L3A_MSG,LINK3A            ! 65
      WRITE(L1A,151) L4A,L4ASTAT,L4A_MSG,LINK4A            ! 66
      WRITE(L1A,151) L4B,L4BSTAT,L4B_MSG,LINK4B            ! 67
      WRITE(L1A,151) L4C,L4CSTAT,L4C_MSG,LINK4C            ! 68
      WRITE(L1A,151) L4D,L4DSTAT,L4D_MSG,LINK4D            ! 69
      WRITE(L1A,151) L5A,L5ASTAT,L5A_MSG,LINK5A            ! 70
      WRITE(L1A,151) L5B,L5BSTAT,L5B_MSG,LINK5B            ! 71
      WRITE(L1A,151) OP2,OP2STAT,OP2_MSG,OP2FIL            ! 72
      DO I=1,MOU4
         WRITE(L1A,151) OU4(I),OU4STAT(I),OU4FIL(I)        ! 72+I
      ENDDO
      DO I=1,MOT4
         WRITE(L1A,151) OT4(I),OT4STAT(I),OT4FIL(I)        ! 72+MOU4+I
      ENDDO

! Write counter info from module SCONTR

      I = 0

      I = I + 1  ;     WRITE(L1A,160) LBAROFF            , 'LBAROFF               (  1)'  !
      I = I + 1  ;     WRITE(L1A,160) LBUSHOFF           , 'LBUSHOFF              (  2)'  !
      I = I + 1  ;     WRITE(L1A,160) LCMASS             , 'LCMASS                (  3)'  !
      I = I + 1  ;     WRITE(L1A,160) LCONM2             , 'LCONM2                (  4)'  !
      I = I + 1  ;     WRITE(L1A,160) LCORD              , 'LCORD                 (  5)'  !
      I = I + 1  ;     WRITE(L1A,160) LDOFG              , 'LDOFG                 (  6)'  !
      I = I + 1  ;     WRITE(L1A,160) LEDAT              , 'LEDAT                 (  7)'  !
      I = I + 1  ;     WRITE(L1A,160) LELE               , 'LELE                  (  8)'  !
      I = I + 1  ;     WRITE(L1A,160) LFORCE             , 'LFORCE                (  9)'  !
      I = I + 1  ;     WRITE(L1A,160) LGRAV              , 'LGRAV                 ( 10)'  !
      I = I + 1  ;     WRITE(L1A,160) LGRID              , 'LGRID                 ( 11)'  !
      I = I + 1  ;     WRITE(L1A,160) LGUSERIN           , 'LGUSERIN              ( 12)'  !
      I = I + 1  ;     WRITE(L1A,160) LIND_GRDS_MPCS     , 'LIND_GRDS_MPCS        ( 13)'  !
      I = I + 1  ;     WRITE(L1A,160) LLOADC             , 'LLOADC                ( 14)'  !
      I = I + 1  ;     WRITE(L1A,160) LLOADR             , 'LLOADR                ( 15)'  !
      I = I + 1  ;     WRITE(L1A,160) LMATANGLE          , 'LMATANGLE             ( 16)'  !
      I = I + 1  ;     WRITE(L1A,160) LMATL              , 'LMATL                 ( 17)'  !
      I = I + 1  ;     WRITE(L1A,160) LMPC               , 'LMPC                  ( 18)'  !
      I = I + 1  ;     WRITE(L1A,160) LMPCADDC           , 'LMPCADDC              ( 19)'  !
      I = I + 1  ;     WRITE(L1A,160) LMPCADDR           , 'LMPCADDR              ( 20)'  !
      I = I + 1  ;     WRITE(L1A,160) LPBAR              , 'LPBAR                 ( 21)'  !
      I = I + 1  ;     WRITE(L1A,160) LPBEAM             , 'LPBEAM                ( 22)'  !
      I = I + 1  ;     WRITE(L1A,160) LPBUSH             , 'LPBUSH                ( 23)'  !
      I = I + 1  ;     WRITE(L1A,160) LPCOMP             , 'LPCOMP                ( 24)'  !
      I = I + 1  ;     WRITE(L1A,160) LPCOMP_PLIES       , 'LPCOMP_PLIES          ( 25)'  !
      I = I + 1  ;     WRITE(L1A,160) LPDAT              , 'LPDAT                 ( 26)'  !
      I = I + 1  ;     WRITE(L1A,160) LPELAS             , 'LPELAS                ( 27)'  !
      I = I + 1  ;     WRITE(L1A,160) LPLATEOFF          , 'LPLATEOFF             ( 28)'  !
      I = I + 1  ;     WRITE(L1A,160) LPLATETHICK        , 'LPLATETHICK           ( 29)'  !
      I = I + 1  ;     WRITE(L1A,160) LPLOAD             , 'LPLOAD                ( 30)'  !
      I = I + 1  ;     WRITE(L1A,160) LPLOTEL            , 'LPLOTEL               ( 31)'  !
      I = I + 1  ;     WRITE(L1A,160) LPMASS             , 'LPMASS                ( 32)'  !
      I = I + 1  ;     WRITE(L1A,160) LPROD              , 'LPROD                 ( 33)'  !
      I = I + 1  ;     WRITE(L1A,160) LPSHEAR            , 'LPSHEAR               ( 34)'  !
      I = I + 1  ;     WRITE(L1A,160) LPSHEL             , 'LPSHEL                ( 35)'  !
      I = I + 1  ;     WRITE(L1A,160) LPSOLID            , 'LPSOLID               ( 36)'  !
      I = I + 1  ;     WRITE(L1A,160) LPUSER1            , 'LPUSER1               ( 37)'  !
      I = I + 1  ;     WRITE(L1A,160) LPUSERIN           , 'LPUSERIN              ( 38)'  !
      I = I + 1  ;     WRITE(L1A,160) LRFORCE            , 'LRFORCE               ( 39)'  !
      I = I + 1  ;     WRITE(L1A,160) LRIGEL             , 'LRIGEL                ( 40)'  !
      I = I + 1  ;     WRITE(L1A,160) LSEQ               , 'LSEQ                  ( 41)'  !
      I = I + 1  ;     WRITE(L1A,160) LSETLN             , 'LSETLN                ( 42)'  !
      I = I + 1  ;     WRITE(L1A,160) LSETS              , 'LSETS                 ( 43)'  !
      I = I + 1  ;     WRITE(L1A,160) LSLOAD             , 'LSLOAD                ( 44)'  !
      I = I + 1  ;     WRITE(L1A,160) LSPC               , 'LSPC                  ( 45)'  !
      I = I + 1  ;     WRITE(L1A,160) LSPC1              , 'LSPC1                 ( 46)'  !
      I = I + 1  ;     WRITE(L1A,160) LSPCADDC           , 'LSPCADDC              ( 47)'  !
      I = I + 1  ;     WRITE(L1A,160) LSPCADDR           , 'LSPCADDR              ( 48)'  !
      I = I + 1  ;     WRITE(L1A,160) LSUB               , 'LSUB                  ( 49)'  !
      I = I + 1  ;     WRITE(L1A,160) LSUSERIN           , 'LSUSERIN              ( 50)'  !
      I = I + 1  ;     WRITE(L1A,160) LTDAT              , 'LTDAT                 ( 51)'  !
      I = I + 1  ;     WRITE(L1A,160) LTERM_KGG          , 'LTERM_KGG             ( 52)'  !
      I = I + 1  ;     WRITE(L1A,160) LTERM_KGGD         , 'LTERM_KGGD            ( 53)'  !
      I = I + 1  ;     WRITE(L1A,160) LTERM_MGGE         , 'LTERM_MGGE            ( 54)'  !
      I = I + 1  ;     WRITE(L1A,160) LVVEC              , 'LVVEC                 ( 55)'  !
      I = I + 1  ;     WRITE(L1A,160) MAX_ELEM_DEGREE    , 'MAX_ELEM_DEGREE       ( 56)'  !
      I = I + 1  ;     WRITE(L1A,160) MAX_GAUSS_POINTS   , 'MAX_GAUSS_POINTS      ( 57)'  !
      I = I + 1  ;     WRITE(L1A,160) MAX_STRESS_POINTS  , 'MAX_STRESS_POINTS     ( 58)'  !
      I = I + 1  ;     WRITE(L1A,160) MDT                , 'MDT                   ( 59)'  !
      I = I + 1  ;     WRITE(L1A,160) MELGP              , 'MELGP                 ( 60)'  !
      I = I + 1  ;     WRITE(L1A,160) MELDOF             , 'MELDOF                ( 61)'  !
      I = I + 1  ;     WRITE(L1A,160) MID1_PCOMP_EQ      , 'MID1_PCOMP_EQ         ( 62)'  !
      I = I + 1  ;     WRITE(L1A,160) MID2_PCOMP_EQ      , 'MID2_PCOMP_EQ         ( 63)'  !
      I = I + 1  ;     WRITE(L1A,160) MID3_PCOMP_EQ      , 'MID3_PCOMP_EQ         ( 64)'  !
      I = I + 1  ;     WRITE(L1A,160) MID4_PCOMP_EQ      , 'MID4_PCOMP_EQ         ( 65)'  !
      I = I + 1  ;     WRITE(L1A,160) MLL_SDIA           , 'MLL_SDIA              ( 66)'  !
      I = I + 1  ;     WRITE(L1A,160) MMPC               , 'MMPC                  ( 67)'  !
      I = I + 1  ;     WRITE(L1A,160) MOFFSET            , 'MOFFSET               ( 68)'  !
      I = I + 1  ;     WRITE(L1A,160) MRBE3              , 'MRBE3                 ( 69)'  !
      I = I + 1  ;     WRITE(L1A,160) MRSPLINE           , 'MRSPLINE              ( 70)'  !
      I = I + 1  ;     WRITE(L1A,160) NAOCARD            , 'NAOCARD               ( 71)'  !
      I = I + 1  ;     WRITE(L1A,160) NBAROFF            , 'NBAROFF               ( 72)'  !
      I = I + 1  ;     WRITE(L1A,160) NBUSHOFF           , 'NBUSHOFF              ( 73)'  !
      I = I + 1  ;     WRITE(L1A,160) NBAROR             , 'NBAROR                ( 74)'  !
      I = I + 1  ;     WRITE(L1A,160) NBEAMOR            , 'NBEAMOR               ( 75)'  !
      I = I + 1  ;     WRITE(L1A,160) NCBAR              , 'NCBAR                 ( 76)'  !
      I = I + 1  ;     WRITE(L1A,160) NCBEAM             , 'NCBEAM                ( 77)'  !
      I = I + 1  ;     WRITE(L1A,160) NCBUSH             , 'NCBUSH                ( 78)'  !
      I = I + 1  ;     WRITE(L1A,160) NCELAS1            , 'NCELAS1               ( 79)'  !
      I = I + 1  ;     WRITE(L1A,160) NCELAS2            , 'NCELAS2               ( 80)'  !
      I = I + 1  ;     WRITE(L1A,160) NCELAS3            , 'NCELAS3               ( 81)'  !
      I = I + 1  ;     WRITE(L1A,160) NC_INFILE          , 'NC_INFILE             ( 82)'  !
      I = I + 1  ;     WRITE(L1A,160) NCELAS4            , 'NCELAS4               ( 83)'  !
      I = I + 1  ;     WRITE(L1A,160) NCHEXA8            , 'NCHEXA8               ( 84)'  !
      I = I + 1  ;     WRITE(L1A,160) NCHEXA20           , 'NCHEXA20              ( 85)'  !
      I = I + 1  ;     WRITE(L1A,160) NCMASS             , 'NCMASS                ( 86)'  !
      I = I + 1  ;     WRITE(L1A,160) NCONM2             , 'NCONM2                ( 87)'  !
      I = I + 1  ;     WRITE(L1A,160) NCORD              , 'NCORD                 ( 88)'  !
      I = I + 1  ;     WRITE(L1A,160) NCORD1             , 'NCORD1                ( 89)'  !
      I = I + 1  ;     WRITE(L1A,160) NCORD2             , 'NCORD2                ( 90)'  !
      I = I + 1  ;     WRITE(L1A,160) NCPENTA6           , 'NCPENTA6              ( 91)'  !
      I = I + 1  ;     WRITE(L1A,160) NCPENTA15          , 'NCPENTA15             ( 92)'  !
      I = I + 1  ;     WRITE(L1A,160) NCQUAD4            , 'NCQUAD4               ( 93)'  !
      I = I + 1  ;     WRITE(L1A,160) NCQUAD4K           , 'NCQUAD4K              ( 94)'  !
      I = I + 1  ;     WRITE(L1A,160) NCROD              , 'NCROD                 ( 95)'  !
      I = I + 1  ;     WRITE(L1A,160) NCSHEAR            , 'NCSHEAR               ( 96)'  !
      I = I + 1  ;     WRITE(L1A,160) NCTETRA4           , 'NCTETRA4              ( 97)'  !
      I = I + 1  ;     WRITE(L1A,160) NCTETRA10          , 'NCTETRA10             ( 98)'  !
      I = I + 1  ;     WRITE(L1A,160) NCTRIA3            , 'NCTRIA3               ( 99)'  !
      I = I + 1  ;     WRITE(L1A,160) NCTRIA3K           , 'NCTRIA3K              (100)'  !
      I = I + 1  ;     WRITE(L1A,160) NCUSER1            , 'NCUSER1               (101)'  !
      I = I + 1  ;     WRITE(L1A,160) NCUSERIN           , 'NCUSERIN              (102)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFA              , 'NDOFA                 (103)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOF_EIG           , 'NDOF_EIG              (104)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFF              , 'NDOFF                 (105)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFG              , 'NDOFG                 (106)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFL              , 'NDOFL                 (107)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFM              , 'NDOFM                 (108)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFN              , 'NDOFN                 (109)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFO              , 'NDOFO                 (110)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFR              , 'NDOFR                 (111)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFS              , 'NDOFS                 (112)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFSA             , 'NDOFSA                (113)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFSB             , 'NDOFSB                (114)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFSE             , 'NDOFSE                (115)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFSG             , 'NDOFSG                (116)'  !
      I = I + 1  ;     WRITE(L1A,160) NDOFSZ             , 'NDOFSZ                (117)'  !
      I = I + 1  ;     WRITE(L1A,160) NEDAT              , 'NEDAT                 (118)'  !
      I = I + 1  ;     WRITE(L1A,160) NELE               , 'NELE                  (119)'  !
      I = I + 1  ;     WRITE(L1A,160) NFORCE             , 'NFORCE                (120)'  !
      I = I + 1  ;     WRITE(L1A,160) NGRAV              , 'NGRAV                 (121)'  !
      I = I + 1  ;     WRITE(L1A,160) NGRDSET            , 'NGRDSET               (122)'  !
      I = I + 1  ;     WRITE(L1A,160) NGRID              , 'NGRID                 (123)'  !
      I = I + 1  ;     WRITE(L1A,160) NIND_GRDS_MPCS     , 'NIND_GRDS_MPCS        (124)'  !
      I = I + 1  ;     WRITE(L1A,160) NLOAD              , 'NLOAD                 (125)'  !
      I = I + 1  ;     WRITE(L1A,160) NMATANGLE          , 'NMATANGLE             (126)'  !
      I = I + 1  ;     WRITE(L1A,160) NMATL              , 'NMATL                 (127)'  !
      I = I + 1  ;     WRITE(L1A,160) NMPC               , 'NMPC                  (128)'  !
      I = I + 1  ;     WRITE(L1A,160) NMPCADD            , 'NMPCADD               (129)'  !
      I = I + 1  ;     WRITE(L1A,160) NPBAR              , 'NPBAR                 (130)'  !
      I = I + 1  ;     WRITE(L1A,160) NPBARL             , 'NPBARL                (131)'  !
      I = I + 1  ;     WRITE(L1A,160) NPBEAM             , 'NPBEAM                (132)'  !
      I = I + 1  ;     WRITE(L1A,160) NPBUSH             , 'NPBUSH                (133)'  !
      I = I + 1  ;     WRITE(L1A,160) NPCARD             , 'NPCARD                (134)'  !
      I = I + 1  ;     WRITE(L1A,160) NPCOMP             , 'NPCOMP                (135)'  !
      I = I + 1  ;     WRITE(L1A,160) NPDAT              , 'NPDAT                 (136)'  !
      I = I + 1  ;     WRITE(L1A,160) NPELAS             , 'NPELAS                (137)'  !
      I = I + 1  ;     WRITE(L1A,160) NPLATEOFF          , 'NPLATEOFF             (138)'  !
      I = I + 1  ;     WRITE(L1A,160) NPLATETHICK        , 'NPLATETHICK           (139)'  !
      I = I + 1  ;     WRITE(L1A,160) NPLOTEL            , 'NPLOTEL               (140)'  !
      I = I + 1  ;     WRITE(L1A,160) NPLOAD             , 'NPLOAD                (141)'  !
      I = I + 1  ;     WRITE(L1A,160) NPLOAD4_3D         , 'NPLOAD4_3D            (142)'  !
      I = I + 1  ;     WRITE(L1A,160) NPMASS             , 'NPMASS                (143)'  !
      I = I + 1  ;     WRITE(L1A,160) NPROD              , 'NPROD                 (144)'  !
      I = I + 1  ;     WRITE(L1A,160) NPSHEAR            , 'NPSHEAR               (145)'  !
      I = I + 1  ;     WRITE(L1A,160) NPSHEL             , 'NPSHEL                (146)'  !
      I = I + 1  ;     WRITE(L1A,160) NPSOLID            , 'NPSOLID               (147)'  !
      I = I + 1  ;     WRITE(L1A,160) NPUSER1            , 'NPUSER1               (148)'  !
      I = I + 1  ;     WRITE(L1A,160) NPUSERIN           , 'NPUSERIN              (149)'  !
      I = I + 1  ;     WRITE(L1A,160) NRBAR              , 'NRBAR                 (150)'  !
      I = I + 1  ;     WRITE(L1A,160) NRBE1              , 'NRBE1                 (151)'  !
      I = I + 1  ;     WRITE(L1A,160) NRBE2              , 'NRBE2                 (152)'  !
      I = I + 1  ;     WRITE(L1A,160) NRFORCE            , 'NRFORCE               (153)'  !
      I = I + 1  ;     WRITE(L1A,160) NRIGEL             , 'NRIGEL                (154)'  !
      I = I + 1  ;     WRITE(L1A,160) NRECARD            , 'NRECARD               (155)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_ACCE     , 'NROWS_OTM_ACCE        (156)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_DISP     , 'NROWS_OTM_DISP        (157)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_MPCF     , 'NROWS_OTM_MPCF        (158)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_SPCF     , 'NROWS_OTM_SPCF        (159)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_ELFE     , 'NROWS_OTM_ELFE        (160)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_ELFN     , 'NROWS_OTM_ELFN        (161)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_STRE     , 'NROWS_OTM_STRE        (162)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_OTM_STRN     , 'NROWS_OTM_STRN        (163)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_ACCE     , 'NROWS_TXT_ACCE        (164)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_DISP     , 'NROWS_TXT_DISP        (165)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_MPCF     , 'NROWS_TXT_MPCF        (166)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_SPCF     , 'NROWS_TXT_SPCF        (167)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_ELFE     , 'NROWS_TXT_ELFE        (168)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_ELFN     , 'NROWS_TXT_ELFN        (169)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_STRE     , 'NROWS_TXT_STRE        (170)'  !
      I = I + 1  ;     WRITE(L1A,160) NROWS_TXT_STRN     , 'NROWS_TXT_STRN        (171)'  !
      I = I + 1  ;     WRITE(L1A,160) NRSPLINE           , 'NRSPLINE              (172)'  !
      I = I + 1  ;     WRITE(L1A,160) NSEQ               , 'NSEQ                  (173)'  !
      I = I + 1  ;     WRITE(L1A,160) NSETS              , 'NSETS                 (174)'  !
      I = I + 1  ;     WRITE(L1A,160) NSLOAD             , 'NSLOAD                (175)'  !
      I = I + 1  ;     WRITE(L1A,160) NSPC               , 'NSPC                  (176)'  !
      I = I + 1  ;     WRITE(L1A,160) NSPC1              , 'NSPC1                 (177)'  !
      I = I + 1  ;     WRITE(L1A,160) NSPCADD            , 'NSPCADD               (178)'  !
      I = I + 1  ;     WRITE(L1A,160) NSPOINT            , 'NSPOINT               (179)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_SPCSIDS        , 'NUM_SPCSIDS           (180)'  !
      I = I + 1  ;     WRITE(L1A,160) NSUB               , 'NSUB                  (181)'  !
      I = I + 1  ;     WRITE(L1A,160) NTCARD             , 'NTCARD                (182)'  !
      I = I + 1  ;     WRITE(L1A,160) NTDAT              , 'NTDAT                 (183)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_ALL          , 'NTERM_ALL             (184)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_CG_LTM       , 'NTERM_CG_LTM          (185)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_DLR          , 'NTERM_DLR             (186)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_GMN          , 'NTERM_GMN             (187)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_GOA          , 'NTERM_GOA             (188)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_HMN          , 'NTERM_HMN             (189)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_IF_LTM       , 'NTERM_IF_LTM          (190)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_IRR          , 'NTERM_IRR             (191)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KAA          , 'NTERM_KAA             (192)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KAAD         , 'NTERM_KAAD            (193)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KAO          , 'NTERM_KAO             (194)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KAOD         , 'NTERM_KAOD            (195)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFF          , 'NTERM_KFF             (196)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFFD         , 'NTERM_KFFD            (197)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFS          , 'NTERM_KFS             (198)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFSD         , 'NTERM_KFSD            (199)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFSe         , 'NTERM_KFSe            (200)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KFSDe        , 'NTERM_KFSDe           (201)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KGG          , 'NTERM_KGG             (202)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KGGD         , 'NTERM_KGGD            (203)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KLL          , 'NTERM_KLL             (204)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KLLD         , 'NTERM_KLLD            (205)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KLLDn        , 'NTERM_KLLDn           (206)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KLLs         , 'NTERM_KLLs            (207)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KLLDs        , 'NTERM_KLLDs           (208)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KMM          , 'NTERM_KMM             (209)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KMMD         , 'NTERM_KMMD            (210)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KMSM         , 'NTERM_KMSM            (211)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KMSMn        , 'NTERM_KMSMn           (212)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KMSMs        , 'NTERM_KMSMs           (213)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KNM          , 'NTERM_KNM             (214)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KNMD         , 'NTERM_KNMD            (215)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KNN          , 'NTERM_KNN             (216)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KNND         , 'NTERM_KNND            (217)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KOO          , 'NTERM_KOO             (218)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KOOD         , 'NTERM_KOOD            (219)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KOOs         , 'NTERM_KOOs            (220)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KOODs        , 'NTERM_KOODs           (221)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRL          , 'NTERM_KRL             (222)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRLD         , 'NTERM_KRLD            (223)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRR          , 'NTERM_KRR             (224)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRRD         , 'NTERM_KRRD            (225)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRRcb        , 'NTERM_KRRcb           (226)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRRcbs       , 'NTERM_KRRcbs          (227)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KRRcbn       , 'NTERM_KRRcbn          (228)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KXX          , 'NTERM_KXX             (229)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KSS          , 'NTERM_KSS             (230)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KSSD         , 'NTERM_KSSD            (231)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KSSe         , 'NTERM_KSSe            (232)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_KSSDe        , 'NTERM_KSSDe           (233)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_LMN          , 'NTERM_LMN             (234)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_LTM          , 'NTERM_LTM             (235)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MAA          , 'NTERM_MAA             (236)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MAO          , 'NTERM_MAO             (237)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MFF          , 'NTERM_MFF             (238)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MFS          , 'NTERM_MFS             (239)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MGG          , 'NTERM_MGG             (240)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MGGC         , 'NTERM_MGGC            (241)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MGGE         , 'NTERM_MGGE            (242)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MGGS         , 'NTERM_MGGS            (243)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MLL          , 'NTERM_MLL             (244)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MLLn         , 'NTERM_MLLn            (245)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MLLs         , 'NTERM_MLLs            (246)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MLR          , 'NTERM_MLR             (247)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MMM          , 'NTERM_MMM             (248)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MNM          , 'NTERM_MNM             (249)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MNN          , 'NTERM_MNN             (250)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MOO          , 'NTERM_MOO             (251)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MPF0         , 'NTERM_MPF0            (252)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MRL          , 'NTERM_MRL             (253)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MRN          , 'NTERM_MRN             (254)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MRR          , 'NTERM_MRR             (255)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MRRcb        , 'NTERM_MRRcb           (256)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MRRcbn       , 'NTERM_MRRcbn          (257)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MXX          , 'NTERM_MXX             (258)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MXXn         , 'NTERM_MXXn            (259)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_MSS          , 'NTERM_MSS             (260)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PA           , 'NTERM_PA              (261)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PF           , 'NTERM_PF              (262)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PFYS         , 'NTERM_PFYS            (263)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PG           , 'NTERM_PG              (264)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIXA        , 'NTERM_PHIXA           (265)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIXG        , 'NTERM_PHIXG           (266)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIZG        , 'NTERM_PHIZG           (267)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIZL        , 'NTERM_PHIZL           (268)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIZL1       , 'NTERM_PHIZL1          (269)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PHIZL2       , 'NTERM_PHIZL2          (270)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PL           , 'NTERM_PL              (271)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PM           , 'NTERM_PM              (272)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PN           , 'NTERM_PN              (273)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PO           , 'NTERM_PO              (274)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PR           , 'NTERM_PR              (275)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_PS           , 'NTERM_PS              (276)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_QM           , 'NTERM_QM              (277)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_QS           , 'NTERM_QS              (278)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_QSYS         , 'NTERM_QSYS            (279)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_RMG          , 'NTERM_RMG             (280)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_RMM          , 'NTERM_RMM             (281)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_RMN          , 'NTERM_RMN             (282)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_ULL          , 'NTERM_ULL             (283)'  !
      I = I + 1  ;     WRITE(L1A,160) NTERM_ULLI         , 'NTERM_ULLI            (284)'  !
      I = I + 1  ;     WRITE(L1A,160) NTSUB              , 'NTSUB                 (285)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_CB_DOFS        , 'NUM_CB_DOFS           (286)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_EIGENS         , 'NUM_EIGENS            (287)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_KLLD_DIAG_ZEROS, 'NUM_KLLD_DIAG_ZEROS   (288)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_MLL_DIAG_ZEROS , 'NUM_MLL_DIAG_ZEROS    (289)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_MPCSIDS        , 'NUM_MPCSIDS           (290)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_PARTVEC_RECORDS, 'NUM_PARTVEC_RECORDS   (291)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_PCHD_SPC1      , 'NUM_PCHD_SPC1         (292)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_SPC_RECORDS    , 'NUM_SPC_RECORDS       (293)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_SPC1_RECORDS   , 'NUM_SPC1_RECORDS      (294)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_SUPT_CARDS     , 'NUM_SUPT_CARDS        (295)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_USET_RECORDS   , 'NUM_USET_RECORDS      (296)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_USETSTR        , 'NUM_USETSTR           (297)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_USET           , 'NUM_USET              (298)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_USET_U1        , 'NUM_USET_U1           (299)'  !
      I = I + 1  ;     WRITE(L1A,160) NUM_USET_U2        , 'NUM_USET_U2           (300)'  !
      I = I + 1  ;     WRITE(L1A,160) NVEC               , 'NVEC                  (301)'  !
      I = I + 1  ;     WRITE(L1A,160) NVVEC              , 'NVVEC                 (302)'  !
      I = I + 1  ;     WRITE(L1A,160) SETLEN             , 'SETLEN                (304)'  !

! Write PARAM's (some need to stay the same in a restart)

      WRITE(L1A,191) ELFORCEN, 'PARAM ELFORCEN (needed for restart)'
      WRITE(L1A,191) HEXAXIS , 'PARAM HEXAXIS  (needed for restart)'
      WRITE(L1A,191) MATSPARS, 'PARAM MATSPARS (needed for restart)'
      WRITE(L1A,191) MIN4TRED, 'PARAM MIN4TRED (needed for restart)'
      WRITE(L1A,191) QUAD4TYP, 'PARAM QUAD4TYP (needed for restart)'
      WRITE(L1A,191) QUADAXIS, 'PARAM QUADAXIS (needed for restart)'
      WRITE(L1A,191) SPARSTOR, 'PARAM SPARSTOR (needed for restart)'
      WRITE(L1A,192) IORQ1M  , 'PARAM IORQ1M   (needed for restart)'
      WRITE(L1A,192) IORQ1S  , 'PARAM IORQ1S   (needed for restart)'
      WRITE(L1A,192) IORQ1B  , 'PARAM IORQ1B   (needed for restart)'
      WRITE(L1A,192) IORQ2B  , 'PARAM IORQ2B   (needed for restart)'
      WRITE(L1A,192) IORQ2T  , 'PARAM IORQ2T   (needed for restart)'
      WRITE(L1A,193) CBMIN3  , 'PARAM CBMIN3   (needed for restart)'
      WRITE(L1A,193) CBMIN4  , 'PARAM CBMIN4   (needed for restart)'

! Write COMM

      WRITE(L1A,103) (COMM(I),I=0,49)

      CALL FILE_CLOSE ( L1A, LINK1A, CLOSE_STAT )



      RETURN

! **********************************************************************************************************************************
  103 FORMAT(1X,50A1)

  110 FORMAT(1X,I14)

  120 FORMAT(1X,A)

  121 FORMAT(1X,1ES13.6)

  140 FORMAT(1X,I14)

  151 FORMAT(1X,I14,1X,A,1X,A,1X,A)

  160 FORMAT(I15,1X,A)

  191 FORMAT(6X,A9,1X,A)

  192 FORMAT(7X,I8,1X,A)

  193 FORMAT(1X,1ES14.6,1X,A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_L1A

   END MODULE FILE_LIFECYCLE
