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

   MODULE CASE_CONTROL_OUTPUTS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CC_ACCE, CC_DISP, CC_ELDA, CC_ELFO, CC_ENFO, CC_GPFO, CC_MPCF, CC_OLOA, CC_SPCF, CC_STRE, CC_STRN, COMPUTE_OUTPUT_TARGETS

   CONTAINS

      SUBROUTINE CC_ACCE ( CARD )

! Processes Case Control ACCE cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  ACCE_F06
      USE MODEL_STUF, ONLY            :  SC_ACCE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_ACCE'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
! CC_OUTPUTS processes all output type Case Control entries (they all have some common code so it is put there)

      CALL CC_OUTPUTS ( CARD, 'ACCE', SETID )

! Check to see if BOTH, ENGR or NODE were in the ELFO request

      CALL COMPUTE_OUTPUT_TARGETS(ACCE_F06)

! Set CASE CONTROL output request variable to SETID

      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_ACCE(I) = SETID
         ENDDO
      ELSE
         SC_ACCE(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_ACCE


      SUBROUTINE CC_DISP ( CARD )

      ! Processes Case Control cards for requests for displacement outputs
      ! - DISP: displacement
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  DISP_F06
      USE MODEL_STUF, ONLY            :  SC_DISP

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_DISP'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'DISP', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(DISP_F06)

      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_DISP(I) = SETID
         ENDDO
      ELSE
         SC_DISP(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_DISP


      SUBROUTINE CC_ELDA ( CARD )

! Processes Case Control ELDATA cards

! NOTE: The coding assumes that ELDATA(i,BOTH) is only valid for i >= IOUTMIN_FIJ (which is 1) and i <= IOUTMAX_FIJ (which is 5)
! ----

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, WARN_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  CCELDT

      USE BDF_SET_SYNTAX, ONLY        :  GET_ANSID, STOKEN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'CC_ELDA'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: ERRTOK            ! Character string that holds part of an error message from subr STOKEN
      CHARACTER( 3*BYTE)              :: EXCEPT            ! An input/output for subr STOKEN
      CHARACTER( 1*BYTE)              :: FIJFIL_WARN = 'N' ! Set to 'Y' if warning message written for FIJFIL option
      CHARACTER( 3*BYTE)              :: THRU              ! An input/output for subr STOKEN
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! Char string output from subr STOKEN, called herein
      CHARACTER(LEN=LEN(CARD))        :: TOKSTR            ! Character string to tokenize
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! Type of the char TOKEN's output from subr STOKEN, called herein

      INTEGER(LONG)                   :: ICOL1       = 0   ! Location, in CARD, where "(" begins
      INTEGER(LONG)                   :: ICOL2       = 0   ! Location, in CARD, where ")" begins
      INTEGER(LONG)                   :: IERROR            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: IOCHK       = 0   ! IOSTAT error number when reading data from a file
      INTEGER(LONG)                   :: IOFF        = 8   ! An index offset into array CCELDT
      INTEGER(LONG)                   :: IOUT              ! Integer number, in ELDATA request, indicating type of output request
      INTEGER(LONG), PARAMETER        :: IOUTMIN_BUG = 0   ! Min val of IOUT (=1,2,3,4,5,6,7,8,9 are the ELDATA print options)
      INTEGER(LONG), PARAMETER        :: IOUTMAX_BUG = 9   ! Max val of IOUT (=1,2,3,4,5,6,7,8,9 are the ELDATA print options)
      INTEGER(LONG), PARAMETER        :: IOUTMIN_FIJ = 2   ! Min val of IOUT (=1,2,3,4,5,6,7,8,9 are the ELDATA file  options)
      INTEGER(LONG), PARAMETER        :: IOUTMAX_FIJ = 6   ! Max val of IOUT (=1,2,3,4,5,6,7,8,9 are the ELDATA file  options)
      INTEGER(LONG)                   :: NTOKEN            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: SETID       = 0   ! Set ID on this Case Control card
      INTEGER(LONG)                   :: STRNG_LEN   = 0   ! Length of character string between "()" in the ELDATA card
      INTEGER(LONG)                   :: TOKEN_BEG   = 0   ! An input to subr STOKEN, called herein


      INTRINSIC INDEX



! **********************************************************************************************************************************
! Process ELDATA cards.

! Processes element debug output or disk file output. Card format is:

!           ELDATA(i,PRINT or FIJFIL or BOTH) = NONE or ALL or SETID

! where i is 0,1,2,3,4 or 5 and either PRINT or FIJFIL or BOTH must be selected. i determines kind of elem output.

! First get i,PRINT, FIJFIL between (). CCELDT(j) array will be set equal to:

!       1)  -1   if "ALL" is requested,
!       2)   0   if "NONE" is requested, or
!       3) SETID if a set ID is requested

! CCELDT( 0) is CC requests for print to BUGFIL of elem geometric data
! CCELDT( 1) is CC requests for print to BUGFIL of elem property and material info
! CCELDT( 2) is CC requests for print to BUGFIL of elem thermal and pressure matrices    : PTE, PPE
! CCELDT( 3) is CC requests for print to BUGFIL of elem mass matrix                      : ME
! CCELDT( 4) is CC requests for print to BUGFIL of elem stiffness matrix                 : KE
! CCELDT( 5) is CC requests for print to BUGFIL of elem stress & strain recovery matrices: SEi, STEi, BEi
! CCELDT( 6) is CC requests for print to BUGFIL of elem displacement and load matrices   : UEL, PEL (all subcases)
! CCELDT( 7) is CC requests for print to BUGFIL of elem shape fcns and Jacobian matrices
! CCELDT( 8) is CC requests for print to BUGFIL of elem strain-displacement matrices
! CCELDT( 9) is CC requests for print to BUGFIL of elem checks on strain-displ matrices for RB motion & constant strain
! CCELDT(10) is CC requests for write to F21FIL unformatted file of element              : PTE, PPE
! CCELDT(11) is CC requests for write to F22FIL unformatted file of element              : ME
! CCELDT(12) is CC requests for write to F23FIL unformatted file of element              : KE
! CCELDT(13) is CC requests for write to F24FIL unformatted file of element              : SEi, STEi, BEi
! CCELDT(14) is CC requests for write to F25FIL unformatted file of element              : UEL, PEL (all subcases)
! CCELDT(15) is CC requests currently not used

! CCELDT will be processed in subroutine SCPRO to fill in output array ELDAT(k) with elements whose output was
! requested by setting bit j-1 in array ELDAT to 1 if CCELDT(j) is nonzero

! Find out if "NONE", "ALL" or SETID

      CALL GET_ANSID ( CARD, SETID )

! Get data in between ()

      ICOL1  = INDEX(CARD(1:),'(')
      ICOL2  = INDEX(CARD(1:),')')
      STRNG_LEN = ICOL2 - ICOL1 - 1
      IF ((ICOL1 == 0) .OR. (ICOL2 == 0) .OR. (ICOL2 < ICOL1)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1201)
         WRITE(F06,1201)
         RETURN
      ELSE IF (STRNG_LEN < 1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1202)
         WRITE(F06,1202)
         RETURN
      ELSE

! Read 1st token after opening paren: "(". Need to make sure that the token is:

!  1) An integer                               : TOKTYP(1) /= 'INTEGER ',
!  2) Integer has <= 8 digits                  : IERROR /= 0, and
!  3) Integer starts before closing paren, ")" :

         TOKSTR(1:STRNG_LEN) = CARD(ICOL1+1:ICOL2-1)
         TOKEN_BEG = 1
         THRU   = 'OFF'
         EXCEPT = 'OFF'
         CALL STOKEN ( SUBR_NAME, TOKSTR, TOKEN_BEG, STRNG_LEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
         IF ((TOKTYP(1) /= 'INTEGER ') .OR. (IERROR /= 0)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1281)
            WRITE(F06,1281)
            RETURN
         ELSE
            READ(TOKEN(1),'(I8)',IOSTAT=IOCHK) IOUT
            IF (IOCHK > 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1264)
               WRITE(F06,1264)
            ENDIF
            IF (TOKEN_BEG <= STRNG_LEN) THEN
               THRU   = 'OFF'
               EXCEPT = 'OFF'

               CALL STOKEN ( SUBR_NAME, TOKSTR, TOKEN_BEG, STRNG_LEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
               IF (IERROR == 1) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1283) 'ELDATA', ERRTOK
                  WRITE(F06,1283) 'ELDATA', ERRTOK
               ELSE
                  IF      (TOKTYP(1) == 'PRINT   ') THEN
                     IF ((IOUT < IOUTMIN_BUG) .OR. (IOUT > IOUTMAX_BUG)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1282) IOUTMIN_BUG, IOUTMAX_BUG, 'PRINT'
                        WRITE(F06,1282) IOUTMIN_BUG, IOUTMAX_BUG, 'PRINT'
                        RETURN
                     ENDIF
                     CCELDT(IOUT) = SETID
                  ELSE IF (TOKTYP(1) == 'FIJFIL  ') THEN
                     IF ((IOUT < IOUTMIN_FIJ) .OR. (IOUT > IOUTMAX_FIJ)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1282) IOUTMIN_FIJ, IOUTMAX_FIJ, 'FIJFIL'
                        WRITE(F06,1282) IOUTMIN_FIJ, IOUTMAX_FIJ, 'FIJFIL'
                        RETURN
                     ENDIF
                     CCELDT(IOUT+IOFF) = SETID
                  ELSE IF (TOKTYP(1) == 'BOTH    ') THEN
                     IF ((IOUT < IOUTMIN_FIJ) .OR. (IOUT > IOUTMAX_FIJ)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1282) IOUTMIN_FIJ, IOUTMAX_FIJ, 'BOTH'
                        WRITE(F06,1282) IOUTMIN_FIJ, IOUTMAX_FIJ, 'BOTH'
                        RETURN
                     ENDIF
                     CCELDT(IOUT)      = SETID
                     CCELDT(IOUT+IOFF) = SETID
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1284) TOKEN(1)
                     WRITE(F06,1284) TOKEN(1)
                  ENDIF
               ENDIF
            ELSE                                        ! TOKEN_BEG > STRNG_LEN, so quit
               CCELDT(IOUT) = SETID                     ! Default for ELDATA output is PRINT
            ENDIF
            IF (TOKEN_BEG <= STRNG_LEN) THEN            ! Error: there shouldn't be more than the integer, PRINT, FIJFIL or BOTH
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1285) CARD,TOKEN(1)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1285) CARD,TOKEN(1)
               ENDIF
               IF (FIJFIL_WARN == 'N') THEN
                  WARN_ERR = WARN_ERR + 1
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! As of 03/07/2020 there is some error in the calculation of the check on strain displ matrices for RB and constant strain modes
! of displacement. Therefore, this check is temporarily suspended

      if (cceldt(9) /= 0) then
         write(f06,12344)
         write(f06,12345)
         write(f06,12344)
         write(f06,*)
      endif



      RETURN

! **********************************************************************************************************************************
 1264 FORMAT(' *ERROR  1264: ERROR READING ELDATA CASE CONTROL ENTRY FOR THE INTEGER DESCRIBING THE TYPE OF OUTPUT. ENTRY IGNORED')

 1201 FORMAT(' *ERROR  1201: COULD NOT FIND MATCHING PARENTHESES () FOLLOWING ELDATA CASE CONTROL COMMAND.'                        &
                    ,/,14X,' THE PARENTHESES ARE FOR THE PURPOSE OF CONTAINING THE SPECIFIC ELDATA OUTPUT REQUEST')

 1202 FORMAT(' *ERROR  1202: THERE IS NO DATA BETWEEN THE MATCHING PARENTHESES FOLLOWING THE ELDATA COMMAND')

 1281 FORMAT(' *ERROR  1281: THE FIRST ENTRY FOLLOWING THE OPENING PARENTHESIS OF THE ELDATA CASE CONTROL ENTRY:'                  &
                    ,/,14X,' (1) MUST BE AN INTEGER,'                                                                              &
                    ,/,14X,' (2) MUST HAVE LESS THAN 8 DIGITS,'                                                                    &
                    ,/,14X,' (3) AND MUST BEGIN BEFORE CLOSING PARENTHESIS')

 1282 FORMAT(' *ERROR  1282: THE INTEGER DESCRIBING THE TYPE OF OUTPUT FOR CASE CONTROL ELDATA ENTRY MUST:'                        &
                    ,/,14X,' HAVE A VALUE >= ',I2,' AND <= ',I2,' FOR OPTION ',A)

 1283 FORMAT(' *ERROR  1283: ENTRIES ON ',A,' CASE CONTROL ENTRY MUST BE <= 8 CHARACTERS. THE FOLLOWING ENTRY EXCEEDS THIS:'       &
                    ,/,15X,A,/)

 1284 FORMAT(' *ERROR  1284: DATA BETWEEN PARENS FOLLOWING INTEGER ON ELFORCE/FORCE CASE CONTROL ENTRY MUST BE "PRINT", "FIJFIL",',&
                           ' OR "BOTH".'                                                                                           &
                    ,/,14X,' ENTRY = ',A8,' NOT RECOGNIZED')

 1285 FORMAT(' *WARNING    : ERROR READING ELDATA CASE CONTROL ENTRY:',/,15X,A,/,15X,'DATA BETWEEN PARENS FOLOWING INTEGER SHOULD',&
                           ' ONLY BE ONE OF: "PRINT", "FIJFIL" OR "BOTH". FIRST ENTRY, ',A,', WILL BE USED')

12344 format('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

12345 format('NOTE: as of 03/07/2020 the check on strain-displacement matrices using Case Control ELDATA(9) is suspended',/,6x,    &
                   'until an error in the calculation is fixed. This can be overridden with DEBUG(202) > 0')

! **********************************************************************************************************************************

      END SUBROUTINE CC_ELDA


      SUBROUTINE CC_ELFO ( CARD )

      ! Processes Case Control ELFO (elforce) entries
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  err
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SC_ELFE, SC_ELFN
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  FORC_F06

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_ELFO'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      CHARACTER( 1*BYTE)              :: FOUND_BOTH        ! CC_CMD_DESCRIBERS has request for "BOTH"
      CHARACTER( 1*BYTE)              :: FOUND_ENGR        ! CC_CMD_DESCRIBERS has request for "ENGR"
      CHARACTER( 1*BYTE)              :: FOUND_NODE        ! CC_CMD_DESCRIBERS has request for "NODE"

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'ELFO', SETID )


      ! Check to see if BOTH, ENGR or NODE were in the ELFO request
      FOUND_BOTH = 'N'
      FOUND_ENGR = 'N'
      FOUND_NODE = 'N'
      DO I=1,NCCCD
         IF (CC_CMD_DESCRIBERS(I)(1:4) == 'BOTH') FOUND_BOTH = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:4) == 'ENGR') FOUND_ENGR = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:4) == 'NODE') FOUND_NODE = 'Y'
      ENDDO

      CALL COMPUTE_OUTPUT_TARGETS(FORC_F06)


      ! Set CASE CONTROL output request variable to SETID
      IF      (FOUND_BOTH == 'Y') THEN

         IF (NSUB == 0) THEN
            DO I = 1,LSUB
               SC_ELFE(I) = SETID
               SC_ELFN(I) = SETID
            ENDDO
         ELSE
            SC_ELFE(NSUB) = SETID
            SC_ELFN(NSUB) = SETID
         ENDIF

      ELSE IF (FOUND_NODE == 'Y') THEN

         IF (NSUB == 0) THEN
            DO I = 1,LSUB
               SC_ELFN(I) = SETID
            ENDDO
         ELSE
            SC_ELFN(NSUB) = SETID
         ENDIF

      ELSE
         ! Default is ENGR
         IF (NSUB == 0) THEN
            DO I = 1,LSUB
               SC_ELFE(I) = SETID
            ENDDO
         ELSE
            SC_ELFE(NSUB) = SETID
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_ELFO


      SUBROUTINE CC_ENFO ( CARD )

! Processes Case Control ENFO (ENFORCED) entries

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ENFFIL, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT
      USE TEXT_FIELD_UTILS, ONLY      :  GET_CHAR_STRING_END

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_ENFO'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: IEND              ! Col where end of data is on CARD1
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error




! **********************************************************************************************************************************
! Process ENFORCED card

      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,8862)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8862)
         ENDIF
      ENDIF

      CALL GET_CHAR_STRING_END ( CARD1, IEND )
      ENFFIL = CARD1(1:IEND)



      RETURN

! **********************************************************************************************************************************
 8862 FORMAT(' *WARNING    : MISSING EQUAL (=) SIGN ON CASE CONTROL CARD: (TITLE, SUBTITLE, LABEL OR ENFORCED).')

! **********************************************************************************************************************************

      END SUBROUTINE CC_ENFO


      SUBROUTINE CC_GPFO ( CARD )

      ! Processes Case Control GPFO cards that define grid point force balance output requests
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  GPFO_F06
      USE MODEL_STUF, ONLY            :  SC_GPFO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_GPFO'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'GPFO', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(GPFO_F06)

      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_GPFO(I) = SETID
         ENDDO
      ELSE
         SC_GPFO(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_GPFO


      SUBROUTINE CC_MPCF ( CARD )

      ! Processes Case Control cards for requests for MPC force output requests
      ! - MPCF: MPC force
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  MPCF_F06
      USE MODEL_STUF, ONLY            :  SC_MPCF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_MPCF'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'MPCF', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(MPCF_F06)


      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_MPCF(I) = SETID
         ENDDO
      ELSE
         SC_MPCF(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_MPCF


      SUBROUTINE CC_OLOA ( CARD )

      ! Processes Case Control cards for requests for applied load requests
      ! - OLOA: applied load
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  OLOA_F06
      USE MODEL_STUF, ONLY            :  SC_OLOA

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_OLOA'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'OLOA', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(OLOA_F06)


      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_OLOA(I) = SETID
         ENDDO
      ELSE
         SC_OLOA(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_OLOA


      SUBROUTINE CC_SPCF ( CARD )

! Processes Case Control SPCF cards for SPC force output requests
      ! - SPCF: SPC force
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  SPCF_F06
      USE MODEL_STUF, ONLY            :  SC_SPCF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SPCF'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'SPCF', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(SPCF_F06)


      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_SPCF(I) = SETID
         ENDDO
      ELSE
         SC_SPCF(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_SPCF


      SUBROUTINE CC_STRE ( CARD )

      ! Processes Case Control STRE cards for element stress output requests
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_F06
      USE MODEL_STUF, ONLY            :  SC_STRE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_STRE'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries (they all
      ! have some common code so it is put there)

      CALL CC_OUTPUTS ( CARD, 'STRE', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(STRE_F06)

      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_STRE(I) = SETID
         ENDDO
      ELSE
         SC_STRE(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_STRE


      SUBROUTINE CC_STRN ( CARD )

      ! Processes Case Control STRN cards for element strain output requests
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_F06
      USE MODEL_STUF, ONLY            :  SC_STRN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_STRN'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)

      CALL CC_OUTPUTS ( CARD, 'STRN', SETID )

      CALL COMPUTE_OUTPUT_TARGETS(STRN_F06)

      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_STRN(I) = SETID
         ENDDO
      ELSE
         SC_STRN(NSUB) = SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_STRN


      SUBROUTINE CC_OUTPUTS ( CARD, WHAT, SETID )

! Process the character string in parens for the following Case Control output request entries (e.g. SORT1, PRINT, etc)

!       ACCE()
!       DISP()
!       ELFO()
!       GPFO()
!       MPCF()
!       OLOA()
!       SPCF()
!       STRE()
!       STRN()


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NCCCD, NSUB
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_SET_SYNTAX, ONLY        :  GET_ANSID
      USE TEXT_FIELD_UTILS, ONLY      :  PARSE_CHAR_STRING

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_OUTPUTS'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT              ! Which CC type output to process (e.g., DISP, SPCF, etc)
      CHARACTER(LEN=LEN(CARD))        :: CHAR_STRING       ! Character string between parens () if it exists

      INTEGER(LONG), INTENT(OUT)      :: SETID             ! Set ID on this Case Control card
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICOL1       = 0   ! Location, in CARD, where "(" begins
      INTEGER(LONG)                   :: ICOL2       = 0   ! Location, in CARD, where ")" begins
      INTEGER(LONG)                   :: IERR        = 0   ! Error designator from subr PARSE_CSV_STRING
      INTEGER(LONG)                   :: NUM_WORDS   = 0   ! Number of words in the string between parens (), if present
      INTEGER(LONG)                   :: STRING_LEN  = 0   ! Length of character string between "()" in the ELDATA card




! **********************************************************************************************************************************
! Initialize

      IERR = 0

      CHAR_STRING(1:) = ' '

      DO I=1,NCCCD
         CC_CMD_DESCRIBERS(I)(1:) = ' '
      ENDDO

! Find out if "NONE", "ALL" or SETID

      CALL GET_ANSID ( CARD, SETID )

! Find out if there is data enclosed between parens (). If so, STRING_LEN will be > 0.

      ICOL1 = INDEX(CARD(1:),'(')
      ICOL2 = INDEX(CARD(1:),')')
      STRING_LEN = ICOL2 - ICOL1 - 1

! If there is a string of data in between the parens then parse it to get all of the words and put them into array CC_CMD_DESCRIBERS

      IF (STRING_LEN > 0) THEN                             ! There is a pair () of parens, so we will have to parse data between ()
         CHAR_STRING(1:STRING_LEN) = CARD(ICOL1+1:ICOL2-1)
         CALL PARSE_CHAR_STRING ( CHAR_STRING, STRING_LEN, NCCCD, LEN(CC_CMD_DESCRIBERS), NUM_WORDS, CC_CMD_DESCRIBERS, IERR )
      ELSE
         RETURN
      ENDIF

! Check that command descriptors (e.g. "SORT1") are allowable

      IF (IERR == 0) THEN
         IF (NUM_WORDS <= NCCCD) THEN
            CALL CHK_CC_CMD_DESCRIBERS ( WHAT, NUM_WORDS )
         ELSE
            CALL CHK_CC_CMD_DESCRIBERS ( WHAT, NCCCD     )
         ENDIF
      ELSE
         RETURN
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_OUTPUTS


      SUBROUTINE CHK_CC_CMD_DESCRIBERS ( WHAT, NUM_WORDS )

     ! Checks Case Control output requests to make sure the descriptors in parens
     ! (e.g. SORT1, PRINT, etc) are valid.
     ! Write warning messages if a descriptor is not valid for MYSTRAN

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG
      USE IOUNT1, ONLY                 :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                 :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, ECHO, FATAL_ERR, WARN_ERR, NSUB
      USE TIMDAT, ONLY                 :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY   :  STRN_LOC, STRN_OPT, STRN_CUR, STRE_LOC, STRE_OPT, FORC_LOC
      USE PARAMS, ONLY                 :  SUPWARN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER         :: NUM_POSS_CCD = 31 ! Number of possible CC command describers (incl all MSC ones as well)
      INTEGER(LONG), PARAMETER         :: NUM_OUT_TYP  =  9 ! Number of OUTPUT_TYPE's

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)) :: SUBR_NAME = 'CHK_CC_CMD_DESCRIBERS'
      CHARACTER(LEN=*), INTENT(IN)     :: WHAT              ! What Case Control output is this call for (e.g. 'DISP')

                                                            ! OUTPUT_TYPE is 'ACCE', 'DISP', etc
      CHARACTER( 4*BYTE)               :: OUTPUT_TYPE(NUM_OUT_TYP)

                                                            ! List of all of the possible values for a MSC NASTRAN CC_CMD_DESCRIBERS
      CHARACTER(LEN(CC_CMD_DESCRIBERS)):: ALLOW_CC_CMD_DESCR(NUM_POSS_CCD,NUM_OUT_TYP)

      CHARACTER( 1*BYTE)               :: FOUND             ! ='Y' if we found something we are looking for in an IF test

      INTEGER(LONG), INTENT(IN)        :: NUM_WORDS         ! Number of words we need to check in CC_CMD_DESCRIBERS
      INTEGER(LONG)                    :: I,J               ! DO loop indices
      INTEGER(LONG)                    :: JCOL              ! Designator of a column in an array

      LOGICAL                          :: IS_PLOT, IS_PRINT, IS_PUNCH


! **********************************************************************************************************************************
      IF      (WHAT == 'ACCE') THEN;   OUTPUT_TYPE( 1) = 'ACCE';   JCOL =  1;
      ELSE IF (WHAT == 'DISP') THEN;   OUTPUT_TYPE( 2) = 'DISP';   JCOL =  2;
      ELSE IF (WHAT == 'ELFO') THEN;   OUTPUT_TYPE( 3) = 'ELFO';   JCOL =  3;
      ELSE IF (WHAT == 'GPFO') THEN;   OUTPUT_TYPE( 4) = 'GPFO';   JCOL =  4;
      ELSE IF (WHAT == 'MPCF') THEN;   OUTPUT_TYPE( 5) = 'MPCF';   JCOL =  5;
      ELSE IF (WHAT == 'OLOA') THEN;   OUTPUT_TYPE( 6) = 'OLOA';   JCOL =  6;
      ELSE IF (WHAT == 'SPCF') THEN;   OUTPUT_TYPE( 7) = 'SPCF';   JCOL =  7;
      ELSE IF (WHAT == 'STRE') THEN;   OUTPUT_TYPE( 8) = 'STRE';   JCOL =  8;
      ELSE IF (WHAT == 'STRN') THEN;   OUTPUT_TYPE( 9) = 'STRN';   JCOL =  9;
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1204) SUBR_NAME, WHAT
         WRITE(F06,1204) SUBR_NAME, WHAT
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IS_PLOT = ((WHAT == 'ACCE') .OR. (WHAT == 'DISP') .OR. (WHAT == 'ELFO') .OR. (WHAT == 'GPFO')  .OR.  &
                 (WHAT == 'MPCF') .OR. (WHAT == 'OLOA') .OR. (WHAT == 'SPCF') .OR.  &
                 (WHAT == 'STRE') .OR. (WHAT == 'STRN'))

      ! same as plot
      IS_PRINT = IS_PLOT

      ! remove GPFO from PUNCH
      IS_PUNCH = ((WHAT == 'ACCE') .OR. (WHAT == 'DISP') .OR. (WHAT == 'ELFO') .OR.  &
                  (WHAT == 'MPCF') .OR. (WHAT == 'OLOA') .OR. (WHAT == 'SPCF') .OR.  &
                  (WHAT == 'STRE') .OR. (WHAT == 'STRN'))

      ! Set all of the allowable values that can be in ALLOW_CC_CMD_DESCR.
      ! These are all of the values from MSC NASTRAN. Not all are implemented in MYSTRAN.

!     =================ACCE=================   =================DISP=================   =================ELFO=================
      ALLOW_CC_CMD_DESCR( 1, 1) = 'SORT1   ' ; ALLOW_CC_CMD_DESCR( 1, 2) = 'SORT1   ' ; ALLOW_CC_CMD_DESCR( 1, 3) = 'SORT1   '
      ALLOW_CC_CMD_DESCR( 2, 1) = 'SORT2   ' ; ALLOW_CC_CMD_DESCR( 2, 2) = 'SORT2   ' ; ALLOW_CC_CMD_DESCR( 2, 3) = 'SORT2   '
      ALLOW_CC_CMD_DESCR( 3, 1) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 2) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 3) = 'PRINT   '
      ALLOW_CC_CMD_DESCR( 4, 1) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 2) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 3) = 'PUNCH   '
      ALLOW_CC_CMD_DESCR( 5, 1) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 2) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 3) = 'PLOT    '
      ALLOW_CC_CMD_DESCR( 6, 1) = 'REAL    ' ; ALLOW_CC_CMD_DESCR( 6, 2) = 'REAL    ' ; ALLOW_CC_CMD_DESCR( 6, 3) = 'REAL    '
      ALLOW_CC_CMD_DESCR( 7, 1) = 'IMAG    ' ; ALLOW_CC_CMD_DESCR( 7, 2) = 'IMAG    ' ; ALLOW_CC_CMD_DESCR( 7, 3) = 'IMAG    '
      ALLOW_CC_CMD_DESCR( 8, 1) = 'MAG     ' ; ALLOW_CC_CMD_DESCR( 8, 2) = 'MAG     ' ; ALLOW_CC_CMD_DESCR( 8, 3) = 'MAG     '
      ALLOW_CC_CMD_DESCR( 9, 1) = 'PHASE   ' ; ALLOW_CC_CMD_DESCR( 9, 2) = 'PHASE   ' ; ALLOW_CC_CMD_DESCR( 9, 3) = 'PHASE   '
      ALLOW_CC_CMD_DESCR(10, 1) = '        ' ; ALLOW_CC_CMD_DESCR(10, 2) = '        ' ; ALLOW_CC_CMD_DESCR(10, 3) = '        '
      ALLOW_CC_CMD_DESCR(11, 1) = '        ' ; ALLOW_CC_CMD_DESCR(11, 2) = '        ' ; ALLOW_CC_CMD_DESCR(11, 3) = '        '
      ALLOW_CC_CMD_DESCR(12, 1) = '        ' ; ALLOW_CC_CMD_DESCR(12, 2) = '        ' ; ALLOW_CC_CMD_DESCR(12, 3) = '        '
      ALLOW_CC_CMD_DESCR(13, 1) = '        ' ; ALLOW_CC_CMD_DESCR(13, 2) = '        ' ; ALLOW_CC_CMD_DESCR(13, 3) = '        '
      ALLOW_CC_CMD_DESCR(14, 1) = '        ' ; ALLOW_CC_CMD_DESCR(14, 2) = '        ' ; ALLOW_CC_CMD_DESCR(14, 3) = '        '
      ALLOW_CC_CMD_DESCR(15, 1) = '        ' ; ALLOW_CC_CMD_DESCR(15, 2) = '        ' ; ALLOW_CC_CMD_DESCR(15, 3) = '        '
      ALLOW_CC_CMD_DESCR(16, 1) = '        ' ; ALLOW_CC_CMD_DESCR(16, 2) = '        ' ; ALLOW_CC_CMD_DESCR(16, 3) = '        '
      ALLOW_CC_CMD_DESCR(17, 1) = '        ' ; ALLOW_CC_CMD_DESCR(17, 2) = '        ' ; ALLOW_CC_CMD_DESCR(17, 3) = '        '
      ALLOW_CC_CMD_DESCR(18, 1) = '        ' ; ALLOW_CC_CMD_DESCR(18, 2) = '        ' ; ALLOW_CC_CMD_DESCR(18, 3) = '        '
      ALLOW_CC_CMD_DESCR(19, 1) = '        ' ; ALLOW_CC_CMD_DESCR(19, 2) = '        ' ; ALLOW_CC_CMD_DESCR(19, 3) = '        '
      ALLOW_CC_CMD_DESCR(20, 1) = 'PSDF    ' ; ALLOW_CC_CMD_DESCR(20, 2) = 'PSDF    ' ; ALLOW_CC_CMD_DESCR(20, 3) = 'PSDF    '
      ALLOW_CC_CMD_DESCR(21, 1) = 'ATOC    ' ; ALLOW_CC_CMD_DESCR(21, 2) = 'ATOC    ' ; ALLOW_CC_CMD_DESCR(21, 3) = 'ATOC    '
      ALLOW_CC_CMD_DESCR(22, 1) = 'CRMS    ' ; ALLOW_CC_CMD_DESCR(22, 2) = 'CRMS    ' ; ALLOW_CC_CMD_DESCR(22, 3) = 'CRMS    '
      ALLOW_CC_CMD_DESCR(23, 1) = 'RALL    ' ; ALLOW_CC_CMD_DESCR(23, 2) = 'RALL    ' ; ALLOW_CC_CMD_DESCR(23, 3) = 'RALL    '
      ALLOW_CC_CMD_DESCR(24, 1) = 'RPRINT  ' ; ALLOW_CC_CMD_DESCR(24, 2) = 'RPRINT  ' ; ALLOW_CC_CMD_DESCR(24, 3) = 'RPRINT  '
      ALLOW_CC_CMD_DESCR(25, 1) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(25, 2) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(25, 3) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(26, 1) = 'NOPRINT ' ; ALLOW_CC_CMD_DESCR(26, 2) = 'NOPRINT ' ; ALLOW_CC_CMD_DESCR(26, 3) = 'NOPRINT '
      ALLOW_CC_CMD_DESCR(27, 1) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(27, 2) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(27, 3) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(28, 1) = 'CID     ' ; ALLOW_CC_CMD_DESCR(28, 2) = 'CID     ' ; ALLOW_CC_CMD_DESCR(28, 3) = 'CID     '
      ALLOW_CC_CMD_DESCR(29, 1) = '        ' ; ALLOW_CC_CMD_DESCR(29, 2) = '        ' ; ALLOW_CC_CMD_DESCR(29, 3) = 'BOTH    '
      ALLOW_CC_CMD_DESCR(30, 1) = '        ' ; ALLOW_CC_CMD_DESCR(30, 2) = '        ' ; ALLOW_CC_CMD_DESCR(30, 3) = 'ENGR    '
      ALLOW_CC_CMD_DESCR(31, 1) = '        ' ; ALLOW_CC_CMD_DESCR(31, 2) = '        ' ; ALLOW_CC_CMD_DESCR(31, 3) = 'NODE    '

!     =================GPFO=================   =================MPCF=================   =================OLOA=================
      ALLOW_CC_CMD_DESCR( 1, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 1, 5) = 'SORT1   ' ; ALLOW_CC_CMD_DESCR( 1, 6) = 'SORT1   '
      ALLOW_CC_CMD_DESCR( 2, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 2, 5) = 'SORT2   ' ; ALLOW_CC_CMD_DESCR( 2, 6) = 'SORT2   '
      ALLOW_CC_CMD_DESCR( 3, 4) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 5) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 6) = 'PRINT   '
      ALLOW_CC_CMD_DESCR( 4, 4) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 5) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 6) = 'PUNCH   '
      ALLOW_CC_CMD_DESCR( 5, 4) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 5) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 6) = 'PLOT    '
      ALLOW_CC_CMD_DESCR( 6, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 6, 5) = 'REAL    ' ; ALLOW_CC_CMD_DESCR( 6, 6) = 'REAL    '
      ALLOW_CC_CMD_DESCR( 7, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 7, 5) = 'IMAG    ' ; ALLOW_CC_CMD_DESCR( 7, 6) = 'IMAG    '
      ALLOW_CC_CMD_DESCR( 8, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 8, 5) = 'MAG     ' ; ALLOW_CC_CMD_DESCR( 8, 6) = 'MAG     '
      ALLOW_CC_CMD_DESCR( 9, 4) = '        ' ; ALLOW_CC_CMD_DESCR( 9, 5) = 'PHASE   ' ; ALLOW_CC_CMD_DESCR( 9, 6) = 'PHASE   '
      ALLOW_CC_CMD_DESCR(10, 4) = '        ' ; ALLOW_CC_CMD_DESCR(10, 5) = '        ' ; ALLOW_CC_CMD_DESCR(10, 6) = '        '
      ALLOW_CC_CMD_DESCR(11, 4) = '        ' ; ALLOW_CC_CMD_DESCR(11, 5) = '        ' ; ALLOW_CC_CMD_DESCR(11, 6) = '        '
      ALLOW_CC_CMD_DESCR(12, 4) = '        ' ; ALLOW_CC_CMD_DESCR(12, 5) = '        ' ; ALLOW_CC_CMD_DESCR(12, 6) = '        '
      ALLOW_CC_CMD_DESCR(13, 4) = '        ' ; ALLOW_CC_CMD_DESCR(13, 5) = '        ' ; ALLOW_CC_CMD_DESCR(13, 6) = '        '
      ALLOW_CC_CMD_DESCR(14, 4) = '        ' ; ALLOW_CC_CMD_DESCR(14, 5) = '        ' ; ALLOW_CC_CMD_DESCR(14, 6) = '        '
      ALLOW_CC_CMD_DESCR(15, 4) = '        ' ; ALLOW_CC_CMD_DESCR(15, 5) = '        ' ; ALLOW_CC_CMD_DESCR(15, 6) = '        '
      ALLOW_CC_CMD_DESCR(16, 4) = '        ' ; ALLOW_CC_CMD_DESCR(16, 5) = '        ' ; ALLOW_CC_CMD_DESCR(16, 6) = '        '
      ALLOW_CC_CMD_DESCR(17, 4) = '        ' ; ALLOW_CC_CMD_DESCR(17, 5) = '        ' ; ALLOW_CC_CMD_DESCR(17, 6) = '        '
      ALLOW_CC_CMD_DESCR(18, 4) = '        ' ; ALLOW_CC_CMD_DESCR(18, 5) = '        ' ; ALLOW_CC_CMD_DESCR(18, 6) = '        '
      ALLOW_CC_CMD_DESCR(19, 4) = '        ' ; ALLOW_CC_CMD_DESCR(19, 5) = '        ' ; ALLOW_CC_CMD_DESCR(19, 6) = '        '
      ALLOW_CC_CMD_DESCR(20, 4) = '        ' ; ALLOW_CC_CMD_DESCR(20, 5) = 'PSDF    ' ; ALLOW_CC_CMD_DESCR(20, 6) = 'PSDF    '
      ALLOW_CC_CMD_DESCR(21, 4) = '        ' ; ALLOW_CC_CMD_DESCR(21, 5) = 'ATOC    ' ; ALLOW_CC_CMD_DESCR(21, 6) = 'ATOC    '
      ALLOW_CC_CMD_DESCR(22, 4) = '        ' ; ALLOW_CC_CMD_DESCR(22, 5) = 'CRMS    ' ; ALLOW_CC_CMD_DESCR(22, 6) = 'CRMS    '
      ALLOW_CC_CMD_DESCR(23, 4) = '        ' ; ALLOW_CC_CMD_DESCR(23, 5) = 'RALL    ' ; ALLOW_CC_CMD_DESCR(23, 6) = 'RALL    '
      ALLOW_CC_CMD_DESCR(24, 4) = '        ' ; ALLOW_CC_CMD_DESCR(24, 5) = 'RPRINT  ' ; ALLOW_CC_CMD_DESCR(24, 6) = 'RPRINT  '
      ALLOW_CC_CMD_DESCR(25, 4) = '        ' ; ALLOW_CC_CMD_DESCR(25, 5) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(25, 6) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(26, 4) = '        ' ; ALLOW_CC_CMD_DESCR(26, 5) = 'NOPRINT ' ; ALLOW_CC_CMD_DESCR(26, 6) = 'NOPRINT '
      ALLOW_CC_CMD_DESCR(27, 4) = '        ' ; ALLOW_CC_CMD_DESCR(27, 5) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(27, 6) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(28, 4) = '        ' ; ALLOW_CC_CMD_DESCR(28, 5) = 'CID     ' ; ALLOW_CC_CMD_DESCR(28, 6) = 'CID     '
      ALLOW_CC_CMD_DESCR(29, 4) = '        ' ; ALLOW_CC_CMD_DESCR(29, 5) = '        ' ; ALLOW_CC_CMD_DESCR(29, 6) = '        '
      ALLOW_CC_CMD_DESCR(30, 4) = '        ' ; ALLOW_CC_CMD_DESCR(30, 5) = '        ' ; ALLOW_CC_CMD_DESCR(30, 6) = '        '
      ALLOW_CC_CMD_DESCR(31, 4) = '        ' ; ALLOW_CC_CMD_DESCR(31, 5) = '        ' ; ALLOW_CC_CMD_DESCR(31, 6) = '        '

!     =================SPCF=================   =================STRE================= ; =================STRN=================
      ALLOW_CC_CMD_DESCR( 1, 7) = 'SORT1   ' ; ALLOW_CC_CMD_DESCR( 1, 8) = 'SORT1   ' ; ALLOW_CC_CMD_DESCR( 1, 9) = 'SORT1   '
      ALLOW_CC_CMD_DESCR( 2, 7) = 'SORT2   ' ; ALLOW_CC_CMD_DESCR( 2, 8) = 'SORT2   ' ; ALLOW_CC_CMD_DESCR( 2, 9) = 'SORT2   '
      ALLOW_CC_CMD_DESCR( 3, 7) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 8) = 'PRINT   ' ; ALLOW_CC_CMD_DESCR( 3, 9) = 'PRINT   '
      ALLOW_CC_CMD_DESCR( 4, 7) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 8) = 'PUNCH   ' ; ALLOW_CC_CMD_DESCR( 4, 9) = 'PUNCH   '
      ALLOW_CC_CMD_DESCR( 5, 7) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 8) = 'PLOT    ' ; ALLOW_CC_CMD_DESCR( 5, 9) = 'PLOT    '
      ALLOW_CC_CMD_DESCR( 6, 7) = 'REAL    ' ; ALLOW_CC_CMD_DESCR( 6, 8) = 'REAL    ' ; ALLOW_CC_CMD_DESCR( 6, 9) = 'REAL    '
      ALLOW_CC_CMD_DESCR( 7, 7) = 'IMAG    ' ; ALLOW_CC_CMD_DESCR( 7, 8) = 'IMAG    ' ; ALLOW_CC_CMD_DESCR( 7, 9) = 'IMAG    '
      ALLOW_CC_CMD_DESCR( 8, 7) = 'MAG     ' ; ALLOW_CC_CMD_DESCR( 8, 8) = 'MAG     ' ; ALLOW_CC_CMD_DESCR( 8, 9) = 'MAG     '
      ALLOW_CC_CMD_DESCR( 9, 7) = 'PHASE   ' ; ALLOW_CC_CMD_DESCR( 9, 8) = 'PHASE   ' ; ALLOW_CC_CMD_DESCR( 9, 9) = 'PHASE   '
      ALLOW_CC_CMD_DESCR(10, 7) = '        ' ; ALLOW_CC_CMD_DESCR(10, 8) = 'VONMISES' ; ALLOW_CC_CMD_DESCR(10, 9) = 'VONMISES'
      ALLOW_CC_CMD_DESCR(11, 7) = '        ' ; ALLOW_CC_CMD_DESCR(11, 8) = 'MAXS    ' ; ALLOW_CC_CMD_DESCR(11, 9) = 'MAXS    '
      ALLOW_CC_CMD_DESCR(12, 7) = '        ' ; ALLOW_CC_CMD_DESCR(12, 8) = 'SHEAR   ' ; ALLOW_CC_CMD_DESCR(12, 9) = 'SHEAR   '
      ALLOW_CC_CMD_DESCR(13, 7) = '        ' ; ALLOW_CC_CMD_DESCR(13, 8) = 'STRCUR  ' ; ALLOW_CC_CMD_DESCR(13, 9) = 'STRCUR  '
      ALLOW_CC_CMD_DESCR(14, 7) = '        ' ; ALLOW_CC_CMD_DESCR(14, 8) = 'FIBER   ' ; ALLOW_CC_CMD_DESCR(14, 9) = 'FIBER   '
      ALLOW_CC_CMD_DESCR(15, 7) = '        ' ; ALLOW_CC_CMD_DESCR(15, 8) = 'CENTER  ' ; ALLOW_CC_CMD_DESCR(15, 9) = 'CENTER  '
      ALLOW_CC_CMD_DESCR(16, 7) = '        ' ; ALLOW_CC_CMD_DESCR(16, 8) = 'CORNER  ' ; ALLOW_CC_CMD_DESCR(16, 9) = 'CORNER  '
      ALLOW_CC_CMD_DESCR(17, 7) = '        ' ; ALLOW_CC_CMD_DESCR(17, 8) = 'BILIN   ' ; ALLOW_CC_CMD_DESCR(17, 9) = 'BILIN   '
      ALLOW_CC_CMD_DESCR(18, 7) = '        ' ; ALLOW_CC_CMD_DESCR(18, 8) = 'SGAGE   ' ; ALLOW_CC_CMD_DESCR(18, 9) = 'SGAGE   '
      ALLOW_CC_CMD_DESCR(19, 7) = '        ' ; ALLOW_CC_CMD_DESCR(19, 8) = 'CUBIC   ' ; ALLOW_CC_CMD_DESCR(19, 9) = 'CUBIC   '
      ALLOW_CC_CMD_DESCR(20, 7) = 'PSDF    ' ; ALLOW_CC_CMD_DESCR(20, 8) = 'PSDF    ' ; ALLOW_CC_CMD_DESCR(20, 9) = 'PSDF    '
      ALLOW_CC_CMD_DESCR(21, 7) = 'ATOC    ' ; ALLOW_CC_CMD_DESCR(21, 8) = 'ATOC    ' ; ALLOW_CC_CMD_DESCR(21, 9) = 'ATOC    '
      ALLOW_CC_CMD_DESCR(22, 7) = 'CRMS    ' ; ALLOW_CC_CMD_DESCR(22, 8) = 'CRMS    ' ; ALLOW_CC_CMD_DESCR(22, 9) = 'CRMS    '
      ALLOW_CC_CMD_DESCR(23, 7) = 'RALL    ' ; ALLOW_CC_CMD_DESCR(23, 8) = 'RALL    ' ; ALLOW_CC_CMD_DESCR(23, 9) = 'RALL    '
      ALLOW_CC_CMD_DESCR(24, 7) = 'RPRINT  ' ; ALLOW_CC_CMD_DESCR(24, 8) = 'RPRINT  ' ; ALLOW_CC_CMD_DESCR(24, 9) = 'RPRINT  '
      ALLOW_CC_CMD_DESCR(25, 7) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(25, 8) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(25, 9) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(26, 7) = 'NOPRINT ' ; ALLOW_CC_CMD_DESCR(26, 8) = 'NOPRINT ' ; ALLOW_CC_CMD_DESCR(26, 9) = 'NOPRINT '
      ALLOW_CC_CMD_DESCR(27, 7) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(27, 8) = 'RPUNCH  ' ; ALLOW_CC_CMD_DESCR(27, 9) = 'RPUNCH  '
      ALLOW_CC_CMD_DESCR(28, 7) = 'CID     ' ; ALLOW_CC_CMD_DESCR(28, 8) = '        ' ; ALLOW_CC_CMD_DESCR(28, 9) = '        '
      ALLOW_CC_CMD_DESCR(29, 7) = '        ' ; ALLOW_CC_CMD_DESCR(29, 8) = '        ' ; ALLOW_CC_CMD_DESCR(29, 9) = '        '
      ALLOW_CC_CMD_DESCR(30, 7) = '        ' ; ALLOW_CC_CMD_DESCR(30, 8) = '        ' ; ALLOW_CC_CMD_DESCR(30, 9) = '        '
      ALLOW_CC_CMD_DESCR(31, 7) = '        ' ; ALLOW_CC_CMD_DESCR(31, 8) = '        ' ; ALLOW_CC_CMD_DESCR(31, 9) = '        '

! Check that every word that is in CC_CMD_DESCRIBERS is a member of ALLOW_CC_CMD_DESCR.

ido_1:DO I=1,NUM_WORDS
jdo_1:   DO J=1,NUM_POSS_CCD
            FOUND = 'N'
            IF (CC_CMD_DESCRIBERS(I) == ALLOW_CC_CMD_DESCR(J,JCOL)) THEN
               FOUND = 'Y'
               EXIT jdo_1
            ENDIF
         ENDDO jdo_1
         IF (FOUND == 'N') THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CC_CMD_DESCRIBERS(I), OUTPUT_TYPE(JCOL)
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,101) CC_CMD_DESCRIBERS(I), OUTPUT_TYPE(JCOL)
               ENDIF
            ENDIF
         ENDIF
      ENDDO ido_1

      DO I=1,NUM_WORDS
         !write(ERR,*) "CC_CMD_DESCRIBERS(I)",I,CC_CMD_DESCRIBERS(I)
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'SORT2') THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,201) WHAT
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,201) WHAT
               ENDIF
            ENDIF
         ENDIF

         IF ((.NOT. IS_PUNCH) .AND. (CC_CMD_DESCRIBERS(I)(1:5) == 'PUNCH')) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,202) WHAT
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,202) WHAT
               ENDIF
            ENDIF
         ENDIF

         IF ((.NOT. IS_PLOT) .AND. (CC_CMD_DESCRIBERS(I)(1:4) == 'PLOT')) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,203) WHAT
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,203) WHAT
               ENDIF
            ENDIF
         ENDIF

      ENDDO

      IF ((WHAT == 'STRE') .OR. (WHAT == 'STRN')) THEN

         DO I=1,NUM_WORDS

            IF (CC_CMD_DESCRIBERS(I)(1:5) == 'STRCUR') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,301) WHAT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,301) WHAT
                  ENDIF
               ENDIF
            ENDIF

            IF (CC_CMD_DESCRIBERS(I)(1:5) == 'SGAGE') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,302) WHAT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,302) WHAT
                  ENDIF
               ENDIF
            ENDIF

            IF (CC_CMD_DESCRIBERS(I)(1:5) == 'CUBIC') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,303) WHAT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,303) WHAT
                  ENDIF
               ENDIF
            ENDIF

         ENDDO

      ENDIF

! Set values for variables in module CC_OUTPUT_DESCRIBERS

      IF (WHAT == 'STRE' ) THEN

         DO I=1,NUM_WORDS

            IF (NSUB <= 1) THEN
               IF      (CC_CMD_DESCRIBERS(I)(1:6) == 'CENTER') THEN
                  STRE_LOC = 'CENTER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:6) == 'CORNER') THEN
                  STRE_LOC = 'CORNER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'BILIN' ) THEN
                  STRE_LOC = 'CORNER'
               ENDIF
            ENDIF

            ! TODO: MISES is valid for VONMISES (test this...)
            IF      (CC_CMD_DESCRIBERS(I)(1:8) == 'VONMISES') THEN
               STRE_OPT = 'VONMISES'
            ELSE IF (CC_CMD_DESCRIBERS(I)(1:4) == 'MAXS'    ) THEN
               STRE_OPT = 'MAXS'
            ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'SHEAR'   ) THEN
               STRE_OPT = 'SHEAR'
            ENDIF

         ENDDO

      ENDIF

      IF (WHAT == 'STRN' ) THEN

         DO I=1,NUM_WORDS

            IF (NSUB <= 1) THEN
               IF      (CC_CMD_DESCRIBERS(I)(1:6) == 'CENTER') THEN
                  STRN_LOC = 'CENTER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:6) == 'CORNER') THEN
                  STRN_LOC = 'CORNER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'BILIN' ) THEN
                  STRN_LOC = 'CORNER'
               ENDIF
            ENDIF

            ! TODO: MISES is valid for VONMISES (test this...)
            IF      (CC_CMD_DESCRIBERS(I)(1:8) == 'VONMISES') THEN
               STRN_OPT = 'VONMISES'
            ELSE IF (CC_CMD_DESCRIBERS(I)(1:4) == 'MAXS'    ) THEN
               STRN_OPT = 'MAXS'
            ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'SHEAR'   ) THEN
               STRN_OPT = 'SHEAR'
            ENDIF

            IF      (CC_CMD_DESCRIBERS(I)(1:6) == 'STRCUR'  ) THEN
               STRN_CUR = 'STRCUR'
            ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'FIBER'   ) THEN
               STRN_CUR = 'FIBER'
            ENDIF

         ENDDO

      ENDIF

      IF (WHAT == 'ELFO' ) THEN

         DO I=1,NUM_WORDS

            IF (NSUB <= 1) THEN
               IF      (CC_CMD_DESCRIBERS(I)(1:6) == 'CENTER') THEN
                  FORC_LOC = 'CENTER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:6) == 'CORNER') THEN
                  FORC_LOC = 'CORNER'
               ELSE IF (CC_CMD_DESCRIBERS(I)(1:5) == 'BILIN' ) THEN
                  FORC_LOC = 'CORNER'
               ENDIF
            ENDIF

         ENDDO

      ENDIF




      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *WARNING    : "',A,'" IS NOT A RECOGNIZED CASE CONTROL COMMAND OPTION FOR OUTPUT TYPE "',A,'". OPTION IGNORED')

  201 FORMAT(' *WARNING    : "SORT2"    IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "SORT1" WILL BE USED')

  202 FORMAT(' *WARNING    : "PUNCH"    IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "PRINT" WILL BE USED')

  203 FORMAT(' *WARNING    : "PLOT"     IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "PLOT"   IGNORED')

  301 FORMAT(' *WARNING    : "STRCUR"   IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "STRCUR" IGNORED')

  302 FORMAT(' *WARNING    : "SGAGE"    IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "SGAGE"  IGNORED')

  303 FORMAT(' *WARNING    : "CUBIC"    IS NOT AN OPTION IN MYSTRAN FOR ',A,'. "CUBIC"  IGNORED')

 1204 FORMAT(' *ERROR  1205: PROGRAMMING ERROR IN SUBROUTINE ',A,'. INVALID VALUE ',A,' FOR VARIABLE "WHAT"')

! **********************************************************************************************************************************

      END SUBROUTINE CHK_CC_CMD_DESCRIBERS

      SUBROUTINE COMPUTE_OUTPUT_TARGETS ( REQUEST_F06 )

         USE PENTIUM_II_KIND, ONLY       :  LONG
         USE SCONTR, ONLY                :  CC_CMD_DESCRIBERS, NCCCD

         IMPLICIT NONE

         LOGICAL, INTENT(OUT)            :: REQUEST_F06
         LOGICAL                         :: FOUND_PRINT
         LOGICAL                         :: FOUND_PLOT
         INTEGER(LONG)                   :: I

         FOUND_PRINT = .FALSE.
         FOUND_PLOT  = .FALSE.
         DO I=1,NCCCD
            IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PRINT') FOUND_PRINT = .TRUE.
            IF (CC_CMD_DESCRIBERS(I)(1:4) == 'PLOT')  FOUND_PLOT  = .TRUE.
         ENDDO

         REQUEST_F06 = FOUND_PRINT .OR. .NOT. FOUND_PLOT

      END SUBROUTINE COMPUTE_OUTPUT_TARGETS


   END MODULE CASE_CONTROL_OUTPUTS
