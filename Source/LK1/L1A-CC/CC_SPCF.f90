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

      SUBROUTINE CC_SPCF ( CARD )

! Processes Case Control SPCF cards for SPC force output requests
      ! - SPCF: SPC force
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  PCHSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  SPCF_OUT
      USE MODEL_STUF, ONLY            :  SC_SPCF

      USE CC_SPCF_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SPCF'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE)              :: FOUND_PRINT       ! CC_CMD_DESCRIBERS has request for "PRINT"
      CHARACTER( 1*BYTE)              :: FOUND_PLOT        ! CC_CMD_DESCRIBERS has request for "PLOT"
      CHARACTER( 1*BYTE)              :: FOUND_PUNCH       ! CC_CMD_DESCRIBERS has request for "PUNCH"
      CHARACTER( 1*BYTE)              :: FOUND_NEU         ! CC_CMD_DESCRIBERS has request for "NEU"
      CHARACTER( 1*BYTE)              :: FOUND_CSV         ! CC_CMD_DESCRIBERS has request for "CSV"

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries
      ! (they all have some common code so it is put there)
      CALL CC_OUTPUTS ( CARD, 'SPCF', SETID )

      ! Check to see if PLOT, PRINT, PUNCH, NEU, CSV were in the request
      FOUND_PRINT = 'N'
      FOUND_PLOT  = 'N'
      FOUND_PUNCH = 'N'
      FOUND_NEU   = 'N'
      FOUND_CSV   = 'N'
      DO I=1,NCCCD
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PRINT') FOUND_PRINT = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:4) == 'PLOT')  FOUND_PLOT  = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PUNCH') FOUND_PUNCH = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'NEU')   FOUND_NEU   = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'CSV')   FOUND_CSV   = 'Y'
      ENDDO
      ! concatenate the strings
      SPCF_OUT = TRIM(FOUND_PRINT) // TRIM(FOUND_PLOT) // TRIM(FOUND_PUNCH) // TRIM(FOUND_NEU) // TRIM(FOUND_CSV)

! --- op2_upgraded begin --- !
      ! For bare "SPCF = ALL" requests, default to PRINT+PLOT so classic
      ! OP2 SPC force output (OQGV1) is emitted without requiring explicit
      ! "(PLOT)" qualifiers on the Case Control entry.
      IF (SPCF_OUT(1:5) == 'NNNNN') THEN
        SPCF_OUT = 'YYNNN'
      ENDIF
! --- op2_upgraded end --- !


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
