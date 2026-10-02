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

   MODULE ELEMENT_RESULT_DISPATCH

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OFP3

   CONTAINS

      SUBROUTINE OFP3 ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Main driver routine for all element node (or engineering force) and stress and strain output requests for one subcase

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, WRT_FIJ
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MFIJ, MOGEL
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ANY_ELFE_OUTPUT, ANY_ELFN_OUTPUT, ANY_STRE_OUTPUT, ANY_STRN_OUTPUT
      USE LINK9_STUFF, ONLY           :  MAXREQ, OGEL

      USE ELEMENT_RESULT_OUTPUT, ONLY: OFP3_ELFN, OFP3_ELFE_1D, OFP3_ELFE_2D, OFP3_STRE_NO_PCOMP, OFP3_STRN_NO_PCOMP
      USE COMPOSITE_RESULT_OUTPUT, ONLY:  OFP3_STRE_PCOMP, OFP3_STRN_PCOMP
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3'

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG), PARAMETER        :: MERROR = 6        ! Number of error indicators used
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR(MERROR)    ! Local error count




! **********************************************************************************************************************************
! Initialize

      OGEL = ZERO

      DO I=1,MERROR
         IERROR(I) = 0
      ENDDO

      DO I=1,MFIJ
         WRT_FIJ(I) = 0
      ENDDO

! Call routines to calc element outputs
! 02/12/2020: Comment out the IF tests. Some of the routines have to be run if there are ELDATA requests to write to file .BUG

      OT4_EROW = 0
!      IF (ANY_ELFN_OUTPUT > 0) THEN
         CALL OFP3_ELFN          ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )
!     ENDIF

      OT4_EROW = 0
!     IF (ANY_ELFE_OUTPUT > 0) THEN
         CALL OFP3_ELFE_1D       ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OEF1X1, OEF1, OEF1X
         CALL OFP3_ELFE_2D       ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OEF1C, OEF1X1, OEF1, OEF1X
!     ENDIF

      OT4_EROW = 0
!     IF (ANY_STRE_OUTPUT > 0) THEN
         CALL OFP3_STRE_PCOMP    ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OES1C
         CALL OFP3_STRE_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OES1X1, OES1, OES1X,
!     ENDIF

      OT4_EROW = 0
!     IF (ANY_STRN_OUTPUT > 0) THEN
         CALL OFP3_STRN_PCOMP    ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OSTR1C
         CALL OFP3_STRN_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )  ! OSTR1X, OSTR1
!     ENDIF

! **********************************************************************************************************************************
      DO I=1,MERROR
         IF (IERROR(I) > 0) THEN
            WRITE(ERR,9996) SUBR_NAME
            WRITE(F06,9996) SUBR_NAME
            CALL OUTA_HERE ( 'Y' )                         ! Quit. (grid ID's not found and/or duplicate definitions in TSET))
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE LISTED ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE OFP3

   END MODULE ELEMENT_RESULT_DISPATCH
