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

   MODULE ROD_BAR_BEAM_CARDS

   IMPLICIT NONE

   PRIVATE

   EXTERNAL :: ELEPRO

   PUBLIC :: BD_BAROR, BD_BAROR0, BD_BEAMOR, BD_BEAMOR0, BD_CBAR, BD_CBAR0, BD_CROD, BD_CONROD, BD_PROD, BD_PBAR, BD_PBARL, BD_PBEAM, BD_PLOTEL

   CONTAINS

      SUBROUTINE BD_BAROR ( CARD )

! Processes BAROR Bulk Data Cards. Reads and checks the property ID, if present, and the V vector,
! if present. The BAROR V vector type (BAROR_VVEC_TYPE) was determined in subr BAROR0

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LVVEC, NBAROR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BAROR_VVEC_TYPE, BAROR_G0, BAROR_VV, BAROR_PID

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_BAROR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: PGM_ERR   = 0     ! A  count of the number of coding errors


      REAL(DOUBLE)                    :: EPS1              ! A small value to compare zero to
      REAL(DOUBLE)                    :: R8INP             ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! BAROR Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    3      Property ID
!    6-8    V-Vector
!
      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read and check data

      IF (JCARD(3)(1:) /= ' ') THEN                        ! Read prop ID
         CALL I4FLD ( JCARD(3), JF(3), I4INP )
         IF (IERRFL(3) == 'N') THEN
            IF (I4INP /= BAROR_PID) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,11851) SUBR_NAME, JF(3), JCARD(1), BAROR_PID, I4INP
               WRITE(F06,11851) SUBR_NAME, JF(3), JCARD(1), BAROR_PID, I4INP
            ENDIF
            IF (I4INP <= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
               WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
            ENDIF
         ENDIF
      ENDIF

      IF (BAROR_VVEC_TYPE == 'VECTOR   ') THEN             ! Get V vec or grid no. Subr BAROR0 determined BAROR_VVEC_TYPE
         LVVEC = LVVEC + 1
         DO J=1,3
            IF (JCARD(J+5)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  IF (DABS(R8INP - BAROR_VV(J)) > EPS1) THEN
                     PGM_ERR = PGM_ERR + 1                 ! Coding error: This value doesn't agree with that read in LOADB0
                     WRITE(ERR,11852) SUBR_NAME, JF(J+5), JCARD(1), BAROR_VV(J), R8INP
                     WRITE(F06,11852) SUBR_NAME, JF(J+5), JCARD(1), BAROR_VV(J), R8INP
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ELSE IF (BAROR_VVEC_TYPE == 'GRID     ') THEN
         CALL I4FLD ( JCARD(6), JF(6), I4INP )
         IF (IERRFL(6) == 'N') THEN
            IF (I4INP /= BAROR_G0) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,11851) SUBR_NAME, JF(6), JCARD(1), BAROR_G0, I4INP
               WRITE(F06,11851) SUBR_NAME, JF(6), JCARD(1), BAROR_G0, I4INP
            ENDIF
            IF (I4INP <= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
               WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
            ENDIF
         ENDIF
      ELSE IF (BAROR_VVEC_TYPE == 'ERROR    ') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1102)
         WRITE(F06,1102)
      ENDIF

      CALL CARD_FLDS_NOT_BLANK ( JCARD,2,0,4,5,0,0,0,9 )   ! Issue warning if fields 2, 4, 5, 9 are not blank
      CALL BD_IMBEDDED_BLANK ( JCARD,0,3,0,0,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 3, 6-8
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF (PGM_ERR > 0) THEN                                ! PGM_ERR /= 0 is a coding error, so quit
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1102 FORMAT(' *ERROR  1102: ERROR IN SPECIFYING V VECTOR ON BAROR CARD. EITHER FIELD 6 MUST BE A POSITIVE INTEGER GRID POINT'     &
                    ,/,14X,' OR FIELDS 6, 7, 8 MUST CONTAIN REAL VECTOR COMPONENTS')

11851 FORMAT(' *ERROR  1185: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE READ FROM FIELD ',I3,' ON ',A,' CARD IN THIS SUBROUTINE'                                        &
                    ,/,14X,' DOES NOT AGREE WITH THAT READ FROM ORIGINAL BULK DATA DECK SCAN.'                                     &
                    ,/,14X,' VALUES ARE: ',I8,' (ORIGINAL BULK DATA SCAN),'                                                        &
                    ,/,14X,'        AND: ',I8,' (HERE)')

11852 FORMAT(' *ERROR  1185: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE READ FROM FIELD ',I3,' ON ',A,' CARD IN THIS SUBROUTINE'                                        &
                    ,/,14X,' DOES NOT AGREE WITH THAT READ FROM ORIGINAL BULK DATA DECK SCAN.'                                     &
                    ,/,14X,' VALUES ARE: ',1ES13.6,' (ORIGINAL BULK DATA SCAN),'                                                   &
                    ,/,14X,'        AND: ',1ES13.6,' (HERE)')

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_BAROR


      SUBROUTINE BD_BAROR0 ( CARD )

! Processes BAROR Bulk Data Card to increment LVVEC if the BAROR card
! has a V vector. Also, determine type of V vector is on this BAROR card.
! When this subr finishes, BAROR_VVEC_TYPE will be either:
!         a) 'VECTOR   ' means this BAROR card had V vector in fields 6-8
!         b) 'GRID     ' means this BAROR card had an integer in field 6
!                        and blank fields 7 and 8 (indicating a grid for V vector).
!         c) 'UNDEFINED' means this BAROR card had blank fields 6, 7 and 8
!         d) 'ERROR    ' means anything but (a), (b), or (c). Subr BD_BAROR will print error

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF, LVVEC, NBAROR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  BAROR_PID, BAROR_G0, BAROR_VV, BAROR_VVEC_TYPE, JBAROR

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_BAROR0'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count


      REAL(DOUBLE)                    :: R8INP             ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Count the number of BAROR cards. The count will be checked to make sure that there is no more than 1 in subr LOADB

      NBAROR = NBAROR + 1

! Only set BAROR values if this is the 1st BAROR card. There should be only 1 BAROR card.
! Reset IERRFL's to 'N', after testing, since CDRERR is not being called until BD_BAROR is called from LOADB

      IF (NBAROR == 1) THEN

         DO J=1,10                                         ! Set JBAROR to JCARD. We need JBAROR(3) (the prop ID) in subr ELEPRO
            JBAROR(J) = JCARD(J)
         ENDDO

         IF (JCARD(3)(1:) /= ' ') THEN
            CALL I4FLD ( JCARD(3), JF(3), I4INP )
            IF (IERRFL(3) == 'N') THEN
               BAROR_PID = I4INP
            ELSE
               JERR = JERR + 1
               IERRFL(3) = 'N'
            ENDIF
         ENDIF

         DO J=1,JCARD_LEN                                  ! See if there is an actual V vector.
            IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
               BAROR_VVEC_TYPE = 'VECTOR   '
               EXIT
            ENDIF
         ENDDO

         IF (BAROR_VVEC_TYPE == 'VECTOR   ') THEN          ! If there was an actual V vector, get components.

            LVVEC = LVVEC + 1
            JERR = 0
            DO J=1,3
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  BAROR_VV(J) = R8INP
               ELSE
                  JERR = JERR + 1
                  IERRFL(J+5) = 'N'                        ! Reset IERRFL - we don't want card field errors written from this subr
               ENDIF
            ENDDO
            IF (JERR /= 0) THEN
               BAROR_VVEC_TYPE = 'ERROR    '               ! Found err in V vector components, so reset BAROR_VVEC_TYPE to 'ERROR'
            ENDIF

         ELSE                                              ! Check to see if there is a grid no. for specifying VVEC

            IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN

               BAROR_VVEC_TYPE = 'GRID     '               ! We will check in subr BD_BAROR for read error
               CALL I4FLD ( JCARD(6), JF(6), I4INP )
               IF (IERRFL(6) == 'N') THEN
                  IF (I4INP > 0) THEN
                     BAROR_G0 = I4INP
                  ELSE
                     BAROR_VVEC_TYPE = 'ERROR    '         ! Found error in field 6, so reset VVEC_TYPE
                     IERRFL(6) = 'N'                       ! Reset IERRFL - we don't want card field errors written from this subr
                  ENDIF
               ENDIF

            ELSE

               IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
                  BAROR_VVEC_TYPE = 'UNDEFINED'
               ELSE
                  BAROR_VVEC_TYPE = 'ERROR    '
               ENDIF

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_BAROR0


      SUBROUTINE BD_BEAMOR ( CARD )

! Processes BEAMOR Bulk Data Cards. Reads and checks the property ID, if present, and the V vector,
! if present. The BEAMOR V vector type (BEAMOR_VVEC_TYPE) was determined in subr BEAMOR0

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LVVEC, NBEAMOR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BEAMOR_VVEC_TYPE, BEAMOR_G0, BEAMOR_VV, BEAMOR_PID

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_BEAMOR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: PGM_ERR   = 0     ! A  count of the number of coding errors


      REAL(DOUBLE)                    :: EPS1              ! A small value to compare zero to
      REAL(DOUBLE)                    :: R8INP             ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! BEAMOR Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    3      Property ID
!    6-8    V-Vector
!
      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read and check data

      IF (JCARD(3)(1:) /= ' ') THEN                        ! Read prop ID
         CALL I4FLD ( JCARD(3), JF(3), I4INP )
         IF (IERRFL(3) == 'N') THEN
            IF (I4INP /= BEAMOR_PID) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,11851) SUBR_NAME, JF(3), JCARD(1), BEAMOR_PID, I4INP
               WRITE(F06,11851) SUBR_NAME, JF(3), JCARD(1), BEAMOR_PID, I4INP
            ENDIF
            IF (I4INP <= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
            ENDIF
         ENDIF
      ENDIF

      IF (BEAMOR_VVEC_TYPE == 'VECTOR   ') THEN             ! Get V vec or grid no. Subr BEAMOR0 determined BEAMOR_VVEC_TYPE
         LVVEC = LVVEC + 1
         DO J=1,3
            IF (JCARD(J+5)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  IF (DABS(R8INP - BEAMOR_VV(J)) > EPS1) THEN
                     PGM_ERR = PGM_ERR + 1                 ! Coding error: This value doesn't agree with that read in LOADB0
                     WRITE(ERR,11852) SUBR_NAME, JF(J+5), JCARD(1), BEAMOR_VV(J), R8INP
                     WRITE(F06,11852) SUBR_NAME, JF(J+5), JCARD(1), BEAMOR_VV(J), R8INP
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ELSE IF (BEAMOR_VVEC_TYPE == 'GRID     ') THEN
         CALL I4FLD ( JCARD(6), JF(6), I4INP )
         IF (IERRFL(6) == 'N') THEN
            IF (I4INP /= BEAMOR_G0) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,11851) SUBR_NAME, JF(6), JCARD(1), BEAMOR_G0, I4INP
               WRITE(F06,11851) SUBR_NAME, JF(6), JCARD(1), BEAMOR_G0, I4INP
            ENDIF
            IF (I4INP <= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
               WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', I4INP
            ENDIF
         ENDIF
      ELSE IF (BEAMOR_VVEC_TYPE == 'ERROR    ') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1102)
         WRITE(F06,1102)
      ENDIF

      CALL CARD_FLDS_NOT_BLANK ( JCARD,2,0,4,5,0,0,0,9 )   ! Issue warning if fields 2, 4, 5, 9 are not blank
      CALL BD_IMBEDDED_BLANK ( JCARD,0,3,0,0,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 3, 6-8
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF (PGM_ERR > 0) THEN                                ! PGM_ERR /= 0 is a coding error, so quit
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1102 FORMAT(' *ERROR  1102: ERROR IN SPECIFYING V VECTOR ON BEAMOR CARD. EITHER FIELD 6 MUST BE A POSITIVE INTEGER GRID POINT'    &
                    ,/,14X,' OR FIELDS 6, 7, 8 MUST CONTAIN REAL VECTOR COMPONENTS')

11851 FORMAT(' *ERROR  1185: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE READ FROM FIELD ',I3,' ON ',A,' CARD IN THIS SUBROUTINE'                                        &
                    ,/,14X,' DOES NOT AGREE WITH THAT READ FROM ORIGINAL BULK DATA DECK SCAN.'                                     &
                    ,/,14X,' VALUES ARE: ',I8,' (ORIGINAL BULK DATA SCAN),'                                                        &
                    ,/,14X,'        AND: ',I8,' (HERE)')

11852 FORMAT(' *ERROR  1185: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE READ FROM FIELD ',I3,' ON ',A,' CARD IN THIS SUBROUTINE'                                        &
                    ,/,14X,' DOES NOT AGREE WITH THAT READ FROM ORIGINAL BULK DATA DECK SCAN.'                                     &
                    ,/,14X,' VALUES ARE: ',1ES13.6,' (ORIGINAL BULK DATA SCAN),'                                                   &
                    ,/,14X,'        AND: ',1ES13.6,' (HERE)')

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_BEAMOR


      SUBROUTINE BD_BEAMOR0 ( CARD )

! Processes BEAMOR Bulk Data Card to increment LVVEC if the BEAMOR card
! has a V vector. Also, determine type of V vector is on this BEAMOR card.
! When this subr finishes, BEAMOR_VVEC_TYPE will be either:
!         a) 'VECTOR   ' means this BEAMOR card had V vector in fields 6-8
!         b) 'GRID     ' means this BEAMOR card had an integer in field 6
!                        and blank fields 7 and 8 (indicating a grid for V vector).
!         c) 'UNDEFINED' means this BEAMOR card had blank fields 6, 7 and 8
!         d) 'ERROR    ' means anything but (a), (b), or (c). Subr BD_BEAMOR will print error

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF, LVVEC, NBEAMOR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  BEAMOR_PID, BEAMOR_G0, BEAMOR_VV, BEAMOR_VVEC_TYPE, JBEAMOR

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_BEAMOR0'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count


      REAL(DOUBLE)                    :: R8INP             ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Count the number of BEAMOR cards. The count will be checked to make sure that there is no more than 1 in subr LOADB

      NBEAMOR = NBEAMOR + 1

! Only set BEAMOR values if this is the 1st BEAMOR card. There should be only 1 BEAMOR card.
! Reset IERRFL's to 'N', after testing, since CDRERR is not being called until BD_BEAMOR is called from LOADB

      IF (NBEAMOR == 1) THEN

         DO J=1,10                                         ! Set JBEAMOR to JCARD. We need JBEAMOR(3) (the prop ID) in subr ELEPRO
            JBEAMOR(J) = JCARD(J)
         ENDDO

         IF (JCARD(3)(1:) /= ' ') THEN
            CALL I4FLD ( JCARD(3), JF(3), I4INP )
            IF (IERRFL(3) == 'N') THEN
               BEAMOR_PID = I4INP
            ELSE
               JERR = JERR + 1
               IERRFL(3) = 'N'
            ENDIF
         ENDIF

         DO J=1,JCARD_LEN                                  ! See if there is an actual V vector.
            IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
               BEAMOR_VVEC_TYPE = 'VECTOR   '
               EXIT
            ENDIF
         ENDDO

         IF (BEAMOR_VVEC_TYPE == 'VECTOR   ') THEN          ! If there was an actual V vector, get components.

            LVVEC = LVVEC + 1
            JERR = 0
            DO J=1,3
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  BEAMOR_VV(J) = R8INP
               ELSE
                  JERR = JERR + 1
                  IERRFL(J+5) = 'N'                        ! Reset IERRFL - we don't want card field errors written from this subr
               ENDIF
            ENDDO
            IF (JERR /= 0) THEN
               BEAMOR_VVEC_TYPE = 'ERROR    '               ! Found err in V vec components, so reset BEAMOR_VVEC_TYPE to 'ERROR'
            ENDIF

         ELSE                                              ! Check to see if there is a grid no. for specifying VVEC

            IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN

               BEAMOR_VVEC_TYPE = 'GRID     '               ! We will check in subr BD_BEAMOR for read error
               CALL I4FLD ( JCARD(6), JF(6), I4INP )
               IF (IERRFL(6) == 'N') THEN
                  IF (I4INP > 0) THEN
                     BEAMOR_G0 = I4INP
                  ELSE
                     BEAMOR_VVEC_TYPE = 'ERROR    '         ! Found error in field 6, so reset VVEC_TYPE
                     IERRFL(6) = 'N'                       ! Reset IERRFL - we don't want card field errors written from this subr
                  ENDIF
               ENDIF

            ELSE

               IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
                  BEAMOR_VVEC_TYPE = 'UNDEFINED'
               ELSE
                  BEAMOR_VVEC_TYPE = 'ERROR    '
               ENDIF

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_BEAMOR0


      SUBROUTINE BD_CBAR ( CARD, LARGE_FLD_INP )

      ! Processes CBAR and CBEAM Bulk Data Cards:
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LBAROFF, LVVEC, MEDAT_CBAR, &
                                         MEDAT_CBEAM, NBAROFF, NBAROR, NBEAMOR, NCBAR, NCBEAM, NEDAT, NELE, NVVEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BAROFF, BAROR_G0, BEAMOR_G0, BAROR_PID, BEAMOR_PID, BAROR_VVEC_TYPE, BEAMOR_VVEC_TYPE,    &
                                         BAROR_VV, BEAMOR_VV, EDAT, ETYPE, JBAROR, JBEAMOR, VVEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK, R8FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CBAR'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN=JCARD_LEN)        :: BAR_OR_BEAM       ! Field 1 of CBAR/CBEAM card
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: ELID              ! Field 2 of CBAR/CBEAM card
      CHARACTER( 1*BYTE)              :: FOUND     = 'N'   ! 'Y' if the V vec is one that is already stored in array VVEC
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER(LEN=JCARD_LEN)        :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 9*BYTE)              :: VVEC_TYPE         ! Type of V vector on this CBAR/CBEAM

      INTEGER(LONG)                   :: G0   = 0          ! Grid specifying V vector for this CBAR/CBEAM, if input
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: VVEC_NUM  = 0     ! V vector number


      REAL(DOUBLE)                    :: VV(3)             ! The 3 components of the V vector for this CBAR/CBEAM elem
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      ! CBAR element Bulk Data Card routine

      !   FIELD   ITEM           ARRAY ELEMENT
      !   -----   ------------   -------------
      !    1      Element type   ETYPE(nele)  =B1 for CBAR
      !    2      Element ID     EDAT(nedat+1)
      !    3      Property ID    EDAT(nedat+2)
      !    4      Grid A         EDAT(nedat+3)
      !    5      Grid B         EDAT(nedat+4)
      !    6-8    V-Vector       (see VVEC explanation below)
      !                          V vector key goes in EDAT(nedat+5)
      ! on optional second card:
      !    2      Pin Flag A     EDAT(nedat+6)
      !    3      Pin Flag B     EDAT(nedat+7)
      !    4-9    Offsets        (see BAROFF explanation below)
      !                          Offset key goes in EDAT(nedat+8)

      ! NOTES:

      ! If fields 3, 6-8 are blank, they are loaded with the data from  the BAROR/BEAMOR entry (these will remain blank if
      ! no BAROR/BEAMOR card exists). If V-vector is specfied via a grid point then EDAT(nedat+5) is set to that grid number.
      ! If V-vector is specified via an actual vector, the vector is loaded into array VVEC(NVVEC,J) (J=1,2,3) unless
      ! a vector equal to it has been put in VVEC. EDAT(nedat+5) is set equal to -NVVEC, where NVVEC is the row number
      ! in array VVEC.

      ! Offsets are in fields 4 - 9 of the first continuation card. If there are any offsets for this element, they are written to
      ! array BAROFF in row NBAROFF and NBAROFF is written in EDAT(nedat+8). If there are no offsets for this element, a zero is entered
      ! in array EDAT(nedat+8).
      EPS1 = EPSIL(1)

      ! Make JCARD from CARD
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      BAR_OR_BEAM = JCARD(1)
      ELID        = JCARD(2)

      ! Set JCARD_EDAT to JCARD
      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

      ! Initialize variables
      VVEC_TYPE = 'UNDEFINED'

      ! Check property ID field. Set to BAROR prop ID, if present, or to this elem ID, if not
      IF (JCARD(3)(1:) == ' ') THEN                        ! Prop ID field is blank, so use one of the following:
         IF (BAR_OR_BEAM(1:4) == 'CBAR') THEN
            IF (BAROR_PID /= 0) THEN                       ! Use BAROR prop ID for this CBAR prop ID
               JCARD_EDAT(3) = JBAROR(3)
            ELSE                                           ! Use CBAR  elem ID for this CBAR prop ID
               JCARD_EDAT(3) = JCARD(2)
            ENDIF
         ELSE
            IF (BEAMOR_PID /= 0) THEN                      ! Use BEAMOR prop ID for this CBEAM prop ID
               JCARD_EDAT(3) = JBEAMOR(3)
            ELSE                                           ! Use CBEAM  elem ID for this CBEAM prop ID
               JCARD_EDAT(3) = JCARD(2)
            ENDIF
         ENDIF
      ENDIF

      ! Call ELEPRO to increment NELE and load some of the connection data into array EDAT
      IF (BAR_OR_BEAM(1:4) == 'CBAR') THEN
         CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CBAR , 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
         NCBAR = NCBAR+1
         ETYPE(NELE) = 'BAR     '
      ELSE
         CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CBEAM, 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
         NCBEAM = NCBEAM+1
         ETYPE(NELE) = 'BEAM    '
      ENDIF

      ! Get the V vector for this CBAR/CBEAM (either from this CBAR/CBEAM or from BAROR/BEAMOR values, if present)
      ! Null all components. Some may be read from CBAR card
      DO J=1,3
         VV(J) = ZERO
      ENDDO

      IF (BAR_OR_BEAM(1:4) == 'CBAR') THEN
         IF (NBAROR > 0) THEN                              ! Set VVEC fields to BAROR values if this CBAR's VVEC fields are blank
            IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
               IF      (BAROR_VVEC_TYPE == 'GRID     ') THEN
                  VVEC_TYPE = 'BAROR_GRD'
                  G0        =  BAROR_G0
               ELSE IF (BAROR_VVEC_TYPE == 'VECTOR   ') THEN
                  VVEC_TYPE = 'BAROR_VEC'
                  DO J=1,3
                     VV(J) = BAROR_VV(J)
                  ENDDO
               ENDIF                                       ! We checked on BAROR_VVEC_TYPE = 'ERROR' in subr BD_BAROR and we will
            ENDIF                                          ! check case where fields 6, 7, 8 are blank below
         ENDIF
      ELSE
         IF (NBEAMOR > 0) THEN
            ! Set VVEC fields to BEAMOR values if this CBEAM's VVEC fields are blank
            IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
               IF      (BEAMOR_VVEC_TYPE == 'GRID     ') THEN
                  VVEC_TYPE = 'BEAMOR_GRD'
                  G0        =  BEAMOR_G0
               ELSE IF (BEAMOR_VVEC_TYPE == 'VECTOR   ') THEN
                  VVEC_TYPE = 'BEAMOR-VEC'
                  DO J=1,3
                     VV(J) = BEAMOR_VV(J)
                  ENDDO
               ENDIF                                       ! We checked on BEAMOR_VVEC_TYPE = 'ERROR' in subr BD_BEAMOR and we will
            ENDIF                                          ! check case where fields 6, 7, 8 are blank below
         ENDIF
      ENDIF

      IF (VVEC_TYPE == 'UNDEFINED') THEN                   ! We did not find a V vector so look for one on this CBAR/CBEAM
         DO J=1,JCARD_LEN                                  ! See if there is an actual V vector.
            IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
               VVEC_TYPE = 'VECTOR   '
               EXIT
            ENDIF
         ENDDO
         IF (VVEC_TYPE == 'VECTOR   ') THEN
            ! If there is an actual V vector, get components
            LVVEC = LVVEC + 1
            JERR = 0
            DO J=1,3
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  VV(J) = R8INP
               ELSE
                  JERR = JERR + 1
               ENDIF
            ENDDO
            IF (JERR /= 0) THEN
               VVEC_TYPE = 'ERROR    '                     ! Found error in V vector components, so reset VVEC_TYPE
            ENDIF
         ELSE                                              ! Check to see if there is a grid no. for specifying VVEC
            IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
               VVEC_TYPE = 'GRID     '
               CALL I4FLD ( JCARD(6), JF(6), I4INP )
               IF (IERRFL(6) == 'N') THEN
                  G0 = I4INP
                  IF (G0 < 0) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1187) JCARD(1), JCARD(2), G0
                     WRITE(F06,1187) JCARD(1), JCARD(2), G0
                  ENDIF
               ELSE
                  VVEC_TYPE = 'ERROR    '             ! Found error in field 6, so reset VVEC_TYPE
               ENDIF
            ELSE
               IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
                  VVEC_TYPE = 'UNDEFINED'
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1188) JCARD(1), JCARD(2)
                  WRITE(F06,1188) JCARD(1), JCARD(2)
               ELSE
                  VVEC_TYPE = 'ERROR    '
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1186) JCARD(1), JCARD(2)
                  WRITE(F06,1186) JCARD(1), JCARD(2)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      ! Load V vector data into EDAT and into VVEC, if not already there
      IF ((VVEC_TYPE == 'GRID     ') .OR. (VVEC_TYPE == 'BAROR_GRD')) THEN

         NEDAT = NEDAT + 1
         EDAT(NEDAT) = G0

      ELSE IF ((VVEC_TYPE == 'VECTOR   ') .OR. (VVEC_TYPE == 'BAROR_VEC')) THEN

         FOUND = 'N'
         DO J=1,NVVEC
            IF ((DABS(VV(1) - VVEC(J,1)) < EPS1) .AND. (DABS(VV(2) - VVEC(J,2)) < EPS1) .AND. (DABS(VV(3) - VVEC(J,3)) < EPS1)) THEN
               VVEC_NUM = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO

         IF (FOUND == 'N') THEN
            NVVEC = NVVEC + 1
            IF (NVVEC > LVVEC) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1132) SUBR_NAME,LVVEC
               WRITE(F06,1132) SUBR_NAME,LVVEC
               CALL OUTA_HERE ( 'Y' )                      ! Coding error, so quit
            ENDIF
            VVEC_NUM = NVVEC
            DO J=1,3
               VVEC(NVVEC,J) = VV(J)
            ENDDO
         ENDIF

         NEDAT = NEDAT + 1
         EDAT(NEDAT) = -VVEC_NUM
      ENDIF

       ! Write warnings and errors if any
      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 2-8
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      ! Optional Second Card:
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (ICONT == 1) THEN
         DO J = 2,3
            ! Get pin flag data, if present
            IF (JCARD(J)(1:) /= ' ') THEN
               CALL IP6CHK ( JCARD(J), JCARDO, IP6TYP, IDUM )
               IF (IP6TYP == 'COMP NOS') THEN
                  CALL I4FLD ( JCARDO, JF(J), I4INP )
                  IF (IERRFL(J) == 'N') THEN
                     NEDAT = NEDAT + 1
                     EDAT(NEDAT) = I4INP
                  ENDIF
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1130) JCARD(J) ,J, JCARD(1) ,ELID
                  WRITE(F06,1130) JCARD(J) ,J, JCARD(1), ELID
               ENDIF
            ELSE
               ! Null EDAT for this pin flag
               NEDAT = NEDAT + 1
               EDAT(NEDAT) = 0
            ENDIF
         ENDDO


         ! Get offsets, if present
         IF ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ') .OR. (JCARD(6)(1:) /= ' ') .OR.                                      &
             (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ') .OR. (JCARD(9)(1:) /= ' ')) THEN
            NBAROFF = NBAROFF + 1
            IF (NBAROFF > LBAROFF) THEN
               ! Coding error, so quit
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1161) SUBR_NAME, JCARD(1), LBAROFF
               WRITE(F06,1161) SUBR_NAME, JCARD(1), LBAROFF
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            NEDAT = NEDAT + 1
            EDAT(NEDAT) = NBAROFF
            DO J=1,6
               CALL R8FLD ( JCARD(J+3), JF(J+3), R8INP )
               IF (IERRFL(J+3) == 'N') THEN
                  BAROFF(NBAROFF,J) = R8INP
               ENDIF
            ENDDO
         ELSE
            ! Null EDAT for the offset flag
            NEDAT = NEDAT + 1
            EDAT(NEDAT) = 0
         ENDIF

         CALL BD_IMBEDDED_BLANK ( JCARD,0,0,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 4-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ELSE                                                 ! Null 2 pin flag, and 1 bar offset, fields in EDAT since no cont card

         DO J=1,3
            NEDAT = NEDAT + 1
            EDAT(NEDAT) = 0
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1130 FORMAT(' *ERROR  1130: INVALID PINFLAG = ',A,' IN FIELD ',I2,' ON CONTINUATION ENTRY OF ',A,' ID = ',A                       &
                    ,/,14X,' PINFLAGS CAN CONTAIN ONLY DIGITS 1-6')

 1132 FORMAT(' *ERROR  1132: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY V VECTORS. LIMIT IS ',I8)

 1161 FORMAT(' *ERROR  1161: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY CBAR/CBEAM OFFSETS. LIMIT IS ',I8)

 1186 FORMAT(' *ERROR  1186: ERROR IN SPECIFYING V VECTOR ON ',A,A,'. EITHER FIELD 6 MUST BE A POSITIVE INTEGER GRID POINT'        &
                    ,/,14X,' OR FIELDS 6, 7, 8 MUST CONTAIN REAL VECTOR COMPONENTS (WITH DECIMAL POINTS)')

 1187 FORMAT(' *ERROR  1187: GRID SPECIFYING V VECTOR ON ',A,A,' MUST BE > 0. VALUE IS = ',I8)

 1188 FORMAT(' *ERROR  1188: NO V VECTOR SPECIFIED FOR ',A,' ELEMENT ID = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CBAR


      SUBROUTINE BD_CBAR0 ( CARD, LARGE_FLD_INP )

! Processes CBAR or CBEAM Bulk Data Cards to increment LVVEC and LBAROFF if the CBAR or CBEAM entry has a V vector or offsets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LBAROFF, LVVEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CBAR0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 1*BYTE)              :: FND_VVEC          ! Indicator of whether there is an actual V vec on this CBAR card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! See if there is an actual V vector (not a grid point). If so, increment LVVEC

      FND_VVEC = 'N'
      DO J=1,JCARD_LEN
         IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
            FND_VVEC = 'Y'
            EXIT
         ENDIF
      ENDDO
      IF (FND_VVEC == 'Y') THEN
         LVVEC = LVVEC + 1
      ENDIF

! Optional Second Card - see if there are any offsets. If so, increment LBAROFF

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         IF ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ') .OR. (JCARD(6)(1:) /= ' ') .OR.                                      &
             (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ') .OR. (JCARD(9)(1:) /= ' ')) THEN
            LBAROFF = LBAROFF + 1
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CBAR0


      SUBROUTINE BD_CROD ( CARD )

! Processes CROD Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, MEDAT_CROD, NCROD, NEDAT, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CROD'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO

      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CROD element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =R1 for CROD
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid A         EDAT(nedat+3)
!    5      Grid B         EDAT(nedat+4)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to element ID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Read and check data

      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CROD, 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
      NCROD = NCROD+1
      ETYPE(NELE) = 'ROD     '

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,0,0,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-5
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )   ! Issue warning if fields 6, 7, 8, 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CROD


      SUBROUTINE BD_CONROD ( CARD )

! Processes CONROD Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF, MEDAT_CROD, NCROD, NELE, NEDAT, NPROD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, PROD, RPROD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CONROD'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD but with fields 5 and 6 switched to get G.P.'s together in EDAT

      INTEGER(LONG)                   :: ELEM_ID           ! Elem ID from field 2
      INTEGER(LONG)                   :: MATL_ID           ! Matl ID from field 5
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! An integer read
      INTEGER(LONG)                   :: IERR              ! Error count


      REAL(DOUBLE)                    :: R8INP             ! A real value read



! **********************************************************************************************************************************
! CONROD element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) = 'ROD     '
!    2      Element ID     EDAT(nedat+1)
!           Property ID    EDAT(nedat+2) (fictitous property PID = -EID)
!    3      Grid A         EDAT(nedat+3)
!    4      Grid B         EDAT(nedat+4)
!           PID            PROD(nprod,1)
!    5      MID            PROD(nprod,2)
!    6      Area          RPROD(nprod,1)
!    7      J             RPROD(nprod,2)
!    8      Tors stress C RPROD(nprod,3)
!    9      NSM           RPROD(nprod,4)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! First, check that fields 2-9 have the proper data type (we are going to have to rearrange the fields prior to calling ELEPRO).
! If any erors, return

      IERR = 0

      DO I=2,5
         CALL I4FLD ( JCARD(I), JF(I), I4INP )
         IF (IERRFL(I) == 'Y') IERR = IERR + 1
      ENDDO

      DO I=6,9
         CALL R8FLD ( JCARD(I), JF(I), R8INP )
         IF (IERRFL(I) == 'Y') IERR = IERR + 1
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )
      IF (IERR > 0) THEN
         RETURN
      ENDIF

! Set JCARD_EDAT to have fields 1-5 like a CROD

      DO I=1,10
         JCARD_EDAT(I)(1:) = ' '
      ENDDO
      JCARD_EDAT(1) = JCARD(1)                             ! Elem type
      JCARD_EDAT(2) = JCARD(2)                             ! Elem ID
      JCARD_EDAT(3) = JCARD(2)                             ! Prop ID = Elem ID for the time being
      JCARD_EDAT(4) = JCARD(3)                             ! Grid 1
      JCARD_EDAT(5) = JCARD(4)                             ! Grid 2

! Now call EDAT to check data as if this were a normal CROD. Here, ELEPRO needs prop ID = elem ID to work properly

      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CROD, 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
      NCROD = NCROD+1
      ETYPE(NELE) = 'ROD     '

! Now change PID in EDAT to -PID so EMG can find the properties

      EDAT(NEDAT-2) = -EDAT(NEDAT-2)

! Put property ID (neg of elem ID) and material ID's into array PROD

      NPROD = NPROD + 1

      ELEM_ID = 0
      CALL I4FLD ( JCARD(2), JF(2), ELEM_ID )
      IF (IERRFL(2) == 'N') THEN
         PROD(NPROD,1) = -ELEM_ID
      ENDIF

      MATL_ID = 0
      CALL I4FLD ( JCARD(5), JF(5), MATL_ID )
      IF (IERRFL(5) == 'N') THEN
         PROD(NPROD,2) = MATL_ID
      ENDIF

! Put real data from CONROD into array RPROD. We already checked that the data in theses fields can be read by R8FLD

      CALL R8FLD ( JCARD(6), JF(6), RPROD(NPROD,1) )
      CALL R8FLD ( JCARD(7), JF(7), RPROD(NPROD,2) )
      CALL R8FLD ( JCARD(8), JF(8), RPROD(NPROD,3) )
      CALL R8FLD ( JCARD(9), JF(9), RPROD(NPROD,4) )



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CONROD


      SUBROUTINE BD_PROD ( CARD )

! Processes PROD Bulk Data Cards. Reads and checks:

!  1) Property ID and material ID and enter them into array PROD
!  2) Area, torsional constant, stress recovery coeff for torsion and nonstructural mass and enter into array RPROD

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPROD, NPROD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PROD, RPROD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PROD'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID   = 0 ! Material ID (field 3 of this property card)
      INTEGER(LONG)                   :: PROPERTY_ID   = 0 ! Property ID (field 2 of this property card)




! **********************************************************************************************************************************
! PROD Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Prop ID         PROD(nprod,1)
!    3      Mat ID          PROD(nprod,2)
!    4      Area           RPROD(nprod,1)
!    5      Torsion J      RPROD(nprod,2)
!    6      Tors stress C  RPROD(nprod,3)
!    7      NSM            RPROD(nprod,4)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NPROD = NPROD+1

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read property ID and enter into array PROD
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPROD-1
            IF (PROPERTY_ID == PROD(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
            ENDIF
         ENDDO
         PROD(NPROD,1) = PROPERTY_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )          ! Read material ID and enter into array PROD
      IF (IERRFL(3) == 'N') THEN
         IF (MATERIAL_ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
         ELSE
            PROD(NPROD,2) = MATERIAL_ID
         ENDIF
      ENDIF

      DO J = 1,4                                           ! Read real property values in fields 4-7
         CALL R8FLD ( JCARD(J+3), JF(J+3), RPROD(NPROD,J) )
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8-9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PROD


      SUBROUTINE BD_PBAR ( CARD, LARGE_FLD_INP )

      ! Processes PBAR Bulk Data Cards. Reads and checks:
      !  1) Prop ID and Material ID and enter into array PBAR
      !  2) Area, moments of inertia, torsional constant ans nonstructural mass and enter into array RPBAR
      !  3) From 1st continuation card (if present): coords of 4 points for stress recovery and enter into array RPBAR
      !  4) From 2nd continuation card (if present): area factors for transverse shear and I12 and enter into array RPBAR

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  EPSIL, SUPINFO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BARTOR, IERRFL, FATAL_ERR, JCARD_LEN, JF, LPBAR, NPBAR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBAR, RPBAR

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE BAR_PROPERTY_VALIDATION, ONLY:  CHECK_BAR_MOIS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   =   'BD_PBAR'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format

      INTEGER(LONG)                   :: ICONT       = 0   ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0   ! Error indicator
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID = 0   ! Material ID (field 3 of this property card)
      INTEGER(LONG)                   :: PROPERTY_ID = 0   ! Property ID (field 2 of this property card)


      REAL(DOUBLE)                    :: I1          = ZERO! Moment of inertia
      REAL(DOUBLE)                    :: I2          = ZERO! Moment of inertia
      REAL(DOUBLE)                    :: I12         = ZERO! Product of inertia
      REAL(DOUBLE)                    :: EPS1              ! A small number



! **********************************************************************************************************************************
      ! PBAR Bulk Data Card routine
      !
      !   FIELD   ITEM                                            ARRAY ELEMENT
      !   -----   ----                                            -------------
      !    2      Property ID                                     PBAR(npbar, 1)
      !    3      Material ID                                     PBAR(npbar, 2)
      !    4      Area                                           RPBAR(npbar, 1)
      !    5      Inertia 1, I1                                  RPBAR(npbar, 2)
      !    6      Inertia 2, I2                                  RPBAR(npbar, 3)
      !    7      Torsional constant, J                          RPBAR(npbar, 4)
      !    8      Non-structural mass                            RPBAR(npbar, 5)
      !
      ! on optional second card:
      !    2      Y1 y coord of 1st point for stress recovery)   RPBAR(npbar, 6)
      !    3      Z1 z coord of 1st point for stress recovery)   RPBAR(npbar, 7)
      !    4      Y2 y coord of 2nd point for stress recovery)   RPBAR(npbar, 8)
      !    5      Z2 z coord of 2nd point for stress recovery)   RPBAR(npbar, 9)
      !    6      Y3 y coord of 3rd point for stress recovery)   RPBAR(npbar,10)
      !    7      Z3 z coord of 3rd point for stress recovery)   RPBAR(npbar,11)
      !    8      Y4 y coord of 4th point for stress recovery)   RPBAR(npbar,12)
      !    9      Z4 z coord of 4th point for stress recovery)   RPBAR(npbar,13)
      ! on optional third card:
      !    2      K1, Plane 1 shear factor                       RPBAR(npbar,14)
      !    3      K2, Plane 2 shear factor                       RPBAR(npbar,15)
      !    4      I12, Product of inertia                        RPBAR(npbar,16)
      !    5      C Torsional stress recovery coefficient        RPBAR(npbar,17)
      EPS1 = EPSIL(1)
      CHILD(1:) = ' '

      ! Make JCARD from CARD
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      ! Check for overflow
      NPBAR = NPBAR+1

     ! Read and check data on parent card
      ID   = JCARD(2)
      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read property ID and enter into array PBAR
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPBAR-1
            IF (PROPERTY_ID == PBAR(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
             ENDIF
         ENDDO
         PBAR(NPBAR,1) = PROPERTY_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )          ! Read material ID and enter into array PBAR
      IF (IERRFL(3) == 'N') THEN
         IF (MATERIAL_ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
         ELSE
            PBAR(NPBAR,2) = MATERIAL_ID
         ENDIF
      ENDIF

      ! Read real property values in fields 4-8
      DO J = 1,5
         CALL R8FLD ( JCARD(J+3), JF(J+3), RPBAR(NPBAR,J) )
      ENDDO
      IF (IERRFL(5)  == 'N') THEN
         I1 = RPBAR(NPBAR,2)
      ENDIF
      IF (IERRFL(6)  == 'N') THEN
         I2 = RPBAR(NPBAR,3)
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 2-8
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      ! Read and check data on optional 2nd, 3rd cards:
      ! Init I12 since the 2nd cont entry, which would read it, may not exist
      I12 = ZERO

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )                 ! Read 2nd card
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (ICONT == 1) THEN
         DO J = 6,13                                       ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-4), JF(J-4), RPBAR(NPBAR,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )              ! Read 3rd card
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         flush(err)
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J = 14,16                                   ! Read real property values in fields 2-4 of 3rd card
               CALL R8FLD ( JCARD(J-12), JF(J-12), RPBAR(NPBAR,J) )
            ENDDO
            IF (IERRFL(4) == 'N') THEN
               I12 = RPBAR(NPBAR,16)
            ENDIF
                                                           ! Read torsional stress coefficient
            CALL R8FLD ( JCARD(5), JF(5), RPBAR(NPBAR,17) )
            IF (DABS(RPBAR(NPBAR,17)) > EPS1) THEN         ! Set BARTOR to 'Y' if any BAR elem has a nonzero torsional stress coeff
               BARTOR = 'Y'
            ENDIF

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,0,0,0)! Make sure that there are no imbedded blanks in fields 2-5
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

         ENDIF
      ENDIF

      ! Call subr to check sensibility of I1, I2, I12 combinations
      CALL CHECK_BAR_MOIs ( 'PBAR', ID, I1, I2, I12, IERR )
      RPBAR(NPBAR, 2) = I1
      RPBAR(NPBAR, 3) = I2
      RPBAR(NPBAR,16) = I12
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE BD_PBAR


      SUBROUTINE BD_PBARL ( CARD, LARGE_FLD_INP, PBARL_TYPE )

! Processes PBARL Bulk Data Cards. Reads and checks:

!  1) Prop ID and Material ID and enter into array PBAR
!  2) Read data on type of crossection and dimensions and calculate section props
!  3) Create an equivalent PBAR entry (data goes into arrays PBAR, RPBAR)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE DERIVED_DATA_TYPES, ONLY    :  CHAR1_INT1
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, IERRFL, FATAL_ERR, JCARD_LEN, JF, LPBAR, NPBAR, NPBARL
      USE PARAMS, ONLY                :  EPSIL, PBARLSHR, SUPINFO
      USE CONSTANTS_1, ONLY           :  PI, ZERO, QUARTER, THIRD, HALF, ONE, TWO, THREE, FOUR, FIVE, SIX, SEVEN, EIGHT, NINE,     &
                                         TEN, TWELVE

      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBAR, RPBAR

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, LEFT_ADJ_BDFLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE BAR_PROPERTY_VALIDATION, ONLY:  CHECK_BAR_MOIS

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: NS        = 25    ! Dimension of array BAR_SHAPE
      TYPE(CHAR1_INT1)                :: BAR_SHAPE(NS)     ! Array with the BAR crossection type and the num of D(i) needed for it

      CHARACTER(LEN(BLNK_SUB_NAM))    :: SUBR_NAME   =   'BD_PBARL'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(OUT)   :: PBARL_TYPE        ! Name of the cross-section (e.g. I, BAR, etc)
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: CS_TYPE           ! Name of the cross-section (e.g. I, BAR, etc)
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 1*BYTE)              :: FOUND             ! 'Y' if a correct name for a cross-section was found in field 5
      CHARACTER(LEN(JCARD))           :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(99*BYTE)              :: MSG               ! Error message written out
      CHARACTER(LEN(JCARD))           :: NAME              ! JCARD(1) from parent entry

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: KD                ! Counter
      INTEGER(LONG)                   :: NUM_D             ! Number of D(i) values to read from the continuation entries
      INTEGER(LONG)                   :: MATL_ID   = 0     ! Material ID (field 3 of this property card)
      INTEGER(LONG)                   :: PROP_ID   = 0     ! Property ID (field 2 of this property card)


      REAL(DOUBLE)                    :: AREA      = ZERO  ! Cross-sectional area
      REAL(DOUBLE)                    :: D(NS)             ! Dimensions of cross-secion of the bar
      REAL(DOUBLE)                    :: I1,I2     = ZERO  ! Moments of inertia
      REAL(DOUBLE)                    :: K1,K2     = ZERO  ! Shear flex factors
      REAL(DOUBLE)                    :: I12       = ZERO  ! Product of inertia
      REAL(DOUBLE)                    :: JTOR      = ZERO  ! Torsional constant
      REAL(DOUBLE)                    :: NSM       = ZERO  ! Non structural mass
      REAL(DOUBLE)                    :: Y(4)      = ZERO  ! Y coords in cross-section for 4 points of data recovery
      REAL(DOUBLE)                    :: Z(4)      = ZERO  ! Z coords in cross-section for 4 points of data recovery
      REAL(DOUBLE)                    :: R8INP             ! A real value read from a field on this PBARL entry



! **********************************************************************************************************************************
! PBAR Bulk Data Card routine

!   FIELD   ITEM                                            ARRAY ELEMENT
!   -----   ----                                            -------------
!    2      Property ID                                     PBAR(npbar, 1)
!    3      Material ID                                     PBAR(npbar, 2)
!    4      GROUP   ***(NOT USED)***
!    5      TYPE (type of cross-section)

! on mandatory second card: (the number of D(i) are set according to TYPE. See BAR_SHAPE array)
!    2      D(1)
!    3      D(2)
!    4      D(3)
!    5      D(4)
!    6      D(5)
!    7      D(6)
!    8      D(7)
!    9      D(8)
! on optional cards:
!    2      D(9)
!    .       .
!    .       .
!    .      NSM (last entry)


      ! Initialize
      AREA = ZERO
      I1   = ZERO
      I2   = ZERO
      I12  = ZERO
      JTOR = ZERO
      Y(1) = ZERO;   Z(1) = ZERO
      Y(2) = ZERO;   Z(2) = ZERO
      Y(3) = ZERO;   Z(3) = ZERO
      Y(4) = ZERO;   Z(4) = ZERO
      K1   = ZERO
      K2   = ZERO

      NSM  = ZERO

      DO I=1,NS
         BAR_SHAPE(I)%Col_1(1:) = ' '
         BAR_SHAPE(I)%Col_2     = 0
         D(I)                   = ZERO
      ENDDO

      ! no idea what the Col_1 and Col_2 part refer to, but it's still
      ! easy to add a new section
      BAR_SHAPE( 1)%Col_1(1:8) = 'BAR     ';  BAR_SHAPE( 1)%Col_2 = 2
      BAR_SHAPE( 2)%Col_1(1:8) = 'BOX     ';  BAR_SHAPE( 2)%Col_2 = 4
      BAR_SHAPE( 3)%Col_1(1:8) = 'BOX1    ';  BAR_SHAPE( 3)%Col_2 = 6
      BAR_SHAPE( 4)%Col_1(1:8) = 'CHAN    ';  BAR_SHAPE( 4)%Col_2 = 4
      BAR_SHAPE( 5)%Col_1(1:8) = 'CHAN1   ';  BAR_SHAPE( 5)%Col_2 = 4
      BAR_SHAPE( 6)%Col_1(1:8) = 'CHAN2   ';  BAR_SHAPE( 6)%Col_2 = 4
      BAR_SHAPE( 7)%Col_1(1:8) = 'CROSS   ';  BAR_SHAPE( 7)%Col_2 = 4
      BAR_SHAPE( 8)%Col_1(1:8) = 'H       ';  BAR_SHAPE( 8)%Col_2 = 4
      BAR_SHAPE( 9)%Col_1(1:8) = 'HAT     ';  BAR_SHAPE( 9)%Col_2 = 4
      BAR_SHAPE(10)%Col_1(1:8) = 'HEXA    ';  BAR_SHAPE(10)%Col_2 = 3
      BAR_SHAPE(11)%Col_1(1:8) = 'I       ';  BAR_SHAPE(11)%Col_2 = 6
      BAR_SHAPE(12)%Col_1(1:8) = 'I1      ';  BAR_SHAPE(12)%Col_2 = 4
      BAR_SHAPE(13)%Col_1(1:8) = 'ROD     ';  BAR_SHAPE(13)%Col_2 = 1
      BAR_SHAPE(14)%Col_1(1:8) = 'T       ';  BAR_SHAPE(14)%Col_2 = 4
      BAR_SHAPE(15)%Col_1(1:8) = 'T1      ';  BAR_SHAPE(15)%Col_2 = 4
      BAR_SHAPE(16)%Col_1(1:8) = 'T2      ';  BAR_SHAPE(16)%Col_2 = 4
      BAR_SHAPE(17)%Col_1(1:8) = 'TUBE    ';  BAR_SHAPE(17)%Col_2 = 2
      BAR_SHAPE(18)%Col_1(1:8) = 'Z       ';  BAR_SHAPE(18)%Col_2 = 4
      BAR_SHAPE(19)%Col_1(1:8) = 'TUBE2   ';  BAR_SHAPE(19)%Col_2 = 2

      NPBAR         = NPBAR  + 1
      PBAR(NPBAR,3) = NPBARL

      ! Make JCARD from CARD
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NAME = JCARD(1)
      ID   = JCARD(2)

      ! Code not written for cross-section dimension data spilling over
      ! to a 2nd cont entry so give error if BAR_SHAPE(i)%Col_2 > 7
      DO I=1,NS
         IF (BAR_SHAPE(I)%Col_2 > 7) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1172) SUBR_NAME, NAME
            WRITE(F06,1172) SUBR_NAME, NAME
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDDO

      ! Read and check data on parent card
      CALL I4FLD ( JCARD(2), JF(2), PROP_ID )              ! Read property ID and enter into array PBAR
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPBAR-1
            IF (PROP_ID == PBAR(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROP_ID
               WRITE(F06,1145) JCARD(1),PROP_ID
               EXIT
             ENDIF
         ENDDO
         PBAR(NPBAR,1) = PROP_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), MATL_ID )              ! Read material ID and enter into array PBAR
      IF (IERRFL(3) == 'N') THEN
         IF (MATL_ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATL_ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATL_ID
         ELSE
            PBAR(NPBAR,2) = MATL_ID
         ENDIF
      ENDIF

      NUM_D = 0
      CALL LEFT_ADJ_BDFLD ( JCARD(5) )                     ! Read cross-section type in field 5
      CALL CHAR_FLD ( JCARD(5), JF(5), CS_TYPE )
      FOUND = 'N'
      DO I=1,NS
         IF (CS_TYPE(1:8) == BAR_SHAPE(I)%Col_1) THEN
            NUM_D = BAR_SHAPE(I)%Col_2
            FOUND = 'Y'
         ENDIF
      ENDDO
      IF (FOUND == 'N') THEN                               ! If a valid CS_TYPE not found give error and return
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1162) CS_TYPE(1:8), JF(5), NAME, JCARD(5)
         WRITE(F06,1162) CS_TYPE(1:8), JF(5), NAME, JCARD(5)
         RETURN
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,5,0,0,0,0 )     ! Make sure that there are no imbedded blanks in fields 2,3,5
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,0,6,7,8,9 )   ! Issue warning if fields 4 and 6-9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      PBARL_TYPE(1:) = ' '
      PBARL_TYPE(1:) = CS_TYPE(1:)

      ! Read D(i) dimension values on continuation entries
      IERR = 0                                             ! Mandatory 1st continuation entry
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         KD = 0
D_do1:   DO J=2,9                                          ! --- Read cross-section dimension data
            IF (KD < NUM_D) THEN
               IF (JCARD(J)(1:) /= ' ') THEN
                  CALL R8FLD ( JCARD(J), JF(J), R8INP )
                  IF (IERRFL(J) == 'N') THEN
                     KD = KD + 1
                     D(KD) = R8INP
                     IF (D(KD) <= ZERO) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1174) NAME, ID, JF(J), D(KD)
                        WRITE(F06,1174) NAME, ID, JF(J), D(KD)
                        CYCLE D_do1
                     ENDIF
                     CYCLE D_do1
                  ENDIF
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1171) NAME, ID, CS_TYPE(1:8), NUM_D
                  WRITE(F06,1171) NAME, ID, CS_TYPE(1:8), NUM_D
                  EXIT D_do1
               ENDIF
            ELSE
               EXIT D_do1
            ENDIF
         ENDDO D_do1
         IF (KD <= 7) THEN
            CALL R8FLD( JCARD(KD+2), JF(KD+2), R8INP )        ! ---  Read NSM
            IF (IERRFL(KD+1) == 'N') THEN
               NSM = R8INP
            ENDIF
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF


      ! Call routines to calc bar cross-section props, depending on CS_TYPE
      IERR = 0
      IF      (CS_TYPE(1:8) == BAR_SHAPE( 1)%Col_1) THEN;  NUM_D = BAR_SHAPE( 1)%Col_2;   CALL SECTION_PROPS_BAR   ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 2)%Col_1) THEN;  NUM_D = BAR_SHAPE( 2)%Col_2;   CALL SECTION_PROPS_BOX   ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 3)%Col_1) THEN;  NUM_D = BAR_SHAPE( 3)%Col_2;   CALL SECTION_PROPS_BOX1  ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 4)%Col_1) THEN;  NUM_D = BAR_SHAPE( 4)%Col_2;   CALL SECTION_PROPS_CHAN  ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 5)%Col_1) THEN;  NUM_D = BAR_SHAPE( 5)%Col_2;   CALL SECTION_PROPS_CHAN1 ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 6)%Col_1) THEN;  NUM_D = BAR_SHAPE( 6)%Col_2;   CALL SECTION_PROPS_CHAN2 ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 7)%Col_1) THEN;  NUM_D = BAR_SHAPE( 7)%Col_2;   CALL SECTION_PROPS_CROSS ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 8)%Col_1) THEN;  NUM_D = BAR_SHAPE( 8)%Col_2;   CALL SECTION_PROPS_H     ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE( 9)%Col_1) THEN;  NUM_D = BAR_SHAPE( 9)%Col_2;   CALL SECTION_PROPS_HAT   ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(10)%Col_1) THEN;  NUM_D = BAR_SHAPE(10)%Col_2;   CALL SECTION_PROPS_HEXA  ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(11)%Col_1) THEN;  NUM_D = BAR_SHAPE(11)%Col_2;   CALL SECTION_PROPS_I     ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(12)%Col_1) THEN;  NUM_D = BAR_SHAPE(12)%Col_2;   CALL SECTION_PROPS_I1    ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(13)%Col_1) THEN;  NUM_D = BAR_SHAPE(13)%Col_2;   CALL SECTION_PROPS_ROD   ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(14)%Col_1) THEN;  NUM_D = BAR_SHAPE(14)%Col_2;   CALL SECTION_PROPS_T     ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(15)%Col_1) THEN;  NUM_D = BAR_SHAPE(15)%Col_2;   CALL SECTION_PROPS_T1    ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(16)%Col_1) THEN;  NUM_D = BAR_SHAPE(16)%Col_2;   CALL SECTION_PROPS_T2    ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(17)%Col_1) THEN;  NUM_D = BAR_SHAPE(17)%Col_2;   CALL SECTION_PROPS_TUBE  ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(18)%Col_1) THEN;  NUM_D = BAR_SHAPE(18)%Col_2;   CALL SECTION_PROPS_Z     ( IERR )
      ELSE IF (CS_TYPE(1:8) == BAR_SHAPE(19)%Col_1) THEN;  NUM_D = BAR_SHAPE(19)%Col_2;   CALL SECTION_PROPS_TUBE2 ( IERR )
      !HAT1
      !DBOX
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1162) JCARD(5), JF(5), NAME, ID
         WRITE(F06,1162) JCARD(5), JF(5), NAME, ID
         RETURN
      ENDIF

      IF (IERR == 0) THEN

         RPBAR(npbar, 1) = AREA
         RPBAR(npbar, 2) = I1
         RPBAR(npbar, 3) = I2
         RPBAR(npbar, 4) = JTOR
         RPBAR(npbar, 5) = NSM
         RPBAR(npbar, 6) = Y(1)
         RPBAR(npbar, 7) = Z(1)
         RPBAR(npbar, 8) = Y(2)
         RPBAR(npbar, 9) = Z(2)
         RPBAR(npbar,10) = Y(3)
         RPBAR(npbar,11) = Z(3)
         RPBAR(npbar,12) = Y(4)
         RPBAR(npbar,13) = Z(4)
         IF (PBARLSHR == 'Y') THEN
            RPBAR(npbar,14) = K1
            RPBAR(npbar,15) = K2
         ELSE
            RPBAR(npbar,14) = ZERO
            RPBAR(npbar,15) = ZERO
         ENDIF
         RPBAR(npbar,16) = I12

      CALL CHECK_BAR_MOIs ( 'PBARL', ID, I1, I2, I12, IERR )     ! Call subr to check sensibility of I1, I2, I12 combinations
      RPBAR(NPBAR, 2) = I1
      RPBAR(NPBAR, 3) = I2
      RPBAR(NPBAR,16) = I12
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      ELSE

         MSG(1:) = ' '
         IF (IERR ==  1) MSG = 'INVALID DIMENSIONS FOR CROSS-SECTION'
         IF (IERR == 99) MSG = 'CODE NOT WRITTEN FOR THIS CROSS-SECTION YET'
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1177) NAME, ID, CS_TYPE(1:8), MSG
         WRITE(F06,1177) NAME, ID, CS_TYPE(1:8), MSG

      ENDIF

! Write PBAR equivalent Bulk Data entries:
!
!     IF (ECHO(1:4) /= 'NONE') THEN
!        CALL WRITE_PBAR_EQUIV
!     ENDIF
!


      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1162 FORMAT(' *ERROR  1162: INVALID NAME = "',A,'" FOR CROSS-SECTION IN FIELD ',I2,' ON ',A,' ENTRY = ',A,                        &
                          '. CHECK USERS MANUAL FOR VALID ENTRIES')

 1171 FORMAT(' *ERROR  1171: CONTINUATION ENTRY  FOR ',A,A,' DOES NOT HAVE ENOUGH DIMENSION DATA FOR CROSS-SECTION "',A,'".',      &
                           ' ENTRY SHOULD HAVE ',I2,' VALUES')

 1172 FORMAT(' *ERROR  1172: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CODE NOT WRITTEN TO READ MORE THAN 7 CROSS-SECTION DIMENSIONS FOR ',A)

 1174 FORMAT(' *ERROR  1174: CROSS SECTION DIMENSIONS ON ',A,A,' MUST BE > 0 BUT FIELD ',I2,' HAS VALUE ',1ES10.2)

 1177 FORMAT(' *ERROR  1177: ',A,A,' FOR "',A,'" CROSS-SECTION: ',A)

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)


! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE WRITE_PBAR_EQUIV

      USE PARAMS, ONLY                :  PBARLDEC, PBARLSHR
      USE SCONTR, ONLY                :  JCARD_LEN

      IMPLICIT NONE

      CHARACTER( 6*BYTE)              :: FMT0              ! Format used in writing data to char array ICARD
      CHARACTER( 8*BYTE)              :: FMT1              ! Format used in writing data to char array ICARD
      CHARACTER(LEN=JCARD_LEN)        :: ICARD(2:9)        ! Char array for fields 2-9 of equiv PBAR B.D. entries

! **********************************************************************************************************************************
      WRITE(F06,101)
      WRITE(F06,102) ID
      WRITE(F06,104)

      FMT0 = '(8X)'

      IF      (PBARLDEC == 0) THEN
         FMT1 = '(F8.0)'
      ELSE IF (PBARLDEC == 1) THEN
         FMT1 = '(F8.1)'
      ELSE IF (PBARLDEC == 2) THEN
         FMT1 = '(F8.2)'
      ELSE IF (PBARLDEC == 3) THEN
         FMT1 = '(F8.3)'
      ELSE IF (PBARLDEC == 4) THEN
         FMT1 = '(F8.4)'
      ELSE IF (PBARLDEC == 5) THEN
         FMT1 = '(F8.5)'
      ELSE IF (PBARLDEC == 6) THEN
         FMT1 = '(F8.6)'
      ELSE
         FMT1 = FMT0                                       ! Write blank fields if PBARLDEC is not caught being out of range 0-6
      ENDIF

      WRITE(ICARD(2),201)  PROP_ID   ;   CALL LEFT_ADJ_BDFLD ( ICARD(2) )
      WRITE(ICARD(3),201)  MATL_ID   ;   CALL LEFT_ADJ_BDFLD ( ICARD(3) )
      WRITE(ICARD(4),FMT1) AREA      ;   CALL LEFT_ADJ_BDFLD ( ICARD(4) )
      WRITE(ICARD(5),FMT1) I1;       ;   CALL LEFT_ADJ_BDFLD ( ICARD(5) )
      WRITE(ICARD(6),FMT1) I2;       ;   CALL LEFT_ADJ_BDFLD ( ICARD(6) )
      WRITE(ICARD(7),FMT1) JTOR;     ;   CALL LEFT_ADJ_BDFLD ( ICARD(7) )
      WRITE(ICARD(8),FMT1) NSM;      ;   CALL LEFT_ADJ_BDFLD ( ICARD(8) )
      WRITE(ICARD(9),FMT0)
      WRITE(F06,301) (ICARD(I),I=2,8)                      ! Write parent entry

      WRITE(ICARD(2),FMT1) Y(1)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(2) )
      WRITE(ICARD(3),FMT1) Z(1)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(3) )
      WRITE(ICARD(4),FMT1) Y(2)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(4) )
      WRITE(ICARD(5),FMT1) Z(2)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(5) )
      WRITE(ICARD(6),FMT1) Y(3)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(6) )
      WRITE(ICARD(7),FMT1) Z(3)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(7) )
      WRITE(ICARD(8),FMT1) Y(4)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(8) )
      WRITE(ICARD(9),FMT1) Z(4)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(9) )
      WRITE(F06,302) (ICARD(I),I=2,9)                      ! Write 1st cont entry

      IF (PBARLSHR == 'Y') THEN                            ! Don't write K1, K2 if param PBARLSHR = 'N'
         WRITE(ICARD(2),FMT1) K1     ;   CALL LEFT_ADJ_BDFLD ( ICARD(2) )
         WRITE(ICARD(3),FMT1) K2     ;   CALL LEFT_ADJ_BDFLD ( ICARD(3) )
      ELSE
         WRITE(ICARD(2),FMT0)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(2) )
         WRITE(ICARD(3),FMT0)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(3) )
      ENDIF
      WRITE(ICARD(4),FMT1) I12       ;   CALL LEFT_ADJ_BDFLD ( ICARD(4) )
      WRITE(ICARD(5),FMT0)
      WRITE(ICARD(6),FMT0)
      WRITE(ICARD(7),FMT0)
      WRITE(ICARD(8),FMT0)
      WRITE(ICARD(9),FMT0)
      IF ((PBARLSHR == 'Y') .OR. (DABS(I12) > 0)) THEN     ! Only write 2nd cont entry if we want K1, K2 included or I12 is /= 0.
         WRITE(F06,302) (ICARD(I),I=2,9)
      ENDIF

      WRITE(F06,101)

! **********************************************************************************************************************************
  101 FORMAT('$*******************************************************************************')

  102 FORMAT('$ PBAR equivalent for PBARL ',A)

  104 FORMAT('$--1---|-------2-------|-------3-------|-------4-------|-------5-------|-------6-------|-------7-------|-------8---',&
                 '----|-------9-------|--10---|')

  201 FORMAT(I8)

  210 FORMAT(F8.0)

  211 FORMAT(F8.1)

  212 FORMAT(F8.2)

  213 FORMAT(F8.3)

  214 FORMAT(F8.4)

  215 FORMAT(F8.5)

  216 FORMAT(F8.6)

  301 FORMAT('PBAR    ',7A8)

  302 FORMAT('        ',8A8)

      END SUBROUTINE WRITE_PBAR_EQUIV

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_BAR ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: A                 ! Short dimension
      REAL(DOUBLE)                    :: B                 ! Longer dimension
      REAL(DOUBLE)                    :: BETA              ! Ratio of dimensions

! **********************************************************************************************************************************
      JERR = 0

      A    = D(1)                                          ! D(1), D(2) /= 0 was checked when D(i) were input
      B    = D(2)

      AREA =  A*B
      I1   =  A*B*B*B/TWELVE
      I2   =  B*A*A*A/TWELVE
      I12  =  ZERO
      IF (A > B) THEN
         BETA = B/A
         JTOR =  A*B*B*B*THIRD*(ONE - 0.630D0*BETA*(ONE - BETA*BETA*BETA*BETA/TWELVE))
      ELSE
         BETA = A/B
         JTOR =  B*A*A*A*THIRD*(ONE - 0.630D0*BETA*(ONE - BETA*BETA*BETA*BETA/TWELVE))
      ENDIF

      K1   =  FIVE/SIX
      K2   =  FIVE/SIX

      Y(1) =  HALF*B   ;   Z(1) =  HALF*A
      Y(2) = -Y(1)     ;   Z(2) =  Z(1)
      Y(3) = -Y(1)     ;   Z(3) = -Z(1)
      Y(4) =  Y(1)     ;   Z(4) = -Z(1)

      END SUBROUTINE SECTION_PROPS_BAR

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_BOX ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: B1,H1             ! Outside dimensions of BOX
      REAL(DOUBLE)                    :: B2,H2             ! Inside  dimensions of BOX
      REAL(DOUBLE)                    :: ASTAR             ! Area enclosed by middle line of wall
      REAL(DOUBLE)                    :: S                 ! Length of middle wall
      REAL(DOUBLE)                    :: T                 ! Thickness of wall

! **********************************************************************************************************************************
      JERR = 0

      B1    = D(1)
      H1    = D(2)
      B2    = D(1) - TWO*D(4)
      H2    = D(2) - TWO*D(3)
      S     = B1 + B2 + H1 + H2
      T     = HALF*(D(3) + D(4))
      ASTAR = QUARTER*(B1 + B2)*(H1 + H2)

      IF (B2 <= ZERO) JERR = 1
      IF (H2 <= ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AREA =  B1*H1 - B2*H2
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      I1   =  (B1*H1*H1*H1 - B2*H2*H2*H2)/TWELVE
      I2   =  (H1*B1*B1*B1 - H2*B2*B2*B2)/TWELVE
      I12  =  ZERO
      JTOR =  FOUR*ASTAR*ASTAR*T/S
      K1   =  H2/D(2)
      K2   =  B2/D(1)
      Y(1) =  HALF*D(2)   ;   Z(1) =  HALF*D(1)
      Y(2) = -HALF*D(2)   ;   Z(1) =  HALF*D(1)
      Y(3) = -HALF*D(2)   ;   Z(3) = -HALF*D(1)
      Y(4) =  HALF*D(2)   ;   Z(4) = -HALF*D(1)

      END SUBROUTINE SECTION_PROPS_BOX

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_BOX1 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ASTAR             ! Area enclosed by middle line of wall
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: H                 ! Height of box
      REAL(DOUBLE)                    :: W                 ! Width  of box
      REAL(DOUBLE)                    :: YCG               ! Y location of C.G. relative to lower left corner
      REAL(DOUBLE)                    :: ZCG               ! Z location of C.G. relative to lower left corner
      REAL(DOUBLE)                    :: TL                ! Thickness of left   web
      REAL(DOUBLE)                    :: TR                ! Thickness of right  web
      REAL(DOUBLE)                    :: TB                ! Thickness of bottom web
      REAL(DOUBLE)                    :: TT                ! Thickness of top    web

! **********************************************************************************************************************************
      JERR = 0

      H  = D(2)
      W  = D(1)
      TL = D(6)
      TR = D(5)
      TB = D(4)
      TT = D(3)

      IF (W - TL - TR <= ZERO) JERR = 1
      IF (H - TT - TB <= ZERO) JERR = 1

      ASTAR = (W - HALF*(TL + TR))*(H - HALF*(TT + TB))
      DEN   = (H - HALF*(TT + TB))*(ONE/TL + ONE/TR) + (W - HALF*(TL + TR))*(ONE/TT + ONE/TB)
      IF (DEN > 0) THEN
         JTOR  =  FOUR*ASTAR*ASTAR/DEN
      ELSE
         JERR = 1
      ENDIF

      AREA = H*(TL + TR) + (W - TL - TR)*(TT + TB)
      IF (AREA <= ZERO) THEN
         JERR = 1
      ENDIF

      IF (JERR > 0) RETURN

      YCG  =  (H*(TL + TR)*HALF*H + (W - TL - TR)*(HALF*TB + (H - HALF*TT)))/AREA
      ZCG  =  (H*TL*HALF*TL + H*TR*(W - HALF*TR) + (W - TL -TR)*(TT + TB)*(TL + HALF*(W - TL - TR)))/AREA
                                                           ! Calc I1,I2 as I for each of parts about own C.G + part area times D^2
      I1   =  (TL + TR)*H*H*H/TWELVE + H*TL*(YCG - HALF*H)*(YCG - HALF*H)                                                          &
            + (W - TL - TR)*(TT*TT*TT + TB*TB*TB)/TWELVE + (W - TL -TR)*TB*(YCG - HALF*TB)*(YCG - HALF*TB)                         &
            +                                              (W - TL -TR)*TT*(H - HALF*TT - YCG)*(H - HALF*TT - YCG)
      I2   =  H*(TL*TL*TL + TR*TR*TR)/TWELVE + H*TL*(ZCG - HALF*TL)*(ZCG - HALF*TL) + H*TR*(W - HALF*TR -ZCG)*(W - HALF*TR -ZCG)   &
            + (TT + TB)*W*W*W/TWELVE + W*(TT + TB)*(HALF*W - ZCG)*(HALF*W - ZCG)
      I12  =  ZERO
      Y(1) =  H - YCG   ;   Z(1) =  W - ZCG
      Y(2) = -YCG       ;   Z(2) =  Z(1)
      Y(3) =  Y(2)      ;   Z(3) = -ZCG
      Y(4) =  Y(1)      ;   Z(4) =  Z(3)
      K1   =  (H - TT - TB)/H
      K2   =  (W - TL - TR)/W

      END SUBROUTINE SECTION_PROPS_BOX1

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_CHAN ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.12D0    ! Constant in calc of JTOR for open cross-sections
      REAL(DOUBLE)                    :: B                 ! Overall width
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: H                 ! Overall height
      REAL(DOUBLE)                    :: TW                ! Web thickness
      REAL(DOUBLE)                    :: TF                ! Flange thickness
      REAL(DOUBLE)                    :: ZSHR              ! Z location of shear center relative to center line of vertical web
      REAL(DOUBLE)                    :: ZCG               ! Z location of C.G. relative to center line of vertical web

! **********************************************************************************************************************************
      JERR = 0

      B  = D(1)
      H  = D(2)
      TW = D(3)
      TF = D(4)

      DEN  = SIX*(B - HALF*TW)*TF + (H - TF)*TW
      IF (DEN > ZERO) THEN
         ZSHR = -THREE*TF*(B - HALF*TW)*(B - HALF*TW)/DEN
         AREA =  (H - TWO*TF)*TW + TWO*B*TF
         ZCG = TWO*TF*(B - TW) * HALF*B / AREA
      ELSE
         JERR = 1
      ENDIF

      IF (H - TWO*TF < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      I1   =  TW*H*H*H/TWELVE + TWO*((B - TW)*TF*TF*TF/TWELVE + (B - TW)*TF*QUARTER*(H - TF)*(H - TF))

      I2   =  H*TW*TW*TW/TWELVE + H*TW*ZCG*ZCG                                                                                     &
            + TF*(B - TW)*(B - TW)*(B - TW)/TWELVE + TF*(B - TW)*(-ZCG + HALF*B)*(-ZCG + HALF*B)                                   &
            + TF*(B - TW)*(B - TW)*(B - TW)/TWELVE + TF*(B - TW)*(-ZCG + HALF*B)*(-ZCG + HALF*B)


      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(H*TW*TW*TW + (B - TW)*TF*TF*TF + (B - TW)*TF*TF*TF)
      Y(1) =  HALF*H   ;   Z(1) = -ZSHR - HALF*TW + B
      Y(2) = -Y(1)     ;   Z(2) =  Z(1)
      Y(3) = -Y(1)     ;   Z(3) = -ZSHR - HALF*TW
      Y(4) =  Y(1)     ;   Z(4) =  Z(3)
      K1   =  (H - TWO*TF)/H
      K2   =  (B - TW)/B

      END SUBROUTINE SECTION_PROPS_CHAN

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_CHAN1 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.12D0    ! Constant in calc of JTOR for open cross-sections
      REAL(DOUBLE)                    :: B                 ! Overall width
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: H                 ! Overall height
      REAL(DOUBLE)                    :: TW                ! Web thickness
      REAL(DOUBLE)                    :: TF                ! Flange thickness
      REAL(DOUBLE)                    :: ZSHR              ! Z location of shear center relative to lower left corner

! **********************************************************************************************************************************
      JERR = 0

      B  = D(1) + D(2)
      H  = D(4)
      TW = D(2)
      TF = HALF*(D(4) - D(3))

      DEN  = SIX*(B - HALF*TW)*TF + (H - TF)*TW
      IF (DEN > ZERO) THEN
         ZSHR = -THREE*TF*(B - HALF*TW)*(B - HALF*TW)/DEN
      ELSE
         JERR = 1
      ENDIF

      IF (TF < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AREA =  (H - TWO*TF)*TW + TWO*B*TF
      I1   =  TW*H*H*H/TWELVE + TWO*((B - TW)*TF*TF*TF/TWELVE + (B - TW)*TF*QUARTER*(H - TF)*(H - TF))

      I2   =  H*TW*TW*TW/TWELVE + H*TW*ZSHR*ZSHR                                                                                   &
            + TF*(B - TW)*(B - TW)*(B - TW) + TWO*TF*(B - TW)*(-ZSHR + HALF*(B + TW))*(-ZSHR + HALF*(B + TW))

      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(H*TW*TW*TW + (B - TF)*TF*TF*TF)
      Y(1) =  HALF*H   ;   Z(1) = -ZSHR - HALF*TW + B
      Y(2) = -Y(1)     ;   Z(2) =  Z(1)
      Y(3) = -Y(1)     ;   Z(3) = -ZSHR - HALF*TW
      Y(4) =  Y(1)     ;   Z(4) =  Z(3)
      K1   =  (H - TWO*TF)/H
      K2   =  (B - TW)/B

      END SUBROUTINE SECTION_PROPS_CHAN1

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_CHAN2 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.12D0    ! Constant in calc of JTOR for open cross-sections
      REAL(DOUBLE)                    :: B                 ! Overall width
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: H                 ! Overall height
      REAL(DOUBLE)                    :: TW                ! Web thickness
      REAL(DOUBLE)                    :: TF                ! Flange thickness
      REAL(DOUBLE)                    :: ZSHR              ! Z location of shear center relative to lower left corner

! **********************************************************************************************************************************
      JERR = 0

      B  = D(3)
      H  = D(4)
      TW = D(2)
      TF = D(1)

      DEN  = SIX*(B - HALF*TW)*TF + (H - TF)*TW
      IF (DEN > ZERO) THEN
         ZSHR = -THREE*TF*(B - HALF*TW)*(B - HALF*TW)/DEN
      ELSE
         JERR = 1
      ENDIF

      IF (TF < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AREA =  (H - TWO*TF)*TW + TWO*B*TF
      I2   =  H*TW*TW*TW/TWELVE + H*TW*ZSHR*ZSHR                                                                                   &
            + TF*(B - TW)*(B - TW)*(B - TW) + TWO*TF*(B - TW)*(-ZSHR + HALF*(B + TW))*(-ZSHR + HALF*(B + TW))

      I1   =  TW*H*H*H/TWELVE + TWO*((B - TW)*TF*TF*TF/TWELVE + (B - TW)*TF*QUARTER*(H - TF)*(H - TF))

      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(H*TW*TW*TW + (B - TF)*TF*TF*TF)
      Y(1) =  HALF*H   ;   Z(1) = -ZSHR - HALF*TW + B
      Y(2) = -Y(1)     ;   Z(2) =  Z(1)
      Y(3) = -Y(1)     ;   Z(3) = -ZSHR - HALF*TW
      Y(4) =  Y(1)     ;   Z(4) =  Z(3)
      K1   =  (H - TWO*TF)/H
      K2   =  (B - TW)/B

      END SUBROUTINE SECTION_PROPS_CHAN2

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_CROSS ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.17 on #10 in table 2.5 but MSC uses 1.0

      REAL(DOUBLE)                    :: AVERT             ! Area of vertical member
      REAL(DOUBLE)                    :: AHORZ             ! Area of horizontal member
      REAL(DOUBLE)                    :: HV                ! Height of vert member
      REAL(DOUBLE)                    :: WH                ! Width of horiz member
      REAL(DOUBLE)                    :: TV                ! Thickness of vertical member
      REAL(DOUBLE)                    :: TH                ! Thickness of horizontal member
      REAL(DOUBLE)                    :: W0                ! Width of left and righ

! **********************************************************************************************************************************
      JERR = 0

      HV = D(3)
      W0 = HALF*D(1)
      WH = D(1) + D(2)
      TV = D(2)
      TH = D(4)

      AVERT = HV*TV
      AHORZ = WH*TH
      AREA = AVERT + AHORZ - TV*TH
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      I1   =  TV*HV*HV*HV/TWELVE + TWO*W0*TH*TH*TH/TWELVE
      I2   =  HV*TV*TV*TV/TWELVE + TWO*(TH*W0*W0*W0/TWELVE + TH*W0*QUARTER*(TV + W0)*(TV+W0))
      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(HV*TV*TV*TV + TWO*W0*TH*TH*TH)
      Y(1) =  HALF*HV  ;   Z(1) =  ZERO
      Y(2) =  ZERO     ;   Z(2) =  HALF*TV + W0
      Y(3) = -Y(1)     ;   Z(3) =  ZERO
      Y(4) =  ZERO     ;   Z(4) = -Z(2)

      K1   =  (FIVE/SIX)*AVERT/AREA
      K2   =  (FIVE/SIX)*AHORZ/AREA

      END SUBROUTINE SECTION_PROPS_CROSS

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_H ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.31, MSC uses 1.0

      REAL(DOUBLE)                    :: AVERT             ! Area of vertical member
      REAL(DOUBLE)                    :: AHORZ             ! Area of horizontal member
      REAL(DOUBLE)                    :: HF                ! Height
      REAL(DOUBLE)                    :: WW                ! Width of center member
      REAL(DOUBLE)                    :: TF                ! Thickness of vertical member
      REAL(DOUBLE)                    :: TW                ! Thickness of center member

! **********************************************************************************************************************************
      JERR = 0

      HF = D(3)
      WW = D(1)
      TF = HALF*D(2)
      TW = D(4)

      AVERT = HF*TF
      AHORZ = WW*TW
      AREA =  TWO*AVERT + AHORZ

      I1   =  TWO*TF*HF*HF*HF/TWELVE + WW*TW*TW*TW/TWELVE
      I2   =  TWO*(HF*TF*TF*TF/TWELVE + HF*TF*QUARTER*(WW + TF)*(WW + TF)) + TW*WW*WW*WW/TWELVE
      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(TWO*HF*TF*TF*TF + WW*TW*TW*TW)
      Y(1) =  HALF*HF  ;   Z(1) =  HALF*WW + TF
      Y(2) = -Y(1)     ;   Z(2) =  Z(1)
      Y(3) =  Y(2)     ;   Z(3) = -Z(1)
      Y(4) =  Y(1)     ;   Z(4) =  Z(3)

      K1 = (FIVE/SIX)*TWO*AVERT/AREA
      K2 = AHORZ/AREA

      END SUBROUTINE SECTION_PROPS_H

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_HAT ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.15D0    ! Constant in calc of JTOR for open cross-sections
!                                                            Avg: Pilkey has 1.12 for channel section and 1.17 for I section
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: NUM               ! Intermediate variable
      REAL(DOUBLE)                    :: H                 ! Height
      REAL(DOUBLE)                    :: T                 ! Thickness of hat walls
      REAL(DOUBLE)                    :: W1                ! Overall width of top of hat
      REAL(DOUBLE)                    :: W2                ! Width of base legs of hat
      REAL(DOUBLE)                    :: W3                ! Inside width of top of hat
      REAL(DOUBLE)                    :: YSHR              ! Location of shear center relative to center of top of hat

! **********************************************************************************************************************************
      JERR = 0

      H  = D(1)
      T  = D(2)
      W1 = D(3)
      W2 = D(4)
      W3 = W1 - TWO*T

      DEN  = W3*W3*W3 + SIX*W3*W3*(H - T) + SIX*(W2 + T)*W3*W3 + EIGHT*(W2 + T)*(W2 + T)*(W2 + T) + TWELVE*(W2 + T)*(W2 + T)*W3
      NUM  = (H - T)*(THREE*(H - T)*W3*W3 - EIGHT*W3*W3*W3)
      IF (DEN > ZERO) THEN
         YSHR = NUM/DEN
      ELSE
         JERR = 1
      ENDIF

      IF (H - T      < ZERO) JERR = 1
      IF (W1 - TWO*T < ZERO) JERR = 1

      AREA =  T*(TWO*W2 + TWO*H + (W1 - TWO*T))
      IF (AREA <= ZERO) THEN
         JERR = 1
      ENDIF

      IF (JERR > 0) RETURN

      I1   =  TWO*W2*T*T*T/TWELVE + TWO*W2*T*(H - T + YSHR)*T*(H - T + YSHR)                                                       &
            + TWO*T*H*H*H/TWELVE + TWO*H*T*(HALF*H + YSHR - HALF*T)*(HALF*H + YSHR - HALF*T)                                       &
            + W3*T*T*T/TWELVE + W3*T*(HALF*T + YSHR)*(HALF*T + YSHR)

      I2   =  TWO*T*W2*W2*W2/TWELVE + W2*T*QUARTER*(W1 + W2)*(W1 + W2)                                                             &
            + TWO*H*T*T*T/TWELVE + TWO*H*T*(W3 + HALF*T)*(W3 + HALF*T)                                                             &
            + T*(W3*W3*W3)/TWELVE

      I12  =  ZERO
      JTOR =  ALPHA*THIRD*( TWO*W2 + TWO*H + (W1 - TWO*T) )*T*T*T
      Y(1) =  HALF*T + YSHR       ; Z(1) =  HALF*W1
      Y(2) = -YSHR + HALF*T - H   ; Z(2) =  HALF*W1 + W2
      Y(3) =  Y(2)                ; Z(3) = -Z(2)
      Y(4) =  Y(1)                ; Z(4) = -Z(1)
      K1   =  (H - TWO*T)/H
      K2   =  (W1 + W3)/(W1 + W2)

      END SUBROUTINE SECTION_PROPS_HAT

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_HEXA ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: BETA              ! Ratio of dimensions
      REAL(DOUBLE)                    :: H0                ! Overall height
      REAL(DOUBLE)                    :: HT                ! Height of triangular portion
      REAL(DOUBLE)                    :: W0                ! Overall width
      REAL(DOUBLE)                    :: W1                ! Width of base
      REAL(DOUBLE)                    :: WAVG              ! Average width
      REAL(DOUBLE)                    :: WT                ! Width of triangular portion

! **********************************************************************************************************************************
      JERR = 0

      H0   = D(3)
      W0   = D(2)
      WT   = D(1)
      W1   = W0 -TWO*WT
      WAVG = HALF*(W0 + W1)
      HT   = HALF*H0

      AREA =  W1*H0 + TWO*WT*HT
      I1   =  W1*H0*H0*H0/TWELVE + FOUR*( WT*HT*HT*HT/(THREE*TWELVE) + HALF*WT*HT*HT*HT/NINE )
      I2   =  H0*W1*W1*W1/TWELVE + FOUR*( HT*WT*WT*WT/(THREE*TWELVE) + HALF*WT*HT*WT*WT/NINE )
      I12  =  ZERO
      IF (WAVG > H0) THEN
         BETA = H0/WAVG
         JTOR =  WAVG*H0*H0*H0*THIRD*(ONE - 0.630D0*BETA*(ONE - BETA*BETA*BETA*BETA/TWELVE))
      ELSE
         BETA = WAVG/H0
         JTOR =  H0*WAVG*WAVG*WAVG*THIRD*(ONE - 0.630D0*BETA*(ONE - BETA*BETA*BETA*BETA/TWELVE))
      ENDIF
      Y(1) =  HT     ;   Z(1) =  ZERO
      Y(2) = -HT     ;   Z(2) =  ZERO
      Y(3) =  ZERO   ;   Z(3) =  HALF*W0
      Y(4) =  ZERO   ;   Z(4) = -Z(3)
      K1   =  FIVE/SIX
      K2   =  FIVE/SIX

      END SUBROUTINE SECTION_PROPS_HEXA

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_I ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.31 in #10 and has 1.0 in #12 on Table 2.5. MSC uses 1.0

      REAL(DOUBLE)                    :: ATOP              ! Area of top flange
      REAL(DOUBLE)                    :: ABOT              ! Area of bottom flange
      REAL(DOUBLE)                    :: AWEB              ! Area of web
      REAL(DOUBLE)                    :: H0                ! Height of I beam
      REAL(DOUBLE)                    :: HW                ! Height of web
      REAL(DOUBLE)                    :: HM                ! Height between C/L of top and bottom
      REAL(DOUBLE)                    :: TB                ! Thickness of top flange
      REAL(DOUBLE)                    :: TT                ! Thickness of bottom flange
      REAL(DOUBLE)                    :: TW                ! Thickness of web
      REAL(DOUBLE)                    :: WB                ! Width of top flange
      REAL(DOUBLE)                    :: WT                ! Width of bottom flange
      REAL(DOUBLE)                    :: YCG               ! Dist from bottom of I section to section C.G.

! **********************************************************************************************************************************
      JERR = 0

      H0 = D(1)
      HW = H0 - D(5) - D(6)
      TB = D(5)
      TT = D(6)
      TW = D(4)
      WB = D(2)
      WT = D(3)
      HM = H0 - HALF*(TT + TB)

      IF (HW < ZERO) JERR = 1
      IF ((WB - TW  < ZERO) .OR. (WT - TW < ZERO)) JERR = 1
      IF (JERR > 0) RETURN

      ATOP = WT*TT
      ABOT = WB*TB
      AWEB = HW*TW
      AREA =  ATOP + ABOT + AWEB

      YCG  =  (WB*TB*(HALF*TB) + HW*TW*(TB + HALF*HW) + WT*TT*(H0 - HALF*TT))/AREA
                                                           ! Calc I1 as I for each of parts about own C.G + part area times D^2
      I1   =  WT*TT*TT*TT/TWELVE + WT*TT*(H0 - HALF*TT - YCG)*(H0 - HALF*TT - YCG)                                                 &
            + WB*TB*TB*TB/TWELVE + WB*TB*(YCG  - HALF*TB)*(YCG  - HALF*TB)                                                         &
            + TW*HW*HW*HW/TWELVE + TW*HW*(TB + HALF*HW - YCG)*(TB + HALF*HW - YCG)
                                                           ! Due to sym about Y only need B*H^3/12 for 2 flanges + web
      I2   =  (TT*WT*WT*WT + HW*TW*TW*TW + TB*WB*WB*WB)/TWELVE

      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(WT*TT*TT*TT + WB*TB*TB*TB + HM*TW*TW*TW)
      Y(1) =  H0 - YCG     ;   Z(1) =  HALF*WT
      Y(2) = -YCG          ;   Z(2) =  HALF*WB
      Y(3) =  Y(2)         ;   Z(3) = -Z(2)
      Y(4) =  Y(1)         ;   Z(4) = -Z(1)

      K1   =  AWEB/AREA                                    ! Assume web carries shear uniformly in plane 1 shear
      K2   =  (FIVE/SIX)*(ATOP + ABOT)/AREA                ! Assume top & bot flanges carry all shear parabolically in plane 2

      END SUBROUTINE SECTION_PROPS_I

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_I1 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.30 in #10 and 1.0 in #12 in table 2.5

      REAL(DOUBLE)                    :: ATOP              ! Area of top flange
      REAL(DOUBLE)                    :: ABOT              ! Area of bottom flange
      REAL(DOUBLE)                    :: AWEB              ! Area of web
      REAL(DOUBLE)                    :: H0                ! Height of I beam
      REAL(DOUBLE)                    :: HW                ! Height of web
      REAL(DOUBLE)                    :: TB                ! Thickness of top flange
      REAL(DOUBLE)                    :: TT                ! Thickness of bottom flange
      REAL(DOUBLE)                    :: TW                ! Thickness of web
      REAL(DOUBLE)                    :: WB                ! Width of top flange
      REAL(DOUBLE)                    :: WT                ! Width of bottom flange
      REAL(DOUBLE)                    :: YCG               ! Dist from bottom of I section to section C.G.

! **********************************************************************************************************************************
      JERR = 0

      H0 = D(4)
      HW = D(3)
      TB = HALF*(D(4) - D(3))
      TT = TB
      TW = D(2)
      WB = D(1) + D(2)
      WT = WB

      IF (HW < ZERO) JERR = 1
      IF ((WB - TW  < ZERO) .OR. (WT - TW < ZERO)) JERR = 1
      IF (JERR > 0) RETURN

      ATOP = WT*TT
      ABOT = WB*TB
      AWEB = HW*TW
      AREA =  ATOP + ABOT + AWEB

      YCG  =  ZERO
                                                           ! Calc I1 as I for each of parts about own C.G + part area times D^2
      I1   =  WT*TT*TT*TT/TWELVE + WT*TT*(HALF*HW + HALF*TT)*(HALF*HW + HALF*TT)                                                   &
            + WB*TB*TB*TB/TWELVE + WB*TB*(HALF*HW + HALF*TB)*(HALF*HW + HALF*TB)                                                   &
            + TW*HW*HW*HW/TWELVE
                                                           ! Due to sym about Y only need B*H^3/12 for 2 flanges + web
      I2   =  (TT*WT*WT*WT + HW*TW*TW*TW + TB*WB*WB*WB)/TWELVE

      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(WB*TB*TB*TB + HW*TW*TW*TW + WT*TT*TT*TT)
      Y(1) =  H0 - YCG     ;   Z(1) =  HALF*WT
      Y(2) = -YCG          ;   Z(2) =  HALF*WB
      Y(3) =  Y(2)         ;   Z(3) = -Z(2)
      Y(4) =  Y(1)         ;   Z(4) = -Z(1)

      K1   =  AWEB/AREA                                    ! Assume web carries shear uniformly in plane 1 shear
      K2   =  (FIVE/SIX)*(ATOP + ABOT)/AREA                ! Assume top & bot flanges carry all shear parabolically in plane 2

      END SUBROUTINE SECTION_PROPS_I1

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_ROD ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: RAD

! **********************************************************************************************************************************
      JERR = 0

      RAD = D(1)

      AREA =  PI*RAD*RAD
      I1   =  QUARTER*PI*RAD*RAD*RAD*RAD
      I2   =  I1
      I12  =  ZERO
      JTOR =  I1 + I2
      Y(1) =  RAD    ;   Z(1) =  ZERO
      Y(2) =  ZERO   ;   Z(2) =  RAD
      Y(3) = -RAD    ;   Z(3) =  ZERO
      Y(4) =  ZERO   ;   Z(4) = -RAD

      K1   =  SIX/SEVEN
      K2   =  SIX/SEVEN

      END SUBROUTINE SECTION_PROPS_ROD

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_T ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.12 in #10 on table 2.5. MSC uses 1.0

      REAL(DOUBLE)                    :: AVERT             ! Area of vertical member
      REAL(DOUBLE)                    :: AHORZ             ! Area of horizontal member
      REAL(DOUBLE)                    :: JCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: KCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: NUM               ! Intermediate variable
      REAL(DOUBLE)                    :: WT                ! Width of top flange
      REAL(DOUBLE)                    :: TT                ! Thickness of top flange
      REAL(DOUBLE)                    :: HW                ! Height of web
      REAL(DOUBLE)                    :: TW                ! Thickness of web
      REAL(DOUBLE)                    :: YCG               ! Dist from bottom of web to C.G.


! **********************************************************************************************************************************
      JERR = 0

      WT = D(1)
      TT = D(3)
      HW = D(2) - D(3)
      TW = D(4)

      IF (HW < ZERO) JERR = 1
      IF (ABS(WT - TW) < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AVERT = HW*TW
      AHORZ = WT*TT
      AREA = AVERT + AHORZ
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      YCG  =  (WT*TT*(HW + HALF*TT) + HW*TW*(HALF*HW))/AREA
                                                           ! Calc I1 as I for each of parts about own C.G + part area times D^2
      I1   =  WT*TT*TT*TT/TWELVE + WT*TT*(HW + HALF*TT - YCG)*(HW + HALF*TT - YCG)                                                &
            + TW*HW*HW*HW/TWELVE + HW*TW*(YCG - HALF*HW)*(YCG - HALF*HW)

      I2   =  (TT*WT*WT*WT + HW*TW*TW*TW)/TWELVE
      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(WT*TT*TT*TT + HW*TW*TW*TW)
      Y(1) =  HW + TT - YCG   ;   Z(1) =  ZERO
      Y(2) =  Y(1)            ;   Z(2) =  HALF*WT
      Y(3) = -YCG             ;   Z(3) =  ZERO
      Y(4) =  Y(1)            ;   Z(4) = -Z(2)

      JCONS = WT*TT/((HW + HALF*TT)*TW)
      KCONS = WT/(HW + HALF*TT)
      NUM   = TEN*(ONE + FOUR*JCONS)*(ONE + FOUR*JCONS)
      DEN   = (12.D0 + 96.D0*JCONS + 276.D0*JCONS*JCONS + 192.D0*JCONS*JCONS*JCONS) + 30.D0*KCONS*KCONS*JCONS*(ONE + JCONS)
      IF (DABS(DEN) > ZERO) THEN
         K1 = NUM/DEN
      ELSE
         K1 = ZERO
      ENDIF
      K2    = (FIVE/SIX)*AHORZ/AREA

      END SUBROUTINE SECTION_PROPS_T

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_T1 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.12 in #10 on table 2.5. MSC uses 1.0

      REAL(DOUBLE)                    :: AVERT             ! Area of vertical member
      REAL(DOUBLE)                    :: AHORZ             ! Area of horizontal member
      REAL(DOUBLE)                    :: JCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: KCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: NUM               ! Intermediate variable
      REAL(DOUBLE)                    :: HV                ! Height of vert flange
      REAL(DOUBLE)                    :: TV                ! Thickness of vert flange
      REAL(DOUBLE)                    :: WH                ! Width of horiz flange
      REAL(DOUBLE)                    :: TH                ! Thickness of horiz flange
      REAL(DOUBLE)                    :: ZCG               ! Dist from bottom of web to C.G.

! **********************************************************************************************************************************
      JERR = 0

      HV = D(1)
      TV = D(3)
      WH = D(2)
      TH = D(4)

      AVERT = HV*TV
      AHORZ = WH*TH
      AREA = AHORZ + AVERT
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      ZCG =  (HV*TV*(WH + HALF*TV) + WH*TH*HALF*WH)/AREA

      I1   =  WH*TH*TH*TH/TWELVE + TV*HV*HV*HV/TWELVE
                                                           ! Calc I2 as I for each of parts about own C.G + part area times D^2
      I2   =  TH*WH*WH*WH/TWELVE + WH*TH*(ZCG - HALF*WH)*(ZCG - HALF*WH)                                                           &
            + HV*TV*TV*TV/TWELVE + HV*TV*(ZCG - WH - HALF*TV)*(ZCG - WH - HALF*TV)
      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(WH*TH*TH*TH + HV*TV*TV*TV)
      Y(1) =  ZERO      ;   Z(1) =  (WH + TV - ZCG)
      Y(2) = -HALF*HV   ;   Z(2) =  Z(1)
      Y(3) =  ZERO      ;   Z(3) = -ZCG
      Y(4) = -Y(2)      ;   Z(4) =  Z(1)

      JCONS = HV*TV/((WH + HALF*TV)*TH)
      KCONS = HV/(WH + HALF*TV)
      NUM   = TEN*(ONE + FOUR*JCONS)*(ONE + FOUR*JCONS)
      DEN   = (12.D0 + 96.D0*JCONS + 276.D0*JCONS*JCONS + 192.D0*JCONS*JCONS*JCONS) + 30.D0*KCONS*KCONS*JCONS*(ONE + JCONS)
      K1    = (FIVE/SIX)*AVERT/AREA
      IF (DABS(DEN) > ZERO) THEN
         K2 = NUM/DEN
      ELSE
         K2 = ZERO
      ENDIF

      END SUBROUTINE SECTION_PROPS_T1

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_T2 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.12 in #10 on table 2.5. MSC uses 1.0

      REAL(DOUBLE)                    :: AVERT             ! Area of vertical member
      REAL(DOUBLE)                    :: AHORZ             ! Area of horizontal member
      REAL(DOUBLE)                    :: JCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: KCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: DEN               ! Intermediate variable
      REAL(DOUBLE)                    :: NUM               ! Intermediate variable
      REAL(DOUBLE)                    :: WB                ! Width of bottom flange
      REAL(DOUBLE)                    :: TB                ! Thickness of bottom flange
      REAL(DOUBLE)                    :: HW                ! Height of web
      REAL(DOUBLE)                    :: TW                ! Thickness of web
      REAL(DOUBLE)                    :: YCG               ! Dist from bottom of I section to section C.G.

! **********************************************************************************************************************************
      JERR = 0

      WB = D(1)
      TB = D(3)
      HW = D(2) - D(3)
      TW = D(4)

      IF (HW <= ZERO) JERR = 1
      IF (ABS(WB - TW) < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AVERT = HW*TW
      AHORZ = WB*TB
      AREA = AVERT + AHORZ
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      YCG  =  (TB*WB*HALF*TB + HW*TW*(TB + HALF*HW))/AREA
                                                           ! Calc I1 as I for each of parts about own C.G + part area times D^2
      I1   =  WB*TB*TB*TB/TWELVE + WB*TB*(YCG - HALF*TB)*(YCG - HALF*TB)                                                           &
            + TW*HW*HW*HW/TWELVE + TW*HW*(TB + HALF*HW - YCG)*(TB + HALF*HW - YCG)
      I2   =  (TB*WB*WB*WB + HW*TW*TW*TW)/TWELVE
      I12  =  ZERO
      JTOR =  ALPHA*THIRD*(WB*TB*TB*TB + HW*TW*TW*TW)
      JTOR =  ALPHA*THIRD*((D(1)-D(4))*D(3)*D(3)*D(3) + D(2)*D(4)*D(4)*D(4))
      Y(1) =  TB + HW - YCG   ;   Z(1) =  HALF*TW
      Y(2) =  -YCG            ;   Z(2) =  HALF*WB
      Y(3) =  Y(2)            ;   Z(3) = -Z(2)
      Y(4) =  Y(1)            ;   Z(4) = -Z(1)

      JCONS = WB*TB/((HW + HALF*TB)*TW)
      KCONS = WB/(HW + HALF*TB)
      NUM   = TEN*(ONE + FOUR*JCONS)*(ONE + FOUR*JCONS)
      DEN   = (12.D0 + 96.D0*JCONS + 276.D0*JCONS*JCONS + 192.D0*JCONS*JCONS*JCONS) + 30.D0*KCONS*KCONS*JCONS*(ONE + JCONS)
      IF (DABS(DEN) > ZERO) THEN
         K1 = NUM/DEN
      ELSE
         K1 = ZERO
      ENDIF
      K2    = (FIVE/SIX)*AHORZ/AREA

      END SUBROUTINE SECTION_PROPS_T2

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_TUBE ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: JCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: NU                ! Poisson ratio
      REAL(DOUBLE)                    :: RAD1,RAD2         ! Radii
      REAL(DOUBLE)                    :: RAD21             ! RAD2/RAD1

! **********************************************************************************************************************************
      JERR = 0

      RAD1 = D(1)
      RAD2 = D(2)

      IF (RAD1 - RAD2 < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AREA =  PI*(RAD1*RAD1 - RAD2*RAD2)
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      I1   =  QUARTER*PI*(RAD1*RAD1*RAD1*RAD1 - RAD2*RAD2*RAD2*RAD2)
      I2   =  I1
      I12  =  ZERO
      JTOR =  TWO*I1
      Y(1) =  RAD1   ;   Z(1) =  ZERO
      Y(2) =  ZERO   ;   Z(2) =  RAD1
      Y(3) = -Y(1)   ;   Z(3) =  ZERO
      Y(4) =  ZERO   ;   Z(4) = -Z(2)

      NU   =  ZERO
      RAD21=  RAD2/RAD1
      JCONS=  (1 + RAD21*RAD21)
      K1   =  SIX*(1 + NU)*JCONS*JCONS/((SEVEN + SIX*NU)*JCONS*JCONS + TWO*(TEN + SIX*NU)*RAD21*RAD21)
      K2   =  K1

      END SUBROUTINE SECTION_PROPS_TUBE

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_TUBE2 ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: JCONS             ! Intermediate variable
      REAL(DOUBLE)                    :: NU                ! Poisson ratio
      REAL(DOUBLE)                    :: RAD1,RAD2         ! Radii
      REAL(DOUBLE)                    :: T                 ! Tube thickness
      REAL(DOUBLE)                    :: RAD21             ! RAD2/RAD1

! **********************************************************************************************************************************
      JERR = 0

      RAD1 = D(1)
      T = D(2)
      RAD2 = RAD1 - T

      IF (T < ZERO) JERR = 1
      IF (JERR > 0) RETURN

      AREA =  PI*(RAD1*RAD1 - RAD2*RAD2)
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      I1   =  QUARTER*PI*(RAD1*RAD1*RAD1*RAD1 - RAD2*RAD2*RAD2*RAD2)
      I2   =  I1
      I12  =  ZERO
      JTOR =  TWO*I1
      Y(1) =  RAD1   ;   Z(1) =  ZERO
      Y(2) =  ZERO   ;   Z(2) =  RAD1
      Y(3) = -Y(1)   ;   Z(3) =  ZERO
      Y(4) =  ZERO   ;   Z(4) = -Z(2)

      NU   =  ZERO
      RAD21=  RAD2/RAD1
      JCONS=  (1 + RAD21*RAD21)
      K1   =  SIX*(1 + NU)*JCONS*JCONS/((SEVEN + SIX*NU)*JCONS*JCONS + TWO*(TEN + SIX*NU)*RAD21*RAD21)
      K2   =  K1

      END SUBROUTINE SECTION_PROPS_TUBE2

! ##################################################################################################################################

      SUBROUTINE SECTION_PROPS_Z ( JERR )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: JERR              ! Error indicator

      REAL(DOUBLE)                    :: ALPHA = 1.0D0     ! Constant in calc of JTOR for open cross-sections
!                                                            Pilkey has 1.12 in #10 on table 2.5. MSC uses 1.0

      REAL(DOUBLE)                    :: ATOP              ! Area of top flange
      REAL(DOUBLE)                    :: ABOT              ! Area of bottom flange
      REAL(DOUBLE)                    :: AWEB              ! Area of web
      REAL(DOUBLE)                    :: WF                ! Width of bottom flange
      REAL(DOUBLE)                    :: TF                ! Thickness of bottom flange
      REAL(DOUBLE)                    :: HW                ! Height of web
      REAL(DOUBLE)                    :: TW                ! Thickness of web

! **********************************************************************************************************************************
      JERR = 0

      WF = D(1) + D(2)
      TF = HALF*(D(4) - D(3))
      HW = D(3)
      TW = D(2)

      IF (HW <= 0) JERR = 1
      IF (JERR > 0) RETURN

      ATOP = WF*TF
      ABOT = ATOP
      AWEB = HW*TW

      AREA =  ATOP + ABOT + AWEB
      IF (AREA <= ZERO) THEN
         JERR = 1
         RETURN
      ENDIF

      I1   =  TW*HW*HW*HW/TWELVE + TWO*( WF*TF*TF*TF/TWELVE + WF*TF*QUARTER*(HW + TF)*(HW + TF) )
      I2   =  HW*TW*TW*TW/TWELVE + TWO*( TF*WF*WF*WF/TWELVE + TF*WF*QUARTER*(WF - TW)*(WF - TW) )
      I12  =  -TF*D(1)*(D(1) + D(2))*(D(3) + TF)/TWO
      JTOR =  ALPHA*THIRD*(TWO*WF*TF*TF*TF + HW*TW*TW*TW)
      Y(1) =  HALF*HW + TF     ;   Z(1) =  HALF*TW
      Y(2) = -Y(1)             ;   Z(2) =  (WF - HALF*TW)
      Y(3) = -Y(1)             ;   Z(3) = -Z(1)
      Y(4) =  Y(1)             ;   Z(4) = -Z(2)
      K1   =  AWEB/AREA
      K2   =  (FIVE/SIX)*(ATOP + ABOT)/AREA


      END SUBROUTINE SECTION_PROPS_Z

      END SUBROUTINE BD_PBARL


      SUBROUTINE BD_PBEAM ( CARD, LARGE_FLD_INP )

! Processes PBEAM Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  EPSIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BEAMTOR, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBEAM, NPBEAM
      USE CONSTANTS_1, ONLY           :  ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBEAM, RPBEAM
      USE PARAMS, ONLY                :  SUPINFO

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, R8FLD
      USE BAR_PROPERTY_VALIDATION, ONLY:  CHECK_BAR_MOIS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   =   'BD_PBEAM'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHRINP             ! Character field read from CARD
      CHARACTER(LEN(JCARD))           :: ID                 ! Property ID for this PBEAM
      CHARACTER(LEN(JCARD))           :: NAME               ! Char name for output error purposes

      INTEGER(LONG)                   :: ICONT       = 0    ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0    ! Error indicator
      INTEGER(LONG)                   :: J                  ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID = 0    ! Material ID (field 3 of this property card)
      INTEGER(LONG)                   :: PROPERTY_ID = 0    ! Property ID (field 2 of this property card)


      REAL(DOUBLE)                    :: AREA_A      = ZERO ! Cross sectional area at end A
      REAL(DOUBLE)                    :: I1_A        = ZERO ! Moment of inertia, plane 1 at end A
      REAL(DOUBLE)                    :: I2_A        = ZERO ! Moment of inertia, plane 2 at end A
      REAL(DOUBLE)                    :: I12_A       = ZERO ! Product of inertia at end A
      REAL(DOUBLE)                    :: JTOR_A      = ZERO ! Torsional constantr at end A
      REAL(DOUBLE)                    :: NSM_A       = ZERO ! Nonstructural mass at end A

      REAL(DOUBLE)                    :: AREA        = ZERO ! Cross sectional area at any location along beam
      REAL(DOUBLE)                    :: I1          = ZERO ! Moment of inertia, plane 1 at any location along beam
      REAL(DOUBLE)                    :: I2          = ZERO ! Moment of inertia, plane 2 at any location along beam
      REAL(DOUBLE)                    :: I12         = ZERO ! Product of inertia at any location along beam
      REAL(DOUBLE)                    :: JTOR        = ZERO ! Torsional constantr at any location along beam
      REAL(DOUBLE)                    :: NSM         = ZERO ! Nonstructural mass at any location along beam



! **********************************************************************************************************************************
! PBEAM Bulk Data Card routine

!  FIELD          ITEM                                                ARRAY ELEMENT
!  -----          ----                                             ------------------
!    2      Property ID                                 PID          PBEAM(npbeam, 1)
!    3      Material ID                                 MID          PBEAM(npbeam, 2)
!    4      Cross sectional area          : end A       A(A)        RPBEAM(npbeam, 1)
!    5      Area moment of inertia 1      :   "         I1(A)       RPBEAM(npbeam, 2)
!    6      Area moment of inertia 2      :   "         I2(A)       RPBEAM(npbeam, 3)
!    7      Area product of inertia 12    :   "         I12(A)      RPBEAM(npbeam, 4)
!    8      Torsion constant              :   "         J(A)        RPBEAM(npbeam, 5)
!    9      Non structural mass           :   "         NSM(A)      RPBEAM(npbeam, 6)

! on mandatory 2nd card:
!    2      Stress coefficient            :   "         C1(A)       RPBEAM(npbeam, 7)
!    3      Stress coefficient            :   "         C2(A)       RPBEAM(npbeam, 8)
!    4      Stress coefficient            :   "         D1(A)       RPBEAM(npbeam, 9)
!    5      Stress coefficient            :   "         D2(A)       RPBEAM(npbeam,10)
!    6      Stress coefficient            :   "         E1(A)       RPBEAM(npbeam,11)
!    7      Stress coefficient            :   "         E2(A)       RPBEAM(npbeam,12)
!    8      Stress coefficient            :   "         F1(A)       RPBEAM(npbeam,13)
!    9      Stress coefficient            :   "         F2(A)       RPBEAM(npbeam,14)

! on mandatory 3rd card:
!    2      Stress output request option                S0           PBEAM(npbeam, 3)
!    3      Loc of next set of data (X/XB)              X/XB        RPBEAM(NPBEAM,15)
!    4      Cross sectional area          : end B       A(B)        RPBEAM(npbeam,16)
!    5      Area moment of inertia 1      :   "         I1(B)       RPBEAM(npbeam,17)
!    6      Area moment of inertia 2      :   "         I2(B)       RPBEAM(npbeam,18)
!    7      Area produc of inertia 12     :   "         I12(B)      RPBEAM(npbeam,19)
!    8      Torsion constant              :   "         J(B)        RPBEAM(npbeam,20)
!    9      Non structural mass           :   "         NSM(B)      RPBEAM(npbeam,21)

! on optional  4th card:
!    2      Stress coefficient            :   "         C1(B)       RPBEAM(npbeam,22)
!    3      Stress coefficient            :   "         C2(B)       RPBEAM(npbeam,23)
!    4      Stress coefficient            :   "         D1(B)       RPBEAM(npbeam,24)
!    5      Stress coefficient            :   "         D2(B)       RPBEAM(npbeam,25)
!    6      Stress coefficient            :   "         E1(B)       RPBEAM(npbeam,26)
!    7      Stress coefficient            :   "         E2(B)       RPBEAM(npbeam,27)
!    8      Stress coefficient            :   "         F1(B)       RPBEAM(npbeam,28)
!    9      Stress coefficient            :   "         F2(B)       RPBEAM(npbeam,29)

! on optional  5th card:
!    2      Shear factor for plane 1                    K1          RPBEAM(npbeam,30)
!    3      Shear factor for plane 2                    K2          RPBEAM(npbeam,31)
!    4      Shear relief coeff due to taper for plane 1 S1          RPBEAM(npbeam,32)
!    5      Shear relief coeff due to taper for plane 2 S2          RPBEAM(npbeam,33)
!    6      NSM MOI/length about NSM C.G. at end A      NSI(A)      RPBEAM(npbeam,34)
!    7      NSM MOI/length about NSM C.G. at end B      NSI(B)      RPBEAM(npbeam,35)
!    8      Warping coefficient for end A               CW(A)       RPBEAM(npbeam,36)
!    9      Warping coefficient for end B               CW(B)       RPBEAM(npbeam,37)

! on optional  6th card:
!    2      y coord of C.G. of NSM at end A             M1(A)       RPBEAM(npbeam,38)
!    3      z coord of C.G. of NSM at end A             M2(A)       RPBEAM(npbeam,39)
!    4      y coord of C.G. of NSM at end B             M1(B)       RPBEAM(npbeam,40)
!    5      z coord of C.G. of NSM at end B             M2(B)       RPBEAM(npbeam,41)
!    6      y coord of neutral axis for end A           N1(A)       RPBEAM(npbeam,42)
!    7      z coord of neutral axis for end A           N2(A)       RPBEAM(npbeam,43)
!    8      y coord of neutral axis for end B           N1(B)       RPBEAM(npbeam,44)
!    9      z coord of neutral axis for end B           N2(B)       RPBEAM(npbeam,45)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Increment NPBEAM

      NPBEAM = NPBEAM + 1

! Read and check data on parent card

      NAME = JCARD(1)
      ID   = JCARD(2)
      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read property ID and enter into array PBEAM
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPBEAM-1
            IF (PROPERTY_ID == PBEAM(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
             ENDIF
         ENDDO
         PBEAM(NPBEAM,1) = PROPERTY_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )          ! Read material ID and enter into array PBEAM
      IF (IERRFL(3) == 'N') THEN
         IF (MATERIAL_ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
         ELSE
            PBEAM(NPBEAM,2) = MATERIAL_ID
         ENDIF
      ENDIF

      DO J = 1,6                                           ! Read real property values in fields 4-8
         CALL R8FLD ( JCARD(J+3), JF(J+3), RPBEAM(NPBEAM,J) )
      ENDDO
      IF (IERRFL(4)  == 'N') THEN
         AREA_A = RPBEAM(NPBEAM,1)
      ENDIF
      IF (IERRFL(5)  == 'N') THEN
         I1_A   = RPBEAM(NPBEAM,2)
      ENDIF
      IF (IERRFL(6)  == 'N') THEN
         I2_A   = RPBEAM(NPBEAM,3)
      ENDIF
      IF (IERRFL(7)  == 'N') THEN
         I12_A  = RPBEAM(NPBEAM,4)
      ENDIF
      IF (IERRFL(8)  == 'N') THEN
         JTOR_A = RPBEAM(NPBEAM,5)
      ENDIF
      IF (IERRFL(9)  == 'N') THEN
         NSM_A  = RPBEAM(NPBEAM,6)
      ENDIF

! Call subr to check sensibility of I1, I2, I12 combinations

      CALL CHECK_BAR_MOIs ( 'PBEAM', ID, I1, I2, I12, IERR )
      RPBEAM(NPBEAM,2) = I1
      RPBEAM(NPBEAM,3) = I2
      RPBEAM(NPBEAM,4) = I12
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! Read and check data on mandatory 2nd card (needed even if all fields are blank since 3rd cont is mandatory):

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         DO J = 7,14                                       ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-5), JF(J-5), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID
      ENDIF

! Read and check data on mandatory 3rd card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         CALL CHAR_FLD ( JCARD(2), JF(2), CHRINP )         ! S0 for 1st set of data after data for end A
         IF      (CHRINP(1:4) == 'YES ') THEN
            PBEAM(NPBEAM,3) = 1
         ELSE IF (CHRINP(1:4) == 'YESA') THEN
            PBEAM(NPBEAM,3) = 1
         ELSE IF (CHRINP(1:4) == 'NO  ') THEN
            PBEAM(NPBEAM,3) = 1
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1167) JF(2), NAME, ID, CHRINP
            WRITE(F06,1167) JF(2), NAME, ID, CHRINP
         ENDIF

         CALL R8FLD ( JCARD(3), JF(2), RPBEAM(NPBEAM,15) ) ! X/XB at next section along BEAM

         IF (JCARD(4)(1:) /= ' ') THEN                  ! AREA at next section along BEAM
            CALL R8FLD ( JCARD(4), JF(4), RPBEAM(NPBEAM,16) )
            IF (IERRFL(4) == 'N') THEN
               AREA = RPBEAM(NPBEAM,16)
            ELSE
               AREA = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,16) = AREA_A
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN                     ! I1   at next section along BEAM
            CALL R8FLD ( JCARD(5), JF(5), RPBEAM(NPBEAM,17) )
            IF (IERRFL(5) == 'N') THEN
               I1 = RPBEAM(NPBEAM,17)
            ELSE
               I1 = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,17) = I1
         ENDIF

         IF (JCARD(6)(1:) /= ' ') THEN                     ! I2   at next section along BEAM
            CALL R8FLD ( JCARD(6), JF(6), RPBEAM(NPBEAM,18) )
            IF (IERRFL(6) == 'N') THEN
               I2 = RPBEAM(NPBEAM,18)
            ELSE
               I2 = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,18) = I2
         ENDIF

         IF (JCARD(7)(1:) /= ' ') THEN                     ! I12  at next section along BEAM
            CALL R8FLD ( JCARD(7), JF(7), RPBEAM(NPBEAM,19) )
            IF (IERRFL(7) == 'N') THEN
               I12 = RPBEAM(NPBEAM,19)
            ELSE
               I12 = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,19) = I12
         ENDIF

         IF (JCARD(8)(1:) /= ' ') THEN                     ! JTOR at next section along BEAM
            CALL R8FLD ( JCARD(8), JF(8), RPBEAM(NPBEAM,20) )
            IF (IERRFL(8) == 'N') THEN
               JTOR = RPBEAM(NPBEAM,20)
            ELSE
               JTOR = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,20) = JTOR
         ENDIF

         IF (JCARD(9)(1:) /= ' ') THEN                     ! NSM  at next section along BEAM
            CALL R8FLD ( JCARD(9), JF(9), RPBEAM(NPBEAM,21) )
            IF (IERRFL(9) == 'N') THEN
               NSM = RPBEAM(NPBEAM,21)
            ELSE
               NSM = ZERO
            ENDIF
         ELSE
            RPBEAM(NPBEAM,21) = NSM
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF

! Read and check data on optional 4th card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         DO J = 22,29                                      ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-20), JF(J-20), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
      ENDIF

! Read and check data on optional 5th card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         DO J = 30,37                                      ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-28), JF(J-28), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
      ENDIF

! Read and check data on optional 6th card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         DO J = 38,45                                      ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-36), JF(J-36), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1167 FORMAT(' *ERROR  1167: INCORRECT VALUE IN FIELD ',I2,' OF ',A,A,' = ',A,' MUST BE "YES", "YESA" OR "NO"')

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PBEAM


      SUBROUTINE BD_PLOTEL ( CARD )

! Processes PLOTEL Bulk Data Cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF, MEDAT_PLOTEL, NELE, NPLOTEL
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PLOTEL'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! An integer read
      INTEGER(LONG)                   :: IERR              ! Error count




! **********************************************************************************************************************************
! PLOTEL element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =R1 for CROD
!    2      Element ID     EDAT(nedat+1)
!    3      Grid A         EDAT(nedat+3)
!    4      Grid B         EDAT(nedat+4)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! First, check that fields 2-4 have the proper data type (we are going to have to rearrange the fields prior to calling ELEPRO).
! If any erors, return

      IERR = 0

      DO I=2,4
         CALL I4FLD ( JCARD(I), JF(I), I4INP )
         IF (IERRFL(I) == 'Y') IERR = IERR + 1
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-4
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )   ! Issue warning if fields 5-9 not blank
      CALL CRDERR ( CARD )
      IF (IERR > 0) THEN
         RETURN
      ENDIF

! Set JCARD_EDAT to have fields 2-5 like CROD (but with prop ID = elem ID, where prop ID is not used anyway)

      DO I=1,10
         JCARD_EDAT(I)(1:) = ' '
      ENDDO
      JCARD_EDAT(1) = JCARD(1)                             ! Elem type
      JCARD_EDAT(2) = JCARD(2)                             ! Elem ID
      JCARD_EDAT(3) = JCARD(2)                             ! Prop ID = elem ID (won't be used)
      JCARD_EDAT(4) = JCARD(3)                             ! G1
      JCARD_EDAT(5) = JCARD(4)                             ! G2

! Read and check data

      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_PLOTEL, 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N' )
      NPLOTEL = NPLOTEL + 1
      ETYPE(NELE) = 'PLOTEL  '

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,0,0,0,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-4
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )   ! Issue warning if fields 5-9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields




      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_PLOTEL

   END MODULE ROD_BAR_BEAM_CARDS
