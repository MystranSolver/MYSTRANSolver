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

   MODULE GRID_COORDINATES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_CORD, BD_GRID, BD_GRDSET0, BD_GRDSET, BD_SEQGP, BD_SPOINT0, BD_SPOINT, BD_SNORM

   CONTAINS

      SUBROUTINE BD_CORD ( CARD, LARGE_FLD_INP )

! Processes CORD1C, CORD1R, CORD1S and CORD2C, CORD2R, CORD2S Bulk Data Cards
!  1) Sets coord type  (0 {2R},1 {2C},2 {2S}) and enters it into array CORD
!  2) Reads coord system ID and reference ID  and enters it into array CORD
!  3) Reads coord data into array RCORD

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, IERRFL, JCARD_LEN, JF, LCORD, NCORD, NCORD1, NCORD2, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CORD, RCORD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CORD'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CORD_CID          ! Field 2 of CORD card (coord sys ID)
      CHARACTER(LEN(JCARD))           :: CORD_NAME         ! Name of coors sys

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! CORD1R Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Cord Type       CORD(ncord,1) 11 is CORD1R, 21 is CORD2R, 22 is CORD2C, 23 is CORD2S
!    2      CID             CORD(ncord,2)
!    3      GA              Temporarily put into CORD(ncord,3)
!    3      GB              Temporarily put into CORD(ncord,4)
!    3      GC              Temporarily put into CORD(ncord,5)
!    3      RID             CORD(ncord,3) ref sys for grid A (will be entered later when GRID array is sorted and we can find GA
!    4      RID             CORD(ncord,4) ref sys for grid B (will be entered later when GRID array is sorted and we can find GB
!    5      RID             CORD(ncord,5) ref sys for grid C (will be entered later when GRID array is sorted and we can find GC


! CORD2C, CORD2R, CORD2S Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
! on first card:
!    1      Cord Type       CORD(ncord,1) =02 {2R},12 {2C},22 {2S}
!    2      CID             CORD(ncord,2)
!    3      RID             CORD(ncord,3)
!    4      A1             RCORD(ncord,1)
!    5      A2             RCORD(ncord,2)
!    6      A3             RCORD(ncord,3)
!    7      B1             RCORD(ncord,4)
!    8      B2             RCORD(ncord,5)
!    9      B3             RCORD(ncord,6)
! on required second card:
!    2      C1             RCORD(ncord,7)
!    3      C2             RCORD(ncord,8)
!    4      C3             RCORD(ncord,9)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      CORD_NAME = JCARD(1)

! ---------------------------------------------------------------------------------------------------------------------------------
      IF (CORD_NAME(1:5) == 'CORD1') THEN
                                                           ! Read data for 1st coord system defined on this entry
         NCORD1 = NCORD1 + 1
         NCORD  = NCORD  + 1

         CORD_CID = JCARD(2)
         IF      (JCARD(1)(1:6) == 'CORD1R') THEN
            CORD(NCORD,1) = 11
         ELSE IF (JCARD(1)(1:6) == 'CORD1C') THEN
            CORD(NCORD,1) = 12
         ELSE IF (JCARD(1)(1:6) == 'CORD1S') THEN
            CORD(NCORD,1) = 13
         ENDIF

         CALL I4FLD ( JCARD(2), JF(2), I4INP )             ! Read CID and make sure it is > 0 (cannot define 0, or basic, system)
         IF (IERRFL(2) == 'N') THEN
            IF (I4INP < 0) THEN                            ! --- CID cannot be negative
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(2), CORD_NAME, JCARD(2), JCARD(2)
               WRITE(F06,1169) JF(2), CORD_NAME, JCARD(2), JCARD(2)
            ELSE IF (I4INP == 0) THEN                      ! --- CID cannot be 0 (can't define basic)
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1170) JF(2), CORD_NAME, JCARD(2), JCARD(2)
               WRITE(F06,1170) JF(2), CORD_NAME, JCARD(2), JCARD(2)
            ELSE                                           ! --- CID is OK
               CORD(NCORD,2) = I4INP
            ENDIF
         ENDIF

         CALL I4FLD ( JCARD(3), JF(3), I4INP )             ! --- Read GA. If > 0 put it, temporarily, into CORD(ncord,3)
         IF (IERRFL(3) == 'N') THEN
            IF (I4INP >= 0) THEN
               CORD(NCORD,3) = I4INP
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(3), CORD_NAME, CORD_CID, JCARD(3)
               WRITE(F06,1169) JF(3), CORD_NAME, CORD_CID, JCARD(3)
            ENDIF
         ENDIF

         CALL I4FLD ( JCARD(4), JF(4), I4INP )             ! --- Read GB. If > 0 put it, temporarily, into CORD(ncord,4)
         IF (IERRFL(4) == 'N') THEN
            IF (I4INP >= 0) THEN
               CORD(NCORD,4) = I4INP
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(4), CORD_NAME, CORD_CID, JCARD(4)
               WRITE(F06,1169) JF(4), CORD_NAME, CORD_CID, JCARD(4)
            ENDIF
         ENDIF

         CALL I4FLD ( JCARD(5), JF(5), I4INP )             ! --- Read GC. If > 0 put it, temporarily, into CORD(ncord,5)
         IF (IERRFL(5) == 'N') THEN
            IF (I4INP >= 0) THEN
               CORD(NCORD,5) = I4INP
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(5), CORD_NAME, CORD_CID, JCARD(5)
               WRITE(F06,1169) JF(5), CORD_NAME, CORD_CID, JCARD(5)
            ENDIF
         ENDIF


         IF (JCARD(6)(1:) == ' ') THEN                     ! There is only 1 CORD1 defined on this entry

            CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,0,0,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )
            CALL CRDERR ( CARD )

         ELSE                                              ! A 2nd CORD1 is on this entry since field 6 is non-lank

            NCORD1 = NCORD1 + 1
            NCORD  = NCORD  + 1

            CORD_CID = JCARD(6)
            CORD(NCORD,1) = 11

            CALL I4FLD ( JCARD(6), JF(6), I4INP )          ! Read CID and make sure it is > 0 (cannot define 0, or basic, system)
            IF (IERRFL(6) == 'N') THEN
               IF (I4INP < 0) THEN                         ! --- CID cannot be negative
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1169) JF(6), CORD_NAME, JCARD(6), JCARD(6)
                  WRITE(F06,1169) JF(6), CORD_NAME, JCARD(6), JCARD(6)
               ELSE IF (I4INP == 0) THEN                   ! --- CID cannot be 0 (can't define basic)
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1170) JF(6), CORD_NAME, JCARD(6), JCARD(6)
                  WRITE(F06,1170) JF(6), CORD_NAME, JCARD(6), JCARD(6)
               ELSE                                        ! --- CID is OK
                  CORD(NCORD,2) = I4INP
               ENDIF
            ENDIF

            CALL I4FLD ( JCARD(7), JF(7), I4INP )          ! --- Read GA. If > 0 put it, temporarily, into CORD(ncord,7)
            IF (IERRFL(7) == 'N') THEN
               IF (I4INP >= 0) THEN
                  CORD(NCORD,3) = I4INP
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1169) JF(7), CORD_NAME, CORD_CID, JCARD(7)
                  WRITE(F06,1169) JF(7), CORD_NAME, CORD_CID, JCARD(7)
               ENDIF
            ENDIF

            CALL I4FLD ( JCARD(8), JF(8), I4INP )          ! --- Read GB. If > 0 put it, temporarily, into CORD(ncord,8)
            IF (IERRFL(8) == 'N') THEN
               IF (I4INP >= 0) THEN
                  CORD(NCORD,4) = I4INP
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1169) JF(8), CORD_NAME, CORD_CID, JCARD(8)
                  WRITE(F06,1169) JF(8), CORD_NAME, CORD_CID, JCARD(8)
               ENDIF
            ENDIF

            CALL I4FLD ( JCARD(9), JF(9), I4INP )          ! --- Read GC. If > 0 put it, temporarily, into CORD(ncord,9)
            IF (IERRFL(9) == 'N') THEN
               IF (I4INP >= 0) THEN
                  CORD(NCORD,5) = I4INP
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1169) JF(9), CORD_NAME, CORD_CID, JCARD(9)
                  WRITE(F06,1169) JF(9), CORD_NAME, CORD_CID, JCARD(9)
               ENDIF
            ENDIF

            CALL BD_IMBEDDED_BLANK   ( JCARD,0,0,0,0,6,7,8,9 )
            CALL CRDERR ( CARD )

         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (CORD_NAME(1:5) == 'CORD2') THEN

         NCORD2 = NCORD2 + 1
         NCORD  = NCORD  + 1

         CORD_CID = JCARD(2)

         IF      (JCARD(1)(1:6) == 'CORD2R') THEN
            CORD(NCORD,1) = 21
         ELSE IF (JCARD(1)(1:6) == 'CORD2C') THEN
            CORD(NCORD,1) = 22
         ELSE IF (JCARD(1)(1:6) == 'CORD2S') THEN
            CORD(NCORD,1) = 23
         ENDIF

         CALL I4FLD ( JCARD(2), JF(2), I4INP )             ! Read CID and make sure it is > 0 (cannot define 0, or basic, system)
         IF (IERRFL(2) == 'N') THEN
            IF (I4INP < 0) THEN                            ! --- CID cannot be negative
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(2), CORD_NAME, JCARD(2), JCARD(2)
               WRITE(F06,1169) JF(2), CORD_NAME, JCARD(2), JCARD(2)
            ELSE IF (I4INP == 0) THEN                      ! --- CID cannot be 0 (can't define basic)
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1170) JF(2), CORD_NAME, JCARD(2), JCARD(2)
               WRITE(F06,1170) JF(2), CORD_NAME, JCARD(2), JCARD(2)
            ELSE                                           ! --- CID is OK
               CORD(NCORD,2) = I4INP
            ENDIF
         ENDIF

         CALL I4FLD ( JCARD(3), JF(3), I4INP )             ! Read RID and make sure it is >= 0
         IF (IERRFL(3) == 'N') THEN
            IF (I4INP >= 0) THEN
               CORD(NCORD,3) = I4INP
            ELSE                                           ! --- RID cannot be negative
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1169) JF(3), CORD_NAME, JCARD(3), JCARD(3)
               WRITE(F06,1169) JF(3), CORD_NAME, JCARD(3), JCARD(3)
            ENDIF
         ENDIF

         DO J = 1,6                                        ! Read real data on parent card
            CALL R8FLD ( JCARD(J+3), JF(J+3), RCORD(NCORD,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
         CALL CRDERR ( CARD )

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )              ! Read 2nd card
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            CALL R8FLD ( JCARD(2), JF(2), RCORD(NCORD,7) )
            CALL R8FLD ( JCARD(3), JF(3), RCORD(NCORD,8) )
            CALL R8FLD ( JCARD(4), JF(4), RCORD(NCORD,9) )

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )
            CALL CRDERR ( CARD )

         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1136) CORD_NAME, CORD_CID
            WRITE(F06,1136) CORD_NAME, CORD_CID
         ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1169 FORMAT(' *ERROR  1169: FIELD ',I3,' ON ',A,' ID ',A,' CANNOT BE NEGATIVE. VALUE IS = ',A)

 1170 FORMAT(' *ERROR  1170: FIELD ',I3,' ON ',A,' ID ',A,' CANNOT BE 0 (CANNOT DEFINE BASIC SYSTEM). VALUE IS = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CORD


      SUBROUTINE BD_GRID ( CARD )

! Processes GRID Bulk Data Cards. Reads and checks:

!  1) Grid ID (field 2) and enters it into array GRID
!  2) Input  coord sys (field 3) and enters it into array GRID (or uses GRDSET3, if input on a GRDSET card)
!  3) Global coord sys (field 7) and enters it into array GRID (or uses GRDSET7, if input on a GRDSET card)
!  4) Permanent SPC's  (field 8) and enters it into array GRID (or uses GRDSET8, if input on a GRDSET card)
!  5) Grid coordinates (fields 4, 5 and 6) and enters tham into array RGRID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LGRID, NGRID, NGRDSET
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID, RGRID, GRDSET3, GRDSET7, GRDSET8

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK, R8FLD
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_GRID'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! The type of TOKEN (looking for 'INTEGER ')

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein




! **********************************************************************************************************************************
! GRID Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Grid number    GRID (ngrid,1)
!    3      Input co-ord   GRID (ngrid,2)
!    4      X1             RGRID(ngrid,1)
!    5      X2             RGRID(ngrid,2)
!    6      X3             RGRID(ngrid,3)
!    7      Disp. co-ord   GRID (ngrid,3)
!    8      PSPC           GRID (ngrid,4)
!   10      LB = 1         GRID (ngrid,5) key to put a line break into F06 file after outputs for this grid
! none      GRID/SPOINT    GRID (ngrid,6) = 6 for a physical grid

! Make JCARD from  ARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NGRID = NGRID+1
!xx   IF (NGRID > LGRID) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LGRID
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LGRID
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), I4INP )                ! Read grid ID
      IF (IERRFL(2) == 'N') THEN
         IF (I4INP > 0) THEN
            GRID(NGRID,1) = I4INP
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1168) JCARD(2)
            WRITE(F06,1168) JCARD(2)
         ENDIF
      ENDIF

      IF (JCARD(3)(1:) == ' ') THEN                        ! Read input coord sys or set to GRDSET value
         IF (NGRDSET > 0) THEN
            GRID(NGRID,2) = GRDSET3
         ELSE
            GRID(NGRID,2) = 0
         ENDIF
      ELSE
         CALL I4FLD ( JCARD(3), JF(3), GRID(NGRID,2) )
      ENDIF

      CALL R8FLD ( JCARD(4), JF(4), RGRID(NGRID,1) )       ! Read 1st component of grid coordinates
      CALL R8FLD ( JCARD(5), JF(5), RGRID(NGRID,2) )       ! Read 2nd component of grid coordinates
      CALL R8FLD ( JCARD(6), JF(6), RGRID(NGRID,3) )       ! Read 3rd component of grid coordinates

      IF (JCARD(7)(1:) == ' ') THEN                        ! Read global coord sys or set to GRDSET value
         IF (NGRDSET > 0) THEN
            GRID(NGRID,3) = GRDSET7
         ELSE
            GRID(NGRID,3) = 0
         ENDIF
      ELSE
         CALL I4FLD ( JCARD(7), JF(7), GRID(NGRID,3) )
      ENDIF

      IF (JCARD(8)(1:) == ' ') THEN                        ! Read perm SPC's or set to GRDSET value
         IF (NGRDSET > 0) THEN
            GRID(NGRID,4) = GRDSET8
         ELSE
            GRID(NGRID,4) = 0
         ENDIF
      ELSE
         CALL IP6CHK ( JCARD(8), JCARDO, IP6TYP, IDUM )
         IF (IP6TYP == 'COMP NOS') THEN
            CALL I4FLD ( JCARDO, JF(8), GRID(NGRID,4) )
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1124) JF(8), JCARD(1), JCARD(2), JF(8), JCARD(8)
            WRITE(F06,1124) JF(8), JCARD(1), JCARD(2), JF(8), JCARD(8)
         ENDIF
      ENDIF

      TOKEN = JCARD(10)(1:8)                               ! Only send the 1st 8 chars of this JCARD. It has been left justified
      CALL TOKCHK ( TOKEN, TOKTYP )                        ! See if field 10 has an integer. If so call it LB
      IF (TOKTYP == 'INTEGER ') THEN
         CALL I4FLD ( JCARD(10), JF(10), I4INP )           ! Read LB (put this many line breaks in F06 after this grid)
         IF (IERRFL(10) == 'N') THEN
            GRID(NGRID,5) = I4INP
         ENDIF
      ENDIF

      GRID(NGRID,6) = 6                                    ! Set GRID(NGRID,6) = 6 to indicate 6 comps for a physical grid

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )     ! Make sure that there are no imb blanks in fields 2-7. Field 8 is PSPC
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1168 FORMAT(' *ERROR  1168: ZERO OR NEGATIVE GRID ID NOT ALLOWED ON GRID CARD. VALUE IS = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_GRID


      SUBROUTINE BD_GRDSET0 ( CARD )

! Processes GRDSET Bulk Data Cards from LOADB0 to set GRDSET3, 7, 8 prior
! to reading bulk data in LOADB. If there are errors reading the GRDSET
! card entries, messages are not printed out until LOADB calls BD_GRDSET

!  1) GRDSET3 is field 3 for a GRID card (the input  coord system)
!  2) GRDSET7 is field 7 for a GRID card (the global coord system)
!  2) GRDSET7 is field 8 for a GRID card (the perm SPC's)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, JF, IERRFL, NGRDSET
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRDSET3, GRDSET7, GRDSET8

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD, IP6CHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_GRDSET0'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein




! **********************************************************************************************************************************
! GRDSET Bulk Data Card routine

!   FIELD         ITEM       ARRAY ELEMENT
!   -----   ---------------- ------------
!    3      Input Coord sys     GRDSET3
!    7      Displ Coord sys     GRDSET7
!    8      PSPC                GRDSET8


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Count the number of GRDSET cards. The count will be checked to make sure that there is no more than 1 in subr LOADB

      NGRDSET = NGRDSET + 1

! Only set GRDSET3, 7, 8 if this is the 1st GRDSET card.
! Reset IERRFL's to 'N', after testing, since CDRERR is not being called until BD_GRDSET is called from LOADB

      IF (NGRDSET == 1) THEN

         IF (JCARD(3)(1:) /= ' ') THEN
            CALL I4FLD ( JCARD(3), JF(3), I4INP )
            IF (IERRFL(3) == 'N') THEN
               GRDSET3 = I4INP
            ELSE
               IERRFL(3) = 'N'
            ENDIF
         ENDIF

         IF (JCARD(7)(1:) /= ' ') THEN
            CALL I4FLD ( JCARD(7), JF(7), I4INP )
            IF (IERRFL(7) == 'N') THEN
               GRDSET7 = I4INP
            ELSE
               IERRFL(7) = 'N'
            ENDIF
         ENDIF

         IF (JCARD(8)(1:) /= ' ') THEN
            CALL IP6CHK ( JCARD(8), JCARDO, IP6TYP, IDUM )
            IF (IP6TYP == 'COMP NOS') THEN
               CALL I4FLD ( JCARDO, JF(8), I4INP )
               IF (IERRFL(8) == 'N') THEN
                  GRDSET8 = I4INP
               ELSE
                  IERRFL(8) = 'N'
               ENDIF
            ENDIF
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_GRDSET0


      SUBROUTINE BD_GRDSET ( CARD )

! Processes GRDSET Bulk Data Cards. Subr GRDSET0 read fields 3, 7 and 8 of the GRDSET card so that those values
! could be used on GRID cards that may have been in the deck before the GRDSET card. This subr reads tham again
! and checks them for error:

!  1) GRDSET3 is field 3 for a GRID card (the input  coord system)
!  2) GRDSET7 is field 7 for a GRID card (the global coord system)
!  2) GRDSET7 is field 8 for a GRID card (the perm SPC's)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NGRDSET
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRDSET3, GRDSET7, GRDSET8

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_GRDSET'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: PGM_ERR   = 0     ! A  count of the number of coding errors




! **********************************************************************************************************************************
! GRDSET Bulk Data Card routine. Values for GRDSET3, 7, 8 have already been
! read when B.D. deck was scanned originally. Here, we read card to detect
! and report read errors and to check for coding error (if values don't
! agree with those read in BD_GRDSET0)

!   FIELD         ITEM       ARRAY ELEMENT
!   -----   ---------------- ------------
!    3      Input Coord sys     GRDSET3
!    7      Displ Coord sys     GRDSET7
!    8      PSPC                GRDSET8


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read and check data

      IF (JCARD(3)(1:) /= ' ') THEN                        ! Read input coord sys ID in field 3 and store in GRDSET3
         CALL I4FLD ( JCARD(3), JF(3), I4INP )
         IF (IERRFL(3) == 'N') THEN
            IF (I4INP /= GRDSET3) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,1185) SUBR_NAME,JF(3),JCARD(1),GRDSET3,I4INP
               WRITE(F06,1185) SUBR_NAME,JF(3),JCARD(1),GRDSET3,I4INP
            ENDIF
         ENDIF
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Read global coord sys ID in field 7 and store in GRDSET7
         CALL I4FLD ( JCARD(7), JF(7), I4INP )
         IF (IERRFL(7) == 'N') THEN
            IF (I4INP /= GRDSET7) THEN
               PGM_ERR = PGM_ERR + 1                       ! Coding error: This value doesn't agree with that read in LOADB0
               WRITE(ERR,1185) SUBR_NAME,JF(7),JCARD(1),GRDSET7,I4INP
               WRITE(F06,1185) SUBR_NAME,JF(7),JCARD(1),GRDSET7,I4INP
            ENDIF
         ENDIF
      ENDIF

      IF (JCARD(8)(1:) /= ' ') THEN                        ! Read perm SPC's in field 8 and store in GRDSET8
         CALL IP6CHK ( JCARD(8), JCARDO, IP6TYP, IDUM )
         IF (IP6TYP == 'COMP NOS') THEN
            CALL I4FLD ( JCARDO, JF(8), I4INP )
            IF (IERRFL(8) == 'N') THEN
               IF (I4INP /= GRDSET8) THEN
                  PGM_ERR = PGM_ERR + 1                    ! Coding error: This value doesn't agree with that read in LOADB0
                  WRITE(ERR,1185) SUBR_NAME,JF(8),JCARD(1),GRDSET8,I4INP
                  WRITE(F06,1185) SUBR_NAME,JF(8),JCARD(1),GRDSET8,I4INP
               ENDIF
            ENDIF
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1123) JF(8),JCARD(1),JF(8),JCARD(8)
            WRITE(F06,1123) JF(8),JCARD(1),JF(8),JCARD(8)
         ENDIF
      ENDIF

      CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,7,0,0 )   ! Make sure that there are no imbed blanks in fields 3,7. Field 8 is PSPC
      CALL CARD_FLDS_NOT_BLANK ( JCARD,2,0,4,5,6,0,0,9 )   ! Issue warning if fields 2, 4, 5, 6, 9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF (PGM_ERR > 0) THEN                                ! PGM_ERR /= 0 is a coding error, so quit
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1123 FORMAT(' *ERROR  1123: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' CARD. MUST BE A COMBINATION OF DIGITS 1-6'                &
                    ,/,14X,'HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1185 FORMAT(' *ERROR  1185: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE READ FROM FIELD ',I3,' ON ',A,' CARD IN THIS SUBROUTINE'                                        &
                    ,/,14X,' DOES NOT AGREE WITH THAT READ FROM ORIGINAL BULK DATA DECK SCAN.'                                     &
                    ,/,14X,' VALUES ARE: ',I8,' (ORIGINAL BULK DATA SCAN),'                                                        &
                    ,/,14X,'        AND: ',I8,' (HERE)')

! **********************************************************************************************************************************

      END SUBROUTINE BD_GRDSET


      SUBROUTINE BD_SEQGP ( CARD )

!  Processes SEQGP Bulk Data Cards. Reads and checks data:

!  1) Grid ID's entered into integer array SEQ1
!  2) Sequence ID's entered into real array SEQ2. If input sequence numbers are integer, they are converted
!     to real before entering them into array SEQ2.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, JF, LSEQ, NSEQ
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  SEQ1, SEQ2

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SEQGP'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: DEC_COL           ! Number that indicates whether entry in a seq field is real or integer
      INTEGER(LONG)                   :: ISEQ              ! An integer sequence number
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JFLD1             ! A field number on the SEQGP card where grid ID's are located
      INTEGER(LONG)                   :: JFLD2             ! A field number on the SEQGP card where sequence numbers are located


      REAL(DOUBLE)                    :: RSEQ              ! A real sequence number

      INTRINSIC INDEX,DBLE



! **********************************************************************************************************************************
!  SEQGP Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!    2,4,6,8 Grid ID's      SEQ1(nseq - nseq+4)
!    3,5,7,9 Seq ID's       SEQ2(nseq - nseq+4) Can be real or integer

!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read and check data

      DO J = 1,4
         JFLD1 = 2*J
         JFLD2 = 2*J + 1
         IF (JCARD(JFLD1)(1:) /= ' ') THEN                 ! Check grid field. If blank CYCLE
            NSEQ = NSEQ+1                                  ! Increment NSEQ if field 2,4,6 or 8 is not blank
            CALL I4FLD(JCARD(JFLD1), JF(JFLD1), SEQ1(NSEQ))! Read grid ID in field 2,4,6 or 8
            IF (JCARD(JFLD2)(1:) /= ' ') THEN              ! If sequence field is not blank, OK. Otherwise, error
               DEC_COL = INDEX(JCARD(JFLD2),'.')
               IF (DEC_COL == 0) THEN                      ! Integer entry for sequence number
                  CALL I4FLD ( JCARD(JFLD2), JF(JFLD2), ISEQ )
                  SEQ2(NSEQ) = DBLE(ISEQ)
               ELSE                                        ! Real entry for sequence number
                  CALL R8FLD ( JCARD(JFLD2), JF(JFLD2), RSEQ )
                  SEQ2(NSEQ) = RSEQ
               ENDIF
               IF (SEQ2(NSEQ) <= ZERO) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1157) JFLD2,JCARD(JFLD2)
                  WRITE(F06,1157) JFLD2,JCARD(JFLD2)
               ENDIF
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1159) JFLD2
               WRITE(F06,1159) JFLD2
            ENDIF
         ELSE
            CYCLE
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1157 FORMAT(' *ERROR  1157: VALUE ON SEQGP ENTRY, FIELD ',I3,', MUST BE >= 0 BUT IS = ',A)

 1159 FORMAT(' *ERROR  1159: FIELD ',I3,' OF SEQGP ENTRY CANNOT BE BLANK. MUST CONTAIN AN INTEGER OR REAL SEQUENCE NUMBER')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_SEQGP


      SUBROUTINE BD_SPOINT0 ( CARD, LARGE_FLD_INP, DELTA_SPOINT )

! Counts the SPOINT's of one SPOINT Bulk Data entry in the first pass over the Bulk Data. The list (fields 2-9 of the parent line and
! fields 2-9 of the continuations, single IDs and ranges "ID1 THRU ID2" in any mix) is read by subr READ_ID_LIST as in BD_SPOINT,
! so the count is the number of SPOINT's that BD_SPOINT enters into array GRID. An entry with an error is not counted.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE BDF_ID_LISTS, ONLY          :  READ_ID_LIST

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card (the parent line; continuations are read)
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      INTEGER(LONG), INTENT(OUT)      :: DELTA_SPOINT      ! Number of SPOINT's defined on this B.D. SPOINT entry
      INTEGER(LONG), ALLOCATABLE      :: RANGES(:,:)       ! The SPOINT IDs: single (S, S) or ranges (S1, S2)
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: NBAD              ! Number of list fields in error
      INTEGER(LONG)                   :: NRANGE            ! Number of items in RANGES

! **********************************************************************************************************************************
      CALL READ_ID_LIST ( CARD, LARGE_FLD_INP, 2_LONG, .TRUE., NRANGE, RANGES, NBAD )

      DELTA_SPOINT = 0
      IF (NBAD == 0) THEN
         DO J=1,NRANGE
            DELTA_SPOINT = DELTA_SPOINT + RANGES(2,J) - RANGES(1,J) + 1
         ENDDO
      ENDIF

      RETURN

      END SUBROUTINE BD_SPOINT0


      SUBROUTINE BD_SPOINT ( CARD, LARGE_FLD_INP )

! Read Bulk Data SPOINT entries. Enter the SPOINT number into array GRID (in col 1) and set GRID(ngrid,6) to 1 to indicate SPOINT.
! The list is fields 2-9 of the parent line and fields 2-9 of the continuations: single IDs and ranges "ID1 THRU ID2" in any mix
! (subr READ_ID_LIST; before, a range was read only as the whole entry "ID1 THRU ID2" and the continuations were not read).

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  NGRID
      USE MODEL_STUF, ONLY            :  GRID
      USE BDF_ID_LISTS, ONLY          :  READ_ID_LIST

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card (the parent line; continuations are read)
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      INTEGER(LONG), ALLOCATABLE      :: RANGES(:,:)       ! The SPOINT IDs: single (S, S) or ranges (S1, S2)
      INTEGER(LONG)                   :: J, K              ! DO loop indices
      INTEGER(LONG)                   :: NBAD              ! Number of list fields in error
      INTEGER(LONG)                   :: NRANGE            ! Number of items in RANGES

! **********************************************************************************************************************************
      CALL READ_ID_LIST ( CARD, LARGE_FLD_INP, 2_LONG, .FALSE., NRANGE, RANGES, NBAD )

      IF (NBAD == 0) THEN                                  ! As counted in BD_SPOINT0 (an entry with an error is not counted)
         DO J=1,NRANGE
            DO K=RANGES(1,J),RANGES(2,J)
               NGRID = NGRID + 1
               GRID(NGRID,1) = K
               GRID(NGRID,6) = 1
            ENDDO
         ENDDO
      ENDIF

      RETURN

      END SUBROUTINE BD_SPOINT


      SUBROUTINE BD_SNORM ( CARD )

! Processes SNORM Bulk Data Cards
!  1) Reads grid point ID and coord system ID and enters it into array SNORM
!  2) Reads normal vector into array RSNORM

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, JF, NSNORM
      USE MODEL_STUF, ONLY            :  SNORM, RSNORM
      USE CONSTANTS_1, ONLY           :  ZERO

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SNORM'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value

      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value


! **********************************************************************************************************************************
! SNORM Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      GID            SNORM(nsnorm,1)
!    3      CID            SNORM(nsnorm,2)
!    4      N1             RSNORM(nsnorm,1)
!    5      N2             RSNORM(nsnorm,2)
!    6      N3             RSNORM(nsnorm,3)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NSNORM  = NSNORM  + 1

      CALL I4FLD ( JCARD(2), JF(2), I4INP )
      SNORM(NSNORM,1) = I4INP

      CALL I4FLD ( JCARD(3), JF(3), I4INP )
      SNORM(NSNORM,2) = I4INP

      CALL R8FLD ( JCARD(4), JF(4), R8INP )
      RSNORM(NSNORM,1) = R8INP

      CALL R8FLD ( JCARD(5), JF(5), R8INP )
      RSNORM(NSNORM,2) = R8INP

      CALL R8FLD ( JCARD(6), JF(6), R8INP )
      RSNORM(NSNORM,3) = R8INP

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,0,0,0 )
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,7,8,9 )
      CALL CRDERR ( CARD )


! **********************************************************************************************************************************

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_SNORM

   END MODULE GRID_COORDINATES
