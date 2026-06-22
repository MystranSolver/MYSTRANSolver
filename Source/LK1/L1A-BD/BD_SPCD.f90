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
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE BD_SPCD ( CARD, CC_LOAD_FND )

! Processes SPCD Bulk Data Cards. Native NX decks use SPCD as a load-set-driven
! enforced displacement, so the set ID must satisfy Case Control LOAD = SID.
! The resulting constraint records are written to LINK1O in the same packed form
! used by SPC/SPC1 so later DOF processing can place them in the SE-set.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1O
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSPC, LSUB, NSPC, NUM_SPC_RECORDS,&
                                         NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SUPWARN
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN
      USE MODEL_STUF, ONLY            :  SPC_SIDS, SUBLOD

      USE BD_SPC_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SPCD'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET
      CHARACTER( 8*BYTE)              :: IP6TYP
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN(JCARD))           :: JCARDO

      INTEGER(LONG)                   :: COMPJ     = 0
      INTEGER(LONG)                   :: GRIDJ     = 0
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IDUM
      INTEGER(LONG)                   :: JERR      = 0
      INTEGER(LONG)                   :: SETID     = 0

      REAL(DOUBLE)                    :: RSPCJ     = ZERO

! Make JCARD from CARD
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow
      NSPC = NSPC + 1

! Read set ID from LOAD namespace
      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
         SPC_SIDS(NSPC) = SETID
      ELSE
         JERR = JERR + 1
      ENDIF

      DO I=1,2
         IF (JCARD(3*I)(1:) /= ' ') THEN

            JERR = 0
            CALL I4FLD ( JCARD(3*I), JF(3*I), GRIDJ )
            CALL IP6CHK ( JCARD(3*I+1), JCARDO, IP6TYP, IDUM )
            IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
               CALL I4FLD ( JCARDO, JF(3*I+1), COMPJ )
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1124) JF(3*I+1),JCARD(1),JCARD(2),JF(3*I+1),JCARD(3*I+1)
               WRITE(F06,1124) JF(3*I+1),JCARD(1),JCARD(2),JF(3*I+1),JCARD(3*I+1)
            ENDIF

            CALL R8FLD ( JCARD(3*I+2), JF(3*I+2), RSPCJ )
            DOFSET = 'SE'

            IF ((JERR == 0 ) .AND. (IERRFL(3*I) == 'N') .AND. (IERRFL(3*I+1) == 'N') .AND. (IERRFL(3*I+2) == 'N')) THEN
               NUM_SPC_RECORDS = NUM_SPC_RECORDS + 1
               WRITE(L1O) SETID,COMPJ,GRIDJ,GRIDJ,RSPCJ,DOFSET
            ENDIF

         ELSE
            IF ((JCARD(3*I+1)(1:) /= ' ') .OR. (JCARD(3*I+2)(1:) /= ' ')) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1144) JCARD(1),JCARD(2),JF(3*I+1),JF(3*I+2),JF(3*I)
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1144) JCARD(1),JCARD(2),JF(3*I+1),JF(3*I+2),JF(3*I)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,0,5,6,0,8,0 )
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )
      CALL CRDERR ( CARD )

      RETURN

  101 FORMAT(A)
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')
 1144 FORMAT(' *WARNING    : ON ',A,' SET ID = ',A,' FIELDS ',I3,' AND ',I3,' ARE IGNORED SINCE GRID FIELD ',I3,' IS BLANK.')

      END SUBROUTINE BD_SPCD
