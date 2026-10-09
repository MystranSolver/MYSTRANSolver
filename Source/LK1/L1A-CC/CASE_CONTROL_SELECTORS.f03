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

   MODULE CASE_CONTROL_SELECTORS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CC_LOAD, CC_METH, CC_MPC, CC_NLPARM, CC_SPC, CC_STATSUB, CC_TEMP

   CONTAINS

      SUBROUTINE CC_LOAD ( CARD )

! Processes Case Control LOAD cards that define load sets to be applied

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SUBLOD

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_LOAD'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
! Process LOAD cards

! Get SETID

      CALL GET_SETID ( CARD, SETID )

! Set CASE CONTROL variable to SETID

      IF (NSUB /= 0) THEN
         SUBLOD(NSUB,1) = SETID
      ELSE
         DO I=1,LSUB
            SUBLOD(I,1) = SETID
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_LOAD


      SUBROUTINE CC_METH ( CARD )

! Processes Case Control eigenvalue METHOD cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, BLNK_SUB_NAM, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, CC_EIGR_SID_SUB, CC_EIGR_SID_DECK, IS_MODES_SUBCASE

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_METH'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
! Process METHOD cards

! Get SETID

      CALL GET_SETID ( CARD, SETID )

! Record per-subcase METHOD assignment so that SOL 103 decks can request a different set of modes for each subcase.
! NSUB is incremented by CC_SUBC at parse time, so:
!   * NSUB == 0  -> this METHOD appears above any SUBCASE card; it is the deck-default that any subcase lacking its own
!                   METHOD inherits during LOADC's post-parse pass.
!   * NSUB >  0  -> this METHOD belongs to the current (most-recently-opened) subcase.

      IF (NSUB == 0) THEN

         IF ((CC_EIGR_SID_DECK /= 0) .AND. (CC_EIGR_SID_DECK /= SETID)) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,8867) CC_EIGR_SID_DECK, SETID
            IF (SUPWARN == 'N') THEN
               WRITE(F06,8867) CC_EIGR_SID_DECK, SETID
            ENDIF
         ENDIF
         CC_EIGR_SID_DECK = SETID

      ELSE

         IF (ALLOCATED(CC_EIGR_SID_SUB)) THEN
            IF ((CC_EIGR_SID_SUB(NSUB) /= 0) .AND. (CC_EIGR_SID_SUB(NSUB) /= SETID)) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,8868) NSUB, CC_EIGR_SID_SUB(NSUB), SETID
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,8868) NSUB, CC_EIGR_SID_SUB(NSUB), SETID
               ENDIF
            ENDIF
            CC_EIGR_SID_SUB(NSUB) = SETID
            IS_MODES_SUBCASE(NSUB) = 'Y'
         ENDIF

      ENDIF

! Maintain the legacy scalar CC_EIGR_SID so existing single-METHOD code paths (BD_EIGR scalar match, WRITE_L1Z, restart
! sanity check, etc.) continue to work unchanged. After LOADC inheritance the scalar reflects the last seen SID.

      CC_EIGR_SID = SETID



      RETURN

! **********************************************************************************************************************************
 8867 FORMAT(' *WARNING    : MORE THAN ONE DECK-LEVEL METHOD ENTRY IN CASE CONTROL. PREVIOUS SET ID = ',I8,', NEW SET ID = ',I8,    &
             '. NEW VALUE WILL BE USED AS THE DECK DEFAULT.')
 8868 FORMAT(' *WARNING    : MORE THAN ONE METHOD ENTRY IN SUBCASE ',I8,'. PREVIOUS SET ID = ',I8,', NEW SET ID = ',I8,             &
             '. NEW VALUE WILL BE USED.')

! **********************************************************************************************************************************

      END SUBROUTINE CC_METH


      SUBROUTINE CC_MPC ( CARD )

! Processes Case Control MPC cards for defining MPC set ID's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  MPCSETS

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_MPC'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! Process MPC cards

! Get SETID

      CALL GET_SETID ( CARD ,SETID )

! Set Case Control variable to SETID

      IF (NSUB /= 0) THEN
         MPCSETS(NSUB) = SETID
      ELSE
         DO I=1,LSUB
            MPCSETS(I) = SETID
         ENDDO
      ENDIF


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_MPC


      SUBROUTINE CC_NLPARM ( CARD )

! Processes Case Control NLPARM cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LSUB, NSUB, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE NONLINEAR_PARAMS, ONLY      :  NL_SID

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_NLPARM'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
! Process NLPARM cards

! Find out if "NONE", "ALL" or SETID

      CALL GET_SETID ( CARD, SETID )

! Set CASE CONTROL variable to SETID

      IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
         IF (NSUB /= 0) THEN
            NL_SID(NSUB) = SETID
         ELSE
            DO I = 1,LSUB
               NL_SID(I) = SETID
           ENDDO
         ENDIF
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1203) 'NLPARM'
         WRITE(F06,1203) 'NLPARM'
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1203 FORMAT(' *ERROR  1203: CASE CONTROL ENTRY ',A,' ONLY VALID IN SOL NLSTATIC (SOL 4)')

! **********************************************************************************************************************************

      END SUBROUTINE CC_NLPARM


      SUBROUTINE CC_SPC ( CARD )

! Processes Case Control SPC cards for defining SPC set ID's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  SPCSETS

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SPC'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! Process SPC cards

! Get SETID

      CALL GET_SETID ( CARD ,SETID )

! Set Case Control variable to SETID

      IF (NSUB /= 0) THEN
         SPCSETS(NSUB) = SETID
      ELSE
         DO I=1,LSUB
            SPCSETS(I) = SETID
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_SPC


      SUBROUTINE CC_STATSUB ( CARD )

! Processes Case Control STATSUB entries for SOL 105 buckling. Syntax accepted:
!   STATSUB           = n
!   STATSUB(PRELOAD)  = n
! The describer (BUCKLING) is recognized and rejected as FATAL since nonlinear preload integration is not implemented.
! "n" is a static SUBCASE id (the external subcase number) whose linear-static solution provides the prestress used to
! build KGGD for the buckling eigenproblem in this subcase. Resolution of "n" against the actual subcase table happens
! later in LOADC, after all SUBCASE cards have been seen.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, BLNK_SUB_NAM, FATAL_ERR, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  CC_STATSUB_DECK, CC_STATSUB_SUB

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_STATSUB'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Case Control card

      INTEGER(LONG)                   :: SETID             ! Integer following '=' on the card (the static subcase id)
      INTEGER(LONG)                   :: LP, RP            ! Positions of '(' and ')' in CARD (0 if absent)
      INTEGER(LONG)                   :: EQ                ! Position of '=' in CARD




! **********************************************************************************************************************************
! Process STATSUB cards

      ! Detect an optional describer ( PRELOAD or BUCKLING ) before the '=' sign.
      EQ = INDEX(CARD,'=')
      LP = INDEX(CARD,'(')
      RP = INDEX(CARD,')')

      IF ((LP > 0) .AND. (RP > LP) .AND. ((EQ == 0) .OR. (LP < EQ))) THEN
         IF (INDEX(CARD(LP+1:RP-1),'BUCKLING') > 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9881)
            WRITE(F06,9881)
            RETURN
         ENDIF
         ! Any describer other than PRELOAD is treated as a warning. Empty describer is silently accepted.
         IF ((INDEX(CARD(LP+1:RP-1),'PRELOAD') == 0) .AND. (LEN_TRIM(CARD(LP+1:RP-1)) > 0)) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9882) CARD(LP:RP)
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9882) CARD(LP:RP)
            ENDIF
         ENDIF
      ENDIF

      ! Pull the integer following '=' (this is a SUBCASE id, not a bulk-data set id, but the parser is identical).
      CALL GET_SETID ( CARD, SETID )

      ! Record the per-subcase or deck-level value. NSUB is incremented by CC_SUBC at parse time, so:
      !   * NSUB == 0  -> this STATSUB appears above any SUBCASE card; it is the deck-default that any subcase lacking
      !                   its own STATSUB inherits during LOADC's post-parse pass.
      !   * NSUB >  0  -> this STATSUB belongs to the current (most-recently-opened) subcase.

      IF (NSUB == 0) THEN

         IF ((CC_STATSUB_DECK /= 0) .AND. (CC_STATSUB_DECK /= SETID)) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9883) CC_STATSUB_DECK, SETID
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9883) CC_STATSUB_DECK, SETID
            ENDIF
         ENDIF
         CC_STATSUB_DECK = SETID

      ELSE

         IF (ALLOCATED(CC_STATSUB_SUB)) THEN
            IF ((CC_STATSUB_SUB(NSUB) /= 0) .AND. (CC_STATSUB_SUB(NSUB) /= SETID)) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,9884) NSUB, CC_STATSUB_SUB(NSUB), SETID
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,9884) NSUB, CC_STATSUB_SUB(NSUB), SETID
               ENDIF
            ENDIF
            CC_STATSUB_SUB(NSUB) = SETID
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 9881 FORMAT(' *ERROR  9881: STATSUB(BUCKLING) IS NOT SUPPORTED. ONLY STATSUB(PRELOAD) (OR THE EQUIVALENT BARE "STATSUB=n")',       &
             ' IS RECOGNIZED.')
 9882 FORMAT(' *WARNING    : UNRECOGNIZED DESCRIBER ',A,' ON STATSUB CASE CONTROL ENTRY. PRELOAD INTERPRETATION ASSUMED.')
 9883 FORMAT(' *WARNING    : MORE THAN ONE DECK-LEVEL STATSUB ENTRY IN CASE CONTROL. PREVIOUS VALUE = ',I8,', NEW VALUE = ',I8,     &
             '. NEW VALUE WILL BE USED AS THE DECK DEFAULT.')
 9884 FORMAT(' *WARNING    : MORE THAN ONE STATSUB ENTRY IN SUBCASE ',I8,'. PREVIOUS VALUE = ',I8,', NEW VALUE = ',I8,              &
             '. NEW VALUE WILL BE USED.')

! **********************************************************************************************************************************

      END SUBROUTINE CC_STATSUB


      SUBROUTINE CC_TEMP ( CARD )

! Processes Case Control TEMP cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  LSUB, NSUB, NTSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SUBLOD

      USE BDF_SET_SYNTAX, ONLY        :  GET_SETID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_TEMP'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
! Process TEMP cards

! Get SETID

      CALL GET_SETID ( CARD, SETID )

! Set Case Control variable to SETID

      NTSUB = NTSUB + 1
      IF (NSUB /= 0) THEN
         SUBLOD(NSUB,2) = SETID
      ELSE
         DO I = 1,LSUB
            SUBLOD(I,2) = SETID
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CC_TEMP

   END MODULE CASE_CONTROL_SELECTORS
