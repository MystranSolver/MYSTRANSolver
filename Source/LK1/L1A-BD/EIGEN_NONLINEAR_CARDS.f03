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

   MODULE EIGEN_NONLINEAR_CARDS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_EIGR, BD_EIGRL, BD_NLPARM

   CONTAINS

      SUBROUTINE BD_EIGR ( CARD, LARGE_FLD_INP, EIGFND )

! Processes EIGR Bulk Data Cards. Reads and checks data and write data to file LINK1M.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1M
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSUB, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, CC_EIGR_SID_SUB, CC_EIGR_SID_DECK, EIG_PARAMS
      USE MODEL_STUF, ONLY            :  EIG_COMP, EIG_CRIT, EIG_CRIT_DEF, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_METH, EIG_MSGLVL,     &
                                         EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, EIG_SID, EIG_SIGMA,    &
                                         EIG_VECS, MAXMIJ, MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, LEFT_ADJ_BDFLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_L1M

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_EIGR'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(1*BYTE), INTENT(INOUT):: EIGFND            ! ='Y' if this EIGR card is the one called for in Case Control
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHRINP            ! Char data in one field of this entry
      CHARACTER( 1*BYTE)              :: USE_THIS_EIG      ! ='Y' if this is the EIGR meth requested in CC

      INTEGER(LONG)                   :: I_SUB             ! DO loop index over subcases
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      LOGICAL                         :: MATCHES_SCALAR    ! True when this card's SID matches the legacy scalar CC_EIGR_SID;
!                                                            in that case we still do the WRITE_L1M / LSUB-hack path so legacy
!                                                            single-METHOD decks continue to behave identically.
      LOGICAL                         :: MATCHES_PER_SUB   ! True when at least one subcase requested this card's SID via its own
!                                                            METHOD card (or via inheritance from the deck-level METHOD default).
      LOGICAL                         :: SUB_WANTS_THIS    ! Per-subcase test inside the EIG_PARAMS population loop




! **********************************************************************************************************************************
! EIGR Bulk Data Card routine

!  Card 1:

!  Field   Item           Description                     Type
!  -----   ------------   -----------                     ----
!   2      EIG_SID        EIGR set ID                     Integer
!   3      EIG_METH       EIGR method                     Char
!   4      EIG_FRQ1       Lower bound on freq in search   Real
!   5      EIG_FRQ2       Upper bound on freq in search   Real, > EIG_FRQ1
!   6      EIG_N1         1st mode number                 Integer, >= 0, default = 1
!   7      EIG_N2         2nd mode number                 Integer, 1 < EIG_N2 <= EIG_N1
!   8      EIG_VECS       Are eigenvecs requested         Char (def = Y)
!   9      EIG_CRIT       Criteria for ortho check        Real >= 0. or blank (used for orthog. check)
!  Required card 2:

!  Field   Item           Description                     Type
!  -----   ------------   -------------                   ----
!   2      EIG_NORM       Type of eigenvec normalization  Char
!   3      EIG_GRID       Grid to normailze on            Integer
!   4      EIG_COMP       DOF comp to normalize on        Integer
!   5      EIG_SIGMA      Shift eigen                     Real


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      JERR = 0
      USE_THIS_EIG    = 'N'
      MATCHES_SCALAR  = .FALSE.
      MATCHES_PER_SUB = .FALSE.
      CALL I4FLD ( JCARD(2), JF(2), EIG_SID )              ! Read set ID and check if it is one requested in Case Control
      IF (IERRFL(2) == 'N') THEN
         IF (EIG_SID == CC_EIGR_SID) THEN
            IF (EIGFND == 'Y') THEN
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
               WRITE(ERR,1117) JCARD(1),JCARD(2)
               WRITE(F06,1117) JCARD(1),JCARD(2)
            ELSE
               EIGFND = 'Y'
               MATCHES_SCALAR = .TRUE.
               USE_THIS_EIG   = 'Y'
            ENDIF
         ENDIF
         ! Also pick up cards requested by any other modes-subcase via its own METHOD entry, or by deck-default
         ! propagation when a subcase has not declared a METHOD of its own. This lets SOL 103 decks define a
         ! distinct set of modes per subcase. We do not bump EIGFND in this branch -- EIGFND is the legacy guard
         ! against duplicate cards for the *scalar* SID only.
         IF (.NOT. MATCHES_SCALAR) THEN
            IF (ALLOCATED(CC_EIGR_SID_SUB)) THEN
               DO I_SUB = 1, NSUB
                  IF (CC_EIGR_SID_SUB(I_SUB) == EIG_SID) THEN
                     MATCHES_PER_SUB = .TRUE.
                     EXIT
                  ENDIF
                  IF ((CC_EIGR_SID_SUB(I_SUB) == 0) .AND. (EIG_SID == CC_EIGR_SID_DECK) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
                     MATCHES_PER_SUB = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
            IF (MATCHES_PER_SUB) THEN
               USE_THIS_EIG = 'Y'
            ELSE
               RETURN
            ENDIF
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

      CALL LEFT_ADJ_BDFLD ( JCARD(3) )
      CALL CHAR_FLD ( JCARD(3), JF(3), CHRINP )            ! Read METHOD and check for valid entry
      IF ((IERRFL(3) == 'N') .AND. (USE_THIS_EIG == 'Y')) THEN
         IF      (CHRINP(1:4) == 'GIV ') THEN
            EIG_METH = 'GIV'
         ELSE IF (CHRINP(1:4) == 'INV ') THEN
            EIG_METH = 'INV'
         ELSE IF (CHRINP(1:4) == 'MGIV') THEN
            EIG_METH = 'MGIV'
         ELSE
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1164) CHRINP
            WRITE(F06,1164) CHRINP
         ENDIF
      ENDIF

      CALL R8FLD ( JCARD(4), JF(4), EIG_FRQ1 )             ! Read lower  frequency of search range
      CALL R8FLD ( JCARD(5), JF(5), EIG_FRQ2 )             ! Read higher frequency of search range
      CALL I4FLD ( JCARD(6), JF(6), EIG_N1 )               ! Read 1st mode number
      CALL I4FLD ( JCARD(7), JF(7), EIG_N2 )               ! Read 2nd mode number
      CALL LEFT_ADJ_BDFLD ( JCARD(8) )                     ! Read indicator of whether eigenvec output is desired
      IF (JCARD(8)(1:1) == 'N') THEN
         EIG_VECS = 'N'
      ELSE
         EIG_VECS = 'Y'
      ENDIF
      IF (JCARD(9)(1:) /= ' ') THEN
         CALL R8FLD ( JCARD(9), JF(9), EIG_CRIT )          ! Read orthogonality check criteria
      ELSE
         EIG_CRIT = EIG_CRIT_DEF
      ENDIF

! Check that the above data read meets requirements.

      CALL EIGR_DATA_CHECK

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF ((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. &! Increment JERR if there were errors reading any of the data fields
          (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. &
          (IERRFL(6) == 'Y') .OR. (IERRFL(7) == 'Y') .OR. &
          (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         JERR = JERR + 1
      ENDIF

! Second Card only required if user wants other than default renormalization of eigenvectors:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         CALL CHAR_FLD ( JCARD(2), JF(2), EIG_NORM )
         IF ((EIG_NORM == 'MASS    ') .OR. (EIG_NORM == 'MAX     ') .OR. (EIG_NORM == 'POINT   ') .OR. (EIG_NORM == 'NONE    '))THEN
            CONTINUE
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1165) EIG_NORM
            WRITE(F06,1165) EIG_NORM
         ENDIF
         IF (EIG_NORM == 'POINT   ') THEN
            CALL I4FLD ( JCARD(3), JF(3), EIG_GRID )
            CALL I4FLD ( JCARD(4), JF(4), EIG_COMP )
            IF(IERRFL(4) == 'N') THEN
               IF ((EIG_COMP < 1) .OR. (EIG_COMP > 6)) THEN
                  JERR = JERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1118) JF(4), EIG_SID, EIG_COMP
                  WRITE(F06,1118) JF(4), EIG_SID, EIG_COMP
               ENDIF
            ENDIF

            IF ((IERRFL(3) == 'Y') .OR. (IERRFL(4) == 'Y')) THEN
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (EIG_METH == 'INV     ') THEN
            CALL R8FLD ( JCARD(5), JF(5), EIG_SIGMA )
         ELSE
            EIG_SIGMA = ZERO
         ENDIF

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF

! Write eigen extract data to file LINK1M if there were no errors and if this is the set ID requested in Case Control

      IF ((JERR == 0) .AND. (USE_THIS_EIG == 'Y')) THEN

         EIG_LAP_MAT_TYPE = '   '
         EIG_MODE         = 0
         EIG_MSGLVL       = 0
         EIG_NCVFACL      = 0

         NUM_FAIL_CRIT    = 0                              ! Following have not been determined yet but write values to L1M anyway
         MAXMIJ           = ZERO
         MIJ_ROW          = 0
         MIJ_COL          = 0

         ! Populate the per-subcase parameter table for every subcase that asked for this SID. The legacy EIG_*
         ! scalars carry the values for the LAST card BD_EIGR processes; LINK4's modes-subcase loop reads back from
         ! EIG_PARAMS(ISUB) into the scalars before each eigensolver invocation, so the scalars' end-of-parsing
         ! state is irrelevant for the multi-METHOD path.
         IF (ALLOCATED(EIG_PARAMS) .AND. ALLOCATED(CC_EIGR_SID_SUB)) THEN
            DO I_SUB = 1, NSUB
               SUB_WANTS_THIS = .FALSE.
               IF (CC_EIGR_SID_SUB(I_SUB) == EIG_SID) SUB_WANTS_THIS = .TRUE.
               IF ((CC_EIGR_SID_SUB(I_SUB) == 0) .AND. (EIG_SID == CC_EIGR_SID_DECK) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
                  SUB_WANTS_THIS = .TRUE.
               ENDIF
               IF (SUB_WANTS_THIS) THEN
                  EIG_PARAMS(I_SUB)%METHOD           = EIG_METH
                  EIG_PARAMS(I_SUB)%NORM             = EIG_NORM
                  EIG_PARAMS(I_SUB)%LAP_MAT_TYPE     = EIG_LAP_MAT_TYPE
                  EIG_PARAMS(I_SUB)%VECS             = EIG_VECS
                  EIG_PARAMS(I_SUB)%SID              = EIG_SID
                  EIG_PARAMS(I_SUB)%N1               = EIG_N1
                  EIG_PARAMS(I_SUB)%N2               = EIG_N2
                  EIG_PARAMS(I_SUB)%COMP             = EIG_COMP
                  EIG_PARAMS(I_SUB)%GRID             = EIG_GRID
                  EIG_PARAMS(I_SUB)%MODE             = EIG_MODE
                  EIG_PARAMS(I_SUB)%MSGLVL           = EIG_MSGLVL
                  EIG_PARAMS(I_SUB)%NCVFACL          = EIG_NCVFACL
                  EIG_PARAMS(I_SUB)%CRIT             = EIG_CRIT
                  EIG_PARAMS(I_SUB)%FRQ1             = EIG_FRQ1
                  EIG_PARAMS(I_SUB)%FRQ2             = EIG_FRQ2
                  EIG_PARAMS(I_SUB)%SIGMA            = EIG_SIGMA
               ENDIF
            ENDDO
         ENDIF

         IF (MATCHES_SCALAR) THEN
            ! to ensure SCNUM is alloc'd right. #subcases = #eigenvecs (legacy hack; needed for single-METHOD path until
            ! LINK9 is fully switched over to MODE_SUBCASE indexing).
            IF (EIG_N2 > LSUB) THEN
               LSUB             = EIG_N2
            ELSE
               ! no idea what the # of eigenvectors should be for now, let's keep
               ! it large for now. this ought to be fixed someday
               LSUB = 1000
            END IF

            CALL WRITE_L1M
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1118 FORMAT(' *ERROR  1118: DOF COMPONENT NUMBER IN FIELD ',I3,' OF EIGR CONTINUATION ENTRY WITH ID = ',I8,' MUST BE A SINGLE',   &
                           ' DIGIT 1-6'                                                                                            &
                    ,/,14X,' BUT VALUE IS = ',I8)

 1164 FORMAT(' *ERROR  1164: METHOD MUST BE GIV, MGIV, OR INV ON EIGR ENTRY. VALUE IS ',A)

 1165 FORMAT(' *ERROR  1165: NORMALIZATION FACTOR ON EIGR CONTINUATION ENTRY MUST BE MASS, MAX, POINT or NONE. VALUE INPUT IS = ',A)

 1117 FORMAT(' *ERROR  1117: ',A,' ENTRY WITH SET ID = ',A,' IS A DUPLICATE SET ID.')

! *********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EIGR_DATA_CHECK

      IMPLICIT NONE

      CHARACTER( 9*BYTE)              :: EIG_SEARCH_CRIT   ! Search criteria (by freg or mode number)

! **********************************************************************************************************************************
      IF ((JCARD(4)(1:) == ' ') .AND. (JCARD(5)(1:) == ' ') .AND. (JCARD(6)(1:) == ' ') .AND.(JCARD(7)(1:) == ' ')) THEN
         WRITE(ERR,1113) JF(4), JF(7), JCARD(1), JCARD(2)
         WRITE(F06,1113) JF(4), JF(7), JCARD(1), JCARD(2)
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      EIG_SEARCH_CRIT(1:) = ' '

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Search criteria must be mode number so check EIG_N1,2
         EIG_SEARCH_CRIT = 'MODE NMBR'
         IF (JCARD(6)(1:) == ' ') THEN
            EIG_N1 = 1
         ENDIF
         IF (EIG_N1 < 0) THEN
            WRITE(ERR,1103) JF(6), JCARD(1), JCARD(2)
            WRITE(F06,1103) JF(6), JCARD(1), JCARD(2)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         IF (EIG_N2 < EIG_N1) THEN
            WRITE(ERR,1105) JF(6), JF(7), JCARD(1), JCARD(2), JF(7), JF(6)
            WRITE(F06,1105) JF(6), JF(7), JCARD(1), JCARD(2), JF(7), JF(6)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ELSE                                                 ! There is no EIG_N2 so make sure that there is no EIG_N1
         IF (JCARD(6)(1:) /= ' ') THEN
            WRITE(ERR,1106) JCARD(1), JCARD(2)
            WRITE(F06,1106) JCARD(1), JCARD(2)
            FATAL_ERR = FATAL_ERR + 1
            IF (EIG_N1 < 0) THEN
               WRITE(ERR,1103) JF(6), JCARD(1), JCARD(2)
               WRITE(F06,1103) JF(6), JCARD(1), JCARD(2)
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF
      ENDIF

      IF (EIG_SEARCH_CRIT == 'MODE NMBR') THEN
         RETURN
      ELSE                                                 ! Search criteria must be frequency so check EIG_FRQ1,2
         IF (JCARD(4)(1:) == ' ') THEN
            EIG_FRQ1 = ZERO
         ENDIF
         IF (EIG_FRQ1 < ZERO) THEN
            WRITE(ERR,1103) JF(4), JCARD(1), JCARD(2)
            WRITE(F06,1103) JF(4), JCARD(1), JCARD(2)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         IF (EIG_FRQ2 < EIG_FRQ1) THEN
            WRITE(ERR,1105) JF(4), JF(5), JCARD(1), JCARD(2), JF(5), JF(4)
            WRITE(F06,1105) JF(4), JF(5), JCARD(1), JCARD(2), JF(5), JF(4)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1103 FORMAT(' *ERROR  1103: ILLEGAL ENTRY IN FIELD ',I2,' OF ',A,A,' ENTRY. ENTRY MUST BE > 0 OR BLANK')

 1105 FORMAT(' *ERROR  1105: ILLEGAL ENTRY IN FIELD ',I2,' OR ',I2,' OF ',A,A,' ENTRY.',                                           &
                           ' FIELD ',I2,' ENTRY MUST BE GREATER THAN FIELD ',I2,' ENTRY')

 1106 FORMAT(' *ERROR  1106: ILLEGAL ENTRY ON ',A,A,' ENTRY. CANNOT HAVE N1 DEFINED IN FIELD 6 WITH FIELD 7 (N2) BLANK')

 1113 FORMAT(' *ERROR  1113: FIELDS ',I2,' -',I2,' OF ',A,A,' ENTRY CANNOT ALL BE BLANK. THERE IS NO DEFINITION FOR A SEARCH RANGE')

! **********************************************************************************************************************************

      END SUBROUTINE EIGR_DATA_CHECK

      END SUBROUTINE BD_EIGR


      SUBROUTINE BD_EIGRL ( CARD, LARGE_FLD_INP, EIGFND )

! Processes EIGRL Bulk Data Cards. Reads and checks data and write data to file LINK1M.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1M
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSUB, NSUB, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPM4
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, CC_EIGR_SID_SUB, CC_EIGR_SID_DECK, EIG_PARAMS,                               &
                                         EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_LANCZOS_NEV_DELT,                   &
                                         EIG_METH, EIG_MSGLVL, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM,  &
                                         EIG_SID, EIG_SIGMA, EIG_VECS, MAXMIJ, MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, R8FLD
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_L1M

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_EIGRL'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(1*BYTE), INTENT(INOUT):: EIGFND            ! ='Y' if this EIGR card is the one called for in Case Control
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 1*BYTE)              :: USE_THIS_EIG      ! ='Y' if this is the EIGR meth requested in CC
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I4INP             ! An integer*4 value read
      INTEGER(LONG)                   :: I_SUB             ! DO loop index over subcases
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      LOGICAL                         :: MATCHES_SCALAR    ! Mirrors BD_EIGR: legacy scalar match drives WRITE_L1M / LSUB hack
      LOGICAL                         :: MATCHES_PER_SUB   ! True when any modes-subcase requested this card
      LOGICAL                         :: SUB_WANTS_THIS    ! Per-subcase test inside the EIG_PARAMS population loop




! **********************************************************************************************************************************
! EIGRL Bulk Data Card routine

!  Card 1:

!  Field   Item                   Description                          Type
!  -----   ------------           -----------                          ----
!   2      EIG_SID               EIGRL set ID                         Integer
!   3      EIG_FRQ1              Lower bound on freq in search        Real   , >= 0.
!   4      EIG_FRQ2              Upper bound on freq in search        Real   , >= 0.
!   5      EIG_N2                Desired number of roots              Integer, >= 0
!   6      EIG_MSGLVL            Message level for subr Lanczos       Integer, >= 0
!   7      EIG_NCVFACL           (see MODEL_STUF for explanation)     Integer, >= 1
!   8      EIG_SIGMA             Lanczos shift eigen                  Real
!   9      EIG_NORM              Renormalization method (MASS or MAX) Char

! Continuation entry:

!   1      EIG_MODE              Lanczos "mode" (dsband)              Integer, 2 or 3
!   2      EIG_LAP_MAT_TYPE      LAPACK matrix type (DGB, DPB)        Char
!   3      EIG_LANCZOS_NEV_DELT  Number to add to est num roots       Integer >= 0

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      JERR = 0
      USE_THIS_EIG    = 'N'
      MATCHES_SCALAR  = .FALSE.
      MATCHES_PER_SUB = .FALSE.

      ! second card deprecated. set defaults:
      !   - ARPACK mode 2 for buckling, 3 for everything else
      !   - DGB matrix type (in case we use the banded solver)
      !   - EIG_LANCZOS_NEV_DELT, previously undocumented, kept default (2)
      IF (SOL_NAME == 'BUCKLING') THEN
         EIG_MODE = 2
      ELSE
         EIG_MODE = 3
      ENDIF
      EIG_LAP_MAT_TYPE = 'DGB     '
      EIG_LANCZOS_NEV_DELT = 2

      CALL I4FLD ( JCARD(2), JF(2), EIG_SID )              ! Read set ID and check if it is one requested in Case Control
      IF (IERRFL(2) == 'N') THEN
         IF (EIG_SID == CC_EIGR_SID) THEN
            IF (EIGFND == 'Y') THEN
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
               WRITE(ERR,1117) JCARD(1),JCARD(2)
               WRITE(F06,1117) JCARD(1),JCARD(2)
            ELSE
               EIGFND = 'Y'
               MATCHES_SCALAR = .TRUE.
               USE_THIS_EIG   = 'Y'
            ENDIF
         ENDIF
         IF (.NOT. MATCHES_SCALAR) THEN
            ! See BD_EIGR for rationale: also pick up cards needed by per-subcase METHOD requests or by
            ! deck-default propagation. Only the scalar match triggers the WRITE_L1M / LSUB-hack path.
            IF (ALLOCATED(CC_EIGR_SID_SUB)) THEN
               DO I_SUB = 1, NSUB
                  IF (CC_EIGR_SID_SUB(I_SUB) == EIG_SID) THEN
                     MATCHES_PER_SUB = .TRUE.
                     EXIT
                  ENDIF
                  IF ((CC_EIGR_SID_SUB(I_SUB) == 0) .AND. (EIG_SID == CC_EIGR_SID_DECK) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
                     MATCHES_PER_SUB = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
            IF (MATCHES_PER_SUB) THEN
               USE_THIS_EIG = 'Y'
            ELSE
               RETURN
            ENDIF
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

      CALL R8FLD ( JCARD(3), JF(3), EIG_FRQ1 )             ! Read field 3: lower  frequency of search range

      CALL R8FLD ( JCARD(4), JF(4), EIG_FRQ2 )             ! Read field 4: higher frequency of search range

      ! default N should be 1, not 0
      EIG_N2 = 1
      IF (JCARD(5)(1:) /= ' ') THEN
         CALL I4FLD ( JCARD(5), JF(5), I4INP )            ! Read field 5: number of desired roots
         IF (IERRFL(5) == 'N') THEN
            EIG_N2 = I4INP
         ENDIF
      ENDIF

      IF (JCARD(6)(1:) /= ' ') THEN                        ! Read field 6: MSGLVL
         CALL I4FLD ( JCARD(6), JF(6), I4INP )
         IF (IERRFL(6) == 'N') THEN
            EIG_MSGLVL = I4INP
         ENDIF
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Read field 7: EIG_NCVFACL
         CALL I4FLD ( JCARD(7), JF(7), I4INP )
         IF (IERRFL(7) == 'N') THEN
            EIG_NCVFACL = I4INP
            IF (JCARD(7)(1:) /= ' ') THEN
               IF (EIG_MODE < 1) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1137) JF(7), EIG_NCVFACL
                  WRITE(F06,1137) JF(7), EIG_NCVFACL
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (JCARD(8)(1:) /= ' ') THEN                        ! Read field 8: Lanczos shift freq
         CALL R8FLD ( JCARD(8), JF(8), EIG_SIGMA )
      ENDIF

      IF (JCARD(9)(1:) /= ' ') THEN                        ! Read field 9: renormalization method if field is not blank
      CALL CHAR_FLD ( JCARD(9), JF(9), EIG_NORM )
         IF ((EIG_NORM == 'MASS    ') .OR. (EIG_NORM == 'MAX     ') .OR. (EIG_NORM == 'NONE    ')) THEN
            CONTINUE
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1107) EIG_NORM
            WRITE(F06,1107) EIG_NORM
         ENDIF
      ENDIF

! Check that the above data read meets requirements.

      CALL EIGRL_DATA_CHECK

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. & ! Increment JERR if there were errors reading any of the data fields
         (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. &
         (IERRFL(6) == 'Y') .OR. (IERRFL(7) == 'Y') .OR. &
         (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         JERR = JERR + 1
      ENDIF

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      EIG_CRIT = ONEPM4                                    ! Use this until code is changed to read a value from the EIGRL entry
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! Write eigen extract data to file LINK1M if there were no errors and if this is the set ID requested in Case Control

      IF ((JERR == 0) .AND. (USE_THIS_EIG == 'Y')) THEN

         EIG_METH      = 'LANCZOS'
         EIG_N1        = 1
         EIG_GRID      = 0
         EIG_COMP      = 0
         EIG_VECS      = 'Y'

         NUM_FAIL_CRIT = 0                                 ! Following have not been determined yet but write values to L1M anyway
         MAXMIJ        = ZERO
         MIJ_ROW       = 0
         MIJ_COL       = 0

         ! Capture per-subcase parameters (see BD_EIGR for the matching rationale).
         IF (ALLOCATED(EIG_PARAMS) .AND. ALLOCATED(CC_EIGR_SID_SUB)) THEN
            DO I_SUB = 1, NSUB
               SUB_WANTS_THIS = .FALSE.
               IF (CC_EIGR_SID_SUB(I_SUB) == EIG_SID) SUB_WANTS_THIS = .TRUE.
               IF ((CC_EIGR_SID_SUB(I_SUB) == 0) .AND. (EIG_SID == CC_EIGR_SID_DECK) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
                  SUB_WANTS_THIS = .TRUE.
               ENDIF
               IF (SUB_WANTS_THIS) THEN
                  EIG_PARAMS(I_SUB)%METHOD            = EIG_METH
                  EIG_PARAMS(I_SUB)%NORM              = EIG_NORM
                  EIG_PARAMS(I_SUB)%LAP_MAT_TYPE      = EIG_LAP_MAT_TYPE
                  EIG_PARAMS(I_SUB)%VECS              = EIG_VECS
                  EIG_PARAMS(I_SUB)%SID               = EIG_SID
                  EIG_PARAMS(I_SUB)%N1                = EIG_N1
                  EIG_PARAMS(I_SUB)%N2                = EIG_N2
                  EIG_PARAMS(I_SUB)%COMP              = EIG_COMP
                  EIG_PARAMS(I_SUB)%GRID              = EIG_GRID
                  EIG_PARAMS(I_SUB)%LANCZOS_NEV_DELT  = EIG_LANCZOS_NEV_DELT
                  EIG_PARAMS(I_SUB)%MODE              = EIG_MODE
                  EIG_PARAMS(I_SUB)%MSGLVL            = EIG_MSGLVL
                  EIG_PARAMS(I_SUB)%NCVFACL           = EIG_NCVFACL
                  EIG_PARAMS(I_SUB)%CRIT              = EIG_CRIT
                  EIG_PARAMS(I_SUB)%FRQ1              = EIG_FRQ1
                  EIG_PARAMS(I_SUB)%FRQ2              = EIG_FRQ2
                  EIG_PARAMS(I_SUB)%SIGMA             = EIG_SIGMA
               ENDIF
            ENDDO
         ENDIF

         IF (MATCHES_SCALAR) THEN
            ! ensure a proper size for SCNUM (legacy single-METHOD path)
            IF (EIG_N2 > LSUB) THEN
               LSUB          = EIG_N2
            ELSE
               ! since we have adaptive lanczos now, we set this to be
               ! INITIAL_NEV*(2**MAX_DOUBLINGS), both being 10 and unlikely to be
               ! changed unless someone *really* wants more than 10k modes AND
               ! doesn't want to specify nmodes manually.
               IF (SOL_NAME /= 'BUCKLING') THEN
                  LSUB = 10240
               END IF
            END IF

            CALL WRITE_L1M
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1107 FORMAT(' *ERROR  1107: NORMALIZATION FACTOR ON EIGRL ENTRY MUST BE MASS, MAX or NONE. VALUE INPUT IS = ',A)

 1108 FORMAT(' *ERROR  1108: LANCZOS SOLUTION "MODE" IN FIELD ',I2,' MUST BE 2 OR 3. VALUE INPUT IS = ',I8)

 1109 FORMAT(' *ERROR  1109: LAPACK MATRIX TYPE SHOULD BE "DGB" OR "DPB". VALUE INPUT IS = ',A)

 1117 FORMAT(' *ERROR  1117: ',A,' ENTRY WITH SET ID = ',A,' IS A DUPLICATE SET ID.')

 1137 FORMAT(' *ERROR  1137: FIELD ',I2,' ON THE EIGRL ENTRY MUST BE >= 1 BUT WAS = ',I8)

! *********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EIGRL_DATA_CHECK

      IMPLICIT NONE

! **********************************************************************************************************************************
      IF ((JCARD(3)(1:) == ' ') .AND. (JCARD(4)(1:) == ' ') .AND. (JCARD(5)(1:) == ' ')) THEN
         WRITE(ERR,1113) JF(3), JF(5), 'EIGRL', EIG_SID
         WRITE(F06,1113) JF(3), JF(5), 'EIGRL', EIG_SID
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      IF (JCARD(5)(1:) /= ' ') THEN                        ! Search criteria must be mode number so check EIG_N1,2
         IF (EIG_N2 < 0) THEN
            WRITE(ERR,1103) JF(5), 'EIGRL', EIG_SID
            WRITE(F06,1103) JF(5), 'EIGRL', EIG_SID
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         RETURN

      ELSE                                                 ! Search criteria must be frequency so check EIG_FRQ1,2

         IF (JCARD(3)(1:) == ' ') THEN
            EIG_FRQ1 = ZERO
         ENDIF
         IF (EIG_FRQ1 < ZERO) THEN
            WRITE(ERR,1103) JF(3), 'EIGRL', EIG_SID
            WRITE(F06,1103) JF(3), 'EIGRL', EIG_SID
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         IF (EIG_FRQ2 < EIG_FRQ1) THEN
            WRITE(ERR,1105) JF(3), JF(4), 'EIGRL', EIG_SID, JF(4), JF(3)
            WRITE(F06,1105) JF(3), JF(4), 'EIGRL', EIG_SID, JF(4), JF(3)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         IF (DABS(EIG_FRQ2 - EIG_FRQ1) == ZERO) THEN
            WRITE(ERR,1111) JF(3), JF(4), 'EIGRL', EIG_SID
            WRITE(F06,1111) JF(3), JF(4), 'EIGRL', EIG_SID
            FATAL_ERR = FATAL_ERR + 1
         ENDIF

      ENDIF

      RETURN

! **********************************************************************************************************************************
 1113 FORMAT(' *ERROR  1113: FIELDS ',I2,' -',I2,' OF ',A,I8,' CANNOT ALL BE BLANK. THERE IS NO DEFINITION FOR AN EIGENVAL ',      &
                            'SEARCH RANGE')

 1103 FORMAT(' *ERROR  1103: ILLEGAL ENTRY IN FIELD ',I2,' OF ',A,I8,'. ENTRY MUST BE > 0 OR BLANK')

 1105 FORMAT(' *ERROR  1105: ILLEGAL ENTRY IN FIELD ',I2,' OR ',I2,' OF ',A,I8,' ENTRY.',                                          &
                           ' FIELD ',I2,' ENTRY MUST BE GREATER THAN FIELD ',I2,' ENTRY')

 1111 FORMAT(' *ERROR  1111: ILLEGAL ENTRY IN FIELD ',I2,' OR ',I2,' OF ',A,I8,' ENTRY.',                                          &
                           ' FREQ RANGE, |F2 - F1| MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE EIGRL_DATA_CHECK

      END SUBROUTINE BD_EIGRL


      SUBROUTINE BD_NLPARM ( CARD, CC_NLSID_FND )

! Processes NLPARM Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSUB
      USE TIMDAT, ONLY                :  TSEC
      USE NONLINEAR_PARAMS, ONLY      :  NL_MAXITER, NL_NUM_LOAD_STEPS, NL_SID

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_NLPARM'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_NLSID_FND(LSUB)! 'Y' if B.D NLPARM card w/ same set ID (SID) as C.C. NLPARM = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CARD_NAME         ! Field 1 of CARD
      CHARACTER(LEN(JCARD))           :: CHAR_SID          ! Field 2 of CARD

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! An integer value read from GRAV entry
      INTEGER(LONG)                   :: SETID             ! NLPARM set id




! **********************************************************************************************************************************
! NLPARM Bulk Data Card routine

!   FIELD   ITEM                VARIABLE
!   -----   ------------        -------------
!     2     NL_SID              Set ID
!     3     NL_NUM_LOAD_STEPS   Number of increments into which the load is to be subdivided
!     7     NL_MAXITER          Maximum number of iterations for convergence in any 1 load step


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Initialize

      NL_NUM_LOAD_STEPS = 1

! Read and check data

      CARD_NAME = JCARD(1)
      CHAR_SID  = JCARD(2)

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,LSUB
            IF (SETID == NL_SID(I)) THEN
               CC_NLSID_FND(I) = 'Y'
            ENDIF
         ENDDO
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), I4INP )
      IF (IERRFL(3) == 'N') THEN
         IF (I4INP >= 0) THEN
            NL_NUM_LOAD_STEPS = I4INP
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1156) JF(3), CARD_NAME, CHAR_SID, I4INP
            WRITE(F06,1156) JF(3), CARD_NAME, CHAR_SID, I4INP
         ENDIF
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN
         CALL I4FLD ( JCARD(7), JF(7), I4INP )
         IF (IERRFL(7) == 'N') THEN
            IF (I4INP >= 0) THEN
               NL_MAXITER = I4INP
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1156) JF(7), CARD_NAME, CHAR_SID, I4INP
               WRITE(F06,1156) JF(7), CARD_NAME, CHAR_SID, I4INP
            ENDIF
         ENDIF
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,0,0,7,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-3
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,0,8,9 )   ! Issue warning if fields 4,5,6,8,9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1156 FORMAT(' *ERROR  1156: ILLEGAL ENTRY IN FIELD ',I2,' OF ',A,A,' CARD. ENTRY MUST BE > 0 BUT WAS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_NLPARM

   END MODULE EIGEN_NONLINEAR_CARDS
