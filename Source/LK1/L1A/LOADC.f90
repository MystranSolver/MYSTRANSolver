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

      SUBROUTINE LOADC

      ! LOADC reads in the CASE CONTROL DECK
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUGOUT, ERR, F06, IN1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, ENFORCED, FATAL_ERR, WARN_ERR, NSUB, NTSUB, PROG_NAME,        &
                                         RESTART, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, MEFFMASS_CALC, MPCSET, MPCSETS, MPFACTOR_CALC, SCNUM, SPCSET, SPCSETS, SUBLOD
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC, FORC_LOC

      USE LOADC_USE_IFs
      USE TO_UPPER_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'LOADC'

      CHARACTER( 1*BYTE)              :: DOLLAR_WARN       ! Indicator of whether there was a $ sign in col 1
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD              ! Case Control card
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'CASE CONTROL'
      CHARACTER(10*BYTE), PARAMETER   :: END_CARD    = 'BEGIN BULK'

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1




! **********************************************************************************************************************************

      ! Process CASE CONTROL DECK

outer:DO
         DOLLAR_WARN = 'N'
         CALL READ_BDF_LINE(IN1, IOCHK, CARD)

         ! Quit if EOF/EOR occurs during read
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs during read.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE outer
         ENDIF

         WRITE(F06,101) CARD

         ! Replace all tab characters with a white space
         CALL REPLACE_TABS_W_BLANKS ( CARD )

         ! Shift card so that it begins in col 1
         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )

         ! Check for CASE CONTROL cards. Exit loop on 'BEGIN BULK'
         IF      (CARD1(1:4) == 'ACCE'    ) THEN
            CALL CC_ACCE ( CARD1 )

         ELSE IF(CARD1(1:10) == 'BEGIN BULK') THEN
            IF (DOLLAR_WARN == 'Y') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1199) CARD1
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1199) CARD1
               ENDIF
            ENDIF
            EXIT outer

         ELSE IF((CARD1(1:4) == 'DISP'    ) .OR.  (CARD1(1:6) == 'VECTOR'  )) THEN
            CALL CC_DISP   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'ECHO'    ) THEN
            CALL CC_ECHO   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'ELDA'    ) THEN
            CALL CC_ELDA   ( CARD1 )
            BUGOUT = 'Y'

         ELSE IF((CARD1(1:7) == 'ELFORCE') .OR. (CARD1(1:4) == 'ELFO').OR. (CARD1(1:5) == 'FORCE')) THEN
            CALL CC_ELFO   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'ENFORCED') THEN
            CALL CC_ENFO   ( CARD1 )
            ENFORCED = 'Y'

         ELSE IF (CARD1(1:4) == 'GPFO'    ) THEN
            CALL CC_GPFO   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'LABE'    ) THEN
            CALL CC_LABE   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'LOAD'    ) THEN
            CALL CC_LOAD   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'MEFFMASS') THEN
            IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
               MEFFMASS_CALC = 'Y'
            ENDIF

         ELSE IF (CARD1(1:4) == 'METH'    ) THEN
            CALL CC_METH   ( CARD1 )

         ELSE IF((CARD1(1:3) == 'MPC'     ) .AND. (CARD1(1:4) /= 'MPCF'    )) THEN
            CALL CC_MPC    ( CARD1 )

         ELSE IF (CARD1(1:4) == 'MPCF'    ) THEN
            CALL CC_MPCF   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'MPFACTOR') THEN
            IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
               MPFACTOR_CALC = 'Y'
            ENDIF

         ELSE IF (CARD1(1:6) == 'NLPARM'  ) THEN
            CALL CC_NLPARM ( CARD1 )

         ELSE IF (CARD1(1:4) == 'OLOA'    ) THEN
            CALL CC_OLOA   ( CARD1 )

         ELSE IF (CARD1(1:6) == 'OUTPUT' ) THEN            ! Normal OUTPUT entry is OK. Ones like OUTPUT(PLOT) we end CC processing.
            IF (INDEX(CARD,"(") > 0) THEN                  ! If we find "(" in an OUTPUT CC entry it indicates, e.g., OUTPUT(PLOT),
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,902) CARD
               WRITE(F06,902) CARD

               ! so read all entries up to BEGIN BULK and then exit loop for CC entries
inner:         DO
                  READ(IN1,101) CARD
                  CARD = TO_UPPER(CARD)
                  IF (CARD(1:10) == 'BEGIN BULK') THEN
                     WRITE(F06,101) CARD
                     EXIT outer
                  ELSE
                     CYCLE inner
                  ENDIF
               ENDDO inner
            ELSE                                           ! The OUTPUT entry was a normal one (with no "(PLOT)", etc delimiter,
               CYCLE outer                                 ! so continue processing CC entries
            ENDIF

         ELSE IF (CARD1(1:3) == 'SET'     ) THEN
            CALL CC_SET    ( CARD1 )

         ELSE IF((CARD1(1:3) == 'SPC'     ) .AND. (CARD1(1:4) /= 'SPCF'    )) THEN
            CALL CC_SPC    ( CARD1 )

         ELSE IF (CARD1(1:4) == 'SPCF'    ) THEN
            CALL CC_SPCF   ( CARD1 )

         ELSE IF((CARD1(1:4) == 'STRA'    ) .OR.  (CARD1(1:4) == 'STRN'    ) .OR.  (CARD1(1:8) == 'ELSTRAIN')) THEN
            CALL CC_STRN   ( CARD1 )

         ELSE IF((CARD1(1:4) == 'STRE'    ) .OR.  (CARD1(1:4) == 'STRS'    ) .OR.  (CARD1(1:8) == 'ELSTRESS')) THEN
            CALL CC_STRE   ( CARD1 )

         ELSE IF (CARD1(1:8) == 'SUBCASE ') THEN
            CALL CC_SUBC   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'SUBT'    ) THEN
            CALL CC_SUBT   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'TEMP'    ) THEN
            CALL CC_TEMP   ( CARD1 )

         ELSE IF (CARD1(1:4) == 'TITL'    ) THEN
            CALL CC_TITL   ( CARD1 )

         ELSE IF (CARD(1:1) == '$') THEN
            DO I=IACHAR('A'),IACHAR('Z')
               IF (CARD(2:2) == ACHAR(I)) THEN
                  DOLLAR_WARN = 'Y'
               ENDIF
            ENDDO

         ELSE IF (CARD(1:CC_ENTRY_LEN) == ' ') THEN
            CYCLE outer

         ELSE                                              ! Card is not recognized.
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,9993) PROG_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9993) PROG_NAME
            ENDIF

         ENDIF

      ENDDO outer

      IF (STRN_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) STRN_LOC, 'STRAIN'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) STRN_LOC, 'STRAIN'
         ENDIF
      ENDIF
      IF (STRE_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) STRE_LOC, 'STRESS'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) STRE_LOC, 'STRESS'
         ENDIF
      ENDIF
      IF (FORC_LOC /= 'CENTER  ') THEN
         WRITE(ERR,1016) FORC_LOC, 'FORCE'
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1016) FORC_LOC, 'FORCE'
         ENDIF
      ENDIF

      IF (NSUB == 0) THEN                                  ! There was no SUBCASE entry in Case Control so set = 1
         NSUB     = 1
         SCNUM(1) = 1
      ENDIF

      ! If SOL is modes or CB or buckilng, then a METH card should have been found in Case Control
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN_CB_MODEL') .OR. (SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF (CC_EIGR_SID == 0) THEN
             WRITE(ERR,1004)
             WRITE(F06,1004)
             FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

     ! If SOL is buckling there should be 2 subcases; the first with a load
     ! and the second with a METH (METH checked above)
     ! If any load, SPC or MPC is found in the 2nd subcase, make sure it is
     ! the same as in subcase 1
      IF (SOL_NAME == 'BUCKLING') THEN
         IF (NSUB /= 2) THEN                               ! Check for 2 subcases
            WRITE(ERR,1101)
            WRITE(F06,1101)
            FATAL_ERR = FATAL_ERR + 1
            IF (NSUB < 2) THEN
               ! further code will crash if we continue with just one subcase
               CALL OUTA_HERE ( 'Y' )
            END IF
         ENDIF
         ! Check that subcase 1 has a mechanical or thermal load
         IF ((SUBLOD(1,1) == 0) .AND. (SUBLOD(1,2) == 0)) THEN
            WRITE(ERR,1101)
            WRITE(F06,1101)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF

         ! Check the 2 subcases for identical loading (if S/C 2 has any stated)
         DO I = 1,2
            IF (SUBLOD(2,I) /= 0) THEN
               IF (SUBLOD(2,I) /= SUBLOD(1,I)) THEN
                  WRITE(ERR,1102)
                  WRITE(F06,1102)
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDIF
         ENDDO

         IF (SPCSETS(2) /= 0) THEN                          ! Check the 2 subcases for identical SPC (if S/C 2 has any stated)
            IF (SPCSETS(2) /= SPCSETS(1)) THEN
               WRITE(ERR,1102)
               WRITE(F06,1102)
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         IF (MPCSETS(2) /= 0) THEN                          ! Check the 2 subcases for identical MPC (if S/C 2 has any stated)
            IF (MPCSETS(2) /= MPCSETS(1)) THEN
               WRITE(ERR,1102)
               WRITE(F06,1102)
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

      ENDIF

      ! Make sure that NTSUB <= NSUB (if user had more TEMP cards in C.C.
      ! than subcases, this will give problems with size of TCASE2)
      IF (NTSUB > NSUB) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1803) NTSUB, NSUB
         WRITE(F06,1803) NTSUB, NSUB
      ENDIF


     ! Make sure that, if there are more than 1 of SPC or MPC set requests
     ! in Case Control, that all (nonzero) are the same. Otherwise issue error.
      SPCSET = SPCSETS(1)
      DO I=2,NSUB
         IF (SPCSETS(I) /= 0) THEN
            IF (SPCSETS(I) /= SPCSET) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1830) 'SPC', (SPCSETS(J),J=1,NSUB)
               WRITE(F06,1830) 'SPC', (SPCSETS(J),J=1,NSUB)
            ENDIF
         ENDIF
      ENDDO

      ! Make sure that, if there are more than 1 of MPC or MPC set requests in
      ! Case Control, that all (nonzero) are the same. Otherwise issue error.
      MPCSET = MPCSETS(1)
      DO I=2,NSUB
         IF (MPCSETS(I) /= 0) THEN
            IF (MPCSETS(I) /= MPCSET) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1830) 'MPC', (MPCSETS(J),J=1,NSUB)
               WRITE(F06,1830) 'MPC', (MPCSETS(J),J=1,NSUB)
            ENDIF
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  901 FORMAT(' *WARNING    : ELDATA ENTRY NOT ALLOWED IN RESTART')

  902 FORMAT(' *WARNING    : CASE CONTROL ENTRY: ',A                                                                               &
                    ,/,14X,' IS NOT ALLOWED. ALL FOLLOWING CASE CONTROL COMMANDS IGNORED')

 1004 FORMAT(' *ERROR  1004: NO METHOD ENTRY FOUND IN CASE CONTROL DECK. CANNOT RUN REAL EIGENVALUE OR LINEAR BUCKLING ANALYSIS.')

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

 1015 FORMAT(' *WARNING    : GRID POINT FORCE BALANCE ONLY ALLOWED IN ',A,' SOLUTION. ABOVE ENTRY IGNORED')

 1016 FORMAT(' *INFORMATION: ALL SUBCASES WILL USE "',A,'" AS THE LOCATION OF ',A,' OUTPUTS FOR PSHELL QUAD4 ELEMENTS'             &
                    ,/,14X,' SINCE THIS WAS THE FIRST REQUEST, OTHER THAN DEFAULT "CENTER  ", DETECTED')

 1028 FORMAT(' *ERROR  1028: THERE MUST BE 2 SUBCASES FOR LINEAR BUCKLING ANALYSES BUT NSUB = ',I8)

 1101 FORMAT(' *ERROR  1101: FOR BUCKLING ANALYSES THERE MUST BE 2 SUBCASES WITH A LOAD (AND/OR TEMP) DEFINED IN SUBCASE 1')

 1102 FORMAT(' *ERROR  1102: FOR BUCKLING ANALYSES ANY LOAD, SPS OR MPC IN 2nd SUBCASE MUST BE THE SAME AS THOSE IN 1st SUBCASE')

 1199 FORMAT(' *WARNING    : BE CAREFUL WITH LINES THAT BEGIN WITH A $ SIGN IN COL 1 FOLLOWED BY AN UPPER CASE LETTER IN EXEC OR', &
                           ' CASE CONTROL.'                                                                                        &
                    ,/,14X,' THE LINE CAN BE MISINTERPRETED AS A DIRECTIVE FOR THE BANDIT GRID RESEQUENCING ALGORITHM.'            &
                    ,/,14X,' SEE THE BANDIT.PDF FILE INSTALLED WHEN YOU RAN SETUP.EXE TO INSTALL MYSTRAN' &
                    ,/,14X, A)

 1803 FORMAT(' *ERROR  1803: CASE CONTROL FOUND NTSUB = ',I8,' TEMP ENTRIES WHEN THERE ARE ONLY NSUB = ',I8,' SUBCASES.',          &
                           ' NTSUB MUST BE <= NSUB')

 1830 FORMAT(' *ERROR  1830: ONLY ONE ',A,' SET ID IS ALLOWED PER RUN. HOWEVER, IN CASE CONTROL THERE WERE THE FOLLOWING SET',     &
                           ' ID''s FOUND: '                                                                                        &
                           ,/,14X,32767(I8,', '))

 9993 FORMAT(' *WARNING    : PRIOR ENTRY NOT PROCESSED BY ',A)

! **********************************************************************************************************************************

      END SUBROUTINE LOADC
