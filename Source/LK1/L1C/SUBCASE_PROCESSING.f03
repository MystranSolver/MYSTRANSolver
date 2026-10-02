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

   MODULE SUBCASE_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SUBCASE_PROC

   CONTAINS

      SUBROUTINE SUBCASE_PROC

! Subcase processor

! Processes all output requests from Case Control. When output requests were read in Case Control, they were temporarily stored into
! the arrays (each array is NSUB x 1):

!         SC_ACCE for displacement output requests
!         SC_DISP for displacement output requests
!         SC_ELFN for elem nodal force output requests
!         SC_ELFE for elem engineering force output requests
!         SC_GPFO for G.P. force balance output requests
!         SC_MPCF for MPC force output requests
!         SC_OLOA for applied load output requests
!         SC_SPCF for SPC force output requests
!         SC_STRE for elem stress output requests
!         SC_STRN for elem strain output requests

! Each of the NSUB entries in each of these arrays is one of the following:
!
!          0   : means no output was requested, or
!         -1   : means output for all (grids or elems) was requested, or
!         SETID: means output for all members of set SETID was requested

! This subr processes this information, along with the set definitions that are contained in array ALL_SETS_ARRAY,
! to create arrays that specify for every grid, element, and subcase the specific output requested. This specification
! is put into arrays OGROUT and GROUT (for grid related outputs), OELOUT and ELOUT (for element related outputs)
! (for element force and stress), and ELDT (for element debug output):

!  (1) ACCE, DISP, GPFO, MPCF, OLOA, SPCF output requests are put into arrays OGROUT and  GROUT
!      (a) OGROUT is NSUB  x 1    with an indicator for every S/C as to whether any grid related outputs were requested
!      (b) GROUT  is NGRID x NSUB with an indicator for every G.P.-S/C for which grid related output was requested

!  (2) ELFN, ELFE, STRE,STRN output requests are put into arrays OELOUT and  ELOUT
!      (a) OELOUT is NSUB x 1    with an indicator for every S/C as to whether any element related output was requested
!      (b) ELOUT  is NELE x NSUB with an indicator for every elem-S/C for which element related outputs were requested

!  (3) Element debug (ELDATA Case Control requests) are put into array ELDT (these requests are above subcase level)
!      (a) ELDT   is NELE x 1    with an indicator for every element for which element debug output was requested


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1D

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, DATA_NAM_LEN, FATAL_ERR, IBIT, WARN_ERR, LSETLN,              &
                                         MELDTS, MELOUTS, METYPE, MGROUTS, NELE, NGRID, NSUB

      USE SCONTR, ONLY                :  GROUT_ACCE_BIT, GROUT_DISP_BIT, GROUT_GPFO_BIT, GROUT_MPCF_BIT, GROUT_OLOA_BIT,           &
                                         GROUT_SPCF_BIT, ELOUT_ELFE_BIT, ELOUT_ELFN_BIT, ELOUT_STRE_BIT, ELOUT_STRN_BIT

      USE PARAMS, ONLY                :  PRTSCP, SUPWARN
      USE TIMDAT, ONLY                :  TSEC

      USE MODEL_STUF, ONLY            :  CCELDT, ONE_SET_ARRAY, SC_ACCE, SC_DISP, SC_ELFN, SC_ELFE, SC_GPFO, SC_MPCF,              &
                                         SC_OLOA, SC_SPCF, SC_STRE, SC_STRN, ELDT, OELDT, ELOUT, OELOUT, GROUT, OGROUT, LABEL,     &
                                         SCNUM, STITLE, TITLE, SUBLOD, GRID, GRID_ID, ESORT1, ETYPE

      USE MODEL_STUF, ONLY            :  ANY_ACCE_OUTPUT, ANY_DISP_OUTPUT, ANY_MPCF_OUTPUT, ANY_SPCF_OUTPUT, ANY_OLOA_OUTPUT,      &
                                         ANY_GPFO_OUTPUT, ANY_ELFE_OUTPUT, ANY_ELFN_OUTPUT, ANY_STRE_OUTPUT, ANY_STRN_OUTPUT
      USE BDF_SET_SYNTAX, ONLY        :  STOKEN, TOKCHK
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SUBCASE_PROC'
      CHARACTER(     1*BYTE)          :: CBIT(16)          ! Flag to indicate if a bit in GROUT, ELOUT, ELDT is turned on
      CHARACTER(20*BYTE)              :: ELM_OUTPUT_REQ    ! Type of elem output request (stress, etc)
      CHARACTER(10*BYTE)              :: ELM_BIT_NAME      ! Char name for output warning purposes
      CHARACTER(LEN(ETYPE))           :: ELM_TYP(METYPE)   ! List of elem types requested for elem output but not allowed
      CHARACTER(LSETLN*BYTE)          :: ERRTOK            ! Error message for printout (returned from subr STOKEN)
      CHARACTER(     3*BYTE)          :: EXCEPT            ! Flag indicating whether EXCEPT is "ON "/"OFF" (INOUT for subr STOKEN)
      CHARACTER(DATA_NAM_LEN*BYTE)    :: DATA_SET_NAME     ! A data set name for output purposes
      CHARACTER(     8*BYTE)          :: WARN_NAME         ! Name for output error message purposes
      CHARACTER(     1*BYTE)          :: PRNTOUT           ! Flag used in deciding what to output if B.D. PARAM PRTSCP = 1
      CHARACTER(     3*BYTE)          :: THRU              ! Flag indicating whether THRU is "ON " or "OFF" (INOUT for subr STOKEN)
      CHARACTER(     8*BYTE)          :: TOKEN(3)          ! Array of 3 char tokens returned from subr STOKEN
      CHARACTER(     8*BYTE)          :: TOKTYP(3)         ! Array of 3 char indicators of what type of tokens are in TOKEN(1-3)
      CHARACTER(LSETLN*BYTE)          :: TOKSTR            ! Character string to tokenize. This is ONE_SET_ARRAY which is generated
!                                                            from ALL_SETS_ARRAY by a call to subr SETPRO (CONTAIN'd herein). We
!                                                            use TOKSTR instead of ONE_SET_ARRAY here because ONE_SET_ARRAY is an
!                                                            LSETLN array of 1 byte characters and we want TOKSTR to be a single
!                                                            character variable of length LSETLN)

      INTEGER(LONG)                   :: AELEM             ! Actual elem ID
      INTEGER(LONG)                   :: AELEM_LO          ! Lower  elem ID (actual) in a range of elem ID's
      INTEGER(LONG)                   :: AELEM_HI          ! Higher elem ID (actual) in a range of elem ID's
      INTEGER(LONG)                   :: AGRID             ! Actual grid ID
      INTEGER(LONG)                   :: AGRID_LO          ! Lower  grid ID (actual) in a range of grid ID's
      INTEGER(LONG)                   :: AGRID_HI          ! Higher grid ID (actual) in a range of grid ID's
      INTEGER(LONG)                   :: ELM_BIT(METYPE)   ! Array used for output warning purposes
      INTEGER(LONG)                   :: ESORT1_ROW_NUM    ! Row number in array ESORT1 where an actual element ID is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: I,J,K,M           ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Intermediate variable when setting/clearing bits in arrays
      INTEGER(LONG)                   :: IERROR        = 0 ! Error number from subr STOKEN
      INTEGER(LONG)                   :: INUM              ! No. of elems or grids found in range xxxxx_LO to xxxxx_HI
      INTEGER(LONG)                   :: IOCHK1            ! IOSTAT error number when reading a file
      INTEGER(LONG)                   :: IOCHK2            ! IOSTAT error number when reading a file
      INTEGER(LONG)                   :: ISTART        = 0 ! Where to start looking for a SET in TOKSTR (INOUT for subr STOKEN)
      INTEGER(LONG)                   :: JERR              ! Local error count
      INTEGER(LONG)                   :: NTOKEN        = 0 ! The number of tokens found in subr STOKEN
      INTEGER(LONG)                   :: NUM_WARN      = 0 ! Count of elem types for which elem output was requested but not allowed
      INTEGER(LONG)                   :: NULSET            ! Output from subr SETPRO. If 0, there were no SET's in ALL_SETS_ARRAY
      INTEGER(LONG)                   :: SETID             ! = 0, -1, or pos integer set ID read from array SC_xxxx
      INTEGER(LONG)                   :: TOKLEN        = 0 ! Length (bytes) of TOKSTR sent to subr STOKEN


      INTRINSIC                       :: IAND,IBCLR,IBSET



! **********************************************************************************************************************************
      JERR = 0

! **********************************************************************************************************************************
! 1st: Process grid output requests to generate OGROUT(i), GROUT(j,i)

! OGROUT(i) is an overall indicator for each S/C for whether any of the types of G.P. output requests were made in
! Case Control.

! GROUT(j,i) is a specific indicator for each S/C for every G.P. of the types of G.P. output requests were made in
! Case Control.

! Both OGROUT and GROUT are set here to be used in LINK9.

! Bit position 0 in OGROUT(i), GROUT(i) is for DISP requests (turned on if SC_DISP(i) > 0 which is set in subr CC_DISP).
! Bit position 1 in OGROUT(i), GROUT(i) is for OLOA requests (turned on if SC_OLOA(i) > 0 which is set in subr CC_OLOA).
! Bit position 2 in OGROUT(i), GROUT(i) is for SPCF requests (turned on if SC_SPCF(i) > 0 which is set in subr CC_SPCF).
! Bit position 3 in OGROUT(i), GROUT(i) is for MPCF requests (turned on if SC_MPCF(i) > 0 which is set in subr CC_MPCF).
! Bit position 4 in OGROUT(i), GROUT(i) is for GPFO requests (turned on if SC_GPFO(i) > 0 which is set in subr CC_GPFO).
! Bit position 5 in OGROUT(i), GROUT(i) is for ACCE requests (turned on if SC_ACCE(i) > 0 which is set in subr CC_ACCE).

! SC_DISP(i), etc. are integer arrays containing SET ID's (-1 for  'ALL', 0 for 'NONE', positive integer no. for
! SETID) of the 4  possible G.P. related outputs requested in Case Control

      WARN_NAME = 'GRID    '

! I ranges over the number of types of G.P. output requests that can be made in Case Control (described above)

      DO I=0,MGROUTS-1

j_loop1: DO J=1,NSUB

            IF      (I == 0) THEN
               SETID = SC_DISP(J)
            ELSE IF (I == 1) THEN
               SETID = SC_OLOA(J)
            ELSE IF (I == 2) THEN
               SETID = SC_SPCF(J)
            ELSE IF (I == 3) THEN
               SETID = SC_MPCF(J)
            ELSE IF (I == 4) THEN
               SETID = SC_GPFO(J)
            ELSE IF (I == 5) THEN
               SETID = SC_ACCE(J)
            ENDIF

            IF (SETID ==  0) THEN                          ! Check for 'NONE'
               CYCLE j_loop1
            ENDIF

            IF (SETID == -1) THEN                          ! Check for 'ALL' (SETID = -1)
               OGROUT(J) = IBSET(OGROUT(J),I)              ! Set bit position I in OGROUT to 1
               DO K = 1,NGRID
                  GROUT(K,J) = IBSET(GROUT(K,J),I)
               ENDDO
               CYCLE j_loop1
            ENDIF

            IF (SETID > 0) THEN                            ! Check for actual set ID

               CALL SETPRO ( SETID, NULSET, TOKLEN )       ! Call SETPRO to get the data which was input for set SETID.
               DO K=1,LSETLN
                  TOKSTR(K:K) = ONE_SET_ARRAY(K)
               ENDDO

               IF (NULSET == 0) THEN                       ! No SET's in array ALL_SETS_ARRAY (error was already written in SETPRO)
                  JERR = JERR + 1
                  CYCLE j_loop1
               ENDIF

               ISTART = 1                                  ! Init ISTART, THRU, EXCEPT for the DO loop looking for tokens in TOKSTR
               THRU   = 'OFF'
               EXCEPT = 'OFF'

token_loop1:   DO                                          ! Call STOKEN in a DO loop until we run out of tokens in TOKSTR

                  CALL STOKEN ( SUBR_NAME, TOKSTR, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

                  IF (IERROR == 1) THEN                    ! Error: too long a token in the set definition
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1421) ERRTOK
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1421) ERRTOK
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop1

                  ELSE IF (IERROR == 2) THEN               ! Error: "I1 THRU I2" is missing I2 in the set definition
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1425) NTOKEN
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1425) NTOKEN
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop1

                  ELSE IF (IERROR == 3) THEN               ! Error: found "EXCEPT" but EXCEPT has already been turned "ON "
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1428)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1428)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop1

                  ELSE IF (IERROR == 4) THEN               ! Error: found "EXCEPT" but "THRU" has not been found yet
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1429)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1429)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop1

                  ELSE IF (NTOKEN == 1) THEN               ! Is the TOKEN a number not connected to a 'THRU'?
                     IF (TOKTYP(1) == 'INTEGER ') THEN
                        OGROUT(J) = IBSET(OGROUT(J),I)
                        READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AGRID
                        IF (IOCHK1 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop1
                        ELSE
                           CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
                           IF (GRID_ID_ROW_NUM == -1) THEN
                              WARN_ERR = WARN_ERR + 1
                              WRITE(ERR,1402) WARN_NAME,AGRID,SETID
                              IF (SUPWARN == 'N') THEN
                                 WRITE(F06,1402) WARN_NAME,AGRID,SETID
                              ENDIF
                           ELSE
                              IF (EXCEPT == 'OFF') THEN
                                 IDUM           = GROUT(GRID_ID_ROW_NUM,J)
                                 GROUT(GRID_ID_ROW_NUM,J) = IBSET(IDUM,I)
                              ENDIF
                              IF (EXCEPT == 'ON ') THEN
                                 IF ((AGRID >= AGRID_LO) .AND. (AGRID <= AGRID_HI)) THEN
                                    IDUM           = GROUT(GRID_ID_ROW_NUM,J)
                                    GROUT(GRID_ID_ROW_NUM,J) = IBCLR(IDUM,I)
                                 ELSE
                                    IDUM           = GROUT(GRID_ID_ROW_NUM,J)
                                    GROUT(GRID_ID_ROW_NUM,J) = IBSET(IDUM,I)
                                    EXCEPT = 'OFF'
                                    THRU   = 'OFF'
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF

                     ELSE IF (TOKTYP(1) == 'EXCEPT  ') THEN
                        IF (THRU /= 'ON ') THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1427)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1427)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop1
                        ENDIF

                     ELSE
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1422) TOKEN(1),TOKTYP(1)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1422) TOKEN(1),TOKTYP(1)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop1
                     ENDIF

                  ELSE IF (NTOKEN == 3) THEN               ! TOKEN's connected by 'THRU'
                     IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
                        OGROUT(J) = IBSET(OGROUT(J),I)
                        READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AGRID_LO
                        IF (IOCHK1 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop1
                        ENDIF
                        READ(TOKEN(3),'(I8)',IOSTAT=IOCHK2) AGRID_HI
                        IF (IOCHK2 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop1
                        ENDIF
                                                           ! Error: "I1 THRU I2" has I1 > I2 (AGRID_LO > AGRID_HI)
                        IF ((IOCHK1 == 0) .AND. (IOCHK2 == 0)) THEN
                           IF (AGRID_LO > AGRID_HI) THEN
                              WRITE(ERR,1410) SETID
                              WRITE(ERR,1411) TOKSTR
                              WRITE(ERR,1423) AGRID_LO,AGRID_HI
                              WRITE(F06,1410) SETID
                              WRITE(F06,1411) TOKSTR
                              WRITE(F06,1423) AGRID_LO,AGRID_HI
                              WRITE(F06,*)
                              FATAL_ERR = FATAL_ERR + 1
                              JERR = JERR + 1
                              EXIT token_loop1
                           ENDIF

                           INUM = 0
                           DO M = AGRID_LO,AGRID_HI
                              CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, M, GRID_ID_ROW_NUM )
                              IF (GRID_ID_ROW_NUM /= -1) THEN
                                 INUM           = INUM + 1
                                 IDUM           = GROUT(GRID_ID_ROW_NUM,J)
                                 GROUT(GRID_ID_ROW_NUM,J) = IBSET(IDUM,I)
                              ENDIF
                           ENDDO
                           IF (INUM == 0) THEN
                              WARN_ERR = WARN_ERR + 1
                              WRITE(ERR,1403) WARN_NAME,AGRID_LO,AGRID_HI,SETID
                              IF (SUPWARN == 'N') THEN
                                 WRITE(F06,1403) WARN_NAME,AGRID_LO,AGRID_HI,SETID
                              ENDIF
                           ENDIF

                        ENDIF

                     ELSE
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop1
                     ENDIF

                  ELSE                                     ! Error: wrong tokens from set definition for some other reason

                     WRITE(ERR,1426) SUBR_NAME,NTOKEN
                     WRITE(F06,1426) SUBR_NAME,NTOKEN
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )                ! Coding error (NTOKEN /= 1 or 3), so quit

                  ENDIF

                  IF (ISTART <= TOKLEN) THEN
                     CYCLE token_loop1
                  ELSE
                     EXIT token_loop1
                  ENDIF

               ENDDO token_loop1

            ENDIF

         ENDDO j_loop1

      ENDDO

! **********************************************************************************************************************************
! 2nd: Process element output requests to generate OELOUT(i),ELOUT(j,i)

! OELOUT(i) is an overall indicator for each S/C for whether any of the types of elem output requests were made in
! Case Control.

! ELOUT(j,i) is a specific indicator for each S/C for every elem of the types of elem output requests that were made in
! Case Control

! Both OELOUT and ELOUT are set here to be used in LINK9.

! Bit position 0 in OELOUT(i), ELOUT(i) is for ELFORCE(NODE) requests (turned on if SC_ELFN(i) > 0 set in subr CC_ELFO).
! Bit position 1 in OELOUT(i), ELOUT(i) is for ELFORCE(ENGR) requests (turned on if SC_ELFE(i) > 0 set in subr CC_ELFO).
! Bit position 2 in OELOUT(i), ELOUT(i) is for STRESS        requests (turned on if SC_STRE(i) > 0 set in subr CC_STRE).
! Bit position 3 in OELOUT(i), ELOUT(i) is for STRAIN        requests (turned on if SC_STRN(i) > 0 set in subr CC_STRN).

! SC_ELFN(i), etc. are integer arrays containing SET ID's (-1 for 'ALL', 0 for 'NONE', positive integer no. for
! SETID) of the 3 possible elem related outputs requested in Case Control

      WARN_NAME = 'ELEMENT '

! I ranges over the number of types of elem output requests that can be made in Case Control (described above)

      DO I=0,MELOUTS-1

         IF (I == 0) ELM_OUTPUT_REQ = 'element nodal forces'
         IF (I == 1) ELM_OUTPUT_REQ = 'element engr forces'
         IF (I == 2) ELM_OUTPUT_REQ = 'element stresses'
         IF (I == 3) ELM_OUTPUT_REQ = 'element strains'

j_loop2: DO J = 1,NSUB

            IF      (I == 0) THEN
               SETID = SC_ELFN(J)
            ELSE IF (I == 1) THEN
               SETID = SC_ELFE(J)
            ELSE IF (I == 2) THEN
               SETID = SC_STRE(J)
            ELSE IF (I == 3) THEN
               SETID = SC_STRN(J)
            ENDIF

            IF (SETID ==  0) THEN                          ! Check for 'NONE'
               CYCLE j_loop2
            ENDIF

            IF (SETID == -1) THEN                          ! Check for 'ALL' (SETID = -1)
               OELOUT(J) = IBSET(OELOUT(J),I)              ! Set bit position I in OELOUT to 1
               DO K = 1,NELE
                  IF (ETYPE(K) /= 'PLOTEL  ') THEN
                     ELOUT(K,J) = IBSET(ELOUT(K,J),I)
                  ELSE
                     WARN_ERR = WARN_ERR + 1
                     WRITE(ERR,1406) ELM_OUTPUT_REQ, ETYPE(K)
                     IF (SUPWARN == 'N') THEN
                        WRITE(F06,1406) ELM_OUTPUT_REQ, ETYPE(K)
                     ENDIF
                  ENDIF
               ENDDO
               CYCLE j_loop2
            ENDIF

            IF (SETID > 0) THEN                            ! Check for actual set ID

               CALL SETPRO ( SETID, NULSET, TOKLEN )       ! Call SETPRO to get the data which was input for set SETID.
               DO K=1,LSETLN
                  TOKSTR(K:K) = ONE_SET_ARRAY(K)
               ENDDO

               IF (NULSET == 0) THEN                       ! Error was already written in SETPRO
                  JERR = JERR + 1
                  CYCLE j_loop2
               ENDIF

               ISTART = 1                                  ! Init ISTART, THRU, EXCEPT for the DO loop looking for tokens in TOKSTR
               THRU   = 'OFF'
               EXCEPT = 'OFF'

token_loop2:   DO                                          ! Call STOKEN in a DO loop until we run out of tokens in TOKSTR

                  CALL STOKEN ( SUBR_NAME, TOKSTR, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

                  IF (IERROR == 1) THEN                    ! Error: too long a token in the set definition
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1421) ERRTOK
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1421) ERRTOK
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop2

                  ELSE IF (IERROR == 2) THEN               ! Error: "I1 THRU I2" is missing I2 in the set definition
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1425) NTOKEN
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1425) NTOKEN
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop2

                  ELSE IF (IERROR == 3) THEN               ! Error: found "EXCEPT" but EXCEPT has already been turned "ON "
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1428)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1428)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop2

                  ELSE IF (IERROR == 4) THEN               ! Error: found "EXCEPT" but "THRU" has not been found yet
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1429)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1429)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop2

                  ELSE IF (NTOKEN == 1) THEN               ! Is the TOKEN a number not connected to a 'THRU'?
                     IF (TOKTYP(1) == 'INTEGER ') THEN
                        OELOUT(J) = IBSET(OELOUT(J),I)
                        READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AELEM
                        IF (IOCHK1 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop2
                        ELSE
                           CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, AELEM, ESORT1_ROW_NUM )
                           IF (ESORT1_ROW_NUM == -1) THEN
                              WARN_ERR = WARN_ERR + 1
                              WRITE(ERR,1402) WARN_NAME,AELEM,SETID
                              IF (SUPWARN == 'N') THEN
                                 WRITE(F06,1402) WARN_NAME,AELEM,SETID
                              ENDIF
                           ELSE
                              IF (EXCEPT == 'OFF') THEN
                                 IDUM        = ELOUT(ESORT1_ROW_NUM,J)
                                 IF (ETYPE(ESORT1_ROW_NUM) /= 'PLOTEL  ') THEN
                                    ELOUT(ESORT1_ROW_NUM,J) = IBSET(IDUM,I)
                                 ELSE
                                    WARN_ERR = WARN_ERR + 1
                                    WRITE(ERR,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                    IF (SUPWARN == 'N') THEN
                                       WRITE(F06,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                    ENDIF
                                 ENDIF
                              ENDIF
                              IF (EXCEPT == 'ON ') THEN
                                 IF ((AELEM >= AELEM_LO) .AND. (AELEM <= AELEM_HI)) THEN
                                    IDUM        = ELOUT(ESORT1_ROW_NUM,J)
                                    IF (ETYPE(ESORT1_ROW_NUM) /= 'PLOTEL  ') THEN
                                       ELOUT(ESORT1_ROW_NUM,J) = IBCLR(IDUM,I)
                                    ELSE
                                       WARN_ERR = WARN_ERR + 1
                                       WRITE(ERR,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                       IF (SUPWARN == 'N') THEN
                                          WRITE(F06,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                       ENDIF
                                    ENDIF
                                 ELSE
                                    IDUM        = ELOUT(ESORT1_ROW_NUM,J)
                                    IF (ETYPE(ESORT1_ROW_NUM) /= 'PLOTEL  ') THEN
                                       ELOUT(ESORT1_ROW_NUM,J) = IBSET(IDUM,I)
                                    ELSE
                                       WARN_ERR = WARN_ERR + 1
                                       WRITE(ERR,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                       IF (SUPWARN == 'N') THEN
                                          WRITE(F06,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                       ENDIF
                                    ENDIF
                                    EXCEPT = 'OFF'
                                    THRU   = 'OFF'
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF

                     ELSE IF (TOKTYP(1) == 'EXCEPT  ') THEN
                        IF (THRU /= 'ON ') THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1427)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1427)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop2
                        ENDIF

                     ELSE
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1422) TOKEN(1),TOKTYP(1)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1422) TOKEN(1),TOKTYP(1)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop2
                     ENDIF

                  ELSE IF (NTOKEN == 3) THEN               ! TOKEN's connected by 'THRU'
                     IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
                        OELOUT(J) = IBSET(OELOUT(J),I)
                        READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AELEM_LO
                        IF (IOCHK1 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop2
                        ENDIF
                        READ(TOKEN(3),'(I8)',IOSTAT=IOCHK2) AELEM_HI
                        IF (IOCHK2 > 0) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1404)
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1404)
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop2
                        ENDIF

                        IF ((IOCHK1 == 0) .AND. (IOCHK2 == 0)) THEN  ! Error: "I1 THRU I2" has I1 > I2 (AELEM_LO > AELEM_HI)
                           IF (AELEM_LO > AELEM_HI) THEN
                              WRITE(ERR,1410) SETID
                              WRITE(ERR,1411) TOKSTR
                              WRITE(ERR,1423) AELEM_LO,AELEM_HI
                              WRITE(F06,1410) SETID
                              WRITE(F06,1411) TOKSTR
                              WRITE(F06,1423) AELEM_LO,AELEM_HI
                              WRITE(F06,*)
                              FATAL_ERR = FATAL_ERR + 1
                              JERR = JERR + 1
                              EXIT token_loop2
                           ENDIF

                           INUM = 0
                           DO M = AELEM_LO,AELEM_HI
                              CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, M, ESORT1_ROW_NUM )
                              IF (ESORT1_ROW_NUM /= -1) THEN
                                 INUM = INUM + 1
                                 IDUM        = ELOUT(ESORT1_ROW_NUM,J)
                                 IF (ETYPE(ESORT1_ROW_NUM) /= 'PLOTEL  ') THEN
                                    ELOUT(ESORT1_ROW_NUM,J) = IBSET(IDUM,I)
                                 ELSE
                                    WARN_ERR = WARN_ERR + 1
                                    WRITE(ERR,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                    IF (SUPWARN == 'N') THEN
                                       WRITE(F06,1406) ELM_OUTPUT_REQ, ETYPE(ESORT1_ROW_NUM)
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDDO
                           IF (INUM == 0) THEN
                              WARN_ERR = WARN_ERR + 1
                              WRITE(ERR,1403) WARN_NAME,AELEM_LO,AELEM_HI,SETID
                              IF (SUPWARN == 'N') THEN
                                 WRITE(F06,1403) WARN_NAME,AELEM_LO,AELEM_HI,SETID
                              ENDIF
                           ENDIF

                        ENDIF

                     ELSE
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop2
                     ENDIF

                  ELSE                                     ! Error: wrong tokens from set definition for some other reason

                     WRITE(ERR,1426) SUBR_NAME,NTOKEN
                     WRITE(F06,1426) SUBR_NAME,NTOKEN
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )                        ! Coding error (NTOKEN /= 1 or 3), so quit

                  ENDIF

                  IF (ISTART <= TOKLEN) THEN
                     CYCLE token_loop2
                  ELSE
                     EXIT token_loop2
                  ENDIF

               ENDDO token_loop2

            ENDIF

         ENDDO j_loop2

      ENDDO

! **********************************************************************************************************************************
! 3rd: Process element output requests to generate OELDT, ELDT(j). The outputs are not subcase dependent. UEL, PEL
! are printed or written to disk for all subcases, if any request is made. CCELDT(1-16) is an integer array containing
! SET ID's (-1 for 'ALL', 0 for 'NONE', positive integer no. for SETID) of the 16 possible outputs requested by the
! Case Control ELDATA command:

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

! CCELDT is used to set bits in arrays OELDT and ELDT(j):
! OELDT tells whether any output requests were made.
! ELDT(j) tells which elements have a request:

! Bit  0 in OELDT, ELDT is for print to BUGFIL of elem geometric data
! Bit  1 in OELDT, ELDT is for print to BUGFIL of elem property and material info
! Bit  2 in OELDT, ELDT is for print to BUGFIL of elem thermal and pressure matrices    : PTE, PPE
! Bit  3 in OELDT, ELDT is for print to BUGFIL of elem mass matrix                      : ME
! Bit  4 in OELDT, ELDT is for print to BUGFIL of elem stiffness matrix                 : KE
! Bit  5 in OELDT, ELDT is for print to BUGFIL of elem stress & strain recovery matrices: SEi, STEi, BEi
! Bit  6 in OELDT, ELDT is for print to BUGFIL of elem displacement and load matrices   : UEL, PEL (all subcases)
! Bit  7 in OELDT, ELDT is for print to BUGFIL of elem shape fcns and Jacobian matrices
! Bit  8 in OELDT, ELDT is for print to BUGFIL of elem strain-displacement matrices
! Bit  9 in OELDT, ELDT is for print to BUGFIL of elem checks on strain-displ matrices for RB motion & constant strain
! Bit 10 in OELDT, ELDT is for write to F21FIL unformatted file of element              : PTE, PPE
! Bit 11 in OELDT, ELDT is for write to F22FIL unformatted file of element              : ME
! Bit 12 in OELDT, ELDT is for write to F23FIL unformatted file of element              : KE
! Bit 13 in OELDT, ELDT is for write to F24FIL unformatted file of element              : SEi, STEi, BEi
! Bit 14 in OELDT, ELDT is for write to F25FIL unformatted file of element              : UEL, PEL (all subcases)
! Bit 15 is not used

! Bit position k is turned on if any of that output is requested by CCELDT(k)

      WARN_NAME = 'ELEMENT '

outer:DO I=0,MELDTS-1

         SETID = CCELDT(I)

         IF (SETID == 0) THEN                              ! Check for 'NONE'
            CYCLE outer
         ENDIF

         IF (SETID == -1) THEN                             ! Check for 'ALL' (SETID = -1)
            OELDT = IBSET(OELDT,I)                         ! Set bit position I in OELDT to 1
            DO K = 1,NELE
               ELDT(K) = IBSET(ELDT(K),I)
            ENDDO
            CYCLE outer
         ENDIF

         IF (SETID > 0) THEN                               ! Check for actual set ID

            CALL SETPRO ( SETID, NULSET, TOKLEN )          ! Call SETPRO to get the data which was input for set SETID.
            DO K=1,LSETLN
               TOKSTR(K:K) = ONE_SET_ARRAY(K)
            ENDDO

            IF (NULSET == 0) THEN                          ! Error was already written in SETPRO
               JERR = JERR + 1
               CYCLE outer
            ENDIF

            ISTART = 1                                     ! Init ISTART, THRU, EXCEPT for the DO loop looking for tokens in TOKSTR
            THRU   = 'OFF'
            EXCEPT = 'OFF'

token_loop3:DO                                             ! Call STOKEN in a DO loop until we run out of tokens in TOKSTR

               CALL STOKEN ( SUBR_NAME, TOKSTR, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

               IF (IERROR == 1) THEN                       ! Error: too long a token in the set definition
                  WRITE(ERR,1410) SETID
                  WRITE(ERR,1411) TOKSTR
                  WRITE(ERR,1421) ERRTOK
                  WRITE(F06,1410) SETID
                  WRITE(F06,1411) TOKSTR
                  WRITE(F06,1421) ERRTOK
                  WRITE(F06,*)
                  FATAL_ERR = FATAL_ERR + 1
                  JERR = JERR + 1
                  EXIT token_loop3

               ELSE IF (IERROR == 2) THEN                  ! Error: "I1 THRU I2" is missing I2 in the set definition
                  WRITE(ERR,1410) SETID
                  WRITE(ERR,1411) TOKSTR
                  WRITE(ERR,1425) NTOKEN
                  WRITE(F06,1410) SETID
                  WRITE(F06,1411) TOKSTR
                  WRITE(F06,1425) NTOKEN
                  WRITE(F06,*)
                  FATAL_ERR = FATAL_ERR + 1
                  JERR = JERR + 1
                  EXIT token_loop3

               ELSE IF (IERROR == 3) THEN                  ! Error: found "EXCEPT" but EXCEPT has already been turned "ON "
                  WRITE(ERR,1410) SETID
                  WRITE(ERR,1411) TOKSTR
                  WRITE(ERR,1428)
                  WRITE(F06,1410) SETID
                  WRITE(F06,1411) TOKSTR
                  WRITE(F06,1428)
                  WRITE(F06,*)
                  FATAL_ERR = FATAL_ERR + 1
                  JERR = JERR + 1
                  EXIT token_loop3

               ELSE IF (IERROR == 4) THEN                  ! Error: found "EXCEPT" but "THRU" has not been found yet
                  WRITE(ERR,1410) SETID
                  WRITE(ERR,1411) TOKSTR
                  WRITE(ERR,1429)
                  WRITE(F06,1410) SETID
                  WRITE(F06,1411) TOKSTR
                  WRITE(F06,1429)
                  WRITE(F06,*)
                  FATAL_ERR = FATAL_ERR + 1
                  JERR = JERR + 1
                  EXIT token_loop3

               ELSE IF (NTOKEN == 1) THEN                  ! Is the TOKEN a number not connected to a 'THRU'?
                  IF (TOKTYP(1) == 'INTEGER ') THEN
                     OELDT = IBSET(OELDT,I)
                     READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AELEM
                     IF (IOCHK1 > 0) THEN
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1404)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1404)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop3
                     ELSE
                        CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, AELEM, ESORT1_ROW_NUM )
                        IF (ESORT1_ROW_NUM == -1) THEN
                           WARN_ERR = WARN_ERR + 1
                           WRITE(ERR,1402) WARN_NAME,AELEM,SETID
                           IF (SUPWARN == 'N') THEN
                              WRITE(F06,1402) WARN_NAME,AELEM,SETID
                           ENDIF
                        ELSE
                           IF (EXCEPT == 'OFF') THEN
                              IDUM     = ELDT(ESORT1_ROW_NUM)
                              ELDT(ESORT1_ROW_NUM) = IBSET(IDUM,I)
                           ENDIF
                           IF (EXCEPT == 'ON ') THEN
                              IF ((AELEM >= AELEM_LO) .AND. (AELEM <= AELEM_HI)) THEN
                                 IDUM     = ELDT(ESORT1_ROW_NUM)
                                 ELDT(ESORT1_ROW_NUM) = IBCLR(IDUM,I)
                              ELSE
                                 IDUM     = ELDT(ESORT1_ROW_NUM)
                                 ELDT(ESORT1_ROW_NUM) = IBSET(IDUM,I)
                                 EXCEPT = 'OFF'
                                 THRU   = 'OFF'
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF

                  ELSE IF (TOKTYP(1) == 'EXCEPT  ') THEN
                     IF (THRU /= 'ON ') THEN
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1427)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1427)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop3
                     ENDIF

                  ELSE
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1422) TOKEN(1),TOKTYP(1)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1422) TOKEN(1),TOKTYP(1)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop3
                  ENDIF

               ELSE IF (NTOKEN == 3) THEN                  ! TOKEN's connected by 'THRU'
                  IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
                     OELDT = IBSET(OELDT,I)
                     READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) AELEM_LO
                     IF (IOCHK1 > 0) THEN
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1404)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1404)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop3
                     ENDIF
                     READ(TOKEN(3),'(I8)',IOSTAT=IOCHK2) AELEM_HI
                     IF (IOCHK2 > 0) THEN
                        WRITE(ERR,1410) SETID
                        WRITE(ERR,1411) TOKSTR
                        WRITE(ERR,1404)
                        WRITE(F06,1410) SETID
                        WRITE(F06,1411) TOKSTR
                        WRITE(F06,1404)
                        WRITE(F06,*)
                        FATAL_ERR = FATAL_ERR + 1
                        JERR = JERR + 1
                        EXIT token_loop3
                     ENDIF

                     IF ((IOCHK1 == 0) .AND. (IOCHK2 == 0)) THEN ! Error: "I1 THRU I2" has I1 > I2 (AELEM_LO > AELEM_HI)
                        IF (AELEM_LO > AELEM_HI) THEN
                           WRITE(ERR,1410) SETID
                           WRITE(ERR,1411) TOKSTR
                           WRITE(ERR,1423) AELEM_LO,AELEM_HI
                           WRITE(F06,1410) SETID
                           WRITE(F06,1411) TOKSTR
                           WRITE(F06,1423) AELEM_LO,AELEM_HI
                           WRITE(F06,*)
                           FATAL_ERR = FATAL_ERR + 1
                           JERR = JERR + 1
                           EXIT token_loop3
                        ENDIF

                        INUM = 0
                        DO M = AELEM_LO,AELEM_HI
                           CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, M, ESORT1_ROW_NUM )
                           IF (ESORT1_ROW_NUM /= -1) THEN
                              INUM = INUM + 1
                              ELDT(ESORT1_ROW_NUM) = IBSET(ELDT(ESORT1_ROW_NUM),I)
                           ENDIF
                        ENDDO
                        IF (INUM == 0) THEN
                           WARN_ERR = WARN_ERR + 1
                           WRITE(ERR,1403) WARN_NAME,AELEM_LO,AELEM_HI,SETID
                           IF (SUPWARN == 'N') THEN
                              WRITE(F06,1403) WARN_NAME,AELEM_LO,AELEM_HI,SETID
                           ENDIF
                        ENDIF

                     ENDIF

                  ELSE
                     WRITE(ERR,1410) SETID
                     WRITE(ERR,1411) TOKSTR
                     WRITE(ERR,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                     WRITE(F06,1410) SETID
                     WRITE(F06,1411) TOKSTR
                     WRITE(F06,1424) TOKEN(1),TOKEN(2),TOKEN(3)
                     WRITE(F06,*)
                     FATAL_ERR = FATAL_ERR + 1
                     JERR = JERR + 1
                     EXIT token_loop3
                  ENDIF

               ELSE                                        ! Error: wrong tokens from set definition for some other reason

                  WRITE(ERR,1426) SUBR_NAME,NTOKEN
                  WRITE(F06,1426) SUBR_NAME,NTOKEN
                  WRITE(F06,*)
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )                   ! Coding error (NTOKEN /= 1 or 3), so quit

               ENDIF

               IF (ISTART <= TOKLEN) THEN
                  CYCLE token_loop3
               ELSE
                  EXIT token_loop3
               ENDIF

            ENDDO token_loop3

         ENDIF

      ENDDO outer

! **********************************************************************************************************************************
! Parameter PRTSCP = 1 requests output of subcase parameters.

      IF (PRTSCP == 1) CALL PARAM_PRTSCP_OUTPUT

! **********************************************************************************************************************************
! Check if any grid output was requested for SPOINT's that do not support these requests. Set WARN_ERR, write message and reset bit.

      NUM_WARN = 0
      DO I=1,NSUB
         DO J=1,NGRID
            IF (GRID(J,6) == 1) THEN
               IF (IAND(GROUT(J,I),IBIT(GROUT_GPFO_BIT)) > 0) THEN
                  NUM_WARN   = NUM_WARN + 1
                  GROUT(J,I) = IBCLR(GROUT(J,I),IBIT(GROUT_GPFO_BIT))
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,*)
      ENDDO

      IF (NUM_WARN > 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9991) 'GRID POINT FORCE BALANCE'
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9991) 'GRID POINT FORCE BALANCE'
         ENDIF
      ENDIF

! Check if any elem output was requested for elements that do not support these requests. Set WARN_ERR, write message and reset bit.

      NUM_WARN = 0
      DO I=1,NSUB
         DO J=1,NELE
            IF (ETYPE(J) == 'USERIN  ') THEN
               DO K=1,MELOUTS
                  IF (IAND(ELOUT(J,I),IBIT(K-1)) > 0) THEN
                     ELOUT(J,I) = IBCLR(ELOUT(J,I),K-1)
                     NUM_WARN   = NUM_WARN + 1
                     ELM_TYP(NUM_WARN) = ETYPE(J)
                     ELM_BIT(NUM_WARN) = IBIT(K-1)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      IF (NUM_WARN > 0) THEN
         DO I=1,NUM_WARN
            WARN_ERR = WARN_ERR + 1
            IF      (ELM_BIT(I) == 0) THEN
               ELM_BIT_NAME = 'NODE FORCE'
            ELSE IF (ELM_BIT(I) == 1) THEN
               ELM_BIT_NAME = 'ENGR FORCE'
            ELSE IF (ELM_BIT(I) == 2) THEN
               ELM_BIT_NAME = 'STRESS    '
            ELSE IF (ELM_BIT(I) == 3) THEN
               ELM_BIT_NAME = 'STRAIN    '
            ENDIF
            WRITE(ERR,9992) ELM_BIT_NAME, ELM_TYP(I)
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9992) ELM_BIT_NAME, ELM_TYP(I)
            ENDIF
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Check for errors and quit if any

      IF (JERR > 0) THEN
         WRITE(ERR,1405)
         WRITE(F06,1405)
         CALL OUTA_HERE ( 'Y' )                            ! Errors reading SET's, so quit
      ENDIF

! **********************************************************************************************************************************
! Set ANY_xxxx_OUTPUT variables to tell whether any output of a particular kind was requested in Case Control

      ANY_ACCE_OUTPUT = 0
      DO J=1,NSUB
         ANY_ACCE_OUTPUT = ANY_ACCE_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_ACCE_BIT))
      ENDDO

      ANY_DISP_OUTPUT = 0
      DO J=1,NSUB
         ANY_DISP_OUTPUT = ANY_DISP_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_DISP_BIT))
      ENDDO

      ANY_GPFO_OUTPUT = 0
      DO J=1,NSUB
         ANY_GPFO_OUTPUT = ANY_GPFO_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_GPFO_BIT))
      ENDDO

      ANY_MPCF_OUTPUT = 0
      DO J=1,NSUB
         ANY_MPCF_OUTPUT = ANY_MPCF_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_MPCF_BIT))
      ENDDO

      ANY_OLOA_OUTPUT = 0
      DO J=1,NSUB
         ANY_OLOA_OUTPUT = ANY_OLOA_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_OLOA_BIT))
      ENDDO

      ANY_SPCF_OUTPUT = 0
      DO J=1,NSUB
         ANY_SPCF_OUTPUT = ANY_SPCF_OUTPUT + IAND(OGROUT(J),IBIT(GROUT_SPCF_BIT))
      ENDDO

      ANY_ELFE_OUTPUT = 0
      DO J=1,NSUB
         ANY_ELFE_OUTPUT = ANY_ELFE_OUTPUT + IAND(OELOUT(J),IBIT(ELOUT_ELFE_BIT))
      ENDDO

      ANY_ELFN_OUTPUT = 0
      DO J=1,NSUB
         ANY_ELFN_OUTPUT = ANY_ELFN_OUTPUT + IAND(OELOUT(J),IBIT(ELOUT_ELFN_BIT))
      ENDDO

      ANY_STRE_OUTPUT = 0
      DO J=1,NSUB
         ANY_STRE_OUTPUT = ANY_STRE_OUTPUT + IAND(OELOUT(J),IBIT(ELOUT_STRE_BIT))
      ENDDO

      ANY_STRN_OUTPUT = 0
      DO J=1,NSUB
         ANY_STRN_OUTPUT = ANY_STRN_OUTPUT + IAND(OELOUT(J),IBIT(ELOUT_STRN_BIT))
      ENDDO

! Now write out subcase data to L1D

      DATA_SET_NAME = 'S/C NUMBERS, TITLING, LOAD SET IDS'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NSUB
      DO I=1,NSUB
         WRITE(L1D) SCNUM(I)
         WRITE(L1D) TITLE(I)
         WRITE(L1D) STITLE(I)
         WRITE(L1D) LABEL(I)
         WRITE(L1D) SUBLOD(I,1)
         WRITE(L1D) SUBLOD(I,2)
      ENDDO

      DATA_SET_NAME = 'OGROUT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NSUB
      DO I = 1,NSUB
         WRITE(L1D) OGROUT(I)
      ENDDO

      DATA_SET_NAME = 'OELOUT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NSUB
      DO I = 1,NSUB
         WRITE(L1D) OELOUT(I)
      ENDDO

      DATA_SET_NAME = 'OELDT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) OELDT

      DATA_SET_NAME = 'GROUT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NGRID
      WRITE(L1D) NSUB
      DO I = 1,NGRID
         DO J = 1,NSUB
            WRITE(L1D) GROUT(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'ELOUT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NELE
      WRITE(L1D) NSUB
      DO I = 1,NELE
         DO J = 1,NSUB
            WRITE(L1D) ELOUT(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'ELDT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) NELE
      DO I=1,NELE
         WRITE(L1D) ELDT(I)
      ENDDO

      DATA_SET_NAME = 'ANY_xxxx_OUTPUT'
      WRITE(L1D) DATA_SET_NAME
      WRITE(L1D) ANY_ACCE_OUTPUT
      WRITE(L1D) ANY_DISP_OUTPUT
      WRITE(L1D) ANY_OLOA_OUTPUT
      WRITE(L1D) ANY_SPCF_OUTPUT
      WRITE(L1D) ANY_MPCF_OUTPUT
      WRITE(L1D) ANY_GPFO_OUTPUT
      WRITE(L1D) ANY_ELFN_OUTPUT
      WRITE(L1D) ANY_ELFE_OUTPUT
      WRITE(L1D) ANY_STRE_OUTPUT
      WRITE(L1D) ANY_STRN_OUTPUT



      RETURN

! **********************************************************************************************************************************
 1402 FORMAT(' *WARNING    : ',A8,1X,I8,' NOT FOUND ON SET ENTRY, ID = ',I8)

 1403 FORMAT(' *WARNING    : NO ',A8,' IN RANGE ',I8,' THRU ',I8,' REQUESTED ON SET ENTRY, ID = ',I8)

 1404 FORMAT('               ERROR ATTEMPTING TO READ INTEGER ON SET ENTRY')

 1405 FORMAT(' PROCESSING TERMINATED IN SUBROUTINE SCPRO DUE TO PRIOR LISTED ERRORS')

 1406 FORMAT(' WARNING     : REQUEST FOR ELEMENT OUTPUT FOR "', A,'", FOR ELEMENT TYPE ',A,' NOT ALLOWED. REQUEST IGNORED')

 1410 FORMAT(' *ERROR  1410: SYNTAX ERROR ON THE CASE CONTROL SET ENTRY ID = ',I8,' LISTED BELOW:')

 1411 FORMAT(15X,A)

 1421 FORMAT(14X,' TOKENS ON SET ENTRIES MUST BE LESS THAN 8 CHARACTERS.'  &
          ,/,14X,' LISTED BELOW IS ONE FOUND THAT EXCEEDS THIS REQUIREMENT:',/,14X,A)

 1422 FORMAT(14X,' EXPECTING TO READ INTEGER TOKEN ON SET ENTRY BUT TOKEN IS: ',A8,' OF TYPE ',A8)

 1423 FORMAT(14X,' WHEN "I1 THRU I2" IS USED ON SET ENTRY, I1 MUST BE GREATER THEN I2. INPUT WAS I1 = ',I8,' AND I2 = ',I8)

 1424 FORMAT(14X,' EXPECTING TO READ "I1 THRU I2" ON SET ENTRY BUT READ: ',3A8)

 1425 FORMAT(14X,' EXPECTING TO READ 3 TOKENS: "I1 THRU I2", BUT ONLY ',I8,' TOKENS WERE FOUND')

 1426 FORMAT(' *ERROR  1426: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NOTKEN SHOULD BE 1 OR 3 BUT IS = ',I8)

 1427 FORMAT(14X,' THE QUALIFIER "EXCEPT" CAN ONLY MODIFY A PREVIOUS "THRU"')

 1428 FORMAT(14X,' READ MORE THAN ONE "EXCEPT" BEFORE NEXT "THRU"')

 1429 FORMAT(14X,' READ "EXCEPT" BEFORE AN "I1 THRU I2" WAS SPECIFIED')

 9991 FORMAT(' *WARNING    : ',A,' NOT ALLOWED FOR SPOINT''s')

 9992 FORMAT(' *WARNING    : ELEMENT ',A,' OUTPUTS REQUESTED IN CASE CONTROL FOR ELEMENT TYPE ',A,' NOT ALLOWED')

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE SETPRO ( SETID_IN, NULSET, SLEN )

! Processes data from Case Control SET cards that have all been loaded into array ALL_SETS_ARRAY (when read by
! subroutine CC_SET) to get one sets' data (which is put into array ONE_SET_ARRAY).

! For example, if the Case Control SET input is:

! SET  1243= 12, 34   ,1358, 8976,
!    201,302, 406
! SET 6    = 1001,   3002

! Then variable ALL_SETS_ARRAY, as output from subroutine LOADC, is:

! SET  1243= 12, 34   ,1358, 8976,201,302, 406SET 6    = 1001,   3002

! There is no limit on the number of SET's in Case Control


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LSETLN, MAX_TOKEN_LEN, SETLEN, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ALL_SETS_ARRAY, ONE_SET_ARRAY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SETPRO'
      CHARACTER(LSETLN*BYTE)          :: SETCHR            ! Equiv to ALL_SETS_ARRAY (used in finding strings in ONE_SET_ARRAY)
      CHARACTER(8*BYTE)               :: TOKENI            ! String that should contain an integer set ID from SETCHR.
      CHARACTER(LSETLN*BYTE)          :: TOKSTR            ! Equiv to ONE_SET_ARRAY (used in setting strings in ONE_SET_ARRAY)
      CHARACTER(8*BYTE)               :: TOKTY1            ! 8 char string from subr TOKSTR (what type of data was in TOKENI)

      INTEGER(LONG), INTENT(IN)       :: SETID_IN          ! The set ID whose data we are looking for
      INTEGER(LONG), INTENT(OUT)      :: NULSET            ! 0 if no set with ID equal to SETID_IN was found, 1 if set was found
      INTEGER(LONG), INTENT(OUT)      :: SLEN              ! Length (in bytes) of the SET description put into ONE_SET_ARRAY
      INTEGER(LONG)                   :: DATA_BEG          ! Position in SETCHR where a sets data begins
      INTEGER(LONG)                   :: DATA_END          ! Position in SETCHR where a sets data ends
      INTEGER(LONG)                   :: ECOL              ! Position in SETCHR where an '=' is located
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: I                 ! DO loop indices
      INTEGER(LONG)                   :: SETID             ! A set ID in a SET in SETCHR
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a file
      INTEGER(LONG)                   :: POSN              ! Pos'n in SETCHR where some character string begins
      INTEGER(LONG)                   :: SET_1_BEG         ! Pos'n in SETCHR where one  'SET' begins (where the chars 'SET' begin)
      INTEGER(LONG)                   :: SET_2_BEG         ! Pos'n in SETCHR where next 'SET' begins (where the chars 'SET' begin)
      INTEGER(LONG)                   :: SID_BEG           ! Pos'n in SETCHR where a sets ID begins
      INTEGER(LONG)                   :: SID_END           ! Pos'n in SETCHR where a sets ID ends


      INTRINSIC INDEX



! **********************************************************************************************************************************
! Initialize outputs

      NULSET = 0
      SLEN   = 0

      DO I=1,LSETLN                                        ! Load SETCHR with data in ALL_SETS_ARRAY and init ONE_SET_ARRAY
         SETCHR(I:I) = ALL_SETS_ARRAY(I)
         ONE_SET_ARRAY(I) = ' '
      ENDDO

      SET_1_BEG = INDEX(SETCHR(1:),'SET')                  ! Are any SET's in SETCHR ?  If not, print error and return
      IF (SET_1_BEG == 0) THEN
         NULSET = 0
         WRITE(ERR,1405) SETID_IN
         WRITE(F06,1405) SETID_IN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! Begin loop on extracting SET's from SETCHR (ALL_SETS_ARRAY)

outer:DO

         TOKSTR(1:) = ' '

         POSN = INDEX(SETCHR(SET_1_BEG+1:),'SET')          ! Find next occurrance of 'SET' so we can find SET definition in between
         IF (POSN == 0) THEN
            SET_2_BEG = SETLEN + 1                         ! There is no other SET, so make "starting" posn > SETLEN
         ELSE
            SET_2_BEG = SET_1_BEG + POSN                   ! This is where the next SET begins
         ENDIF

         POSN = INDEX(SETCHR(SET_1_BEG:),'=')              ! Find the col where "=" is for current SET. The set ID (SETID) should
         ECOL = SET_1_BEG + POSN - 1                       ! be between SET and "=". The "=" is in col ECOL.

         DATA_BEG = ECOL  + 1                              ! Find col where the SET data is (following "=" and skipping delimiters).
         DATA_END = SET_2_BEG - 1

i_loop1: DO I=DATA_BEG,DATA_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) THEN
               DATA_BEG = DATA_BEG + 1
            ELSE
               EXIT i_loop1
            ENDIF
         ENDDO i_loop1

         TOKENI(1:) = ' '                                  ! Now find the set ID and load it into TOKENI
         SID_BEG    = SET_1_BEG + 3
         SID_END    = ECOL - 1
         SETID      = 0
         K          = 0

i_loop2: DO I=SID_BEG,SID_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) THEN
               CYCLE i_loop2                               ! We haven't gotten to the set ID yet, so cycle
            ENDIF
            K = K + 1
            IF (K > MAX_TOKEN_LEN) THEN                    ! Token too long. Set NULSET, SET_2_BEG and cycle outer
               NULSET = 0                                  ! so we can continue to look for
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
            TOKENI(K:K) = SETCHR(I:I)                      ! Load set ID chars into TOKENI
         ENDDO i_loop2

         CALL TOKCHK (TOKENI, TOKTY1 )                     ! TOKCHK will determine type of token in TOKENI
         IF (TOKTY1 == 'INTEGER ') THEN                    ! If token is an integer, read it into SETID
            READ(TOKENI,'(I8)',IOSTAT=IOCHK) SETID
            IF (IOCHK > 0) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1407)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1407)
               ENDIF
               NULSET = 0
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ELSE                                              ! Token wasn't an integer, so set NULSET, SET_1_BEG and cycle outer
            NULSET = 0
            SET_1_BEG = SET_2_BEG
            CYCLE outer
         ENDIF

         IF (SETID == SETID_IN) THEN                       ! Now get the data for set ID = SETID
            NULSET = 1
            SLEN = DATA_END - DATA_BEG + 1
            TOKSTR(1:) = SETCHR(DATA_BEG:DATA_END)
            EXIT                                           ! We have our set data for set ID = SETID_IN, so exit outer
         ELSE
            NULSET = 0
            SLEN   = 0
            IF (SET_2_BEG > SETLEN) THEN
               IF (NULSET == 0) THEN
                  WRITE(ERR,1405) SETID_IN
                  WRITE(F06,1405) SETID_IN
                  FATAL_ERR = FATAL_ERR + 1
                  RETURN
               ENDIF
            ELSE
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ENDIF

      ENDDO outer

      IF (NULSET /= 0) THEN
         DO I=1,LSETLN                                     ! Load TOKSTR data into ONE_SET_ARRAY
            ONE_SET_ARRAY(I) = TOKSTR(I:I)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1405 FORMAT(' *ERROR  1405: SET ID ',I8,' NOT FOUND')

 1407 FORMAT(' *WARNING    : ERROR READING SET ID ON CASE CONTROL SET ENTRY. ENTRY IGNORED')

! **********************************************************************************************************************************

      END SUBROUTINE SETPRO

! ##################################################################################################################################

      SUBROUTINE PARAM_PRTSCP_OUTPUT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      IMPLICIT NONE

      CHARACTER(1*BYTE)               :: NOT_ALLOWED_MSG   ! Message indicator

      INTEGER(LONG)                   :: II, JJ, KK        ! Local DO loop indices


! **********************************************************************************************************************************
      WRITE(F06,1001)

! Print grid DISP/FORCE/GPFB output requests

      PRNTOUT = 'N'
      NOT_ALLOWED_MSG = 'N'
      DO II=1,NSUB
         IF (OGROUT(II) > 0) THEN
            PRNTOUT = 'Y'
         ENDIF
      ENDDO
      IF (PRNTOUT == 'Y') THEN
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1013)
         DO II=1,NSUB
            DO JJ=1,NGRID
               PRNTOUT = 'N'
               DO KK=1,MGROUTS
                  IF (IAND(GROUT(JJ,II),IBIT(KK-1)) > 0) THEN
                     PRNTOUT = 'Y'
                     IF ((GRID(JJ,6) == 1) .AND. (KK-1 == GROUT_GPFO_BIT)) THEN
                        CBIT(KK) = '*'
                        NOT_ALLOWED_MSG = 'Y'
                     ELSE
                        CBIT(KK) = 'Y'
                     ENDIF
                  ELSE
                     CBIT(KK) = '-'
                  ENDIF
               ENDDO
               IF (PRNTOUT  == 'Y') THEN
                  IF (GRID(JJ,6) == 6) THEN
                     WRITE(F06,1014) GRID_ID(JJ),SCNUM(II),(CBIT(KK),KK=1,MGROUTS)
                  ELSE
                     WRITE(F06,1015) GRID_ID(JJ),SCNUM(II),(CBIT(KK),KK=1,MGROUTS)
                  ENDIF
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF (NOT_ALLOWED_MSG == 'Y') THEN
            WRITE(F06,*) '* indicates that output request was made but it is not allowed for scalar points'
         ENDIF
         WRITE(F06,1016)
      ENDIF

! Print element ELFi/STR output requests

      PRNTOUT = 'N'
      NOT_ALLOWED_MSG = 'N'
      DO II=1,NSUB
         IF (OELOUT(II) > 0) THEN
            PRNTOUT = 'Y'
         ENDIF
      ENDDO
      IF (PRNTOUT == 'Y') THEN
         WRITE(F06,1011)
         WRITE(F06,1021)
         WRITE(F06,1022)
         DO II=1,NSUB
            DO JJ=1,NELE
               PRNTOUT = 'N'
               DO KK=1,MELOUTS
                  IF (IAND(ELOUT(JJ,II),IBIT(KK-1)) > 0) THEN
                     PRNTOUT = 'Y'
                     IF (ETYPE(JJ) == 'USERIN  ') THEN
                        CBIT(KK) = '*'                     ! Output for this elem type requested but not allowed
                        NOT_ALLOWED_MSG = 'Y'
                     ELSE
                        CBIT(KK) = 'Y'                     ! Output requested and allowed
                     ENDIF
                  ELSE
                     CBIT(KK) = '-'                        ! Output not requested
                  ENDIF
               ENDDO
               IF (PRNTOUT == 'Y') THEN
                  WRITE(F06,1023) ESORT1(JJ),ETYPE(JJ),SCNUM(II),(CBIT(KK),KK=1,MELOUTS)
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF (NOT_ALLOWED_MSG == 'Y') THEN
            WRITE(F06,*) '* indicates that output request was made but it is not allowed for that element type'
         ENDIF
         WRITE(F06,1024)
      ENDIF

! Print ELDATA output requests

      IF (OELDT /= 0) THEN
         WRITE(F06,1011)
         WRITE(F06,1031)
         WRITE(F06,1032)
         PRNTOUT = 'Y'
      ENDIF
      IF (PRNTOUT == 'Y') THEN
         DO JJ=1,NELE
            PRNTOUT = 'N'
            DO KK=1,MELDTS
               IF (IAND(ELDT(JJ),IBIT(KK-1)) > 0) THEN
                  CBIT(KK) = 'Y'
                  PRNTOUT = 'Y'
               ELSE
                  CBIT(KK) = '-'
               ENDIF
            ENDDO
            IF (PRNTOUT == 'Y') THEN
               WRITE(F06,1033) ESORT1(JJ),ETYPE(JJ),(CBIT(KK),KK=1,MELDTS)
            ENDIF
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1034)
         WRITE(F06,*)
      ENDIF


      WRITE(F06,1035)

! **********************************************************************************************************************************
 1001 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::::START PARAM PRTSCP OUTPUT FROM SUBROUTINE SUBCASE_PROC::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1011 FORMAT(' ******************************************************************************************************************',&
              '*****************')

 1012 FORMAT(32X,'CASE CONTROL GRID POINT OUTPUT REQUESTS',/)

 1013 FORMAT(76X,'D  O  S  M  G  A',/,76X,'I  L  P  P  P  C',/,76X,'S  O  C  C  F  C' ,/,76X,'P  A  F  F  O  E',/)

 1014 FORMAT(10X,' Grid   point ',I8,' output requests for subcase ',I8,' are:',8(2X,A1))

 1015 FORMAT(10X,' Scalar point ',I8,' output requests for subcase ',I8,' are:',8(2X,A1))

 1016 FORMAT(1X                                                                                                                 ,/,&
             ' DISP indicates grid point displacement output requests'                                                          ,/,&
             ' OLOA indicates grid point applied load output requests'                                                          ,/,&
             ' SPCF indicates grid point forces of single point constraint output requests'                                     ,/,&
             ' MPCF indicates grid point forces of multi  point constraint output requests'                                     ,/,&
             ' GPFO indicates grid point force balance output requests'                                                         ,/,&
             ' ACCE indicates grid point acceleration output requests',/)

 1021 FORMAT(24X,'CASE CONTROL ELEMENT OUTPUT REQUESTS',/)


 1022 FORMAT(76X,'E  E  S  S  ',/,76X,'L  L  T  T  ',/,76X,'F  F  R  R  ',/,76X,'N  E  E  N',/)

 1023 FORMAT(' Element ',I8,',type ',A,', output requests for subcase ',I8,' are:',8(2X,A1))

 1024 FORMAT(1X                                                                                                                 ,/,&
             ' ELFN indicates element nodal load output requests'                                                               ,/,&
             ' ELFE indicates element engineering force output requests'                                                        ,/,&
             ' STRE indicates element stress output requests'                                                                   ,/,&
             ' STRN indicates element strain output requests',/)

 1031 FORMAT(82X,'CASE CONTROL ELDATA OUTPUT REQUESTS',/)

 1032 FORMAT(88X,'B I T   P O S I T I O N',/,77X,'0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15',/)

 1033 FORMAT(1X,'ELEMENT ',I8,',TYPE ',A,', ELDATA OUTPUT REQUESTS (ALL SUBCASES) ARE:',16(2X,A1))

 1034 FORMAT(1X                                                                                                                 ,/,&
             ' Bit position  0 is for printed output of element grid geometric data'                                            ,/,&
             ' Bit position  1 is for printed output of element property, material and other data'                              ,/,&
             ' Bit position  2 is for printed output of element nodal mass.'                                                    ,/,&
             ' Bit position  3 is for printed output of element thermal and pressure loads in local element coords.'            ,/,&
             ' Bit position  4 is for printed output of element stiffness matrix in local element coords.'                      ,/,&
             ' Bit position  5 is for printed output of element stress recovery matrices in local element coords.'              ,/,&
             ' Bit position  6 is for printed output of element displacements and nodal loads.'                                 ,/,&
             ' Bit position  7 is for printed output of elem shape fcns and Jacobian matrices from subrs SHPxyz, JAC2D etc.'    ,/,&
             ' Bit position  8 is for printed output of strain-displacement matrices in subroutines BMQMEM, etc.'               ,/,&
             ' Bit position  9 is for printed output of elem rigid body and constant strain checks (subroutine BCHECK).'        ,/,&
             ' Bit position 10 is for unformatted file output of element mass.'                                                 ,/,&
             ' Bit position 11 is for unformatted file output of element thermal and pressure loads.'                           ,/,&
             ' Bit position 12 is for unformatted file output of element stiffness matrix.'                                     ,/,&
             ' Bit position 13 is for unformatted file output of element stress recovery matrices.'                             ,/,&
             ' Bit position 14 is for unformatted file output of element displs and loads.'                                     ,/,&
             ' Bit position 15 is not used.')

 1035 FORMAT(' :::::::::::::::::::::::::::::::::::::::END PARAM PRTSCP OUTPUT FROM SUBROUTINE SUBCASE_PROC:::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE PARAM_PRTSCP_OUTPUT

      END SUBROUTINE SUBCASE_PROC

   END MODULE SUBCASE_PROCESSING
