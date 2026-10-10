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

   MODULE INPUT_STAGE_DISPATCH

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELEPRO, LOADB, LOADB0, LOADB_RESTART, LOADC, LOADC0, LOADE, LOADE0, READ_BDF_LINE

   CONTAINS

      SUBROUTINE ELEPRO ( INCR_NELE, JCARD, NFIELD, NMORE,                                                                         &
                          CHK_FLD2, CHK_FLD3, CHK_FLD4, CHK_FLD5, CHK_FLD6, CHK_FLD7, CHK_FLD8, CHK_FLD9 )

! Element connection data processor.

!  1) Increments count on number of elements and checks that the number does not exceed with the total count made by subr LOADB0
!  2) Checks to make sure that, when all integer data is put into EDAT for this element, that NEDAT won't exceed the max value that
!     was counted in LOADB0. NMORE
!  3) Verifies that all integer ID's are > 0 (since it is generally reading only elem ID, prop ID and grid connections)
!  4) Loads NFIELD fields of data for this element into integer array EDAT (the calling routine may add more data)

! If there are no errors:

!  3) Reads element connection data into EDAT
!  4) Resets pointer array, EPNT.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR, JF, LEDAT, LELE, NEDAT, NELE, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, EPNT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELEPRO'
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD2          ! If 'N', then if field 2 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD3          ! If 'N', then if field 3 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD4          ! If 'N', then if field 4 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD5          ! If 'N', then if field 5 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD6          ! If 'N', then if field 6 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD7          ! If 'N', then if field 7 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD8          ! If 'N', then if field 8 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: CHK_FLD9          ! If 'N', then if field 9 is blank it will not be checked for > 0
      CHARACTER(LEN=*), INTENT(IN)    :: INCR_NELE         ! If 'Y', increment NELE. Otherwise do not increment NELE
      CHARACTER(LEN=*), INTENT(IN)    :: JCARD(10)         ! The 10 fields of a Bulk Data card
      CHARACTER( 9*BYTE)              :: NAME = '         '! Name for output error purposes
      CHARACTER(LEN=LEN(CHK_FLD2))    :: CHK_FLD_ARRAY(2:9)! Array containing CHK_FLDi's

      INTEGER(LONG), INTENT(IN)       :: NFIELD            ! Number of card fields to read from JCARD (start w/ field 2)
      INTEGER(LONG), INTENT(IN)       :: NMORE             ! Number of terms that have to be written to EDAT for this element
!                                                            in total (not just here). There may be multiple calls to this subr for
!                                                            one element in which case NMORE is the additional amount of data to be
!                                                            written in this call. In addition, some data fields for this element
!                                                            may be written to EDAT in the subr that called this subr in which case
!                                                            NMORE may be greater than NFIELD so that a check on whether that
!                                                            additional data will fit into EDAT can be made here.
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: I,J               ! DO loop indices




! **********************************************************************************************************************************
      CHK_FLD_ARRAY(2) = CHK_FLD2
      CHK_FLD_ARRAY(3) = CHK_FLD3
      CHK_FLD_ARRAY(4) = CHK_FLD4
      CHK_FLD_ARRAY(5) = CHK_FLD5
      CHK_FLD_ARRAY(6) = CHK_FLD6
      CHK_FLD_ARRAY(7) = CHK_FLD7
      CHK_FLD_ARRAY(8) = CHK_FLD8
      CHK_FLD_ARRAY(9) = CHK_FLD9

      IF      (INCR_NELE == 'Y') THEN
         NELE = NELE + 1                                   ! Increment element counter, 'NELE'
      ELSE IF (INCR_NELE /= 'N') THEN
         NAME = 'INCR_NELE'
         WRITE(ERR,1013) NAME, INCR_NELE
         WRITE(F06,1013) NAME, INCR_NELE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO I=2,9
         IF ((CHK_FLD_ARRAY(I) /= 'Y') .AND. (CHK_FLD_ARRAY(I) /= 'N')) THEN
            NAME = 'CHK_FLD'
            WRITE(ERR,1013) NAME, CHK_FLD_ARRAY(I)
            WRITE(F06,1013) NAME, CHK_FLD_ARRAY(I)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDDO

      IF (NELE > LELE) THEN
         WRITE(ERR,1000) SUBR_NAME, LELE
         WRITE(F06,1000) SUBR_NAME, LELE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (too many elems), so quit
      ENDIF

      IF ((NEDAT + NMORE) > LEDAT) THEN
         WRITE(ERR,1001) SUBR_NAME, LEDAT
         WRITE(F06,1001) SUBR_NAME, LEDAT
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (too much EDAT data), so quit
      ENDIF

      IF (INCR_NELE == 'Y') THEN
         EPNT(NELE) = NEDAT+1                              ! Set element pointer EPNT to start for new element
      ENDIF

jdo:  DO J=1,NFIELD                                        ! Load element data into array EDAT

! Testing for (1:6) == 'CBUSH ' or 'CBUSH*' instead of just (1:5) == 'CBUSH' may help protect against
! matching keywords such as CBUSH1D and CBUSH2D that could be added in the future. It's also
! used in BD_CQUAD, BD_CTRIA, BD_PLOAD2, and LOADB_RESTART.

         IF ((J == 3) .AND. ((JCARD(1)(1:6) == 'CBUSH ') .OR. (JCARD(1)(1:6) == 'CBUSH*'))) THEN
            NEDAT = NEDAT + 1
            IF (JCARD(5)(1:) == ' ') THEN                  ! CBUSH has G2 blank so BUSH elem has only G1 grid. Other end is grounded
               EDAT(NEDAT) = 1
            ELSE                                           ! CBUSH G2 not blank so BUSH elem has 2 grids
               EDAT(NEDAT) = 2
            ENDIF
         ENDIF

         CALL I4FLD ( JCARD(J+1), JF(J+1), I4INP )
         IF (IERRFL(J+1) == 'N') THEN
            IF (CHK_FLD_ARRAY(J+1) == 'Y') THEN
               IF (I4INP <= 0) THEN                        ! No error if BUSH and G2 field is blank (BUSH grounded at end G2)
                  IF (((JCARD(1)(1:8) == 'CBUSH   ') .OR. (JCARD(1)(1:8) == 'CBUSH*  ')) .AND. (JCARD(5)(1:) == ' ')) THEN
                     EXIT jdo
                  ENDIF
                  WRITE(ERR,1021) JCARD(1), JCARD(2), I4INP, JF(J+1)
                  WRITE(F06,1021) JCARD(1), JCARD(2), I4INP, JF(J+1)
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDIF
            NEDAT = NEDAT + 1
            EDAT(NEDAT) = I4INP
         ENDIF

      ENDDO jdo



      RETURN

! **********************************************************************************************************************************
 1000 FORMAT(' *ERROR  1000: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ELEMENTS. LIMIT = ',I8)

 1001 FORMAT(' *ERROR  1001: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MUCH ELEMENT DATA IN ARRAY EDAT. LIMIT = ',I8)

 1013 FORMAT(' *ERROR  1013: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' SUBROUTINE CALLED WITH ILLEGAL VALUE FOR ',A,' = ',A,'. MUST BE EITHER Y OR N')

 1021 FORMAT(' *ERROR  1021: ',A,A,' HAS INTEGER = ',I8,' IN FIELD ',I3,'. MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE ELEPRO

      SUBROUTINE LOADB

      ! LOADB reads in the Bulk Data deck.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, IMB_BLANK, JF, LIND_GRDS_MPCS,               &
                                         LSUB, LLOADC, LMPCADDC, LSPCADDC, MDT, MTDAT_TEMPP1, MTDAT_TEMPRB,                        &
                                         MAX_GAUSS_POINTS, MAX_STRESS_POINTS,                                                      &
                                         MELGP, MELDOF, MMPC, MOFFSET, NBAROR, NBEAMOR, NFORCE,NGRAV, NGRDSET, NGRID, NLOAD, NMPC, &
                                         NMPCADD, NPCOMP, NRBAR, NRBE1, NRBE2, NRFORCE, NRSPLINE, NSLOAD, NSPOINT, NSPC, NSPC1,    &
                                         NSPCADD, NPBAR, NPBARL, NPLOAD, NSUB, NUM_MPCSIDS, NUM_PARTVEC_RECORDS, PROG_NAME,        &
                                         SOL_NAME, NCBAR, NCBEAM, NCBUSH, NCHEXA20, NCHEXA8, NCPENTA15, NCPENTA6, NCQUAD4,         &
                                         NCQUAD4K, NCQUAD8, NCROD, NCSHEAR, NCTETRA10, NCTETRA4, NCTRIA3, NCTRIA3K, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  GRIDSEQ, IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T, QUADAXIS, SUPINFO, SUPWARN
      USE OUTPUT4_MATRICES, ONLY      :  NUM_PARTN_REQUESTS
      USE MODEL_STUF, ONLY            :  FORMOM_SIDS, GRAV_SIDS, IOR3D_MAX, LOAD_SIDS,                                             &
                                         MPCSET, MPC_SIDS, MPCSIDS, MPCADD_SIDS, PBAR, RPCOMP, PRESS_SIDS, RFORCE_SIDS,            &
                                         RPBAR, SLOAD_SIDS, SPC_SIDS, SPC1_SIDS, SPCADD_SIDS, SPCSET, CC_EIGR_SID, SCNUM, SUBLOD


      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  FFIELD, FFIELD2
      USE DOF_SETS, ONLY              :  BD_ASET, BD_ASET1, BD_USET, BD_USET1, BD_SUPORT
      USE ROD_BAR_BEAM_CARDS, ONLY    :  BD_BAROR, BD_BEAMOR, BD_CBAR, BD_CROD, BD_CONROD, BD_PROD, BD_PBAR, BD_PBARL, BD_PBEAM, BD_PLOTEL
      USE SPRING_BUSH_MASS, ONLY      :  BD_CELAS1, BD_CELAS2, BD_CELAS3, BD_CELAS4, BD_PELAS, BD_CBUSH, BD_PBUSH, BD_CMASS1, BD_CMASS2, BD_CMASS3, BD_CMASS4, BD_PMASS, BD_CONM2
      USE SOLID_CARDS, ONLY           :  BD_CHEXA, BD_CPENTA, BD_CTETRA, BD_PSOLID
      USE GRID_COORDINATES, ONLY      :  BD_CORD, BD_GRID, BD_GRDSET, BD_SEQGP, BD_SPOINT, BD_SNORM
      USE SHELL_COMPOSITE_CARDS, ONLY :  BD_CQUAD, BD_CQUAD8, BD_CTRIA, BD_CSHEAR, BD_PSHEAR, BD_PSHEL, BD_PCOMP, BD_PCOMP1
      USE USER_ELEMENTS, ONLY         :  BD_CUSER1, BD_CUSERIN, BD_PUSER1, BD_PUSERIN
      USE DEBUG_CARDS, ONLY           :  BD_DEBUG
      USE EIGEN_NONLINEAR_CARDS, ONLY :  BD_EIGR, BD_EIGRL, BD_NLPARM
      USE BULK_DATA_LOADS, ONLY       :  BD_LOAD, BD_FORMOM, BD_GRAV, BD_PLOAD2, BD_PLOAD4, BD_RFORCE, BD_SLOAD
      USE MATERIAL_CARDS, ONLY        :  BD_MAT1, BD_MAT2, BD_MAT8, BD_MAT9
      USE CONSTRAINT_CARDS, ONLY      :  BD_SPC, BD_SPC1, BD_SPCADD, BD_MPC, BD_MPCADD
      USE PARAM_CARDS, ONLY           :  BD_PARAM
      USE PARTITION_VECTORS, ONLY     :  BD_PARVEC, BD_PARVEC1
      USE RIGID_ELEMENTS, ONLY        :  BD_RBAR, BD_RBE1, BD_RBE2, BD_RBE3, BD_RSPLINE
      USE BULK_DATA_TEMPERATURES, ONLY:  BD_TEMP, BD_TEMPD, BD_TEMPRP
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE SORTING, ONLY               :  SORT_INT1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'LOADB'
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD1              ! BD card (a small field card or the 1st half of a large field card)
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD2              ! 2nd half of a large field card
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD               ! 16 col field card (either CARD1 if small field or CARD1 + CARD2 if
!                                                             a large field card. This is output from subr FFIELD
      CHARACTER( 1*BYTE)              :: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER( 1*BYTE)              :: CC_MPC_FND         ! 'Y' if B.D. MPC card w/ same set ID (SID) as C.C. MPC = SID
      CHARACTER( 1*BYTE)              :: CC_NLSID_FND(LSUB) ! 'Y' if B.D NLPARM card w/ same set ID (SID) as C.C. NLPARM = SID
      CHARACTER( 1*BYTE)              :: CC_SPC_FND         ! 'Y' if B.D. SPC/SPC1 card w/ same set ID (SID) as C.C. SPC = SID
      CHARACTER( 3*BYTE)              :: CONSTR_TYPE = '   '! Used for output error message (='SPC' or 'MPC')
      CHARACTER( 9*BYTE)              :: DECK_NAME   = 'BULK DATA'
      CHARACTER( 1*BYTE)              :: EIGFND      = 'N'  ! 'Y' if B.D. EIGR card w/ same set ID (SID) as C.C METHOD = SID
      CHARACTER( 1*BYTE)              :: FOUND       = 'N'  ! Used to indicate if we found something we are looking for
      CHARACTER( 1*BYTE)              :: LARGE_FLD_INP      ! If 'Y', card is in large field format

                                                            ! PBARL section names
      CHARACTER(132*BYTE)             :: MESSAG1            ! Messag to print out
      CHARACTER(132*BYTE)             :: MESSAG2            ! Messag to print out
      CHARACTER( 8*BYTE)              :: PBARL_SEC_TYPES(NPBARL)

      CHARACTER( 8*BYTE)              :: SEC_TYPE           ! Cross-section name for a PBARL
      CHARACTER( 1*BYTE)              :: SID_ON_MPCADD_FND  ! 'Y' if B.D. MPC card w/ set ID on MPCADD card found
      CHARACTER( 1*BYTE)              :: SID_ON_SPCADD_FND  ! 'Y' if B.D. SPC/SPC1 card w/ set ID on SPCADD card found
      CHARACTER( 1*BYTE)              :: SID_ON_LOAD_FND    ! 'Y' if B.D. FORCE/MOMENT/GRAV/PLOAD card w/ set ID on SPCADD card fnd
      CHARACTER( 7*BYTE), PARAMETER   :: END_CARD    = 'ENDDATA'

      INTEGER(LONG)                   :: COMMENT_COL        ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: I,J,K,L            ! DO loop indices
      INTEGER(LONG)                   :: IERR               ! Error indicator from subr FFIELD
      INTEGER(LONG)                   :: IOCHK              ! IOSTAT error number when reading Bulk Data cards from unit IN1
      INTEGER(LONG)                   :: IOR3D              ! Integration order from a PSOLID entry
      INTEGER(LONG)                   :: IPBARL             ! Count of number of PBARL entries as they are read
      INTEGER(LONG)                   :: ELEM_NUM_GRDS      ! Number of grids for an elem
      INTEGER(LONG)                   :: ELEM_NUM_DOFS      ! Number of DOF   for an elem
      INTEGER(LONG)                   :: NG                 ! Actual num grids on CUSERIN (not incl SPOINT's)
      INTEGER(LONG)                   :: NS                 ! Actual num SPOINT'ss on CUSERIN
      INTEGER(LONG)                   :: NUM_QUADS          ! Number of quadrilateral elements




! **********************************************************************************************************************************
      ! Initialize
      IOR3D_MAX = 0

      IF (SPCSET == 0) THEN
         CC_SPC_FND = '0'
      ELSE
         CC_SPC_FND = 'N'
      ENDIF

      IF (MPCSET == 0) THEN
         CC_MPC_FND = '0'
      ELSE
         CC_MPC_FND = 'N'
      ENDIF

      ! Initialize CC_LOAD_FND

      !  CC_LOAD_FND = Character array containing 'Y' or 'N' indicating whether the load
      !           or temperature load requested in Case Control was found in B.D.deck.
      !           Each row has the following indicators for 1 subcase
      !            (1) Col 1: indicators for Case Control LOAD
      !            (2) Col 2: indicators for Case Control TEMP
      !            (a) CC_LOAD_FND(I,1) = '0' if there were no LOAD requests in Case Control
      !                            = 'N' if there were LOAD requests in Case Control
      !            (b) CC_LOAD_FND(I,2) = '0' if there were no TEMP requests in Case Control
      !                            = 'N' if there were TEMP requests in Case Control
      !           Subroutines BD_FORMOM, BD_GRAV, BD_LOAD, BD_PLOAD2, BD_TEMP, BD_TEMPD, BD_TEMPRP
      !           reset CC_LOAD_FND to 'Y' if a load (or temp) with set ID matching one in a
      !           Case Control request is found.
      DO I=1,NSUB
         IF (SUBLOD(I,1) == 0) THEN
            CC_LOAD_FND(I,1) = '0'
         ELSE
            CC_LOAD_FND(I,1) = 'N'
         ENDIF
         IF (SUBLOD(I,2) == 0) THEN
            CC_LOAD_FND(I,2) = '0'
         ELSE
            CC_LOAD_FND(I,2) = 'N'
         ENDIF
      ENDDO

! **********************************************************************************************************************************
      CARD(1:) = ' '

      IPBARL    = 0

      ELEM_NUM_GRDS = 0
      NUM_QUADS     = 0
      MELGP         = 2   ! This max num grids/elem DOF's covers all 2 node/6 comp per node elems
      MELDOF        = 12  ! (other elems will be checked and MELGP, MELDOF reset if necessary)

      ! Process Bulk Data cards in a large loop that runs until either an
      ! ENDDATA card is found or when an error or EOF/EOR occurs
bdf:  DO

         CALL READ_BDF_LINE(IN1, IOCHK, CARD1)

         ! Must have this since CARD goes to BD_xxxx, not CARD1.
         ! This will get reset if CARD1 is a large field format
         CARD(1:) = CARD1(1:)

         ! Quit if EOF/EOR occurs.
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         ! Write out BULK DATA card.
         IF (ECHO(1:4) /= 'NONE') THEN
            WRITE(F06,101) CARD1
         ENDIF

         ! Determine if the card is large or small format
         LARGE_FLD_INP = 'N'
         DO I=1,8
            IF (CARD1(I:I) == '*') THEN
               LARGE_FLD_INP = 'Y'
            ENDIF
         ENDDO

         ! FFIELD converts free-field card to fixed field and left justifies
         ! the data in fields 2-9 and outputs a 10 field, 16 col/field CARD
         IF ((CARD1(1:1) /= '$') .AND. (CARD1(1:) /= ' ')) THEN

            IF (LARGE_FLD_INP == 'N') THEN
               CALL FFIELD ( CARD1, IERR )
               CARD(1:) = CARD1(1:)
            ELSE
               ! Read 2nd physical entry for a large field parent B.D. entry
               CALL READ_BDF_LINE(IN1, IOCHK, CARD2)

               IF (IOCHK < 0) THEN
                  WRITE(ERR,1011) END_CARD
                  WRITE(F06,1011) END_CARD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               IF (IOCHK > 0) THEN
                  WRITE(ERR,1010) DECK_NAME
                  WRITE(F06,1010) DECK_NAME
                  WRITE(F06,'(A)') CARD2
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ENDIF

               IF (ECHO(1:4) /= 'NONE') THEN
                  WRITE(F06,101) CARD2
               ENDIF

               IF      (CARD2( 1: 8) == CARD1(73:80)) THEN
                  CONTINUE
               ELSE IF ((CARD2( 1: 1) == '*') .AND. (CARD1(73:73) == ' ') .AND. (CARD2(2:8) == CARD1(74:80))) THEN
                  CONTINUE
               ELSE IF ((CARD2( 1: 1) == ' ') .AND. (CARD1(73:73) == '*') .AND. (CARD2(2:8) == CARD1(74:80))) THEN
                  CONTINUE
               ELSE
                  ! CARD2 is not a continuation of CARD1 so backspace IN1
                  BACKSPACE(IN1)
                  CARD2(1:) = ' '
                  CARD2(1:8) = CARD1(73:80)
               ENDIF

               CALL FFIELD2 ( CARD1, CARD2, CARD, IERR )
            ENDIF

            IF (IERR /= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
!xx            WRITE(ERR,101) CARD
!xx            WRITE(ERR,1003)
!xx            IF (ECHO == 'NONE  ') THEN
!xx               WRITE(F06,101) CARD
!xx            ENDIF
!xx            WRITE(F06,1003)
               CYCLE
            ENDIF

         ENDIF

         ! Process Bulk Data card
         !write(err,*) 'cardname', card
         IF((     CARD(1:5) == 'ASET '   ) .OR. (CARD(1:5) == 'ASET*'   ) .OR.                                                     &
                 (CARD(1:5) == 'OMIT '   ) .OR. (CARD(1:5) == 'OMIT*'   )) THEN
            CALL BD_ASET    ( CARD )

         ELSE IF((CARD(1:5) == 'ASET1'   ) .OR. (CARD(1:5) == 'OMIT1'   )) THEN
            CALL BD_ASET1   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'BAROR'   )  THEN
            CALL BD_BAROR   ( CARD )

         ELSE IF (CARD(1:6) == 'BEAMOR'  )  THEN
            CALL BD_BEAMOR  ( CARD )

         ELSE IF((CARD(1:4) == 'CBAR'    ) .OR. (CARD(1:5) == 'CBEAM'   ))  THEN
            CALL BD_CBAR    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'CBUSH'   )  THEN
            CALL BD_CBUSH   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CELAS1'  )  THEN
            CALL BD_CELAS1  ( CARD )

         ELSE IF (CARD(1:6) == 'CELAS2'  )  THEN
            CALL BD_CELAS2  ( CARD )

         ELSE IF (CARD(1:6) == 'CELAS3'  )  THEN
            CALL BD_CELAS3  ( CARD )

         ELSE IF (CARD(1:6) == 'CELAS4'  )  THEN
            CALL BD_CELAS4  ( CARD )

         ELSE IF (CARD(1:5) == 'CHEXA'   ) THEN
            CALL BD_CHEXA   ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CMASS1'  )  THEN
            CALL BD_CMASS1  ( CARD )

         ELSE IF (CARD(1:6) == 'CMASS2'  )  THEN
            CALL BD_CMASS2  ( CARD )

         ELSE IF (CARD(1:6) == 'CMASS3'  )  THEN
            CALL BD_CMASS3  ( CARD )

         ELSE IF (CARD(1:6) == 'CMASS4'  )  THEN
            CALL BD_CMASS4  ( CARD )

         ELSE IF (CARD(1:6) == 'CONROD'  )  THEN
            CALL BD_CONROD  ( CARD )

         ELSE IF (CARD(1:5) == 'CONM2'   )  THEN
            CALL BD_CONM2   ( CARD, LARGE_FLD_INP )

         ELSE IF ((CARD(1:6) == 'CORD1C'  ) .OR. (CARD(1:6) == 'CORD1R'  ) .OR. (CARD(1:6) == 'CORD1S'  ) .OR.                     &
                  (CARD(1:6) == 'CORD2C'  ) .OR. (CARD(1:6) == 'CORD2R'  ) .OR. (CARD(1:6) == 'CORD2S'  )) THEN
            CALL BD_CORD    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CPENTA'  ) THEN
            CALL BD_CPENTA  ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CQUAD4'  ) THEN
            NUM_QUADS = NUM_QUADS + 1
            CALL BD_CQUAD   ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP   = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CQUAD8'  ) THEN
            CALL BD_CQUAD8   ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP   = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:4) == 'CROD'    )  THEN
            CALL BD_CROD    ( CARD )

         ELSE IF (CARD(1:6) == 'CSHEAR'  )  THEN
            NUM_QUADS = NUM_QUADS + 1
            CALL BD_CSHEAR  ( CARD, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP   = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CTETRA'  ) THEN
            CALL BD_CTETRA  ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CTRIA3'  ) THEN
            CALL BD_CTRIA   ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP   = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:6) == 'CUSER1'  )  THEN
            CALL BD_CUSER1  ( CARD, LARGE_FLD_INP, ELEM_NUM_GRDS )
            ELEM_NUM_DOFS = 6*ELEM_NUM_GRDS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:7) == 'CUSERIN' )  THEN
            CALL BD_CUSERIN ( CARD, LARGE_FLD_INP, NG, NS )
            ELEM_NUM_GRDS =   NG + NS
            ELEM_NUM_DOFS = 6*NG + NS
            IF (MELGP < ELEM_NUM_GRDS) THEN
               MELGP = ELEM_NUM_GRDS
            ENDIF
            IF (MELDOF < ELEM_NUM_DOFS) THEN
               MELDOF = ELEM_NUM_DOFS
            ENDIF

         ELSE IF (CARD(1:5) == 'DEBUG'   )  THEN
            CALL BD_DEBUG   ( CARD )

         ELSE IF((CARD(1:5) == 'EIGR '   ) .OR. (CARD(1:5) == 'EIGR*'   ))  THEN
            CALL BD_EIGR    ( CARD, LARGE_FLD_INP, EIGFND )

         ELSE IF (CARD(1:5) == 'EIGRL'   )  THEN
            CALL BD_EIGRL   ( CARD, LARGE_FLD_INP, EIGFND )

         ELSE IF((CARD(1:5) == 'FORCE'   ) .OR. (CARD(1:6) == 'MOMENT'  )) THEN
            CALL BD_FORMOM  ( CARD, CC_LOAD_FND )

         ELSE IF (CARD(1:4) == 'GRAV'    )  THEN
            CALL BD_GRAV    ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

         ELSE IF (CARD(1:6) == 'GRDSET'  )  THEN
            CALL BD_GRDSET  ( CARD )

         ELSE IF (CARD(1:4) == 'GRID'    )  THEN
            CALL BD_GRID    ( CARD )

         ELSE IF (CARD(1:4) == 'LOAD'    )  THEN
            CALL BD_LOAD    ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

         ELSE IF (CARD(1:4) == 'MAT1'    )  THEN
            CALL BD_MAT1    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'MAT2'    )  THEN
            CALL BD_MAT2    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'MAT8'    )  THEN
            CALL BD_MAT8    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'MAT9'    )  THEN
            CALL BD_MAT9    ( CARD, LARGE_FLD_INP )

         ELSE IF((CARD(1:4) == 'MPC '    ) .OR. (CARD(1:4) == 'MPC*'    ))  THEN
            CALL BD_MPC     ( CARD, LARGE_FLD_INP, CC_MPC_FND )

         ELSE IF (CARD(1:6) == 'MPCADD'  )  THEN
            CALL BD_MPCADD  ( CARD, LARGE_FLD_INP, CC_MPC_FND )

         ELSE IF (CARD(1:6) == 'NLPARM'  )  THEN
            IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
               CALL BD_NLPARM  ( CARD, CC_NLSID_FND )
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,9993) PROG_NAME
               IF (SUPWARN == 'N') THEN
                  IF (ECHO(1:4) /= 'NONE') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,9993) PROG_NAME
               ENDIF
            ENDIF

         ELSE IF (CARD(1:5) == 'PARAM'   )  THEN
            CALL BD_PARAM   ( CARD )

         ELSE IF((CARD(1:7) == 'PARVEC ' ) .OR. (CARD(1:7) == 'PARVEC*'  )) THEN
            CALL BD_PARVEC   ( CARD )

         ELSE IF((CARD(1:8) == 'PARVEC1 ')  .OR. (CARD(1:8) == 'PARVEC1*'))THEN
            CALL BD_PARVEC1  ( CARD, LARGE_FLD_INP )

         ELSE IF((CARD(1:5) == 'PBAR '   ) .OR. (CARD(1:5) == 'PBAR*'   ))  THEN
            CALL BD_PBAR    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'PBARL'   )  THEN
            CALL BD_PBARL    ( CARD, LARGE_FLD_INP, SEC_TYPE )
            IPBARL = IPBARL + 1
            PBARL_SEC_TYPES(IPBARL) = SEC_TYPE

         ELSE IF (CARD(1:5) == 'PBEAM'   )  THEN
            CALL BD_PBEAM   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'PBUSH'   )  THEN
            CALL BD_PBUSH   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'PCOMP'   )  THEN
            CALL BD_PCOMP   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'PCOMP1'  )  THEN
            CALL BD_PCOMP1  ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'PELAS'   )  THEN
            CALL BD_PELAS   ( CARD )

         ELSE IF (CARD(1:6) == 'PLOAD4'  )  THEN
            CALL BD_PLOAD4  ( CARD, CC_LOAD_FND )

         ELSE IF (CARD(1:6) == 'PLOAD2'  )  THEN
            CALL BD_PLOAD2  ( CARD, CC_LOAD_FND )

         ELSE IF (CARD(1:6) == 'PLOTEL'  )  THEN
            CALL BD_PLOTEL  ( CARD )

         ELSE IF (CARD(1:5) == 'PMASS'   )  THEN
            CALL BD_PMASS   ( CARD )

         ELSE IF (CARD(1:4) == 'PROD'    )  THEN
            CALL BD_PROD    ( CARD )

         ELSE IF (CARD(1:6) == 'PSHEAR'  )  THEN
            CALL BD_PSHEAR  ( CARD )

         ELSE IF (CARD(1:6) == 'PSHELL'  )  THEN
            CALL BD_PSHEL   ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'PSOLID'  )  THEN
            CALL BD_PSOLID  ( CARD, IOR3D )
            IF (IOR3D > IOR3D_MAX) THEN
               IOR3D_MAX = IOR3D
            ENDIF

         ELSE IF (CARD(1:6) == 'PUSER1'  )  THEN
            CALL BD_PUSER1  ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:7) == 'PUSERIN' )  THEN
            CALL BD_PUSERIN ( CARD )

         ELSE IF (CARD(1:4) == 'RBAR'    )  THEN
            CALL BD_RBAR    ( CARD )

         ELSE IF (CARD(1:4) == 'RBE1'    )  THEN
            CALL BD_RBE1    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'RBE2'    )  THEN
            CALL BD_RBE2    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'RBE3'    )  THEN
            CALL BD_RBE3    ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'RFORCE'  )  THEN
            CALL BD_RFORCE  ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

         ELSE IF (CARD(1:7) == 'RSPLINE' )  THEN
            CALL BD_RSPLINE ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'SEQGP'   )  THEN
            IF (GRIDSEQ(1:6) == 'BANDIT') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1021) GRIDSEQ
               IF (SUPWARN == 'N') THEN
                  IF (ECHO(1:4) /= 'NONE') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1021) GRIDSEQ
               ENDIF
               CYCLE
            ELSE
               CALL BD_SEQGP( CARD )
            ENDIF

         ELSE IF (CARD(1:5) == 'SLOAD'   )  THEN
            CALL BD_SLOAD   ( CARD, CC_LOAD_FND )

         ELSE IF (CARD(1:5) == 'SNORM'   )  THEN
            CALL BD_SNORM   ( CARD )

         ELSE IF((CARD(1:4) == 'SPC '    ) .OR. (CARD(1:4) == 'SPC*'    )) THEN
            CALL BD_SPC     ( CARD, CC_SPC_FND )

         ELSE IF (CARD(1:4) == 'SPC1'    )  THEN
            CALL BD_SPC1    ( CARD, LARGE_FLD_INP, CC_SPC_FND )

         ELSE IF (CARD(1:6) == 'SPCADD'  )  THEN
            CALL BD_SPCADD  ( CARD, LARGE_FLD_INP, CC_SPC_FND )

         ELSE IF (CARD(1:6) == 'SPOINT'  )  THEN
            CALL BD_SPOINT  ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'SUPORT'  )  THEN
            CALL BD_SUPORT  ( CARD )

         ELSE IF((CARD(1:5) == 'TEMP '   ) .OR. (CARD(1:5) == 'TEMP*'   ))  THEN
            CALL BD_TEMP    ( CARD, CC_LOAD_FND )

         ELSE IF (CARD(1:5) == 'TEMPD'   )  THEN
            CALL BD_TEMPD   ( CARD, CC_LOAD_FND )

         ELSE IF((CARD(1:6) == 'TEMPRB'  ) .OR. (CARD(1:6) == 'TEMPP1'  )) THEN
            CALL BD_TEMPRP  ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

         ELSE IF((CARD(1:5) == 'USET '   ) .OR. (CARD(1:5) == 'USET*'    )) THEN
            CALL BD_USET     ( CARD )

         ELSE IF((CARD(1:6) == 'USET1 '  )  .OR. (CARD(1:6) == 'USET1*'  ))THEN
            CALL BD_USET1    ( CARD, LARGE_FLD_INP )

         ELSE IF ((CARD(1:1) == '$') .OR. (CARD(1:BD_ENTRY_LEN) == ' ')) THEN
            CYCLE

         ELSE IF (CARD(1:7) == 'ENDDATA' )  THEN
            EXIT
         ELSE IF ((CARD(1:1) == ' ') .OR. (CARD(1:1) == '+') .OR. (CARD(1:1) == '*'))  THEN
            !WRITE(ERR,*) 'FAILED WHEN FINDING A CONTINUATION'
            !WRITE(F06,*) 'FAILED WHEN FINDING A CONTINUATION'
            ! only defined when it's a large field continuation
            !WRITE(ERR,'(A)') CARD2
            !WRITE(F06,'(A)') CARD2

         ELSE                                              ! CARD not processed by MYSTRAN
            WARN_ERR = WARN_ERR + 1
            write(err,*) 'card name not found...'
            WRITE(ERR,101) CARD
            WRITE(ERR,9993) PROG_NAME
            flush(err)
            IF (SUPWARN == 'N') THEN
               IF (ECHO(1:4) /= 'NONE') THEN
                  WRITE(F06,101) CARD
               ENDIF
               WRITE(F06,9993) PROG_NAME
            ENDIF

         ENDIF

      ENDDO bdf

! **********************************************************************************************************************************
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ! Check to make sure there are no elemenmts not coded for DIFFEREN or BUCKLING solutions
      !  If there are, tell user that they have to recognize this.
      !  If they want to continue put DEBUG 201 w/val /= 0 entry in the BDF

      IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF((NCSHEAR   > 0)) THEN
            MESSAG1= ' *WARNING: BUCKLING and DIFFERN SOL are only coded for the BAR, HEXA, PENTA, TETRA, TRIA, and QUAD elements'
            MESSAG2= '           Either remove all other elements or include a Bulk Data entry: DEBUG   201, with value /= 0'
            IF (DEBUG(201) == 0) THEN
               WRITE(F06,*) MESSAG1
               WRITE(F06,*) MESSAG2
               WRITE(ERR,*) MESSAG1
               WRITE(ERR,*) MESSAG2
               CALL OUTA_HERE ( 'Y' )
            ELSE
               WRITE(F06,9991)
               WRITE(ERR,9991)
            ENDIF
         ENDIF
      ENDIF
 9991 FORMAT( /,' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'&
            ,'+++++++++++++++++++',/                                                                                               &
            ,' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'&
            ,'++++++++++++++++',/,' ++',127X,'++',/                                                                                &
             ' ++                                                          ___    _____    ____                     ' ,28X,'++',/  &
             ' ++                            * * * * * * * * * *  |\  |   /   \     |     |      * * * * * * * * * *' ,28X,'++',/  &
             ' ++                            * * * * * * * * * *  | \ |   | O |     |     |---   * * * * * * * * * *' ,28X,'++',/  &
             ' ++                            * * * * * * * * * *  |  \|   \___/     |     |____  * * * * * * * * * *' ,28X,'++',/  &
             ' ++',127X,'++',/,                                                                                                    &
             ' ++',127X,'++',/,                                                                                                    &
   ' ++', 6X,' B U C K L I N G   &   D I F F E R E N   S O L   a r e   o n l y   c o d e d   f o r   e l e m e n t s:',18X,'++',/  &
             ' ++',127X,'++',/,                                                                                                    &
             ' ++                                BAR,  HEXA,  PENTA,  TETRA,  TRIA3,  TRIA3K,  QUAD4,  QUAD4K       ' ,28X,'++',/  &
             ' ++',127X,'++',/                                                                                                     &
   ' ++', 6X,'   A l l   o t h e r   e l e m s   i g n o r e d   i n   t h e  d i f f e r e n t i a l   s t i f f   c a l c s',    &
  10X,'++',/,' ++',127X,'++',/                                                                                                     &
             ' ++',127X,'++',/                                                                                                     &
            ,' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'&
            ,'++++++++++++++++',/                                                                                                  &
            ,' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'&
            ,'++++++++++++++++')
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      ! Check that, if there were PARTN requests in Exec Control
      ! (i.e. NUM_PARTN_REQUESTS > 0) that PARVEC,1 records were in Bulk Data
      IF ((NUM_PARTN_REQUESTS > 0) .AND. (NUM_PARTVEC_RECORDS == 0)) THEN
         WRITE(ERR,1805) NUM_PARTN_REQUESTS
         WRITE(F06,1805) NUM_PARTN_REQUESTS
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      ! Write PBAR equivalent props for the PBARL entries
      IF (NPBARL > 0) THEN

         IF (DEBUG(113) > 0) THEN
            IPBARL = 0
            WRITE(F06,*)
            WRITE(F06,1028)
            DO I=1,NPBAR
               IF (PBAR(I,3) > 0) THEN
                  IPBARL = IPBARL + 1
                  WRITE(F06,1029) PBARL_SEC_TYPES(IPBARL), PBAR(I,1)  ,PBAR(I,2)  ,(RPBAR(I,J),J=1,16)
               ENDIF
            ENDDO
         ENDIF
         WRITE(F06,*)

         IF (ECHO(1:4) /= 'NONE') THEN
            IF (DEBUG(113) == 0) THEN
               WRITE(F06,*)
               WRITE(F06,1025)
               IPBARL = 0
               DO I=1,NPBAR
                  IF (PBAR(I,3) > 0) THEN
                     IPBARL = IPBARL + 1
                     WRITE(F06,1026) PBAR(I,1)  ,PBAR(I,2)  ,RPBAR(I, 1),RPBAR(I, 2),RPBAR(I, 3),RPBAR(I, 4),RPBAR(I, 5),          &
                                     PBARL_SEC_TYPES(IPBARL)
                     WRITE(F06,1027) RPBAR(I, 6),RPBAR(I, 7),RPBAR(I, 8),RPBAR(I, 9),RPBAR(I,10),RPBAR(I,11),RPBAR(I,12),RPBAR(I,13)
                     WRITE(F06,1027) RPBAR(I,14),RPBAR(I,15),RPBAR(I,16)
                     WRITE(F06,*)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

      ENDIF

      ! Set MOFFSET based on the element type that requires the most offset points
      IF ((NCBAR   > 0) .OR. (NCBEAM  > 0) .OR. (NCBUSH   > 0) .OR. (NCROD    > 0)) MOFFSET = 2
      IF ((NCTRIA3 > 0) .OR. (NCTRIA3K > 0)) MOFFSET = 3
      IF ((NCQUAD4 > 0) .OR. (NCQUAD4K > 0)) MOFFSET = 4
      IF  (NCQUAD8 > 0) MOFFSET = 8

      ! Determine max num of indep grids on MPC's and rigid elems
      ! (for dimensioning array MPC_IND_GRIDS)
      LIND_GRDS_MPCS = NMPC*MMPC + 2*NRBAR + 6*NRBE1 + NRBE2 + 2*NRSPLINE

      ! Determine the max number of BEi, SEi strain/stress recovery matrices
      ! will be needed for this execution
      CALL CALC_MAX_STRESS_POINTS

      ! Determine the max number of Gauss points that will be used in this execution
      CALL CALC_MAX_GAUSS_POINTS

      ! Give warning if there are any SPOINT's
      IF (NSPOINT > 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9994) NSPOINT
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9994) NSPOINT
         ENDIF
      ENDIF

      ! Give error if more than 1 BAROR, BEAMOR or GRDSET card was in
      ! Bulk Data (these counted by BD_BAROR0, BD_BEAMOR0, BD_GRDSET0)
      IF (NBAROR > 1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1022) 'BAROR'
         WRITE(F06,1022) 'BAROR'
      ENDIF

      IF (NBEAMOR > 1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1022) 'BEAMOR'
         WRITE(F06,1022) 'BEAMOR'
      ENDIF

      IF (NGRDSET > 1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1023)
         WRITE(F06,1023)
      ENDIF

! If there are quad elements give message on how the local element x axis will be determined

      IF (NUM_QUADS > 0) THEN
         IF (QUADAXIS == 'SPLITD') THEN
            WRITE(ERR,1101) QUADAXIS
            IF (SUPINFO == 'N') THEN
               WRITE(F06,1101) QUADAXIS
            ENDIF
         ENDIF
         IF (QUADAXIS == 'SIDE12') THEN
            WRITE(ERR,1102) QUADAXIS
            IF (SUPINFO == 'N') THEN
               WRITE(F06,1102) QUADAXIS
            ENDIF
         ENDIF
      ENDIF

! **********************************************************************************************************************************
      ! If SOL_NAME is modes or CB then an EIGR card should have been found with SID matching a SID in Case Control.
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         IF (EIGFND == 'N') THEN
            WRITE(ERR,1005) CC_EIGR_SID
            WRITE(F06,1005) CC_EIGR_SID
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

      ! Check SPC entries
      ! -----------------

      ! (a) If SPC was requested in CASE CONTROL, we should have found a Bulk Data SPC, SPC1 or MPCADD with SID matching C.C SPC SID
      IF (CC_SPC_FND == 'N') THEN
         CONSTR_TYPE = 'SPC'
         WRITE(ERR,1006) CONSTR_TYPE,SPCSET
         WRITE(F06,1006) CONSTR_TYPE,SPCSET
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      ! (b) If there are SPCADD's check to make sure that all SPC, SPC1 Bulk Data requested on them are in the deck
      DO I=1,NSPCADD
         IF (SPCADD_SIDS(I,1) == SPCSET) THEN
            DO K=2,LSPCADDC
               IF (SPCADD_SIDS(I,K) == 0) CYCLE
               SID_ON_SPCADD_FND = 'N'
               DO L=1,NSPC
                  IF (SPCADD_SIDS(I,K) == SPC_SIDS(L)) THEN
                     SID_ON_SPCADD_FND = 'Y'
                  ENDIF
               ENDDO
               DO L=1,NSPC1
                  IF (SPCADD_SIDS(I,K) == SPC1_SIDS(L)) THEN
                     SID_ON_SPCADD_FND = 'Y'
                  ENDIF
               ENDDO
               IF (SID_ON_SPCADD_FND == 'N') THEN
                  WRITE(ERR,1018) SPCADD_SIDS(I,K),SPCADD_SIDS(I,1)
                  WRITE(F06,1018) SPCADD_SIDS(I,K),SPCADD_SIDS(I,1)
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDDO
         ENDIF
      ENDDO

! Check MPC entries
! -----------------

! (a) If MPC was requested in CASE CONTROL, we should have found a Bulk Data MPC or MPCADD with a SID matching Case Control MPC SID

      IF (CC_MPC_FND == 'N') THEN
         CONSTR_TYPE = 'MPC'
         WRITE(ERR,1006) CONSTR_TYPE,MPCSET
         WRITE(F06,1006) CONSTR_TYPE,MPCSET
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! (b) If there are MPCADD's check to make sure that all MPC Bulk Data requested on them are in the deck

      DO I=1,NMPCADD
         IF (MPCADD_SIDS(I,1) == MPCSET) THEN
            DO K=2,LMPCADDC
               IF (MPCADD_SIDS(I,K) == 0) CYCLE
               SID_ON_MPCADD_FND = 'N'
               DO L=1,NMPC
                  IF (MPCADD_SIDS(I,K) == MPC_SIDS(L)) THEN
                     SID_ON_MPCADD_FND = 'Y'
                  ENDIF
               ENDDO
               IF (SID_ON_MPCADD_FND == 'N') THEN
                  WRITE(ERR,1015) MPCADD_SIDS(I,K),MPCADD_SIDS(I,1)
                  WRITE(F06,1015) MPCADD_SIDS(I,K),MPCADD_SIDS(I,1)
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDDO
         ENDIF
      ENDDO

! (c) If there are MPC's, create array MPCSIDS (which has 1 entry if there is any B.D. MPC and > 1 entry if there are MPCADD's.
!     The first one in the array is the MPC set ID called for in Case Control. Subsequent ones are the set ID's from any
!     MPCADD cards that have the set ID called for in Case Control. So, if there are 2 MPCADD's that have the same ID as
!     called for in Case Control and there are (for example) 12 individual MPC set ID's identified on these
!     MPCADD cards, then MPCSIDS will contain 13 entries.

      IF ((NMPC > 0) .OR. (NMPCADD > 0)) THEN

         FOUND = 'N'                                       ! FOUND will indicate if there is are MPCADD's with set ID = MPCSET
         NUM_MPCSIDS = 1                                   ! NUM_MPCSIDS will be 1 if there are any MPC's. It will be > 1 if there
!                                                            are MPCADD entries that match the MPC set requested in Case Control
         DO I=1,NMPCADD
            IF (MPCADD_SIDS(I,1) == MPCSET) THEN
               FOUND = 'Y'
j_do1:         DO J=2,LMPCADDC
                  IF (MPCADD_SIDS(I,J) /= 0) THEN
                     NUM_MPCSIDS = NUM_MPCSIDS + 1
                     CYCLE j_do1
                  ELSE
                     EXIT j_do1
                  ENDIF
               ENDDO j_do1
            ENDIF
         ENDDO

         CALL ALLOCATE_MODEL_STUF ( 'MPCSIDS', SUBR_NAME ) ! Allocate enough memory for array MPCSIDS
         MPCSIDS(1) = MPCSET

         IF (FOUND == 'Y') THEN                            ! There are MPCADD entries so count the MPC set ID's on those MPCADD's
!                                                            whose set ID matches Case Control request
            K = 1                                          ! Fill array MPCSIDS, as described above, with SID's that will be used
            DO I=1,NMPCADD
               IF (MPCADD_SIDS(I,1) == MPCSET) THEN
j_do2:            DO J=2,LMPCADDC
                     IF (MPCADD_SIDS(I,J) /= 0) THEN
                        K = K + 1
                        MPCSIDS(K) = MPCADD_SIDS(I,J)
                        CYCLE j_do2
                     ELSE
                        EXIT j_do2
                     ENDIF
                  ENDDO j_do2
               ENDIF
            ENDDO


            ! (d) Check for duplicate MPC set ID's on MPCADD. There was a
            !     check on each MPCADD entry read in subr BD_MPCADD but no
            !     check was made across duplicate MPCADD entries (i.e.
            !     there may be more than 1 MPCADD entry with the same set ID
            !     and we need to make sure that an MPC set ID on the 2nd, etc.,
            !     is not the same as one on the 1st, etc)
            CALL SORT_INT1 ( SUBR_NAME, 'MPCSIDS', NUM_MPCSIDS, MPCSIDS )
            DO I=2,NUM_MPCSIDS
               IF (MPCSIDS(I) == MPCSIDS(I-1)) THEN
                  WRITE(ERR,1024) MPCSIDS(I), MPCSET, MPCSET
                  WRITE(F06,1024) MPCSIDS(I), MPCSET, MPCSET
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDDO

         ENDIF

      ENDIF

      ! Check if load (LOAD, FORCE, MOMENT, GRAV, TEMP, etc.) Bulk Data cards
      ! with SID matching Case Control were found
      DO I=1,NSUB
        IF (CC_LOAD_FND(I,1) == 'N') THEN
           WRITE(ERR,1007) SUBLOD(I,1)
           WRITE(F06,1007) SUBLOD(I,1)
           FATAL_ERR = FATAL_ERR + 1
        ENDIF
        IF (CC_LOAD_FND(I,2) == 'N') THEN
           WRITE(ERR,1008) SUBLOD(I,2)
           WRITE(F06,1008) SUBLOD(I,2)
           FATAL_ERR = FATAL_ERR + 1
        ENDIF
      ENDDO

      ! Check to make sure that all forces, moments and GRAV cards requested
      ! on LOAD Bulk Data cards are in the deck
      DO I=1,NLOAD
         DO J=1,NSUB
            IF (LOAD_SIDS(I,1) == SUBLOD(J,1)) THEN
               DO K=2,LLOADC
                  IF (LOAD_SIDS(I,K) == 0) CYCLE
                  SID_ON_LOAD_FND = 'N'
                  DO L=1,NFORCE
                     IF (LOAD_SIDS(I,K) == FORMOM_SIDS(L)) THEN
                        SID_ON_LOAD_FND = 'Y'
                     ENDIF
                  ENDDO
                  DO L=1,NGRAV
                     IF (LOAD_SIDS(I,K) == GRAV_SIDS(L)) THEN
                        SID_ON_LOAD_FND = 'Y'
                     ENDIF
                  ENDDO
                  DO L=1,NPLOAD
                     IF (LOAD_SIDS(I,K) == PRESS_SIDS(L)) THEN
                        SID_ON_LOAD_FND = 'Y'
                     ENDIF
                  ENDDO
                  DO L=1,NRFORCE
                     IF (LOAD_SIDS(I,K) == RFORCE_SIDS(L)) THEN
                        SID_ON_LOAD_FND = 'Y'
                     ENDIF
                  ENDDO
                  DO L=1,NSLOAD
                     IF (LOAD_SIDS(I,K) == SLOAD_SIDS(L)) THEN
                        SID_ON_LOAD_FND = 'Y'
                     ENDIF
                  ENDDO
                  IF (SID_ON_LOAD_FND == 'N') THEN
                     WRITE(ERR,1009) LOAD_SIDS(I,K),LOAD_SIDS(I,1),SCNUM(J)
                     WRITE(F06,1009) LOAD_SIDS(I,K),LOAD_SIDS(I,1),SCNUM(J)
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      ! Write message regarding max number of elem grids/DOF's
      MDT = MAX(MTDAT_TEMPP1, MTDAT_TEMPRB, MELGP+3, 5)
      WRITE(F06,*)
      WRITE(ERR,1197) MELGP, MELDOF
      IF (SUPINFO == 'N') THEN
         WRITE(F06,1197) MELGP, MELDOF
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT('ECHO: loadb ', A)

 1003 FORMAT(' *ERROR  1003: ALL FIELDS ON THE ABOVE ENTRY MUST BE NO LONGER THAN 8 CHARACTERS')

 1005 FORMAT(' *ERROR  1005: NO EIGENVALUE ENTRY WAS FOUND IN BULK DATA DECK MATCHING SID = ',I8,' REQUESTED IN CASE CONTROL')

 1006 FORMAT(' *ERROR  1006: NO ',A3,' ENTRY WAS FOUND IN BULK DATA DECK MATCHING SID = ',I8,' REQUESTED IN CASE CONTROL')

 1007 FORMAT(' *ERROR  1007: NO LOADING ENTRY WAS FOUND IN BULK DATA DECK MATCHING SID = ',I8,' REQUESTED IN CASE CONTROL')

 1008 FORMAT(' *ERROR  1008: NO TEMPERATURE ENTRY WAS FOUND IN BULK DATA DECK MATCHING SID = ',I8,' REQUESTED IN CASE CONTROL')

 1009 FORMAT(' *ERROR  1009: MISSING FORCE, MOMENT, GRAV OR PLOAD ENTRY WITH SID =',I8,' ON LOAD ENTRY WITH SID = ',I8,            &
                           ' FOR SUBCASE ',I8)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

 1012 FORMAT(' *ERROR  1012: IMBEDDED BLANKS FOUND IN FIELD ',I2,' NOT ALLOWED')

 1015 FORMAT(' *ERROR  1015: MISSING MPC OR MPC1 ENTRY WITH SID =',I8,' REQUESTED ON MPCADD ENTRY SID = ',I8)

 1018 FORMAT(' *ERROR  1018: MISSING SPC OR SPC1 ENTRY WITH SID =',I8,' REQUESTED ON SPCADD ENTRY SID = ',I8)

 1019 FORMAT(' *ERROR  1019: THE GRID TO BE USED AS REFERENCE IN THE GRID POINT EQUIL CHECK, ',I8,', CANNOT BE AN SPOINT')

 1021 FORMAT(' *WARNING    : WHEN PARAM GRIDSEQ = ',A,' SEQGP ENTRIES IN THE BULK DATA DECK ARE NOT ALLOWED. ENTRY IGNORED.')

 1022 FORMAT(' *ERROR  1022: ONLY ONE ',A,' ENTRY ALLOWED IN DATA DECK')

 1023 FORMAT(' *ERROR  1023: ONLY ONE GRDSET ENTRY ALLOWED IN DATA DECK.')

 1024 FORMAT(' *ERROR  1024: MPC SET NUMBER ',I8,' ON MPCADD ',I8,' IS A DUPLICATE AMONG ALL SET IDs ON MPCADD ',I8,' ENTRIES')

 1025 FORMAT(1X,'THE PBARL BULK DATA ENTRIES IN THE BULK DATA DECK ARE REPLACED BY THE FOLLOWING PBAR ENTRIES:',/,                 &
             1x,'--------------------------------------------------------------------------------------------')

 1026 FORMAT(' PBAR    ',2I12,5(1ES12.4),'                cross-section type: ',A)

 1027 FORMAT(9X,8(1ES12.4))

 1028 FORMAT('------------------------------------------------------------------------------------------------',                   &
             'P B A R L   P R O P E R T I E S',                                                                                    &
             '-------------------------------------------------------------------------------------------------',/,                &
             'Sec Type     Prop ID     Matl ID     Area         I1          I2          J          NSM          Y1          Z1    '&
             , '      Y2          Z2          Y3          Z3          Y4          Z4          K1          K2         I12')

 1029 FORMAT(A8,2I12,16(1ES12.4))

 1101 FORMAT(' *INFORMATION: BASED ON PARAM QUADAXIS = ',A,' THE LOCAL X AXIS OF QUAD4 AND SHEAR ELEMENTS WILL BE DETERMINED BY'   &
                    ,/,14X,' SPLITTING THE ANGLE BETWEEN THE 2 DIAGONALS OF THE ELEMENT')

 1102 FORMAT(' *INFORMATION: BASED ON PARAM QUADAXIS = ',A,' THE LOCAL X AXIS OF QUAD4 AND SHEAR ELEMENTS WILL BE DETERMINED BY'   &
                    ,/,14X,' THE LINE FROM THE 1ST TO 2ND NODES ON THE ELEMENT CONNECTION ENTRY')

 1197 FORMAT(' *INFORMATION: NUMBER OF GRID''s (INCL SPOINT''s) FOR ANY ELEMENT IN THIS MODEL IS     <= ',I12,/,                   &
             ' *INFORMATION: NUMBER OF DOF''s  FOR ANY ELEMENT IN THIS MODEL IS                     <= ',I12,/)
 1198 format(32767(i14))

 1199 format(32767(1es14.6))

 1805 FORMAT(' *ERROR  1805: THERE WERE ',I8,' PARTN REQUEST(S) FOR PARTITIONING OUTPUT4 MATRICES IN EXEC CONTROL BUT NO',         &
                           ' PARTITIONING'                                                                                         &
                    ,/,14X,' VECTORS (BULK DATA PARVEC OR PARVEC1 ENTRIES) WERE FOUND IN THE BULK DATA DECK')

 9993 FORMAT(' *LOADB-WARNING    : PRIOR ENTRY NOT PROCESSED BY ',A)

 9994 FORMAT(' *WARNING    : Due to the presence of ',I8,' scalar points (SPOINT''s) the user should be aware of the following:'   &
                    ,/,14X,'    a) They have no geometry; however their displ, forces, etc are reported in F06 as T1 components'   &
                    ,/,14X,'    b) As a consequence of (a), displ components 2-6 are undefined for SPOINT''s'                      &
                    ,/,14X,'    c) SPOINT''s will be treated as if their basic and global coordinates are all zero'                &
                    ,/,14X,'    d) SPOINT''s will not have loads due to GRAV, RFORCE, or TEMP Bulk Data entries'                   &
                    ,/,14X,'    e) All elements except ELASi (including rigid elements) connected to SPOINT''s will result in a',  &
                                 ' fatal error'                                                                                    &
                    ,/,14X,'    f) If CONMi''s are connected to SPOINT''s the offset and moment of inertia terms are ignored'      &
                    ,/,14X,'    g) Terms in matrix RBGLOBAL (rigid body displ matrix) for SPOINT''s are zero which will probably'  &
                    ,/,14x,'       result in incorrect stiffness matrix equilibrium checks for models with SPOINT''s')





! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CALC_MAX_STRESS_POINTS

      USE SCONTR, ONLY                :  METYPE
      USE MODEL_STUF, ONLY            :  NUM_SEi

      IMPLICIT NONE

! **********************************************************************************************************************************

      MAX_STRESS_POINTS = 1
      DO I=1,METYPE
         IF (NUM_SEi(I) > MAX_STRESS_POINTS) THEN
            MAX_STRESS_POINTS = NUM_SEi(I)
         ENDIF
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE CALC_MAX_STRESS_POINTS

! ##################################################################################################################################

      SUBROUTINE CALC_MAX_GAUSS_POINTS

      USE SCONTR, ONLY                :  METYPE
      USE MODEL_STUF, ONLY            :  NUM_SEi

      IMPLICIT NONE

! **********************************************************************************************************************************

      MAX_GAUSS_POINTS = MAX ( IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T, IOR3D_MAX )

! **********************************************************************************************************************************

      END SUBROUTINE CALC_MAX_GAUSS_POINTS

      END SUBROUTINE LOADB

      !---------------------------------------------------------
      SUBROUTINE READ_BDF_LINE(IN1, IOCHK, LINE)

      USE IOUNT1, ONLY                :  ERR, INFILE, F06 !
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE TEXT_FIELD_UTILS, ONLY      :  TO_UPPER

      USE FILE_LIFECYCLE, ONLY        :  READERR
      IMPLICIT NONE

      CHARACTER(256), INTENT (INOUT)  :: LINE
      CHARACTER(256)                  :: TRIM_LINE
      INTEGER, INTENT (IN)            :: IN1
      INTEGER, INTENT (INOUT)         :: IOCHK
      INTEGER                         :: I, FATAL_ERR
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06
      MESSAG = 'BULK DATA CARD          '
      FATAL_ERR = 0

      READ(IN1,101,IOSTAT=IOCHK) LINE

      TRIM_LINE = TRIM(LINE)
      DO WHILE(TRIM_LINE(1:1) == '$')
         IF (IOCHK /= 0) THEN
           REC_NO = -99
           CALL READERR (IOCHK, INFILE, MESSAG, REC_NO, OUNT )
           FATAL_ERR = FATAL_ERR + 1
         ENDIF
         READ(IN1,101,IOSTAT=IOCHK) LINE
         TRIM_LINE = TRIM(LINE)
      ENDDO

      DO I = 1,256
          ! remove $
          IF (LINE(I:I) == '$') THEN
              LINE(I:) = ' '
              EXIT
          ENDIF
      ENDDO

      LINE = TO_UPPER(LINE)


  101 FORMAT(A)
      END SUBROUTINE READ_BDF_LINE


      SUBROUTINE LOADB0

      ! Preliminary reading of the Bulk Data to count several data sizes so
      ! that arrays may be allocated prior to the final reading of the Bulk Data.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, LCMASS, LDOFG, LELE,                               &
                                         LEDAT, LFORCE, LCONM2, LCORD, LGRAV, LGRID, LGUSERIN, LLOADC, LLOADR,                     &
                                         LMATL, LMPC, LMPCADDC, LMPCADDR, LPBAR, LPBEAM, LPBUSH, LPCOMP, LPCOMP_PLIES, LPDAT,      &
                                         LPELAS, LPLOAD, LPMASS, LPROD, LPSHEL, LPSHEAR, LPSOLID, LPUSER1, LPUSERIN, LRFORCE,      &
                                         LRIGEL, LSEQ, LSLOAD, LSNORM, LSPC, LSPC1, LSPCADDC, LSPCADDR, LSUSERIN, LTDAT,           &
                                         MEDAT_CBAR, MEDAT_CBEAM, MEDAT_CBUSH,                                                     &
                                         MEDAT_CELAS1, MEDAT_CELAS2, MEDAT_CELAS3, MEDAT_CELAS4,                                   &
                                         MEDAT_CQUAD, MEDAT_CQUAD8, MEDAT_CROD, MEDAT_CSHEAR, MEDAT_CTRIA, MEDAT_CUSER1,           &
                                         MEDAT0_CUSERIN, MMPC,                                                                     &
                                         MPDAT_PLOAD2, MPDAT_PLOAD4, MEDAT_PLOTEL, MRBE3, MRSPLINE, MTDAT_TEMPRB, MTDAT_TEMPP1,    &
                                         NPBARL, NSPOINT, PROG_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRDSET3, GRDSET7, GRDSET8
      USE PARAMS, ONLY                :  GRIDSEQ

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  FFIELD, FFIELD2
      USE ROD_BAR_BEAM_CARDS, ONLY    :  BD_BAROR0, BD_BEAMOR0, BD_CBAR0
      USE SPRING_BUSH_MASS, ONLY      :  BD_CBUSH0
      USE SOLID_CARDS, ONLY           :  BD_CHEXA0, BD_CPENTA0, BD_CTETRA0
      USE SHELL_COMPOSITE_CARDS, ONLY :  BD_CQUAD0, BD_CQUAD80, BD_CTRIA0, BD_PCOMP0, BD_PCOMP10
      USE USER_ELEMENTS, ONLY         :  BD_CUSERIN0
      USE DEBUG_CARDS, ONLY           :  BD_DEBUG0
      USE GRID_COORDINATES, ONLY      :  BD_GRDSET0, BD_SPOINT0
      USE BULK_DATA_LOADS, ONLY       :  BD_LOAD0, BD_SLOAD0
      USE CONSTRAINT_CARDS, ONLY      :  BD_SPCADD0, BD_MPC0, BD_MPCADD0
      USE PARAM_CARDS, ONLY           :  BD_PARAM0
      USE RIGID_ELEMENTS, ONLY        :  BD_RBE30, BD_RSPLINE0

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'LOADB0'
      CHARACTER( 7*BYTE), PARAMETER   :: END_CARD    = 'ENDDATA'

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! 16 col field card (either CARD1 if small field or CARD1 + CARD2 if
!                                                            a large field card. This is output from subr FFIELD

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD1             ! BD card (a small field card or the 1st half of a large field card)
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD2             ! 2nd half of a large field card
      CHARACTER( 9*BYTE)              :: DECK_NAME   = 'BULK DATA'
      CHARACTER( 1*BYTE)              :: LARGE_FLD_INP     ! If 'Y', card is in large field format

      INTEGER(LONG)                   :: COMMENT_COL       ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: DELTA_LEDAT       ! Delta number of words to add to LEDAT for an element (for elements
!                                                            where the delta cannot be discerned from the BD card name, such as
!                                                            CHEXA which can have 8 or 20 nodes
      INTEGER(LONG)                   :: DELTA_SLOAD       ! Delta number of words to add to LSLOAD for SLOAD's
      INTEGER(LONG)                   :: DELTA_SPOINT      ! Delta number of words to add to LGRID for SPOINT's
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICNT              ! Card count
      INTEGER(LONG)                   :: IERR              ! Error indicatro from subr FFIELD
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a file
      INTEGER(LONG)                   :: ILOAD             ! Number of loads defined on 1 logical B.D. LOAD card
      INTEGER(LONG)                   :: IMPC              ! Number of grid/comp/coeff triplets defined on 1 logical B.D. MPC  card
      INTEGER(LONG)                   :: IMPCADD           ! Number of MPC set ID's defined on 1 B.D. MPCADD card
      INTEGER(LONG)                   :: IRBE3             ! Number of grid/comp/coeff triplets defined on 1 logical B.D. RBE3 card
      INTEGER(LONG)                   :: IRSPLINE          ! Number of fields that can have a grid or comp num on an RSPLINE entry
      INTEGER(LONG)                   :: ISPCADD           ! Number of SPC set ID's defined on 1 B.D. SPCADD card
      INTEGER(LONG)                   :: IPLIES            ! Number of composite layers on 1 B.D. PCOMP card
      INTEGER(LONG)                   :: NG_USERIN         ! Number of grids found on USERIN elems (not incl SPOINT's)
      INTEGER(LONG)                   :: NS_USERIN         ! Number of SPOINT's found on USERIN elems




! **********************************************************************************************************************************
      CARD(1:) = ' '

! Set MRBE3, MRSPLINE so that we will make sure that they will be at least 1

      MRBE3    = 1
      MRSPLINE = 1

! Process Bulk Data cards in a loop that runs until either an ENDDATA card is found or when an error or EOF/EOR occurs

      ICNT = 0
!xx   REWIND (IN1)
      DO

         ICNT = ICNT + 1
         CALL READ_BDF_LINE(IN1, IOCHK, CARD1)
         CARD(1:)  = CARD1(1:)

         ! Quit if EOF/EOR occurs.
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         ! Remove any comments within the CARD1 by deleting everything
         ! from $ on (after col 1)
         COMMENT_COL = 1
         DO I=2,BD_ENTRY_LEN
            IF (CARD1(I:I) == '$') THEN
               COMMENT_COL = I
               EXIT
            ENDIF
         ENDDO

         IF (COMMENT_COL > 1) THEN
            CARD1(COMMENT_COL:) = ' '
         ENDIF

         ! Determine if the card is large or small format
         LARGE_FLD_INP = 'N'
         DO I=1,8
            IF (CARD1(I:I) == '*') THEN
               LARGE_FLD_INP = 'Y'
            ENDIF
         ENDDO

         ! FFIELD converts free-field card to fixed field and left justifies
         ! data in fields 2-9 and outputs a 10 field, 16 col/field CARD1
          IF ((CARD1(1:1) /= '$')  .AND. (CARD1(1:) /= ' ')) THEN

            IF (LARGE_FLD_INP == 'N') THEN

               CALL FFIELD ( CARD1, IERR )
               CARD(1:) = CARD1(1:)

            ELSE

               ! Read 2nd physical entry for a large field parent B.D. entry
               CALL READ_BDF_LINE(IN1, IOCHK, CARD2)

               IF (IOCHK < 0) THEN
                  WRITE(ERR,1011) END_CARD
                  WRITE(F06,1011) END_CARD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               IF (IOCHK > 0) THEN
                  WRITE(ERR,1010) DECK_NAME
                  WRITE(F06,1010) DECK_NAME
                  WRITE(F06,'(A)') CARD2
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ENDIF

               COMMENT_COL = 1
               DO I=2,BD_ENTRY_LEN
                  IF (CARD2(I:I) == '$') THEN
                     COMMENT_COL = I
                     EXIT
                  ENDIF
               ENDDO

               IF (COMMENT_COL > 1) THEN
                  CARD2(COMMENT_COL:) = ' '
               ENDIF

!              IF (CARD2(1:8) /= CARD1(73:80)) THEN
!                 BACKSPACE(IN1)
!                 CARD2(1:) = ' '
!                 CARD2(1:8) = CARD1(73:80)
!              ENDIF

               IF      (CARD2( 1: 8) == CARD1(73:80)) THEN
                  CONTINUE     !ICONT = 1
               ELSE IF ((CARD2( 1: 8) == '*       ') .AND. (CARD1(73:80) == '        ')) THEN
                  CONTINUE     !ICONT = 1
               ELSE IF ((CARD2( 1: 8) == '        ') .AND. (CARD1(73:80) == '*       ')) THEN
                  CONTINUE     !ICONT = 1
               ELSE
                  BACKSPACE(IN1)
                  CARD2(1:) = ' '
                  CARD2(1:8) = CARD1(73:80)
               ENDIF

               CALL FFIELD2 ( CARD1, CARD2, CARD, IERR )

            ENDIF

            IF (IERR /= 0) THEN
               CYCLE
            ENDIF

         ENDIF

          ! No errors, so process Bulk Data card. No need to check for
          ! imbedded blanks found when FFIELD was run - this will be
          ! checked when LOADB reads the bulk data
          IF      (CARD(1:5) == 'BAROR'   )  THEN
            CALL BD_BAROR0 ( CARD )

         ELSE IF (CARD(1:6) == 'BEAMOR'  )  THEN
            CALL BD_BEAMOR0( CARD )

         ELSE IF ((CARD(1:4) == 'CBAR'    ) .OR. (CARD(1:5) == 'CBEAM'   ))  THEN
            LELE  = LELE + 1
            IF (CARD(1:4) == 'CBAR'    ) THEN
               LEDAT = LEDAT + MEDAT_CBAR
            ELSE
               LEDAT = LEDAT + MEDAT_CBEAM
            ENDIF
            CALL BD_CBAR0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'CBUSH'   )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CBUSH
            CALL BD_CBUSH0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CELAS1'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CELAS1

         ELSE IF (CARD(1:6) == 'CELAS2'  )  THEN
            LELE   = LELE + 1
            LEDAT  = LEDAT + MEDAT_CELAS2
            LPELAS = LPELAS + 1                            ! CELAS2 has props on conn entry and we create a PELAS for them

         ELSE IF (CARD(1:6) == 'CELAS3'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CELAS3

         ELSE IF (CARD(1:6) == 'CELAS4'  )  THEN
            LELE   = LELE + 1
            LEDAT  = LEDAT + MEDAT_CELAS4
            LPELAS = LPELAS + 1                            ! CELAS4 has props on conn entry and we create a PELAS for them

         ELSE IF (CARD(1:5) == 'CHEXA'   ) THEN
            LELE  = LELE + 1
            CALL BD_CHEXA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CMASS1'  )  THEN
            LCMASS = LCMASS + 1

         ELSE IF (CARD(1:6) == 'CMASS2'  )  THEN
            LCMASS = LCMASS + 1
            LPMASS = LPMASS + 1

         ELSE IF (CARD(1:6) == 'CMASS3'  )  THEN
            LCMASS = LCMASS + 1

         ELSE IF (CARD(1:6) == 'CMASS4'  )  THEN
            LCMASS = LCMASS + 1
            LPMASS = LPMASS + 1

         ELSE IF (CARD(1:5) == 'CONM2'   )  THEN
            LCONM2 = LCONM2 + 1

         ELSE IF (CARD(1:6) == 'CONROD'  )  THEN
            LELE  = LELE  + 1
            LPROD = LPROD + 1
            LEDAT = LEDAT + MEDAT_CROD

         ELSE IF((CARD(1:6) == 'CORD1C'  ) .OR. (CARD(1:6) == 'CORD1R'  ) .OR. (CARD(1:6) == 'CORD1S'  )) THEN
            LCORD = LCORD + 1
            IF (CARD(41:48) /= '        ') LCORD = LCORD + 1

         ELSE IF((CARD(1:6) == 'CORD2C'  ) .OR. (CARD(1:6) == 'CORD2R'  ) .OR. (CARD(1:6) == 'CORD2S'  )) THEN
            LCORD = LCORD + 1

         ELSE IF (CARD(1:6) == 'CPENTA'  ) THEN
            LELE  = LELE + 1
            CALL BD_CPENTA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CQUAD4'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CQUAD
            CALL BD_CQUAD0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CQUAD8'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CQUAD8
            CALL BD_CQUAD80 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:4) == 'CROD'    )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CROD

         ELSE IF (CARD(1:6) == 'CSHEAR'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CSHEAR

         ELSE IF (CARD(1:6) == 'CTETRA'  ) THEN
            LELE  = LELE + 1
            CALL BD_CTETRA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CTRIA3'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CTRIA
            CALL BD_CTRIA0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CUSER1'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CUSER1

         ELSE IF (CARD(1:7) == 'CUSERIN' )  THEN
            LELE  = LELE + 1
            CALL BD_CUSERIN0 ( CARD, NG_USERIN, NS_USERIN )
            IF (NG_USERIN > LGUSERIN) THEN
               LGUSERIN = NG_USERIN
            ENDIF
            IF (NS_USERIN > LSUSERIN) THEN
               LSUSERIN = NS_USERIN
            ENDIF
                                                           ! LEDAT has "+ 1" term since last record is NUM_BDY_DOF not in MEDAT0
            LEDAT = LEDAT + MEDAT0_CUSERIN + 2*NG_USERIN + NS_USERIN +  1

         ELSE IF (CARD(1:5) == 'DEBUG'   )  THEN
            CALL BD_DEBUG0 ( CARD )

         ELSE IF((CARD(1:5) == 'FORCE'   ) .OR. (CARD(1:6) == 'MOMENT'  )) THEN
            LFORCE = LFORCE +1

         ELSE IF (CARD(1:4) == 'GRAV'    )  THEN
            LGRAV = LGRAV + 1

         ELSE IF (CARD(1:6) == 'GRDSET'  )  THEN
            CALL BD_GRDSET0 ( CARD )

         ELSE IF (CARD(1:4) == 'GRID'    )  THEN
            LGRID = LGRID + 1
            LDOFG = LDOFG + 6

         ELSE IF (CARD(1:4) == 'LOAD'    )  THEN
            LLOADR = LLOADR + 1
            CALL BD_LOAD0 ( CARD, LARGE_FLD_INP, ILOAD )
            IF (ILOAD > LLOADC) THEN
               LLOADC = ILOAD
            ENDIF

         ELSE IF((CARD(1:4) == 'MAT1'    ) .OR. (CARD(1:4) == 'MAT2'    )  .OR.                                                    &
                 (CARD(1:4) == 'MAT8'    ) .OR. (CARD(1:4) == 'MAT9'    )) THEN
            LMATL = LMATL + 1

         ELSE IF((CARD(1:4) == 'MPC '    ) .OR. (CARD(1:4) == 'MPC*'    )) THEN
            LMPC = LMPC + 1
            CALL BD_MPC0 ( CARD, LARGE_FLD_INP, IMPC )
            IF (IMPC > MMPC) THEN
               MMPC = IMPC
            ENDIF

         ELSE IF (CARD(1:6) == 'MPCADD'  )  THEN
            LMPCADDR = LMPCADDR + 1
            CALL BD_MPCADD0 ( CARD, LARGE_FLD_INP, IMPCADD )
            IF (IMPCADD > LMPCADDC) THEN
               LMPCADDC = IMPCADD
            ENDIF

         ELSE IF (CARD(1:6) == 'PARAM '  )  THEN
            CALL BD_PARAM0 ( CARD )

         ELSE IF((CARD(1:5) == 'PBAR '   ) .OR. (CARD(1:5) == 'PBAR*'   ))  THEN
            LPBAR = LPBAR + 1

         ELSE IF (CARD(1:5) == 'PBARL'   )  THEN
            LPBAR  = LPBAR  + 1
            NPBARL = NPBARL + 1

         ELSE IF (CARD(1:5) == 'PBEAM'   )  THEN
            LPBEAM = LPBEAM + 1

         ELSE IF (CARD(1:5) == 'PBUSH'   )  THEN
            LPBUSH = LPBUSH + 1

         ELSE IF (CARD(1:5) == 'PCOMP'   )  THEN
            LPCOMP = LPCOMP + 1
            CALL BD_PCOMP0 ( CARD, LARGE_FLD_INP, IPLIES )
            IF (IPLIES > LPCOMP_PLIES) THEN
               LPCOMP_PLIES = IPLIES
            ENDIF

         ELSE IF (CARD(1:6) == 'PCOMP1'  )  THEN
            LPCOMP = LPCOMP + 1
            CALL BD_PCOMP10 ( CARD, LARGE_FLD_INP, IPLIES )
            IF (IPLIES > LPCOMP_PLIES) THEN
               LPCOMP_PLIES = IPLIES
            ENDIF

         ELSE IF (CARD(1:5) == 'PELAS'   )  THEN
            LPELAS = LPELAS + 1

         ELSE IF (CARD(1:6) == 'PLOAD2'  )  THEN
            LPDAT  = LPDAT  + MPDAT_PLOAD2
            LPLOAD = LPLOAD + 1

         ELSE IF (CARD(1:6) == 'PLOAD4'  )  THEN
            LPDAT  = LPDAT  + MPDAT_PLOAD4
            LPLOAD = LPLOAD + 1

         ELSE IF (CARD(1:6) == 'PLOTEL'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_PLOTEL

         ELSE IF (CARD(1:5) == 'PMASS'   )  THEN
            LPMASS = LPMASS + 4

         ELSE IF (CARD(1:4) == 'PROD'    )  THEN
            LPROD = LPROD + 1

         ELSE IF (CARD(1:6) == 'PSHEAR'  )  THEN
            LPSHEAR = LPSHEAR + 1

         ELSE IF (CARD(1:6) == 'PSHELL'  )  THEN
            LPSHEL = LPSHEL + 1

         ELSE IF (CARD(1:6) == 'PSOLID'  )  THEN
            LPSOLID = LPSOLID + 1

         ELSE IF (CARD(1:6) == 'PUSER1'  )  THEN
            LPUSER1 = LPUSER1 + 1

         ELSE IF (CARD(1:7) == 'PUSERIN' )  THEN
            LPUSERIN = LPUSERIN + 1

         ELSE IF (CARD(1:4) == 'RBAR'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE1'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE2'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE3'    )  THEN
            LRIGEL  = LRIGEL + 1
            CALL BD_RBE30 ( CARD, LARGE_FLD_INP, IRBE3 )
            IF (IRBE3 > MRBE3) THEN
               MRBE3 = IRBE3
            ENDIF

         ELSE IF (CARD(1:6) == 'RFORCE'  )  THEN
            LRFORCE = LRFORCE + 1

         ELSE IF (CARD(1:7) == 'RSPLINE' )  THEN
            CALL BD_RSPLINE0 ( CARD, LARGE_FLD_INP, IRSPLINE )
            LRIGEL = LRIGEL + 1
            IF (IRSPLINE > MRSPLINE) THEN
               MRSPLINE = IRSPLINE
            ENDIF

         ELSE IF (CARD(1:5) == 'SLOAD'   )  THEN
            CALL BD_SLOAD0 ( CARD, DELTA_SLOAD )
            LSLOAD = LSLOAD + DELTA_SLOAD

         ELSE IF (CARD(1:5) == 'SEQGP'   )  THEN
            LSEQ = LSEQ + 4                                ! Conservative estimate. There can only be 4 entries per card

         ELSE IF (CARD(1:5) == 'SNORM'   )  THEN
            LSNORM = LSNORM + 1

         ELSE IF((CARD(1:4) == 'SPC '    ) .OR. (CARD(1:4) == 'SPC*'    )) THEN
            LSPC = LSPC + 1

         ELSE IF (CARD(1:4) == 'SPC1'    )  THEN
            LSPC1 = LSPC1 + 1

         ELSE IF (CARD(1:6) == 'SPCADD'  )  THEN
            LSPCADDR = LSPCADDR + 1
            CALL BD_SPCADD0 ( CARD, LARGE_FLD_INP, ISPCADD )
            IF (ISPCADD > LSPCADDC) THEN
               LSPCADDC = ISPCADD
            ENDIF

         ELSE IF (CARD(1:6) == 'SPOINT'  )  THEN
            CALL BD_SPOINT0 ( CARD, LARGE_FLD_INP, DELTA_SPOINT )
            NSPOINT = NSPOINT + DELTA_SPOINT               ! DELTA_SPOINT = number of SPOINTS defined on this SPOINT Bulk Data entry
            LGRID   = LGRID   + DELTA_SPOINT               ! Each SPOINT counts as 1 in the number of grids
            LDOFG   = LDOFG   + DELTA_SPOINT

         ELSE IF (CARD(1:6) == 'TEMPRB'  )  THEN
            LTDAT = LTDAT + MTDAT_TEMPRB

         ELSE IF (CARD(1:6) == 'TEMPP1'  )  THEN
            LTDAT = LTDAT + MTDAT_TEMPP1

         ELSE IF (CARD(1:7) == 'ENDDATA' )  THEN
            EXIT

         ENDIF

      ENDDO

! ! Reset FATAL_ERR to 0. Errors that incremented FATAL_ERR in LOADBO (and routines it calls) will be recorded and
! ! reported when LOADB runs.
! !
! !   FATAL_ERR = 0


      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

! **********************************************************************************************************************************

      END SUBROUTINE LOADB0


      SUBROUTINE LOADB_RESTART

      ! LOADB_RESTART reads in some entries in the Bulk Data deck
      ! (e.g., DEBUG, PARAM) for a RESTART run
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN, JF, PROG_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  FFIELD, FFIELD2
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKCARD, MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD
      USE PARAM_CARDS, ONLY           :  BD_PARAM

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: NUM_PARMS = 25      ! Number of PARAM entries allowed in RESTART
      INTEGER(LONG), PARAMETER        :: NUM_DEB   = 28      ! Number of DEBUG entries allowed in RESTART

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'LOADB_RESTART'

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD                ! 16 col field card (either CARD1 if small field or CARD1 + CARD2 if
!                                                              a large field card. This is output from subr FFIELD

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD1               ! BD card (a small field card or the 1st half of a large field card)
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD2               ! 2nd half of a large field card
      CHARACTER( 9*BYTE)              :: DECK_NAME   = 'BULK DATA'
      CHARACTER( 7*BYTE), PARAMETER   :: END_CARD    = 'ENDDATA'
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of 8 characters making up CARD
      CHARACTER( 1*BYTE)              :: LARGE_FLD_INP       ! If 'Y', card is in large field format
      CHARACTER( 8*BYTE)              :: PARM_NAME(NUM_PARMS)! Names of PARAM entries allowed in RESTART

      INTEGER(LONG)                   :: COMMENT_COL         ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: DEB_NUM(NUM_DEB)    ! Allowable DEBUG indices in restart
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: INDEX               ! Index
      INTEGER(LONG)                   :: INT_VAL             ! Integer value read fron a card field
      INTEGER(LONG)                   :: IERR                ! Error indicator from subr FFIELD
      INTEGER(LONG)                   :: IOCHK               ! IOSTAT error number when reading Bulk Data cards from unit IN1




! **********************************************************************************************************************************
      DO I=1,NUM_PARMS
         PARM_NAME(I) = '        '
      ENDDO

      PARM_NAME( 1) = 'AUTOSPC '
      PARM_NAME( 2) = 'ELFORCEN'
      PARM_NAME( 3) = 'EQCHECK '
      PARM_NAME( 4) = 'GRDPNT  '
      PARM_NAME( 5) = 'POST    '
      PARM_NAME( 6) = 'PRTBASIC'
      PARM_NAME( 7) = 'PRTCORD '
      PARM_NAME( 8) = 'PRTDISP'
      PARM_NAME( 9) = 'PRTDOF  '
      PARM_NAME(10) = 'PRTFOR  '
      PARM_NAME(11) = 'PRTGMN  '
      PARM_NAME(12) = 'PRTGOA  '
      PARM_NAME(13) = 'PRTjMN  '
      PARM_NAME(14) = 'PRTMASS '
      PARM_NAME(15) = 'PRTQSYS '
      PARM_NAME(16) = 'PRTRMG  '
      PARM_NAME(17) = 'PRTSCP  '
      PARM_NAME(18) = 'PRTTSET '
      PARM_NAME(19) = 'PRTSTIFF'
      PARM_NAME(20) = 'PRTSTIFD'
      PARM_NAME(21) = 'PRTUO0  '
      PARM_NAME(22) = 'PRTYS   '
      PARM_NAME(23) = 'RELINK3 '
      PARM_NAME(24) = 'SORT_MAX'
      PARM_NAME(25) = 'TINY    '

      DEB_NUM( 1) =   1
      DEB_NUM( 2) =   2
      DEB_NUM( 3) =   3
      DEB_NUM( 4) =   5
      DEB_NUM( 5) =   6
      DEB_NUM( 6) =   9
      DEB_NUM( 7) =  11
      DEB_NUM( 8) =  13
      DEB_NUM( 9) =  15
      DEB_NUM(10) =  18
      DEB_NUM(11) =  19
      DEB_NUM(12) =  21
      DEB_NUM(13) =  22
      DEB_NUM(14) =  81
      DEB_NUM(15) =  82
      DEB_NUM(16) =  83
      DEB_NUM(17) =  84
      DEB_NUM(18) =  85
      DEB_NUM(19) =  86
      DEB_NUM(20) =  91
      DEB_NUM(21) =  92
      DEB_NUM(22) = 101
      DEB_NUM(23) = 195
      DEB_NUM(24) = 196
      DEB_NUM(25) = 197
      DEB_NUM(26) = 198
      DEB_NUM(27) = 199
      DEB_NUM(28) = 200

      WRITE(F06,100)

      ! Process Bulk Data cards in a large loop that runs until either an
      ! ENDDATA card is found or when an error or EOF/EOR occurs
      DO
         CALL READ_BDF_LINE(IN1, IOCHK, CARD1)
         CARD(1:) = CARD1(1:)

         ! Quit if EOF/EOR occurs.
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         ! Remove any comments within the CARD1 by deleting everything
         ! from $ on (after col 1)
         COMMENT_COL = 1
         DO I=2,BD_ENTRY_LEN
            IF (CARD1(I:I) == '$') THEN
               COMMENT_COL = I
               EXIT
            ENDIF
         ENDDO

         IF (COMMENT_COL > 1) THEN
            CARD1(COMMENT_COL:) = ' '
         ENDIF

         ! Determine if the card is large or small format
         LARGE_FLD_INP = 'N'
         DO I=1,8
            IF (CARD1(I:I) == '*') THEN
               LARGE_FLD_INP = 'Y'
            ENDIF
         ENDDO

         ! FFIELD converts free-field card to fixed field and left justifies
         ! data in fields 2-9 and outputs a 10 field, 16 col/field CARD1
         IF ((CARD1(1:1) /= '$')  .AND. (CARD1(1:) /= ' ')) THEN

            IF (LARGE_FLD_INP == 'N') THEN

               CALL FFIELD ( CARD1, IERR )
               CARD(1:) = CARD1(1:)

            ELSE
               ! Read 2nd physical entry for a large field parent B.D. entry
               CALL READ_BDF_LINE(IN1, IOCHK, CARD2)

               IF (IOCHK < 0) THEN
                  WRITE(ERR,1011) END_CARD
                  WRITE(F06,1011) END_CARD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               IF (IOCHK > 0) THEN
                  WRITE(ERR,1010) DECK_NAME
                  WRITE(F06,1010) DECK_NAME
                  WRITE(F06,'(A)') CARD2
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ENDIF

               IF (ECHO /= 'NONE  ') THEN
                  WRITE(F06,101) CARD2
               ENDIF

               COMMENT_COL = 1
               DO I=2,BD_ENTRY_LEN
                  IF (CARD2(I:I) == '$') THEN
                     COMMENT_COL = I
                     EXIT
                  ENDIF
               ENDDO

               IF (COMMENT_COL > 1) THEN
                  CARD2(COMMENT_COL:) = ' '
               ENDIF

               CALL FFIELD2 ( CARD1, CARD2, CARD, IERR )

            ENDIF

            IF (IERR /= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1003)
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,101) CARD
               ENDIF
               WRITE(F06,1003)
               CYCLE
            ENDIF

         ENDIF

         ! Process Bulk Data card
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

         IF      ((JCARD(1)(1:6) == 'DEBUG ') .OR. (JCARD(1)(1:6) == 'DEBUG*'))  THEN
            READ(JCARD(2),'(I8)') INDEX
            DO I=1,NUM_DEB
               CALL I4FLD ( JCARD(2), JF(2), INDEX )
               IF (INDEX == DEB_NUM(I)) THEN
                  CALL I4FLD ( JCARD(3), JF(3), INT_VAL )
                  DEBUG(INDEX) = INT_VAL
                  WRITE(F06,101) CARD
                  EXIT
               ENDIF
            ENDDO

         ELSE IF ((JCARD(1)(1:6) == 'PARAM ') .OR. (JCARD(1)(1:6) == 'PARAM*'))  THEN
            IF (JCARD(2) == 'AUTOSPC ') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,103)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,103)
               ENDIF
               JCARD(4)(1:) = ' '
               JCARD(5)(1:) = ' '
               JCARD(6)(1:) = ' '
               CALL MKCARD ( JCARD, CARD )
            ENDIF
            DO I=1,NUM_PARMS
               IF (JCARD(2) == PARM_NAME(I)) THEN
                  WRITE(F06,101) CARD
                  CALL BD_PARAM  ( CARD )
                  EXIT
               ENDIF
            ENDDO

         ELSE IF (CARD(1:8) == 'ENDDATA ')  THEN
            WRITE(F06,101) CARD
            EXIT

         ELSE IF ((CARD(1:1) == '$') .OR. (CARD(1:BD_ENTRY_LEN) == ' ')) THEN
            CYCLE

         ELSE                                              ! CARD not processed by MYSTRAN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,9993) PROG_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) CARD
               WRITE(F06,9993) PROG_NAME
            ENDIF

         ENDIF

      ENDDO



      RETURN

! **********************************************************************************************************************************
  100 FORMAT('$'                                                                                                                ,/,&
             '$ Following are the only Bulk Data entries that were processed for this restart. All other model data will be the'  ,&
              ' same as from the'                                                                                               ,/,&
             '$ original run except that any PARAM or DEBUG Bulk Data entries not in this restart deck will revert to their'      ,&
              ' default values.'                                                                                                ,/,&
             '$ Some PARAM and DEBUG entries may not be processed in this restart. See the MYSTRAN documentaton discussion on'    ,&
              ' restarts.'                                                                                                      ,/,&
             '$')


  101 FORMAT(A)

  103 FORMAT(' *WARNING    : THE ONLY FIELD FROM THE PARAM AUTOSPC THAT WILL BE USED IN THIS RESTART IS FIELD 7 FOR SPC FORCE INFO')

 1003 FORMAT(' *ERROR  1003: ALL FIELDS ON THE ABOVE ENTRY MUST BE NO LONGER THAN 8 CHARACTERS')

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

 9993 FORMAT(' *WARNING    : PRIOR ENTRY NOT PROCESSED BY ',A)

! **********************************************************************************************************************************

      END SUBROUTINE LOADB_RESTART


      SUBROUTINE LOADC

      ! LOADC reads in the CASE CONTROL DECK
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUGOUT, ERR, F06, IN1, WRT_ERR
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, ENFORCED, FATAL_ERR, WARN_ERR, NSUB, NTSUB,                                 &
                                         NUM_BUCKLING_SUBS, PROG_NAME, RESTART, SOL_NAME, CC_CMD_DESCRIBERS
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, CC_EIGR_SID_SUB, CC_EIGR_SID_DECK, CC_STATSUB_DECK, CC_STATSUB_SUB,          &
                                         IS_BUCKLING_SUBCASE, IS_MODES_SUBCASE,                                                    &
                                         MEFFMASS_CALC, MPCSET, MPCSETS, MPFACTOR_CALC, SCNUM, SPCSET, SPCSETS, SUBLOD,            &
                                         SC_STRE, SC_STRN, SC_ELFE, SC_ELFN
      USE MODEL_STUF, ONLY            :  EIG_PARAMS,                                                                               &
                                         EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_LANCZOS_NEV_DELT, EIG_METH,         &
                                         EIG_MSGLVL, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, EIG_SID,   &
                                         EIG_SIGMA, EIG_VECS
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC, FORC_LOC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, REPLACE_TABS_W_BLANKS
      USE CASE_CONTROL_OUTPUTS, ONLY  :  CC_ACCE, CC_DISP, CC_ELDA, CC_ELFO, CC_ENFO, CC_GPFO, CC_MPCF, CC_OLOA, CC_SPCF, CC_STRE, CC_STRN
      USE CASE_CONTROL_METADATA, ONLY :  CC_ECHO, CC_LABE, CC_SUBT, CC_TITL
      USE CASE_CONTROL_SELECTORS, ONLY:  CC_LOAD, CC_METH, CC_MPC, CC_NLPARM, CC_SPC, CC_STATSUB, CC_TEMP
      USE CASE_CONTROL_SETS, ONLY     :  CC_SET, CC_SUBC
      USE TEXT_FIELD_UTILS, ONLY      :  TO_UPPER

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: DOLLAR_WARN       ! Indicator of whether there was a $ sign in col 1
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD              ! Case Control card
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'CASE CONTROL'
      CHARACTER(10*BYTE), PARAMETER   :: END_CARD    = 'BEGIN BULK'

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: QUAD4_LOC


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

         ELSE IF (CARD1(1:7) == 'STATSUB') THEN
            CALL CC_STATSUB ( CARD1 )

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

                                                           ! Assign the same CENTER/CORNER location to all of FORCE,
                                                           ! STRESS, and STRAIN outputs according to the priority rules.
                                                           ! Because NONE is stored as 0 in SC_STRE, etc., this treats
                                                           ! STRESS(CORNER) = NONE as if STRESS wasn't defined at
                                                           ! all and defaults to CENTER.
      IF       (SC_STRE(1) /= 0) THEN
         QUAD4_LOC = STRE_LOC
      ELSE IF  (SC_STRN(1) /= 0) THEN
         QUAD4_LOC = STRN_LOC
      ELSE IF  (SC_ELFE(1) /= 0 .OR. SC_ELFN(1) /= 0) THEN
         QUAD4_LOC = FORC_LOC
      ELSE
         QUAD4_LOC = 'CENTER'
      ENDIF
      STRE_LOC = QUAD4_LOC
      STRN_LOC = QUAD4_LOC
      FORC_LOC = QUAD4_LOC
                                                           ! From here on, STRE_LOC = STRN_LOC = FORC_LOC.

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

      ! Propagate the deck-level METHOD default down to every subcase that does not have its own METHOD.
      ! This must happen *after* the NSUB==0 -> NSUB=1 fixup above, so that an "implicit" single subcase still picks
      ! up the METHOD declared at the deck level. The deck default is also recorded as the SID for any subcase that
      ! omitted METHOD. For SOL types that do not consume eigenvalue extraction the arrays simply stay zero.
      IF (ALLOCATED(CC_EIGR_SID_SUB)) THEN
         DO I=1,NSUB
            IF ((CC_EIGR_SID_SUB(I) == 0) .AND. (CC_EIGR_SID_DECK /= 0)) THEN
               CC_EIGR_SID_SUB(I)  = CC_EIGR_SID_DECK
               IS_MODES_SUBCASE(I) = 'Y'
            ENDIF
         ENDDO
         ! Keep the legacy scalar in sync: prefer the deck default; otherwise the first non-zero per-subcase SID.
         IF (CC_EIGR_SID_DECK /= 0) THEN
            CC_EIGR_SID = CC_EIGR_SID_DECK
         ELSE
            DO I=1,NSUB
               IF (CC_EIGR_SID_SUB(I) /= 0) THEN
                  CC_EIGR_SID = CC_EIGR_SID_SUB(I)
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         ! Final fallback: for any subcase whose EIG_PARAMS slot is still empty (could happen when NSUB was 0 at the
         ! time BD_EIGR/BD_EIGRL ran -- e.g. a deck with a deck-level METHOD and no SUBCASE cards), copy the values
         ! from the legacy EIG_* scalars. After BD_EIGR/BD_EIGRL has WRITTEN L1M for the canonical (scalar-match) card,
         ! those scalars hold that card's data, which is exactly what an inherited subcase should use.
         IF (ALLOCATED(EIG_PARAMS)) THEN
            DO I=1,NSUB
               IF ((CC_EIGR_SID_SUB(I) /= 0) .AND. (EIG_PARAMS(I)%SID == 0)) THEN
                  EIG_PARAMS(I)%METHOD            = EIG_METH
                  EIG_PARAMS(I)%NORM              = EIG_NORM
                  EIG_PARAMS(I)%LAP_MAT_TYPE      = EIG_LAP_MAT_TYPE
                  EIG_PARAMS(I)%VECS              = EIG_VECS
                  EIG_PARAMS(I)%SID               = CC_EIGR_SID_SUB(I)
                  EIG_PARAMS(I)%N1                = EIG_N1
                  EIG_PARAMS(I)%N2                = EIG_N2
                  EIG_PARAMS(I)%COMP              = EIG_COMP
                  EIG_PARAMS(I)%GRID              = EIG_GRID
                  EIG_PARAMS(I)%LANCZOS_NEV_DELT  = EIG_LANCZOS_NEV_DELT
                  EIG_PARAMS(I)%MODE              = EIG_MODE
                  EIG_PARAMS(I)%MSGLVL            = EIG_MSGLVL
                  EIG_PARAMS(I)%NCVFACL           = EIG_NCVFACL
                  EIG_PARAMS(I)%CRIT              = EIG_CRIT
                  EIG_PARAMS(I)%FRQ1              = EIG_FRQ1
                  EIG_PARAMS(I)%FRQ2              = EIG_FRQ2
                  EIG_PARAMS(I)%SIGMA             = EIG_SIGMA
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      ! If SOL is modes or CB or buckilng, then a METH card should have been found in Case Control
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN_CB_MODEL') .OR. (SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF (CC_EIGR_SID == 0) THEN
             WRITE(ERR,1004)
             WRITE(F06,1004)
             FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

     ! For SOL 105 buckling, resolve STATSUB(PRELOAD) references for each buckling subcase. A buckling subcase is one
     ! that has a resolved METHOD (and therefore IS_MODES_SUBCASE(I)=='Y'). The remaining subcases supply linear-static
     ! preloads. STATSUB(PRELOAD)=n names the external SUBCASE id of the static subcase whose displacement field is
     ! used to assemble KGGD for the buckling eigenproblem. If no STATSUB is given anywhere, fall back to the legacy
     ! behavior: the first non-buckling subcase that carries a LOAD or TEMP request acts as the preload source.
      IF (SOL_NAME == 'BUCKLING') THEN

         ! Tag every modes-subcase as a buckling-subcase, and count them.
         NUM_BUCKLING_SUBS = 0
         IF (ALLOCATED(IS_BUCKLING_SUBCASE) .AND. ALLOCATED(IS_MODES_SUBCASE)) THEN
            DO I = 1, NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  IS_BUCKLING_SUBCASE(I) = 'Y'
                  NUM_BUCKLING_SUBS = NUM_BUCKLING_SUBS + 1
               ENDIF
            ENDDO
         ENDIF

         IF (NUM_BUCKLING_SUBS == 0) THEN                  ! no buckling subcase -> nothing to solve
            WRITE(ERR,1101)
            WRITE(F06,1101)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Confirm at least one non-buckling subcase carries a load (mechanical or thermal). Without one there is no
         ! linear-static preload to drive KGGD.
         BLOCK
            LOGICAL :: HAS_STATIC_LOAD
            INTEGER(LONG) :: STATSUB_REQ, RESOLVED_IDX, K
            HAS_STATIC_LOAD = .FALSE.
            DO I = 1, NSUB
               IF ((IS_BUCKLING_SUBCASE(I) == 'N') .AND. ((SUBLOD(I,1) /= 0) .OR. (SUBLOD(I,2) /= 0))) THEN
                  HAS_STATIC_LOAD = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT. HAS_STATIC_LOAD) THEN
               WRITE(ERR,1101)
               WRITE(F06,1101)
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

            ! Resolve STATSUB(PRELOAD) for each buckling subcase. Per-subcase value beats deck-level which beats the
            ! legacy fallback (first static subcase that carries a load).
            DO I = 1, NSUB
               IF (IS_BUCKLING_SUBCASE(I) /= 'Y') CYCLE
               STATSUB_REQ = 0
               IF (ALLOCATED(CC_STATSUB_SUB)) STATSUB_REQ = CC_STATSUB_SUB(I)
               IF (STATSUB_REQ == 0) STATSUB_REQ = CC_STATSUB_DECK

               RESOLVED_IDX = 0
               IF (STATSUB_REQ == 0) THEN
                  ! Legacy fallback: first static subcase with LOAD or TEMP.
                  DO K = 1, NSUB
                     IF ((IS_BUCKLING_SUBCASE(K) == 'N') .AND. ((SUBLOD(K,1) /= 0) .OR. (SUBLOD(K,2) /= 0))) THEN
                        RESOLVED_IDX = K
                        EXIT
                     ENDIF
                  ENDDO
                  IF (RESOLVED_IDX == 0) THEN
                     WRITE(ERR,1103) I
                     WRITE(F06,1103) I
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ELSE
                  ! Look up external subcase id (STATSUB_REQ) in SCNUM(:) to find its internal index.
                  DO K = 1, NSUB
                     IF (SCNUM(K) == STATSUB_REQ) THEN
                        RESOLVED_IDX = K
                        EXIT
                     ENDIF
                  ENDDO
                  IF (RESOLVED_IDX == 0) THEN
                     WRITE(ERR,1104) STATSUB_REQ, I
                     WRITE(F06,1104) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                  ELSE IF (IS_BUCKLING_SUBCASE(RESOLVED_IDX) == 'Y') THEN
                     WRITE(ERR,1105) STATSUB_REQ, I
                     WRITE(F06,1105) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                     RESOLVED_IDX = 0
                  ELSE IF ((SUBLOD(RESOLVED_IDX,1) == 0) .AND. (SUBLOD(RESOLVED_IDX,2) == 0)) THEN
                     WRITE(ERR,1106) STATSUB_REQ, I
                     WRITE(F06,1106) STATSUB_REQ, I
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF

               IF (ALLOCATED(EIG_PARAMS) .AND. (RESOLVED_IDX > 0)) THEN
                  EIG_PARAMS(I)%STATSUB_REF = RESOLVED_IDX
               ENDIF
            ENDDO
         END BLOCK

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

 1016 FORMAT(' *INFORMATION: ALL SUBCASES WILL USE "',A,'" AS THE LOCATION OF ',A,' OUTPUTS FOR PSHELL QUAD4 ELEMENTS')

 1028 FORMAT(' *ERROR  1028: THERE MUST BE 2 SUBCASES FOR LINEAR BUCKLING ANALYSES BUT NSUB = ',I8)

 1101 FORMAT(' *ERROR  1101: FOR BUCKLING ANALYSES THERE MUST BE AT LEAST ONE SUBCASE WITH A METHOD (BUCKLING EIGENPROBLEM)',     &
                           ' AND AT LEAST ONE OTHER SUBCASE THAT CARRIES A LOAD (AND/OR TEMP) TO ACT AS THE STATSUB PRELOAD.')

 1102 FORMAT(' *ERROR  1102: FOR BUCKLING ANALYSES ANY LOAD, SPS OR MPC IN 2nd SUBCASE MUST BE THE SAME AS THOSE IN 1st SUBCASE')

 1103 FORMAT(' *ERROR  1103: BUCKLING SUBCASE (INTERNAL #',I8,') HAS NO STATSUB AND NO STATIC SUBCASE WITH A LOAD WAS FOUND',     &
                           ' TO ACT AS THE LEGACY DEFAULT PRELOAD.')

 1104 FORMAT(' *ERROR  1104: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES A SUBCASE THAT DOES NOT',  &
                           ' EXIST IN THIS CASE CONTROL DECK.')

 1105 FORMAT(' *ERROR  1105: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES ANOTHER BUCKLING SUBCASE.',&
                           ' STATSUB MUST POINT TO A LINEAR-STATIC SUBCASE.')

 1106 FORMAT(' *ERROR  1106: STATSUB(PRELOAD)=',I8,' ON BUCKLING SUBCASE (INTERNAL #',I8,') REFERENCES A SUBCASE THAT CARRIES',   &
                           ' NEITHER LOAD NOR TEMP.')

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


      SUBROUTINE LOADC0

      ! Preliminary reading of the Case Control to count several data sizes
      ! so that arrays may be allocated prior to the final reading of the
      ! Case Control.
      !
      ! LOADC0 reads in the CASE CONTROL DECK and:
      !   1) Counts the number of subcases and increments LSUB
      !   2) Counts the number of SET cards and calls CC_SET0 to count the
      !      number of characters in SET's to determine LSETLN

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, FATAL_ERR, LSETS, LSUB
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, REPLACE_TABS_W_BLANKS
      USE CASE_CONTROL_SETS, ONLY     :  CC_SET0

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LOADC0'
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD              ! Case Control card
      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'CASE CONTROL'
      CHARACTER(10*BYTE)              :: END_CARD

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: JERR              ! Error count




! **********************************************************************************************************************************
      ! Initialize
      JERR     = 0
      END_CARD = 'CEND      '

      ! Process CASE CONTROL DECK
      JERR = 0
      END_CARD = 'BEGIN BULK'

      DO
         ! Read an input deck card
         CALL READ_BDF_LINE(IN1, IOCHK, CARD)

         IF (IOCHK < 0) THEN                               ! Quit if EOF/EOR occurs during read
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (IOCHK > 0) THEN                               ! If error occurs during read, write message & CYCLE back to read again
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         CALL REPLACE_TABS_W_BLANKS ( CARD )               ! Replace all tab characters with a white space

         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )  ! Shift CARD so it begins in col 1

         ! Check for Case Control cards. Exit loop on 'BEGIN BULK'
         IF      (CARD1(1: 4) == 'SUBC'      ) THEN
            LSUB = LSUB + 1

         ELSE IF (CARD1(1: 3) == 'SET'       ) THEN
            LSETS = LSETS + 1
            CALL CC_SET0 ( CARD1 )

         ELSE IF (CARD1(1:10) == 'BEGIN BULK') THEN
            EXIT
         ENDIF

      ENDDO

      ! If there were no explicit subcases, then default LSUB to 1
      IF (LSUB == 0) THEN
         LSUB = 1
      ENDIF

      IF (JERR > 0) THEN
         WRITE(ERR,10141)
         WRITE(F06,10141)
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

10141 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE LOADC0


      SUBROUTINE LOADE

      ! LOADE reads in the EXEC CONTROL DECK
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN, CHKPNT, FATAL_ERR, WARN_ERR, JCARD_LEN, JF,     &
                                         PROG_NAME, SOL_NAME, RESTART
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      use debug_parameters, only      :  debug

      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES, ALLOW_OU4_MYSTRAN_NAMES,                     &
                                         ALLOW_OU4_OUTPUT_NAMES, OU4_PART_MAT_NAMES, OU4_PART_VEC_NAMES, NUM_OU4_VALID_NAMES


      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, REPLACE_TABS_W_BLANKS
      USE EXECUTIVE_CONTROL_IO, ONLY  :  EC_DEBUG, EC_IN4FIL, EC_OUTPUT4, EC_PARTN
      USE BDF_SET_SYNTAX, ONLY        :  STOKEN
      USE BDF_FIELD_VALIDATION, ONLY  :  CRDERR, I4FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LOADE'
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD              ! Exec Control deck card
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(LEN=JCARD_LEN)        :: CHARFLD           ! Character field used when suvr I4FLD is called
      CHARACTER(12*BYTE)              :: DECK_NAME = 'EXEC CONTROL'
      CHARACTER( 1*BYTE)              :: DOLLAR_WARN       ! Indicator of whether there was a $ sign in col 1
      CHARACTER( 4*BYTE), PARAMETER   :: END_CARD  = 'CEND'
      CHARACTER(LEN=EC_ENTRY_LEN)     :: ERRTOK            ! An error message that may be returned from subr STOKEN
      CHARACTER( 3*BYTE)              :: EXCEPT    = 'OFF' ! An input/output variable for subr STOKEN, called herein


      CHARACTER( 1*BYTE)              :: PRT_OU4_VALID_NAMES! If 'Y', print valid OUTPUT4 matrix names

      CHARACTER( 3*BYTE)              :: THRU      = 'OFF' ! An input/output variable for subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! Character tokens returned from subr STOKEN
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! Description of the TOKEN's returned from subr STOKEN

                                                           ! Proper SOL number
      CHARACTER( 10*BYTE)             :: SOL_NUM_SHOULD_BE = '1, 3 or 31'
      CHARACTER( 1*BYTE)              :: ANY_OU4_NAME_BAD  ! 'Y'/'N' if requested OUTPUT4 matrix name is valid

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: EC_OUTPUT4_ERR = 0! Count of errors when readig OUTPUT4 entries
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR           = 0! Error indicator.
      INTEGER(LONG)                   :: JERR           = 0! Error indicator.
      INTEGER(LONG)                   :: IERROR            ! An error number returned from subr STOKEN
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: ISTART            ! An input/output for subr STOKEN (where a token begins in CARD)
      INTEGER(LONG)                   :: NTOKEN            ! An output from subr STOKEN (how many tokens were read)
      INTEGER(LONG)                   :: SOL_INT           ! Integer value read from an Exec Control SOL entry
      INTEGER(LONG)                   :: TOKLEN            ! Length of character string sent to subr STOKEN (= LEN(CARD))




! **********************************************************************************************************************************
      DOLLAR_WARN = 'N'

      EC_OUTPUT4_ERR      = 0
      ANY_OU4_NAME_BAD    = 'Y'
      PRT_OU4_VALID_NAMES = 'N'

! Initialize arrays used in OUTPUT4 processing

      DO I=1, NUM_OU4_VALID_NAMES
         ACT_OU4_MYSTRAN_NAMES(I)(1:) = ' '
         ACT_OU4_OUTPUT_NAMES(I)(1:)  = ' '
         OU4_PART_MAT_NAMES(I,1)(1:)    = ' '
         OU4_PART_VEC_NAMES(I,1)(1:)    = ' '
         OU4_PART_VEC_NAMES(I,2)(1:)    = ' '
      ENDDO

! Process EXECUTIVE CONTROL DECK

      DO
         CALL READ_BDF_LINE(IN1, IOCHK, CARD)

         IF (IOCHK < 0) THEN                               ! Quit if EOF/EOR occurs during read
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (IOCHK > 0) THEN                               ! Check if error occurs during read.
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         WRITE(F06,101) CARD

         CALL REPLACE_TABS_W_BLANKS ( CARD )               ! Replace all tab characters with a white space

         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )  ! Shift CARD so it begins in col 1

         IF (CARD1(1:1) == '$') THEN
            DO I=IACHAR('A'),IACHAR('Z')
               IF (CARD1(2:2) == ACHAR(I)) THEN
                  DOLLAR_WARN = 'Y'
               ENDIF
            ENDDO

         ELSE IF (CARD1(1:EC_ENTRY_LEN) == ' ') THEN
            CYCLE

         ELSE IF (CARD1(1: 3) == 'APP'       ) THEN
            CONTINUE

         ELSE IF (CARD1(1: 4) == 'CEND'      ) THEN
            EXIT

         ELSE IF (CARD1(1: 6) == 'CHKPNT'    ) THEN
            CHKPNT = 'Y'

         ELSE IF (CARD1(1: 5) == 'DEBUG'     ) THEN
            CALL EC_DEBUG ( CARD1 )

         ELSE IF (CARD1(1: 2) == 'ID'        ) THEN
            CONTINUE

         ELSE IF (CARD1(1: 3) == 'IN4'       ) THEN
            CALL EC_IN4FIL ( CARD1 )

         ELSE IF (CARD1(1: 7) == 'OUTPUT4'   ) THEN
            CALL EC_OUTPUT4 ( CARD1, JERR, ANY_OU4_NAME_BAD )
            EC_OUTPUT4_ERR = EC_OUTPUT4_ERR + JERR
            IF (ANY_OU4_NAME_BAD == 'Y') THEN
               PRT_OU4_VALID_NAMES = 'Y'
            ENDIF

         ELSE IF (CARD1(1: 5) == 'PARTN'     ) THEN
            CALL EC_PARTN ( CARD1, JERR )
            EC_OUTPUT4_ERR = EC_OUTPUT4_ERR + JERR

         ELSE IF (CARD1(1: 7) == 'RESTART'   ) THEN
            RESTART = 'Y'

         ELSE IF (CARD1(1: 3) == 'SOL'       ) THEN

            SOL_NAME(1:LEN(SOL_NAME)) = ' '
            ISTART = 4
            TOKLEN = LEN(CARD1)
            CALL STOKEN ( SUBR_NAME, CARD1, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

            IF (NTOKEN > 1) THEN

               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1027)
               WRITE(F06,1027)
               CYCLE

            ELSE

               IF (TOKTYP(1) == 'UNKNOWN ') THEN           ! Maybe TOKEN(1) is a name (i.e. SOL_NAME)
                                                           ! Code removed that was here to check IERROR = 1. See explanation above
                  IF      ((TOKEN(1) == 'STATICS ') .OR. (TOKEN(1) == 'SESTATIC')) THEN
                     SOL_NAME(1:7)  =  'STATICS'
                     SOL_INT        = 1

                  ELSE IF ((TOKEN(1) == 'MODAL   ') .OR. (TOKEN(1) == 'MODES   ') .OR. (TOKEN(1) == 'NORMAL M').OR.                &
                           (TOKEN(1) == 'SEMODES ')) THEN
                     SOL_NAME(1:5)  =  'MODES'
                     SOL_INT        = 3

                  ELSE IF ((TOKEN(1) == 'DIFFEREN') .OR. (TOKEN(1) == 'DIFF STI')) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:8)  =  'DIFFEREN'
                     SOL_INT        = 4

                  ELSE IF ((TOKEN(1) == 'BUCKLING') .OR. (TOKEN(1) == 'SEBUCKL')) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:8)  =  'BUCKLING'
                     SOL_INT        = 5

                  ELSE IF ((TOKEN(1) == 'NLSTATIC') .OR. (TOKEN(1) == 'GNOLIN  ') .OR. (TOKEN(1) == 'DIFFEREN')) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:8)  =  'NLSTATIC'
                     SOL_INT        = 106

                  ELSE IF (TOKEN(1)(1:3) == 'GEN') THEN
                     CALL STOKEN ( SUBR_NAME, CARD1, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
                     IF   (TOKEN(1) == 'CB      ') THEN ! This is TOKEN(1) instead of (2) since ISTART is after 'NORMAL'
                        CALL STOKEN ( SUBR_NAME, CARD1, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
                        IF   (TOKEN(1) == 'MODEL   ') THEN ! This is TOKEN(1) instead of (2) since ISTART is after 'NORMAL'
                           SOL_NAME(1:LEN(SOL_NAME)) = ' '
                           SOL_NAME(1:12)= 'GEN CB MODEL'
                        ELSE
                           FATAL_ERR = FATAL_ERR + 1
                           WRITE(ERR,1030)
                           WRITE(F06,1030)
                           CYCLE
                        ENDIF
                     ENDIF
                  ELSE
                     WRITE(ERR,1017) TOKEN(1)
                     WRITE(F06,1017) TOKEN(1)
                     FATAL_ERR = FATAL_ERR + 1
                     CYCLE
                  ENDIF

               ELSE IF (TOKTYP(1) == 'INTEGER ') THEN      ! TOKEN(1) is the integer SOL_INT value

                  CHARFLD(1:) = ' '
                  CHARFLD(1:) = TOKEN(1)(1:)
                  CALL I4FLD ( CHARFLD, JF(2), SOL_INT )
                  CALL CRDERR ( CARD1 )
                  IF      ((SOL_INT == 1 ) .OR. (SOL_INT == 101)) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:7) = 'STATICS'
                  ELSE IF ((SOL_INT == 3 ) .OR. (SOL_INT == 103))  THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:5) = 'MODES'
                  ELSE IF ((SOL_INT == 4 ) .OR. (SOL_INT == 104)) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:8) = 'DIFFEREN'
                  ELSE IF ((SOL_INT == 5 ) .OR. (SOL_INT == 105)) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME(1:8) = 'BUCKLING'
                  ELSE IF (SOL_INT == 31) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME      = 'GEN CB MODEL'
                  ELSE IF ((SOL_INT == 66) .OR. (SOL_INT == 106)) THEN
                     SOL_NAME(1:LEN(SOL_NAME)) = ' '
                     SOL_NAME      = 'NLSTATIC'
                  ELSE
                     WRITE(ERR,999) SOL_NUM_SHOULD_BE, SOL_INT
                     WRITE(F06,999) SOL_NUM_SHOULD_BE, SOL_INT
                     FATAL_ERR = FATAL_ERR + 1
                     CYCLE
                  ENDIF

               ELSE
                  WRITE(ERR,1017) TOKEN(1)
                  WRITE(F06,1017) TOKEN(1)
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ENDIF
            ENDIF

         ELSE IF (CARD1(1: 4) == 'TIME'      ) THEN
            CONTINUE

         ELSE                                              ! Card is not recognized.
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,9993) PROG_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9993) PROG_NAME
            ENDIF
         ENDIF

      ENDDO

      IF (DOLLAR_WARN == 'Y') THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1101) CARD
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1101) CARD
         ENDIF
      ENDIF

! Check to make sure that a SOL card was read in EXEC CONTROL.

      IF (SOL_NAME(1:) == ' ') THEN
        WRITE(ERR,1016)
        WRITE(F06,1016)
        FATAL_ERR = FATAL_ERR + 1
      ENDIF

! Make sure user has not requested RESTART in nonlinear analyses

      IF ((RESTART == 'Y') .AND. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1102) SOL_NAME
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1102) SOL_NAME
         ENDIF
         RESTART = 'N'
      ENDIF

! List valid OUTPUT4 names if needed

      IF (PRT_OU4_VALID_NAMES == 'Y') THEN
         WRITE(F06,*)
         WRITE(F06,1028)
         DO I=1,NUM_OU4_VALID_NAMES
            IF (ALLOW_OU4_MYSTRAN_NAMES(I) == ALLOW_OU4_OUTPUT_NAMES(I)) THEN
               WRITE(F06,996) I, ALLOW_OU4_MYSTRAN_NAMES(I)
            ELSE
               WRITE(F06,997) I, ALLOW_OU4_MYSTRAN_NAMES(I), ALLOW_OU4_OUTPUT_NAMES(I)
            ENDIF
         ENDDO
         WRITE(F06,*)
      ENDIF

! NLSTATIC not completed yet so give error and quit

      if (SOL_NAME == 'NLSTATIC') then
         if (debug(202) == 0) then                         ! Use  non-advertised debug(202) /= 0 to continue testing NLSTATIC
            Write(err,1103)
            Write(f06,1103)
            ec_output4_err = ec_output4_err + 1
         endif
      endif

! If EC_OUTPUT4_ERR > 0 quit

      IF (EC_OUTPUT4_ERR > 0) THEN
         WRITE(ERR,998) EC_OUTPUT4_ERR
         WRITE(F06,998) EC_OUTPUT4_ERR
         CALL OUTA_HERE ( 'Y' )
      ENDIF


      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  996 FORMAT(' (',I3,')',5X,A)

  997 FORMAT(' (',I3,')',5X,A,' or its NASTRAN alias NAME: ',A)

  998 FORMAT(' PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS')

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE ',A,' BUT IS SOL = ',I8)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

 1016 FORMAT(' *ERROR  1016: A CORRECT SOL ENTRY WAS NOT FOUND IN THE EXEC CONTROL DECK')

 1017 FORMAT(' *ERROR  1017: INCORRECT ENTRY : ',A8,' ON SOL EXECUTIVE CONTROL ENTRY')

 1027 FORMAT(' *ERROR  1027: SOL EXEC CONTROL ENTRY MUST HAVE ONLY ONE INTEGER NUMBER FOLLOWING SOL.',                             &
                           ' ABOVE ENTRY DOES NOT MEET THIS REQUIREMENT')

 1028 FORMAT(' The list of valid names of matrices for OUTPUT4 are:',/)

 1030 FORMAT(' *ERROR  1030: THE ENTRY "GEN CB MODEL" ON SOL ENTRY FOR CRAIG-BAMPTON MODEL GENERATION MAY HAVE BEEN MISSPELLED')

 1101 FORMAT(' *WARNING    : BE CAREFUL WITH LINES THAT BEGIN WITH A $ SIGN IN COL 1 FOLLOWED BY AN UPPER CASE LETTER IN EXEC OR', &
                           ' CASE CONTROL.'                                                                                        &
                    ,/,14X,' THE LINE CAN BE MISINTERPRETED AS A DIRECTIVE FOR THE BANDIT GRID RESEQUENCING ALGORITHM.'            &
                    ,/,14X,' SEE THE BANDIT.PDF FILE INSTALLED WHEN YOU RAN SETUP.EXE TO INSTALL MYSTRAN' &
                    ,/,14X, A)

 1102 FORMAT(' *WARNING    : REQUEST FOR RESTART IS NOT ALLOWED IN SOLUTION ',A)

 1103 format(' *INFORMATION: NLSTATIC code not completed so quiting')

 9901 FORMAT(' *WARNING:   : EXEC CONTROL ENTRY "SOL" CANNOT BE "',A,'". SOL ',A,' WILL BE USED')

 9993 FORMAT(' *WARNING    : PRIOR ENTRY NOT PROCESSED BY ',A)


! **********************************************************************************************************************************

      END SUBROUTINE LOADE


      SUBROUTINE LOADE0

      ! LOADE0 does a preliminary read of the EXEC CONTROL DECK to find
      ! if there is a RESTART entry
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, FILE_NAM_MAXLEN, IN0, IN1, INC, LEN_INPUT_FNAME, INFILE,           &
                                         LEN_RESTART_FNAME, LNUM_IN4_FILES, RESTART_FILNAM, SCR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN, FATAL_ERR, RESTART
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT, REPLACE_TABS_W_BLANKS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LOADE0'
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD              ! Exec Control deck card
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'EXEC CONTROL'
      CHARACTER( 4*BYTE), PARAMETER   :: END_CARD  = 'CEND'
      CHARACTER(LEN(INFILE))          :: FILNAM

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IBEG              ! Col where FILNAM begins after leading blanks
      INTEGER(LONG)                   :: IEND              ! Col where FILNAM ends after trailing blanks
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator.
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1




! **********************************************************************************************************************************
      ! Initialize
      LNUM_IN4_FILES = 0
      RESTART        = 'N'

!xx   REWIND (IN1)
main: DO
         CALL READ_BDF_LINE(IN1, IOCHK, CARD)

         IF (IOCHK < 0) THEN                               ! Quit if EOF/EOR occurs during read
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (IOCHK > 0) THEN                               ! Check if error occurs during read.
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         CALL REPLACE_TABS_W_BLANKS ( CARD )               ! Replace all tab characters with a white space

         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )  ! Shift CARD so it begins in col 1

         IF (CARD1(1:7) == 'RESTART'   ) THEN              ! No errors, so look for RESTART
            RESTART = 'Y'
            FILNAM(1:) = CARD1(8:)

            IBEG = 1
i_do2:      DO I=1,FILE_NAM_MAXLEN
               IF (FILNAM(I:I) == ' ') THEN
                  CYCLE
               ELSE
                  IBEG = I
                  EXIT i_do2
               ENDIF
            ENDDO i_do2

            IEND = 1
i_do3:      DO I=FILE_NAM_MAXLEN,1,-1
               IF (FILNAM(I:I) == ' ') THEN
                  CYCLE
               ELSE
                  IEND = I
                  EXIT i_do3
               ENDIF
            ENDDO i_do3

            IF (IEND == IBEG) THEN
               LEN_RESTART_FNAME = LEN_INPUT_FNAME
               RESTART_FILNAM = INFILE(1:LEN_INPUT_FNAME)
            ELSE
               LEN_RESTART_FNAME = IEND - IBEG + 2
               RESTART_FILNAM(1:) = FILNAM(IBEG:IEND) // '.'
            ENDIF

         ELSE IF (CARD1(1:3) == 'IN4') THEN
            LNUM_IN4_FILES = LNUM_IN4_FILES + 1

         ELSE IF (CARD1(1:4) == 'CEND' ) THEN              ! Check for CEND card
            EXIT main

         ENDIF

      ENDDO main



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

! **********************************************************************************************************************************

      END SUBROUTINE LOADE0

   END MODULE INPUT_STAGE_DISPATCH

      SUBROUTINE ELEPRO ( INCR_NELE, JCARD, NFIELD, NMORE, CHK_FLD2, CHK_FLD3, CHK_FLD4, CHK_FLD5, CHK_FLD6, CHK_FLD7, CHK_FLD8, CHK_FLD9 )
      USE PENTIUM_II_KIND, ONLY: LONG
      USE INPUT_STAGE_DISPATCH, ONLY: ELEPRO_MOD => ELEPRO
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: INCR_NELE, CHK_FLD2, CHK_FLD3, CHK_FLD4, CHK_FLD5, CHK_FLD6, CHK_FLD7, CHK_FLD8, CHK_FLD9
      CHARACTER(LEN=*), INTENT(IN) :: JCARD(10)
      INTEGER(LONG), INTENT(IN) :: NFIELD, NMORE
      CALL ELEPRO_MOD ( INCR_NELE, JCARD, NFIELD, NMORE, CHK_FLD2, CHK_FLD3, CHK_FLD4, CHK_FLD5, CHK_FLD6, CHK_FLD7, CHK_FLD8, CHK_FLD9 )
      END SUBROUTINE ELEPRO

      SUBROUTINE READ_BDF_LINE ( IN1, IOCHK, LINE )
      USE INPUT_STAGE_DISPATCH, ONLY: READ_BDF_LINE_MOD => READ_BDF_LINE
      IMPLICIT NONE
      CHARACTER(256), INTENT(INOUT) :: LINE
      INTEGER, INTENT(IN) :: IN1
      INTEGER, INTENT(INOUT) :: IOCHK
      CALL READ_BDF_LINE_MOD ( IN1, IOCHK, LINE )
      END SUBROUTINE READ_BDF_LINE
