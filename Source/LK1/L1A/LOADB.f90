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


      USE LOADB_USE_IFs

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
            CALL BD_SPOINT  ( CARD )

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
      USE TO_UPPER_Interface

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
