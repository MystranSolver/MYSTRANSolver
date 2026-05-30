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

      SUBROUTINE BD_PARAM ( CARD )

! Processes PARAM Bulk Data Cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEPSIL, MPBARLU, NUM_USETSTR,       &
                                         WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MACHINE_PARAMS, ONLY        :  MACH_PREC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE PARAMS, ONLY                :  ARP_TOL         , ART_KED         , ART_ROT_KED     , ART_TRAN_KED    ,                   &
                                         ART_MASS        , ART_ROT_MASS    , ART_TRAN_MASS   , AUTOSPC         , AUTOSPC_NSET    , &
                                         AUTOSPC_RAT     , AUTOSPC_INFO    , AUTOSPC_SPCF    , BAILOUT         , BANDEDOPT       , &
                                         CRS_CCS         , &
                                         CBMIN3          , CBMIN4          , CBMIN4T         , CHKGRDS         ,                   &
                                         CUSERIN         , CUSERIN_EID     , CUSERIN_IN4     , CUSERIN_PID     , CUSERIN_SPNT_ID , &
                                         CUSERIN_XSET    , CUSERIN_COMPTYP , DARPACK         ,                                     &
                                         DELBAN          , EIGESTL         , EIGNORM2        , ELFORCEN        , EPSERR          , &
                                         EQCHK_REF_GRID  , EQCHK_NORM      , EQCHK_OUTPUT    , EQCHK_TINY      ,                   &
                                         EPSIL           , EMP0_PAUSE      , ESP0_PAUSE      , F06_COL_START   ,                   &
                                         GRDPNT          , GRDPNT_IN       , GRIDSEQ         , HEXAXIS         ,                   &
                                         IORQ1M          , IORQ1S          , IORQ1B          , IORQ2B          , IORQ2T          , &
                                         ITMAX           , KLLRAT          , KOORAT          , LANCMETH        , MATSPARS        , &
                                         MEMAFAC         , MIN4TRED        , MXALLOCA        , MAXRATIO        ,                   &
                                         MEFMCORD        , MEFMLOC         , MEFMGRID        ,                                     &
                                         MPFOUT          , MXITERI         , MXITERL         , OTMSKIP         , POST            , &
                                         PBARLDEC        , PBARLSHR        , PCOMPEQ         , PCHSPC1         , PCMPTSTM        , &
                                         PRTBASIC        , PRTCONN         , PRTCORD         , PRTDISP         , PRTDLR          , &
                                         PRTDOF          , PRTFOR          , PRTHMN          , PRTGMN          , PRTGOA          , &
                                         PRTCGLTM        , PRTPHIZL        , PRTIFLTM        , PRTKXX          , PRTMXX          , &
                                         PRTOU4          , PRTPHIXA        , PRTMASS         , PRTMASSD        , PRTRMG          , &
                                         PRTSCP          , PRTPSET         , PRTTSET         , PRTUSET         ,                   &
                                         PRTSTIFD        , PRTSTIFF        , PRTUO0          ,                                     &
                                         PRTYS           , PRTQSYS         ,                                                       &
                                         Q4SURFIT        , QUADAXIS        , QUAD4TYP        , RCONDK          , RELINK3         , &
                                         SEQPRT          , SEQQUIT         , SETLKTM         , SETLKTK         , SHRFXFAC        , &
                                         SKIPMGG         , SOLLIB          , SPARSE_FLAVOR   , SPARSTOR        ,                   &
                                         SPC1QUIT        , SORT_MAX        , SPC1SID         , STR_CID                           , &
                                         SUPINFO         , SUPWARN         , NOCOUNTS                                            , &
                                         THRESHK         , THRESHK_LAP     , TINY            ,                                     &
                                         TSTM_DEF        , USR_JCT         , USR_LTERM_KGG   , USR_LTERM_MGG   , WINAMEM         , &
                                         WTMASS          , K6ROT,                                                                  &
                                         PRTALL          , PRTF06          , PRTNEU          , PRTOP2          ,                   &
                                         SPIENV6         , SPIENV7         , SPIENV8         , SLU_NTHR

      USE BD_PARAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PARAM'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHRPARM           ! A char value read from a PARAM Bulk Data card
      CHARACTER(15*BYTE)              :: PARNAM            ! The name of a parameter
      CHARACTER(8*BYTE)               :: PARAM_NAME        ! The name of a parameter

      INTEGER(LONG)                   :: I4PARM    = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: I4PARMI(5)        ! Values read from input file that should be an integer value
      INTEGER(LONG)                   :: LOWER             ! Lower allowable value for an integer parameter
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR      = 0     ! Local error indicator
      INTEGER(LONG)                   :: II                ! An index in array EPSIL
      INTEGER(LONG)                   :: UPPER             ! Upper allowable value for an integer parameter


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: R8PARM            ! A value read from input file that should be a real value

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! PARAM Bulk Data Card routine

      EPS1 = EPSIL(1)

! Check for PARAM cards. If PARAM name is not legal, write warning

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      PARAM_NAME = JCARD(2)


      ! ARP_TOL is the variable TOL used in ARPACK routines to decide convergence.
      ! NOTE: if a value of -1. is input on a PARAM ARP_TOL entry, then the Lanczos algorithm will use machine precision for ARP_TOL
      IF      (JCARD(2)(1:7) == 'ARP_TOL') THEN
         PARNAM = 'ARP_TOL  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - ARP_TOL) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,ARP_TOL,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,ARP_TOL,R8PARM
               ENDIF
               ARP_TOL = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! ART_KED is used to add artificial differential stiffness terms to translation and/or rotation DOF's in the element KED matrix
! (to make sure there are no zero's on diag of element KED matrices)

      ELSE IF (JCARD(2)(1:8) == 'ART_KED ') THEN

         PARNAM = 'ART_KED'
         ART_KED = ' '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM(1:1) == 'Y') THEN
               ART_KED = 'Y'
            ELSE IF (CHRPARM(1:1) == 'N') THEN
               ART_KED = 'N'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM, 'Y OR N', CHRPARM, ART_KED
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM, 'Y OR N', CHRPARM, ART_KED
               ENDIF
            ENDIF
         ENDIF

         IF (ART_KED == 'Y') THEN

            IF (JCARD(4)(1:) /= ' ') THEN
               PARNAM = 'ART_TRAN_KED'
               CALL R8FLD ( JCARD(4), JF(4), R8PARM )
               IF (IERRFL(4) == 'N') THEN
                  WRITE(ERR,1147) PARNAM, ART_TRAN_KED, R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM, ART_TRAN_KED, R8PARM
                  ENDIF
                  ART_TRAN_KED = R8PARM
               ENDIF
            ENDIF

            IF (JCARD(5)(1:) /= ' ') THEN
               PARNAM = 'ART_ROT_KED'
               CALL R8FLD ( JCARD(5), JF(5), R8PARM )
               IF (IERRFL(5) == 'N') THEN
                  WRITE(ERR,1147) PARNAM, ART_ROT_KED, R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM, ART_ROT_KED, R8PARM
                  ENDIF
                  ART_ROT_KED = R8PARM
               ENDIF
            ENDIF

         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! ART_MASS is used to add artificial mass terms to transl and/or rot DOF's (to make sure there are no zero's on diag of MGG)

      ELSE IF (JCARD(2)(1:8) == 'ART_MASS ') THEN

         PARNAM = 'ART_MASS'
         ART_MASS = ' '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM(1:1) == 'Y') THEN
               ART_MASS = 'Y'
            ELSE IF (CHRPARM(1:1) == 'N') THEN
               ART_MASS = 'N'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM, 'Y OR N', CHRPARM, ART_MASS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM, 'Y OR N', CHRPARM, ART_MASS
               ENDIF
            ENDIF
         ENDIF

         IF (ART_MASS == 'Y') THEN

            IF (JCARD(4)(1:) /= ' ') THEN
               PARNAM = 'ART_TRAN_MASS'
               CALL R8FLD ( JCARD(4), JF(4), R8PARM )
               IF (IERRFL(4) == 'N') THEN
                  WRITE(ERR,1147) PARNAM, ART_TRAN_MASS, R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM, ART_TRAN_MASS, R8PARM
                  ENDIF
                  ART_TRAN_MASS = R8PARM
               ENDIF
            ENDIF

            IF (JCARD(5)(1:) /= ' ') THEN
               PARNAM = 'ART_ROT_MASS'
               CALL R8FLD ( JCARD(5), JF(5), R8PARM )
               IF (IERRFL(5) == 'N') THEN
                  WRITE(ERR,1147) PARNAM, ART_ROT_MASS, R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM, ART_ROT_MASS, R8PARM
                  ENDIF
                  ART_ROT_MASS = R8PARM
               ENDIF
            ENDIF

         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! AUTOSPC = 'Y' causes  MYSTRAN to SPC DOF's

      ELSE IF (JCARD(2)(1:8) == 'AUTOSPC ') THEN
         PARNAM = 'AUTOSPC '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM(1:1) == 'Y') THEN
               AUTOSPC = 'Y'
            ELSE IF (CHRPARM(1:1) == 'N') THEN
               AUTOSPC = 'N'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(4)(1:) /= ' ') THEN
            PARNAM = 'AUTOSPC_RAT'
            CALL R8FLD ( JCARD(4), JF(4), R8PARM )
            IF (IERRFL(4) == 'N') THEN
               IF (DABS(R8PARM - AUTOSPC_RAT) > EPS1) THEN
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,1147) PARNAM, AUTOSPC_RAT, R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM, AUTOSPC_RAT, R8PARM
                  ENDIF
                  AUTOSPC_RAT = R8PARM
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN
            PARNAM = 'AUTOSPC_NSET'
            CALL I4FLD ( JCARD(5), JF(5), I4PARM )
            IF (IERRFL(5) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
                  AUTOSPC_NSET = I4PARM
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(5),AUTOSPC_NSET
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(5),AUTOSPC_NSET
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(6)(1:) /= ' ') THEN
            PARNAM = 'AUTOSPC_INFO'
            CALL CHAR_FLD ( JCARD(6), JF(6), CHRPARM )
            IF (IERRFL(6) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  AUTOSPC_INFO = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  AUTOSPC_INFO = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_INFO
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_INFO
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(7)(1:) /= ' ') THEN
            PARNAM = 'AUTOSPC_SPCF'
            CALL CHAR_FLD ( JCARD(7), JF(7), CHRPARM )
            IF (IERRFL(7) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  AUTOSPC_SPCF = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  AUTOSPC_SPCF = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_SPCF
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_SPCF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in fields 3-6
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 7-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! BAILOUT defines whether to quit if a singularity is found in matrix decomp

      ELSE IF (JCARD(2)(1:8) == 'BAILOUT ') THEN
         PARNAM = 'BAILOUT '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= BAILOUT) THEN
               WRITE(ERR,1146) PARNAM,BAILOUT,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,BAILOUT,I4PARM
               ENDIF
               BAILOUT = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! BANDEDOPT enables experimental banded-order optimization path

      ELSE IF (JCARD(2)(1:8) == 'BANDEDOP') THEN
         PARNAM = 'BANDEDOPT'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, BANDEDOPT)

! CBMIN3 is a parameter for the Mindlin (thick) triangular plate element (CTRIA3).
!   It is used in calculating PHISQ, a scalar multiple of the transverse shear stiff

      ELSE IF (JCARD(2)(1:8) == 'CBMIN3  ') THEN
         PARNAM = 'CBMIN3  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - CBMIN3) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,CBMIN3,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,CBMIN3,R8PARM
               ENDIF
               CBMIN3 = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! CBMIN4 is a parameter for the Mindlin (thick) quad plate element (CQUAD4).
! It is used in calculating PHISQ, a scalar multiple of the transverse shear stiff

      ELSE IF (JCARD(2)(1:8) == 'CBMIN4  ') THEN
         PARNAM = 'CBMIN4  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - CBMIN4) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,CBMIN4,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,CBMIN4,R8PARM
               ENDIF
               CBMIN4 = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! CBMIN4T is a parameter for the Mindlin (thick) quad plate element (CQUAD4).
      ! It is used in calculating PHISQ, a scalar multiple of the transverse shear stiff
      ELSE IF (JCARD(2)(1:8) == 'CBMIN4T ') THEN
         PARNAM = 'CBMIN4T '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - CBMIN4T) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,CBMIN4T,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,CBMIN4T,R8PARM
               ENDIF
               CBMIN4T = R8PARM
            ENDIF
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! CHKGRDS tells whether to call GET_ELEM_AGRID_BGRID to make sure
      ! all grids on elem connection entries are defined
      ELSE IF (JCARD(2)(1:8) == 'CHKGRDS ') THEN
         PARNAM = 'CHKGRDS  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, CHKGRDS)


      ! CRS_CCS specifies the storage method (CRS or CCS)
      ELSE IF (JCARD(2)(1:8) == 'CRS_CCS ') THEN
         PARNAM = 'CRS_CCS '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'CRS     ') THEN
               CRS_CCS  = 'CRS'
            ELSE IF (CHRPARM == 'CCS     ') THEN
               CRS_CCS  = 'CCS'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'CRS or CCS',CHRPARM,CRS_CCS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'CRS or CCS',CHRPARM,CRS_CCS
               ENDIF
            ENDIF
         ENDIF

! CUSERIN = 'Y' causes  MYSTRAN to write CUSERIN, PUSERIN B.D. card images for a substructure to F06

      ELSE IF (JCARD(2)(1:8) == 'CUSERIN ') THEN
         CUSERIN = 'Y'
         PARNAM = 'CUSERIN '
         IF (JCARD(3)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_EID'
            CALL I4FLD ( JCARD(3), JF(3), I4PARM )
            IF (IERRFL(3) == 'N') THEN
               WRITE(ERR,1146) PARNAM, CUSERIN_EID, I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM, CUSERIN_EID, I4PARM
               ENDIF
               CUSERIN_EID = I4PARM
            ENDIF
         ENDIF

         IF (JCARD(4)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_PID'
            CALL I4FLD ( JCARD(4), JF(4), I4PARM )
            IF (IERRFL(4) == 'N') THEN
               WRITE(ERR,1146) PARNAM, CUSERIN_PID, I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM, CUSERIN_PID, I4PARM
               ENDIF
               CUSERIN_PID = I4PARM
            ENDIF
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_SPNT_ID'
            CALL I4FLD ( JCARD(5), JF(5), I4PARM )
            IF (IERRFL(5) == 'N') THEN
               WRITE(ERR,1146) PARNAM, CUSERIN_SPNT_ID, I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM, CUSERIN_SPNT_ID, I4PARM
               ENDIF
               CUSERIN_SPNT_ID = I4PARM
            ENDIF
         ENDIF

         IF (JCARD(6)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_IN4'
            CALL I4FLD ( JCARD(6), JF(6), I4PARM )
            IF (IERRFL(6) == 'N') THEN
               WRITE(ERR,1146) PARNAM, CUSERIN_IN4, I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM, CUSERIN_IN4, I4PARM
               ENDIF
               CUSERIN_IN4 = I4PARM
            ENDIF
         ENDIF

         IF (JCARD(7)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_XSET'
            CALL CHAR_FLD ( JCARD(7), JF(7), CHRPARM )
            IF (IERRFL(7) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               CUSERIN_XSET = CHRPARM(1:LEN(TSET_CHR_LEN))
            ENDIF
         ENDIF

         IF (JCARD(8)(1:) /= ' ') THEN
            PARNAM = 'CUSERIN_COMPTYP'
            CALL I4FLD ( JCARD(8), JF(8), I4PARM )
            IF (IERRFL(8) == 'N') THEN
               CUSERIN_COMPTYP = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in fields 3-6
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! Delta to add to EIG_N2 so that ARPACK will find a few more eigen
      ! than user requested (since higher ones seem to be a little bad)
      ELSE IF (JCARD(2)(1:8) == 'DARPACK ') THEN
         PARNAM = 'DARPACK '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM >= 0) THEN
               DARPACK = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1172) PARNAM,'>= 0',I4PARM,DARPACK
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= 0',I4PARM,DARPACK
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! DELBAN, where if it is 1, delete the bandit files left over, if 0, don't.
      ELSE IF (JCARD(2)(1:8) == 'DELBAN  ') THEN
         PARNAM = 'DELBAN  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= DELBAN) THEN
               WRITE(ERR,1146) PARNAM,DELBAN,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,DELBAN,I4PARM
               ENDIF
               DELBAN = I4PARM
            ENDIF
         ENDIF



         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! EIGESTL defines upper limit on NDOFL for running code to estimate the number of eigens < EIG_FRQ2
      ELSE IF (JCARD(2)(1:8) == 'EIGESTL ') THEN
         PARNAM = 'EIGESTL '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= EIGESTL) THEN
               WRITE(ERR,1146) PARNAM,EIGESTL,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,EIGESTL,I4PARM
               ENDIF
               EIGESTL = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! EIGNORM2 tells whether to rescale eigenvecors by a set of user specified +/-1 scale factors

      ELSE IF (JCARD(2)(1:8) == 'EIGNORM2') THEN
         PARNAM = 'EIGNORM2'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, EIGNORM2)

! ELFORCEN changes default from global to local or basic coords for calculating element nodal forces in LINK9

      ELSE IF (JCARD(2)(1:8) == 'ELFORCEN') THEN
         PARNAM = 'ELFORCEN'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'LOCAL') THEN
               ELFORCEN = 'LOCAL'
            ELSE IF (CHRPARM == 'BASIC   ') THEN
               ELFORCEN = 'BASIC'
            ELSE IF (CHRPARM == 'GLOBAL  ') THEN
               ELFORCEN = 'GLOBAL'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'BASIC OR GLOBAL',CHRPARM,ELFORCEN
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'BASIC OR GLOBAL',CHRPARM,ELFORCEN
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! EPSERR tells whether to calculate the NASTRAN like epsilon error estimate
      ELSE IF (JCARD(2)(1:8) == 'EPSERR  ') THEN
         PARNAM = 'EPSERR  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, EPSERR)


! EPSIL are roundoff numbers used in comparing computed values to zero

      ELSE IF (JCARD(2)(1:8) == 'EPSIL   ') THEN
         PARNAM = 'EPSIL   '
         CALL I4FLD ( JCARD(3), JF(3), II )
         IF (IERRFL(3) == 'N') THEN
            IF ((II > 0) .AND. (II <= MEPSIL)) THEN
               CALL R8FLD ( JCARD(4), JF(4), R8PARM )
               IF (R8PARM /= EPSIL(II)) THEN
                  WRITE(ERR,1148) PARNAM,II,EPSIL(II),R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1148) PARNAM,II,EPSIL(II),R8PARM
                  ENDIF
                  EPSIL(II) = R8PARM
               ENDIF
            ELSE
               WRITE(ERR,101) CARD
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1149) MEPSIL
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1149) MEPSIL
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )! Issue warning if fields 5-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! EQCHECK requests that an equilibrium check be made of the stiffness matrices

      ELSE IF (JCARD(2)(1:8) == 'EQCHECK ') THEN
         EQCHK_REF_GRID = 0                                ! Reset default value to zero from what it was set in module PARAMS so we
!                                                            can set DEFAULT grid to basic origin if user asks for equil check
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         PARNAM = 'EQCHK_REF_GRID'
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM >= -1) THEN
               EQCHK_REF_GRID = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'>= -1',I4PARM,EQCHK_REF_GRID
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= -1',I4PARM,EQCHK_REF_GRID
               ENDIF
            ENDIF
         ENDIF

         DO I=1,5
            EQCHK_OUTPUT(I) = 0
            IF (JCARD(I+3)(1:) /= ' ') THEN
               IF      (I == 1) THEN
                  PARNAM = 'EQCHK_OUTPUT(1)'
               ELSE IF (I == 2 ) THEN
                  PARNAM = 'EQCHK_OUTPUT(2)'
               ELSE IF (I == 3 ) THEN
                  PARNAM = 'EQCHK_OUTPUT(3)'
               ELSE IF (I == 4 ) THEN
                  PARNAM = 'EQCHK_OUTPUT(4)'
               ELSE IF (I == 5 ) THEN
                  PARNAM = 'EQCHK_OUTPUT(5)'
               ENDIF
               CALL I4FLD ( JCARD(I+3), JF(I+3), I4PARMI(I) )
               IF (IERRFL(I+3) == 'N') THEN
                  LOWER = 0
                  UPPER = 3
                  IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                     EQCHK_OUTPUT(I) = I4PARMI(I)
                  ELSE
                     WARN_ERR = WARN_ERR + 1
                     WRITE(ERR,101) CARD
                     WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+3),EQCHK_OUTPUT(I)
                     IF (SUPWARN == 'N') THEN
                        IF (ECHO == 'NONE  ') THEN
                           WRITE(F06,101) CARD
                        ENDIF
                        WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+3),EQCHK_OUTPUT(I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         IF (JCARD(9)(1:) /= ' ') THEN
            PARNAM = 'EQCHK_TINY'
            CALL R8FLD ( JCARD(9), JF(9), R8PARM )
            IF (IERRFL(9) == 'N') THEN
               IF (DABS(R8PARM - EQCHK_TINY) > EPS1) THEN
                  WRITE(ERR,1147) PARNAM,EQCHK_TINY,R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM,EQCHK_TINY,R8PARM
                  ENDIF
                  EQCHK_TINY = R8PARM
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(10)(1:) /= ' ') THEN
            PARNAM = 'EQCHK_NORM'
            CALL CHAR_FLD ( JCARD(10), JF(10), CHRPARM )
            IF (IERRFL(10) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  EQCHK_NORM = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  EQCHK_NORM = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,EQCHK_NORM
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,EQCHK_NORM
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,8,9 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! F06COL is the 1st col in F06 file for output data to begin. If it is not > 2, then
! output will be written with each main header centered on one another

!zzzz ELSE IF (JCARD(2)(1:8) == 'F06COL  ') THEN
!zzzz    PARNAM = 'F06_COL_START'
!zzzz    CALL I4FLD ( JCARD(3), JF(3), I4PARM )
!zzzz    IF (IERRFL(3) == 'N') THEN
!zzzz       LOWER = 2
!zzzz       UPPER = 99999999
!zzzz       IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
!zzzz          F06_COL_START = I4PARM
!zzzz       ELSE
!zzzz          WARN_ERR = WARN_ERR + 1
!zzzz          WRITE(ERR,101) CARD
!zzzz          WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),F06_COL_START
!zzzz          IF (SUPWARN == 'N') THEN
!zzzz             IF (ECHO == 'NONE  ') THEN
!zzzz                WRITE(F06,101) CARD
!zzzz             ENDIF
!zzzz             WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),F06_COL_START
!zzzz          ENDIF
!zzzz       ENDIF
!zzzz    ENDIF

!zzzz    CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
!zzzz    CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
!zzzz    CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTALL writes all outputs for all output files regardless of other flags
      ELSE IF ((PARAM_NAME(1:8) == 'PRTALL  ') .OR. (PARAM_NAME(1:8) == 'FILES   ')) THEN
         PARNAM = 'PRTALL  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, PRTALL)
         IF (PRTALL == 'Y') THEN
             PRTF06 = 'Y'
             PRTNEU = 'Y'
             PRTOP2 = 'Y'
         ENDIF

      ! PRTOP2 writes all outputs for the f06 file regardless of other flags besides PRTALL
      ELSE IF ((PARAM_NAME(1:8) == 'PRTF06  ') .OR. (PARAM_NAME(1:8) == 'F06     ')) THEN
         PARNAM = 'PRTF06  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, PRTF06)

      ! PRTNEU writes all outputs for the neu file regardless of other flags besides PRTALL
      ELSE IF ((PARAM_NAME(1:8) == 'PRTNEU  ') .OR. (PARAM_NAME(1:8) == 'NEU     ')) THEN
         PARNAM = 'PRTNEU  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, PRTNEU)

      ! PRTOP2 writes all outputs for the op2 file regardless of other flags besides PRTALL
      ELSE IF ((PARAM_NAME(1:8) == 'PRTOP2  ') .OR. (PARAM_NAME(1:8) == 'OP2     ')) THEN
         PARNAM = 'PRTOP2  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, PRTOP2)

      ! GRDPNT causes the grid point weight generator to be run to calculate mass of the model relative to G.P defined by PARAM GRDPNT.
      ELSE IF (JCARD(2)(1:8) == 'GRDPNT  ') THEN
         ! GRDPNT_IN is used to decide if GPWG is run in LINK0
         PARNAM = 'GRDPNT  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM >= 0) THEN
               GRDPNT    = I4PARM
               GRDPNT_IN = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'>= 0',I4PARM,GRDPNT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= 0',I4PARM,GRDPNT
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! GRIDSEQ specifies the grid point sequencing method

      ELSE IF (JCARD(2)(1:8) == 'GRIDSEQ ') THEN
         PARNAM = 'GRIDSEQ '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'BANDIT  ') THEN
               GRIDSEQ = 'BANDIT  '
            ELSE IF (CHRPARM == 'GRID    ') THEN
               GRIDSEQ = 'GRID    '
            ELSE IF (CHRPARM == 'INPUT   ') THEN
               GRIDSEQ = 'INPUT   '
            ELSE IF (CHRPARM(1:3) == 'RCM') THEN
               GRIDSEQ = 'RCM     '
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'BANDIT, GRID, INPUT or RCM',CHRPARM,GRIDSEQ
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'BANDIT, GRID, INPUT or RCM',CHRPARM,GRIDSEQ
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(4)(1:) /= ' ') THEN
            PARNAM = 'SEQQUIT '
            CALL CHAR_FLD ( JCARD(4), JF(4), CHRPARM )
            IF (IERRFL(4) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  SEQQUIT = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  SEQQUIT = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,SEQQUIT
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,SEQQUIT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN
            PARNAM = 'SEQPRT  '
            CALL CHAR_FLD ( JCARD(5), JF(5), CHRPARM )
            IF (IERRFL(5) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  SEQPRT = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  SEQPRT = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,SEQPRT
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,SEQPRT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! HEXAXIS sets the method for how the local x axis of quad elements is determined

      ELSE IF (JCARD(2)(1:8) == 'HEXAXIS ') THEN
         PARNAM = 'HEXAXIS'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'SPLITD  ') THEN
               HEXAXIS = 'SPLITD'
            ELSE IF (CHRPARM == 'SIDE12  ') THEN
               HEXAXIS = 'SIDE12'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'SPLITD or SIDE12',CHRPARM,HEXAXIS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'SPLITD or SIDE12',CHRPARM,HEXAXIS
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! IORQ1M defines the Gaussian integration order for in-plane direct stresses for isoparametric quad plate elements in subr QMEM1

      ELSE IF (JCARD(2)(1:8) == 'IORQ1M  ') THEN
         PARNAM = 'IORQ1M  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= IORQ1M) THEN
               WRITE(ERR,1146) PARNAM,IORQ1M,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,IORQ1M,I4PARM
               ENDIF
               IORQ1M = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! IORQ1S defines the Gaussian integration order for in-plane shear stresses for isoparametric quad plate elems in subr QMEM1

      ELSE IF (JCARD(2)(1:8) == 'IORQ1S  ') THEN
         PARNAM = 'IORQ1S  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= IORQ1S) THEN
               WRITE(ERR,1146) PARNAM,IORQ1S,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,IORQ1S,I4PARM
               ENDIF
               IORQ1S = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! IORQ1B defines the Gaussian integration order for bending stresses for Kirchoff (thin) isoparam quad plate elems in subr QPLT1

      ELSE IF (JCARD(2)(1:8) == 'IORQ1B  ') THEN
         PARNAM = 'IORQ1B  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= IORQ1B) THEN
               WRITE(ERR,1146) PARNAM,IORQ1B,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,IORQ1B,I4PARM
               ENDIF
               IORQ1B = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! IORQ2B defines the Gaussian integration order for bending stresses for Mindlin (thick) isoparam quad plate elems in subr QPLT2

      ELSE IF (JCARD(2)(1:8) == 'IORQ2B  ') THEN
         PARNAM = 'IORQ2B  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= IORQ2B) THEN
               WRITE(ERR,1146) PARNAM,IORQ2B,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,IORQ2B,I4PARM
               ENDIF
               IORQ2B = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! IORQ2T defines the Gaussian integ order for transv shear stresses for Mindlin (thick) isoparam quad plate elems in subr QPLT2

      ELSE IF (JCARD(2)(1:8) == 'IORQ2T  ') THEN
         PARNAM = 'IORQ2T  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= IORQ2T) THEN
               WRITE(ERR,1146) PARNAM,IORQ2T,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,IORQ2T,I4PARM
               ENDIF
               IORQ2T = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! ITMAX is max iterations used in refining the soln when LAPACK is used in LINK8

      ELSE IF (JCARD(2)(1:8) == 'ITMAX   ') THEN
         PARNAM = 'ITMAX   '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM > 0) THEN
               ITMAX = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'> 0',I4PARM,ITMAX
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'> 0',I4PARM,ITMAX
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! KLLRAT tells whether to calculate max ratio of matrix diagonal to factor diagonal

      ELSE IF (JCARD(2)(1:8) == 'KLLRAT  ') THEN
         PARNAM = 'KLLRAT  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, KLLRAT)

! KOORAT tells whether to calculate max ratio of matrix diagonal to factor diagonal

      ELSE IF (JCARD(2)(1:8) == 'KOORAT  ') THEN
         PARNAM = 'KOORAT  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, KOORAT)

! MATSPARS Sets whether to use sparse matrix add and multiply routines for situations where the matrices are in any sparse format.
! If 'Y', use sparse routines. If 'N', use full matrix routines in place of all sparse routines

      ELSE IF (JCARD(2)(1:8) == 'MATSPARS') THEN
         PARNAM = 'MATSPARS'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, MATSPARS)

! MXALLOCA defines the max allowable attempts when attempting to allocate an array before quiting with a fatal error

      ELSE IF (JCARD(2)(1:8) == 'MXALLOCA') THEN
         PARNAM = 'MXALLOCA '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= MXALLOCA) THEN
               WRITE(ERR,1146) PARNAM,MXALLOCA,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,MXALLOCA,I4PARM
               ENDIF
               MXALLOCA = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MAXRATIO sets the value for comparing ratio of matrix diagonal to factor diagonal before action is taken

      ELSE IF (JCARD(2)(1:8) == 'MAXRATIO') THEN
         PARNAM = 'MAXRATIO  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - MAXRATIO) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,MAXRATIO,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,MAXRATIO,R8PARM
               ENDIF
               MAXRATIO = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MEFMCORD is the coordinate system for output of modal effective masses

      ELSE IF (JCARD(2)(1:8) == 'MEFMCORD') THEN
         PARNAM = 'MEFMCORD'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 99999999
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               MEFMCORD = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),MEFMCORD
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),MEFMCORD
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MEFMLOC tells what ref loaction to use for reducing modal effective mass to a 6 DOF location.

      ELSE IF (JCARD(2)(1:8) == 'MEFMLOC ') THEN
         PARNAM = 'MEFMLOC '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'CG      ') THEN
               MEFMLOC = 'CG    '
            ELSE IF (CHRPARM == 'GRID    ') THEN
               MEFMLOC = 'GRID  '
            ELSE IF (CHRPARM == 'GRDPNT  ') THEN
               MEFMLOC = 'GRDPNT'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'CG, GRID, GRDSET',CHRPARM,MEFMLOC
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'CG, GRID, GRDSET',CHRPARM,MEFMLOC
               ENDIF
            ENDIF
         ENDIF

         IF (MEFMLOC == 'GRID') THEN
            PARNAM = 'MEFMGRID'
            CALL I4FLD ( JCARD(4), JF(4), I4PARM )
            IF (IERRFL(4) == 'N') THEN
               IF (I4PARM >= 0) THEN
                  MEFMGRID = I4PARM
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1172) PARNAM,'>= 0',I4PARM,MEFMLOC
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1172) PARNAM,'>= 0',I4PARM,MEFMLOC
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (MEFMLOC(1:2) == 'CG') THEN
            CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )
         ELSE
            CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,0,0,0,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )
            CALL CRDERR ( CARD )
         ENDIF

! MEMAFAC is a factor to multiply the size request of memory to be allocated when looping to find an allowable amount of
! memory to allocate. Used when the initial request for memory (in subrs ESP or EMP) cannot be met and we know that the request
! is prpbably on the high side.

      ELSE IF (JCARD(2)(1:8) == 'MEMAFAC ') THEN
         PARNAM = 'MEMAFAC  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - MEMAFAC) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,MEMAFAC,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,MEMAFAC,R8PARM
               ENDIF
               MEMAFAC = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MIN4TRED sets the method for how the 5th node of the MIN4T element is reduced out (to get a 4 node quad element)

      ELSE IF (JCARD(2)(1:8) == 'MIN4TRED') THEN
         PARNAM = 'MIN4TRED'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'B54     ') THEN
               MIN4TRED = 'B54'
            ELSE IF (CHRPARM == 'STC     ') THEN
               MIN4TRED = 'STC'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'B54 or STC',CHRPARM,MIN4TRED
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'B54 or STC',CHRPARM,MIN4TRED
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MPFOUT

      ELSE IF (JCARD(2)(1:8) == 'MPFOUT  ') THEN
         PARNAM = 'MPFOUT  '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == '6       ') THEN
               MPFOUT = '6'
            ELSE IF (CHRPARM == 'R       ') THEN
               MPFOUT = 'R'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'6 or R',CHRPARM,MPFOUT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'6 or R',CHRPARM,MPFOUT
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MXITERI is the max number of iterations in the Inverse Power eigenvalue algorithm

      ELSE IF (JCARD(2)(1:8) == 'MXITERI ') THEN
         PARNAM = 'MXITERI '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM >= 1) THEN
               IF (I4PARM /= MXITERI) THEN
                  WRITE(ERR,1146) PARNAM,MXITERI,I4PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1146) PARNAM,MXITERI,I4PARM
                  ENDIF
                  MXITERI = I4PARM
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'>= 1',I4PARM,MXITERI
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= 1',I4PARM,MXITERI
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! MXITERL is the max number of iterations in the Lanczos eigenvalue algorithm

      ELSE IF (JCARD(2)(1:8) == 'MXITERL ') THEN
         PARNAM = 'MXITERL '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM >= 1) THEN
               MXITERL = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'>= 1',I4PARM,MXITERL
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= 1',I4PARM,MXITERL
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! --- chase_feast_add --- begin !
! LANCMETH is a deprecated alias for selecting the EIGRL extract backend when no EIGRL continuation is present

      ELSE IF (JCARD(2)(1:8) == 'LANCMETH') THEN
         PARNAM = 'LANCMETH'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM(1:6) == 'ARPACK') THEN
               LANCMETH = 'ARPACK  '
            ELSE IF (CHRPARM(1:5) == 'CHASE') THEN
               LANCMETH = 'CHASE   '
            ELSE IF (CHRPARM(1:5) == 'FEAST') THEN
               LANCMETH = 'FEAST   '
            ELSE IF (CHRPARM(1:5) == 'SUBSP') THEN
               LANCMETH = 'SUBSP   '
            ELSE IF (CHRPARM(1:5) == 'DENSE') THEN
               LANCMETH = 'DENSE   '
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'ARPACK/CHASE/FEAST/SUBSP/DENSE',CHRPARM,LANCMETH
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'ARPACK/CHASE/FEAST/SUBSP/DENSE',CHRPARM,LANCMETH
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )
         CALL CRDERR ( CARD )

! OTMSKIP defines whether tp quit if a singularity is found in matrix decomp

      ELSE IF (JCARD(2)(1:8) == 'OTMSKIP ') THEN
         PARNAM = 'OTMSKIP '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM /= OTMSKIP) THEN
               WRITE(ERR,1146) PARNAM,OTMSKIP,I4PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1146) PARNAM,OTMSKIP,I4PARM
               ENDIF
               OTMSKIP = I4PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PBARLDEC is the the number of decimal digits to use when writing the PBAR equivalents for PBARL entry real data

      ELSE IF (JCARD(2)(1:8) == 'PBARLDEC') THEN
         PARNAM = 'PBARLDEC'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF ((I4PARM >= 0) .AND. (I4PARM <= MPBARLU)) THEN
               PBARLDEC = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PBARLDEC
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PBARLDEC
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PBARLSHR indicator of whether to include shear flexibility (k1, k2) in equivalent PBAR for PBARL Bulk Data entries

      ELSE IF (JCARD(2)(1:8) == 'PBARLSHR') THEN
         PARNAM = 'PBARLSHR'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, PBARLSHR)

! PCOMPEQ Indicator to write equiv PSHELL, MAT2 to F06 for PCOMP's

      ELSE IF (JCARD(2)(1:8) == 'PCOMPEQ ') THEN
         PARNAM = 'PCOMPEQ'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 2
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PCOMPEQ = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PCOMPEQ
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PCOMPEQ
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PCMPTSTM is used for specifying PCOMP ratio shear thickness to total plate thickness

      ELSE IF (JCARD(2)(1:8) == 'PCMPTSTM') THEN
         PARNAM = 'PCMPTSTM'
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - PCMPTSTM) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,PCMPTSTM,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,PCMPTSTM,R8PARM
               ENDIF
               PCMPTSTM = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PCHSPC1 decides whether to write text file of SPC1's based on the G.P. singularity processor results.
! SPC1SID is the set ID (if not entered, then default value in module PARAMS will be used).

      ELSE IF (JCARD(2)(1:8) == 'PCHSPC1 ') THEN
         PARNAM = 'PCHSPC1 '
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM(1:1) == 'Y') THEN
               PCHSPC1 = 'Y'
            ELSE IF (CHRPARM(1:1) == 'N') THEN
               PCHSPC1 = 'N'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,PCHSPC1
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,PCHSPC1
               ENDIF
            ENDIF
         ENDIF

         IF (JCARD(4)(1:) /= ' ') THEN
            PARNAM = 'SPC1SID'
            CALL I4FLD ( JCARD(4), JF(4), I4PARM )
            IF (IERRFL(4) == 'N') THEN
               SPC1SID = I4PARM
            ENDIF
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN
            PARNAM = 'SPC1QUIT'
            CALL CHAR_FLD ( JCARD(5), JF(5), CHRPARM )
            IF (IERRFL(5) == 'N') THEN
               CALL LEFT_ADJ_BDFLD ( CHRPARM )
               IF      (CHRPARM(1:1) == 'Y') THEN
                  SPC1QUIT = 'Y'
               ELSE IF (CHRPARM(1:1) == 'N') THEN
                  SPC1QUIT = 'N'
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,SPC1QUIT
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,SPC1QUIT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! POST writes FEMAP data to file NEU
      ELSE IF (JCARD(2)(1:8) == 'POST    ') THEN
         PARNAM = 'POST    '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM <= 0) THEN
               POST = I4PARM
               PRTNEU = 'Y'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'>= 0',I4PARM,POST
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'>= 0',I4PARM,POST
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTBASIC prints the grid coordinates in the basic coord system
      ELSE IF (JCARD(2)(1:8) == 'PRTBASIC') THEN
         PARNAM = 'PRTBASIC'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTBASIC = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTBASIC
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTBASIC
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTCONN prints table of elements connected to each grid

      ELSE IF (JCARD(2)(1:8) == 'PRTCONN ') THEN
         PARNAM = 'PRTCONN '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM > 0) THEN
               PRTCONN = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'> 0',I4PARM,PRTCONN
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'> 0',I4PARM,PRTCONN
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTCGLTM prints the CB matrix

      ELSE IF (JCARD(2)(1:8) == 'PRTCGLTM') THEN
         PARNAM = 'PRTCGLTM'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTCGLTM = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTCGLTM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTCGLTM
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTCORD prints the grid coordinates in the basic coord system

      ELSE IF (JCARD(2)(1:8) == 'PRTCORD ') THEN
         PARNAM = 'PRTCORD '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 2
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTCORD = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTCORD
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTCORD
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTDISP  prints G, N, F, A and/or L set displ in LINK5

      ELSE IF (JCARD(2)(1:8) == 'PRTDISP ') THEN
         PARNAM = 'PRTDISP '
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTDISP(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTDISP(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTDISP(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,0,0,0 )! Make sure that there are no imbedded blanks in field 3-6
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,7,8,9 )! Issue warning if fields 7-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTDLR  prints DLR constraint matrix

      ELSE IF (JCARD(2)(1:8) == 'PRTDLR  ') THEN
         PARNAM = 'PRTDLR  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTDLR = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDLR
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDLR
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTDOF prints the DOF tables (TDOF and/or TDOFI). TDOF is sorted by external grid, TDOFI by resequenced DOF number
      ELSE IF (JCARD(2)(1:8) == 'PRTDOF  ') THEN
         PARNAM = 'PRTDOF  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 3
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTDOF = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTFOR prints G, N, F, A and/or L set stiffness matrices in LINK2
      ELSE IF (JCARD(2)(1:8) == 'PRTFOR  ') THEN
         PARNAM = 'PRTFOR  '
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTFOR(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTFOR(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTFOR(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTGMN  prints GMN constraint matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTGMN  ') THEN
         PARNAM = 'PRTGMN  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTGMN = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTGMN
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTGMN
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! PRTGOA  prints GOA constraint matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTGOA  ') THEN
         PARNAM = 'PRTGOA  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTGOA = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTGOA
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTGOA
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTHMN  prints HMN constraint matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTHMN  ') THEN
         PARNAM = 'PRTHMN  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTHMN = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTHMN
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTHMN
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTIFLTM prints the CB matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTIFLTM') THEN
         PARNAM = 'PRTIFLTM '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTIFLTM = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTIFLTM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTIFLTM
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTKXX prints the CB matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTKXX  ') THEN
         PARNAM = 'PRTKXX  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTKXX = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTKXX
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTKXX
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTMASS  prints G, N, F, A and/or L set stiffness matrices in LINK2
      ELSE IF (JCARD(2)(1:8) == 'PRTMASS ') THEN
         PARNAM = 'PRTMASS '
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTMASS(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTMASS(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTMASS(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTMASSD prints G, N, F, A and/or L set mass matrix diagonals in LINK2
      ELSE IF (JCARD(2)(1:8) == 'PRTMASSD') THEN
         PARNAM = 'PRTMASSD'
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTMASSD(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTMASSD(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTMASSD(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTMXX prints the CB matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTMXX  ') THEN
         PARNAM = 'PRTMXX  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTMXX = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTMXX
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTMXX
               ENDIF
            ENDIF
         ENDIF


      ! PRTOU4 prints the CB matrix
      ELSE IF (JCARD(2)(1:8) == 'PRTOU4  ') THEN
         PARNAM = 'PRTOU4  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTOU4 = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTOU4
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTOU4
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! PRTPHIXA prints the CB matrix PHIXA
      ELSE IF (JCARD(2)(1:8) == 'PRTPHIXA') THEN
         PARNAM = 'PRTPHIXA'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTPHIXA = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTPHIXA
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTPHIXA
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTPHIZL prints the CB matrix

      ELSE IF (JCARD(2)(1:8) == 'PRTPHIZL') THEN
         PARNAM = 'PRTPHIZL'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTPHIZL = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTPHIZL
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTPHIZL
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTQSYS  prints QSYS enforced displs

      ELSE IF (JCARD(2)(1:8) == 'PRTQSYS ') THEN
         PARNAM = 'PRTQSYS  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTQSYS = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTQSYS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTQSYS
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTRMG  prints RMG constraint matrix and it's partitions

      ELSE IF (JCARD(2)(1:8) == 'PRTRMG  ') THEN
         PARNAM = 'PRTRMG  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 3
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTRMG = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTRMG
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTRMG
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTSCP prints subcase information in the subcase processor subroutine

      ELSE IF (JCARD(2)(1:8) == 'PRTSCP  ') THEN
         PARNAM = 'PRTSCP  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTSCP = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTSCP
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTSCP
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTPSET prints the OUTPUT4 matrix partitioning vector sets (PSET).

      ELSE IF (JCARD(2)(1:8) == 'PRTPSET ') THEN
         PARNAM = 'PRTPSET '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 999999999
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTPSET = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTTSET prints the DOF set table (TSET).

      ELSE IF (JCARD(2)(1:8) == 'PRTTSET ') THEN
         PARNAM = 'PRTTSET '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 999999999
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTTSET = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTSTIFD prints G, N, F, A and/or L set stiffness matrix diagonals in LINK2

      ELSE IF (JCARD(2)(1:8) == 'PRTSTIFD') THEN
         PARNAM = 'PRTSTIFD'
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTSTIFD(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTSTIFD(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTSTIFD(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTSTIFF prints G, N, F, A and/or L set stiffness matrices in LINK2

      ELSE IF (JCARD(2)(1:8) =='PRTSTIFF') THEN
         PARNAM = 'PRTSTIFF'
         DO I=1,5
            CALL I4FLD ( JCARD(I+2), JF(I+2), I4PARMI(I) )
            IF (IERRFL(I+2) == 'N') THEN
               LOWER = 0
               UPPER = 3
               IF ((I4PARMI(I) >= LOWER) .AND. (I4PARMI(I) <= UPPER)) THEN
                  PRTSTIFF(I) = I4PARMI(I)
               ELSE
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTSTIFF(I)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(I+2),PRTSTIFF(I)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in field 3-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTUO0  prints UO0 enforced displs

      ELSE IF (JCARD(2)(1:8) == 'PRTUO0  ') THEN
         PARNAM = 'PRTUO0  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTUO0 = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTUO0
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTUO0
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTUSET states whether to print the USET table (DOF's defined on USET Bulk Data entries that belong to a
! user defined set named either "U1" or "U2").

      ELSE IF (JCARD(2)(1:8) == 'PRTUSET ') THEN
         PARNAM = 'PRTUSET '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 999999999
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTUSET = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTDOF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! PRTYS  prints YS enforced displs

      ELSE IF (JCARD(2)(1:8) == 'PRTYS   ') THEN
         PARNAM = 'PRTYS  '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 0
            UPPER = 1
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               PRTYS = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTYS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),PRTYS
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! Q4SURFIT sets the polynomial order for the surface fit of QUAD4 stress/strain when stresses are requested for other than corner

      ELSE IF (JCARD(2)(1:8) == 'Q4SURFIT') THEN
         PARNAM = 'Q4SURFIT'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            LOWER = 1
            UPPER = 6
            IF ((I4PARM >= LOWER) .AND. (I4PARM <= UPPER)) THEN
               Q4SURFIT = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1196) PARNAM,LOWER,UPPER,JCARD(3),Q4SURFIT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1196) PARNAM,LOWER,UPPER,JCARD(3),Q4SURFIT
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! QUADAXIS sets the method for how the local x axis of quad elements is determined

      ELSE IF (JCARD(2)(1:8) == 'QUADAXIS') THEN
         PARNAM = 'QUADAXIS'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'SPLITD  ') THEN
               QUADAXIS = 'SPLITD'
            ELSE IF (CHRPARM == 'SIDE12  ') THEN
               QUADAXIS = 'SIDE12'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'SPLITD or SIDE12',CHRPARM,QUADAXIS
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'SPLITD or SIDE12',CHRPARM,QUADAXIS
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! QUAD4TYP tells which element to use for the MYSTRAN QUAD4 element

      ELSE IF (JCARD(2)(1:8) == 'QUAD4TYP') THEN
         PARNAM = 'QUAD4TYP'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'MIN4T   ') THEN
               QUAD4TYP = 'MIN4T '
            ELSE IF (CHRPARM == 'MIN4    ') THEN
               QUAD4TYP = 'MIN4  '
            ELSE IF (CHRPARM == 'MITC4   ') THEN
               QUAD4TYP = 'MITC4 '
            ELSE IF (CHRPARM == 'MITC4+  ') THEN
               QUAD4TYP = 'MITC4+'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'MIN4T, MIN4, MITC4 or MITC4+',CHRPARM,QUAD4TYP
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'MIN4T, MIN4, MITC4 or MITC4+',CHRPARM,QUAD4TYP
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! RCONDK = 'Y' executes LAPACK code in LINK3 to calc the recriprocal of the condition number, RCOND, of a matrix to be decomposed
      ELSE IF (JCARD(2)(1:8) == 'RCONDK  ') THEN
         PARNAM = 'RCONDK  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, RCONDK)

      ! RELINK3 = 'Y' causes LINK3 (and therefore also LINK5) to be rerun in a RESTART. This is only used in linear statics
      ELSE IF (JCARD(2)(1:8) == 'RELINK3 ') THEN
         PARNAM = 'RELINK3  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, RELINK3)

! SETLKTK parameter gives the option number for how LTERM_KGG is calculated

      ELSE IF (JCARD(2)(1:8) == 'SETLKTK ') THEN
         PARNAM = 'SETLKTK '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF ((I4PARM == 0) .OR. (I4PARM == 1) .OR. (I4PARM == 2) .OR. (I4PARM == 3) .OR. (I4PARM == 4)) THEN
               SETLKTK = I4PARM
            ELSE                                           ! If not 0,1,2,3,4 leave SETLKTK = default value
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'0,1,2,3 OR 4',I4PARM,SETLKTM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'0,1,2,3 OR 4',I4PARM,SETLKTM
               ENDIF
            ENDIF                                          ! If 0,1,2,3 see if user wants to pause to input value at execution time
            IF (SETLKTK /= 4) THEN
               CHRPARM(1:) = ' '
               PARNAM = 'ESP0_PAUSE'
               CALL CHAR_FLD ( JCARD(4), JF(4), CHRPARM )
               IF (IERRFL(4) == 'N') THEN
                  IF (CHRPARM(1:) /= ' ') THEN
                     CALL LEFT_ADJ_BDFLD ( CHRPARM )
                     IF (CHRPARM == 'PAUSE   ') THEN
                        ESP0_PAUSE = 'Y'
                     ELSE
                        WARN_ERR = WARN_ERR + 1
                        WRITE(ERR,1178) CHRPARM,PARNAM,JF(4),'ESP0'
                        IF (SUPWARN == 'N') THEN
                           WRITE(F06,1178) CHRPARM,PARNAM,JF(4),'ESP0'
                        ENDIF
                        ESP0_PAUSE = 'Y'
                     ENDIF
                  ELSE
                     ESP0_PAUSE = 'N'
                  ENDIF
               ENDIF
            ELSE                                           ! Read user value for LTERM_MGGE from field 5 when SETLKTM = 4
               I4PARM = 0
               PARNAM = 'USR_LTERM_KGG'
               CALL I4FLD ( JCARD(5), JF(5), I4PARM )
               IF (IERRFL(5) == 'N') THEN
                  IF (I4PARM > 0) THEN
                     USR_LTERM_KGG = I4PARM
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1179) I4PARM,PARNAM,JF(5)
                     WRITE(F06,1179) I4PARM,PARNAM,JF(5)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! SETLKTM parameter gives the option number for how LTERM_MGG is calculated

      ELSE IF (JCARD(2)(1:8) == 'SETLKTM ') THEN
         PARNAM = 'SETLKTM '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )            ! Field 3 has SETLKTM value
         IF (IERRFL(3) == 'N') THEN
            IF ((I4PARM == 0) .OR. (I4PARM == 3) .OR. (I4PARM == 4)) THEN
               SETLKTM = I4PARM
            ELSE                                           ! If not 0,3,4 leave SETLKTM = default value
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'0, 3 OR 4',I4PARM,SETLKTM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'0, 3 OR 4',I4PARM,SETLKTM
               ENDIF
            ENDIF
            IF (SETLKTM /= 4) THEN                         ! Field 4 has ESP0_PAUSE (only used if SETLKTM /= 4)
               CHRPARM(1:) = ' '
               PARNAM = 'EMP0_PAUSE'
               CALL CHAR_FLD ( JCARD(4), JF(4), CHRPARM )
               IF (IERRFL(4) == 'N') THEN
                  IF (CHRPARM(1:) /= ' ') THEN
                     CALL LEFT_ADJ_BDFLD ( CHRPARM )
                     IF (CHRPARM == 'PAUSE   ') THEN
                        EMP0_PAUSE = 'Y'
                     ELSE
                        WARN_ERR = WARN_ERR + 1
                        WRITE(ERR,1178) CHRPARM,PARNAM,JF(4),'EMP0'
                        IF (SUPWARN == 'N') THEN
                           WRITE(F06,1178) CHRPARM,PARNAM,JF(4),'EMP0'
                        ENDIF
                        EMP0_PAUSE = 'Y'
                     ENDIF
                  ELSE
                     EMP0_PAUSE = 'N'
                  ENDIF
               ENDIF
            ELSE                                           ! Read user value for LTERM_MGG from field 5 when SETLKTK = 4
               I4PARM = 0
               PARNAM = 'USR_LTERM_MGG'
               CALL I4FLD ( JCARD(5), JF(5), I4PARM )
               IF (IERRFL(5) == 'N') THEN
                  IF (I4PARM > 0) THEN
                     USR_LTERM_MGG = I4PARM
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1179) I4PARM,PARNAM,JF(5)
                     WRITE(F06,1179) I4PARM,PARNAM,JF(5)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,4,5,0,0,0,0 )! Make sure that there are no imbedded blanks in fields 3-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! SHRFXFAC Factor used to adjust transverse shear stiffness when user has indicated zero shear flexibility for shell elements.
! The shear stiffness will be reset from infinite (zero flexibility) to SHRFXFAC times the average of the bending stiffnesses
! in the 2 planes

      ELSE IF (JCARD(2)(1:8) == 'SHRFXFAC') THEN
         PARNAM = 'SHRFXFAC '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (DABS(R8PARM - SHRFXFAC) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,SHRFXFAC,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,SHRFXFAC,R8PARM
               ENDIF
               SHRFXFAC = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! SKIPMGG 'Y', 'N' indicator to say whether to skip calculation of MGG
      ELSE IF (JCARD(2)(1:7) == 'SKIPMGG') THEN
         PARNAM = 'SKIPMGG '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, SKIPMGG)

! SOLLIB sets the method for solving equations (BANDED, SPARSE)

      ELSE IF (JCARD(2)(1:8) == 'SOLLIB  ') THEN
         PARNAM = 'SOLLIB'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF (CHRPARM == 'BANDED  ') THEN
               SOLLIB = 'BANDED  '
            ELSE IF (CHRPARM == 'SPARSE  ') THEN
               SOLLIB = 'SPARSE  '
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'BANDED or SPARSE',CHRPARM,SOLLIB
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'BANDED or SPARSE',CHRPARM,SOLLIB
               ENDIF
            ENDIF
         ENDIF

         IF (SOLLIB == 'SPARSE  ') THEN
            IF (JCARD(4)(1:) /= ' ') THEN
               PARNAM = 'SPARSE_FLAVOR'
               CALL CHAR_FLD ( JCARD(4), JF(4), CHRPARM )
               IF (IERRFL(4) == 'N') THEN
                  CALL LEFT_ADJ_BDFLD ( CHRPARM )
                  IF      (CHRPARM(1:7) == 'SUPERLU') THEN
                     SPARSE_FLAVOR = 'SUPERLU '
                  ELSE
                     WARN_ERR = WARN_ERR + 1
                     WRITE(ERR,101) CARD
                     WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_INFO
                     IF (SUPWARN == 'N') THEN
                        IF (ECHO == 'NONE  ') THEN
                           WRITE(F06,101) CARD
                        ENDIF
                        WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,AUTOSPC_INFO
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! SORT_MAX sets the max number of times SORT algorithms are executed

      ELSE IF (JCARD(2)(1:8) == 'SORT_MAX') THEN
         PARNAM = 'SORT_MAX'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM > 0) THEN
               SORT_MAX = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'> 0',I4PARM,SORT_MAX
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'> 0',I4PARM,SORT_MAX
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! SPARSTOR sets whether compressed storage matrices are stored as symmetric (only nonzeros on and above the diagonal are stored)
! or are stored as nonsymmetric (all nonzeros are stored). This affects ALL compressed storage matrices their matrix operations

      ELSE IF (JCARD(2)(1:8) == 'SPARSTOR') THEN
         PARNAM = 'SPARSTOR'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (IERRFL(3) == 'N') THEN
            CALL LEFT_ADJ_BDFLD ( CHRPARM )
            IF      (CHRPARM == 'SYM     ') THEN
               SPARSTOR = 'SYM   '
            ELSE IF (CHRPARM == 'NONSYM  ') THEN
               SPARSTOR = 'NONSYM'
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1189) PARNAM,'SYM OR NONSYM',CHRPARM,SPARSTOR
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1189) PARNAM,'SYM OR NONSYM',CHRPARM,SPARSTOR
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! STR_CID defines output coord sys for elem stress/strain/engr force (if it is not to be elem local sys).
! STR_CID can be 0 (basic coord system) or an actual coord system ID (i.e. one defined on a Bulk data CORDij entry)

      ELSE IF (JCARD(2)(1:8) == 'STR_CID ') THEN
         PARNAM = 'STR_CID'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            STR_CID = I4PARM
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ! SUPINFO determines whether warning messages are suppressed in the output file
      ELSE IF (JCARD(2)(1:8) == 'SUPINFO ') THEN
         PARNAM = 'SUPINFO  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, SUPINFO)


      ! SUPWARN determines whether warning messages are suppressed in the output file
      ELSE IF (JCARD(2)(1:8) == 'SUPWARN ') THEN
         PARNAM = 'SUPWARN  '
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, SUPWARN)

      ! NOCOUNT suppresses all "counter", progress-indicating, non-advance writes to
      ! standard output. It's good when debugging, and basically required if stdout
      ! has no rewrite capabilities (like a file). We should probably check if the
      ! output is not rewrite-friendly and set it automatically at some point.
      ELSE IF (JCARD(2)(1:8) == 'NOCOUNTS') THEN
         PARNAM = 'NOCOUNTS'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, NOCOUNTS)

      ! THRESHK is used in determining if equilibration is performed (see LAPACK subr DLAQSB in MODULE LAPACK_BLAS_AUX_1)
      ELSE IF (JCARD(2)(1:8) == 'THRESHK ') THEN
         PARNAM = 'THRESHK '
         CALL R8FLD ( JCARD(3), JF(3), THRESHK )
         IF (IERRFL(3) == 'N') THEN
            IF (DABS(THRESHK - THRESHK_LAP) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,THRESHK_LAP,THRESHK
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,THRESHK_LAP,THRESHK
               ENDIF
            ENDIF
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! TINY acts as a filter for small terms in matrix print
      ELSE IF (JCARD(2)(1:8) == 'TINY    ') THEN
         PARNAM = 'TINY'
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         IF (CHRPARM == 'MACHPREC') THEN
            TINY = MACH_PREC
         ELSE
            CALL R8FLD ( JCARD(3), JF(3), R8PARM )
            IF (IERRFL(3) == 'N') THEN
               IF (DABS(R8PARM - TINY) > EPS1) THEN
                  WRITE(ERR,1147) PARNAM,TINY,R8PARM
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,1147) PARNAM,TINY,R8PARM
                  ENDIF
                  TINY = R8PARM
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! TSTM_DEF changes default value of TS/TM (ratio of shear to membrane thickness) for PSHELL's

      ELSE IF (JCARD(2)(1:8) == 'TSTM_DEF') THEN
         PARNAM = 'TSTM_DEF'
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (DABS(R8PARM - TSTM_DEF) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,TSTM_DEF,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,TSTM_DEF,R8PARM
               ENDIF
               TSTM_DEF = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! USETPRT is not currently a MYSTRAN PARAM but MYSTRAN will automatically respond with USETPRT = 0 if there are any USETSTR entries
! USETSEL is not currently a MYSTRAN PARAM but MYSTRAN will automatically respond with USETSEL = 0 if there are any USETSTR entries

      ELSE IF (JCARD(2)(1:8) == 'USETPRT ') THEN
         WRITE(ERR,102) JCARD(2), JCARD(2), '0'
         WRITE(F06,102) JCARD(2), JCARD(2), '0'

      ELSE IF (JCARD(2)(1:8) == 'USETSEL ') THEN
         WRITE(ERR,102) JCARD(2), JCARD(2), ' 0 OR -1'
         WRITE(F06,102) JCARD(2), JCARD(2), ' 0 OR -1'

! USETSTR defines DOF sets to be written out in the F06 file with grid/comp and DOF number

      ELSE IF (JCARD(2)(1:8) == 'USETSTR ') THEN
         CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
         CALL LEFT_ADJ_BDFLD ( CHRPARM )
do_i:    DO I=1,JCARD_LEN
            IF (CHRPARM(I:I) == ':') THEN
               IERR = IERR + 1
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,104) CHRPARM
               WRITE(F06,104) CHRPARM
               EXIT do_i
            ENDIF
         ENDDO do_i
         IF (IERR == 0) THEN
            CALL LOAD_USETSTR_TABLE ( CHRPARM(1:2), JCARD(2) )
         ENDIF

! USR_JCT is a user supplied value for JCT - used in sort routines

      ELSE IF (JCARD(2)(1:8) == 'USR_JCT ') THEN
         PARNAM = 'USR_JCT '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM)
         IF (IERRFL(3) == 'N') THEN
            IF (I4PARM > 0) THEN
               USR_JCT = I4PARM
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1172) PARNAM,'> 0',I4PARM,USR_JCT
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1172) PARNAM,'> 0',I4PARM,USR_JCT
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! --- BANDED_optimizisation -begin-- !
! WINAMEM is an optional per-array MB cap. Historical MYSTRAN used this for the Windows XP 2 GB address-space limit.
! In modern 64-bit builds, WINAMEM <= 0 disables this legacy cap and lets ALLOCATE/STAT report memory failure.
! --- BANDED_optimizisation -end-- !

      ELSE IF (JCARD(2)(1:8) == 'WINAMEM ') THEN
         PARNAM = 'WINAMEM  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(2) == 'N') THEN
            IF (DABS(R8PARM - WINAMEM) > EPS1) THEN
               WRITE(ERR,1147) PARNAM,WINAMEM,R8PARM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,1147) PARNAM,WINAMEM,R8PARM
               ENDIF
               WINAMEM = R8PARM
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! WTMASS multiplies the mass matrix (after the grid pt wgt generator) by the real WTMASS value

      ELSE IF (JCARD(2)(1:8) == 'WTMASS  ') THEN
         PARNAM = 'WTMASS  '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (R8PARM > ZERO) THEN
               WTMASS = R8PARM
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1173) PARNAM,'> 0.D0',R8PARM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1173) PARNAM,'> 0.D0',R8PARM
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


! K6ROT is a small stiffness added to the drilling dof of 5-dof shells.

      ELSE IF (JCARD(2)(1:8) == 'K6ROT   ') THEN
         PARNAM = 'K6ROT   '
         CALL R8FLD ( JCARD(3), JF(3), R8PARM )
         IF (IERRFL(3) == 'N') THEN
            IF (R8PARM >= 0) THEN
               K6ROT = R8PARM
               IF (K6ROT > 0) THEN
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,101) CARD
                  WRITE(ERR,1198)
                  IF (SUPWARN == 'N') THEN
                     IF (ECHO == 'NONE  ') THEN
                        WRITE(F06,101) CARD
                     ENDIF
                     WRITE(F06,1198)
                  ENDIF
               ENDIF
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1173) PARNAM,'>= 0.D0',R8PARM
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1173) PARNAM,'>= 0.D0',R8PARM
               ENDIF
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         !WRITE(*,*) 'K6ROT = ', K6ROT

! This guard is necessary since the SLU_SET_* subroutines are only present when
! linking against SuperLU_MT. That's why this file is .F90 and not .f90. Sorry!
#ifdef USE_SUPERLU_MT
      ! SuperLU_MT memory growth factors
      ELSE IF (JCARD(2)(1:7) == 'SPIENV6') THEN
         PARNAM = 'SPIENV6 '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            SPIENV6 = I4PARM
            CALL SLU_SET_MEMFAC(6, SPIENV6)
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! SuperLU_MT memory growth factors
      ELSE IF (JCARD(2)(1:7) == 'SPIENV7') THEN
         PARNAM = 'SPIENV7 '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            SPIENV7 = I4PARM
            CALL SLU_SET_MEMFAC(7, SPIENV7)
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! SuperLU_MT memory growth factors
      ELSE IF (JCARD(2)(1:7) == 'SPIENV8') THEN
         PARNAM = 'SPIENV8 '
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            SPIENV8 = I4PARM
            CALL SLU_SET_MEMFAC(8, SPIENV8)
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ! SuperLU_MT number of threads
      ELSE IF (JCARD(2)(1:8) == 'SLU_NTHR') THEN
         PARNAM = 'SLU_NTHR'
         CALL I4FLD ( JCARD(3), JF(3), I4PARM )
         IF (IERRFL(3) == 'N') THEN
            SLU_NTHR = I4PARM
            CALL SLU_SET_NTHR(SLU_NTHR)
         ENDIF
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
#endif

! PARAM parameter name not recognized

      ELSE

         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,101) CARD
         WRITE(ERR,1171) JCARD(2)
         IF (SUPWARN == 'N') THEN
            IF (ECHO == 'NONE  ') THEN
               WRITE(F06,101) CARD
            ENDIF
            WRITE(F06,1171) JCARD(2)
         ENDIF
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF
! --- chase_feast_add --- end !



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  102 FORMAT(' *INFORMATION: PARAMETER ',A8,' IS NOT USED IN MYSTRAN. PRINTOUT OF USETSTR REQUESTS WILL BE HONORED',&
                           ' AS IF MSC NASTRAN ',A8,' = ',A)

  104 FORMAT(' *WARNING    : USE OF ":" SEPERATOR FOR THE PARAM USETSTR BULK DATA ENTRY NOT ALLOWED. ENTRY HAD "',A,'" IN FIELD 3')

  106 FORMAT(" *WARNING    : PARAMETER K6ROT SHOULD NOT BE SET TO ZERO SAVE FOR DEBUGGING.")

 1110 FORMAT(' *ERROR  1110: PARAMETER NAMED ',A,' MUST BE >= ',I2,' AND <= ',I2,' BUT INPUT VALUE IS: ',A)

 1146 FORMAT(' *INFORMATION: PARAMETER ',A,' CHANGED FROM ',I8,' TO ',I8,/)

 1147 FORMAT(' *INFORMATION: PARAMETER ',A,' CHANGED FROM ',1ES13.6,' TO ',1ES13.6,/)

 1148 FORMAT(' *INFORMATION: PARAMETER ',A,I3,' CHANGED FROM ',1ES13.6,' TO ',1ES13.6,/)

 1149 FORMAT(' *WARNING    : ONLY EPSIL( 1) THRU EPSIL(',I2,') CAN BE CHANGED ON PARAM ENTRY')

 1171 FORMAT(' *WARNING    : PARAMETER NAMED ',A,' NOT RECOGNIZED. ENTRY IGNORED')

 1172 FORMAT(' *WARNING    : PARAMETER NAMED ',A,' MUST BE ',A,' BUT INPUT VALUE IS: ',I8,'. DEFAULT VALUE ',I8,' WILL BE USED')

 1173 FORMAT(' *ERROR  1173: PARAMETER NAMED ',A,' MUST BE ',A,' BUT INPUT VALUE IS ',1ES13.6)

 1178 FORMAT(' *WARNING    : VALUE = ',A,' FOR PARAMETER ',A,' INPUT IN FIELD ',I3,' IS INCORRECT. SHOULD BE "PAUSE".'             &
                    ,/,14X,' MYSTRAN WLL BE PAUSE''d AFTER SUBROUTINE ',A,' FOR USER INPUT OF LTERM_KGG')

 1179 FORMAT(' *ERROR  1179: VALUE = ',I8,' FOR DEBUG PARAMETER ',A,' INPUT IN FIELD ',I3,' MUST BE >= 0')

 1189 FORMAT(' *WARNING    : PARAMETER NAMED ',A,' MUST BE ',A,' BUT INPUT VALUE IS: ',A,'. DEFAULT VALUE = ',A,' WILL BE USED')

 1196 FORMAT(' *WARNING    : PARAMETER NAMED ',A,' MUST BE >= ',I2,' AND <= ',I2,' BUT INPUT VALUE IS: ',A,'. DEFAULT VALUE = ',   &
                             I8,' WILL BE USED')

 1197 FORMAT(' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
             '>>>>>>>>>>>>>>>>>',/,                                                                                             &
             ' *WARNING    : PARAMETER ',A15,' SHOULD NOT BE CHANGED FROM "SYM" TO "NONSYM". WRONG ANSWERS WILL PROBABLY RESULT',/,&
             ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
             '>>>>>>>>>>>>>>>>>')

 1198 FORMAT("*WARNING    : THE K6ROT PARAMETER IS AN EXPERIMENTAL FEATURE. WHEN IN DOUBT, SET IT TO ZERO.",/)
! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, VALUE)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)         :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN), INTENT(IN) :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(15*BYTE), INTENT(IN)       :: PARNAM            ! The name  of the parameter
      !
      CHARACTER(LEN(JCARD)), INTENT(INOUT) :: CHRPARM           ! A char value read from a PARAM Bulk Data card
      CHARACTER(1*BYTE), INTENT(INOUT)     :: VALUE             ! The vakue of the parameter

      CALL CHAR_FLD ( JCARD(3), JF(3), CHRPARM )
      IF (IERRFL(3) == 'N') THEN
         CALL LEFT_ADJ_BDFLD ( CHRPARM )
         IF      (CHRPARM(1:1) == 'Y') THEN
            VALUE = 'Y'
         ELSE IF (CHRPARM(1:1) == 'N') THEN
            VALUE = 'N'
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,1189) PARNAM,'Y OR N',CHRPARM,VALUE
            IF (VALUE == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,101) CARD
               ENDIF
               WRITE(F06,1189) PARNAM,'Y OR N',CHRPARM,VALUE
            ENDIF
         ENDIF
      ENDIF
      CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,0,0,0,0,0 )! Make sure that there are no imbedded blanks in field 3
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-9 not blank
      CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      RETURN
  101 FORMAT(A)
 1189 FORMAT(' *WARNING    : PARAMETER NAMED ',A,' MUST BE ',A,' BUT INPUT VALUE IS: ',A,'. DEFAULT VALUE = ',A,' WILL BE USED')
      END SUBROUTINE YES_NO_CHECK

! ##################################################################################################################################
      SUBROUTINE LOAD_USETSTR_TABLE ( DOF_SET_NAME, PARNAM )
     ! Formulates a table that shows which of the DOF sets (G, M, N, SA, SB, SG, SZ, SE, S, F, O, S, R, L, U1, U2) that have been
     ! requested for output. The table has 2 columns. Col 1 has the DOF set names and col 2 has a 1 if that set was requested for
     ! output or 0 otherwise.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE DOF_TABLES, ONLY            :  USETSTR_TABLE

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: DOF_SET_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: PARNAM
      CHARACTER( 1*BYTE)              :: DOF_SET_NAME_FOUND

      INTEGER(LONG)                   :: I

! **********************************************************************************************************************************
      DOF_SET_NAME_FOUND = 'N'
      DO I=1,16
         IF (DOF_SET_NAME(1:2) == USETSTR_TABLE(I,1)(1:2)) THEN
            USETSTR_TABLE(I,2) = '1 '
            DOF_SET_NAME_FOUND = 'Y'
         ENDIF
      ENDDO

      IF (DOF_SET_NAME_FOUND == 'Y') THEN
         NUM_USETSTR = NUM_USETSTR + 1
      ELSE
         WRITE(ERR,101) DOF_SET_NAME, PARNAM, (USETSTR_TABLE(I,1),I=1,16)
         WRITE(F06,101) DOF_SET_NAME, PARNAM, (USETSTR_TABLE(I,1),I=1,16)
         RETURN
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *WARNING    : SNAME = "',A,'" ON BULK DATA PARAM "',A,'" NOT RECOGNIZED.'                                           &
                    ,/,14x,' ALLOWABLE SNAME VALUES ARE: ',16(A2,2X))

! **********************************************************************************************************************************

      END SUBROUTINE LOAD_USETSTR_TABLE

      END SUBROUTINE BD_PARAM
