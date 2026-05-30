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

      SUBROUTINE SEQ_PROC

! Generates the grid point sequence order.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     SEQ,     L1B
      USE IOUNT1, ONLY                :  WRT_ERR, SEQFIL
      USE IOUNT1, ONLY                :  WRT_ERR, SEQSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, FATAL_ERR, NELE, NGRID, NSEQ, PROG_NAME, WARN_ERR
      USE PARAMS, ONLY                :  BANDEDOPT, EPSIL, GRIDSEQ, SOLLIB
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  ETYPE, GRID_ID, GRID_SEQ, INV_GRID_SEQ, SEQ1, SEQ2
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SEQ_PROC_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SEQ_PROC'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: NSEQ_SAVE         ! Save of incoming SEQ count (for diagnostics)
      INTEGER(LONG), ALLOCATABLE      :: TMP_GRID_ID(:)    ! Set to array GRID_ID for aid in sorting GRID_SEQ
      INTEGER(LONG), ALLOCATABLE      :: TMP_GRD_SEQ(:)    ! Set to array GRID_SEQ so we can sort it and get array INV_GRID_SEQ
!                                                            without disturbing GRID_SEQ sequence


      REAL(DOUBLE), ALLOCATABLE       :: R_GSEQ(:)         ! Real sequence numbers (since SEQGP cards can have real no's). In the
!                                                            end, the sequence array that will be used is integer array GRID_SEQ

      INTRINSIC                       :: DBLE
      LOGICAL                         :: DO_RCM            ! Activate in-core RCM sequencing add-on
      LOGICAL                         :: RCM_OK            ! True if RCM graph build/order succeeded



! **********************************************************************************************************************************
! Coming in to this subr, GRID_SEQ is in the order of the grids as read in the input data deck.

      ALLOCATE ( R_GSEQ(NGRID), TMP_GRID_ID(NGRID), TMP_GRD_SEQ(NGRID) )
      NSEQ_SAVE = NSEQ
      DO_RCM = .FALSE.
      RCM_OK = .FALSE.

      IF ((GRIDSEQ(1:3) == 'RCM') .OR. ((BANDEDOPT == 'Y') .AND. (SOLLIB == 'BANDED  ') .AND. (GRIDSEQ(1:6) /= 'BANDIT'))) THEN
         DO_RCM = .TRUE.
      ENDIF

! Generate initial R_GSEQ based on the GRID_SEQ value. R_GSEQ(I) is  the (real) sequence number for Grid Point GRID_ID(I).
! If there are no SEQGP sequencing cards, then this will be the final grid point sequence order (as a real number).
! It will be converted to an array of consecutive integers in GRID_SEQ later.

! NOTE: when Bandit sequencing is requested we use Bandit to generate a set of SEQGP cards (if it runs correctly and finds
!       resequencing is necessary) and then proceed with that complete set of SEQGP cards. Need to do it this way because we have to
!       get INV_GRID_SEQ (later) and to write the seq arrays to L1B.

      IF (DO_RCM) THEN
         CALL RCM_SEQ_PROC ( R_GSEQ, RCM_OK )
         IF (RCM_OK) THEN
            NSEQ = 0
            WRITE(ERR,102) NGRID, NSEQ_SAVE, GRIDSEQ, BANDEDOPT, SOLLIB
            IF (SUPINFO == 'N') THEN
               WRITE(F06,102) NGRID, NSEQ_SAVE, GRIDSEQ, BANDEDOPT, SOLLIB
            ENDIF
         ELSE
            WRITE(ERR,103) GRIDSEQ, BANDEDOPT, SOLLIB
            IF (SUPINFO == 'N') THEN
               WRITE(F06,103) GRIDSEQ, BANDEDOPT, SOLLIB
            ENDIF
         ENDIF
      ENDIF

      IF (.NOT. RCM_OK) THEN
         IF       (GRIDSEQ(1:6) == 'BANDIT') THEN             ! Call subr AUTO_SEQ_PROC to generate SEQ1, SEQ2 from SEQGP card images
            CALL AUTO_SEQ_PROC
            IF (NSEQ == NGRID) THEN                           ! Bandit did reseq grids. Set R_GSEQ to the SEQ2 from Bandit SEQGP cards
               DO I=1,NGRID
!                 R_GSEQ(I) = DBLE(SEQ2(I))                   ! Shouldn't need this, SEQ2 is REAL(DOUBLE)
                  R_GSEQ(I) = SEQ2(I)
               ENDDO
            ELSE                                              ! AUTO_SEQ_PROC didn't reseq all grids so reset GRIDSEQ = 'GRID'
               WRITE(ERR,101) NGRID,NSEQ,PROG_NAME
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,101) NGRID,NSEQ,PROG_NAME
               ENDIF
               GRIDSEQ = 'GRID    '                           ! Need this to cover case where AUTO_SEQ_PROC returned without completing
            ENDIF
         ENDIF

         IF (GRIDSEQ(1:4) == 'GRID'  ) THEN                   ! Sequence grids in numerical order, but include SEQGP entries (later)
            DO I=1,NGRID
               R_GSEQ(I) = DBLE(I)
            ENDDO
         ELSE IF ((GRIDSEQ(1:5) == 'INPUT') .OR. (GRIDSEQ(1:3) == 'RCM'))  THEN
            DO I=1,NGRID
              R_GSEQ(I) = DBLE(GRID_SEQ(I))
            ENDDO
         ENDIF
      ENDIF

! Check to make sure that all grid points on SEQGP cards are defined

      IERROR = 0
      DO I = 1,NSEQ
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SEQ1(I), IGRID )
         IF (IGRID == -1) THEN
            WRITE(ERR,1361) 'GRID', SEQ1(I), 'SEQGP BULK DATA ENTRY'
            WRITE(F06,1361) 'GRID', SEQ1(I), 'SEQGP BULK DATA ENTRY'
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO
      IF (IERROR > 0 ) THEN
         WRITE(ERR,9999) SUBR_NAME, IERROR
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9999) SUBR_NAME, IERROR
         ENDIF
         CALL OUTA_HERE ( 'Y' )                            ! Some grid ID's on SEQGP cards not defined
      ENDIF

! Reset GRID_SEQ to consecutive integers. R_GSEQ will be the grid sequence number for the time being. At the end,
! the integer sequence numbers will be back in GRID_SEQ

      DO I=1,NGRID
         GRID_SEQ(I) = I
      ENDDO

! If there are SEQGP cards in the data deck, or if auto sequencing has produced them, then arrays SEQ1, SEQ2 exist.
! Sort SEQ1 (grid numbers) with SEQ2 (sequence numbers) so that SEQ1 is in grid point numerical order (like GRID_ID) and then
! check to make sure that no sequence numbers are duplicated.

      IF (NSEQ > 1) THEN

         CALL SORT_INT1_REAL1 ( SUBR_NAME, 'SEQ1, SEQ2', NSEQ, SEQ1, SEQ2 )

         DO I = 1,NSEQ-1
            IF (SEQ1(I) == SEQ1(I+1)) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1399) SEQ1(I+1),SEQ2(I+1)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1399) SEQ1(I+1),SEQ2(I+1)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

! Change sequence numbers in R_GSEQ to values in SEQ2 (from SEQGP Bulk Data cards). Note: SEQ1 is in grid numerical order
! so R_GSEQ(1) will be for grid GRID_ID(1) and so on.

      IF (NSEQ == NGRID) THEN
         DO I=1,NGRID
            R_GSEQ(I) = SEQ2(I)
         ENDDO
      ELSE
         DO I=1,NGRID
            DO J=1,NSEQ
               IF (SEQ1(J) == GRID_ID(I)) THEN
                  R_GSEQ(I) = SEQ2(J)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! Set TMP_GRID_ID to GRID_ID and then sort TMP_GRID_ID with R_GSEQ so that R_GSEQ is in numerically increasing order.

      DO I=1,NGRID
         TMP_GRID_ID(I) = GRID_ID(I)
      ENDDO


      CALL SORT_REAL1_INT1 ( SUBR_NAME, 'R_GSEQ, TMP_GRID_ID', NGRID, R_GSEQ, TMP_GRID_ID )

! Check to make sure that there are no redundant R_GSEQ sequence numbers


      IERROR = 0
      DO I=1,NGRID-1
         IF ((DABS(R_GSEQ(I+1)) - DABS(R_GSEQ(I))) < EPSIL(1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1360) TMP_GRID_ID(I), TMP_GRID_ID(I+1), R_GSEQ(I)
            WRITE(F06,1360) TMP_GRID_ID(I), TMP_GRID_ID(I+1), R_GSEQ(I)
         ENDIF
      ENDDO

      IF (IERROR > 0) THEN
         WRITE(ERR,9999) SUBR_NAME, IERROR
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9999) SUBR_NAME, IERROR
         ENDIF
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Sort TMP_GRID_ID with R_GSEQ and with GRID_SEQ so that TMP_GRID_ID is in numerically increasing order.

      CALL SORT_INT2_REAL1 ( SUBR_NAME, 'TMP_GRID_ID, GRID_SEQ, R_GSEQ', NGRID, TMP_GRID_ID, GRID_SEQ, R_GSEQ )

! Generate array INV_GRID_SEQ

      DO I=1,NGRID
         TMP_GRD_SEQ(I)  = GRID_SEQ(I)
         INV_GRID_SEQ(I) = I
      ENDDO

      CALL SORT_INT2 ( SUBR_NAME, 'TMP_GRD_SEQ, INV_GRID_SEQ', NGRID, TMP_GRD_SEQ, INV_GRID_SEQ )

! Print table showing R_GSEQ, GRID_SEQ and INV_GRID_SEQ

      IF (DEBUG(13) == 1) THEN
         WRITE(F06,111)
         DO I=1,NGRID
            WRITE(F06,112) GRID_ID(I), I, R_GSEQ(I), GRID_SEQ(I), INV_GRID_SEQ(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,113)
         WRITE(F06,*)
      ENDIF

! **********************************************************************************************************************************
! Now write sequence and GRID_SEQ data to file LINK1B.

      DATA_SET_NAME = 'GRID_SEQ, INV_GRID_SEQ'
      WRITE(L1B) DATA_SET_NAME
      WRITE(L1B) NGRID
      DO I=1,NGRID
         WRITE(L1B) GRID_SEQ(I), INV_GRID_SEQ(I)
      ENDDO

      DATA_SET_NAME = 'SEQ1, SEQ2'
      WRITE(L1B) DATA_SET_NAME
      WRITE(L1B) NSEQ
      DO I=1,NSEQ
         WRITE(L1B) SEQ1(I),SEQ2(I)
      ENDDO



      IF (ALLOCATED(R_GSEQ)) THEN
         DEALLOCATE ( R_GSEQ )
      ENDIF
      IF (ALLOCATED(TMP_GRID_ID)) THEN
         DEALLOCATE ( TMP_GRID_ID )
      ENDIF
      IF (ALLOCATED(TMP_GRD_SEQ)) THEN
         DEALLOCATE ( TMP_GRD_SEQ )
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: SUBR AUTO_SEQ_PROC DID NOT SEQUENCE ALL OF THE ',I8,' GRIDS. ONLY ',I8,' GRIDS WERE SEQUENCED.'       &
                  ,/,15X,A,' WILL DEFAULT TO A SEQUENCE THAT IS IN GRID NUMERICAL ORDER',/)
! --- BANDED_optimizisation -begin-- !
  102 FORMAT(' *INFORMATION: IN-CORE RCM GRID SEQUENCING APPLIED TO ',I8,' GRIDS (IGNORED ',I8,' INPUT/AUTO SEQGP ENTRY(IES)).'    &
                  ,/,15X,' PARAM GRIDSEQ=',A8,' BANDEDOPT=',A1,' SOLLIB=',A8                                                    &
                  ,/,15X,' RCM UPDATES GRID_SEQ/INV_GRID_SEQ BEFORE DOF NUMBERING, SO MASS AND STIFFNESS MATRICES SHARE ORDER.')
! --- BANDED_optimizisation -end-- !
  103 FORMAT(' *WARNING    : IN-CORE RCM SEQUENCING REQUESTED BUT A VALID CONNECTIVITY GRAPH COULD NOT BE BUILT.'                  &
                  ,/,15X,' FALLING BACK TO EXISTING GRIDSEQ FLOW. PARAM GRIDSEQ=',A8,' BANDEDOPT=',A1,' SOLLIB=',A8)

  111 FORMAT(56X,'GRID SEQUENCE DATA',//16X,'GRID ID                  I                   R_GSEQ(I)               GRID_SEQ(I)   ', &
                 '     INV_GRID_SEQ(I)',/,12X,'(Actual grid ID)    (Internal grid ID)     (Grid seq - real num)',                  &
                 '  (Grid seq - integer num)  (* - see below)',/)

  112 FORMAT(3X,2I20,1ES28.6,2I20)

  113 FORMAT(15X,'* INV_GRID_SEQ(I) = internal grid ID that is sequenced I-th')

 1399 FORMAT(' *WARNING    : REDUNDANT VALUE IN G.P. SEQUENCE ARRAY.'                                                              &
                    ,/,14X,' GRID POINT ',I8,' WILL USE SEQUENCE NUMBER ',1ES13.6)

 1360 FORMAT(' *ERROR  1360: SEQUENCE NUMBERS FOR GRIDS ',I8,' AND ',I8,' ARE BOTH ',1ES13.6,'. SEQUENCE NUMBERS MUST BE UNIQUE')


 1361 FORMAT(' *ERROR  1361: UNDEFINED ',A,I8,' ON ',A)

 9999 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE RCM_SEQ_PROC ( R_GSEQ, RCM_OK )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, medat0_cuserin, MELGP, NELE, NGRID
      USE MODEL_STUF, ONLY            :  EDAT, ELGP, EPNT, ETYPE, GRID_ID
      USE PARAMS, ONLY                :  EPSIL
      USE SEQ_PROC_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME_RCM = 'RCM_SEQ_PROC'
      CHARACTER(8*BYTE)               :: TYPE_E

      REAL(DOUBLE), INTENT(OUT)       :: R_GSEQ(NGRID)
      LOGICAL, INTENT(OUT)            :: RCM_OK

      INTEGER(LONG)                   :: ASTAT
      INTEGER(LONG)                   :: DELTA
      INTEGER(LONG)                   :: E
      INTEGER(LONG)                   :: EPNTK
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IDX
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: MIN_DEG
      INTEGER(LONG)                   :: NG
      INTEGER(LONG)                   :: NODE
      INTEGER(LONG)                   :: ORDER_LEN
      INTEGER(LONG)                   :: P
      INTEGER(LONG)                   :: START_NODE
      INTEGER(LONG)                   :: STAMP
      INTEGER(LONG)                   :: TAIL
      INTEGER(LONG)                   :: HEAD
      INTEGER(LONG)                   :: TOTAL_INC
      INTEGER(LONG)                   :: GRID_ROW_NUM
      INTEGER(LONG)                   :: NBR_LEN
      INTEGER(LONG)                   :: COMP_LEN
      INTEGER(LONG)                   :: NBR_NODE
      INTEGER(LONG)                   :: NSEQ_LOC
      INTEGER(LONG)                   :: TMP
      INTEGER(LONG)                   :: TMP_DEG
      INTEGER(LONG)                   :: ELROW

      INTEGER(LONG), ALLOCATABLE      :: ELEM_NGP(:)
      INTEGER(LONG), ALLOCATABLE      :: ELEM_BGRID(:,:)
      INTEGER(LONG), ALLOCATABLE      :: NODE_ELEM_CNT(:)
      INTEGER(LONG), ALLOCATABLE      :: NODE_ELEM_PTR(:)
      INTEGER(LONG), ALLOCATABLE      :: NODE_ELEM_WRK(:)
      INTEGER(LONG), ALLOCATABLE      :: NODE_ELEM_LIST(:)
      INTEGER(LONG), ALLOCATABLE      :: DEGREE(:)
      INTEGER(LONG), ALLOCATABLE      :: MARK(:)
      INTEGER(LONG), ALLOCATABLE      :: ORDER(:)
      INTEGER(LONG), ALLOCATABLE      :: QUEUE(:)
      INTEGER(LONG), ALLOCATABLE      :: COMP_ORDER(:)
      INTEGER(LONG), ALLOCATABLE      :: NEIGH(:)
      LOGICAL, ALLOCATABLE            :: VISITED(:)

      RCM_OK = .FALSE.

      ALLOCATE(ELEM_NGP(NELE), ELEM_BGRID(MELGP,NELE), NODE_ELEM_CNT(NGRID), NODE_ELEM_PTR(NGRID+1), STAT=ASTAT)
      IF (ASTAT /= 0) GOTO 900
      ELEM_NGP      = 0
      ELEM_BGRID    = 0
      NODE_ELEM_CNT = 0
      NODE_ELEM_PTR = 0

      DO E=1,NELE
         CALL GET_ELGP ( E )
         NG = ELGP
         ELEM_NGP(E) = NG
         EPNTK = EPNT(E)
         TYPE_E = ETYPE(E)
         DELTA = 1
         IF (TYPE_E == 'BUSH    ') DELTA = 2
         IF (TYPE_E == 'USERIN  ') DELTA = medat0_cuserin - 1

         NSEQ_LOC = 0
         DO I=1,NG
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME_RCM, NGRID, GRID_ID, EDAT(EPNTK+I+DELTA), GRID_ROW_NUM )
            IF (GRID_ROW_NUM > 0) THEN
               NSEQ_LOC = NSEQ_LOC + 1
               ELEM_BGRID(NSEQ_LOC,E) = GRID_ROW_NUM
               NODE_ELEM_CNT(GRID_ROW_NUM) = NODE_ELEM_CNT(GRID_ROW_NUM) + 1
            ENDIF
         ENDDO
         ELEM_NGP(E) = NSEQ_LOC
      ENDDO

      TOTAL_INC = 0
      DO I=1,NGRID
         NODE_ELEM_PTR(I) = TOTAL_INC + 1
         TOTAL_INC = TOTAL_INC + NODE_ELEM_CNT(I)
      ENDDO
      NODE_ELEM_PTR(NGRID+1) = TOTAL_INC + 1

      ALLOCATE(NODE_ELEM_LIST(MAX(TOTAL_INC,1)), NODE_ELEM_WRK(NGRID), DEGREE(NGRID), MARK(NGRID), ORDER(NGRID),                   &
               QUEUE(NGRID), COMP_ORDER(NGRID), NEIGH(NGRID), VISITED(NGRID), STAT=ASTAT)
      IF (ASTAT /= 0) GOTO 900

      NODE_ELEM_LIST = 0
      NODE_ELEM_WRK  = NODE_ELEM_PTR(1:NGRID)
      DEGREE         = 0
      MARK           = 0
      ORDER          = 0
      QUEUE          = 0
      COMP_ORDER     = 0
      NEIGH          = 0
      VISITED        = .FALSE.

      DO E=1,NELE
         NG = ELEM_NGP(E)
         IF (NG < 2) CYCLE
         DO I=1,NG
            NODE = ELEM_BGRID(I,E)
            IF (NODE <= 0) CYCLE
            IDX = NODE_ELEM_WRK(NODE)
            NODE_ELEM_LIST(IDX) = E
            NODE_ELEM_WRK(NODE) = IDX + 1
         ENDDO
      ENDDO

      STAMP = 0
      DO I=1,NGRID
         STAMP = STAMP + 1
         DO P=NODE_ELEM_PTR(I),NODE_ELEM_PTR(I+1)-1
            ELROW = NODE_ELEM_LIST(P)
            IF (ELROW <= 0) CYCLE
            NG = ELEM_NGP(ELROW)
            DO J=1,NG
               NODE = ELEM_BGRID(J,ELROW)
               IF ((NODE > 0) .AND. (NODE /= I)) THEN
                  IF (MARK(NODE) /= STAMP) THEN
                     MARK(NODE) = STAMP
                     DEGREE(I) = DEGREE(I) + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ORDER_LEN = 0
      DO
         START_NODE = 0
         MIN_DEG = HUGE(MIN_DEG)
         DO I=1,NGRID
            IF (.NOT. VISITED(I)) THEN
               IF (DEGREE(I) < MIN_DEG) THEN
                  MIN_DEG = DEGREE(I)
                  START_NODE = I
               ENDIF
            ENDIF
         ENDDO
         IF (START_NODE == 0) EXIT

         HEAD = 1
         TAIL = 1
         QUEUE(1) = START_NODE
         VISITED(START_NODE) = .TRUE.
         COMP_LEN = 0

         DO WHILE (HEAD <= TAIL)
            NODE = QUEUE(HEAD)
            HEAD = HEAD + 1
            COMP_LEN = COMP_LEN + 1
            COMP_ORDER(COMP_LEN) = NODE

            STAMP = STAMP + 1
            NBR_LEN = 0
            DO P=NODE_ELEM_PTR(NODE),NODE_ELEM_PTR(NODE+1)-1
               ELROW = NODE_ELEM_LIST(P)
               IF (ELROW <= 0) CYCLE
               NG = ELEM_NGP(ELROW)
               DO J=1,NG
                  NBR_NODE = ELEM_BGRID(J,ELROW)
                  IF ((NBR_NODE > 0) .AND. (.NOT. VISITED(NBR_NODE)) .AND. (NBR_NODE /= NODE)) THEN
                     IF (MARK(NBR_NODE) /= STAMP) THEN
                        MARK(NBR_NODE) = STAMP
                        NBR_LEN = NBR_LEN + 1
                        NEIGH(NBR_LEN) = NBR_NODE
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO

            IF (NBR_LEN > 1) THEN
               DO J=2,NBR_LEN
                  TMP = NEIGH(J)
                  TMP_DEG = DEGREE(TMP)
                  K = J - 1
                  DO WHILE (K >= 1)
                     IF ((DEGREE(NEIGH(K)) > TMP_DEG) .OR. ((DEGREE(NEIGH(K)) == TMP_DEG) .AND. (NEIGH(K) > TMP))) THEN
                        NEIGH(K+1) = NEIGH(K)
                        K = K - 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO
                  NEIGH(K+1) = TMP
               ENDDO
            ENDIF

            DO J=1,NBR_LEN
               NBR_NODE = NEIGH(J)
               IF (.NOT. VISITED(NBR_NODE)) THEN
                  VISITED(NBR_NODE) = .TRUE.
                  TAIL = TAIL + 1
                  QUEUE(TAIL) = NBR_NODE
               ENDIF
            ENDDO
         ENDDO

         DO J=COMP_LEN,1,-1
            ORDER_LEN = ORDER_LEN + 1
            ORDER(ORDER_LEN) = COMP_ORDER(J)
         ENDDO
      ENDDO

      IF (ORDER_LEN /= NGRID) GOTO 900

      DO I=1,NGRID
         R_GSEQ(ORDER(I)) = DBLE(I)
      ENDDO

      RCM_OK = .TRUE.

  900 CONTINUE
      IF (ALLOCATED(ELEM_NGP)) DEALLOCATE(ELEM_NGP)
      IF (ALLOCATED(ELEM_BGRID)) DEALLOCATE(ELEM_BGRID)
      IF (ALLOCATED(NODE_ELEM_CNT)) DEALLOCATE(NODE_ELEM_CNT)
      IF (ALLOCATED(NODE_ELEM_PTR)) DEALLOCATE(NODE_ELEM_PTR)
      IF (ALLOCATED(NODE_ELEM_WRK)) DEALLOCATE(NODE_ELEM_WRK)
      IF (ALLOCATED(NODE_ELEM_LIST)) DEALLOCATE(NODE_ELEM_LIST)
      IF (ALLOCATED(DEGREE)) DEALLOCATE(DEGREE)
      IF (ALLOCATED(MARK)) DEALLOCATE(MARK)
      IF (ALLOCATED(ORDER)) DEALLOCATE(ORDER)
      IF (ALLOCATED(QUEUE)) DEALLOCATE(QUEUE)
      IF (ALLOCATED(COMP_ORDER)) DEALLOCATE(COMP_ORDER)
      IF (ALLOCATED(NEIGH)) DEALLOCATE(NEIGH)
      IF (ALLOCATED(VISITED)) DEALLOCATE(VISITED)

      RETURN

      END SUBROUTINE RCM_SEQ_PROC

! ##################################################################################################################################

      SUBROUTINE AUTO_SEQ_PROC

! Reads SEQGP card images from bandit output file (filename.SEQ) using subr BD_SEQGP which creates SEQ1, SEQ2 arrays from SEQGP info

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, sc1
      USE SCONTR, ONLY                :  BANDIT_ERR, BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, LSEQ, NGRID, NSEQ,          &
                                         PROG_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE PARAMS, ONLY                :  GRIDSEQ, SEQPRT, SEQQUIT, SUPINFO

      IMPLICIT NONE

      LOGICAL                         :: LEXIST              ! T/F depending on whether file bandit.f07 exists

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'AUTO_SEQ_PROC'
      CHARACTER( 1*BYTE)              :: IS_GRID_SEQD(NGRID) ! 'Y'/'N' indicator of whether a grid was sequenced in file SEQ
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD                ! Card image read from SEQ file (which was generated by Bandit)
      CHARACTER(85*BYTE)              :: VEC_DESCR           ! Char name of array which has the grids that Bandit did not sequence

      INTEGER(LONG)                   :: BANDIT_MAX_SEQ_NUM  ! Max sequence number from SEQ2 that is generated by Bandit sequencing
      INTEGER(LONG)                   :: BANDIT_NSEQ         ! Number of grids sequenced by Bandit
      INTEGER(LONG)                   :: GRID_NOT_SEQD(NGRID)! Array of grid numbers not sequenced by Bandit
      INTEGER(LONG)                   :: GRID_ROW_NUM        ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: IERR0       = 0     ! Local error count
      INTEGER(LONG)                   :: IERR1       = 0     ! Local error count
      INTEGER(LONG)                   :: INT_SEQ2            ! Integer value for real SEQ2 sequence number
      INTEGER(LONG)                   :: IOCHK       = 0     ! IOSTAT error number when opening/reading a file
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD
      INTEGER(LONG)                   :: K                   ! Counter
      INTEGER(LONG)                   :: NUM_GRIDS_NOT_SEQD=0! Number of grids not sequenced by Bandit
      INTEGER(LONG)                   :: NUM_LEFT            ! Number used in writing 4 pairs of SEQ1,2 to SEQ file
      INTEGER(LONG)                   :: NUM_SEQ_FILE_LINES  ! Number of lines in the SEQ file (after STIME)
      INTEGER(LONG)                   :: OUNT(2)             ! File units to write messages to
      INTEGER(LONG)                   :: REC_NO      = 0     ! Indicator of record number when error encountered reading file
      INTEGER(LONG)                   :: XTIME               ! Time stamp read from a file


      INTRINSIC                       :: DBLE, INT



! **********************************************************************************************************************************
      IF (BANDIT_ERR /= 0) THEN
         WRITE(ERR,8888) BANDIT_ERR
         WARN_ERR = WARN_ERR + 1
         WRITE(F06,8888) BANDIT_ERR
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8888) BANDIT_ERR
         ENDIF
         CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, SEQSTAT )
         RETURN
      ENDIF

      NUM_SEQ_FILE_LINES = 0

      OUNT(1) = ERR
      OUNT(2) = F06

      INQUIRE ( FILE=SEQFIL, EXIST=LEXIST )
exist:IF (LEXIST) THEN                                     ! SEQFIL does exist, so open and read data

         IERR0 = 0                                         ! Open SEQ file, read and check STIME
         OPEN (SEQ,FILE=SEQFIL,STATUS='OLD',IOSTAT=IOCHK)
         IF (IOCHK /= 0) THEN
            CALL OPNERR ( IOCHK, SEQFIL, OUNT )
            IERR0 = IERR0 +1
            IF (IOCHK < 0) THEN                            ! File cannot be opened
               WRITE(ERR,9991) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9991) SEQFIL, GRIDSEQ
                  WRITE(F06,9998)
               ENDIF
            ELSE                                           ! ERR reading file
               WRITE(ERR,9999) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               WRITE(F06,9999) SEQFIL, GRIDSEQ
               WRITE(F06,9998)
            ENDIF
            CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, 'KEEP' )
            RETURN
         ELSE                                              ! No OPEN error, so read and check STIME
            READ(SEQ,'(1X,I11)',IOSTAT=IOCHK) XTIME
            IF (IOCHK /= 0) THEN
               REC_NO = 1
               CALL READERR ( IOCHK, SEQFIL, 'STIME', REC_NO, OUNT )
               IERR0 = IERR0 + 1
            ELSE                                           ! No error reading XTIME, so check XTIME = STIME
               IF (XTIME /= STIME) THEN
                  CALL STMERR ( XTIME, SEQFIL, OUNT )
                  IERR0 = IERR0 + 1
               ENDIF
            ENDIF
         ENDIF

 err0:   IF (IERR0 == 0) THEN                              ! No problem opening SEQ file and reading STIME, so process SEQGP entries
            CALL DEALLOCATE_MODEL_STUF ( 'SEQ1,2' )
            LSEQ = NGRID
            NSEQ = 0
            CALL ALLOCATE_MODEL_STUF ( 'SEQ1,2', SUBR_NAME )

            REC_NO = 0
            IERR1  = 0
i_do1:      DO                                             ! Loop reading SEQGP records until end of file
               READ (SEQ,'(A)',IOSTAT=IOCHK) CARD
               REC_NO = REC_NO + 1
               IF      (IOCHK <  0) THEN                   ! EOF/EOR so exit
                  EXIT
               ELSE IF (IOCHK >  0) THEN                   ! Error reading a SEQGP card
                  CALL READERR ( IOCHK, SEQFIL, 'SEQGP cards', REC_NO, OUNT )
                  IERR1 = IERR1 + 1
               ELSE                                        ! READ was OK so process record
                  IF (CARD(1:5) == 'SEQGP   ') THEN        ! If SEQGP image, call BD_SEQGP to read it and to add to SEQ1,2 arrays
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                     CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
                     DO J=1,10
                        CALL LEFT_ADJ_BDFLD ( JCARD(J) )
                     ENDDO
                     CALL MKCARD ( JCARD, CARD )
                     CALL BD_SEQGP ( CARD )
                  ELSE
                     CYCLE i_do1
                  ENDIF
               ENDIF
            ENDDO i_do1
            BANDIT_NSEQ = NSEQ

            IF (IERR1 > 0) THEN
               WRITE(ERR,*)
               WRITE(ERR,9992) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,*)
                  WRITE(F06,9992) SEQFIL, GRIDSEQ
                  WRITE(F06,9998)
               ENDIF
               CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, 'KEEP' )
            ENDIF

            IF (NUM_SEQ_FILE_LINES > 0) THEN
               WRITE(ERR,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
               WRITE(ERR,*)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
                  WRITE(F06,*)
               ENDIF
            ELSE
               WRITE(ERR,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
               WRITE(ERR,9994)
               WRITE(ERR,9998)
               WRITE(ERR,*)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
                  WRITE(F06,9994)
                  WRITE(F06,9998)
                  WRITE(F06,*)
               ENDIF
            ENDIF

! Finished reading the SEQGP entries in the Bandit SEQ file. Now find out if all grids were sequenced. If not, sequence them
! after the last Bandit sequence number, add them to the SEQ file, and tell user.

            NUM_GRIDS_NOT_SEQD = NGRID - NSEQ

            IF (NUM_GRIDS_NOT_SEQD == 0) THEN                ! Bandit sequenced all grids

               WRITE(ERR,10000) NGRID, GRIDSEQ
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,10000) NGRID, GRIDSEQ
               ENDIF

            ELSE                                             ! Some grids were not sequenced by Bandit

               WRITE(ERR,9996) SEQFIL, GRIDSEQ, NGRID, NSEQ
               WRITE(F06,9996) SEQFIL, GRIDSEQ, NGRID, NSEQ

               DO I=1,NGRID                                  ! Array IS_GRID_SEQD will show which grids were not seq'd by Bandit
                  IS_GRID_SEQD(I) = 'N'
               ENDDO
               DO I=1,NSEQ
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SEQ1(I), GRID_ROW_NUM )
                  IF (GRID_ROW_NUM > 0) THEN
                     IS_GRID_SEQD(GRID_ROW_NUM) = 'Y'
                  ENDIF
               ENDDO

               K = 0                                       ! Array GRID_NOT_SEQD has grid no's for grids not seq'd by Bandit
               DO I=1,NGRID
                  IF (IS_GRID_SEQD(I) == 'N') THEN
                     K = K + 1
                     GRID_NOT_SEQD(K) = GRID_ID(I)
                  ENDIF
               ENDDO
                                                           ! Write list of grids not sequenced by Bandit
               CALL SORT_INT1 ( SUBR_NAME, 'GRID_NOT_SEQD', NUM_GRIDS_NOT_SEQD, GRID_NOT_SEQD )
               VEC_DESCR(1:) = 'GRIDS WHICH BANDIT DID NOT SEQUENCE (PERHAPS THEY HAD NO ELEMENTS CONNECTED TO THEM)'
               IF (NUM_GRIDS_NOT_SEQD > 10) THEN
                  CALL WRITE_INTEGER_VEC ( VEC_DESCR, GRID_NOT_SEQD, NUM_GRIDS_NOT_SEQD )
               ELSE
                  WRITE(F06,10003) VEC_DESCR
                  DO I=1,NUM_GRIDS_NOT_SEQD
                     WRITE(F06,10004) GRID_NOT_SEQD(I)
                  ENDDO
               ENDIF
               WRITE(F06,*)
                                                           ! Sequence the grids not sequenced by Bandit (in grid numerical order)
               BANDIT_MAX_SEQ_NUM = 0                      !  (a) Get max sequence number generated by Bandit
               DO I=1,NSEQ
                  INT_SEQ2 = INT(SEQ2(I))
                  IF (INT_SEQ2 > BANDIT_MAX_SEQ_NUM) THEN
                     BANDIT_MAX_SEQ_NUM = INT_SEQ2
                  ENDIF
               ENDDO
               IF ((NSEQ + NUM_GRIDS_NOT_SEQD) > LSEQ) THEN!  (b) Check to make sure we won't exceed size of allocated arrays SEQ1,2
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1301) SUBR_NAME,NSEQ,NUM_GRIDS_NOT_SEQD,LSEQ
                  WRITE(F06,1301) SUBR_NAME,NSEQ,NUM_GRIDS_NOT_SEQD,LSEQ
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               DO I=1,NUM_GRIDS_NOT_SEQD                   !  (c) Sequence the grids starting with previous NSEQ (from Bandit)
                  NSEQ = NSEQ + 1
                  SEQ1(NSEQ) = GRID_NOT_SEQD(I)            !      NOTE: GRID_NOT_SEQD was sorted above so it is in grid num order
                  SEQ2(NSEQ) = DBLE(BANDIT_MAX_SEQ_NUM + I)
               ENDDO

               NUM_LEFT = NUM_GRIDS_NOT_SEQD               ! Write data for grids not seq'd by Bandit to the SEQ file
               BACKSPACE (SEQ)
               WRITE(SEQ,300) PROG_NAME, NUM_LEFT
               NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 2
               DO I=1,NUM_GRIDS_NOT_SEQD,4
                  IF (NUM_LEFT >= 4) THEN
                     WRITE(SEQ,301) (SEQ1(K),INT(SEQ2(K)),K=I+BANDIT_NSEQ,I+BANDIT_NSEQ+3)
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                     NUM_LEFT = NUM_LEFT - 4
                  ELSE
                     WRITE(SEQ,301) (SEQ1(K),INT(SEQ2(K)),K=I+BANDIT_NSEQ,I+BANDIT_NSEQ+NUM_LEFT-1)
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                  ENDIF
               ENDDO

               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,9998)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,9998)
               ENDIF

            ENDIF

! Write SEQ cards to F06, if user requested it

            IF (SEQPRT == 'Y') THEN                        ! See if user wants to write SEQGP card images to F06

               IERR1 = 0
               REWIND (SEQ)
               READ(SEQ,'(1X,I11)',IOSTAT=IOCHK) XTIME
               IF (IOCHK /= 0) THEN
                  REC_NO = 1
                  CALL READERR ( IOCHK, SEQFIL, 'STIME', REC_NO, OUNT )
                  IERR1 = IERR1 + 1
               ELSE                                        ! No error reading XTIME, so check XTIME = STIME
                  IF (XTIME /= STIME) THEN
                     CALL STMERR ( XTIME, SEQFIL, OUNT )
                     IERR1 = IERR1 + 1
                  ENDIF
               ENDIF

               IF (IERR1 == 0) THEN
                  REC_NO = 0
                  WRITE(F06,9997) SEQFIL
i_do2:            DO I=1,NUM_SEQ_FILE_LINES
                     READ (SEQ,'(A)',IOSTAT=IOCHK) CARD
                     REC_NO = REC_NO + 1
                     IF (IOCHK /= 0) THEN
                        WRITE(F06,10001) I
                     ELSE
                        WRITE(F06,10002) I, CARD
                     ENDIF
                  ENDDO i_do2
                  WRITE(F06,*)
               ENDIF

            ENDIF

            CALL FILE_CLOSE ( SEQ, SEQFIL, SEQSTAT )

         ELSE

            CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, SEQSTAT )

         ENDIF err0

      ELSE                                                 ! bandit.f07 does not exist. Write message and quit if SEQQUIT = Y

         WRITE(ERR,9990) SEQFIL, GRIDSEQ
         WRITE(ERR,9994)
         WRITE(ERR,9998)

         IF (SUPINFO == 'N') THEN
            WRITE(F06,9990) SEQFIL, GRIDSEQ
            WRITE(F06,9994)
            WRITE(F06,9998)
         ENDIF

         RETURN

      ENDIF exist



      RETURN

! **********************************************************************************************************************************
  300 FORMAT('$ The following SEQGP entries were generated by ',A,' to sequence',/,                                                &
             '$ the ',I8,' grids not sequenced by Bandit.')

  301 FORMAT('SEQGP   ',8I8)

 1301 FORMAT(' *ERROR  1301: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' (NSEQ + NUM_GRIDS_NOT_SEQD) = ',I8,' MUST BE < LSEQ = ',I8,' SO ALLOCATED MEM FOR ARRAYS SEQ1,2',     &
                           ' WON''T BE EXCEEDED')


 8888 FORMAT(' *WARNING    : BANDIT DID NOT RUN SUCCESSFULLY. THIS IS A NON-CRITICAL OPTIMIZATION. IT QUIT WITH ERROR = ',I8       &
                    ,/,14X,' CHECK FILE BANDIT.OUT IN THE DIRECTORY WHERE MYSTRAN.EXE RESIDES')

 9990 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), DOES NOT EXIST',/)

 9991 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), CANNOT BE OPENED.',/)

 9999 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), CANNOT BE READ.',/)

 9992 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), HAS HAD ABOVE LISTED READ ERROR(S).',/)

 9993 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), HAS BEEN READ. THERE WERE ',I8,' SEQGP ENTRY IMAGES READ.')

 9994 FORMAT(14X,' IT MAY BE THAT BANDIT FOUND THAT NO RESEQUENCING WAS NEEDED OR DUE TO ERROR IN RUNNING BANDIT.',/)

 9996 FORMAT(' *WARNING    : FILE ',A                                                                                              &
                    ,/,14x,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,') DOES NOT HAVE ALL GRIDS DEFINED ON SEQGP ENTRIES.'                  &
                    ,/,14X,' THE NUMBER OF GRIDS = ',I8,' BUT THE NUMBER OF GRIDS SEQUENCED = ',I8,/)

 9997 FORMAT(' *INFORMATION: THE FOLLOWING IS A LISTING OF THE SEQGP RECORDS READ FROM FILE ',A,':',/)

 9998 FORMAT(14X,' MAKE SURE BANDIT HAS RUN SUCCESSFULLY (CHECK FILE bandit.out IN THE DIRECTORY WHERE MYSTRAN.EXE RESIDES).'/)

10000 FORMAT(' *INFORMATION: ALL BANDIT SEQGP ENTRY IMAGES HAVE BEEN USED TO SUCCESSFULLY RESEQUENCE ',I8,' GRIDS BASED ON PARAM', &
                           ' GRIDSEQ ',A,/)

10001 FORMAT(15X,'Card image ',I7,': could not be written')

10002 FORMAT(15X,'Card image ',I7,': ',A)

10003 FORMAT(24X,A,/)

10004 FORMAT(62X,I8)

! **********************************************************************************************************************************

      END SUBROUTINE AUTO_SEQ_PROC

! ##################################################################################################################################

      SUBROUTINE AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, CLOSE_STAT )

! Writes message for subr AUTO_SEQ_PROC

      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  SEQQUIT

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: SUBR_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! Status for closing SEQFIL

! **********************************************************************************************************************************
      CALL FILE_CLOSE ( SEQ, SEQFIL, CLOSE_STAT )

      IF (SEQQUIT == 'Y') THEN
         WRITE(ERR,8881) SUBR_NAME, SEQQUIT
         WRITE(F06,8881) SUBR_NAME, SEQQUIT
         CALL OUTA_HERE ( 'Y' )
      ELSE
         WRITE(ERR,8886) SEQQUIT
         WRITE(F06,8886) SEQQUIT
      ENDIF

      RETURN

! **********************************************************************************************************************************
 8881 FORMAT(14X,' PROCESSING TERMINATED IN SUBR ',A,' BASED ON ABOVE PROBLEM AND FIELD 4 OF PARAM GRIDSEQ = ',A,'.'               &
          ,/,14X,' (THIS WAS EITHER ENTERED ON A BULK DATA PARAM GRIDSEQ ENTRY OR IS THE DEFAULT. SEE MYSTRAN DOCUMENTATION)',/)

 8886 FORMAT(14X,' SINCE FIELD 4 OF PARAM GRIDSEQ = ',A,', MYSTRAN WILL DEFAULT TO GRID SEQUENCING BASED ON GRID NUMERICAL ORDER.' &
          ,/,14X,' (THIS WAS EITHER ENTERED ON A BULK DATA PARAM GRIDSEQ ENTRY OR IS THE DEFAULT. SEE MYSTRAN DOCUMENTATION)',/)

! **********************************************************************************************************************************

      END SUBROUTINE AUTO_SEQ_PROC_WRAPUP

      END SUBROUTINE SEQ_PROC
