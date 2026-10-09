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

   MODULE LINK1_LOAD_SUPPORT

   USE DOF_LOOKUP_UTILS, ONLY :  TDOF_COL_NUM

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EPTL, GET_GRID_6X6_MASS, PRESSURE_DATA_PROC, SLOAD_PROC, TEMPERATURE_DATA_PROC, YS_ARRAY

   CONTAINS

      SUBROUTINE EPTL

! Element pressure and thermal loads processor

! Processes the elems to generate the elem pressure and thermal loads using the EMG set of routines.  The thermal
! and pressure loads are inserted into the system loads array SYS_LOAD.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, F21, F21FIL, F21_MSG, SC1, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_P_T_BIT, ELDT_F21_P_T_BIT, IBIT, LINKNO, MBUG, MELDOF, NCORD,      &
                                         NELE, NGRID, NSUB, NTSUB
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  ELDOF, ELDT, GRID, GRID_ID, CORD, AGRID, ELGP, NUM_EMG_FATAL_ERRS, OELDT, PLY_NUM, PPE,   &
                                         PTE, SYS_LOAD, TYPE, SUBLOD

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE EMG_MOD, ONLY               :  EMG
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_FIJFIL
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EPTL'
      CHARACTER( 1*BYTE)              :: OPT(6)

      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable used in setting WRT_BUG(3) and OUT10
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: L                 ! Counter
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: ROW_NUM           ! Row no. in array TDOF corresponding to an elem DOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! Row no. in array TDOF where data begins for AGRID
      INTEGER(LONG)                   :: TCASE2(NSUB)      ! TCASE2(I) gives the internal subcase no. for internal thermal case I
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE2(1-5) = 3, 5, 0, 0, 0
                                                           ! Indicator for output of elem data to BUG file


      REAL(DOUBLE)                    :: DZE(MELDOF,MELDOF)! A dummy array for the call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      INTRINSIC IAND

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

!! Initialize

      EPS1 = EPSIL(1)

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Null dummy array DZE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,MELDOF
            DZE(I,J) = ZERO
         ENDDO
      ENDDO

! Thermal load flag (OPT(2)). Need TCASE2 - table relating jth thermal case to ith subcase

      IF (NTSUB > 0) THEN
         OPT(2) = 'Y'                                       ! OPT(2) is for calc of PTE
         DO I = 1,NSUB
            TCASE2(I) = 0
         ENDDO
         J = 0
         DO I = 1,NSUB
            IF (SUBLOD(I,2) /= 0) THEN
               J = J + 1
               TCASE2(J) = I
            ENDIF
         ENDDO
      ELSE
         OPT(2) = 'N'
      ENDIF

! Set element pressure calc flag

      OPT(5) = 'Y'                                          ! OPT(5) is for calc of PPE

! Process the elements:

      IERROR = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      CALL COUNTER_INIT('Calculating load matrix for element', NELE)
      DO I = 1,NELE

         DO J=0,MBUG-1
            WRT_BUG(J) = 0
         ENDDO

         IF (LINKNO == 1) THEN                             ! Only want element PTE, PPE, written to BUG file in LINK 1
            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_P_T_BIT))         ! WRT_BUG(2) is for printed output of PTE, PPE
            IF (I1 > 0) THEN
               WRT_BUG(2) = 1
            ENDIF
         ENDIF


         OPT(1) = 'N'                                      ! OPT(1) is for calc of ME
         OPT(3) = 'N'                                      ! OPT(3) is for calc of SEi, STEi
         OPT(4) = 'N'                                      ! OPT(4) is for calc of KE
         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'Y' )       ! 'Y' means write to BUG file

         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE
         ENDIF

         I1 = IAND(OELDT,IBIT(ELDT_F21_P_T_BIT))           ! Do we need to write elem thermal load matrices to F21 files
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 1, 0 )
         ENDIF

         L = 0                                             ! Generate element DOF'S
         DO J = 1,ELGP
!xx         CALL CALC_TDOF_ROW_NUM ( AGRID(J), ROW_NUM_START, 'N' )
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(J), IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
               ROW_NUM = ROW_NUM_START + K - 1
               L       = L + 1
               EDOF(L) = TDOF(ROW_NUM,G_SET_COL_NUM)
            ENDDO
         ENDDO

         IF (OPT(2) == 'Y') THEN                           ! Process the element thermal loads - PTE array:
                                                           ! Transform PTE from local-basic-global
            IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN
               CALL ELEM_TRANSFORM_LBG ( 'PTE', DZE, PTE )
            ENDIF

            DO J = 1,ELDOF                                 ! Load PTE into SYS_LOAD array.
k_do113:       DO K = 1,NTSUB
                  IF (DABS(PTE(J,K)) < EPS1) CYCLE k_do113
                  SYS_LOAD(EDOF(J),TCASE2(K)) = SYS_LOAD(EDOF(J),TCASE2(K)) + PTE(J,K)
               ENDDO k_do113
            ENDDO

         ENDIF

         IF (OPT(5) == 'Y') THEN                           ! Process element pressure loads - PPE array.
                                                           ! Transform PTE from local-basic-global
            IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN
               CALL ELEM_TRANSFORM_LBG ( 'PPE', DZE, PPE )
            ENDIF

            DO J = 1,ELDOF                                 ! Load PPE into SYS_LOAD array.
k_do123:       DO K = 1,NSUB
                  IF (DABS(PPE(J,K)) < EPS1) CYCLE k_do123
                  SYS_LOAD(EDOF(J),K) = SYS_LOAD(EDOF(J),K) + PPE(J,K)
               ENDDO k_do123
            ENDDO

         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO

      WRITE(SC1,*) CR13

! Reset option values:

      OPT(2) = 'N'
      OPT(5) = 'N'

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! Errors from subr EMG, so quit
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

98712 format('J, PPE(J) = ',i8,10(1es15.6))

! **********************************************************************************************************************************

      END SUBROUTINE EPTL


      SUBROUTINE GET_GRID_6X6_MASS ( AGRID, IGRID, FOUND, GRID_MGG )

! Gets a 6 x 6 mass matrix for 1 grid point from the MGG mass matrix (which contains block diagonal 6 x 6 grid mass matrices).
! THis subr was not coded for SPOINT's so check if AGRID is an SPOINT and give program error and quit if it is

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NGRID, NTERM_MGG
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOF
      USE SPARSE_MATRICES, ONLY       :  I2_MGG, J_MGG, MGG
      USE MODEL_STUF, ONLY            :  GRID_SEQ

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_GRID_6X6_MASS'
      CHARACTER( 1*BYTE), INTENT(OUT) :: FOUND             ! 'Y' if there is a mass matrix for this grid and 'N' otherwise

      INTEGER(LONG), INTENT(IN)       :: AGRID             ! Actual grid number of grid for which we want the 6 x 6 mass matrix
      INTEGER(LONG), INTENT(IN)       :: IGRID             ! Internal grid number of grid for which we want the 6 x 6 mass matrix
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: I1,J1             ! Indices
      INTEGER(LONG)                   :: IGRID_DOF_NUM     ! G-set DOF number for IGRID
      INTEGER(LONG)                   :: NUM_COMPS         ! No. displ components (1 for SPOINT, 6 for actual grid)


      REAL(DOUBLE), INTENT(OUT)       :: GRID_MGG(6,6)     ! 6 x 6 mass matrix for internal grid IGRID

      INTRINSIC                       :: MODULO



! **********************************************************************************************************************************
! If AGRID is an SPOINT give error and quit

      CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
      IF (NUM_COMPS /= 6) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1502) AGRID, NUM_COMPS
         WRITE(F06,1502) AGRID, NUM_COMPS
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Initialize outputs

      FOUND = 'N'

      DO I=1,6
         DO J=1,6
            GRID_MGG(I,J) = ZERO
         ENDDO
      ENDDO

      IGRID_DOF_NUM = 6*(GRID_SEQ(IGRID) - 1) + 1
k_do: DO K=1,NTERM_MGG

         IF ((I2_MGG(K) >= IGRID_DOF_NUM) .AND. (I2_MGG(K) <= IGRID_DOF_NUM + 5)) THEN

            I1 = MODULO(I2_MGG(K),6)
            IF (I1 > 0) THEN
               I = I1
            ELSE
               I = 6
            ENDIF

            J1 = MODULO( J_MGG(K),6)
            IF (J1 > 0) THEN
               J = J1
            ELSE
               J = 6
            ENDIF

            GRID_MGG(I,J) = MGG(K)
            FOUND = 'Y'

         ELSE

            IF (FOUND == 'Y') THEN
               EXIT k_do
            ELSE
               CYCLE k_do
            ENDIF

         ENDIF

      ENDDO k_do

      DO I=1,6
         DO J=1,I-1
            GRID_MGG(I,J) = GRID_MGG(J,I)
         ENDDO
      ENDDO




      RETURN

! **********************************************************************************************************************************
 1500 FORMAT(' *ERROR  1500: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT FIND INTERNAL GRID ID, IN ARRAY TDOF, FOR ACTUAL GRID ',I8)

 1501 FORMAT(' *ERROR  1501: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INTERNAL GRID ID, FOUND IN ARRAY TDOF, FOR ACTUAL GRID ',I8,' MUST BE > 0 BUT IS ',I8)

 1502 FORMAT(' *ERROR  1515: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' SUBR NOT PROGRAMMED FOR ANYTHING BUT 6 COMP GRIDS BUT WAS CALLED FOR GRID ',I8,' WHICH HAS ',I3,      &
                           ' NUMBER OF COMPONENTS')


! **********************************************************************************************************************************

      END SUBROUTINE GET_GRID_6X6_MASS


      SUBROUTINE PRESSURE_DATA_PROC

! Processes element pressure data. The element overall pressure and/or element G.P. pressure data was written to file L1Q when the
! element pressure Bulk Data entries were read. In this subr, records are read from L1Q and the data processed into arrays PPNT and
! PDATA (arrays used in the element generation routines)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     L1Q
      USE IOUNT1, ONLY                :  WRT_ERR,                            LINK1Q
      USE IOUNT1, ONLY                :  WRT_ERR,                            L1Q_MSG
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, DATA_NAM_LEN, FATAL_ERR, JCARD_LEN, LPDAT, LLOADC,            &
                                         MPDAT_PLOAD1, MPDAT_PLOAD2, MPDAT_PLOAD4, MPLOAD4_3D_DATA, NELE, NLOAD, NPCARD,           &
                                         NPLOAD4_3D, NPDAT, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  LOAD_SIDS, LOAD_FACS, SUBLOD, PDATA, PPNT, PLOAD4_3D_DATA, PTYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN, READERR
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PRESSURE_DATA_PROC'
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! B.D elem pressure card image read from file LINK1Q
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes
      CHARACTER( 1*BYTE)              :: EFLAG             ! Flag to decide whether a situation in subr EPPUT is an error or not
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters in CARD
      CHARACTER( 8*BYTE)              :: NAME              ! Card name (PLOAD1,2 or 4)
      CHARACTER( 6*BYTE)              :: PLATE_OR_SOLID    ! 'PLATE' or 'SOLID' element designation from PLOAD4 entries
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Variable to test whether "THRU" option was used on B.D. PLOAD2 card
      CHARACTER( 8*BYTE)              :: THRU              ! ='Y' if THRU option used on TEMPRB, TEMPP1 continuation card

      INTEGER(LONG)                   :: EID               ! Actual element ID
      INTEGER(LONG)                   :: EID1,EID2         ! The 2 actual elem ID's in "EID1 THRU EID2" on elem press B.D. card
      INTEGER(LONG)                   :: EL_PRES_ERR       ! Count of error messages when elements have redundant pressures
      INTEGER(LONG)                   :: EL_REDUNDANT_PRES ! Count of warning messages when elements have redundant pressures
      INTEGER(LONG)                   :: G1,G34            ! Grid numbers from a PLOAD4 entry for solid elements
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IELEM             ! Internal element number for the actual element ID, EID
      INTEGER(LONG)                   :: IERROR    = 0     ! Cum. count of errors as we read, and check cards from file LINK1Q
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: IPLOAD4_3D        ! Count of 3D elements that have PLOAD4 face pressure definition
      INTEGER(LONG)                   :: IPPNT             ! Index in array PPNT (a pointer to where press data for EID starts)
      INTEGER(LONG)                   :: LSID(LLOADC+1)    ! Array of load SID's, from a LOAD Bulk Data card, for one S/C
      INTEGER(LONG)                   :: NFIELD            ! No. fields on PLOAD1, PLOAD2 B.D. cards that have pressure data
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SETID             ! Pressure load set ID read from an elem pressure B.D. card
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file


      REAL(DOUBLE)                    :: SCALE             ! Scale factor from a LOAD Bulk Data card
      REAL(DOUBLE)                    :: RPDAT             ! Real pressure value read from file LINK1Q
      REAL(DOUBLE)                    :: RPDAT1            ! Real pressure value read from file LINK1Q
      REAL(DOUBLE)                    :: RSID(LLOADC+1)    ! Array of load magnitudes (for LSID set ID's) needed for one S/C



! **********************************************************************************************************************************
      EL_REDUNDANT_PRES  =  0                              ! Initialize warn, err indicators for redundant elem press definition
      EL_PRES_ERR        =  0
      IPLOAD4_3D         =  0
      PLATE_OR_SOLID(1:) = ' '

      OUNT(1) = ERR                                        ! Set units for writing errors (for subr READERR)
      OUNT(2) = F06

      DO I=1,LLOADC                                         ! Initialize LSID, RSID arrays
         LSID(I) = 0
         RSID(I) = ZERO
      ENDDO

isubc:DO I=1,NSUB                                          ! Loop through the S/C's

         NPDAT = 0                                         ! 09/21/21: Init NPDAT before each S/C. Otherwise can get error 1523

         IF (SUBLOD(I,1) == 0) THEN                        ! If no load for this S/C, CYCLE
            CYCLE isubc
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(I,1)                             ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the elem pressure data in file LINK1Q that matches SUBLOD(I,1)
j_do1:   DO J = 1,NLOAD                                    ! Then, the actual mag. will come from RSID(1)
            IF (SUBLOD(I,1) == LOAD_SIDS(J,1)) THEN
k_do1:         DO K = 2,LLOADC
                  IF (LOAD_SIDS(J,K) == 0) THEN
                     EXIT k_do1
                  ELSE
                     NSID = K                              ! Note: NSID will not get larger than LLOADC
                     RSID(K) = LOAD_FACS(J,1)*LOAD_FACS(J,K)
                     LSID(K) = LOAD_SIDS(J,K)
                   ENDIF
                ENDDO k_do1
            ENDIF
         ENDDO j_do1

pcards:  DO J=1,NPCARD                                     ! Process elem pressure card info
                                                           ! (2-b-  i) Read a record from scratch - forces are in global coords
            READ(L1Q,IOSTAT=IOCHK) CARD                    ! Read element pressure CARD from LINK1Q
            IF (IOCHK /= 0) THEN
               REC_NO = J + 1
               CALL READERR ( IOCHK, LINK1Q, L1Q_MSG, REC_NO, OUNT )
               IERROR = IERROR + 1
               CYCLE isubc
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )        ! Make 10 fields of 8 chars from CARD

            NAME = CARD(1:8)                               ! Set NFIELD based on parent CARD type
            IF      ((NAME(1:7) == 'PLOAD1 ') .OR. (NAME(1:7) == 'PLOAD1*')) THEN
               NFIELD = MPDAT_PLOAD1
            ELSE IF ((NAME(1:7) == 'PLOAD2 ') .OR. (NAME(1:7) == 'PLOAD2*')) THEN
               NFIELD = MPDAT_PLOAD2
            ELSE IF ((NAME(1:7) == 'PLOAD4 ') .OR. (NAME(1:7) == 'PLOAD4*')) THEN
               NFIELD = MPDAT_PLOAD4
            ELSE
               CYCLE pcards                                 ! Ignore record if not for PLOAD1 or PLOAD2
            ENDIF

            IPPNT = NPDAT + 1                              ! Set index for pointer array, PPNT

            READ(JCARD(2),'(I8)') SETID                    ! Get pressure load SID

            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
k_do2:      DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                  EXIT k_do2
               ENDIF
            ENDDO k_do2

            IF (FOUND == 'N') THEN                         ! Cycle back on J loop and read another elem pressure card
               CYCLE pcards
            ENDIF

! Put pressure data into PDATA
            IF ((NPDAT + NFIELD) > LPDAT) THEN             ! Check for overflow in PDATA
               WRITE(ERR,1523) SUBR_NAME,LPDAT
               WRITE(F06,1523) SUBR_NAME,LPDAT
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (dim of array PDATA too small), so quit
            ENDIF

                                                           ! Set the field number where we expect to find the 1st pressure value
            IF      ((NAME(1:7) == 'PLOAD1 ') .OR. (NAME(1:7) == 'PLOAD1*')) THEN
               Write(err,'(a,a)') ' *ERROR     : Code not written for PLOAD1 processing in subr ', subr_name
               Write(f06,'(a,a)') ' *ERROR     : Code not written for PLOAD1 processing in subr ', subr_name
               fatal_err = fatal_err + 1
               call outa_here ( 'Y' )
            ELSE IF ((NAME(1:7) == 'PLOAD2 ') .OR. (NAME(1:7) == 'PLOAD2*')) THEN
               NPDAT = NPDAT + 1
               READ(JCARD(3),'(F16.0)') RPDAT
               PDATA(NPDAT) = SCALE*RPDAT
            ELSE IF ((NAME(1:7) == 'PLOAD4 ') .OR. (NAME(1:7) == 'PLOAD4*')) THEN
               NPDAT = NPDAT + 1
               READ(JCARD(4),'(F16.0)') RPDAT1
               PDATA(NPDAT) = SCALE*RPDAT1
               DO K = 5,7
                  NPDAT = NPDAT + 1
                  IF (JCARD(K)(1:) /= ' ') THEN
                     READ(JCARD(K),'(F16.0)') RPDAT
                     PDATA(NPDAT) = SCALE*RPDAT
                  ELSE
                     PDATA(NPDAT) = SCALE*RPDAT1
                  ENDIF
               ENDDO
            ENDIF

! Process EID's. First check for the 2 options on specifying elem data. For PLOAD2, either all data are EID's or THRU option is used
!                For PLOAD4, either "THRU" is used or there is only 1 EID

            IF      (NAME(1:6) == 'PLOAD2') THEN
               TOKEN = JCARD(5)(1:8)                       ! Only send the 1st 8 chars of this JCARD. It has been left justified
            ELSE IF (NAME(1:6) == 'PLOAD4') THEN
               TOKEN = JCARD(8)(1:8)
            ENDIF
            CALL TOKCHK ( TOKEN, TOKTYP )

            THRU = 'N'
            IF (TOKTYP == 'THRU    ') THEN
               THRU = 'Y'
            ENDIF

            IF (NAME(1:6) == 'PLOAD2') THEN

               IF (THRU == 'N') THEN
                  EFLAG = 'Y'
k_do4:            DO K=4,9
                     IF (JCARD(K) == '        ') THEN      ! See if there is any data in other JCARD fields
                        CYCLE k_do4
                     ELSE
                        READ(JCARD(K),'(I8)') EID
                        CALL EPPUT ( SETID, EID, I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                     ENDIF
                  ENDDO k_do4
               ELSE
                  READ(JCARD(4),'(I8)') EID1
                  READ(JCARD(6),'(I8)') EID2
                  EFLAG = 'N'
k_do5:            DO K=EID1,EID2
                     CALL EPPUT ( SETID, K,   I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  ENDDO k_do5
               ENDIF

            ELSE IF (NAME(1:6) == 'PLOAD4') THEN

               IF (THRU == 'N') THEN
                  READ(JCARD(3),'(I8)') EID
                  CALL EPPUT ( SETID, EID, I, IPPNT, NAME, 'Y', IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  IF (IELEM == -1) CYCLE pcards
                  IF ((JCARD(8)(1:) /= ' ') .AND. (JCARD(9)(1:) /= ' ')) THEN
                     READ(JCARD(8),'(I8)') G1
                     READ(JCARD(9),'(I9)') G34
                     PLATE_OR_SOLID = 'SOLID'
                     IPLOAD4_3D = IPLOAD4_3D + 1
                     PLOAD4_3D_DATA(IPLOAD4_3D,1) = EID
                     PLOAD4_3D_DATA(IPLOAD4_3D,2) = IELEM
                     PLOAD4_3D_DATA(IPLOAD4_3D,3) = I      ! I is the internal subcase number (loop ID is: isubc)
                     PLOAD4_3D_DATA(IPLOAD4_3D,4) = G1
                     PLOAD4_3D_DATA(IPLOAD4_3D,5) = G34
                  ELSE
                     PLATE_OR_SOLID = 'PLATE'
                  ENDIF
               ELSE
                  READ(JCARD(3),'(I8)') EID1
                  READ(JCARD(9),'(I8)') EID2
                  EFLAG = 'N'
k_do6:            DO K=EID1,EID2
                     CALL EPPUT ( SETID, K,   I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  ENDDO k_do6
               ENDIF

            ENDIF

         ENDDO pcards

         IF (EL_REDUNDANT_PRES > 0) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9999)
            IF (SUPWARN == 'Y') THEN
               WRITE(F06,9999)
            ENDIF
         ENDIF

         REWIND (L1Q)                                       ! Need to read all of the elem pressure records again for the next S/C
         READ(L1Q,IOSTAT=IOCHK) XTIME
         IF (IOCHK /= 0) THEN
            REC_NO = 1
            CALL READERR ( IOCHK, LINK1Q, L1Q_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Cannot read STIME from temperature data file, so quit
         ENDIF

      ENDDO isubc

      IF ((IERROR /= 0) .OR. (EL_PRES_ERR > 0)) THEN
         WRITE(ERR,9998) IERROR+EL_PRES_ERR
         WRITE(F06,9998) IERROR+EL_PRES_ERR
         CALL OUTA_HERE ( 'Y' )                            ! Quit due to errors reading elem press file
      ENDIF

! **********************************************************************************************************************************
! Now finally write processed pressure data to L1Q (destroying the PLOAD card data on the records)

! First close and delete L1Q file

      CALL FILE_CLOSE ( L1Q, LINK1Q, 'DELETE' )

! Open L1Q for write:

      CALL FILE_OPEN ( L1Q, LINK1Q, OUNT, 'REPLACE', L1Q_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

      DATA_SET_NAME = 'PPNT'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NELE
      WRITE(L1Q) NSUB
      DO I=1,NELE
         DO J=1,NSUB
            WRITE(L1Q) PPNT(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PDATA'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NPDAT
      DO I=1,NPDAT
         WRITE(L1Q) PDATA(I)
      ENDDO

      DATA_SET_NAME = 'PTYPE'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NELE
      DO I=1,NELE
         WRITE(L1Q) PTYPE(I)
      ENDDO

      DATA_SET_NAME = 'PLOAD4_3D_DATA'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NPLOAD4_3D
      DO I=1,NPLOAD4_3D
         WRITE(L1Q) (PLOAD4_3D_DATA(I,J),J=1,MPLOAD4_3D_DATA)
      ENDDO



      RETURN

! **********************************************************************************************************************************




 9998 FORMAT('PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS')

 9999 FORMAT(' *WARNING    : CHECK ERR OUTPUT FILE FOR WARNING MESSAGES REGARDING REDUNDANT ELEM PRESSURE DEFINITION')

 1523 FORMAT(' *ERROR  1523: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MUCH ELEMENT PRESSURE DATA. MAX IS LPDAT = ',I8)

! ##################################################################################################################################

      CONTAINS

! Change log (changes following completion of Version 1.02 on 05/01/03)

! 08/13/04: (1) Add EL_PRES_ERR to arg list (like EL_REDUNDANT_PRES it gives actual error count as opposed to warning count).
!               Then remove CALL OUTA_HERE and replace with updating EL_PRES_ERR and RETURN
!           (2) Change EL_REDUNDANT_PRES to INOUT arg (was only OUT so total amount of warn errors was not correct)
!           (3) Change arg CARD to SETID and change *ERROR  1320

! ##################################################################################################################################

      SUBROUTINE EPPUT ( SETID, EID, JSUB, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )

! Element pressure routine - generates the PPNT(i,J) array and PTYPE(i) array

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NELE, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ESORT1, ETYPE, SUBLOD, PPNT, PTYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EPPUT'
      CHARACTER(LEN=*), INTENT(IN)    :: EFLAG             ! Flag to decide whether a situation in subr EPPUT is an error or not
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! 'PLOAD1', 'PLOAD2' or 'PLOAD4'

      INTEGER(LONG), INTENT(IN)       :: EID               ! Actual elem ID
      INTEGER(LONG), INTENT(IN)       :: IPPNT             ! Index in array PPNT (a pointer to where press data for EID starts)
      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Internal subcase number
      INTEGER(LONG), INTENT(IN)       :: SETID             ! Actual elem ID
      INTEGER(LONG), INTENT(OUT)      :: IELEM    ! Internal elem ID for actual elem ID EID
      INTEGER(LONG), INTENT(INOUT)    :: EL_REDUNDANT_PRES      ! Count of warning messages when elements have redundant pressures
      INTEGER(LONG), INTENT(INOUT)    :: EL_PRES_ERR       ! Count of errors where elem ID is wrong (*ERROR  1320)




! **********************************************************************************************************************************
! Initialize outputs

      IELEM = 0

! First convert EID to interal value

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, EID, IELEM )

! If EID not found then: error if EFLAG = Y, or return if EFLAG /= Y

      IF (IELEM == -1) THEN
         IF (EFLAG == 'Y') THEN
            WRITE(ERR,1520) EID, SETID
            WRITE(F06,1520) EID, SETID
            FATAL_ERR = FATAL_ERR + 1
            EL_PRES_ERR = EL_PRES_ERR + 1
         ENDIF
         RETURN
      ENDIF

! No error, so put pointer into PPNT(i,j) and give warning if the elem has had a pressure defined previously

      IF (PPNT(IELEM,JSUB) /= 0) THEN
         EL_REDUNDANT_PRES = EL_REDUNDANT_PRES + 1
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1521) JSUB,ESORT1(IELEM)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1521) JSUB,ESORT1(IELEM)
         ENDIF
      ENDIF

      PPNT(IELEM,JSUB) = IPPNT

! Set PTYPE for this element

      IF      (NAME(1:6) == 'PLOAD1') THEN
         PTYPE(IELEM) = '2'
      ELSE IF (NAME(1:6) == 'PLOAD2') THEN
         PTYPE(IELEM) = '1'
      ELSE IF (NAME(1:6) == 'PLOAD4') THEN
         IF      ((ETYPE(IELEM)(1:4) == 'TRIA') .OR. (ETYPE(IELEM)(1:5) == 'TETRA') .OR. (ETYPE(IELEM)(1:5) == 'PENTA')) THEN
            PTYPE(IELEM) = '3'
         ELSE IF ((ETYPE(IELEM)(1:4) == 'QUAD') .OR. (ETYPE(IELEM)(1:4) == 'HEXA')) THEN
            PTYPE(IELEM) = '4'
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1520 FORMAT(' *ERROR  1520: ELEMENT ',I8,' ON PLOAD2 ',I8,' DOES NOT EXIST OR IS OF WRONG TYPE FOR THE PRESSURE CARD:')

 1521 FORMAT(' *WARNING    : FOR INTERNAL SUBCASE NUMBER ',I8,' ELEMENT ',I8,' HAS PRESSURE DEFINED MORE THAN ONCE.',              &
                           ' LAST VALUE IN INPUT DECK WILL BE USED')

 1522 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE PRESSURE DATA ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE EPPUT

      END SUBROUTINE PRESSURE_DATA_PROC


      SUBROUTINE SLOAD_PROC

! Scalar load processor

! Transforms the input B.D. scalar load data to system force data in the SYS_LOAD array for dof's in the G set.
! File LINK1W was written when SLOAD B.D entries were read, with one record for each such pair of scalar point/load values.
! There are NSLOAD total no. records written to file LINK1W with each record containing:

!               SETID  = Load set ID on the SLOAD card
!               SPOINT = Scalar point where load acts
!               FMAG   = Scalar point load magnitude

! The process in creating array SYS_LOAD from this information is as follows:

!  (1) For each subcase (1 to NSUB):

!      (a) Generate LSID, RSID tables of load set ID's/scale factors for load sets for this subcase:

!          (  i) LLOADC is the max number of pairs of scale factors/load set ID's over all LOAD Bulk Data cards
!                including the pair defined by the set ID and overall scale factor on the LOAD Bulk Data card.
!                LSID and RSID are dimensioned 1 to LLOADC and:
!                   LSID(1) is always the load set ID requested in Case Control for this subcase.
!                   RSID(1) is always 1.0

!          ( ii) If the load set requested in Case Control is not a set ID from a LOAD Bulk data card then the
!                load must be on a separate FORCE/MOMENT in which case LSID(1) and RSID(1) are all that is needed
!                to define the load set contribution due to the FORCE/MOMENT cards for this subcase.

!          (iii) If the load set requested in Case Control is a set ID from a LOAD Bulk data card then the ramainder
!                (K = 2,LLOADC) of entries into LSID and RSID will be the pairs of load set ID's/scale factors from
!                that LOAD Bulk Data card (with RSID also multiplied by the overall scale factor on the LOAD Bulk data
!                card. The load set ID's are in array LOAD_SIDS(i,j) created when LOAD Bulk Data cards were read.
!                The scale factors are in array LOAD_FACS(i,j) also created when LOAD Bulk Data cards were read.
!                Note, there may not be as many as LLOADC pairs of set ID's/scale factors on a given LOAD Bulk Data
!                card since LLOADC is the max, from all LOAD Bulk Data cards, of pairs.
!                Thus, the entries in LSID from the last entry (for a given LOAD card) to LLOADC will be zero (LSID
!                was initialized to zero). This fact is used in a DO loop to EXIT when LSID(K) = 0

!      (b) For each record in LINK1W (1 to NSLOAD)

!          (  i) Read a record from file: (SETID, SPOINT, FMAG)

!          ( ii) Scan LSID and RSID to get the scale factor for the SLOAD components in SETID, if this
!                SLOAD's set ID is in LSID

!          (iii) Load these force values values into the system load array, SYS_LOAD


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06, L1W, LINK1W, L1W_MSG, L1WSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LLOADC, NGRID, NLOAD, NSLOAD, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  LOAD_SIDS, LOAD_FACS, SYS_LOAD, SUBLOD, GRID, GRID_ID

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  FILERR, FILE_CLOSE, READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SLOAD_PROC'
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for

      INTEGER(LONG)                   :: GDOF              ! G-set DOF no. for actual grid AGRID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: SETID             ! Load set ID read from record in file LINK1I
      INTEGER(LONG)                   :: LSID(LLOADC+1)     ! Array of load SID's, from a LOAD Bulk Data card, for one S/C
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM           ! Row no. in array TDOF corresponding to GDOF
      INTEGER(LONG)                   :: SPOINT            ! Scalra point read from a record of L1W (point where force acts)
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from file


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: FMAG              ! Force magnitude read from a L1W record (force on the SPOINT)
      REAL(DOUBLE)                    :: RSID(LLOADC+1)     ! Array of load magnitudes (for LSID set ID's) needed for one S/C
      REAL(DOUBLE)                    :: SCALE             ! Scale factor for a load

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! Now process forces into SYS_LOAD for each subcase.

      DO I=1,LLOADC                                         ! Initialize LSID, RSID arrays
         LSID(I) = 0
         RSID(I) = ZERO
      ENDDO

i_do2:DO I=1,NSUB                                          ! Loop through the S/C's

         IF (SUBLOD(I,1) == 0) THEN                        ! If no load for this S/C, CYCLE
            CYCLE i_do2
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(I,1)                             ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the FORCE or MOMENT card in file LINK1I that matches SUBLOD(I,1)
         DO J = 1,NLOAD                                    ! Then, the actual mag. will come from RSID(1) & FMAG
            IF (SUBLOD(I,1) == LOAD_SIDS(J,1)) THEN
k_do21:        DO K = 2,LLOADC
                  IF (LOAD_SIDS(J,K) == 0) THEN
                     EXIT k_do21
                  ELSE
                     NSID = K                              ! Note: NSID will not get larger than LLOADC
                     RSID(K) = LOAD_FACS(J,1)*LOAD_FACS(J,K)
                     LSID(K) = LOAD_SIDS(J,K)
                   ENDIF
                ENDDO k_do21
            ENDIF
         ENDDO

j_do22:  DO J=1,NSLOAD                                     ! Process SLOAD card info
                                                           ! (2-b-  i) Read a record from scratch - forces are in global coords
            READ(L1W,IOSTAT=IOCHK) SETID, SPOINT, FMAG
            IF (IOCHK /= 0) THEN
               REC_NO = J
               CALL READERR ( IOCHK, LINK1W, 'SLOAD FILE', REC_NO, OUNT )
               CALL FILE_CLOSE ( L1W, LINK1W, L1WSTAT )
               CALL OUTA_HERE ( 'Y' )                              ! Error reading scratch file, so quit
            ENDIF

            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
k_do221:    DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                  EXIT k_do221
               ENDIF
            ENDDO k_do221

            IF (FOUND == 'N') THEN                         ! Cycle back on J loop and read another force/moment card
               CYCLE j_do22
            ENDIF
                                                           ! Get GRID_ID_ROW_NUM, we checked it's existence earlier
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SPOINT, GRID_ID_ROW_NUM )

            IF (DABS(FMAG) < EPS1) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1513) 'SLOAD', SETID
               IF (SUPWARN == 'N') THEN                     ! Issue warning if all force components zero
                  WRITE(F06,1513) 'SLOAD', SETID
               ENDIF
            ENDIF

            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SPOINT, IGRID )
            ROW_NUM = TDOF_ROW_START(IGRID)

            CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )      ! (2-b-iii). Put forces and moments into SYS_LOAD array
            GDOF = TDOF(ROW_NUM,G_SET_COL_NUM)
            SYS_LOAD(GDOF,I) = SYS_LOAD(GDOF,I) + SCALE*FMAG

         ENDDO j_do22

         REWIND (L1W)                                       ! Need to read all of the FORCE/MOMENT records again for the next S/C
         READ(L1W,IOSTAT=IOCHK) XTIME
         IF (IOCHK /= 0) THEN
            REC_NO = 1
            CALL READERR ( IOCHK, LINK1W, 'SLOAD FILE', REC_NO, OUNT )
            CALL FILERR ( OUNT )
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDDO i_do2



      RETURN

! **********************************************************************************************************************************






 1513 FORMAT(' *WARNING    : ',A8,1X,I8,' HAS ZERO FORCE MAGNITUDE')

! **********************************************************************************************************************************

      END SUBROUTINE SLOAD_PROC


      SUBROUTINE TEMPERATURE_DATA_PROC

! Grid and element temperature data processor

! Processes model, grid and elem temperature data. The model (TEMPD), grid (TEMP), or element temperature (TEMPRB, TEMPP1) entries
! were written to file LINK1K when the Bulk Data was read. Here, that data is read and processed into arrays needed for temp load
! calculation.

! There are NTCARD records in file LINK1K. There is
!   1 record  for any model TEMPD         B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   1 record  for any grid  TEMP          B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   n records for any elem  TEMPRB/TEMPP1 B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!     where n are the no. of physical cards (parent + continuations)

! The process involves 3 passes through the data in file LINK1K and is as follows:

! Pass 1: Process model temperature cards: TEMPD          (see subr TEMPD_DATA_PROC     for details)
! Pass 2: Process grid  temperature cards: TEMP           (see subr GRID_TEMP_DATA_PROC for details)
! Pass 3: Process elem  temperature cards: TEMPRB, TEMPP1 (see subr ELEM_TEMP_DATA_PROC for details)

! The three passes are done in the order shown so that TEMPD cards will first set all grids to a value that can be overridden
! by a TEMP card. Finally, element temperature cards can override TEMP cards for the grids belonging to that element.
! During these 3 passes through LINK1K, the grid and element temperature data arrays are constructed:

!   GTEMP : Array of temperatures for all grids and subcases
!   ETEMP : Array of avg bulk temperatures for all elems and subcases
!   TPNT  : Pointer array for where, in array TDATA, that elem temperature data begins for an element
!   TDATA : Array of element temperature data from TEMPRB, TEMPP1 cards
!   CGTEMP: Char array that indicates whether or not a grid has a temperature defined, for a subcase
!   CETEMP: Char array that indicates whether or not a elem has a temperature defined, for a subcase

! At the end, array CETEMP is scanned to make sure that every elem has had a temperature defined
! either on a TEMPRB or TEMPP1 card directly, or indirectly on TEMPD/TEMP cards for each grid belonging to the element.
! If all elems have temperature defined, file LINK1K is closed and reopened to write the above temperature arrays
! for further processing in the element subr ELMDAT.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     L1K
      USE IOUNT1, ONLY                :  WRT_ERR,                            LINK1K
      USE IOUNT1, ONLY                :  WRT_ERR,                            L1K_MSG
      USE SCONTR, ONLY                :  DATA_NAM_LEN, NELE, NGRID, NTDAT, NTSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CETEMP, CETEMP_ERR, CGTEMP, CGTEMP_ERR, ETEMP, GTEMP, TDATA, TPNT, GRID_ID, ESORT1, ETYPE,&
                                         SCNUM, SUBLOD, eid
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN, READERR
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELGP

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TEMPERATURE_DATA_PROC'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes
      CHARACTER(  1*BYTE)             :: NOTE              ! Used to indicate whether or not to print out a message
      CHARACTER(  1*BYTE)             :: TEMP_ELM          ! Descriptor of how elem temp is defined

      INTEGER(LONG)                   :: ELID_ERR          ! Output from subr ELEM_TEMP_DATA_PROC (elem ID's undefined).
      INTEGER(LONG)                   :: GID_ERR           ! Output from subr GRID_TEMP_DATA_PROC (grid ID's undefined).
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERRT             ! Cum. count of grids with no temp defined after all data processed
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: READ_ERR  = 0     ! Count of read errors when temperature data file is read
      INTEGER(LONG)                   :: TCASE1(NSUB)      ! TCASE1(I) gives the internal thermal case number for internal S/C I
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE1(1-5) = 0, 0, 1, 0, 2
      INTEGER(LONG)                   :: TCASE2(NSUB)      ! TCASE2(I) gives the internal subcase no. for internal thermal case I
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE2(1-5) = 3, 5, 0, 0, 0




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Generate TCASE1, TCASE2 arrays.

      J = 0
      DO I=1,NSUB
         TCASE1(I) = 0
      ENDDO
      DO I = 1,NSUB
         IF (SUBLOD(I,2) > 0) THEN
            J = J + 1
            TCASE1(I) = J
         ENDIF
      ENDDO

      DO I = 1,NSUB
         TCASE2(I) = 0
      ENDDO
      J = 0
      DO I = 1,NSUB
         IF (SUBLOD(I,2) /= 0) THEN
            J = J + 1
            TCASE2(J) = I
         ENDIF
      ENDDO

! **********************************************************************************************************************************
! (1) Pass 1: Read and process TEMPD cards from L1K and set G.P. temperatures

      CALL TEMPD_DATA_PROC ( TCASE1, OUNT, READ_ERR )

! **********************************************************************************************************************************
! (2) Pass 2: Read and process TEMP cards from L1K and set G.P. temperatures

      CALL GRID_TEMP_DATA_PROC ( TCASE1, OUNT, READ_ERR, GID_ERR )

! **********************************************************************************************************************************
! (3) Pass 3: Read and process the element temperature cards, TEMPRB, TEMPP1 and set elem temperatures

      CALL ELEM_TEMP_DATA_PROC ( TCASE1, OUNT, READ_ERR, ELID_ERR )

! **********************************************************************************************************************************
! (4) If no errors, check thermal data for completeness

      IF ((READ_ERR /= 0) .OR. (GID_ERR /= 0) .OR. (ELID_ERR /= 0)) THEN

         IERRT = READ_ERR + GID_ERR + ELID_ERR
         WRITE(ERR,9996) SUBR_NAME,IERRT
         WRITE(F06,9996) SUBR_NAME,IERRT
         CALL OUTA_HERE ( 'Y' )                                    ! Quit due to errors reading/processing temperature data

      ELSE                                                 ! Check that all elems have a temperature defined

         IERRT = 0
         DO J=1,NTSUB
            DO I=1,NELE
               IF ((ETYPE(I)(1:4) == 'ELAS') .OR. (ETYPE(I)(1:4) == 'BUSH') .OR. (ETYPE(I)(1:6) == 'PLOTEL')) THEN
                  CETEMP(I,J) = 'N'
               ELSE
                  IF (CETEMP(I,J) /= 'E') THEN
                     CALL ELEM_TEMP_CHK ( TCASE2(J), J, I, TEMP_ELM )
                     IF (TEMP_ELM /= CETEMP_ERR) THEN
                        CETEMP(I,J) = TEMP_ELM
                     ELSE
                        IERRT = IERRT + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

! Make sure that CETEMP has no error terms

! Debug output of grid and element temperature data

      IF ((DEBUG(8) == 1) .OR. (DEBUG(8) == 3)) THEN
         NOTE = 'N'
         WRITE(F06,99971)
         DO J=1,NTSUB
            DO I=1,NGRID
               IF (CGTEMP(I,J) == CGTEMP_ERR) THEN
                  WRITE(F06,99972) SCNUM(TCASE2(J)),GRID_ID(I),CGTEMP(I,J)
                  NOTE = 'Y'
               ELSE
                  WRITE(F06,99973) SCNUM(TCASE2(J)),GRID_ID(I),GTEMP(I,J),CGTEMP(I,J)
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF (NOTE == 'Y') THEN
            WRITE(F06,99974) CGTEMP_ERR
         ENDIF
         WRITE(F06,*)
      ENDIF

      IF ((DEBUG(8) == 2) .OR. (DEBUG(8) == 3)) THEN
         NOTE = 'N'
         WRITE(F06,99991)
         DO J=1,NTSUB
            DO I=1,NELE
               IF (CETEMP(I,J) == CETEMP_ERR) THEN
                  WRITE(F06,99972) SCNUM(TCASE2(J)),ESORT1(I),CETEMP(I,J)
                  NOTE = 'Y'
               ELSE
                  WRITE(F06,99973) SCNUM(TCASE2(J)),ESORT1(I),ETEMP(I,J),CETEMP(I,J)
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF (NOTE == 'Y') THEN
            WRITE(F06,99994) CETEMP_ERR
         ENDIF
         WRITE(F06,*)
      ENDIF

! Quit if there were errors setting element temperatures

      IF (IERRT > 0) THEN
         WRITE(ERR,9998) IERRT
         WRITE(F06,9998) IERRT
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! (5) Write processed temperature data to L1K

! First close and delete L1K file

      CALL FILE_CLOSE ( L1K, LINK1K, 'DELETE' )

! Open L1K for write:

      CALL FILE_OPEN ( L1K, LINK1K, OUNT, 'REPLACE', L1K_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

      DATA_SET_NAME = 'TPNT'
      WRITE(L1K) DATA_SET_NAME
      WRITE(L1K) NELE
      WRITE(L1K) NTSUB
      DO I = 1,NELE
         DO J = 1,NTSUB
            WRITE(L1K) TPNT(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'TDATA'
      WRITE(L1K) DATA_SET_NAME
      WRITE(L1K) NTDAT
      DO I = 1,NTDAT
         WRITE(L1K) TDATA(I)
      ENDDO

      DATA_SET_NAME = 'GTEMP'
      WRITE(L1K) DATA_SET_NAME
      WRITE(L1K) NGRID
      WRITE(L1K) NTSUB
      DO I = 1,NGRID
         DO J = 1,NTSUB
            WRITE(L1K) GTEMP(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'CGTEMP'
      WRITE(L1K) DATA_SET_NAME
      WRITE(L1K) NGRID
      WRITE(L1K) NTSUB
      DO I = 1,NGRID
         DO J = 1,NTSUB
            WRITE(L1K) CGTEMP(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'CETEMP'
      WRITE(L1K) DATA_SET_NAME
      WRITE(L1K) NELE
      WRITE(L1K) NTSUB
      DO I = 1,NGRID
         DO J = 1,NTSUB
            WRITE(L1K) CGTEMP(I,J)
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
 9998 FORMAT(/,' Processing terminated due to',I8,' element temperatures undefined. see output file for details')

 9996 FORMAT(/,' Processing terminated in subroutine ',A,' due to above ',I8,' ERRORS')

99971 FORMAT(/,' OUTPUT FROM SUBROUTINE TEMPERATURE_DATA_PROC:',//,'   SUBCASE      GRID         TEMPERATURE    CGTEMP')

99972 FORMAT(1X,2I10,20X    ,9X,A1)

99973 FORMAT(1X,2I10,1ES20.6,9X,A1)

99974 FORMAT(1X,A,' INDICATES NO TEMPERATURE DEFINED FOR THIS SUBCASE FOR THIS GRID')

99991 FORMAT(/,' OUTPUT FROM SUBROUTINE TEMPERATURE_DATA_PROC:',//,'   SUBCASE   ELEMENT         TEMPERATURE    CETEMP')

99994 FORMAT(1X,A,' INDICATES NO TEMPERATURE DEFINED FOR THIS SUBCASE FOR THIS ELEMENT')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE TEMPD_DATA_PROC ( TCASE1, OUNT, IERR )

! Overall model temperature data processor

! Processes data on TEMPD B.D. cards. These card images were written to file LINK1K when the Bulk Data was read.
! Here, that data is read and processed into arrays needed for temp load calculation.

! There are NTCARD records in file LINK1K. There is
!   1 record  for any model TEMPD         B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   1 record  for any grid  TEMP          B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   n records for any elem  TEMPRB/TEMPP1 B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!     where n are the no. of physical cards (parent + continuations)

! This 1st of 3 passes through file LINK1K processes TEMPD cards as follows:

!  (1) Read a record from file LINK1K. If it is TEMPD and the set ID (SID) matches a subcase request, then
!      ( i) Load the temperature value into array GTEMP for all grids for this subcase
!      (ii) Write 'D' to array CGTEMP for all grids for this subcase (grid temperature defined) to indicate that the
!           temperature is defined via a default (hence 'D') TEMPD card.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, NGRID, NSUB, NTCARD, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CGTEMP, GTEMP, SUBLOD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TEMPD_DATA_PROC'
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters in CARD
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! B.D elem temp card image read from file LINK1K

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), INTENT(IN)       :: TCASE1(NSUB)      ! TCASE1(I) gives the internal thermal case number for internal S/C I
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE1(1-5) = 0, 0, 1, 0, 2
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG), INTENT(INOUT)    :: IERR              ! Cum. count of errors as we read, and check cards from file LINK1K
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SID               ! Thermal load set ID read from an elem temperature B.D. card


      REAL(DOUBLE)                    :: RTEMP             ! Real value of a temperature on a TEMPD or TEMP B.D. card



! **********************************************************************************************************************************
      DO I=1,NTCARD
         READ(L1K,IOSTAT=IOCHK) CARD
         IF (IOCHK /= 0) THEN
            REC_NO = I
            CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
            IERR = IERR + 1
            CYCLE
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF ((CARD(1:6) == 'TEMPD ') .OR. (CARD(1:6) == 'TEMPD*')) THEN
            DO J = 1,4
               IF (JCARD(2*J) == '        ') THEN
                  CYCLE
               ELSE
                  READ(JCARD(2*J)  ,'(I8  )') SID
                  READ(JCARD(2*J+1),'(F16.0)') RTEMP
                  DO K = 1,NSUB
                     IF (SUBLOD(K,2) == SID) THEN
                        DO L = 1,NGRID
                            GTEMP(L,TCASE1(K)) = RTEMP
                           CGTEMP(L,TCASE1(K)) = 'D'
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE TEMPD_DATA_PROC

! ##################################################################################################################################

      SUBROUTINE GRID_TEMP_DATA_PROC ( TCASE1, OUNT, READ_ERR, GID_ERR )

! Grid temperature data processor

! Processes data on TEMP B.D. cards. These card images were written to file LINK1K when the Bulk Data was read.
! Here, that data is read and processed into arrays needed for temp load calculation.

! There are NTCARD records in file LINK1K. There is
!   1 record  for any model TEMPD         B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   1 record  for any grid  TEMP          B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   n records for any elem  TEMPRB/TEMPP1 B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!     where n are the no. of physical cards (parent + continuations)

! This 2nd of 3 passes through file LINK1K processes TEMP cards as follows (after REWINDing file LINK1K:

!  (1) Read a record from file LINK1K. If it is TEMP and the set ID (SID) matches a subcase request, then
!      ( i) Load the temperature value into array GTEMP for the grids on this TEMP card for this subcase
!      (ii) Write 'G' to array CGTEMP for these grids for this subcase (grid temperature defined) to indicate that the
!           temperature is defined via a grid (hence 'G') TEMP card. If that grid already has a 'G' issue a warning.



      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, NGRID, NSUB, NTCARD, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  CGTEMP, GTEMP, GRID_ID, SUBLOD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GRID_TEMP_DATA_PROC'
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters in CARD
      CHARACTER( 8*BYTE), PARAMETER   :: NAME = 'GRID    ' ! name for output error purposes
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! B.D elem temp card image read from file LINK1K

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), INTENT(IN)       :: TCASE1(NSUB)      ! TCASE1(I) gives the internal thermal case number for internal S/C I
      INTEGER(LONG), INTENT(INOUT)    :: READ_ERR          ! Count of read errors when temperature data file is read
      INTEGER(LONG), INTENT(OUT)      :: GID_ERR           ! Count of errors due to a grid ID not defined
      INTEGER(LONG)                   :: AGRID             ! Actual grid ID
      INTEGER(LONG)                   :: GD_REDUNDANT_TEMP ! Count of warning messages when elements have redundant temperatures
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SID               ! Thermal load set ID read from an elem temperature B.D. card
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE1(1-5) = 0, 0, 1, 0, 2


      REAL(DOUBLE)                    :: RTEMP             ! Real value of a temperature on a TEMPD or TEMP B.D. card



! **********************************************************************************************************************************
      GID_ERR = 0
      GD_REDUNDANT_TEMP = 0                                ! Initialize warning indicator for redundant grid temp definition

      REWIND (L1K)
      READ(L1K,IOSTAT=IOCHK) XTIME
      IF (IOCHK /= 0) THEN
         REC_NO = 1
         CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
         CALL OUTA_HERE ( 'Y' )                            ! Cannot read STIME from temperature data file, so quit
      ENDIF

      DO I=1,NTCARD
         READ(L1K,IOSTAT=IOCHK) CARD
         IF (IOCHK /= 0) THEN
            REC_NO = I+1
            CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
            READ_ERR = READ_ERR + 1
            CYCLE
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF ((CARD(1:5) == 'TEMP ') .OR. (CARD(1:5) == 'TEMP*')) THEN
            READ(JCARD(2),'(I8)') SID
            DO J = 1,3
               IF (JCARD(J*2+1)(1:) == ' ') CYCLE
               READ(JCARD(J*2+1),'(I8)'  ) AGRID
               READ(JCARD(J*2+2),'(F16.0)') RTEMP
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  WRITE(ERR,1822) 'GRID ', AGRID, 'TEMP ', SID
                  WRITE(F06,1822) 'GRID ', AGRID, 'TEMP ', SID
                  GID_ERR   = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ELSE
                  DO K = 1,NSUB
                    IF (SID == SUBLOD(K,2)) THEN
                        GTEMP(GRID_ID_ROW_NUM,TCASE1(K)) = RTEMP
                        IF (CGTEMP(GRID_ID_ROW_NUM,TCASE1(K)) == 'G') THEN
                           GD_REDUNDANT_TEMP = GD_REDUNDANT_TEMP + 1
                           WARN_ERR = WARN_ERR + 1
                           WRITE(ERR,1531) NAME,AGRID
                           IF (SUPWARN == 'N') THEN
                              WRITE(F06,1531) NAME,AGRID
                           ENDIF
                        ENDIF
                        CGTEMP(GRID_ID_ROW_NUM,TCASE1(K)) = 'G'
                    ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (GD_REDUNDANT_TEMP > 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1527)
         IF (SUPWARN == 'U') THEN
            WRITE(F06,1527)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1527 FORMAT(' *WARNING    : CHECK ERR OUTPUT FILE FOR WARNING MESSAGES REGARDING REDUNDANT TEMPERATURE DEFINITION')

 1531 FORMAT(' *WARNING    : ',A8,1X,I8,' HAS TEMPERATURE DEFINED MORE THAN ONCE. LAST VALUE IN INPUT DECK WILL BE USED')

! **********************************************************************************************************************************

      END SUBROUTINE GRID_TEMP_DATA_PROC

! ##################################################################################################################################

      SUBROUTINE ELEM_TEMP_DATA_PROC ( TCASE1, OUNT, READ_ERR, ELID_ERR )

! Element temperature data processor

! Processes data on TEMPRB, TEMPP1 B.D. cards. These card images were written to file LINK1K when the Bulk Data was read.
! Here, that data is read and processed into arrays needed for temp load calculation.

! There are NTCARD records in file LINK1K. There is
!   1 record  for any model TEMPD         B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   1 record  for any grid  TEMP          B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!   n records for any elem  TEMPRB/TEMPP1 B.D. card whose set ID (SID) matches that of a subcase TEMP = SID request
!     where n are the no. of physical cards (parent + continuations)

! This 3rd of 3 passes through file LINK1K processes TEMPRB, TEMPP1 cards as follows (after REWINDing file LINK1K:

!  (1) In an outer DO loop, cycle through records from file LINK1K until a TEMPRB or TEMPP1 parent card is found. If one
!      is found, load the temperature data into arrays TDATA and get subr ETPUT to write data to array TPNT

!  (2) In an inner DO loop, read continuation cards for the above parent. If any are found, they will contain additional
!      element ID's in one of several formats (checked when B.D. was read). Call ETPUT to load add these elements to
!      array TPNT.

!      The calls to ETPUT involve up to 3 different sets of circumstances:
!      ( i ) The additional element ID's on the cont. card are all separate element ID's, or
!      (iia) The additional element ID's on the cont. card are specified in an EID1 THRU EID2 format in fields 2-4
!      (iib) The additional element ID's on the cont. card are specified in an EID1 THRU EID2 format in fields 5-7
!      In case ( i ) an error is given if any elem ID's are undefined.
!      In cases (iia) and (iib) the elem ID's are specified by a range and some, or all, of these do not have to exist.
!      Thus, no error is given for these cases. If some elems wind up with no temperature defined at all, them this will
!      be detected in the calling routine


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, LTDAT, MTDAT_TEMPRB, MTDAT_TEMPP1,      &
                                         NTCARD, NTDAT,  &
                                         NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  TWO
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  TDATA

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELEM_TEMP_DATA_PROC'
      CHARACTER( 1*BYTE)              :: IELEM_ERR         ! Flag to decide whether to CYCLE a loop
      CHARACTER( 1*BYTE)              :: EFLAG             ! Flag to decide whether a situation in subr ETPUT is an error or not
      CHARACTER(LEN=JCARD_LEN)        :: CARD_NAME         ! The 1st field of the parent card of an entry in file L1K
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters in CARD
      CHARACTER(LEN=JCARD_LEN)        :: OLDTAG            ! Field 10 (continuation field) on CARD
      CHARACTER( 8*BYTE)              :: THRU              ! ='Y' if THRU option used on TEMPRB, TEMPP1 continuation card
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Variable to test whether "THRU" option was used on B.D. temp card
      CHARACTER(BD_ENTRY_LEN)         :: CARD              ! B.D elem temp card image read from file LINK1K

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), INTENT(IN)       :: TCASE1(NSUB)      ! TCASE1(I) gives the internal thermal case number for internal S/C I
      INTEGER(LONG), INTENT(INOUT)    :: READ_ERR          ! Count of read errors when temperature data file is read
      INTEGER(LONG), INTENT(OUT)      :: ELID_ERR          ! Count of errors due to an elem ID not defined
      INTEGER(LONG)                   :: EID               ! Actual element ID
      INTEGER(LONG)                   :: EID1,EID2         ! The 2 actual elem ID's in "EID1 THRU EID2" on elem temp B.D. card
      INTEGER(LONG)                   :: EL_REDUNDANT_TEMP ! Count of warning messages when elements have redundant temperatures
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: ICRD              ! Count of the records (cards) read from file LINK1K
      INTEGER(LONG)                   :: IELEM             ! Internal element number for the actual element ID, EID
      INTEGER(LONG)                   :: IFIELD, JFIELD    ! DO loop limits when reading temperature data from a B.D. temp card
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: ITPNT             ! Index in array TPNT (a pointer to where temp data for EID starts)
      INTEGER(LONG)                   :: NFIELD            ! No. fields on TEMPRB, TEMPP1 B.D. cards that have temperature data
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SID               ! Thermal load set ID read from an elem temperature B.D. card
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE1(1-5) = 0, 0, 1, 0, 2


      REAL(DOUBLE)                    :: TB1               ! Bulk temperature from TEMPRB card
      REAL(DOUBLE)                    :: TB2               ! Bulk temperature from TEMPRB card
      REAL(DOUBLE)                    :: TE_BULK           ! Element bulk temperature



! **********************************************************************************************************************************
      ELID_ERR = 0
      EL_REDUNDANT_TEMP = 0                                ! Initialize warning indicator for redundant elem temp definition

      REWIND (L1K)                                          ! L1K has been read from by subr GRID_TEMP_DATA_PROC earlier
      READ(L1K,IOSTAT=IOCHK) XTIME
      IF (IOCHK /= 0) THEN
         REC_NO = 1
         CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
         CALL OUTA_HERE ( 'Y' )                            ! Cannot read STIME from temperature data file, so quit
      ENDIF

      ICRD = 0                                             ! ICRD will count records read from file LINK1K

ntcrd:DO                                                   ! Top of loop for reading parent CARD

         ICRD = ICRD + 1
         IF (ICRD <= NTCARD) THEN                          ! Read and process CARD as long as card count <= NTCARD
            READ(L1K,IOSTAT=IOCHK) CARD
            CARD_NAME = CARD(1:8)
            IF (IOCHK /= 0) THEN
               REC_NO = ICRD
               CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
               READ_ERR = READ_ERR + 1
               CYCLE ntcrd
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
                                                           ! If we have TEMPRB or TEMPP1, set NFIELD
            IF      ((CARD(1:7) == 'TEMPRB ') .OR. (CARD(1:7) == 'TEMPRB*')) THEN
               NFIELD = MTDAT_TEMPRB
               READ(JCARD(4),'(F16.0)') TB1
               READ(JCARD(5),'(F16.0)') TB2
               TE_BULK = (TB1 + TB2)/TWO
            ELSE IF ((CARD(1:7) == 'TEMPP1 ') .OR. (CARD(1:7) == 'TEMPP1*')) THEN
               NFIELD = MTDAT_TEMPP1
               READ(JCARD(4),'(F16.0)') TE_BULK
            ELSE                                           ! otherwise, CYCLE back to read another parent
               CYCLE ntcrd
            ENDIF

            ITPNT = NTDAT + 1                              ! Set index for pointer array, TPNT

            READ(JCARD(2),'(I8)') SID                      ! Get thermal load SID, EID
            READ(JCARD(3),'(I8)') EID

            IF ((NTDAT + NFIELD) > LTDAT) THEN             ! Check for overflow in TDATA
               WRITE(ERR,1525) SUBR_NAME,LTDAT
               WRITE(F06,1525) SUBR_NAME,LTDAT
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Coding error (dim of array TDATA too small), so quit
            ENDIF

            IFIELD = 4                                     ! Put temperature data into array TDATA
            JFIELD = IFIELD - 1 + NFIELD
            DO J = IFIELD,JFIELD
               NTDAT = NTDAT + 1
               READ(JCARD(J),'(F16.0)') TDATA(NTDAT)
            ENDDO

            EFLAG = 'Y'                                    ! First process EID on parent card into TPNT array
            CALL ETPUT ( CARD_NAME, EFLAG, EID, SID, TCASE1, ITPNT, TE_BULK, IELEM, EL_REDUNDANT_TEMP )
            IF (IELEM == -1) THEN
               ELID_ERR = ELID_ERR + 1
            ENDIF

! Process continuation cards, if any. If we get read error on CARD, go back to top of main loop to look for new parent card.

cont_cards: DO                                             ! Top of loop for reading optional cont. CARD's for the parent CARD
               IF (ICRD == NTCARD) EXIT cont_cards
               OLDTAG = JCARD(10)
               READ(L1K,IOSTAT=IOCHK) CARD
               CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
               IF (IOCHK /= 0) THEN                        ! Error reading temp. data file, so set error & cycle to read another
                  REC_NO = ICRD
                  CALL READERR ( IOCHK, LINK1K, L1K_MSG, REC_NO, OUNT )
                  READ_ERR = READ_ERR + 1
                  CYCLE ntcrd
               ENDIF
               IF (OLDTAG == JCARD(1)) THEN
                  ICRD = ICRD + 1
               ELSE                                        ! Next card was not a continuation, so cycle back to read another parent
                  BACKSPACE(L1K)
                  CYCLE ntcrd
               ENDIF

! First check for the 2 options on specifying data on continuation card. Either all data are EID's or the THRU option is used
! in which case field 3 and/or 6 will have "THRU". Notice that EFLAG = 'N' when ETPUT is called for the case where EID1 THRU EID2
! is on temp card. This is so that all elements in the range EID1 thru EID2 do not have to exist to avoid an error being flagged.

               TOKEN = JCARD(3)(1:8)                       ! Only send the 1st 8 chars of this JCARD. It has been left justified
               CALL TOKCHK ( TOKEN, TOKTYP )
               THRU = 'N'
               IF (TOKTYP == 'THRU    ') THEN
                  THRU = 'Y'
               ENDIF

               IF (THRU == 'N') THEN
                  IELEM_ERR = 'N'
                  EFLAG = 'Y'                              ! Any elem ID's read in fields 2-9 must be valid
                  DO J=2,9
                     IF (JCARD(J) == '        ') EXIT
                     READ(JCARD(J),'(I8)') EID             ! Process individual EID's on continuation card into TPNT array
                     CALL ETPUT ( CARD_NAME, EFLAG, EID, SID, TCASE1, ITPNT, TE_BULK, IELEM, EL_REDUNDANT_TEMP )
                     IF (IELEM == -1) THEN                 ! ETPUT didn't find elem ID = EID
                        ELID_ERR = ELID_ERR + 1
                        IELEM_ERR = 'Y'
                     ENDIF
                  ENDDO
                  IF (IELEM_ERR == 'Y') THEN
                     CYCLE cont_cards
                  ENDIF
               ELSE
                  READ(JCARD(2),'(I8)') EID1
                  READ(JCARD(4),'(I8)') EID2
                  EFLAG = 'N'                              ! Some, or all, elem ID's can be missing from range EID1 -> EID2
                  DO J=EID1,EID2                           ! Process "EID1 THRU EID2" (fields 2-4) on cont card into TPNT array
                     CALL ETPUT ( CARD_NAME, EFLAG, J, SID, TCASE1, ITPNT, TE_BULK, IELEM, EL_REDUNDANT_TEMP )
                  ENDDO

                  TOKEN = JCARD(6)(1:8)                      ! Only send the 1st 8 chars of this JCARD. It has been left justified
                  CALL TOKCHK ( TOKEN, TOKTYP )
                  IF (TOKTYP == 'THRU    ') THEN
                     READ(JCARD(5),'(I8)') EID1
                     READ(JCARD(7),'(I8)') EID2
                     EFLAG = 'N'                           ! Some, or all, elem ID's can be missing from range EID1 -> EID2
                     DO J=EID1,EID2                        ! Process "EID1 THRU EID2" (fields 5-7) on cont card into TPNT array
                        CALL ETPUT ( CARD_NAME, EFLAG, J, SID, TCASE1, ITPNT, TE_BULK, IELEM, EL_REDUNDANT_TEMP )
                     ENDDO
                  ENDIF
               ENDIF

            ENDDO cont_cards

         ELSE

            EXIT                                           ! EXIT loop when ICRD exceeds NTCARD

         ENDIF

      ENDDO ntcrd

      IF (EL_REDUNDANT_TEMP > 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1527)
         IF (SUPWARN == 'Y') THEN
            WRITE(F06,1527)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1525 FORMAT(' *ERROR  1525: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MUCH ELEMENT TEMPERATURE DATA. LIMIT = ',I12)

 1527 FORMAT(' *WARNING    : CHECK ERR OUTPUT FILE FOR WARNING MESSAGES REGARDING REDUNDANT TEMPERATURE DEFINITION')

! **********************************************************************************************************************************

      END SUBROUTINE ELEM_TEMP_DATA_PROC

! ##################################################################################################################################

      SUBROUTINE ELEM_TEMP_CHK( ISCNO, JTCOL, IELEM, TEMP_ELM )

! Called by subr TEMPERATURE_PROC for all elements that do not have their temperature defined via a B.D.
! element temperature card (TEMPRB, TEMPP1) to determine if the element's grids have a temperature
! defined via a TEMPD or TEMP entry. If it does, then an avg bulk temperature is calculated from the grid
! tamperature data. This avg bulk temp is loaded into array ETEMP. Output variable TEMP_ELM is a character
! which gives info about how the elem temp was arrived at.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MELGP, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ETEMP, CGTEMP, CETEMP_ERR, GTEMP, EDAT, EPNT, ETYPE, GRID_ID, SCNUM,                      &
                                         AGRID, ELGP, EID, TYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELEM_TEMP_CHK'

      CHARACTER( 1*BYTE), INTENT(OUT) :: TEMP_ELM          ! Descriptor of how temperature is defined for elem IELEM:
!                                                             = 'D' if all grids have temperature defined on TEMPD
!                                                             = 'G' if all grids have temperature defined on TEMP
!                                                             = 'B' if grids have temperature defined on TEMPD and TEMP
!                                                             = CETEMP_ERR if one or more grids have no temperature defined
      CHARACTER( 1*BYTE)              :: TEMP_DEF          ! = 'Y' if any grid for elem IELEM has TEMPD definition
      CHARACTER( 1*BYTE)              :: TEMP_GRD          ! = 'Y' if any grid for elem IELEM has TEMP  definition
      CHARACTER( 1*BYTE)              :: TEMP_ERR          ! = 'Y' if any grid for elem IELEM has no temperature defined

      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where AGRID(I) is found
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG), INTENT(IN)       :: IELEM             ! Internal element number for a specific actual element ID
      INTEGER(LONG), INTENT(IN)       :: ISCNO             ! Internal subcase number
      INTEGER(LONG), INTENT(IN)       :: JTCOL             ! Col in thermal array CGTEMP, CETEMP for internal subcase no. ISCNO




! **********************************************************************************************************************************
! Initialize outputs

      TEMP_ELM = ' '

      TEMP_DEF = 'N'
      TEMP_GRD = 'N'
      TEMP_ERR = 'N'

      EID  = EDAT(EPNT(IELEM))
      TYPE = ETYPE(IELEM)

      CALL GET_ELGP ( IELEM )

      ETEMP(IELEM,JTCOL) = ZERO
      DO I=1,ELGP
         AGRID(I) = EDAT(EPNT(IELEM)+I+1)
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(I), GRID_ID_ROW_NUM )
         IF      (CGTEMP(GRID_ID_ROW_NUM,JTCOL) == 'D') THEN
            TEMP_DEF = 'Y'
         ELSE IF (CGTEMP(GRID_ID_ROW_NUM,JTCOL) == 'G') THEN
            TEMP_GRD = 'Y'
         ELSE
            TEMP_ERR = 'Y'
            WRITE(ERR,1526) EID, SCNUM(ISCNO), AGRID(I)
            WRITE(F06,1526) EID, SCNUM(ISCNO), AGRID(I)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         ETEMP(IELEM,JTCOL) = ETEMP(IELEM,JTCOL) + GTEMP(GRID_ID_ROW_NUM,JTCOL)
      ENDDO
      ETEMP(IELEM,JTCOL) = ETEMP(IELEM,JTCOL)/ELGP

      TEMP_ELM = ' '
      IF (TEMP_ERR == 'Y') THEN
         TEMP_ELM = CETEMP_ERR
      ELSE IF ((TEMP_DEF == 'Y') .AND. (TEMP_GRD == 'N')) THEN
         TEMP_ELM = 'D'
      ELSE IF ((TEMP_DEF == 'N') .AND. (TEMP_GRD == 'Y')) THEN
         TEMP_ELM = 'G'
      ELSE IF ((TEMP_DEF == 'Y') .AND. (TEMP_GRD == 'Y')) THEN
         TEMP_ELM = 'B'
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1526 FORMAT(' *ERROR  1526: CANNOT CALC TEMPERATURE FOR ELEM ',I8,' SUBCASE ',I8,'. GRID ',I8,' HAS NO TEMPERATURE DEFINED')

! **********************************************************************************************************************************

      END SUBROUTINE ELEM_TEMP_CHK

! ##################################################################################################################################

      SUBROUTINE ETPUT( CARD_NAME, EFLAG, EID, SID, TCASE1, ITPNT, TE_BULK, IELEM, EL_REDUNDANT_TEMP )

! Elem temperature routine - generates the TPNT(i,j) array that points into array TDATA to indicate where temperature
! data for an element begins for a thermal subcase. Also loads avg bulk temp for this elem into array ETEMP and sets
! array CETEMP to state that this elem temp was specified on an elem temp card

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NELE, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETEMP, ESORT1, CETEMP, TPNT, TYPE, SUBLOD
      USE PARAMS, ONLY                :  SUPWARN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ETPUT'
      CHARACTER( 1*BYTE), INTENT(IN)  :: EFLAG             ! Flag used to decide whether it is an error if an elem in the range
      CHARACTER( 8*BYTE), PARAMETER   :: NAME = 'ELEMENT ' ! name for output error purposes
                                                           ! IEL_LO to IEL_HI does not exist
      CHARACTER(LEN=*), INTENT(IN)    :: CARD_NAME         ! The 1st field of the parent card of an entry in file L1K

      INTEGER(LONG), INTENT(IN)       :: ITPNT             ! Index in array TPNT (a pointer to where temp data for EID starts)
      INTEGER(LONG), INTENT(IN)       :: EID               ! Actual element ID
      INTEGER(LONG), INTENT(IN)       :: SID               ! Set ID from the elem temperature card being processed
      INTEGER(LONG), INTENT(OUT)      :: IELEM             ! Internal element number for the actual element ID, EID
      INTEGER(LONG), INTENT(INOUT)    :: EL_REDUNDANT_TEMP ! Count of warning messages when elements have redundant temperatures
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG), INTENT(IN)       :: TCASE1(NSUB)      ! TCASE1(I) gives the internal thermal case number for internal S/C I
!                                                            If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                            thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                            TCASE1(1-5) = 0, 0, 1, 0, 2


      REAL(DOUBLE) , INTENT(IN)       :: TE_BULK           ! Bulk temperature from element temperature B.D. card



! **********************************************************************************************************************************
! Initialize outputs

      IELEM = 0

! First convert EID to interal value

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, EID, IELEM )

! If EID not found then: error if EFLAG = 'Y', or return if EFLAG /= 'Y'

      IF (IELEM == -1) THEN
         IF (EFLAG == 'Y') THEN
            WRITE(ERR,1530) EID, CARD_NAME, SID, TYPE
            WRITE(F06,1530) EID, CARD_NAME, SID, TYPE
            FATAL_ERR = FATAL_ERR + 1
         ENDIF

         RETURN
      ENDIF

! No error, so put pointer into TPNT(i,j) and give warning if the
! element has had a temperature defined previously

      DO I = 1,NSUB
         IF (SID == SUBLOD(I,2)) THEN
            IF (TPNT(IELEM,TCASE1(I)) /= 0) THEN
               EL_REDUNDANT_TEMP = EL_REDUNDANT_TEMP + 1
!xx            FATAL_ERR = FATAL_ERR + 1
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1531) NAME,ESORT1(IELEM)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1531) NAME,ESORT1(IELEM)
               ENDIF
            ENDIF
            TPNT(IELEM,TCASE1(I))   = ITPNT
            ETEMP(IELEM,TCASE1(I))  = TE_BULK
            CETEMP(IELEM,TCASE1(I)) = 'E'
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1530 FORMAT(' *ERROR  1530: ELEMENT ',I8,' ON ',A,I8,' DOES NOT EXIST OR IS OF WRONG TYPE = "',A,'" FOR THE TEMPERATURE CARD')

 1531 FORMAT(' *WARNING    : ',A8,1X,I8,' HAS TEMPERATURE DEFINED MORE THAN ONCE. LAST VALUE IN INPUT DECK WILL BE USED')

! **********************************************************************************************************************************

      END SUBROUTINE ETPUT

      END SUBROUTINE TEMPERATURE_DATA_PROC


      SUBROUTINE YS_ARRAY

! Process enforced displacement data in file LINK1H and write the enforced displacement array, YSe, for use in subsequent LINK's.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,    L1H
      USE IOUNT1, ONLY                :  WRT_ERR,                           LINK1H
      USE IOUNT1, ONLY                :  WRT_ERR,                           L1H_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFSE, NGRID
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE COL_VECS, ONLY              :  YSe

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  FILERR, READERR
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'YS_ARRAY'

      INTEGER(LONG)                   :: COMP              ! Displ component number read from file LINK1H
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row num in array GRID_ID of the actual grid that the YS displ is at
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR1H            ! Count of errors as file LINK1H is read
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM_START     ! Row no. in array TDOF where data begins for GRID_ID(GRID_ID_ROW_NUM)
      INTEGER(LONG)                   :: SE_SET_COL_NUM    ! Col no., in TDOF array, of the SE-set DOF list
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row num in array TDOF for DOF corresponding to GRID_ID_ROW_NUM, COMP
      INTEGER(LONG)                   :: YSDOF             ! SE-set DOF number for the DOF corresponding to GRID_ID_ROW_NUM, COMP


      REAL(DOUBLE)                    :: YSV               ! Enforced displ value read from file LINK1H



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Get YSV values from file L1H array This subr is only called if NDOFSE > 0 so we can assume there is a SE-set

      IERR1H = 0
      DO I = 1,NDOFSE
         READ(L1H,IOSTAT=IOCHK) GRID_ID_ROW_NUM,COMP,YSV
         IF (IOCHK /= 0) THEN
            REC_NO = I
            CALL READERR ( IOCHK, LINK1H, L1H_MSG, REC_NO, OUNT )
            IERR1H = IERR1H + 1
         ELSE
            CALL TDOF_COL_NUM ( 'SE', SE_SET_COL_NUM )
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(GRID_ID_ROW_NUM), IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            TDOF_ROW_NUM = ROW_NUM_START + COMP - 1
            YSDOF = TDOF(TDOF_ROW_NUM,SE_SET_COL_NUM)
            YSe(YSDOF) = YSV
         ENDIF
      ENDDO

! If there were any errors based on reading above file, quit.

      IF (IERR1H > 0) THEN
         CALL FILERR ( OUNT )
         CALL OUTA_HERE ( 'Y' )                                    ! Errors reading YSe file, so quit
      ENDIF

! Rewind L1H and write out ordered YSe array

      REWIND (L1H)
      WRITE(L1H) STIME
      DO I = 1,NDOFSE
         WRITE(L1H) YSe(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE YS_ARRAY

   END MODULE LINK1_LOAD_SUPPORT
