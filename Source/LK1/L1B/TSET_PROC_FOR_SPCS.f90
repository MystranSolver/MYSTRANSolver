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


      SUBROUTINE TSET_PROC_FOR_SPCS ( IERRT )

! DOF Processor for SPC's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1H, L1O, L1O_MSG, LINK1O
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ENFORCED, FATAL_ERR, LSPCADDC, NDOFSB, NDOFSE, NDOFSG, NGRID, NSPCADD,      &
                                         NUM_SPC_RECORDS, NUM_SPC1_RECORDS, NUM_SPCSIDS, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, SPCADD_SIDS, SPCSET, SPCSIDS, SUBLOD

      USE TSET_PROC_FOR_SPCS_USE_IFs                       ! Added 2019/07/14

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC_FOR_SPCS'
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! The name of a DOF set (e.g. 'SB', 'A ', etc)

      INTEGER(LONG), INTENT(INOUT)    :: IERRT             ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: GRID1             ! An actual grid ID
      INTEGER(LONG)                   :: GRID2             ! An actual grid ID
      INTEGER(LONG)                   :: DOF_ERR   = 0     ! Count of errors that result from setting displ sets in TSET
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: I,J,K,GRID_NUM    ! DO loop indices
      INTEGER(LONG)                   :: ICOMP             ! DOF components read from file LINK1O (SPC's) or LINK1N (ASET/OMIT's)
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO    = 0     ! Record number when reading a file
      INTEGER(LONG)                   :: SETID             ! An SPC set ID read from file LINK1O
      LOGICAL                         :: PROCESS_SETID     ! True when this SPC/SPCD record is active for current run


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: RSPC              ! SPC displ value (nonzero's are enforced displ's)



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

      EPS1 = EPSIL(1)

! **********************************************************************************************************************************
! Process the permanent SPC DOF data on the GRID card. The PSPC's are in GRID(I,4) for each GRID PT., I

      GID_ERR = 0
      DOF_ERR = 0
      DOFSET = 'SG'

      DO I=1,NGRID
         IF (GRID(I,4) /= 0) THEN
            CALL RDOF ( GRID(I,4), CDOF )
            CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               IF (CDOF(K) == '1') THEN
                  IF ((TSET(I,K) == '  ') .OR. (TSET(I,K) == 'SB')) THEN
                     TSET(I,K) = DOFSET
                     NDOFSG = NDOFSG + 1
                  ELSE
                     DOF_ERR = DOF_ERR + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1331) GRID(I,1), K, DOFSET, TSET(I,K)
                     WRITE(F06,1331) GRID(I,1), K, DOFSET, TSET(I,K)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IERRT = IERRT + GID_ERR + DOF_ERR

! ----------------------------------------------------------------------------------------------------------------------------------
! Process SPC data from file L1O (data written when SPC and SPC1 Bulk Data cards were read)

! Make a list of all of the individual SPC sets that make up what was called for in Case Control (name this list array SPCSIDS).
! The first one in the list is the SPC set ID called for in Case Control. Subsequent ones are the set ID's from any
! SPCADD cards that have the set ID called for in Case Control. So, if there are 2 SPCADD's that have the same ID as
! called for in Case Control and there are (for example) 25 individual SPC/SPC1 set ID's identified on these
! SPCADD cards, then SPCSIDS will contain 26 entries.

      GID_ERR = 0
      DOF_ERR = 0
      REC_NO  = 0

      NUM_SPCSIDS = 1                                      ! First, count the number, NUM_SPCSIDS, of SPC SID's that will be used
i_do4:DO I=1,NSPCADD
         IF (SPCADD_SIDS(I,1) == SPCSET) THEN
j_do4:      DO J=2,LSPCADDC
               IF (SPCADD_SIDS(I,J) /= 0) THEN
                  NUM_SPCSIDS = NUM_SPCSIDS + 1
                  CYCLE
               ELSE
                  EXIT j_do4
               ENDIF
            ENDDO j_do4
         ENDIF
      ENDDO i_do4

      CALL ALLOCATE_MODEL_STUF ( 'SPCSIDS', SUBR_NAME )    ! Allocate enough memory for array SPCSIDS

      K = 1                                                ! Fill array SPCSIDS, as described above, with SID's that will be used
      SPCSIDS(K) = SPCSET
i_do5:DO I=1,NSPCADD
         IF (SPCADD_SIDS(I,1) == SPCSET) THEN
j_do5:      DO J=2,LSPCADDC
               IF (SPCADD_SIDS(I,J) /= 0) THEN
                  K = K + 1
                  SPCSIDS(K) = SPCADD_SIDS(I,J)
                  CYCLE
               ELSE
                  EXIT j_do5
               ENDIF
            ENDDO j_do5
         ENDIF
      ENDDO i_do5

i_do6:DO I=1,NUM_SPC_RECORDS + NUM_SPC1_RECORDS

         READ(L1O,IOSTAT=IOCHK) SETID, ICOMP, GRID1, GRID2, RSPC, DOFSET
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1O, L1O_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Error reading SPC file . No sense continuing
         ENDIF

         IF ((DOFSET == 'SB') .OR. (DOFSET == 'SE')) THEN  ! Make sure that DOFSET = 'SB' or 'SE'
            CONTINUE
         ELSE
            WRITE(ERR,1369) SUBR_NAME, LINK1O, DOFSET
            WRITE(F06,1369) SUBR_NAME, LINK1O, DOFSET
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Pgm error (data in file LINK1O must be for sets 'SB' or 'SE')
         ENDIF

! No error, so processes data. First, make sure SETID is the one called for in Case Control:

         PROCESS_SETID = .FALSE.
         DO J=1,NUM_SPCSIDS
            IF (SETID == SPCSIDS(J)) THEN
               PROCESS_SETID = .TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (.NOT. PROCESS_SETID) THEN
            DO J=1,NSUB
               IF (SETID == SUBLOD(J,1)) THEN
                  PROCESS_SETID = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         IF (PROCESS_SETID) THEN
                                                        ! Check for existence of grid pt
            CALL RDOF ( ICOMP, CDOF )                   ! Convert ICOMP to CDOF char form for use below

            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID1, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', GRID1, 'SPC OR SPC1/SPCD', SETID
               WRITE(F06,1822) 'GRID ', GRID1, 'SPC OR SPC1/SPCD', SETID
            ENDIF
            IF (GRID2 /= GRID1) THEN
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID2, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', GRID2, 'SPC OR SPC1/SPCD', SETID
                  WRITE(F06,1822) 'GRID ', GRID2, 'SPC OR SPC1/SPCD', SETID
               ENDIF
            ENDIF

            IF (GID_ERR == 0) THEN                      ! Put CDOF data in TSET for GRID1 thru GRID2
               DO GRID_NUM=GRID1,GRID2                  ! GRID2 >= GRID1 was checked in subr BD_SPC1
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, GRID_ID_ROW_NUM )
                  IF (GRID_ID_ROW_NUM /= -1) THEN
                     CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM, NUM_COMPS, SUBR_NAME )
                     DO K = 1,NUM_COMPS                 ! Put data in TSET and write enforced displ to L1H.
                        IF (CDOF(K) == '1') THEN
                           IF      (DOFSET == 'SE') THEN
! --- cbeam_005f_enforced_spc begin --- !
                              IF ((TSET(GRID_ID_ROW_NUM,K) == '  ') .OR. (TSET(GRID_ID_ROW_NUM,K) == 'SB')) THEN
                                 IF (TSET(GRID_ID_ROW_NUM,K) == 'SB') THEN
                                    NDOFSB = NDOFSB - 1
                                 ENDIF
                                 TSET(GRID_ID_ROW_NUM,K) = 'SE'
                                 NDOFSE = NDOFSE + 1
                                                           ! Write enforced displs to L1H. If ENFORCED = 'Y' write all terms
                                 IF (ENFORCED == 'Y') THEN
                                    IF (DOFSET == 'SE') THEN
                                       WRITE(L1H) GRID_ID_ROW_NUM, K, RSPC
                                    ENDIF
                                 ELSE                   ! If ENFORCED = 'N' write only nonzero terms
                                    IF ((DOFSET == 'SE') .AND. (ABS(RSPC) > EPS1)) THEN
                                       WRITE(L1H) GRID_ID_ROW_NUM, K, RSPC
                                    ENDIF
                                 ENDIF
                              ELSE
                                 DOF_ERR = DOF_ERR + 1
                                 FATAL_ERR = FATAL_ERR + 1
                                 WRITE(ERR,1332) SETID, GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                                 WRITE(F06,1332) SETID, GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                              ENDIF
! --- cbeam_005f_enforced_spc end --- !
                           ELSE IF (DOFSET == 'SB') THEN
                              IF ((TSET(GRID_ID_ROW_NUM,K) == '  ') .OR. (TSET(GRID_ID_ROW_NUM,K) == 'SG') .OR.                 &
                                  (TSET(GRID_ID_ROW_NUM,K) == 'SB')) THEN
                                 TSET(GRID_ID_ROW_NUM,K) = 'SB'
                                 NDOFSB = NDOFSB + 1
                              ELSE IF (TSET(GRID_ID_ROW_NUM,K) == 'SE') THEN
                                 CONTINUE
!xx-DOFSET CAN'T BE 'SE' HERE       IF ((DOFSET == 'SE') .AND. (DABS(RSPC) > EPS1)) THEN
!xx-DUE TO TEST 4 LINES ABOVE          WRITE(L1H) GRID_ID_ROW_NUM, K, RSPC
!xx-MOD ON 12/28/06                 ENDIF
                              ELSE
                                 DOF_ERR = DOF_ERR + 1
                                 FATAL_ERR = FATAL_ERR + 1
                                 WRITE(ERR,1332) SETID, GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                                 WRITE(F06,1332) SETID, GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

      ENDDO i_do6

      CALL DEALLOCATE_MODEL_STUF ( 'SPCSIDS' )

      IERRT = IERRT + GID_ERR + DOF_ERR



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1331 FORMAT(' *ERROR  1331: GRID POINT ',I8,' HAS COMPONENT ',I2,' IN THE ',A2,' DISPL SET. (PERM SPC ON GRID ENTRY)',            &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')

 1332 FORMAT(' *ERROR  1332: SPC SET ',I8,' HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',                      &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')

 1369 FORMAT(' *ERROR  1369: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RECORD READ FROM FILE: ',A                                                                            &
                    ,/,14X,' INDICATES THAT AN SPCd DOF BELONGS TO THE "',A2,'" SET. MUST BE EITHER "SE" OR "SB"')

! **********************************************************************************************************************************

      END SUBROUTINE TSET_PROC_FOR_SPCS
