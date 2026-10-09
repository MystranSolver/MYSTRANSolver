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

   MODULE DOF_SET_CONSTRUCTION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TSET_PROC, RDOF

   CONTAINS

      SUBROUTINE TSET_PROC

! Generate TSET table. The displ set definitions are:

!   G---------->N + M (M = DOF's identified as dependent on MPC's or rigid elements)
!               |
!               |
!               |
!               |
!               |---------->F + S (S is SZ (SPC'd to 0 displ) or SE (enforced displ))
!                           |           SZ can be either: SB (defined on SPC/SPC1 cards) or
!                           |                             SG (defined on GRID cards) or
!                           |                             SA (based on AUTOSPC)
!                           |
!                           |---------->A + O (O = DOF's omitted via Guyan reduction)
!                                       |
!                                       |
!                                       |
!                                       |---------->L + R (R = DOF on SUPORT B.D. entries))

! Every DOF is in either: M (MPC'd via MPC eqns or rigid elements), or
!                         S (SPC'd),or
!                         O (Omitted via Guyan Reduction),or
!                         R (SUPORT'd), or
!                         L (Left over)

! TSET is a table that defines which displ set each grid/component pair belongs to. It has: R, A, O, (SG, SB, SE) or M for each DOF.
! The table has NGRID rows and 6 columns (1 col for each of the 6 components of displ at a grid). An example of a TSET table written
! to the F06 file is shown below:

!                                      DEGREE OF FREEDOM SET TABLE (TSET)

!             GRID   GRID_SEQ   INV_GRID_SEQ*      T1       T2       T3       R1       R2       R3

!             101           4              7       L        L        L        SG       SG       SG
!             102           5              5       M        M        M        L        L        L
!             103           3              3       L        L        L        L        L        L
!             104           6              1       L        L        L        L        L        L
!             105           2              2       L        L        L        L        L        L
!             106           7              4       L        L        L        O        O        O
!             107           1              6       SE       SG       SG       SG       SG       SG

!              * INV_GRID_SEQ = J meams that the J-th entry under GRID is sequenced GRID_SEQ(J)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NGRID, NAOCARD, NUM_SUPT_CARDS,                                  &
                                         NDOFL, NDOFM, NDOFO, NDOFR, NDOFS, NDOFSA, NDOFSG, NDOFSB, NDOFSE, NDOFSZ
      USE PARAMS, ONLY                :  PRTTSET
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET
      USE MODEL_STUF, ONLY            :  GRID

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_TSET

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC'

      INTEGER(LONG)                   :: I,K               ! DO loop indices
      INTEGER(LONG)                   :: IERRT     = 0     ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)




! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

! ----------------------------------------------------------------------------------------------------------------------------------
! TSET for cols 2-6 are not defined for SPOINT's, so set = '--'

      WRITE(SC1,12345,ADVANCE='NO') '       Initializing           ', CR13
      DO I=1,NGRID
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS == 1) THEN
            DO K=2,6
               TSET(I,K) = '--'
            ENDDO
         ENDIF
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Process M-set (MPC's and rigid elements)

      NDOFM = 0
      WRITE(SC1,12345,ADVANCE='NO') '       Process MPCs           ', CR13
      CALL TSET_PROC_FOR_MPCS ( IERRT )                    ! Process MPC data from file L1S

      WRITE(SC1,12345,ADVANCE='NO') '       Process Rigid Elements', CR13
      CALL TSET_PROC_FOR_RIGELS ( IERRT )                  ! Process the Rigid Element data from file L1F

! ----------------------------------------------------------------------------------------------------------------------------------
! Process S-set (SPC's)

      NDOFSG  = 0
      NDOFSB  = 0
      NDOFSE  = 0
      WRITE(SC1,12345,ADVANCE='NO') '       Process SPCs           ', CR13
      CALL TSET_PROC_FOR_SPCS ( IERRT )                    ! Process SPC data from file L1O
      NDOFSZ = NDOFSA + NDOFSB + NDOFSG
      NDOFS  = NDOFSZ + NDOFSE

! ----------------------------------------------------------------------------------------------------------------------------------
! Process O-set (OMIT's/ASET's)

      NDOFO = 0
      IF (NAOCARD == 0) THEN                               ! If no ASET,1/OMIT,1 cards, then, for time being, set all remaining DOF
         DO I = 1,NGRID                                    ! to A-set  (if there are SUPORT's some will get changed to R set below)
            CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               IF (TSET(I,K) == '  ') THEN
                  TSET(I,K) = 'A '
               ENDIF
            ENDDO
         ENDDO
      ELSE
         WRITE(SC1,12345,ADVANCE='NO') '       Process OMITs          ', CR13
         CALL TSET_PROC_FOR_OMITS ( IERRT )
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Check that all grid/components have been set

      DO I=1,NGRID
         DO K=1,6
            IF (TSET(I,K) == '  ') THEN
               WRITE(ERR,1368) SUBR_NAME,GRID(I,1),K
               WRITE(F06,1368) SUBR_NAME,GRID(I,1),K
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ENDDO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Process R-set (SUPORT's)

      NDOFR = 0
      IF (NUM_SUPT_CARDS > 0) THEN
         WRITE(SC1,12345,ADVANCE='NO') '       Process SUPORTs        ', CR13
         CALL TSET_PROC_FOR_SUPORTS ( IERRT )
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! All remaining A set are also L set

      NDOFL = 0
      DO I = 1,NGRID
         DO K = 1,6
            IF (TSET(I,K) == 'A ') THEN
               TSET(I,K) = 'L '
               NDOFL = NDOFL + 1
            ENDIF
         ENDDO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Print the TSET table, if requested

      IF (PRTTSET == 1) THEN
         CALL WRITE_TSET
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
      WRITE(SC1,*) CR13

! Quit if there were errors

      IF (IERRT > 0) THEN
         WRITE(ERR,9996) SUBR_NAME,IERRT
         WRITE(F06,9996) SUBR_NAME,IERRT
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! *********************************************************************************************************************************
 1368 FORMAT(' *ERROR  1368: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' GRID POINT ',I8,' COMPONENT ',I2,' IS NOT DEFINED IN ANY DISPL SET')

 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

12345 FORMAT(A, A)

! *********************************************************************************************************************************

      END SUBROUTINE TSET_PROC


      SUBROUTINE TSET_PROC_FOR_MPCS ( IERRT )

! DOF Processor for MPC's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1S, L1S_MSG, LINK1S
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LIND_GRDS_MPCS, LMPCADDC, NDOFM, NGRID, NIND_GRDS_MPCS, NMPC,    &
                                         NMPCADD, NTERM_RMG, NUM_MPCSIDS
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID_ID, MPC_IND_GRIDS, MPCSET, MPCSIDS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_ARRAY_ROW_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC_FOR_MPCS'
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! The name of a DOF set (e.g. 'SB', 'A ', etc)
      CHARACTER( 1*BYTE)              :: MPC_SET_USED      ! 'Y'/'N' indicator if an MPC set in B.D. is used

      INTEGER(LONG), INTENT(INOUT)    :: IERRT             ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: AGRID_D           ! Dep grid for an MPC read from file LINK1S
      INTEGER(LONG)                   :: AGRID_I           ! An independent grid in an MPC equation

      INTEGER(LONG)                   :: AGRID_I_PREV = 0  ! Previous value of AGRID_I read from file L1F. (each RBE2 element has a
!                                                            separate record in file L1F for each dep  grid, so, if there are 10
!                                                            dep grids on one RBE2 then there will be 10 records.The 10 AGRID_I's
!                                                            read from these 10 records are all identical and we only want to write
!                                                            the AGRID_I once to array MPC_IND_GRIDS.

      INTEGER(LONG)                   :: COMPS_D           ! Comp number for the dep grid for an MPC eqn (read from file LINK1S)
      INTEGER(LONG)                   :: COMPS_D_TSET      ! COMPS_D converted (e.g. if COMPS_D read from L1S is 0, change to 1)
      INTEGER(LONG)                   :: DOF_ERR   = 0     ! Count of errors that result from setting displ sets in TSET
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: IJUNK             ! A grid read from file LINK1S that we do not need
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: NUM_TRIPLES       ! Counter on number of pairs of grid/comp/coeff triplets on this MPC
!                                                            card. Must be <= MMPC which was counted in subr BD_MPC0
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO    = 0     ! Record number when reading a file
      INTEGER(LONG)                   :: SETID             ! An SPC set ID read from file LINK1O


      REAL(DOUBLE)                    :: RJUNK             ! An MPC coeff value read from file LINK1S that we do not need



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! Process MPC data from file L1S (data written when MPC Bulk Data cards were read)

! File LINK1S contains data from the NMPC number of logical MPC cards in the input B.D. deck. For each logical MPC card, LINK1S has:
!     1st record          for 1st MPC: the MPC set ID
!     2nd record          for 1st MPC: the num of triplets of grid/comp/coeff (incl ones for the dependent DOF) on this logical MPC
!     3rd record          for 1st MPC: grid/comp/coeff for the dependent DOF on the this MPC logical card
!     4th record          for 1st MPC: grid/comp/coeff for the 1st independent DOF on the this MPC logical card
!     5th record, and on, for 1st MPC: grid/comp/coeff for the 2nd, and on, independent DOF's (if any) on the this MPC logical card

! The above record structure is repeated for each MPC logical card in the data deck (in the order in which they were read from the
! B.D. deck). All logical MPC cards are included, not only the ones that may be used in a particular execution of MYSTRAN

      GID_ERR = 0
      DOF_ERR = 0

      REC_NO = 0
i_do3:DO I=1,NMPC                                          ! Process data from file LINK1S (contains all info from the NMPC MPC's)

         READ(L1S,IOSTAT=IOCHK) SETID                      ! Read the SETID for the i-th logical MPC
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         READ(L1S,IOSTAT=IOCHK) NUM_TRIPLES                ! Read the number of triplets of grid/comp/coeff for the i-th logical MPC
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )
         ENDIF


         MPC_SET_USED = 'N'

j_do3:   DO J=1,NUM_MPCSIDS                                ! NUM_MPCSIDS will be 1 if there are any MPC's. It will be > 1 if there
!                                                            are MPCADD entries that match the MPC set requested in Case Control
            IF (SETID == MPCSIDS(J)) THEN

               MPC_SET_USED = 'Y'
               NTERM_RMG = NTERM_RMG + NUM_TRIPLES
                                                           ! Read dependent grid/comp/coeff from for the i-th logical MPC
               READ(L1S,IOSTAT=IOCHK) AGRID_D, COMPS_D, RJUNK

               IF (COMPS_D == 0) THEN                      ! If SPOINT change COMPS_D from 0 to 1
                  COMPS_D_TSET = 1
               ELSE
                  COMPS_D_TSET = COMPS_D
               ENDIF

               REC_NO = REC_NO + 1
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
                                                           ! Get the row number, in array GRID_ID, for dependent grid, AGRID_D
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', AGRID_D, 'MPC ', MPCSIDS(J)
                  WRITE(F06,1822) 'GRID ', AGRID_D, 'MPC ', MPCSIDS(J)
               ENDIF

               DO K=1,NUM_TRIPLES-1                        ! Read the indep grid/comp/coeff values (only need indep grid ID's)
                  READ(L1S,IOSTAT=IOCHK) AGRID_I, IJUNK, RJUNK
                  REC_NO = REC_NO + 1
                  IF (IOCHK /= 0) THEN
                     CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  IF (AGRID_I /= AGRID_I_PREV) THEN
                     NIND_GRDS_MPCS = NIND_GRDS_MPCS + 1      ! Up count on num indep grids & store grid in array MPC_IND_GRIDS
                     IF (NIND_GRDS_MPCS > LIND_GRDS_MPCS) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, LIND_GRDS_MPCS, 'MPC_IND_GRIDS')
                     MPC_IND_GRIDS(NIND_GRDS_MPCS) = AGRID_I
                     AGRID_I_PREV = AGRID_I
                  ENDIF
               ENDDO

               IF (GID_ERR > 0) THEN
                  CYCLE j_do3
               ENDIF

               DOFSET = 'M '
               IF(TSET(GRID_ID_ROW_NUM,COMPS_D_TSET) == '  ') THEN
                  TSET(GRID_ID_ROW_NUM,COMPS_D_TSET) = DOFSET
                  NDOFM = NDOFM + 1
               ELSE
                  DOF_ERR = DOF_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1333) SETID, AGRID_D, COMPS_D_TSET, DOFSET, TSET(GRID_ID_ROW_NUM, COMPS_D_TSET)
                  WRITE(F06,1333) SETID, AGRID_D, COMPS_D_TSET, DOFSET, TSET(GRID_ID_ROW_NUM, COMPS_D_TSET)
               ENDIF

               EXIT j_do3                                  ! We found and processed an MPC set, so exit this loop

            ELSE

               MPC_SET_USED = 'N'
               CYCLE j_do3

            ENDIF

         ENDDO j_do3

         IF (MPC_SET_USED == 'N') THEN                     ! This MPC set is not to be used, so skip all grid/comp/coeff records
            DO K=1,NUM_TRIPLES
               READ(L1S,IOSTAT=IOCHK) IJUNK, IJUNK, RJUNK
               REC_NO = REC_NO + 1
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDDO
         ENDIF

      ENDDO i_do3

      IERRT = IERRT + GID_ERR + DOF_ERR



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1333 FORMAT(' *ERROR  1333: MPC SET ',I8,' HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',                      &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')

! **********************************************************************************************************************************

      END SUBROUTINE TSET_PROC_FOR_MPCS


      SUBROUTINE TSET_PROC_FOR_OMITS ( IERRT )

! DOF Processor for MPC's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1N, L1N_MSG, LINK1N
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NAOCARD, NDOFO, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC_FOR_OMITS'
      CHARACTER( 1*BYTE)              :: ASET_FND          ! 'Y' if there are ASET/ASET1 data in file LINK1N
      CHARACTER( 1*BYTE)              :: CDOF1(6)          ! An output from subr RDOF
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! The name of a DOF set (e.g. 'SB', 'A ', etc)
      CHARACTER( 1*BYTE)              :: OMIT_FND          ! 'Y' if there are OMIT/OMIT1 data in file LINK1N

      INTEGER(LONG), INTENT(INOUT)    :: IERRT             ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: DOF_ERR   = 0     ! Count of errors that result from setting displ sets in TSET
      INTEGER(LONG)                   :: GID1              ! 1st grid on RBAR element read from file LINK1F
      INTEGER(LONG)                   :: GID2              ! 2nd grid on RBAR element read from file LINK1F
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: ICOMP             ! DOF components read from file LINK1O (SPC's) or LINK1N (ASET/OMIT's)
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO    = 0     ! Record number when reading a file




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! Process the ASET / OMIT data in file L1N
! If there is only ASET/ASET1 data in file L1N, then all remaining DOF will be O-set
! If there is only OMIT/OMIT1 data in file L1N, then all remaining DOF will be A-set
! If there is both ASET/ASET1 and OMIT/OMIT1 then they must define all remaining DOF's, or error

      GID_ERR = 0
      DOF_ERR = 0
      REC_NO  = 0

      IF (NAOCARD > 0) THEN

         ASET_FND = 'N'
         OMIT_FND = 'N'

         DO I=1,NAOCARD

            READ(L1N,IOSTAT=IOCHK) ICOMP,GID1,GID2,DOFSET  ! Read a record from L1N
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1N, L1N_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading ASET/OMIT file. No sense continuing
            ENDIF

            IF      (DOFSET == 'A ') THEN                  ! Make sure DOFSET is "A " or "O "
               ASET_FND = 'Y'
            ELSE IF (DOFSET == 'O ') THEN
               OMIT_FND = 'Y'
            ELSE
               WRITE(ERR,1363) SUBR_NAME,LINK1N,DOFSET
               WRITE(F06,1363) SUBR_NAME,LINK1N,DOFSET
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Pgm error (DOFSET is not A-set or O-set), so quit
            ENDIF

            CALL RDOF ( ICOMP, CDOF1 )                     ! Convert ICOMP to CDOF1 char forn for use below
                                                           ! Make sure grid GID1 & GID2 exist
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID1, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1365) GID1
               WRITE(F06,1365) GID1
            ENDIF
            IF (GID2 /= GID1) THEN
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID2, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1365) GID2
                  WRITE(F06,1365) GID2
               ENDIF
            ENDIF

            IF ((GID_ERR == 0) .AND. (GID2 >= GID1)) THEN  ! Put CDOF1 data in TSET for GID1 thru GID2
               DO J=GID1,GID2                              ! GID2 > GID1 was checked when ASET/OMIT B.D. cards were read
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, J, GRID_ID_ROW_NUM )
                  IF (GRID_ID_ROW_NUM /= -1) THEN
                     CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM, NUM_COMPS, SUBR_NAME )
                     DO K = 1,NUM_COMPS
                        IF (CDOF1(K) == '1') THEN
                           IF ((TSET(GRID_ID_ROW_NUM,K) == '  ') .OR. (TSET(GRID_ID_ROW_NUM,K) == DOFSET)) THEN
                              TSET(GRID_ID_ROW_NUM,K) = DOFSET
                              IF (DOFSET == 'O ') THEN
                                 NDOFO = NDOFO + 1
                              ENDIF
                           ELSE
                              DOF_ERR = DOF_ERR + 1
                              FATAL_ERR = FATAL_ERR + 1
                              WRITE(ERR,1366) J,K,DOFSET,TSET(GRID_ID_ROW_NUM,K)
                              WRITE(F06,1366) J,K,DOFSET,TSET(GRID_ID_ROW_NUM,K)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

         ENDDO                                             ! DO loop on I = 1 to NAOCARD

         IF ((ASET_FND == 'Y').AND.(OMIT_FND == 'N')) THEN ! Make all DOF's not as yet set be O-set
            DO I=1,NGRID
               CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
               DO K=1,NUM_COMPS
                  IF (TSET(I,K) == '  ') THEN
                     TSET(I,K) = 'O '
                     NDOFO = NDOFO + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         IF ((ASET_FND == 'N').AND.(OMIT_FND == 'Y')) THEN ! Make all DOF's not as yet set be A-set
            DO I=1,NGRID
               CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
               DO K=1,NUM_COMPS
                  IF (TSET(I,K) == '  ') THEN
                     TSET(I,K) = 'A '
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

      ENDIF                                                ! IF (NAOCARD > 0)

      IERRT = IERRT + GID_ERR + DOF_ERR



      RETURN

! **********************************************************************************************************************************
 1363 FORMAT(' *ERROR  1363: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RECORD READ FROM FILE: '                                                                              &
                    ,/,14X,' HAS INCORRECT SET DEFINITION = ',A2,'. SHOULD BE EITHER "A " OR "O "')

 1365 FORMAT(' *ERROR  1365: GRID ',I8,' SPECIFIED ON ASET/ASET1 OR OMIT/OMIT1 ENTRY DOES NOT EXIST')

 1366 FORMAT(' *ERROR  1366: GRID POINT ',I8,' HAS COMPONENT ',I2,' IN THE ',A2,' DISPL SET',                                      &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')


! **********************************************************************************************************************************

      END SUBROUTINE TSET_PROC_FOR_OMITS


      SUBROUTINE TSET_PROC_FOR_RIGELS ( IERRT )

! DOF Processor for rigid elements (incl RBE3 and RSPLINE)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1F, L1F_MSG, LINK1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LIND_GRDS_MPCS, NDOFM, NGRID, NIND_GRDS_MPCS, NRECARD
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, MPC_IND_GRIDS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_ARRAY_ROW_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC_FOR_RIGELS'
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF
      CHARACTER( 1*BYTE)              :: CDOF1(6)          ! An output from subr RDOF
      CHARACTER( 1*BYTE)              :: CDOF2(6)          ! An output from subr RDOF
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! The name of a DOF set (e.g. 'SB', 'A ', etc)
      CHARACTER( 8*BYTE)              :: RTYPE             ! Rigid elem type (RBAR, RBE1, RBE2)

      INTEGER(LONG), INTENT(INOUT)    :: IERRT             ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: COMPS_D           ! Dep DOF's on RBE1 or RBE2 elem read from file LINK1F
      INTEGER(LONG)                   :: DDOF1             ! Dep DOF's at GID1 on RBAR elem read from file LINK1F
      INTEGER(LONG)                   :: DDOF2             ! Dep DOF's at GID2 on RBAR elem read from file LINK1F
      INTEGER(LONG)                   :: AGRID_D           ! Dep grid for RBE1 or RBE2 elem read from file LINK1F/LINK1S
      INTEGER(LONG)                   :: DOF_ERR    = 0    ! Count of errors that result from setting displ sets in TSET
      INTEGER(LONG)                   :: GID1              ! A grid number
      INTEGER(LONG)                   :: GID2              ! A grid number
      INTEGER(LONG)                   :: GID_ERR    = 0    ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_I ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IDUM(6)           ! Integer values read that are not used
      INTEGER(LONG)                   :: AGRIDI_I(6)       ! Up to 6 indep grids for RBE1 elem read from file LINK1F
      INTEGER(LONG)                   :: AGRID_I           ! Indep grid on RBE elem read from file LINK1F
      INTEGER(LONG)                   :: AGRID_I1          ! Indep grid on RBE elem read from file LINK1F
      INTEGER(LONG)                   :: AGRID_I2          ! Indep grid on RBE elem read from file LINK1F
      INTEGER(LONG)                   :: AGRID_I_PREV      ! Previous value of AGRID_I read from file L1F. (each RBE2 element has a
!                                                            separate record in file L1F for each dep  grid, so, if there are 10
!                                                            dep grids on one RBE2 then there will be 10 records.The 10 IGID0's
!                                                            read from these 10 records are all identical and we only want to write
!                                                            the IGID0 once to array MPC_IND_GRIDS.
      INTEGER(LONG)                   :: AGRID_I1_PREV     ! Same general description as for AGRID_I_PREV
      INTEGER(LONG)                   :: AGRID_I2_PREV     ! Same general description as for AGRID_I_PREV
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: IRBE3             ! Number of triplets of grid/comp/weight on L1F for RBE# elems
      INTEGER(LONG)                   :: NUMI              ! Noumber of pairs of grid/indep DOF's read from file LINK1F for RBE1's
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO     = 0    ! Record number when reading a file
      INTEGER(LONG)                   :: REFC              ! Dependent components on RBE3
      INTEGER(LONG)                   :: REFGRID           ! Dependent grid on RBE3
      INTEGER(LONG)                   :: REID              ! Rigid elem ID read from file LINK1F


      REAL(DOUBLE)                    :: RDUM              ! Real value read that is not used



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! Process the Rigid Element data from file L1F (data written when RBAR, RBE1, RBE2 Bulk Data cards were read)

      AGRID_I_PREV  = 0
      AGRID_I1_PREV = 0
      AGRID_I2_PREV = 0

      DO I=1,NRECARD

         READ(L1F,IOSTAT=IOCHK) RTYPE
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Error reading rigid elem file. No sense continuing
         ENDIF

! RBAR rigid element

         IF      (RTYPE == 'RBAR    ') THEN

            READ(L1F,IOSTAT=IOCHK) REID,GID1,IDUM(1),DDOF1,GID2,IDUM(2),DDOF2
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading rigid elem file. No sense continuing
            ENDIF
                                                           ! Check for existence of grid pt. GID1
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID1, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', GID1, RTYPE, REID
               WRITE(F06,1822) 'GRID ', GID1, RTYPE, REID
            ENDIF
                                                           ! Check for existence of grid pt. GID2
            IF (DDOF1 > 0) THEN                            ! Convert DDOF1 to CDOF1 char forn for use below
               CALL RDOF ( DDOF1, CDOF1 )
            ELSE
               DO K=1,6
                  CDOF1(K) = '0'
               ENDDO
            ENDIF

            DOFSET = 'M '
            IF (GRID_ID_ROW_NUM > 0) THEN                  ! Put CDOF1 data in TSET
               DO K = 1,6
                  IF (CDOF1(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, GID1, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, GID1, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, GID1, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                           WRITE(F06,1330) RTYPE, REID, GID1, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID2, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', GID2, RTYPE, REID
               WRITE(F06,1822) 'GRID ', GID2, RTYPE, REID
            ENDIF

            IF (DDOF2 > 0) THEN                            ! Convert DDOF2 to CDOF2 char forn for use below
               CALL RDOF ( DDOF2, CDOF2 )
            ELSE
               DO K=1,6
                  CDOF2(K) = '0'
               ENDDO
            ENDIF

            IF (GRID_ID_ROW_NUM > 0) THEN                  ! Put CDOF2 data in TSET
               DO K = 1,6
                  IF (CDOF2(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, GID2, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, GID2, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, GID2, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                           WRITE(F06,1330) RTYPE, REID, GID2, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

! RBE1 rigid element

         ELSE IF (RTYPE == 'RBE1    ') THEN
                                                           ! NUMI was checked to be <= 6 when RBE1 was read in Bulk Data
            READ(L1F,IOSTAT=IOCHK) REID,AGRID_D,COMPS_D,NUMI,(AGRIDI_I(J),IDUM(J),J=1,NUMI)
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading rigid elem file. No sense continuing
            ENDIF

            DO J=1,NUMI                                    ! Check for existence of independent grids, AGRIDI_I(1-6)
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRIDI_I(J), GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', AGRIDI_I(J), RTYPE, REID
                  WRITE(F06,1822) 'GRID ', AGRIDI_I(J), RTYPE, REID
               ENDIF
            ENDDO
                                                            ! Check for existence of dependent grid, AGRID_D
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_D, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_D, RTYPE, REID
            ENDIF

            CALL RDOF ( COMPS_D, CDOF )                    ! Convert COMPS_D to CDOF char forn for use below

            DOFSET = 'M '                                  ! Put data in TSET
            IF (GRID_ID_ROW_NUM > 0) THEN
               DO K = 1,6
                  IF (CDOF(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, AGRID_D, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, AGRID_D, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                           WRITE(F06,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

! RBE2 rigid element

         ELSE IF (RTYPE == 'RBE2    ') THEN
            READ(L1F,IOSTAT=IOCHK) REID, AGRID_D, COMPS_D, AGRID_I
            IF (AGRID_I /= AGRID_I_PREV) THEN
               NIND_GRDS_MPCS = NIND_GRDS_MPCS + 1         ! Up count on num indep grids & store grid in array MPC_IND_GRIDS
               IF (NIND_GRDS_MPCS > LIND_GRDS_MPCS) CALL ARRAY_SIZE_ERROR_1  ( SUBR_NAME, LIND_GRDS_MPCS, 'MPC_IND_GRIDS' )
               MPC_IND_GRIDS(NIND_GRDS_MPCS) = AGRID_I
               AGRID_I_PREV = AGRID_I
            ENDIF
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading rigid elem file. No sense continuing
            ENDIF
                                                           ! Check for existence of dependent grid, AGRID_D
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )
            IF (GRID_ID_ROW_NUM_D == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_D, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_D, RTYPE, REID
            ENDIF
                                                           ! Check for existence of independent grid, AGRID_I
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I, GRID_ID_ROW_NUM_I )
            IF (GRID_ID_ROW_NUM_I == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_I, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_I, RTYPE, REID
            ENDIF

            CALL RDOF ( COMPS_D, CDOF )                    ! Convert COMPS_D to CDOF char forn for use below

            DOFSET = 'M '                                  ! Put data in TSET
            IF (GRID_ID_ROW_NUM_D > 0) THEN
               DO K = 1,6
                  IF (CDOF(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM_D,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM_D,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM_D,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, AGRID_D, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, AGRID_D, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                           WRITE(F06,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

! RBE3 rigid element

         ELSE IF (RTYPE == 'RBE3    ') THEN
            READ(L1F,IOSTAT=IOCHK) REID, REFGRID, REFC, IRBE3, RDUM
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading rigid elem file. No sense continuing
            ENDIF
                                                           ! Check for existence of dependent grid, REFGRID
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, REFGRID, GRID_ID_ROW_NUM_D )
            IF (GRID_ID_ROW_NUM_D == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', REFGRID, RTYPE, REID
               WRITE(F06,1822) 'GRID ', REFGRID, RTYPE, REID
            ENDIF
                                                           ! Check for existence of independent grids
            DO J=1,IRBE3
               READ(L1F) AGRID_I, IDUM(1), RDUM
               REC_NO = REC_NO + 1
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
                  CALL OUTA_HERE ( 'Y' )                   ! Error reading rigid elem file. No sense continuing
               ENDIF
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I, GRID_ID_ROW_NUM_I )
               IF (GRID_ID_ROW_NUM_I == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', AGRID_I, RTYPE, REID
                  WRITE(F06,1822) 'GRID ', AGRID_I, RTYPE, REID
               ENDIF
            ENDDO

            CALL RDOF ( REFC, CDOF )                       ! Convert REFC to CDOF char forn for use below

            DOFSET = 'M '                                  ! Put data in TSET
            IF (GRID_ID_ROW_NUM_D > 0) THEN
               DO K = 1,6
                  IF (CDOF(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM_D,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM_D,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM_D,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, REFGRID, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, REFGRID, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, REFGRID, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                           WRITE(F06,1330) RTYPE, REID, REFGRID, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

! RSPLINE rigid element

         ELSE IF (RTYPE == 'RSPLINE ') THEN

            READ(L1F,IOSTAT=IOCHK) REID, IDUM(1), IDUM(2), AGRID_I1, AGRID_I2, AGRID_D, COMPS_D, RDUM

            IF (AGRID_I1 /= AGRID_I1_PREV) THEN
               NIND_GRDS_MPCS = NIND_GRDS_MPCS + 1         ! Up count on num indep grids & store grid in array MPC_IND_GRIDS
               IF (NIND_GRDS_MPCS > LIND_GRDS_MPCS) CALL ARRAY_SIZE_ERROR_1  ( SUBR_NAME, LIND_GRDS_MPCS, 'MPC_IND_GRIDS' )
               MPC_IND_GRIDS(NIND_GRDS_MPCS) = AGRID_I1
               AGRID_I1_PREV = AGRID_I1
            ENDIF

            IF (AGRID_I2 /= AGRID_I2_PREV) THEN
               NIND_GRDS_MPCS = NIND_GRDS_MPCS + 1         ! Up count on num indep grids & store grid in array MPC_IND_GRIDS
               IF (NIND_GRDS_MPCS > LIND_GRDS_MPCS) CALL ARRAY_SIZE_ERROR_1  ( SUBR_NAME, LIND_GRDS_MPCS, 'MPC_IND_GRIDS' )
               MPC_IND_GRIDS(NIND_GRDS_MPCS) = AGRID_I2
               AGRID_I2_PREV = AGRID_I2
            ENDIF

            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading rigid elem file. No sense continuing
            ENDIF
                                                           ! Check for existence of dependent grid, AGRID_D
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )
            IF (GRID_ID_ROW_NUM_D == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_D, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_D, RTYPE, REID
            ENDIF
                                                           ! Check for existence of independent grid, AGRID_I1
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I1, GRID_ID_ROW_NUM_I )
            IF (GRID_ID_ROW_NUM_I == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_I1, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_I1, RTYPE, REID
            ENDIF
                                                           ! Check for existence of independent grid, AGRID_I2
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I2, GRID_ID_ROW_NUM_I )
            IF (GRID_ID_ROW_NUM_I == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', AGRID_I2, RTYPE, REID
               WRITE(F06,1822) 'GRID ', AGRID_I2, RTYPE, REID
            ENDIF

            CALL RDOF ( COMPS_D, CDOF )                    ! Convert COMPS_D to CDOF char forn for use below

            DOFSET = 'M '                                  ! Put data in TSET
            IF (GRID_ID_ROW_NUM_D > 0) THEN
               DO K = 1,6
                  IF (CDOF(K) == '1') THEN
                     IF (TSET(GRID_ID_ROW_NUM_D,K) == '  ') THEN
                        TSET(GRID_ID_ROW_NUM_D,K) = DOFSET
                        NDOFM = NDOFM + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        IF (TSET(GRID_ID_ROW_NUM_D,K) == '--') THEN
                           WRITE(ERR,1313) RTYPE, REID, AGRID_D, K, DOFSET
                           WRITE(F06,1313) RTYPE, REID, AGRID_D, K, DOFSET
                        ELSE
                           WRITE(ERR,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                           WRITE(F06,1330) RTYPE, REID, AGRID_D, K, DOFSET, TSET(GRID_ID_ROW_NUM_D,K)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

         ELSE

            WRITE(ERR,1329) SUBR_NAME,RTYPE
            WRITE(F06,1329) SUBR_NAME,RTYPE
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ENDDO

      IERRT = IERRT + GID_ERR + DOF_ERR



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1329 FORMAT(' *ERROR  1329: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RIGID ELEMENT TYPE MUST BE "RBAR", "RBE1", OR "RBE2" BUT IS "',A8,'"')

 1313 FORMAT(' *ERROR  1330: ',A8,' RIGID ELEMENT NUMBER ',I8,' HAS SPOINT POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',&
                           ' HOWEVER THIS COMPONENT IS UNDEFINED FOR SPOINT''s')

 1330 FORMAT(' *ERROR  1330: ',A8,' RIGID ELEMENT NUMBER ',I8,' HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',  &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')





! **********************************************************************************************************************************

      END SUBROUTINE TSET_PROC_FOR_RIGELS



      SUBROUTINE TSET_PROC_FOR_SPCS ( IERRT )

! DOF Processor for SPC's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1H, L1O, L1O_MSG, LINK1O
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ENFORCED, FATAL_ERR, LSPCADDC, NDOFSB, NDOFSE, NDOFSG, NGRID, NSPCADD,      &
                                         NUM_SPC_RECORDS, NUM_SPC1_RECORDS, NUM_SPCSIDS
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, SPCADD_SIDS, SPCSET, SPCSIDS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF

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

j_do6:   DO J=1,NUM_SPCSIDS

            IF (SETID == SPCSIDS(J)) THEN
                                                           ! Check for existence of grid pt
               CALL RDOF ( ICOMP, CDOF )                   ! Convert ICOMP to CDOF char form for use below

               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID1, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', GRID1, 'SPC OR SPC1', SPCSIDS(J)
                  WRITE(F06,1822) 'GRID ', GRID1, 'SPC OR SPC1', SPCSIDS(J)
               ENDIF
               IF (GRID2 /= GRID1) THEN
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID2, GRID_ID_ROW_NUM )
                  IF (GRID_ID_ROW_NUM == -1) THEN
                     GID_ERR = GID_ERR + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1822) 'GRID ', GRID2, 'SPC OR SPC1', SPCSIDS(J)
                     WRITE(F06,1822) 'GRID ', GRID2, 'SPC OR SPC1', SPCSIDS(J)
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
                                 IF (TSET(GRID_ID_ROW_NUM,K) == '  ') THEN
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
                              ELSE IF (DOFSET == 'SB') THEN
                                 IF ((TSET(GRID_ID_ROW_NUM,K) == '  ') .OR. (TSET(GRID_ID_ROW_NUM,K) == 'SG') .OR.                 &
                                     (TSET(GRID_ID_ROW_NUM,K) == 'SB')) THEN
                                    TSET(GRID_ID_ROW_NUM,K) = 'SB'
                                    NDOFSB = NDOFSB + 1
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

               EXIT j_do6                                     ! We found and processed the SPC set, so quit

            ENDIF

         ENDDO j_do6

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


      SUBROUTINE TSET_PROC_FOR_SUPORTS ( IERRT )

! DOF Processor for SUPORT's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1T, L1T_MSG, LINK1T
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFR, NGRID, NUM_SUPT_CARDS
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TSET_PROC_FOR_SUPTS'
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! The name of a DOF set (e.g. 'SB', 'A ', etc)

      INTEGER(LONG), INTENT(INOUT)    :: IERRT             ! Sum of all grid and DOF errors
      INTEGER(LONG)                   :: GRID_NUM          ! An actual grid ID
      INTEGER(LONG)                   :: DOF_ERR   = 0     ! Count of errors that result from setting displ sets in TSET
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: I,K               ! DO loop indices
      INTEGER(LONG)                   :: ICOMP             ! DOF components read from file LINK1T (SUPORT's)
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO    = 0     ! Record number when reading a file




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! Process SPC data from file L1T (data written when SUPORT Bulk Data cards were read)

      GID_ERR = 0
      DOF_ERR = 0
      DOFSET = 'R '
      REC_NO = 0

i_do6:DO I=1,NUM_SUPT_CARDS

         READ(L1T,IOSTAT=IOCHK) GRID_NUM, ICOMP
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1T, L1T_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Error reading SUPORT data file . No sense continuing
         ENDIF

         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, GRID_ID_ROW_NUM )
         IF (GRID_ID_ROW_NUM == -1) THEN
            GID_ERR = GID_ERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1310) 'GRID ', GRID_NUM, 'SUPORT'
            WRITE(F06,1310) 'GRID ', GRID_NUM, 'SUPORT'
         ENDIF

         IF (GID_ERR == 0) THEN
            CALL RDOF ( ICOMP, CDOF )                      ! Convert ICOMP to CDOF char form for use below
            DO K = 1,6                                     ! Put data in TSET
               IF (CDOF(K) == '1') THEN
                  IF      (DOFSET == 'R ') THEN
                     IF (TSET(GRID_ID_ROW_NUM,K) == 'A ') THEN
                        TSET(GRID_ID_ROW_NUM,K) = 'R '
                        NDOFR = NDOFR + 1
                     ELSE
                        DOF_ERR = DOF_ERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1334) GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                        WRITE(F06,1334) GRID_NUM, K, DOFSET, TSET(GRID_ID_ROW_NUM,K)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

      ENDDO i_do6

      IERRT = IERRT + GID_ERR + DOF_ERR



      RETURN

! **********************************************************************************************************************************
 1310 FORMAT(' *ERROR  1310: ',A,I8,' ON ',A,' BULK DATA ENTRY IS UNDEFINED')

 1331 FORMAT(' *ERROR  1331: GRID POINT ',I8,' HAS COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',                                     &
                           ' HOWEVER THIS GRID/COMP IS ALREADY IN THE ',A2,' DISPL SET')

 1334 FORMAT(' *ERROR  1334: SUPORT B.D. ENTRY HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' DISPL SET.',                   &
                           ' HOWEVER THIS GRID/COMP IS ALREADY IN THE ',A2,' DISPL SET')

 1369 FORMAT(' *ERROR  1369: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RECORD READ FROM FILE: ',A                                                                            &
                    ,/,14X,' INDICATES THAT AN SPCd DOF BELONGS TO THE "',A2,'" SET. MUST BE EITHER "SE" OR "SB"')

! **********************************************************************************************************************************

      END SUBROUTINE TSET_PROC_FOR_SUPORTS


      SUBROUTINE RDOF ( INTDOF, CDOF )

! Convert DOF's from integer to char flag form. This version assumes that the data in INTDOF contains only digits 0-6.
! This subr is only used in situations where we know that this is the case. That is true since all DOF fields on any
! Bulk Data card are checked by subr IP6CHK for validity when the bulk data was read.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RDOF'
      CHARACTER( 8*BYTE)              :: CINT              ! 8 char field of integers in INTDOF
      CHARACTER( 1*BYTE), INTENT(OUT) :: CDOF(6)           ! Contains 1 in each of the 6 pos'ns corresponding to a DOF from INTDOF
!                                                            and zero otherwise. For example, if INTDOF = 135 then CDOF = '101010'

      INTEGER(LONG), INTENT(IN)       :: INTDOF            ! Integer field which should contain only the digits 1 - 6
      INTEGER(LONG)                   :: I




! **********************************************************************************************************************************
! Initialize CDOF

      DO I=1,6
         CDOF(I) = '0'
      ENDDO

! Write INTDOF to CINT

      WRITE(CINT,'(I8)') INTDOF

! Process CINT to form CDOF

      DO I = 1,8
         IF      (CINT(I:I) == '0') THEN                   ! CINT = 0 is SPOINT so make it comp 1 same as T1 for physical grid
            CDOF(1) = '1'
         ELSE IF (CINT(I:I) == '1') THEN
            CDOF(1) = '1'
         ELSE IF (CINT(I:I) == '2') THEN
            CDOF(2) = '1'
         ELSE IF (CINT(I:I) == '3') THEN
            CDOF(3) = '1'
         ELSE IF (CINT(I:I) == '4') THEN
            CDOF(4) = '1'
         ELSE IF (CINT(I:I) == '5') THEN
            CDOF(5) = '1'
         ELSE IF (CINT(I:I) == '6') THEN
            CDOF(6) = '1'
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE RDOF

   END MODULE DOF_SET_CONSTRUCTION
