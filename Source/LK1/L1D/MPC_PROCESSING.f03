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

   MODULE MPC_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MPC_PROC

   CONTAINS

      SUBROUTINE MPC_PROC

! Processes MPC equations to get terms for the RMG constraint matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     L1J,     L1S, LINK1S, L1S_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LMPCADDC, NGRID, NMPC, NMPCADD, NUM_MPCSIDS, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID_ID, MPCSET, MPCSIDS
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MPC_PROC'
      CHARACTER( 1*BYTE)              :: MPC_SET_USED      ! 'Y'/'N' indicator if an MPC set in B.D. is used

      INTEGER(LONG)                   :: COMP              ! Comp number for a grid for an MPC eqn (read from file LINK1S)
      INTEGER(LONG)                   :: COMP_JUNK         ! A displ component read from file LINK1S that we do not need
      INTEGER(LONG)                   :: GID               ! Grid for an MPC read from file LINK1S
      INTEGER(LONG)                   :: GID_JUNK          ! An grid read from file LINK1S that we do not need
!xx   INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no., in TDOF array, of the G-set DOF list
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: M_SET_COL_NUM     ! Col no., in TDOF array, of the M-set DOF list
      INTEGER(LONG)                   :: NUM_TRIPLES       ! Counter on number of pairs of grid/comp/coeff triplets on an MPC
!                                                            logical card. Must be <= MMPC which was counted in subr BD_MPC0
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: RMG_COL_NUM       ! Col no. of a term in array RMG
      INTEGER(LONG)                   :: RMG_ROW_NUM       ! Row no. of a term in array RMG
      INTEGER(LONG)                   :: ROW_NUM           ! A row number in array TDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: SETID             ! An SPC set ID read from file LINK1O


      REAL(DOUBLE)                    :: COEFF             ! An MPC coeff value read from file LINK1S that we do not need
      REAL(DOUBLE)                    :: COEFF_JUNK        ! An MPC coeff value read from file LINK1S that we do not need



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

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )

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
j_do3:   DO J=1,NUM_MPCSIDS

            IF (SETID == MPCSIDS(J)) THEN

               MPC_SET_USED = 'Y'

               READ(L1S,IOSTAT=IOCHK) GID,COMP,COEFF       ! Read dependent grid/comp/coeff
               REC_NO = REC_NO + 1
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
                                                           ! Get row num (in GRID_ID) corresponding to grid GID (we know GID exists)
!xx            CALL CALC_TDOF_ROW_NUM ( GID, ROW_NUM_START, 'N' )
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID, IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
                                                           ! Determine the row and col for COEFF in RMG and write to L1J
               ROW_NUM = ROW_NUM_START + COMP - 1
               RMG_ROW_NUM = TDOF(ROW_NUM,M_SET_COL_NUM)
               RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)
               IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
                  WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,COEFF
               ELSE IF (RMG_ROW_NUM  <= 0) THEN
                  WRITE(ERR,1517) SUBR_NAME,GID,COMP,RMG_ROW_NUM
                  WRITE(F06,1517) SUBR_NAME,GID,COMP,RMG_ROW_NUM
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ELSE IF (RMG_COL_NUM <= 0) THEN
                  WRITE(ERR,1513) SUBR_NAME,GID,COMP,RMG_COL_NUM
                  WRITE(F06,1513) SUBR_NAME,GID,COMP,RMG_COL_NUM
                  FATAL_ERR = FATAL_ERR + 1
                 CALL OUTA_HERE ( 'Y' )
               ENDIF

               DO K=1,NUM_TRIPLES-1                        ! Cycle over other terms in MPC eqn. Use above RMG_ROW_NUM & write to L1J

                  READ(L1S,IOSTAT=IOCHK) GID,COMP,COEFF

                  REC_NO = REC_NO + 1
                  IF (IOCHK /= 0) THEN
                     CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF

!xx05/05/07       CALL GET_ARRAY_ROW_NUM ( '', SUBR_NAME, NGRID, GRID_ID, GID, GRID_ID_ROW_NUM )
!xx05/05/07          IF (GRID_ID_ROW_NUM == -1) THEN
!xx05/05/07             WRITE(ERR,1822) 'GRID', GID, 'MPC', MPCSIDS(J)
!xx05/05/07             WRITE(F06,1822) 'GRID', GID, 'MPC', MPCSIDS(J)
!xx05/05/07             FATAL_ERR = FATAL_ERR + 1
!xx05/05/07             CYCLE j_do3
!xx05/05/07          ENDIF
!xx               CALL CALC_TDOF_ROW_NUM ( GID, ROW_NUM_START, 'N' )
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GID, IGRID )
                  ROW_NUM_START = TDOF_ROW_START(IGRID)

!xx05/05/07       ROW_NUM = 6*(GRID_ID_ROW_NUM - 1) + COMP
                  ROW_NUM = ROW_NUM_START + COMP - 1
                  RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)
                  IF (RMG_COL_NUM > 0) THEN
                     WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,COEFF
                  ELSE
                     WRITE(ERR,1513) SUBR_NAME,GID,COMP,RMG_COL_NUM
                     WRITE(F06,1513) SUBR_NAME,GID,COMP,RMG_COL_NUM
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF

               ENDDO

               EXIT j_do3

            ELSE

               MPC_SET_USED = 'N'
               CYCLE j_do3

            ENDIF

         ENDDO j_do3

         IF (MPC_SET_USED == 'N') THEN                     ! This MPC set is not to be used, so skip all grid/comp/coeff records
            DO K=1,NUM_TRIPLES
               READ(L1S,IOSTAT=IOCHK) GID_JUNK,COMP_JUNK,COEFF_JUNK
               REC_NO = REC_NO + 1
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1S, L1S_MSG, REC_NO, OUNT )
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDDO
         ENDIF

      ENDDO i_do3

      ! do not deallocate if we're on the first step of a buckling sol
      IF ((SOL_NAME(1:8) /= 'BUCKLING') .OR. (LOAD_ISTEP /= 1)) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'MPCSIDS' )
      END IF



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1517 FORMAT(' *ERROR  1517: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ROW NUMBER IN ARRAY RMG CALCULATED FOR GRID ',I8,', COMPONENT ',I2,' MUST BE > 0 BUT IS = ',I8)

 1513 FORMAT(' *ERROR  1513: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' COL NUMBER IN ARRAY RMG CALCULATED FOR GRID ',I8,', COMPONENT ',I2,' MUST BE > 0 BUT IS = ',I8)





! **********************************************************************************************************************************

      END SUBROUTINE MPC_PROC

   END MODULE MPC_PROCESSING
