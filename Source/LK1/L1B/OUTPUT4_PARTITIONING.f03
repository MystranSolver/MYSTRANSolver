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

   MODULE OUTPUT4_PARTITIONING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OU4_PARTVEC_PROC

   CONTAINS

      SUBROUTINE OU4_PARTVEC_PROC ( INDEX, OU4_MAT_NAME, NROWS_F, NCOLS_F, ROW_SET, COL_SET, CAN_PARTN, NROWS_P, NCOLS_P,          &
                                    VAL_ROWS, VAL_COLS )

! Calculates partitioning vectors OU4_PARTVEC_ROW and OU4_PARTVEC_COL used to partition OUTPUT4 matrices.

! The Exec Control may have requests for OUTPUT4 matrices and there may be PARTN requests there to partition particular OUTPUT4
! matrices to output a subset of the rows and/or cols of that matrix.

! The PARTN requests in Exec Control have partitioning vector names with the specific grid/component partitioning data defined on
! Bulk Data PARVEC and/or PARVEC1 entries.

! The partitioning vector data from the Bulk Data PARVEC,1 entries is written to file LINK1V.

! This subr puts the data from file LINK1V into a table (PSET) for each OUTPUT4 matrix to be partitioned. This PSET table has
! 6 cols (1 for each of the 6 components of motion: T1, T2, T3, R1, R2, R3) and as many rows as there are grids in the model.

! Each OUTPUT4 matrix to be partitioned has a PSET table associated with it. The PSET table is then used to create the row and col
! partitioning vectors OU4_PARTVEC_ROW and OU4_PARTVEC_COL which are used in the calling program as input to the subr which does
! the actual partitioning (subr PARTITION_SS)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, L1V, L1V_MSG, LINK1V
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MTSET, NDOFG, NGRID, NUM_PARTVEC_RECORDS, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, TDOF, TDOFI, TDOF_ROW_START
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, OU4_PART_VEC_NAMES, OU4_PARTVEC_COL, OU4_PARTVEC_ROW,              &
                                         OU4_MAT_ROW_GRD_COMP, OU4_MAT_COL_GRD_COMP

      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, INV_GRID_SEQ
      USE PARAMS, ONLY                :  PRTPSET, SUPINFO

      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN, READERR
      USE DOF_SET_CONSTRUCTION, ONLY: RDOF
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY: TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN(BLNK_SUB_NAM))                     :: SUBR_NAME = 'OU4_PARTVEC_PROC'

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES)),INTENT(IN) :: OU4_MAT_NAME
                                                            ! Name of OU4 matrix to be partitioned

      CHARACTER(LEN(TSET_CHR_LEN))         ,INTENT(IN) :: COL_SET
                                                            ! Displ set for the cols of matrix whose name is OU4_MAT_NAME

      CHARACTER(LEN(TSET_CHR_LEN))         ,INTENT(IN) :: ROW_SET
                                                            ! Displ set for the rows of matrix whose name is OU4_MAT_NAME

      CHARACTER( 1*BYTE)                   ,INTENT(OUT):: CAN_PARTN
                                                            ! 'Y' if matrix should be part'd (rows, cols or both)

      CHARACTER( 1*BYTE)                               :: CDOF(6)
                                                            ! An output from subr RDOF

      CHARACTER( 3*BYTE), PARAMETER                    :: COL_OR_ROW(2) = (/'COL', 'ROW'/)
                                                            ! Designator for PSET table heading

      CHARACTER(LEN(TSET_CHR_LEN))                     :: COL_AND_ROW_SET(2)
                                                            ! COL_SET and ROW_SET in an array

      CHARACTER( 1*BYTE)                               :: FOUND_PART_VEC(2)
                                                            ! 'Y' if we found the col and row partitioning vectors

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES))            :: PART_VEC_NAME(2)
                                                            ! The name of the col and row  partitioning vectors

      CHARACTER(LEN(TSET_CHR_LEN)), PARAMETER          :: PSET_CHAR(2)  = (/'CP', 'RP'/)
                                                            ! Char to enter into PSET table for either col or row


      CHARACTER(LEN=MTSET)                             :: PSET(NGRID,MTSET,2)
                                                            ! Tables, like TSET and USET but for partitioning vectors and only
!                                                             in the current pass of this subr (i.e. it will get regenerated each
!                                                             time this subr is called with diff vals, depending on the row/col
!                                                             partitioning vecs being processed. The entries in this table will be
!                                                             either "CP" (for the col partitioning DOF's) or "RP" for the rows

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES))            :: VNAME
                                                            ! The name of a part vector from a PARTN E.C. entry


      INTEGER(LONG), INTENT(IN)       :: INDEX              ! Index into array ACT_OU4_MYSTRAN_NAMES where OU4_MAT_NAME exists
      INTEGER(LONG), INTENT(IN)       :: NCOLS_F            ! Number of cols in the complete input OUTPUT4 matrix OU4_MAT_NAME
      INTEGER(LONG), INTENT(IN)       :: NROWS_F            ! Number of cols in the complete input OUTPUT4 matrix OU4_MAT_NAME
      INTEGER(LONG), INTENT(OUT)      :: NCOLS_P            ! Number of cols that will be in the OUTPUT4 matrix when partitioned
      INTEGER(LONG), INTENT(OUT)      :: NROWS_P            ! Number of cols that will be in the OUTPUT4 matrix when partitioned

      INTEGER(LONG)                   :: GRID1              ! An actual grid ID
      INTEGER(LONG)                   :: GRID2              ! An actual grid ID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM    ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GID_ERR   = 0      ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: I,J,K,GRID_NUM     ! DO loop indices
      INTEGER(LONG)                   :: ICOMP              ! DOF components read from file LINK1V (SPC's) or LINK1N (ASET/OMIT's)
      INTEGER(LONG)                   :: IDOF               ! Row in vectors OU4_PARTVEC_COL, OU4_PARTVEC_ROW
      INTEGER(LONG)                   :: IERR_ALLOC(2)      ! Indicator of errors allocating arrays OU4_PARTVEC_COL, OU4_PARTVEC_ROW
      INTEGER(LONG)                   :: IERR_CHK           ! Number of errors found when checking sensibility of partitioning vecs

      INTEGER(LONG)                   :: IERR_PART(2)       ! Error indicator if matrix does not have cols or rows that are not grid
!                                                             displ set oriented

      INTEGER(LONG)                   :: IERRT              ! Count of total number of errors found here
      INTEGER(LONG)                   :: IGRID              ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK              ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: IROW               ! Row number in array TDOF
      INTEGER(LONG)                   :: NDIM_F(2)          ! Array of NCOLS_F, NROWS_F
      INTEGER(LONG)                   :: NUM_PSET_ENTRIES   ! Number of entries made into the PSET table
      INTEGER(LONG)                   :: NUM_COMPS          ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: OUNT(2)            ! File units to write messages to.
      INTEGER(LONG)                   :: PSET_ERR   = 0     ! Count of errors that result from setting displ sets in PSET
      INTEGER(LONG)                   :: REC_NO     = 0     ! Record number when reading a file
      INTEGER(LONG)                   :: TDOF_COL(2)        ! Col number in TDOF/TDOFI where CHAR_SET set data exists
      INTEGER(LONG), INTENT(OUT)      :: VAL_COLS           ! Number to enter into PARTVEC_COL for a col that is to be partitioned
      INTEGER(LONG), INTENT(OUT)      :: VAL_ROWS           ! Number to enter into PARTVEC_ROW for a row that is to be partitioned




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! ----------------------------------------------------------------------------------------------------------------------------------
! Initialize

      CALL DEALLOCATE_PART_VECS
      CALL DEALLOCATE_COL_VEC ( 'OU4_MAT_COL_GRD_COMP' )   ! Do this here instead of at end so calc'd array will be returned
      CALL DEALLOCATE_COL_VEC ( 'OU4_MAT_ROW_GRD_COMP' )   ! Do this here instead of at end so calc'd array will be returned

      VAL_COLS = 2
      VAL_ROWS = 1

      CAN_PARTN = 'Y'

      GID_ERR  = 0
      IERRT    = 0
      PSET_ERR = 0
      REC_NO   = 0

      COL_AND_ROW_SET(1) = COL_SET
      COL_AND_ROW_SET(2) = ROW_SET

      NDIM_F(1) = NCOLS_F
      NDIM_F(2) = NROWS_F

      NCOLS_P   = 0
      NROWS_P   = 0

! If both COL_SET and ROW_SET are not displ sets then we cannot partition matrix OU4_MAT_NAME using this subr. Therefore construct
! partitioning vectors that have all nonzero terms. These will be "partitioning" vectors that include all cols and rows in the
! partition. When we return from this subr this matrix will not actually be partitioned into itself but we will, at least, have a
! set of partitioning vectors to return. These null vectors will not be used since CAN_PARTN = 'N' is returned also

      IF ((COL_SET == '--') .AND. (ROW_SET == '--')) THEN  ! If both are null we are finished generating vectors for OU4_MATRIX_NAME
         CALL ALLOCATE_PART_VECS ( IERR_ALLOC, IERR_PART, TDOF_COL )
         WRITE(ERR,101) OU4_MAT_NAME
         WRITE(F06,101) OU4_MAT_NAME
         CAN_PARTN = 'N'
         NCOLS_P = NCOLS_F
         NROWS_P = NROWS_F
         CALL CHECK_PART_VECS ( IERR_CHK )
         RETURN
      ENDIF

! If either set is not partitionable, set the output rows/cols to the input ones

      IF      (COL_SET == '--') THEN
         NCOLS_P = NCOLS_F
      ELSE IF (ROW_SET == '  ') THEN
         NROWS_P = NROWS_F
      ENDIF

! Check consistency of INDEX and OU4_MAT_NAME

      IF (OU4_MAT_NAME /= ACT_OU4_MYSTRAN_NAMES(INDEX)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1335) SUBR_NAME, INDEX, ACT_OU4_MYSTRAN_NAMES(INDEX), OU4_MAT_NAME
         WRITE(F06,1335) SUBR_NAME, INDEX, ACT_OU4_MYSTRAN_NAMES(INDEX), OU4_MAT_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Get the names of the column and row partitioning vectors for matrix OU4_MAT_NAME.

      PART_VEC_NAME(1) = OU4_PART_VEC_NAMES(INDEX,1)       ! Name of col partitioning vector for this OU4_MAT_NAME
      PART_VEC_NAME(2) = OU4_PART_VEC_NAMES(INDEX,2)       ! Name of row partitioning vector for this OU4_MAT_NAME

! Allocate memory for, and initialize, arrays OU4_PARTVEC_COL, OU4_PARTVEC_ROW. Initialization to 1's will mean
! that the partitioning vector, unless modified below, will "partition" all rows and cols of matrix OU4_MAT_NAME

      CALL ALLOCATE_PART_VECS ( IERR_ALLOC, IERR_PART, TDOF_COL )

      DO I=1,2
         IF (IERR_ALLOC(I) > 0) THEN
            WRITE(ERR,1336) PART_VEC_NAME(I), SUBR_NAME
            WRITE(F06,1336) PART_VEC_NAME(I), SUBR_NAME
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO
      IF ((IERR_ALLOC(1) > 0) .OR. (IERR_ALLOC(2) > 0)) THEN
         WRITE(ERR,9999) SUBR_NAME
         WRITE(F06,9999) SUBR_NAME
         CAN_PARTN = 'N'
         CALL DEALLOCATE_PART_VECS
         RETURN
      ENDIF

! Process PARTVEC data from file L1V (data written when PARTVEC and PARTVEC1 Bulk Data cards were read)
! Scan file L1V to find the records that have the CP_NAME and RP_NAME partitioning vector data

i_do1:DO I=1,2                                             ! J=1 is for the col part. vector and J=2 is for the row part. vector

         IF (PART_VEC_NAME(I)(1:) == ' ') CYCLE i_do1

         NUM_PSET_ENTRIES = 0

         DO J=1,NGRID
            DO K=1,MTSET
               PSET(J,K,I)(1:) = '--'
            ENDDO
         ENDDO

         FOUND_PART_VEC(I) = 'N'

         CALL FILE_OPEN ( L1V, LINK1V, OUNT, 'OLD', L1V_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
j_do1:   DO J=1,NUM_PARTVEC_RECORDS

            READ(L1V,IOSTAT=IOCHK) VNAME, ICOMP, GRID1, GRID2
            REC_NO = REC_NO + 1
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK1V, L1V_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading PARTVEC file . No sense continuing
            ENDIF

            IF (VNAME == PART_VEC_NAME(I)) THEN

               FOUND_PART_VEC(I) = 'Y'

               CALL RDOF ( ICOMP, CDOF )                   ! Convert ICOMP to CDOF char form for use below

               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID1, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM == -1) THEN
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', GRID1, 'PARTVEC OR PARTVEC1', VNAME
                  WRITE(F06,1822) 'GRID ', GRID1, 'PARTVEC OR PARTVEC1', VNAME
               ENDIF
               IF (GRID2 /= GRID1) THEN
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID2, GRID_ID_ROW_NUM )
                  IF (GRID_ID_ROW_NUM == -1) THEN
                     GID_ERR = GID_ERR + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1822) 'GRID ', GRID1, 'PARTVEC OR PARTVEC1', VNAME
                     WRITE(F06,1822) 'GRID ', GRID1, 'PARTVEC OR PARTVEC1', VNAME
                  ENDIF
               ENDIF

               IF (GID_ERR == 0) THEN                      ! Put CDOF data in PSET for GRID1 thru GRID2
                  DO GRID_NUM=GRID1,GRID2                  ! GRID2 >= GRID1 was checked in subr BD_SPC1
                     CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, GRID_ID_ROW_NUM )
                     IF (GRID_ID_ROW_NUM /= -1) THEN
                        CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM, NUM_COMPS, SUBR_NAME )
                        DO K = 1,NUM_COMPS                 ! Put data in PSET
                           IF (CDOF(K) == '1') THEN
                              IF ((PSET(GRID_ID_ROW_NUM,K,I) == '--') .OR. (PSET(GRID_ID_ROW_NUM,K,I) == PSET_CHAR(I))) THEN
                                 PSET(GRID_ID_ROW_NUM,K,I) = PSET_CHAR(I)
                                 NUM_PSET_ENTRIES = NUM_PSET_ENTRIES + 1
                              ELSE
                                 PSET_ERR  = PSET_ERR  + 1
                                 FATAL_ERR = FATAL_ERR + 1
                                 WRITE(ERR,1332) VNAME, GRID_NUM, K, PSET_CHAR(I), PSET(GRID_ID_ROW_NUM,K,I)
                                 WRITE(F06,1332) VNAME, GRID_NUM, K, PSET_CHAR(I), PSET(GRID_ID_ROW_NUM,K,I)
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ELSE

                  CYCLE j_do1

               ENDIF

            ENDIF

         ENDDO j_do1
         CALL FILE_CLOSE ( L1V, LINK1V, 'KEEP' )

         IF (FOUND_PART_VEC(I) == 'N') THEN

            WRITE(ERR,102) PART_VEC_NAME(I), OU4_MAT_NAME, COL_OR_ROW(I)
            WRITE(F06,102) PART_VEC_NAME(I), OU4_MAT_NAME, COL_OR_ROW(I)

         ELSE

            IF (IERR_PART(I) == 0) THEN
               IF (COL_AND_ROW_SET(I) /= '--') THEN        ! Calc part vec if set is not '--'. If '--', it was done in the ALLOCATE
                  IF (PRTPSET > 0) THEN
                     WRITE(F06,56) COL_OR_ROW(I), OU4_MAT_NAME, NUM_PSET_ENTRIES
                     WRITE(F06,57)
                     DO K=1,NGRID
                        WRITE(F06,58) GRID(K,1), (PSET(K,J,I),J=1,MTSET)
                     ENDDO
                     WRITE(F06,*)
                  ENDIF

                  IF (NDIM_F(I) > 0) THEN
                     DO K=1,NGRID
                        IGRID = INV_GRID_SEQ(K)
                        CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
                        DO J=1,NUM_COMPS
                           IF (PSET(IGRID,J,I) == PSET_CHAR(I)) THEN
                              IROW = TDOF_ROW_START(IGRID) + J - 1
                              CALL TDOF_COL_NUM ( COL_AND_ROW_SET(I), TDOF_COL(I) )
                              IDOF = TDOF(IROW,TDOF_COL(I))
                              IF (IDOF > 0) THEN
                                 IF (I == 1) THEN
                                    OU4_PARTVEC_COL(IDOF) = VAL_COLS
                                    NCOLS_P = NCOLS_P + 1
                                 ELSE
                                    OU4_PARTVEC_ROW(IDOF) = VAL_ROWS
                                    NROWS_P = NROWS_P + 1
                                 ENDIF
                              ELSE
                                 WARN_ERR = WARN_ERR + 1
                                 WRITE(ERR,103) PART_VEC_NAME(I), OU4_MAT_NAME, GRID_ID(IGRID), ICOMP, COL_AND_ROW_SET(I)
                                 WRITE(F06,103) PART_VEC_NAME(I), OU4_MAT_NAME, GRID_ID(IGRID), ICOMP, COL_AND_ROW_SET(I)
                                 CAN_PARTN = 'N'
                                 CALL DEALLOCATE_PART_VECS
                                 RETURN
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF

               ELSE
                  WRITE(ERR,104) OU4_MAT_NAME, COL_OR_ROW(I), COL_OR_ROW(I)
                  WRITE(F06,104) OU4_MAT_NAME, COL_OR_ROW(I), COL_OR_ROW(I)
                  IF (I == 1) THEN
                     NCOLS_P = NCOLS_F
                  ELSE
                     NROWS_P = NROWS_F
                  ENDIF

               ENDIF

            ENDIF

         ENDIF

      ENDDO i_do1

      IF (DEBUG(114) > 0) CALL OU4_PARTVEC_PROC_DEB ( '1' )

! Check NCOLS_P, NROWS_P vs the number of nonzero entries in the partitioning vectors

      CALL CHECK_PART_VECS ( IERR_CHK )

! Check for other errors

      IERRT = GID_ERR + PSET_ERR + IERR_ALLOC(1) + IERR_ALLOC(2) + IERR_CHK
      IF (IERRT > 0) THEN
         WRITE(ERR,9999) IERRT, SUBR_NAME
         WRITE(F06,9999) IERRT, SUBR_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Calc arrays that will have the grid/comp headers for rows and cols of the partitioned matrix

      CALL GET_OU4_PART_GRD_COMP



      RETURN

! **********************************************************************************************************************************
   56 FORMAT(36X,'PSET FOR THE ',A,' PARTITIONING OF MATRIX "',A,'" WITH',I9,' ENTRIES'                                            &
          ,/,36X,'--------------------------------------------------------------------------------',/)

   57 FORMAT(42x,'     Grid       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(42X,I9,32767(7X,A2))

  101 FORMAT(' *WARNING    : OUTPUT4 MATRIX "',A,'" IS NOT AVAILABLE FOR PARTITIONING. IT DOES NOT HAVE EITHER COLS OR ROWS THAT', &
                           ' ARE GRID/COMP SET ORIENTED. IT WILL NOT BE PARTITIONED')

  102 FORMAT(' *WARNING    : PARVEC VECTOR "',A,'" FOR PARTITIONING OUTPUT4 MATRIX "',A,'" COULD NOT BE FOUND. MATRIX ',A,         &
                           'S WILL NOT BE PARTITIONED')

  103 FORMAT(' *WARNING    : PARTITIONING VECTOR "',A,'" FOR PARTITIONING OUTPUT4 MATRIX "',A,'" REFERENCES GRID/COMPS ',          &
       I8,'/',I6,'.',/,14X,' AT LEAST ONE OF THESE COMPONENTS IS NOT IN THE "',A,'" SET. THIS MATRIX WILL NOT BE PARTITIONED')

  104 FORMAT(' *INFORMATION: OUTPUT4 MATRIX "',A,'" DOES NOT HAVE ',A,'S THAT BELONG TO A GRID ORIENTED DISPL SET.',               &
                           ' THE ',A,'S OF THIS MATRIX WILL NOT BE PARTITIONED')

 1322 FORMAT(' *ERROR  1322: PROCESSING STOPPED DUE TO THE ABOVE LISTED ',I8,' ERRORS IN SUBR ',A)

 1331 FORMAT(' *ERROR  1331: GRID POINT ',I8,' HAS COMPONENT ',I2,' IN THE ',A2,' DISPL SET. (PERM SPC ON GRID CARD)',             &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')

 1332 FORMAT(' *ERROR  1332: PARTITIONING VECTOR ',A,' HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE "',A2,'" PSET.',              &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE "',A2,'" PSET')

 1335 FORMAT(' *ERROR  1335: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' FOR "INDEX" =',I3,' THE NAME OF THE OU4 MATRIX TO BE PARTITIONED SHOULD BE "',A,'" BUT IS "',A,'"')

 1336 FORMAT(' *ERROR  1336: ERROR ALLOCATING MEMORY TO PARTITION VECTOR "',A,'" IN SUBR ',A)

 1369 FORMAT(' *ERROR  1369: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RECORD READ FROM FILE: ',A                                                                            &
                    ,/,14X,' INDICATES THAT AN SPCd DOF BELONGS TO THE "',A2,'" SET. MUST BE EITHER "SE" OR "SB"')

 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,2X,A,' IS UNDEFINED')

 9999 FORMAT(' PROCESSING STOPPED DUE TO ABOVE ',I8,' FATAL ERRORS IN SUBR ',A)



! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ALLOCATE_PART_VECS ( IERR_ALLOC, IERR_PART, TDOF_COL )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: IERR_ALLOC(2)     ! Tot count of errors allocating arrays OU4_PARTVEC_COL, OU4_PARTVEC_ROW
      INTEGER(LONG), INTENT(OUT)      :: IERR_PART(2)      ! Error indicator if matrix does not have cols or rows that are not grid
!                                                            displ set oriented

      INTEGER(LONG), INTENT(OUT)      :: TDOF_COL(2)       ! Col number in TDOF/TDOFI where CHAR_SET set data exists
      INTEGER(LONG)                   :: II                ! DO loop index

! **********************************************************************************************************************************
      IERR_PART(1)   = 0   ;   IERR_PART(2)  = 0
      IERR_ALLOC(1)  = 0   ;   IERR_ALLOC(2) = 0

      TDOF_COL(1) = -1
      ALLOCATE ( OU4_PARTVEC_COL(NDIM_F(1)), STAT=IERR_ALLOC(1))
      DO II=1,NDIM_F(1)
         IF (COL_SET == '--') THEN
            OU4_PARTVEC_COL(II) = VAL_COLS                 ! All cols will be "partitioned" if COL_SET = '--'
         ELSE
            OU4_PARTVEC_COL(II) = 0                        ! If COL_SET /= '--', initialize to 0. Col part vec will be calc'd later
         ENDIF
      ENDDO

      TDOF_COL(2) = -1
      ALLOCATE ( OU4_PARTVEC_ROW(NDIM_F(2)), STAT=IERR_ALLOC(2))
      DO II=1,NDIM_F(2)
         IF (ROW_SET == '--') THEN
            OU4_PARTVEC_ROW(II) = VAL_ROWS                 ! All rows will be "partitioned" if ROW_SET = '--'
         ELSE
            OU4_PARTVEC_ROW(II) = 0                        ! If ROW_SET /= '--', initialize to 0. Row part vec will be calc'd later
         ENDIF
      ENDDO

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_PART_VECS

! ##################################################################################################################################

      SUBROUTINE DEALLOCATE_PART_VECS

      USE PENTIUM_II_KIND, ONLY       :  LONG

      IMPLICIT NONE

      INTEGER(LONG)                   :: IERR_DEALLOC      ! Error indicator in any 1 of the allocate statements

! **********************************************************************************************************************************
      IF (ALLOCATED ( OU4_PARTVEC_COL )) THEN
         DEALLOCATE ( OU4_PARTVEC_COL, STAT=IERR_DEALLOC )
         IF (IERR_DEALLOC > 0) THEN
            WRITE(F06,1337) IERR_DEALLOC, 'OU4_PARTVEC_COL', 'DEALLOCATE_PART_VECS'
            IF (SUPINFO /= 'N') THEN
               WRITE(F06,1337) IERR_DEALLOC, 'OU4_PARTVEC_COL', 'DEALLOCATE_PART_VECS'
            ENDIF
         ENDIF
      ENDIF

      IF (ALLOCATED ( OU4_PARTVEC_ROW )) THEN
         DEALLOCATE ( OU4_PARTVEC_ROW, STAT=IERR_DEALLOC )
         IF (IERR_DEALLOC > 0) THEN
            WRITE(F06,1337) IERR_DEALLOC, 'OU4_PARTVEC_ROW', 'DEALLOCATE_PART_VECS'
            IF (SUPINFO /= 'N') THEN
               WRITE(F06,1337) IERR_DEALLOC, 'OU4_PARTVEC_ROW', 'DEALLOCATE_PART_VECS'
            ENDIF
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1337 FORMAT(' *ERROR  1337: ERROR = ',I8,' DEALLOCATING ARRAY ',A,' IN SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_PART_VECS

! ##################################################################################################################################

      SUBROUTINE GET_OU4_PART_GRD_COMP

! Calcs the grid and comp numbers that represent the rows and cols of the partitioned matrix. This will be used when the partitioned
! matrix is written in subrs WRITE_OU4_FULL_MAT, WRITE_OU4_SPARSE_MAT as info for the user.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE DOF_TABLES, ONLY            :  TDOFI

      IMPLICIT NONE

      INTEGER(LONG)                   :: II,KK,LL


! **********************************************************************************************************************************

      IF (COL_SET /= '--') THEN

         CALL ALLOCATE_COL_VEC ( 'OU4_MAT_COL_GRD_COMP', NCOLS_P, SUBR_NAME )

         CALL TDOF_COL_NUM ( COL_SET, TDOF_COL(1) )
         KK = 0
         LL = 0
         DO II=1,NDOFG
            IF (TDOFI(II,TDOF_COL(1)) > 0) THEN
               KK = KK + 1
               IF (OU4_PARTVEC_COL(KK) == VAL_COLS) THEN
                  LL = LL + 1
                  OU4_MAT_COL_GRD_COMP(LL,1) = TDOFI(II,1)
                  OU4_MAT_COL_GRD_COMP(LL,2) = TDOFI(II,2)
               ENDIF
            ENDIF
         ENDDO
                                                           ! Alloc and null array (but won't be used)


      ELSE

         CALL ALLOCATE_COL_VEC ( 'OU4_MAT_COL_GRD_COMP', NCOLS_F, SUBR_NAME )

      ENDIF

      IF (ROW_SET /= '--') THEN

         CALL ALLOCATE_COL_VEC ( 'OU4_MAT_ROW_GRD_COMP', NROWS_P, SUBR_NAME )

         CALL TDOF_COL_NUM ( ROW_SET, TDOF_COL(1) )
         KK = 0
         LL = 0
         DO II=1,NDOFG
            IF (TDOFI(II,TDOF_COL(1)) > 0) THEN
               KK = KK + 1
               IF (OU4_PARTVEC_ROW(KK) == VAL_ROWS) THEN
                  LL = LL + 1
                  OU4_MAT_ROW_GRD_COMP(LL,1) = TDOFI(II,1)
                  OU4_MAT_ROW_GRD_COMP(LL,2) = TDOFI(II,2)
               ENDIF
            ENDIF
         ENDDO
                                                           ! Alloc and null array (but won't be used)
      ELSE

         CALL ALLOCATE_COL_VEC ( 'OU4_MAT_ROW_GRD_COMP', NROWS_F, SUBR_NAME )

      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_OU4_PART_GRD_COMP

! ##################################################################################################################################

      SUBROUTINE CHECK_PART_VECS ( IIERRT )

! Checks the partitioning vectors generated in calling subr for sensibility

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      IMPLICIT NONE

      INTEGER(LONG)                   :: II                ! DO loop index
      INTEGER(LONG)                   :: IICNT             ! Counter
      INTEGER(LONG)                   :: IIERR             ! Local error count
      INTEGER(LONG)                   :: IIERRT            ! Total error count

! **********************************************************************************************************************************
      IIERRT = 0

! Make sure OU4_PARTVEC_COL has only 0's or VAL_COLS entries and that the number of nonzeros is NCOLS_P

      IIERR = 0
      IICNT = 0
      DO II=1,NCOLS_F
         IF ((OU4_PARTVEC_COL(II) /= 0) .AND. (OU4_PARTVEC_COL(II) /= VAL_COLS)) THEN
            IIERR  = IIERR  + 1
            IIERRT = IIERRT + 1
         ENDIF
         IF (OU4_PARTVEC_COL(II) /= 0) THEN
            IICNT = IICNT + 1
         ENDIF
      ENDDO

      IF (IIERR > 0 ) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1338) SUBR_NAME, 'OU4_PARTVEC_COL', OU4_MAT_NAME, IIERR, VAL_COLS
         WRITE(F06,1338) SUBR_NAME, 'OU4_PARTVEC_COL', OU4_MAT_NAME, IIERR, VAL_COLS
      ENDIF

      IF (IICNT /= NCOLS_P) THEN
         IIERRT = IIERRT + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1339) SUBR_NAME, 'OU4_PARTVEC_COL', IICNT, OU4_MAT_NAME, 'NCOLS_P', NCOLS_P
         WRITE(F06,1339) SUBR_NAME, 'OU4_PARTVEC_COL', IICNT, OU4_MAT_NAME, 'NCOLS_P', NCOLS_P
      ENDIF

! Make sure OU4_PARTVEC_ROW has only 0's or VAL_ROWS entries and that the number of nonzeros is NROWS_P

      IIERR = 0
      IICNT = 0
      DO II=1,NROWS_F
         IF ((OU4_PARTVEC_ROW(II) /= 0) .AND. (OU4_PARTVEC_ROW(II) /= VAL_ROWS)) THEN
            IIERR  = IIERR  + 1
            IIERRT = IIERRT + 1
         ENDIF
         IF (OU4_PARTVEC_ROW(II) /= 0) THEN
            IICNT = IICNT + 1
         ENDIF
      ENDDO

      IF (IIERR > 0 ) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1338) SUBR_NAME, 'OU4_PARTVEC_ROW', OU4_MAT_NAME, IIERR, VAL_ROWS
         WRITE(F06,1338) SUBR_NAME, 'OU4_PARTVEC_ROW', OU4_MAT_NAME, IIERR, VAL_ROWS
      ENDIF

      IF (IICNT /= NROWS_P) THEN
         IIERRT = IIERRT + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1339) SUBR_NAME, 'OU4_PARTVEC_ROW', IICNT, OU4_MAT_NAME, 'NROWS_P', NROWS_P
         WRITE(F06,1339) SUBR_NAME, 'OU4_PARTVEC_ROW', IICNT, OU4_MAT_NAME, 'NROWS_P', NROWS_P
      ENDIF

! **********************************************************************************************************************************
 1338 FORMAT(' *ERROR  1138: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARTITIONING VECTOR ',A,' FOR OUTPUT4 MATRIX "',A,'" HAS ',I8,' TERMS THAT ARE NOT EITHER 0 OR',      &
                           ' VAL_COLS = ',I8,'. CANNOT CONTINUE WITH CORRUPT PARTITIONING VECTOR')

 1339 FORMAT(' *ERROR  1339: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZEROS IN PARTITIONING VECTOR ',A,' = ',I8,' FOR OUTPUT4 MATRIX "',A,                 &
                         '". IT SHOULD BE ',A,' = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE CHECK_PART_VECS

! ##################################################################################################################################

      SUBROUTINE OU4_PARTVEC_PROC_DEB ( WHICH )

      CHARACTER( 1*BYTE)              :: WHICH             ! Decides what to print out for this call to this subr
      CHARACTER( 3*BYTE)              :: DASH = '---'

! **********************************************************************************************************************************
      IF (WHICH == '1') THEN

         WRITE(F06,87001) OU4_MAT_NAME,NROWS_F,NCOLS_F

         WRITE(F06,*)
         WRITE(F06,87874) (I,I=1,NROWS_F)
         WRITE(F06,87873) (DASH,I=1,NROWS_F)
         WRITE(F06,87877) (OU4_PARTVEC_ROW(I),I=1,NROWS_F)
         WRITE(F06,*)
         WRITE(F06,87874) (I,I=1,NCOLS_F)
         WRITE(F06,87873) (DASH,I=1,NCOLS_F)
         WRITE(F06,87878) (OU4_PARTVEC_COL(I),I=1,NCOLS_F)
         WRITE(F06,*)

         WRITE(F06,87999) NROWS_P, NCOLS_P

      ENDIF

! **********************************************************************************************************************************
87001 FORMAT(' Partitioning vectors for matrix "',a,'" with NROWS_F =',i5,' rows and NCOLS_F =',i5,' cols',/,                      &
             ' ---------------------------------------------------------------------------------------------------')

87873 FORMAT(20x,32767a3)

87874 FORMAT(' Index, I       =  ',32767I3)

87877 FORMAT(' PARTVEC_ROW(I) =  ',32767I3)

87878 FORMAT(' PARTVEC_COL(I) =  ',32767I3)

87999 FORMAT(' The partitioned matrix will have NROWS_P =',i8,' rows and NCOLS_P =',i8,' cols')

! **********************************************************************************************************************************

      END SUBROUTINE OU4_PARTVEC_PROC_DEB

      END SUBROUTINE OU4_PARTVEC_PROC


   END MODULE OUTPUT4_PARTITIONING
