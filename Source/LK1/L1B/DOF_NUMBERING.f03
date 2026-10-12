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

   MODULE DOF_NUMBERING

   USE DOF_LOOKUP_UTILS, ONLY : TDOF_COL_NUM

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: DOF_PROC, TDOF_PROC, TDOF_COL_NUM

   CONTAINS

      SUBROUTINE DOF_PROC ( TDOF_MSG )

! DOF Processor

! Part 1: Generate TSET table (see subr TSET_PROC for explanation)
! ------
!    TSET is a table that the DOF set (e.g. "G ", "N ", etc) for each of the 6 components for every grid)

! Part 2: Generate USET table (see subr TSET_PROC for explanation)
! ------
!    USET is a table that the user defined set set ("U1" or "U2") for each of the 6 components for every grid)

! Part 3: Generate TDOF table from TSET and USET
! ------
!    TDOF is a table that has the DOF number for every DOF and every DOF set

! Part 4: Check for errors
! ------

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFSE, NUM_USETSTR, SOL_NAME
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_SET_CONSTRUCTION, ONLY: TSET_PROC
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DOF_PROC'
      CHARACTER(LEN=*), INTENT(IN)    :: TDOF_MSG          ! Message to be printed out regarding at what point in the run the TDOF,I
!                                                            tables are printed out
      CHARACTER(43*BYTE)              :: MODNAM            ! Name to write to screen to describe module being run





! **********************************************************************************************************************************
! Part 1:  Generate TSET table
! ------
      CALL OURTIM
      MODNAM = ' DOF Set Table                              '
      WRITE(SC1,1093) MODNAM, HOUR, MINUTE, SEC, SFRAC
      CALL TSET_PROC

! Part 2:  Generate USET table if there are any USETSTR entries entries in the Bulk Data.
! ------
      IF (NUM_USETSTR > 0) THEN
         CALL USET_PROC
      ENDIF

! Part 3: Generate TDOF table from TSET
! ------
      CALL OURTIM
      MODNAM = ' DOF Number Table                           '
      WRITE(SC1,1093) MODNAM, HOUR, MINUTE, SEC, SFRAC
      CALL TDOF_PROC ( TDOF_MSG )

! Part 4: Make sure that NDOFSE /= 0 only in statics
! ------
      IF (NDOFSE > 0) THEN
         IF ((SOL_NAME(1:7) /= 'STATICS') .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
            WRITE(ERR,1323) SOL_NAME
            WRITE(F06,1323) SOL_NAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1093 FORMAT(5X,A,18X,2X,I2,':',I2,':',I2,'.',I3)

 1323 FORMAT(' *ERROR  1323: ENFORCED DISPLACEMENTS ONLY ALLOWED IN STATICS SOLUTION. HOWEVER, SOL = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DOF_PROC


      SUBROUTINE TDOF_PROC ( TDOF_MSG )

! TDOF table generation. TDOF is a table that has the DOF number for every DOF. The table has NDOFG rows and MTDOF columns where:

!   NDOFG = the total number of degrees of freedom in the model (6 times the number of grids plus the number of SPOINT's)
!   MTDOF = the number of different displ sets (currently 14 plus 2 USET user's sets) plus 4 cols for the internal grid/component
!           and external grid/component

! The TDOF table is sorted in grid numerical order. It's sister table, TDOFI, is TDOF sorted in internal G-set order.
! Cols 5 thru MTDOF of these tables give the degree of freedom (DOF) number for each of the 16 different displ sets/user's sets.
! An example of a TDOF table written out to an F06 file is shown below:

!                                              DEGREE OF FREEDOM TABLE SORTED ON GRID NUMBER (TDOF)

! EXTERNAL  INTERNAL                                     DOF NUMBER FOR DISPLACEMENT SET:
! GRD-COMP  GRD-COMP ---------------------------------------------------------------------------------------------------------------
! NUMBER    NUMBER         G      M      N     SA     SB     SG     SZ     SE      S      F      O      A     R      L     U1     U2

!    1022-1       1-1      1      0      1      0      0      0      0      0      0      1      0      1     0      1      0      0
!        -2        -2      2      0      2      0      0      0      0      0      0      2      0      2     0      2      0      1
!        -3        -3      3      0      3      0      0      0      0      0      0      3      0      3     0      3      0      0
!        -4        -4      4      0      4      0      0      0      0      0      0      4      0      4     0      4      0      2
!        -5        -5      5      0      5      0      0      0      0      0      0      5      0      5     0      5      0      0
!        -6        -6      6      0      6      0      0      0      0      0      0      6      0      6     0      6      0      3

!    1043-1       2-1      7      0      7      0      1      0      1      0      1      0      0      0     0      0      1      0
!        -2        -2      8      0      8      0      2      0      2      0      2      0      0      0     0      0      0      0
!        -3        -3      9      0      9      0      0      0      0      0      0      7      0      7     0      7      2      0
!        -4        -4     10      0     10      0      0      0      0      0      0      8      0      8     0      8      0      0
!        -5        -5     11      0     11      0      0      0      0      0      0      9      0      9     0      9      3      0
!        -6        -6     12      0     12      0      0      0      0      0      0     10      0     10     0     10      0      0


!                     ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ----- ------ ------ ------
! TOTAL NUMBER OF DOF:    12      0     12      0      2      6      2      0      2     12      0     10     0     10      3      3


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LDOFG, MTDOF, NDOFA, NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, NDOFO,   &
                                         NDOFR, NDOFS, NDOFSA, NDOFSB, NDOFSE, NDOFSG, NDOFSZ, NGRID, NUM_USET_U1, NUM_USET_U2,    &
                                         SOL_NAME, WARN_ERR
      USE PARAMS, ONLY                :  EIGESTL, PRTDOF
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET, TDOF, TDOFI, TDOF_ROW_START, USET
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EIG_N2, GRID, GRID_ID, GRID_SEQ, INV_GRID_SEQ

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, CALC_TDOF_ROW_START, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SORTING, ONLY               :  SORT_TDOF
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_TDOF
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TDOF_PROC'
      CHARACTER(LEN=*), INTENT(IN)    :: TDOF_MSG          ! Message to be printed out regarding at what point in the run the TDOF,I
!                                                            tables are printed out
      CHARACTER(  5*BYTE)             :: SET_NAME          ! A data set name for output purposes

      INTEGER(LONG)                   ::  A_SET_COL        ! Col no. in array TDOF where the  A-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  F_SET_COL        ! Col no. in array TDOF where the  F-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  G_SET_COL        ! Col no. in array TDOF where the  G-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  L_SET_COL        ! Col no. in array TDOF where the  L-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  M_SET_COL        ! Col no. in array TDOF where the  M-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  N_SET_COL        ! Col no. in array TDOF where the  N-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  O_SET_COL        ! Col no. in array TDOF where the  O-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  R_SET_COL        ! Col no. in array TDOF where the  R-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   ::  S_SET_COL        ! Col no. in array TDOF where the  S-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: SA_SET_COL        ! Col no. in array TDOF where the SA-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: SB_SET_COL        ! Col no. in array TDOF where the SB-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: SE_SET_COL        ! Col no. in array TDOF where the SE-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: SG_SET_COL        ! Col no. in array TDOF where the SG-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: SZ_SET_COL        ! Col no. in array TDOF where the SZ-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: U1_SET_COL        ! Col no. in array TDOF where the U1-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: U2_SET_COL        ! Col no. in array TDOF where the U2-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I_USET_U1         ! Counter for USET U1
      INTEGER(LONG)                   :: I_USET_U2         ! Counter for USET U2
      INTEGER(LONG)                   :: IGRID             ! Internal grid number
      INTEGER(LONG)                   :: IROW              ! Row number in array TDOF or TDOFI
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)




! **********************************************************************************************************************************
      WRITE(SC1, * ) '     TDOF PROC'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

! Call routine to calc what row in TDOF each grid's data begins

      CALL CALC_TDOF_ROW_START ( 'Y' )

! First, set NDOFG = LDOFG. It will be counted later.

      NDOFG = LDOFG

! Get column numbers for array TDOF different displ sets

      CALL TDOF_COL_NUM ( 'G ',  G_SET_COL )
      CALL TDOF_COL_NUM ( 'N ',  N_SET_COL )
      CALL TDOF_COL_NUM ( 'M ',  M_SET_COL )
      CALL TDOF_COL_NUM ( 'F ',  F_SET_COL )
      CALL TDOF_COL_NUM ( 'S ',  S_SET_COL )
      CALL TDOF_COL_NUM ( 'SA', SA_SET_COL )
      CALL TDOF_COL_NUM ( 'SB', SB_SET_COL )
      CALL TDOF_COL_NUM ( 'SG', SG_SET_COL )
      CALL TDOF_COL_NUM ( 'SZ', SZ_SET_COL )
      CALL TDOF_COL_NUM ( 'SE', SE_SET_COL )
      CALL TDOF_COL_NUM ( 'O ',  O_SET_COL )
      CALL TDOF_COL_NUM ( 'A ',  A_SET_COL )
      CALL TDOF_COL_NUM ( 'R ',  R_SET_COL )
      CALL TDOF_COL_NUM ( 'L ',  L_SET_COL )
      CALL TDOF_COL_NUM ( 'U1', U1_SET_COL )
      CALL TDOF_COL_NUM ( 'U2', U2_SET_COL )

! Set 1st 4 cols of TDOF (actual grid ID - component number, internal grid ID - component number). The rows are in GRID_ID order,
! TDOF_ROW_START(I) + J - 1 for component J of point I (subr CALC_TDOF_ROW_START), so each point is written with its own number of
! components. (It took the number of components of INV_GRID_SEQ(I), the point at sequence position I: with scalar points and grids
! sequenced in another order than their IDs, the rows were labelled with the wrong point and component, and subr MGGS_MASS_MATRIX,
! which finds each point's first row by its component label 1, put the CMASS masses at other DOFs: ERROR 989 (MLL singular), or
! a wrong eigenvalue without a message.)

      IROW = 0
      CALL COUNTER_INIT('       Process col 1-4 of TDOF', NGRID)
      DO I = 1,NGRID
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
         DO J = 1,NUM_COMPS
            IROW = IROW + 1
            TDOF(IROW,1) = GRID_ID(I)
            TDOF(IROW,2) = J
            TDOF(IROW,3) = GRID_SEQ(I)
            TDOF(IROW,4) = J
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Calc TDOF for G-set (col 5 in TDOF). We can do this at this point since all components go in G-set.

      NDOFG = 0
      CALL COUNTER_INIT('       Process G -set         ', NGRID)

      DO I=1,NGRID
         IGRID = INV_GRID_SEQ(I)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            IROW = TDOF_ROW_START(IGRID) + J - 1
            NDOFG = NDOFG + 1
            IF (NDOFG > LDOFG) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, LDOFG, 'TDOF' )
            TDOF(IROW,G_SET_COL) = NDOFG
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Make sure NDOFG exactly equals LDOFG (otherwise coding error)

      IF (NDOFG /= LDOFG) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1311) SUBR_NAME, NDOFG, LDOFG
         WRITE(F06,1311) SUBR_NAME, NDOFG, LDOFG
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Put M-set DOF numbers into TDOF (at M_SET_COL)

      IF (NDOFM > 0) THEN
         NDOFM = 0
         CALL COUNTER_INIT('       Process M -set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'M ') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFM = NDOFM + 1
                  TDOF(IROW,M_SET_COL) = NDOFM
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put SA-set DOF numbers into TDOF (at SA_SET_COL)

      IF (NDOFSA > 0) THEN
         NDOFSA = 0
         CALL COUNTER_INIT('       Process SA-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'SA') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFSA = NDOFSA + 1
                  TDOF(IROW,SA_SET_COL) = NDOFSA
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put SB-set DOF numbers into TDOF (at SB_SET_COL)

      IF (NDOFSB > 0) THEN
         NDOFSB = 0
         CALL COUNTER_INIT('       Process SB-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'SB') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFSB = NDOFSB + 1
                  TDOF(IROW,SB_SET_COL) = NDOFSB
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put SG-set DOF numbers into TDOF (at SG_SET_COL)

      IF (NDOFSG > 0) THEN
         NDOFSG = 0
         CALL COUNTER_INIT('       Process SG-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'SG') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFSG = NDOFSG + 1
                  TDOF(IROW,SG_SET_COL) = NDOFSG
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put SE-set DOF numbers into TDOF (at SE_SET_COL)

      IF (NDOFSE > 0) THEN
         NDOFSE = 0
         CALL COUNTER_INIT('       Process SE-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'SE') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFSE = NDOFSE + 1
                  TDOF(IROW,SE_SET_COL) = NDOFSE
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put O-set DOF numbers into TDOF (at O_SET_COL)

      IF (NDOFO > 0) THEN
         NDOFO = 0
         CALL COUNTER_INIT('       Process O -set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'O ') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFO = NDOFO + 1
                  TDOF(IROW,O_SET_COL) = NDOFO
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Put R-set DOF numbers into TDOF (at R_SET_COL)

      IF (NDOFR > 0) THEN
         NDOFR = 0
         CALL COUNTER_INIT('       Process R -set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (TSET(IGRID,J) == 'R ') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  NDOFR = NDOFR + 1
                  TDOF(IROW,R_SET_COL) = NDOFR
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Calc TDOF for N-set based on G-set minus M-set = S-set + O-set + R-set + L-set

      NDOFN = 0
      CALL COUNTER_INIT('       Process N -set         ', NGRID)
      DO I=1,NGRID
         IGRID = INV_GRID_SEQ(I)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            IROW  = TDOF_ROW_START(IGRID) + J - 1
            IF ((TDOF(IROW,G_SET_COL) > 0) .AND. (TDOF(IROW,M_SET_COL) == 0)) THEN
               NDOFN = NDOFN + 1
               TDOF(IROW,N_SET_COL) = NDOFN
            ELSE
               TDOF(IROW,N_SET_COL) = 0
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Calc DOF'S in SZ-set (all zero SPC's) based on SA + SB + SG

      IF ((NDOFSA > 0) .OR. (NDOFSB > 0) .OR. (NDOFSG > 0)) THEN
         NDOFSZ = 0
         CALL COUNTER_INIT('       Process SZ-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IROW  = TDOF_ROW_START(IGRID) + J - 1
               IF ((TDOF(IROW,SA_SET_COL) > 0) .OR. (TDOF(IROW,SB_SET_COL) > 0) .OR. (TDOF(IROW,SG_SET_COL) > 0)) THEN
                  NDOFSZ = NDOFSZ + 1
                  TDOF(IROW,SZ_SET_COL) = NDOFSZ
               ELSE
                  TDOF(IROW,SZ_SET_COL) = 0
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Calc DOF'S in S-set based on SZ + SE

      IF ((NDOFSZ > 0) .OR. (NDOFSE > 0)) THEN
         NDOFS = 0
         CALL COUNTER_INIT('       Process S -set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IROW  = TDOF_ROW_START(IGRID) + J - 1
               IF ((TDOF(IROW,SZ_SET_COL) > 0) .OR. (TDOF(IROW,SE_SET_COL) > 0)) THEN
                  NDOFS = NDOFS + 1
                  TDOF(IROW,S_SET_COL) = NDOFS
               ELSE
                  TDOF(IROW,S_SET_COL) = 0
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Calc TDOF for F-set based on N-set minus S-set

      NDOFF = 0
      CALL COUNTER_INIT('       Process F -set         ', NGRID)
      DO I=1,NGRID
         IGRID = INV_GRID_SEQ(I)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            IROW  = TDOF_ROW_START(IGRID) + J - 1
            IF ((TDOF(IROW,N_SET_COL) > 0) .AND. (TDOF(IROW,S_SET_COL) == 0)) THEN
               NDOFF = NDOFF + 1
               TDOF(IROW,F_SET_COL) = NDOFF
            ELSE
               TDOF(IROW,F_SET_COL) = 0
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Calc TDOF for A-set based on F-set minus O-set

      NDOFA = 0
      CALL COUNTER_INIT('       Process A -set         ', NGRID)
      DO I=1,NGRID
         IGRID = INV_GRID_SEQ(I)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            IROW  = TDOF_ROW_START(IGRID) + J - 1
            IF ((TDOF(IROW,F_SET_COL) > 0) .AND. (TDOF(IROW,O_SET_COL) == 0)) THEN
               NDOFA = NDOFA + 1
               TDOF(IROW,A_SET_COL) = NDOFA
            ELSE
               TDOF(IROW,A_SET_COL) = 0
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Calc TDOF for L-set based on A-set minus R-set

      NDOFL = 0
      CALL COUNTER_INIT('       Process L -set         ', NGRID)
      DO I=1,NGRID
         IGRID = INV_GRID_SEQ(I)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            IROW  = TDOF_ROW_START(IGRID) + J - 1
            IF ((TDOF(IROW,A_SET_COL) > 0) .AND. (TDOF(IROW,R_SET_COL) == 0)) THEN
               NDOFL = NDOFL + 1
               TDOF(IROW,L_SET_COL) = NDOFL
            ELSE
               TDOF(IROW,L_SET_COL) = 0
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO

! Calc TDOF for USET U1-set and put DOF numbers into TDOF at U1_SET_COL

      IF (NUM_USET_U1 > 0) THEN
         I_USET_U1 = 0
         CALL COUNTER_INIT('       Process U1-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (USET(IGRID,J) == 'U1') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  I_USET_U1 = I_USET_U1 + 1
                  TDOF(IROW,U1_SET_COL) = I_USET_U1
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Calc TDOF for USET U2-set and put DOF numbers into TDOF at U2_SET_COL

      IF (NUM_USET_U2 > 0) THEN
         I_USET_U2 = 0
         CALL COUNTER_INIT('       Process U2-set         ', NGRID)
         DO I=1,NGRID
            IGRID = INV_GRID_SEQ(I)
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               IF (USET(IGRID,J) == 'U2') THEN
                  IROW = TDOF_ROW_START(IGRID) + J - 1
                  I_USET_U2 = I_USET_U2 + 1
                  TDOF(IROW,U2_SET_COL) = I_USET_U2
               ENDIF
            ENDDO
            CALL COUNTER_PROGRESS(I)
         ENDDO
      ENDIF

! Sort TDOF so that G-set DOF's are in numerical order


      WRITE(SC1,12345,ADVANCE='NO') '       Setting up to get TDOFI', CR13
      DO I=1,NDOFG
         DO J=1,MTDOF
            TDOFI(I,J) = TDOF(I,J)
         ENDDO
      ENDDO

      WRITE(SC1,12345,ADVANCE='NO') '       Sort TDOF to get TDOFI ', CR13
      CALL SORT_TDOF ( SUBR_NAME, 'TDOF', NDOFG, TDOFI, G_SET_COL )

! Table TDOF is printed in the F06 file if B.D. PARAM PRTDOF = 1 or 3

      IF (PRTDOF > 0) THEN
         CALL WRITE_TDOF ( TDOF_MSG )
      ENDIF

! Make sure that an R-set has been defined if this is CB soln

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         IF (NDOFR < 6) THEN
            SET_NAME = 'R-SET'
            WRITE(ERR,1312) SET_NAME, NDOFR
            WRITE(F06,1312) SET_NAME, NDOFR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! Make sure that if the R set exists it has at least 6 DOF's defined

      IF (NDOFR > 0) THEN
         IF (NDOFR < 6) THEN
            SET_NAME = 'R-SET'
            WRITE(ERR,1304) SET_NAME, NDOFR
            WRITE(F06,1304) SET_NAME, NDOFR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! If this is an eigen analysis (or CB) make sure that some spec for the number of modes to be found exists

      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         IF ((NDOFL > EIGESTL) .AND. (EIG_N2 == 0)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1313) NDOFL, EIGESTL
            WRITE(F06,1313) NDOFL, EIGESTL
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,*) CR13



      RETURN

! **********************************************************************************************************************************
 1304 FORMAT(' *ERROR  1304: AN ',A,' HAS BEEN DEFINED VIA SUPORT BULK DATA ENTRIES AND THEREFORE MUST HAVE AT LEAST 6 DOF IN IT.' &
                        ,/,14X,' HOWEVER, THERE WERE ONLY ',I2,' DOF DEFINED IN THIS SET')

 1308 FORMAT('  Array TDOF_ROW_START gives the row in array TDOF where the DOF data begins for the GRID listed below',//,          &
             '                  I     GRID   TDOF_ROW_START',/)

 1311 FORMAT(' *ERROR  1311: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NDOFG = ',I12,' CALCULATED IN THIS SUBR MUST AGREE WITH LDOFG = ',I12,'  DETERMINED EARLIER')

 1312 FORMAT(' *ERROR  1312: FOR SOL = ''GEN CB MODEL'' THERE MUST BE AN ',A,'-SET WITH AT LEAST NDOFR = 6 DOF''s.'                &
                    ,/,14X,' HOWEVER ONLY ',I1,' DOF''s WERE DEFINED ON BULK DATA SUPORT ENTRIES')

 1313 format(' *ERROR  1313: FOR SOL = "MODES" OR "GEN CB MODEL" THE EIGRL ENTRY MUST HAVE THE NUMBER OF DESIRED MODES > 0 OR THE',&
                           ' PROBLEM DOF SIZE'                                                                                     &
                    ,/,14X,' (NDOFL = ',I8,') MUST BE LESS THAT PARAM EIGESTL = ',I8,' (OR USE LARGER VALUE FOR PARAM EIGESTL)')

12345 FORMAT(A, A)

! **********************************************************************************************************************************

      END SUBROUTINE TDOF_PROC
      SUBROUTINE USET_PROC

! USET is a table that specifies which grid/component pairs are ones defined by the user on Bulk Data USET or USET1 entries.
! The table has NGRID rows and 6 columns (1 col for each of the 6 components of displ at a grid). The table can have entries that
! are either 'U1' or 'U2' (i.e. user set U1 or U2). The table is constructed like tha TSET table (see explanation in subr TSET_PROC)
! An example of a USET table written to the F06 file is shown below:


!                        DEGREES OF FREEDOM DEFINED ON USET BULK DATA ENTRIES
!                        ----------------------------------------------------

!                     Grid       T1       T2       T3       R1       R2       R3

!                     1012       U1       U2       U1       --       U1       --
!                     1022       --       U2       --       U2       --       U2
!                     1031       U2       U1       U1       U2       --       --
!                     1032       U1       --       U1       --       U1       --
!                     1033       --       U1       --       U1       --       --
!                     1041       U1       --       U1       U2       U1       U2


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, L1X, L1X_MSG, LINK1X
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ENFORCED, FATAL_ERR, NGRID, NUM_USET_RECORDS, NUM_USET_U1, NUM_USET_U2
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN, USET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID

      USE DOF_TABLE_LIFECYCLE, ONLY   :  ALLOCATE_DOF_TABLES
      USE FILE_LIFECYCLE, ONLY        :  FILE_OPEN, READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_SET_CONSTRUCTION, ONLY: RDOF
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_USET

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'USET_PROC'
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: SNAME             ! The name of a set ('U1' or 'U2')

      INTEGER(LONG)                   :: USET_ERR   = 0    ! Count of errors that result from setting displ sets in USET
      INTEGER(LONG)                   :: GRID1             ! An actual grid ID
      INTEGER(LONG)                   :: GRID2             ! An actual grid ID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number, in array GRID_ID, where an actual grid ID resides
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of errors that result from undefined grid ID's
      INTEGER(LONG)                   :: I,J,GRID_NUM      ! DO loop indices
      INTEGER(LONG)                   :: ICOMP             ! DOF components read from file LINK1X (SPC's) or LINK1N (ASET/OMIT's)
      INTEGER(LONG)                   :: IERRT             ! Total number of errors found here
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: REC_NO    = 0     ! Record number when reading a file




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! ----------------------------------------------------------------------------------------------------------------------------------
! Initialize

      GID_ERR     = 0
      IERRT       = 0
      USET_ERR    = 0
      NUM_USET_U1 = 0
      NUM_USET_U2 = 0
      REC_NO      = 0

! Allocate memory to USET and initialize

      CALL ALLOCATE_DOF_TABLES ( 'USET', SUBR_NAME )
      DO I=1,NGRID
         DO J=1,6
            USET(I,J) = '--'
         ENDDO
      ENDDO

! Process USET data from file L1X (data written when USET and USET1 Bulk Data entries were read)

      CALL FILE_OPEN ( L1X, LINK1X, OUNT, 'OLD', L1X_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

i_do6:DO I=1,NUM_USET_RECORDS

         READ(L1X,IOSTAT=IOCHK) SNAME, ICOMP, GRID1, GRID2
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1X, L1X_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Error reading SPC file . No sense continuing
         ENDIF

         IF ((SNAME /= 'U1') .AND. (SNAME /= 'U2')) THEN   ! Make sure that SNAME = 'U1' or 'U2'
            WRITE(ERR,1369) SUBR_NAME, LINK1X, SNAME
            WRITE(F06,1369) SUBR_NAME, LINK1X, SNAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Pgm error (data in file LINK1X must be for sets 'U1' or 'U2')
         ENDIF

! No error, so processes data.

         CALL RDOF ( ICOMP, CDOF )                   ! Convert ICOMP to CDOF char form for use below

         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID1, GRID_ID_ROW_NUM )
         IF (GRID_ID_ROW_NUM == -1) THEN
            GID_ERR = GID_ERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1822) 'GRID ', GRID1, 'USET OR USET1', SNAME
            WRITE(F06,1822) 'GRID ', GRID1, 'USET OR USET1', SNAME
         ENDIF
         IF (GRID2 /= GRID1) THEN
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID2, GRID_ID_ROW_NUM )
            IF (GRID_ID_ROW_NUM == -1) THEN
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'GRID ', GRID1, 'USET OR USET1', SNAME
               WRITE(F06,1822) 'GRID ', GRID1, 'USET OR USET1', SNAME
            ENDIF
         ENDIF

         IF (GID_ERR == 0) THEN                      ! Put CDOF data in USET for GRID1 thru GRID2
            DO GRID_NUM=GRID1,GRID2                  ! GRID2 >= GRID1 was checked in subr BD_SPC1
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM /= -1) THEN
                  CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM, NUM_COMPS, SUBR_NAME )
                  DO J = 1,NUM_COMPS                 ! Put data in USET and write enforced displ to L1H.
                     IF (CDOF(J) == '1') THEN
                        IF      (SNAME == 'U1') THEN
                           IF ((USET(GRID_ID_ROW_NUM,J) == '--') .OR. (USET(GRID_ID_ROW_NUM,J) == 'U1')) THEN
                              IF (USET(GRID_ID_ROW_NUM,J) == '--') THEN
                                 NUM_USET_U1 = NUM_USET_U1 + 1
                              ENDIF
                              USET(GRID_ID_ROW_NUM,J) = 'U1'
                           ELSE
                              USET_ERR = USET_ERR + 1
                              FATAL_ERR = FATAL_ERR + 1
                              WRITE(ERR,1332) SNAME, GRID_NUM, J, SNAME, USET(GRID_ID_ROW_NUM,J)
                              WRITE(F06,1332) SNAME, GRID_NUM, J, SNAME, USET(GRID_ID_ROW_NUM,J)
                           ENDIF
                        ELSE IF (SNAME == 'U2') THEN
                           IF ((USET(GRID_ID_ROW_NUM,J) == '--') .OR. (USET(GRID_ID_ROW_NUM,J) == 'U2')) THEN
                              IF (USET(GRID_ID_ROW_NUM,J) == '--') THEN
                                 NUM_USET_U2 = NUM_USET_U2 + 1
                              ENDIF
                              USET(GRID_ID_ROW_NUM,J) = 'U2'
                           ELSE
                              USET_ERR = USET_ERR + 1
                              FATAL_ERR = FATAL_ERR + 1
                              WRITE(ERR,1332) SNAME, GRID_NUM, J, SNAME, USET(GRID_ID_ROW_NUM,J)
                              WRITE(F06,1332) SNAME, GRID_NUM, J, SNAME, USET(GRID_ID_ROW_NUM,J)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

      ENDDO i_do6

      IERRT = GID_ERR + USET_ERR
      IF (IERRT > 0) THEN
         WRITE(ERR,1322) IERRT, SUBR_NAME
         WRITE(F06,1322) IERRT, SUBR_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL WRITE_USET



      RETURN

! **********************************************************************************************************************************
 1322 FORMAT(' *ERROR  1322: PROCESSING STOPPED DUE TO THE ABOVE LISTED ',I8,' ERRORS IN SUBR ',A)

 1331 FORMAT(' *ERROR  1331: GRID POINT ',I8,' HAS COMPONENT ',I2,' IN THE ',A2,' DISPL SET. (PERM SPC ON GRID CARD)',             &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' DISPL SET')

 1332 FORMAT(' *ERROR  1332: USET SET ',A,' HAS GRID POINT ',I8,' COMPONENT ',I2,' IN THE ',A2,' USET.',                           &
                           ' HOWEVER THIS GRID/COMPONENT IS ALREADY IN THE ',A2,' USET')

 1369 FORMAT(' *ERROR  1369: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' RECORD READ FROM FILE: ',A                                                                            &
                    ,/,14X,' INDICATES THAT AN SPCd DOF BELONGS TO THE "',A2,'" SET. MUST BE EITHER "SE" OR "SB"')

 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,2X,A,' IS UNDEFINED')

! **********************************************************************************************************************************

      END SUBROUTINE USET_PROC


   END MODULE DOF_NUMBERING
