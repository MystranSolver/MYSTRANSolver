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

   MODULE TEMP_FILE_WRITERS

   USE FILE_LIFECYCLE, ONLY :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN
   USE DOF_ARRAY_INDEXING, ONLY :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
   USE DIAGNOSTICS_MEMORY_REPORTING, ONLY :  LINK_MESSAGE
   USE DOF_LOOKUP_UTILS, ONLY :  TDOF_COL_NUM

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_L1M, WRITE_L1Z, WRITE_EDAT, WRITE_ELM_OT4, WRITE_FIJFIL, WRITE_GRD_OT4, WRITE_PARTND_MAT_HDRS, WRITE_DOF_TABLES, WRITE_TDOF, WRITE_TSET, WRITE_USET, WRITE_USETSTR, WRITE_USERIN_BD_CARDS

   CONTAINS

      SUBROUTINE WRITE_L1M

! Writes data from file LINK1M of eigenvalue extraction data (actual eigenvalues/vectors are not in this file)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE SCONTR, ONLY                :  LINKNO, NUM_EIGENS
      USE IOUNT1, ONLY                :  ERR, F06, L1M, L1M_MSG, L1MSTAT, LINK1M, SC1, WRT_ERR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL, GEN_MASS, MODE_NUM

      USE MODEL_STUF, ONLY            :  EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_METH, EIG_MSGLVL, EIG_LAP_MAT_TYPE, &
                                         EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, EIG_SID, EIG_SIGMA, EIG_VECS, MAXMIJ,    &
                                         MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT

      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  LINK_MESSAGE

      IMPLICIT NONE

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to

! **********************************************************************************************************************************
! Make units for writing errors the screen and output file

      OUNT(1) = ERR
      OUNT(2) = F06

!xx   STATUS = 'OLD    '
!xx   RW     = 'WRITE'
      CALL FILE_OPEN ( L1M, LINK1M, OUNT, 'REPLACE', L1M_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

      CALL LINK_MESSAGE('WRITE EIGENVALUE DATA FROM PRIOR LINK')

      WRITE(L1M) EIG_SID
      WRITE(L1M) EIG_METH
      WRITE(L1M) EIG_FRQ1
      WRITE(L1M) EIG_FRQ2
      WRITE(L1M) EIG_N1
      WRITE(L1M) EIG_N2
      WRITE(L1M) EIG_VECS
      WRITE(L1M) EIG_CRIT
      WRITE(L1M) EIG_NORM
      WRITE(L1M) EIG_GRID
      WRITE(L1M) EIG_COMP
      WRITE(L1M) EIG_MODE
      WRITE(L1M) EIG_SIGMA
      WRITE(L1M) EIG_LAP_MAT_TYPE
      WRITE(L1M) EIG_MSGLVL
      WRITE(L1M) EIG_NCVFACL
      WRITE(L1M) NUM_FAIL_CRIT
      WRITE(L1M) MAXMIJ
      WRITE(L1M) MIJ_ROW
      WRITE(L1M) MIJ_COL

      DO I=1,NUM_EIGENS
         WRITE(L1M) MODE_NUM(I), EIGEN_VAL(I), GEN_MASS(I)
      ENDDO

      CALL FILE_CLOSE ( L1M, LINK1M, 'KEEP' )


! **********************************************************************************************************************************

      END SUBROUTINE WRITE_L1M


      SUBROUTINE WRITE_L1Z

! Writes file LINK1Z of some of the data needed in a restart if user has CHKPNT in the Exec Control

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06, L1Z, LINK1Z, L1Z_MSG, L1ZSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSUB, SOL_NAME
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, MPCSET, SPCSET, SUBLOD

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_L1Z'

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN




! **********************************************************************************************************************************
      CALL FILE_OPEN ( L1Z, LINK1Z, OUNT, 'REPLACE', L1Z_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      WRITE(L1Z) SOL_NAME
      WRITE(L1Z) NSUB
      WRITE(L1Z) MPCSET
      WRITE(L1Z) SPCSET
      DO I=1,NSUB
         WRITE(L1Z) SUBLOD(I,1), SUBLOD(I,2)
      ENDDO
      WRITE(L1Z) CC_EIGR_SID
      CALL FILE_CLOSE ( L1Z, LINK1Z, L1ZSTAT )

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_L1Z


      SUBROUTINE WRITE_EDAT

! This subr write array EDAT based on user input Bulk Data DEBUG (see module DEBUG_PARAMS). EDAT is an array that contains all of
! the information read fron element connection entries in the Bulk Data

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06

      USE SCONTR, ONLY                :  BLNK_SUB_NAM  , LGUSERIN      , LSUSERIN      , NELE          , NCUSERIN      , WARN_ERR, &
                                         MEDAT_CBAR    , MEDAT_CBEAM   , MEDAT_CBUSH   , MEDAT_CELAS1  , MEDAT_CELAS2  ,           &
                                         MEDAT_CELAS3  , MEDAT_CELAS4  , MEDAT_CHEXA8  , MEDAT_CHEXA20 , MEDAT_CPENTA6 ,           &
                                         MEDAT_CPENTA15, MEDAT_PLOTEL  , MEDAT_CQUAD   , MEDAT_CROD    ,                           &
                                         MEDAT_CSHEAR  , MEDAT_CTETRA4 , MEDAT_CTETRA10, MEDAT_CTRIA   ,                           &
                                         MEDAT_CUSER1  , MEDAT0_CUSERIN, METYPE

      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, ETYPE
      USE PARAMS, ONLY                :  SUPWARN

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_EDAT'
      CHARACTER(10*BYTE)              :: ICHAR             !

                                                           ! Descriptors for items in EDAT for each element in the Bulk Data deck
      CHARACTER(16*BYTE)              :: EDAT_DESCR(0:MAX(2*LGUSERIN+LSUSERIN+6,25),METYPE)
      CHARACTER(56*BYTE)              :: EDAT_EXPLAIN_CID  = ' (-99 is a placeholder to indicate that field was blank)'
      CHARACTER(61*BYTE)              :: EDAT_EXPLAIN_VVEC = ' (neg number -i indicates v-vector is in row i in array VVEC)'

      INTEGER(LONG)                   :: EPNTK             ! Value in array EPNT where data begins for an element
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1,L              ! Counters
      INTEGER(LONG)                   :: MEDAT             ! Number of terms in EDAT for a specific element type
      INTEGER(LONG)                   :: NG                ! Number of grids defined on a CUSERIN entry
      INTEGER(LONG)                   :: NS                ! Number of scalar points defined on a CUSERIN entry




! **********************************************************************************************************************************
! Initialize

      DO I=0,MAX(2*LGUSERIN+LSUSERIN+6,25)
         DO j=1,METYPE
            EDAT_DESCR(I,J)(1:) = ' '
         ENDDO
      ENDDO

! Set values for EDAT_DESCR array

      EDAT_DESCR( 0, 1) = 'BAR             '
      EDAT_DESCR( 1, 1) = 'Elem ID         '
      EDAT_DESCR( 2, 1) = 'Prop ID         '
      EDAT_DESCR( 3, 1) = 'Grid A          '
      EDAT_DESCR( 4, 1) = 'Grid B          '
      EDAT_DESCR( 5, 1) = 'V-vector key    '
      EDAT_DESCR( 6, 1) = 'Pin Flag A      '
      EDAT_DESCR( 7, 1) = 'Pin Flag B      '
      EDAT_DESCR( 8, 1) = 'Offset key      '

      EDAT_DESCR( 0, 2) = 'BEAM            '
      EDAT_DESCR( 1, 2) = 'Elem ID         '
      EDAT_DESCR( 2, 2) = 'Prop ID         '
      EDAT_DESCR( 3, 2) = 'Grid A          '
      EDAT_DESCR( 4, 2) = 'Grid B          '
      EDAT_DESCR( 5, 2) = 'V-vector key    '
      EDAT_DESCR( 6, 2) = 'Pin Flag A      '
      EDAT_DESCR( 7, 2) = 'Pin Flag B      '
      EDAT_DESCR( 8, 2) = 'Offset key      '

      EDAT_DESCR( 0, 3) = 'BUSH            '
      EDAT_DESCR( 1, 3) = 'Elem ID         '
      EDAT_DESCR( 2, 3) = 'Prop ID         '
      EDAT_DESCR( 3, 3) = 'Number of grids '               ! If CBUSH G2 is blank BUSH has only 1 grid. Otherwise 2 grids
      EDAT_DESCR( 4, 3) = 'Grid A          '
      EDAT_DESCR( 5, 3) = 'Grid B          '
      EDAT_DESCR( 6, 3) = 'V-vector key    '
      EDAT_DESCR( 7, 3) = 'CID             '
      EDAT_DESCR( 8, 3) = 'OCID            '
      EDAT_DESCR( 9, 3) = 'Offset key      '

      EDAT_DESCR( 0, 4) = 'ELAS1           '
      EDAT_DESCR( 1, 4) = 'Elem ID         '
      EDAT_DESCR( 2, 4) = 'Prop ID         '
      EDAT_DESCR( 3, 4) = 'Grid A          '
      EDAT_DESCR( 4, 4) = 'Comps at Grid A '
      EDAT_DESCR( 5, 4) = 'Grid B          '
      EDAT_DESCR( 6, 4) = 'Comps at Grid B '

      EDAT_DESCR( 0, 5) = 'ELAS2           '
      EDAT_DESCR( 1, 5) = 'Elem ID         '
      EDAT_DESCR( 2, 5) = 'Prop ID         '
      EDAT_DESCR( 3, 5) = 'Grid A          '
      EDAT_DESCR( 4, 5) = 'Comps at Grid A '
      EDAT_DESCR( 5, 5) = 'Grid B          '
      EDAT_DESCR( 6, 5) = 'Comps at Grid B '

      EDAT_DESCR( 0, 6) = 'ELAS3           '
      EDAT_DESCR( 1, 6) = 'Elem ID         '
      EDAT_DESCR( 2, 6) = 'Prop ID         '
      EDAT_DESCR( 3, 6) = 'Scalar point A  '
      EDAT_DESCR( 4, 6) = 'Scalar point B  '

      EDAT_DESCR( 0, 7) = 'ELAS4           '
      EDAT_DESCR( 1, 7) = 'Elem ID         '
      EDAT_DESCR( 2, 7) = 'Prop ID         '
      EDAT_DESCR( 3, 7) = 'Scalar point A  '
      EDAT_DESCR( 4, 7) = 'Scalar point B  '

      EDAT_DESCR( 0, 8) = 'HEXA8           '
      EDAT_DESCR( 1, 8) = 'Elem ID         '
      EDAT_DESCR( 2, 8) = 'Prop ID         '
      EDAT_DESCR( 3, 8) = 'Grid 1          '
      EDAT_DESCR( 4, 8) = 'Grid 2          '
      EDAT_DESCR( 5, 8) = 'Grid 3          '
      EDAT_DESCR( 6, 8) = 'Grid 4          '
      EDAT_DESCR( 7, 8) = 'Grid 5          '
      EDAT_DESCR( 8, 8) = 'Grid 6          '
      EDAT_DESCR( 9, 8) = 'Grid 7          '
      EDAT_DESCR(10,  8) = 'Grid 8          '

      EDAT_DESCR( 0, 9) = 'HEXA20          '
      EDAT_DESCR( 1, 9) = 'Elem ID         '
      EDAT_DESCR( 2, 9) = 'Prop ID         '
      EDAT_DESCR( 3, 9) = 'Grid  1         '
      EDAT_DESCR( 4, 9) = 'Grid  2         '
      EDAT_DESCR( 5, 9) = 'Grid  3         '
      EDAT_DESCR( 6, 9) = 'Grid  4         '
      EDAT_DESCR( 7, 9) = 'Grid  5         '
      EDAT_DESCR( 8, 9) = 'Grid  6         '
      EDAT_DESCR( 9, 9) = 'Grid  7         '
      EDAT_DESCR(10, 9) = 'Grid  8         '
      EDAT_DESCR(11, 9) = 'Grid  9         '
      EDAT_DESCR(12, 9) = 'Grid 10         '
      EDAT_DESCR(13, 9) = 'Grid 11         '
      EDAT_DESCR(14, 9) = 'Grid 12         '
      EDAT_DESCR(15, 9) = 'Grid 13         '
      EDAT_DESCR(16, 9) = 'Grid 14         '
      EDAT_DESCR(17, 9) = 'Grid 15         '
      EDAT_DESCR(18, 9) = 'Grid 16         '
      EDAT_DESCR(19, 9) = 'Grid 17         '
      EDAT_DESCR(20, 9) = 'Grid 18         '
      EDAT_DESCR(21, 9) = 'Grid 19         '
      EDAT_DESCR(22, 9) = 'Grid 20         '

      EDAT_DESCR( 0,10) = 'PENTA6          '
      EDAT_DESCR( 1,10) = 'Elem ID         '
      EDAT_DESCR( 2,10) = 'Prop ID         '
      EDAT_DESCR( 3,10) = 'Grid 1          '
      EDAT_DESCR( 5,10) = 'Grid 2          '
      EDAT_DESCR( 5,10) = 'Grid 3          '
      EDAT_DESCR( 6,10) = 'Grid 4          '
      EDAT_DESCR( 7,10) = 'Grid 5          '
      EDAT_DESCR( 8,10) = 'Grid 6          '

      EDAT_DESCR( 0,11) = 'PENTA15         '
      EDAT_DESCR( 1,11) = 'Elem ID         '
      EDAT_DESCR( 2,11) = 'Prop ID         '
      EDAT_DESCR( 3,11) = 'Grid  1         '
      EDAT_DESCR( 4,11) = 'Grid  2         '
      EDAT_DESCR( 5,11) = 'Grid  3         '
      EDAT_DESCR( 6,11) = 'Grid  4         '
      EDAT_DESCR( 7,11) = 'Grid  5         '
      EDAT_DESCR( 8,11) = 'Grid  6         '
      EDAT_DESCR( 9,11) = 'Grid  7         '
      EDAT_DESCR(10,11) = 'Grid  8         '
      EDAT_DESCR(11,11) = 'Grid  9         '
      EDAT_DESCR(12,11) = 'Grid 10         '
      EDAT_DESCR(13,11) = 'Grid 11         '
      EDAT_DESCR(14,11) = 'Grid 12         '
      EDAT_DESCR(15,11) = 'Grid 13         '
      EDAT_DESCR(16,11) = 'Grid 14         '
      EDAT_DESCR(17,11) = 'Grid 15         '

      EDAT_DESCR( 0,12) = 'PLOTEL          '
      EDAT_DESCR( 1,12) = 'Elem ID         '
      EDAT_DESCR( 2,12) = 'Elem ID         '
      EDAT_DESCR( 3,12) = 'Grid  1         '
      EDAT_DESCR( 4,12) = 'Grid  2         '

      EDAT_DESCR( 0,13) = 'QUAD4           '
      EDAT_DESCR( 1,13) = 'Elem ID         '
      EDAT_DESCR( 2,13) = 'Prop ID         '
      EDAT_DESCR( 3,13) = 'Grid A          '
      EDAT_DESCR( 4,13) = 'Grid B          '
      EDAT_DESCR( 5,13) = 'Grid C          '
      EDAT_DESCR( 6,13) = 'Grid D          '
      EDAT_DESCR( 7,13) = 'Mtrl angle key  '
      EDAT_DESCR( 8,13) = 'Basic CID or not'
      EDAT_DESCR( 9,13) = 'Plate offset key'
      EDAT_DESCR(10,13) = 'PSHELL or PCOMP '
      EDAT_DESCR(11,13) = 'Plate thick key '

      EDAT_DESCR( 0,14) = 'QUAD4K          '
      EDAT_DESCR( 1,14) = 'Elem ID         '
      EDAT_DESCR( 2,14) = 'Prop ID         '
      EDAT_DESCR( 3,14) = 'Grid A          '
      EDAT_DESCR( 4,14) = 'Grid B          '
      EDAT_DESCR( 5,14) = 'Grid C          '
      EDAT_DESCR( 6,14) = 'Grid D          '
      EDAT_DESCR( 7,14) = 'Mtrl angle key  '
      EDAT_DESCR( 8,14) = 'Basic CID or not'
      EDAT_DESCR( 9,14) = 'Plate offset key'
      EDAT_DESCR(10,14) = 'PSHELL or PCOMP '
      EDAT_DESCR(11,14) = 'Plate thick key '

      EDAT_DESCR( 0,15) = 'ROD             '
      EDAT_DESCR( 1,15) = 'Elem ID         '
      EDAT_DESCR( 2,15) = 'Prop ID         '
      EDAT_DESCR( 3,15) = 'Grid A          '
      EDAT_DESCR( 4,15) = 'Grid B          '

      EDAT_DESCR( 0,16) = 'SHEAR           '
      EDAT_DESCR( 1,16) = 'Elem ID         '
      EDAT_DESCR( 2,16) = 'Prop ID         '
      EDAT_DESCR( 3,16) = 'Grid A          '
      EDAT_DESCR( 4,16) = 'Grid B          '
      EDAT_DESCR( 5,16) = 'Grid C          '
      EDAT_DESCR( 6,16) = 'Grid D          '

      EDAT_DESCR( 0,17) = 'TETRA4          '
      EDAT_DESCR( 1,17) = 'Elem ID         '
      EDAT_DESCR( 2,17) = 'Prop ID         '
      EDAT_DESCR( 3,17) = 'Grid 1          '
      EDAT_DESCR( 4,17) = 'Grid 2          '
      EDAT_DESCR( 5,17) = 'Grid 3          '
      EDAT_DESCR( 6,17) = 'Grid 4          '

      EDAT_DESCR( 0,18) = 'TETRA10         '
      EDAT_DESCR( 1,18) = 'Elem ID         '
      EDAT_DESCR( 2,18) = 'Prop ID         '
      EDAT_DESCR( 3,18) = 'Grid  1         '
      EDAT_DESCR( 4,18) = 'Grid  2         '
      EDAT_DESCR( 5,18) = 'Grid  3         '
      EDAT_DESCR( 6,18) = 'Grid  4         '
      EDAT_DESCR( 7,18) = 'Grid  5         '
      EDAT_DESCR( 8,18) = 'Grid  6         '
      EDAT_DESCR( 9,18) = 'Grid  7         '
      EDAT_DESCR(10,18) = 'Grid  8         '
      EDAT_DESCR(11,18) = 'Grid  9         '
      EDAT_DESCR(12,18) = 'Grid 10         '

      EDAT_DESCR( 0,19) = 'TRIA3           '
      EDAT_DESCR( 1,19) = 'Elem ID         '
      EDAT_DESCR( 2,19) = 'Prop ID         '
      EDAT_DESCR( 3,19) = 'Grid A          '
      EDAT_DESCR( 4,19) = 'Grid B          '
      EDAT_DESCR( 5,19) = 'Grid C          '
      EDAT_DESCR( 6,19) = 'Mtrl angle key  '
      EDAT_DESCR( 7,19) = 'Basic CID or not'
      EDAT_DESCR( 8,19) = 'Plate offset key'
      EDAT_DESCR( 9,19) = 'PSHELL or PCOMP '
      EDAT_DESCR(10,19) = 'Plate thick key '

      EDAT_DESCR( 0,20) = 'TRIA3K          '
      EDAT_DESCR( 1,20) = 'Elem ID         '
      EDAT_DESCR( 2,20) = 'Prop ID         '
      EDAT_DESCR( 3,20) = 'Grid A          '
      EDAT_DESCR( 4,20) = 'Grid B          '
      EDAT_DESCR( 5,20) = 'Grid C          '
      EDAT_DESCR( 6,20) = 'Mtrl angle key  '
      EDAT_DESCR( 7,20) = 'Basic CID or not'
      EDAT_DESCR( 8,20) = 'Plate offset key'
      EDAT_DESCR( 9,20) = 'PSHELL or PCOMP '
      EDAT_DESCR(10,20) = 'Plate thick key '

      EDAT_DESCR( 0,21) = 'USER1           '
      EDAT_DESCR( 1,21) = 'Elem ID         '
      EDAT_DESCR( 2,21) = 'Prop ID         '
      EDAT_DESCR( 3,21) = 'Grid A          '
      EDAT_DESCR( 4,21) = 'Grid B          '
      EDAT_DESCR( 5,21) = 'Grid C          '
      EDAT_DESCR( 6,21) = 'Grid D          '
      EDAT_DESCR( 7,21) = 'Grid V          '
      EDAT_DESCR( 8,21) = 'Pin Flag A      '
      EDAT_DESCR( 9,21) = 'Pin Flag B      '
      EDAT_DESCR(10,21) = 'Pin Flag C      '
      EDAT_DESCR(11,21) = 'Pin Flag D      '

      EDAT_DESCR( 0,22) = 'USERIN          '
      EDAT_DESCR( 1,22) = 'Elem ID         '
      EDAT_DESCR( 2,22) = 'Prop ID         '
      EDAT_DESCR( 3,22) = 'NG              '
      EDAT_DESCR( 4,22) = 'NS              '
      EDAT_DESCR( 5,22) = 'CID0            '

! ----------------------------------------------------------------------------------------------------------------------------------
! Write EDAT table

      WRITE(F06,2000)
      WRITE(F06,2001)
      WRITE(F06,2000)
      WRITE(F06,*)
      WRITE(F06,2002)

      I1 = 1
do_1: DO K=1,NELE

         IF      (ETYPE(K) == 'BAR     ') THEN   ;   MEDAT = MEDAT_CBAR
         ELSE IF (ETYPE(K) == 'BEAM    ') THEN   ;   MEDAT = MEDAT_CBEAM
         ELSE IF (ETYPE(K) == 'BUSH    ') THEN   ;   MEDAT = MEDAT_CBUSH
         ELSE IF (ETYPE(K) == 'ELAS1   ') THEN   ;   MEDAT = MEDAT_CELAS1
         ELSE IF (ETYPE(K) == 'ELAS2   ') THEN   ;   MEDAT = MEDAT_CELAS2
         ELSE IF (ETYPE(K) == 'ELAS3   ') THEN   ;   MEDAT = MEDAT_CELAS3
         ELSE IF (ETYPE(K) == 'ELAS4   ') THEN   ;   MEDAT = MEDAT_CELAS4
         ELSE IF (ETYPE(K) == 'HEXA8   ') THEN   ;   MEDAT = MEDAT_CHEXA8
         ELSE IF (ETYPE(K) == 'HEXA20  ') THEN   ;   MEDAT = MEDAT_CHEXA20
         ELSE IF (ETYPE(K) == 'PENTA6  ') THEN   ;   MEDAT = MEDAT_CPENTA6
         ELSE IF (ETYPE(K) == 'PENTA15 ') THEN   ;   MEDAT = MEDAT_CPENTA15
         ELSE IF (ETYPE(K) == 'PLOTEL  ') THEN   ;   MEDAT = MEDAT_PLOTEL
         ELSE IF (ETYPE(K) == 'QUAD4   ') THEN   ;   MEDAT = MEDAT_CQUAD
         ELSE IF (ETYPE(K) == 'QUAD4K  ') THEN   ;   MEDAT = MEDAT_CQUAD
         ELSE IF (ETYPE(K) == 'ROD     ') THEN   ;   MEDAT = MEDAT_CROD
         ELSE IF (ETYPE(K) == 'SHEAR   ') THEN   ;   MEDAT = MEDAT_CSHEAR
         ELSE IF (ETYPE(K) == 'TETRA4  ') THEN   ;   MEDAT = MEDAT_CTETRA4
         ELSE IF (ETYPE(K) == 'TETRA10 ') THEN   ;   MEDAT = MEDAT_CTETRA10
         ELSE IF (ETYPE(K) == 'TRIA3   ') THEN   ;   MEDAT = MEDAT_CTRIA
         ELSE IF (ETYPE(K) == 'TRIA3K  ') THEN   ;   MEDAT = MEDAT_CTRIA
         ELSE IF (ETYPE(K) == 'USER1   ') THEN   ;   MEDAT = MEDAT_CUSER1
         ELSE IF (ETYPE(K) == 'USERIN  ') THEN
            EPNTK = EPNT(K)
            NG    = EDAT(EPNTK+2)
            NS    = EDAT(EPNTK+3)
            MEDAT = MEDAT0_CUSERIN + 2*NG + NS + 1
            DO J=1,NG
               WRITE(ICHAR,'(I9)') J
               EDAT_DESCR(J+5,22) =  'Grid   ' // ICHAR
            ENDDO
            DO J=1,NS
               WRITE(ICHAR,'(I9)') J
               EDAT_DESCR(J+5+NG,22) = 'SPOINT ' // ICHAR
            ENDDO
            DO J=1,NG
               WRITE(ICHAR,'(I9)') J
               EDAT_DESCR(J+5+NG+NS,22) = 'Comps  ' // ICHAR
            ENDDO
            EDAT_DESCR(2*NG+NS+6,22) = '# boundary DOF  '
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,1001) ETYPE(K), SUBR_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1001) ETYPE(K), SUBR_NAME
            ENDIF
            EXIT do_1
         ENDIF

         DO J=1,METYPE
            IF (ETYPE(K) == EDAT_DESCR(0,J)(1:8)) THEN
               L = 1
               DO I=I1,I1+MEDAT-1
                  IF (I == I1) THEN
                     WRITE (F06,2003) K, EPNT(K), ETYPE(K), I, 'EPNTK', EDAT(I), EDAT_DESCR(L,J)
                  ELSE
                     IF ( L-1 < 10) THEN
                        IF (EDAT(I) == -99) THEN
                           WRITE (F06,2004) I, '"  +', L-1, EDAT(I), TRIM(EDAT_DESCR(L,J)) // EDAT_EXPLAIN_CID
                        ELSE IF ((EDAT(I) < 0) .AND. (TRIM(EDAT_DESCR(L,J)) == 'V-vector key')) THEN
                           WRITE (F06,2004) I, '"  +', L-1, EDAT(I), TRIM(EDAT_DESCR(L,J)) // EDAT_EXPLAIN_VVEC
                        ELSE
                           WRITE (F06,2004) I, '"  +', L-1, EDAT(I), EDAT_DESCR(L,J)
                        ENDIF
                     ELSE
                        WRITE (F06,2010) I, '"  +', L-1, EDAT(I), EDAT_DESCR(L,J)
                     ENDIF
                  ENDIF
                  L = L + 1
               ENDDO
               I1 = I1 + MEDAT
               WRITE(F06,*)
            ENDIF
         ENDDO

      ENDDO do_1

      IF (NCUSERIN > 0) THEN
         WRITE(F06,*)
         WRITE(F06,2100)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1001 FORMAT(' *WARNING    : ELEMENT TYPE ',A,' IS NOT INCLUDED IN SUBR ',A                                                        &
                   ,/,14X, ' CANNOT PRINT EDAT TABLE AS REQUESTED BY DEBUG(108)')

 2000 FORMAT( 6X,'-------------------------------------------------------------------------------------------------')

 2001 FORMAT( 5X,'|  A R R A Y S   E P N T   A N D   E D A T   O F   E L E M E N T   C O N N E C T I O N   D A T A  |')

 2002 FORMAT(24X,'EDAT is an array of element connection data for all elements',/,                                                 &
             14X,'EPNT is an array that gives the starting position in EDAT for internal element K',//,                            &
             17X,'K         EPNT(K)   Elem type        I     Loction in EDAT   EDAT(I)   Description',/,                           &
              9X,'Internal elem ID  (EPNTK)                         relative to EPNTK',/)

 2003 FORMAT(10X,I8,I13,7X,A8,I9,9X,A,I13,7X,A)

 2004 FORMAT(47X,I8,11X,A,I1,I11,7X,A)

 2010 FORMAT(47X,I8,11X,A,I2,I10,7X,A)

 2100 FORMAT(8X,'NOTE: In the above table for USERIN elements :',/,                                                                &
             8X,'      (1) "Comps  i" are the boundary displacement components for "Grid i"',/,                                    &
             8X,'      (2) "SPOINT i" are the scalar points for the modal DOF for the USERIN element',/)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_EDAT



      SUBROUTINE WRITE_ELM_OT4 ( MAT_NAME, NROWS_MAT, NROWS_TXT, NCOLS, TXT, UNT )

! Writes CB OTM text file that describes the rows of element related OTM matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NROWS_TXT         ! Number of rows in TXT

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_ELM_OT4'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name of the OTM that MAT describes
      CHARACTER(LEN=*), INTENT(IN)    :: TXT(NROWS_TXT)    ! Lines of this array describe the rows of CB OTM of MAT_NAME

      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols in MAT
      INTEGER(LONG), INTENT(IN)       :: NROWS_MAT         ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number where to write matrix
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
      WRITE(UNT,11) NROWS_MAT, NCOLS, MAT_NAME
      IF      (MAT_NAME(1:8) == 'OTM_ELFE') THEN
         WRITE(UNT,21)
      ELSE IF (MAT_NAME(1:8) == 'OTM_ELFN') THEN
         WRITE(UNT,22)
      ELSE IF (MAT_NAME(1:8) == 'OTM_STRE') THEN
         IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
            WRITE(UNT,231)
         ELSE
            WRITE(UNT,232)
         ENDIF
      ELSE IF (MAT_NAME(1:8) == 'OTM_STRN') THEN
         IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
            WRITE(UNT,241)
         ELSE
            WRITE(UNT,242)
         ENDIF
      ENDIF

      DO I=1,NROWS_TXT
         WRITE(UNT,31) TXT(I)
      ENDDO

      WRITE(UNT, *)



      RETURN

! **********************************************************************************************************************************
   11 FORMAT('-------------------------------------------------------------------------------------------------------------------',&
             '----------------',/,                                                                                                 &
             '       Explanation of rows of ',I8,' row by ',I8,' col matrix ',A,/)

   21 FORMAT('   ROW            DESCRIPTION             TYPE     EID              ITEM',/,                                         &
             ' ------- ------------------------------ -------- -------    --------------------')

   22 FORMAT('   ROW            DESCRIPTION             TYPE     EID     GRID     COMP',/,                                         &
             ' ------- ------------------------------ -------- ------- -------    ----')

  231 FORMAT('   ROW            DESCRIPTION             TYPE     EID        LOCATION            ITEM',/,                           &
             ' ------- ------------------------------ -------- -------    ------------   --------------------')

  232 FORMAT('   ROW            DESCRIPTION             TYPE     EID              ITEM',/,                                         &
             ' ------- ------------------------------ -------- -------    --------------------')

  241 FORMAT('   ROW            DESCRIPTION             TYPE     EID        LOCATION            ITEM',/,                           &
             ' ------- ------------------------------ -------- -------    ------------   --------------------')

  242 FORMAT('   ROW            DESCRIPTION             TYPE     EID              ITEM',/,                                         &
             ' ------- ------------------------------ -------- -------    --------------------')

   31 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_ELM_OT4


      SUBROUTINE WRITE_FIJFIL ( WHICH, JVEC )

! Writes elem matrices to unformatted disk files if disk file output for elem data is requested.
! User must have Case Control entries ELDATA in order to get these files written

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06, F21, F22, F23, F24, F25, F21_MSG, F22_MSG, F23_MSG, F24_MSG, F25_MSG
      USE DEBUG_PARAMETERS
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MAX_STRESS_POINTS, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EID, TYPE, ELGP, ELDOF, KE, ME, PEB, PEG, PEL, PPE, PTE,                                  &
                                         SE1, SE2, SE3, STE1, STE2, STE3, UEB, UEG, UEL
      USE PARAMS, ONLY                :  ELFORCEN

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: FILE_INQUIRE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FIJFIL'

      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Internal subcase or vector number for data to be written
      INTEGER(LONG), INTENT(IN)       :: WHICH             ! Which F2j file to write to
      INTEGER(LONG)                   :: I,J, K            ! DO loop indices




! **********************************************************************************************************************************

      IF (( DEBUG(193) == 1) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'In WRITE_FIJFIL' )
      ENDIF

! Write element thermal, pressure, loads to disk file if requested

      IF ((WHICH == 1) .OR. (WHICH == 9999)) THEN
         WRITE(F21) F21_MSG
         WRITE(F21) EID
         WRITE(F21) TYPE
         WRITE(F21) ELDOF
         WRITE(F21) NTSUB
         WRITE(F21) NSUB
         DO I=1,NTSUB
            DO J=1,ELDOF
               WRITE(F21) PTE(J,I)
            ENDDO
         ENDDO
         DO I=1,NSUB
            DO J=1,ELDOF
               WRITE(F21) PPE(J,I)
            ENDDO
         ENDDO
!xx      WRITE(F21) 'FINISHED'
      ENDIF

! Write element mass matrix to disk file if requested

      IF ((WHICH == 2) .OR. (WHICH == 9999)) THEN
         WRITE(F22) F22_MSG
         WRITE(F22) EID
         WRITE(F22) TYPE
         WRITE(F22) ELDOF
         DO I=1,ELDOF
            DO J=I,ELDOF
               WRITE(F22) ME(I,J)
            ENDDO
         ENDDO
!xx      WRITE(F22) 'FINISHED'
      ENDIF

! Write element stiffness matrix to disk file if requested

      IF ((WHICH == 3) .OR. (WHICH == 9999)) THEN
         WRITE(F23) F23_MSG
         WRITE(F23) EID
         WRITE(F23) TYPE
         WRITE(F23) ELDOF
         DO I=1,ELDOF
            DO J=I,ELDOF
               WRITE(F23) KE(I,J)
            ENDDO
         ENDDO
!xx      WRITE(F23) 'FINISHED'
      ENDIF

! Write element stress recovery matrices to disk file if requested

      IF ((WHICH == 4) .OR. (WHICH == 9999)) THEN

         WRITE(F24) F24_MSG
         WRITE(F24) EID
         WRITE(F24) TYPE
         WRITE(F24) ELDOF
         WRITE(F24) NTSUB
         WRITE(F24) MAX_STRESS_POINTS

         DO K=1,MAX_STRESS_POINTS+1
             DO I=1,3
                DO J=1,ELDOF
                   WRITE(F24) SE1(I,J,K)
                ENDDO
             ENDDO

             DO I=1,3
                DO J=1,ELDOF
                   WRITE(F24) SE2(I,J,K)
                ENDDO
             ENDDO

             DO I=1,3
                DO J=1,ELDOF
                   WRITE(F24) SE3(I,J,K)
                ENDDO
             ENDDO

             DO J=1,NTSUB
                DO I=1,3
                   WRITE(F24) STE1(I,J,K)
                ENDDO
             ENDDO

             DO J=1,NTSUB
                DO I=1,3
                   WRITE(F24) STE2(I,J,K)
                ENDDO
             ENDDO

             DO J=1,NTSUB
                DO I=1,3
                   WRITE(F24) STE3(I,J,K)
                ENDDO
             ENDDO

         ENDDO

!xx      WRITE(F24) 'FINISHED'
      ENDIF

! Write element loads, displ's to disk file if requested

      IF ((WHICH == 5) .OR. (WHICH == 9999)) THEN
         WRITE(F25) F25_MSG
         WRITE(F25) 'Displs and forces are in coord system: ', ELFORCEN
         WRITE(F25) EID
         WRITE(F25) TYPE
         WRITE(F25) ELDOF
         WRITE(F25) JVEC
         IF      (ELFORCEN == 'LOCAL') THEN
            DO I=1,ELDOF
               WRITE(F25) UEL(I),PEL(I)
            ENDDO
         ELSE IF (ELFORCEN == 'BASIC') THEN
            DO I=1,ELDOF
               WRITE(F25) UEB(I),PEB(I)
            ENDDO
         ELSE IF (ELFORCEN == 'GLOBAL') THEN
            DO I=1,ELDOF
               WRITE(F25) UEG(I),PEG(I)
            ENDDO
         ENDIF
!xx      WRITE(F25) 'FINISHED'
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_FIJFIL


      SUBROUTINE WRITE_GRD_OT4 ( MAT_NAME, NROWS_MAT, NROWS_TXT, NCOLS, TXT, UNT )

! Writes CB OTM text file that describes the rows of grid related OTM matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NROWS_TXT         ! Number of rows in TXT

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_GRD_OT4'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name of the OTM that MAT describes
      CHARACTER(LEN=*), INTENT(IN)    :: TXT(NROWS_TXT)    ! Lines of this array describe the rows of CB OTM of MAT_NAME

      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols in MAT
      INTEGER(LONG), INTENT(IN)       :: NROWS_MAT         ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number where to write matrix
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
      WRITE(UNT, 1) NROWS_MAT, NCOLS, MAT_NAME
      WRITE(UNT, 2)

      DO I=1,NROWS_TXT
         WRITE(UNT,3) TXT(I)
      ENDDO

      WRITE(UNT, *)



      RETURN

! **********************************************************************************************************************************
    1 FORMAT('-------------------------------------------------------------------------------------------------------------------',&
             '----------------',/,                                                                                                 &
             'Explanation of rows of ',I8,' row by ',I8,' col matrix ',A,/)

    2 FORMAT('   ROW            DESCRIPTION              GRID     COMP',/,                                                         &
             ' ------- ------------------------------  -------    ----')

    3 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_GRD_OT4


      SUBROUTINE WRITE_PARTNd_MAT_HDRS ( MAT_NAME, ROW_SET, COL_SET, NROWS, NCOLS )

! Writes the grid/comp pairs corresponding to the cols and rows of a partitioned matrix. Used for OUTPUT4 matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MTDOF, NDOFA, NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, NDOFO, NDOFR,   &
                                         NDOFS, NDOFSA, NDOFSB, NDOFSE, NDOFSG, NDOFSZ, NUM_USET_U1, NUM_USET_U2, TSET_CHR_LEN
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOFI
      USE OUTPUT4_MATRICES, ONLY      :  OU4_MAT_COL_GRD_COMP, OU4_MAT_ROW_GRD_COMP

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_PARTNd_MAT_HDRS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME            ! Name of the partitioned matrix whose row/col headings will be output
      CHARACTER(LEN=*), INTENT(IN)    :: COL_SET             ! Set name that the cols of MAT_NAME belong to (e.g. 'G ', 'l ', etc)
      CHARACTER(LEN=*), INTENT(IN)    :: ROW_SET             ! Set name that the rows of MAT_NAME belong to (e.g. 'G ', 'l ', etc)

      INTEGER(LONG), INTENT(IN)       :: NCOLS               ! Number of cols in the partitioned matrix MAT_NAME
      INTEGER(LONG), INTENT(IN)       :: NROWS               ! Number of rows in the partitioned matrix MAT_NAME
      INTEGER(LONG)                   :: J,K                 ! DO loop indices or counters
      INTEGER(LONG)                   :: OUTPUT_1(10)        ! Part of a line of output
      INTEGER(LONG)                   :: OUTPUT_2(10)        ! Part of a line of output
      INTEGER(LONG)                   :: NUM_LEFT            ! Used when printing a line of 10 values in the set




! **********************************************************************************************************************************
! Write the matrix col headers out to F06

      IF (COL_SET /= '--') THEN

         WRITE(F06,101) 'C O L S', MAT_NAME
         WRITE(F06,201) 'COLS', COL_SET

         NUM_LEFT = NCOLS
         DO J=1,NCOLS,10
            IF (NUM_LEFT >= 10) THEN
               DO K=1,10
                  OUTPUT_1(K) = OU4_MAT_COL_GRD_COMP(J+K-1,1)
                  OUTPUT_2(K) = OU4_MAT_COL_GRD_COMP(J+K-1,2)
               ENDDO
               WRITE(F06,401) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,10), J+9
            ELSE
               DO K=1,NUM_LEFT
                  OUTPUT_1(K) = OU4_MAT_COL_GRD_COMP(J+K-1,1)
                  OUTPUT_2(K) = OU4_MAT_COL_GRD_COMP(J+K-1,2)
               ENDDO
               IF (NUM_LEFT == 1) WRITE(F06,301) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 2) WRITE(F06,302) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 3) WRITE(F06,303) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 4) WRITE(F06,304) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 5) WRITE(F06,305) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 6) WRITE(F06,306) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 7) WRITE(F06,307) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 8) WRITE(F06,308) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 9) WRITE(F06,309) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
            ENDIF
            NUM_LEFT = NUM_LEFT - 10
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ENDIF

! Write the matrix row headers out to F06

      IF (ROW_SET /= '--') THEN

         WRITE(F06,101) 'R O W S', MAT_NAME
         WRITE(F06,201) 'ROWS', ROW_SET

         NUM_LEFT = NROWS
         DO J=1,NROWS,10
            IF (NUM_LEFT >= 10) THEN
               DO K=1,10
                  OUTPUT_1(K) = OU4_MAT_ROW_GRD_COMP(J+K-1,1)
                  OUTPUT_2(K) = OU4_MAT_ROW_GRD_COMP(J+K-1,2)
               ENDDO
               WRITE(F06,401) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,10), J+9
            ELSE
               DO K=1,NUM_LEFT
                  OUTPUT_1(K) = OU4_MAT_ROW_GRD_COMP(J+K-1,1)
                  OUTPUT_2(K) = OU4_MAT_ROW_GRD_COMP(J+K-1,2)
               ENDDO
               IF (NUM_LEFT == 1) WRITE(F06,301) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 2) WRITE(F06,302) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 3) WRITE(F06,303) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 4) WRITE(F06,304) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 5) WRITE(F06,305) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 6) WRITE(F06,306) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 7) WRITE(F06,307) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 8) WRITE(F06,308) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
               IF (NUM_LEFT == 9) WRITE(F06,309) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
            ENDIF
            NUM_LEFT = NUM_LEFT - 10
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ENDIF


      RETURN

! **********************************************************************************************************************************
  101 FORMAT(5X,'G R I D - C O M P   N U M B E R S   F O R   T H E   ',A,'   O F   P A R T I T I O N E D   M A T R I X  "',A,'"')

  201 FORMAT(47X,A,' ARE FROM THE "',A2,'"',' DISPLACEMENT SET',/,                                                                 &
               16X,'-1-        -2-        -3-        -4-        -5-        -6-        -7-        -8-        -9-       -10-',/)

  301 FORMAT(1X,I6,'=', 1(I9,'-',I1),9(11X))

  302 FORMAT(1X,I6,'=', 2(I9,'-',I1),8(11X))

  303 FORMAT(1X,I6,'=', 3(I9,'-',I1),7(11X))

  304 FORMAT(1X,I6,'=', 4(I9,'-',I1),6(11X))

  305 FORMAT(1X,I6,'=', 5(I9,'-',I1),5(11X))

  306 FORMAT(1X,I6,'=', 6(I9,'-',I1),4(11X))

  307 FORMAT(1X,I6,'=', 7(I9,'-',I1),3(11X))

  308 FORMAT(1X,I6,'=', 8(I9,'-',I1),2(11X))

  309 FORMAT(1X,I6,'=', 9(I9,'-',I1),1(11X))

  401 FORMAT(1X,I6,'=',10(I9,'-',I1),'      =',I6)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_PARTNd_MAT_HDRS


      SUBROUTINE WRITE_DOF_TABLES

! Writess DOF table data to file LINK1C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  L1C, LINK1C, L1C_MSG, ERR, F06
      USE SCONTR, ONLY                :  DATA_NAM_LEN, MTDOF, NDOFG, NGRID, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOFI, TDOF, TSET
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_DOF_TABLES'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME      ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters

      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN



! **********************************************************************************************************************************
! Write TSET, TDOF, TDOFI tables to file L1C

      ! Ensure L1C is open. if this is a two-step run, we need to reopen L1C.
      IF (LOAD_ISTEP > 1) THEN
         OUNT(1) = ERR
         OUNT(2) = F06
         CALL FILE_OPEN ( L1C, LINK1C, OUNT, 'REPLACE', L1C_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      DATA_SET_NAME = 'TSET'
      WRITE(L1C) DATA_SET_NAME
      WRITE(L1C) NGRID
      WRITE(L1C) ((TSET(I,J),J=1,6),I=1,NGRID)            ! One record per table (one per entry took about 2 s for 200k DOF)
      DATA_SET_NAME = 'TDOFI'
      WRITE(L1C) DATA_SET_NAME
      WRITE(L1C) NDOFG
      WRITE(L1C) MTDOF
      WRITE(L1C) ((TDOFI(I,J),J=1,MTDOF),I=1,NDOFG)
      DATA_SET_NAME = 'TDOF'
      WRITE(L1C) DATA_SET_NAME
      WRITE(L1C) NDOFG
      WRITE(L1C) MTDOF
      WRITE(L1C) ((TDOF(I,J),J=1,MTDOF),I=1,NDOFG)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_DOF_TABLES


      SUBROUTINE WRITE_TDOF ( TDOF_MSG )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MTDOF, NDOFG, NDOFM, NDOFN, NDOFSA, NDOFSB, NDOFSG, NDOFSZ, NDOFSE, NDOFS,  &
                                         NDOFF, NDOFO, NDOFA, NDOFR, NDOFL, NGRID, NUM_USET_U1, NUM_USET_U2
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI
      USE PARAMS, ONLY                :  PRTDOF
      USE MODEL_STUF, ONLY            :  GRID_ID, INV_GRID_SEQ

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_TDOF'
      CHARACTER(LEN=*), INTENT(IN)    :: TDOF_MSG          ! Message to be printed out regarding at what point in the run the TDOF,I
!                                                            tables are printed out
      CHARACTER( 10*BYTE)             :: NAME1             ! A data set name for output purposes
      CHARACTER(  7*BYTE)             :: NAME2             ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IROW              ! Row number in array TDOF or TDOFI
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)




! **********************************************************************************************************************************
! Table TDOF is printed in the F06 file if B.D. PARAM PRTDOF = 1 or 3

      IF ((PRTDOF == 1) .OR. (PRTDOF == 3)) THEN
         NAME1 = 'GRID POINT'
         NAME2 = '(TDOF)'
         IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
            WRITE(F06,102) NAME1, NAME2, TDOF_MSG
         ELSE
            WRITE(F06,101) NAME1, NAME2, TDOF_MSG
         ENDIF
         IROW = 0
         DO I = 1,NGRID
            CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
            DO J = 1,NUM_COMPS
               IROW = IROW + 1
               IF (J == 1) THEN
                  IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
                     WRITE(F06,202) (TDOF(IROW,K),K = 1,MTDOF)
                  ELSE
                     WRITE(F06,201) (TDOF(IROW,K),K = 1,MTDOF-2)
                  ENDIF
               ELSE
                  IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
                     WRITE(F06,302) TDOF(IROW,2),TDOF(IROW,4),(TDOF(IROW,K),K = 5,MTDOF)
                  ELSE
                     WRITE(F06,301) TDOF(IROW,2),TDOF(IROW,4),(TDOF(IROW,K),K = 5,MTDOF-2)
                  ENDIF
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
            WRITE(F06,402) NDOFG, NDOFM, NDOFN, NDOFSA, NDOFSB, NDOFSG, NDOFSZ, NDOFSE, NDOFS, NDOFF, NDOFO, NDOFA, NDOFR, NDOFL,  &
                           NUM_USET_U1, NUM_USET_U2
         ELSE
            WRITE(F06,401) NDOFG, NDOFM, NDOFN, NDOFSA, NDOFSB, NDOFSG, NDOFSZ, NDOFSE, NDOFS, NDOFF, NDOFO, NDOFA, NDOFR, NDOFL
         ENDIF
         WRITE(F06,'(//)')
      ENDIF

! Table TDOFI is printed in the F06 file if B.D. PARAM PRTDOF = 2 or 3

      IF ((PRTDOF == 2) .OR. (PRTDOF == 3)) THEN
         NAME1 = 'DOF NUMBER'
         NAME2 = '(TDOFI)'
         IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
            WRITE(F06,102) NAME1, NAME2, TDOF_MSG
         ELSE
            WRITE(F06,101) NAME1, NAME2, TDOF_MSG
         ENDIF
         IROW = 0
         DO I = 1,NGRID
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J = 1,NUM_COMPS
               IROW = IROW + 1
               IF (J == 1) THEN
                  IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
                     WRITE(F06,202) (TDOFI(IROW,K),K = 1,MTDOF)
                  ELSE
                     WRITE(F06,201) (TDOFI(IROW,K),K = 1,MTDOF-2)
                  ENDIF
               ELSE
                  IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
                     WRITE(F06,302) TDOFI(IROW,2),TDOFI(IROW,4),(TDOFI(IROW,K),K = 5,MTDOF)
                  ELSE
                     WRITE(F06,301) TDOFI(IROW,2),TDOFI(IROW,4),(TDOFI(IROW,K),K = 5,MTDOF-2)
                  ENDIF
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         IF ((NUM_USET_U1 > 0) .OR. (NUM_USET_U2 > 0)) THEN
            WRITE(F06,402) NDOFG, NDOFM, NDOFN, NDOFSA, NDOFSB, NDOFSG, NDOFSZ, NDOFSE, NDOFS, NDOFF, NDOFO, NDOFA, NDOFR, NDOFL,  &
                           NUM_USET_U1, NUM_USET_U2
         ELSE
            WRITE(F06,401) NDOFG, NDOFM, NDOFN, NDOFSA, NDOFSB, NDOFSG, NDOFSZ, NDOFSE, NDOFS, NDOFF, NDOFO, NDOFA, NDOFR, NDOFL
         ENDIF
         WRITE(F06,'(//)')
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(6X,'                                         DEGREE OF FREEDOM TABLE SORTED ON ',A,1X,A,/,7X,A                        &
          ,/                                                                                                                       &
          ,/,1X,                                                                                                                   &
'EXTERNAL  INTERNAL                                     DOF NUMBER FOR DISPLACEMENT SET:'                                          &
          ,/,1X,                                                                                                                   &
'GRD-COMP  GRD-COMP --------------------------------------------------------------------------------------------------------',     &
                   '--------'                                                                                                      &
          ,/,1X,                                                                                                                   &
'NUMBER    NUMBER          G       M       N      SA      SB      SG      SZ      SE       S       F       O       A       ',      &
                          'R       L',/)
  102 FORMAT(6X,'                                         DEGREE OF FREEDOM TABLE SORTED ON ',A,1X,A,/,7X,A                        &
          ,/                                                                                                                       &
          ,/,1X,                                                                                                                   &
'EXTERNAL  INTERNAL                                     DOF NUMBER FOR DISPLACEMENT SET:'                                          &
          ,/,1X,                                                                                                                   &
'GRD-COMP  GRD-COMP --------------------------------------------------------------------------------------------------------',     &
                   '------------------------'                                                                                      &
          ,/,1X,                                                                                                                   &
'NUMBER    NUMBER          G       M       N      SA      SB      SG      SZ      SE       S       F       O       A       ',      &
                          'R       L      U1      U2',/)

  201 FORMAT(2(I8,'-',I1),14(I8))

  202 FORMAT(2(I8,'-',I1),16(I8))

  301 FORMAT(2(8X,'-',I1),14(I8))

  302 FORMAT(2(8X,'-',I1),16(I8))

  401 FORMAT(8X,'             ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------',    &
                ' ------- -------'                                                                                              ,/,&
             'TOTAL NUMBER OF DOF:',14(I8)                                                                                 ,//,    &
             27X,        'G       M       N      SA      SB      SG      SZ      SE       S       F       O       A       ',       &
                         'R       L',/)

  402 FORMAT(8X,'             ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------',    &
                ' ------- ------- ------- -------'                                                                              ,/,&
             'TOTAL NUMBER OF DOF:',16(I8)                                                                                 ,//,    &
             27X,        'G       M       N      SA      SB      SG      SZ      SE       S       F       O       A       ',       &
                         'R       L      U1      U2',/)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_TDOF


      SUBROUTINE WRITE_TSET

! Writes the NGRID x 6 TSET degree of freedom table to the F06 file based on user supplied Bulk Data Param PRTTSET

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MTSET, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID, GRID_SEQ, INV_GRID_SEQ
      USE DOF_TABLES, ONLY            :  TSET

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_TSET'

      INTEGER(LONG)                   :: I,J               ! DO loop indices




! *********************************************************************************************************************************
      WRITE(F06,56)
      WRITE(F06,57)
      DO I=1,NGRID
         WRITE(F06,58) GRID(I,1), GRID_SEQ(I), INV_GRID_SEQ(I), (TSET(I,J),J = 1,MTSET)
      ENDDO
      WRITE(F06,59)



      RETURN

! *********************************************************************************************************************************
   56 FORMAT(54X,'DEGREE OF FREEDOM SET TABLE (TSET)',/)

   57 FORMAT(28X,'GRID    GRID_SEQ   INV_GRID_SEQ*      T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(24X,I8,I12,I15,6(7X,A2))

   59 FORMAT(1X,/,29X,'* INV_GRID_SEQ = J meams that the J-th entry under GRID is sequenced GRID_SEQ(J)',//)

! *********************************************************************************************************************************

      END SUBROUTINE WRITE_TSET


      SUBROUTINE WRITE_USET

! Writes the NGRID x 6 USET table to the F06 file based on user supplied Bulk Data Param PRTUSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MTSET, NDOFG, NGRID, NUM_USET_U1, NUM_USET_U2
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID, GRID_SEQ, INV_GRID_SEQ
      USE PARAMS, ONLY                :  PRTUSET
      USE DOF_TABLES, ONLY            :  TDOF, USET, USETSTR_TABLE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_USET'

      INTEGER(LONG)                   :: I,J               ! DO loop indices




! **********************************************************************************************************************************
! Write the USET table

      IF (PRTUSET > 0) THEN

         WRITE(F06,56)
         WRITE(F06,57)

i_do:    DO I=1,NGRID
            IF ((USET(I,1)(1:1) == 'U') .OR. (USET(I,2)(1:1) == 'U') .OR. (USET(I,3)(1:1) == 'U') .OR.                             &
               (USET(I,4)(1:1) == 'U') .OR. (USET(I,5)(1:1) == 'U') .OR. (USET(I,6)(1:1) == 'U')) THEN
               WRITE(F06,58) GRID(I,1), (USET(I,J),J = 1,MTSET)
            ELSE
               CYCLE i_do
            ENDIF
         ENDDO i_do
         WRITE(F06,'(//)')

      ENDIF



      RETURN

! **********************************************************************************************************************************
   56 FORMAT(50X,'DEGREES OF FREEDOM DEFINED ON USET BULK DATA ENTRIES'                                                            &
          ,/,50X,'----------------------------------------------------',/)

   57 FORMAT(42x,'     Grid       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(42X,I9,32767(7X,A2))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_USET


      SUBROUTINE WRITE_USETSTR

! Write USET definition tables for any of the MYSTRAN displ sets (G-set, M-set, etc). The user must have had a Bulk Data
! PARAM USETSTR entry with a dilpl set (like G) in field 3. Any number of these PARAM entries for valid displ sets will cause that
! displ set table to be written

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MTDOF, NDOFA, NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, NDOFO, NDOFR,   &
                                         NDOFS, NDOFSA, NDOFSB, NDOFSE, NDOFSG, NDOFSZ, NUM_USET_U1, NUM_USET_U2, TSET_CHR_LEN
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOFI, USETSTR_TABLE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_LOOKUP_UTILS, ONLY      :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_USETSTR'
      CHARACTER(LEN(TSET_CHR_LEN))    :: CHAR_SET            ! Name of a DOF set: G, N, M, etc
      CHARACTER(LEN(TSET_CHR_LEN))    :: NULL_SET_NAME(MTDOF)! Array of names of sets that are null (for output message purposes)
      CHARACTER( 1*BYTE)              :: USETSTR_OUTPUT      ! If 'Y' then output of USET tables is requested

      INTEGER(LONG)                   :: COL_NUM             ! Column number in TDOF where a DOF set exists
      INTEGER(LONG)                   :: GRID_NUM(NDOFG)     ! Array of grid numbers for members of a DOF set requested in USETSTR
      INTEGER(LONG)                   :: COMP_NUM(NDOFG)     ! Array of comp numbers for members of a DOF set requested in USETSTR
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_LEFT            ! Used when printing a line of 10 values in the set
      INTEGER(LONG)                   :: NUM_IN_SET          ! A set length (e.g. NDOFM for the M-set
      INTEGER(LONG)                   :: NUM_NULL            ! Number of sets that have been requested for output that are null
      INTEGER(LONG)                   :: OUTPUT_1(10)        ! Part of a line of output
      INTEGER(LONG)                   :: OUTPUT_2(10)        ! Part of a line of output




! **********************************************************************************************************************************
! Initialize

      NUM_NULL = 0

! Scan table USETSTR to see if there are output requests

      USETSTR_OUTPUT = 'N'
      DO I=1,16
         IF (USETSTR_TABLE(I,2) /= '0 ') THEN
            USETSTR_OUTPUT = 'Y'
         ENDIF
      ENDDO

! If there were requests, process them and output them to F06

      IF (USETSTR_OUTPUT == 'Y') THEN

         WRITE(F06,200)
         DO I=1,MTDOF-4

            IF (USETSTR_TABLE(I,2) /= '0 ') THEN           ! If /= 0 then a request for that DOF set was made via B.D. PARAM USETSTR

               CHAR_SET = USETSTR_TABLE(I,1)

               IF      (CHAR_SET == 'G ') THEN   ;   NUM_IN_SET = NDOFG
               ELSE IF (CHAR_SET == 'M ') THEN   ;   NUM_IN_SET = NDOFM
               ELSE IF (CHAR_SET == 'N ') THEN   ;   NUM_IN_SET = NDOFN
               ELSE IF (CHAR_SET == 'SA') THEN   ;   NUM_IN_SET = NDOFSA
               ELSE IF (CHAR_SET == 'SB') THEN   ;   NUM_IN_SET = NDOFSB
               ELSE IF (CHAR_SET == 'SG') THEN   ;   NUM_IN_SET = NDOFSG
               ELSE IF (CHAR_SET == 'SZ') THEN   ;   NUM_IN_SET = NDOFSZ
               ELSE IF (CHAR_SET == 'SE') THEN   ;   NUM_IN_SET = NDOFSE
               ELSE IF (CHAR_SET == 'S ') THEN   ;   NUM_IN_SET = NDOFS
               ELSE IF (CHAR_SET == 'F ') THEN   ;   NUM_IN_SET = NDOFF
               ELSE IF (CHAR_SET == 'O ') THEN   ;   NUM_IN_SET = NDOFO
               ELSE IF (CHAR_SET == 'A ') THEN   ;   NUM_IN_SET = NDOFA
               ELSE IF (CHAR_SET == 'R ') THEN   ;   NUM_IN_SET = NDOFR
               ELSE IF (CHAR_SET == 'L ') THEN   ;   NUM_IN_SET = NDOFL
               ELSE IF (CHAR_SET == 'U1') THEN   ;   NUM_IN_SET = NUM_USET_U1
               ELSE IF (CHAR_SET == 'U2') THEN   ;   NUM_IN_SET = NUM_USET_U2
               ELSE
                  WRITE(ERR,1327) SUBR_NAME,CHAR_SET
                  WRITE(F06,1327) SUBR_NAME,CHAR_SET
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )                   ! Coding error, so quit
               ENDIF

               IF (NUM_IN_SET > 0) THEN

                  WRITE(F06,201) CHAR_SET
                  DO J=1,NDOFG
                     GRID_NUM(J) = 0
                     COMP_NUM(J) = 0
                  ENDDO
                  CALL TDOF_COL_NUM ( CHAR_SET, COL_NUM )
                  K = 0
                  DO J=1,NDOFG
                     IF (TDOFI(J,COL_NUM) /= 0) THEN
                        K = K + 1
                        GRID_NUM(K) = TDOFI(J,1)
                        COMP_NUM(K) = TDOFI(J,2)
                     ENDIF
                  ENDDO

                  NUM_LEFT = NUM_IN_SET
                  DO J=1,NUM_IN_SET,10
                     IF (NUM_LEFT >= 10) THEN
                        DO K=1,10
                           OUTPUT_1(K) = GRID_NUM(J+K-1)
                           OUTPUT_2(K) = COMP_NUM(J+K-1)
                        ENDDO
                        WRITE(F06,401) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,10), J+9
                     ELSE
                        DO K=1,NUM_LEFT
                           OUTPUT_1(K) = GRID_NUM(J+K-1)
                           OUTPUT_2(K) = COMP_NUM(J+K-1)
                        ENDDO
                        IF (NUM_LEFT == 1) WRITE(F06,301) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 2) WRITE(F06,302) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 3) WRITE(F06,303) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 4) WRITE(F06,304) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 5) WRITE(F06,305) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 6) WRITE(F06,306) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 7) WRITE(F06,307) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 8) WRITE(F06,308) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                        IF (NUM_LEFT == 9) WRITE(F06,309) J, (OUTPUT_1(K), OUTPUT_2(K), K=1,NUM_LEFT)
                     ENDIF
                     NUM_LEFT = NUM_LEFT - 10
                  ENDDO
                  WRITE(F06,*)
                  WRITE(F06,*)

               ELSE

                  NUM_NULL = NUM_NULL + 1
                  NULL_SET_NAME(NUM_NULL) = CHAR_SET

               ENDIF

            ENDIF

         ENDDO

         IF (NUM_NULL > 0) THEN
            DO I=1,NUM_NULL
               WRITE(ERR,901) NULL_SET_NAME(I)
               WRITE(F06,901) NULL_SET_NAME(I)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF




      RETURN

! **********************************************************************************************************************************
  200 FORMAT(16X,'U S E T   D E F I N I T I O N   T A B L E   ( I N T E R N A L   S E Q U E N C E ,   R O W   S O R T )')

  201 FORMAT(55X,A2,7X,'DISPLACEMENT SET',/,                                                                                       &
               16X,'-1-        -2-        -3-        -4-        -5-        -6-        -7-        -8-        -9-       -10-',/)

  301 FORMAT(1X,I6,'=', 1(I9,'-',I1),9(11X))

  302 FORMAT(1X,I6,'=', 2(I9,'-',I1),8(11X))

  303 FORMAT(1X,I6,'=', 3(I9,'-',I1),7(11X))

  304 FORMAT(1X,I6,'=', 4(I9,'-',I1),6(11X))

  305 FORMAT(1X,I6,'=', 5(I9,'-',I1),5(11X))

  306 FORMAT(1X,I6,'=', 6(I9,'-',I1),4(11X))

  307 FORMAT(1X,I6,'=', 7(I9,'-',I1),3(11X))

  308 FORMAT(1X,I6,'=', 8(I9,'-',I1),2(11X))

  309 FORMAT(1X,I6,'=', 9(I9,'-',I1),1(11X))

  401 FORMAT(1X,I6,'=',10(I9,'-',I1),'      =',I6)

  901 FORMAT(' *INFORMATION: BULK DATA PARAM USETSTR HAD REQUEST FOR USET DEFINITION TABLE FOR THE "',A,'" SET.',                  &
                           ' HOWEVER, THAT SET IS NULL')

 1327 FORMAT(' *ERROR  1327: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INPUT VARIABLE CHAR_SET = "',A2,'" IS NOT ONE OF THE CORRECT DESIGNATIONS FOR A DISPL SET')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_USETSTR


      SUBROUTINE WRITE_USERIN_BD_CARDS ( NROWS, X_SET )

! Write out B.D. CUSERIN card images, created in MYSTRAN for use in INPUT4 for a substructure, if user has requested it via Bulk
! Data PARAM CUSERIN,

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, F06FIL, MOU4, OU4, OU4FIL
      USE SCONTR, ONLY                :  JCARD_LEN, NCORD, NDOFG, NGRID, NVEC, WARN_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  START_YEAR, START_MONTH, START_DAY, START_HOUR, START_MINUTE, START_SEC, START_SFRAC
      USE DOF_TABLES, ONLY            :  TDOFI
      USE PARAMS, ONLY                :  CUSERIN_EID, CUSERIN_IN4, CUSERIN_PID, CUSERIN_SPNT_ID, CUSERIN_XSET,                     &
                                         CUSERIN_COMPTYP, SUPWARN
      USE MODEL_STUF, ONLY            :  CORD, RCORD, GRID_ID, GRID, RGRID
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_OUTPUT_NAMES, NUM_OU4_REQUESTS, OU4_FILE_UNITS
      USE DOF_LOOKUP_UTILS, ONLY      :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE BDF_FIELD_VALIDATION, ONLY  :  LEFT_ADJ_BDFLD

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NROWS                ! Size that is at least as large as any of the arrays below need

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_USERIN_BD_ENTRIES'
      CHARACTER(LEN=*), INTENT(IN)    :: X_SET                ! Displ set that the USERIN element is connected to
      CHARACTER( 1*BYTE)              :: COMPS(6)             ! Array of displ comps for 1 grid in the USERIN_GRIDS array
      CHARACTER( 1*BYTE)              :: CORD_FND             ! Indicator of whether we found a coord sys in array CORD

      CHARACTER( 8*BYTE)              :: CORD_MSG = '"CID0"  '! Field 6 of CUSERIN and field 3 of CORD2R to indicate to user that
!                                                               this should be replaced by act coord sys ID which defines USERIN
!                                                               element basic coord sys relative to the OA model basic coord sys

      CHARACTER(LEN=JCARD_LEN)        :: ICARD(NROWS,9)       ! Char array for fields 2-9 of CUSERIN or PUSERIN B.D. entries
      CHARACTER( 8*BYTE)              :: USERIN_COMPS(NROWS)  ! Char array of all of the displ comps for 1 USERIN_GRID

      INTEGER(LONG)                   :: CORD_ID              ! ID of a coord system on a grid for the global coord sys ID
      INTEGER(LONG)                   :: FIRST                ! First of the USERIN_CORDS that is nonzero
      INTEGER(LONG)                   :: I,J                  ! DO loop indices
      INTEGER(LONG)                   :: IROW                 ! Row number in array GRID_ID
      INTEGER(LONG)                   :: ICORD                ! Internal coord sys ID
      INTEGER(LONG)                   :: K                    ! Counter
      INTEGER(LONG)                   :: NUM_LEFT             ! Counter
      INTEGER(LONG)                   :: NUM_USERIN_GRIDS     ! Number of different R set grids
      INTEGER(LONG)                   :: X_SET_COMPS(NROWS)   ! Array of displ components in the R set
      INTEGER(LONG)                   :: X_SET_GRIDS(NROWS)   ! Array of grids in the R set
      INTEGER(LONG)                   :: X_SET_COL            ! Col in TDOF where R set exists
      INTEGER(LONG)                   :: USERIN_CORDS(NROWS,2)! Array of global coord ID's from X_SET_GRIDS (act, int ID's)
      INTEGER(LONG)                   :: USERIN_GRIDS(NROWS)  ! Array of different grid values from X_SET_GRIDS

! **********************************************************************************************************************************
      CALL TDOF_COL_NUM ( X_SET, X_SET_COL )

! Initialize

      DO I=1,NROWS
         X_SET_GRIDS(I)    = 0
         X_SET_COMPS(I)    = 0
         USERIN_COMPS(I)   = '        '
         USERIN_GRIDS(I)   = 0
         USERIN_CORDS(I,1) = 0   ;   USERIN_CORDS(I,2) = 0
      ENDDO

! Get array of X set grids and displ comps

      K = 0
      DO I=1, NDOFG
         IF (TDOFI(I,X_SET_COL) > 0) THEN
            K = K + 1
            X_SET_GRIDS(K) = TDOFI(I,1)
            X_SET_COMPS(K) = TDOFI(I,2)
         ENDIF
      ENDDO


! Determine number of different X_SET_GRIDS

      NUM_USERIN_GRIDS = 1
      USERIN_GRIDS(1) = X_SET_GRIDS(1)
      DO I=2,NROWS
         IF (X_SET_GRIDS(I) /= X_SET_GRIDS(I-1)) THEN
            NUM_USERIN_GRIDS = NUM_USERIN_GRIDS + 1
            USERIN_GRIDS(NUM_USERIN_GRIDS) = X_SET_GRIDS(I)
         ENDIF
      ENDDO


! Get list of displ comps for each of the different grids in USERIN_GRIDS

      DO I=1,NUM_USERIN_GRIDS

         USERIN_COMPS(I) = '        '

         DO J=1,6
            COMPS(J) = ' '
         ENDDO

         DO J=1,NROWS
            IF (X_SET_GRIDS(J) == USERIN_GRIDS(I)) THEN
               IF      (X_SET_COMPS(J) == 1) THEN
                  COMPS(1) = '1'
               ELSE IF (X_SET_COMPS(J) == 2) THEN
                  COMPS(2) = '2'
               ELSE IF (X_SET_COMPS(J) == 3) THEN
                  COMPS(3) = '3'
               ELSE IF (X_SET_COMPS(J) == 4) THEN
                  COMPS(4) = '4'
               ELSE IF (X_SET_COMPS(J) == 5) THEN
                  COMPS(5) = '5'
               ELSE IF (X_SET_COMPS(J) == 6) THEN
                  COMPS(6) = '6'
               ENDIF
            ENDIF
         ENDDO

         IF (CUSERIN_COMPTYP == 0) THEN
            K = 0
            DO J=1,6
               IF (COMPS(J) /= ' ') THEN
                  K = K + 1
                  USERIN_COMPS(I)(K:K) = COMPS(J)
               ENDIF
            ENDDO
         ELSE
            DO J=1,6
               USERIN_COMPS(I)(J:J) = COMPS(J)
            ENDDO
         ENDIF

      ENDDO


! Get list of different global coord systems (NOTE: we don't need input coord systems since the GRID B.D. entries written out
! below have grid coords in basic system; which is the way RGRID is after subr GRID_PROC runs)

      DO I=1,NUM_USERIN_GRIDS
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, USERIN_GRIDS(I), IROW )
         CORD_ID = GRID(IROW,3)
         USERIN_CORDS(i,1) = CORD_ID
      ENDDO

! Now get the internal ID for the actual coord sys ID's in USERIN_CORDS(I,1)

      DO I=1,NUM_USERIN_GRIDS
         IF (USERIN_CORDS(I,1) /= 0) THEN
            CORD_FND = 'N'
               DO J=1,NCORD
               IF (USERIN_CORDS(I,1) == CORD(J,2)) THEN
                  CORD_FND = 'Y'
                  USERIN_CORDS(I,2) = J
                  EXIT
               ENDIF
            ENDDO
            IF (CORD_FND == 'N') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,9988) USERIN_CORDS(I,1), SUBR_NAME
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,9988) USERIN_CORDS(I,1), SUBR_NAME
               ENDIF
            ENDIF
         ENDIF
      ENDDO


      WRITE(F06,*)
      WRITE(F06,101)
      WRITE(F06,103) CUSERIN_EID, CUSERIN_XSET, F06FIL, START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE,      &
                     START_SEC, START_SFRAC
      WRITE(F06,104)

      DO I=1,NROWS
         ICARD(I,1) = '        '
         ICARD(I,2) = '        '
      ENDDO

! Write grid entries

      WRITE(F06,*)
      DO I=1,NUM_USERIN_GRIDS

         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, USERIN_GRIDS(I), IROW )
                                                           ! Grid ID
         WRITE(ICARD(I,2),201) USERIN_GRIDS(I)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,2) )

         WRITE(ICARD(I,3),203) CORD_MSG

         DO J=1,3                                          ! Grid coordinates
            WRITE(ICARD(I,J+3),204) RGRID(I,J)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,J+3) )
         ENDDO

         IF (GRID(I,3) /= 0) THEN                          ! Global coord sys
            WRITE(ICARD(I,7),201) GRID(I,3)           ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,7) )
         ELSE
            WRITE(ICARD(I,7),203) '        '
         ENDIF

         WRITE(F06,205) (ICARD(I,J),J=2,7)

      ENDDO
      WRITE(F06,*)

! Write CORD entries (for grid global coord systems)

      FIRST = 1
      DO I=1,NUM_USERIN_GRIDS
         IF (USERIN_CORDS(I,1) /= 0) THEN
            FIRST = I
            EXIT
         ENDIF
      ENDDO

      IF (USERIN_CORDS(FIRST,1) /= 0) THEN

         WRITE(ICARD(FIRST,2),201) USERIN_CORDS(FIRST,1)      ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,2) )

         WRITE(ICARD(FIRST,3),203) CORD_MSG

         ICORD = USERIN_CORDS(FIRST,2)
         DO J=1,3
            WRITE(ICARD(FIRST,J+3),204) RCORD(ICORD,J)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,J+3) )
         ENDDO

         WRITE(ICARD(FIRST,7),204) RCORD(ICORD, 6)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,7) )
         WRITE(ICARD(FIRST,8),204) RCORD(ICORD, 9)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,8) )
         WRITE(ICARD(FIRST,9),204) RCORD(ICORD,12)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,9) )

         WRITE(F06,206) (ICARD(FIRST,J),J=2,9)

         WRITE(ICARD(FIRST,2),204) RCORD(ICORD, 4)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,2) )
         WRITE(ICARD(FIRST,3),204) RCORD(ICORD, 7)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,3) )
         WRITE(ICARD(FIRST,4),204) RCORD(ICORD,10)            ;   CALL LEFT_ADJ_BDFLD ( ICARD(FIRST,4) )

         WRITE(F06,207) (ICARD(FIRST,J),J=2,4)

      ENDIF

      DO I=FIRST+1,NUM_USERIN_GRIDS

         IF ((USERIN_CORDS(I,1) /= USERIN_CORDS(I-1,1)) .AND. (USERIN_CORDS(I,1) /= 0)) THEN

            WRITE(ICARD(I,2),201) USERIN_CORDS(I,1)        ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,2) )

            WRITE(ICARD(I,3),203) CORD_MSG

            ICORD = USERIN_CORDS(I,2)
            DO J=1,3
               WRITE(ICARD(I,J+3),203) '0.      '
            ENDDO

            DO J=1,3
               WRITE(ICARD(I,J+6),204) RCORD(ICORD,J+9)    ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,J+6) )
            ENDDO

            WRITE(F06,206) (ICARD(I,J),J=2,9)

            DO J=1,3
               WRITE(ICARD(I,J+1),204) RCORD(ICORD,J+3)    ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,J+3) )
            ENDDO

            WRITE(F06,207) (ICARD(I,J),J=2,4)

         ENDIF

      ENDDO
      WRITE(F06,*)

! Write SPOINT entry

      IF (NVEC > 0) THEN
         WRITE(ICARD(2,1),201) CUSERIN_SPNT_ID             ;   CALL LEFT_ADJ_BDFLD ( ICARD(2,1) )
         WRITE(ICARD(4,1),201) CUSERIN_SPNT_ID+NVEC-1      ;   CALL LEFT_ADJ_BDFLD ( ICARD(4,1) )
         WRITE(F06,211) ICARD(2,1), ICARD(4,1)
      ENDIF

! Write CUSERIN entry

      WRITE(ICARD(2,1),201) CUSERIN_EID                    ;   CALL LEFT_ADJ_BDFLD ( ICARD(2,1) )
      WRITE(ICARD(3,1),201) CUSERIN_PID                    ;   CALL LEFT_ADJ_BDFLD ( ICARD(3,1) )
      WRITE(ICARD(4,1),201) NUM_USERIN_GRIDS               ;   CALL LEFT_ADJ_BDFLD ( ICARD(4,1) )
      WRITE(ICARD(5,1),201) NVEC                           ;   CALL LEFT_ADJ_BDFLD ( ICARD(5,1) )
      WRITE(ICARD(6,1),203) CORD_MSG
      WRITE(F06,202) (ICARD(I,1),I=2,6)

      DO I=1,NUM_USERIN_GRIDS                              ! Write CUSERIN cont cards for grids to internal file (ICARD)
         WRITE(ICARD(I,1),201) USERIN_GRIDS(I)             ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,1) )
         WRITE(ICARD(I,2),203) USERIN_COMPS(I)             ;   CALL LEFT_ADJ_BDFLD ( ICARD(I,2) )
      ENDDO

      NUM_LEFT = NUM_USERIN_GRIDS                          ! Write grids/comps (from internal file ICARD) to F06
      DO I=1,NUM_USERIN_GRIDS,4

         IF      (NUM_LEFT >=  4) THEN
            WRITE(F06,224) ICARD(I,1)  , ICARD(I,2)  , ICARD(I+1,1), ICARD(I+1,2),                                                 &
                           ICARD(I+2,1), ICARD(I+2,2), ICARD(I+3,1), ICARD(I+3,2)
            NUM_LEFT = NUM_LEFT - 4

         ELSE IF (NUM_LEFT == 3) THEN
            WRITE(F06,223) ICARD(I,1)  , ICARD(I,2)  , ICARD(I+1,1), ICARD(I+1,2),                                                 &
                           ICARD(I+2,1), ICARD(I+2,2)

         ELSE IF (NUM_LEFT == 2) THEN
            WRITE(F06,222) ICARD(I,1)  , ICARD(I,2)  , ICARD(I+1,1), ICARD(I+1,2)

         ELSE IF (NUM_LEFT == 1) THEN
            WRITE(F06,222) ICARD(I,1)  , ICARD(I,2)

         ENDIF

      ENDDO

! Write PUSERIN entry

      WRITE(F06,*)
      WRITE(ICARD(2,1),201) CUSERIN_PID                    ;   CALL LEFT_ADJ_BDFLD ( ICARD(2,1) )
      WRITE(ICARD(3,1),201) CUSERIN_IN4                    ;   CALL LEFT_ADJ_BDFLD ( ICARD(3,1) )
      WRITE(F06,241) ICARD(2,1), ICARD(3,1)

! Write message regarding CORD_MSG

      WRITE(F06,298)
      WRITE(F06,299) CORD_MSG
      WRITE(F06,*)
      WRITE(F06,399)
      WRITE(F06,101)

! **********************************************************************************************************************************
  101 FORMAT('$*******************************************************************************')

  103 FORMAT('$ B U L K   D A T A   E N T R I E S   F O R   C U S E R I N   E L E M',I8      ,/,                                   &
             '$      (to be used to define a substructure in an overall systems model)'      ,//,                                  &
             '$ The GRID and CUSERIN entries below are for the ',A,'-set from file:',/,                                            &
             '$ ',A,/,                                                                                                             &
             '$ run on ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3)

  104 FORMAT('$--1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|--10---|')

  201 FORMAT(I8)

  202 FORMAT('CUSERIN ',5A8)

  203 FORMAT(A8)

  204 FORMAT(G8.1)

  205 FORMAT('GRID    ',6A8)

  206 FORMAT('CORD2R  ',8A8)

  207 FORMAT('        ',3A8)

  211 FORMAT('SPOINT  ',A8,'THRU    ',A8,/)

  221 FORMAT(8X,2A8)

  222 FORMAT(8X,4A8)

  223 FORMAT(8X,6A8)

  224 FORMAT(8X,8A8)

  241 FORMAT('PUSERIN ',2A8,'<-"mat 1 name"-><-"mat 2 name"-><-"mat 3 name"->',/)

  298 FORMAT('$ NOTES:',/,'$ -----')

  299 FORMAT('$ ',A,' is to be replaced with the coord sys ID that is used to define the'    ,/,                                   &
             '$ basic coord sys of this USERIN elem rel to the system model basic coord     ',/,                                   &
             '$ system in the system model run'                                              ,//,                                  &
             '$ If the above grid entries are used, and are different than the corresponding',/,                                   &
             '$ grids in the system model, RBE2''s should be included to connect them to the',/,                                   &
             '$ corresponding grids in the system model.'                                    ,/,                                   &
             '$ **NOTE: If RBE2''s are NOT used, it is imperative that the grids in the     ',/,                                   &
             '$ system model that this USERIN element is connected to have the same global  ',/,                                   &
             '$ coordinate system as was used in generating this substructure.')

  399 FORMAT('$ name 1 is to be replaced with the stiffness matrix name:'                    ,/,                                   &
             '$        For CB model generation, KXX or its alias, KRRGN'                     ,/,                                   &
             '$        For statics, KGG, KAA, etc'                                           ,//,                                  &
             '$ name 2 is to be replaced with the mass matrix name (if one is input):'       ,/,                                   &
             '$        For CB model generation, MXX or its alias, MRRGN'                     ,/,                                   &
             '$        For statics, MGG, MAA, etc'                                           ,//,                                  &
             '$ name 3 is to be replaced with:'                                              ,/,                                   &
             '$        For CB model generation, RBM0 (not required)'                         ,/,                                   &
             '$        For statics, load matrix PG, PA, etc'                                 ,//,                                  &
             '$ The matrices whose names are "name i" must have been requested to be'        ,/,                                   &
             '$ written to binary files via Exec Control OUTPUT4 statement(s) in this run'   ,//,                                  &
             '$ Finally, make sure that the real numbers above have enough decimal places to',/,                                   &
             '$ accurately represent the quantities. Otherwise replace them before using them')


 9988 FORMAT(' *WARNING    : CANNOT FIND COORD SYSTEM ',I8,' IN SUBR ',A)





! **********************************************************************************************************************************

      END SUBROUTINE WRITE_USERIN_BD_CARDS


   END MODULE TEMP_FILE_WRITERS
