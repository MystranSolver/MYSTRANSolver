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

      MODULE MODEL_STUF

! This module contains many of the allocatable arrays that get allocated in subr ALLOCATE_MODEL_STUF as well as some other
! non-allocateble model variables. Some of the dimensions of the arrays described here are defined in module SCONTR

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ONEPM4, ZERO, TEN, ONE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, JCARD_LEN, MELDTS, MEMATR, MEMATC, MEPROP, METYPE,            &
                                         MPSOLID, MEFE, MEFEI, MEFER, MEWE, MEWEI, MEWER, MAX_STRESS_POINTS, MPBEAM_STATIONS

      IMPLICIT NONE

      INTEGER, PRIVATE                :: I, J

      SAVE

! **********************************************************************************************************************************
! Concentrated mass data
! ----------------------

      INTEGER(LONG), ALLOCATABLE      :: CONM2(:,:)          ! NCONM2 x 3 integer array of CONM2 data.
      INTEGER(LONG), ALLOCATABLE      :: CMASS(:,:)          ! NPMASS x 7 integer array of scalar mass connection data.
      INTEGER(LONG), ALLOCATABLE      :: PMASS(:,:)          ! NPMASS x 1 integer array of scalar mass property data.

      REAL(DOUBLE) , ALLOCATABLE      :: RCONM2(:,:)         ! NCONM2 x 10 real array of CONM2 data.
      REAL(DOUBLE) , ALLOCATABLE      :: RPMASS(:,:)         ! NPMASS x  1 real array of scalar mass data.

! Each row of CONM2 has:

!             Col. 1: CONM2 ID
!             Col. 2: GID = grid point where CONM2 is located
!             Col. 3: CID = ID of coord. sys. that CONM2 is described in

! Each row of RCONM2 has:

!           After Bulk Data has been read:
!             Col.  1: mass
!             Col.  2: x offset in coord system CID at the mass
!             Col.  3: y offset in coord system CID at the mass
!             Col.  4: z offset in coord system CID at the mass
!             Col.  5: MOI-11   in coord system CID at the mass
!             Col.  6: MOI-21   in coord system CID at the mass
!             Col.  7: MOI-22   in coord system CID at the mass
!             Col.  8: MOI-31   in coord system CID at the mass
!             Col.  9: MOI-32   in coord system CID at the mass
!             Col. 10: MOI-33   in coord system CID at the mass

!           After running subr CONM2_PROC_1:
!             Col.  1: mass
!             Col.  2: x offset in basic coords at the mass
!             Col.  3: y offset in basic coords at the mass
!             Col.  4: z offset in basic coords at the mass
!             Col.  5: MOI-11   in basic coords at the mass
!             Col.  6: MOI-21   in basic coords at the mass
!             Col.  7: MOI-22   in basic coords at the mass
!             Col.  8: MOI-31   in basic coords at the mass
!             Col.  9: MOI-32   in basic coords at the mass
!             Col. 10: MOI-33   in basic coords at the mass

!           After running subr CONM2_PROC_2:
!             Col.  1: mass
!             Col.  2: x offset in basic coords at grid GID
!             Col.  3: y offset in basic coords at grid GID
!             Col.  4: z offset in basic coords at grid GID
!             Col.  5: MOI-11   in basic coords at grid GID
!             Col.  6: MOI-21   in basic coords at grid GID
!             Col.  7: MOI-22   in basic coords at grid GID
!             Col.  8: MOI-31   in basic coords at grid GID
!             Col.  9: MOI-32   in basic coords at grid GID
!             Col. 10: MOI-33   in basic coords at grid GID

! CONM2 and RCONM2 are sorted in the order that the CONM2's were input in the Bulk Data. They are initially created
! in subroutine BD_CONM2. Array RCONM2 is modified in subroutines CONM2_PROC_1 and CONM2_PROC_2 as described above.

! **********************************************************************************************************************************
! Model total mass data (used in modal efective mass output to get percentages)
! ----------------------

      REAL(DOUBLE)                    :: MCG(6,6)            ! CG mass matrix (not incl in any other module)
      REAL(DOUBLE)                    :: MEFM_RB_MASS(6,6)   ! Rigid body mass rel to MEFMGRIF

      REAL(DOUBLE)                    :: MODEL_XCG           ! Model X-CG (from subr GPWG)
      REAL(DOUBLE)                    :: MODEL_YCG           ! Model Y-CG (from subr GPWG)
      REAL(DOUBLE)                    :: MODEL_ZCG           ! Model Z-CG (from subr GPWG)

! The following mass, inertia terms are from MCG

      REAL(DOUBLE)                    :: MODEL_MASS          ! Model total mass or weight (from subr GPWG): MCG(i,i), i=1,2,3
      REAL(DOUBLE)                    :: MODEL_IXX           ! Model moment of inertia about x (from subr GPWG): MCG(4,4)
      REAL(DOUBLE)                    :: MODEL_IYY           ! Model moment of inertia about y (from subr GPWG): MCG(5,5)
      REAL(DOUBLE)                    :: MODEL_IZZ           ! Model moment of inertia about z (from subr GPWG): MCG(6,6)

! **********************************************************************************************************************************
! Coordinate system data
! ----------------------

      INTEGER(LONG), ALLOCATABLE      :: CORD(:,:)           ! NCORD x MCORD integer array of coord system data.

      REAL(DOUBLE) , ALLOCATABLE      :: RCORD(:,:)          ! NCORD x MRCORD real array of coord system data.
      REAL(DOUBLE) , ALLOCATABLE      :: TN(:,:,:)           ! 3 x 3 x LCORD real array of coord transformation matrices
!                                                              (transforms a vector to its ref sys, RID, from its master sys, CID)

! Each row of CORD, for a CORD1R,C,S has:

!          Col. 1: Coord type
!                  11 for CORD1R
!                  12 for CORD1C
!                  13 for CORD1S
!          Col. 2: CID = Coord system ID
!          Col. 3: RID = Reference coord system ID for Grid A
!          Col. 4: RID = Reference coord system ID for Grid B
!          Col. 5: RID = Reference coord system ID for Grid C

! Each row of CORD, for a CORD2R,C,S has the same as for CORD1 above but only needs 1 reference system:

!          Col. 1: Coord type
!                  21 for CORD2R
!                  22 for CORD2C
!                  23 for CORD2S
!          Col. 2: CID = Coord system ID
!          Col. 3: RID = Reference coord system ID for this CORD2R,C,S

! Each row of RCORD has:

!          After Bulk Data has been read:

!            Col.  1: X coord of Pt A in RID coord system
!            Col.  2: Y coord of Pt A in RID coord system
!            Col.  3: Z coord of Pt A in RID coord system
!            Col.  4: X coord of Pt B in RID coord system
!            Col.  5: Y coord of Pt B in RID coord system
!            Col.  6: Z coord of Pt B in RID coord system
!            Col.  7: X coord of Pt C in RID coord system
!            Col.  8: Y coord of Pt C in RID coord system
!            Col.  9: Z coord of Pt C in RID coord system
!            Col. 10: 0.
!            Col. 11: 0.
!            Col. 12: 0.

!          After subr CORD_PROC has run:

!            Col.  1: Location of CID origin in basic coord sys X dir
!            Col.  2: Location of CID origin in basic coord sys Y dir
!            Col.  3: Location of CID origin in basic coord sys Z dir

!            Col.  4: TN(1,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
!            Col.  5: TN(1,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 1st row of TN in cols  4- 6 of RCORD
!            Col.  6: TN(1,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]

!            Col.  7: TN(2,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
!            Col.  8: TN(2,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 2nd row of TN in cols  7- 9 of RCORD
!            Col.  9: TN(2,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]

!            Col. 10: TN(3,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
!            Col. 11: TN(3,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 3rd row of TN in cols 10-12 of RCORD
!            Col. 12: TN(3,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
!            where:
!                      | TN(1,1,M)  TN(1,2,M)  TN(1,3,M) |
!              TN(M) = | TN(2,1,M)  TN(2,2,M)  TN(2,3,M) |
!                      | TN(3,1,M)  TN(3,2,M)  TN(3,3,M) |

!            is the coord system transformation which will transform a vector
!            in internal coordinate system ID = M (of coord system principal axes)
!            to a vector in basic coord system

! CORD and RCORD are sorted in the order that the CORDij Bulk Data entries were input in the Bulk Data. Subroutine
! BD_CORD initially creates arrays CORD and RCORD. Subr CORD_PROC changes the entries in RCORD as explained above.

! **********************************************************************************************************************************
! Grid sequence data
! ------------------

      INTEGER(LONG), ALLOCATABLE      :: SEQ1(:)             ! NSEQ integer array of Grid ID's on SEQGP Bulk Data entries

      REAL(DOUBLE),  ALLOCATABLE      :: SEQ2(:)             ! NSEQ real array of sequence no's on SEQGP Bulk Data entries
!                                                            SEQ2(I) is the sequence number for grid SEQ1(I)

! After Bulk Data is read, SEQ1 and SEQ2 are in the order in which they were encountered in the Bulk Data. They are
! initially created in subr BD_SEQP. After subr SEQ_PROC has run, they are in an order in which the grids in SEQ1 are
! in numerical order (and this is the order they are in when they are written to filename.L1B)

! **********************************************************************************************************************************
! Grid data
! ---------

      INTEGER(LONG), ALLOCATABLE      :: GRID(:,:)           ! Array of int data from GRID Bulk Data entries (see comments below)
      INTEGER(LONG)                   :: GRDSET3  = 0        ! Input coord system defined in field 3 of a GRDSET entries, if present
      INTEGER(LONG)                   :: GRDSET7  = 0        ! Displ coord system defined in field 7 of a GRDSET entries, if present
      INTEGER(LONG)                   :: GRDSET8  = 0        ! Permanent SPC's defined in field 8 of a GRDSET entries, if present
      INTEGER(LONG), ALLOCATABLE      :: GRID_ID (:)         ! Array of grid ID's in numerical order
      INTEGER(LONG), ALLOCATABLE      :: GRID_SEQ(:)         ! GRID_SEQ(i) is the sequence number for grid GRID_ID(i) (see below)

                                                             ! Array that shows the elements connected to each grid
      INTEGER(LONG), ALLOCATABLE      :: GRID_ELEM_CONN_ARRAY(:,:)

      INTEGER(LONG), ALLOCATABLE      :: INV_GRID_SEQ(:)     ! INV_GRID_SEQ(i) = internal grid ID that is sequenced i-th (see below)

      REAL(DOUBLE) , ALLOCATABLE      :: RGRID(:,:)          ! Array of real data from GRID Bulk Data entries (see comments below)

! Each row of GRID is for one grid point and contains:

!              (1) Grid point number    in col 1
!              (2) Input coord system   in col 2
!              (3) Global coord system  in col 3
!              (4) Permanent SPC's      in col 4
!              (5) Line break indicator in col 5 (put this many line breaks in F06 after this grid number)
!              (6) Num of comps         in col 6 (1 for SPOINT, 6 for actual grid)

!            The array is sorted in the following order:
!              (1) After the B.D. deck is read it is in GRID input order
!              (2) After subr GRID_PROC it is in grid point numerical order

! Each row of RGRID is for one grid point and contains:

!          After Bulk Data has been read: the 3 coords of the grid in the coord sys defined in col 2 of array GRID for this G.P.
!          After subr GRID_PROC has run : the 3 coords of the grid in the basic (0) coord sys of the model

! The following example explains GRID_ID(I), GRID_SEQ(I) and INV_GRID_SEQ(I). The model has 5 grid points and they are sequenced
! in the order 401, 201, 501, 301, 101. Array GRID_ID has the grids in numerical order. GRID_SEQ(I) is the sequence number for
! GRID_ID(I). INV_GRID_SEQ(I) does the following: the 4 for INV_GRID_SEQ(1) means that the 4th entry in GRID_ID (grid 401) is
! sequenced as GRID_SEQ(4) or first (as stated in the example).

!                       I  GRID_ID(I)  GRID_SEQ(I)  INV_GRID_SEQ(I)
!
!                       1     101          5              4 (i.e. the 4th entry in GRID_ID, grid 401, is sequenced 1st)
!                       2     201          2              2 ( "    "  2nd   "   "     "       "  201  "     "      2nd)
!                       3     301          4              5 ( "    "  5th   "   "     "       "  501  "     "      3rd)
!                       4     401          1              3 ( "    "  3rd   "   "     "       "  301  "     "      4th)
!                       5     501          3              1 ( "    "  1st   "   "     "       "  101  "     "      5th)

! N O T E:  A R R A Y S   G R I D   A N D   G R I D_I D   M U S T   B E   S O R T E D   T H E   S A M E   A F T E R   S U B R
!           G R I D_P R O C   H A S   R U N

! GRID_ELEM_CONN_ARRAY has NGRID rows, one for each grid i (in numerical order).
! It has a number of cols = 2 + number of elems connected to grid i. A typical array is:

!                     Table of elements connected to each grid

!         Grid    Num elems     ID's of elements connected to this grid -->

!         1011            3       11     1121     1141
!         1012            2       11       12
!         1013            3       12     1323     3143
!         1021            5       11       21       22     1121     2131

! **********************************************************************************************************************************
! Grid data, con't
! ----------------

      INTEGER(LONG), ALLOCATABLE      :: MPC_IND_GRIDS(:)    ! Array of actual grid numbers that are referenced as independent on
!                                                              MPC's and rigid elements. KGG singularity processor needs to skip
!                                                              these grids since they may have no stiffness in KGG but do in KNN.

! **********************************************************************************************************************************
! SNORM data (grid point surface normals)
! ----------

      INTEGER(LONG), ALLOCATABLE      :: SNORM(:,:)          ! Array SNORM GID and CID.
                                                             ! One row for each SNORM entry in the file.
                                                             ! Only used in LINK0 to read the entries.
      REAL(DOUBLE), ALLOCATABLE       :: RSNORM(:,:)         ! Array SNORM vectors in the input (CID) coordinate systems.
                                                             ! One row for each SNORM entry in the file, corresponding to SNORM.
                                                             ! Only used in LINK0 to read the entries.
      REAL(DOUBLE), ALLOCATABLE       :: GRID_SNORM(:,:)     ! NGRID x MGRID_SNORM array normalized surface normal vectors in
                                                             ! the basic coordinate system indexed by internal grid point number
                                                             ! (BGRID), the same as the GRID and RGRID arrays. Where all 3
                                                             ! components are exactly zero, there is no normal for that grid point.

! **********************************************************************************************************************************
! Force and load set ID's
! -----------------------

      INTEGER(LONG), ALLOCATABLE      :: LOAD_SIDS(:,:)      ! Array of integer data from Bulk Data LOAD entries.
      INTEGER(LONG), ALLOCATABLE      :: FORMOM_SIDS(:)      ! NFORCE  x 1 array of all set ID's from FORCE and MOMENT B.D. entries
      INTEGER(LONG), ALLOCATABLE      :: GRAV_SIDS(:)        ! NGRAV   x 1 array of all set ID's from GRAV   B.D. entries
      INTEGER(LONG), ALLOCATABLE      :: PRESS_SIDS(:)       ! NPLOAD  x 1 array of all set ID's from PLOAD  B.D. entries
      INTEGER(LONG), ALLOCATABLE      :: RFORCE_SIDS(:)      ! NRFORCE x 1 array of all set ID's from RFORCE B.D. entries
      INTEGER(LONG), ALLOCATABLE      :: SLOAD_SIDS(:)       ! NSLOAD  x 1 array of all set ID's from SLOAD  B.D. entries

      REAL(DOUBLE) , ALLOCATABLE      :: LOAD_FACS(:,:)      ! Array of real data from Bulk Data LOAD entries

! Arrays LOAD_SIDS and LOAD_FACS are built in subr BD_LOAD from data on B.D LOAD entries. Each B.D. LOAD entry adds a new row to
! arrays LOAD_SIDS and LOAD_FACS.

!  LOAD_SIDS: The cols of LOAD_SIDS contain:
!                    (1) Col 1: set ID of the LOAD B.D. entry
!                    (2) Col 2: set ID of the 1st load being combined
!                     .    .     .    .   .  .  .   .    .    .      .
!
!                    (i) Col i: load set ID of the i-1 load being combined

!  LOAD_FACS: The cols of LOAD_FACS contain:
!                    (1) Col 1: The overall scale factor for the LOAD combination
!                    (2) Col 2: The scale factor for the 1st load being combined
!                    (2) Col 3: The scale factor for the 2nd load being combined
!                     .    .     .    .   .  .  .   .    .    .      .
!
!                    (i) Col i: The scale factor for the i-1 load being combined

! **********************************************************************************************************************************
! SPC, MPC set ID's
! -----------------

      INTEGER(LONG), ALLOCATABLE      :: MPC_SIDS(:)         ! NMPC  x 1 array of set ID's from all MPC  Bulk data entries
      INTEGER(LONG), ALLOCATABLE      :: MPCSIDS(:)          ! NUM_MPCSIDS  x 1 array of set ID's from MPC, MPCADD for 1 run
      INTEGER(LONG), ALLOCATABLE      :: MPCADD_SIDS(:,:)    ! Array of integer data from Bulk Data MPCADD entries.
      INTEGER(LONG), ALLOCATABLE      :: SPC_SIDS(:)         ! NSPC  x 1 array of set ID's from all SPC  Bulk data entries
      INTEGER(LONG), ALLOCATABLE      :: SPC1_SIDS(:)        ! NSPC1 x 1 array of set ID's from all SPC1 Bulk data entries
      INTEGER(LONG), ALLOCATABLE      :: SPCSIDS(:)          ! NUM_SPCSIDS  x 1 array of set ID's from SPC, SPCADD for 1 run
      INTEGER(LONG), ALLOCATABLE      :: SPCADD_SIDS(:,:)    ! Array of integer data from Bulk Data SPCADD entries.

!  SPCADD_SIDS: Each SPCADD entry adds a new row to array SPCADD_SIDS (below). The cols of SPCADD_SIDS are:
!                    (1) Col 1: SPCADD set ID of the SPCADD B.D. entry
!                    (2) Col 2: SPC/SCP1 set ID of the 1st SPC requested on the SPCADD entry
!                    (2) Col 3: SPC/SCP1 set ID of the 2nd SPC requested on the SPCADD entry
!                     .    .     .    .   .  .  .   .    .    .      .
!
!                    (i) Col i: SPC/SCP1 set ID of the i-1 SPC requested on the SPCADD entry

!  MPCADD_SIDS: Each MPCADD entry adds a new row to array MPCADD_SIDS (below). The cols of MPCADD_SIDS are:
!                    (1) Col 1: MPCADD set ID of the MPCADD B.D. entry
!                    (2) Col 2: MPC set ID of the 1st MPC requested on the MPCADD entry
!                    (2) Col 3: MPC set ID of the 2nd MPC requested on the MPCADD entry
!                     .    .     .    .   .  .  .   .    .    .      .
!
!                    (i) Col i: MPC set ID of the i-1 MPC requested on the MPCADD entry

! **********************************************************************************************************************************
! Element and grid temperature data
! ---------------------------------

      CHARACTER( 1*BYTE), ALLOCATABLE :: CETEMP(:,:)         ! NELE  x NTSUB array of indicators of whether a elem has a
!                                                              temperature defined for a subcase
      CHARACTER( 1*BYTE), ALLOCATABLE :: CGTEMP(:,:)         ! NGRID x NTSUB array of indicators of whether a grid has a
!                                                              temperature defined for a subcase
      CHARACTER( 1*BYTE), PARAMETER   :: CETEMP_ERR = '*'    ! Value used to indicate an elem has no temperature defined for a S/C.
!                                                              Also, value used in initializing array CETEMP
      CHARACTER( 1*BYTE), PARAMETER   :: CGTEMP_ERR = '*'    ! Value used to indicate a grid has no temperature defined for a S/C.
!                                                              Also, value used in initializing array CGTEMP

      INTEGER(LONG)    , ALLOCATABLE  :: TPNT(:,:)           ! NELE x 1 array of pointers to where in TDATA the temperature data
!                                                              begins for an element

      REAL(DOUBLE)     , ALLOCATABLE  :: ETEMP(:,:)          ! NELE  x NSUB array of element temperatures
      REAL(DOUBLE)     , ALLOCATABLE  :: GTEMP(:,:)          ! NGRID x NSUB array of grid point temperatures from TEMP, TEMPD
      REAL(DOUBLE)                    :: ETEMP_INIT = ZERO   ! Value used in initializing array ETEMP
      REAL(DOUBLE)                    :: GTEMP_INIT = ZERO   ! Value used in initializing array GTEMP
      REAL(DOUBLE)     , ALLOCATABLE  :: TDATA(:)            ! NTDAT x 1 array of element temperature data from TEMPRB and TEMPP1

! **********************************************************************************************************************************
! Element pressure data
! ----------------------

      CHARACTER( 1*BYTE), ALLOCATABLE :: PTYPE(:)            ! NELE x 1 array of how many pressures there are in PDATA for an elem.
!                                                              ('6' for PLOAD1, '1' for PLOAD2, and '3' or '4' for PLOAD4)

      INTEGER(LONG), ALLOCATABLE      :: PPNT(:,:)           ! NELE x 1 array of pointers to where in PDATA the elem pressure data
!                                                              begins for an element

      INTEGER(LONG), ALLOCATABLE      :: PLOAD4_3D_DATA(:,:) ! NPLOAD4_3D x 1 array of elem ID and 3 or 4 grids for the face on
!                                                              a solid element where a pressure is defined (via PLOAD4 B.D. entry)
!                                                              Col 1 is the actual element ID
!                                                              Col 2 is the internal element ID
!                                                              Col 3 is the internal subcase number for the PLOAD4 entry
!                                                              Col 4 is G1  (see PLOAD4 B.D. entry description)
!                                                              Col 5 is G34 (see PLOAD4 B.D. entry description)

      REAL(DOUBLE) , ALLOCATABLE      :: PDATA(:)            ! NPDAT x 1 array of elem press data from PLOADi B.D. entries

! **********************************************************************************************************************************
! Subcase items
! -------------

      CHARACTER(LEN=CC_ENTRY_LEN)                           &! Character array of Case Control TITLE  info
                        , ALLOCATABLE :: TITLE(:)

      CHARACTER(LEN=CC_ENTRY_LEN)                           &! Character array of Case Control STITLE info
                        , ALLOCATABLE :: STITLE(:)

      CHARACTER(LEN=CC_ENTRY_LEN)                           &! Character array of Case Control LABEL  info
                        , ALLOCATABLE :: LABEL(:)


      CHARACTER( 1*BYTE), ALLOCATABLE :: ALL_SETS_ARRAY(:)   ! Array containing all C.C. SET entries (complete entry incl "SET id=")
      CHARACTER( 1*BYTE), ALLOCATABLE :: ONE_SET_ARRAY(:)    ! Array containing one SET (after "=" sign) from ALL_SETS_ARRAY
      CHARACTER( 1*BYTE)              :: MEFFMASS_CALC = 'N' ! If 'Y', calc modal effective masses
      CHARACTER( 1*BYTE)              :: MPFACTOR_CALC = 'N' ! If 'Y', calc modal participation factors

      INTEGER(LONG)                   :: CCELDT(0:MELDTS-1) = (/(0,I=0,MELDTS-1)/)
                                                             ! See description below

      INTEGER(LONG)                   :: ANY_ACCE_OUTPUT     ! > 0 if requests for output of accels in a any S/C
      INTEGER(LONG)                   :: ANY_DISP_OUTPUT     ! > 0 if requests for output of displs in a any S/C
      INTEGER(LONG)                   :: ANY_OLOA_OUTPUT     ! > 0 if requests for output of applied loads in a any S/C
      INTEGER(LONG)                   :: ANY_SPCF_OUTPUT     ! > 0 if requests for output of SPC forces in a any S/C
      INTEGER(LONG)                   :: ANY_MPCF_OUTPUT     ! > 0 if requests for output of MPC forces in a any S/C
      INTEGER(LONG)                   :: ANY_GPFO_OUTPUT     ! > 0 if requests for output of G.P. force balance in a any S/C
      INTEGER(LONG)                   :: ANY_ELFN_OUTPUT     ! > 0 if requests for output of elem node forces in a any S/C
      INTEGER(LONG)                   :: ANY_ELFE_OUTPUT     ! > 0 if requests for output of elem engr forces in a any S/C
      INTEGER(LONG)                   :: ANY_STRE_OUTPUT     ! > 0 if requests for output of elem stresses in a any S/C
      INTEGER(LONG)                   :: ANY_STRN_OUTPUT     ! > 0 if requests for output of elem strains in a any S/C

      INTEGER(LONG)     , ALLOCATABLE :: OGROUT(:)           ! See description below
      INTEGER(LONG)     , ALLOCATABLE :: GROUT(:,:)          ! See description below

      INTEGER(LONG)     , ALLOCATABLE :: ELOUT(:,:)          ! See description below
      INTEGER(LONG)     , ALLOCATABLE :: OELOUT(:)           ! See description below

      INTEGER(LONG)                   :: OELDT               ! See description below
      INTEGER(LONG)     , ALLOCATABLE :: ELDT(:)             ! See description below

      INTEGER(LONG)     , ALLOCATABLE :: SETS_IDS(:)         ! Array of set ID's from all C.C SET entries
      INTEGER(LONG)     , ALLOCATABLE :: SCNUM(:)            ! Array of subcase numbers
      INTEGER(LONG)     , ALLOCATABLE :: SUBLOD(:,:)         ! LSUB x 2 array. Row i (for internal S/C i) gives ID's of C.C. LOAD
!                                                            (in col 1) and C.C. TEMP (in col 2) found for that S/C in C.C.

      INTEGER(LONG)                   :: CC_EIGR_SID         = 0
                                                             ! Eigenvalue extraction method requested by the Case Control METH entry

      INTEGER(LONG)    , ALLOCATABLE  :: SC_ACCE(:)          ! C.C. ACCEL request info.         Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_DISP(:)          ! C.C. DISP request info.          Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_ELFN(:)          ! C.C. ELFORCE(NODE) request info. Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_ELFE(:)          ! C.C. ELFORCE(engr) request info. Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_GPFO(:)          ! C.C. GPFORCE request info.       Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_MPCF(:)          ! C.C. MPCFORCE request info.      Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_OLOA(:)          ! C.C. OLOAD request info.         Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_SPCF(:)          ! C.C. SPCFORCE request info.      Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_STRE(:)          ! C.C. STRESS request info.        Indicator for "ALL", "NONE", "SETID"
      INTEGER(LONG)    , ALLOCATABLE  :: SC_STRN(:)          ! C.C. STRAIN request info.        Indicator for "ALL", "NONE", "SETID"

      INTEGER(LONG)    , ALLOCATABLE  :: MPCSETS(:)          ! Set ID's of MPC sets requested in Case Control
      INTEGER(LONG)    , ALLOCATABLE  :: SPCSETS(:)          ! Set ID's of SPC sets requested in Case Control

      INTEGER(LONG)                   :: MPCSET = 0          ! Set ID of MPC set used
      INTEGER(LONG)                   :: SPCSET = 0          ! Set ID of SPC set used

! CCELDT(1-16) is an integer array containing SET ID's (-1 for 'ALL', 0 for 'NONE', positive integer no. for SETID) of the 16
! possible outputs requested by the Case Control ELDATA command. All output is in the element local coordinate system (the one
! in which the elements are developed):

!     CCELDT( 0) is CC requests for print to BUGFIL of elem geometric data
!     CCELDT( 1) is CC requests for print to BUGFIL of elem property and material info
!     CCELDT( 2) is CC requests for print to BUGFIL of elem thermal and pressure matrices    : PTE, PPE
!     CCELDT( 3) is CC requests for print to BUGFIL of elem mass matrix                      : ME
!     CCELDT( 4) is CC requests for print to BUGFIL of elem stiffness matrix                 : KE
!     CCELDT( 5) is CC requests for print to BUGFIL of elem stress & strain recovery matrices: SEi, STEi, BEi
!     CCELDT( 6) is CC requests for print to BUGFIL of elem displacement and load matrices   : UEL, PEL (all subcases)
!     CCELDT( 7) is CC requests for print to BUGFIL of elem shape fcns and Jacobian matrices
!     CCELDT( 8) is CC requests for print to BUGFIL of elem strain-displacement matrices
!     CCELDT( 9) is CC requests for print to BUGFIL of elem checks on strain-displ matrices for RB motion & constant strain
!     CCELDT(10) is CC requests for write to F21FIL unformatted file of element              : PTE, PPE
!     CCELDT(11) is CC requests for write to F22FIL unformatted file of element              : ME
!     CCELDT(12) is CC requests for write to F23FIL unformatted file of element              : KE
!     CCELDT(13) is CC requests for write to F24FIL unformatted file of element              : SEi, STEi, BEi
!     CCELDT(14) is CC requests for write to F25FIL unformatted file of element              : UEL, PEL (all subcases)
!     CCELDT(15) is CC requests currently not used

!  The SC_xxxx arrays all use the same convention to indicate "ALL", "NONE" or "SETID":

!     (1) "ALL"   is indicated  by a value of -1
!     (2) "NONE"  is indicated  by a value of  0
!     (3) "SETID" is indicated  by a positive integer ID of a C.C. SET entry

!  OGROUT: One row for each S/C. Bits indicate if C.C. requests were made for G.P. related outputs for that S/C.

!     Bit  0 in a row of OGROUT indicates some Case Control DISP  output requests were found for that S/C
!     Bit  1 in a row of OGROUT indicates some Case Control OLOAD output requests were found for that S/C
!     Bit  2 in a row of OGROUT indicates some Case Control SPCF  output requests were found for that S/C
!     Bit  3 in a row of OGROUT indicates some Case Control GPFO  output requests were found for that S/C

!     Array OGROUT is set in subr SCPRO based on the values read from SC_DISP, SC_OLOA, SC_SPCF, SC_GPFO

!  GROUT : Each row is for 1 G.P. and each col is for 1 S/C. If, for a subcase:

!     SC_DISP = -1   , bit 0 in OGROUT and in all rows of GROUT are set
!     SC_DISP =  0   , bit 0 in OGROUT and in all rows of GROUT are left cleared
!     SC_DISP = setid, bit 0 in OGROUT and in the rows of GROUT for the setid grids are set

!     SC_OLOA = -1   , bit 1 in OGROUT and in all rows of GROUT are set
!     SC_OLOA =  0   , bit 1 in OGROUT and in all rows of GROUT are left cleared
!     SC_OLOA = setid, bit 1 in OGROUT and in the rows of GROUT for the setid grids are set

!     SC_SPCF = -1   , bit 2 in OGROUT and in all rows of GROUT are set
!     SC_SPCF =  0   , bit 2 in OGROUT and in all rows of GROUT are left cleared
!     SC_SPCF = setid, bit 2 in OGROUT and in the rows of GROUT for the setid grids are set

!     SC_GPFO = -1   , bit 3 in OGROUT and in all rows of GROUT are set
!     SC_GPFO =  0   , bit 3 in OGROUT and in all rows of GROUT are left cleared
!     SC_GPFO = setid, bit 3 in OGROUT and in the rows of GROUT for the setid grids are set

!     Array GROUT is set in subr SCPRO based on the values read from SC_DISP, SC_OLOA, SC_SPCF, SC_GPFO

!  OELOUT: One row for each S/C. Bits indicate if C.C. requests were made for elem related outputs for that S/C.

!     Bit  0 in a row of OELOUT indicates some Case Control ELFORCE(NODE) output requests were found
!     Bit  1 in a row of OELOUT indicates some Case Control ELFORCE(ENGR) output requests were found
!     Bit  2 in a row of OELOUT indicates some Case Control STRESS output requests were found
!     Array OELOUT is set in subr SCPRO based on the values read from SC_ELFN, SC_ELFE, SC_STRE, SC_STRN

!  ELOUT : Each row is for 1 elem and each col is for 1 S/C. If, for a subcase:

!     SC_ELFN = -1   , bit 0 in OELOUT and in all rows of ELOUT are set
!     SC_ELFN =  0   , bit 0 in OELOUT and in all rows of ELOUT are left cleared
!     SC_ELFN = setid, bit 0 in OELOUT and in the rows of ELOUT for the setid ELids are set

!     SC_ELFE = -1   , bit 1 in OELOUT and in all rows of ELOUT are set
!     SC_ELFE =  0   , bit 1 in OELOUT and in all rows of ELOUT are left cleared
!     SC_ELFE = setid, bit 1 in OELOUT and in the rows of ELOUT for the setid ELids are set

!     SC_STRE = -1   , bit 2 in OELOUT and in all rows of ELOUT are set
!     SC_STRE =  0   , bit 2 in OELOUT and in all rows of ELOUT are left cleared
!     SC_STRE = setid, bit 2 in OELOUT and in the rows of ELOUT for the setid grids are set

!     SC_STRN = -1   , bit 3 in OELOUT and in all rows of ELOUT are set
!     SC_STRN =  0   , bit 3 in OELOUT and in all rows of ELOUT are left cleared
!     SC_STRN = setid, bit 3 in OELOUT and in the rows of ELOUT for the setid grids are set

!     Array OELOUT is set in subr SCPRO based on the values read from SC_ELFN, SC_ELFE, SC_STRE, SC_STRN

!  OELDT tells whether any element has ELDATA requests:

!        Bit  0 in OELDT indicates some requests for elem geometric data
!        Bit  1 in OELDT indicates some requests for elem property and material info
!        Bit  2 in OELDT indicates some requests for elem thermal and pressure matrices       : PTE, PPE
!        Bit  3 in OELDT indicates some requests for elem mass matrix                         : ME
!        Bit  4 in OELDT indicates some requests for elem stiffness matrix                    : KE
!        Bit  5 in OELDT indicates some requests for elem stress & strain recovery matrices   : SEi, STEi, BEi
!        Bit  6 in OELDT indicates some requests for elem displacement and load matrices      : UEL, PEL (all subcases)
!        Bit  7 in OELDT indicates some requests for elem shape fcns and Jacobian matrices
!        Bit  8 in OELDT indicates some requests for elem strain-displacement matrices
!        Bit  9 in OELDT indicates some requests for elem checks on strain-displ matrices for RB motion & constant strain
!        Bit 10 in OELDT indicates some requests for unformatted file of element              : PTE, PPE
!        Bit 11 in OELDT indicates some requests for unformatted file of element              : ME
!        Bit 12 in OELDT indicates some requests for unformatted file of element              : KE
!        Bit 13 in OELDT indicates some requests for unformatted file of element              : SEi, STEi, BEi
!        Bit 14 in OELDT indicates some requests for unformatted file of element              : UEL, PEL (all subcases)
!        Bit 15 in OELDT is currently not used

!     Array OELDT is set in subr SCPRO based on the values read from CCELDT

!  ELDT: Array of 16 bits for each element to indicate whether any Case Control ELDATA requests were made for that element
!  ELDT(j) tells whether internal element j has requests:

!        Bit  0 in ELDT is for print to BUGFIL of elem geometric data
!        Bit  1 in ELDT is for print to BUGFIL of elem property and material info
!        Bit  2 in ELDT is for print to BUGFIL of elem thermal and pressure matrices    : PTE, PPE
!        Bit  3 in ELDT is for print to BUGFIL of elem mass matrix                      : ME
!        Bit  4 in ELDT is for print to BUGFIL of elem stiffness matrix                 : KE
!        Bit  5 in ELDT is for print to BUGFIL of elem stress & strain recovery matrices: SEi, STEi, BEi
!        Bit  6 in ELDT is for print to BUGFIL of elem displacement and load matrices   : UEL, PEL (all subcases)
!        Bit  7 in ELDT is for print to BUGFIL of elem shape fcns and Jacobian matrices
!        Bit  8 in ELDT is for print to BUGFIL of elem strain-displacement matrices
!        Bit  9 in ELDT is for print to BUGFIL of elem checks on strain-displ matrices for RB motion & constant strain
!        Bit 10 in ELDT is for write to F21FIL unformatted file of element              : PTE, PPE
!        Bit 11 in ELDT is for write to F22FIL unformatted file of element              : ME
!        Bit 12 in ELDT is for write to F23FIL unformatted file of element              : KE
!        Bit 13 in ELDT is for write to F24FIL unformatted file of element              : SEi, STEi, BEi
!        Bit 14 in ELDT is for write to F25FIL unformatted file of element              : UEL, PEL (all subcases)
!        Bit 15 in ELDT is currently not used
!  Array ELDT is set in subr SCPRO based on the values read from CCELDT

!  SUBLOD: Array with one row for each DS/C and 2 cols. Row i (for internal S/C i) gives the ID's of C.C. LOAD (col 1) and
!         C.C. TEMP (col 2) found for that S/C in C.C.

! **********************************************************************************************************************************
! Full system load array
! ----------------------

      REAL(DOUBLE) , ALLOCATABLE      :: SYS_LOAD(:,:)       ! Array of applied (incl thermal) loads on all G-set DOF's and subcases

! **********************************************************************************************************************************
! Element property and material data
! ----------------------------------

      INTEGER(LONG), ALLOCATABLE      :: MATL   (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PBAR   (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PBEAM  (:,:)         ! See description below
! --- CBEAM_standard begin --- !
      INTEGER(LONG), ALLOCATABLE      :: PBEAM_NSTATIONS(:)   ! Number of stored station x/L values for each PBEAM property
! --- CBEAM_standard end --- !
      INTEGER(LONG), ALLOCATABLE      :: PBUSH  (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PCOMP  (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PELAS  (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PROD   (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PSHEAR (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PSHEL  (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PSOLID (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PUSER1 (:,:)         ! See description below
      INTEGER(LONG), ALLOCATABLE      :: PUSERIN(:,:)         ! See description below

      REAL(DOUBLE) , ALLOCATABLE      :: RMATL  (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPBAR  (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPBEAM (:,:)         ! See description below
! --- CBEAM_standard begin --- !
      REAL(DOUBLE) , ALLOCATABLE      :: PBEAM_XL(:,:)        ! Stored station x/L values for each PBEAM property
      REAL(DOUBLE) , ALLOCATABLE      :: PBEAM_RPROPS(:,:,:)  ! Stored [A,I1,I2,I12,J,NSM] for each PBEAM station
! --- CBEAM_standard end --- !
      REAL(DOUBLE) , ALLOCATABLE      :: RPBUSH (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPCOMP (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPELAS (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPROD  (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPSHEAR(:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPSHEL (:,:)         ! See description below
      REAL(DOUBLE) , ALLOCATABLE      :: RPUSER1(:,:)         ! See description below

!  PBAR   = Array of integer data from PBAR  Bulk Data entries. Each row is for one PBAR  entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID

!  RPBAR  = Array of real data from PBAR  Bulk Data entries. Each row is for one PBAR  entry read in B.D. and contains:
!             ( 1) Col  1: Cross sectional area, A              , (parent        entry, field 4)
!             ( 2) Col  2: Area moment of inertia 1, I1         , (parent        entry, field 5)
!             ( 3) Col  3: Area moment of inertia 2, I2         , (parent        entry, field 6)
!             ( 4) Col  4: Torsion constant, J                  , (parent        entry, field 7)
!             ( 5) Col  5: Non structural mass, NSM             , (parent        entry, field 8)
!             ( 6) Col  6: Stress coefficient                   , (optional  2nd entry, field 2)
!             ( 7) Col  7: Stress coefficient                   , (optional  2nd entry, field 3)
!             ( 8) Col  8: Stress coefficient                   , (optional  2nd entry, field 4)
!             ( 9) Col  9: Stress coefficient                   , (optional  2nd entry, field 5)
!             (10) Col 10: Stress coefficient                   , (optional  2nd entry, field 6)
!             (11) Col 11: Stress coefficient                   , (optional  2nd entry, field 7)
!             (12) Col 12: Stress coefficient                   , (optional  2nd entry, field 8)
!             (13) Col 13: Stress coefficient                   , (optional  2nd entry, field 9)
!             (14) Col 14: Shear factor for plane 1, K1         , (optional  3rd entry, field 2)
!             (15) Col 15: Shear factor for plane 2, K2         , (optional  3rd entry, field 3)
!             (16) Col 16: Area MOI, I12                        , (optional  3rd entry, field 4)
!             (17) Col 17: Torsional stress recovery coefficient, (optional  3rd entry, field 5)

!  PBEAM  = Array of integer data from PBEAM Bulk Data entries. Each row is for one PBEAM entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID
!
!  PBEAM_NSTATIONS = Number of station x/L values stored for each PBEAM property. For NX-oriented beam work this includes:
!                    station 1 = 0.0 at end A, followed by each continuation station found on the PBEAM entry.

!  RPBEAM = Array of real data from PBEAM Bulk Data entries. Each row is for one PBEAM entry read in B.D. and contains:
!             ( 1) Col  1: Cross sectional area, A       : end A          , (parent        entry, field 4)
!             ( 2) Col  2: Area moment of inertia 1, I1  :   "            , (parent        entry, field 5)
!             ( 3) Col  3: Area moment of inertia 2, I2  :   "            , (parent        entry, field 6)
!             ( 4) Col  4: Area MOI, I12                 :   "            , (parent        entry, field 7)
!             ( 5) Col  5: Torsion constant, J           :   "            , (parent        entry, field 8)
!             ( 6) Col  6: Non structural mass, NSM      :   "            , (parent        entry, field 9)
!             ( 7) Col  7: Stress coefficient            :   "            , (mandatory 2nd entry, field 2) (since 3rd is mandatory)
!             ( 8) Col  8: Stress coefficient            :   "            , (mandatory 2nd entry, field 3)
!             ( 9) Col  9: Stress coefficient            :   "            , (mandatory 2nd entry, field 4)
!             (10) Col 10: Stress coefficient            :   "            , (mandatory 2nd entry, field 5)
!             (11) Col 11: Stress coefficient            :   "            , (mandatory 2nd entry, field 6)
!             (12) Col 12: Stress coefficient            :   "            , (mandatory 2nd entry, field 7)
!             (13) Col 13: Stress coefficient            :   "            , (mandatory 2nd entry, field 8)
!             (14) Col 14: Stress coefficient            :   "            , (mandatory 2nd entry, field 9)
!             (15) Col  1: Cross sectional area, A       : end B          , (mandatory 3rd entry, field 4)
!             (16) Col  2: Area moment of inertia 1, I1  :   "            , (mandatory 3rd entry, field 5)
!             (17) Col  3: Area moment of inertia 2, I2  :   "            , (mandatory 3rd entry, field 6)
!             (18) Col  4: Area MOI, I12                 :   "            , (mandatory 3rd entry, field 7)
!             (19) Col  5: Torsion constant, J           :   "            , (mandatory 3rd entry, field 8)
!             (20) Col  6: Non structural mass, NSM      :   "            , (mandatory 3rd entry, field 9)
!             (21) Col  7: Stress coefficient            :   "            , (optional  4th entry, field 2)
!             (22) Col  8: Stress coefficient            :   "            , (optional  4th entry, field 3)
!             (23) Col  9: Stress coefficient            :   "            , (optional  4th entry, field 4)
!             (24) Col 10: Stress coefficient            :   "            , (optional  4th entry, field 5)
!             (25) Col 11: Stress coefficient            :   "            , (optional  4th entry, field 6)
!             (26) Col 12: Stress coefficient            :   "            , (optional  4th entry, field 7)
!             (27) Col 13: Stress coefficient            :   "            , (optional  4th entry, field 8)
!             (28) Col 14: Stress coefficient            :   "            , (optional  4th entry, field 9)
!             (29) Col 29: Shear factor for plane 1, K1                   , (optional  5th entry, field 2)
!             (30) Col 30: Shear factor for plane 2, K2                   , (optional  5th entry, field 3)
!             (31) Col 31: Shear relief coeff due to taper for plane 1, S1, (optional  5th entry, field 4)
!             (32) Col 32: Shear relief coeff due to taper for plane 2, S2, (optional  5th entry, field 5)
!             (33) Col 33: NSM MOI/length about NSM C.G. at end A, NSI(A) , (optional  5th entry, field 6)
!             (34) Col 34: NSM MOI/length about NSM C.G. at end B, NSI(B) , (optional  5th entry, field 7)
!             (35) Col 35: Warping coefficient for end A, CW(A)           , (optional  5th entry, field 8)
!             (36) Col 36: Warping coefficient for end B, CW(B)           , (optional  5th entry, field 9)
!             (37) Col 37: y coord of C.G. of NSM at end A, M1(A)         , (optional  6th entry, field 2)
!             (38) Col 38: z coord of C.G. of NSM at end A, M2(A)         , (optional  6th entry, field 3)
!             (39) Col 39: y coord of C.G. of NSM at end B, M1(B)         , (optional  6th entry, field 4)
!             (40) Col 40: z coord of C.G. of NSM at end B, M2(B)         , (optional  6th entry, field 5)
!             (41) Col 41: y coord of neutral axis for end A, N1(A)       , (optional  6th entry, field 6)
!             (42) Col 42: z coord of neutral axis for end A, N2(A)       , (optional  6th entry, field 7)
!             (43) Col 43: y coord of neutral axis for end B, N1(B)       , (optional  6th entry, field 8)
!             (44) Col 44: z coord of neutral axis for end B, N2(B)       , (optional  6th entry, field 9)
!
!  PBEAM_XL = Stored x/L station values for each PBEAM property. Phase-1 beam redevelopment stores the continuation chain explicitly
!             so that CBEAM output can later be made station-aware independently of legacy CBAR end-only output.
!
!  PBEAM_RPROPS = Stored real section properties for each PBEAM station. The 3rd index stores:
!             (1) A
!             (2) I1
!             (3) I2
!             (4) I12
!             (5) J
!             (6) NSM

!  PBUSH  = Array of integer data from PBUSH Bulk Data entries
!             ( 1) Col  1: PID          Prop ID

!  RPBUSH = Array of real data from PBUSH Bulk Data entries
!             ( 1) Col   : K1           Stiffness in element direction
!             ( 2) Col   : K2           Stiffness in element direction
!             ( 3) Col   : K3           Stiffness in element direction
!             ( 4) Col   : K4           Stiffness in element direction
!             ( 5) Col   : K5           Stiffness in element direction
!             ( 6) Col   : K6           Stiffness in element direction
!             ( 7) Col   : B1           Viscous damping coeff in element direction
!             ( 8) Col   : B2           Viscous damping coeff in element direction
!             ( 9) Col   : B3           Viscous damping coeff in element direction
!             (10) Col   : B4           Viscous damping coeff in element direction
!             (11) Col   : B5           Viscous damping coeff in element direction
!             (12) Col   : B6           Viscous damping coeff in element direction
!             (13) Col   : GE           Structural damping coefficient
!             (14) Col   : RCV1         Stress recovery coeff for translation
!             (15) Col   : RCV2         Stress recovery coeff for rotation
!             (16) Col   : RCV3         Strain recovery coeff for translation
!             (17) Col   : RCV4         Strain recovery coeff for rotation

!  PCOMP  = Array of integer data from PCOMP, PCOMP1 Bulk Data entry
!             ( 1) Col  1: PID          Prop ID
!             ( 2) Col  2: PID_TYPE     Type of PCOMP (0 for PCOMP, 1 for PCOMP1)
!             ( 3) Col  3: FT           Failure theory (1=HILL, 2=HOFF, 3=TSAI, 4=STRN)
!             ( 4) Col  4: LAM          Symm lamination opt
!             ( 5) Col  5: PCOMP_PLIES  Number of plies for this PCOMP
!             ( 6) Col  6:              Indicator of whether this PCOMP's equiv PSHELL, MAT2's have been written to F06
!             ( 7) Col  7: MIDi         Material ID of ply i.
!             ( 8) Col  8: SOUTi        Stress or strain output request (1=YES or 0=NO) for ply i
!             ( the last 2 are repeated for each ply)

!  RPCOMP = Array of real data from PCOMP Bulk Data entries
!             ( 1) Col  1: Z0      Dist from ref plane to bottom surface (default = -0.5 times layer thickness)
!             ( 2) Col  2: NSM     Non structural mass per unit area
!             ( 3) Col  3: SB      Allowable interlaminar shear stress. Required if FT is specified
!             ( 4) Col  4: TREF    Ref temperature
!             ( 5) Col  5: GE      Damping coeff
!             ( 6) Col  6: TT      Total plate thickness
!             ( 7) Col  7: Ti      Thickness of ply i (real or blank but T1 -1st ply - must be specified)
!             ( 7) Col  8: THETAi  Orientation angle of longitudinal dir of ply i wrt material axis for the composite elem
!             ( 7) Col  9: Zi      z coord (+/-) from ref plane to center of ply
!             ( the last 3 are repeated for each ply)

!  PELAS  = Array of integer data from PELAS Bulk Data entries. Each row is for one PELAS entry read in B.D. and contains:
!             ( 1) Col  1: Property ID

!  RPELAS = Array of real data from PELAS Bulk Data entries. Each row is for one PELAS entry read in B.D. and contains:
!             ( 1) Col  1: Spring rate, K
!             ( 2) Col  2: Damping coefficient, GE
!             ( 3) Col  3: Stress recovery coefficient, S

!  PROD   = Array of integer data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID

!  RPROD  = Array of real data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
!             ( 1) Col  1: Cross sectional area, A
!             ( 2) Col  2: Torsion constant, J
!             ( 3) Col  3: Torsion stress recov. coeff, C
!             ( 4) Col  4: Non structural mass, NSM

!  PSHEAR = Array of integer data from PSHEL Bulk Data entries. Each row is for one PSHELL entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID1

!  RPSHEAR= Array of real data from PSHEL Bulk Data entries. Each row is for one PSHEL entry read in B.D. and contains:
!             ( 1) Col  1: Membrane thickness, TM
!             ( 4) Col  4: Non structural mass, NSM
!             ( 5) Col  5: F1, effectiveness factor for extensional stiffness along edges 1-2 and 3-4
!             ( 6) Col  6: F2, effectiveness factor for extensional stiffness along edges 2-3 and 1-4

!  PSHEL  = Array of integer data from PSHEL Bulk Data entries. Each row is for one PSHELL entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID 1
!             ( 3) Col  3: Material ID 2
!             ( 4) Col  4: Material ID 3
!             ( 5) Col  5: Material ID 4

!  RPSHEL = Array of real data from PSHEL Bulk Data entries. Each row is for one PSHEL entry read in B.D. and contains:
!             ( 1) Col  1: Membrane thickness, TM
!             ( 2) Col  2: Normalized area MOI, I'=12I/TM^3
!             ( 3) Col  3: Normalized shear thick., TS/TM
!             ( 4) Col  4: Non structural mass, NSM
!             ( 5) Col  5: Stress recovery dist. Z1
!             ( 6) Col  6: Stress recovery dist. Z2

!  PSOLID = Array of integer data from PSOLIDL Bulk Data entries. Each row is for one PSOLID entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID 1

!  PUSER1 = Array of integer data from PUSER1 Bulk Data entries. Each row is for one PUSER1 entry read in B.D. and contains:
!             ( 1) Col  1: Property ID
!             ( 2) Col  2: Material ID

!  RPUSER1= Array of real data from PUSER1 Bulk Data entries. Each row is for one PUSER1 entry read in B.D. and contains:
!             ( 1) Col  1: User value  1
!             ( 2) Col  2: User value  2
!             ( 3) Col  3: User value  3
!             ( 4) Col  4: User value  4
!             ( 5) Col  5: User value  5
!             ( 6) Col  6: User value  6
!             ( 7) Col  7: User value  7
!             ( 8) Col  8: User value  8
!             ( 9) Col  9: User value  9
!             (10) Col 10: User value 10

!  MATL   = Array of integer data from Material Bulk Data entries. Each row is for one Material entry read in B.D. and contains:
!             ( 1) Col  1: Material ID
!             ( 2) Col  2: Material type: 1, in MATL array for MAT1 entry (isotropic)
!                                         2, in MATL array for MAT2 entry (anisotropic, 2D elements)
!                                         8, in MATL array for MAT8 entry (orthotropic)
!                                         9, in MATL array for MAT9 entry (anisotropic, 3D elements)

!  RMATL  = Array of real data from material Bulk Data entries. Each row is for one material entry read in B.D. and contains:
!
!           MAT1 (isotropic):
!           ----------------
!                ( 1) Col  1: E     Young's modulus, E
!                ( 2) Col  2: G     Shear modulus, G
!                ( 3) Col  3: NU    Poisson's ratio, NU
!                ( 4) Col  4: RHO   Mass density, RHO
!                ( 5) Col  5: ALPHA Thermal expansion coeff, ALPHA
!                ( 6) Col  6: TREF  Reference temperature, TREF
!                ( 7) Col  7: GE    Damping coefficient, GE
!                ( 8) Col  8: ST    Tension limit, ST
!                ( 9) Col  9: SC    Compression limit, SC
!                (10) Col 10: SS    Shear limit, SS

!           MAT2 (2D anisotropic):
!           ---------------------
!                ( 1) Col  1: G11   11 term in stress/strain matrix
!                ( 2) Col  2: G12   12 term in stress/strain matrix
!                ( 3) Col  3: G13   13 term in stress/strain matrix
!                ( 4) Col  4: G22   22 term in stress/strain matrix
!                ( 5) Col  5: G23   23 term in stress/strain matrix
!                ( 6) Col  6: G33   33 term in stress/strain matrix
!                ( 7) Col  7: RHO   Mass density
!                ( 8) Col  8: A1    Thermal expansion coeff in material x direction
!                ( 9) Col  9: A2    Thermal expansion coeff in material y direction
!                (10) Col 10: A3    Thermal expansion coeff in shear
!                (11) Col 11: TREF  Reference temperature for thermal expansion
!                (12) Col 12: GE    Damping coeff
!                (13) Col 13: ST    Tension limit, ST
!                (14) Col 14: SC    Compression limit, SC
!                (15) Col 15: SS    Shear limit, SS

!           MAT8 (orthotropic):
!           ------------------
!                ( 1) Col  1: E1    Modulus of elasticity in the longitudinal direction
!                ( 2) Col  2: E2    Modulus of elasticity in the lateral direction
!                ( 3) Col  3: NU12  Poisson's ratio
!                ( 4) Col  4: G12   In-plane shear modulus
!                ( 5) Col  5: G1Z   Transverse shear modulus in 1-Z plane
!                ( 6) Col  6: G2Z   Transverse shear modulus in 2-Z plane
!                ( 7) Col  7: RHO   Mass density
!                ( 8) Col  8: A1    Thermal expansion coeff in the longitudinal direction
!                ( 9) Col  9: A2    Thermal expansion coeff in the lateral direction
!                (10) Col 10: TREF  Reference temperature for thermal expansion
!                (11) Col 11: Xt    Allowable stress or strain in tension in the longitudinal directio
!                (12) Col 12: Xc    Allowable stress or strain in compression in the longitudinal direction
!                (13) Col 13: Yt    Allowable stress or strain in tension in the lateral direction
!                (14) Col 14: Yc    Allowable stress or strain in compression in the lateral direction
!                (15) Col 15: S     Allowable stress or strain for in-plane shear
!                (16) Col 16: GE    Damping coeff
!                (17) Col 17: F12   Interaction term in the tensor polynomial theory of failure
!                (18) Col 18: STRN  For max strain failure theory only. Indicates whether allowables are stress or strain allowables

! **********************************************************************************************************************************
! Element numbering and connection data
! -------------------------------------

      CHARACTER( 8*BYTE), ALLOCATABLE :: ETYPE(:)             ! NELE  x 1 array of elem types
      CHARACTER( 1*BYTE), ALLOCATABLE :: EOFF(:)              ! NELE  x 1 array of 'Y' for elem offsets or 'N' if not

      INTEGER(LONG)     , ALLOCATABLE :: EDAT(:)              ! NEDAT x 1 array of elem connection data
      INTEGER(LONG)     , ALLOCATABLE :: EPNT(:)              ! NELE  x 1 array of pointers to EDAT where data begins for an elem
      INTEGER(LONG)     , ALLOCATABLE :: ESORT1(:)            ! NELE  x 1 array of actual elem ID's for each elem
      INTEGER(LONG)     , ALLOCATABLE :: ESORT2(:)            ! NELE  x 1 array of internal elem numbers

!  ETYPE contains on 8 byte character value for each element in the B.D. deck:
! (ETYPE is set in the subrs that read Bulk Data elem connection entries: BD_CBAR, etc)

!        1)  'BAR     ' for BAR     element
!        2)  'BEAM    ' for BEAM    element
!        3)  'BUSH    ' for BUSH    element
!        4)  'ELAS1   ' for ELAS1   element
!        5)  'ELAS2   ' for ELAS2   element
!        6)  'ELAS3   ' for ELAS3   element
!        7)  'ELAS4   ' for ELAS4   element
!        8)  'HEXA8   ' for HEXA8   element
!        9)  'HEXA20  ' for HEXA20  element
!       10)  'PENTA6  ' for PENTA6  element
!       11)  'PENTA15 ' for PENTA15 element
!       12)  'PLOTEL  ' for PLOTEL  element
!       13)  'QUAD4   ' for QUAD4   element
!       14)  'QUAD4K  ' for QUAD4K  element
!       15)  'ROD     ' for ROD     element
!       16)  'SHEAR   ' for SHEAR   element
!       17)  'TETRA4  ' for TETRA4  element
!       18)  'TETRA10 ' for TETRA10 element
!       19)  'TRIA3K  ' for TRIA3K  element
!       20)  'TRIA3   ' for TRIA3   element
!       21)  'USER1   ' for USER1   element
!       22)  'USERIN  ' for USERIN  element

!  EDAT is an array of elem connection data containing the following info for:

!      1-2)   BAR, BEAM  8 words: (1st 4 read by call to subr ELEPO from subr BD_CBAR/BEAM, next 4 read by subr BD_CBAR/BEAM)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) V-vector key
!                                 6) Pin Flag A
!                                 7) Pin Flag B
!                                 8) Offset key

!                                 The V vector key is either:
!                                   1) a grid ID if the V vector is defined using a grid, or
!                                   2) -NVVEC where NVVEC is row number in array VVEC where the 3 comps of the V vector can be found

!                                 The Offset key is either:
!                                   1) 0 if there are no offsets for this CBAR, or
!                                   2) Offset key (row in array BAROFF where the 3 components at each of the 2 grids can be found)

!        3)   BUSH       9 words: (1st 5 read by call to subr ELEPRO, next 4 by subr BD_CBUSH)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Number of grids (2 if both G1 and G2 >0. If CBUSH G2 is blank, then number of grids = 1)
!                                 4) Grid A
!                                 5) Grid B (set to a zero if CBUSH field 4 is blank)
!                                 6) V-vector key
!                                 7) CID : Coord system ID in field 9 to define VVEC if VVEC not defined in fields 6-8.
!                                          If -99 then CID field was blank so V vec is defined by G0 or X1,2,3 on CBUSH
!                                 8) OCID: coord system ID for offset
!                                          If -1 then use elem local system with S (not S1,2,3) as offset
!                                 9) Offset key

!                                 The V vector key is either:
!                                   1) a grid ID if the V vector is defined using a grid, or
!                                   2) -NVVEC where NVVEC is row number in array VVEC where the 3 comps of the V vector can be found

!                                 The Offset key is either:
!                                   1) 0 if there are no offsets for this CBUSH, or
!                                   2) Offset key (row in array BUSHOFF where the 3 components at each of the 2 grids can be found)

!        4)   ELAS1      6 words: (all read by call to subr ELEPO from subr BD_CELAS1)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Components at Grid A
!                                 6) Components at Grid B

!        5)   ELAS2      6 words: (all read by call to subr ELEPO from subr BD_CELAS2)
!                                 1) Elem ID
!                                 2) Prop ID which is set to -EID since real props are on the CELAS2 entry
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Components at Grid A
!                                 6) Components at Grid B

!        6)   ELAS3      4 words: (all read by call to subr ELEPO from subr BD_CELAS3)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Scalar point A
!                                 4) Scalar point B

!        7)   ELAS4      4 words: (all read by call to subr ELEPO from subr BD_CELAS4)
!                                 1) Elem ID
!                                 2) Prop ID which is set to -EID since real props are on the CELAS2 entry
!                                 3) Scalar point A
!                                 4) Scalar point B

!        8)   HEXA8     10 words: (read by 1 or more calls to subr ELEPO from subr BD_HEXA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-8

!        9)   HEXA20    22 words: (read by 1 or more calls to subr ELEPO from subr BD_HEXA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-20

!       10)   PENTA6     8 words: (read by 1 or more calls to subr ELEPO from subr BD_PENTA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-6

!       11)   PENTA15   17 words: (read by 1 or more calls to subr ELEPO from subr BD_PENTA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-15

!       12)   PLOTEL

!    13-14)   QUAD4,K   10 words: (1st 6 read by call to subr ELEPO from subr BD_CQUAD, next 2 read by subr BD_CQUAD)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Grid C
!                                 6) Grid D
!                                 7) Material property angle key:
!                                     (a) For actual angle specified. the row in array MATANGLE where the angle is stored
!                                     (b) For angle defined by a coord sys ID, the negative of that ID
!                                 8) If matl angle is defined by a coord sys ID then this key is:
!                                     (a) 1 if the coord sys is basic, or
!                                     (b) 2 if the coord sys is anything other than basic
!                                 9) Plate offset key (row number in real array PLATEOFF where the angle is stored)
!                                10) Key for PSHELL (1) or PCOMP (2) property definition
!                                11) Plate thickness key (row number in real array PLATETHICK where thicknesses are stored if
!                                    defined on the connection entry and not on the PSHELL entry

!       15)   ROD        4 words: (all read by call to subr ELEPO from subr BD_ROD1)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B

!       16)   SHEAR               1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Grid C
!                                 6) Grid D

!       17)   TETRA4     6 words: (read by 1 or more calls to subr ELEPO from subr BD_TETRA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-4

!       18)   TETRA10   12 words: (read by 1 or more calls to subr ELEPO from subr BD_TETRA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) etc, Grids 1-10

!       19)   TRIA3,K    9 words: (1st 5 read by call to subr ELEPO from subr BD_CTRIA, next 2 read by subr BD_CTRIA)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Grid C
!                                 6) Material property angle key:
!                                     (a) For actual angle specified. the row in array MATANGLE where the angle is stored
!                                     (b) For angle defined by a coord sys ID, the negative of that ID
!                                 7) If matl angle is defined by a coord sys ID then this key is:
!                                     (a) 1 if the coord sys is basic, or
!                                     (b) 2 if the coord sys is anything other than basic
!                                 8) Plate offset key (row number in real array PLATEOFF where the angle is stored)
!                                 9) Key for PSHELL (1) or PCOMP (2) property definition
!                                10) Plate thickness key (row number in real array PLATETHICK where thicknesses are stored if
!                                    defined on the connection entry and not on the PSHELL entry

!       21)   USER1     11 words: (1st 7 read by call to subr ELEPO from subr BD_CUSER1, next 4 read by subr BD_CUSER1)
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) Grid A
!                                 4) Grid B
!                                 5) Grid C
!                                 6) Grid D
!                                 7) Grid V
!                                 8) - 11) Pin Flags (A, B, C, D)

!       21)   USERIN
!                                 1) Elem ID
!                                 2) Prop ID
!                                 3) NG, number of actual grids that the element is connected to (not incl SPOINT's)
!                                 4) NS, number of SPOINT's the element is connected to
!                                 5) CID0, actual coord system ID that defines the basic coord system for this elem relative to the
!                                    basic coord sys of the model in which this USERIN element is used
!                                 6) etc,

! Sorting of arrays:
! ------------------
!  1) After the Bulk Data has been read, all of the above arrays are in the order in which the elems were read from the B.D. deck.
!     Thus, at this point, ESORT2 is an array of consecutive integers starting with 1 and ending with NELE (the number of elements).

!  2) After subr ELESORT has run, arrays ESORT1, ESORT2, EOFF, EPNT, ETYPE are sorted together such that the elems in ESORT1 are in
!     numerically increasing order. EDAT remains as it was constructed when the Bulk Data was read - in the order in which elems
!     were placed in the Bulk Data Deck.

! **********************************************************************************************************************************
! BAR, BEAM, BUSH specific data
! -----------------------------

! BAR element specific data
! -------------------------

      CHARACTER( 9*BYTE)              :: BAROR_VVEC_TYPE     = ' '
                                                             ! Indicator of type of V vec on BAROR B.D. entry (grid or vector)

      CHARACTER(LEN=JCARD_LEN)        :: JBAROR(10)          = (/(' ', I=1,10)/)
                                                             ! JBAROR: copy of the 10 fields of a BAROR B.D. entry, if present

      INTEGER(LONG)                   :: BAROR_PID           = 0
                                                             ! Prop ID from BAROR Bulk Data entry

      INTEGER(LONG)                   :: BAROR_G0            = 0
                                                             ! V vector grid from BAROR Bulk Data entry

      REAL(DOUBLE)                    :: BAROR_VV(3)         = (/ZERO,ZERO,ZERO/)
                                                             ! BAROR_VV holds V vector comps on BAROR B.D. entry, if present

! **********************************************************************************************************************************
! BEAM element specific data
! --------------------------

! --- cbeam_stations begin --- !
      INTEGER(LONG)                   :: CBEAM_ACTIVE_NSTATIONS = 0
                                                             ! Number of active x/L stations copied into the current BEAM runtime state

      REAL(DOUBLE)                    :: CBEAM_ACTIVE_XL(MPBEAM_STATIONS) = (/ (ZERO, I=1,MPBEAM_STATIONS) /)
                                                             ! Active x/L station positions for the current BEAM runtime state
      REAL(DOUBLE)                    :: CBEAM_ACTIVE_RPROPS(MPBEAM_STATIONS,6) = ZERO
                                                             ! Active [A,I1,I2,I12,J,NSM] station properties for the current BEAM runtime state
      REAL(DOUBLE)                    :: CBEAM_ACTIVE_AREA_SCALE = ONE
                                                             ! Active stiffness-only area scale for the current BEAM runtime state
      REAL(DOUBLE)                    :: CBEAM_FORCE_B1(3,6) = ZERO
                                                             ! Beam section-force-to-stress map at the reference side
      REAL(DOUBLE)                    :: CBEAM_FORCE_B2(3,6) = ZERO
                                                             ! Beam section-force-to-stress map used in station interpolation
! --- cbeam_stations end --- !

      CHARACTER( 9*BYTE)              :: BEAMOR_VVEC_TYPE    = '         '
                                                             ! Indicator of type of V vec on BEAMOR B.D. entry (grid or vector)

      CHARACTER(LEN=JCARD_LEN)        :: JBEAMOR(10)         = (/('        ', I=1,10)/)
                                                             ! JBEAMOR: copy of the 10 fields of a BEAMOR B.D. entry, if present

      INTEGER(LONG)                   :: BEAMOR_PID          = 0
                                                             ! Prop ID from BEAMOR Bulk Data entry

      INTEGER(LONG)                   :: BEAMOR_G0           = 0
                                                             ! V vector grid from BEAMOR Bulk Data entry

      REAL(DOUBLE)                    :: BEAMOR_VV(3)        = (/ZERO,ZERO,ZERO/)
                                                             ! BAROR_VV holds V vector comps on BAROR B.D. entry, if present

! Data for BAR or BEAM
! --------------------

      REAL(DOUBLE), ALLOCATABLE       :: BAROFF(:,:)         ! BAR/BEAM elem offsets (3 offsets at each of 2 grids)

      REAL(DOUBLE), ALLOCATABLE       :: VVEC(:,:)           ! BAR elem V vector components for all distinct BAR elem V vecs

! **********************************************************************************************************************************
! BUSH element specific data
! --------------------------

      CHARACTER(1*BYTE)               :: BUSH_VVEC_OR_CID    = 'N'
                                                             ! Indicator if either one of BUSH_VVEC or BUSH_CID is present

      INTEGER(LONG)                   :: BUSH_CID            = -99
                                                             ! BUSH elem coord system ID (defines elem coord system if not blank on
!                                                              the CBUSH entry. Val of -99 is default and is used to indicate that
!                                                              the CID field on CBUSH was blank. BUSH_CID >= 0 indicates the coord
!                                                              system number of BUSH axes (and x doesn't have to be along grid line)

      INTEGER(LONG)                   :: BUSH_G2             ! Grid 2 grid number. If 0 then BUSH is grounded there

      INTEGER(LONG)                   :: BUSH_OCID           = -1
                                                             ! BUSH elem offset coord system ID (default -1 is local elem coord sys)

      INTEGER(LONG)                   :: BUSH_VVEC           = 0
!                                                            ! If > 0 it is the grid number used to define the BUSB v-vecror
!                                                            ! If < 0 is row number in array VVEC where v vector components are found
!                                                            ! If = 0 no v-vector was specified

      REAL(DOUBLE), ALLOCATABLE       :: BUSHOFF(:,:)        ! BAR/BEAM elem offsets (3 offsets at each of 2 grids)

      REAL(DOUBLE)                    :: BUSH_DXA            ! Offset dist of BUSH elem from end A in elem Xe dir.
      REAL(DOUBLE)                    :: BUSH_DXB            ! Offset dist of BUSH elem from end B in elem Xe dir.=L - DXA
      REAL(DOUBLE)                    :: BUSH_DY             ! Offset distance of BUSH elem from  in elem Ye direction
      REAL(DOUBLE)                    :: BUSH_DZ             ! Offset distance of BUSH elem from  in elem Ze direction

!  BAROFF  = NBAROFF  x 3 real array of offsets for CBAR, CBEAM

!  BUSHOFF = NBUSHOFF x 3 real array of offsets for CBUSH

!  VVEC   = NVVEC x 3 real array of V vector components for CBAR, CBEAM, CBUSH elems. Only the distinct V vectors are put into VVEC.
!           If a CBAR, CBEAM or CBUSH elem has a V vector that is already loaded into array VVEC, then that V vector is used.

! **********************************************************************************************************************************
! Plate element specific data
! ---------------------------

      REAL(DOUBLE), ALLOCATABLE       :: MATANGLE(:)         ! Array of material property angles (one value per nonzero THETA on the
!                                                              plate element connection entries)
      REAL(DOUBLE), ALLOCATABLE       :: PLATEOFF(:)         ! Array of plate offsets (one value per nonzero ZOFFS on the
!                                                              plate element connection entries)

      REAL(DOUBLE), ALLOCATABLE       :: PLATETHICK(:)       ! Array of plate thicknesses for quads, trias that have thickness
!                                                              defined on the connection entry

! **********************************************************************************************************************************
! Overall element data
! --------------------

                                                             ! Character variables that describe the various element types
      CHARACTER( 8*BYTE)              :: ELMTYP(METYPE)      = (/'BAR     ',      & !          1
                                                                 'BEAM    ',      & !          2
                                                                 'BUSH    ',      & !          3
                                                                 'ELAS1   ',      & !          4
                                                                 'ELAS2   ',      & !          5
                                                                 'ELAS3   ',      & !          6
                                                                 'ELAS4   ',      & !          7
                                                                 'HEXA8   ',      & !          8
                                                                 'HEXA20  ',      & !          9
                                                                 'PENTA6  ',      & !         10
                                                                 'PENTA15 ',      & !         11
                                                                 'PLOTEL  ',      & !         12
                                                                 'QUAD4   ',      & !         13
                                                                 'QUAD4K  ',      & !         14
                                                                 'QUAD8   ',      & !         15
                                                                 'ROD     ',      & !         16
                                                                 'SHEAR   ',      & !         17
                                                                 'TETRA4  ',      & !         18
                                                                 'TETRA10 ',      & !         19
                                                                 'TRIA3K  ',      & !         20
                                                                 'TRIA3   ',      & !         21
                                                                 'USER1   ',      & !         22
                                                                 'USERIN  '/)       !         23

                                                             ! Character name for output purposed in LINK9 WRTELi subr's
      CHARACTER(13*BYTE)              :: ELEM_ONAME(METYPE)  = (/'B A R        ', & !          1
                                                                 'B E A M      ', & !          2
                                                                 'B U S H      ', & !          3
                                                                 'E L A S 1    ', & !          4
                                                                 'E L A S 2    ', & !          5
                                                                 'E L A S 3    ', & !          6
                                                                 'E L A S 4    ', & !          7
                                                                 'H E X A  8   ', & !          8
                                                                 'H E X A  20  ', & !          9
                                                                 'P E N T A  6 ', & !         10
                                                                 'P E N T A  15', & !         11
                                                                 'P L O T E L  ', & !         12
                                                                 'Q U A D 4    ', & !         13
                                                                 'Q U A D 4 K  ', & !         14
                                                                 'Q U A D 8    ', & !         15
                                                                 'R O D        ', & !         16
                                                                 'S H E A R    ', & !         17
                                                                 'T E T R A  4 ', & !         18
                                                                 'T E T R A  10', & !         19
                                                                 'T R I A 3 K  ', & !         20
                                                                 'T R I A 3    ', & !         21
                                                                 'U S E R 1    ', & !         22
                                                                 'U S E R I N  '/)  !         23

                                                             ! Array of number of grid points for the various element types
      INTEGER(LONG)                   :: NELGP(METYPE)       =  (/ 2,             & ! BAR      1
                                                                   2,             & ! BEAM     2
                                                                   2,             & ! BUSH     3
                                                                   2,             & ! ELAS1    4
                                                                   2,             & ! ELAS2    5
                                                                   2,             & ! ELAS3    6
                                                                   2,             & ! ELAS4    7
                                                                   8,             & ! HEXA8    8
                                                                  20,             & ! HEXA20   9
                                                                   6,             & ! PENTA6  10
                                                                  15,             & ! PENTA15 11
                                                                   2,             & ! PLOTEL  12
                                                                   4,             & ! QUAD4   13
                                                                   4,             & ! QUAD4K  14
                                                                   8,             & ! QUAD8   15
                                                                   2,             & ! ROD     16
                                                                   4,             & ! SHEAR   17
                                                                   4,             & ! TETRA4  18
                                                                  10,             & ! TETRA10 19
                                                                   3,             & ! TRIA3K  20
                                                                   3,             & ! TRIA3   21
                                                                   4,             & ! USER1   22
                                                                   0/)              ! USERIN  23

                                                             ! Array of number of stress recovery points for various elem types
      INTEGER(LONG)                   :: NUM_SEi(METYPE)     =  (/ 1,             & ! BAR      1
! --- cbeam_stations begin --- !
                                                                  10,             & ! BEAM     2
! --- cbeam_stations end --- !
                                                                   1,             & ! BUSH     3
                                                                   1,             & ! ELAS1    4
                                                                   1,             & ! ELAS2    5
                                                                   1,             & ! ELAS3    6
                                                                   1,             & ! ELAS4    7
                                                                   9,             & ! HEXA8    8
                                                                   9,             & ! HEXA20   9
                                                                   7,             & ! PENTA6  10
                                                                   7,             & ! PENTA15 11
                                                                   0,             & ! PLOTEL  12
                                                                   5,             & ! QUAD4   13
                                                                   5,             & ! QUAD4K  14
                                                                   5,             & ! QUAD8   15
                                                                   1,             & ! ROD     16
                                                                   1,             & ! SHEAR   17
                                                                   5,             & ! TETRA4  18
                                                                   5,             & !,TETRA10 19
                                                                   1,             & ! TRIA3K  20
                                                                   1,             & ! TRIA3   21
                                                                   1,             & ! USER1   22
                                                                   0/)              ! USERIN  23

! **********************************************************************************************************************************
! Individual element data generated one element at a time (in subr EMG)
! ---------------------------------------------------------------------

! Fixed size arrays

      CHARACTER(4*BYTE)               :: ANY_FAILURE_THEORY  = 'N'
                                                             ! Indicates type of failure theory

      CHARACTER(1*BYTE)               :: CAN_ELEM_TYPE_OFFSET= 'N'
                                                             ! Indicates if element type can be offset


      CHARACTER(LEN(BLNK_SUB_NAM))    :: ERR_SUB_NAM(MEFE)   = ' '
                                                             ! Subrs where errors occur. Used for debug Space GASS err array outputs

      CHARACTER(4*BYTE)               :: FAILURE_THEORY      = 'NONE'
                                                             ! Indicates type of failure theory

      CHARACTER(8*BYTE)               :: PCOMP_PROPS         = 'N'
                                                             ! If 'Y' then element uses PCOMP properties

      CHARACTER(3*BYTE)               :: PCOMP_LAM           = '   '
                                                             ! If 'SYM' then PCOMP is a symmetric layup, otherwise it is nonsym

      CHARACTER(1*BYTE)               :: SHELL_T_MOD         = 'N'
                                                             ! Reset to 'Y' if SHELL_T is modified (as for example when
                                                             ! INTL_MID(3) = 0 or G1Z, G2Z transverse shear modulii = 0 (in order
                                                             ! to get subrs QPLT2,3 and TPLT2 to work when transv shear flex = 0

      CHARACTER(8*BYTE)               :: TE_IDENT            = 'N'
                                                             ! If 'Y' then TE element transformation matrix is an identity matrix

      CHARACTER(8*BYTE)               :: TYPE                = '        '
                                                             ! The type of the specific elem being processed (value from ELMTYP)

      INTEGER(LONG)                   :: EID                 = 0
                                                             ! Current elem ID (actual, not internal)

      INTEGER(LONG)                   :: ELAS_COMP(2)        ! Displ comp nos (1 thru 6), at each elem G.P. for a CELAS elem

      INTEGER(LONG)                   :: ELDOF               = 0
                                                             ! The specific number of DOF's for the current elem

      INTEGER(LONG)                   :: ELGP                = 0
                                                             ! The number of grid points that the the current elem connects to

      INTEGER(LONG)                   :: ELNO                = 0
                                                             ! The internal elem ID ( 1 to NELE) of the current element

      INTEGER(LONG)                   :: EMG_IFE(MEFE,MEFEI) = RESHAPE ( (/(ZERO, I=1,MEFE*MEFEI)/), (/MEFE,MEFEI/) )
                                                              ! Array of integer data for EMG errors

      INTEGER(LONG)                   :: EMG_IWE(MEWE,MEWEI) = RESHAPE ( (/(ZERO, I=1,MEWE*MEWEI)/), (/MEWE,MEWEI/) )
                                                             ! Array of integer data for EMG warnings

      INTEGER(LONG)                   :: NUM_EMG_FATAL_ERRS  = 0
                                                             ! The number of errors found in one execution of subr EMG

      INTEGER(LONG)                   :: INTL_MID(MEMATC)    = (/(0, I=1,MEMATC)/)
                                                             ! Internal material ID's for the current elem

      INTEGER(LONG)                   :: INTL_PID            = 0
                                                             ! Internal property ID for the current elem

      INTEGER(LONG)                   :: IOR3D_MAX           = 0
                                                             ! Max integration order for 3D solid elements in any one execution


      INTEGER(LONG)                   :: ISOLID(MPSOLID)     = (/0, 0, -1, 0, 0, 0/)
                                                             ! Data from PSOLID B.D. entry on integration for solid elements
!                                                              ISOLID(1) = PSOLID ID
!                                                              ISOLID(2) = Matl ID
!                                                              ISOLID(3) = Matl coord system ID
!                                                              ISOLID(4) = Integration order
!                                                              ISOLID(5) = Stress locations
!                                                              ISOLID(6) = Integration scheme

      INTEGER(LONG)                   :: MTRL_TYPE(MEMATC)   = (/(0, I=1,MEMATC)/)
                                                             ! Material types for this element (1 for MAT1, 8 for MAT8, etc). There
                                                             ! can be diff types for membrane, bending, shear and mem/bend coupling
                                                             ! Currently the 1st 4 col's are for:
                                                             ! 1) membrane, 2) bending, 3) trans shear, 4) membrane/bending coupling

      INTEGER(LONG)                   :: NUM_PLIES           = 1
                                                             ! Number of plies in a shell using PCOMP

      INTEGER(LONG)                   :: PLY_NUM             = 1
                                                             ! The ply number in a composite layup defined by a PCOMP

      REAL(DOUBLE)                    :: ALPVEC(6,MEMATC)    = RESHAPE ( (/(ZERO, I=1,6*MEMATC)/), (/6,MEMATC/) )
                                                             ! Vector of CTE's. x,y,z,xy,yz,zx in material coordinates initially
                                                             ! then transformed to element coordinates for most element types.
                                                             ! For shells, column i is for material MIDi. i=1..4.
                                                             ! For solids, only the column 1 is used.
                                                             ! They are engineering strain per unit temperature change
                                                             ! rather than tensor components so the shear components
                                                             ! must be divided by 2 before transforming and used as-is
                                                             ! for multiplying by elasticity.

      REAL(DOUBLE)                    :: BENSUM              = ZERO
                                                             ! Sum of diag stiffness from KB for rotation DOF's

      REAL(DOUBLE)                    :: BMEANT(8,12)        = RESHAPE ( (/(ZERO, I=1,8*12)/), (/8,12/) )
                                                             ! Transpose of strain-displ matrix for the quad memb with warp (HBAR)

      REAL(DOUBLE)                    :: EB(3,3)             = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Plane stress material matrix for bending

      REAL(DOUBLE)                    :: EM(3,3)             = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Plane stress material matrix for in-plane stress

      REAL(DOUBLE)                    :: EBM(3,3)            = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Bending/membrane coupling material matrix

      REAL(DOUBLE)                    :: ES(6,6)             = RESHAPE ( (/(ZERO, I=1,6*6)/), (/6,6/) )
                                                             ! 3D stress material matrix

      REAL(DOUBLE)                    :: ET(2,2)             = RESHAPE ( (/(ZERO, I=1,2*2)/), (/2,2/) )
                                                             ! 2D transverse shear material matrix

      REAL(DOUBLE)                    :: ELEM_LEN_AB         = ZERO
                                                             ! Length of elem between elem ends A and B including effects of offsets

      REAL(DOUBLE)                    :: ELEM_LEN_12         = ZERO
                                                             ! Length of elem between grids 1 and 2 including effects of offsets

      REAL(DOUBLE)                    :: EMG_RFE(MEFE,MEFER) = RESHAPE ( (/(ZERO, I=1,MEFE*MEFER)/), (/MEFE,MEFER/) )
                                                             ! Array of real data for EMG errors

      REAL(DOUBLE)                    :: EMG_RWE(MEWE,MEWER) = RESHAPE ( (/(ZERO, I=1,MEWE*MEWER)/), (/MEWE,MEWER/) )
                                                             ! Array of real data for EMG warnings

      REAL(DOUBLE)                    :: EMAT(MEMATR,MEMATC) = RESHAPE ( (/(ZERO, I=1,MEMATR*MEMATC)/), (/MEMATR,MEMATC/) )
                                                             ! Array of matl props for the elem

      REAL(DOUBLE)                    :: EPROP(MEPROP)       = (/(ZERO, I=1,MEPROP)/)
                                                             ! Array of elem properties (area, area MOI;s, etc) for the elem

      REAL(DOUBLE)                    :: FCONV(3)            = (/(ZERO, I=1,3)/)
                                                             ! Array of constants to convert stresses to engr forces for elem

      REAL(DOUBLE)                    :: HBAR                = ZERO
                                                             ! For the quad elem, the dist from the mean plane to the G.P.'s

      REAL(DOUBLE)                    :: HEXA_DELTA          = ZERO
                                                             ! Angle to rotate TE to go from local x axis being parallel to side 1-2
!                                                              to local x axis splitting the angle between the 2 diagonals

      REAL(DOUBLE)                    :: HEXA_GAMMA          = ZERO
                                                             ! Angle between side 1-2 and diag from G.P 2 to G.P.4 of a HEXA element

      REAL(DOUBLE)                    :: HEXA_THETA          = ZERO
                                                             ! Angle between side 1-2 and diag from G.P 1 to G.P.3 of a HEXA element

      REAL(DOUBLE)                    :: MASS_PER_UNIT_AREA  = ZERO
                                                             ! Shell element mass/area (incl NSM and mat'l density times thickness)

      REAL(DOUBLE)                    :: MXWARP              = ZERO
                                                             ! Current quad elem maximum warp before warning message written
!                                                              A finite value for this is calcd in sub ELMGM2 based on the elem size

      REAL(DOUBLE)                    :: PHI_SQ              = ZERO
                                                             ! Finite elem shear correction factor

      REAL(DOUBLE)                    :: PSI_HAT             = ZERO
                                                             ! BENSUM/SHRSUM (Note: this goes to 0 as KS goes to infinity)

      REAL(DOUBLE)                    :: QUAD_DELTA          = ZERO
                                                             ! Angle to rotate TE to go from local x axis being parallel to side 1-2
!                                                              to local x axis splitting the angle between the 2 diagonals

      REAL(DOUBLE)                    :: QUAD_GAMMA          = ZERO
                                                             ! Angle between side 1-2 and diag from G.P 2 to G.P.4 of a QUAD element

      REAL(DOUBLE)                    :: QUAD_THETA          = ZERO
                                                             ! Angle between side 1-2 and diag from G.P 1 to G.P.3 of a QUAD element

      REAL(DOUBLE)                    :: RHO(MEMATC)         = (/(ZERO, I=1,MEMATC)/)
                                                             ! Vector of material densities

      REAL(DOUBLE)                    :: SHELL_A(3,3)        = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Membrane force resultant/strain matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_B(3,3)        = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Membrane/bend coupling force resultant/strain matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_D(3,3)        = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Bending force resultant/strain matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_ALP(6,MEMATC) = RESHAPE ( (/(ZERO, I=1,6*MEMATC)/), (/6,MEMATC/) )
                                                             ! Effective CTE matrix for shell elems (used for MAT2 output on PCOMP)

      REAL(DOUBLE)                    :: SHELL_AALP(3)       = (/(ZERO, I=1,3)/)
                                                             ! Membrane matl matrix times CTE matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_BALP(3)       = (/(ZERO, I=1,3)/)
                                                             ! Mem/bend coupling matl matrix times CTE matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_DALP(3)       = (/(ZERO, I=1,3)/)
                                                             ! Bending matl matrix times CTE matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_TALP(2)       = (/(ZERO, I=1,2)/)
                                                             ! Transverse shear matl matrix times CTE matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_T(2,2)        = RESHAPE ( (/(ZERO, I=1,2*2)/), (/2,2/) )
                                                             ! Transverse shear force resultant/strain matrix for shell elems

      REAL(DOUBLE)                    :: SHELL_PROP_ALP(3)   = (/(ZERO, I=1,3)/)
                                                             ! matrix resulting from material matrix times coeff of thermal
!                                                              expansion vector times a property (thickness or bending MOI)

      REAL(DOUBLE)                    :: SHELL_STR_ANGLE(5) = (/(ZERO, I=1,5)/)
                                                             ! Angle of element coordinate system from cartesial local coordinate system
                                                             ! at each stress output point

      REAL(DOUBLE)                    :: SHRSUM              = ZERO
                                                             ! Sum of diag stiffness from KS for rotation DOF's

      REAL(DOUBLE)                    :: ULT_STRE(9,MEMATC)  = RESHAPE ( (/(ZERO, I=1,9*MEMATC)/), (/9,MEMATC/) )
                                                             ! Array of material stress allowables for the elem
                                                             ! tens/compr dir 1, tens/compr dir 2, tens/compr dir 3, shear 12,23,31

      REAL(DOUBLE)                    :: ULT_STRN(9,MEMATC)  = RESHAPE ( (/(ZERO, I=1,9*MEMATC)/), (/9,MEMATC/) )
                                                             ! Array of material strain allowables for the elem

      REAL(DOUBLE)                    :: STRAIN(9)           = (/(ZERO, I=1,9)/)
                                                             ! Array of elem strains for one S/C

      REAL(DOUBLE)                    :: STRESS(9)           = (/(ZERO, I=1,9)/)
                                                             ! Array of elem stresses for one S/C

      REAL(DOUBLE)                    :: TE(3,3)             = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Coord system transformation matrix (3 x 3) such that UEL = TE*UEB

      REAL(DOUBLE)                    :: TE_GA_GB(3,3)       = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Coord system transformation matrix (3 x 3) for BUSH

      REAL(DOUBLE)                    :: T1P(6,6)            = RESHAPE ( (/(ZERO, I=1,6*6)/), (/6,6/) )
                                                             ! Coord transf matrix for 6x6 ply material matrices (6 stresses)

      REAL(DOUBLE)                    :: T1M(3,3)            = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Portion of T1P for shell membrane and bending stresses

      REAL(DOUBLE)                    :: T1T(2,2)            = RESHAPE ( (/(ZERO, I=1,2*2)/), (/2,2/) )
                                                             ! Portion of T1P for shell transverse shear stresses

      REAL(DOUBLE)                    :: T2P(6,6)            = RESHAPE ( (/(ZERO, I=1,6*6)/), (/6,6/) )
                                                             ! Inverse of T1P

      REAL(DOUBLE)                    :: T2M(3,3)            = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Portion of T2P for shell membrane and bending stresses

      REAL(DOUBLE)                    :: T2T(2,2)            = RESHAPE ( (/(ZERO, I=1,2*2)/), (/2,2/) )
                                                             ! Portion of T2P for shell transverse shear stresses

      REAL(DOUBLE)                    :: THETAM              = ZERO
                                                             ! Material orientation angle relative to side 1-2 of element

      REAL(DOUBLE)                    :: THETA_PLY           = ZERO
                                                             ! Angle of ply longitudinal axis relative to THETAM

      REAL(DOUBLE)                    :: TPLY                = ZERO
                                                             ! Ply thickness

      REAL(DOUBLE)                    :: TREF(MEMATC)        = (/(ZERO, I=1,MEMATC)/)
                                                             ! Vector of ref temps (membrane, bending, transverse shear)

      REAL(DOUBLE)                    :: WARP_WARN           = ZERO
                                                             ! Amount of quad warp before a warning msg is written. It is calculated
!                                                              in subr ELMGM2 as EPS4 times the average diagonal length

      REAL(DOUBLE)                    :: XTB(3,3)            = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Basic elem coords of the 3 grids of 1 of the 4 TPLT2 tria elems that
!                                                              form a MIN4T QUAD4

      REAL(DOUBLE)                    :: XTL(3,3)            = RESHAPE ( (/(ZERO, I=1,3*3)/), (/3,3/) )
                                                             ! Local elem coords of the 3 grids of 1 of the 4 TPLT2 tria elems that
!                                                              form a MIN4T QUAD4

      REAL(DOUBLE)                    :: ZPLY                = ZERO
                                                             ! Coord of mid plane of a ply rel to mid plane of element

      REAL(DOUBLE)                    :: ZS(9)               = (/(ZERO, I=1,9)/)
                                                             ! Array of stress recovery coefficients for the current elem

      REAL(DOUBLE)                    :: ZOFFS               = ZERO
                                                             ! Plate offset

! Allocatable arrays

      CHARACTER(1*BYTE), ALLOCATABLE  :: OFFSET(:)           ! 'Y' or 'N' to indicate whether elem has offsets at any of it's G.P.'s

      INTEGER(LONG), ALLOCATABLE      :: AGRID(:)            ! Actual grid points for the current elem

      INTEGER(LONG), ALLOCATABLE      :: BGRID(:)            ! Row numbers in array GRID_ID where current elems actual grids exist

      INTEGER(LONG), ALLOCATABLE      :: DOFPIN(:)           ! Integers 1-6 (or 0) to indicate displ comps pinned at each elem grid

      REAL(DOUBLE) , ALLOCATABLE      :: BE1(:,:,:)          ! Array of elem membrane strain per unit displ for the elem

      REAL(DOUBLE) , ALLOCATABLE      :: BE2(:,:,:)          ! Array of elem bending strain per unit displ for the elem

      REAL(DOUBLE) , ALLOCATABLE      :: BE3(:,:,:)          ! Array of elem transv shear strain per unit displ for elem

      REAL(DOUBLE) , ALLOCATABLE      :: DT(:,:)             ! Temp and temp grad data for the current elem for all subcases

      REAL(DOUBLE) , ALLOCATABLE      :: KE(:,:)             ! Current elem stiffness matrix

      REAL(DOUBLE)                    :: KEO_BUSH(12,12)     ! BUSH elem stiffness matrix in local elem coords with offsets

      REAL(DOUBLE) , ALLOCATABLE      :: KED(:,:)            ! Current elem linear differential stiffness matrix

      REAL(DOUBLE) , ALLOCATABLE      :: KEG(:,:)            ! Current elem stiffness matrix in global coords with offsets

      REAL(DOUBLE) , ALLOCATABLE      :: KEM(:,:)            ! Current elem material stiffness matrix

      REAL(DOUBLE) , ALLOCATABLE      :: ME(:,:)             ! Current elem mass matrix

      REAL(DOUBLE) , ALLOCATABLE      :: OFFDIS(:,:)         ! Elem offset distances used in offset calcs in subr ELMOFF
!                                                              Set equal to OFFDIS_G in subr ELMGM1 since offsets are processed
!                                                              in global coords. So, any reference to OFFDIS after ELMGM1 will
!                                                              mean offsets in global coords.

      REAL(DOUBLE)                    :: OFFDIS_L(2,3)       ! Elem offsets for BUSH in local, CID, coord system

      REAL(DOUBLE) , ALLOCATABLE      :: OFFDIS_O(:,:)       ! Elem offsets based on input values (saved for posterity)

      REAL(DOUBLE) , ALLOCATABLE      :: OFFDIS_B(:,:)       ! Elem offsets in basic coord system

      REAL(DOUBLE) , ALLOCATABLE      :: OFFDIS_G(:,:)       ! Elem offsets in global coord system

      REAL(DOUBLE)                    :: OFFDIS_GA_GB(2,3)   ! Elem offsets for BUSH in coord sys along & normal to BUSH axis

      REAL(DOUBLE) , ALLOCATABLE      :: PEB(:)              ! Array of current elem nodal forces in basic coords for one S/C

      REAL(DOUBLE) , ALLOCATABLE      :: PEG(:)              ! Array of current elem nodal forces in global coords for one S/C

      REAL(DOUBLE)                    :: PE_GA_GB(12)        ! BUSH element nodal loads in coord sys with x along the line GA-GB

      REAL(DOUBLE) , ALLOCATABLE      :: PEL(:)              ! Array of current elem nodal forces in local elem coords for one S/C
!                                                              PEL is always calculated, and stays in, the local element coord system
!                                                              If PARAM ELFORCEN is BASIC or GLOBAL then PEB or PEG is used
!                                                              for output to the E L E M  N O D A L   F O R C E S

      REAL(DOUBLE) , ALLOCATABLE      :: PRESS(:,:)          ! Elem pressure data for the current elem for all subcases

      REAL(DOUBLE) , ALLOCATABLE      :: PTE(:,:)            ! Elem thermal loads for the current elem for all subcases

      REAL(DOUBLE) , ALLOCATABLE      :: PPE(:,:)            ! Elem pressure loads for the current elem for all subcases

      REAL(DOUBLE) , ALLOCATABLE      :: SE1(:,:,:)          ! Array of elem membrane stress per unit displ for the elem

      REAL(DOUBLE) , ALLOCATABLE      :: SE2(:,:,:)          ! Array of elem bending stress per unit displ for the elem

      REAL(DOUBLE) , ALLOCATABLE      :: SE3(:,:,:)          ! Array of elem transv shear strain per unit displ for elem

      REAL(DOUBLE) , ALLOCATABLE      :: STE1(:,:,:)         ! Elem thermal stress corrections for membrane stress

      REAL(DOUBLE) , ALLOCATABLE      :: STE2(:,:,:)         ! Elem thermal stress corrections for bending stress

      REAL(DOUBLE) , ALLOCATABLE      :: STE3(:,:,:)         ! Elem thermal stress corrections for transv shear stress

      REAL(DOUBLE) , ALLOCATABLE      :: UGG(:)              ! Array of displs at grids of one elem in global coords for one S/C
                                                             ! UGG is a partition of the G-set global displs solved for in MYSTRAN
                                                             ! for one elements grids

      REAL(DOUBLE) , ALLOCATABLE      :: UEG(:)              ! Array of displs at nodes* of one elem in global coords for one S/C

      REAL(DOUBLE) , ALLOCATABLE      :: UEB(:)              ! Array of displs at nodes* of one elem in basic coords for one S/C

      REAL(DOUBLE) , ALLOCATABLE      :: UEL(:)              ! Array of displs at nodes* of one elem in local elem coords for 1 S/C

      REAL(DOUBLE) , ALLOCATABLE      :: XEB(:,:)            ! Array of basic coordinates of the grids of the elem + V vector

      REAL(DOUBLE) , ALLOCATABLE      :: XEL(:,:)            ! Array of local elem coordinates of the grids of the elem

      REAL(DOUBLE) , ALLOCATABLE      :: XGL(:,:)            ! Array of local elem coordinates of the Gauss points of the elem


! Data particular to the USERIN element

      CHARACTER( 8*BYTE), ALLOCATABLE :: USERIN_MAT_NAMES(:,:)
!                                                            ! USERIN_MAT_NAMES(i,1) is the USERIN elem i stiff matrix name
!                                                            ! USERIN_MAT_NAMES(i,2) is the USERIN elem i mass  matrix name
!                                                            ! USERIN_MAT_NAMES(i,3) is the USERIN elem i load  matrix name
!                                                            ! USERIN_MAT_NAMES(i,4) is the USERIN elem i RBM0  matrix name

      CHARACTER( 8*BYTE)              :: USERIN_LOAD_MAT_NAME! Load  matrix name from USERIN_MAT_NAMES

      CHARACTER( 8*BYTE)              :: USERIN_MASS_MAT_NAME! Mass  matrix name from USERIN_MAT_NAMES

      CHARACTER( 8*BYTE)              :: USERIN_RBM0_MAT_NAME! RBM0  matrix name from USERIN_MAT_NAMES

      CHARACTER( 8*BYTE)              :: USERIN_STIF_MAT_NAME! Stiff matrix name from USERIN_MAT_NAMES

      INTEGER(LONG), ALLOCATABLE      :: USERIN_ACT_GRIDS(:) ! Array of actual grids (not incl SPOINT's) for a USERIN elem

      INTEGER(LONG), ALLOCATABLE      :: USERIN_ACT_COMPS(:) ! Array of integers 1-6 of grid components for USERIN_ACT_GRIDS. Each
!                                                              row has a list of all comps for that grid, e.g. 125 for comps 1,2,5

      INTEGER(LONG)                   :: USERIN_CID0         ! ID that defines basic coord sys of USERIN elem rel to model basic sys

      INTEGER(LONG)                   :: USERIN_IN4_INDEX  = 0
                                                             ! Index number (1 to NUM_IN4_files) of an IN4fil

      INTEGER(LONG)                   :: USERIN_NUM_BDY_DOF  = 0
                                                             ! Number of physical boundary DOF's for a USERIN elem

      INTEGER(LONG)                   :: USERIN_NUM_ACT_GRDS = 0
                                                             ! Number of actual grids for a USERIN elem (not counting SPOINT's)

      INTEGER(LONG)                   :: USERIN_NUM_SPOINTS  = 0
                                                             ! Number of SPOINT's for a USERIN elem

      REAL(DOUBLE)                    :: USERIN_RBM0(6,6)    = RESHAPE ( (/(ZERO, I=1,6*6)/), (/6,6/) )
                                                             ! 6x6 rigid body mass matrix for a USERIN element rel to GRDPNT loc

! * "nodes" refers to an element ends which may be offset from the grid points to which the element is attached

! Notes:
! ------

!  EMAT     : Each col has real matl props for one of the aspects (memb, bend, transverse shear, memb-bend coupling) for the elem

!  MXWARP   : See subr ELMGM2 for calculation. It is a small number times the avg diagonal.

!  OFFDIS   : The 3 cols are the  x, y, z (elem coords) of the offsets at each of the elem grids

!  PEB      : PEB = TE' x PEL

!  PEG      : PEG = T0G' x PEB (T0G transformation matrix is generated, and used, in subr ELMDIS)

!  PEL      : PEL = KE*UEL - PTE

!  ULT_STRE : Each col has stress allowables for one of the aspects (memb, bend, transverse shear, memb-bend coupling) for the elem

!  ULT_STRN : Each col has strain allowables for one of the aspects (memb, bend, transverse shear, memb-bend coupling) for the elem

!  STEi     : These are the corrections that must be made to stresses calculated based on total element strains to subtract out the
!             stress due to the free thermal expansion part of the total strain (i.e. they are E*ALPHA*DT types of terms)

!  STRESS   : No's 1 - 3 are membrane stresses, no's 4 - 6 are bending stresses and no's 7, 8 are transverse shear stresses.
!             STRESS(9) not currently used. The transverse shear stresses are only used for calculating elem engineering forces
!             (so transverse shear elem forces can be output, but the transverse stress itself cannot be output)

!  XEB      : XEB has 1 more row than needed for the number of grids. The last row in XEB is a vector, measured from grid a, in the
!             basic coord directions, along the v vector.

!  ZS       : e.g. the Y, Z coords of 4 points on a CBAR cross-section at which to calc stresses


! **********************************************************************************************************************************
! Eigenvalue extraction info - from a B.D. entry EIGR or EIGRL
! ------------------------------------------------------------

      CHARACTER(LEN=JCARD_LEN)        :: EIG_METH            = ' '
                                                             ! Eigenvalue extraction method

      CHARACTER(LEN=JCARD_LEN)        :: EIG_NORM            = 'MASS'
                                                             ! Eigenvector renorm nethod (e.g. MAX, POINT)

      CHARACTER(LEN=JCARD_LEN)        :: EIG_LAP_MAT_TYPE    = 'DPB'
                                                             ! Type of LAPACK matrices to use in LANCZOS.
!                                                              DGB is the type used in the original ARPACK subr dsband
!                                                              DPB uses less disk storage but may not work for free-free eigens

      CHARACTER(1*BYTE)               :: EIG_VECS            = 'Y'
                                                             ! Indicator of whether to calc eigenvecs

      INTEGER(LONG)                   :: EIG_SID             = 0
                                                             ! Eigenvalue extraction method set ID

      INTEGER(LONG)                   :: EIG_N1              = 0
                                                             ! For GIV, MGIV the lowest  mode number for vector calc

      INTEGER(LONG)                   :: EIG_N2              = 0
                                                             ! The highest mode number for vector calc

      INTEGER(LONG)                   :: EIG_COMP            = 0
                                                             ! For GIV, MGIV the displ component value for vector renorm

      INTEGER(LONG)                   :: EIG_GRID            = 0
                                                             ! For GIV, MGIV the grid for eigenvector renormalization

      INTEGER(LONG)                   :: EIG_LANCZOS_NEV_DELT= 2
                                                             ! Number to add to est num eigens when search is on freq range

      INTEGER(LONG)                   :: EIG_MODE            = 2
                                                             ! For Lanczos, the "mode" (see IPARAM(7) in ARPACK subr dsband)

      INTEGER(LONG)                   :: EIG_MSGLVL          = 0
                                                             ! For Lanczos, the message level in ARPACK subr dsband)

      INTEGER(LONG)                   :: EIG_NCVFACL         = 3
                                                             ! For Lanczos, the factor in ARPACK subr dsband to multiply the number
!                                                              of requested eigenvalues/vectors by to get NCV. This is used to
!                                                              dim several arrays in the ARPACK version of the Lanczos eigenval
!                                                              extraction method.(see subr DSBAND in module ARPACK_LANCZOS_1).
!                                                              It must be > 1

      INTEGER(LONG)                   :: MIJ_ROW             = 0
                                                             ! Row no. of largest off-diag gen. mass term.

      INTEGER(LONG)                   :: MIJ_COL             = 0
                                                             ! Col no. of largest off-diag gen. mass term.

      INTEGER(LONG)                   :: NUM_FAIL_CRIT       = 0
                                                             ! Number of terms in the gen mass matrix failing orthog criteria.

      REAL(DOUBLE)                    :: EIG_CRIT            = 0
                                                             ! For GIV, MGIV the criteria for off-diag gen mass terms

      REAL(DOUBLE)                    :: EIG_CRIT_DEF        = ONEPM4
                                                             ! For GIV, MGIV the criteria for off-diag gen mass terms

      REAL(DOUBLE)                    :: EIG_FRQ1            = ZERO
                                                             ! Lower  frequency in the search range

      REAL(DOUBLE)                    :: EIG_FRQ2            = ZERO
                                                             ! Higher frequency in the search range

      REAL(DOUBLE)                    :: EIG_SIGMA           = -ONE
                                                             ! For Lanczos, the shift frequency

      REAL(DOUBLE)                    :: MAXMIJ              = ZERO
                                                             ! Largest off-diag term in generalized mass matrix.                                                            
! **********************************************************************************************************************************
! Per-subcase eigenvalue extraction parameters (Feature A: SOL 103 multi-METHOD).
!
! The legacy EIG_* scalars above remain as the "currently active" view used by the eigensolvers (EIG_LANCZOS_ARPACK*,
! EIG_GIV_MGIV, EIG_INV_PWR, etc.). When more than one modes-subcase exists with distinct METHOD requests, LINK4 loops
! over the modes-subcases and copies the per-subcase parameters from EIG_PARAMS(isub) into those scalars at the head of
! each iteration. Results are captured back into the per-subcase arrays at the bottom of the iteration.
!
! For legacy single-METHOD decks, EIG_PARAMS(1) holds the same data as the scalars and the loop runs exactly once.

      TYPE EIG_PARAMS_TYPE
         CHARACTER(LEN=JCARD_LEN)     :: METHOD              = ' '
         CHARACTER(LEN=JCARD_LEN)     :: NORM                = 'MASS'
         CHARACTER(LEN=JCARD_LEN)     :: LAP_MAT_TYPE        = 'DPB'
         CHARACTER(1*BYTE)            :: VECS                = 'Y'
         INTEGER(LONG)                :: SID                 = 0       ! EIGR/EIGRL set ID resolved from CC METHOD
         INTEGER(LONG)                :: N1                  = 0
         INTEGER(LONG)                :: N2                  = 0
         INTEGER(LONG)                :: COMP                = 0
         INTEGER(LONG)                :: GRID                = 0
         INTEGER(LONG)                :: LANCZOS_NEV_DELT    = 2
         INTEGER(LONG)                :: MODE                = 2
         INTEGER(LONG)                :: MSGLVL              = 0
         INTEGER(LONG)                :: NCVFACL             = 3
         REAL(DOUBLE)                 :: CRIT                = ZERO
         REAL(DOUBLE)                 :: FRQ1                = ZERO
         REAL(DOUBLE)                 :: FRQ2                = ZERO
         REAL(DOUBLE)                 :: SIGMA               = -ONE
         ! Results captured from the eigensolver run for this subcase:
         INTEGER(LONG)                :: NUM_EIGENS          = 0       ! No. of eigenvalues actually extracted
         INTEGER(LONG)                :: NVEC                = 0       ! No. of eigenvectors retained
         INTEGER(LONG)                :: NUM_FAIL_CRIT       = 0
         INTEGER(LONG)                :: MIJ_ROW             = 0
         INTEGER(LONG)                :: MIJ_COL             = 0
         REAL(DOUBLE)                 :: MAXMIJ              = ZERO
         ! Per-subcase result buffers, populated by LINK4 inside the modes-subcase loop. They are sized exactly to
         ! NUM_EIGENS at solve time, so memory is not wasted on subcases that ask for fewer modes than the deck maximum.
         REAL(DOUBLE)   , ALLOCATABLE :: EIGEN_VAL(:)        ! (NUM_EIGENS)
         INTEGER(LONG)  , ALLOCATABLE :: MODE_NUM(:)         ! (NUM_EIGENS) extraction order index for each mode
         REAL(DOUBLE)   , ALLOCATABLE :: GEN_MASS(:)         ! (NUM_EIGENS) generalized masses
         REAL(DOUBLE)   , ALLOCATABLE :: EIGEN_VEC(:,:)      ! (NDOFL, NVEC)
         ! Reserved for STATSUB (feature B): static subcase id whose stiffness is used as the prestress source.
         ! Default 0 = legacy BUCKLING single-prestress mechanism.
         INTEGER(LONG)                :: STATSUB_REF         = 0
      END TYPE EIG_PARAMS_TYPE

      TYPE(EIG_PARAMS_TYPE), ALLOCATABLE :: EIG_PARAMS(:)    ! Sized LSUB; index by internal subcase number
      INTEGER(LONG)        , ALLOCATABLE :: CC_EIGR_SID_SUB(:) ! Per-subcase METHOD set ID (0 if subcase has no METHOD).
!                                                              Populated by CC_METH; for subcases lacking a METHOD card,
!                                                              the deck-default (METHOD declared above any SUBCASE) is
!                                                              propagated by the post-parse inheritance pass.
      INTEGER(LONG)                      :: CC_EIGR_SID_DECK = 0
!                                                              METHOD SID declared above the first SUBCASE card (deck-default)
      CHARACTER(1*BYTE) , ALLOCATABLE    :: IS_MODES_SUBCASE(:)
!                                                              'Y' for each modes-subcase (has a resolved METHOD); 'N' otherwise.

      INTEGER(LONG)        , ALLOCATABLE :: CC_STATSUB_SUB(:)   ! Per-subcase STATSUB(PRELOAD) value (0 if subcase has no STATSUB).
!                                                              Populated by CC_STATSUB; for buckling subcases lacking a STATSUB
!                                                              card, the deck-default (STATSUB declared above any SUBCASE) is
!                                                              propagated by the post-parse inheritance pass in LOADC.
      INTEGER(LONG)                      :: CC_STATSUB_DECK  = 0
!                                                              STATSUB(PRELOAD) declared above the first SUBCASE card (deck-default).
      CHARACTER(1*BYTE) , ALLOCATABLE    :: IS_BUCKLING_SUBCASE(:)
!                                                              'Y' for each buckling-subcase (resolved by LOADC for SOL 105
!                                                              decks); 'N' otherwise.

! Per-subcase eigen result tally (allocated alongside EIG_PARAMS by ALLOCATE_MODEL_STUF). The full per-subcase result
! buffers live as allocatable components inside EIG_PARAMS(ISUB).
      INTEGER(LONG)     , ALLOCATABLE :: NUM_EIGENS_SUB(:)   ! No. of eigenvalues extracted per modes-subcase
                                                             
! --- feast_add --- begin !
      CHARACTER(LEN=JCARD_LEN)        :: EIG_EXTRACT_METHOD  = 'ARPACK'
                                                             ! Extract backend selected for EIGRL/LANCZOS family.

      CHARACTER(LEN=JCARD_LEN)        :: EIG_EXTRACT_MODE    = ' '
                                                             ! Optional method-specific mode selector from EIGRL continuation.

      CHARACTER(LEN=JCARD_LEN)        :: EIG_EXTRACT_SOURCE  = 'DEFAULT'
                                                             ! DEFAULT/EIGRL/PARAM source for extract-method selection.

      INTEGER(LONG)                   :: EIG_FEAST_M0        = 48
                                                             ! FEAST search subspace size.

      INTEGER(LONG)                   :: EIG_FEAST_TOL_DIGITS= 8
                                                             ! FEAST convergence digits (fpm(3)).

      INTEGER(LONG)                   :: EIG_FEAST_MAX_LOOP  = 60
                                                             ! FEAST iteration limit (fpm(4)).

      INTEGER(LONG)                   :: EIG_FEAST_N_CONTOUR = 8
                                                             ! FEAST contour integration points (fpm(8)).

      INTEGER(LONG)                   :: EIG_SUBSPACE_NSUB   = 24
                                                             ! Dense inverse-subspace working subspace dimension.

      INTEGER(LONG)                   :: EIG_SUBSPACE_MAX_ITER = 40
                                                             ! Dense inverse-subspace iteration limit.

      INTEGER(LONG)                   :: EIG_DENSE_NEX       = 64
                                                             ! Reserved dense oversampling/workspace knob for parity with other methods.
                                                             
      REAL(DOUBLE)                    :: EIG_FEAST_SEARCH_SCALE = 1.10D0
                                                             ! FEAST range expansion factor when upper frequency is requested.

      REAL(DOUBLE)                    :: EIG_SUBSPACE_TOL    = 1.0D-06
                                                             ! Dense inverse-subspace convergence tolerance.
                                                             
! --- feast_add --- end !

! **********************************************************************************************************************************
! Rigid element ID's

      INTEGER(LONG), ALLOCATABLE      :: RIGID_ELEM_IDS(:)   ! Rigid element ID's

! **********************************************************************************************************************************

      END MODULE MODEL_STUF
