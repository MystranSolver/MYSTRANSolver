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

      MODULE DEBUG_PARAMETERS

! DEBUG is an array and its values are specified on Bulk Data DEBUG entries with field 2 containing the index (DEBUG number) and
! field 3 the value of the DEBUG index. These are used, obviously, for debug output. See the MYSTRAN User's Reference Manual for
! the Bulk Data entry DEBUG.

      USE PENTIUM_II_KIND, ONLY       :  LONG

      IMPLICIT NONE

      INTEGER(LONG), PRIVATE          :: I                      ! Only use is in implied DO loop to initialize DEBUG
      INTEGER(LONG), PARAMETER        :: NDEBUG        =  250   ! Size of DEBUG array
      INTEGER(LONG)                   :: DEBUG(NDEBUG) = (/(0, I=1,NDEBUG)/)

      SAVE

! DEBUG parameter actions:
!        DEBUG   1  = 1 print KIND parameters defined in module PENTIUM_II_KIND to F06 file
!        DEBUG   2  = 1 print constants (parameters) defined in module CONSTANTS_1
!        DEBUG   3  = 1 print machine parameters as determined by LAPACK function DLAMCH
!        DEBUG   4  = 1 do not use BMEANT in QMEM1 calcs
!        DEBUG   5  = 1 print Gauss quadriture formula abscissa's and weights calculated in subr ORDER
!        DEBUG   6  = 1 print some quad element data to BUG file (over and above what is printed with C.C ELDATA command)
!                   = 2 print some hexa element data to BUG file (over and above what is printed with C.C ELDATA command)
!        DEBUG   7  = 1 print arrays ESORT1, ESORT2, EPNT and ETYPE in subr ELESORT before/after sorting elements
!        DEBUG   8  = 1 print grid temperature data in subr TEMPERATURE_DATA_PROC
!                   = 2 print elem temperature data in subr TEMPERATURE_DATA_PROC
!                   = 3 print both grid and elem temperature data in subr TEMPERATURE_DATA_PROC
!        DEBUG   9  > 0 print debug info in subr PINFLG (pin flag processing for elements)

!        DEBUG  10  = 11 or 33 prints data on algorithm to create arrays STFKEY, STFPNT, STFCOL, STF at end of subr ESP
!                   = 12 or 32 prints detailed data on algorithm to create arrays STFKEY, STFPNT, STFCOL, STF in subr SPARSE
!                   = 13 or 33 also prints template of nonzero terms in K in subr ESP if PARAM SETLKTK = 1 or 2
!                   = 21 or 33 prints data on algorithm to create arrays EMSKEY, EMSPNT, EMSCOL, EMS at end of subr ESP
!                   = 22 or 32 prints detailed data on algorithm to create arrays EMSKEY, EMSPNT, EMSCOL, EMS in subr SPARSE

!        DEBUG  11  = 1 print individual 6x6 r.b. displ matrices in basic and global coords for each grid (subr RB_DISP_MATRIX_PROC)
!                   = 2 print ngrid by 6 rigid body displ matrix in global coords for the model
!                   = 3 print both

!        DEBUG  12  = 1 Use BAR K1 and K2 shear factors in calculating stiffness matrix when I12 /= 0 as well as when I12 = 0
!        DEBUG  13  = 1 print GRID_SEQ and INV_GRID_SEQ table
!        DEBUG  14  = 1 print matrices generated in the rigid element generation subr's
!        DEBUG  15  = 1 print concentrated mass data in subrs CONM2_PROC_1 and CONM2_PROC_2
!        DEBUG  16  = 1 use static equivalent instead of work equivalent pressure loads in subrs QPLT1, TPLT2

!        DEBUG  17  > 0 print some diagnostic info in subr KGG_SINGULARITY_PROC for grids that have AUTOSPC'd components
!                   > 1 do above for all grids (not just ones that have AUTOSPC's)

!        DEBUG  18  > 0 print diagnostics in subr QMEM1 regarding checks on the BMEAN matrix satisfying R.B. motion
!        DEBUG  19  = 1 print output from subr STOKEN (which, among other uses, processes Case Control SET cards)

!        DEBUG  20  = 1 bypass the simple solution for GMN when RMM is diagonal. Use subr SOLVE_GMN instead

!        DEBUG  21  = 0 In subr STIFF_MAT_EQUIL_CHK use MATMULT_SFF to multiply stiffness matrix times rigid body displs
!                   = 1 In subr STIFF_MAT_EQUIL_CHK use DSBMV       to multiply stiffness matrix times rigid body displs

!        DEBUG  22  = 1 print RBMAT in subr STIFF_MAT_EQUIL_CHK
!        DEBUG  23  > 0 do equil checks on stiffness matrix even though model has SPOINT's
!        DEBUG  24  = 1 or 3 print KFSe matrix, 2 or 3 print KSSe matrix in subrs REDUCE_N_FS, REDUCE_PN_TO_PF
!        DEBUG  25  = 1 or 3 print PFYS matrix, 2 or 3 print QSYS matrix in subrs REDUCE_N_FS, REDUCE_PN_TO_PF
!        DEBUG  26  = 1 print YSe matrix (S-set enforcorced displs) in REDUCE_N_FS
!        DEBUG  27    Not used
!        DEBUG  28    Not used
!        DEBUG  29    Not used

!        DEBUG  30    Not used
!        DEBUG  31  = 1 print KLL stiff matrix in subr LINK3
!        DEBUG  32  = 1 print PL load matrix in subr LINK3
!        DEBUG  33  = 1 print UL displ matrix
!        DEBUG  34  = 1 or 3 print ABAND matrix (KLL in band form) before equilibrate
!                     2 or 3 print ABAND equil'd in LINK3
!        DEBUG  35    Not used
!        DEBUG  36  = 1 print grid 6x6 mass for every grid in LINK2
!        DEBUG  37    Not used
!        DEBUG  38    Not used
!        DEBUG  39    Not used

!        DEBUG  40  = 1 or 3 print banded stiffness matrix ABAND in subr EIG_GIV_MGIV
!                   = 1 print RFAC = KLL - sigma*MLL             in subr EIG_INV
!                   = 1 print RFAC = KLL - sigma*MLL             in subr EIG_LANCZOS
!                   = 2 or 3 print banded mass      matrix ABAND in subr EIG_GIV_MGIV

!        DEBUG  41  = 1 print KLL stiff matrix in subr LINK4
!        DEBUG  42  = 1 print MLL  mass matrix in subr LINK4
!                   = 2 print MLLn mass matrix in subr LINK4
!        DEBUG  43  = 1 print eigenvectors in LINK4
!        DEBUG  44    Not used
!        DEBUG  45    Not used
!        DEBUG  46  = 1 print debug info for Inverse Power eigen extraction
!        DEBUG  47  = 1 print eigenvalue estimate at each iteration in Lanczos
!        DEBUG  48  = 1 do not calculate off-diag terms in generalized mass matrix in LINK 4
!        DEBUG  49  = 1 print KMSMn in EIG_LANCZOS_ARPACK
!        DEBUG  50  = 1 print debug info in EIG_LANCZOS_ARPACK

!        DEBUG  51  = 1 print debug info in ROT_AXES_MATL_TO_LOC

!        DEBUG  52  = 1 print debug info in FILE_OPEN

!        DEBUG  53  = 1 print debug info in TPLT2

!        DEBUG  54  = 1 print element stiffness matrix to f06 file

!        DEBUG  55  = 1 Write PHIXG in full format in EXPAND_PHIXA_TO_PHIXG
!                   = 2 Write PHIZG in full format in LINK5
!                   = 3 Write both

!        DEBUG  56  = 1 print PEL in subr CALC_ELEM_NODE_FORCES.f03

!        DEBUG  57  > 0 print KGG full in LINK2.f03

!        DEBUG  64    Not used
!        DEBUG  65    Not used
!        DEBUG  66    Not used
!        DEBUG  67    Not used
!        DEBUG  80  > 0 print LAPACK_S scale factors, in subr EQUILIBRATE, used to equilibrate the stiffness matrices

!        DEBUG  82  = 1 print data from subr MATMULT_SFF

!        DEBUG  83  = 1 print data from subr MATMULT_SFS_NTERM
!                     2 print data from subr MATMULT_SFS
!                     3 print data from both subrs

!        DEBUG  84  = 1 print data from subr MATMULT_SSS_NTERM
!                     2 print data from subr MATMULT_SSS
!                     3 print data from both subrs

!        DEBUG  85  = 1 print data from subr MATTRNSP_SS

!        DEBUG  86  = 1 print data from subr PARTITION_SS_NTERM
!                     2 print data from subr PARTITION_SS
!                     3 print data from both subrs

!        DEBUG  87  = 1 print data on algorithm to convert sparse CRS matrix to sparse CCS in subr SPARSE_CRS_SPARSE_CCS
!        DEBUG  88  = 1 do not write separator line between grids several places(matrix diagonal output, equil check)
!        DEBUG  89  = 1 write row numbers where there are zero diag terms in subr SPARSE_MAT_DIAG_ZEROS

!        DEBUG  91  = 1 print info on how the number of output requests are counted to get MAXREQ (used in dimensioning array OGEL)
!        DEBUG  92  = 1 print OLOAD, SPCF, MPCF totals even if ALL_SAME_CID = 'N'

!        DEBUG 100  > 0 check allocation status of allocatable arrays
!                   > 1 also write memory allocated to all arrays to F06 file

!        DEBUG 101  > 0 write sparse I_MATOUT array in subr READ_MATRIX_1
!                   > 1 call subr to check I_MATOUT array to make sure that terms are nondecreasing

!        DEBUG 102  > 0 print debug info in subr MERGE_MAT_COLS_SSS

!        DEBUG 103  > 0 Do not use MRL (or MLR) in calc of modal participation factors and effective mass

!        DEBUG 104  > 0 check if KRRcb is singular - only done with a BANDED solution (but KRRcb is generally small)

!        DEBUG 105  > 0 write KLLs matrix to unformatted file

!        DEBUG 106  > 0 write info on all files in subr WRITE_ALLOC_MEM_TABLE (if 0 only write for those arrays that have memory
!                       allocated to them

!        DEBUG 108  > 0 write EDAT table

!        DEBUG 109  > 0 write debug info in subr ELMDIS

!        DEBUG 110  > 0 write debug info for BUSH elem in subrs ELMDAT1, ELMGM1

!        DEBUG 111  > 0 write some debug info on RSPLINE

!        DEBUG 112  > 0 write THETAM (plate element material angle) and the location in subr EMG where it was calculated

!        DEBUG 113  > 0 write PBARL entries in a special format that has 1 line per PBAR entry

!        DEBUG 114  > 0 write debug info in subr OU4_PARTVEC_PROC

!        DEBUG 115  > 0 write debug info in subr READ_INCLUDE_FILNAM

!        DEBUG 172  > 0 calc PHI_SQ for the MIN4T based on area weighting of the TRIA3's. Otherwise, use simple average

!        DEBUG 173  = 1 write some debug info in subr PARSE_CSV_STRING
!                   = 2 write some more detailed data

!        DEBUG 174  > 0 print MPFACTOR, MEFFMASS values with 2 decimal places of accuracy rather than 6

!        DEBUG 175  > 0 Debug output from subr SURFACE_FIT

!        DEBUG 176  > 0 calc stresses using SEi and elem displacements instead of from STRAIN and material matrix

!        DEBUG 177  > 0 print BAR, ROD margins of safety whether or not they would otherwise be

!        DEBUG 178  > 0 print G-set differential stiffness matrix

!        DEBUG 179  = 1 print blank space at beg of lines of output for CUSERIN entries in the F06 file

!        DEBUG 180  > 0 write debug info to F06 for USERIN elements

!        DEBUG 181  = 1 include USERIN RB mass in subr GPWG even though user did not input 3rd matrix (RBM0) on IN4FIL

!        DEBUG 182  = 1 print debug data in subr MGGS_MASS_MATRIX for scalar mass matrix

!        DEBUG 183  = 1 write some debug data for generating TDOF array

!        DEBUG 184  > 0 write L1M data to F06

!        DEBUG 185  > 0 Let eigen routines find and process all eigenval, vecs found even if NVEC > NDOFL - NUM_MLL_DIAG_ZEROS

!        DEBUG 186  > 0 print debug info for pressure loads on faces of solid elements

!        DEBUG 187  > 0 print list of number of the various elastic elements

!        DEBUG 188  > 0 Do not abort in QPLT3 if KOO is reported to be singular

!        DEBUG 189  = 1 print messages in subr ESP for KE in local coords if element diagonal stiffness < 0
!                   = 2 print these messages in subr ESP after transformation to global
!                   = 3 do both

!        DEBUG 190  > 0 do not round off FAILURE_INDEX to 0 in subr POLY_FAILURE_INDEX

!        DEBUG 191  = 0 Use temperatures at Gauss points for thermal loads in solid elements

!        DEBUG 192  > 0 print some summary info for max abs value of GP force balance for each soln vector

!        DEBUG 193  = 1  , call FILE_INQUIRE at end of LINK1
!                   = 2  , call FILE_INQUIRE at end of LINK2
!                   = 3  , call FILE_INQUIRE at end of LINK3
!                   = 4  , call FILE_INQUIRE at end of LINK4
!                   = 5  , call FILE_INQUIRE at end of LINK5
!                   = 6  , call FILE_INQUIRE at end of LINK6
!                   = 9  , call FILE_INQUIRE at end of LINK9
!                   = 100, call FILE_INQUIRE at end of MAIN
!                   = 999, do all of the above

!        DEBUG 194  = 1 or 3 skip check on CW/CCW numbering of QUAD's
!                   = 2 or 3 skip check on QUAD interior angles < 180 deg
!                   = 3 skip both

!        DEBUG 195  > 0 Write OTM matrices to F06 at end of LINK9

!        DEBUG 196  = 0 Matrix output filter SMALL = EPSIL(1)
!                   > 0 Matrix output filter SMALL = TINY (param defined by user with default = 0.D0)

!        DEBUG 197  > 0 print debug info in subr EC_ENTRY_OUTPUT4

!        DEBUG 198  > 0 write debug info in subr QPLT3

!        DEBUG 199  > 0 check matrix times its inverse = identity matrix in

!        DEBUG 201 /= 0 allow SOL = BUCKLING or DIFFEREN to run even if some elements are not coded for these soln's

!        DEBUG 202  > 0 calculate RB and constant strain sanity checks on strain-displacement matrices

!        DEBUG 203  > 0 print debug info in subr BAR1

!        DEBUG 204  = 1, debug BUSH for calculation of the offsets in local element coords
!                   = 2, debug BUSH for calculation of the stiffness and stress/strain rcovery matrices
!                   = 9, do all of above in BUSH

!        DEBUG 248  > 0 override fatal error and continue with orthotropic material properties for MIN4T QUAD4

!        DEBUG 249  > 0 In BREL1 call code for Timoshenko BART instead of BAR1

      END MODULE DEBUG_PARAMETERS
