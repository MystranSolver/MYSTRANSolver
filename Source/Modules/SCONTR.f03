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

      MODULE SCONTR

! System control variables

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, SINGLE, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), PRIVATE   :: I                         ! Only use is in implied DO loops below to initialize variables

      INTEGER(LONG), PARAMETER :: NCCCD                  = 31                ! Max number of words in array CC_CMD_DESCRIBERS

! System control / counters

      CHARACTER(  7*BYTE), PARAMETER :: PROG_NAME        = 'MYSTRAN'         ! This computer program's name
      CHARACTER( 31*BYTE), PARAMETER :: BLNK_SUB_NAM     = ' '               ! Blank subr name

      CHARACTER(  1*BYTE)            :: BANDIT_IMPROVE   = 'N'               ! Set to 'Y' if Bandit improves bandwidth
      CHARACTER(  1*BYTE)            :: BARTOR           = 'N'               ! Set to 'Y' if any BAR elem has tor stress coeff
      CHARACTER(  1*BYTE)            :: BEAMTOR          = 'N'               ! Set to 'Y' if any BEAM elem has tor stress coeff
      CHARACTER(  8*BYTE)            :: CC_CMD_DESCRIBERS(NCCCD)             ! Array of words from data in between () on CC entry
      CHARACTER(  1*BYTE)            :: COMM(0:49)       = (/('0', I=0,49)/) ! Flag (0 or C) to indicate a LINK has been completed
      CHARACTER(  1*BYTE)            :: CHKPNT           = 'N'               ! 'Y' if run is to be check pointed
      CHARACTER(  6*BYTE)            :: ECHO             = 'UNSORT'          ! Case Control ECHO (either UNSORT, SORT or NONE)
      CHARACTER(  1*BYTE)            :: ELESORT_RUN      = 'N'               ! 'Y'/'N' indicator of whether subr ELESORT has run
      CHARACTER(  1*BYTE)            :: EPSIL1_SET       = 'N'               ! 'Y'/'N' indicator of whether EPSIL(1) is in B.D. deck
      CHARACTER(  6*BYTE)            :: ENFORCED         = 'N'               ! If 'Y' then this is a run where all DOF's are SE-set
      CHARACTER(  5*BYTE)            :: FACTORED_MATRIX  = 'none '           ! Name of matrix that has been decomposed into
!                                                                              triangular factors (e.g.'KLL  ')
      CHARACTER(  1*BYTE)            :: IERRFL(10)       = (/('N', I=1,10)/) ! 'Y', 'N' indicates errs in Bulk Data card fields
      CHARACTER(  1*BYTE)            :: IMB_BLANK(2:9)   = (/('N', I=2, 9)/) ! 'Y', 'N' indicates imbedded blanks in B.D. field
      CHARACTER(  1*BYTE)            :: PRINTENV         = 'N'               ! 'Y' if Software Passport env vars are to be printed
      CHARACTER(  1*BYTE)            :: RESTART          = 'N'               ! 'Y' if run is a restart
! --- rsa_nastran begin --- !
      CHARACTER(  1*BYTE)            :: RSA_NX_SEMODES   = 'N'               ! 'Y' if SOL SEMODES alias was requested
! --- rsa_nastran end --- !
      CHARACTER( 16*BYTE)            :: SOL_NAME         = '                '! Name for the solution (e.g. 'STATICS')
      CHARACTER( 2*BYTE)             :: TSET_CHR_LEN     = '  '              ! Char len of entries in TSET

      CHARACTER(  1*BYTE)            :: UNLOCK           = 'N'               ! Y or N indicator for whether program will run
!                                                                              unlimited size problems
      INTEGER(LONG)            :: IBIT(0:15)          = (/(2**I, I=0,15)/)
      INTEGER(LONG)            :: JF(10)              = (/(   I, I=1,10)/)

      INTEGER(LONG), PARAMETER :: BD_ENTRY_LEN        = 256      ! Length, in chars, of a Bulk Data deck entry
      INTEGER(LONG), PARAMETER :: CC_ENTRY_LEN        = 256      ! Length, in chars, of a Case Control deck entry
      INTEGER(LONG), PARAMETER :: DATA_NAM_LEN        = 256      ! Length, in chars, of a data set name
      INTEGER(LONG), PARAMETER :: EC_ENTRY_LEN        = 256      ! Length, in chars, of a Executive Control deck entry
      INTEGER(LONG), PARAMETER :: INI_ENTRY_LEN       = 256      ! Length, in chars, of an entry in file Program.INI
      INTEGER(LONG), PARAMETER :: JCARD_LEN           =  16      ! Length, in chars, of an entry in the B.D. file (after expanding)
      INTEGER(LONG), PARAMETER :: MAX_INTEGER_LEN     =  10      ! Max number of digits in any integer to keep integers <= 4 bytes
      INTEGER(LONG), PARAMETER :: NUM_TRIA_ORDERS     =   3      ! Num of triangular integration orders for isoparametric elems

      INTEGER(LONG)            :: BANDIT_ERR          =   0      ! Error report from Bandit
      INTEGER(LONG)            :: DEMO_GRID_LIMIT     =  10      ! Max number of grids allowed in demo ver of the program
      INTEGER(LONG)            :: FATAL_ERR           =   0      ! Fatal err count from input (Exec & Case Control, Bulk Data)
      INTEGER(LONG)            :: WARN_ERR            =   0      ! Warn err count from input (Exec & Case Control, Bulk Data)

      INTEGER(LONG)            :: INT_SC_NUM          =   0      ! Internal subcase number
      INTEGER(LONG)            :: INT_EIG_NUM         =   0      ! Per-subcase local eigenvector index (resets to 1 at each subcase)
      INTEGER(LONG)            :: JTSUB               =   0      ! Internal thermal array col no. corresponding to INT_SC_NUM

      INTEGER(LONG)            :: KLL_SDIA            =   0      ! Number of super-diagonals in matrix KLL
      INTEGER(LONG)            :: KLLD_SDIA           =   0      ! Number of super-diagonals in matrix KLLD
      INTEGER(LONG)            :: KMSM_SDIA           =   0      ! Number of super-diagonals in matrix KMSM
      INTEGER(LONG)            :: KOO_SDIA            =   0      ! Number of super-diagonals in matrix KOO
      INTEGER(LONG)            :: KRRcb_SDIA          =   0      ! Number of super-diagonals in matrix KRRcb

      INTEGER(LONG)            :: KMAT_BW             =   0      ! 6 times grid BW returned from subr BANDIT called in subr LINK1

      INTEGER(LONG)            :: LINKNO              =   0      ! LINK num being run
      INTEGER(LONG)            :: LINKNO_L1A          =   0      ! LINK num MYSTRAN was executing when file LINK1A was written
      INTEGER(LONG)            :: LINKNO_START        =   0      ! LINK num to start MYSTRAN with. Normally this is LINK 1
!                                                                  however, if the program was terminated abnormally right as
!                                                                  a new LINK is beginning, this allows a "restart"

      INTEGER(LONG)            :: LBAROFF             =   0      ! Max allow num of CBAR, CBEAM offset vectors      (see note (1))
      INTEGER(LONG)            :: LBUSHOFF            =   0      ! Max allow num of CBUSH       offset vectors      (see note (1))
      INTEGER(LONG)            :: LCMASS              =   0      ! Max allow num of CMASSi Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LCONM2              =   0      ! Max allow num of CONM2  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LCORD               =   0      ! Max allow num of CORD   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LDOFG               =   0      ! Max allow num of G-set DOF's (checked against NDOFG)
      INTEGER(LONG)            :: LEDAT               =   0      ! Max allow num of element connection data  (see note (1))
      INTEGER(LONG)            :: LELE                =   0      ! Max allow num of elems                 (see note (1))
      INTEGER(LONG)            :: LFORCE              =   0      ! Max allow num of FORCE   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LGRAV               =   0      ! Max allow num of GRAV    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LGRID               =   0      ! Max allow num of GRID    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LGUSERIN            =   0      ! Max no. of grids found on any USERIN entry (not conting SPOINT's)
      INTEGER(LONG)            :: LIND_GRDS_MPCS      =   0      ! Max allow num of independent grids on MPC's and rigid elems
      INTEGER(LONG)            :: LLOADC              =   0      ! Max no. of pairs of (load fac/load mag) over all LOAD B.D. cards
!                                                                  incl the pair defined by the LOAD card set ID and overall scale
!                                                                  factor. Counted  by subr LOADB0 based on info from subr BD_LOAD0
      INTEGER(LONG)            :: LLOADR              =   0      ! Max allow num of LOAD    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LMATANGLE           =   0      ! Max allow num of plate elem matl angles    (see note (1))
      INTEGER(LONG)            :: LMATL               =   0      ! Max allow num of matl    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LMPC                =   0      ! Max allow num of MPC's                     (see note (1))
      INTEGER(LONG)            :: LMPCADDC            =   0      ! Max no. of pairs of set ID's over all MPCADD B.D. cards incl
!                                                                  the MPCADD set ID on the MPCADD card
      INTEGER(LONG)            :: LMPCADDR            =   0      ! Max allow num of MPCADD  Bulk Data cards   (see Note (1))

      INTEGER(LONG)            :: LPBAR               =   0      ! Max allow num of PBAR    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPBEAM              =   0      ! Max allow num of PBEAM   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPBUSH              =   0      ! Max allow num of PBUSH   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPCOMP              =   0      ! Max allow num of PCOMP   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPCOMP_PLIES        =   0      ! Max allow num of plies for any PCOMP entry
      INTEGER(LONG)            :: LPDAT               =   0      ! Max allow num of rows for array PDATA
      INTEGER(LONG)            :: LPELAS              =   0      ! Max allow num of PELAS   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPLATEOFF           =   0      ! Max allow num of plate elem offsets        (see note (1))
      INTEGER(LONG)            :: LPLATETHICK         =   0      ! Max allow num of plate elem thicknesses    (see note (1))
      INTEGER(LONG)            :: LPLOAD              =   0      ! Max allow num of PLOADi  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPLOTEL             =   0      ! Max allow num of PLOTEL  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPMASS              =   0      ! Max allow num of PMASS   Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPROD               =   0      ! Max allow num of PROD    Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPSHEAR             =   0      ! Max allow num of PSHEAR  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPSHEL              =   0      ! Max allow num of PSHELL  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPSOLID             =   0      ! Max allow num of PSOLID  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPUSER1             =   0      ! Max allow num of PUSER1  Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LPUSERIN            =   0      ! Max allow num of PUSERIN Bulk Data cards   (see note (1))
      INTEGER(LONG)            :: LRFORCE             =   0      ! Max allow num of RFORCE  Bulk data entries (see note (1))
      INTEGER(LONG)            :: LRIGEL              =   0      ! Max allow num of rigid elems           (see note (1))
      INTEGER(LONG)            :: LSEQ                =   0      ! Max allow num of rows for arrays SEQ1, SEQ2
      INTEGER(LONG)            :: LSETLN              =   0      ! Max allow num of chars in SET definitions (see Note (1))
      INTEGER(LONG)            :: LSETS               =   0      ! Max allow num of Case Control SET's       (see Note (1))
      INTEGER(LONG)            :: LSLOAD              =   0      ! Max allow num of SLOAD pairs              (see Note (1))
      INTEGER(LONG)            :: LSNORM              =   0      ! Max allow num of SNORM Bulk Data cards    (see Note (1))
      INTEGER(LONG)            :: LSPC                =   0      ! Max allow num of SPC    Bulk Data cards   (see Note (1))
      INTEGER(LONG)            :: LSPC1               =   0      ! Max allow num of SPC1   Bulk Data cards   (see Note (1))
      INTEGER(LONG)            :: LSPCADDC            =   0      ! Max no. of pairs of set ID's over all SPCADD B.D. cards incl
!                                                                  the SPCADD set ID on the SPCADD card
      INTEGER(LONG)            :: LSPCADDR            =   0      ! Max allow num of SPCADD Bulk Data cards   (see Note (1))
      INTEGER(LONG)            :: LSUB                =   0      ! Max allow num of subcases                 (see Note (1))
      INTEGER(LONG)            :: LSUSERIN            =   0      ! Max no. of SPOINT's found on any USERIN entry
      INTEGER(LONG)            :: LTDAT               =   0      ! Max allow num of rows for arrays TDATA
      INTEGER(LONG)            :: LTERM_KGG           =   0      ! Max allow num of terms in matrix KGG      (see note (1))
      INTEGER(LONG)            :: LTERM_KGGD          =   0      ! Max allow num of terms in matrix KGGD     (see note (1))
      INTEGER(LONG)            :: LTERM_MGGE          =   0      ! Max allow num of terms in matrix MGG for elems (see note (1))
      INTEGER(LONG)            :: LVVEC               =   0      ! Max allow num of v vectors (CBAR)         (see Note (1))

      INTEGER(LONG)            :: MAX_ELEM_DEGREE     =   0      ! Max number of elements connected to any one grid
      INTEGER(LONG)            :: MAX_GAUSS_POINTS    =   0      ! Max number of Gauss pts for any element
      INTEGER(LONG)            :: MAX_STRESS_POINTS   =   0      ! Max number of pts for stress recovery within 1 elem (for output
!                                                                  of STRESS at CORNER, GAUSS, CENTER options)
      INTEGER(LONG)            :: MDT                 =   0      ! Max number of rows in dimensioning element array DT
      INTEGER(LONG)            :: MELGP               =   0      ! Max no. grid points for any finite element
      INTEGER(LONG)            :: MELDOF              =   0      ! Max no. DOF's for any finite element
      INTEGER(LONG)            :: MID1_PCOMP_EQ       = 1000000  ! MAT2 MID1 matl ID base num for PSHELL equiv mtrl's for PCOMP
      INTEGER(LONG)            :: MID2_PCOMP_EQ       = 2000000  ! MAT2 MID2 matl ID base num for PSHELL equiv mtrl's for PCOMP
      INTEGER(LONG)            :: MID3_PCOMP_EQ       = 3000000  ! MAT2 MID3 matl ID base num for PSHELL equiv mtrl's for PCOMP
      INTEGER(LONG)            :: MID4_PCOMP_EQ       = 4000000  ! MAT2 MID4 matl ID base num for PSHELL equiv mtrl's for PCOMP
      INTEGER(LONG)            :: MIN4T_QUAD4_TRIA_NO =   0      ! The num (1 thru 4) of the TRIA3 elem that makes up a MIN4T QUAD4
      INTEGER(LONG)            :: MLL_SDIA            =   0      ! Number of super-diagonals in matrix MLL
      INTEGER(LONG)            :: MMPC                =   0      ! Max no. of grid/comp/coeff triplets on any one MPC  logical card
      INTEGER(LONG)            :: MOFFSET             =   0      ! No. rows allowed in dim array OFFDIS, OFFSET. Set in subr LOADB
      INTEGER(LONG)            :: MRBE3               =   0      ! Max no. of grid/comp/coeff triplets on any one RBE3 logical card
      INTEGER(LONG)            :: MRSPLINE            =   0      ! Max no. of non-blank fields on any one RSPLINE logical card

      INTEGER(LONG)            :: NAOCARD             =   0      ! Count of no. of ASET/OMIT cards written to filename.L1N
      INTEGER(LONG)            :: NBAROFF             =   0      ! Count of CBAR  offset vectors as Bulk Data is read
      INTEGER(LONG)            :: NBUSHOFF            =   0      ! Count of CBUSH offset vectors as Bulk Data is read
      INTEGER(LONG)            :: NBAROR              =   0      ! Count of no. of BAROR  Bulk Data cards
      INTEGER(LONG)            :: NBEAMOR             =   0      ! Count of no. of BEAMOR Bulk Data cards
      INTEGER(LONG)            :: NCBAR               =   0      ! Count of no. of CBAR    elems
      INTEGER(LONG)            :: NCBEAM              =   0      ! Count of no. of CBEAM   elems
      INTEGER(LONG)            :: NCBUSH              =   0      ! Count of no. of CBUSH   elems
      INTEGER(LONG)            :: NCELAS1             =   0      ! Count of no. of CELAS1  elems
      INTEGER(LONG)            :: NCELAS2             =   0      ! Count of no. of CELAS2  elems
      INTEGER(LONG)            :: NCELAS3             =   0      ! Count of no. of CELAS3  elems
      INTEGER(LONG)            :: NCELAS4             =   0      ! Count of no. of CELAS4  elems
      INTEGER(LONG)            :: NC_INFILE           =   0      ! Length of INFILE, incl directory
      INTEGER(LONG)            :: NCHEXA8             =   0      ! Count of no. of CHEXA   elems with  8 nodes
      INTEGER(LONG)            :: NCHEXA20            =   0      ! Count of no. of CHEXA   elems with 20 nodes
      INTEGER(LONG)            :: NCMASS              =   0      ! Count of no. of CMASSi  elems
      INTEGER(LONG)            :: NCONM2              =   0      ! Count of no. of CONM2   elems
      INTEGER(LONG)            :: NCORD               =   0      ! Count of no. of CORD cards
      INTEGER(LONG)            :: NCORD1              =   0      ! Count of no. of CORD1R
      INTEGER(LONG)            :: NCORD2              =   0      ! Count of no. of CORD2C, 2R or 2S cards
      INTEGER(LONG)            :: NCPENTA6            =   0      ! Count of no. of CPENTA  elems with  6 nodes
      INTEGER(LONG)            :: NCPENTA15           =   0      ! Count of no. of CPENTA  elems with 15 nodes
      INTEGER(LONG)            :: NCPYRA5             =   0      ! Count of no. of CPYRA   elems with  5 nodes
      INTEGER(LONG)            :: NCPYRA14            =   0      ! Count of no. of CPYRA   elems with 14 nodes
      INTEGER(LONG)            :: NCQUAD4             =   0      ! Count of no. of CQUAD4  elems
      INTEGER(LONG)            :: NCQUAD4K            =   0      ! Count of no. of CQUAD4K elems
      INTEGER(LONG)            :: NCQUAD8             =   0      ! Count of no. of CQUAD8  elems
      INTEGER(LONG)            :: NCQUADR             =   0      ! Count of no. of CQUADR  elems
      INTEGER(LONG)            :: NCROD               =   0      ! Count of no. of CROD    elems
      INTEGER(LONG)            :: NCSHEAR             =   0      ! Count of no. of CSHEAR  elems
      INTEGER(LONG)            :: NCTETRA4            =   0      ! Count of no. of CTETRA  elems with  4 nodes
      INTEGER(LONG)            :: NCTETRA10           =   0      ! Count of no. of CTETRA  elems with 10 nodes
      INTEGER(LONG)            :: NCTRIA3             =   0      ! Count of no. of CTRIA3  elems
      INTEGER(LONG)            :: NCTRIA3K            =   0      ! Count of no. of CTRIA3K elems
      INTEGER(LONG)            :: NCUSER1             =   0      ! Count of no. of CUSER1  elems
      INTEGER(LONG)            :: NCUSERIN            =   0      ! Count of no. of CUSERIN elems

      INTEGER(LONG)            :: NDOFA               =   0      ! Count of no. of DOF's in the A-set
      INTEGER(LONG)            :: NDOF_EIG            =   0      ! DOF size to replace NDOFL for eigprob to avoid zero mass modes
      INTEGER(LONG)            :: NDOFF               =   0      ! Count of no. of DOF's in the F-set
      INTEGER(LONG)            :: NDOFG               =   0      ! Count of no. of DOF's in the G-set
      INTEGER(LONG)            :: NDOFL               =   0      ! Count of no. of DOF's in the L-set
      INTEGER(LONG)            :: NDOFM               =   0      ! Count of no. of DOF's in the M-set
      INTEGER(LONG)            :: NDOFN               =   0      ! Count of no. of DOF's in the N-set
      INTEGER(LONG)            :: NDOFO               =   0      ! Count of no. of DOF's in the O-set
      INTEGER(LONG)            :: NDOFR               =   0      ! Count of no. of DOF's in the R-set
      INTEGER(LONG)            :: NDOFS               =   0      ! Count of no. of DOF's in the S-set
      INTEGER(LONG)            :: NDOFSA              =   0      ! Count of no. of DOF's in the S-set constr to 0 via AUTOSPC
      INTEGER(LONG)            :: NDOFSB              =   0      ! Count of no. of DOF's in the S-set constr to 0 on SPC/SPC1 cards
      INTEGER(LONG)            :: NDOFSE              =   0      ! Count of no. of DOF's in the S-set with enforced displacement
      INTEGER(LONG)            :: NDOFSG              =   0      ! Count of no. of DOF's in the S-set perm SPC'd on GRID cards
      INTEGER(LONG)            :: NDOFSZ              =   0      ! NDOFSB + NDOFSG Count of no. of DOF's in the S-set with 0 displ

      INTEGER(LONG)            :: NEDAT               =   0      ! Count of no. of terms in EDAT array
      INTEGER(LONG)            :: NELE                =   0      ! Count of no. of elastic elems
      INTEGER(LONG)            :: NFORCE              =   0      ! Count of no. of FORCE/MOMENT Bulk Data cards
      INTEGER(LONG)            :: NGRAV               =   0      ! Count of no. of GRAV    Bulk Data cards
      INTEGER(LONG)            :: NGRDSET             =   0      ! Count of no. of GRDSET  Bulk Data cards
      INTEGER(LONG)            :: NGRID               =   0      ! Count of no. of GRID    Bulk Data cards
      INTEGER(LONG)            :: NIND_GRDS_MPCS      =   0      ! Count of no. of independent grids on MPC's and rigid elems
      INTEGER(LONG)            :: NLOAD               =   0      ! Count of no. of LOAD    Bulk Data cards
      INTEGER(LONG)            :: NMATANGLE           =   0      ! Count of no. of matl property angles on plate elem conn entries
      INTEGER(LONG)            :: NMATL               =   0      ! Count of no. of matl    Bulk Data cards
      INTEGER(LONG)            :: NMPC                =   0      ! Count of no. of MPC     Bulk Data cards
      INTEGER(LONG)            :: NMPCADD             =   0      ! Count of no. of MPCADD cards
      INTEGER(LONG)            :: NPBAR               =   0      ! Count of no. of PBAR    Bulk Data cards
      INTEGER(LONG)            :: NPBARL              =   0      ! Count of no. of PBARL   Bulk Data cards
      INTEGER(LONG)            :: NPBEAM              =   0      ! Count of no. of PBEAM   Bulk Data cards
      INTEGER(LONG)            :: NPBUSH              =   0      ! Count of no. of PBUSH   Bulk Data cards
      INTEGER(LONG)            :: NPCARD              =   0      ! Count of no. of PLOAD1, PLOAD2 cards written to filename.L1Q
      INTEGER(LONG)            :: NPCOMP              =   0      ! Count of no. of PCOMP   Bulk Data cards
      INTEGER(LONG)            :: NPDAT               =   0      ! Count of no. of rows that go into array PDATA
      INTEGER(LONG)            :: NPELAS              =   0      ! Count of no. of PELAS  Bulk Data cards
      INTEGER(LONG)            :: NPLATEOFF           =   0      ! Count of no. of plate element offsets on plate elem conn entries
      INTEGER(LONG)            :: NPLATETHICK         =   0      ! Count of no. of plate thicknesses on plate elem conn entries
      INTEGER(LONG)            :: NPLOTEL             =   0      ! Count of no. of PLOTEL
      INTEGER(LONG)            :: NPLOAD              =   0      ! Count of no. of PLOADi  Bulk Data cards
      INTEGER(LONG)            :: NPLOAD4_3D          =   0      ! Count of no. of PLOAD4 entries that are for solid elements
      INTEGER(LONG)            :: NPMASS              =   0      ! Count of no. of PMASS  Bulk Data cards
      INTEGER(LONG)            :: NPROD               =   0      ! Count of no. of PROD   Bulk Data cards
      INTEGER(LONG)            :: NPSHEAR             =   0      ! Count of no. of PSHEAR Bulk Data cards
      INTEGER(LONG)            :: NPSHEL              =   0      ! Count of no. of PSHELL Bulk Data cards
      INTEGER(LONG)            :: NPSOLID             =   0      ! Count of no. of PSOLID Bulk Data cards
      INTEGER(LONG)            :: NPUSER1             =   0      ! Count of no. of PUSER1 Bulk Data cards
      INTEGER(LONG)            :: NPUSERIN            =   0      ! Count of no. of PUSER1 Bulk Data cards
      INTEGER(LONG)            :: NRBAR               =   0      ! Count of no. of RBAR rigid elems
      INTEGER(LONG)            :: NRBE1               =   0      ! Count of no. of RBE1 rigid elems
      INTEGER(LONG)            :: NRBE2               =   0      ! Count of no. of RBE2 rigid elems
      INTEGER(LONG)            :: NRFORCE             =   0      ! Count of no. of RFORCE  Bulk Data cards
      INTEGER(LONG)            :: NRIGEL              =   0      ! Count of no. of rigid elems (RBAR+RBE1+RBE2+RBE3)
      INTEGER(LONG)            :: NRECARD             =   0      ! Count of no. of rigid element records written to file L1F
      INTEGER(LONG)            :: NROWS_OTM_ACCE      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_DISP      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_MPCF      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_SPCF      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_ELFE      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_ELFN      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_STRE      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_OTM_STRN      =   0      ! Num of rows for OTM
      INTEGER(LONG)            :: NROWS_TXT_ACCE      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_DISP      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_MPCF      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_SPCF      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_ELFE      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_ELFN      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_STRE      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NROWS_TXT_STRN      =   0      ! Num of rows for TXT
      INTEGER(LONG)            :: NRSPLINE            =   0      ! Count of no. of RSPLINE elements
      INTEGER(LONG)            :: NSEQ                =   0      ! Count of no. of grid points that are on SEQGP Bulk Data cards
      INTEGER(LONG)            :: NSETS               =   0      ! Count of no. of SET's in Case Control
      INTEGER(LONG)            :: NSLOAD              =   0      ! Count of no. of SLOAD pairs (point/mag)
      INTEGER(LONG)            :: NSNORM              =   0      ! Count of no. of SNORM Bulk Data cards
      INTEGER(LONG)            :: NSPC                =   0      ! Count of no. of SPC  cards written to filename.L1O
      INTEGER(LONG)            :: NSPC1               =   0      ! Count of no. of SPC1 cards written to filename.L1O
      INTEGER(LONG)            :: NSPCADD             =   0      ! Count of no. of SPCADD cards
      INTEGER(LONG)            :: NSPOINT             =   0      ! Count of no. of SPOINT's
      INTEGER(LONG)            :: NUM_SPCSIDS         =   0      ! The number of SPC set ID's called for in an execution
      INTEGER(LONG)            :: NSUB                =   0      ! Count of no. of subcases
      INTEGER(LONG)            :: NUM_MODES_SUBS      =   0      ! Count of subcases that have a resolved METHOD (modes-subcases)
      INTEGER(LONG)            :: NUM_BUCKLING_SUBS   =   0      ! Count of subcases tagged as buckling-subcases (SOL 105 only)
      INTEGER(LONG)            :: TOTAL_MODES         =   0      ! sum(NUM_EIGENS_SUB(:)) = total eigenvectors across all modes-subs
      INTEGER(LONG), ALLOCATABLE :: MODE_SUBCASE(:)              ! Length = TOTAL_MODES. Maps each global mode index (the JVEC
!                                                                  iteration in LINK5/LINK9) to its owning internal subcase index.
!                                                                  For legacy single-METHOD SOL 103 decks this is just
!                                                                  MODE_SUBCASE(:) = 1.
      INTEGER(LONG)            :: NTCARD              =   0      ! Count of no. of TEMP/TEMPRB/TEMPP1 cards written to filename.L1K
      INTEGER(LONG)            :: NTDAT               =   0      ! Count of no. of rows that go into array TDATA
      INTEGER(LONG)            :: NTERM_ALL           =   0      ! Count of no. of terms in ALL             matrix
      INTEGER(LONG)            :: NTERM_CG_LTM        =   0      ! Count of no. of terms in CG_LDS_LTMcb    matrix
      INTEGER(LONG)            :: NTERM_DLR           =   0      ! Count of no. of terms in DLR             matrix
      INTEGER(LONG)            :: NTERM_GMN           =   0      ! Count of no. of terms in GMN             matrix
      INTEGER(LONG)            :: NTERM_GOA           =   0      ! Count of no. of terms in GOA             matrix
      INTEGER(LONG)            :: NTERM_HMN           =   0      ! Count of no. of terms in HMN             matrix
      INTEGER(LONG)            :: NTERM_IF_LTM        =   0      ! Count of no. of terms in IF_LTM          matrix
      INTEGER(LONG)            :: NTERM_IRR           =   0      ! Count of no. of terms in IRR             matrix
      INTEGER(LONG)            :: NTERM_KAA           =   0      ! Count of no. of terms in KAA             matrix
      INTEGER(LONG)            :: NTERM_KAAD          =   0      ! Count of no. of terms in KAAD            matrix
      INTEGER(LONG)            :: NTERM_KAO           =   0      ! Count of no. of terms in KAO             matrix
      INTEGER(LONG)            :: NTERM_KAOD          =   0      ! Count of no. of terms in KAOD            matrix
      INTEGER(LONG)            :: NTERM_KFF           =   0      ! Count of no. of terms in KFF             matrix
      INTEGER(LONG)            :: NTERM_KFFD          =   0      ! Count of no. of terms in KFFD            matrix
      INTEGER(LONG)            :: NTERM_KFS           =   0      ! Count of no. of terms in KFS             matrix
      INTEGER(LONG)            :: NTERM_KFSD          =   0      ! Count of no. of terms in KFSD            matrix
      INTEGER(LONG)            :: NTERM_KFSe          =   0      ! Count of no. of terms in KFSe            matrix
      INTEGER(LONG)            :: NTERM_KFSDe         =   0      ! Count of no. of terms in KFSDe           matrix
      INTEGER(LONG)            :: NTERM_KGG           =   0      ! Count of no. of terms in KGG             matrix
      INTEGER(LONG)            :: NTERM_KGGD          =   0      ! Count of no. of terms in KGGD            matrix
      INTEGER(LONG)            :: NTERM_KLL           =   0      ! Count of no. of terms in KLL             matrix
      INTEGER(LONG)            :: NTERM_KLLD          =   0      ! Count of no. of terms in KLLD            matrix
      INTEGER(LONG)            :: NTERM_KLLDn         =   0      ! Count of no. of terms in KLLDn           matrix
      INTEGER(LONG)            :: NTERM_KLLs          =   0      ! Count of no. of terms in KLLs            matrix
      INTEGER(LONG)            :: NTERM_KLLDs         =   0      ! Count of no. of terms in KLLDs           matrix
      INTEGER(LONG)            :: NTERM_KMM           =   0      ! Count of no. of terms in KMM             matrix
      INTEGER(LONG)            :: NTERM_KMMD          =   0      ! Count of no. of terms in KMMD            matrix
      INTEGER(LONG)            :: NTERM_KMSM          =   0      ! Count of no. of terms in KMSM            matrix
      INTEGER(LONG)            :: NTERM_KMSMn         =   0      ! Count of no. of terms in KMSMn           matrix
      INTEGER(LONG)            :: NTERM_KMSMs         =   0      ! Count of no. of terms in KMSMs           matrix
      INTEGER(LONG)            :: NTERM_KNM           =   0      ! Count of no. of terms in KNM             matrix
      INTEGER(LONG)            :: NTERM_KNMD          =   0      ! Count of no. of terms in KNMD            matrix
      INTEGER(LONG)            :: NTERM_KNN           =   0      ! Count of no. of terms in KNN             matrix
      INTEGER(LONG)            :: NTERM_KNND          =   0      ! Count of no. of terms in KNND            matrix
      INTEGER(LONG)            :: NTERM_KOO           =   0      ! Count of no. of terms in KOO             matrix
      INTEGER(LONG)            :: NTERM_KOOD          =   0      ! Count of no. of terms in KOOD            matrix
      INTEGER(LONG)            :: NTERM_KOOs          =   0      ! Count of no. of terms in KOOs            matrix
      INTEGER(LONG)            :: NTERM_KOODs         =   0      ! Count of no. of terms in KOODs           matrix
      INTEGER(LONG)            :: NTERM_KRL           =   0      ! Count of no. of terms in KRL             matrix
      INTEGER(LONG)            :: NTERM_KRLD          =   0      ! Count of no. of terms in KRLD            matrix
      INTEGER(LONG)            :: NTERM_KRR           =   0      ! Count of no. of terms in KRR             matrix
      INTEGER(LONG)            :: NTERM_KRRD          =   0      ! Count of no. of terms in KRRD            matrix
      INTEGER(LONG)            :: NTERM_KRRcb         =   0      ! Count of no. of terms in KRRcb           matrix
      INTEGER(LONG)            :: NTERM_KRRcbs        =   0      ! Count of no. of terms in KRRcbs          matrix
      INTEGER(LONG)            :: NTERM_KRRcbn        =   0      ! Count of no. of terms in KRRcbn          matrix
      INTEGER(LONG)            :: NTERM_KXX           =   0      ! Count of no. of terms in KRRcb           matrix
      INTEGER(LONG)            :: NTERM_KSS           =   0      ! Count of no. of terms in KSS             matrix
      INTEGER(LONG)            :: NTERM_KSSD          =   0      ! Count of no. of terms in KSSD            matrix
      INTEGER(LONG)            :: NTERM_KSSe          =   0      ! Count of no. of terms in KSSe            matrix
      INTEGER(LONG)            :: NTERM_KSSDe         =   0      ! Count of no. of terms in KSSDe           matrix
      INTEGER(LONG)            :: NTERM_LMN           =   0      ! Count of no. of terms in LMN             matrix
      INTEGER(LONG)            :: NTERM_LTM           =   0      ! Count of no. of terms in LTM             matrix
      INTEGER(LONG)            :: NTERM_MAA           =   0      ! Count of no. of terms in MAA             matrix
      INTEGER(LONG)            :: NTERM_MAO           =   0      ! Count of no. of terms in MAO             matrix
      INTEGER(LONG)            :: NTERM_MFF           =   0      ! Count of no. of terms in MPF0            matrix
      INTEGER(LONG)            :: NTERM_MFS           =   0      ! Count of no. of terms in MFS             matrix
      INTEGER(LONG)            :: NTERM_MGG           =   0      ! Count of no. of terms in MGG             matrix
      INTEGER(LONG)            :: NTERM_MGGC          =   0      ! Count of no. of terms in MGG             matrix for CONM2's
      INTEGER(LONG)            :: NTERM_MGGE          =   0      ! Count of no. of terms in MGG             matrix for elems
      INTEGER(LONG)            :: NTERM_MGGS          =   0      ! Count of no. of terms in MGG             matrix for scalar mass
      INTEGER(LONG)            :: NTERM_MLL           =   0      ! Count of no. of terms in MLL             matrix
      INTEGER(LONG)            :: NTERM_MLLn          =   0      ! Count of no. of terms in MLLn            matrix
      INTEGER(LONG)            :: NTERM_MLLs          =   0      ! Count of no. of terms in MLLs            matrix
      INTEGER(LONG)            :: NTERM_MLR           =   0      ! Count of no. of terms in MLR             matrix
      INTEGER(LONG)            :: NTERM_MMM           =   0      ! Count of no. of terms in MMM             matrix
      INTEGER(LONG)            :: NTERM_MNM           =   0      ! Count of no. of terms in MNM             matrix
      INTEGER(LONG)            :: NTERM_MNN           =   0      ! Count of no. of terms in MNN             matrix
      INTEGER(LONG)            :: NTERM_MOO           =   0      ! Count of no. of terms in MOO             matrix
      INTEGER(LONG)            :: NTERM_MPF0          =   0      ! Count of no. of terms in MFF             matrix
      INTEGER(LONG)            :: NTERM_MRL           =   0      ! Count of no. of terms in MRL             matrix
      INTEGER(LONG)            :: NTERM_MRN           =   0      ! Count of no. of terms in MRN             matrix
      INTEGER(LONG)            :: NTERM_MRR           =   0      ! Count of no. of terms in MRR             matrix
      INTEGER(LONG)            :: NTERM_MRRcb         =   0      ! Count of no. of terms in MRRcb           matrix
      INTEGER(LONG)            :: NTERM_MRRcbn        =   0      ! Count of no. of terms in MRRcbn          matrix
      INTEGER(LONG)            :: NTERM_MXX           =   0      ! Count of no. of terms in MRRcb           matrix
      INTEGER(LONG)            :: NTERM_MXXn          =   0      ! Count of no. of terms in MRRcb           matrix
      INTEGER(LONG)            :: NTERM_MSS           =   0      ! Count of no. of terms in MSS             matrix
      INTEGER(LONG)            :: NTERM_PA            =   0      ! Count of no. of terms in load            matrix PA
      INTEGER(LONG)            :: NTERM_PF            =   0      ! Count of no. of terms in load            matrix PF
      INTEGER(LONG)            :: NTERM_PFYS          =   0      ! Count of no. of terms in load            matrix PFYS
      INTEGER(LONG)            :: NTERM_PG            =   0      ! Count of no. of terms in load            matrix PG
      INTEGER(LONG)            :: NTERM_PHIXA         =   0      ! Count of no. of terms in CB              matrix PHIXA
      INTEGER(LONG)            :: NTERM_PHIXG         =   0      ! Count of no. of terms in CB              matrix PHIXG
      INTEGER(LONG)            :: NTERM_PHIZG         =   0      ! Count of no. of terms in CB              matrix PHIZG
      INTEGER(LONG)            :: NTERM_PHIZL         =   0      ! Count of no. of terms in PHIZL           matrix
      INTEGER(LONG)            :: NTERM_PHIZL1        =   0      ! Count of no. of terms in PHIZL1          matrix
      INTEGER(LONG)            :: NTERM_PHIZL2        =   0      ! Count of no. of terms in PHIZL2          matrix
      INTEGER(LONG)            :: NTERM_PL            =   0      ! Count of no. of terms in PL              matrix
      INTEGER(LONG)            :: NTERM_PM            =   0      ! Count of no. of terms in load            matrix PM
      INTEGER(LONG)            :: NTERM_PN            =   0      ! Count of no. of terms in load            matrix PN
      INTEGER(LONG)            :: NTERM_PO            =   0      ! Count of no. of terms in load            matrix PO
      INTEGER(LONG)            :: NTERM_PR            =   0      ! Count of no. of terms in PR              matrix
      INTEGER(LONG)            :: NTERM_PS            =   0      ! Count of no. of terms in load            matrix PS
      INTEGER(LONG)            :: NTERM_QM            =   0      ! Count of no. of terms in load            matrix QM
      INTEGER(LONG)            :: NTERM_QS            =   0      ! Count of no. of terms in load            matrix QS
      INTEGER(LONG)            :: NTERM_QSYS          =   0      ! Count of no. of terms in load            matrix QSYS
      INTEGER(LONG)            :: NTERM_RMG           =   0      ! Count of no. of terms in RMG             matrix
      INTEGER(LONG)            :: NTERM_RMM           =   0      ! Count of no. of terms in RMM             matrix
      INTEGER(LONG)            :: NTERM_RMN           =   0      ! Count of no. of terms in RMN             matrix
      INTEGER(LONG)            :: NTERM_ULL           =   0      ! Count of no. of terms in ULL             matrix
      INTEGER(LONG)            :: NTERM_ULLI          =   0      ! Count of no. of terms in ULLI            matrix
      INTEGER(LONG)            :: NTSUB               =   0      ! Count of no. of subcases with a thermal load
      INTEGER(LONG)            :: NUM_CB_DOFS         =   0      ! Num of Craig-Bampton DOF's (2*NDOFR+NVEC)
      INTEGER(LONG)            :: NUM_EIGENS          =   0      ! Num of eigenvals calc'd (may be > NVEC)
      INTEGER(LONG)            :: NUM_KLLD_DIAG_ZEROS =   0      ! Num of zeros on the diagonal of the KLLD stiff matrix
      INTEGER(LONG)            :: NUM_MLL_DIAG_ZEROS  =   0      ! Num of zeros on the diagonal of the MLL mass matrix
      INTEGER(LONG)            :: NUM_MPCSIDS         =   0      ! The number of MPC set ID's called for in an execution
      INTEGER(LONG)            :: NUM_PARTVEC_RECORDS =   0      ! Num of PARTVEC, PARTVEC1 records written to file LINK1V
      INTEGER(LONG)            :: NUM_PCHD_SPC1       =   0      ! Num of SPC1 records written to file SPC if req by param PCHSPC1
      INTEGER(LONG)            :: NUM_SPC_RECORDS     =   0      ! Num of SPC  cards written to file LINK1O.
      INTEGER(LONG)            :: NUM_SPC1_RECORDS    =   0      ! Num of SPC1 cards written to file LINK1O.
      INTEGER(LONG)            :: NUM_SUPT_CARDS      =   0      ! Num of SUPORT cards written to file LINK1T.
      INTEGER(LONG)            :: NUM_USET_RECORDS    =   0      ! Num of USET, USET1 records written to file LINK1X.
      INTEGER(LONG)            :: NUM_USETSTR         =   0      ! Num of PARAM, USETSTR entries in the Bulk Data
      INTEGER(LONG)            :: NUM_USET            =   0      ! Num of GRID/COMP pairs in array all USET's
      INTEGER(LONG)            :: NUM_USET_U1         =   0      ! Num of GRID/COMP pairs in array USET_U1 (1 entry for each COMP)
      INTEGER(LONG)            :: NUM_USET_U2         =   0      ! Num of GRID/COMP pairs in array USET_U2 (1 entry for each COMP)
      INTEGER(LONG)            :: NVEC                =   0      ! Num of eigenvecs requested for output (may be < NUM_EIGENS)
      INTEGER(LONG)            :: NVVEC               =   0      ! Count of the no. of v vectors (CBAR)
      INTEGER(LONG)            :: PCH_LINE_NUM        =   0      ! Line number in PCH ("punch") file
      INTEGER(LONG)            :: SETLEN              =   0      ! Count of the number characters in SET definitions
      INTEGER(LONG)            :: TRIA_ORDER_NUMS(NUM_TRIA_ORDERS) =    (/1,3,7/)
                                                                 ! Triangular isoparametric integration orders

      INTEGER(LONG), PARAMETER :: DEDAT_Q4_MATANG_KEY =   6      ! Delta in EDAT for QUAD4 to get from EID to the matl angle key
      INTEGER(LONG), PARAMETER :: DEDAT_Q4_POFFS_KEY  =   8      ! Delta in EDAT for QUAD4 to get from EID to the offset key
      INTEGER(LONG), PARAMETER :: DEDAT_Q4_SHELL_KEY  =   9      ! Delta in EDAT for QUAD4 to get from EID to the shell/pcomp key
      INTEGER(LONG), PARAMETER :: DEDAT_Q4_THICK_KEY  =  10      ! Delta in EDAT for QUAD4 to get from EID to the thickness key
      INTEGER(LONG), PARAMETER :: DEDAT_T3_MATANG_KEY =   5      ! Delta in EDAT for TRIA3 to get from EID to the matl angle key
      INTEGER(LONG), PARAMETER :: DEDAT_T3_POFFS_KEY  =   7      ! Delta in EDAT for TRIA3 to get from EID to the offset key
      INTEGER(LONG), PARAMETER :: DEDAT_T3_SHELL_KEY  =   8      ! Delta in EDAT for TRIA3 to get from EID to the shell/pcomp key
      INTEGER(LONG), PARAMETER :: DEDAT_T3_THICK_KEY  =   9      ! Delta in EDAT for TRIA3 to get from EID to the thickness key
      INTEGER(LONG), PARAMETER :: DEDAT_Q8_POFFS_KEY  =  12      ! Delta in EDAT for QUAD8 to get from EID to the offset key
      INTEGER(LONG), PARAMETER :: DEDAT_Q8_SHELL_KEY  =  13      ! Delta in EDAT for QUAD8 to get from EID to the shell/pcomp key
      INTEGER(LONG), PARAMETER :: DEDAT_Q8_THICK_KEY  =  14      ! Delta in EDAT for QUAD8 to get from EID to the thickness key
      INTEGER(LONG), PARAMETER :: MAX_FEMAP_COLS      =  24      ! Max number of columns for array FEMAP_ELEM_VECS
      INTEGER(LONG), PARAMETER :: MAX_NUM_STR         =   9      ! Number of different stresses/strains
      INTEGER(LONG), PARAMETER :: MAX_ORDER_GAUSS     =  10      ! Max order that can be used when subr ORDER_GAUSS is called
      INTEGER(LONG), PARAMETER :: MAX_ORDER_TETRA     =   4      ! Max order that can be used when subr ORDER_TETRA is called
      INTEGER(LONG), PARAMETER :: MAX_ORDER_TRIA      =   7      ! Max order that can be used when subr ORDER_TRIA  is called
      INTEGER(LONG), PARAMETER :: MAX_TOKEN_LEN       =   8      ! Max length (chars) of tokens in SET's
      INTEGER(LONG), PARAMETER :: MBUG                =  10      ! No. of kinds BUG outputs (dimension of WRT_BUG)
      INTEGER(LONG), PARAMETER :: MCMASS              =   7      ! No. cols allowed in dimensioning array CMASS
      INTEGER(LONG), PARAMETER :: MCONM2              =   3      ! No. cols allowed in dimensioning array CONM2
      INTEGER(LONG), PARAMETER :: MCORD               =   5      ! No. cols allowed in dimensioning array CORD
      INTEGER(LONG), PARAMETER :: MEDAT_CBAR          =   8      ! No. terms that go into EDAT array for CBAR   elems
      INTEGER(LONG), PARAMETER :: MEDAT_CBEAM         =   8      ! No. terms that go into EDAT array for CBEAM  elems
      INTEGER(LONG), PARAMETER :: MEDAT_CBUSH         =   9      ! No. terms that go into EDAT array for CBEAM  elems
      INTEGER(LONG), PARAMETER :: MEDAT_CELAS1        =   6      ! No. terms that go into EDAT array for CELAS1 elems
      INTEGER(LONG), PARAMETER :: MEDAT_CELAS2        =   6      ! No. terms that go into EDAT array for CELAS3 elems
      INTEGER(LONG), PARAMETER :: MEDAT_CELAS3        =   4      ! No. terms that go into EDAT array for CELAS3 elems
      INTEGER(LONG), PARAMETER :: MEDAT_CELAS4        =   4      ! No. terms that go into EDAT array for CELAS3 elems
      INTEGER(LONG), PARAMETER :: MEDAT_CHEXA8        =  10      ! No. terms that go into EDAT array for CHEXA  elems with  8 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CHEXA20       =  22      ! No. terms that go into EDAT array for CHEXA  elems with 20 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CPENTA6       =   8      ! No. terms that go into EDAT array for CPENTA elems with  6 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CPENTA15      =  17      ! No. terms that go into EDAT array for CPENTA elems with 15 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CPYRA5        =   7      ! No. terms that go into EDAT array for CPYRA  elems with  5 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CPYRA14       =  16      ! No. terms that go into EDAT array for CPYRA  elems with 14 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CQUAD         =  11      ! No. terms that go into EDAT array for CQUAD  elems
      INTEGER(LONG), PARAMETER :: MEDAT_CQUAD8        =  15      ! No. terms that go into EDAT array for CQUAD8 elems
      INTEGER(LONG), PARAMETER :: MEDAT_CROD          =   4      ! No. terms that go into EDAT array for CROD   elems
      INTEGER(LONG), PARAMETER :: MEDAT_CSHEAR        =   6      ! No. terms that go into EDAT array for CSHEAR elems with  4 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CTETRA4       =   6      ! No. terms that go into EDAT array for CTETRA elems with  4 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CTETRA10      =  12      ! No. terms that go into EDAT array for CTETRA elems with 10 nodes
      INTEGER(LONG), PARAMETER :: MEDAT_CTRIA         =  10      ! No. terms that go into EDAT array for CTRIA  elems
      INTEGER(LONG), PARAMETER :: MEDAT_CUSER1        =  11      ! No. terms that go into EDAT array for CUSER1 elems
      INTEGER(LONG), PARAMETER :: MEDAT0_CUSERIN      =   5      ! No. terms that go into EDAT array for USERIN not incl grids/comps
      INTEGER(LONG), PARAMETER :: MEDAT_PLOTEL        =   4      ! No. terms that go into EDAT array for PLOTEL elems
      INTEGER(LONG), PARAMETER :: MEFE                = 1000     ! Max number of messages to write to EMG error arrays
      INTEGER(LONG), PARAMETER :: MEFEI               =   3      ! Max number of cols in int  array for EMG error numbers
      INTEGER(LONG), PARAMETER :: MEFER               =   2      ! Max number of cols in real array for EMG error numbers
      INTEGER(LONG), PARAMETER :: MEWE                =   2      ! Max number of messages to write to EMG error arrays
      INTEGER(LONG), PARAMETER :: MEWEI               =   1      ! Max number of cols in int  array for EMG error numbers
      INTEGER(LONG), PARAMETER :: MEWER               =   2      ! Max number of cols in real array for EMG error numbers
      INTEGER(LONG), PARAMETER :: MELDTS              =  16      ! No. of types of ELDATA related outputs
      INTEGER(LONG), PARAMETER :: MELOUTS             =   4      ! No. of types of elem related outputs (ELFN, ELFE, STRE)
      INTEGER(LONG), PARAMETER :: MEOFIL              =   4      ! Max no. elem disk debug output files
      INTEGER(LONG), PARAMETER :: MEPROP              =  50      ! Max no. element properties that cab be stored in array EPROP
      INTEGER(LONG), PARAMETER :: MEPSIL              =   6      ! Max no. variables in EPSIL array
      INTEGER(LONG), PARAMETER :: METYPE              =  23      ! Max number of element types ! 
      INTEGER(LONG), PARAMETER :: MFIJ                =   5      ! Max number of disk files for WRT_FIJ (F21, F22, etc files)
      INTEGER(LONG), PARAMETER :: MGRID               =   6      ! No. cols allowed in dimensioning array GRID
      INTEGER(LONG), PARAMETER :: MGROUTS             =   7      ! No. of types of grid related outputs
!                                                                  (ACCE, DISP, OLOA, SPCF, GPFO, MPCF)
      INTEGER(LONG), PARAMETER :: MMATL               =   2      ! No. cols allowed in dimensioning array MATL
      INTEGER(LONG), PARAMETER :: MMSPRNT             =   3      ! No. cols allowed in dimensioning array MSPRNT
      INTEGER(LONG), PARAMETER :: MOGEL               =  12      ! No. cols allowed in dimensioning array OGEL
! --- cbeam_add begin --- !
      INTEGER(LONG), PARAMETER :: MPDAT_PLOAD1        =  26      ! Beam/bar PLOAD1 data slots (6 comps x [P1,P2,X1,X2]) + 2 legacy dir flags
      INTEGER(LONG), PARAMETER :: MPDAT_PLOAD2        =   1      ! No. pressures on PLOAD2 Bulk Data card
      INTEGER(LONG), PARAMETER :: MPDAT_PLOAD4        =   4      ! No. pressuresa on PLOAD4 Bulk Data card
      INTEGER(LONG), PARAMETER :: MPBAR               =   3      ! No. cols allowed in dimensioning array PBAR
      INTEGER(LONG), PARAMETER :: MPBARLU             =   6      ! Max num of dec places in format for writing PBAR equivs of PBARL
      INTEGER(LONG), PARAMETER :: MPBEAM              =   4      ! No. cols allowed in dimensioning array PBEAM
      INTEGER(LONG), PARAMETER :: MPBEAM_STATIONS     =  11      ! Max stored x/L stations for NX-style PBEAM/CBEAM output
      INTEGER(LONG), PARAMETER :: MPBUSH              =   2      ! No. cols allowed in dimensioning array PBUSH
      INTEGER(LONG), PARAMETER :: MPLOAD4_3D_DATA     =   5      ! No. cols allowed for array PLOAD4_3D_DATA
      INTEGER(LONG), PARAMETER :: MPCOMP0             =   6      ! No. integer data on PCOMP parent entry (PID,FT,LAM) + NUM_LAYERS
      INTEGER(LONG), PARAMETER :: MPCOMP_PLIES        =   2      ! No. integer data for each layer on PCOMP (MIDi, SOUTi)
      INTEGER(LONG), PARAMETER :: MPELAS              =   1      ! No. cols allowed in dimensioning array PELAS
      INTEGER(LONG), PARAMETER :: MPMASS              =   1      ! No. cols allowed in dimensioning array PMASS
      INTEGER(LONG), PARAMETER :: MPRESS              =  26      ! No. rows allowed in dimensioning array PRESS
! --- cbeam_add end --- !
      INTEGER(LONG), PARAMETER :: MPROD               =   2      ! No. cols allowed in dimensioning array PROD
      INTEGER(LONG), PARAMETER :: MPSHEAR             =   2      ! No. cols allowed in dimensioning array PSHEAR
      INTEGER(LONG), PARAMETER :: MPSHEL              =   6      ! No. cols allowed in dimensioning array PSHEL
      INTEGER(LONG), PARAMETER :: MPSOLID             =   6      ! No. cols allowed in dimensioning array PSOLID
      INTEGER(LONG), PARAMETER :: MPUSER1             =   2      ! No. cols allowed in dimensioning array PUSER1
      INTEGER(LONG), PARAMETER :: MPUSERIN            =   2      ! No. cols allowed in dimensioning array PUSERIN
      INTEGER(LONG), PARAMETER :: MRCONM2             =  10      ! No. cols allowed in dimensioning array RCONM2
      INTEGER(LONG), PARAMETER :: MRCORD              =  12      ! No. cols allowed in dimensioning array RCORD
      INTEGER(LONG), PARAMETER :: MRGRID              =   3      ! No. cols allowed in dimensioning array RGRID
      INTEGER(LONG), PARAMETER :: MRMATLC             =  30      ! No. cols allowed in dimensioning array RMATL
      INTEGER(LONG), PARAMETER :: MRSNORM             =   3      ! No. cols allowed in dimensioning array RSNORM
      INTEGER(LONG), PARAMETER :: MGRID_SNORM         =   3      ! No. cols allowed in dimensioning array GRID_SNORM
      INTEGER(LONG), PARAMETER :: MSNORM              =   2      ! No. cols allowed in dimensioning array SNORM
      INTEGER(LONG), PARAMETER :: MEMATR              = MRMATLC+2! No. rows allowed in dimensioning array EMAT
      INTEGER(LONG), PARAMETER :: MEMATC              =   4      ! No. cols allowed in dimensioning array EMAT
      INTEGER(LONG), PARAMETER :: MRPBAR              =  17      ! No. cols allowed in dimensioning array RPBAR
      INTEGER(LONG), PARAMETER :: MRPBEAM             =  45      ! No. cols allowed in dimensioning array RPBEAM
      INTEGER(LONG), PARAMETER :: MRPBUSH             =  23      ! No. cols allowed in dimensioning array RPBUSH
      INTEGER(LONG), PARAMETER :: MRPCOMP0            =   6      ! No. real data on PCOMP parent entry (Z0, NSM, SB, TREF, GE)
      INTEGER(LONG), PARAMETER :: MRPCOMP_PLIES       =   3      ! No. real data for each layer on PCOMP (Ti, THETAi) plus
!                                                                  calculated Di = +/- dist from ref plane to center of ply
      INTEGER(LONG), PARAMETER :: MRPELAS             =   3      ! No. cols allowed in dimensioning array RPELAS
      INTEGER(LONG), PARAMETER :: MRPMASS             =   1      ! No. cols allowed in dimensioning array RPMASS
      INTEGER(LONG), PARAMETER :: MRPROD              =   4      ! No. cols allowed in dimensioning array RPROD
      INTEGER(LONG), PARAMETER :: MRPSHEAR            =   4      ! No. cols allowed in dimensioning array RPSHEAR
      INTEGER(LONG), PARAMETER :: MRPSHEL             =   7      ! No. cols allowed in dimensioning array RPSHEL
      INTEGER(LONG), PARAMETER :: MRPUSER1            =  10      ! No. cols allowed in dimensioning array RPUSER1
      INTEGER(LONG), PARAMETER :: MTDOF               =  20      ! No. cols allowed in dimensioning array TDOF
      INTEGER(LONG), PARAMETER :: MTDAT_TEMPP1        =   2      ! No. fields of real data on TEMPP1 Bulk Data card
      INTEGER(LONG), PARAMETER :: MTDAT_TEMPRB        =   6      ! No. fields of real data on TEMPRB Bulk Data card
      INTEGER(LONG), PARAMETER :: MTSET               =   6      ! No. cols allowed in dimensioning array TSET
      INTEGER(LONG), PARAMETER :: MUSERIN_MAT_NAMES   =   4      ! No. cols allowed in dimensioning array USERIN_MAT_NAMES

! Parameters used to denote bit positions in subcase arrays GROUT, ELOUT and ELDT. Note: if any of these change here, the changes
! must be reflected in subr SUBCASE_PROC which assumes the following assignments since it uses IBIT(J), J = 0, 1, 2... and not the
! named variables below. GROUT, ELOUT and ELDT all have 16 bits with the following ones used:

      INTEGER(LONG), PARAMETER :: GROUT_ACCE_BIT      =   5      ! Bit pos in OGROUT, GROUT for G.P. accel print requests
      INTEGER(LONG), PARAMETER :: GROUT_DISP_BIT      =   0      ! Bit pos in OGROUT, GROUT for displ print requests
      INTEGER(LONG), PARAMETER :: GROUT_GPFO_BIT      =   4      ! Bit pos in OGROUT, GROUT for G.P. force balance print requests
      INTEGER(LONG), PARAMETER :: GROUT_MPCF_BIT      =   3      ! Bit pos in OGROUT, GROUT for MPC force print requests
      INTEGER(LONG), PARAMETER :: GROUT_OLOA_BIT      =   1      ! Bit pos in OGROUT, GROUT for OLOAD print requests
      INTEGER(LONG), PARAMETER :: GROUT_SPCF_BIT      =   2      ! Bit pos in OGROUT, GROUT for SPC force print requests

      INTEGER(LONG), PARAMETER :: ELOUT_ELFE_BIT      =   1      ! Bit pos in OELOUT, ELOUT: print requests for elem engr forces
      INTEGER(LONG), PARAMETER :: ELOUT_ELFN_BIT      =   0      ! Bit pos in OELOUT, ELOUT: print requests for elem node forces
      INTEGER(LONG), PARAMETER :: ELOUT_STRE_BIT      =   2      ! Bit pos in OELOUT, ELOUT: print requests for elem stress
      INTEGER(LONG), PARAMETER :: ELOUT_STRN_BIT      =   3      ! Bit pos in OELOUT, ELOUT: print requests for elem strain

      INTEGER(LONG), PARAMETER :: ELDT_BUG_DAT1_BIT   =   0      ! Bit pos in OELDT, ELDT for elem data printed requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_DAT2_BIT   =   1      ! Bit pos in OELDT, ELDT for elem data printed requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_P_T_BIT    =   2      ! Bit pos in OELDT, ELDT for elem pres/therm load print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_ME_BIT     =   3      ! Bit pos in OELDT, ELDT for elem mass print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_KE_BIT     =   4      ! Bit pos in OELDT, ELDT for elem stiff matrix print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_SE_BIT     =   5      ! Bit pos in OELDT, ELDT for elem stress rec mat print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_U_P_BIT    =   6      ! Bit pos in OELDT, ELDT for elem displ/load print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_SHPJ_BIT   =   7      ! Bit pos in OELDT, ELDT for elem displ/load print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_BMAT_BIT   =   8      ! Bit pos in OELDT, ELDT for elem displ/load print requests
      INTEGER(LONG), PARAMETER :: ELDT_BUG_BCHK_BIT   =   9      ! Bit pos in OELDT, ELDT for elem displ/load print requests

      INTEGER(LONG), PARAMETER :: ELDT_F21_P_T_BIT    =  10      ! Bit pos in OELDT, ELDT for elem mass file requests
      INTEGER(LONG), PARAMETER :: ELDT_F22_ME_BIT     =  11      ! Bit pos in OELDT, ELDT for elem press/therm load file requests
      INTEGER(LONG), PARAMETER :: ELDT_F23_KE_BIT     =  12      ! Bit pos in OELDT, ELDT for elem stiff matrix file requests
      INTEGER(LONG), PARAMETER :: ELDT_F24_SE_BIT     =  13      ! Bit pos in OELDT, ELDT for elem stress rec mat file requests
      INTEGER(LONG), PARAMETER :: ELDT_F25_U_P_BIT    =  14      ! Bit pos in OELDT, ELDT for elem displ/load file requests


      REAL(DOUBLE) , PARAMETER :: FEMAP_VERSION       = 8.2      ! FEMAP   Version number
      REAL(SINGLE)             :: KMAT_DEN            = ZERO     ! Stiffness matrix density from BANDIT
      REAL(DOUBLE)             :: TOT_MB_MEM_ALLOC    = ZERO     ! The total amount of memory, in MB, allocated at any time

      INTEGER(LONG)            :: COUNTER_VALUE       = ZERO     ! Current value of the counter
      INTEGER(LONG)            :: COUNTER_PERC        = ZERO     ! Current percentage of the counter.
      INTEGER(LONG)            :: COUNTER_TOTAL       = ZERO     ! Max value of the counter
      INTEGER(LONG)            :: COUNTER_STARTED     = ZERO     ! Timestamp of counter start
      INTEGER(LONG)            :: COUNTER_UPDATED     = ZERO     ! Timestamp of counter update
      INTEGER(LONG)            :: COUNTER_LIMITER     = ZERO     ! Counter updates this second
      CHARACTER(:), ALLOCATABLE :: COUNTER_PREFIX
      CHARACTER(LEN=20) :: COUNTER_FMT
      INTEGER(LONG)            :: NUM_SUBC_CARDS      = ZERO     ! Number of subcase CARDS.

!   Notes:
!    (1) The max limits are only used to check that the dimensions of arrays have been set high enough
!        They are initialized here to 0 and are counted when the input deck is read in subrs LOADC0, LOADB0

      END MODULE SCONTR
