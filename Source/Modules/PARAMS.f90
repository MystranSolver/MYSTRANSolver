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

      MODULE PARAMS

! Variables that can be set by Bulk Data PARAM cards. Most are actual PARAM names; however, a few are variables set on a PARAM card.

      USE PENTIUM_II_KIND, ONLY :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY     :  ONEPM5, ONEPM6, ONEPM8, ONEPM15, ONEPM14, ZERO, QUARTER, TENTH, ONE, TWO, THREEP6, FIVE, SIX,   &
                                   ONEPP6, ONEPP7, ONE_THOUSAND 
      USE SCONTR, ONLY          :  MEPSIL, TSET_CHR_LEN

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), PRIVATE   :: I                          ! Only use is in implied DO loops below to initialize variables

! **********************************************************************************************************************************

! Parameters used ------------------------------------------------------------------------------------------------------------------
! >>>>>>>>>>>>>>>

      REAL(DOUBLE)             :: ARP_TOL        =   1.0E-9    ! Input to ARPACK subr dsband to decide convergence in subr dsconv.
!                                                              NOTE: if a value of -1. is input on a PARAM ARP_TOL entry, then the
!                                                                    Lanczos algorithm will use machine precision for ARP_TOL !ONEPM6 

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: ART_KED        =    'N'    ! Indicates whether to add artificial differ stiff terms to KED
      REAL(DOUBLE)             :: ART_TRAN_KED   =   ONEPM6  ! Artificial differ stiff for translational DOF's
      REAL(DOUBLE)             :: ART_ROT_KED    =   ONEPM6  ! Artificial differ stiff for rotational DOF's

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: ART_MASS       =    'N'    ! Indicates whether to add artificial mass terms to transl or rot DOF's
      REAL(DOUBLE)             :: ART_TRAN_MASS  =   ONEPM6  ! Artificial mass for translational DOF's
      REAL(DOUBLE)             :: ART_ROT_MASS   =   ONEPM6  ! Artificial mass for rotational DOF's

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: AUTOSPC        =    'Y'    ! 'Y'/'N' indicates whether to use automatic SPC for singular DOF's
      REAL(DOUBLE)             :: AUTOSPC_RAT    =   ONEPM8  
                                                             !*For each of the 2 3x3 (K33) stiffness matrices for a grid, the ratio
!                                                              of each of the 3 eigenvalues of K33 to the max eigenvalue is calc'd.
!                                                              For any ratio <= AUTOSPC_RAT, a component will be AUTOSPC'd. The
!                                                              grid comp that is AUTOSPC'd is the one whose eigenvector value for
!                                                              that eigenvalue is maximum.  
      INTEGER(LONG)            :: AUTOSPC_NSET   =     1     !*If = 1, SPC KNN for null rows. If 2, SPC KNN for small diag terms
!                                                              If = 3, do both                               
      CHARACTER(  1*BYTE)      :: AUTOSPC_INFO   =    'N'    ! If 'Y' print information on AUTOSPC's
      CHARACTER(  1*BYTE)      :: AUTOSPC_SPCF   =    'N'    ! If 'Y' write SPC forces on SA DOF's in subr OFP2

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: BAILOUT        =     0     ! If >= 0 quit if a singularity in decomposing a matrix is detected

! ----------------------------------------------------------------------------------------------------------------------------------
! --- banded_optimization begin --- !
      CHARACTER(  1*BYTE)      :: BANDEDOPT      =    'N'    ! 'Y'/'N' enable experimental banded-order optimization path RCM
! --- banded_optimization end --- !
! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: CBMIN3         =  TWO      ! Trans shear factor for MIN3  triangle elems (TRIA3)

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: CBMIN4         =  THREEP6  ! Trans shear factor for MIN4  quad    elems (QUAD4, QUAD4TYP='MIN4  ')

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: CBMIN4T        =  THREEP6  ! Trans shear factor for MIN4T quad    elems (QUAD4, QUAD4TYP='MIN4T ')

! ----------------------------------------------------------------------------------------------------------------------------------
! --- cbeam_add begin --- !
      REAL(DOUBLE)             :: CBEAMSHR       =  ONE      ! Global multiplier on beam shear factors K1/K2 before CBAR/CBEAM
!                                                              stiffness assembly. 0.0 gives Bernoulli/no-shear behavior.
      INTEGER(LONG), PARAMETER :: MCBEAMSHR_PID  =   200     ! Max per-property beam shear overrides from PARAM,CBEAMSHR,value,PID
      INTEGER(LONG)            :: NCBEAMSHR_PID  =     0     ! Number of active per-property beam shear overrides
      INTEGER(LONG)            :: CBEAMSHR_PID(MCBEAMSHR_PID) = 0
      REAL(DOUBLE)             :: CBEAMSHR_VAL(MCBEAMSHR_PID) = ONE
      REAL(DOUBLE)             :: CBEAMAREA      =  ONE      ! Global multiplier on beam stiffness area before CBAR/CBEAM
!                                                              stiffness assembly. Mass and self-weight continue to use geometric area.
      INTEGER(LONG), PARAMETER :: MCBEAMAREA_PID =   200     ! Max per-property beam area overrides from PARAM,CBEAMAREA,value,PID
      INTEGER(LONG)            :: NCBEAMAREA_PID =     0     ! Number of active per-property beam area overrides
      INTEGER(LONG)            :: CBEAMAREA_PID(MCBEAMAREA_PID) = 0
      REAL(DOUBLE)             :: CBEAMAREA_VAL(MCBEAMAREA_PID) = ONE

      REAL(DOUBLE)             :: BEAMAMO       =  ONE      ! Global multiplier on beam axial stiffness EA before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMAMO_PID  =   200     ! Max per-property beam axial overrides from PARAM,BEAMAMO,value,PID
      INTEGER(LONG)            :: NBEAMAMO_PID  =     0     ! Number of active per-property beam axial overrides
      INTEGER(LONG)            :: BEAMAMO_PID(MBEAMAMO_PID) = 0
      REAL(DOUBLE)             :: BEAMAMO_VAL(MBEAMAMO_PID) = ONE

      REAL(DOUBLE)             :: BEAMV1MO       =  ONE      ! Global multiplier on beam plane-1 stiffness before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMV1MO_PID  =   200     ! Max per-property beam plane-1 overrides from PARAM,BEAMV1MO,value,PID
      INTEGER(LONG)            :: NBEAMV1MO_PID  =     0     ! Number of active per-property beam plane-1 overrides
      INTEGER(LONG)            :: BEAMV1MO_PID(MBEAMV1MO_PID) = 0
      REAL(DOUBLE)             :: BEAMV1MO_VAL(MBEAMV1MO_PID) = ONE

      REAL(DOUBLE)             :: BEAMV2MO       =  ONE      ! Global multiplier on beam plane-2 stiffness before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMV2MO_PID  =   200     ! Max per-property beam plane-2 overrides from PARAM,BEAMV2MO,value,PID
      INTEGER(LONG)            :: NBEAMV2MO_PID  =     0     ! Number of active per-property beam plane-2 overrides
      INTEGER(LONG)            :: BEAMV2MO_PID(MBEAMV2MO_PID) = 0
      REAL(DOUBLE)             :: BEAMV2MO_VAL(MBEAMV2MO_PID) = ONE

      REAL(DOUBLE)             :: BEAMM1MO       =  ONE      ! Global multiplier on beam major-axis bending stiffness EI1 before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMM1MO_PID  =   200     ! Max per-property beam major-bending overrides from PARAM,BEAMM1MO,value,PID
      INTEGER(LONG)            :: NBEAMM1MO_PID  =     0     ! Number of active per-property beam major-bending overrides
      INTEGER(LONG)            :: BEAMM1MO_PID(MBEAMM1MO_PID) = 0
      REAL(DOUBLE)             :: BEAMM1MO_VAL(MBEAMM1MO_PID) = ONE

      REAL(DOUBLE)             :: BEAMM2MO       =  ONE      ! Global multiplier on beam minor-axis bending stiffness EI2 before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMM2MO_PID  =   200     ! Max per-property beam minor-bending overrides from PARAM,BEAMM2MO,value,PID
      INTEGER(LONG)            :: NBEAMM2MO_PID  =     0     ! Number of active per-property beam minor-bending overrides
      INTEGER(LONG)            :: BEAMM2MO_PID(MBEAMM2MO_PID) = 0
      REAL(DOUBLE)             :: BEAMM2MO_VAL(MBEAMM2MO_PID) = ONE

      REAL(DOUBLE)             :: BEAMTMO       =  ONE      ! Global multiplier on beam torsional stiffness GJ before CBAR/CBEAM
      INTEGER(LONG), PARAMETER :: MBEAMTMO_PID  =   200     ! Max per-property beam torsion overrides from PARAM,BEAMTMO,value,PID
      INTEGER(LONG)            :: NBEAMTMO_PID  =     0     ! Number of active per-property beam torsion overrides
      INTEGER(LONG)            :: BEAMTMO_PID(MBEAMTMO_PID) = 0
      REAL(DOUBLE)             :: BEAMTMO_VAL(MBEAMTMO_PID) = ONE
! --- cbeam_add end --- !
! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: CHKGRDS        =    'Y'    ! If 'Y' call GET_ELEM_AGRID_BGRID to check all grids on elems exist

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  3*BYTE)      :: CRS_CCS        =  'CRS'    ! If CRS matrices are stored in compressed row storage
!                                                              If CCS matrices are stored in compressed col storage

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: CUSERIN        =    'N'    ! If 'Y' write CUSERIN B.D. card images
      INTEGER(LONG)            :: CUSERIN_EID    =  9999999  ! USERIN elem ID if CUSERIN B.D. card images written to F06
      INTEGER(LONG)            :: CUSERIN_PID    =  9999999  ! USERIN prop ID if CUSERIN B.D. card images written to F06
      INTEGER(LONG)            :: CUSERIN_SPNT_ID=    1001   ! USERIN SPOINT beginning ID if CUSERIN B.D. card images written to F06
      INTEGER(LONG)            :: CUSERIN_IN4    =  9999999  ! USERIN IN4 index ID if CUSERIN B.D. card images written to F06
      CHARACTER(LEN(TSET_CHR_LEN))                                                                                                 &
                               :: CUSERIN_XSET   =   '  '    ! Displ set to use when gen tables of grids/comps for a USERIN elem
      INTEGER(LONG)            :: CUSERIN_COMPTYP=     0     ! If 0, write components on CUSERIN B.D. entries in compacted form,
!                                                              otherwise write out in expanded form of 6 characters
!                                                              (i.e 1346 would be written as 1 34 6)

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: DARPACK        =     2     ! Delta to add to EIG_N2 so that ARPACK will find a few more eigens 
!                                                              than user requested (since higher ones seem to be a little bad)

      INTEGER(LONG)            :: DELBAN         =     1     ! Delete the Bandit files left over, 0 is not to. 1 is to delete them.

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: EIGESTL        =  5000     ! Upper limit on NDOFL for running code to est # eigens below EIG_FRQ2

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: EIGNORM2       =    'N'    ! If 'Y' then eigenvectors will be renormed a last time by multiplying
!                                                              by a set of scale factors (1 per eigenvector) supplied in a file with
!                                                              the same name as the input file and extension 'EIN' (if it exists)

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: ELFORCEN       = 'GLOBAL'  ! Coord sys for calc elem nodal forces ('GLOBAL', 'LOCAL', 'BASIC')

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: EPSERR         =    'Y'    ! 'Y'/'N' indicates whether to calc epsilon error estimate

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: EPSIL(MEPSIL)  =  (/ONEPM15,ONEPM14,ONEPM5,TENTH,ONEPM6,ONEPM15/)
                                                             ! Roundoff parameters:
!                                                              (1) EPSIL(1) - will change to LAPACK MACH_PREC in LINK1
!                                                              (2) EPSIL(2) - used in subr DGETF2 for testing convergence
!                                                              (3) EPSIL(3) - used in subr EIG_INV_PWR for testing convergence
!                                                              (4) EPSIL(4) - used in QUAD element warping warning message
!                                                              (5) EPSIL(5) - used in BAR and ROD element margin of safety calcs
!                                                              (6) EPSIL(6) - used in BAR margin of safety calcs

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: EQCHK_REF_GRID =     0     ! Ref grid for stiffness matrix equilibrium check
      INTEGER(LONG)            :: EQCHK_OUTPUT(5)=     (/0,0,0,0,0/)
!                                                              I = 1 is G-set, 2 is N-set, 3 is F-set and 4 is A-set and 5 is L-set
!                                                              EQCHK_OUTPUT(I) = 1 print loads due to rigid body displ RBGLOBAL
!                                                                              = 2 print strain energy due to RBGLOBAL
!                                                                              = 3 print both
      REAL(DOUBLE)             :: EQCHK_TINY     =  ONEPM5   ! When equil check is done, forces smaller than this won't be printed
      CHARACTER(  1*BYTE)      :: EQCHK_NORM     =    'N'    ! 'Y'/'N' indicates whether to norm the equil check on stiff mat diags

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: PRTALL         =    'N'    ! 'Y', 'N' flag to write all output for all files regardless of other flags
      CHARACTER(  1*BYTE)      :: PRTF06         =    'N'    ! 'Y', 'N' flag to write all f06 outputs regardless of other flags besides PRTALL
      CHARACTER(  1*BYTE)      :: PRTNEU         =    'N'    ! 'Y', 'N' flag to write all neu outputs regardless of other flags besides PRTALL
      CHARACTER(  1*BYTE)      :: PRTOP2         =    'N'    ! 'Y', 'N' flag to write all op2 outputs regardless of other flags besides PRTALL

! case 1: PRTALL=Y, PRTOP2=N -> all op2 output will be created and all neu output 
! case 2: PRTALL=N, PRTOP2=Y -> all op2 output will be created
! case 3: PRTALL=N, PRTOP2=N -> do whatever the case control says
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: GRDPNT         =    -1     ! Ref grid for GPWG (neg integer so 0 can be interp as basic origin)

      INTEGER(LONG)            :: GRDPNT_IN      =    -1     ! Value of GRDPNT read in the Bulk Data File

! ----------------------------------------------------------------------------------------------------------------------------------
! --- banded_optimization begin --- !
      CHARACTER(  8*BYTE)      :: GRIDSEQ        = 'INPUT   '! Method for sequencing grids:
!                                                              BANDIT for bandit auto grid swquencing
!                                                              GRID for grid numerical order
!                                                              INPUT for grid input order
!                                                              RCM for reserved/add-on RCM path (currently mapped to INPUT flow)
! --- banded_optimization end --- !
      CHARACTER(  1*BYTE)      :: SEQQUIT        =    'N'    !*'Y', 'N' indicator to stop processing if G.P. auto sequencing failed
!                                                               (goes in field 4 of PARAM GRIDSEQ entry)
      CHARACTER(  1*BYTE)      :: SEQPRT         =    'N'    !*'Y', 'N' indicator to print SEQGP card images from bandit
!                                                               (goes in field 5 of PARAM GRIDSEQ entry)

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: HEXAXIS        = 'SPLITD'  ! 'SIDE12', use side 1-2 as the local elem x axis.
!                                                              'SPLITD', use angle that splits the 2 diags to define the elem x axis

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: IORQ1M         =     2     ! Integration order for membrane         strains for QUAD4,4K    

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: IORQ1S         =     1     ! Integration order for in-plane shear   strains for QUAD4,4K

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: IORQ1B         =     2     ! Integration order for bending          strains for QUAD4K

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: IORQ2B         =     2     ! Integration order for bending          strains for QUAD4

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: IORQ2T         =     3     ! Integration order for transverse shear strains for QUAD4

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: ITMAX          =     5     ! Max iter's used in refining the soln when LAPACK is used in LINK3

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: KLLRAT         =    'Y'    ! 'Y', 'N' to tell whether to calc ratio of max/min KLL diagonal terms

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: KOORAT         =    'Y'    ! 'Y', 'N' to tell whether to calc ratio of max/min KOO diagonal terms

! ----------------------------------------------------------------------------------------------------------------------------------
! --- dense_subspace_feast_add  begin --- !
      CHARACTER(  6*BYTE)      :: LANCMETH       = 'ARPACK'  ! Lanczos backend:
!                                                              ARPACK (implemented)
!                                                              FEAST / CHASE (accepted; currently fallback to ARPACK in LINK4)
! --- dense_subspace_feast_add  end --- !

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: MATSPARS       =    'Y'    ! 'Y' for use of sparse SFF, SFS, SSS or 'N' for full matrix add/mult

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: MAXRATIO       =  ONEPP7   ! Max value of matrix diagonal to factor diagonal before messages are 
!                                                              written and BAILOUT tested for aborting run

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  3*BYTE)      :: MIN4TRED       =   'STC'   ! How to reduce the 5 node MIN4T element to 4 nodes
!                                                              'B54': Use Tessler's reduction with Kirchoff constraints
!                                                              'STC': Use static (Guyan) condensation

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: MEFMCORD       =     0     ! Coordinate system for output of modal effective masses

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: MEFMLOC        = 'GRDPNT'  ! Ref loc for calc modal effective mass (can be 'CG', 'GRID', 'GRDPNT'
      INTEGER(LONG)            :: MEFMGRID       =     0     ! Ref grid for calc modal effec mass (if MEFMLOC = 'GRID' pr 'GRDPNT'

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: MEMAFAC        = 0.9       ! Factor to multiply the size request of memory to be allocated when
!                                                              looping to find an allowable amount of memory to allocate. Used when
!                                                              the initial request for memory (in subrs ESP or EMP) cannot be met
!                                                              and we know that the request is conservative.

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: MPFOUT         =    '6'    ! MPFACTOR's can be output for 6 components (rel to MEFMGRID)
!                                                              or for all NDOFR (MPFOUT = 'R') DOF's

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: MXALLOCA       =    10     ! Max number of attempts to allow when trying to allocate memory

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: MXITERI        =    50     ! Max number of iterations in Inverse Power eigenvalue method

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: MXITERL        =    50     ! Max number of iterations in Lanczos eigenvalue method

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: OTMSKIP        =     0     ! Number of lines to skip between segments of OTM text file output

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PBARLDEC       =     5     ! Num of dec digits when writing PBAR equivs for PBARL entry real data

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: PBARLSHR       = 'Y'       ! Include K1, K2 for PBAR equiv to PBARL BAR props

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: PCHSPC1        =    'N'    ! 'Y', 'N' indicator to write file with SPC's from G.P. sing. proc.
      INTEGER(LONG)            :: SPC1SID        =  9999999  !*Set ID to use on SPC1's from G.P. singularity processor
!                                                              (goes in field 4 of PARAM PCHSPC1 entry)
      CHARACTER(  1*BYTE)      :: SPC1QUIT       =    'N'    ! 'Y', 'N' indicator to stop processing if G.P. singularities found

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PCOMPEQ        =     0     ! Indicator to write equiv PSHELL, MAT2 to F06 for PCOMP's. If > 0,
!                                                              write the equivalent PSHELL amd MAT2 Bulk Data entries for the PCOMP.
!                                                              If > 1 also write the data in a format with a greater number of
!                                                              digits of accuracy

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: PCMPTSTM       = FIVE/SIX  ! Ratio: (shear thick/total plate thick) for a plate using PCOMP props.
!                                                              A value of (pi^2)/12 is more accurate for Mindlin plates but 5/6
!                                                              seems to be a more popular value in other FEA programs

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: POST           =     0     ! If = -1 write FEMAP Neutral File (NEU)

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTBASIC       =     0     ! If = 1, print grids in the basic coordinate system

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTCGLTM       =     0     ! If = 1, print CB matrix CG_LTM  

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTCONN        =     0     ! If = 1, print table of elements connected to each grid

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTCORD        =     0     ! If = 1, print coord sys transform data. If 2 print more detailed info

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTDISP(5)     =     (/0,0,0,0,0/)
!                                                                   PRTDISP (1) = 1 prints UG
!                                                                   PRTDISP (2) = 1 or 3 prints sparse  UN, 2 or 3 prints  UM
!                                                                   PRTDISP (3) = 1 or 3 prints sparse  UF, 2 or 3 prints  US
!                                                                   PRTDISP (4) = 1 or 3 prints sparse  UA, 2 or 3 prints  UO
!                                                                   PRTDISP (5) = 1 or 3 prints sparse  UL, 2 or 3 prints  UR
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTDOF         =     0     ! If = 1, print TDOF. If 2, print TDOFI. If 3, print both

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTDLR         =     0     ! PRTDLR  = 1 prints DLR

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTFOR(5)      =     (/0,0,0,0,0/)
!                                                                   PRTFOR  (1) = 1 prints PG
!                                                                   PRTFOR  (2) = 1 or 3 prints sparse  PN, 2 or 3 prints  PM
!                                                                   PRTFOR  (3) = 1 or 3 prints sparse  PF, 2 or 3 prints  PS
!                                                                   PRTFOR  (4) = 1 or 3 prints sparse  PA, 2 or 3 prints  PO
!                                                                   PRTFOR  (5) = 1 or 3 prints sparse  PL, 2 or 3 prints  PR
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTGMN         =     0     ! PRTGMN  = 1 prints GMN

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTGOA         =     0     ! PRTGOA  = 1 prints GOA

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTHMN         =     0     ! If = 1, print HMN constraint matrix

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTIFLTM       =     0     ! If = 1, print CB matrix IF_LTM  

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTKXX         =     0     ! If = 1, print CB matrix KXX  

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTMASS(5)     =     (/0,0,0,0,0/)
!                                                                   PRTMASS(1) = 1 prints sparse MGG, 2 also MGGC, MGGE, MGGS
!                                                                   PRTMASS(2) = 1 or 3 prints sparse MNN, 2 or 3 prints MNM, MMM
!                                                                   PRTMASS(3) = 1 or 3 prints sparse MFF, 2 or 3 prints MFS, MSS
!                                                                   PRTMASS(4) = 1 or 3 prints sparse MAA, 2 or 3 prints MAO, MOO
!                                                                   PRTMASS(5) = 1 or 3 prints sparse MLL, 2 or 3 prints MRL, MRR

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTMASSD(5)    =     (/0,0,0,0,0/)
!                                                                   PRTMASSD(1) = 1 prints diag of MGG
!                                                                   PRTMASSD(2) = 1 prints diag of MNN
!                                                                   PRTMASSD(3) = 1 prints diag of MFF
!                                                                   PRTMASSD(4) = 1 prints diag of MAA
!                                                                   PRTMASSD(5) = 1 prints diag of MLL
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTMXX         =     0     ! If = 1, print CB matrix MXX  

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTOU4         =     0     ! If > 0, print all OU4 matrices written to OPi files for a run 
      INTEGER(LONG)            :: PRTOU4_FMT     =     0     ! Format to write OU4 matrices

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTPHIXA       =     0     ! If = 1, print CB matrix PHIXA

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTPHIZL       =     0     ! If = 1, print CB matrix PHIZL  

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTQSYS        =     0     ! Print QSYS matrix

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTRMG         =     0     ! Print RMG matrix and/or it's partitions
!                                                                    PRTRMG = 1 prints RMG
!                                                                    PRTRMG = 2 prints partitions RMM, RMN of RMG
!                                                                    PRTRMG = 3 prints both RMG and it's partitions RMM, RMN

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTSCP         =     0     ! If = 1, print data generated in the subcase processor

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTPSET        =     0     ! If > 0, print PSET tables

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTTSET        =     0     ! If > 0, print TSET table

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTUSET        =     0     ! If > 0, print USET tables

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTSTIFD(5)    =     (/0,0,0,0,0/)
!                                                                   PRTSTIFD(1) = 1 prints diag of KGG
!                                                                   PRTSTIFD(2) = 1 prints diag of KNN
!                                                                   PRTSTIFD(3) = 1 prints diag of KFF
!                                                                   PRTSTIFD(4) = 1 prints diag of KAA
!                                                                   PRTSTIFD(5) = 1 prints diag of KLL
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTSTIFF(5)    =     (/0,0,0,0,0/)
!                                                                   PRTSTIFF(1) = 1 prints sparse KGG
!                                                                   PRTSTIFF(2) = 1 or 3 prints sparse KNN, 2 or 3 prints KNM, KMM
!                                                                   PRTSTIFF(3) = 1 or 3 prints sparse KFF, 2 or 3 prints KFS, KSS
!                                                                   PRTSTIFF(4) = 1 or 3 prints sparse KAA, 2 or 3 prints KAO, KOO
!                                                                   PRTSTIFF(5) = 1 or 3 prints sparse KLL, 2 or 3 prints KRL, KRR
! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTUO0         =     0     ! Print UO0  matrix

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: PRTYS          =     0     ! Print YS   matrix of enforced displs

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: Q4SURFIT       =     3     ! Number of terms in the polynomial surface fit for extrapolating
!                                                              QUAD4 stress/strain to element corners.
!                                                                    Q4SURFIT = 1: Constant fit
!                                                                    Q4SURFIT = 3: Bi-linear function without a cross term
!                                                                    Q4SURFIT = 6: Bi-quadratic function

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: QUADAXIS       = 'SPLITD'  ! 'SIDE12', use side 1-2 as the local elem x axis.
!                                                              'SPLITD', use angle that splits the 2 diags to define the elem x axis

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: QUAD4TYP       =  'MIN4  ' ! Which element to use in MYSTRAN as the QUAD4 element
!                                                              'MIN4T ': Use Tessler's MIN4T element made up of 4 MIN3 triangles
!                                                              'MIN4  ': Use Tessler's MIN4 element 
!                                                              'MITC4 ': Use Bathe's MITC4 (1983)
!                                                              'MITC4+': Use Ko/Bathe's 2nd MITC4+ (Nov 2016)

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: RELINK3        =    'N'    ! 'Y', 'N' indicator to redo LINK3,5 on a restart

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: RCONDK         =    'N'    ! 'Y', 'N' indicator to calc KLL LAPACK cond no. (needed for error est)

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: SETLKTK        =     0     ! Integer option for how to calculate LTERM for G-set stiff matrices
!                                                              If = 0 then estimate of LTERM is based on full elem KE
!                                                                     matrices unconnected (i.e. connected DOF's recounted)
!                                                              If = 1 then estimate of LTERM is based on matrix bandwidth
!                                                                     from BANDIT times number of rows in stiff matrix.  
!                                                              If = 2 then estimate of LTERM is based on actual elem KE
!                                                                     matrices unconnected.
!                                                              If = 3 then the value in field 4 of the PARAM SETLKTK card is
!                                                                     used as the estimate for LTERM
      CHARACTER(  1*BYTE)      :: ESP0_PAUSE     =    'N'    !*Flag to indicate "PAUSE" LINK1 after subr ESP0
!                                                              (goes in field 4 of PARAM SETLKTK entry)
      INTEGER(LONG)            :: USR_LTERM_KGG  =     0     ! User supplied value for LTERM_KGG (via PARAM SETLKT B.D. card)

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: SETLKTM         =    0     ! Integer option number for how to calculate LTERM_MGGE
!                                                              If = 0 then estimate of LTERM_MGGE is based on full elem ME
!                                                                     matrices unconnected (i.e. connected DOF's recounted)
!                                                              If = 3 then estimate of LTERM_MGGE is based on actual elem ME
!                                                                     matrices unconnected.  
      CHARACTER(  1*BYTE)      :: EMP0_PAUSE     =    'N'    !*Flag to indicate "PAUSE" LINK1 after subr EMP0
!                                                             (goes in field 4 of PARAM SETLKTM entry)

      INTEGER(LONG)            :: USR_LTERM_MGG  =     0     ! User supplied value for LTERM_MGG (via PARAM SETLKT B.D. card)

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: SHRFXFAC       =  ONEPP6   ! Factor used to adjust transverse shear stiffness when user has
!                                                              indicated zero shear flexibility for shell elements. The shear
!                                                              stiffness will be reset from infinite (zero flexibility) to
!                                                              SHRFXFAC times the average of the bending stiffnesses in the 2 planes
! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: SKIPMGG        =    'N'    ! 'Y', 'N' indicator to say whether to skip calculation of MGG
!                                                               in which case MGG will be read from previously generated,
!                                                               and saved, files (LINK1R for MGG)

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  8*BYTE)      :: SOLLIB         = 'SPARSE  '! If 'BANDED  ', use LAPACK and ARPACK for eqn soln and eigens.
!                                                              If 'SPARSE  ', use value determined by parameter SPARSE_FLAVOR
!                                                              defined in field 4 of the PARAM, SOLLIB entry
                                                             
      CHARACTER(  8*BYTE)      :: SPARSE_FLAVOR  = 'SUPERLU '! This denotes which SPARSE SOLLIB to use. Currently SuperLU is the
!                                                              only option

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: SORT_MAX       =     5     ! Max number of times to run sort algorithm before stopping with error.
!                                                              SORT_MAX can be changed with a B.D. PARAM SORT_MAX card

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  6*BYTE)      :: SPARSTOR       = 'NONSYM'  ! If 'SYM', compressed storage for matrices is done in symmetric form
!                                                              where only nonzero terms on and above the diagonal are stored.
!                                                              If 'NONSYM', all nonzero terms are stored

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: STR_CID        =    -1     ! Indicator for the coordinate system to use ID for elem stress, strain
!                                                              and emgineering force output:
!                                                              -2 is material coordinate system for solids and local element
!                                                                 coord system for other element types
!                                                              -1 is local element coord system (default)
!                                                               0 is basic coord system
!                                                               j (any other integer) is a defined coord system for output

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: SUPINFO        =    'Y'    ! 'Y', 'N' indicator to supress info msg's in the F06 file

! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: SUPWARN        =    'Y'    ! 'Y', 'N' indicator to supress warn msg's in the F06 file

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: THRESHK        =  TENTH    !*Actual value used for threshold in deciding whether to equilibrate 
!                                                              stiffness matrix in LAPACK subr DLAQSB in module LAPACK_BLAS_AUX_1

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: TINY           =  ZERO     ! Filter for small terms in matrix print

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: TSTM_DEF       =  FIVE/SIX ! Default value for TS/TM on PSHEL cards. Can be reset with PARAM card 

! ----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG)            :: USR_JCT        =     0     ! User supplied (PARAM B.D. card) value for JCT - used in sort subr's

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: WINAMEM        =2147.483647! Max MB of memory that Windows XP allows for arrays

! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: WTMASS         =   ONE     ! Value from PARAM WTMASS Bulk Data card
! ----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE)             :: K6ROT          =   0.0     ! Virtual rotational stiffness
! ----------------------------------------------------------------------------------------------------------------------------------
      CHARACTER(  1*BYTE)      :: NOCOUNTS        =    'N'   ! 'Y' supresses the "counter" non-advance writes to stdout.
                                                             ! Good for debugging or when writing to a log file.

! Fixed values ---------------------------------------------------------------------------------------------------------------------
! >>>>>>>>>>>>

      REAL(DOUBLE) , PARAMETER :: THRESHK_LAP    =  TENTH    ! LAPACK suggested value for threshold used in equilibrate decision for
!                                                              stiffness matrix in LAPACK subr DLAQSB in module LAPACK_BLAS_AUX_1

! Obsolete, or not used ------------------------------------------------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>>>

      INTEGER(LONG)            :: F06_COL_START  =     0     ! 1st col in F06 file for output data to begin. If it is not > 2, then
!                                                              output will be written with each main header centered on one another

      INTEGER(LONG)            :: SPIENV6        =   50     ! Memory growth factor for SuperLU_MT -- corresponds to -sp_ienv(6)
      INTEGER(LONG)            :: SPIENV7        =   50     ! Memory growth factor for SuperLU_MT -- corresponds to -sp_ienv(7)
      INTEGER(LONG)            :: SPIENV8        =   30     ! Memory growth factor for SuperLU_MT -- corresponds to -sp_ienv(8)

      INTEGER(LONG)            :: SLU_NTHR        =    0     ! Number of threads for SuperLU_MT.
                                                             ! Using 0 will use as many threads as the system has.
      END MODULE PARAMS

!---  cbeam add --- end!
