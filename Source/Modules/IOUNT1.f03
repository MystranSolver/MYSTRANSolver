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

      MODULE IOUNT1

! Input/output logical units, file names, status and file description (msg) variables

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  MBUG, MFIJ

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: FILE_NAM_MAXLEN = 256 ! Max length (chars) in any file name (incl folder length & ext)
      INTEGER(LONG)                   :: LEN_INPUT_FNAME       ! Length of input   dat filename up to and including the decimal pt
      INTEGER(LONG)                   :: LEN_RESTART_FNAME     ! Length of restart dat filename up to and including the decimal pt
      INTEGER(LONG), PRIVATE          :: I                     ! Only use is in implied DO loops below to initialize variables
      INTEGER(LONG), ALLOCATABLE      :: IN4FIL_NUM(:)         ! File number in E.C entry IN4
      INTEGER(LONG), PARAMETER        :: MOT4            =  9  ! Number of OT4 files allowable
      INTEGER(LONG), PARAMETER        :: MOU4            =  9  ! Number of OU4 files allowable
      INTEGER(LONG)                   :: NUM_IN4_FILES         ! Count of IN4 files specified by user
      INTEGER(LONG)                   :: LNUM_IN4_FILES        ! Max number of IN4 files allowed
      INTEGER(LONG), PARAMETER        :: NUM_OU4_FILES   = 20  ! Number of OU4 files allowable
      INTEGER(LONG), PARAMETER        :: MAX_FIL         = 72  ! Number of files (except OU4, SCR)


      SAVE

      CHARACTER(  1*BYTE)             :: BUGOUT  = 'N'         ! Y/N indicator if anything has been written to BUG file
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: DEFDIR                ! User inputs dir of where to find input (and to write output) here
      CHARACTER( 3*BYTE), PARAMETER   :: OU4_EXT(MOU4) = (/'OU1','OU2','OU3','OU4','OU5','OU6','OU7','OU8','OU9'/)
      CHARACTER( 3*BYTE), PARAMETER   :: OT4_EXT(MOT4) = (/'OT1','OT2','OT3','OT4','OT5','OT6','OT7','OT8','OT9'/)

! Following are the variable names for all files (except for units SC1, SCR, which do not have file names) used by Program

      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: BUGFIL                ! (filename.BUG) Debug file: ELDATA C.C. request elem debug info
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: EINFIL                ! (filename.EIN) Eigenvector scale factors (used to change signs)
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: ENFFIL                ! (filename.ENF) Enforced displs - file for DOF's that are enforced
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: ERRFIL                ! (filename.ERR) Error file: error messages written (and to F06)
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F06FIL                ! (filename.F06) Output file
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: IN0FIL                ! (filename.F06) Input file with all INCLUDE files
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: INFILE                ! (filename.DAT) Input file
      CHARACTER(FILE_NAM_MAXLEN*BYTE),                                                                                             &
                          ALLOCATABLE :: IN4FIL(:)             ! Names of IN4 files containing USERIN matrices
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: INCFIL                ! Names of B.D. INCLUDE files
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: INIFIL                ! Program initialization (formatted text) file (program.INI)
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1A                ! (filename.L1A) Formatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: NEUFIL                ! (filename.NEU) FEMAP neutral file
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SEQFIL                ! (filename.SEQ) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SPCFIL                ! (filename.SPC) Formatted file , see descr. below

      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F21FIL                ! (filename.F21) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F22FIL                ! (filename.F22) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F23FIL                ! (filename.F23) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F24FIL                ! (filename.F24) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: F25FIL                ! (filename.F25) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1B                ! (filename.L1B) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1C                ! (filename.L1C) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1D                ! (filename.L1D) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1E                ! (filename.L1E) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1F                ! (filename.L1F) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1G                ! (filename.L1G) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1H                ! (filename.L1H) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1I                ! (filename.L1I) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1J                ! (filename.L1J) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1K                ! (filename.L1K) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1L                ! (filename.L1L) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1M                ! (filename.L1M) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1N                ! (filename.L1N) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1O                ! (filename.L1O) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1P                ! (filename.L1P) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1Q                ! (filename.L1Q) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1R                ! (filename.L1R) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1S                ! (filename.L1S) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1T                ! (filename.L1T) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1U                ! (filename.L1U) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1V                ! (filename.L1V) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1W                ! (filename.L1W) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1X                ! (filename.L1X) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1Y                ! (filename.L1Y) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK1Z                ! (filename.L1Z) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2A                ! (filename.L2A) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2B                ! (filename.L2B) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2C                ! (filename.L2C) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2D                ! (filename.L2D) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2E                ! (filename.L2E) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2F                ! (filename.L2F) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2G                ! (filename.L2G) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2H                ! (filename.L2H) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2I                ! (filename.L2I) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2J                ! (filename.L2J) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2K                ! (filename.L2K) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2L                ! (filename.L2L) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2M                ! (filename.L2M) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2N                ! (filename.L2N) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2O                ! (filename.L2O) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2P                ! (filename.L2P) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2Q                ! (filename.L2Q) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2R                ! (filename.L2R) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2S                ! (filename.L2S) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK2T                ! (filename.L2T) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK3A                ! (filename.L3A) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK4A                ! (filename.L4A) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK4B                ! (filename.L4B) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK4C                ! (filename.L4C) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK4D                ! (filename.L4D) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK5A                ! (filename.L5A) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: LINK5B                ! (filename.L5B) Unformatted file , see descr. below

      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: OP2FIL                ! (filename.OP2) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: OU4FIL(MOU4)          ! (filename.OUi) Unformatted file , see descr. below
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: OT4FIL(MOT4)          ! (filename.OTi) Formatted file   , see descr. below

      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: RESTART_FILNAM        ! Name of file to restart

! The following STAT variables are for deciding what to do with the files when they are closed. The user can override the
! initial values in file Program.INI. Not included are status for file units SC1, SCR. Also, note that several of the file status'
! are PARAMETER. These are files that the user would not want to be deleted (BUGFIL would want to be kept but only if something
! has been written to it

      CHARACTER(  8*BYTE)             :: BUGSTAT       = 'DELETE  '  ! close status for file BUGFIL
      CHARACTER(  8*BYTE)             :: EINSTAT       = 'KEEP    '  ! close status for file EINFIL
      CHARACTER(  8*BYTE)             :: ENFSTAT       = 'KEEP    '  ! close status for file ENFFIL
      CHARACTER(  8*BYTE)             :: ERRSTAT       = 'KEEP    '  ! close status for file ERRFIL
      CHARACTER(  8*BYTE)             :: F06STAT       = 'KEEP    '  ! close status for file F06FIL
      CHARACTER(  8*BYTE)             :: IN0STAT       = 'KEEP    '  ! close status for file INFILE plus all INCLUDE files
      CHARACTER(  8*BYTE)             :: IN1STAT       = 'KEEP    '  ! close status for file INFILE
      CHARACTER(  8*BYTE)             :: IN4STAT       = 'KEEP    '  ! close status for file IN4FIL
      CHARACTER(  8*BYTE)             :: INCSTAT       = 'KEEP    '  ! close status for file INCFIL
      CHARACTER(  8*BYTE)             :: INISTAT       = 'KEEP    '  ! close status for file INIFIL
      CHARACTER(  8*BYTE)             :: L1ASTAT       = 'DELETE  '  ! close status for file LINK1A
      CHARACTER(  8*BYTE)             :: NEUSTAT       = 'DELETE  '  ! close status for file NEUFIL
      CHARACTER(  8*BYTE)             :: SEQSTAT       = 'DELETE  '  ! close status for file SEQFIL
      CHARACTER(  8*BYTE)             :: SPCSTAT       = 'DELETE  '  ! close status for file SPCFIL

      CHARACTER(  8*BYTE)             :: F21STAT       = 'DELETE  '  ! close status for file F21FIL
      CHARACTER(  8*BYTE)             :: F22STAT       = 'DELETE  '  ! close status for file F22FIL
      CHARACTER(  8*BYTE)             :: F23STAT       = 'DELETE  '  ! close status for file F23FIL
      CHARACTER(  8*BYTE)             :: F24STAT       = 'DELETE  '  ! close status for file F24FIL
      CHARACTER(  8*BYTE)             :: F25STAT       = 'DELETE  '  ! close status for file F25FIL
      CHARACTER(  8*BYTE)             :: L1BSTAT       = 'DELETE  '  ! close status for file LINK1B
      CHARACTER(  8*BYTE)             :: L1CSTAT       = 'DELETE  '  ! close status for file LINK1C
      CHARACTER(  8*BYTE)             :: L1DSTAT       = 'DELETE  '  ! close status for file LINK1D
      CHARACTER(  8*BYTE)             :: L1ESTAT       = 'DELETE  '  ! close status for file LINK1E
      CHARACTER(  8*BYTE)             :: L1FSTAT       = 'DELETE  '  ! close status for file LINK1F
      CHARACTER(  8*BYTE)             :: L1GSTAT       = 'DELETE  '  ! close status for file LINK1G
      CHARACTER(  8*BYTE)             :: L1HSTAT       = 'DELETE  '  ! close status for file LINK1H
      CHARACTER(  8*BYTE)             :: L1ISTAT       = 'DELETE  '  ! close status for file LINK1I
      CHARACTER(  8*BYTE)             :: L1JSTAT       = 'DELETE  '  ! close status for file LINK1J
      CHARACTER(  8*BYTE)             :: L1KSTAT       = 'DELETE  '  ! close status for file LINK1K
      CHARACTER(  8*BYTE)             :: L1LSTAT       = 'DELETE  '  ! close status for file LINK1L
      CHARACTER(  8*BYTE)             :: L1MSTAT       = 'DELETE  '  ! close status for file LINK1M
      CHARACTER(  8*BYTE)             :: L1NSTAT       = 'DELETE  '  ! close status for file LINK1N
      CHARACTER(  8*BYTE)             :: L1OSTAT       = 'DELETE  '  ! close status for file LINK1O
      CHARACTER(  8*BYTE)             :: L1PSTAT       = 'DELETE  '  ! close status for file LINK1P
      CHARACTER(  8*BYTE)             :: L1QSTAT       = 'DELETE  '  ! close status for file LINK1Q
      CHARACTER(  8*BYTE)             :: L1RSTAT       = 'DELETE  '  ! close status for file LINK1R
      CHARACTER(  8*BYTE)             :: L1SSTAT       = 'DELETE  '  ! close status for file LINK1S
      CHARACTER(  8*BYTE)             :: L1TSTAT       = 'DELETE  '  ! close status for file LINK1T
      CHARACTER(  8*BYTE)             :: L1USTAT       = 'DELETE  '  ! close status for file LINK1U
      CHARACTER(  8*BYTE)             :: L1VSTAT       = 'DELETE  '  ! close status for file LINK1V
      CHARACTER(  8*BYTE)             :: L1WSTAT       = 'DELETE  '  ! close status for file LINK1W
      CHARACTER(  8*BYTE)             :: L1XSTAT       = 'DELETE  '  ! close status for file LINK1X
      CHARACTER(  8*BYTE)             :: L1YSTAT       = 'DELETE  '  ! close status for file LINK1Y
      CHARACTER(  8*BYTE)             :: L1ZSTAT       = 'DELETE  '  ! close status for file LINK1Z
      CHARACTER(  8*BYTE)             :: L2ASTAT       = 'DELETE  '  ! close status for file LINK2A
      CHARACTER(  8*BYTE)             :: L2BSTAT       = 'DELETE  '  ! close status for file LINK2B
      CHARACTER(  8*BYTE)             :: L2CSTAT       = 'DELETE  '  ! close status for file LINK2C
      CHARACTER(  8*BYTE)             :: L2DSTAT       = 'DELETE  '  ! close status for file LINK2D
      CHARACTER(  8*BYTE)             :: L2ESTAT       = 'DELETE  '  ! close status for file LINK2E
      CHARACTER(  8*BYTE)             :: L2FSTAT       = 'DELETE  '  ! close status for file LINK2F
      CHARACTER(  8*BYTE)             :: L2GSTAT       = 'DELETE  '  ! close status for file LINK2G
      CHARACTER(  8*BYTE)             :: L2HSTAT       = 'DELETE  '  ! close status for file LINK2H
      CHARACTER(  8*BYTE)             :: L2ISTAT       = 'DELETE  '  ! close status for file LINK2I
      CHARACTER(  8*BYTE)             :: L2JSTAT       = 'DELETE  '  ! close status for file LINK2J
      CHARACTER(  8*BYTE)             :: L2KSTAT       = 'DELETE  '  ! close status for file LINK2K
      CHARACTER(  8*BYTE)             :: L2LSTAT       = 'DELETE  '  ! close status for file LINK2L
      CHARACTER(  8*BYTE)             :: L2MSTAT       = 'DELETE  '  ! close status for file LINK2M
      CHARACTER(  8*BYTE)             :: L2NSTAT       = 'DELETE  '  ! close status for file LINK2N
      CHARACTER(  8*BYTE)             :: L2OSTAT       = 'DELETE  '  ! close status for file LINK2O
      CHARACTER(  8*BYTE)             :: L2PSTAT       = 'DELETE  '  ! close status for file LINK2P
      CHARACTER(  8*BYTE)             :: L2QSTAT       = 'DELETE  '  ! close status for file LINK2Q
      CHARACTER(  8*BYTE)             :: L2RSTAT       = 'DELETE  '  ! close status for file LINK2R
      CHARACTER(  8*BYTE)             :: L2SSTAT       = 'DELETE  '  ! close status for file LINK2S
      CHARACTER(  8*BYTE)             :: L2TSTAT       = 'DELETE  '  ! close status for file LINK2T
      CHARACTER(  8*BYTE)             :: L3ASTAT       = 'DELETE  '  ! close status for file LINK3A
      CHARACTER(  8*BYTE)             :: L4ASTAT       = 'DELETE  '  ! close status for file LINK4A
      CHARACTER(  8*BYTE)             :: L4BSTAT       = 'DELETE  '  ! close status for file LINK4B
      CHARACTER(  8*BYTE)             :: L4CSTAT       = 'DELETE  '  ! close status for file LINK4C
      CHARACTER(  8*BYTE)             :: L4DSTAT       = 'DELETE  '  ! close status for file LINK4D
      CHARACTER(  8*BYTE)             :: L5ASTAT       = 'DELETE  '  ! close status for file LINK5A
      CHARACTER(  8*BYTE)             :: L5BSTAT       = 'DELETE  '  ! close status for file LINK5B

      CHARACTER(  8*BYTE)             :: OP2STAT       = 'KEEP    '  ! close status for file OP2FIL

                                                                     ! close status for file OU4FIL's
      CHARACTER(  8*BYTE)             :: OU4STAT(MOU4) = (/('DELETE', I=1,MOU4)/)
      CHARACTER(  8*BYTE)             :: OT4STAT(MOT4) = (/('DELETE', I=1,MOT4)/)

      CHARACTER(  8*BYTE)             :: BUGSTAT_OLD   = 'DELETE  '  ! close status for file BUGFIL for use in restart
      CHARACTER(  8*BYTE)             :: ERRSTAT_OLD   = 'DELETE  '  ! close status for file ERRFIL for use in restart

! The following are messages that describe what the files are (no message for SC1)

      CHARACTER( 64*BYTE)             :: BUG_MSG       = 'ELEMENT DEBUG OUTPUT FILE'
      CHARACTER( 64*BYTE)             :: EIN_MSG       = 'EIGENVEC NUMBERS FOR SIGN CHANGE'
      CHARACTER( 64*BYTE)             :: ENF_MSG       = 'ENFORCED DISPL FOR ALL DOFs FILE'
      CHARACTER( 64*BYTE)             :: ERR_MSG       = 'ERROR FILE'
      CHARACTER( 64*BYTE)             :: F06_MSG       = 'PROGRAM OUTPUT DATA FILE'
      CHARACTER( 64*BYTE)             :: IN0_MSG       = 'PROGRAM INPUT DATA FILE WITH ALL INCLUDE FILES'
      CHARACTER( 64*BYTE)             :: IN1_MSG       = 'PROGRAM INPUT DATA FILE'
      CHARACTER( 64*BYTE)             :: IN4_MSG       = 'IN4 INPUT MATRICES FILE'
      CHARACTER( 64*BYTE)             :: INC_MSG       = 'BULK DATA INCLUDE FILES'
      CHARACTER( 64*BYTE)             :: INI_MSG       = 'PROGRAM INITIALIZATION FILE'
      CHARACTER( 64*BYTE)             :: L1A_MSG       = 'LINK1A TEXT FILE (I/O UNITS, COUNTERS, ETC)'
      CHARACTER( 64*BYTE)             :: NEU_MSG       = 'FEMAP NEUTRAL FILE'
      CHARACTER( 64*BYTE)             :: SEQ_MSG       = 'BANDIT SEQGP CARD IMAGES'
      CHARACTER( 64*BYTE)             :: SPC_MSG       = 'SPC1 TEXT FILE'

      CHARACTER( 64*BYTE)             :: F21_MSG       = 'ELEM THERM & PRESS LOADS DISK FILE - IN LOCAL ELEM COORDS'
      CHARACTER( 64*BYTE)             :: F22_MSG       = 'ELEM MASS DISK FILE - IN LOCAL ELEM COORDS'
      CHARACTER( 64*BYTE)             :: F23_MSG       = 'ELEM STIFF MATRIX DISK FILE - IN LOCAL ELEM COORDS'
      CHARACTER( 64*BYTE)             :: F24_MSG       = 'ELEM STRESS RECOVERY MATRICES DISK FILE - IN LOCAL ELEM COORDS'
      CHARACTER( 64*BYTE)             :: F25_MSG       = 'ELEM NODAL DISPL & LOADS DISK FILE - IN LOCAL ELEM COORDS'
      CHARACTER( 64*BYTE)             :: L1B_MSG       = 'GRID, GRID SEQ, AND COORD DATA'
      CHARACTER( 64*BYTE)             :: L1C_MSG       = 'DOF TABLES'
      CHARACTER( 64*BYTE)             :: L1D_MSG       = 'SUBCASE DATA'
      CHARACTER( 64*BYTE)             :: L1E_MSG       = 'PG LOAD MATRIX'
      CHARACTER( 64*BYTE)             :: L1F_MSG       = 'RIGID ELEMENT DATA'
      CHARACTER( 64*BYTE)             :: L1G_MSG       = 'ELEM CONN, PROP, MATL, ETC'
      CHARACTER( 64*BYTE)             :: L1H_MSG       = 'ENFORCED DISPLACEMENTS'
      CHARACTER( 64*BYTE)             :: L1I_MSG       = 'FORCE, MOMENT DATA'
      CHARACTER( 64*BYTE)             :: L1J_MSG       = 'RMG MATRIX'
      CHARACTER( 64*BYTE)             :: L1K_MSG       = 'ELEM TEMPERATURE DATA'
      CHARACTER( 64*BYTE)             :: L1L_MSG       = 'KGG STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L1M_MSG       = 'EIGENVALUE EXTRACTION DATA'
      CHARACTER( 64*BYTE)             :: L1N_MSG       = 'ASET, OMIT DATA'
      CHARACTER( 64*BYTE)             :: L1O_MSG       = 'SPC, SPC1 DATA'
      CHARACTER( 64*BYTE)             :: L1P_MSG       = 'GRAV LOAD DATA'
      CHARACTER( 64*BYTE)             :: L1Q_MSG       = 'ELEM PRESSURE DATA'
      CHARACTER( 64*BYTE)             :: L1R_MSG       = 'MGG MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L1S_MSG       = 'MPC GRIDS, COMPONENTS, EQN COEFFICIENTS'
      CHARACTER( 64*BYTE)             :: L1T_MSG       = 'SUPORT DATA'
      CHARACTER( 64*BYTE)             :: L1U_MSG       = 'RFORCE LOAD DATA'
      CHARACTER( 64*BYTE)             :: L1V_MSG       = 'PARTVEC, PARTVEC1 DATA'
      CHARACTER( 64*BYTE)             :: L1X_MSG       = 'USET, USET1 DATA'
      CHARACTER( 64*BYTE)             :: L1W_MSG       = 'SCALAR LOAD DATA'
      CHARACTER( 64*BYTE)             :: L1Y_MSG       = 'CONCENTRATED AND SCALAR MASS DATA'
      CHARACTER( 64*BYTE)             :: L1Z_MSG       = 'CHKPNT DATA'
      CHARACTER( 64*BYTE)             :: L2A_MSG       = 'GMN CONSTRAINT MATRIX'
      CHARACTER( 64*BYTE)             :: L2B_MSG       = 'KSF STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L2C_MSG       = 'QSYS ENFORCED DISPL EQUIVALENT LOADS'
      CHARACTER( 64*BYTE)             :: L2D_MSG       = 'PS LOAD MATRIX'
      CHARACTER( 64*BYTE)             :: L2E_MSG       = 'GOA MATRIX'
      CHARACTER( 64*BYTE)             :: L2F_MSG       = 'UO0 MATRIX'
      CHARACTER( 64*BYTE)             :: L2G_MSG       = 'KLL STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L2H_MSG       = 'PL LOAD MATRIX'
      CHARACTER( 64*BYTE)             :: L2I_MSG       = 'MLL MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L2J_MSG       = 'MPC FORCE RECOVERY STIFF MATRIX, HMN'
      CHARACTER( 64*BYTE)             :: L2K_MSG       = 'KRL STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L2L_MSG       = 'KRR STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L2M_MSG       = 'MRL MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L2N_MSG       = 'MRR MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L2O_MSG       = 'KAA STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L2P_MSG       = 'MAA MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L2Q_MSG       = 'PA  LOAD MATRIX'
      CHARACTER( 64*BYTE)             :: L2R_MSG       = 'MPC FORCE RECOVERY MASS MATRIX, LMN'
      CHARACTER( 64*BYTE)             :: L2S_MSG       = 'MSF MASS MATRIX'
      CHARACTER( 64*BYTE)             :: L2T_MSG       = 'KLLs STIFFNESS MATRIX'
      CHARACTER( 64*BYTE)             :: L3A_MSG       = 'UL DISPL VECTORS'
      CHARACTER( 64*BYTE)             :: L4A_MSG       = 'ULLI - inverse of KMSM triangular factor ULL'
      CHARACTER( 64*BYTE)             :: L4B_MSG       = 'ALL L-set matrix for the standard eigenproblem'
      CHARACTER( 64*BYTE)             :: L4C_MSG       = 'LINK4C'
      CHARACTER( 64*BYTE)             :: L4D_MSG       = 'LINK4D'
      CHARACTER( 64*BYTE)             :: L5A_MSG       = 'UG DISPL VECTORS'
      CHARACTER( 64*BYTE)             :: L5B_MSG       = 'COLS OF PHIXA EXPANDED TO G-SET'

      CHARACTER( 64*BYTE)             :: OP2_MSG       = 'OP2 PROGRAM OUTPUT FILE'
      CHARACTER( 64*BYTE)             :: OU4_MSG(MOU4) = (/'OUTPUT4 FILE 21                                                 ',     &
                                                           'OUTPUT4 FILE 22                                                 ',     &
                                                           'OUTPUT4 FILE 23                                                 ',     &
                                                           'OUTPUT4 FILE 24                                                 ',     &
                                                           'OUTPUT4 FILE 25                                                 ',     &
                                                           'OUTPUT4 FILE 26                                                 ',     &
                                                           'OUTPUT4 FILE 27                                                 ',     &
                                                           'OUTPUT4 FILE 28                                                 ',     &
                                                           'OUTPUT4 FILE 29                                                 '/)
      CHARACTER( 64*BYTE)             :: OT4_MSG(MOT4) = (/'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 21        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 22        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 23        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 24        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 25        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 26        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR OUTPUT4 MATRICES IN OU4 FILE 27        ',     &
                                                           'TEXT FILE DESCRIPTOR FOR GRID RELATED OTMs                      ',     &
                                                           'TEXT FILE DESCRIPTOR FOR ELEM RELATED OTMs                      '/)

! The following are unit numbers for all files used by the program

      INTEGER(LONG)                   :: SC1           =    6 ! Unit no. for screen

      INTEGER(LONG)                   :: BUG           =    2 ! Unit no. for debug output file
      INTEGER(LONG)                   :: EIN           = 1001 ! Unit no. for text file w/ eigenvec scale facs (if supplied)
      INTEGER(LONG)                   :: ENF           = 1002 ! Unit no. for text file w/ enforced displ for all grids/comps
      INTEGER(LONG)                   :: ERR           =    3 ! Unit no. for error file
      INTEGER(LONG)                   :: F06           =    7 ! Unit no. for output file
      INTEGER(LONG)                   :: IN0           = 1003 ! Unit no. for input file with all INCLUDE files
      INTEGER(LONG)                   :: IN1           =    8 ! Unit no. for input file
      INTEGER(LONG)                   :: IN4           = 1004 ! Unit no. for IN4file
      INTEGER(LONG)                   :: INC           = 1005 ! Unit no. for Bulk Data INCLUDE files

! BANDIT file unit nos: IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20 defined in DATA stmt
! in BANDIT_BLOCK_DATA. Values are 1006-1020

      INTEGER(LONG)                   :: INI           =    9 ! Unit no. for MYSTRAN.INI file
      INTEGER(LONG)                   :: L1A           =  101 ! Unit no. for file LINK1A
      INTEGER(LONG)                   :: NEU           =   10 ! Unit no. for FEMPA neutral file
      INTEGER(LONG)                   :: SEQ           =   12 ! Unit no. for SEQGP card images from BANDIT
      INTEGER(LONG)                   :: SPC           =   13 ! Unit no. for SPC1 text file

      INTEGER(LONG)                   :: F21           =   91 ! Unit no. for file F21FIL
      INTEGER(LONG)                   :: F22           =   92 ! Unit no. for file F22FIL
      INTEGER(LONG)                   :: F23           =   93 ! Unit no. for file F23FIL
      INTEGER(LONG)                   :: F24           =   94 ! Unit no. for file F24FIL
      INTEGER(LONG)                   :: F25           =   95 ! Unit no. for file F24FIL
      INTEGER(LONG)                   :: L1B           =  102 ! Unit no. for file LINK1B
      INTEGER(LONG)                   :: L1C           =  103 ! Unit no. for file LINK1C
      INTEGER(LONG)                   :: L1D           =  104 ! Unit no. for file LINK1D
      INTEGER(LONG)                   :: L1E           =  105 ! Unit no. for file LINK1E
      INTEGER(LONG)                   :: L1F           =  106 ! Unit no. for file LINK1F
      INTEGER(LONG)                   :: L1G           =  107 ! Unit no. for file LINK1G
      INTEGER(LONG)                   :: L1H           =  108 ! Unit no. for file LINK1H
      INTEGER(LONG)                   :: L1I           =  109 ! Unit no. for file LINK1I
      INTEGER(LONG)                   :: L1J           =  110 ! Unit no. for file LINK1J
      INTEGER(LONG)                   :: L1K           =  111 ! Unit no. for file LINK1K
      INTEGER(LONG)                   :: L1L           =  112 ! Unit no. for file LINK1L
      INTEGER(LONG)                   :: L1M           =  113 ! Unit no. for file LINK1M
      INTEGER(LONG)                   :: L1N           =  114 ! Unit no. for file LINK1N
      INTEGER(LONG)                   :: L1O           =  115 ! Unit no. for file LINK1O
      INTEGER(LONG)                   :: L1P           =  116 ! Unit no. for file LINK1P
      INTEGER(LONG)                   :: L1Q           =  117 ! Unit no. for file LINK1Q
      INTEGER(LONG)                   :: L1R           =  118 ! Unit no. for file LINK1R
      INTEGER(LONG)                   :: L1S           =  119 ! Unit no. for file LINK1S
      INTEGER(LONG)                   :: L1T           =  120 ! Unit no. for file LINK1S
      INTEGER(LONG)                   :: L1U           =  121 ! Unit no. for file LINK1U
      INTEGER(LONG)                   :: L1V           =  122 ! Unit no. for file LINK1V
      INTEGER(LONG)                   :: L1W           =  123 ! Unit no. for file LINK1W
      INTEGER(LONG)                   :: L1X           =  124 ! Unit no. for file LINK1X
      INTEGER(LONG)                   :: L1Y           =  125 ! Unit no. for file LINK1Y
      INTEGER(LONG)                   :: L1Z           =  126 ! Unit no. for file LINK1Z
      INTEGER(LONG)                   :: L2A           =  201 ! Unit no. for file LINK2A
      INTEGER(LONG)                   :: L2B           =  202 ! Unit no. for file LINK2B
      INTEGER(LONG)                   :: L2C           =  203 ! Unit no. for file LINK2C
      INTEGER(LONG)                   :: L2D           =  204 ! Unit no. for file LINK2D
      INTEGER(LONG)                   :: L2E           =  205 ! Unit no. for file LINK2E
      INTEGER(LONG)                   :: L2F           =  206 ! Unit no. for file LINK2F
      INTEGER(LONG)                   :: L2G           =  207 ! Unit no. for file LINK2G
      INTEGER(LONG)                   :: L2H           =  208 ! Unit no. for file LINK2H
      INTEGER(LONG)                   :: L2I           =  209 ! Unit no. for file LINK2I
      INTEGER(LONG)                   :: L2J           =  210 ! Unit no. for file LINK2J
      INTEGER(LONG)                   :: L2K           =  211 ! Unit no. for file LINK2K
      INTEGER(LONG)                   :: L2L           =  212 ! Unit no. for file LINK2L
      INTEGER(LONG)                   :: L2M           =  213 ! Unit no. for file LINK2M
      INTEGER(LONG)                   :: L2N           =  214 ! Unit no. for file LINK2N
      INTEGER(LONG)                   :: L2O           =  215 ! Unit no. for file LINK2O
      INTEGER(LONG)                   :: L2P           =  216 ! Unit no. for file LINK2P
      INTEGER(LONG)                   :: L2Q           =  217 ! Unit no. for file LINK2Q
      INTEGER(LONG)                   :: L2R           =  218 ! Unit no. for file LINK2R
      INTEGER(LONG)                   :: L2S           =  219 ! Unit no. for file LINK2S
      INTEGER(LONG)                   :: L2T           =  220 ! Unit no. for file LINK2T
      INTEGER(LONG)                   :: L3A           =  301 ! Unit no. for file LINK3A
      INTEGER(LONG)                   :: L4A           =  401 ! Unit no. for file LINK4A
      INTEGER(LONG)                   :: L4B           =  402 ! Unit no. for file LINK4B
      INTEGER(LONG)                   :: L4C           =  403 ! Unit no. for file LINK4B
      INTEGER(LONG)                   :: L4D           =  404 ! Unit no. for file LINK4B
      INTEGER(LONG)                   :: L5A           =  501 ! Unit no. for file LINK5A
      INTEGER(LONG)                   :: L5B           =  502 ! Unit no. for file LINK5B

      INTEGER(LONG)                   :: OP2           =   14 ! Unit no. for file OP2FIL

                                                              ! Unit no's for OU4FIL files
      INTEGER(LONG)                   :: OU4(MOU4)     = (/21,22,23,24,25,26,27,28,29/)
      INTEGER(LONG)                   :: OT4(MOT4)     = (/31,32,33,34,35,36,37,38,39/)
      INTEGER(LONG)                   :: OT4_ELM_OTM   =   38 ! OT4 unit above reserved for elem related CB OTM text msgs
      INTEGER(LONG)                   :: OT4_GRD_OTM   =   39 ! OT4 unit above reserved for grid related CB OTM text msgs
      INTEGER(LONG)                   :: OU4_ELM_OTM   =   28 ! OU4 unit above reserved for elem related CB OTM values
      INTEGER(LONG)                   :: OU4_GRD_OTM   =   29 ! OU4 unit above reserved for grid related CB OTM values

                                                              ! Unit no's. for  scratch files
      INTEGER(LONG)                   :: SCR(9)        = (/991,992,993,994,995,996,997,998,999/)

! The following are indicators of whether to write to BUG, ERR

      INTEGER(LONG)                   :: WRT_BUG(0:MBUG-1)  = (/(0, I=0,MBUG-1)/)
                                                               ! WRT_BUG specifies what to write to the BUG file. Set by C.C. ELDATA

      INTEGER(LONG)                   :: WRT_FIJ(MFIJ)      = (/(0, I=1,MFIJ)/)
                                                               ! WRT_BUG specifies what to write to the BUG file. Set by C.C. ELDATA

      INTEGER(LONG)                   :: WRT_ERR = 1           ! WRT_ERR says whether to write ERR file or not


! Description of files:
! ---------------------

! BUGFIL is a formatted file containing element data written if ELDATA Case Control requests are made

! EINFIL is a formatted file containing scale factors to apply to eigenvectors (after renormalization according to EIGR/EIGRL
!        Bulk Data entry). This file will only be used if Bulk Data PARAM EIGNORM2 is present with value 'Y'

! ENFFIL is a text file containing enforced displs for all grids/components (used when Case Control ENFORCED is present)

! ERRFIL is a formatted file containing all warning and error messages (also written to F06FIL)

! F06FIL is a formatted file containing the normal output from MYSTRAN

! OP2FIL is an unformatted file containing:
!      the MSC Nastran formatted OP2 (basically the binary version of the F06)

! F21FIL is an unformatted file containing:
!      Array ME (element mass) for elements requested in ELDATA Case Control command
!      See module MODEL_STUF for description. Written by subr WEOFIL (called by subr's ELEM_MASS & EPTL in LINK1)

! F22FIL is an unformatted file containing:
!      ArrayS PTE, PPE (element thermal & pressure loads) for elements requested in ELDATA Case Control command
!      See module MODEL_STUF for description. Written by subr WEOFIL (called by subr's ELEM_MASS & EPTL in LINK1)

! F23FIL is an unformatted file containing:
!      Array KE (element stiffness) for elements requested in ELDATA Case Control command
!      See module MODEL_STUF for description. Written by subr WEOFIL (called by subr ESP in LINK1)

! F24FIL is an unformatted file containing:
!      Arrays SE1, SE2, SE3, STE1, STE2, STE3 (elem stress recovery matrices) for elems requested in ELDATA Case Control command
!      See module MODEL_STUF for description. Written by subr WEOFIL (called by subr ESP in LINK1)

! F25FIL is an unformatted file containing:
!      Arrays UEL, PEL (elem displ's and nodal loads) for elements requested in ELDATA Case Control command
!      See module MODEL_STUF for description. Written by subr WEOFIL (called by subr OFP3 in LINK9)

! INFILE is a formatted file containing the normal input

! INIFIL is a formatted file containing initialization variables for the program

! LINK1A is a formatted file first created in the LINK1. It is read at the beginning of each subsequent LINK and rewritten at the
! end of each link. This file MUST be present to start any other LINK in the Program execution after the 1st one has run.

!    LINK number
!    SOL name                                    See module SCONTR for description.
!    I/O unit numbers, names and close status'   See this   module for description
!    Array sizes, etc                            See module SCONTR for description.
!    some PARAMS
!    COMM

! LINK1B is an unformatted file containing:
!    Arrays GRID, RGRID                          See module MODEL_STUF for description.  Written in subr GRID_PROC
!    Arrays CORD, RCORD                          See module MODEL_STUF for description.  Written in subr GRID_PROC
!    Array  GRID_SEQ                             See module MODEL_STUF for description.  Written in subr SEQ_PROC
!    Arrays SEQ1, SEQ2                           See module MODEL_STUF for description.  Written in subr SEQ_PROC

! LINK1C is an unformatted file containing:
!    Arrays TSET, TDOFI, TDOF                    See module DOF_TABLES for description.  Written in subr DOF_PROC

! LINK1D is an unformatted file containing:
!    Arrays SCNUM, TITLE, STITLE, LABEL          See module MODEL_STUF for description.  Written in subr SUBCASE_PROC
!    Array  SUBLOD                               See module MODEL_STUF for description.  Written in subr SUBCASE_PROC
!    Arrays OGROUT, GROUT                        See module MODEL_STUF for description.  Written in subr SUBCASE_PROC
!    Arrays OELOUT, ELOUT                        See module MODEL_STUF for description.  Written in subr SUBCASE_PROC
!    Arrays OELDT, ELDT                          See module MODEL_STUF for description.  Written in subr SUBCASE_PROC

! LINK1E is an unformatted file containing G-set loads for all subcases:
!    Array SYS_LOAD  (NDOFG x NSUB)              See module SPARSE_PG for description. Written in LINK1

! LINK1F is an unformatted file containing rigid element data:
!    After the B.D. deck is read, L1F contains 2 records for each rigid element:
!    The 1st record is the element type ('RBAR    ', 'RBE1    ', or 'RBE2    ')
!    The 2nd record contains the rigid element ID and all other element data (see subr's BD_RBAR, BD_RBE1, BD_RBE2)
!

! LINK1G is an unformatted file containing:
!    Element data:                               See module MODEL_STUF  for description. Written in subr ELSAVE
!      Arrays ETYPE, EPNT, ESORT1, ESORT2, EOFF
!      Array  EDAT
!      Elem integration order and thick plate constants (MIN3, MIN4)
!      Array  VVEC of element V vectors
!      Array  BAROFF
!      Array  PLATEOFF
!      Array  PLATETHICK
!      Arrays PBAR  , RPBAR
!      Arrays PBEAM , RPBEAM
!      Arrays PROD  , RPROD
!      Arrays PELAS , RPELAS
!      Arrays PSHEL , RPSHEL
!      Arrays PCOMP , RPCOMP
!      Array  PSOLID
!      Arrays PUSER1, RPUSER1
!      Array  PUSERIN
!      Array  USERIN_MAT_NAMES
!      Arrays MAT1  , RMATL
!      Array  MATANGLE

! LINK1H is an unformatted file containing enforced displacement data:
!      After subr DOF_PROC, LINK1H contains enforced SPC data: grid no., component (1,2,3,4,5 or 6), SPC value
!      After subr YS_ARRAY, LINK1H has array YS. See module SPARSE_MATRICES for description.  Written in subr YS_ARRAY

! LINK1I is an unformatted file containing FORCE/MOMENT B.D. card data:
!      When Bulk Data deck is read, any FORCE or MOMENT card data is writtrn to LINK1I for later processing

! LINK1J is an unformatted file containing:
!      After subr RIGID_ELEM_PROC, RMG constraint matrix coefficients.

! LINK1K is an unformatted file containing element and grid point temperature data:
!      After the B.D. deck is read, LINK1K has TEMP, TEMPD, TEMPRB, TEMPP1 card data
!      After subr TEMPERATURE_PROC has run, LINK1K has:
!         Arrays TPNT, TDATA                     See module MODEL_STUF  for description. Written in subr TEMPERATURE_PROC
!         Arrays GTEMP, CGTEMP                   See module MODEL_STUF  for description. Written in subr TEMPERATURE_PROC

! LINK1L is an unformatted file containing the G-set stiffness matrix:
!      For each row: row num, col num, KGG val.  See module STF_ARRAYS  for description. Written in subr SPARSE_KGG

! LINK1M is an unformatted file that has eigenvalue extraction data read in LINK1 and generalized masses from LINK4.
!      After the B.D. deck is read, LINK1M has the 10 items read from a B.D. EIGR card:
!         EIG_SID, EIG_METH, EIG_FRQ1, EIG_FRQ2, EIG_N1, EIG_N2, EIG_CRIT, EIG_NORM, EIG_GRID, EIG_COMP
!      After LINK4 has run, file LINK1M has:
!         One record for each of: EIG_SID, EIG_METH, EIG_FRQ1, EIG_FRQ2, EIG_N1, EIG_N2, EIG_CRIT, EIG_NORM, EIG_GRID, EIG_COMP
!         MAXMIJ = largest off-diag term in the gen mass matrix
!         One record for each vector containing: mode number, extraction order, eigenvalue, radian freq, cycles, gen mass, gen stiff

! LINK1N is an unformatted file containing ASET/OMIT data:
!      After the B.D. deck is read, LINK1N has ASET, ASET1, OMIT, OMIT1 information (DOF, GRID range, A/O set)

! LINK1O is an unformatted file containing SPC, SPC1 data. This file is read in subr DOF_PROC:
!      After the B.D. deck is read, LINK1O has Set ID, DOF, GRID no., enforced displ (0. or value)

! LINK1P is an unformatted file containing gravity load data that is read in subr GRAV_PROC:
!      After the B.D. deck is read, LINK1P has:
!               SETID         = Load set ID
!               ACID_L        = Local coord sys ID that gravity load is given in
!               GRAV_GRID     = ID of grid that rotational (components 4, 5, 6) grav accels are about
!               ACCEL_RB(1-6) = 6 components of gravity load

! LINK1Q is an unformatted file containing element pressure load data:
!      After the B.D. deck is read, LINK1K has PLOAD1 and PLOAD2 card data
!      After subr PRESSURE_DATA_PROC has run, LINK1K has:
!      Arrays PPNT, PDATA, PTYPE                 See module MODEL_STUF for description. Written in subr PRESSURE_DATA_PROC

! LINK1R is an unformatted file containing mass data:
!      Arrays I_MGG, J_MGG, MGG                  See module SPARSE_MATRICES  for description. Written in LINK1

! LINK1S is an unformatted file containing MPC data:
!      After the B.D. deck is read, LINK1S has the following records for each MPC logical card:
!         MPC set ID
!         Number of triplets of grid/component/MPC coeff for this MPC set
!         One record containing grid, component, MPC coefficient for each triplet in the MPC set

! LINK1T is an unformatted file containing SUPORT data:
!      After the B.D. deck is read, LINK1T has grid, component values for SUPORT'd DOF (one pair/record):

! LINK1U is an unformatted file containing RFORCE load data that is read in subr RFORCE_PROC:
!      After the B.D. deck is read, LINK1P has:
!               SETID         = Load set ID
!               ACID_L        = Local coord sys ID that RFORCE load is given in
!               GID           = ID of grid that the angular velocity and/or angular accel is about
!               SCALE_AV      = Scale factor for angular velocity
!               SCALE_AA      = Scale factor for angular accel
!               ACCEL_RB(1-3) = 3 components of RFORCE vector

! LINK1V is an unformatted file containing PARTVEC, PARTVEC1 data. This file is read in subr PARTVEC_PROC:
!      After the B.D. deck is read, LINK1V has Set name, DOF, GRID no.

! LINK1W is an unformatted file containing SLOAD load data that is read in subr SLOAD_PROC:
!      After the B.D. deck is read, LINK1W has (for each scalar load):
!               SETID         = Load set ID
!               SPOINT        = Scalar point at which the scalar load acts
!               SLOAD         = Scalar load magnitude

! LINK1X is an unformatted file containing USET, USET1 data. This file is read in subr DOF_PROC:
!      After the B.D. deck is read, LINK1X has Set name, DOF, GRID no.

! LINK1Y is an unformatted file containing concentrated and scalar mass data
!      After subr CONM2_PROC_1 has run, LINK1Y has:
!

! LINK1Z is an unformatted file containing CHKPNT data:
!      After the C.C deck is read, LINK1Z has the following records
!         NSUB
!         MPCSET
!         SPCSET
!         SUBLOD(I,1), SUBLOD(I,2) one record for each I=1,NSUB

! LINK2A is an unformatted file containing the GMN constraint matrix. Written in LINK2.

! LINK2B is an unformatted file containing the KSF stiffness matrix partition.

! LINK2C is an unformatted file containing the QSYS equivalent loads due to enforced displ's. QSYS = KFS(transpose) x YS

! LINK2D is an unformatted file containing the PS loads on SPC'd-DOf's

! LINK2E is an unformatted file containing the GOA matrix

! LINK2F is an unformatted file containing UO0 matrix UO0 = KOO(-1) x KAO(transpose)

! LINK2G is an unformatted file containing the reduced KLL stiffness matrix for the L-set

! LINK2H is an unformatted file containing the reduced PL load matrix for the L-set

! LINK2I is an unformatted file containing the reduced MLL mass matrix for the L-set

! LINK2J is an unformatted file containing the MPC force recovery matrix HMN = (KNMt + KMM*GMN)

! LINK2K is an unformatted file containing the reduced KRL stiffness matrix

! LINK2L is an unformatted file containing the reduced KRR stiffness matrix

! LINK2M is an unformatted file containing the reduced MRL mass matrix

! LINK2N is an unformatted file containing the reduced MRR mass matrix

! LINK2O is an unformatted file containing the reduced KAA stiffness matrix

! LINK2P is an unformatted file containing the reduced MAA mass matrix

! LINK2Q is an unformatted file containing the reduced PA  load matrix

! LINK3A is an unformatted file containing A-set displacement solution data from:
!      LINK3 - for SOL=1, the A-set displacements
!      LINK4 - for SOL=3, the A-set eigenvectors

! LINK4A is an unformated file containing ULLI, the inverse of the triangular factor of KMSM

! LINK4B is an unformated file containing ALL, the symmetric standard eigenproblem matrix

! LINK5A is an unformated file containing G-set displacement solution data from LINK5

! LINK5B is an unformated file containing cols of PHIXA expanded to G-set size from LINK5

! NEUFIL is a formatted file containing data for FEMAP

! IN4FIL is an unformatted file of matrix headers, matrix data in NASTRAN INPUT4 format

! OU4FIL is an unformatted file of matrix headers, matrix data in NASTRAN OUTPUT4 format

! SEQFIL is a formatted file containing SEQGP card images generated in subr BANDIT (automatic grid sequencer)

! SPCFIL is a text file containing the B.D. SPC1 card images from the grid point singularity processor

      END MODULE IOUNT1
