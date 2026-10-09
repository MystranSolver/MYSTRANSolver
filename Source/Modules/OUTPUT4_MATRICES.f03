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

      MODULE OUTPUT4_MATRICES

! This module defines the characteristics for matrices that are allowed for OUTPUT4 requests in Case Control. The definition of
! the actual matrices (with dimensions) appear in a variety of other modules (with the names listed in ALLOW_OU4_MYSTRAN_NAMES
! below. In addition, all of the output transformation matrices related to Craig-Bampton analyses are defined here. Finally, the
! allowable OUTPUT4 matrix RBM0 is defined here as a 6x6 rigid body mass matrix relative to the basic origin for the 6 DOF's at the
! boundary of a Craig-Bampton model

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                 :  BLNK_SUB_NAM, TSET_CHR_LEN

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), PARAMETER         :: NUM_OU4_VALID_NAMES = 33                   ! Num of OUTPUT4 matrices defined here
      INTEGER(LONG)                    :: NUM_OU4_REQUESTS    = 0                    ! Number of OUTPUT4 requests in a run
      INTEGER(LONG)                    :: NUM_PARTN_REQUESTS  = 0                    ! Number of PARTN requests in a run
      INTEGER(LONG), PRIVATE           :: II

! **********************************************************************************************************************************
! V A R I A B L E S   F O R   O U T P U T 4   M A T R I X   D E F I N I T I O N
! -----------------------------------------------------------------------------

      INTEGER(LONG), PARAMETER         :: TAPE_ACTION_MAX_VAL =  0                   ! See NASTRAN OUTPUT4 ITAPE allow values
      INTEGER(LONG), PARAMETER         :: TAPE_ACTION_MIN_VAL = -3                   ! See NASTRAN OUTPUT4 ITAPE allow values

      INTEGER(LONG)                    :: OU4_FILE_UNITS(NUM_OU4_VALID_NAMES)        ! File unit no's user has req for the mats
      INTEGER(LONG)                    :: OU4_TAPE_ACTION(NUM_OU4_VALID_NAMES)       ! Tape action (0 None, -1 rewind before write
!                                                                                      -2 EOF and rewind after write, -3 both)
      CHARACTER(16*BYTE)               :: ACT_OU4_MYSTRAN_NAMES(NUM_OU4_VALID_NAMES) ! MYSTRAN names: act OUTPUT4 requested mats
      CHARACTER(16*BYTE)               :: ACT_OU4_OUTPUT_NAMES(NUM_OU4_VALID_NAMES)  ! Output  names: act OUTPUT4 requested mats

                                                           ! Following are MYSTRAN names that can be requested in Exec Control for
                                                           ! OUTPUT4. Corresponding NASTRAN names (incl for compatibility between
                                                           ! MYSTRAN and NASTRAN) are the ALLOW_OU4_OUTPUT_NAMES listed later
      CHARACTER(16*BYTE), PARAMETER    :: ALLOW_OU4_MYSTRAN_NAMES(NUM_OU4_VALID_NAMES)=                      &
                                                                                       (/'CG_LTM          ', &! (  1)
                                                                                         'DLR             ', &! (  2)
                                                                                         'EIGEN_VAL       ', &! (  3)
                                                                                         'EIGEN_VEC       ', &! (  4)
                                                                                         'GEN_MASS        ', &! (  5)
                                                                                         'IF_LTM          ', &! (  6)
                                                                                         'KAA             ', &! (  7)
                                                                                         'KGG             ', &! (  8)
                                                                                         'KLL             ', &! (  9)
                                                                                         'KRL             ', &! ( 10)
                                                                                         'KRR             ', &! ( 11)
                                                                                         'KRRcb           ', &! ( 12)
                                                                                         'KXX             ', &! ( 13)
                                                                                         'LTM             ', &! ( 14)
                                                                                         'MCG             ', &! ( 15)
                                                                                         'MEFFMASS        ', &! ( 16)
                                                                                         'MPFACTOR        ', &! ( 17)
                                                                                         'MAA             ', &! ( 18)
                                                                                         'MGG             ', &! ( 19)
                                                                                         'MLL             ', &! ( 20)
                                                                                         'MRL             ', &! ( 21)
                                                                                         'MRN             ', &! ( 22)
                                                                                         'MRR             ', &! ( 23)
                                                                                         'MRRcb           ', &! ( 24)
                                                                                         'MXX             ', &! ( 25)
                                                                                         'PA              ', &! ( 26)
                                                                                         'PG              ', &! ( 27)
                                                                                         'PL              ', &! ( 28)
                                                                                         'PHIXG           ', &! ( 29)
                                                                                         'PHIZG           ', &! ( 30)
                                                                                         'RBM0            ', &! ( 31)
                                                                                         'TR6_0           ', &! ( 32)
                                                                                         'TR6_CG          ' /)! ( 33)

                                                           ! The following are allowable names Exec Control OUTPUT4 requests
                                                           ! in order to have compatibility with NASTRAN. The corresponding
                                                           ! MYSTRN names are the ALLOW_OU4_MYSTRAN_NAMES, above
      CHARACTER(16*BYTE), PARAMETER    :: ALLOW_OU4_OUTPUT_NAMES(NUM_OU4_VALID_NAMES)=                       &
                                                                                       (/'CG_LTM          ', &! (  1) for CG_LTM
                                                                                         'DM              ', &! (  2) for DLR
                                                                                         'LAMA            ', &! (  3) for EIGEN_VAL
                                                                                         'PHIG            ', &! (  4) for EIGEN_VEC
                                                                                         'MI              ', &! (  5) for GEN_MASS
                                                                                         'IF_LTM          ', &! (  6) for IF_LTM
                                                                                         'KAA             ', &! (  7) for KAA
                                                                                         'KGG             ', &! (  8) for KGG
                                                                                         'KLL             ', &! (  9) for KLL
                                                                                         'KLR(t)          ', &! ( 10) for KRL
                                                                                         'KRR             ', &! ( 11) for KRR
                                                                                         'KBB             ', &! ( 12) for KRRcb
                                                                                         'KRRGN           ', &! ( 13) for KXX
                                                                                         'LTM             ', &! ( 14) for LTM
                                                                                         'RBMCG           ', &! ( 15) for MCG
                                                                                         'MEFFMASS        ', &! ( 16) for MEFFMASS
                                                                                         'MPFACTOR        ', &! ( 17) for MPFACTOR
                                                                                         'MAA             ', &! ( 18) for MAA
                                                                                         'MGG             ', &! ( 19) for MGG
                                                                                         'MLL             ', &! ( 20) for MLL
                                                                                         'MRL             ', &! ( 21) for MRL
                                                                                         'MRN             ', &! ( 22) for MRN
                                                                                         'MRR             ', &! ( 23) for MRR
                                                                                         'MR              ', &! ( 24) for MRRcb
                                                                                         'MRRGN           ', &! ( 25) for MXX
                                                                                         'PA              ', &! ( 26) for PA
                                                                                         'PG              ', &! ( 27) for PG
                                                                                         'PL              ', &! ( 28) for PL
                                                                                         'PHIXG           ', &! ( 29) for PHIXG
                                                                                         'PHIZG           ', &! ( 30) for PHIZG
                                                                                         'RBM0            ', &! ( 31) for RBM0
                                                                                         'TR6_0           ', &! ( 32) for TR6_0
                                                                                         'RBRCG           ' /)! ( 33) for TR6_CG

      CHARACTER(LEN(BLNK_SUB_NAM)), PARAMETER    :: SUBR_WHEN_TO_WRITE_OU4_MATS(NUM_OU4_VALID_NAMES) =       &
                                                                                       (/'LINK6           ', &! (  1) for CG_LTM
                                                                                         'LINK6           ', &! (  2) for DLR
                                                                                         'LINK4           ', &! (  3) for EIGEN_VAL
                                                                                         'LINK5           ', &! (  4) for EIGEN_VEC
                                                                                         'LINK4           ', &! (  5) for GEN_MASS
                                                                                         'LINK6           ', &! (  6) for IF_LTM
                                                                                         'LINK2           ', &! (  7) for KAA
                                                                                         'LINK2           ', &! (  8) for KGG
                                                                                         'LINK2           ', &! (  9) for KLL
                                                                                         'LINK2           ', &! ( 10) for KRL
                                                                                         'LINK2           ', &! ( 11) for KRR
                                                                                         'LINK6           ', &! ( 12) for KRRcb
                                                                                         'LINK6           ', &! ( 13) for KXX
                                                                                         'LINK6           ', &! ( 14) for LTM
                                                                                         'LINK6           ', &! ( 15) for MCG
                                                                                         'LINK9           ', &! ( 16) for MEFFMASS
                                                                                         'LINK9           ', &! ( 17) for MPFACTOR
                                                                                         'LINK2           ', &! ( 18) for MAA
                                                                                         'LINK2           ', &! ( 19) for MGG
                                                                                         'LINK2           ', &! ( 20) for MLL
                                                                                         'LINK2           ', &! ( 21) for MRL
                                                                                         'LINK6           ', &! ( 22) for MRN
                                                                                         'LINK2           ', &! ( 23) for MRR
                                                                                         'LINK6           ', &! ( 24) for MRRcb
                                                                                         'LINK6           ', &! ( 25) for MXX
                                                                                         'LINK2           ', &! ( 26) for PA
                                                                                         'LINK2           ', &! ( 27) for PG
                                                                                         'LINK2           ', &! ( 28) for PL
                                                                                         'LINK5           ', &! ( 29) for PHIXG
                                                                                         'LINK5           ', &! ( 30) for PHIZG
                                                                                         'LINK6           ', &! ( 31) for RBM0
                                                                                         'LINK6           ', &! ( 32) for TR6_0
                                                                                         'LINK6           ' /)! ( 33) for TR6_CG

      CHARACTER(LEN(TSET_CHR_LEN)), PARAMETER                                                                                      &
                                       :: OU4_MAT_COL_SETS(NUM_OU4_VALID_NAMES) =     (/'--'              , &! (  1) for CG_LTM
                                                                                        'R '              , &! (  2) for DLR
                                                                                        '--'              , &! (  3) for EIGEN_VAL
                                                                                        '--'              , &! (  4) for EIGEN_VEC
                                                                                        '--'              , &! (  5) for GEN_MASS
                                                                                        '--'              , &! (  6) for IF_LTM
                                                                                        'A '              , &! (  7) for KAA
                                                                                        'G '              , &! (  8) for KGG
                                                                                        'L '              , &! (  9) for KLL
                                                                                        'L '              , &! ( 10) for KRL
                                                                                        'R '              , &! ( 11) for KRR
                                                                                        'R '              , &! ( 12) for KRRcb
                                                                                        '--'              , &! ( 13) for KXX
                                                                                        '--'              , &! ( 14) for LTM
                                                                                        '--'              , &! ( 15) for MCG
                                                                                        '--'              , &! ( 16) for MEFFMASS
                                                                                        'R '              , &! ( 17) for MPFACTOR
                                                                                        'A '              , &! ( 18) for MAA
                                                                                        'G '              , &! ( 19) for MGG
                                                                                        'L '              , &! ( 20) for MLL
                                                                                        'L '              , &! ( 21) for MRL
                                                                                        '--'              , &! ( 22) for MRN
                                                                                        'R '              , &! ( 23) for MRR
                                                                                        'R '              , &! ( 24) for MRRcb
                                                                                        '--'              , &! ( 25) for MXX
                                                                                        '--'              , &! ( 26) for PA
                                                                                        '--'              , &! ( 27) for PG
                                                                                        '--'              , &! ( 28) for PL
                                                                                        '--'              , &! ( 29) for PHIXG
                                                                                        '--'              , &! ( 30) for PHIZG
                                                                                        '--'              , &! ( 31) for RBM0
                                                                                        '--'              , &! ( 32) for TR6_0
                                                                                        '--'               /)! ( 33) for TR6_CG

      CHARACTER(LEN(TSET_CHR_LEN)), PARAMETER                                                                                      &
                                       :: OU4_MAT_ROW_SETS(NUM_OU4_VALID_NAMES) =     (/'--'              , &! (  1) for CG_LTM
                                                                                        'L '              , &! (  2) for DLR
                                                                                        '--'              , &! (  3) for EIGEN_VAL
                                                                                        'G '              , &! (  4) for EIGEN_VEC
                                                                                        '--'              , &! (  5) for GEN_MASS
                                                                                        'R '              , &! (  6) for IF_LTM
                                                                                        'A '              , &! (  7) for KAA
                                                                                        'G '              , &! (  8) for KGG
                                                                                        'L '              , &! (  9) for KLL
                                                                                        'R '              , &! ( 10) for KRL
                                                                                        'R '              , &! ( 11) for KRR
                                                                                        'R '              , &! ( 12) for KRRcb
                                                                                        '--'              , &! ( 13) for KXX
                                                                                        '--'              , &! ( 14) for LTM
                                                                                        '--'              , &! ( 15) for MCG
                                                                                        '--'              , &! ( 16) for MEFFMASS
                                                                                        '--'              , &! ( 17) for MPFACTOR
                                                                                        'A '              , &! ( 18) for MAA
                                                                                        'G '              , &! ( 19) for MGG
                                                                                        'L '              , &! ( 20) for MLL
                                                                                        'R '              , &! ( 21) for MRL
                                                                                        'R '              , &! ( 22) for MRN
                                                                                        'R '              , &! ( 23) for MRR
                                                                                        'R '              , &! ( 24) for MRRcb
                                                                                        '--'              , &! ( 25) for MXX
                                                                                        'A '              , &! ( 26) for PA
                                                                                        'G '              , &! ( 27) for PG
                                                                                        'L '              , &! ( 28) for PL
                                                                                        'G '              , &! ( 29) for PHIXG
                                                                                        'G '              , &! ( 30) for PHIZG
                                                                                        '--'              , &! ( 31) for RBM0
                                                                                        'R '              , &! ( 32) for TR6_0
                                                                                        'R '               /)! ( 33) for TR6_CG

      INTEGER(LONG), ALLOCATABLE       :: OU4_MAT_ROW_GRD_COMP(:,:)                 ! Grid/comp nos for rows in the part'd matrix
      INTEGER(LONG), ALLOCATABLE       :: OU4_MAT_COL_GRD_COMP(:,:)                 ! Grid/comp nos for cols in the part'd matrix

! 1st col of HAS_OU4_MAT_BEEN_PROCESSED is whether OU4 output request was processed. 2nd, 3rd col for whether partition request
! was processed for cols, rows

      CHARACTER(1*BYTE)                :: HAS_OU4_MAT_BEEN_PROCESSED(NUM_OU4_VALID_NAMES,3) =                                      &
                                                    RESHAPE ( (/('N', II=1,3*NUM_OU4_VALID_NAMES)/), (/NUM_OU4_VALID_NAMES,3/) )

! **********************************************************************************************************************************
! V A R I A B L E S   F O R   P A R T I T I O N I N G   O U T P U T 4   M A T R I C E S
! -------------------------------------------------------------------------------------

      CHARACTER(16*BYTE)               :: OU4_PART_MAT_NAMES(NUM_OU4_VALID_NAMES,5)  ! Names of matrices (from PARTN entry)

      CHARACTER(16*BYTE)               :: OU4_PART_VEC_NAMES(NUM_OU4_VALID_NAMES,2)  ! Row/col partition vec names for OU4 matrices

      INTEGER(LONG), ALLOCATABLE       :: OU4_PARTVEC_COL(:)                         ! Col partition vector for an OU4 matrix
      INTEGER(LONG), ALLOCATABLE       :: OU4_PARTVEC_ROW(:)                         ! Row partition vector for an OU4 matrix

! **********************************************************************************************************************************
! V A R I A B L E S   F O R   O T M's   F O R   O U T P U T 4   M A T R I C E S

! In a Craig-Bampton model generation run, Case Control output requests for displacements, element forces, stresses, etc are
! interpreted as requests for the Output Transformation Matrix (OTM) for that quantity. See Appendix D to the MYSTRAN User's
! Reference Manual for a discussion of Craig-Bampton model OTM's

      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_ACCE(:)      ! CB accel     text descriptor for rows of OTM_ACCE
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_DISP(:)      ! CB disp      text descriptor for rows of OTM_DISP
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_MPCF(:)      ! CB MPC force text descriptor for rows of OTM_MPCF
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_SPCF(:)      ! CB SPC force text descriptor for rows of OTM_SPCF

      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_ELFE(:)      ! CB elm engr force text descriptor for rows of OTM_ELFE
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_ELFN(:)      ! CB elm node force text descriptor for rows of OTM_ELFN
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_STRE(:)      ! CB elm stress     text descriptor for rows of OTM_STRE
      CHARACTER(132*BYTE), ALLOCATABLE :: TXT_STRN(:)      ! CB elm strain     text descriptor for rows of OTM_STRE

      REAL(DOUBLE), ALLOCATABLE        :: OTM_ACCE(:,:)    ! CB accel     OTM for grids req in Case Control ACCE = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_DISP(:,:)    ! CB disp      OTM for grids req in Case Control DISP = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_MPCF(:,:)    ! CB MPC force OTM for grids req in Case Control MPCF = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_SPCF(:,:)    ! CB SPC force OTM for grids req in Case Control SPCF = sid entry

      REAL(DOUBLE), ALLOCATABLE        :: OTM_ELFE(:,:)    ! CB elm engr force OTM for elems req in CC ELFORCE(ENGR) = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_ELFN(:,:)    ! CB elm node force OTM for elems req in CC ELFORCE(NODE) = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_STRE(:,:)    ! CB elm stress     OTM for elems req in CC STRESS        = sid entry
      REAL(DOUBLE), ALLOCATABLE        :: OTM_STRN(:,:)    ! CB elm strain     OTM for elems req in CC STRESS        = sid entry

! **********************************************************************************************************************************
! O T H E R

      REAL(DOUBLE)                     :: RBM0(6,6)        ! Rigid body mass matrix rel to basic origin for the 6 DOF's at the
!                                                            boundary of a Craig-Bampton model

      END MODULE OUTPUT4_MATRICES
