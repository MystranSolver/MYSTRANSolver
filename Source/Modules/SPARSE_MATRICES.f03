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

      MODULE SPARSE_MATRICES

! Arrays for sparse stiffness, mass, load, constraint and some other matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

! **********************************************************************************************************************************
! G-set arrays:

      CHARACTER(  1*BYTE)             :: SYM_KGG     = 'x' ! Matrix KGG     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KGGD    = 'x' ! Matrix KGGD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MGG     = 'x' ! Matrix MGG     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MGGC    = 'x' ! Matrix MGGC    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MGGE    = 'x' ! Matrix MGGE    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MGGS    = 'x' ! Matrix MGGE    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PG      = 'x' ! Matrix PG      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_RMG     = 'x' ! Matrix RMG     sparse storage method (SYM or NONSYM)

      INTEGER(LONG), ALLOCATABLE      :: I_KGG(:)          ! Row indicators for nonzero terms in stiff matrix KGG
      INTEGER(LONG), ALLOCATABLE      :: I_KGGD(:)         ! Row indicators for nonzero terms in stiff matrix KGGD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_MGG(:)          ! Row indicators for nonzero terms in mass matrix MGG
      INTEGER(LONG), ALLOCATABLE      :: I_MGGE(:)         ! Row indicators for nonzero terms in mass matrix MGGE
      INTEGER(LONG), ALLOCATABLE      :: I_MGGC(:)         ! Row indicators for nonzero terms in mass matrix MGGE
      INTEGER(LONG), ALLOCATABLE      :: I_MGGS(:)         ! Row indicators for nonzero terms in mass matrix MGGS
      INTEGER(LONG), ALLOCATABLE      :: I_PG(:)           ! Row indicators for nonzero terms in load matrix PG
      INTEGER(LONG), ALLOCATABLE      :: I_RMG(:)          ! Row indicators for nonzero terms in stiff matrix RMG
      INTEGER(LONG), ALLOCATABLE      :: I2_MGG(:)         ! Array of row numbers for all nonzero terms in matrix MGG

      INTEGER(LONG), ALLOCATABLE      :: J_KGG(:)          ! Col numbers for nonzero terms in stiff matrix KGG
      INTEGER(LONG), ALLOCATABLE      :: J_KGGD(:)         ! Col numbers for nonzero terms in stiff matrix KGG-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_MGG(:)          ! Col numbers for nonzero terms in mass matrix MGG
      INTEGER(LONG), ALLOCATABLE      :: J_MGGC(:)         ! Col numbers for nonzero terms in mass matrix MGGE
      INTEGER(LONG), ALLOCATABLE      :: J_MGGE(:)         ! Col numbers for nonzero terms in mass matrix MGGE
      INTEGER(LONG), ALLOCATABLE      :: J_MGGS(:)         ! Col numbers for nonzero terms in mass matrix MGGS
      INTEGER(LONG), ALLOCATABLE      :: J_PG(:)           ! Col numbers for nonzero terms in load matrix PG
      INTEGER(LONG), ALLOCATABLE      :: J_RMG(:)          ! Col numbers for nonzero terms in stiff matrix RMG

      REAL(DOUBLE) , ALLOCATABLE      :: KGG(:)            ! G-set nonzero stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KGGD(:)           ! G-set nonzero stiff terms-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: MGG(:)            ! G-set nonzero mass terms for elem and conc masses
      REAL(DOUBLE) , ALLOCATABLE      :: MGGC(:)           ! G-set nonzero mass terms due to concentrated masses
      REAL(DOUBLE) , ALLOCATABLE      :: MGGE(:)           ! G-set nonzero mass terms due to element mass
      REAL(DOUBLE) , ALLOCATABLE      :: MGGS(:)           ! G-set nonzero mass terms due to scalar masses
      REAL(DOUBLE) , ALLOCATABLE      :: PG(:)             ! G-set applied loads
      REAL(DOUBLE) , ALLOCATABLE      :: RMG(:)            ! Nonzero terms in constraint matrix (from rigid elems and MPC's)

! **********************************************************************************************************************************
! N, M-set arrays

      CHARACTER(  1*BYTE)             :: SYM_KNN     = 'x' ! Matrix KNN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KNM     = 'x' ! Matrix KNM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMM     = 'x' ! Matrix KMM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMN     = 'x' ! Matrix KMN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KNND    = 'x' ! Matrix KNND    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KNMD    = 'x' ! Matrix KNMD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMMD    = 'x' ! Matrix KMMD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMND    = 'x' ! Matrix KMND    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MNN     = 'x' ! Matrix MNN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MNM     = 'x' ! Matrix MNM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MMN     = 'x' ! Matrix MMN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MMM     = 'x' ! Matrix MMM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PN      = 'x' ! Matrix PN      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PM      = 'x' ! Matrix PM      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_RMN     = 'x' ! Matrix RMN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_RMM     = 'x' ! Matrix RMM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_GMN     = 'x' ! Matrix GMN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_GMNt    = 'x' ! Matrix GMNt    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_HMN     = 'x' ! Matrix HMN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_LMN     = 'x' ! Matrix KMN     sparse storage method (SYM or NONSYM)

      INTEGER(LONG), ALLOCATABLE      :: I_KNN(:)          ! Row indicators for nonzeros in stiff matrix KNN
      INTEGER(LONG), ALLOCATABLE      :: I_KNM(:)          ! Row indicators for nonzeros in stiff matrix KNM
      INTEGER(LONG), ALLOCATABLE      :: I_KMM(:)          ! Row indicators for nonzeros in stiff matrix KMM
      INTEGER(LONG), ALLOCATABLE      :: I_KMN(:)          ! Row indicators for nonzeros in stiff matrix KMN

      INTEGER(LONG), ALLOCATABLE      :: I_KNND(:)         ! Row indicators for nonzeros in stiff matrix KNND-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KNMD(:)         ! Row indicators for nonzeros in stiff matrix KNMD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KMMD(:)         ! Row indicators for nonzeros in stiff matrix KMMD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KMND(:)         ! Row indicators for nonzeros in stiff matrix KMND-diff stiff

      INTEGER(LONG), ALLOCATABLE      :: I_MNN(:)          ! Row indicators for nonzeros in mass matrix MNN
      INTEGER(LONG), ALLOCATABLE      :: I_MNM(:)          ! Row indicators for nonzeros in mass matrix MNM
      INTEGER(LONG), ALLOCATABLE      :: I_MMN(:)          ! Row indicators for nonzeros in mass matrix MMN
      INTEGER(LONG), ALLOCATABLE      :: I_MMM(:)          ! Row indicators for nonzeros in mass matrix MMM

      INTEGER(LONG), ALLOCATABLE      :: I_PN(:)           ! Row indicators for nonzeros in load matrix PN
      INTEGER(LONG), ALLOCATABLE      :: I_PM(:)           ! Row indicators for nonzeros in load matrix PM

      INTEGER(LONG), ALLOCATABLE      :: I_RMN(:)          ! Row indicators for nonzeros in constraint matrix RMN
      INTEGER(LONG), ALLOCATABLE      :: I_RMM(:)          ! Row indicators for nonzeros in constraint matrix RMM

      INTEGER(LONG), ALLOCATABLE      :: I_GMN(:)          ! Row indicators for nonzeros in matrix GMN = RMM(-1)*RMN
      INTEGER(LONG), ALLOCATABLE      :: I_GMNt(:)         ! Row indicators for nonzeros in GMNt (transp of GMN)

      INTEGER(LONG), ALLOCATABLE      :: I_HMN(:)          ! Row indicators for nonzeros in HMN MPC force rec matrix
      INTEGER(LONG), ALLOCATABLE      :: I_LMN(:)          ! Row indicators for nonzeros in LMN MPC force rec matrix

      INTEGER(LONG), ALLOCATABLE      :: I2_GMN(:)         ! Array of row numbers for all nonzero terms in matrix GMN

      INTEGER(LONG), ALLOCATABLE      :: J_KNN(:)          ! Col numbers for nonzeros in stiff matrix KNN
      INTEGER(LONG), ALLOCATABLE      :: J_KNM(:)          ! Col numbers for nonzeros in stiff matrix KNM
      INTEGER(LONG), ALLOCATABLE      :: J_KMM(:)          ! Col numbers for nonzeros in stiff matrix KMM
      INTEGER(LONG), ALLOCATABLE      :: J_KMN(:)          ! Col numbers for nonzeros in stiff matrix KMN

      INTEGER(LONG), ALLOCATABLE      :: J_KNND(:)         ! Col numbers for nonzeros in stiff matrix KNND-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KNMD(:)         ! Col numbers for nonzeros in stiff matrix KNMD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KMMD(:)         ! Col numbers for nonzeros in stiff matrix KMMD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KMND(:)         ! Col numbers for nonzeros in stiff matrix KMND-diff stiff

      INTEGER(LONG), ALLOCATABLE      :: J_MNN(:)          ! Col numbers for nonzeros in mass matrix MNN
      INTEGER(LONG), ALLOCATABLE      :: J_MNM(:)          ! Col numbers for nonzeros in mass matrix MNM
      INTEGER(LONG), ALLOCATABLE      :: J_MMN(:)          ! Col numbers for nonzeros in mass matrix MMN
      INTEGER(LONG), ALLOCATABLE      :: J_MMM(:)          ! Col numbers for nonzeros in mass matrix MMM

      INTEGER(LONG), ALLOCATABLE      :: J_PN(:)           ! Col numbers for nonzeros in load matrix PN
      INTEGER(LONG), ALLOCATABLE      :: J_PM(:)           ! Col numbers for nonzeros in load matrix PM

      INTEGER(LONG), ALLOCATABLE      :: J_RMN(:)          ! Col numbers for nonzeros in constraint matrix RMN
      INTEGER(LONG), ALLOCATABLE      :: J_RMM(:)          ! Col numbers for nonzeros in constraint matrix RMM

      INTEGER(LONG), ALLOCATABLE      :: J_GMN(:)          ! Col numbers for nonzeros in matrix GMN = RMM(-1)*RMN
      INTEGER(LONG), ALLOCATABLE      :: J_GMNt(:)         ! Col numbers for nonzeros in GMNt (transp of GMN)

      INTEGER(LONG), ALLOCATABLE      :: J_HMN(:)          ! Col numbers for nonzeros in HMN MPC force rec matrix
      INTEGER(LONG), ALLOCATABLE      :: J_LMN(:)          ! Col numbers for nonzeros in LMN MPC force rec matrix

      REAL(DOUBLE) , ALLOCATABLE      :: KNN(:)            ! Initially, N-set partition of KGG. Finally, N-set red stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KNM(:)            ! N-set row and M-set col part. of KGG
      REAL(DOUBLE) , ALLOCATABLE      :: KMM(:)            ! M-set row and M-set col part. of KGG
      REAL(DOUBLE) , ALLOCATABLE      :: KMN(:)            ! M-set row and N-set col part. of KGG

      REAL(DOUBLE) , ALLOCATABLE      :: KNND(:)           ! Initially, N-set partition of KGG. Finally, N-set red stiff-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KNMD(:)           ! N-set row and M-set col part. of KGGD-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KMMD(:)           ! M-set row and M-set col part. of KGGD-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KMND(:)           ! M-set row and N-set col part. of KGGDS-diff stiff

      REAL(DOUBLE) , ALLOCATABLE      :: MNN(:)            ! Initially, N-set part. of MGG. Finally, N-set reduced mass
      REAL(DOUBLE) , ALLOCATABLE      :: MNM(:)            ! N-set row and M-set col part. of MGG
      REAL(DOUBLE) , ALLOCATABLE      :: MMN(:)            ! M-set row and N-set col part. of MGG
      REAL(DOUBLE) , ALLOCATABLE      :: MMM(:)            ! M-set row and M-set col part. of MGG

      REAL(DOUBLE) , ALLOCATABLE      :: PN(:)             ! Initially, N-set part. of PG. Finally, N-set reduced load
      REAL(DOUBLE) , ALLOCATABLE      :: PM(:)             ! M-set row part. of PG

      REAL(DOUBLE) , ALLOCATABLE      :: RMN(:)            ! N-set cols of RMG
      REAL(DOUBLE) , ALLOCATABLE      :: RMM(:)            ! M-set cols of RMG

      REAL(DOUBLE) , ALLOCATABLE      :: GMN(:)            ! Nonzero terms in matrix GMN = -RMM(-1) x RMN
      REAL(DOUBLE) , ALLOCATABLE      :: GMNt(:)           ! Nonzero terms in matrix GMNt

      REAL(DOUBLE) , ALLOCATABLE      :: HMN(:)            ! Nonzero terms in matrix HMN
      REAL(DOUBLE) , ALLOCATABLE      :: LMN(:)            ! Nonzero terms in matrix LMN

! **********************************************************************************************************************************
! F, S-set arrays

      CHARACTER(  1*BYTE)             :: SYM_KFF     = 'x' ! Matrix KFF     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSF     = 'x' ! Matrix KSF     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KFS     = 'x' ! Matrix KFS     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSS     = 'x' ! Matrix KSS     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KFSe    = 'x' ! Matrix KFSe    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSSe    = 'x' ! Matrix KSSe    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KFFD    = 'x' ! Matrix KFFD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSFD    = 'x' ! Matrix KSFD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KFSD    = 'x' ! Matrix KFSD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSSD    = 'x' ! Matrix KSSD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KFSDe   = 'x' ! Matrix KFSDe   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KSSDe   = 'x' ! Matrix KSSDe   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MFF     = 'x' ! Matrix MFF     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MFS     = 'x' ! Matrix MFS     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MSF     = 'x' ! Matrix MFS     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MSS     = 'x' ! Matrix MSS     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PF      = 'x' ! Matrix PF      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PF_TMP  = 'x' ! Matrix PF_TMP  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PFYS    = 'x' ! Matrix PFYS    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PFYS1   = 'x' ! Matrix PFYS1   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PS      = 'x' ! Matrix PS      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_QSYS    = 'x' ! Matrix QSYS    sparse storage method (SYM or NONSYM)

      INTEGER(LONG), ALLOCATABLE      :: I_KFF(:)          ! Row indicators for nonzero terms in stiff matrix KNN
      INTEGER(LONG), ALLOCATABLE      :: I_KFS(:)          ! Row indicators for nonzero terms in stiff matrix KFS
      INTEGER(LONG), ALLOCATABLE      :: I_KSF(:)          ! Row indicators for nonzero terms in stiff matrix KSF
      INTEGER(LONG), ALLOCATABLE      :: I_KSS(:)          ! Row indicators for nonzero terms in stiff matrix KSS
      INTEGER(LONG), ALLOCATABLE      :: I_KSSe(:)         ! Row indicators for nonzero terms in stiff matrix KSSe
      INTEGER(LONG), ALLOCATABLE      :: I_KFSe(:)         ! Row indicators for nonzero terms in stiff matrix KFSe

      INTEGER(LONG), ALLOCATABLE      :: I_KFFD(:)         ! Row indicators for nonzero terms in stiff matrix KNND
      INTEGER(LONG), ALLOCATABLE      :: I_KFSD(:)         ! Row indicators for nonzero terms in stiff matrix KFSD
      INTEGER(LONG), ALLOCATABLE      :: I_KSFD(:)         ! Row indicators for nonzero terms in stiff matrix KSFD
      INTEGER(LONG), ALLOCATABLE      :: I_KSSD(:)         ! Row indicators for nonzero terms in stiff matrix KSS
      INTEGER(LONG), ALLOCATABLE      :: I_KSSDe(:)        ! Row indicators for nonzero terms in stiff matrix KSSe
      INTEGER(LONG), ALLOCATABLE      :: I_KFSDe(:)        ! Row indicators for nonzero terms in stiff matrix KFSe

      INTEGER(LONG), ALLOCATABLE      :: I_MFF(:)          ! Row indicators for nonzero terms in mass matrix MNN
      INTEGER(LONG), ALLOCATABLE      :: I_MFS(:)          ! Row indicators for nonzero terms in mass matrix MFS
      INTEGER(LONG), ALLOCATABLE      :: I_MSF(:)          ! Row indicators for nonzero terms in mass matrix MSF
      INTEGER(LONG), ALLOCATABLE      :: I_MSS(:)          ! Row indicators for nonzero terms in mass matrix MSS

      INTEGER(LONG), ALLOCATABLE      :: I_PF(:)           ! Row indicators for nonzero terms in load matrix PF
      INTEGER(LONG), ALLOCATABLE      :: I_PF_TMP(:)       ! Row indicators for nonzero terms in load matrix PF_TMP
      INTEGER(LONG), ALLOCATABLE      :: I_PS(:)           ! Row indicators for nonzero terms in load matrix PS
      INTEGER(LONG), ALLOCATABLE      :: I_PFYS(:)         ! Row indicators for nonzero terms in load matrix PFYS
      INTEGER(LONG), ALLOCATABLE      :: I_PFYS1(:)        ! Row indicators for nonzero terms in load matrix PFYS1
      INTEGER(LONG), ALLOCATABLE      :: I_QSYS(:)         ! Row indicators for nonzeros in constr load matrix QSYS

      INTEGER(LONG), ALLOCATABLE      :: J_KFF(:)          ! Col numbers for nonzero terms in stiff matrix KNN
      INTEGER(LONG), ALLOCATABLE      :: J_KFS(:)          ! Col numbers for nonzero terms in stiff matrix KFS
      INTEGER(LONG), ALLOCATABLE      :: J_KSF(:)          ! Col numbers for nonzero terms in stiff matrix KSF
      INTEGER(LONG), ALLOCATABLE      :: J_KSS(:)          ! Col numbers for nonzero terms in stiff matrix KSS
      INTEGER(LONG), ALLOCATABLE      :: J_KSSe(:)         ! Col numbers for nonzero terms in stiff matrix KSSe
      INTEGER(LONG), ALLOCATABLE      :: J_KFSe(:)         ! Col numbers for nonzero terms in stiff matrix KFSe

      INTEGER(LONG), ALLOCATABLE      :: J_KFFD(:)         ! Col numbers for nonzero terms in stiff matrix KNND-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KFSD(:)         ! Col numbers for nonzero terms in stiff matrix KFSD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KSFD(:)         ! Col numbers for nonzero terms in stiff matrix KSFD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KSSD(:)         ! Col numbers for nonzero terms in stiff matrix KSSD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KSSDe(:)        ! Col numbers for nonzero terms in stiff matrix KSSDe-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KFSDe(:)        ! Col numbers for nonzero terms in stiff matrix KFSDe-diff stiff

      INTEGER(LONG), ALLOCATABLE      :: J_MFF(:)          ! Col numbers for nonzero terms in mass matrix MNN
      INTEGER(LONG), ALLOCATABLE      :: J_MFS(:)          ! Col numbers for nonzero terms in mass matrix MFS
      INTEGER(LONG), ALLOCATABLE      :: J_MSF(:)          ! Col numbers for nonzero terms in mass matrix MSF
      INTEGER(LONG), ALLOCATABLE      :: J_MSS(:)          ! Col numbers for nonzero terms in mass matrix MSS

      INTEGER(LONG), ALLOCATABLE      :: J_PF(:)           ! Col numbers for nonzero terms in load matrix PF
      INTEGER(LONG), ALLOCATABLE      :: J_PF_TMP(:)       ! Col numbers for nonzero terms in load matrix PF_TMP
      INTEGER(LONG), ALLOCATABLE      :: J_PS(:)           ! Col numbers for nonzero terms in load matrix PS
      INTEGER(LONG), ALLOCATABLE      :: J_PFYS(:)         ! Col numbers for nonzero terms in load matrix PFYS
      INTEGER(LONG), ALLOCATABLE      :: J_PFYS1(:)        ! Col numbers for nonzero terms in load matrix PFYS1
      INTEGER(LONG), ALLOCATABLE      :: J_QSYS(:)         ! Col numbers for nonzeros in constr load matrix QSYS

      REAL(DOUBLE) , ALLOCATABLE      :: KFF(:)            ! Initially, F-set part. of KNN . Finally, F-set reduced stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KFS(:)            ! F-set row and S-set col part. of KNN
      REAL(DOUBLE) , ALLOCATABLE      :: KSF(:)            ! S-set row and F-set col part. of KNN
      REAL(DOUBLE) , ALLOCATABLE      :: KSS(:)            ! S-set row and S-set col part. of KNN
      REAL(DOUBLE) , ALLOCATABLE      :: KFSe(:)           ! Cols of KFS for SE-set (enforced displ DOF's)
      REAL(DOUBLE) , ALLOCATABLE      :: KSSe(:)           ! Cols of KSS for SE-set (enforced displ DOF's)

      REAL(DOUBLE) , ALLOCATABLE      :: KFFD(:)           ! Initially, F-set part. of KNND. Finally, F-set reduced stiff-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KFSD(:)           ! F-set row and S-set col part. of KNND
      REAL(DOUBLE) , ALLOCATABLE      :: KSFD(:)           ! S-set row and F-set col part. of KNND
      REAL(DOUBLE) , ALLOCATABLE      :: KSSD(:)           ! S-set row and S-set col part. of KNND
      REAL(DOUBLE) , ALLOCATABLE      :: KFSDe(:)          ! Cols of KFSDe for SE-set (enforced displ DOF's)
      REAL(DOUBLE) , ALLOCATABLE      :: KSSDe(:)          ! Cols of KSSDe for SE-set (enforced displ DOF's)

      REAL(DOUBLE) , ALLOCATABLE      :: MFF(:)            ! Initially, F-set part. of MNN. Finally, F-set reduced mass
      REAL(DOUBLE) , ALLOCATABLE      :: MFS(:)            ! F-set row and S-set col part. of MNN
      REAL(DOUBLE) , ALLOCATABLE      :: MSF(:)            ! F-set row and S-set col part. of MNN
      REAL(DOUBLE) , ALLOCATABLE      :: MSS(:)            ! S-set row and S-set col part. of MNN

      REAL(DOUBLE) , ALLOCATABLE      :: PF(:)             ! Initially, F-set part. of PN. Finally, F-set reduced load
      REAL(DOUBLE) , ALLOCATABLE      :: PF_TMP(:)         ! Temp array for PF while it is being formed
      REAL(DOUBLE) , ALLOCATABLE      :: PS(:)             ! S-set row part. of PN
      REAL(DOUBLE) , ALLOCATABLE      :: PFYS1(:)          ! = -KFS x YSe
      REAL(DOUBLE) , ALLOCATABLE      :: PFYS(:)           ! NSUB columns, each = -KFSe*YSe
      REAL(DOUBLE) , ALLOCATABLE      :: QSYS(:)           ! Nonzero terms in matrix QSYS = KSSe x YSe

! **********************************************************************************************************************************
! A, O-set arrays

      CHARACTER(  1*BYTE)             :: SYM_KAA     = 'x' ! Matrix KAA     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KAO     = 'x' ! Matrix KAO     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KOO     = 'x' ! Matrix KOO     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KAAD    = 'x' ! Matrix KAAD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KAOD    = 'x' ! Matrix KAOD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KOOD    = 'x' ! Matrix KOOD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MAA     = 'x' ! Matrix MAA     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MAO     = 'x' ! Matrix MAO     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MOO     = 'x' ! Matrix MOO     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PA      = 'x' ! Matrix PA      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PO      = 'x' ! Matrix PO      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_GOA     = 'x' ! Matrix GOA     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_GOAt    = 'x' ! Matrix GOAt    sparse storage method (SYM or NONSYM)

      INTEGER(LONG), ALLOCATABLE      :: I_KAA(:)          ! Row indicators for nonzero terms in stiff matrix KAA
      INTEGER(LONG), ALLOCATABLE      :: I_KAO(:)          ! Row indicators for nonzero terms in stiff matrix KAO
      INTEGER(LONG), ALLOCATABLE      :: I_KOO(:)          ! Row indicators for nonzero terms in stiff matrix KOO
      INTEGER(LONG), ALLOCATABLE      :: I2_KOO(:)         ! Row values     for nonzero terms in stiff matrix KOO
      INTEGER(LONG), ALLOCATABLE      :: I_KOOs(:)         ! Row indicators for nonzero terms in stiff matrix KOOs
      INTEGER(LONG), ALLOCATABLE      :: I2_KOOs(:)        ! Row values     for nonzero terms in stiff matrix KOOs

      INTEGER(LONG), ALLOCATABLE      :: I_KAAD(:)         ! Row indicators for nonzero terms in stiff matrix KAAD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KAOD(:)         ! Row indicators for nonzero terms in stiff matrix KAOD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KOOD(:)         ! Row indicators for nonzero terms in stiff matrix KOOD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I2_KOOD(:)        ! Row values     for nonzero terms in stiff matrix KOOD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I_KOODs(:)        ! Row indicators for nonzero terms in stiff matrix KOODs-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: I2_KOODs(:)       ! Row values     for nonzero terms in stiff matrix KOODs-diff stiff

      INTEGER(LONG), ALLOCATABLE      :: I_MAA(:)          ! Row indicators for nonzeros in mass matrix MAA
      INTEGER(LONG), ALLOCATABLE      :: I_MAO(:)          ! Row indicators for nonzero terms in mass matrix MAO
      INTEGER(LONG), ALLOCATABLE      :: I_MOO(:)          ! Row indicators for nonzero terms in mass matrix MOO

      INTEGER(LONG), ALLOCATABLE      :: I_PA(:)           ! Row indicators for nonzero terms in load matrix PA
      INTEGER(LONG), ALLOCATABLE      :: I_PO(:)           ! Row indicators for nonzero terms in load matrix PO

      INTEGER(LONG), ALLOCATABLE      :: I_GOA(:)          ! Row indicators for nonzeros in O-set reduction matrix GOA
      INTEGER(LONG), ALLOCATABLE      :: I2_GOA(:)         ! Array of row numbers for all nonzero terms in matrix GOA
      INTEGER(LONG), ALLOCATABLE      :: I_GOAt(:)         ! Row indicators for nonzeros in matrix GOAt (transp of GOA)

      INTEGER(LONG), ALLOCATABLE      :: J_KAA(:)          ! Col numbers for nonzero terms in stiff matrix KAA
      INTEGER(LONG), ALLOCATABLE      :: J_KAO(:)          ! Col numbers for nonzero terms in stiff matrix KAO
      INTEGER(LONG), ALLOCATABLE      :: J_KOO(:)          ! Col numbers for nonzero terms in stiff matrix KOO
      INTEGER(LONG), ALLOCATABLE      :: J_KOOs(:)         ! Col numbers for nonzero terms in stiff matrix KOOs

      INTEGER(LONG), ALLOCATABLE      :: J_KAAD(:)         ! Col numbers for nonzero terms in stiff matrix KAAD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KAOD(:)         ! Col numbers for nonzero terms in stiff matrix KAOD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KOOD(:)         ! Col numbers for nonzero terms in stiff matrix KOOD-diff stiff
      INTEGER(LONG), ALLOCATABLE      :: J_KOODs(:)        ! Col numbers for nonzero terms in stiff matrix KOODs-diff stiff

      INTEGER(LONG), ALLOCATABLE      :: J_MAA(:)          ! Col numbers for nonzeros in mass matrix MAA
      INTEGER(LONG), ALLOCATABLE      :: J_MAO(:)          ! Col numbers for nonzero terms in mass matrix MAO
      INTEGER(LONG), ALLOCATABLE      :: J_MOO(:)          ! Col numbers for nonzero terms in mass matrix MOO

      INTEGER(LONG), ALLOCATABLE      :: J_PA(:)           ! Col numbers for nonzero terms in load matrix PA
      INTEGER(LONG), ALLOCATABLE      :: J_PO(:)           ! Col numbers for nonzero terms in load matrix PO

      INTEGER(LONG), ALLOCATABLE      :: J_GOA(:)          ! Col numbers for nonzeros in O-set reduction matrix GOA
      INTEGER(LONG), ALLOCATABLE      :: J_GOAt(:)         ! Col numbers for nonzeros in matrix GOAt (transp of GOA)

      REAL(DOUBLE) , ALLOCATABLE      :: KAA(:)            ! Init, A-set part. of KFF. Finally, A-set red stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KAO(:)            ! A-set row and O-set col part. of KFF
      REAL(DOUBLE) , ALLOCATABLE      :: KOO(:)            ! O-set row and O-set col part. of KFF
      REAL(DOUBLE) , ALLOCATABLE      :: KOOs(:)           ! O-set row and O-set col part. of KFF

      REAL(DOUBLE) , ALLOCATABLE      :: KAAD(:)           ! Init, A-set part. of KFF. Finally, A-set red stiff terms-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KAOD(:)           ! A-set row and O-set col part. of KFFD-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KOOD(:)           ! O-set row and O-set col part. of KFFD-diff stiff
      REAL(DOUBLE) , ALLOCATABLE      :: KOODs(:)          ! O-set row and O-set col part. of KFFD-diff stiff

      REAL(DOUBLE) , ALLOCATABLE      :: MAA(:)            ! Initially, A-set part. of MFF. Finally, A-set reduced mass terms
!                                                            stored symmetrically
      REAL(DOUBLE) , ALLOCATABLE      :: MAO(:)            ! A-set row and O-set col part. of MFF
      REAL(DOUBLE) , ALLOCATABLE      :: MOO(:)            ! O-set row and O-set col part. of MFF

      REAL(DOUBLE) , ALLOCATABLE      :: PA(:)             ! Initially, A-set part. of PF. Finally, A-set reduced load terms
      REAL(DOUBLE) , ALLOCATABLE      :: PO(:)             ! O-set row part. of PF

      REAL(DOUBLE) , ALLOCATABLE      :: GOA(:)            ! Nonzero terms in matrix GOA = = -KOO(-1) x KAOt
      REAL(DOUBLE) , ALLOCATABLE      :: GOAt(:)           ! Nonzero terms in matrix GOAt

! **********************************************************************************************************************************
! L, R-set arrays

      CHARACTER(  1*BYTE)             :: SYM_CG_LTM  = 'x' ! Matrix CG_LTM  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_DLR     = 'x' ! Matrix DLR     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_DLRt    = 'x' ! Matrix DLRt    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIZL   = 'x' ! Matrix PHIZL   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIZL1  = 'x' ! Matrix PHIZL1  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIZL1t = 'x' ! Matrix PHIZL2  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIZL2  = 'x' ! Matrix PHIZL2  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_IF_LTM  = 'x' ! Matrix IF_LTM  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_IRR     = 'x' ! Matrix IRR     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KLL     = 'x' ! Matrix KLL     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KLLs    = 'x' ! Matrix KLLs    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRL     = 'x' ! Matrix KRL     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRR     = 'x' ! Matrix KRRD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KLLD    = 'x' ! Matrix KLLD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KLLDs   = 'x' ! Matrix KLLDs   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRLD    = 'x' ! Matrix KRLD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRRD    = 'x' ! Matrix KRRD    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRRcbn  = 'x' ! Matrix KRRcb   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRRcbs  = 'x' ! Matrix KRRcb   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KRRcb   = 'x' ! Matrix KRRcb   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KXX     = 'x' ! Matrix KXX     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMSM    = 'x' ! Matrix KMSM    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_LTM     = 'x' ! Matrix LTM     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_KMSMn   = 'x' ! Matrix MMSMn   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MPF0    = 'x' ! Matrix MPF0    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MLL     = 'x' ! Matrix MLL     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MLLn    = 'x' ! Matrix MLLn    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MLLs    = 'x' ! Matrix MLLs    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MLR     = 'x' ! Matrix MLR     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MRL     = 'x' ! Matrix MRL     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MRR     = 'x' ! Matrix MRR     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MRN     = 'x' ! Matrix MRN     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MRRcb   = 'x' ! Matrix MRRcb   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MRRcbn  = 'x' ! Matrix MRRcbn  sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MXX     = 'x' ! Matrix MXX     sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_MXXn    = 'x' ! Matrix MXXn    sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIXA   = 'x' ! Matrix PHIXA   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PHIXG   = 'x' ! Matrix PHIXG   sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PL      = 'x' ! Matrix PL      sparse storage method (SYM or NONSYM)
      CHARACTER(  1*BYTE)             :: SYM_PR      = 'x' ! Matrix PR      sparse storage method (SYM or NONSYM)

      INTEGER(LONG), ALLOCATABLE      :: I_CG_LTM  (:)     ! Row indicators for nonzeros in matrix CG_LTM

      INTEGER(LONG), ALLOCATABLE      :: I_DLR(:)          ! Row indicators for nonzeros in matrix DLR
      INTEGER(LONG), ALLOCATABLE      :: I_DLRt(:)         ! Row indicators for nonzeros in matrix DLRt
      INTEGER(LONG), ALLOCATABLE      :: I2_DLR(:)         ! Array of row numbers for all nonzero terms in matrix DLR
      INTEGER(LONG), ALLOCATABLE      :: I2_DLRt(:)        ! Array of row numbers for all nonzero terms in matrix DLRt

      INTEGER(LONG), ALLOCATABLE      :: I_IF_LTM(:)       ! Row indicators for nonzero terms in matrix IF_LTM
      INTEGER(LONG), ALLOCATABLE      :: I_IRR(:)          ! Row indicators for nonzero terms in identity matrix IRR

      INTEGER(LONG), ALLOCATABLE      :: I_KLL(:)          ! Row indicators for nonzero terms in stiff matrix KLL
      INTEGER(LONG), ALLOCATABLE      :: I2_KLL(:)         ! Row numbers for nonzero terms in stiff matrix KLL
      INTEGER(LONG), ALLOCATABLE      :: I_KLLs(:)         ! Row indicators for nonzero terms in stiff matrix KLLs
      INTEGER(LONG), ALLOCATABLE      :: I2_KLLs(:)        ! Row numbers for nonzero terms in stiff matrix KLL
      INTEGER(LONG), ALLOCATABLE      :: I_KRL(:)          ! Row indicators for nonzero terms in stiff matrix KRL
      INTEGER(LONG), ALLOCATABLE      :: I_KRR(:)          ! Row indicators for nonzero terms in stiff matrix KRR
      INTEGER(LONG), ALLOCATABLE      :: I_KRRcb(:)        ! Row indicators for nonzero terms in stiff matrix KRRcb
      INTEGER(LONG), ALLOCATABLE      :: I_KRRcbn(:)       ! Row indicators for nonzero terms in stiff matrix KRRcbn
      INTEGER(LONG), ALLOCATABLE      :: I_KRRcbs(:)       ! Row indicators for nonzero terms in stiff matrix KRRcbs
      INTEGER(LONG), ALLOCATABLE      :: I_KXX(:)          ! Row indicators for nonzero terms in stiff matrix KXX

      INTEGER(LONG), ALLOCATABLE      :: I_KLLD(:)         ! Row indicators for nonzero terms in stiff matrix KLLD
      INTEGER(LONG), ALLOCATABLE      :: I_KLLDn(:)        ! Row indicators for nonzero terms in stiff matrix KLLDn
      INTEGER(LONG), ALLOCATABLE      :: I2_KLLD(:)        ! Row numbers for nonzero terms in stiff matrix KLLD
      INTEGER(LONG), ALLOCATABLE      :: I_KLLDs(:)        ! Row indicators for nonzero terms in stiff matrix KLLDs
      INTEGER(LONG), ALLOCATABLE      :: I2_KLLDs(:)       ! Row numbers for nonzero terms in stiff matrix KLLDs
      INTEGER(LONG), ALLOCATABLE      :: I_KRLD(:)         ! Row indicators for nonzero terms in stiff matrix KRLD
      INTEGER(LONG), ALLOCATABLE      :: I_KRRD(:)         ! Row indicators for nonzero terms in stiff matrix KRRD

      INTEGER(LONG), ALLOCATABLE      :: I_KMSM(:)         ! Row ind for nonzero terms in KLL - sigma*MLL
      INTEGER(LONG), ALLOCATABLE      :: I2_KMSM(:)        ! Row values for nonzero terms in KLL - sigma*MLL
      INTEGER(LONG), ALLOCATABLE      :: I_KMSMn(:)        ! Row ind for nonzero terms in KLL - sigma*MLL (nonsym)
      INTEGER(LONG), ALLOCATABLE      :: I_KMSMs(:)        ! Row ind for nonzero terms in KLL - sigma*MLL  (sym)
      INTEGER(LONG), ALLOCATABLE      :: I2_KMSMs(:)       ! Row values for nonzero terms in KLL - sigma*MLL (sym)

      INTEGER(LONG), ALLOCATABLE      :: I_MPF0 (:)        ! Row indicators for nonzeros in temp matrix MPF0
      INTEGER(LONG), ALLOCATABLE      :: I_MLL(:)          ! Row indicators for nonzeros in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: I2_MLL(:)         ! Row numbers for nonzero terms in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: I_MLLn(:)         ! Row indicators for nonzeros in mass matrix MLLn
      INTEGER(LONG), ALLOCATABLE      :: I_MLLs(:)         ! Row indicators for nonzeros in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: I2_MLLs(:)        ! Row numbers for nonzero terms in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: I_MLR(:)          ! Row indicators for nonzero terms in mass matrix MLR
      INTEGER(LONG), ALLOCATABLE      :: I_MRL(:)          ! Row indicators for nonzero terms in mass matrix MRL
      INTEGER(LONG), ALLOCATABLE      :: I_MRR(:)          ! Row indicators for nonzero terms in mass matrix MRR
      INTEGER(LONG), ALLOCATABLE      :: I_MRN(:)          ! Row indicators for nonzero terms in mass matrix MRN
      INTEGER(LONG), ALLOCATABLE      :: I_MRRcb(:)        ! Row indicators for nonzero terms in mass matrix MRRcb
      INTEGER(LONG), ALLOCATABLE      :: I_MRRcbn(:)       ! Row indicators for nonzero terms in mass matrix MRRcb
      INTEGER(LONG), ALLOCATABLE      :: I_MXX(:)          ! Row indicators for nonzero terms in mass matrix MRRcb
      INTEGER(LONG), ALLOCATABLE      :: I_MXXn(:)         ! Row indicators for nonzero terms in mass matrix MRRcbn

      INTEGER(LONG), ALLOCATABLE      :: I_LTM(:)          ! Row indicators for nonzero terms in matrix LTM

      INTEGER(LONG), ALLOCATABLE      :: I_PHIXA(:)        ! Row indicators for nonzero terms in CB matrix PHIXA
      INTEGER(LONG), ALLOCATABLE      :: I_PHIXG(:)        ! Row indicators for nonzero terms in CB matrix PHIXG

      INTEGER(LONG), ALLOCATABLE      :: I_PHIZL(:)        ! Row indicators for nonzeros in matrix PHIZL
      INTEGER(LONG), ALLOCATABLE      :: I_PHIZL1(:)       ! Row indicators for nonzeros in matrix PHIZL1
      INTEGER(LONG), ALLOCATABLE      :: I_PHIZL1t(:)      ! Row indicators for nonzeros in matrix PHIZL1
      INTEGER(LONG), ALLOCATABLE      :: I_PHIZL2(:)       ! Row indicators for nonzeros in matrix PHIZL2
      INTEGER(LONG), ALLOCATABLE      :: I2_PHIZL1(:)      ! Array of row numbers for all nonzero terms in matrix PHIZL1
      INTEGER(LONG), ALLOCATABLE      :: I2_PHIZL1t(:)     ! Array of row numbers for all nonzero terms in matrix PHIZL1t

      INTEGER(LONG), ALLOCATABLE      :: I_PL(:)           ! Row indicators for nonzero terms in load matrix PL
      INTEGER(LONG), ALLOCATABLE      :: I_PR(:)           ! Row indicators for nonzero terms in load matrix PR

      INTEGER(LONG), ALLOCATABLE      :: J_CG_LTM(:)       ! Col numbers for nonzeros in matrix CG_LTM

      INTEGER(LONG), ALLOCATABLE      :: J_DLR(:)          ! Col numbers for nonzeros in matrix DLR
      INTEGER(LONG), ALLOCATABLE      :: J_DLRt(:)         ! Col numbers for nonzeros in matrix DLRt

      INTEGER(LONG), ALLOCATABLE      :: J_IF_LTM(:)       ! Col numbers for nonzero term in matrix IF_LTM
      INTEGER(LONG), ALLOCATABLE      :: J_IRR(:)          ! Col numbers for nonzero term in identity matrix IRR

      INTEGER(LONG), ALLOCATABLE      :: J_KLL(:)          ! Col numbers for nonzero terms in stiff matrix KLL
      INTEGER(LONG), ALLOCATABLE      :: J_KLLs(:)         ! Col numbers for nonzero terms in stiff matrix KLL
      INTEGER(LONG), ALLOCATABLE      :: J_KRL(:)          ! Col numbers for nonzero terms in stiff matrix KRL
      INTEGER(LONG), ALLOCATABLE      :: J_KRR(:)          ! Col numbers for nonzero terms in stiff matrix KRR
      INTEGER(LONG), ALLOCATABLE      :: J_KRRcb(:)        ! Col numbers for nonzero terms in stiff matrix KRRcb
      INTEGER(LONG), ALLOCATABLE      :: J_KRRcbn(:)       ! Col numbers for nonzero terms in stiff matrix KRRcbn
      INTEGER(LONG), ALLOCATABLE      :: J_KRRcbs(:)       ! Col numbers for nonzero terms in stiff matrix KRRcbs
      INTEGER(LONG), ALLOCATABLE      :: J_KXX  (:)        ! Col numbers for nonzero terms in stiff matrix KXX

      INTEGER(LONG), ALLOCATABLE      :: J_KLLD(:)         ! Col numbers for nonzero terms in stiff matrix KLLD
      INTEGER(LONG), ALLOCATABLE      :: J_KLLDn(:)        ! Col numbers for nonzero terms in stiff matrix KLLDn
      INTEGER(LONG), ALLOCATABLE      :: J_KLLDs(:)        ! Col numbers for nonzero terms in stiff matrix KLLD
      INTEGER(LONG), ALLOCATABLE      :: J_KRLD(:)         ! Col numbers for nonzero terms in stiff matrix KRLD
      INTEGER(LONG), ALLOCATABLE      :: J_KRRD(:)         ! Col numbers for nonzero terms in stiff matrix KRRD

      INTEGER(LONG), ALLOCATABLE      :: J_KMSM(:)         ! Col ind for nonzero terms in KLL - sigma*MLL
      INTEGER(LONG), ALLOCATABLE      :: J_KMSMn(:)        ! Col ind for nonzero terms in KLL - sigma*MLL (nonsym)
      INTEGER(LONG), ALLOCATABLE      :: J_KMSMs(:)        ! Col ind for nonzero terms in KLL - sigma*MLL (sym)

      INTEGER(LONG), ALLOCATABLE      :: J_LTM(:)          ! Col numbers for nonzero term in matrix LTM

      INTEGER(LONG), ALLOCATABLE      :: J_MPF0 (:)        ! Col numbers for nonzeros in temp matrix MPF0
      INTEGER(LONG), ALLOCATABLE      :: J_MLL(:)          ! Col numbers for nonzeros in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: J_MLLn(:)         ! Col numbers for nonzeros in mass matrix MLLn
      INTEGER(LONG), ALLOCATABLE      :: J_MLLs(:)         ! Col numbers for nonzeros in mass matrix MLL
      INTEGER(LONG), ALLOCATABLE      :: J_MLR(:)          ! Col numbers for nonzero terms in mass matrix MLR
      INTEGER(LONG), ALLOCATABLE      :: J_MRL(:)          ! Col numbers for nonzero terms in mass matrix MRL
      INTEGER(LONG), ALLOCATABLE      :: J_MRR(:)          ! Col numbers for nonzero terms in mass matrix MRR
      INTEGER(LONG), ALLOCATABLE      :: J_MRN  (:)        ! Col numbers for nonzero terms in mass matrix MRN
      INTEGER(LONG), ALLOCATABLE      :: J_MRRcb(:)        ! Col numbers for nonzero terms in mass matrix MRRcb
      INTEGER(LONG), ALLOCATABLE      :: J_MRRcbn(:)       ! Col numbers for nonzero terms in mass matrix MRRcbn
      INTEGER(LONG), ALLOCATABLE      :: J_MXX(:)          ! Col numbers for nonzero terms in mass matrix MRRcb
      INTEGER(LONG), ALLOCATABLE      :: J_MXXn(:)         ! Col numbers for nonzero terms in mass matrix MRRcbn

      INTEGER(LONG), ALLOCATABLE      :: J_PHIXA(:)        ! Col numbers for nonzero terms in CB matrix PHIXA
      INTEGER(LONG), ALLOCATABLE      :: J_PHIXG(:)        ! Col numbers for nonzero terms in CB matrix PHIXG

      INTEGER(LONG), ALLOCATABLE      :: J_PHIZL(:)        ! Col numbers for nonzeros in matrix PHIZL
      INTEGER(LONG), ALLOCATABLE      :: J_PHIZL1(:)       ! Col numbers for nonzeros in matrix PHIZL1
      INTEGER(LONG), ALLOCATABLE      :: J_PHIZL1t(:)      ! Col numbers for nonzeros in matrix PHIZL1t
      INTEGER(LONG), ALLOCATABLE      :: J_PHIZL2(:)       ! Col numbers for nonzeros in matrix PHIZL2

      INTEGER(LONG), ALLOCATABLE      :: J_PL(:)           ! Col numbers for nonzero terms in load matrix PL
      INTEGER(LONG), ALLOCATABLE      :: J_PR(:)           ! Col numbers for nonzero terms in load matrix PR

      REAL(DOUBLE) , ALLOCATABLE      :: CG_LTM(:)         ! Real values for matrix CG_LTM

      REAL(DOUBLE) , ALLOCATABLE      :: DLR(:)            ! -KLL(-1)*KRL
      REAL(DOUBLE) , ALLOCATABLE      :: DLRt(:)           ! DLR'
      REAL(DOUBLE) , ALLOCATABLE      :: PHIZL(:)          ! CB Displ Transformation Matrix for the L-set

      REAL(DOUBLE) , ALLOCATABLE      :: IF_LTM(:)         ! Real terms in matrix IF_LTM
      REAL(DOUBLE) , ALLOCATABLE      :: IRR(:)            ! R-set row and R-set col identity matrix

      REAL(DOUBLE) , ALLOCATABLE      :: KLL(:)            ! Initially, L-set part. of KAA.  Finally, A-set reduced stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KLLs(:)           ! Initially, L-set part. of KAA. Finally, A-set reduced stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KRL(:)            ! L-set row and R-set col part. of KAA
      REAL(DOUBLE) , ALLOCATABLE      :: KRR(:)            ! R-set row and R-set col part. of KAA
      REAL(DOUBLE) , ALLOCATABLE      :: KRRcb(:)          ! R-set row and R-set col part. of CB stiffness matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KRRcbn(:)         ! R-set row and R-set col part. of CB stiffness matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KRRcbs(:)         ! R-set row and R-set col part. of CB stiffness matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KXX  (:)          ! R-set row and R-set col part. of CB stiffness matrix

      REAL(DOUBLE) , ALLOCATABLE      :: KLLD(:)           ! Initially, L-set part. of KAAD. Finally, A-set reduced stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KLLDn(:)          ! KLLD in stored nonsym format
      REAL(DOUBLE) , ALLOCATABLE      :: KLLDs(:)          ! Initially, L-set part. of KAAD. Finally, A-set reduced stiff terms
      REAL(DOUBLE) , ALLOCATABLE      :: KRLD(:)           ! L-set row and R-set col part. of KAAD
      REAL(DOUBLE) , ALLOCATABLE      :: KRRD(:)           ! R-set row and R-set col part. of KAAD

      REAL(DOUBLE) , ALLOCATABLE      :: KMSM(:)           ! KAA - sigma*MAA stored symmetrically
      REAL(DOUBLE) , ALLOCATABLE      :: KMSMn(:)          ! KAA - sigma*MAA stored nonsymmetrically
      REAL(DOUBLE) , ALLOCATABLE      :: KMSMs(:)          ! KAA - sigma*MAA stored symmetrically

      REAL(DOUBLE) , ALLOCATABLE      :: LTM(:)            ! Real terms in matrix LTM

      REAL(DOUBLE) , ALLOCATABLE      :: MPF0 (:)          ! Temp array for eff mass calc
      REAL(DOUBLE) , ALLOCATABLE      :: MLL(:)            ! Initially, L-set part. of MAA. Finally, L-set reduced mass terms
      REAL(DOUBLE) , ALLOCATABLE      :: MLLn(:)           ! Initially, L-set part. of MAA. Finally, L-set reduced mass terms
      REAL(DOUBLE) , ALLOCATABLE      :: MLLs(:)           ! Initially, L-set part. of MAA. Finally, L-set reduced mass terms
      REAL(DOUBLE) , ALLOCATABLE      :: MLR(:)            ! L-set row and R-set col part. of MLR
      REAL(DOUBLE) , ALLOCATABLE      :: MRL(:)            ! R-set row and L-set col part. of MRL
      REAL(DOUBLE) , ALLOCATABLE      :: MRR(:)            ! R-set row and R-set col part. of MRR
      REAL(DOUBLE) , ALLOCATABLE      :: MRN  (:)          ! R-set row and NVEC-set col part. of CB mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MRRcb(:)          ! R-set row and R-set col part. of CB mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MRRcbn(:)         ! R-set row and R-set col part. of CB mass matrix in nonsym format
      REAL(DOUBLE) , ALLOCATABLE      :: MXX  (:)          ! R-set row and R-set col part. of CB mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MXXn(:)           ! R-set row and R-set col part. of CB mass matrix

      REAL(DOUBLE) , ALLOCATABLE      :: PHIXA(:)          ! CB matrix PHIXA
      REAL(DOUBLE) , ALLOCATABLE      :: PHIXG(:)          ! CB matrix PHIXG

      REAL(DOUBLE) , ALLOCATABLE      :: PHIZL1(:)         ! -KLL(-1)*(MLR + MLL*DLR) part of PHIZL
      REAL(DOUBLE) , ALLOCATABLE      :: PHIZL1t(:)        ! PHIZL1 transposed
      REAL(DOUBLE) , ALLOCATABLE      :: PHIZL2(:)         ! -EIGEN_VEC*diag(EIGEN_VAL) part of PHIZL

      REAL(DOUBLE) , ALLOCATABLE      :: PL(:)             ! Initially, L-set part. of PA. Finally, L-set reduced load terms
      REAL(DOUBLE) , ALLOCATABLE      :: PR(:)             ! R-set row part. of PA

! * NOTE:
!   In the sparse (compressed row storage) format:

!     (1) Array I_KGG (for example) is dimensioned NDOFG+1 and contains integers such that: row I of the G-set stiffness matrix will
!         have I_KGG(I+1) - I_KGG(I) nonzero terms in sparse array KGG

!     (2) Array J_KGG is dimensioned NTERM_KGG (the number of nonzeros in the G-set stiff matrix). Term j gives the row number the
!         G-set stiffness matrix for the j-th term in sparse array KGG

      END MODULE SPARSE_MATRICES
