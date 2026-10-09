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

      MODULE FULL_MATRICES

! Arrays for full matrices. In the early days of MYSTRN, the capability for using full 2-D arrays for matrices (KGG stiffness, etc)
! was allowed. As MYSTRAN grew to employ sparse matrix technology, this full matrix capability was only used to make sure that the
! sparse algorithms coded were working properly. As of approx 2005, the full capability can be considered obsolete.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      REAL(DOUBLE) , ALLOCATABLE      :: KGG_FULL(:,:)     ! Full representation of KGG stiff matrix

      REAL(DOUBLE) , ALLOCATABLE      :: KNN_FULL(:,:)     ! Full representation of KNN stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KNM_FULL(:,:)     ! Full representation of KNM stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KMM_FULL(:,:)     ! Full representation of KMM stiff matrix

      REAL(DOUBLE) , ALLOCATABLE      :: MNN_FULL(:,:)     ! Full representation of MNN mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MNM_FULL(:,:)     ! Full representation of MNM mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MMM_FULL(:,:)     ! Full representation of MMM mass matrix

      REAL(DOUBLE) , ALLOCATABLE      :: PN_FULL(:,:)      ! Full representation of PN load vectors
      REAL(DOUBLE) , ALLOCATABLE      :: PM_FULL(:,:)      ! Full representation of PM load vectors

      REAL(DOUBLE) , ALLOCATABLE      :: KFF_FULL(:,:)     ! Full representation of KFF stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KFS_FULL(:,:)     ! Full representation of KFS stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KFSe_FULL(:,:)    ! Full representation of KFSe stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KSS_FULL(:,:)     ! Full representation of KSS stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KSSe_FULL(:,:)    ! Full representation of KSSe stiff matrix

      REAL(DOUBLE) , ALLOCATABLE      :: MFF_FULL(:,:)     ! Full representation of MFF mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MFS_FULL(:,:)     ! Full representation of MFS mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MSS_FULL(:,:)     ! Full representation of MSS mass matrix

      REAL(DOUBLE) , ALLOCATABLE      :: PF_FULL(:,:)      ! Full representation of PF load vectors
      REAL(DOUBLE) , ALLOCATABLE      :: PS_FULL(:,:)      ! Full representation of PS load vectors

      REAL(DOUBLE) , ALLOCATABLE      :: KAA_FULL(:,:)     ! Full representation of KAA stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KAO_FULL(:,:)     ! Full representation of KAO stiff matrix
      REAL(DOUBLE) , ALLOCATABLE      :: KOO_FULL(:,:)     ! Full representation of KOO stiff matrix

      REAL(DOUBLE) , ALLOCATABLE      :: MAA_FULL(:,:)     ! Full representation of MAA mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MAO_FULL(:,:)     ! Full representation of MAO mass matrix
      REAL(DOUBLE) , ALLOCATABLE      :: MOO_FULL(:,:)     ! Full representation of MOO mass matrix

      REAL(DOUBLE) , ALLOCATABLE      :: PA_FULL(:,:)      ! Full representation of PA load vectors
      REAL(DOUBLE) , ALLOCATABLE      :: PO_FULL(:,:)      ! Full representation of PA load vectors

      REAL(DOUBLE) , ALLOCATABLE      :: PFYS_FULL(:,:)    ! Full representation of PFYS load vector
      REAL(DOUBLE) , ALLOCATABLE      :: QSYS_FULL(:,:)    ! Full representation of QSYS load vector

      REAL(DOUBLE) , ALLOCATABLE      :: RMM_FULL(:,:)     ! Full representation of RMM constraint matrix
      REAL(DOUBLE) , ALLOCATABLE      :: GMN_FULL(:,:)     ! Full representation of GMN constraint matrix
      REAL(DOUBLE) , ALLOCATABLE      :: GMNt_FULL(:,:)    ! Full representation of transpose of GMN

      REAL(DOUBLE) , ALLOCATABLE      :: HMN_FULL(:,:)     ! Full representation of HMN MPC force recovery matrix
      REAL(DOUBLE) , ALLOCATABLE      :: JMN_FULL(:,:)     ! Full representation of JMN MPC force recovery matrix

      REAL(DOUBLE) , ALLOCATABLE      :: GOA_FULL(:,:)     ! Full representation of GOA constraint matrix
      REAL(DOUBLE) , ALLOCATABLE      :: GOAt_FULL(:,:)    ! Full representation of transpose of GOA

      REAL(DOUBLE) , ALLOCATABLE      :: DUM1(:,:)         ! Dummy matrix used in calculating KNN_FULL, MNN_FULL
      REAL(DOUBLE) , ALLOCATABLE      :: DUM2(:,:)         ! Dummy matrix used in calculating KNN_FULL, MNN_FULL
      REAL(DOUBLE) , ALLOCATABLE      :: DUM3(:,:)         ! Dummy matrix used in calculating KNN_FULL, MNN_FULL

      REAL(DOUBLE) , ALLOCATABLE      :: PHIZG_FULL(:,:)   ! Full representation of CB matrix PHIZG (NDOFG by NDOFR+NVEC)

      END MODULE FULL_MATRICES
