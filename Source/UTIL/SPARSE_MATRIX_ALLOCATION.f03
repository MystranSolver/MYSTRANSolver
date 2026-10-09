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

   MODULE SPARSE_MATRIX_ALLOCATION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ALLOCATE_SPARSE_MAT

   CONTAINS

      SUBROUTINE ALLOCATE_SPARSE_MAT ( NAME, NROWS, NTERMS, CALLING_SUBR )

! Allocate arrays for MYSTRAN sparse matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFM, NDOFO, NDOFS, NDOFR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SPARSE_MATRICES , ONLY      :  I_KGG   , J_KGG   , KGG   , I_MGG   , J_MGG   , MGG   , I_PG    , J_PG    , PG    ,       &
                                         I_KGGD  , J_KGGD  , KGGD  ,                                                               &
                                         I_KNN   , J_KNN   , KNN   , I_KNM   , J_KNM   , KNM   , I_KMM   , J_KMM   , KMM   ,       &
                                         I_KMN   , J_KMN   , KMN   ,                                                               &
                                         I_KNND  , J_KNND  , KNND  , I_KNMD  , J_KNMD  , KNMD  , I_KMMD  , J_KMMD  , KMMD  ,       &
                                         I_KMND  , J_KMND  , KMND  ,                                                               &
                                         I_MNN   , J_MNN   , MNN   , I_MNM   , J_MNM   , MNM   , I_MMN   , J_MMN   , MMN   ,       &
                                         I_MMM   , J_MMM   , MMM   , I_PN    , J_PN    , PN    , I_PM    , J_PM    , PM    ,       &
                                         I_GMNt  , J_GMNt  , GMNt  , I_HMN   , J_HMN   , HMN   , I_LMN   , J_LMN   , LMN   ,       &
                                         I_RMN   , J_RMN   , RMN   , I_RMM   , J_RMM   , RMM   , I_GMN   , J_GMN   , GMN   ,       &
                                         I_RMG   , J_RMG   , RMG

      USE SPARSE_MATRICES , ONLY      :  I_KFF   , J_KFF   , KFF   , I_KFS   , J_KFS   , KFS   , I_KSF   , J_KSF   , KSF   ,       &
                                         I_KSS   , J_KSS   , KSS   , I_KFSe  , J_KFSe  , KFSe  , I_KSSe  , J_KSSe  , KSSe  ,       &
                                         I_KFFD  , J_KFFD  , KFFD  , I_KFSD  , J_KFSD  , KFSD  , I_KSFD  , J_KSFD  , KSFD  ,       &
                                         I_KSSD  , J_KSSD  , KSSD  , I_KFSDe , J_KFSDe , KFSDe , I_KSSDe , J_KSSDe , KSSDe ,       &
                                         I_MFF   , J_MFF   , MFF   , I_MFS   , J_MFS   , MFS   , I_MSF   , J_MSF   , MSF   ,       &
                                         I_MSS   , J_MSS   , MSS   , I_PF    , J_PF    , PF    , I_PS    , J_PS    , PS    ,       &
                                         I_PF_TMP, J_PF_TMP, PF_TMP, I_PFYS  , J_PFYS  , PFYS  , I_PFYS1 , J_PFYS1 , PFYS1 ,       &
                                         I_QSYS  , J_QSYS  , QSYS

      USE SPARSE_MATRICES , ONLY      :  I_KAA   , J_KAA   , KAA   , I_KAO   , J_KAO   , KAO   ,                                   &
                                         I_KOO   , J_KOO   , KOO   , I_KOOs  , J_KOOs  , KOOs  ,                                   &
                                         I_KAAD  , J_KAAD  , KAAD  , I_KAOD  , J_KAOD  , KAOD  ,                                   &
                                         I_KOOD  , J_KOOD  , KOOD  , I_KOODs , J_KOODs , KOODs ,                                   &
                                         I_MAA   , J_MAA   , MAA   , I_MAO   , J_MAO   , MAO   , I_MOO   , J_MOO   , MOO   ,       &
                                         I_PA    , J_PA    , PA    , I_PO    , J_PO    , PO    ,                                   &
                                         I_GOA   , J_GOA   , GOA   , I_GOAt  , J_GOAt  , GOAt

      USE SPARSE_MATRICES , ONLY      :  I_KLL   , J_KLL   , KLL   , I_KLLs  , J_KLLs  , KLLs  , I_KRL   , J_KRL   , KRL   ,       &
                                         I_KRR   , J_KRR   , KRR   ,                                                               &
                                         I_KLLD  , J_KLLD  , KLLD  , I_KLLDs , J_KLLDs , KLLDs , I_KRLD  , J_KRLD  , KRLD  ,       &
                                         I_KRRD  , J_KRRD  , KRRD  , I_KLLDn , J_KLLDn , KLLDn ,                                   &
                                         I_MLL   , J_MLL   , MLL   , I_MLLn  , J_MLLn  , MLLn  , I_MLLs  , J_MLLs  , MLLs  ,       &
                                         I_MLR   , J_MLR   , MLR   ,                                                               &
                                         I_MRL   , J_MRL   , MRL   , I_MRR   , J_MRR   , MRR   , I_MPF0  , J_MPF0  , MPF0  ,       &
                                         I_KMSM  , J_KMSM  , KMSM  , I_KMSMn , J_KMSMn , KMSMn , I_KMSMs , J_KMSMs , KMSMs ,       &
                                         I_DLR   , J_DLR   , DLR   , I_DLRt  , J_DLRt  , DLRt  , I_IRR   , J_IRR   , IRR   ,       &
                                         I_PHIXA , J_PHIXA , PHIXA , I_PHIXG , J_PHIXG , PHIXG ,                                   &
                                         I_KRRcb , J_KRRcb , KRRcb , I_KRRcbn, J_KRRcbn, KRRcbn, I_KRRcbs, J_KRRcbs, KRRcbs,       &
                                         I_KXX   , J_KXX   , KXX   , I_LTM   , J_LTM   , LTM   ,                                   &
                                         I_MRN   , J_MRN   , MRN   , I_MRRcb , J_MRRcb , MRRcb , I_MRRcbn, J_MRRcbn, MRRcbn,       &
                                         I_MXX   , J_MXX   , MXX   , I_MXXn  , J_MXXn  , MXXn  ,                                   &
                                         I_IF_LTM, J_IF_LTM, IF_LTM, I_CG_LTM, J_CG_LTM, CG_LTM,                                   &
                                         I_PHIZL , J_PHIZL , PHIZL , I_PHIZL1, J_PHIZL1, PHIZL1, I_PHIZL1t, J_PHIZL1t, PHIZL1t,    &
                                         I_PHIZL2, J_PHIZL2, PHIZL2, I_PL    , J_PL    , PL    , I_PR     , J_PR     , PR

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_SPARSE_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name of the matrix to be allocated in sparse format
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(6*BYTE)               :: NAME1             ! Array name (used for output error message)
      CHARACTER(6*BYTE)               :: NAME2             ! Array name (used for output error message)

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows for matrix NAME
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of nonzero terms that will be in matrix NAME
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      MB_ALLOCATED = (REAL(LONG)*REAL(NROWS + 1 + NTERMS) + REAL(DOUBLE)*REAL(NTERMS))/ONEPP6
      JERR = 0

! Write blank line to screen

      IF ((NAME == 'KNN') .OR. (NAME == 'MNN') .OR. (NAME == 'PN')) THEN
         IF (ALLOCATED(I_KGG)) THEN
            IF (NDOFM == 0) THEN
         !xx   WRITE(SC1, * )
            ENDIF
         ENDIF
      ENDIF

      IF ((NAME == 'KFF') .OR. (NAME == 'MFF') .OR. (NAME == 'PF')) THEN
         IF (ALLOCATED(I_KNN)) THEN
            IF (NDOFS == 0) THEN
         !xx   WRITE(SC1, * )
            ENDIF
         ENDIF
      ENDIF

      IF ((NAME == 'KAA') .OR. (NAME == 'MAA') .OR. (NAME == 'PA')) THEN
         IF (ALLOCATED(I_KFF)) THEN
            IF (NDOFO == 0) THEN
         !xx   WRITE(SC1, * )
            ENDIF
         ENDIF
      ENDIF

      IF ((NAME == 'KLL') .OR. (NAME == 'MLL') .OR. (NAME == 'PL')) THEN
         IF (ALLOCATED(I_KAA)) THEN
            IF (NDOFR == 0) THEN
         !xx   WRITE(SC1, * )
            ENDIF
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! G-set arrays

      IF (NAME == 'KGG') THEN                              ! Allocate arrays for KGG

         NAME1 = 'I_KGG'
         IF (ALLOCATED(I_KGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KGG(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KGG(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KGG'
         IF (ALLOCATED(J_KGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KGG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KGG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KGG'
         IF (ALLOCATED(KGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KGG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KGG(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KGGD') THEN

         NAME1 = 'I_KGGD'
         IF (ALLOCATED(I_KGGD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KGGD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KGGD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KGGD'
         IF (ALLOCATED(J_KGGD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KGGD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KGGD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KGGD'
         IF (ALLOCATED(KGGD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KGGD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KGGD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MGG') THEN                         ! Allocate arrays for MGG

         NAME1 = 'I_MGG'
         IF (ALLOCATED(I_MGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MGG(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MGG(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MGG'
         IF (ALLOCATED(J_MGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MGG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MGG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MGG'
         IF (ALLOCATED(MGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MGG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MGG(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PG') THEN                          ! Allocate arrays for PG

         NAME1 = 'I_PG'
         IF (ALLOCATED(I_PG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PG(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PG(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PG'
         IF (ALLOCATED(J_PG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PG'
         IF (ALLOCATED(PG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PG(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RMG') THEN                         ! Allocate arrays for RMG

         NAME1 = 'I_RMG'
         IF (ALLOCATED(I_RMG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_RMG(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_RMG(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_RMG'
         IF (ALLOCATED(J_RMG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_RMG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_RMG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'RMG'
         IF (ALLOCATED(RMG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RMG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  RMG(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! N, M-set arrays

      ELSE IF (NAME == 'KNN') THEN                         ! Allocate arrays for KNN

         NAME1 = 'I_KNN'
         IF (ALLOCATED(I_KNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KNN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(I_KGG))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KNN(I) = I_KGG(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KNN(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KNN'
         IF (ALLOCATED(J_KNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KNN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(J_KGG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KNN(I) = J_KGG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KNN(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KNN'
         IF (ALLOCATED(KNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KNN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(KGG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KNN(I) = KGG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KNN(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KNM') THEN                         ! Allocate arrays for KNM

         NAME1 = 'I_KNM'
         IF (ALLOCATED(I_KNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KNM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KNM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KNM'
         IF (ALLOCATED(J_KNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KNM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KNM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KNM'
         IF (ALLOCATED(KNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KNM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KNM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMN') THEN                         ! Allocate arrays for KMN

         NAME1 = 'I_KMN'
         IF (ALLOCATED(I_KMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMN'
         IF (ALLOCATED(J_KMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMN'
         IF (ALLOCATED(KMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMM') THEN                         ! Allocate arrays for KMM

         NAME1 = 'I_KMM'
         IF (ALLOCATED(I_KMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMM'
         IF (ALLOCATED(J_KMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMM'
         IF (ALLOCATED(KMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KNND') THEN                         ! Allocate arrays for KNND
         NAME1 = 'I_KNND'
         IF (ALLOCATED(I_KNND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KNND(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(I_KGGD))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KNND(I) = I_KGGD(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KNND(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KNND'
         IF (ALLOCATED(J_KNND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KNND(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(J_KGGD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KNND(I) = J_KGGD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KNND(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KNND'
         IF (ALLOCATED(KNND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KNND(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(KGGD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KNND(I) = KGGD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KNND(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KNMD') THEN                         ! Allocate arrays for KNMD

         NAME1 = 'I_KNMD'
         IF (ALLOCATED(I_KNMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KNMD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KNMD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KNMD'
         IF (ALLOCATED(J_KNMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KNMD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KNMD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KNMD'
         IF (ALLOCATED(KNMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KNMD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KNMD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMND') THEN                         ! Allocate arrays for KMND

         NAME1 = 'I_KMND'
         IF (ALLOCATED(I_KMND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMND(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMND(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMND'
         IF (ALLOCATED(J_KMND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMND(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMND(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMND'
         IF (ALLOCATED(KMND)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMND(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMND(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMMD') THEN                         ! Allocate arrays for KMMD

         NAME1 = 'I_KMMD'
         IF (ALLOCATED(I_KMMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMMD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMMD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMMD'
         IF (ALLOCATED(J_KMMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMMD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMMD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMMD'
         IF (ALLOCATED(KMMD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMMD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMMD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MNN') THEN                         ! Allocate arrays for MNN

         NAME1 = 'I_MNN'
         IF (ALLOCATED(I_MNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MNN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(I_MGG))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MNN(I) = I_MGG(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MNN(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MNN'
         IF (ALLOCATED(J_MNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MNN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(J_MGG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MNN(I) = J_MGG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MNN(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MNN'
         IF (ALLOCATED(MNN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MNN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(MGG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MNN(I) = MGG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MNN(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MNM') THEN                         ! Allocate arrays for MNM

         NAME1 = 'I_MNM'
         IF (ALLOCATED(I_MNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MNM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MNM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MNM'
         IF (ALLOCATED(J_MNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MNM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MNM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MNM'
         IF (ALLOCATED(MNM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MNM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MNM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MMN') THEN                         ! Allocate arrays for MMN

         NAME1 = 'I_MMN'
         IF (ALLOCATED(I_MMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MMN'
         IF (ALLOCATED(J_MMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MMN'
         IF (ALLOCATED(MMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MMM') THEN                         ! Allocate arrays for MMM

         NAME1 = 'I_MMM'
         IF (ALLOCATED(I_MMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MMM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MMM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MMM'
         IF (ALLOCATED(J_MMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MMM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MMM'
         IF (ALLOCATED(MMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MMM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PN') THEN                          ! Allocate arrays for PN

         NAME1 = 'I_PN'
         IF (ALLOCATED(I_PN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(I_PG))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_PN(I) = I_PG(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_PN(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PN'
         IF (ALLOCATED(J_PN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(J_PG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_PN(I) = J_PG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_PN(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PN'
         IF (ALLOCATED(PN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFM == 0) .AND. (ALLOCATED(PG))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     PN(I) = PG(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     PN(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PM') THEN                          ! Allocate arrays for PM

         NAME1 = 'I_PM'
         IF (ALLOCATED(I_PM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PM'
         IF (ALLOCATED(J_PM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PM'
         IF (ALLOCATED(PM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RMN') THEN                         ! Allocate arrays for RMN

         NAME1 = 'I_RMN'
         IF (ALLOCATED(I_RMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_RMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_RMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_RMN'
         IF (ALLOCATED(J_RMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_RMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_RMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'RMN'
         IF (ALLOCATED(RMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  RMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'RMM') THEN                         ! Allocate arrays for RMM

         NAME1 = 'I_RMM'
         IF (ALLOCATED(I_RMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_RMM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_RMM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_RMM'
         IF (ALLOCATED(J_RMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_RMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_RMM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'RMM'
         IF (ALLOCATED(RMM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RMM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  RMM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'GMN') THEN                         ! Allocate arrays for GMN

         NAME1 = 'I_GMN'
         IF (ALLOCATED(I_GMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_GMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_GMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_GMN'
         IF (ALLOCATED(J_GMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_GMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_GMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'GMN'
         IF (ALLOCATED(GMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  GMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'GMNt') THEN                        ! Allocate arrays for GMNt

         NAME1 = 'I_GMNt'
         IF (ALLOCATED(I_GMNt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_GMNt(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_GMNt(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_GMNt'
         IF (ALLOCATED(J_GMNt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_GMNt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_GMNt(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'GMNt'
         IF (ALLOCATED(GMNt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GMNt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  GMNt(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'HMN') THEN                         ! Allocate arrays for HMN

         NAME1 = 'I_HMN'
         IF (ALLOCATED(I_HMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_HMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_HMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_HMN'
         IF (ALLOCATED(J_HMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_HMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_HMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'HMN'
         IF (ALLOCATED(HMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (HMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  HMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'LMN') THEN                         ! Allocate arrays for LMN

         NAME1 = 'I_LMN'
         IF (ALLOCATED(I_LMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_LMN(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_LMN(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_LMN'
         IF (ALLOCATED(J_LMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_LMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_LMN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'LMN'
         IF (ALLOCATED(LMN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (LMN(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  LMN(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! F, S-set arrays

      ELSE IF (NAME == 'KFF') THEN                         ! Allocate arrays for KFF

         NAME1 = 'I_KFF'
         IF (ALLOCATED(I_KFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFF(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(I_KNN))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KFF(I) = I_KNN(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KFF(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFF'
         IF (ALLOCATED(J_KFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(J_KNN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KFF(I) = J_KNN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KFF(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFF'
         IF (ALLOCATED(KFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(KNN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KFF(I) = KNN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KFF(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KFS') THEN                         ! Allocate arrays for KFS

         NAME1 = 'I_KFS'
         IF (ALLOCATED(I_KFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KFS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFS'
         IF (ALLOCATED(J_KFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KFS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFS'
         IF (ALLOCATED(KFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KFS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSF') THEN                         ! Allocate arrays for KSF

         NAME1 = 'I_KSF'
         IF (ALLOCATED(I_KSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSF(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSF(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSF'
         IF (ALLOCATED(J_KSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSF(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSF'
         IF (ALLOCATED(KSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSF(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSS') THEN                         ! Allocate arrays for KSS

         NAME1 = 'I_KSS'
         IF (ALLOCATED(I_KSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSS'
         IF (ALLOCATED(J_KSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSS'
         IF (ALLOCATED(KSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KFSe') THEN                        ! Allocate arrays for KFSe

         NAME1 = 'I_KFSe'
         IF (ALLOCATED(I_KFSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFSe(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KFSe(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFSe'
         IF (ALLOCATED(J_KFSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFSe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KFSe(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFSe'
         IF (ALLOCATED(KFSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFSe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KFSe(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSSe') THEN                        ! Allocate arrays for KSSe

         NAME1 = 'I_KSSe'
         IF (ALLOCATED(I_KSSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSSe(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSSe(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSSe'
         IF (ALLOCATED(J_KSSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSSe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSSe(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSSe'
         IF (ALLOCATED(KSSe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSSe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSSe(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KFFD') THEN                         ! Allocate arrays for KFFD

         NAME1 = 'I_KFFD'
         IF (ALLOCATED(I_KFFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFFD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(I_KNND))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KFFD(I) = I_KNND(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KFFD(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFFD'
         IF (ALLOCATED(J_KFFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFFD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(J_KNND))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KFFD(I) = J_KNND(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KFFD(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFFD'
         IF (ALLOCATED(KFFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFFD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(KNND))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KFFD(I) = KNND(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KFFD(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KFSD') THEN                         ! Allocate arrays for KFSD

         NAME1 = 'I_KFSD'
         IF (ALLOCATED(I_KFSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFSD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KFSD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFSD'
         IF (ALLOCATED(J_KFSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFSD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KFSD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFSD'
         IF (ALLOCATED(KFSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFSD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KFSD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSFD') THEN                         ! Allocate arrays for KSFD

         NAME1 = 'I_KSFD'
         IF (ALLOCATED(I_KSFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSFD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSFD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSFD'
         IF (ALLOCATED(J_KSFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSFD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSFD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSFD'
         IF (ALLOCATED(KSFD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSFD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSFD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSSD') THEN                         ! Allocate arrays for KSSD

         NAME1 = 'I_KSSD'
         IF (ALLOCATED(I_KSSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSSD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSSD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSSD'
         IF (ALLOCATED(J_KSSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSSD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSSD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSSD'
         IF (ALLOCATED(KSSD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSSD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSSD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KFSDe') THEN                        ! Allocate arrays for KFSDe

         NAME1 = 'I_KFSDe'
         IF (ALLOCATED(I_KFSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KFSDe(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KFSDe(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KFSDe'
         IF (ALLOCATED(J_KFSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KFSDe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KFSDe(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KFSDe'
         IF (ALLOCATED(KFSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KFSDe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KFSDe(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KSSDe') THEN                        ! Allocate arrays for KSSDe

         NAME1 = 'I_KSSDe'
         IF (ALLOCATED(I_KSSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KSSDe(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KSSDe(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KSSDe'
         IF (ALLOCATED(J_KSSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KSSDe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KSSDe(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KSSDe'
         IF (ALLOCATED(KSSDe)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KSSDe(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KSSDe(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MFF') THEN                         ! Allocate arrays for MFF

         NAME1 = 'I_MFF'
         IF (ALLOCATED(I_MFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MFF(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(I_MNN))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MFF(I) = I_MNN(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MFF(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MFF'
         IF (ALLOCATED(J_MFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MFF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(J_MNN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MFF(I) = J_MNN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MFF(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MFF'
         IF (ALLOCATED(MFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MFF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(MNN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MFF(I) = MNN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MFF(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MFS') THEN                         ! Allocate arrays for MFS

         NAME1 = 'I_MFS'
         IF (ALLOCATED(I_MFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MFS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MFS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MFS'
         IF (ALLOCATED(J_MFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MFS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MFS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MFS'
         IF (ALLOCATED(MFS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MFS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MFS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MSF') THEN                         ! Allocate arrays for MSF

         NAME1 = 'I_MSF'
         IF (ALLOCATED(I_MSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MSF(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MSF(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MSF'
         IF (ALLOCATED(J_MSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MSF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MSF(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MSF'
         IF (ALLOCATED(MSF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MSF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MSF(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MSS') THEN                         ! Allocate arrays for MSS

         NAME1 = 'I_MSS'
         IF (ALLOCATED(I_MSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MSS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MSS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MSS'
         IF (ALLOCATED(J_MSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MSS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MSS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MSS'
         IF (ALLOCATED(MSS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MSS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MSS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PF') THEN                          ! Allocate arrays for PF

         NAME1 = 'I_PF'
         IF (ALLOCATED(I_PF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PF(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(I_PN))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_PF(I) = I_PN(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_PF(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PF'
         IF (ALLOCATED(J_PF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(J_PN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_PF(I) = J_PN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_PF(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PF'
         IF (ALLOCATED(PF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PF(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFS == 0) .AND. (ALLOCATED(PN))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     PF(I) = PN(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     PF(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PS') THEN                          ! Allocate arrays for PS

         NAME1 = 'I_PS'
         IF (ALLOCATED(I_PS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PS'
         IF (ALLOCATED(J_PS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PS'
         IF (ALLOCATED(PS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PF_TMP') THEN                      ! Allocate arrays for PF_TMP

         NAME1 = 'I_PF_TMP'
         IF (ALLOCATED(I_PF_TMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PF_TMP(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PF_TMP(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PF_TMP'
         IF (ALLOCATED(J_PF_TMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PF_TMP(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PF_TMP(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PF_TMP'
         IF (ALLOCATED(PF_TMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PF_TMP(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PF_TMP(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PFYS') THEN                        ! Allocate arrays for PFYS

         NAME1 = 'I_PFYS'
         IF (ALLOCATED(I_PFYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PFYS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PFYS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PFYS'
         IF (ALLOCATED(J_PFYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PFYS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PFYS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PFYS'
         IF (ALLOCATED(PFYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PFYS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PFYS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PFYS1') THEN                       ! Allocate arrays for PFYS1

         NAME1 = 'I_PFYS1'
         IF (ALLOCATED(I_PFYS1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PFYS1(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PFYS1(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PFYS1'
         IF (ALLOCATED(J_PFYS1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PFYS1(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PFYS1(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PFYS1'
         IF (ALLOCATED(PFYS1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PFYS1(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PFYS1(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'QSYS') THEN                        ! Allocate arrays for QSYS

         NAME1 = 'I_QSYS'
         IF (ALLOCATED(I_QSYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_QSYS(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_QSYS(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_QSYS'
         IF (ALLOCATED(J_QSYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_QSYS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_QSYS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'QSYS'
         IF (ALLOCATED(QSYS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (QSYS(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  QSYS(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! A, O-set arrays

      ELSE IF (NAME == 'KAA') THEN                         ! Allocate arrays for KAA

         NAME1 = 'I_KAA'
         IF (ALLOCATED(I_KAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KAA(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(I_KFF))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KAA(I) = I_KFF(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KAA(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KAA'
         IF (ALLOCATED(J_KAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KAA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(J_KFF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KAA(I) = J_KFF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KAA(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KAA'
         IF (ALLOCATED(KAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KAA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(KFF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KAA(I) = KFF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KAA(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KAO') THEN                         ! Allocate arrays for KAO

         NAME1 = 'I_KAO'
         IF (ALLOCATED(I_KAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KAO(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KAO(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KAO'
         IF (ALLOCATED(J_KAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KAO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KAO(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KAO'
         IF (ALLOCATED(KAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KAO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KAO(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KOO') THEN                         ! Allocate arrays for KOO

         NAME1 = 'I_KOO'
         IF (ALLOCATED(I_KOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KOO(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KOO(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KOO'
         IF (ALLOCATED(J_KOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KOO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KOO(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KOO'
         IF (ALLOCATED(KOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KOO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KOO(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KOOs') THEN                        ! Allocate arrays for KOOs

         NAME1 = 'I_KOOs'
         IF (ALLOCATED(I_KOOs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KOOs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KOOs(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KOOs'
         IF (ALLOCATED(J_KOOs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KOOs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KOOs(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KOOs'
         IF (ALLOCATED(KOOs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KOOs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KOOs(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KAAD') THEN                         ! Allocate arrays for KAAD

         NAME1 = 'I_KAAD'
         IF (ALLOCATED(I_KAAD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KAAD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(I_KFFD))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KAAD(I) = I_KFFD(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KAAD(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KAAD'
         IF (ALLOCATED(J_KAAD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KAAD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(J_KFFD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KAAD(I) = J_KFFD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KAAD(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KAAD'
         IF (ALLOCATED(KAAD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KAAD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(KFFD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KAAD(I) = KFFD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KAAD(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KAOD') THEN                         ! Allocate arrays for KAOD

         NAME1 = 'I_KAOD'
         IF (ALLOCATED(I_KAOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KAOD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KAOD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KAOD'
         IF (ALLOCATED(J_KAOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KAOD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KAOD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KAOD'
         IF (ALLOCATED(KAOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KAOD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KAOD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KOOD') THEN                         ! Allocate arrays for KOOD

         NAME1 = 'I_KOOD'
         IF (ALLOCATED(I_KOOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KOOD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KOOD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KOOD'
         IF (ALLOCATED(J_KOOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KOOD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KOOD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KOOD'
         IF (ALLOCATED(KOOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KOOD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KOOD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KOODs') THEN                        ! Allocate arrays for KOODs

         NAME1 = 'I_KOODs'
         IF (ALLOCATED(I_KOODs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KOODs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KOODs(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KOODs'
         IF (ALLOCATED(J_KOODs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KOODs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KOODs(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KOODs'
         IF (ALLOCATED(KOODs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KOODs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KOODs(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MAA') THEN                         ! Allocate arrays for MAA

         NAME1 = 'I_MAA'
         IF (ALLOCATED(I_MAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MAA(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(I_MFF))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MAA(I) = I_MFF(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MAA(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MAA'
         IF (ALLOCATED(J_MAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MAA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(J_MFF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MAA(I) = J_MFF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MAA(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MAA'
         IF (ALLOCATED(MAA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MAA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(MFF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MAA(I) = MFF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MAA(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MAO') THEN                         ! Allocate arrays for MAO

         NAME1 = 'I_MAO'
         IF (ALLOCATED(I_MAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MAO(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MAO(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MAO'
         IF (ALLOCATED(J_MAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MAO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MAO(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MAO'
         IF (ALLOCATED(MAO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MAO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MAO(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MOO') THEN                         ! Allocate arrays for MOO

         NAME1 = 'I_MOO'
         IF (ALLOCATED(I_MOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MOO(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MOO(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MOO'
         IF (ALLOCATED(J_MOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MOO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MOO(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MOO'
         IF (ALLOCATED(MOO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MOO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MOO(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PA') THEN                          ! Allocate arrays for PA

         NAME1 = 'I_PA'
         IF (ALLOCATED(I_PA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PA(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(I_PF))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_PA(I) = I_PF(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_PA(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PA'
         IF (ALLOCATED(J_PA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(PF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_PA(I) = J_PF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_PA(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PA'
         IF (ALLOCATED(PA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFO == 0) .AND. (ALLOCATED(PF))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     PA(I) = PF(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     PA(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PO') THEN                          ! Allocate arrays for PO

         NAME1 = 'I_PO'
         IF (ALLOCATED(I_PO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PO(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PO(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PO'
         IF (ALLOCATED(J_PO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PO(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PO'
         IF (ALLOCATED(PO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PO(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PO(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'GOA') THEN                         ! ALLOCATE ARRAYS FOR GOA

         NAME1 = 'I_GOA'
         IF (ALLOCATED(I_GOA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_GOA(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_GOA(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_GOA'
         IF (ALLOCATED(J_GOA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_GOA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_GOA(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'GOA'
         IF (ALLOCATED(GOA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GOA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  GOA(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'GOAt') THEN                        ! Allocate arrays for GOAt

         NAME1 = 'I_GOAt'
         IF (ALLOCATED(I_GOAt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_GOAt(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_GOAt(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_GOAt'
         IF (ALLOCATED(J_GOAt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_GOAt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_GOAt(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'GOAt'
         IF (ALLOCATED(GOAt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GOAt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  GOAt(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! L, R-set arrays

      ELSE IF (NAME == 'KLL') THEN                         ! Allocate arrays for KLL

         NAME1 = 'I_KLL'
         IF (ALLOCATED(I_KLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KLL(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_KAA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KLL(I) = I_KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KLL(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KLL'
         IF (ALLOCATED(J_KLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KLL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_KAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KLL(I) = J_KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KLL(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KLL'
         IF (ALLOCATED(KLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KLL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(KAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KLL(I) = KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KLL(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KLLs') THEN                         ! Allocate arrays for KLLs

         NAME1 = 'I_KLLs'
         IF (ALLOCATED(I_KLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KLLs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_KAA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KLLs(I) = I_KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KLLs(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KLLs'
         IF (ALLOCATED(J_KLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KLLs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_KAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KLLs(I) = J_KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KLLs(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KLLs'
         IF (ALLOCATED(KLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KLLs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(KAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KLLs(I) = KAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KLLs(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRL') THEN                         ! Allocate arrays for KRL

         NAME1 = 'I_KRL'
         IF (ALLOCATED(I_KRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRL(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRL(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRL'
         IF (ALLOCATED(J_KRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRL(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRL'
         IF (ALLOCATED(KRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRL(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRR') THEN                         ! Allocate arrays for KRR

         NAME1 = 'I_KRR'
         IF (ALLOCATED(I_KRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRR'
         IF (ALLOCATED(J_KRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRR'
         IF (ALLOCATED(KRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KLLD') THEN                        ! Allocate arrays for KLLD

         NAME1 = 'I_KLLD'
         IF (ALLOCATED(I_KLLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KLLD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_KAAD))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_KLLD(I) = I_KAAD(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_KLLD(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KLLD'
         IF (ALLOCATED(J_KLLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KLLD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_KAAD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_KLLD(I) = J_KAAD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_KLLD(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KLLD'
         IF (ALLOCATED(KLLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KLLD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(KAAD))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     KLLD(I) = KAAD(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     KLLD(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KLLDn') THEN                         ! Allocate arrays for KLLDn

         NAME1 = 'I_KLLDn'
         IF (ALLOCATED(I_KLLDn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KLLDn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KLLDn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KLLDn'
         IF (ALLOCATED(J_KLLDn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KLLDn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KLLDn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KLLDn'
         IF (ALLOCATED(KLLDn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KLLDn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KLLDn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KLLDs') THEN                         ! Allocate arrays for KLLDs

         NAME1 = 'I_KLLDs'
         IF (ALLOCATED(I_KLLDs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KLLDs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KLLDs(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KLLDs'
         IF (ALLOCATED(J_KLLDs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KLLDs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KLLDs(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KLLDs'
         IF (ALLOCATED(KLLDs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KLLDs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KLLDs(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRLD') THEN                         ! Allocate arrays for KRLD

         NAME1 = 'I_KRLD'
         IF (ALLOCATED(I_KRLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRLD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRLD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRLD'
         IF (ALLOCATED(J_KRLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRLD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRLD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRLD'
         IF (ALLOCATED(KRLD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRLD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRLD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRRD') THEN                         ! Allocate arrays for KRRD

         NAME1 = 'I_KRRD'
         IF (ALLOCATED(I_KRRD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRRD(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRRD(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRRD'
         IF (ALLOCATED(J_KRRD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRRD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRRD(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRRD'
         IF (ALLOCATED(KRRD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRRD(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRRD(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MLL') THEN                         ! Allocate arrays for MLL

         NAME1 = 'I_MLL'
         IF (ALLOCATED(I_MLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MLL(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_MAA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MLL(I) = I_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MLL(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MLL'
         IF (ALLOCATED(J_MLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MLL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MLL(I) = J_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MLL(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MLL'
         IF (ALLOCATED(MLL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MLL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MLL(I) = MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MLL(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MLLn') THEN                        ! Allocate arrays for MLLn

         NAME1 = 'I_MLLn'
         IF (ALLOCATED(I_MLLn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MLLn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MLLn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MLLn'
         IF (ALLOCATED(J_MLLn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MLLn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MLLn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MLLn'
         IF (ALLOCATED(MLLn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MLLn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MLLn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MLLs') THEN                        ! Allocate arrays for MLLs

         NAME1 = 'I_MLLs'
         IF (ALLOCATED(I_MLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MLLs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_MAA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MLLs(I) = I_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MLLs(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MLLs'
         IF (ALLOCATED(J_MLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MLLs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MLLs(I) = J_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MLLs(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MLLs'
         IF (ALLOCATED(MLLs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MLLs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MLLs(I) = MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MLLs(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MPF0 ') THEN                       ! Allocate arrays for MPF0

         NAME1 = 'I_MPF0 '
         IF (ALLOCATED(I_MPF0 )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MPF0 (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_MAA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_MPF0 (I) = I_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_MPF0 (I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MPF0 '
         IF (ALLOCATED(J_MPF0 )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MPF0 (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(J_MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_MPF0 (I) = J_MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_MPF0 (I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MPF0 '
         IF (ALLOCATED(MPF0 )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPF0 (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(MAA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     MPF0 (I) = MAA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     MPF0 (I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MRL') THEN                         ! Allocate arrays for MRL

         NAME1 = 'I_MRL'
         IF (ALLOCATED(I_MRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MRL(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MRL(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MRL'
         IF (ALLOCATED(J_MRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MRL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MRL(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MRL'
         IF (ALLOCATED(MRL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MRL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MRL(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MLR') THEN                         ! Allocate arrays for MLR

         NAME1 = 'I_MLR'
         IF (ALLOCATED(I_MLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MLR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MLR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MLR'
         IF (ALLOCATED(J_MLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MLR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MLR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MLR'
         IF (ALLOCATED(MLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MLR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MLR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MRR') THEN                         ! Allocate arrays for MRR

         NAME1 = 'I_MRR'
         IF (ALLOCATED(I_MRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MRR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MRR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MRR'
         IF (ALLOCATED(J_MRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MRR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MRR'
         IF (ALLOCATED(MRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MRR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PL') THEN                          ! Allocate arrays for PL

         NAME1 = 'I_PL'
         IF (ALLOCATED(I_PL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PL(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(I_PA))) THEN
                  WRITE(SC1,12345,ADVANCE='NO') NAME1, NROWS+1, CR13
                  DO I=1,NROWS+1
                     I_PL(I) = I_PA(I)
                  ENDDO
               ELSE
                  DO I=1,NROWS+1
                     I_PL(I) = 1
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PL'
         IF (ALLOCATED(J_PL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(PA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     J_PL(I) = J_PA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     J_PL(I) = 0
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PL'
         IF (ALLOCATED(PL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PL(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               IF ((NDOFR == 0) .AND. (ALLOCATED(PA))) THEN
                  WRITE(SC1,22345,ADVANCE='NO') NAME1, NTERMS, CR13
                  DO I=1,NTERMS
                     PL(I) = PA(I)
                  ENDDO
               ELSE
                  DO I=1,NTERMS
                     PL(I) = ZERO
                  ENDDO
               ENDIF
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PR') THEN                          ! Allocate arrays for PR

         NAME1 = 'I_PR'
         IF (ALLOCATED(I_PR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PR'
         IF (ALLOCATED(J_PR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PR'
         IF (ALLOCATED(PR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMSM') THEN                        ! Allocate arrays for KMSM

         NAME1 = 'I_KMSM'
         IF (ALLOCATED(I_KMSM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMSM(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMSM(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMSM'
         IF (ALLOCATED(J_KMSM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMSM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMSM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMSM'
         IF (ALLOCATED(KMSM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMSM(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMSM(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMSMn') THEN                       ! Allocate arrays for KMSMn

         NAME1 = 'I_KMSMn'
         NAME1 = ''
         IF (ALLOCATED(I_KMSMn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMSMn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMSMn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMSMn'
         IF (ALLOCATED(J_KMSMn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMSMn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMSMn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMSMn'
         IF (ALLOCATED(KMSMn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMSMn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMSMn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KMSMs') THEN                       ! Allocate arrays for KMSMs

         NAME1 = 'I_KMSMs'
         NAME1 = ''
         IF (ALLOCATED(I_KMSMs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KMSMs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KMSMs(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KMSMs'
         IF (ALLOCATED(J_KMSMs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KMSMs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KMSMs(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KMSMs'
         IF (ALLOCATED(KMSMs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KMSMs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KMSMs(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'DLR') THEN                         ! Allocate arrays for DLR

         NAME1 = 'I_DLR'
         IF (ALLOCATED(I_DLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_DLR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_DLR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_DLR'
         IF (ALLOCATED(J_DLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_DLR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_DLR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'DLR'
         IF (ALLOCATED(DLR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (DLR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  DLR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'DLRt') THEN                        ! Allocate arrays for DLRt

         NAME1 = 'I_DLRt'
         IF (ALLOCATED(I_DLRt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_DLRt(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_DLRt(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_DLRt'
         IF (ALLOCATED(J_DLRt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_DLRt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_DLRt(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'DLRt'
         IF (ALLOCATED(DLRt)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (DLRt(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  DLRt(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'CG_LTM  ') THEN                    ! Allocate arrays for CG_LTM

         NAME1 = 'I_CG_LTM  '
         IF (ALLOCATED(I_CG_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_CG_LTM  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_CG_LTM  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_CG_LTM  '
         IF (ALLOCATED(J_CG_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_CG_LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_CG_LTM  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'CG_LTM  '
         IF (ALLOCATED(CG_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CG_LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  CG_LTM  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIZL  ') THEN                     ! Allocate arrays for PHIZL

         NAME1 = 'I_PHIZL  '
         IF (ALLOCATED(I_PHIZL  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIZL  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIZL  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIZL  '
         IF (ALLOCATED(J_PHIZL  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIZL  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIZL  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIZL  '
         IF (ALLOCATED(PHIZL  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIZL  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIZL  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIZL1  ') THEN                    ! Allocate arrays for PHIZL1

         NAME1 = 'I_PHIZL1  '
         IF (ALLOCATED(I_PHIZL1  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIZL1  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIZL1  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIZL1  '
         IF (ALLOCATED(J_PHIZL1  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIZL1  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIZL1  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIZL1  '
         IF (ALLOCATED(PHIZL1  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIZL1  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIZL1  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIZL1t  ') THEN                   ! Allocate arrays for PHIZL1t

         NAME1 = 'I_PHIZL1t  '
         IF (ALLOCATED(I_PHIZL1t  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIZL1t  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIZL1t  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIZL1t  '
         IF (ALLOCATED(J_PHIZL1t  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIZL1t  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIZL1t  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIZL1t  '
         IF (ALLOCATED(PHIZL1t  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIZL1t  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIZL1t  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIZL2  ') THEN                    ! Allocate arrays for PHIZL2

         NAME1 = 'I_PHIZL2  '
         IF (ALLOCATED(I_PHIZL2  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIZL2  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIZL2  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIZL2  '
         IF (ALLOCATED(J_PHIZL2  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIZL2  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIZL2  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIZL2  '
         IF (ALLOCATED(PHIZL2  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIZL2  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIZL2  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'IF_LTM  ') THEN                    ! Allocate arrays for IF_LTM

         NAME1 = 'I_IF_LTM  '
         IF (ALLOCATED(I_IF_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_IF_LTM  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_IF_LTM  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_IF_LTM  '
         IF (ALLOCATED(J_IF_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_IF_LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_IF_LTM  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'IF_LTM  '
         IF (ALLOCATED(IF_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (IF_LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  IF_LTM  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'IRR') THEN                         ! Allocate arrays for IRR

         NAME1 = 'I_IRR'
         IF (ALLOCATED(I_IRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_IRR(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_IRR(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_IRR'
         IF (ALLOCATED(J_IRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_IRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_IRR(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'IRR'
         IF (ALLOCATED(IRR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (IRR(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  IRR(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIXA') THEN                       ! Allocate arrays for PHIXA

         NAME1 = 'I_PHIXA'
         IF (ALLOCATED(I_PHIXA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIXA(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIXA(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIXA'
         IF (ALLOCATED(J_PHIXA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIXA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIXA(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIXA'
         IF (ALLOCATED(PHIXA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIXA(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIXA(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'PHIXG') THEN                       ! Allocate arrays for PHIXG

         NAME1 = 'I_PHIXG'
         IF (ALLOCATED(I_PHIXG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_PHIXG(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_PHIXG(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_PHIXG'
         IF (ALLOCATED(J_PHIXG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_PHIXG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_PHIXG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'PHIXG'
         IF (ALLOCATED(PHIXG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PHIXG(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  PHIXG(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRRcb') THEN                       ! Allocate arrays for KRRcb

         NAME1 = 'I_KRRcb'
         IF (ALLOCATED(I_KRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRRcb(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRRcb(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRRcb'
         IF (ALLOCATED(J_KRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRRcb(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRRcb(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRRcb'
         IF (ALLOCATED(KRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRRcb(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRRcb(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRRcbn') THEN                      ! Allocate arrays for KRRcbn

         NAME1 = 'I_KRRcbn'
         IF (ALLOCATED(I_KRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRRcbn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRRcbn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRRcbn'
         IF (ALLOCATED(J_KRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRRcbn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRRcbn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRRcbn'
         IF (ALLOCATED(KRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRRcbn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRRcbn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KRRcbs') THEN                      ! Allocate arrays for KRRcbs

         NAME1 = 'I_KRRcbs'
         IF (ALLOCATED(I_KRRcbs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KRRcbs(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KRRcbs(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KRRcbs'
         IF (ALLOCATED(J_KRRcbs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KRRcbs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KRRcbs(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KRRcbs'
         IF (ALLOCATED(KRRcbs)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KRRcbs(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KRRcbs(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'KXX  ') THEN                       ! Allocate arrays for KXX

         NAME1 = 'I_KXX  '
         IF (ALLOCATED(I_KXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_KXX  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_KXX  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_KXX  '
         IF (ALLOCATED(J_KXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_KXX  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_KXX  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'KXX  '
         IF (ALLOCATED(KXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KXX  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  KXX  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'LTM  ') THEN                       ! Allocate arrays for LTM

         NAME1 = 'I_LTM  '
         IF (ALLOCATED(I_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_LTM  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_LTM  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_LTM  '
         IF (ALLOCATED(J_LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_LTM  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'LTM  '
         IF (ALLOCATED(LTM  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (LTM  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  LTM  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MRN  ') THEN                       ! Allocate arrays for MRN

         NAME1 = 'I_MRN  '
         IF (ALLOCATED(I_MRN  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MRN  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MRN  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MRN  '
         IF (ALLOCATED(J_MRN  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MRN  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MRN  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MRN  '
         IF (ALLOCATED(MRN  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MRN  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MRN  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MRRcb') THEN                       ! Allocate arrays for MRRcb

         NAME1 = 'I_MRRcb'
         IF (ALLOCATED(I_MRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MRRcb(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MRRcb(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MRRcb'
         IF (ALLOCATED(J_MRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MRRcb(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MRRcb(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MRRcb'
         IF (ALLOCATED(MRRcb)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MRRcb(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MRRcb(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MRRcbn') THEN                      ! Allocate arrays for MRRcbn

         NAME1 = 'I_MRRcbn'
         IF (ALLOCATED(I_MRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MRRcbn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MRRcbn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MRRcbn'
         IF (ALLOCATED(J_MRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MRRcbn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MRRcbn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MRRcbn'
         IF (ALLOCATED(MRRcbn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MRRcbn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MRRcbn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MXX  ') THEN                       ! Allocate arrays for MXX

         NAME1 = 'I_MXX  '
         IF (ALLOCATED(I_MXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MXX  (NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MXX  (I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MXX  '
         IF (ALLOCATED(J_MXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MXX  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MXX  (I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MXX  '
         IF (ALLOCATED(MXX  )) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MXX  (NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MXX  (I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'MXXn') THEN                        ! Allocate arrays for MXXn

         NAME1 = 'I_MXXn'
         IF (ALLOCATED(I_MXXn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (I_MXXn(NROWS+1),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NROWS+1
                  I_MXXn(I) = 1
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'J_MXXn'
         IF (ALLOCATED(J_MXXn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (J_MXXn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  J_MXXn(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME1 = 'MXXn'
         IF (ALLOCATED(MXXn)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME1
            WRITE(F06,990) SUBR_NAME, NAME1
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MXXn(NTERMS),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NTERMS
                  MXXn(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               WRITE(F06,991) MB_ALLOCATED,NAME1,SUBR_NAME,IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF
      !WRITE(SC1,*) CR13

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ALLOCATION ERRORS WERE FOUND (',A,')'                                         &
                    ,/,14X,' WAS CALLED BY SUBR ',A)


12345 FORMAT(7X,'Equate ',A5,' with ',I8,' rows ',A)

22345 FORMAT(7X,'Equate ',A5,' with ',I8,' terms',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_SPARSE_MAT

   END MODULE SPARSE_MATRIX_ALLOCATION
