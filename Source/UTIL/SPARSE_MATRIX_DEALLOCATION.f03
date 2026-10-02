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

   MODULE SPARSE_MATRIX_DEALLOCATION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: DEALLOCATE_SPARSE_MAT

   CONTAINS

      SUBROUTINE DEALLOCATE_SPARSE_MAT ( NAME_IN )

! Deallocate arrays for MYSTRAN sparse matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS

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
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_SPARSE_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Array name (used for output error message)
      CHARACTER(6*BYTE)               :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! **********************************************************************************************************************************
! G-set arrays

      IF (NAME_IN == 'KGG') THEN                           ! Deallocate arrays for KGG

         NAME = 'I_KGG'
         IF (ALLOCATED(I_KGG)) THEN
            DEALLOCATE (I_KGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KGG'
         IF (ALLOCATED(J_KGG)) THEN
            DEALLOCATE (J_KGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KGG'
         IF (ALLOCATED(KGG)) THEN
            DEALLOCATE (KGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KGGD') THEN                      ! Deallocate arrays for KGGD

         NAME = 'I_KGGD'
         IF (ALLOCATED(I_KGGD)) THEN
            DEALLOCATE (I_KGGD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KGGD'
         IF (ALLOCATED(J_KGGD)) THEN
            DEALLOCATE (J_KGGD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KGGD'
         IF (ALLOCATED(KGGD)) THEN
            DEALLOCATE (KGGD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MGG') THEN                      ! Deallocate arrays for MGG

         NAME = 'I_MGG'
         IF (ALLOCATED(I_MGG)) THEN
            DEALLOCATE (I_MGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MGG'
         IF (ALLOCATED(J_MGG)) THEN
            DEALLOCATE (J_MGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MGG'
         IF (ALLOCATED(MGG)) THEN
            DEALLOCATE (MGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PG') THEN                       ! Deallocate arrays for PG

         NAME = 'I_PG'
         IF (ALLOCATED(I_PG)) THEN
            DEALLOCATE (I_PG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PG'
         IF (ALLOCATED(J_PG)) THEN
            DEALLOCATE (J_PG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PG'
         IF (ALLOCATED(PG)) THEN
            DEALLOCATE (PG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RMG') THEN                      ! Deallocate arrays for RMG

         NAME = 'I_RMG'
         IF (ALLOCATED(I_RMG)) THEN
            DEALLOCATE (I_RMG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_RMG'
         IF (ALLOCATED(J_RMG)) THEN
            DEALLOCATE (J_RMG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RMG'
         IF (ALLOCATED(RMG)) THEN
            DEALLOCATE (RMG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! N, M-set arrays

      ELSE IF (NAME_IN == 'KNN') THEN                      ! Deallocate arrays for KNN

         NAME = 'I_KNN'
         IF (ALLOCATED(I_KNN)) THEN
            DEALLOCATE (I_KNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KNN'
         IF (ALLOCATED(J_KNN)) THEN
            DEALLOCATE (J_KNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KNN'
         IF (ALLOCATED(KNN)) THEN
            DEALLOCATE (KNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KNM') THEN                      ! Deallocate arrays for KNM

         NAME = 'I_KNM'
         IF (ALLOCATED(I_KNM)) THEN
            DEALLOCATE (I_KNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KNM'
         IF (ALLOCATED(J_KNM)) THEN
            DEALLOCATE (J_KNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KNM'
         IF (ALLOCATED(KNM)) THEN
            DEALLOCATE (KNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMN') THEN                      ! Deallocate arrays for KMN

         NAME = 'I_KMN'
         IF (ALLOCATED(I_KMN)) THEN
            DEALLOCATE (I_KMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMN'
         IF (ALLOCATED(J_KMN)) THEN
            DEALLOCATE (J_KMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMN'
         IF (ALLOCATED(KMN)) THEN
            DEALLOCATE (KMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMM') THEN                      ! Deallocate arrays for KMM

         NAME = 'I_KMM'
         IF (ALLOCATED(I_KMM)) THEN
            DEALLOCATE (I_KMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMM'
         IF (ALLOCATED(J_KMM)) THEN
            DEALLOCATE (J_KMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMM'
         IF (ALLOCATED(KMM)) THEN
            DEALLOCATE (KMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KNND') THEN                      ! Deallocate arrays for KNND

         NAME = 'I_KNND'
         IF (ALLOCATED(I_KNND)) THEN
            DEALLOCATE (I_KNND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KNND'
         IF (ALLOCATED(J_KNND)) THEN
            DEALLOCATE (J_KNND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KNND'
         IF (ALLOCATED(KNND)) THEN
            DEALLOCATE (KNND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KNMD') THEN                      ! Deallocate arrays for KNMD

         NAME = 'I_KNMD'
         IF (ALLOCATED(I_KNMD)) THEN
            DEALLOCATE (I_KNMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KNMD'
         IF (ALLOCATED(J_KNMD)) THEN
            DEALLOCATE (J_KNMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KNMD'
         IF (ALLOCATED(KNMD)) THEN
            DEALLOCATE (KNMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMND') THEN                      ! Deallocate arrays for KMND

         NAME = 'I_KMND'
         IF (ALLOCATED(I_KMND)) THEN
            DEALLOCATE (I_KMND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMND'
         IF (ALLOCATED(J_KMND)) THEN
            DEALLOCATE (J_KMND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMND'
         IF (ALLOCATED(KMND)) THEN
            DEALLOCATE (KMND,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMMD') THEN                      ! Deallocate arrays for KMMD

         NAME = 'I_KMMD'
         IF (ALLOCATED(I_KMMD)) THEN
            DEALLOCATE (I_KMMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMMD'
         IF (ALLOCATED(J_KMMD)) THEN
            DEALLOCATE (J_KMMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMMD'
         IF (ALLOCATED(KMMD)) THEN
            DEALLOCATE (KMMD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MNN') THEN                      ! Deallocate arrays for MNN

         NAME = 'I_MNN'
         IF (ALLOCATED(I_MNN)) THEN
            DEALLOCATE (I_MNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MNN'
         IF (ALLOCATED(J_MNN)) THEN
            DEALLOCATE (J_MNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MNN'
         IF (ALLOCATED(MNN)) THEN
            DEALLOCATE (MNN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MNM') THEN                      ! Deallocate arrays for MNM

         NAME = 'I_MNM'
         IF (ALLOCATED(I_MNM)) THEN
            DEALLOCATE (I_MNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MNM'
         IF (ALLOCATED(J_MNM)) THEN
            DEALLOCATE (J_MNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MNM'
         IF (ALLOCATED(MNM)) THEN
            DEALLOCATE (MNM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MMN') THEN                      ! Deallocate arrays for MMN

         NAME = 'I_MMN'
         IF (ALLOCATED(I_MMN)) THEN
            DEALLOCATE (I_MMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MMN'
         IF (ALLOCATED(J_MMN)) THEN
            DEALLOCATE (J_MMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MMN'
         IF (ALLOCATED(MMN)) THEN
            DEALLOCATE (MMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MMM') THEN                      ! Deallocate arrays for MMM

         NAME = 'I_MMM'
         IF (ALLOCATED(I_MMM)) THEN
            DEALLOCATE (I_MMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MMM'
         IF (ALLOCATED(J_MMM)) THEN
            DEALLOCATE (J_MMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MMM'
         IF (ALLOCATED(MMM)) THEN
            DEALLOCATE (MMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PN') THEN                       ! Deallocate arrays for PN

         NAME = 'I_PN'
         IF (ALLOCATED(I_PN)) THEN
            DEALLOCATE (I_PN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PN'
         IF (ALLOCATED(J_PN)) THEN
            DEALLOCATE (J_PN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PN'
         IF (ALLOCATED(PN)) THEN
            DEALLOCATE (PN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PM') THEN                       ! Deallocate arrays for PM

         NAME = 'I_PM'
         IF (ALLOCATED(I_PM)) THEN
            DEALLOCATE (I_PM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PM'
         IF (ALLOCATED(J_PM)) THEN
            DEALLOCATE (J_PM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PM'
         IF (ALLOCATED(PM)) THEN
            DEALLOCATE (PM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RMN') THEN                      ! Deallocate arrays for RMN

         NAME = 'I_RMN'
         IF (ALLOCATED(I_RMN)) THEN
            DEALLOCATE (I_RMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_RMN'
         IF (ALLOCATED(J_RMN)) THEN
            DEALLOCATE (J_RMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RMN'
         IF (ALLOCATED(RMN)) THEN
            DEALLOCATE (RMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RMM') THEN                      ! Deallocate arrays for RMM

         NAME = 'I_RMM'
         IF (ALLOCATED(I_RMM)) THEN
            DEALLOCATE (I_RMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_RMM'
         IF (ALLOCATED(J_RMM)) THEN
            DEALLOCATE (J_RMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RMM'
         IF (ALLOCATED(RMM)) THEN
            DEALLOCATE (RMM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GMN') THEN                      ! Deallocate arrays for GMN

         NAME = 'I_GMN'
         IF (ALLOCATED(I_GMN)) THEN
            DEALLOCATE (I_GMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_GMN'
         IF (ALLOCATED(J_GMN)) THEN
            DEALLOCATE (J_GMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GMN'
         IF (ALLOCATED(GMN)) THEN
            DEALLOCATE (GMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GMNt') THEN                     ! Deallocate arrays for GMNt

         NAME = 'I_GMNt'
         IF (ALLOCATED(I_GMNt)) THEN
            DEALLOCATE (I_GMNt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_GMNt'
         IF (ALLOCATED(J_GMNt)) THEN
            DEALLOCATE (J_GMNt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GMNt'
         IF (ALLOCATED(GMNt)) THEN
            DEALLOCATE (GMNt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'HMN') THEN                      ! Deallocate arrays for HMN

         NAME = 'I_HMN'
         IF (ALLOCATED(I_HMN)) THEN
            DEALLOCATE (I_HMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_HMN'
         IF (ALLOCATED(J_HMN)) THEN
            DEALLOCATE (J_HMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'HMN'
         IF (ALLOCATED(HMN)) THEN
            DEALLOCATE (HMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'LMN') THEN                      ! Deallocate arrays for LMN

         NAME = 'I_LMN'
         IF (ALLOCATED(I_LMN)) THEN
            DEALLOCATE (I_LMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_LMN'
         IF (ALLOCATED(J_LMN)) THEN
            DEALLOCATE (J_LMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LMN'
         IF (ALLOCATED(LMN)) THEN
            DEALLOCATE (LMN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! F, S-set arrays

      ELSE IF (NAME_IN == 'KFF') THEN                      ! Deallocate arrays for KFF

         NAME = 'I_KFF'
         IF (ALLOCATED(I_KFF)) THEN
            DEALLOCATE (I_KFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFF'
         IF (ALLOCATED(J_KFF)) THEN
            DEALLOCATE (J_KFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFF'
         IF (ALLOCATED(KFF)) THEN
            DEALLOCATE (KFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KFS') THEN                      ! Deallocate arrays for KFS

         NAME = 'I_KFS'
         IF (ALLOCATED(I_KFS)) THEN
            DEALLOCATE (I_KFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFS'
         IF (ALLOCATED(J_KFS)) THEN
            DEALLOCATE (J_KFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFS'
         IF (ALLOCATED(KFS)) THEN
            DEALLOCATE (KFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSF') THEN                      ! Deallocate arrays for KSF

         NAME = 'I_KSF'
         IF (ALLOCATED(I_KSF)) THEN
            DEALLOCATE (I_KSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSF'
         IF (ALLOCATED(J_KSF)) THEN
            DEALLOCATE (J_KSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSF'
         IF (ALLOCATED(KSF)) THEN
            DEALLOCATE (KSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSS') THEN                      ! Deallocate arrays for KSS

         NAME = 'I_KSS'
         IF (ALLOCATED(I_KSS)) THEN
            DEALLOCATE (I_KSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSS'
         IF (ALLOCATED(J_KSS)) THEN
            DEALLOCATE (J_KSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSS'
         IF (ALLOCATED(KSS)) THEN
            DEALLOCATE (KSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KFSe') THEN                     ! Deallocate arrays for KFSe

         NAME = 'I_KFSe'
         IF (ALLOCATED(I_KFSe)) THEN
            DEALLOCATE (I_KFSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFSe'
         IF (ALLOCATED(J_KFSe)) THEN
            DEALLOCATE (J_KFSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFSe'
         IF (ALLOCATED(KFSe)) THEN
            DEALLOCATE (KFSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSSe') THEN                     ! Deallocate arrays for KSSe

         NAME = 'I_KSSe'
         IF (ALLOCATED(I_KSSe)) THEN
            DEALLOCATE (I_KSSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSSe'
         IF (ALLOCATED(J_KSSe)) THEN
            DEALLOCATE (J_KSSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSSe'
         IF (ALLOCATED(KSSe)) THEN
            DEALLOCATE (KSSe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KFFD') THEN                      ! Deallocate arrays for KFFD

         NAME = 'I_KFFD'
         IF (ALLOCATED(I_KFFD)) THEN
            DEALLOCATE (I_KFFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFFD'
         IF (ALLOCATED(J_KFFD)) THEN
            DEALLOCATE (J_KFFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFFD'
         IF (ALLOCATED(KFFD)) THEN
            DEALLOCATE (KFFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KFSD') THEN                      ! Deallocate arrays for KFSD

         NAME = 'I_KFSD'
         IF (ALLOCATED(I_KFSD)) THEN
            DEALLOCATE (I_KFSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFSD'
         IF (ALLOCATED(J_KFSD)) THEN
            DEALLOCATE (J_KFSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFSD'
         IF (ALLOCATED(KFSD)) THEN
            DEALLOCATE (KFSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSFD') THEN                      ! Deallocate arrays for KSFD

         NAME = 'I_KSFD'
         IF (ALLOCATED(I_KSFD)) THEN
            DEALLOCATE (I_KSFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSFD'
         IF (ALLOCATED(J_KSFD)) THEN
            DEALLOCATE (J_KSFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSFD'
         IF (ALLOCATED(KSFD)) THEN
            DEALLOCATE (KSFD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSSD') THEN                      ! Deallocate arrays for KSSD

         NAME = 'I_KSSD'
         IF (ALLOCATED(I_KSSD)) THEN
            DEALLOCATE (I_KSSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSSD'
         IF (ALLOCATED(J_KSSD)) THEN
            DEALLOCATE (J_KSSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSSD'
         IF (ALLOCATED(KSSD)) THEN
            DEALLOCATE (KSSD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KFSDe') THEN                     ! Deallocate arrays for KFSDe

         NAME = 'I_KFSDe'
         IF (ALLOCATED(I_KFSDe)) THEN
            DEALLOCATE (I_KFSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KFSDe'
         IF (ALLOCATED(J_KFSDe)) THEN
            DEALLOCATE (J_KFSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KFSDe'
         IF (ALLOCATED(KFSDe)) THEN
            DEALLOCATE (KFSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KSSDe') THEN                     ! Deallocate arrays for KSSDe

         NAME = 'I_KSSDe'
         IF (ALLOCATED(I_KSSDe)) THEN
            DEALLOCATE (I_KSSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KSSDe'
         IF (ALLOCATED(J_KSSDe)) THEN
            DEALLOCATE (J_KSSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KSSDe'
         IF (ALLOCATED(KSSDe)) THEN
            DEALLOCATE (KSSDe,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MFF') THEN                      ! Deallocate arrays for MFF

         NAME = 'I_MFF'
         IF (ALLOCATED(I_MFF)) THEN
            DEALLOCATE (I_MFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MFF'
         IF (ALLOCATED(J_MFF)) THEN
            DEALLOCATE (J_MFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MFF'
         IF (ALLOCATED(MFF)) THEN
            DEALLOCATE (MFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MFS') THEN                      ! Deallocate arrays for MFS

         NAME = 'I_MFS'
         IF (ALLOCATED(I_MFS)) THEN
            DEALLOCATE (I_MFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MFS'
         IF (ALLOCATED(J_MFS)) THEN
            DEALLOCATE (J_MFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MFS'
         IF (ALLOCATED(MFS)) THEN
            DEALLOCATE (MFS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MSF') THEN                      ! Deallocate arrays for MSF

         NAME = 'I_MSF'
         IF (ALLOCATED(I_MSF)) THEN
            DEALLOCATE (I_MSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MSF'
         IF (ALLOCATED(J_MSF)) THEN
            DEALLOCATE (J_MSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MSF'
         IF (ALLOCATED(MSF)) THEN
            DEALLOCATE (MSF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MSS') THEN                      ! Deallocate arrays for MSS

         NAME = 'I_MSS'
         IF (ALLOCATED(I_MSS)) THEN
            DEALLOCATE (I_MSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MSS'
         IF (ALLOCATED(J_MSS)) THEN
            DEALLOCATE (J_MSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MSS'
         IF (ALLOCATED(MSS)) THEN
            DEALLOCATE (MSS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PF') THEN                       ! Deallocate arrays for PF

         NAME = 'I_PF'
         IF (ALLOCATED(I_PF)) THEN
            DEALLOCATE (I_PF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PF'
         IF (ALLOCATED(J_PF)) THEN
            DEALLOCATE (J_PF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PF'
         IF (ALLOCATED(PF)) THEN
            DEALLOCATE (PF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PS') THEN                       ! Deallocate arrays for PS

         NAME = 'I_PS'
         IF (ALLOCATED(I_PS)) THEN
            DEALLOCATE (I_PS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PS'
         IF (ALLOCATED(J_PS)) THEN
            DEALLOCATE (J_PS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PS'
         IF (ALLOCATED(PS)) THEN
            DEALLOCATE (PS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PF_TMP') THEN                      ! Deallocate arrays for PF_TMP

         NAME = 'I_PF_TMP'
         IF (ALLOCATED(I_PF_TMP)) THEN
            DEALLOCATE (I_PF_TMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PF_TMP'
         IF (ALLOCATED(J_PF_TMP)) THEN
            DEALLOCATE (J_PF_TMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PF_TMP'
         IF (ALLOCATED(PF_TMP)) THEN
            DEALLOCATE (PF_TMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PFYS') THEN                     ! Deallocate arrays for PFYS

         NAME = 'I_PFYS'
         IF (ALLOCATED(I_PFYS)) THEN
            DEALLOCATE (I_PFYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PFYS'
         IF (ALLOCATED(J_PFYS)) THEN
            DEALLOCATE (J_PFYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PFYS'
         IF (ALLOCATED(PFYS)) THEN
            DEALLOCATE (PFYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PFYS1') THEN                    ! Deallocate arrays for PFYS1

         NAME = 'I_PFYS1'
         IF (ALLOCATED(I_PFYS1)) THEN
            DEALLOCATE (I_PFYS1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PFYS1'
         IF (ALLOCATED(J_PFYS1)) THEN
            DEALLOCATE (J_PFYS1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PFYS1'
         IF (ALLOCATED(PFYS1)) THEN
            DEALLOCATE (PFYS1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'QSYS') THEN                     ! Deallocate arrays for QSYS

         NAME = 'I_QSYS'
         IF (ALLOCATED(I_QSYS)) THEN
            DEALLOCATE (I_QSYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_QSYS'
         IF (ALLOCATED(J_QSYS)) THEN
            DEALLOCATE (J_QSYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'QSYS'
         IF (ALLOCATED(QSYS)) THEN
            DEALLOCATE (QSYS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! A, O-set arrays

      ELSE IF (NAME_IN == 'KAA') THEN                      ! Deallocate arrays for KAA

         NAME = 'I_KAA'
         IF (ALLOCATED(I_KAA)) THEN
            DEALLOCATE (I_KAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KAA'
         IF (ALLOCATED(J_KAA)) THEN
            DEALLOCATE (J_KAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KAA'
         IF (ALLOCATED(KAA)) THEN
            DEALLOCATE (KAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KAO') THEN                      ! Deallocate arrays for KAO

         NAME = 'I_KAO'
         IF (ALLOCATED(I_KAO)) THEN
            DEALLOCATE (I_KAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KAO'
         IF (ALLOCATED(J_KAO)) THEN
            DEALLOCATE (J_KAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KAO'
         IF (ALLOCATED(KAO)) THEN
            DEALLOCATE (KAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KOO') THEN                      ! Deallocate arrays for KOO

         NAME = 'I_KOO'
         IF (ALLOCATED(I_KOO)) THEN
            DEALLOCATE (I_KOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KOO'
         IF (ALLOCATED(J_KOO)) THEN
            DEALLOCATE (J_KOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KOO'
         IF (ALLOCATED(KOO)) THEN
            DEALLOCATE (KOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KOOs') THEN                     ! Deallocate arrays for KOOs

         NAME = 'I_KOOs'
         IF (ALLOCATED(I_KOOs)) THEN
            DEALLOCATE (I_KOOs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KOOs'
         IF (ALLOCATED(J_KOOs)) THEN
            DEALLOCATE (J_KOOs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KOOs'
         IF (ALLOCATED(KOOs)) THEN
            DEALLOCATE (KOOs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KAAD') THEN                      ! Deallocate arrays for KAAD

         NAME = 'I_KAAD'
         IF (ALLOCATED(I_KAAD)) THEN
            DEALLOCATE (I_KAAD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KAAD'
         IF (ALLOCATED(J_KAAD)) THEN
            DEALLOCATE (J_KAAD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KAAD'
         IF (ALLOCATED(KAAD)) THEN
            DEALLOCATE (KAAD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KAOD') THEN                      ! Deallocate arrays for KAOD

         NAME = 'I_KAOD'
         IF (ALLOCATED(I_KAOD)) THEN
            DEALLOCATE (I_KAOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KAOD'
         IF (ALLOCATED(J_KAOD)) THEN
            DEALLOCATE (J_KAOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KAOD'
         IF (ALLOCATED(KAOD)) THEN
            DEALLOCATE (KAOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KOOD') THEN                      ! Deallocate arrays for KOOD

         NAME = 'I_KOOD'
         IF (ALLOCATED(I_KOOD)) THEN
            DEALLOCATE (I_KOOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KOOD'
         IF (ALLOCATED(J_KOOD)) THEN
            DEALLOCATE (J_KOOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KOOD'
         IF (ALLOCATED(KOOD)) THEN
            DEALLOCATE (KOOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KOODs') THEN                     ! Deallocate arrays for KOODs

         NAME = 'I_KOODs'
         IF (ALLOCATED(I_KOODs)) THEN
            DEALLOCATE (I_KOODs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KOODs'
         IF (ALLOCATED(J_KOODs)) THEN
            DEALLOCATE (J_KOODs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KOODs'
         IF (ALLOCATED(KOODs)) THEN
            DEALLOCATE (KOODs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MAA') THEN                      ! Deallocate arrays for MAA

         NAME = 'I_MAA'
         IF (ALLOCATED(I_MAA)) THEN
            DEALLOCATE (I_MAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MAA'
         IF (ALLOCATED(J_MAA)) THEN
            DEALLOCATE (J_MAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MAA'
         IF (ALLOCATED(MAA)) THEN
            DEALLOCATE (MAA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MAO') THEN                      ! Deallocate arrays for MAO

         NAME = 'I_MAO'
         IF (ALLOCATED(I_MAO)) THEN
            DEALLOCATE (I_MAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MAO'
         IF (ALLOCATED(J_MAO)) THEN
            DEALLOCATE (J_MAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MAO'
         IF (ALLOCATED(MAO)) THEN
            DEALLOCATE (MAO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MOO') THEN                      ! Deallocate arrays for MOO

         NAME = 'I_MOO'
         IF (ALLOCATED(I_MOO)) THEN
            DEALLOCATE (I_MOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MOO'
         IF (ALLOCATED(J_MOO)) THEN
            DEALLOCATE (J_MOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MOO'
         IF (ALLOCATED(MOO)) THEN
            DEALLOCATE (MOO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PA') THEN                       ! Deallocate arrays for PA

         NAME = 'I_PA'
         IF (ALLOCATED(I_PA)) THEN
            DEALLOCATE (I_PA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PA'
         IF (ALLOCATED(J_PA)) THEN
            DEALLOCATE (J_PA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PA'
         IF (ALLOCATED(PA)) THEN
            DEALLOCATE (PA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PO') THEN                       ! Deallocate arrays for PO

         NAME = 'I_PO'
         IF (ALLOCATED(I_PO)) THEN
            DEALLOCATE (I_PO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PO'
         IF (ALLOCATED(J_PO)) THEN
            DEALLOCATE (J_PO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PO'
         IF (ALLOCATED(PO)) THEN
            DEALLOCATE (PO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GOA') THEN                      ! Deallocate arrays for GOA

         NAME = 'I_GOA'
         IF (ALLOCATED(I_GOA)) THEN
            DEALLOCATE (I_GOA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_GOA'
         IF (ALLOCATED(J_GOA)) THEN
            DEALLOCATE (J_GOA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GOA'
         IF (ALLOCATED(GOA)) THEN
            DEALLOCATE (GOA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GOAt') THEN                     ! Deallocate arrays for GOAt

         NAME = 'I_GOAt'
         IF (ALLOCATED(I_GOAt)) THEN
            DEALLOCATE (I_GOAt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_GOAt'
         IF (ALLOCATED(J_GOAt)) THEN
            DEALLOCATE (J_GOAt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GOAt'
         IF (ALLOCATED(GOAt)) THEN
            DEALLOCATE (GOAt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************
! L, R-set arrays

      ELSE IF (NAME_IN == 'KLL') THEN                      ! Deallocate arrays for KLL

         NAME = 'I_KLL'
         IF (ALLOCATED(I_KLL)) THEN
            DEALLOCATE (I_KLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KLL'
         IF (ALLOCATED(J_KLL)) THEN
            DEALLOCATE (J_KLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KLL'
         IF (ALLOCATED(KLL)) THEN
            DEALLOCATE (KLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KLLs') THEN                     ! Deallocate arrays for KLLs

         NAME = 'I_KLLs'
         IF (ALLOCATED(I_KLLs)) THEN
            DEALLOCATE (I_KLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KLLs'
         IF (ALLOCATED(J_KLLs)) THEN
            DEALLOCATE (J_KLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KLLs'
         IF (ALLOCATED(KLLs)) THEN
            DEALLOCATE (KLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRL') THEN                      ! Deallocate arrays for KRL

         NAME = 'I_KRL'
         IF (ALLOCATED(I_KRL)) THEN
            DEALLOCATE (I_KRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRL'
         IF (ALLOCATED(J_KRL)) THEN
            DEALLOCATE (J_KRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRL'
         IF (ALLOCATED(KRL)) THEN
            DEALLOCATE (KRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRR') THEN                      ! Deallocate arrays for KRR

         NAME = 'I_KRR'
         IF (ALLOCATED(I_KRR)) THEN
            DEALLOCATE (I_KRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRR'
         IF (ALLOCATED(J_KRR)) THEN
            DEALLOCATE (J_KRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRR'
         IF (ALLOCATED(KRR)) THEN
            DEALLOCATE (KRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KLLD') THEN                     ! Deallocate arrays for KLLD

         NAME = 'I_KLLD'
         IF (ALLOCATED(I_KLLD)) THEN
            DEALLOCATE (I_KLLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KLLD'
         IF (ALLOCATED(J_KLLD)) THEN
            DEALLOCATE (J_KLLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KLLD'
         IF (ALLOCATED(KLLD)) THEN
            DEALLOCATE (KLLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KLLDn') THEN                    ! Deallocate arrays for KLLDn

         NAME = 'I_KLLDn'
         IF (ALLOCATED(I_KLLDn)) THEN
            DEALLOCATE (I_KLLDn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KLLDn'
         IF (ALLOCATED(J_KLLDn)) THEN
            DEALLOCATE (J_KLLDn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KLLDn'
         IF (ALLOCATED(KLLDn)) THEN
            DEALLOCATE (KLLDn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KLLDs') THEN                    ! Deallocate arrays for KLLDs

         NAME = 'I_KLLDs'
         IF (ALLOCATED(I_KLLDs)) THEN
            DEALLOCATE (I_KLLDs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KLLDs'
         IF (ALLOCATED(J_KLLDs)) THEN
            DEALLOCATE (J_KLLDs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KLLDs'
         IF (ALLOCATED(KLLDs)) THEN
            DEALLOCATE (KLLDs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRLD') THEN                      ! Deallocate arrays for KRLD

         NAME = 'I_KRLD'
         IF (ALLOCATED(I_KRLD)) THEN
            DEALLOCATE (I_KRLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRLD'
         IF (ALLOCATED(J_KRLD)) THEN
            DEALLOCATE (J_KRLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRLD'
         IF (ALLOCATED(KRLD)) THEN
            DEALLOCATE (KRLD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRRD') THEN                      ! Deallocate arrays for KRRD

         NAME = 'I_KRRD'
         IF (ALLOCATED(I_KRRD)) THEN
            DEALLOCATE (I_KRRD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRRD'
         IF (ALLOCATED(J_KRRD)) THEN
            DEALLOCATE (J_KRRD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRRD'
         IF (ALLOCATED(KRRD)) THEN
            DEALLOCATE (KRRD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MLL') THEN                      ! Deallocate arrays for MLL

         NAME = 'I_MLL'
         IF (ALLOCATED(I_MLL)) THEN
            DEALLOCATE (I_MLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MLL'
         IF (ALLOCATED(J_MLL)) THEN
            DEALLOCATE (J_MLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MLL'
         IF (ALLOCATED(MLL)) THEN
            DEALLOCATE (MLL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MLLn') THEN                     ! Deallocate arrays for MLLn

         NAME = 'I_MLLn'
         IF (ALLOCATED(I_MLLn)) THEN
            DEALLOCATE (I_MLLn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MLLn'
         IF (ALLOCATED(J_MLLn)) THEN
            DEALLOCATE (J_MLLn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MLLn'
         IF (ALLOCATED(MLLn)) THEN
            DEALLOCATE (MLLn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MLLs') THEN                     ! Deallocate arrays for MLLs

         NAME = 'I_MLLs'
         IF (ALLOCATED(I_MLLs)) THEN
            DEALLOCATE (I_MLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MLLs'
         IF (ALLOCATED(J_MLLs)) THEN
            DEALLOCATE (J_MLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MLLs'
         IF (ALLOCATED(MLLs)) THEN
            DEALLOCATE (MLLs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MLR') THEN                      ! Deallocate arrays for MLR

         NAME = 'I_MLR'
         IF (ALLOCATED(I_MLR)) THEN
            DEALLOCATE (I_MLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MLR'
         IF (ALLOCATED(J_MLR)) THEN
            DEALLOCATE (J_MLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MLR'
         IF (ALLOCATED(MLR)) THEN
            DEALLOCATE (MLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPF0 ') THEN                    ! Deallocate arrays for MPF0

         NAME = 'I_MPF0 '
         IF (ALLOCATED(I_MPF0 )) THEN
            DEALLOCATE (I_MPF0 ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MPF0 '
         IF (ALLOCATED(J_MPF0 )) THEN
            DEALLOCATE (J_MPF0 ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MPF0 '
         IF (ALLOCATED(MPF0 )) THEN
            DEALLOCATE (MPF0 ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MRL') THEN                      ! Deallocate arrays for MRL

         NAME = 'I_MRL'
         IF (ALLOCATED(I_MRL)) THEN
            DEALLOCATE (I_MRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MRL'
         IF (ALLOCATED(J_MRL)) THEN
            DEALLOCATE (J_MRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MRL'
         IF (ALLOCATED(MRL)) THEN
            DEALLOCATE (MRL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MRR') THEN                      ! Deallocate arrays for MRR

         NAME = 'I_MRR'
         IF (ALLOCATED(I_MRR)) THEN
            DEALLOCATE (I_MRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MRR'
         IF (ALLOCATED(J_MRR)) THEN
            DEALLOCATE (J_MRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MRR'
         IF (ALLOCATED(MRR)) THEN
            DEALLOCATE (MRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PL') THEN                       ! Deallocate arrays for PL

         NAME = 'I_PL'
         IF (ALLOCATED(I_PL)) THEN
            DEALLOCATE (I_PL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PL'
         IF (ALLOCATED(J_PL)) THEN
            DEALLOCATE (J_PL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PL'
         IF (ALLOCATED(PL)) THEN
            DEALLOCATE (PL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PR') THEN                       ! Deallocate arrays for PR

         NAME = 'I_PR'
         IF (ALLOCATED(I_PR)) THEN
            DEALLOCATE (I_PR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PR'
         IF (ALLOCATED(J_PR)) THEN
            DEALLOCATE (J_PR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PR'
         IF (ALLOCATED(PR)) THEN
            DEALLOCATE (PR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMSM') THEN                     ! Deallocate arrays for KMSM

         NAME = 'I_KMSM'
         IF (ALLOCATED(I_KMSM)) THEN
            DEALLOCATE (I_KMSM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMSM'
         IF (ALLOCATED(J_KMSM)) THEN
            DEALLOCATE (J_KMSM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMSM'
         IF (ALLOCATED(KMSM)) THEN
            DEALLOCATE (KMSM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMSMn') THEN                    ! Deallocate arrays for KMSMn

         NAME = 'I_KMSMn'
         IF (ALLOCATED(I_KMSMn)) THEN
            DEALLOCATE (I_KMSMn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMSMn'
         IF (ALLOCATED(J_KMSMn)) THEN
            DEALLOCATE (J_KMSMn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMSMn'
         IF (ALLOCATED(KMSMn)) THEN
            DEALLOCATE (KMSMn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KMSMs') THEN                    ! Deallocate arrays for KMSMs

         NAME = 'I_KMSMs'
         IF (ALLOCATED(I_KMSMs)) THEN
            DEALLOCATE (I_KMSMs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KMSMs'
         IF (ALLOCATED(J_KMSMs)) THEN
            DEALLOCATE (J_KMSMs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KMSMs'
         IF (ALLOCATED(KMSMs)) THEN
            DEALLOCATE (KMSMs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'DLR') THEN                      ! Deallocate arrays for DLR

         NAME = 'I_DLR'
         IF (ALLOCATED(I_DLR)) THEN
            DEALLOCATE (I_DLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_DLR'
         IF (ALLOCATED(J_DLR)) THEN
            DEALLOCATE (J_DLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DLR'
         IF (ALLOCATED(DLR)) THEN
            DEALLOCATE (DLR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'DLRt') THEN                     ! Deallocate arrays for DLRt

         NAME = 'I_DLRt'
         IF (ALLOCATED(I_DLRt)) THEN
            DEALLOCATE (I_DLRt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_DLRt'
         IF (ALLOCATED(J_DLRt)) THEN
            DEALLOCATE (J_DLRt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DLRt'
         IF (ALLOCATED(DLRt)) THEN
            DEALLOCATE (DLRt,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CG_LTM  ') THEN                 ! Deallocate arrays for CG_LTM

         NAME = 'I_CG_LTM  '
         IF (ALLOCATED(I_CG_LTM  )) THEN
            DEALLOCATE (I_CG_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_CG_LTM  '
         IF (ALLOCATED(J_CG_LTM  )) THEN
            DEALLOCATE (J_CG_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'CG_LTM  '
         IF (ALLOCATED(CG_LTM  )) THEN
            DEALLOCATE (CG_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIZL  ') THEN                  ! Deallocate arrays for PHIZL

         NAME = 'I_PHIZL  '
         IF (ALLOCATED(I_PHIZL  )) THEN
            DEALLOCATE (I_PHIZL  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIZL  '
         IF (ALLOCATED(J_PHIZL  )) THEN
            DEALLOCATE (J_PHIZL  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIZL  '
         IF (ALLOCATED(PHIZL  )) THEN
            DEALLOCATE (PHIZL  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIZL1  ') THEN                 ! Deallocate arrays for PHIZL1

         NAME = 'I_PHIZL1  '
         IF (ALLOCATED(I_PHIZL1  )) THEN
            DEALLOCATE (I_PHIZL1  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIZL1  '
         IF (ALLOCATED(J_PHIZL1  )) THEN
            DEALLOCATE (J_PHIZL1  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIZL1  '
         IF (ALLOCATED(PHIZL1  )) THEN
            DEALLOCATE (PHIZL1  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIZL1t  ') THEN                ! Deallocate arrays for PHIZL1t

         NAME = 'I_PHIZL1t  '
         IF (ALLOCATED(I_PHIZL1t  )) THEN
            DEALLOCATE (I_PHIZL1t  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIZL1t  '
         IF (ALLOCATED(J_PHIZL1t  )) THEN
            DEALLOCATE (J_PHIZL1t  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIZL1t  '
         IF (ALLOCATED(PHIZL1t  )) THEN
            DEALLOCATE (PHIZL1t  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIZL2  ') THEN                 ! Deallocate arrays for PHIZL2

         NAME = 'I_PHIZL2  '
         IF (ALLOCATED(I_PHIZL2  )) THEN
            DEALLOCATE (I_PHIZL2  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIZL2  '
         IF (ALLOCATED(J_PHIZL2  )) THEN
            DEALLOCATE (J_PHIZL2  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIZL2  '
         IF (ALLOCATED(PHIZL2  )) THEN
            DEALLOCATE (PHIZL2  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'IF_LTM  ') THEN                 ! Deallocate arrays for IF_LTM

         NAME = 'I_IF_LTM  '
         IF (ALLOCATED(I_IF_LTM  )) THEN
            DEALLOCATE (I_IF_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_IF_LTM  '
         IF (ALLOCATED(J_IF_LTM  )) THEN
            DEALLOCATE (J_IF_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'IF_LTM  '
         IF (ALLOCATED(IF_LTM  )) THEN
            DEALLOCATE (IF_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'IRR') THEN                      ! Deallocate arrays for IRR

         NAME = 'I_IRR'
         IF (ALLOCATED(I_IRR)) THEN
            DEALLOCATE (I_IRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_IRR'
         IF (ALLOCATED(J_IRR)) THEN
            DEALLOCATE (J_IRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'IRR'
         IF (ALLOCATED(IRR)) THEN
            DEALLOCATE (IRR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIXA') THEN                    ! Deallocate arrays for PHIXA

         NAME = 'I_PHIXA'
         IF (ALLOCATED(I_PHIXA)) THEN
            DEALLOCATE (I_PHIXA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIXA'
         IF (ALLOCATED(J_PHIXA)) THEN
            DEALLOCATE (J_PHIXA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIXA'
         IF (ALLOCATED(PHIXA)) THEN
            DEALLOCATE (PHIXA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PHIXG') THEN                    ! Deallocate arrays for PHIXG

         NAME = 'I_PHIXG'
         IF (ALLOCATED(I_PHIXG)) THEN
            DEALLOCATE (I_PHIXG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_PHIXG'
         IF (ALLOCATED(J_PHIXG)) THEN
            DEALLOCATE (J_PHIXG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PHIXG'
         IF (ALLOCATED(PHIXG)) THEN
            DEALLOCATE (PHIXG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRRcb') THEN                    ! Deallocate arrays for KRRcb

         NAME = 'I_KRRcb'
         IF (ALLOCATED(I_KRRcb)) THEN
            DEALLOCATE (I_KRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRRcb'
         IF (ALLOCATED(J_KRRcb)) THEN
            DEALLOCATE (J_KRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRRcb'
         IF (ALLOCATED(KRRcb)) THEN
            DEALLOCATE (KRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRRcbn') THEN                   ! Deallocate arrays for KRRcbn

         NAME = 'I_KRRcbn'
         IF (ALLOCATED(I_KRRcbn)) THEN
            DEALLOCATE (I_KRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRRcbn'
         IF (ALLOCATED(J_KRRcbn)) THEN
            DEALLOCATE (J_KRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRRcbn'
         IF (ALLOCATED(KRRcbn)) THEN
            DEALLOCATE (KRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KRRcbs') THEN                   ! Deallocate arrays for KRRcbs

         NAME = 'I_KRRcbs'
         IF (ALLOCATED(I_KRRcbs)) THEN
            DEALLOCATE (I_KRRcbs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KRRcbs'
         IF (ALLOCATED(J_KRRcbs)) THEN
            DEALLOCATE (J_KRRcbs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KRRcbs'
         IF (ALLOCATED(KRRcbs)) THEN
            DEALLOCATE (KRRcbs,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'KXX  ') THEN                    ! Deallocate arrays for KXX

         NAME = 'I_KXX  '
         IF (ALLOCATED(I_KXX  )) THEN
            DEALLOCATE (I_KXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_KXX  '
         IF (ALLOCATED(J_KXX  )) THEN
            DEALLOCATE (J_KXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KXX  '
         IF (ALLOCATED(KXX  )) THEN
            DEALLOCATE (KXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'LTM  ') THEN                    ! Deallocate arrays for LTM

         NAME = 'I_LTM  '
         IF (ALLOCATED(I_LTM  )) THEN
            DEALLOCATE (I_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_LTM  '
         IF (ALLOCATED(J_LTM  )) THEN
            DEALLOCATE (J_LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LTM  '
         IF (ALLOCATED(LTM  )) THEN
            DEALLOCATE (LTM  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MRN  ') THEN                    ! Deallocate arrays for MRN

         NAME = 'I_MRN  '
         IF (ALLOCATED(I_MRN  )) THEN
            DEALLOCATE (I_MRN  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MRN  '
         IF (ALLOCATED(J_MRN  )) THEN
            DEALLOCATE (J_MRN  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MRN  '
         IF (ALLOCATED(MRN  )) THEN
            DEALLOCATE (MRN  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MRRcb') THEN                    ! Deallocate arrays for MRRcb

         NAME = 'I_MRRcb'
         IF (ALLOCATED(I_MRRcb)) THEN
            DEALLOCATE (I_MRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MRRcb'
         IF (ALLOCATED(J_MRRcb)) THEN
            DEALLOCATE (J_MRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MRRcb'
         IF (ALLOCATED(MRRcb)) THEN
            DEALLOCATE (MRRcb,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MRRcbn') THEN                   ! Deallocate arrays for MRRcbn

         NAME = 'I_MRRcbn'
         IF (ALLOCATED(I_MRRcbn)) THEN
            DEALLOCATE (I_MRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MRRcbn'
         IF (ALLOCATED(J_MRRcbn)) THEN
            DEALLOCATE (J_MRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MRRcbn'
         IF (ALLOCATED(MRRcbn)) THEN
            DEALLOCATE (MRRcbn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MXX  ') THEN                    ! Deallocate arrays for MXX

         NAME = 'I_MXX  '
         IF (ALLOCATED(I_MXX  )) THEN
            DEALLOCATE (I_MXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MXX  '
         IF (ALLOCATED(J_MXX  )) THEN
            DEALLOCATE (J_MXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MXX  '
         IF (ALLOCATED(MXX  )) THEN
            DEALLOCATE (MXX  ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MXXn') THEN                     ! Deallocate arrays for MXXn

         NAME = 'I_MXXn'
         IF (ALLOCATED(I_MXXn)) THEN
            DEALLOCATE (I_MXXn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'J_MXXn'
         IF (ALLOCATED(J_MXXn)) THEN
            DEALLOCATE (J_MXXn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MXXn'
         IF (ALLOCATED(MXXn)) THEN
            DEALLOCATE (MXXn,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! **********************************************************************************************************************************

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME_IN
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_SPARSE_MAT

   END MODULE SPARSE_MATRIX_DEALLOCATION
