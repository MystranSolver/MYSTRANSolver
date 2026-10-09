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

   MODULE CB_EFFECTIVE_PROPERTIES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_CB_MEFM_MPF

   CONTAINS

      SUBROUTINE CALC_CB_MEFM_MPF

! Calculates the modal participation factors and modal mass via the technique described in my document MPF.doc:

!                                  MPF = EIGEN_VEC'(MLL*DLR + MRL)/GEN_MASS

! MPF has a number of rows equal to the number of eigenvectors and a number of columns equal to the number of L-set DOF's.
! In this subr the transpose of MPF (MPFt) is first calculated which is then transposed to get MPF, where MPFt is

!                                  MPFt = (MLL*DLR + MRL)'*EIGEN_VEC//GEN_MASS

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR,                                                                  &
                                         NDOFL, NDOFR, NTERM_MPF0 , NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTs_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  MPFOUT
      USE RIGID_BODY_DISP_MATS, ONLY  :  TR6_MEFM
      USE SPARSE_MATRICES, ONLY       :  I_MPF0 , J_MPF0 , MPF0 , SYM_MPF0
      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VEC, GEN_MASS, MEFFMASS, MPFACTOR_NR, MPFACTOR_N6

      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_TO_FULL
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_CB_MEFM_MPF'
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG)                   :: I,J,K                   ! DO loop indices


      REAL(DOUBLE)                    :: DUM1(NDOFR,6)           ! Intermediate matrix
      REAL(DOUBLE)                    :: MEFW_MAT_RR(NDOFR,NDOFR)! Modal eff wgt for 1 mode for all R DOF's
      REAL(DOUBLE)                    :: MEFW_MAT_66(6,6)        ! Modal eff wgt for 1 mode transformed to CG via TR6_CG
      REAL(DOUBLE)                    :: MEFW_DIAG_NR(NVEC,NDOFR)! Matrix whose i-th row is the diagonal from MEFW for mode i
      REAL(DOUBLE)                    :: MPFt(NDOFR,NVEC)        ! Transpose of MPF
      REAL(DOUBLE)                    :: MPFi(1,NDOFR)           ! i-th row of MPF
      REAL(DOUBLE)                    :: MPFit(NDOFR,1)          ! MPFi'



! **********************************************************************************************************************************
      CALL ALLOCATE_EIGEN1_MAT ( 'MEFFMASS', NVEC, 6    , SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'MPFACTOR_NR', NVEC, NDOFR, SUBR_NAME )

! Convert sparse NDOFRxNVEC matrix MPF0 = (MRL + DLRt*MLL)*EIGEN_VEC to full matrix. MPF0 is gen mass times transpose of MPFACTOR

      CALL SPARSE_CRS_TO_FULL  ( 'MPF0', NTERM_MPF0, NDOFR, NVEC, SYM_MPF0, I_MPF0, J_MPF0, MPF0, MPFt )

      DO I=1,NVEC                                          ! Scale by gen mass. Transpose is MPFACTOR matrix
         DO J=1,NDOFR
            MPFACTOR_NR(I,J) = MPFt(J,I)/GEN_MASS(I)
         ENDDO
      ENDDO

      IF (MPFOUT == '6') THEN
         CALL ALLOCATE_EIGEN1_MAT ( 'MPFACTOR_N6', NVEC, 6, SUBR_NAME )
         CALL MATMULT_FFF ( MPFACTOR_NR, TR6_MEFM, NVEC, NDOFR, 6, MPFACTOR_N6 )
      ENDIF

! Calculate modal effective mass from MPF's

      DO I=1,NVEC
                                                           ! MPF is an NVEC x NDOFR matrix of modal participation factors
         DO J=1,NDOFR
            MPFit(J,1) = MPFACTOR_NR(I,J)                  ! MPFit is a col vec equal to the I-th row of MPF (MPF's for mode i)
         ENDDO

         DO J=1,NDOFR
            MPFi(1,J)  = MPFACTOR_NR(I,J)                  ! MPFi is a row vec of the MPF's for mode i (MPFit = MPFi')
         ENDDO
                                                           ! MEFW_MAT_RR = MPFi'xMPFi ( a square matrix of NDOFR x NDOFR)
         CALL MATMULT_FFF ( MPFit, MPFi, NDOFR, 1, NDOFR, MEFW_MAT_RR )
         DO J=1,NDOFR
            DO K=1,NDOFR                                   ! Scale MEFW_MAT_RR by GEN_MASS
               MEFW_MAT_RR(J,K) = GEN_MASS(I)*MEFW_MAT_RR(J,K)
            ENDDO
         ENDDO

         DO J=1,NDOFR                                      ! MEFW_DIAG_NR are diagonals from MEFW_MAT_RR (NDOFR x NDOFR)
            MEFW_DIAG_NR(I,J) = MEFW_MAT_RR(J,J)
         ENDDO
                                                           ! Reduce NDOFR x NDOFR MEFW_MAT_RR to 6 x 6 MEFW_MAT_66
         CALL MATMULT_FFF   ( MEFW_MAT_RR, TR6_MEFM, NDOFR, NDOFR, 6, DUM1 )
         CALL MATMULT_FFF_T ( TR6_MEFM, DUM1, NDOFR, 6, 6, MEFW_MAT_66 )

         DO J=1,6                                          ! Get diags from MEFW_MAT_66 and put into the ith row of MEFFMASS
            MEFFMASS(I,J) = MEFW_MAT_66(J,J)
         ENDDO

      ENDDO



      RETURN

! **********************************************************************************************************************************


97865 format(32767(1es14.6))

! **********************************************************************************************************************************

      END SUBROUTINE CALC_CB_MEFM_MPF

   END MODULE CB_EFFECTIVE_PROPERTIES
