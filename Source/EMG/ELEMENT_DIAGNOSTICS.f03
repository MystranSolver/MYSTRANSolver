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

   MODULE ELEMENT_DIAGNOSTICS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELMOFF, ELMOUT, ELMTLB

   CONTAINS

      SUBROUTINE ELMOFF ( OPT, WRITE_WARN )

! Processes element mass, stiffness, thermal load, pressure load, stress recovery matrices if there are any offsets of the element
! at any grid points. This is a general routine which can be used by any of the elements as long as the
! element has no more than 4 grid points.
! ======================================

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, NSUB, NTSUB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  CAN_ELEM_TYPE_OFFSET, ELGP, EID, KE, ME, NUM_EMG_FATAL_ERRS,                &
                                         OFFDIS, OFFSET, PPE, PTE, SE1, SE2, SE3, TYPE
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMOFF'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: I,J,K,L,M,N       ! DO loop indices
      INTEGER(LONG)                   :: II,JJ             ! Computed indices
      INTEGER(LONG)                   :: JBEG              ! Index
      INTEGER(LONG)                   :: KBEG              ! Index
      INTEGER(LONG)                   :: ROW               ! A computed row number in the elem stiff matrix
      INTEGER(LONG)                   :: COL               ! A computed col number in the elem stiff matrix
      INTEGER(LONG)                   :: NCOL              ! An input to subr MULT_OFFSET, called herein
      INTEGER(LONG)                   :: METH              ! An input to subr MULT_OFFSET, called herein


      REAL(DOUBLE)                    :: DUM3(3,3)         ! An intermediate result when calculating offset SEi
      REAL(DOUBLE)                    :: DUM4(3,3)         ! An intermediate result when calculating offset SEi
      REAL(DOUBLE)                    :: DUM11(3,3)        ! An intermediate result when calculating offset KE
      REAL(DOUBLE)                    :: DUM12(3,3)        ! An intermediate result when calculating offset KE
      REAL(DOUBLE)                    :: DUM21(3,3)        ! An intermediate result when calculating offset KE
      REAL(DOUBLE)                    :: DUM22(3,3)        ! An intermediate result when calculating offset KE
      REAL(DOUBLE)                    :: DXI               ! An offset distance in direction 1
      REAL(DOUBLE)                    :: DYI               ! An offset distance in direction 2
      REAL(DOUBLE)                    :: DZI               ! An offset distance in direction 3
      REAL(DOUBLE)                    :: DXJ               ! An offset distance in direction 1
      REAL(DOUBLE)                    :: DYJ               ! An offset distance in direction 2
      REAL(DOUBLE)                    :: DZJ               ! An offset distance in direction 3
      REAL(DOUBLE)                    :: PDUM1(3,NSUB)     ! An intermediate result when calculating offset PTE, PPE
      REAL(DOUBLE)                    :: PDUM2(3,NSUB)     ! An intermediate result when calculating offset PTE, PPE

      REAL(DOUBLE)                    :: DUM_KE(6*ELGP,6*ELGP)
      REAL(DOUBLE)                    :: E(6*ELGP,6*ELGP)
      REAL(DOUBLE)                    :: Ei(ELGP,6,6)
      REAL(DOUBLE)                    :: KE1(6*ELGP,6*ELGP)


! **********************************************************************************************************************************
! Make sure we are not here for an element that does not support offsets

      IF (CAN_ELEM_TYPE_OFFSET /= 'Y') THEN
         WRITE(ERR,1955) SUBR_NAME, TYPE, EID
         WRITE(F06,1955) SUBR_NAME, TYPE, EID
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! The processing is done in 1 major loop over the  number of G.P.'s in which PTE, SEi are processed by pre or post multiplying by
! the offset matrix (or it's transpose).

! A minor loop within the major loop takes care of ME and KE which gets post-multiplied by the offset matrix and pre-multiplied by
! it's transpose. The offset matrices for each G.P. are a 6 x 6 matrix which is an identity matrix plus a small 3 x 3 submatrix
! containing only 3 independent terms, and the processing takes advantage of this and simplifies the matrix multiplications.
! The offset matrix is called E for each G.P. but is never written out as a 6 x 6 matrix.

! The general form of E for one grid point is:

!                             | 1  0  0 |  0    DZ  -DY |
!                             | 0  1  0 | -DZ   0    DX |
!                             | 0  0  1 |  DY  -DX   0  |
!                         E = |---------|---------------|
!                             | 0  0  0 |  1    0    0  |
!                             | 0  0  0 |  0    1    0  |
!                             | 0  0  0 |  0    0    1  |

! where DX, DY and DZ are the 3 components of the offset of the element at a grid and are in global coords

! With this E matrix, the transformed element matrices are (prime indicates matrix transposition):

!                                MEg = E'* MEe * E

!                                KEg = E'* KEe * E

!                               PTEg = E'* PTEe

!                                SEg = SEe * E

! where MEe, KEe, PTEe and SEe are the mass, stiffness, thermal loads and stress recovery matrices developed in local element
! coordinates at the element nodes and MEg, KEg, PTEg, SEg are the same matrices but in terms of degrees of freedom at the grids and
! in global coordinates.

! Initialize

      E  = ZERO
      Ei = ZERO

      DO I=1,ELGP

         II  = 6*(I-1)
         DXI = ZERO
         DYI = ZERO
         DZI = ZERO
         IF (OFFSET(I) == 'Y') THEN
            DXI = OFFDIS(I,1)
            DYI = OFFDIS(I,2)
            DZI = OFFDIS(I,3)
         ENDIF
         DO J=1,6
            Ei(I,J,J) = ONE
         ENDDO

         Ei(I,1,5) =  DZI
         Ei(I,1,6) = -DYI
         Ei(I,2,4) = -DZI
         Ei(I,2,6) =  DXI
         Ei(I,3,4) =  DYI
         Ei(I,3,5) = -DXI

      ENDDO

! Set E matrix

      DO I=1,ELGP
         JBEG = 6*(I-1)
         DO J=1,6
            KBEG = 6*(I-1)
            DO K=1,6
               E(JBEG+J,KBEG+K) = Ei(I,J,K)
            ENDDO
         ENDDO
      ENDDO

! Apply offset to stiffness matrix KE

      IF (OPT(4) == 'Y') THEN
         KE1 = KE(:6*ELGP, :6*ELGP)
         DUM_KE = MATMUL(KE1,E)
         KE1 = MATMUL(TRANSPOSE(E),DUM_KE)
         KE(:6*ELGP, :6*ELGP) = KE1
      ENDIF

! Process offsets

      DO I=1,ELGP
         II  = 6*(I-1)
         DXI = ZERO
         DYI = ZERO
         DZI = ZERO
         IF (OFFSET(I) == 'Y') THEN
            DXI = OFFDIS(I,1)
            DYI = OFFDIS(I,2)
            DZI = OFFDIS(I,3)
            IF (OPT(2) == 'Y') THEN                        ! Process PTE. Generate E'* PTE
               DO J=1,3
                  DO K=1,NTSUB
                     PDUM1(J,K) = PTE(II+J,K)
                  ENDDO
               ENDDO
               NCOL = NTSUB
               METH = 2
               CALL MULT_OFFSET ( PDUM1, DXI, DYI, DZI, NCOL, METH, PDUM2 )
               DO J=1,3
                  DO K=1,NTSUB
                     PTE(II+J+3,K) = PTE(II+J+3,K) + PDUM2(J,K)
                  ENDDO
               ENDDO
            ENDIF

            IF (OPT(1) == 'Y') THEN                        ! Process ME. Generate E(transp)*ME*E.
               DO J=I,ELGP
                  JJ = 6*(J-1)
                  IF (OFFSET(J) == 'Y') THEN
                     DXJ = OFFDIS(J,1)
                     DYJ = OFFDIS(J,2)
                     DZJ = OFFDIS(J,3)

                     DO K=1,3                              ! Partition ME for this grid point pair (i,j) into 4-3x3 matrices
                        DO L=1,3
                           DUM11(K,L) = ME(II+K,JJ+L)
                           DUM12(K,L) = ME(II+K,JJ+L+3)
                           DUM21(K,L) = ME(II+K+3,JJ+L)
                           DUM22(K,L) = ME(II+K+3,JJ+L+3)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify upper right 3x3 partition of ME
                     METH = 1
                     CALL MULT_OFFSET ( DUM11, DXJ, DYJ, DZJ, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           ME(II+K,JJ+L+3) = DUM12(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify lower left 3x3 partition of ME
                     METH = 2
                     CALL MULT_OFFSET ( DUM11, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           ME(II+K+3,JJ+L) = DUM21(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify lower right 3x3 partition of ME
                     METH = 1
                     CALL MULT_OFFSET ( DUM21, DXJ, DYJ, DZJ, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           ME(II+K+3,JJ+L+3) = DUM22(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3
                     METH = 2
                     CALL MULT_OFFSET ( DUM12, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           ME(II+K+3,JJ+L+3) = ME(II+K+3,JJ+L+3) +DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3
                     METH = 1
                     CALL MULT_OFFSET ( DUM11, DXJ, DYJ, DZJ, NCOL, METH, DUM4 )
                     NCOL = 3
                     METH = 2
                     CALL MULT_OFFSET ( DUM4, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           ME(II+K+3,JJ+L+3) = ME(II+K+3,JJ+L+3) +DUM3(K,L)
                        ENDDO
                     ENDDO

                  ENDIF

               ENDDO

               DO K=2,ELGP                                 ! Generate the remaining Mij using symmetry.
                  DO L=1,K-1
                     DO M=1,6
                        DO N=1,6
                           ROW = 6*(K-1) + M
                           COL = 6*(L-1) + N
                           ME(ROW,COL) = ME(COL,ROW)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

            IF (OPT(3) == 'Y') THEN                        ! Process SEi. Generate SEi*E

               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                DUM3(J,K) = SE1(L,J,II+K)
                        DUM3(J,K) = SE1(J,II+K,L)
                     ENDDO
                  ENDDO
               ENDDO
               NCOL = 3
               METH = 1
               CALL MULT_OFFSET ( DUM3, DXI, DYI, DZI, NCOL, METH, DUM4 )
               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                SE1(L,J,II+K+3) = SE1(L,J,II+K+3) + DUM4(J,K)
                        SE1(J,II+K+3,L) = SE1(J,II+K+3,L) + DUM4(J,K)
                     ENDDO
                  ENDDO
               ENDDO

               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                DUM3(J,K) = SE2(L,J,II+K)
                        DUM3(J,K) = SE2(J,II+K,L)
                     ENDDO
                  ENDDO
               ENDDO
               NCOL = 3
               METH = 1
               CALL MULT_OFFSET ( DUM3, DXI, DYI, DZI, NCOL, METH, DUM4 )
               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                SE2(L,J,II+K+3) = SE2(L,J,II+K+3) + DUM4(J,K)
                        SE2(J,II+K+3,L) = SE2(J,II+K+3,L) + DUM4(J,K)
                     ENDDO
                  ENDDO
               ENDDO

               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                DUM3(J,K) = SE3(L,J,II+K)
                        DUM3(J,K) = SE3(J,II+K,L)
                     ENDDO
                  ENDDO
               ENDDO
               NCOL = 3
               METH = 1
               CALL MULT_OFFSET ( DUM3, DXI, DYI, DZI, NCOL, METH, DUM4 )
               DO L=1,MAX_STRESS_POINTS+1
                  DO J=1,3
                     DO K=1,3
!xxError                SE3(L,J,II+K+3) = SE3(L,J,II+K+3) + DUM4(J,K)
                        SE3(J,II+K+3,L) = SE3(J,II+K+3,L) + DUM4(J,K)
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

            IF (OPT(4) == 'Z') THEN                        ! Process KE. Generate E(transp)*KE*E.
!xx         IF (OPT(4) == 'Y') THEN                        ! Process KE. Generate E(transp)*KE*E.
               DO J=I,ELGP
                  JJ = 6*(J-1)
                  IF (OFFSET(J) == 'Y') THEN
                     DXJ = OFFDIS(J,1)
                     DYJ = OFFDIS(J,2)
                     DZJ = OFFDIS(J,3)

                     DO K=1,3                              ! Partition KE for this grid point pair (i,j) into 4-3x3 matrices
                        DO L=1,3
                           DUM11(K,L) = KE(II+K,JJ+L)
                           DUM12(K,L) = KE(II+K,JJ+L+3)
                           DUM21(K,L) = KE(II+K+3,JJ+L)
                           DUM22(K,L) = KE(II+K+3,JJ+L+3)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify upper right 3x3 partition of KE
                     METH = 1
                     CALL MULT_OFFSET ( DUM11, DXJ, DYJ, DZJ, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           KE(II+K,JJ+L+3) = DUM12(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify lower left 3x3 partition of KE
                     METH = 2
                     CALL MULT_OFFSET ( DUM11, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           KE(II+K+3,JJ+L) = DUM21(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3                              ! Modify lower right 3x3 partition of KE
                     METH = 1
                     CALL MULT_OFFSET ( DUM21, DXJ, DYJ, DZJ, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           KE(II+K+3,JJ+L+3) = DUM22(K,L) + DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3
                     METH = 2
                     CALL MULT_OFFSET ( DUM12, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           KE(II+K+3,JJ+L+3) = KE(II+K+3,JJ+L+3) +DUM3(K,L)
                        ENDDO
                     ENDDO

                     NCOL = 3
                     METH = 1
                     CALL MULT_OFFSET ( DUM11, DXJ, DYJ, DZJ, NCOL, METH, DUM4 )
                     NCOL = 3
                     METH = 2
                     CALL MULT_OFFSET ( DUM4, DXI, DYI, DZI, NCOL, METH, DUM3 )
                     DO K=1,3
                        DO L=1,3
                           KE(II+K+3,JJ+L+3) = KE(II+K+3,JJ+L+3) +DUM3(K,L)
                        ENDDO
                     ENDDO

                  ENDIF

               ENDDO

               DO K=2,ELGP                                 ! Generate the remaining Kij using symmetry.
                  DO L=1,K-1
                     DO M=1,6
                        DO N=1,6
                           ROW = 6*(K-1) + M
                           COL = 6*(L-1) + N
                           KE(ROW,COL) = KE(COL,ROW)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

            IF (OPT(5) == 'Y') THEN                        ! Process PPE. Generate E(transp.)*PPE
               DO J=1,3
                  DO K=1,NSUB
                     PDUM1(J,K) = PPE(II+J,K)
                  ENDDO
               ENDDO
               NCOL = NSUB
               METH = 2
               CALL MULT_OFFSET ( PDUM1, DXI, DYI, DZI, NCOL, METH, PDUM2 )
               DO J=1,3
                  DO K=1,NSUB
                     PPE(II+J+3,K) = PPE(II+J+3,K) + PDUM2(J,K)
                  ENDDO
               ENDDO

            ENDIF

         ENDIF

      ENDDO





      RETURN

! **********************************************************************************************************************************
 3006 format(6(1es14.6))

 1955 FORMAT(' *ERROR  1955: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE ',A,' DOES NOT SUPPORT OFFSETS. ERROR OCCURRED FOR ELEMENT NUMBER ',I8)

1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,', HAS ZERO OR NEGATIVE ',A,' = ',1ES9.1)

1948 FORMAT(' *ERROR  1948: ',A,I8,' MUST HAVE INTEGRATION ORDERS FOR PARAMS ',A,' = ',I3,' IF THE ELEMENT IS A PCOMP'            &
                             ,/,14X,' WITH SYM LAYUP. HOWEVER, THE TWO INTEGRATION ORDERS WERE: ',A,' = ',I3,' AND ',A,' = ',I3)

98761 FORMAT('   X1-X2 = ',1ES14.6,'   Y1-Y2 = ',1ES14.6)

98762 FORMAT('   X2-X3 = ',1ES14.6,'   Y2-Y3 = ',1ES14.6)

98763 FORMAT('   X3-X4 = ',1ES14.6,'   Y3-Y4 = ',1ES14.6)

98764 FORMAT('   X4-X1 = ',1ES14.6,'   Y4-Y1 = ',1ES14.6)
! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MULT_OFFSET ( A, DX, DY, DZ, NCOLA, METH, B )

! Perform matrix multiply to get A*E or E(transp)*A for elem offsets. Matrix E, the offset matrix, is a simple form.
! It is an identity  6 x 6 plus a 3 x 3 in the upper right corner containing the 3 offset distances. Due to this
! simplicity, A*E or E(transp)*A is calculated explicitly

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEFE
      USE MODEL_STUF, ONLY            :  EMG_IFE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MULT_OFFSET'

      INTEGER(LONG)                   :: IK,JK             ! DO loop indices
      INTEGER(LONG), INTENT(IN)       :: METH              ! = 1 if A*E is to be calculated
                                                           ! = 2 if E(transp)*A is to be calculated
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! Number of cols in matrix A


      REAL(DOUBLE) , INTENT(IN)       :: A(3,NCOLA)        ! Matrix to either post-multiply E by or pre-multiply E(transp) by
      REAL(DOUBLE) , INTENT(IN)       :: DX                ! Offset distance in direction 1
      REAL(DOUBLE) , INTENT(IN)       :: DY                ! Offset distance in direction 2
      REAL(DOUBLE) , INTENT(IN)       :: DZ                ! Offset distance in direction 3
      REAL(DOUBLE) , INTENT(INOUT)    :: B(3,NCOLA)        ! Result matrix of either A*E or E(transp)*A



! **********************************************************************************************************************************
! Do not initialize B. It is an output in one call and maybe is input back as A in next call

      IF      (METH == 1) THEN
         DO IK=1,3
            B(IK,1) = -A(IK,2)*DZ + A(IK,3)*DY
            B(IK,2) =  A(IK,1)*DZ - A(IK,3)*DX
            B(IK,3) = -A(IK,1)*DY + A(IK,2)*DX
         ENDDO
      ELSE IF (METH == 2) THEN
         DO JK=1,NCOLA
            B(1,JK) = -A(2,JK)*DZ + A(3,JK)*DY
            B(2,JK) =  A(1,JK)*DZ - A(3,JK)*DX
            B(3,JK) = -A(1,JK)*DY + A(2,JK)*DX
         ENDDO
      ELSE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR /= 0) THEN
            WRITE(ERR,1917) SUBR_NAME, METH
            WRITE(F06,1917) SUBR_NAME, METH
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1917
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = METH
            ENDIF
         ENDIF
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1917 FORMAT(' *ERROR  1917: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' METHOD INDICATOR MUST BE 1 OR 2. VALUE IS ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE MULT_OFFSET

      END SUBROUTINE ELMOFF


      SUBROUTINE ELMOUT ( INT_ELEM_ID, DUM_BUG, CASE_NUM, OPT )

! Prints elem related data (controlled by Case Control ELDATA requests and situational variable WRT_BUG(i) ).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_DAT1_BIT, ELDT_BUG_DAT2_BIT, ELDT_BUG_ME_BIT, ELDT_BUG_P_T_BIT,    &
                                         ELDT_BUG_SE_BIT, ELDT_BUG_KE_BIT, ELDT_BUG_U_P_BIT, MBUG, MDT, MELGP, METYPE,             &
                                         MEMATR, MEMATC, MEPROP, MPRESS, NSUB, NTSUB, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  CONV_RAD_DEG, ZERO
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, ELFORCEN, QUADAXIS, QUAD4TYP
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, BE1, BE2, BE3, BENSUM, BMEANT, CAN_ELEM_TYPE_OFFSET, DOFPIN, DT, ELAS_COMP, &
                                         EID, EB, EM, ES, ET, ELEM_LEN_AB, ELDOF, ELMTYP, ELGP, EMAT, EPROP, FCONV, HBAR, KE, KED, &
                                         ME, MXWARP, NUM_PLIES, NUM_SEi, OFFDIS, OFFSET, PCOMP_PROPS, PEB, PEG, PEL, PHI_SQ,       &
                                         PPE, PRESS, PSI_HAT, PTE, QUAD_DELTA, QUAD_GAMMA, QUAD_THETA, SE1, SE2, SE3,              &
                                         SHELL_T, SHRSUM, STE1, STE2, STE3, THETAM, TE, TYPE, UEB, UEG, UEL, XEB, XEL, SCNUM,      &
                                         SUBLOD, ULT_STRE, ULT_STRN
      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMOUT'
      CHARACTER(LEN=*), INTENT(IN)    :: OPT(6)              ! Array of EMG option indicators explained above
      CHARACTER( 1*BYTE)              :: FOUND               ! Used in determining if we found something we were looking for
      CHARACTER(12*BYTE)              :: GRID_TYPE(ELGP)     ! Type of grid: scalar point or act grid point
      CHARACTER(60*BYTE)              :: NAME1               ! Text used for output print purposes
      CHARACTER(21*BYTE)              :: NAME2               ! Text used for output print purposes
      CHARACTER(12*BYTE)              :: NAME3               ! Text used for output print purposes

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID         ! Internal element ID for which
      INTEGER(LONG), INTENT(IN)       :: CASE_NUM            ! Can be subcase number (e.g. for UEL, PEL output)
      INTEGER(LONG), INTENT(IN)       :: DUM_BUG(0:MBUG-1)   ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: I2                  ! Counter
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices
      INTEGER(LONG)                   :: NUM_COMPS           ! No. displ components (1 for SPOINT, 6 for actual grid)
      INTEGER(LONG)                   :: NUM_STRESS_MATS     ! Number of SEi/BEi matrices for this element
      INTEGER(LONG)                   :: TCASE2(NSUB)        ! TCASE2(I) gives the internal subcase no. for internal thermal case I
!                                                              If there are 5 subcases and internal S/C 3 is the 1-st S/C to have
!                                                              thermal load and internal S/C 5 is the 2-nd to have thermal load:
!                                                              TCASE2(1-5) = 3, 5, 0, 0, 0


      REAL(DOUBLE)                    :: OEL(6)              ! Temp array for holding elem displ, node loads
      REAL(DOUBLE)                    :: SHELL_T_avg         ! Average of the diag terms from transverse shear matrix SHELL_T

      INTRINSIC                       :: ABS


! **********************************************************************************************************************************
! Set GRID_TYPE

      DO I=1,ELGP
         GRID_TYPE(I) = 'undefined   '
      ENDDO
      DO I=1,ELGP
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         IF      (NUM_COMPS == 1) THEN
            GRID_TYPE(I) = 'scalar point'
         ELSE IF (NUM_COMPS == 6) THEN
            GRID_TYPE(I) = 'grid   point'
         ENDIF
      ENDDO

! Write output from ELMDAT subroutine.

      IF (DUM_BUG(0) > 0) THEN

         WRITE(BUG,1000)

         WRITE(BUG,1001) ELDT_BUG_DAT1_BIT, TYPE, EID
         WRITE(BUG,*)

         IF ((TYPE == 'QDPLT2   ') .OR. (TYPE == 'QUAD4   ')) THEN
            WRITE(BUG,*) '  Bending portion of QUAD4 is based on QUAD4TYP formulation = ',QUAD4TYP
            WRITE(BUG,*)
         ENDIF

         WRITE(BUG,*) '  Internal element number,       INT_ELEM_ID  = ' ,INT_ELEM_ID
         WRITE(BUG,*) '  Number of grids elem is connected to, ELGP  = ' ,ELGP
         WRITE(BUG,*) '  Number of DOFs  elem is connected to, ELDOF = ' ,ELDOF
         WRITE(BUG,*)

         WRITE(BUG,*) '  Actual & Internal G.P.s and basic coordinates (AGRID, BGRID, I, XEB)'
         WRITE(BUG,*) '  --------------------------------------------------------------------'
         WRITE(BUG,*) '         Grid Number                   Basic System Coordinates'
         WRITE(BUG,*) '  Actual  Internal  Element        X              Y              Z'
         DO I=1,ELGP
            WRITE(BUG,1004)  AGRID(I), BGRID(I),I,(XEB(I,J),J=1,3)
         ENDDO
         IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'USER1   ')) THEN
            WRITE(BUG,1113) (XEB(ELGP+1,J),J=1,3)
         ENDIF
         WRITE(BUG,*)

         IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'USER1   ')) THEN
            WRITE(BUG,*) '  Basic coordinate directions of v vector'
            WRITE(BUG,*) '  ---------------------------------------'
            WRITE(BUG,3003) (XEB(ELGP+1,J),J=1,3)
            WRITE(BUG,*)
         ENDIF

      ENDIF

      IF (DUM_BUG(1) > 0) THEN

!zzzz    IF (PCOMP_PROPS == 'N') THEN

            WRITE(BUG,1000)

            WRITE(BUG,1001) ELDT_BUG_DAT2_BIT, TYPE, EID

            WRITE(BUG,*) '  EPROP, Array of element property data'
            WRITE(BUG,*) '  -------------------------------------'
            WRITE(BUG,3005) (EPROP(I),I=1,MEPROP)
            WRITE(BUG,*)

            WRITE(BUG,*) '  EMAT, Array of element matl data'
            WRITE(BUG,*) '  --------------------------------'
            WRITE(BUG,*) '    Col 1         Col 2         Col 3         Col 4'
            DO I=1,MEMATR
               WRITE(BUG,3004) (EMAT(I,J),J=1,MEMATC)
            ENDDO
            WRITE(BUG,*)

            WRITE(BUG,*) '  ULT_STRE, material stress allowables'
            WRITE(BUG,*) '  ------------------------------------'
            WRITE(BUG,*) '    Col 1         Col 2         Col 3         Col 4'
            WRITE(BUG,3024) (ULT_STRE(1,J),J=1,MEMATC),'   Tension ultimate'
            WRITE(BUG,3024) (ULT_STRE(2,J),J=1,MEMATC),'   Compr   ultimate'
            WRITE(BUG,3024) (ULT_STRE(3,J),J=1,MEMATC),'   Shear   ultimate'
            WRITE(BUG,*)

            WRITE(BUG,*) '  ULT_STRN, material strain allowables'
            WRITE(BUG,*) '  ------------------------------------'
            WRITE(BUG,*) '    Col 1         Col 2         Col 3         Col 4'
            WRITE(BUG,3024) (ULT_STRN(1,J),J=1,MEMATC),'   Tension ultimate'
            WRITE(BUG,3024) (ULT_STRN(2,J),J=1,MEMATC),'   Compr   ultimate'
            WRITE(BUG,3024) (ULT_STRN(3,J),J=1,MEMATC),'   Shear   ultimate'
            WRITE(BUG,*)

            FOUND = 'N'
            DO I=1,ELDOF
               IF (DOFPIN(I) /= 0) THEN
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'Y') THEN
               WRITE(BUG,*) '  DOFPIN, Array of element DOFs pinned'
               WRITE(BUG,*) '  ------------------------------------'
               WRITE(BUG,2024) (DOFPIN(I),I=1,ELDOF)
               WRITE(BUG,*)
            ENDIF

            IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
               WRITE(BUG,*) '  OFFDIS, Array of element offsets'
               WRITE(BUG,*) '  --------------------------------'
               DO I=1,ELGP
                  WRITE(BUG,3003) (OFFDIS(I,J),J=1,3)
               ENDDO
               WRITE(BUG,*)
            ENDIF

            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
               WRITE(BUG,5001) ELEM_LEN_AB
            ENDIF

            IF(TYPE(1:4) == 'ELAS') THEN                   ! For ELAS, write displ comps, at elem grids, that elem connects to
               WRITE(BUG,4001) (ELAS_COMP(I), GRID_TYPE(I), AGRID(I), I=1,2)
            ENDIF

! Write output from ELMGMi

            IF     ((TYPE == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR.                                  &
                    (TYPE == 'TRMEM   ') .OR. (TYPE == 'TRPLT1  ') .OR. (TYPE == 'TRPLT2  ') .OR.                                  &
                    (TYPE == 'TRIA3K  ') .OR. (TYPE == 'TRIA3   ')) THEN
                  WRITE(BUG,*) '  TE coord transformation matrix from subr ELMGM1'

            ELSE IF ((TYPE == 'QDMEM   ') .OR. (TYPE == 'QDPLT1  ') .OR. (TYPE == 'QDPLT2   ') .OR.                                &
                     (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD4   ')) THEN
                  WRITE(BUG,*) '  TE coord transformation matrix from subr ELMGM2 with QUADAXIS = ',QUADAXIS

            ELSE IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ') .OR.                                                           &
                     (TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                           &
                     (TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN
                  WRITE(BUG,*) '  TE coord transformation matrix from subr ELMGM3'

            ENDIF
            WRITE(BUG,*) '  ----------------------------------------------------------------------'

            IF ((TYPE(1:1) == 'Q') .OR. (TYPE(1:1) == 'H')) THEN
               IF (QUADAXIS == 'SPLITD') THEN
                  WRITE(BUG,*) '  (UEL = TE x UEB and local x axis splits angle between 2 diagonals)'
               ELSE
                  WRITE(BUG,*) '  (UEL = TE x UEB and local x axis is along side 1-2 of the element)'
               ENDIF
            ELSE
               WRITE(BUG,*) '  (UEL = TE x UEB and local x axis is along side 1-2 of the element)'
            ENDIF
            DO I=1,3
               WRITE(BUG,3003) (TE(I,J),J=1,3)
            ENDDO
            WRITE(BUG,*)

            WRITE(BUG,*) '  Actual & Internal G.P.s and local coordinates (AGRID, BGRID, I, XEL)'
            WRITE(BUG,*) '  --------------------------------------------------------------------'
            WRITE(BUG,*) '    Grid Number     Element          Local Element System Coords'
            WRITE(BUG,*) '  Actual Internal     No.          X              Y              Z'
            DO I=1,ELGP
               WRITE(BUG,1004)  AGRID(I), BGRID(I),I,(XEL(I,J),J=1,3)
            ENDDO
            WRITE(BUG,*)

            IF ((TYPE == 'QDMEM   ') .OR. (TYPE == 'QDPLT1  ') .OR. (TYPE == 'QDPLT2   ') .OR.                                     &
                (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD4   ')) THEN
               WRITE(BUG,5002) HBAR, MXWARP
               WRITE(BUG,5003) CONV_RAD_DEG*QUAD_THETA
               WRITE(BUG,5004) CONV_RAD_DEG*QUAD_GAMMA
               WRITE(BUG,5005) CONV_RAD_DEG*QUAD_DELTA
               WRITE(BUG,*)
               IF (ABS(HBAR) > MXWARP) THEN
                  WRITE(BUG,*) '  Matrix BMEAN correction to strain displ matrix used for this elem since HBAR > MXWARP'
                  WRITE(BUG,*) '  Rows of BMEAN (transpose of BMEANT):'
                  WRITE(BUG,*) '  -----------------------------------'
                  DO J=1,4
                     DO K=1,3
                        WRITE(BUG,1005) (BMEANT(I,3*(J-1)+K),I=1,8)
                     ENDDO
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDIF

            ENDIF

            WRITE(BUG,5006) CONV_RAD_DEG*THETAM               ! Write material orientation axes angle from elem x-axis

! Write material matrices for all but 1D elements

            IF((TYPE == 'TRMEM   ') .OR. (TYPE == 'TRPLT1  ') .OR. (TYPE == 'TRPLT2  ') .OR.                                       &
               (TYPE == 'TRIA3K  ') .OR. (TYPE == 'TRIA3   ') .OR.                                                                 &
               (TYPE == 'QDMEM   ') .OR. (TYPE == 'QDPLT1  ') .OR. (TYPE == 'QDPLT2   ') .OR.                                      &
               (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD4   ')) THEN

               WRITE(BUG,*) '  EM material matrix, in local element coordinate system, for membrane stresses'
               WRITE(BUG,*) '  -----------------------------------------------------------------------------'
               DO I=1,3
                  WRITE(BUG,3003) (EM(I,J),J=1,3)
               ENDDO
               WRITE(BUG,*)

               WRITE(BUG,*) '  EB material matrix, in local element coordinate system, for bending stresses'
               WRITE(BUG,*) '  ----------------------------------------------------------------------------'
               DO I=1,3
                  WRITE(BUG,3003) (EB(I,J),J=1,3)
               ENDDO
               WRITE(BUG,*)

               WRITE(BUG,*) '  ET material matrix, in local element coordinate system, for transverse shear stresses'
               WRITE(BUG,*) '  -------------------------------------------------------------------------------------'
               DO I=1,2
                  WRITE(BUG,3002) (ET(I,J),J=1,2)
               ENDDO
               WRITE(BUG,*)

            ELSE IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ') .OR.                                                           &
                     (TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                           &
                     (TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN

               WRITE(BUG,*) '  ES material matrix, in local element coordinate system, for 3D solid elements'
               WRITE(BUG,*) '  -----------------------------------------------------------------------------'
               DO I=1,6
                  WRITE(BUG,3006) (ES(I,J),J=1,6)
               ENDDO
               WRITE(BUG,*)

            ENDIF

!zzzz    ENDIF

      ENDIF

! **********************************************************************************************************************************
! Element thermal and pressure load matrices

      IF (DUM_BUG(2) > 0) THEN

         WRITE(BUG,1000)

         WRITE(BUG,1001) ELDT_BUG_P_T_BIT, TYPE, EID

         IF (NTSUB > 0) THEN

            DO I = 1,NSUB
               TCASE2(I) = 0
            ENDDO
            J = 0
            DO I = 1,NSUB
               IF (SUBLOD(I,2) /= 0) THEN
                  J = J + 1
                  TCASE2(J) = I
               ENDIF
            ENDDO

            WRITE(BUG,*) '  DT, Array of element temperature data (one row for each thermal subcase):'
            WRITE(BUG,*) '  ------------------------------------------------------------------------'
            DO J=1,NTSUB
               WRITE(BUG,1006) J, NTSUB
               WRITE(BUG,3026) (DT(I,J),I=1,MDT)
               WRITE(BUG,*)
            ENDDO
            WRITE(BUG,*)

            DO J=1,NTSUB
               WRITE(BUG,2001) SCNUM(TCASE2(J))
               IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
                  WRITE(BUG,1007)
               ELSE
                  IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                     WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
                  ENDIF
               ENDIF
               WRITE(BUG,3006) (PTE(I,J),I=1,ELDOF)
               WRITE(BUG,*)
            ENDDO
         ELSE
            WRITE(BUG,*) '  There is no thermal load data'
         ENDIF
         WRITE(BUG,*)

         WRITE(BUG,*) '  PRESS, Array of element pressure data (one row for each subcase):'
         WRITE(BUG,*) '  ----------------------------------------------------------------'
         DO J=1,NSUB
            WRITE(BUG,1006) J, NSUB
            WRITE(BUG,3026) (PRESS(I,J),I=1,MPRESS)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

         DO J=1,NSUB
            WRITE(BUG,2003) SCNUM(J)
            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
               WRITE(BUG,1007)
            ELSE
               IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                  WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
               ENDIF
            ENDIF
            WRITE(BUG,3006) (PPE(I,J),I=1,ELDOF)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

! **********************************************************************************************************************************
! Element mass

      IF (DUM_BUG(3) > 0) THEN

         WRITE(BUG,1000)

         WRITE(BUG,1001) ELDT_BUG_ME_BIT, TYPE, EID

         WRITE(BUG,*) '  ME element mass matrix in local element coordinate system'
         WRITE(BUG,*) '  ---------------------------------------------------------'
         IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
            WRITE(BUG,1007)
         ELSE
            IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
               WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
            ENDIF
         ENDIF
         DO I=1,ELDOF
            WRITE(BUG,*) '  Row',I
            WRITE(BUG,3006) (ME(I,J),J=1,ELDOF)
            WRITE(BUG,*)
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Element stiffness matrix

      IF (DUM_BUG(4) > 0) THEN

         WRITE(BUG,1000)

         WRITE(BUG,1001) ELDT_BUG_KE_BIT, TYPE, EID

         IF (TYPE(1:5) == 'TRIA3') THEN
            SHELL_T_avg = 0.5*( SHELL_T(1,1) + SHELL_T(2,2) )
            WRITE(BUG,*) '  TRIA3 plate element parameters used in calculating transverse shear stiffness'
            WRITE(BUG,*) '  -----------------------------------------------------------------------------'
            NAME1 = '      CBMIN3                                             = ' ; WRITE(BUG,5007) NAME1, CBMIN3
            NAME1 = '      SHELL_T_avg = 0.5*[(SHELL_T(1,1) + SHELL_T(2,2)]   = ' ; WRITE(BUG,5007) NAME1, SHELL_T_avg
            NAME1 = '      BENSUM                                             = ' ; WRITE(BUG,5007) NAME1, BENSUM
            NAME1 = '      SHRSUM                                             = ' ; WRITE(BUG,5007) NAME1, SHRSUM
            NAME1 = '      PSI_HAT = BENSUM/SHRSUM                            = ' ; WRITE(BUG,5007) NAME1, PSI_HAT
            NAME1 = '      PHI_SQ  = CBMIN3*PSI_HAT/(1 + CBMIN3*PSI_HAT)      = ' ; WRITE(BUG,5007) NAME1, PHI_SQ
            NAME1 = '      SHELL_T_avg*PHI_SQ                                 = ' ; WRITE(BUG,5007) NAME1, SHELL_T_avg*PHI_SQ
            WRITE(BUG,*)
         ENDIF

         IF (TYPE(1:5) == 'QUAD4') THEN
            SHELL_T_avg = 0.5*( SHELL_T(1,1) + SHELL_T(2,2) )
            WRITE(BUG,*) '  QUAD4 plate element parameters used in calculating transverse shear stiffness'
            WRITE(BUG,*) '  -----------------------------------------------------------------------------'
            NAME1 = '      CBMIN4                                             = ' ; WRITE(BUG,5007) NAME1, CBMIN4
            NAME1 = '      SHELL_T_avg = 0.5*[(SHELL_T(1,1) + SHELL_T(2,2)]   = ' ; WRITE(BUG,5007) NAME1, SHELL_T_avg
            NAME1 = '      BENSUM                                             = ' ; WRITE(BUG,5007) NAME1, BENSUM
            NAME1 = '      SHRSUM                                             = ' ; WRITE(BUG,5007) NAME1, SHRSUM
            NAME1 = '      PSI_HAT = BENSUM/SHRSUM                            = ' ; WRITE(BUG,5007) NAME1, PSI_HAT
            NAME1 = '      PHI_SQ  = CBMIN4*PSI_HAT/(1 + CBMIN4*PSI_HAT)      = ' ; WRITE(BUG,5007) NAME1, PHI_SQ
            NAME1 = '      SHELL_T_avg*PHI_SQ                                 = ' ; WRITE(BUG,5007) NAME1, SHELL_T_avg*PHI_SQ
            WRITE(BUG,*)
         ENDIF

         IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            WRITE(BUG,*) '  KED element stiffness matrix in local element coordinate system'
            WRITE(BUG,*) '  ---------------------------------------------------------------'
            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
               WRITE(BUG,1007)
            ELSE
               IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                  WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
               ENDIF
            ENDIF
            DO I=1,ELDOF
               WRITE(BUG,*) '  Row',I
               WRITE(BUG,3006) (KED(I,J),J=1,ELDOF)
               WRITE(BUG,*)
            ENDDO
         ELSE
            WRITE(BUG,*) '  KE element stiffness matrix in local element coordinate system'
            WRITE(BUG,*) '  --------------------------------------------------------------'
            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
               WRITE(BUG,1007)
            ELSE
               IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                  WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
               ENDIF
            ENDIF
            DO I=1,ELDOF
               WRITE(BUG,*) '  Row',I
               WRITE(BUG,3006) (KE(I,J),J=1,ELDOF)
               WRITE(BUG,*)
            ENDDO
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! SEi, STEi stress recovery matrices

      IF (DUM_BUG(5) > 0) THEN

         WRITE(BUG,1000)

!zzzz    IF (PCOMP_PROPS == 'N') THEN

            WRITE(BUG,1001) ELDT_BUG_SE_BIT, TYPE, EID

!xx         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE(1:5) /= 'TRIA3') .AND. (TYPE(1:5) /= 'QUAD4') .AND                               &
!xx             (TYPE(1:4) /= 'HEXA') .AND. (TYPE(1:5) /= 'PENTA') .AND. (TYPE(1:5) /= 'TETRA')) THEN

            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
               IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                  WRITE(BUG,1003) (OFFSET(I),I=1,ELGP)
               ENDIF
            ELSE
               WRITE(BUG,*) '  SEi, STEi matrices do not include the effects of offsets (this is done later in the code)'
               WRITE(BUG,*)
            ENDIF

            DO I=1,METYPE
               IF (TYPE == ELMTYP(I)) THEN
                  NUM_STRESS_MATS = NUM_SEi(I)
               ENDIF
            ENDDO
            WRITE(BUG,*) '  This element has ',NUM_STRESS_MATS,' stress recovery points. The first is at the center of the element.'
            WRITE(BUG,*) '  Subsequent stress recovery points are at the element corners or at a mesh of Gauss points.'
            WRITE(BUG,*)
            WRITE(BUG,*) '  For the MIN4T QUAD4 element the 2nd-5th recovery points are values at the center of each of the'
            WRITE(BUG,*) '  4 non-overlapping TRIA elements making up the QUAD4 and the 1st recovery point matrices are'
            WRITE(BUG,*) '  the average of the ones for points 2-5'
            WRITE(BUG,*)

            DO K=1,NUM_STRESS_MATS
               WRITE(BUG,*) '  SE1 matrix for stress data recovery point ',k
               WRITE(BUG,*) '  ---------------------------------------------'
               DO I=1,3
                  WRITE(BUG,*) '  Row',I
                  WRITE(BUG,3006) (SE1(I,J,K),J=1,ELDOF)
                  WRITE(BUG,*)
               ENDDO
               WRITE(BUG,*)
            ENDDO

            DO K=1,NUM_STRESS_MATS
               WRITE(BUG,*) '  SE2 matrix for stress data recovery point ',k
               WRITE(BUG,*) '  ---------------------------------------------'
               DO I=1,3
                  WRITE(BUG,*) '  Row',I
                  WRITE(BUG,3006) (SE2(I,J,K),J=1,ELDOF)
                  WRITE(BUG,*)
               ENDDO
               WRITE(BUG,*)
            ENDDO

            DO K=1,NUM_STRESS_MATS
               WRITE(BUG,*) '  SE3 matrix for stress data recovery point ',k
               WRITE(BUG,*) '  ---------------------------------------------'
               DO I=1,3
                  WRITE(BUG,*) '  Row',I
                  WRITE(BUG,3006) (SE3(I,J,K),J=1,ELDOF)
                  WRITE(BUG,*)
               ENDDO
               WRITE(BUG,*)
            ENDDO

            IF (NTSUB > 0) THEN

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  STE1 matrix for thermal stress/strain effects for stress data recovery point ',k
                  WRITE(BUG,*) '  --------------------------------------------------------------------------------'
                  DO J=1,NTSUB
                     WRITE(BUG,*) '  Column',J
                     WRITE(BUG,3003) (STE1(I,J,K),I=1,3)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  STE2 matrix for thermal stress/strain effects for stress data recovery point ',k
                  WRITE(BUG,*) '  --------------------------------------------------------------------------------'
                  DO J=1,NTSUB
                     WRITE(BUG,*) '  Column',J
                     WRITE(BUG,3003) (STE2(I,J,K),I=1,3)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  STE3 matrix for thermal stress/strain effects for stress data recovery point ',k
                  WRITE(BUG,*) '  --------------------------------------------------------------------------------'
                  DO J=1,NTSUB
                     WRITE(BUG,*) '  Column',J
                     WRITE(BUG,3003) (STE3(I,J,K),I=1,3)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

            ENDIF

            WRITE(BUG,*) '  Factors that convert engineering stresses to enginering forces:'
            WRITE(BUG,*) '  --------------------------------------------------------------'
            WRITE(BUG,'(A,1ES14.6)') '  Factor for membrane     stresses, FCONV(1) = ', FCONV(1)
            WRITE(BUG,'(A,1ES14.6)') '  Factor for bending      stresses, FCONV(2) = ', FCONV(2)
            WRITE(BUG,'(A,1ES14.6)') '  Factor for transv shear stresses, FCONV(3) = ', FCONV(3)
            WRITE(BUG,*)
            WRITE(BUG,*)

            IF ((TYPE(1:4) == 'BUSH' ) .OR. (TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:6) == 'USERIN') .OR.   &
                (TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  BE1 matrix for strain data recovery point ',k
                  WRITE(BUG,*) '  ---------------------------------------------'
                  DO I=1,3
                     WRITE(BUG,*) '  Row',I
                     WRITE(BUG,3006) (BE1(I,J,K),J=1,ELDOF)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  BE2 matrix for strain data recovery point ',k
                  WRITE(BUG,*) '  ---------------------------------------------'
                  DO I=1,3
                     WRITE(BUG,*) '  Row',I
                     WRITE(BUG,3006) (BE2(I,J,K),J=1,ELDOF)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

               DO K=1,NUM_STRESS_MATS
                  WRITE(BUG,*) '  BE3 matrix for strain data recovery point ',k
                  WRITE(BUG,*) '  ---------------------------------------------'
                  DO I=1,3
                     WRITE(BUG,*) '  Row',I
                     WRITE(BUG,3006) (BE3(I,J,K),J=1,ELDOF)
                     WRITE(BUG,*)
                  ENDDO
                  WRITE(BUG,*)
               ENDDO

            ENDIF

!zzzz    ENDIF

      ENDIF

! **********************************************************************************************************************************
! Write output from LINK9 calculation of element displs, nodal forces

      IF (DUM_BUG(6) > 0) THEN

         WRITE(BUG,1000)

         WRITE(BUG,1001) ELDT_BUG_U_P_BIT, TYPE, EID

!        IF     (ELFORCEN == 'LOCAL') THEN

            NAME2 = 'AGRID, UEL displs of ' ; NAME3 = ' local elem ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = UEL(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

            NAME2 = 'AGRID, PEL loads at  ' ; NAME3 = ' local elem ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = PEL(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

!        ELSE IF (ELFORCEN == 'BASIC') THEN

            NAME2 = 'AGRID, UEB displs of ' ; NAME3 = ' basic      ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = UEB(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

            NAME2 = 'AGRID, PEB loads at  ' ; NAME3 = ' basic      ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = PEB(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

!        ELSE IF (ELFORCEN == 'GLOBAL') THEN

            NAME2 = 'AGRID, UEG displs of ' ; NAME3 = ' global     ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = UEG(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

            NAME2 = 'AGRID, PEG loads at  ' ; NAME3 = ' global     ' ; WRITE(BUG,2010) NAME2, NAME3, CASE_NUM
            I2 = 0
            DO I=1,ELGP
               OEL = ZERO
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  OEL(J) = PEG(I2)
               ENDDO
               WRITE(BUG,3007) AGRID(I), (OEL(J),J=1,6)
            ENDDO
            WRITE(BUG,*)
            WRITE(BUG,*)

!        ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1000 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------')

 1001 FORMAT(' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,/,                                                      &
             ' ==============================================================',/)

 1003 FORMAT('   Grid point offset indicators: ',4(1X,A1))

 1004 FORMAT(1X,I8,1X,I8,1X,I6,2X,3(1ES15.6))

 1113 FORMAT(3X,'End of v vector:',3(1ES15.6))

 1005 FORMAT(8(1ES14.6))

 1006 FORMAT('   Row ',I8,' of ',I8,/,'   ------------------------')

 1007 FORMAT(' (does not include effects of offsets for BAR, BEAM or ROD)')

 2001 FORMAT('   Thermal  load matrix for subcase number ',  I8  ,' (in local element coordinates)',/,    &
             '   -------------------------------------------------------------------------------------------------------')

 2003 FORMAT('   Pressure load matrix for subcase number ',  I8  ,' (in local element coordinates)',/,    &
             '   ---------------------------------------------------------------------------------------------------------')

 2010 FORMAT(3X,A,' elem nodes in ',A,' coordinate system for subcase ',I8                                                     ,//,&
             '    Grid           T1            T2            T3            R1            R2            R3'                      ,/,&
             '    ----           --            --            --            --            --            --')
 2024 FORMAT(24I3)

 3002 FORMAT(2(1ES14.6))

 3003 FORMAT(3(1ES14.6))

 3004 FORMAT(4(1ES14.6))

 3005 FORMAT(5(1ES14.6))

 3006 FORMAT(6(1ES14.6))

 3007 FORMAT(I8,4X,6(1ES14.6))

 3024 FORMAT(4(1ES14.6),A)

 3026 FORMAT(10(1ES14.6))

 4001 FORMAT('   ELAS_COMP - Displacement components that ELAS elem connects to: comp',I2,' at ',A,I8,' and comp',I2,' at ',A,I8,/,&
             '   ----------------------------------------------------------------------------------------------------------------',&
             '-----------------',/)

 5001 FORMAT('   Length of element between grids 1 and 2 (including effects of offsets) = ',1ES13.6,/,                             &
             '   ----------------------------------------------------------------------',/)

 5002 FORMAT('   HBAR warp of quadrilateral element = ',1ES13.6,'. BMEAN will be used if HBAR > MXWARP = ',                        &
                 1ES13.6/,                                                                                                         &
             '   ---------------------------------------------------------------------------------------------------------------',/)

 5003 FORMAT('   QUAD_THETA angle between side 1-2 and diagonal 1-3                        = ',1ES13.6,' (degrees)')

 5004 FORMAT('   QUAD_GAMMA angle between side 1-2 and diagonal 2-4                        = ',1ES13.6,' (degrees)')

 5005 FORMAT('   QUAD_DELTA angle to rotate from side 1-2 to split angle between diagonals = ',1ES13.6,' (degrees) = ',            &
                 '(QUAD_THETA - QUAD_GAMMA)/2'                                                                                  ,/,&
             '   -------------------------------------------------------------------------')

 5006 FORMAT('   THETAM angle to rotate from local x axis to material orientation axis = ',1ES13.6,' (degrees)',/,                 &
             '   ---------------------------------------------------------------------',/)

 5007 FORMAT(A,1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE ELMOUT


      SUBROUTINE ELMTLB ( OPT )

! Transforms element matrices from local to basic coordinates. Matrices transformed are: ME, KE, KED, PTE, PPE, using elem coord
! transformation matrix TE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MELDOF, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ELDOF, ELGP, KE, KED, ME, PTE, PPE, TE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATERIAL_MATRIX_TRANSFER, ONLY:  MATGET, MATPUT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMTLB'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)

      INTEGER(LONG)                   :: BEG_COL           ! Beginning col of matrix to get partition from
      INTEGER(LONG)                   :: BEG_ROW           ! Beginning row of matrix to get partition from
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: NCOL              ! No. cols to get/put for subrs MATGET/MATPUT, called herein
      INTEGER(LONG), PARAMETER        :: NCOLA     = 3     ! No. cols in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOLB             ! No. cols in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG), PARAMETER        :: NROW      = 3     ! No. rows to get/put for subrs MATGET/MATPUT, called herein
      INTEGER(LONG), PARAMETER        :: NROWA     = 3     ! No. rows in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein


      REAL(DOUBLE)                    :: DUM11(3,3)        ! An intermediate result when calculating transformed KE
      REAL(DOUBLE)                    :: DUM12(3,3)        ! An intermediate result when calculating transformed KE
      REAL(DOUBLE)                    :: PDUM1(3,NSUB)     ! An intermediate result when calculating transformed PTE, PPE
      REAL(DOUBLE)                    :: PDUM2(3,NSUB)     ! An intermediate result when calculating transformed PTE, PPE



! **********************************************************************************************************************************
      IF (OPT(1) == 'Y') THEN                              ! Transform ME to TE' x ME x TE
         NCOL  = 3
         NCOLB = 3
         DO I=1,2*ELGP
            BEG_ROW = 3*I - 2
            DO J=I,2*ELGP
               BEG_COL = 3*J - 2
               CALL MATGET ( ME, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, DUM11 )
               CALL MATMULT_FFF   ( DUM11, TE, NROWA, NCOLA, NCOLB, DUM12 )
               CALL MATMULT_FFF_T ( TE, DUM12, NROWA, NCOLA, NCOLB, DUM11 )
               CALL MATPUT ( DUM11, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, ME )
            ENDDO
         ENDDO

         DO I=1,ELDOF                                      ! Set lower portion of ME using symmetry.
            DO J=1,I-1
               ME(I,J) = ME(J,I)
            ENDDO
         ENDDO

      ENDIF

      IF ((OPT(2) == 'Y') .AND. (NTSUB > 0)) THEN          ! Transform PTE to TE' x PTE
         NCOL  = NTSUB
         NCOLB = NTSUB
         DO I=1,2*ELGP
            BEG_ROW = 3*I - 2
            BEG_COL = 1
            CALL MATGET ( PTE, MELDOF, NTSUB, BEG_ROW, BEG_COL, NROW, NCOL, PDUM1 )
            CALL MATMULT_FFF_T ( TE, PDUM1, NROWA, NCOLA, NCOLB, PDUM2 )
            CALL MATPUT ( PDUM2, MELDOF, NTSUB, BEG_ROW, BEG_COL, NROW, NTSUB, PTE )
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN                              ! Transform KE to TE' x KE x TE
         NCOL  = 3
         NCOLB = 3
         DO I=1,2*ELGP
            BEG_ROW = 3*I - 2
            DO J=I,2*ELGP
               BEG_COL = 3*J - 2
               CALL MATGET ( KE, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, DUM11 )
               CALL MATMULT_FFF   ( DUM11, TE, NROWA, NCOLA, NCOLB, DUM12 )
               CALL MATMULT_FFF_T ( TE, DUM12, NROWA, NCOLA, NCOLB, DUM11 )
               CALL MATPUT ( DUM11, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, KE )
            ENDDO
         ENDDO



         DO I=1,ELDOF                                      ! Set lower portion of KE using symmetry.
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF

      IF (OPT(5) == 'Y') THEN                              ! Transform PPE to TE' x PPE
         NCOL  = NSUB
         NCOLB = NSUB
         DO I=1,2*ELGP
            BEG_ROW = 3*I - 2
            BEG_COL = 1
            CALL MATGET ( PPE, MELDOF, NSUB, BEG_ROW, BEG_COL, NROW, NCOL, PDUM1 )
            CALL MATMULT_FFF_T ( TE, PDUM1, NROWA, NCOLA, NCOLB, PDUM2 )
            CALL MATPUT ( PDUM2, MELDOF, NSUB, BEG_ROW, BEG_COL, NROW, NSUB, PPE )
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN                              ! Transform KED to TE' x KED x TE
         NCOL  = 3
         NCOLB = 3
         DO I=1,2*ELGP
            BEG_ROW = 3*I - 2
            DO J=I,2*ELGP
               BEG_COL = 3*J - 2
               CALL MATGET ( KED, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, DUM11 )
               CALL MATMULT_FFF   ( DUM11, TE, NROWA, NCOLA, NCOLB, DUM12 )
               CALL MATMULT_FFF_T ( TE, DUM12, NROWA, NCOLA, NCOLB, DUM11 )
               CALL MATPUT ( DUM11, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, KED )
            ENDDO
         ENDDO

         DO I=1,ELDOF                                      ! Set lower portion of KED using symmetry.
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ELMTLB

   END MODULE ELEMENT_DIAGNOSTICS
