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

   MODULE MATERIAL_TRANSFORMATIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MATL_TRANSFORM_MATRIX, ROT_AXES_MATL_TO_LOC

   CONTAINS

      SUBROUTINE MATL_TRANSFORM_MATRIX ( T21, TS )

! Formulates 6x6 transformation matrix TS which can do a 2D plane coord transformation for stresses, strains and material matrices

! TS is formed from the coefficients in 3x3 coord transformation matrix T21 where T21 is the coord transformation matrix that
! transforms a unit vector in coord sys 1 (U1) to a unit vector in coord sys 2 (U2):

!                                                  U2 = T21*U1                                                         (1)

! For a plane coord transformation (2D, not 3D)

!                                                       |  cos(THETA)  sin(THETA)  0 |
!                                                 T21 = | -sin(THETA)  cos(THETA)  0 |                                 (2)
!                                                       |      0           0       1 |

! The stresses in any coord sys is a tensor (call it S) that is written as:
!                                                    ~
!                                                       | Sxx  Sxy  Sxz |
!                                                   S = | Syx  Syy  Syz |                                              (3)
!                                                   ~   | Szx  Szy  Szz |

! Matrix T21 is used to transform stress tensors from coord sys 1 to coord sys 2 via the equation (NOTE: T21' indicates
! the transpose of matrix T21 which is also T12):


!                                                 S22 = T21'*S11*T21 = T12*S11*T21                                     (4)
!                                                 ~          ~             ~

! The terms in this triple matrix product can be used to rewrite the stress transformation with the 6 independent stresses
! in a vector (instead of as the 3x3 stress tensor). To this end, define the vector:

!                                                       | Sxx |
!                                                       | Syy |
!                                                   S = | Szz |                                                        (5)
!                                                       | Sxy |
!                                                       | Syz |
!                                                       | Szx |

! Then the relationship between S in coord sys 2 and that in coord sys 1 can be written as:

!                                                --------------
!                                               |              |
!                                               |  S2 = TS*S1  |                                                       (6)
!                                               |              |
!                                                --------------

! where TS is a 6x6 matrix of coefficients from T12 (eqn (2)). The 9 terms in T21 are denoted as Aij and the 36 terms in TS,
! developed from the T21 using eqn (4), are developed in the code below. Again, they are developed by expanding the triple
! matrix product (symbolically) in eqn (4) and then rewriting the 6 independent of the resulting 9 equations in the 6x6
! matrix form of (6).

! Matrix TS can also be used to rotate the material matrices from coord sys 1 to coord sys 2 as follows.
! The stress-strain relations in coord sys 1 and 2 are:

!                                                  S1 = E1*e1
!                                                                                                                      (7)
!                                                  S2 = E2*e2

! where ei are strains and Ei is the material matrix in coord sys 1 and 2. Eqn 6 shows how to transform stresses. It is
! the same for strains so that eqns (7) can be written:

!                                                  S2 = TS*S1 = TS*E1*e1                                               (8)

! but e1 in terms of e2, using the same logic in eqn(6) is

!                                                  e1 = TS'*e2                                                         (9)

! Substituting 9 into 8 yields:

!                                                  S2 = (Ts*E1*TS')*e2 = E2*e2                                         (10)

! so that we conclude:

!                                                ------------------
!                                               |                  |
!                                               |  E2 = TS*E1*TS'  |                                                   (11)
!                                               |                  |
!                                                ------------------

! ----------------------------------------------------------------------------------------------------------------------------------

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATL_TRANSFORM_MATRIX'



      REAL(DOUBLE), INTENT(IN)        :: T21(3,3)          ! 3x3 matrix that transforms a vector in coord sys 1 to coord sys 2
      REAL(DOUBLE), INTENT(OUT)       :: TS(6,6)           ! 6x6 stress transformation matrix
      REAL(DOUBLE)                    :: A11,A12,A13       ! Coefficients from matrix TME
      REAL(DOUBLE)                    :: A21,A22,A23       ! Coefficients from matrix TME
      REAL(DOUBLE)                    :: A31,A32,A33       ! Coefficients from matrix TME



! **********************************************************************************************************************************
! Formulate 6x6 matrix TS from terms in 3x3 matrix T12

      A11 = T21(1,1)     ;     A12 = T21(1,2)     ;     A13 = T21(1,3)
      A21 = T21(2,1)     ;     A22 = T21(2,2)     ;     A23 = T21(2,3)
      A31 = T21(3,1)     ;     A32 = T21(3,2)     ;     A33 = T21(3,3)

! Matrix TS is formed from the coefficients in matrix T12 (denoted as Aij coefficients)

! Row 1 of TS:

      TS(1,1)= A11*A11           ;   TS(1,2)= A21*A21           ;   TS(1,3)= A31*A31;
      TS(1,4)= 2*A11*A21         ;   TS(1,5)= 2*A21*A31         ;   TS(1,6)= 2*A11*A31

! Row 2 of TS:

      TS(2,1)= A12*A12           ;   TS(2,2)= A22*A22           ;   TS(2,3)= A32*A32;
      TS(2,4)= 2*A12*A22         ;   TS(2,5)= 2*A22*A32         ;   TS(2,6)= 2*A12*A32

! Row 3 of TS:

      TS(3,1)= A13*A13           ;   TS(3,2)= A23*A23           ;   TS(3,3)= A33*A33;
      TS(3,4)= 2*A13*A23         ;   TS(3,5)= 2*A23*A33         ;   TS(3,6)= 2*A13*A33

! Row 4 of TS:

      TS(4,1)= A11*A12           ;   TS(4,2)= A21*A22           ;   TS(4,3)= A31*A32;
      TS(4,4)= A11*A22+A21*A12   ;   TS(4,5)= A21*A32+A31*A22   ;   TS(4,6)= A11*A32+A31*A12

! Row 5 of TS:

      TS(5,1)= A12*A13           ;   TS(5,2)= A22*A23           ;   TS(5,3)= A32*A33;
      TS(5,4)= A12*A23+A22*A13   ;   TS(5,5)= A22*A33+A32*A23   ;   TS(5,6)= A12*A33+A32*A13

! Row 6 of TS:

      TS(6,1)= A13*A11           ;   TS(6,2)= A23*A21           ;   TS(6,3)= A33*A31;
      TS(6,4)= A13*A21+A23*A11   ;   TS(6,5)= A23*A31+A33*A21   ;   TS(6,6)= A13*A31+A33*A11



      RETURN

! **********************************************************************************************************************************
99664 format('  Transformation matrix TS: transforms 6 stress and 6x6 matl matrices from matl to elem axes')

99667 format(3(1es14.6),4x,3(1es14.6))

! **********************************************************************************************************************************

      END SUBROUTINE MATL_TRANSFORM_MATRIX



      SUBROUTINE ROT_AXES_MATL_TO_LOC ( WRITE_WARN )

! Rotates material and CTE matrices from the material axes (specified on connection entry) to element local axes

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEMATC, NCORD
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE MODEL_STUF, ONLY            :  ALPVEC, CORD, ISOLID, EB, EBM, EM, ES, ET, MTRL_TYPE, NUM_EMG_FATAL_ERRS, QUAD_DELTA,     &
                                         RCORD, TE, THETAM, TYPE
      USE PARAMS, ONLY                :  EPSIL
      USE DEBUG_PARAMETERS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  PLANE_COORD_TRANS_21
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ROT_AXES_MATL_TO_LOC'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y' write warning messages, otherwise do not
      CHARACTER( 1*BYTE)              :: FOUND             ! If 'Y' we found something we were looking for


      INTEGER(LONG)                   :: CORDM             ! Actual coord system ID (CORDM on PSOLID Bulk Data entry)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ICORD             ! Internal coord system ID for CORDM
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE)                    :: DUM22(2,2)        ! Intermediate matrix in calculating outputs
      REAL(DOUBLE)                    :: DUM33(3,3)        ! Intermediate matrix in calculating outputs
      REAL(DOUBLE)                    :: DUM66(6,6)        ! Intermediate matrix in calculating outputs
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: MATL_AXES_ROTATE  ! Angle to rotate material axes to coincide with local elem x axis
      REAL(DOUBLE)                    :: T0M(3,3)          ! Transform matrix from material coords to basic coords
      REAL(DOUBLE)                    :: T1(6,6)           ! Transforms 6 stress and 6x6 material mats from material to element axes
!                                                            (T1' transforms strains)
      REAL(DOUBLE)                    :: T1_MB(3,3)        ! Portion of T1: transforms 3x3 EM, EB, EBM from material to elem axes
      REAL(DOUBLE)                    :: T1_TS(2,2)        ! Portion of T1: transforms 3x3 ET from material to elem axes
      REAL(DOUBLE)                    :: TEM(3,3)          ! TE*T0M
      REAL(DOUBLE)                    :: TME(3,3)          ! Coord transf matrix which will rotate a vector in local element coord
!                                                            system to a vector in the material coord system (Um = TME*Ue)

      REAL(DOUBLE)                    :: EB0(3,3)          ! Plane stress matl matrix for bending before coord transformation
      REAL(DOUBLE)                    :: EM0(3,3)          ! Plane stress matl matrix for in-plane stress bef coord transformation
      REAL(DOUBLE)                    :: EBM0(3,3)         ! Bend/membr coupling matl matrix before coord transformation
      REAL(DOUBLE)                    :: ES0(6,6)          ! 3D stress matl matrix before coord transformation
      REAL(DOUBLE)                    :: ET0(2,2)          ! 2D transverse shear matl matrix before coord transformation



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Set initial values before coord transformation (NOTE: this is only needed for use in DEBUG(51))

      DO I=1,3
         DO J=1,3
            EB0(I,J)  = EB(I,J)
            EM0(I,J)  = EM(I,J)
            EBM0(I,J) = EBM(I,J)
         ENDDO
      ENDDO

      DO I=1,6
         DO J=1,6
            ES0(I,J)  = ES(I,J)
         ENDDO
      ENDDO

      DO I=1,2
         DO J=1,2
            ET0(I,J)  = ET(I,J)
         ENDDO
      ENDDO

! Now do coord transformation

      IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'TRIA3')) THEN

         IF      (TYPE(1:5) == 'QUAD4') THEN

            MATL_AXES_ROTATE = THETAM - QUAD_DELTA         ! We need angle which would rotate local elem x axis to mat'l axis,
!                                                            NOT the opposite (since we want a transform matrix that will take
!                                                            a vector in elem coords and convert it into a vector in mat'l coords
         ELSE IF (TYPE(1:5) == 'TRIA3') THEN

            MATL_AXES_ROTATE = THETAM

         ENDIF

         IF (DABS(MATL_AXES_ROTATE) > EPS1) THEN

            CALL PLANE_COORD_TRANS_21 ( MATL_AXES_ROTATE, TME, SUBR_NAME )
            CALL MATL_TRANSFORM_MATRIX ( TME, T1 )
            T1 = TRANSPOSE(T1)

                                                           ! T1_MB is for Sxx, Syy, Sxy which are rows and cols 1,2,4 from T1
            T1_MB(1,1) = T1(1,1)     ;     T1_MB(1,2) = T1(1,2)     ;     T1_MB(1,3) = T1(1,4)
            T1_MB(2,1) = T1(2,1)     ;     T1_MB(2,2) = T1(2,2)     ;     T1_MB(2,3) = T1(2,4)
            T1_MB(3,1) = T1(4,1)     ;     T1_MB(3,2) = T1(4,2)     ;     T1_MB(3,3) = T1(4,4)
                                                          ! T1_ET is for Syz, Szx which are rows and cols 5,6 from T1
            T1_TS(1,1) = T1(5,5)     ;     T1_TS(1,2) = T1(5,6)
            T1_TS(2,1) = T1(6,5)     ;     T1_TS(2,2) = T1(6,6)

            IF ((MTRL_TYPE(1) == 2) .OR. (MTRL_TYPE(1) == 8)) THEN   ! Transform membrane material matrix
               CALL MATMULT_FFF   ( EM , T1_MB  , 3, 3, 3, DUM33 )
               CALL MATMULT_FFF_T ( T1_MB, DUM33, 3, 3, 3, EM    )
            ENDIF

            IF ((MTRL_TYPE(2) == 2) .OR. (MTRL_TYPE(2) == 8)) THEN   ! Transform bending material matrix
               CALL MATMULT_FFF   ( EB , T1_MB  , 3, 3, 3, DUM33 )
               CALL MATMULT_FFF_T ( T1_MB, DUM33, 3, 3, 3, EB    )
            ENDIF

            IF ((MTRL_TYPE(3) == 2) .OR. (MTRL_TYPE(3) == 8)) THEN   ! Transform transverse shear material matrix
               CALL MATMULT_FFF   ( ET , T1_TS  , 2, 2, 2, DUM22 )
               CALL MATMULT_FFF_T ( T1_TS, DUM22, 2, 2, 2, ET    )
            ENDIF

            IF ((MTRL_TYPE(4) == 2) .OR. (MTRL_TYPE(4) == 8)) THEN   ! Transform coupled bending/membrane material matrix
               CALL MATMULT_FFF   ( EBM, T1_MB  , 3, 3, 3, DUM33 )
               CALL MATMULT_FFF_T ( T1_MB, DUM33, 3, 3, 3, EBM   )
            ENDIF

            ALPVEC(4:6,:) = ALPVEC(4:6,:) / TWO                      ! Remove shear factor of 2 to transform.
            ALPVEC = MATMUL(TRANSPOSE(T1), ALPVEC)
            ALPVEC(4:6,:) = ALPVEC(4:6,:) * TWO                      ! Reinstate shear factor of 2 after transform.

         ENDIF

         IF (DEBUG(51) > 0) THEN
            CALL DEBUG_ROT_AXES_1
         ENDIF

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

         IF (ISOLID(3) /= -1) THEN                         ! If -1 ES already in elem coords, else transf ES from basic to elem axes

            CORDM = ISOLID(3)                              ! Get transf matrix T0M. CORDM mtrl coord sys from PSOLID Bulk data entry

            T0M(:,:) = ZERO

            IF (CORDM == 0) THEN                           ! CORDM was basic so set T0M to the identity matrix

               DO I=1,3
                  T0M(I,I) = ONE
               ENDDO

            ELSE                                           ! Need to transform ES mat'l matrix to basic coords

               FOUND = 'N'
               DO I=1,NCORD                                ! Get the internal coord system ID for CORDM
                  IF (CORDM == CORD(I,2)) THEN
                     ICORD = I
                     FOUND = 'Y'
                     EXIT
                  ENDIF
               ENDDO

               IF (FOUND == 'Y') THEN
                  DO I=1,3
                     DO J=1,3
                        K = 3 + 3*(I-1) + J
                        T0M(I,J) = RCORD(ICORD,K)
                     ENDDO
                  ENDDO
               ELSE
                  WRITE(ERR,1822) 'COORD SYSTEM ', CORDM, 'PSOLID', ISOLID(1)
                  WRITE(F06,1822) 'COORD SYSTEM ', CORDM, 'PSOLID', ISOLID(1)
                  FATAL_ERR = FATAL_ERR + 1
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               ENDIF

            ENDIF

            CALL MATMULT_FFF ( TE , T0M  , 3, 3, 3, TEM )
            TME = TRANSPOSE(TEM)

            CALL MATL_TRANSFORM_MATRIX ( TME, T1 )
            T1 = TRANSPOSE(T1)

            CALL MATMULT_FFF   ( ES , T1   , 6, 6, 6, DUM66 )
            CALL MATMULT_FFF_T ( T1 , DUM66, 6, 6, 6, ES    )

            ALPVEC(4:6,:) = ALPVEC(4:6,:) / TWO                      ! Remove shear factor of 2 to transform.
            ALPVEC = MATMUL(TRANSPOSE(T1), ALPVEC)
            ALPVEC(4:6,:) = ALPVEC(4:6,:) * TWO                      ! Reinstate shear factor of 2 after transform.

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_ROT_AXES_1

      IMPLICIT NONE

      CHARACTER(LEN=100)              :: MI(9)             ! Messages to be written out

! **********************************************************************************************************************************
      DO I=1,9
         MI(I)(1:) = ' '
      ENDDO

      WRITE(F06,98720)
      IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'TRIA3')) THEN

         IF (DABS(MATL_AXES_ROTATE) > EPS1) THEN

            MI(1) = ' Transforms 6 stress and 6x6 matl matrices from matl to elem axes (T1 transpose transforms strains)'
            WRITE(F06,99664) 'T1', MI(1)
            DO I=1,3
               WRITE(F06,99667) (T1(I,J),J=1,6)
            ENDDO
            WRITE(F06,*)
            DO I=4,6
               WRITE(F06,99667) (T1(I,J),J=1,6)
            ENDDO
            WRITE(F06,*)

            MI(2) = ' Portion of T1: transforms 3x3 EM, EB, EBM from matl to elem axes'
            WRITE(F06,99664) 'T1_MB', MI(2)
            DO I=1,3
               WRITE(F06,99668) (T1_MB(I,J),J=1,3)
            ENDDO
            WRITE(F06,*)

            MI(3) = '  Portion of T1: transforms 2x2 ET from material to elem axes'
            WRITE(F06,99664) 'T1_TS', MI(3)
            DO I=1,2
               WRITE(F06,99669) (T1_TS(I,J),J=1,2)
            ENDDO
            WRITE(F06,*)

            MI(4) = ' Coord transf matrix which will rotate a vector in local elem'
            MI(5) = ' coord system to a vector in the material coord system (Um = TME*Ue)'
            WRITE(F06,99664) 'TME', MI(4)   ;   WRITE(F06,99663) MI(5)
            DO I=1,3
               WRITE(F06,99668) (TME(I,J),J=1,3)
            ENDDO
            WRITE(F06,*)

mem:        IF      ((MTRL_TYPE(1) == 2) .OR. (MTRL_TYPE(1) == 8)) THEN

               WRITE(F06,99665) 'EM '
               DO I=1,3
                  WRITE(F06,99668) (EM0(I,J),J=1,3), (EM(I,J),J=1,3)
               ENDDO
               WRITE(F06,*)   ;   WRITE(F06,*)

               WRITE(F06,99665) 'EB '

            ENDIF mem

bend:       IF ((MTRL_TYPE(2) == 2) .OR. (MTRL_TYPE(2) == 8)) THEN

               DO I=1,3
                  WRITE(F06,99668) (EB0(I,J),J=1,3), (EB(I,J),J=1,3)
               ENDDO
               WRITE(F06,*)   ;   WRITE(F06,*)

            ENDIF bend

tr_shr:     IF ((MTRL_TYPE(3) == 2) .OR. (MTRL_TYPE(3) == 8)) THEN

               WRITE(F06,99665) 'ET '
               DO I=1,2
                  WRITE(F06,99669) (ET0(I,J),J=1,2), (ET(I,J),J=1,2)
               ENDDO
               WRITE(F06,*)   ;   WRITE(F06,*)

            ENDIF tr_shr

mem_bend:   IF ((MTRL_TYPE(4) == 2) .OR. (MTRL_TYPE(4) == 8)) THEN

               WRITE(F06,99665) 'EBM'
               DO I=1,3
                  WRITE(F06,99668) (EBM0(I,J),J=1,3), (EBM(I,J),J=1,3)
               ENDDO
               WRITE(F06,*)   ;   WRITE(F06,*)

            ENDIF mem_bend

         ENDIF

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

         WRITE(F06,'()') 'ISOLID(3) = ',ISOLID(3)
         WRITE(F06,'()') 'MATL MATRIX ES BEFORE COORD TRANSFORMATION = '
         DO I=1,6
            WRITE(F06,67549) (ES0(I,J),J=1,6)
         ENDDO
         WRITE(F06,*)

         IF (ISOLID(3) /= -1) THEN                         ! If -1 ES already in elem coords, else transf ES from basic to elem axes
            WRITE(F06,'()') 'MATL MATRIX ES BEFORE COORD TRANSFORMATION = '
            DO I=1,6
               WRITE(F06,67549) (ES(I,J),J=1,6)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

      WRITE(F06,*)

      WRITE(F06,98799)

      WRITE(F06,*)

! **********************************************************************************************************************************
67549 format(6(1es14.6))

99663 format(1x,a)

99664 format('  Transformation matrix ',2a)

99665 format(20x,'Material matrix ',a3,' before/after T1 transformation',/,                                                        &
             19x,'before',40x,'after',/,'  ----------------------------------------      ----------------------------------------')

99667 format(3(1es14.6), 4x,3(1es14.6))

99668 format(3(1es14.6), 4x,3(1es14.6))

99669 format(7x,2(1es14.6),18x,2(1es14.6))

98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(51) OUTPUT FROM SUBROUTINE ROT_AXES_MATL_TO_LOC:::::::::::::::::::',&
             ':::::::::::::::::',/)

98799 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(51) OUTPUT FROM SUBROUTINE ROT_AXES_MATL_TO_LOC::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_ROT_AXES_1

      END SUBROUTINE ROT_AXES_MATL_TO_LOC




   END MODULE MATERIAL_TRANSFORMATIONS
