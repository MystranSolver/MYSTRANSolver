! #################################################################################################################################
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

   MODULE MITC_KERNELS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MITC_ADD_TO_B, MITC_CONTRAVARIANT_BASIS, MITC_COVARIANT_BASIS, MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION, MITC_DETJ, MITC_ELASTICITY, MITC_INITIALIZE, MITC_SHAPE_FUNCTIONS, MITC_TRANSFORM_B, MITC_TRANSFORM_CONTRAVARIANT_TO_LOCAL, MITC4_CARTESIAN_LOCAL_BASIS

   CONTAINS
      SUBROUTINE MITC_ADD_TO_B ( B, POINT, COL, SCALAR, TENSOR )

! Add the UT of the 3x3 tensor times the scalar to a column of the B for sampling point POINT.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(INOUT)    :: B(6,6*8,4)
      REAL(DOUBLE) , INTENT(IN)       :: SCALAR
      REAL(DOUBLE) , INTENT(IN)       :: TENSOR(3,3)

      INTEGER(LONG), INTENT(IN)       :: POINT
      INTEGER(LONG), INTENT(IN)       :: COL

! **********************************************************************************************************************************

      B(1, COL, POINT) = B(1, COL, POINT) + SCALAR * TENSOR(1,1)
      B(2, COL, POINT) = B(2, COL, POINT) + SCALAR * TENSOR(2,2)
      B(3, COL, POINT) = B(3, COL, POINT) + SCALAR * TENSOR(3,3)
      B(4, COL, POINT) = B(4, COL, POINT) + SCALAR * TENSOR(1,2)
      B(5, COL, POINT) = B(5, COL, POINT) + SCALAR * TENSOR(2,3)
      B(6, COL, POINT) = B(6, COL, POINT) + SCALAR * TENSOR(1,3)

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_ADD_TO_B

      SUBROUTINE MITC_CONTRAVARIANT_BASIS ( G, G_CONTRA )

! Calculates the contravariant basis vectors g^1, g^2, g^3 to the specified covariant basis vectors
! g_r, g_s, g_t. This is the inverse Jacobian matrix where g_r/s/t are the columns of the Jacobian matrix.
! G(:,1) is g_r, etc.
! G_CONTRA(:,1) is g^1, etc.

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: G(3,3)            ! Covariant basis vectors
      REAL(DOUBLE) , INTENT(OUT)      :: G_CONTRA(3,3)     ! Contravariant basis vectors
      REAL(DOUBLE)                    :: DUM1(3)
      REAL(DOUBLE)                    :: DETJ


! **********************************************************************************************************************************

                                                           ! DET(J) = g_r . (g_s x g_t)
      CALL CROSS(G(:,2), G(:,3), DUM1)
      DETJ = G(1,1)*DUM1(1) + G(2,1)*DUM1(2) + G(3,1)*DUM1(3)

      CALL CROSS(G(:,2), G(:,3), G_CONTRA(:,1))
      CALL CROSS(G(:,3), G(:,1), G_CONTRA(:,2))
      CALL CROSS(G(:,1), G(:,2), G_CONTRA(:,3))

      G_CONTRA = G_CONTRA / DETJ

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_CONTRAVARIANT_BASIS

      SUBROUTINE MITC_COVARIANT_BASIS ( R, S, T, G )

! Calculates g_r, g_s, g_t in XEL element coordinates.
! These are also the columns of the Jacobian matrix.
! G(:,1) is g_r, etc.

! Ref [3] SesamX blog https://www.sesamx.io/blog/shell_finite_element/

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, XEL
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS


      IMPLICIT NONE

      INTEGER(LONG)                   :: GP                ! Element grid point number

      REAL(DOUBLE) , INTENT(IN)       :: R, S, T           ! Isoparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: G(3,3)            ! basis vector in basic coordinates
      REAL(DOUBLE)                    :: PSH(ELGP)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to xi and eta.

! **********************************************************************************************************************************

      CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

      G(:,:) = ZERO

      ! Interpolate from the values at nodes
      DO GP=1,ELGP
         ! g_r(r, s, t) = dX/dr = d/dr X + t/2 * d/dr (hv)
         !     = sum over nodes[ dN/dr X + t/2 * dN/dr (hv) ]
         G(:,1) = G(:,1) + XEL(GP,:) * DPSHG(1,GP) + DIRECTOR(:,GP) * T/TWO * DPSHG(1,GP) * DIR_THICKNESS(GP)
         G(:,2) = G(:,2) + XEL(GP,:) * DPSHG(2,GP) + DIRECTOR(:,GP) * T/TWO * DPSHG(2,GP) * DIR_THICKNESS(GP)
         ! Interpolate director vector * thickness.
         G(:,3) = G(:,3) + DIRECTOR(:,GP) * DIR_THICKNESS(GP) / TWO * PSH(GP)
      ENDDO

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_COVARIANT_BASIS

      SUBROUTINE MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION ( R, S, T, ROW_FROM, ROW_TO, B )

! Covariant strain-displacement components at point (R,S,T) directly evaluated from the displacement and rotation interpolations
!
!        Grid point 1        Grid point 2      ...
!      ux uy uz rx ry rz   ux uy uz rx ry rz
! 11 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]
! 22 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]
! 33 [ 0  0  0  0  0  0  | 0  0  0  0  0  0  |     ]
! 12 [ #  #  #  #  #  #  | #  #  #  #  #  #  | ... ]
! 23 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]
! 13 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]


      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, TYPE
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, FOUR
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS


      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: R,S,T             ! Isparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: B(6, 6*ELGP)      ! Strain-displacement matrix.
      REAL(DOUBLE)                    :: PSH(ELGP)         ! Shape functions
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: DPSHG3(3,ELGP)    ! Derivatives of shape functions with respect to R, S, T.
      REAL(DOUBLE)                    :: G(3,3)            ! Covariant basis vectors (Jacobian matrix) in basic coordinates

      INTEGER(LONG), INTENT(IN)       :: ROW_FROM          ! First row of B to generate. Strain component index 1-6.
      INTEGER(LONG), INTENT(IN)       :: ROW_TO            ! Last row of B to generate. Strain component index 1-6.
      INTEGER(LONG)                   :: GP                ! Element grid point number
      INTEGER(LONG)                   :: I,J,K,L           ! Loop and tensor indices
      INTEGER(LONG)                   :: ROW

! **********************************************************************************************************************************

! Reference [2]:
!  MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!     by Dvorkin and Bathe

                                                           ! Shape function derivatives at R,S
      CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

                                                           ! Extend shape function derivates to include
                                                           ! thickness direction (0) for convenience.
      DO GP=1,ELGP
        DO J=1,2
          DPSHG3(J,GP) = DPSHG(J,GP)
        ENDDO
        DPSHG3(3,GP) = ZERO
      ENDDO


! Jacobian matrix
      CALL MITC_COVARIANT_BASIS( R, S, T, G )

! Build strain-displacement matrix
! From "A continuum mechanics based four-node shell element for general nonlinear analysis" by Dvorkin and Bathe
! equations (21a), (22a), (23a), (24a)

      DO ROW=1,6
        DO GP=1,ELGP

          K = (GP-1) * 6

          IF ((ROW >= ROW_FROM) .AND. (ROW <= ROW_TO) .AND. ROW /= 3) THEN

                                                           ! Tensor indices for the row
            SELECT CASE (ROW)
              CASE (1); I=1; J=1                           ! In-layer normal strain
              CASE (2); I=2; J=2                           ! In-layer normal strain
              CASE (4); I=1; J=2                           ! In-layer shear strain
              CASE (5); I=2; J=3                           ! Transverse shear strain
              CASE (6); I=1; J=3                           ! Transverse shear strain
            END SELECT

                                                ! e_ij = sum over nodes [ 1/2 d/di u_0 dot g_j  +  1/2 d/dj u_0 dot g_i  + ...
            B(ROW, K+1) = (DPSHG3(I,GP) * G(1,J) + DPSHG3(J,GP) * G(1,I)) / TWO
            B(ROW, K+2) = (DPSHG3(I,GP) * G(2,J) + DPSHG3(J,GP) * G(2,I)) / TWO
            B(ROW, K+3) = (DPSHG3(I,GP) * G(3,J) + DPSHG3(J,GP) * G(3,I)) / TWO
                                                           !... 1/4 d/di (t h phi) dot g_j  + ...
            B(ROW, K+4) = DIR_THICKNESS(GP) * T * DPSHG3(I,GP) * (G(3,J) * DIRECTOR(2,GP) - G(2,J) * DIRECTOR(3,GP)) / FOUR
            B(ROW, K+5) = DIR_THICKNESS(GP) * T * DPSHG3(I,GP) * (G(1,J) * DIRECTOR(3,GP) - G(3,J) * DIRECTOR(1,GP)) / FOUR
            B(ROW, K+6) = DIR_THICKNESS(GP) * T * DPSHG3(I,GP) * (G(2,J) * DIRECTOR(1,GP) - G(1,J) * DIRECTOR(2,GP)) / FOUR
                                                           !... 1/4 d/dj (t h phi) dot g_i  ]
            IF (J == 3) THEN
                                                           ! Transverse shear rows are special. Eqn. (23a), (24a)
                                                           ! 1/4 d/dt (t h phi) = 1/4 h phi
              B(ROW, K+4) = B(ROW, K+4) + DIR_THICKNESS(GP) * PSH(GP) * (G(3,I) * DIRECTOR(2,GP) - G(2,I) * DIRECTOR(3,GP)) / FOUR
              B(ROW, K+5) = B(ROW, K+5) + DIR_THICKNESS(GP) * PSH(GP) * (G(1,I) * DIRECTOR(3,GP) - G(3,I) * DIRECTOR(1,GP)) / FOUR
              B(ROW, K+6) = B(ROW, K+6) + DIR_THICKNESS(GP) * PSH(GP) * (G(2,I) * DIRECTOR(1,GP) - G(1,I) * DIRECTOR(2,GP)) / FOUR
            ELSE
                                                           ! 1/4 d/dr (t h phi) = 1/4 t h dN/dr phi
                                                           ! 1/4 d/ds (t h phi) = 1/4 t h dN/ds phi
              B(ROW, K+4) = B(ROW, K+4) +                                                                                          &
                DIR_THICKNESS(GP) * T * DPSHG3(J,GP) * (G(3,I) * DIRECTOR(2,GP) - G(2,I) * DIRECTOR(3,GP)) / FOUR
              B(ROW, K+5) = B(ROW, K+5) +                                                                                          &
                DIR_THICKNESS(GP) * T * DPSHG3(J,GP) * (G(1,I) * DIRECTOR(3,GP) - G(3,I) * DIRECTOR(1,GP)) / FOUR
              B(ROW, K+6) = B(ROW, K+6) +                                                                                          &
                DIR_THICKNESS(GP) * T * DPSHG3(J,GP) * (G(2,I) * DIRECTOR(1,GP) - G(1,I) * DIRECTOR(2,GP)) / FOUR
            ENDIF


          ELSE
                                                           ! Zero row
            DO L=1,6
              B(ROW, K+L) = 0
            ENDDO

          ENDIF

        ENDDO
      ENDDO



      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION

      FUNCTION MITC_DETJ ( R, S, T )

! Calculates the Jacobian determinant in basic coordinates at a point in isoparametric coordinates.

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      REAL(DOUBLE)                    :: MITC_DETJ
      REAL(DOUBLE) , INTENT(IN)       :: R
      REAL(DOUBLE) , INTENT(IN)       :: S
      REAL(DOUBLE) , INTENT(IN)       :: T
      REAL(DOUBLE)                    :: G(3,3)            ! covariant basis vectors in basic coordinates
      REAL(DOUBLE)                    :: DUM1(3)

! **********************************************************************************************************************************

      CALL MITC_COVARIANT_BASIS( R, S, T, G )

      !DET(J) = G_R . (G_S x G_T)
      CALL CROSS(G(:,2), G(:,3), DUM1)
      MITC_DETJ = G(1,1)*DUM1(1) + G(2,1)*DUM1(2) + G(3,1)*DUM1(3)

      RETURN


! **********************************************************************************************************************************

      END FUNCTION MITC_DETJ

      FUNCTION MITC_ELASTICITY ()

! Returns the 6x6 elasticity matrix.
! In the element coordinate system for QUAD4.
! In the material coordinate system for QUAD8.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  EM, ET, EPROP
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE)                    :: MITC_ELASTICITY(6,6)


! **********************************************************************************************************************************

! Victor todo this is not very useful because different matrices are used in different places:
! MITC4+ stiffness:     membrane, bending+shear
! MITC4+ thermal load:  membrane
! MITC8 stiffness:      membrane+shear. Could be bending+shear because QUAD8 requires MID1=MID2.

                                                           ! Convert 2D material elasticity matrices to 3D.
      MITC_ELASTICITY(:,:) = ZERO
      MITC_ELASTICITY(1,1) = EM(1,1)
      MITC_ELASTICITY(1,2) = EM(1,2)
      MITC_ELASTICITY(1,4) = EM(1,3)
      MITC_ELASTICITY(2,2) = EM(2,2)
      MITC_ELASTICITY(2,4) = EM(2,3)
      MITC_ELASTICITY(4,4) = EM(3,3)
      MITC_ELASTICITY(5,5) = ET(2,2) * EPROP(3)
      MITC_ELASTICITY(5,6) = ET(2,1) * EPROP(3)
      MITC_ELASTICITY(6,6) = ET(1,1) * EPROP(3)

      DO I=2,6                                             ! Copy UT to LT because it's symmetric.
         DO J=1,I-1
            MITC_ELASTICITY(I,J) = MITC_ELASTICITY(J,I)
         ENDDO
      ENDDO

      RETURN

! **********************************************************************************************************************************

      END FUNCTION MITC_ELASTICITY

      SUBROUTINE MITC_INITIALIZE ()

! Initialize element variables in MITC_STUF.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, EPROP, XEL, BGRID, GRID_SNORM, TYPE, TE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS, GP_RS
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE

      INTEGER(LONG)                   :: GP                ! Element grid point number
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE)                    :: PSH(ELGP)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: TANGENT_R(3)
      REAL(DOUBLE)                    :: TANGENT_S(3)
      REAL(DOUBLE)                    :: NORMAL(3)         ! Intermediate value of some normal.

                                                           ! Normal defined by node positions at each node
      REAL(DOUBLE)                    :: MIDSURFACE_NORMAL(3,ELGP)
      REAL(DOUBLE)                    :: THICKNESS_FACTOR

! **********************************************************************************************************************************

! ----------------------------------------------------------------------------------------------------------------------------------
! R, S coordiantes of each node

      IF (TYPE(1:5) == 'QUAD4') THEN

         ! Bathe's node numbering relationship to R,S coordinates.
         GP_RS(1,1) =  ONE
         GP_RS(1,2) = -ONE
         GP_RS(1,3) = -ONE
         GP_RS(1,4) =  ONE

         GP_RS(2,1) =  ONE
         GP_RS(2,2) =  ONE
         GP_RS(2,3) = -ONE
         GP_RS(2,4) = -ONE

      ELSEIF (TYPE(1:5) == 'QUAD8') THEN

         GP_RS(1,1) = -ONE
         GP_RS(1,2) =  ONE
         GP_RS(1,3) =  ONE
         GP_RS(1,4) = -ONE
         GP_RS(1,5) =  ZERO
         GP_RS(1,6) =  ONE
         GP_RS(1,7) =  ZERO
         GP_RS(1,8) = -ONE

         GP_RS(2,1) = -ONE
         GP_RS(2,2) = -ONE
         GP_RS(2,3) =  ONE
         GP_RS(2,4) =  ONE
         GP_RS(2,5) = -ONE
         GP_RS(2,6) =  ZERO
         GP_RS(2,7) =  ONE
         GP_RS(2,8) =  ZERO

      ELSE

         WRITE(ERR,*) ' *ERROR: INCORRECT ELEMENT TYPE ', TYPE
         WRITE(F06,*) ' *ERROR: INCORRECT ELEMENT TYPE ', TYPE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF


! ----------------------------------------------------------------------------------------------------------------------------------
! Director vector at each node

      DO GP=1,ELGP

         CALL MITC_SHAPE_FUNCTIONS(GP_RS(1,GP), GP_RS(2,GP), PSH, DPSHG)

         TANGENT_R(:)=ZERO
         TANGENT_S(:)=ZERO

         ! TANGENT_R(r, s) = dX/dR = d/dR X = sum over nodes[ dN/dR X ]
         ! TANGENT_S(r, s) = dX/dS = d/dS X = sum over nodes[ dN/dS X ]
         DO I=1,ELGP
            DO J=1,3
               TANGENT_R(J) = TANGENT_R(J) + XEL(I,J) * DPSHG(1,I)
               TANGENT_S(J) = TANGENT_S(J) + XEL(I,J) * DPSHG(2,I)
            ENDDO
         ENDDO

         CALL CROSS(TANGENT_R, TANGENT_S, NORMAL)

         NORMAL = NORMAL / DSQRT(DOT_PRODUCT(NORMAL, NORMAL))

         MIDSURFACE_NORMAL(:,GP) = NORMAL

         NORMAL = GRID_SNORM(BGRID(GP),:)

                                                           ! Use the midsurface normal unless SNORM exists and it's
                                                           ! a linear element.
         IF (ANY(NORMAL /= ZERO) .AND. (TYPE(1:5) == 'QUAD4')) THEN
                                                           ! Transform SNORM from basic to XEL element coordinates.
            CALL MATMULT_FFF(TE, NORMAL, 3, 3, 1, DIRECTOR(:,GP))
         ELSE
            DIRECTOR(:,GP) = MIDSURFACE_NORMAL(:,GP)
         ENDIF

      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Thickness in the direction of the director vector at each node

      ! Thickness is treated as uniform.
      ! To allow grid point thicknesses, this should be interpolated to midside nodes.
      DO GP=1,ELGP

         ! If SNORM is used, the director vector may be different from the midsurface normal where thickness is defined
         ! by the user. Scale the thickness to make it in the direction of the director vector.
         THICKNESS_FACTOR = DOT_PRODUCT(DIRECTOR(:,GP), MIDSURFACE_NORMAL(:,GP))

                                                           ! Error if the angle betwen them is greater than ~=89 degrees.
                                                           ! 90 degrees would divide by zero and >90 degrees would be
                                                           ! inside-out.
         IF(THICKNESS_FACTOR < 0.01) THEN
            WRITE(ERR,*) ' *ERROR: SNORM IS TOO FAR FROM MIDSURFACE NORMAL'
            WRITE(F06,*) ' *ERROR: SNORM IS TOO FAR FROM MIDSURFACE NORMAL'
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DIR_THICKNESS(GP)   = EPROP(1) / THICKNESS_FACTOR

      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_INITIALIZE

      SUBROUTINE MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )


      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, TYPE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR

      USE QUADRILATERAL_SHAPE_FUNCTIONS, ONLY:  SHP2DQ
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: R,S               ! Isoparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: PSH(ELGP)         ! Shape functions
      REAL(DOUBLE) , INTENT(OUT)      :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: DUM(2,ELGP)
      REAL(DOUBLE)                    :: DUM2(ELGP)


! **********************************************************************************************************************************

                                                           ! Shape function derivatives at R,S
      IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8')) THEN

         CALL SHP2DQ ( 0, 0, ELGP, 'MITC_SHAPE_FUNCTIONS', '', 0, R, S, 'N', PSH, DPSHG )

                                                           ! Change node numbering to match Bathe's MITC4+ paper.
         IF (TYPE(1:5) == 'QUAD4') THEN
            DUM = DPSHG
            DPSHG(:,1) = DUM(:,3)
            DPSHG(:,2) = DUM(:,4)
            DPSHG(:,3) = DUM(:,1)
            DPSHG(:,4) = DUM(:,2)
            DUM2 = PSH
            PSH(1) = DUM2(3)
            PSH(2) = DUM2(4)
            PSH(3) = DUM2(1)
            PSH(4) = DUM2(2)
         ENDIF

      ELSE

        WRITE(ERR,*) ' *ERROR: INCORRECT ELEMENT TYPE ', TYPE
        WRITE(F06,*) ' *ERROR: INCORRECT ELEMENT TYPE ', TYPE
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_SHAPE_FUNCTIONS

      SUBROUTINE MITC_TRANSFORM_B ( TRANSFORM, B )

! Transform the strain-displacement matrix (tensor components) to a different basis.
! Equivalent to the sum
!    B_kl = A_ki Alj B~_ij

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      REAL(DOUBLE),  INTENT(INOUT)    :: B(6,6*ELGP)
      REAL(DOUBLE),  INTENT(IN)       :: TRANSFORM(3,3)
      REAL(DOUBLE)                    :: B_TRANSFORMED(6,6*ELGP)

      INTEGER(LONG)                   :: I,J,K,L           ! Tensor indices
      INTEGER(LONG)                   :: INDEX1(6)         ! Mapping of 6x1 vector index to 3x3 tensor first index.
      INTEGER(LONG)                   :: INDEX2(6)         ! Mapping of 6x1 vector index to 3x3 tensor second index.
      INTEGER(LONG)                   :: ROWS(3,3)         ! Mapping of 3x3 tensor indices to 6x1 vector index
      INTEGER(LONG)                   :: ROW
      INTEGER(LONG)                   :: IJ_ROW

      INTRINSIC                       :: DSQRT


! **********************************************************************************************************************************

      INDEX1 = (/ 1, 2, 3, 1, 2, 1 /)
      INDEX2 = (/ 1, 2, 3, 2, 3, 3 /)

      ROWS = RESHAPE((/ 1, 4, 6, 4, 2, 5, 6, 5, 3 /), SHAPE(ROWS))

      B_TRANSFORMED(:,:) = ZERO

      DO ROW=1,6
        K = INDEX1(ROW)
        L = INDEX2(ROW)
        DO I=1,3
          DO J=1,3
            IJ_ROW = ROWS(I,J)
            B_TRANSFORMED(ROW,:) = B_TRANSFORMED(ROW,:) + TRANSFORM(K,I) * TRANSFORM(L,J) * B(IJ_ROW,:)
          ENDDO
        ENDDO
      ENDDO

      B(:,:) = B_TRANSFORMED(:,:)

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC_TRANSFORM_B

      SUBROUTINE MITC_TRANSFORM_CONTRAVARIANT_TO_LOCAL ( R, S, T, B )

! Transform covariant strain components from the contravariant basis to the cartesian local basis.
!   ε_kl = ε~_ij (g^i dot e_k)(g^j dot e_l)

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP
      USE CONSTANTS_1, ONLY           :  ZERO, QUARTER, ONE, TWO, THREE

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T


      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: R, S, T           ! Isoparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: B(6, 6*ELGP)      ! Strain-displacement matrix
      REAL(DOUBLE)                    :: TRANSFORM(3,3)    ! Transformation matrix
      REAL(DOUBLE)                    :: G(3,3)            ! Array of 3 covariant basis vectors in basic coordinates
      REAL(DOUBLE)                    :: G_CONTRA(3,3)     ! Array of 3 contravariant basis vectors in basic coordinates
      REAL(DOUBLE)                    :: E(3,3)            ! Basis vectors

! **********************************************************************************************************************************


      CALL MITC_COVARIANT_BASIS( R, S, T, G )
      CALL MITC_CONTRAVARIANT_BASIS( G, G_CONTRA )
      E = MITC4_CARTESIAN_LOCAL_BASIS(R, S, T)

      ! The 3x3 transformation matrix A is defined by
      !    A_xy = g^y dot e_x
      ! or two transformation matrices multiplied
      CALL MATMULT_FFF_T (E, G_CONTRA, 3, 3, 3, TRANSFORM)

      CALL MITC_TRANSFORM_B( TRANSFORM, B)


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MITC_TRANSFORM_CONTRAVARIANT_TO_LOCAL

      FUNCTION MITC4_CARTESIAN_LOCAL_BASIS ( R, S, T )

! Reference [2]:
!  MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!     by Dvorkin and Bathe

! First index of the result (row) is a vector component in basic coordinates (x,y,z)
! Second index of the result (column) is basis vector (x_l, y_l, normal)

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      REAL(DOUBLE)                    :: MITC4_CARTESIAN_LOCAL_BASIS(3,3)
      REAL(DOUBLE) , INTENT(IN)       :: R
      REAL(DOUBLE) , INTENT(IN)       :: S
      REAL(DOUBLE) , INTENT(IN)       :: T
      REAL(DOUBLE)                    :: G(3,3)
      REAL(DOUBLE)                    :: X(3)
      REAL(DOUBLE)                    :: Y(3)
      REAL(DOUBLE)                    :: Z(3)

! **********************************************************************************************************************************

      CALL MITC_COVARIANT_BASIS( R, S, T, G )

      Z = G(:,3)
      Z = Z / DSQRT(DOT_PRODUCT(Z, Z))

      CALL CROSS(G(:,2), Z, X)
      X = X / DSQRT(DOT_PRODUCT(X, X))

      CALL CROSS(Z, X, Y)

      MITC4_CARTESIAN_LOCAL_BASIS(:,1) = X
      MITC4_CARTESIAN_LOCAL_BASIS(:,2) = Y
      MITC4_CARTESIAN_LOCAL_BASIS(:,3) = Z

      RETURN

! **********************************************************************************************************************************

      END FUNCTION MITC4_CARTESIAN_LOCAL_BASIS

   END MODULE MITC_KERNELS
