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

   MODULE VECTOR_GEOMETRY

   USE FULL_MATRIX_ALGEBRA, ONLY :  MATMULT_FFF, MATMULT_FFF_T
   USE DOF_ARRAY_INDEXING, ONLY :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
   USE QUADRILATERAL_SHAPE_FUNCTIONS, ONLY :  SHP2DQ

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CONVERT_VEC_COORD_SYS, CROSS, GEN_T0L, PARAM_CORDS_ACT_CORDS, PLANE_COORD_TRANS_21, PROJ_VEC_ONTO_PLANE, RIGID_BODY_DISP_MAT, SURFACE_FIT

   CONTAINS

      SUBROUTINE CONVERT_VEC_COORD_SYS ( MESSAG, INPUT_VEC, OUTPUT_VEC, NCID )

! Convert coordinate system of a G-set vector. The input vector is in global coords and the output vector is in one system for all
! grids. The input vector is first converted to basic coords and then to the final output system, if that system is not the basic
! system.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NCORD, NDOFG, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CORD, RCORD, GRID, GRID_ID, INV_GRID_SEQ

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CONVERT_VEC_COORD_SYS'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Text description of INPUT_VEC in case of undefined NCID
      CHARACTER( 1*BYTE)              :: CORD_FND          ! 'Y' if the internal coordinate system number for GCID, NCID is found

      INTEGER(LONG), INTENT(IN)       :: NCID              ! Actual coord system number. INPUT_VEC is to be transformed to this sys.
      INTEGER(LONG)                   :: AGRID             ! Actual grid number for the 6 components of vector being transformed
      INTEGER(LONG)                   :: GCID              ! Global actual coord system number for AGRID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where AGRID exists
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: JCORD             ! Internal coord system number for either GCID or NCID
      INTEGER(LONG)                   :: JFLD              ! Used in error message to indicate a coord sys ID undefined
      INTEGER(LONG)                   :: NUM_COMPS         ! No. displ components (1 for SPOINT, 6 for actual grid)


      REAL(DOUBLE), INTENT(IN)        :: INPUT_VEC(NDOFG)  ! G-set input vector to be transformed from global to NCID
      REAL(DOUBLE), INTENT(OUT)       :: OUTPUT_VEC(NDOFG) ! Transformed output vector
      REAL(DOUBLE)                    :: DUM_GVEC_T(6)     ! Vector for the 3 translation components of 1 grid
      REAL(DOUBLE)                    :: DUM_GVEC_R(6)     ! Vector for the 3 rotation    components of 1 grid
      REAL(DOUBLE)                    :: INPUT_GVEC_T(6)   ! Vector for the 3 translation components of 1 grid
      REAL(DOUBLE)                    :: INPUT_GVEC_R(6)   ! Vector for the 3 rotation    components of 1 grid
      REAL(DOUBLE)                    :: OUTPUT_GVEC_T(6)  ! Vector for the 3 translation components of 1 grid
      REAL(DOUBLE)                    :: OUTPUT_GVEC_R(6)  ! Vector for the 3 rotation    components of 1 grid
      REAL(DOUBLE)                    :: PHID              ! Dummy arg for subr GEN_T0L that is not used here
      REAL(DOUBLE)                    :: THETAD            ! Dummy arg for subr GEN_T0L that is not used here
      REAL(DOUBLE)                    :: T_0_GCID(3,3)     ! Coord transformation matrix from basic to GCID system
      REAL(DOUBLE)                    :: T_0_NCID(3,3)     ! Coord transformation matrix from basic to NCID system



! **********************************************************************************************************************************
! Set OUTPUT_VEC = INPUT_VEC in case no transformation is done (e.g. all grids have basic global and transformed system is basic)


      DO I=1,NDOFG
         OUTPUT_VEC(I) = INPUT_VEC(I)
      ENDDO

! Transform each grid of INPUT_VEC to basic coords in OUTPUT_VEC (unless it is in basic already)

      DO K=1,NGRID

         AGRID = GRID_ID(INV_GRID_SEQ(K))
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(K), NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS == 6) THEN                          ! Only 6 comp grids need transforming

            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
            GCID = GRID(GRID_ID_ROW_NUM,3)

            IF (GCID > 0) THEN

               CORD_FND = 'N'
j_do_1:        DO J=1,NCORD                                ! GCID should be a valid coord sys no. It was checked in subr GRID_PROC
                  IF (GCID == CORD(J,2)) THEN
                     JCORD = J
                     CORD_FND = 'Y'
                     EXIT j_do_1
                  ENDIF
               ENDDO j_do_1

               IF (CORD_FND == 'Y') THEN

                  DO J=1,3
                     INPUT_GVEC_T(J) = INPUT_VEC(6*(K-1)+J)
                     INPUT_GVEC_R(J) = INPUT_VEC(6*(K-1)+J+3)
                  ENDDO

                  CALL GEN_T0L ( GRID_ID_ROW_NUM, JCORD, THETAD, PHID, T_0_GCID )

                  CALL MATMULT_FFF ( T_0_GCID, INPUT_GVEC_T, 3, 3, 1, OUTPUT_GVEC_T )
                  CALL MATMULT_FFF ( T_0_GCID, INPUT_GVEC_R, 3, 3, 1, OUTPUT_GVEC_R )

                  DO J=1,3
                     OUTPUT_VEC(6*(K-1)+J)   = OUTPUT_GVEC_T(J)
                     OUTPUT_VEC(6*(K-1)+J+3) = OUTPUT_GVEC_R(J)
                  ENDDO


               ELSE                                        ! Should not get here (or pgm error) since GCID has been checked to exist

                  JFLD = 7
                  WRITE(ERR,910) GCID, JFLD, GRID(K,1), JFLD
                  WRITE(F06,910) GCID, JFLD, GRID(K,1), JFLD
                  CALL OUTA_HERE ( 'Y' )

               ENDIF

            ENDIF

         ENDIF

      ENDDO

! Transform OUTPUT_VEC to NCID if it is not basic

      IF (NCID /= 0) THEN

         CORD_FND = 'N'
j_do_2:  DO J=1,NCORD
           IF (NCID == CORD(J,2)) THEN
               JCORD = J
               CORD_FND = 'Y'
               EXIT j_do_2
            ENDIF
         ENDDO j_do_2

         IF (CORD_FND == 'Y') THEN

            DO K=1,NGRID

               AGRID = GRID_ID(INV_GRID_SEQ(K))
               CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(K), NUM_COMPS, SUBR_NAME )

               IF (NUM_COMPS == 6) THEN                          ! Only 6 comp grids need transforming

                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )


                  DO J=1,3
                     OUTPUT_GVEC_T(J) = OUTPUT_VEC(6*(K-1)+J)
                     OUTPUT_GVEC_R(J) = OUTPUT_VEC(6*(K-1)+J+3)
                  ENDDO

                  CALL GEN_T0L ( GRID_ID_ROW_NUM, JCORD, THETAD, PHID, T_0_NCID )

                  CALL MATMULT_FFF_T ( T_0_NCID, OUTPUT_GVEC_T, 3, 3, 1, DUM_GVEC_T )
                  CALL MATMULT_FFF_T ( T_0_NCID, OUTPUT_GVEC_R, 3, 3, 1, DUM_GVEC_R )

                  DO J=1,3
                     OUTPUT_VEC(6*(K-1)+J)   = DUM_GVEC_T(J)
                     OUTPUT_VEC(6*(K-1)+J+3) = DUM_GVEC_R(J)
                  ENDDO

               ENDIF

            ENDDO

         ELSE

            WRITE(ERR,912) NCID, MESSAG
            WRITE(F06,912) NCID, MESSAG
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ENDIF





      RETURN

! **********************************************************************************************************************************
  912 FORMAT(' *ERROR   912: COORDINATE SYSTEM ',I8,' DOES NOT EXIST. CANNOT TRANSFORM ',A                                         &
                    ,/,14X,' VECTOR FROM BASIC SYSTEM TO IT, SO VECTOR WILL BE LEFT IN BASIC COORD SYSTEM')

  910 FORMAT(' *ERROR   910: COORD SYSTEM',I8,' IN FIELD ',I2,' ON GRID ',I8,' (OR FROM GRDSET ENTRY FIELD ',I2,') IS UNDEFINED')

86954 format(' K, AGRID, J, INPUT_VEC(J), OUTPUT_VEC(J) = ',3i8,2(1es14.6))




! **********************************************************************************************************************************

      END SUBROUTINE CONVERT_VEC_COORD_SYS


      SUBROUTINE CROSS ( A, B, C )

! Cross product of 3x1 vectors: C = A (x) B

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CROSS'



      REAL(DOUBLE), INTENT(IN)        :: A(3)              ! Components of input  vector A
      REAL(DOUBLE), INTENT(IN)        :: B(3)              ! Components of input  vector B
      REAL(DOUBLE), INTENT(OUT)       :: C(3)              ! Components of output vector C



! **********************************************************************************************************************************
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CROSS


      SUBROUTINE GEN_T0L (RGRID_ROW, ICORD, THETAD, PHID, T0L )

! Generates 3x3 transformation matrix, T0L, to basic (0) coordinate system from another coordinate system (L) with internal coord
! system number ICORD at the grid whose RGRID data is at row RGRID_ROW in array RGRID. Thus, if U0 is a basic system vector
! and UL is a vector in coord system L, then:

!                       U0 = T0L*UL                      (1)

! In the RCORD database, we have transformations, T0P, that transform a vector to basic (0) from the principal
! directions (P) of some other coordinate system. If the local (L) system is rectangular then T0P is T0L and all we
! have to do is get T0P from RCORD. However, if the local system is cylindrical or spherical, then we have to also
! transform from the point at R, THETA, Z (cylindrical) or R, THETA, PHI (spherical) to the respective principal axes
! of that cylindrical or spherical system. Thus, we can write eqn (1) as:

!                       U0 = T0P*UP                      (2)
! and
!                       UP = TPL*UL                      (3)

! Substituting (3) into (2) and comparing with (1) it is seen that:

!                      T0L = T0P*TPL

! For a rectangular system, TPL is the identity matrix.

! For a cylindrical system TPL is the matrix in the relationship:

!              |UPx|   |cos(THETA) -sin(THETA) 0 ||Ur    |
!              |UPy| = |sin(THETA)  cos(THETA) 0 ||Utheta|
!              |UPz|   |     0           0     1 ||Uz    |

! where UPx, UPy, UPz are the displacements in the principal directiona of the cylindrical system and Ur, Utheta, Uz
! are the cylindrical system displacements at some cordinate location R, THETA, Z in the cylindrical system.

! For a spherical system TPL is the matrix in the relationship:

!  |UPx|   |sin(THETA)*cos(PHI) cos(THETA)*cos(PHI) -sin(PHI) ||Ur    |
!  |UPy| = |sin(THETA)*sin(PHI) cos(THETA)*sin(PHI)  cos(PHI) ||Utheta|
!  |UPz|   |    cos(THETA)          -sin(THETA)         0     ||Uphi  |

! where UPx, UPy, UPz are the displacements in the principal directions of the spherical system and Ur, Utheta, Uphi
! are the spherical system displacements at some coordinate location R, THETA, PHI in the spherical system.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, ONE80, PI
      USE IOUNT1, ONLY                :  WRT_ERR, f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE PARAMS, ONLY                :  EPSIL
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  RGRID, CORD, RCORD

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GEN_T0L'

      INTEGER(LONG), INTENT(IN)       :: RGRID_ROW         ! Row number in array RGRID where the RGRID data is stored for the grid
!                                                            point whose coord transformation we seek. Since RGRID (at this point)
!                                                            is sorted in grid numerical order, RGRID_ROW will be the same as what
!                                                            we would get from subr GET_ARRAY_ROW_NUM with actual grid equal to
!                                                            the grid whose transformation we seek.
      INTEGER(LONG), INTENT(IN)       :: ICORD             ! Internal coord ID for coord sys L
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices


      REAL(DOUBLE),  INTENT(OUT)      :: THETAD,PHID       ! Azimuth and elevation angles (deg) for cylindrical/spherical coord sys
      REAL(DOUBLE),  INTENT(OUT)      :: T0L(3,3)          ! 3 x 3 coord transformation matrix described above
      REAL(DOUBLE)                    :: CT,ST,CP,SP       ! SIN or COS of THETA or PHI
      REAL(DOUBLE)                    :: EPS1              ! EPSIL(1), a small number for comparison to zero
      REAL(DOUBLE)                    :: RAD               ! Rad from origin of sph or cyl coord. sys. ICORD to grid point RGRID_ROW
      REAL(DOUBLE)                    :: RADXY             ! Radius for calculating azimuth angle
      REAL(DOUBLE)                    :: RDEG              ! Number of degrees in a radian
      REAL(DOUBLE)                    :: THETA,PHI         ! Radian values for THETAD, PHID
      REAL(DOUBLE)                    :: T0P(3,3)          ! 3 x 3 transformation matrix described above
      REAL(DOUBLE)                    :: TPL(3,3)          ! 3 x 3 transformation matrix described above
      REAL(DOUBLE)                    :: XR0(3)            ! Array of relative coords of RGRID_ROW and origin of coord system ICORD
!                                                            in basic coords.
      REAL(DOUBLE)                    :: XRP(3)            ! Array of relative coords of RGRID_ROW and origin of coord system ICORD
!                                                            in local coord system L.

      INTRINSIC                       :: DASIN, DATAN2, DSIN, DCOS



! **********************************************************************************************************************************
! Initialize outputs

      THETAD = ZERO
      PHID   = ZERO

      DO I=1,3
         DO J=1,3
            T0L(I,J) = ZERO
         ENDDO
      ENDDO

! If ICORD is rectangular, then no more transformation is needed. If it is cylindrical or spherical, more transformation is needed
! and depends on the location of the grid at RGRID_ROW and the location of the origin of system ICORD (both locations in basic sys).
! The location of the grid at RGRID_ROW in basic is in RGRID(RGRID_ROW,1-3). The location of the origin of system ICORD in the basic
! system is in RCORD (first 3 words of the row ICORD).

! Set elements of T0P to the coord. sys. transf. matrix found.
! If rectangular, we are finished.

      EPS1   = EPSIL(1)
      RDEG   = ONE80/PI

      DO I=1,3
         DO J=1,3
            K = 3 + 3*(I-1) + J
            T0P(I,J) = RCORD(ICORD,K)
         ENDDO
      ENDDO
                                                           ! If 11 or 21 then rectangular
      IF ((CORD(ICORD,1) == 11) .OR. (CORD(ICORD,1) == 21)) THEN
         DO I=1,3
            DO J=1,3
               T0L(I,J) = T0P(I,J)
            ENDDO
         ENDDO

      ELSE                                                 ! Not Rectangular

! Get relative coordinates of RGRID_ROW and origin of ICORD in basic coordinate system.

         XR0(1) = RGRID(RGRID_ROW,1) - RCORD(ICORD,1)
         XR0(2) = RGRID(RGRID_ROW,2) - RCORD(ICORD,2)
         XR0(3) = RGRID(RGRID_ROW,3) - RCORD(ICORD,3)

! Transform relative coords from basic to the local coordinate system. Relative coords (3 of them) in local system are in
! XRP. Relative coords in basic system are in XR0. The transformation is: XRP = T0P(transpose)*XR0. (Note: XR0 = T0P*XRP)

         CALL MATMULT_FFF_T ( T0P, XR0, 3, 3, 1, XRP )

! Initialize TPL

         DO I=1,3
            DO J=1,3
               TPL(I,J) = ZERO
            ENDDO
         ENDDO

! **********************************************************************************************************************************
! Calculate T0L = T0P*TPL for cylindrical and spherical systems. Note that if the radius to the grid point (from the
! coord system) is zero, we define THETA = 0

         IF (CORD(ICORD,1) == 22) THEN                     ! Cylindrical
            RAD = DSQRT(XRP(1)*XRP(1) + XRP(2)*XRP(2))
            IF (DABS(RAD) < EPS1) THEN
               DO I=1,3
                  DO J=1,3
                     T0L(I,J) = T0P(I,J)
                     tpl(i,j) = zero
                     tpl(i,i) = one
                  ENDDO
               ENDDO
            ELSE
               THETA    = DATAN2(XRP(2),XRP(1))
               THETAD   = RDEG*THETA
               CT       = DCOS(THETA)
               ST       = DSIN(THETA)
               TPL(1,1) =  CT
               TPL(1,2) = -ST
               TPL(1,3) =  ZERO
               TPL(2,1) =  ST
               TPL(2,2) =  CT
               TPL(2,3) =  ZERO
               TPL(3,1) =  ZERO
               TPL(3,2) =  ZERO
               TPL(3,3) =  ONE
               CALL MATMULT_FFF ( T0P, TPL, 3, 3, 3, T0L )
            ENDIF
         ELSE IF (CORD(ICORD,1) == 23) THEN                ! Spherical
            RAD   = DSQRT(XRP(1)*XRP(1) + XRP(2)*XRP(2) + XRP(3)*XRP(3))
            RADXY = DSQRT(XRP(1)*XRP(1) + XRP(2)*XRP(2))
            IF ((DABS(RAD) < EPS1) .OR. (DABS(RADXY) < EPS1)) THEN
               THETA  = ZERO
               PHI    = ZERO
               CT     = ONE
               ST     = ZERO
               CP     = ONE
               SP     = ZERO
               DO I=1,3
                  DO J=1,3
                     T0L(I,J) = T0P(I,J)
                  ENDDO
               ENDDO
            ELSE
               THETA    =  DASIN(RADXY/RAD)
               PHI      =  DATAN2(XRP(2),XRP(1))
               CT       =  DCOS(THETA)
               ST       =  DSIN(THETA)
               CP       =  DCOS(PHI)
               SP       =  DSIN(PHI)
               PHID     =  RDEG*PHI
               THETAD   =  RDEG*THETA
               TPL(1,1) =  ST*CP
               TPL(1,2) =  CT*CP
               TPL(1,3) = -SP
               TPL(2,1) =  ST*SP
               TPL(2,2) =  CT*SP
               TPL(2,3) =  CP
               TPL(3,1) =  CT
               TPL(3,2) = -ST
               TPL(3,3) =  ZERO
               CALL MATMULT_FFF ( T0P, TPL, 3, 3, 3, T0L )
            ENDIF
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE GEN_T0L


      SUBROUTINE PARAM_CORDS_ACT_CORDS ( NROW, IORD, XEP, XEA )

! Converts element isoparametric coordinates to actual local element coordinates using bilinear shape functions.
! This is used in subr POLYNOM_FIT_STRE_STRN for extrapolating stress/strain values at the points at which the stress/strain
! matrices were calculated to element corner nodes

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  TYPE, XEL
      USE FILE_LIFECYCLE, ONLY        :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PARAM_CORDS_ACT_CORDS'

      INTEGER(LONG), INTENT(IN)       :: IORD              ! Gaussian integration order to be used in obtaining the PSH shape fcns
      INTEGER(LONG), INTENT(IN)       :: NROW              ! Number of rows in XEP, XEA


      REAL(DOUBLE), INTENT(IN)        :: XEP(NROW,3)       ! Parametric coords of NCOL points
      REAL(DOUBLE), INTENT(OUT)       :: XEA(NROW,3)       ! Actual local element coords corresponding to XEP



! **********************************************************************************************************************************
      IF     ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8')) THEN
         CALL GET_QUAD_COORDS
      ELSE
         Write(err,*) ' *ERROR      : Code not written in subr PARAM_CORDS_ACT_CORDS for', type
         Write(f06,*) ' *ERROR      : Code not written in subr PARAM_CORDS_ACT_CORDS for', type
       fatal_err = fatal_err + 1
         call outa_here ( 'y' )
      ENDIF



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################


      SUBROUTINE GET_QUAD_COORDS

! Parametric coords of points are in array XEP. They are obtained from the corner node coords in array XEL from:

!                                   XEA = PSH_MAT^T * XEL

! The terms in PSH_MAT are the shape functions from the PSH rows from subr SHP2DQ for each of the 4 XEP points.

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  WRT_BUG
      USE MODEL_STUF, ONLY            :  XEL, ELGP

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T
      USE QUADRILATERAL_SHAPE_FUNCTIONS, ONLY:  SHP2DQ
      IMPLICIT NONE

      INTEGER(LONG)                   :: J                   ! DO loop index

      REAL(DOUBLE)                    :: DPSHG(2,ELGP)       ! Derivatives of PSH wrt elem isopar coords (not used here).

                                                             ! 4x4 matrix used to calc Gauss pt coords from node coords
      REAL(DOUBLE)                    :: PSH_MAT(ELGP,IORD*IORD)

! **********************************************************************************************************************************

! The PSH_MAT columns are from subr SHP2DQ for each of the 4 XEP parametric coord points for the element.
! We want the XEA orderd in the same fashion as the element node coords in XEL (namely 1-2-3-4 clockwise around the element).

      CALL SHP2DQ ( 1, 1, ELGP, SUBR_NAME, ' ', IORD, XEP(1,1), XEP(1,2), 'Y', PSH_MAT(:,1), DPSHG )
      CALL SHP2DQ ( 2, 1, ELGP, SUBR_NAME, ' ', IORD, XEP(2,1), XEP(2,2), 'Y', PSH_MAT(:,2), DPSHG )
      CALL SHP2DQ ( 2, 2, ELGP, SUBR_NAME, ' ', IORD, XEP(3,1), XEP(3,2), 'Y', PSH_MAT(:,3), DPSHG )
      CALL SHP2DQ ( 1, 2, ELGP, SUBR_NAME, ' ', IORD, XEP(4,1), XEP(4,2), 'Y', PSH_MAT(:,4), DPSHG )


! Multiply shape functions by grid point coordinates to get Gauss point coordinates
! Only the first ELGP rows of XEL are used because it may have additional unused rows.
      CALL MATMULT_FFF_T ( PSH_MAT, XEL(1:ELGP,:), ELGP, IORD*IORD, 3, XEA )


! Debug output

      IF (WRT_BUG(5) == 1) THEN
         WRITE(BUG,*)
         WRITE(BUG,*) '  Parametric and element local coords. (Pt no., parametric coords (XEP), actual coords (XEA)'
         WRITE(BUG,*) '        Pt No.         Parametric Coords, XEP          XEA, coords in local elem coord system'
         WRITE(BUG,*) '                          XI            ET                X             Y               Z'
         WRITE(BUG,1001) '1',(XEP(1,J),J=1,3) , (XEA(1,J),J=1,3)
         WRITE(BUG,1001) '2',(XEP(2,J),J=1,3) , (XEA(2,J),J=1,3)
         WRITE(BUG,1001) '3',(XEP(3,J),J=1,3) , (XEA(3,J),J=1,3)
         WRITE(BUG,1001) '4',(XEP(4,J),J=1,3) , (XEA(4,J),J=1,3)
         WRITE(BUG,*)
      ENDIF

! **********************************************************************************************************************************
 1001 FORMAT(8X,A,10X,2(1ES15.6),2X,3(1ES15.6))

! **********************************************************************************************************************************

      END SUBROUTINE GET_QUAD_COORDS


      END SUBROUTINE PARAM_CORDS_ACT_CORDS


      SUBROUTINE PLANE_COORD_TRANS_21 ( THETA, T21, CALLING_SUBR )

! Creates a coordinate transformation matrix for a plane rotation of a vector in coordinate system 1, through an angle THETA, to
! a vector in coordinate system 2

!                             | U2 |   |  cos(THETA)  sin(THETA)  0 | | U1 |
!                             | V2 | = | -sin(THETA)  cos(THETA)  0 | | V1 |
!                             | W2 |   |      0           0       1 | | W1 |

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE, ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PLANE_COORD_TRANS_21'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this one



      REAL(DOUBLE), INTENT(IN)        :: THETA             ! Angle from x axis of system 1 to x axis of system 2
      REAL(DOUBLE), INTENT(OUT)       :: T21(3,3)          ! Transformation matrix which will transform a vector, U1, in coord sys
!                                                            1 to a vector, U2, in coord sys 2 (i.e. U2 = T21*U1)

      INTRINSIC                       :: DSIN, DCOS



! **********************************************************************************************************************************
! Row 1

      T21(1,1) =  DCOS( THETA )
      T21(1,2) =  DSIN( THETA )
      T21(1,3) =  ZERO

! Row 2

      T21(2,1) = -T21(1,2)
      T21(2,2) =  T21(1,1)
      T21(2,3) =  ZERO

! Row 3

      T21(3,1) =  ZERO
      T21(3,2) =  ZERO
      T21(3,3) =  ONE




      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE PLANE_COORD_TRANS_21


      SUBROUTINE PROJ_VEC_ONTO_PLANE ( VEC_A, VEC_B, VEC_C )

! Calcs the projection of a vector onto a plane given a normal to the plane (VEC_B) and the vector to be procected (VEC_A).
! Result is VEC_C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PROJ_VEC_ONTO_PLANE'

      INTEGER(LONG)                   :: I                 ! DO loop index


      REAL(DOUBLE) , INTENT(IN)       :: VEC_A(3)           ! Vector to be projected
      REAL(DOUBLE) , INTENT(IN)       :: VEC_B(3)           ! Vector normal to the plane onto which VEC_A is to be projected
      REAL(DOUBLE) , INTENT(OUT)      :: VEC_C(3)           ! Vector projection of VEC_A onto plane to which VEC_B is normal
      REAL(DOUBLE)                    :: VEC_DUM(3)         ! Dummy vector in the calc of VEC_C



! **********************************************************************************************************************************
      CALL CROSS ( VEC_A, VEC_B, VEC_DUM )
      CALL CROSS ( VEC_B, VEC_DUM, VEC_C )




      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE PROJ_VEC_ONTO_PLANE


       SUBROUTINE RIGID_BODY_DISP_MAT ( GRD_COORDS, REF_COORDS, RB_DISP )

! Generates a set of 6 rigid body displacement vectors for the 6 displacement components for one grid. The rigid body displacements
! are relative to REF_GRID and are in basic coords

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RIGID_BODY_DISP_MAT'

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE) , INTENT(IN)       :: GRD_COORDS(3)     ! Coords of grid point for which the RB matrix is to be formulated
      REAL(DOUBLE) , INTENT(IN)       :: REF_COORDS(3)     ! Coords of reference grid (grid about which the RB disps occur)
      REAL(DOUBLE) , INTENT(OUT)      :: RB_DISP(6,6)      ! The set of 6 RB displ vectors for the 6 disp components for GRID_NUM
      REAL(DOUBLE)                    :: XBAR              ! Basic X coordinate of GRID_NUM relative to REF_GRID
      REAL(DOUBLE)                    :: YBAR              ! Basic Y coordinate of GRID_NUM relative to REF_GRID
      REAL(DOUBLE)                    :: ZBAR              ! Basic Z coordinate of GRID_NUM relative to REF_GRID



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,6
         DO J=1,6
            RB_DISP(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      XBAR = GRD_COORDS(1) - REF_COORDS(1)
      YBAR = GRD_COORDS(2) - REF_COORDS(2)
      ZBAR = GRD_COORDS(3) - REF_COORDS(3)

! Calc 6 x 6 RB matrix

      DO I=1,6
         DO J=1,6
            RB_DISP(I,J) = ZERO
         ENDDO
         RB_DISP(I,I) = ONE
      ENDDO

      RB_DISP(1,5) =  ZBAR
      RB_DISP(1,6) = -YBAR

      RB_DISP(2,4) = -ZBAR
      RB_DISP(2,6) =  XBAR

      RB_DISP(3,4) =  YBAR
      RB_DISP(3,5) = -XBAR



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE RIGID_BODY_DISP_MAT


      SUBROUTINE SURFACE_FIT ( NUM_FITS, NUM_COEFFS, XI, YI, WI, XO, YO, WO, DEB, MESSAGE, OUNT, POLY_PCT_ERR, PCT_ERR_MAX, IERR )

! Fits a 2D surface polynomial to a set of input values. The subroutine is coded to fit a surface of up to 2nd order:

!              WF(X,Y) = B(0) + B(1)*X + B(2)*Y + B(3)*XY + B(4)*X^2 + B(5)*Y^2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE LSQ_MYSTRAN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SURFACE_FIT'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAGE               ! Message printed if debug output is printed

      INTEGER(LONG), INTENT(IN)       :: DEB                   ! If > 0 then DEB_SURFACE_FIT will be run
      INTEGER(LONG), INTENT(IN)       :: NUM_FITS              ! Number of data points to fit
      INTEGER(LONG), INTENT(IN)       :: NUM_COEFFS            ! Number of coefficients in the fitting polynomial
      INTEGER(LONG), INTENT(IN)       :: OUNT(2)               ! Output units for SURFACE_FIT
      INTEGER(LONG), INTENT(OUT)      :: IERR                  ! Error indicator
      INTEGER(LONG)                   :: I,J                   ! DO loop indices
      INTEGER(LONG)                   :: IFAULT                ! Return code from subr REGCF
      INTEGER(LONG), PARAMETER        :: MAX_COEFFS = 6        ! Maximum number of coefficients coded for ther polynomial fit


      LOGICAL                         :: LINDEP(MAX_COEFFS)

      REAL(DOUBLE), INTENT(IN)        :: WI(NUM_FITS)          ! Values of the function to fit at the input data points
      REAL(DOUBLE), INTENT(IN)        :: XI(NUM_FITS)          ! X coords of the input  data points
      REAL(DOUBLE), INTENT(IN)        :: YI(NUM_FITS)          ! Y coords of the input  data points
      REAL(DOUBLE), INTENT(IN)        :: XO(NUM_FITS)          ! X coords of the output data points
      REAL(DOUBLE), INTENT(IN)        :: YO(NUM_FITS)          ! Y coords of the output data points
      REAL(DOUBLE), INTENT(OUT)       :: POLY_PCT_ERR(NUM_FITS)! % difference between fitted & actual data (norm'd to max act data)
      REAL(DOUBLE), INTENT(OUT)       :: PCT_ERR_MAX           ! Max value from array POLY_PCT_ERR (max error at any input data pt)
      REAL(DOUBLE), INTENT(OUT)       :: WO(NUM_FITS)          ! Values of the function to fit at the output data points
      REAL(DOUBLE)                    :: DEN                   ! Intermediate value in a calculation
      REAL(DOUBLE)                    :: B(0:MAX_COEFFS-1)     ! Coefficients in the polynomial fit
      REAL(DOUBLE)                    :: RES(NUM_FITS)         ! Actual minus fitted data
      REAL(DOUBLE)                    :: RES_ABS               ! Max abs of the residual
      REAL(DOUBLE)                    :: WF(NUM_FITS)          ! Values at the input XI, YI coords of the fitted polynomial
      REAL(DOUBLE)                    :: WI_MAX                ! Max abs value of WI


      REAL(DOUBLE)                    :: XROW(0:MAX_COEFFS-1)! Array of XYI values for one data point

                                                             ! Values of XI, YI and their products to go into the poly fit eqn
      REAL(DOUBLE)                    :: XYI(NUM_FITS,0:MAX_COEFFS-1)

                                                             ! Values of XO, YO and their products to go into the poly fit eqn
      REAL(DOUBLE)                    :: XYO(NUM_FITS,0:MAX_COEFFS-1)

      REAL(DOUBLE), PARAMETER         :: WT = ONE            ! Parameter



! **********************************************************************************************************************************
      IERR = 0

! Make sure that we have not requested an impossible situation with NUM_COEFFS

      IF ((NUM_COEFFS < 1) .OR. (NUM_COEFFS > MAX_COEFFS)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(OUNT(1),958) SUBR_NAME, NUM_COEFFS, MAX_COEFFS
         IF (OUNT(2) /= OUNT(1)) THEN
            WRITE(OUNT(2),958) SUBR_NAME, NUM_COEFFS, MAX_COEFFS
         ENDIF
         IERR = 1
         RETURN
      ENDIF

! Initialize

      DO I=1,NUM_FITS
         DO J=1,MAX_COEFFS
            XYI(I,J-1) = ZERO
            XROW(J-1)  = ZERO
            LINDEP(J) = .FALSE.
         ENDDO
      ENDDO

! Calc polynomial fit to the input WI at its coords, XI, YI

      CALL STARTUP (NUM_COEFFS-1, .TRUE.)

      WI_MAX = ZERO
      DO I=1,NUM_FITS
         IF (DABS(WI(I)) > WI_MAX) WI_MAX = WI(I)
         IF (NUM_COEFFS >= 1) THEN
            XYI(I,0) = ONE
         ENDIF
         IF (NUM_COEFFS >= 2) THEN
            XYI(I,1) = XI(I)
         ENDIF
         IF (NUM_COEFFS >= 3) THEN
            XYI(I,2) = YI(I)
         ENDIF
         IF (NUM_COEFFS >= 4) THEN
            XYI(I,3) = XI(I)*XI(I)
         ENDIF
         IF (NUM_COEFFS >= 5) THEN
            XYI(I,4) = XI(I)*YI(I)
         ENDIF
         IF (NUM_COEFFS >= 6) THEN
            XYI(I,5) = YI(I)*YI(I)
         ENDIF
         XROW(1:NUM_COEFFS-1) = XYI(I,1:NUM_COEFFS-1)
         XROW(0) = ONE
         CALL INCLUD ( WT, XROW, WI(I) )
      ENDDO

      CALL TOLSET ()

      CALL SING ( LINDEP, IFAULT )

      CALL REGCF ( B, NUM_COEFFS, IFAULT )

! Calc percent errors in polynomial fit to input data

      RES_ABS = ZERO
      DO I=1,NUM_FITS
         WF(I) = ZERO
         DO J=0,NUM_COEFFS-1
            WF(I) = WF(I) + B(J)*XYI(I,J)
         ENDDO
         RES(I) = WI(I) - WF(I)
         IF (DABS(RES(I)) > RES_ABS) THEN
            RES_ABS = DABS(RES(I))
         ENDIF
      ENDDO

      PCT_ERR_MAX = ZERO
      DO I=1,NUM_FITS
         IF (DABS(WI_MAX) > 0.0) THEN
            DEN = DABS(WI_MAX)
         ELSE
            DEN = ONE
         ENDIF
         POLY_PCT_ERR(I) = 100*RES(I)/DEN
         IF (DABS(PCT_ERR_MAX) < DABS(POLY_PCT_ERR(I))) THEN
            PCT_ERR_MAX = POLY_PCT_ERR(I)
         ENDIF
      ENDDO

! Calc fit to output coords, XO, YO

      DO I=1,NUM_FITS
         IF (NUM_COEFFS >= 1) THEN
            XYO(I,0)   = ONE
         ENDIF
         IF (NUM_COEFFS >= 2) THEN
            XYO(I,1)   = XO(I)
         ENDIF
         IF (NUM_COEFFS >= 3) THEN
            XYO(I,2)   = YO(I)
         ENDIF
         IF (NUM_COEFFS >= 4) THEN
            XYO(I,3)   = XO(I)*XO(I)
         ENDIF
         IF (NUM_COEFFS >= 5) THEN
            XYO(I,4)   = XO(I)*YO(I)
         ENDIF
         IF (NUM_COEFFS >= 6) THEN
            XYO(I,5)   = YO(I)*YO(I)
         ENDIF
         WO(I) = ZERO
         DO J=0,NUM_COEFFS-1
            WO(I) = WO(I) + B(J)*XYO(I,J)
         ENDDO
      ENDDO

      IF (DEB > 0) CALL DEB_SURFACE_FIT



      RETURN

! **********************************************************************************************************************************
  958 FORMAT(' *ERROR   958: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE VX VECTOR FOR ',A,' ELEMENT ',I8,' WAS LEFT UNSORTED. IT MUST BE SORTED TO DETERMINE VY, VZ')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEB_SURFACE_FIT

      IMPLICIT NONE

      INTEGER(LONG)                   :: II                ! DO loop index

! **********************************************************************************************************************************
      WRITE(OUNT(2),2000)
      WRITE(OUNT(2),'(A)') ' DEBUG 175 results from subr SURFACE_FIT'//MESSAGE//':'
      WRITE(OUNT(2),*)     ' --------------------------------------'
      WRITE(OUNT(2),2001) IFAULT
      WRITE(OUNT(2),'(A,1ES10.2,A)') ' PCT_ERR_MAX = ',PCT_ERR_MAX,'%'
      WRITE(OUNT(2),*)

      IF      (NUM_COEFFS == 1) THEN
         WRITE(OUNT(2),2011) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2021) (B(I), I=0,NUM_COEFFS-1)
      ELSE IF (NUM_COEFFS == 2) THEN
         WRITE(OUNT(2),2012) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2022) (B(I), I=0,NUM_COEFFS-1)
      ELSE IF (NUM_COEFFS == 3) THEN
         WRITE(OUNT(2),2013) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2023) (B(I), I=0,NUM_COEFFS-1)
      ELSE IF (NUM_COEFFS == 4) THEN
         WRITE(OUNT(2),2014) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2024) (B(I), I=0,NUM_COEFFS-1)
      ELSE IF (NUM_COEFFS == 5) THEN
         WRITE(OUNT(2),2014) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2025) (B(I), I=0,NUM_COEFFS-1)
      ELSE IF (NUM_COEFFS == 6) THEN
         WRITE(OUNT(2),2014) NUM_FITS, NUM_COEFFS
         WRITE(OUNT(2),2026) (B(I), I=0,NUM_COEFFS-1)
      ENDIF
      WRITE(OUNT(2),*)

      IF (RES_ABS > 0) THEN
         WRITE(OUNT(2),2031)
      ELSE
         WRITE(OUNT(2),2032)
      ENDIF

      DO II=1,NUM_FITS
         IF (DABS(WI_MAX) > 0.0) THEN
            WRITE(OUNT(2),2051) II, XI(II), YI(II), WI(II), WF(II), RES(II), POLY_PCT_ERR(II)
         ELSE
            WRITE(OUNT(2),2052) II, XI(II), YI(II), WI(II), WF(II), RES(II)
         ENDIF
      ENDDO

      WRITE(OUNT(2),*)
      WRITE(OUNT(2),2061)
      DO II=1,NUM_FITS
         WRITE(OUNT(2),2062) II, XO(II), YO(II), WO(II)
      ENDDO

! **********************************************************************************************************************************
 2000 FORMAT(' ******************************************************************************************************************',&
             '******************')

 2001 FORMAT(' Return code from polynomial fit: IFAULT = ',I3,/)

 2011 FORMAT(' Fitted quadratic surface for ',I2,' data points using a ',I2,'st order polynomial:',/,                              &
             ' -------------------------------------------------------------------------')

 2012 FORMAT(' Fitted quadratic surface for ',I2,' data points using a ',I2,'nd order polynomial:',/,                              &
             ' -------------------------------------------------------------------------')

 2013 FORMAT(' Fitted quadratic surface for ',I2,' data points using a ',I2,'rd order polynomial:',/,                              &
             ' -------------------------------------------------------------------------')

 2014 FORMAT(' Fitted quadratic surface for ',I2,' data points using a ',I2,'th order polynomial:',/,                              &
             ' -------------------------------------------------------------------------')

 2021 FORMAT('  Y = B0                                              ',//,                                                          &
             '      where: B0 = ',1ES14.6,/)

 2022 FORMAT('  Y = B0 + B1*X1                                      ',//,                                                          &
             '      where: B0 = ',1ES14.6,/,                                                                                       &
             '             B1 = ',1ES14.6,/)

 2023 FORMAT('  Y = B0 + B1*X1 + B2*X2                              ',//,                                                          &
             '      where: B0 = ',1ES14.6,/,                                                                                       &
             '             B1 = ',1ES14.6,/,                                                                                       &
             '             B2 = ',1ES14.6,/)

 2024 FORMAT('  Y = B0 + B1*X1 + B2*X2 + B3*X1^2                    ',//,                                                          &
             '      where: B0 = ',1ES14.6,/,                                                                                       &
             '             B1 = ',1ES14.6,/,                                                                                       &
             '             B2 = ',1ES14.6,/,                                                                                       &
             '             B3 = ',1ES14.6,/)

 2025 FORMAT('  Y = B0 + B1*X1 + B2*X2 + B3*X1^2 +B4*X1*Y1          ',//,                                                          &
             '      where: B0 = ',1ES14.6,/,                                                                                       &
             '             B1 = ',1ES14.6,/,                                                                                       &
             '             B2 = ',1ES14.6,/,                                                                                       &
             '             B3 = ',1ES14.6,/,                                                                                       &
             '             B4 = ',1ES14.6,/)

 2026 FORMAT('  Y = B0 + B1*X1 + B2*X2 +B3*X1^2 + B4*X1*Y1 + B5*Y2^2',//,                                                          &
             '      where: B0 = ',1ES14.6,/,                                                                                       &
             '             B1 = ',1ES14.6,/,                                                                                       &
             '             B2 = ',1ES14.6,/,                                                                                       &
             '             B3 = ',1ES14.6,/,                                                                                       &
             '             B4 = ',1ES14.6,/,                                                                                       &
             '             B5 = ',1ES14.6,/)


 2031 FORMAT(' Comparison of actual data to fitted curve:',/,' -----------------------------------------',/,                       &
             '     I      XI(I)         YI(I)         WI(I)         WF(I)       WI - WF       Err')

 2032 FORMAT(' Comparison of actual data to fitted curve:',/,' -----------------------------------------',/,                       &
             '     I      XI(I)         YI(I)         WI(I)         WF(I)       WI - WF')

 2051 FORMAT(I6,5(1ES14.6),F8.3,'%')

 2052 FORMAT(I6,5(1ES14.6))

 2061 FORMAT(' Polynomial values at output data points:',/,' -----------------------------------------',/,                       &
             '     I      XO(I)         YO(I)         WO(I)')

 2062 FORMAT(I6,3(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE DEB_SURFACE_FIT

      END SUBROUTINE SURFACE_FIT


   END MODULE VECTOR_GEOMETRY
