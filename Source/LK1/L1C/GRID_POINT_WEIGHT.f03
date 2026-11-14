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

! End MIT license text

   MODULE GRID_POINT_WEIGHT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GPWG, GPWG_USERIN, RB_DISP_MATRIX_PROC

   CONTAINS

      SUBROUTINE GPWG ( WHICH )

! Generates rigid body mass properties for the finite element model

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2, SC1, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_ME_BIT, IBIT, MBUG, NCONM2, NCORD, NELE, NGRID, SOL_NAME, WARN_ERR
      USE PARAMS, ONLY                :  EPSIL, GRDPNT, MEFMGRID, MEFMLOC, SUPWARN, WTMASS
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, CONM2, CORD, CAN_ELEM_TYPE_OFFSET, ELDT, ELGP, NUM_EMG_FATAL_ERRS,          &
                                         GRID, GRID_ID, MCG, ME, MEFFMASS_CALC, MEFM_RB_MASS,                                      &
                                         MODEL_MASS, MODEL_IXX, MODEL_IYY, MODEL_IZZ, MODEL_XCG, MODEL_YCG, MODEL_ZCG,             &
                                         OFFDIS, OFFSET, PLY_NUM, RCONM2, RGRID, TYPE, USERIN_RBM0

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE EMG_MOD, ONLY               :  EMG
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE, WRITE_ITABLE, WRITE_TABLE_HEADER
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GPWG'
      CHARACTER(12*BYTE), INTENT(IN)  :: WHICH             ! Whether to get mass props for
!                                                            (1) OA model
!                                                            (2) residual str, or
!                                                            (3) USERIN elems

      CHARACTER(1*BYTE)               :: OPT(6)            ! Option flags for what to calculate when subr EMG is called

      INTEGER(LONG)                   :: ACID_G            ! Actual global coord sys ID for a grid
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable resulting from an IAND operation
      INTEGER(LONG)                   :: ICID_G            ! Internal coord sys ID corresponding to actual coord sys ACID_G
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: JDOF              ! Array index used in getting mass terms from the elem mass matrix, ME
      INTEGER(LONG)                   :: GRID_NUM          ! An actual grid ID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: INFO        = 0   ! An output from subr GPWG_PMOI, called herein
      INTEGER(LONG)                   :: NUM_COMPS         ! Either 6 or 1 depending on whether grid is a physical grid or a SPOINT
      INTEGER(LONG)                   :: REFPNT            ! Reference point for GPWG calc (either GRDPNT of MEFMGRID)
      INTEGER(LONG)                   :: REFPNT_DEF        ! Default value of GRDPNT


      REAL(DOUBLE)                    :: BASIC_OFF(3)      ! Offsets of an element at a grid in basic coords
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: DX                ! X offset of a mass from its' grid
      REAL(DOUBLE)                    :: DY                ! Y offset of a mass from its' grid
      REAL(DOUBLE)                    :: DZ                ! Z offset of a mass from its' grid
      REAL(DOUBLE)                    :: GLOBAL_OFF(3)     ! Offsets of an element at a grid in global coords
      REAL(DOUBLE)                    :: MASS              ! Total mass
      REAL(DOUBLE)                    :: M0                ! An intermediate variable used in calc model mass props
      REAL(DOUBLE)                    :: MOI1(3,3)         ! Moments of inertia (several diff interps during exec of this subr)
      REAL(DOUBLE)                    :: S(3,3)            ! Tranformation matrix from basic to principal axes of inertia
      REAL(DOUBLE)                    :: IS(3,3)           ! Moments of Inertia relative to the CG
      REAL(DOUBLE)                    :: MX                ! First moment of MASS about X axis
      REAL(DOUBLE)                    :: MY                ! First moment of MASS about Y axis
      REAL(DOUBLE)                    :: MZ                ! First moment of MASS about Z axis
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: Q(3,3)            ! Transformation between S and Q axes
                                                           ! Output from subr GPWG_PMOI, called herein (transform to princ dir's)
      REAL(DOUBLE)                    :: RB_MASS_BASIC(6,6)! MO: 6x6 rigid body mass matrix about ref point in basic coords
      REAL(DOUBLE)                    :: T0G(3,3)          ! Coord transformation matrix from basic to global for a grid
      REAL(DOUBLE)                    :: TRANS(3,3)        ! Transfer terms when MOI's about c.g. ard calc'd from MOI's about XREF
      REAL(DOUBLE)                    :: XB(3)             ! Basic coord diffs bet c.g. and XREF in X, Y, Z directions
      REAL(DOUBLE)                    :: XD(3)             ! Basic coord diffs bet a mass (at it's c.m.) and XREF in X, Y, Z dirs
      REAL(DOUBLE)                    :: XREF(3)           ! GRDPNT basic coords (or origin of basic sys if GRDPNT doesn't exist)
      REAL(DOUBLE)                    :: X2                ! XD(1)*XD(1)
      REAL(DOUBLE)                    :: Y2                ! XD(2)*XD(2)
      REAL(DOUBLE)                    :: Z2                ! XD(3)*XD(3)
      REAL(DOUBLE)                    :: XY                ! XD(1)*XD(2)
      REAL(DOUBLE)                    :: XZ                ! XD(1)*XD(3)
      REAL(DOUBLE)                    :: YZ                ! XD(2)*XD(3)

      ! op2
      INTEGER(LONG)                   :: ITABLE            ! the op2 subtable counter
      INTEGER(LONG)                   :: ANALYSIS_CODE     ! the result type
      INTEGER(LONG), PARAMETER        :: TABLE_CODE = 0    ! is this right?
      CHARACTER(LEN=8)                :: TABLE_NAME         ! Name of the op2 table that we're writing

      INTRINSIC                       :: DABS
      INTRINSIC                       :: IAND



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      EPS1 = EPSIL(1)

      ! Set defaults in case GRDPNT grid cannot be found or is input as 0 (basic origin)
      REFPNT_DEF = 0
      XREF(1)    = ZERO
      XREF(2)    = ZERO
      XREF(3)    = ZERO

      IF ((SOL_NAME(1:5) == 'MODES') .AND. (MEFFMASS_CALC == 'Y')) THEN
         REFPNT = MEFMGRID
      ELSE
         REFPNT = GRDPNT
      ENDIF

      ! Get reference point coordinates in basic system for the reference point
      IF (REFPNT /= -1) THEN
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, REFPNT, GRID_ID_ROW_NUM )
         IF (GRID_ID_ROW_NUM /= -1) THEN                   ! REFPNT is a grid point in the model, so get its basic coords
            XREF(1) = RGRID(GRID_ID_ROW_NUM,1)
            XREF(2) = RGRID(GRID_ID_ROW_NUM,2)
            XREF(3) = RGRID(GRID_ID_ROW_NUM,3)
         ELSE                                              ! REFPNT was not 0 and is not a grid number in the model
            IF (REFPNT /= 0) THEN
               WARN_ERR  = WARN_ERR + 1
               WRITE(ERR,1402) REFPNT
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1402) REFPNT
               ENDIF
            ENDIF
            REFPNT = REFPNT_DEF
         ENDIF
      ENDIF

! Generate total mass, first and second moments by summing up mass terms. XD(i) are components of vector from
! ref point to a mass point. At this time, mass units are input units without PARAM WTMASS which is what we want for
! the grid point weight generator. Later the mass will be converted by multiplying by WTMASS.

      MASS  = ZERO
      MX    = ZERO
      MY    = ZERO
      MZ    = ZERO
      DO I=1,3
         DO J=1,3
            MOI1(I,J) = ZERO
         ENDDO
      ENDDO
      XB(1) = ZERO
      XB(2) = ZERO
      XB(3) = ZERO

      ! First process element mass terms
      OPT(1) = 'Y'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      IERROR = 0
      CALL COUNTER_INIT('    Working on element', NELE)
elems:DO I = 1,NELE
         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'N' )       ! 'N' means do not write to BUG file

         IF (NUM_EMG_FATAL_ERRS == 0) THEN

            IF (TYPE /= 'USERIN  ') THEN

res_or_oa:     IF ((WHICH(1:8) == 'OA MODEL') .OR. (WHICH == 'RESIDUAL STR')) THEN

elgp_do:          DO J=1,ELGP

                     BASIC_OFF(1) = ZERO
                     BASIC_OFF(2) = ZERO
                     BASIC_OFF(3) = ZERO
                     IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
                        IF (OFFSET(J) == 'Y') THEN               ! Calc (or set to 0) elem offset in basic (BASIC_OFF) at this grid
                           GLOBAL_OFF(1) = OFFDIS(J,1)
                           GLOBAL_OFF(2) = OFFDIS(J,2)
                           GLOBAL_OFF(3) = OFFDIS(J,3)
                           CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(J), GRID_ID_ROW_NUM )
                           ACID_G = GRID(GRID_ID_ROW_NUM,3)
                           IF (ACID_G /= 0) THEN
k_do111:                      DO K=1,NCORD
                                 IF (ACID_G == CORD(K,2)) THEN
                                    ICID_G = K
                                    EXIT k_do111
                                 ENDIF
                              ENDDO k_do111
                              CALL GEN_T0L ( GRID_ID_ROW_NUM, ICID_G, THETAD, PHID, T0G )
                              CALL MATMULT_FFF ( T0G, GLOBAL_OFF, 3, 3, 1, BASIC_OFF )
                           ELSE
                              BASIC_OFF(1) = GLOBAL_OFF(1)
                              BASIC_OFF(2) = GLOBAL_OFF(2)
                              BASIC_OFF(3) = GLOBAL_OFF(3)
                           ENDIF
                        ENDIF
                     ENDIF

                     CALL GET_GRID_NUM_COMPS ( BGRID(J), NUM_COMPS, SUBR_NAME )
                     JDOF  = NUM_COMPS*(J - 1) + 1
                     M0    = ME(JDOF,JDOF)
                     MASS  = MASS + M0
                     XD(1) = RGRID(BGRID(J),1) + BASIC_OFF(1) - XREF(1)
                     XD(2) = RGRID(BGRID(J),2) + BASIC_OFF(2) - XREF(2)
                     XD(3) = RGRID(BGRID(J),3) + BASIC_OFF(3) - XREF(3)
                     MX    = MX + M0*XD(1)
                     MY    = MY + M0*XD(2)
                     MZ    = MZ + M0*XD(3)
                     X2    = XD(1)*XD(1)
                     Y2    = XD(2)*XD(2)
                     Z2    = XD(3)*XD(3)
                     XY    = XD(1)*XD(2)
                     XZ    = XD(1)*XD(3)
                     YZ    = XD(2)*XD(3)
                     MOI1(1,1) = MOI1(1,1) + M0*(Y2 + Z2)
                     MOI1(2,2) = MOI1(2,2) + M0*(X2 + Z2)
                     MOI1(3,3) = MOI1(3,3) + M0*(X2 + Y2)
                     MOI1(2,1) = MOI1(2,1) - M0*XY
                     MOI1(3,1) = MOI1(3,1) - M0*XZ
                     MOI1(3,2) = MOI1(3,2) - M0*YZ

                  ENDDO elgp_do

               ENDIF res_or_oa

            ELSE                                           ! TYPE is USERIN.

userin:        IF ((WHICH(1:8) == 'OA MODEL') .OR. (WHICH(1:6) == 'USERIN')) THEN

                  M0   = USERIN_RBM0(1,1)
                  MASS = MASS + M0
                  IF (DABS(M0) > EPS1) THEN
                     XD(1) = USERIN_RBM0(2,6)/M0 !- XREF(1)
                     XD(2) = USERIN_RBM0(3,4)/M0 !- XREF(2)
                     XD(3) = USERIN_RBM0(1,5)/M0 !- XREF(3)
                     MX    = MX + M0*XD(1)
                     MY    = MY + M0*XD(2)
                     MZ    = MZ + M0*XD(3)

                     ! NOTE: USERIN_RBM0 was calc'd rel to GRDPNT grid point
                     MOI1(1,1) = MOI1(1,1) + USERIN_RBM0(4,4)
                     MOI1(2,2) = MOI1(2,2) + USERIN_RBM0(5,5)
                     MOI1(3,3) = MOI1(3,3) + USERIN_RBM0(6,6)
                     MOI1(2,1) = MOI1(2,1) + USERIN_RBM0(5,4)
                     MOI1(3,1) = MOI1(3,1) + USERIN_RBM0(6,4)
                     MOI1(3,2) = MOI1(3,2) + USERIN_RBM0(6,5)
                  ELSE
                     XD(1) = ZERO
                     XD(2) = ZERO
                     XD(3) = ZERO
                  ENDIF

               ENDIF userin

            ENDIF

         ELSE

            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE

         ENDIF
         CALL COUNTER_PROGRESS(I)

      ENDDO elems

      WRITE(SC1,*) CR13

      OPT(1) = 'N'                                         ! Reset subr EMG option value:

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! Errors from subr EMG, so quit
      ENDIF

! Now process mass in CONM's and calc MOI's about XREF in basic coords. The CONM2 data must be values at the mass in basic coords
! so that this subr must be run after subr CONM2_PROC_1 and before subr CONM2_PROC_2.
! We process all CONM2's, and add values for repeated grids (i.e. if model has 2 diff CONM2's at the same grid the mass props
! for that grid will be the sum of the mass props for each grid)

      DO I=1,NCONM2

         GRID_NUM = CONM2(I,2)
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, GRID_ID_ROW_NUM )
         MASS  = MASS + RCONM2(I,1)                        ! Add all masses

         ! DX, DY, DZ are offsets from the grid to the i-th mass in basic coords
         DX        = RCONM2(I,2)
         DY        = RCONM2(I,3)
         DZ        = RCONM2(I,4)

         ! XD(i) are distances from XREF to the i-th mass
         XD(1)     = RGRID(GRID_ID_ROW_NUM,1) + DX - XREF(1)
         XD(2)     = RGRID(GRID_ID_ROW_NUM,2) + DY - XREF(2)
         XD(3)     = RGRID(GRID_ID_ROW_NUM,3) + DZ - XREF(3)

         M0        = RCONM2(I,1)

         ! First moments about XREF in basic coords
         MX        = MX + M0*XD(1)
         MY        = MY + M0*XD(2)
         MZ        = MZ + M0*XD(3)

         ! Terms needed below for MOI calc
         X2        = XD(1)*XD(1)
         Y2        = XD(2)*XD(2)
         Z2        = XD(3)*XD(3)
         XY        = XD(1)*XD(2)
         XZ        = XD(1)*XD(3)
         YZ        = XD(2)*XD(3)

         ! Following are MOI's about XREF in basic coords
         MOI1(1,1) = MOI1(1,1) + RCONM2(I, 5) + M0*(Y2 + Z2)
         MOI1(2,2) = MOI1(2,2) + RCONM2(I, 7) + M0*(X2 + Z2)
         MOI1(3,3) = MOI1(3,3) + RCONM2(I,10) + M0*(X2 + Y2)
         MOI1(2,1) = MOI1(2,1) - RCONM2(I, 6) - M0*XY      ! RCONM2 has the products of inertia I21, I31, I32 of the CONM2
         MOI1(3,1) = MOI1(3,1) - RCONM2(I, 8) - M0*XZ      ! entry; the tensor terms are their negatives
         MOI1(3,2) = MOI1(3,2) - RCONM2(I, 9) - M0*YZ

      ENDDO

      ! Set other 3 products of inertia based on symmetry
      MOI1(1,2) = MOI1(2,1)
      MOI1(1,3) = MOI1(3,1)
      MOI1(2,3) = MOI1(3,2)

      ! backup MOI1 to create MO
      !MO(1,1) = MOI1(1,1)
      !MO(2,2) = MOI1(2,2)
      !MO(3,3) = MOI1(3,3)
      !
      !MO(1,2) = MOI1(1,2)
      !MO(1,3) = MOI1(1,3)
      !MO(2,3) = MOI1(2,3)
      !
      !MO(2,1) = MOI1(1,2)
      !MO(3,1) = MOI1(1,3)
      !MO(3,2) = MOI1(2,3)

      ! XB(I) are components of distance from reference point, XREF, to c.g.
      IF (DABS(MASS) > EPS1) THEN
         XB(1) = MX/MASS
         XB(2) = MY/MASS
         XB(3) = MZ/MASS
      ENDIF

      ! M0 - RB_MASS_BASIC: is 6x6 rigid body mass matrix about ref point in basic coords
      DO I=1,6                                             ! Init
         DO J=1,6
            RB_MASS_BASIC(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,3
         RB_MASS_BASIC(I,I) = MASS
      ENDDO

      RB_MASS_BASIC(1,4) =  ZERO         ;   RB_MASS_BASIC(1,5) =  MASS*XB(3)   ;   RB_MASS_BASIC(1,6) = -MASS*XB(2)
      RB_MASS_BASIC(2,4) = -MASS*XB(3)   ;   RB_MASS_BASIC(2,5) =  ZERO         ;   RB_MASS_BASIC(2,6) =  MASS*XB(1)
      RB_MASS_BASIC(3,4) =  MASS*XB(2)   ;   RB_MASS_BASIC(3,5) = -MASS*XB(1)   ;   RB_MASS_BASIC(3,6) =  ZERO

      DO I=4,6
         DO J=4,6
            RB_MASS_BASIC(I,J) = MOI1(I-3,J-3)
         ENDDO
      ENDDO

      ! Calc MEFM_RB_MASS here for MODES. It is calc'd in CALC_MRRcb for CB
      IF ((SOL_NAME(1:5) == 'MODES') .AND. (MEFFMASS_CALC == 'Y')) THEN
         DO I=1,6
            DO J=1,6
               MEFM_RB_MASS(I,J) = RB_MASS_BASIC(I,J)
            ENDDO
         ENDDO
      ENDIF

      ! Output results so far
      IF (REFPNT >= 0) THEN
         IF      (WHICH(1:12) == 'RESIDUAL STR' ) THEN
            WRITE(F06,902)
         ELSE IF (WHICH(1: 8) == 'OA MODEL'     ) THEN
            WRITE(F06,903)
         ENDIF

         IF (REFPNT == 0) THEN
            WRITE(F06,1001)
         ELSE IF (REFPNT > 0) THEN
            WRITE(F06,1002) REFPNT
            IF ((SOL_NAME(1:5) == 'MODES') .AND. (MEFFMASS_CALC == 'Y')) THEN
               IF (GRDPNT /= MEFMGRID) THEN
                  WRITE(F06,1013) MEFMGRID
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (REFPNT >= 0) THEN
         WRITE(F06,1004) MASS
         WRITE(F06,1005) (XB(I),I=1,3)

         WRITE(F06,1021)
         WRITE(F06,1901)
         DO I=1,3
            WRITE(F06,1022) (RB_MASS_BASIC(I,J),J=1,6)
         ENDDO
         WRITE(F06,1902)
         DO I=4,6
            WRITE(F06,1022) (RB_MASS_BASIC(I,J),J=1,6)
         ENDDO
         WRITE(F06,1901)
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,1006)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (MOI1(I,J),J=1,3)              ! MOI1 now are MOI's about ref pt in basic
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)

      ENDIF

      ! S, TRANS: Generate moments of inertia about c.g. in basic coord. system
      TRANS(1,1) =  MASS*(XB(2)*XB(2) + XB(3)*XB(3))
      TRANS(2,2) =  MASS*(XB(1)*XB(1) + XB(3)*XB(3))
      TRANS(3,3) =  MASS*(XB(1)*XB(1) + XB(2)*XB(2))
      TRANS(1,2) = -MASS*XB(1)*XB(2)
      TRANS(1,3) = -MASS*XB(1)*XB(3)
      TRANS(2,3) = -MASS*XB(2)*XB(3)
      TRANS(2,1) =  TRANS(1,2)
      TRANS(3,1) =  TRANS(1,3)
      TRANS(3,2) =  TRANS(2,3)
      DO I=1,3
         DO J=1,3
            MOI1(I,J) = MOI1(I,J) - TRANS(I,J)
         ENDDO
      ENDDO

      ! backup MOI1 to create IS
      IS(1,1) = MOI1(1,1)
      IS(2,2) = MOI1(2,2)
      IS(3,3) = MOI1(3,3)

      IS(1,2) = MOI1(1,2)
      IS(1,3) = MOI1(1,3)
      IS(2,3) = MOI1(2,3)

      IS(2,1) = MOI1(1,2)
      IS(3,1) = MOI1(1,3)
      IS(3,2) = MOI1(2,3)

      ! Set model total mass parameters about c.g.
      MODEL_MASS = WTMASS*MASS
      MODEL_IXX  = WTMASS*MOI1(1,1)
      MODEL_IYY  = WTMASS*MOI1(2,2)
      MODEL_IZZ  = WTMASS*MOI1(3,3)
      MODEL_XCG  = XREF(1) + XB(1)
      MODEL_YCG  = XREF(2) + XB(2)
      MODEL_ZCG  = XREF(3) + XB(3)
      DO I=1,6
         DO J=1,6
            MCG(I,J) = ZERO
         ENDDO
      ENDDO
      MCG(1,1) = MODEL_MASS
      MCG(2,2) = MODEL_MASS
      MCG(3,3) = MODEL_MASS
      MCG(4,4) = MODEL_IXX
      MCG(5,5) = MODEL_IYY
      MCG(6,6) = MODEL_IZZ

      ! Output MOI's about c.g.
      IF (REFPNT >= 0) THEN
         WRITE(F06,1007)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (MOI1(I,J),J=1,3)              ! MOI1 now are MOI's about cg in basic
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)
      ENDIF

      ! Q - Get principal MOI's and S transformation matrix (eigenvectors of MOI1)
      CALL GPWG_PMOI ( MOI1, Q, INFO )

      ! Write out princ MOI's and coord transf. Otherwise errors were written in subr GPWG_PMOI
      IF ((INFO == 0) .AND. (REFPNT >= 0)) THEN

         IF(REFPNT > -1) THEN
           ITABLE = -3
 1         FORMAT("WRITE OGPWG OP2: ",A)
 2         FORMAT("* DEBUG OGPWG ITABLE=",i4)
           WRITE(ERR,1) "START"
           WRITE(ERR,2) ITABLE

           TABLE_NAME = "OGPWG   "
           CALL WRITE_TABLE_HEADER(TABLE_NAME)
           ANALYSIS_CODE = 1   ! TODO: this is probably wrong, but is weird for this table

           WRITE(ERR,2) ITABLE
           CALL WRITE_ITABLE(ITABLE)
           WRITE(ERR,2) ITABLE

           WRITE(ERR,1) "WRITE_OPGWG_TABLE3"
           CALL WRITE_OPGWG_TABLE3(ITABLE, ANALYSIS_CODE, TABLE_CODE, REFPNT)
           WRITE(ERR,2) ITABLE
           CALL WRITE_ITABLE(ITABLE)
           ITABLE = ITABLE - 1

           WRITE(ERR,1) "WRITE_OPGWG_TABLE4"
           WRITE(ERR,2) ITABLE

           !data = (self.MO.ravel().tolist() + self.S.ravel().tolist() +
           !        mcg.ravel().tolist() + self.IS.ravel().tolist() + self.IQ.ravel().tolist() +
           !        self.Q.ravel().tolist())

           ! not verified
           ! mass shouldn't be a scalar (it's a vector),
           ! but for physical structure they're all the same
           WRITE(OP2) 78   ! the number of values we're going to write
           WRITE(OP2) ((REAL(RB_MASS_BASIC(I,J), 4), J=1,6), I=1,6),              & ! (6,6) MO - 36
                      ((REAL(TRANS(I,J), 4), J=1,3), I=1,3),                      & ! (3,3) S - 9
                      ! mass should be (3,1) instead of (1,)
                      ! CG should be (3,3) instead of (3,1)
                      ! faking...
                      REAL(MASS, 4), (REAL(XB(I), 4),I=1,3),                      & ! mass - cg
                      REAL(MASS, 4), (REAL(XB(I), 4),I=1,3),                      & ! mass - cg
                      REAL(MASS, 4), (REAL(XB(I), 4),I=1,3),                      & ! mass - cg - 12

                      ((REAL(IS(I,J), 4), J=1,3), I=1,3),                         & ! (3,3) IS
                      REAL(MODEL_IXX, 4), REAL(MODEL_IYY, 4), REAL(MODEL_IZZ, 4), & ! (3,)  IQ
                      ((REAL(Q(I,J), 4), J=1,3), I=1,3)                             ! (3,3) Q - 21
           !CALL WRITE_ITABLE(ITABLE)
           !ITABLE = ITABLE - 1
           CALL END_OP2_TABLE(ITABLE)
           WRITE(ERR,2) ITABLE
           WRITE(ERR,1) "WRITE_OPGWG_TABLE END"
         ENDIF


         WRITE(F06,1008)
         WRITE(F06,1900)

         ! I(Q) - MOI1 is now principal MOI matrix
         DO I=1,3
            WRITE(F06,1101) (MOI1(I,J),J=1,3)
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,1009)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (Q(I,J),J=1,3)
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
  902 FORMAT(/,'  O U T P U T   F R O M   T H E   G R I D   P O I N T   W E I G H T   G E N E R A T O R   F O R   R E S I D U A L' &
             ,'   S T R U C T U R E')

  903 FORMAT(/,7X,'O U T P U T   F R O M   T H E   G R I D   P O I N T   W E I G H T   G E N E R A T O R   F O R   O V E R A L L', &
               '   M O D E L')

 1001 FORMAT(26X,'                   (reference point is basic coord system origin)'                                               &
            ,/)

 1002 FORMAT(26X,'                      (reference point is grid point '    ,I8,')'                                                &
            ,/)

 1004 FORMAT(  26X,'                            Total mass = '    ,1ES13.6                                                         &
            ,/)

 1005 FORMAT(61X,                                    'X             Y             Z'                                               &
          ,/,26X,'             C.G. location :',3(1ES14.6),                                                                        &
           /,26X,'              (relative to reference point in basic coordinate system)'                                          &
          ,//)

 1006 FORMAT(36X,'M.O.I. matrix - about reference point in basic coordinate system')

 1007 FORMAT(34X,'M.O.I. matrix - about above c.g. location in basic coordinate system')

 1008 FORMAT(37X,'M.O.I. matrix - about above c.g. location in principal directions')

 1009 FORMAT(37X,'Transformation from basic coordinates to principal directions')

 1013 FORMAT(  20X,' N O T E: Since user requested modal effective masses to be calculated the reference point',/,                 &
               20X,'          for this mass calc has been changed to the PARAM MEFMASS grid, ',I8,/)

 1021 FORMAT(29X,'6x6 Rigid body mass matrix - about reference point in basic coordinate system')

 1022 FORMAT(22X,'*',3(1ES14.6),'  *',3(1ES14.6),'  *')

 1901 FORMAT('                      ***                                                                                     ***')

 1902 FORMAT('                      *  ************  ************  ************  *  ************  ************  ************  *')

 1101 FORMAT(45X,'*',3(1ES14.6),'  *')

 1900 FORMAT(45X,'***',40X,'***')

 1402 FORMAT(' *WARNING    : PARAM GRDPNT (OR PARAM MEFMGRID) REFERENCES NONEXISTENT GRID POINT ',I8,'. BASIC ORIGIN WILL BE USED')

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')


! **********************************************************************************************************************************

      END SUBROUTINE GPWG

!-------------------------------------------------------------------------------------------------------------
      SUBROUTINE WRITE_OPGWG_TABLE3(ITABLE, ANALYSIS_CODE, TABLE_CODE, REFERENCE_POINT)

     USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY               :  ERR, F06, OP2

      USE ROD_BAR_OUTPUT, ONLY       :  WRITE_ROD

      IMPLICIT NONE

      !INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      CHARACTER(LEN=128)              :: TITLE             ! the model TITLE
      CHARACTER(LEN=128)              :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABEL             ! the subcase LABEL
      INTEGER(LONG), INTENT(IN) :: ANALYSIS_CODE, TABLE_CODE, REFERENCE_POINT

      INTEGER(LONG) :: DEVICE_CODE, APPROACH_CODE
      INTEGER(LONG), PARAMETER    :: ISUBCASE = 0     ! the current subcase
      INTEGER(LONG), PARAMETER    :: NUM_WIDE = 79     ! the number of "words" for an element

      ! dummy
      TITLE = ""
      SUBTITLE = ""
      LABEL = ""
      DEVICE_CODE = 1  ! plot

      APPROACH_CODE = ANALYSIS_CODE*10 + DEVICE_CODE

      ! 584 bytes
      WRITE(OP2) 146
      WRITE(OP2) APPROACH_CODE, TABLE_CODE, REFERENCE_POINT, ISUBCASE, 0,   &
            0, 0, 0, 0, NUM_WIDE, &
            0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, &
            TITLE, SUBTITLE, LABEL

      ITABLE = ITABLE - 1        ! flip it to -4, -6, ... so we don't have to do this later
      END SUBROUTINE WRITE_OPGWG_TABLE3


      SUBROUTINE GPWG_PMOI (MOI1, Q, INFO )

! Jacobi solution for 3x3 eigenvalue problem used in finding principal moments of inertia for the Grid Point Weight Generator (GPWG)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  SUPWARN, WTMASS

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GPWG_PMOI'
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER( 1*BYTE), PARAMETER   :: JOBZ      = 'V'   ! Indicates to solve for eigenvalues and vectors in LAPACK subr DSYEV
      CHARACTER( 1*BYTE)              :: TRANSA            ! Transpose indicator for subr DGEMM
      CHARACTER( 1*BYTE)              :: TRANSB            ! Transpose indicator for subr DGEMM
      CHARACTER( 1*BYTE), PARAMETER   :: UPLO      = 'U'   ! Indicates array KBAND is the upper triangular part of KAA

      INTEGER(LONG), INTENT(OUT)      :: INFO              ! = 0:  successful exit
!                                                            < 0:  if INFO = -i, the i-th argument had an illegal value
!                                                            > 0:  if INFO = i, the algorithm failed to converge; i off_diag
!                                                            elems of an intermediate tridiag form did not converge to zero
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: N         = 3     ! Order of matrix MOI1
      INTEGER(LONG), PARAMETER        :: LWORK     = 3*N-1 ! Size of array WORK


      REAL(DOUBLE) , INTENT(INOUT)    :: MOI1(3,3)         ! On entry, the MOI's about c.g. in basic coords
!                                                            On exit , the principal MOI's in basic coords (if INFO = 0)
      REAL(DOUBLE) , INTENT(OUT)      :: Q(3,3)            ! Transformation from basic to principal directions
!                                                            Q = Z'. That is;
!                                                            U(principal) = Q*U(basic), where U is a displ vector
      REAL(DOUBLE)                    :: Z(3,3)            ! When subr DSYEV is called, Z = input MOI1
!                                                            When subr DSYEV returns, Z = eigenvecs: PMOI = Z'*MOI1*Z
      REAL(DOUBLE)                    :: DUM(N,N)          ! Intermediate result in calculating Z'*MOI1*Z
      REAL(DOUBLE) , PARAMETER        :: ALPHA     = ONE   ! Scalar multiplier for subr DGEMM
      REAL(DOUBLE) , PARAMETER        :: BETA      = ZERO  ! Scalar multiplier for subr DGEMM
      REAL(DOUBLE)                    :: PMOI(3)           ! Principal MOI's in a 1-D array
      REAL(DOUBLE)                    :: WORK(LWORK)       ! Workspace

      EXTERNAL                        :: DGEMM



! **********************************************************************************************************************************
      ! Initialize outputs
      INFO = 0

      DO I=1,N
         DO J=1,N
            Q(I,J) = ZERO
         ENDDO
      ENDDO

      ! Use LAPACK driver DSYEV to get all eigenvalues (the principal MOI's) and eigenvectors (returned in Z)
      DO I=1,N
         DO J=1,N
            Z(I,J) = MOI1(I,J)
            Q(I,J) = ZERO
         ENDDO
      ENDDO

      CALL DSYEV ( JOBZ, UPLO, N, Z, N, PMOI, WORK, LWORK, INFO )

      ! Set MOI1 matrix to zero for off-diag terms and to PMOI for diag terms.
      ! Set Q = Z' (' = transpose)

      CALLED_SUBR = 'DSTEQR'
      IF      (INFO < 0) THEN
         ! LAPACK subr XERBLA should have reported error on an illegal argument
         ! in a called LAPACK subr, so we should not have gotten here
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO > 0) THEN
        ! No convergence in subr DSYEV
        WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1001) CALLED_SUBR, SUBR_NAME
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1001) CALLED_SUBR, SUBR_NAME
         ENDIF
         RETURN

      ELSE
         ! INFO=0, so no error

         ! Set Q = Z'
         DO I=1,N
            DO J=1,N
               Q(J,I) = Z(I,J)
            ENDDO
         ENDDO

         ! MOI1 <-- Z'*MOI1*Z, diags should be PMOI's, off-diags should be zero
         TRANSA = 'N'
         TRANSB = 'N'
         CALL DGEMM ( TRANSA, TRANSB, N, N, N, ALPHA, MOI1, N, Z,   N, BETA, DUM,  N )
         TRANSA = 'T'
         TRANSB = 'N'
         CALL DGEMM ( TRANSA, TRANSB, N, N, N, ALPHA, Z,    N, DUM, N, BETA, MOI1, N )

      ENDIF




      RETURN

! **********************************************************************************************************************************
  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')

 1001 FORMAT(' *WARNING    : LAPACK SUBR ',A8,' CALLED BY SUBROUTINE ',A                                                           &
                    ,/,14X,' CANNOT CONVERGE IN ATTEMPTING TO FIND PRINCIPAL MOIs'                                                 &
                    ,/,14X,' THE ALGORITHM HAS FAILED TO FIND ALL THE EIGENVALUES (PRINCIPAL MOIs) IN 90 ITERATIONS'               &
                    ,/,14X,' PRINCIPAL MOIs AND PRINCIPAL DIRECTIONS CANNOT BE FOUND')

! **********************************************************************************************************************************

      END SUBROUTINE GPWG_PMOI


      SUBROUTINE GPWG_USERIN ( IEID )

! Generates rigid body mass properties for one USERIN element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NGRID, SOL_NAME, WARN_ERR
      USE PARAMS, ONLY                :  EPSIL, GRDPNT, MEFMGRID, SUPWARN, WTMASS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS, EID, GRID_ID, ME, PLY_NUM, RGRID, USERIN_RBM0


      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE EMG_MOD, ONLY               :  EMG
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GPWG_USERIN'
      CHARACTER(1*BYTE)               :: OPT(6)            ! Option flags for what to calculate when subr EMG is called

      INTEGER(LONG), INTENT(IN)       :: IEID              ! Internal element ID for the USERIN element to process
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: INFO        = 0   ! An output from subr GPWG_PMOI, called herein
      INTEGER(LONG)                   :: GRDPNT_DEF        ! Default value of GRDPNT


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: M0                ! An intermediate variable used in calc model mass props
      REAL(DOUBLE)                    :: MOI1(3,3)         ! Moments of inertia (several diff interps during exec of this subr)
      REAL(DOUBLE)                    :: MX                ! First moment of MASS about X axis
      REAL(DOUBLE)                    :: MY                ! First moment of MASS about Y axis
      REAL(DOUBLE)                    :: MZ                ! First moment of MASS about Z axis
      REAL(DOUBLE)                    :: Q(3,3)            ! Output from subr GPWG_PMOI, called herein (transform to princ dir's)
      REAL(DOUBLE)                    :: TRANS(3,3)        ! Transfer terms when MOI's about c.g. ard calc'd from MOI's about XREF
      REAL(DOUBLE)                    :: XB(3)             ! Basic coord diffs bet c.g. and XREF in X, Y, Z directions
      REAL(DOUBLE)                    :: XD(3)             ! Basic coord diffs bet a mass (at it's c.m.) and XREF in X, Y, Z dirs
      REAL(DOUBLE)                    :: XREF(3)           ! GRDPNT basic coords (or origin of basic sys if GRDPNT doesn't exist)

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Set defaults in case GRDPNT grid cannot be found or is input as 0 (basic origin)

      GRDPNT_DEF = 0
      XREF(1)    = ZERO
      XREF(2)    = ZERO
      XREF(3)    = ZERO

! Get reference point coordinates in basic system for the reference point

      IF (GRDPNT /= -1) THEN
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRDPNT, GRID_ID_ROW_NUM )
         IF (GRID_ID_ROW_NUM /= -1) THEN                   ! GRDPNT is a grid point in the model, so get its basic coords
            XREF(1) = RGRID(GRID_ID_ROW_NUM,1)
            XREF(2) = RGRID(GRID_ID_ROW_NUM,2)
            XREF(3) = RGRID(GRID_ID_ROW_NUM,3)
         ELSE                                              ! GRDPNT was not 0 and is not a grid number in the model
            IF (GRDPNT /= 0) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1402) GRDPNT
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1402) GRDPNT
               ENDIF
            ENDIF
            GRDPNT = GRDPNT_DEF
         ENDIF
      ENDIF

! Generate total mass, first and second moments by summing up mass terms. XD(i) are components of vector from
! ref point to a mass point. At this time, mass units are input units without PARAM WTMASS which is what we want for
! the grid point weight generator. Later the mass will be converted by multiplying by WTMASS.

      M0  = ZERO
      MX    = ZERO
      MY    = ZERO
      MZ    = ZERO

      XB(1) = ZERO
      XB(2) = ZERO
      XB(3) = ZERO

! Process this USERIN element mass terms

      OPT(1) = 'Y'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      IERROR  = 0
      PLY_NUM = 0
      CALL EMG ( IEID, OPT, 'N', SUBR_NAME, 'N' )          ! 'N' means do not write to BUG file

      IF (NUM_EMG_FATAL_ERRS == 0) THEN

         M0   = USERIN_RBM0(1,1)
         IF (DABS(M0) > EPS1) THEN
            XD(1) = USERIN_RBM0(2,6)/M0 !- XREF(1)
            XD(2) = USERIN_RBM0(3,4)/M0 !- XREF(2)
            XD(3) = USERIN_RBM0(1,5)/M0 !- XREF(3)
            MX    = MX + M0*XD(1)
            MY    = MY + M0*XD(2)
            MZ    = MZ + M0*XD(3)
         ELSE
            XD(1) = ZERO
            XD(2) = ZERO
            XD(3) = ZERO
         ENDIF


      ENDIF

      OPT(1) = 'N'                                         ! Reset subr EMG option value:

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! Errors from subr EMG, so quit
      ENDIF

! XB(I) are components of distance from reference point, XREF, to c.g.

      IF (DABS(M0) > EPS1) THEN
         XB(1) = MX/M0
         XB(2) = MY/M0
         XB(3) = MZ/M0
      ENDIF

! Output results so far

      IF (GRDPNT >= 0) THEN
         WRITE(F06,1000) EID
         IF (GRDPNT == 0) THEN
            WRITE(F06,1001)
         ELSE IF (GRDPNT > 0) THEN
            WRITE(F06,1002) GRDPNT
         ENDIF
      ENDIF

      IF (GRDPNT >= 0) THEN

         WRITE(F06,1004) M0

         WRITE(F06,1005) (XB(I),I=1,3)

         WRITE(F06,1021)
         WRITE(F06,1901)
         DO I=1,3
            WRITE(F06,1022) (USERIN_RBM0(I,J),J=1,6)
         ENDDO
         WRITE(F06,1902)
         DO I=4,6
            WRITE(F06,1022) (USERIN_RBM0(I,J),J=1,6)
         ENDDO
         WRITE(F06,1901)
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,1006)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (USERIN_RBM0(I+3,J+3),J=1,3)
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)

      ENDIF

! Generate moments of inertia about c.g. in basic coord. system

      TRANS(1,1) =  M0*(XB(2)*XB(2) + XB(3)*XB(3))
      TRANS(2,2) =  M0*(XB(1)*XB(1) + XB(3)*XB(3))
      TRANS(3,3) =  M0*(XB(1)*XB(1) + XB(2)*XB(2))
      TRANS(1,2) = -M0*XB(1)*XB(2)
      TRANS(1,3) = -M0*XB(1)*XB(3)
      TRANS(2,3) = -M0*XB(2)*XB(3)
      TRANS(2,1) =  TRANS(1,2)
      TRANS(3,1) =  TRANS(1,3)
      TRANS(3,2) =  TRANS(2,3)
      DO I=1,3
         DO J=1,3
            MOI1(I,J) = USERIN_RBM0(I+3,J+3) - TRANS(I,J)
         ENDDO
      ENDDO

! Output MOI's about c.g.

      IF (GRDPNT >= 0) THEN
         WRITE(F06,1007)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (MOI1(I,J),J=1,3)              ! MOI1 now are MOI's about cg in basic
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)
      ENDIF

! Get principal MOI's and transformation matrix (eigenvectors of MOI1)

      CALL GPWG_PMOI ( MOI1, Q, INFO )

! Write out princ MOI's and coord transf. Otherwise errors were written in subr GPWG_PMOI

      IF ((INFO == 0) .AND. (GRDPNT >= 0)) THEN
         WRITE(F06,1008)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (MOI1(I,J),J=1,3)              ! MOI1 is now principal MOI matrix
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,1009)
         WRITE(F06,1900)
         DO I=1,3
            WRITE(F06,1101) (Q(I,J),J=1,3)
         ENDDO
         WRITE(F06,1900)
         WRITE(F06,*)
         WRITE(F06,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1000 FORMAT(/,'     O U T P U T   F R O M   T H E   G R I D   P O I N T   W E I G H T   G E N E R A T O R   F O R   U S E R I N', &
               '   E L E M',I8)

 1001 FORMAT(26X,'                   (reference point is basic coord system origin)'                                               &
            ,//)

 1002 FORMAT(26X,'                      (reference point is grid point '    ,I8,')'                                                &
            ,//)

 1004 FORMAT(  26X,'                            Total mass = '    ,1ES13.6                                                         &
            ,//)

 1005 FORMAT(61X,                                    'X             Y             Z'                                               &
          ,/,26X,'             C.G. location :',3(1ES14.6),                                                                        &
           /,26X,'              (relative to reference point in basic coordinate system)'                                          &
          ,//)

 1006 FORMAT(36X,'M.O.I. matrix - about reference point in basic coordinate system')

 1007 FORMAT(34X,'M.O.I. matrix - about above c.g. location in basic coordinate system')

 1008 FORMAT(37X,'M.O.I. matrix - about above c.g. location in principal directions')

 1009 FORMAT(37X,'Transformation from basic coordinates to principal directions')

 1021 FORMAT(29X,'6x6 Rigid body mass matrix - about reference point in basic coordinate system')

 1022 FORMAT(22X,'*',3(1ES14.6),'  *',3(1ES14.6),'  *')

 1901 FORMAT('                      ***                                                                                     ***')

 1902 FORMAT('                      *  ************  ************  ************  *  ************  ************  ************  *')

 1101 FORMAT(45X,'*',3(1ES14.6),'  *')

 1900 FORMAT(45X,'***',40X,'***')

 1402 FORMAT(' *WARNING    : PARAM GRDPNT (OR PARAM MEFMGRID) REFERENCES NONEXISTENT GRID POINT ',I8,'. BASIC ORIGIN WILL BE USED')

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE GPWG_USERIN


      SUBROUTINE RB_DISP_MATRIX_PROC ( REF_PT_TXT, REF_PT )

! Generates a 6 x 6 rigid body displacement matrix in global coords for one grid point via the following procedure:
!    1) for each grid generate a 6 x 6 rigid body displ matrix in basic coords relative to the reference grid which can be:
!         a) a grid defined by param (REF_PT_TXT = EQCHK_REF_GRID)
!         b) basic systen origin (REF_PT_TXT = 'BASIC_ORIGIN')
!         c) model CG (REF_PT_TXT = 'CG')
!         d) an arbitrary grid (REF_PT_TXT = 'GRID' and REF_PT = the grid ID)
!    2) transform to global coords at the grid
!    3) assemble the 6 x 6 rigid body displ matrices, in global coords, for all grids

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD, NGRID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI, TDOF_ROW_START
      USE PARAMS, ONLY                :  EQCHK_REF_GRID, SUPWARN
      USE MODEL_STUF, ONLY            :  CORD, GRID, RGRID, GRID_ID, INV_GRID_SEQ, MODEL_XCG, MODEL_YCG, MODEL_ZCG
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RB_DISP_MATRIX_PROC'
      CHARACTER(LEN=*), INTENT(IN)    :: REF_PT_TXT        ! Reference point used in calculating the 6 rigid body displ vectors
      CHARACTER( 2*BYTE)              :: COMP(6)           ! Text reference to the 6 components of displ (T1, T2, etc)

      INTEGER(LONG), INTENT(IN)       :: REF_PT            ! An actual grid ID (only used if REF_PT_TXT = 'GRID')
      INTEGER(LONG)                   :: ECORD_K     = 0   ! Global coord ID (actual) for grid AGRID_I
      INTEGER(LONG)                   :: ECORD_R     = 0   ! Global coord ID (actual) for grid AGRID_R
      INTEGER(LONG)                   :: AGRID             ! An actual grid number
      INTEGER(LONG)                   :: AGRID_K           ! Grid ID (actual) for a grid whose RB matrix is being generated
      INTEGER(LONG)                   :: AGRID_R           ! Grid ID (actual) for the reference grid (EQCHK_REF_GRID)

      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_R ! Row number in array GRID_ID where AGRID_R is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_K ! Row number in array GRID_ID where AGRID_K is found
      INTEGER(LONG)                   :: G_SET_COL         ! Col no. in array TDOF where the  G-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: G_SET_DOF_NUM     ! DOF number from TDOF for a G-set DOF
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices or counters
      INTEGER(LONG)                   :: ICORD_K           ! Internal coord ID corresponding to ECORD_I
      INTEGER(LONG)                   :: ICORD_R           ! Internal coord ID corresponding to ECORD_R
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid


      REAL(DOUBLE)                    :: DUM1(6,6)         ! Intermediate result in obtaining RB_GRID_GLOBL
      REAL(DOUBLE)                    :: DX0               ! X coord difference between grid I and ref grid
      REAL(DOUBLE)                    :: DY0               ! Y coord difference between grid I and ref grid
      REAL(DOUBLE)                    :: DZ0               ! Z coord difference between grid I and ref grid
      REAL(DOUBLE)                    :: PHID, THETAD      ! Angles output from subr GEN_T0L, called herein but not needed here
      REAL(DOUBLE)                    :: RB_GRID_BASIC(6,6)! Rigid body displ matrix for grid AGRID_I in basic  coords
      REAL(DOUBLE)                    :: RB_GRID_GLOBL(6,6)! Rigid body displ matrix for grid AGRID_I in global coords
      REAL(DOUBLE)                    :: T0G_K(3,3)        ! Transforms a vector in basic coords to global coords for grid AGRID_I
      REAL(DOUBLE)                    :: T0G_R(3,3)        ! Transforms a vector in basic coords to global coords for grid AGRID_R
      REAL(DOUBLE)                    :: TTR_K(6,6)        ! 6x6 matrix with 2 TOG_K(3,3) matrices (one for trans, one for rot)
      REAL(DOUBLE)                    :: TTR_R(6,6)        ! 6x6 matrix with 2 TOG_R(3,3) matrices (one for trans, one for rot)
      REAL(DOUBLE)                    :: X0_R              ! Basic X coord of AGRID_R (ref grid)
      REAL(DOUBLE)                    :: Y0_R              ! Basic Y coord of AGRID_R (ref grid)
      REAL(DOUBLE)                    :: Z0_R              ! Basic Z coord of AGRID_R (ref grid)
      REAL(DOUBLE)                    :: X0_K              ! Basic X coord of AGRID_I
      REAL(DOUBLE)                    :: Y0_K              ! Basic Y coord of AGRID_I
      REAL(DOUBLE)                    :: Z0_K              ! Basic Z coord of AGRID_I



! **********************************************************************************************************************************
      CALL TDOF_COL_NUM ( 'G ',  G_SET_COL )

! Initialize

      DO I=1,6
         DO J=1,6
            TTR_R(I,J) = ZERO
         ENDDO
      ENDDO

! Get basic coords of the ref grid

      IF ((REF_PT_TXT == 'EQCHK REF GRID') .OR. (REF_PT_TXT == 'GRID')) THEN

         IF (REF_PT_TXT == 'EQCHK REF GRID') THEN
            AGRID_R = EQCHK_REF_GRID
         ELSE
            AGRID_R = REF_PT
         ENDIF

         IF (AGRID_R == 0) THEN
            ECORD_R = 0
            X0_R = ZERO
            Y0_R = ZERO
            Z0_R = ZERO
            DO I=1,6
               TTR_R(I,I) = ONE
            ENDDO
         ELSE
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_R, GRID_ID_ROW_NUM_R )
            IF (GRID_ID_ROW_NUM_R == -1) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1408) SUBR_NAME, AGRID_R, 'GRID_ID'
               WRITE(F06,1408) SUBR_NAME, AGRID_R, 'GRID_ID'
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_R, NUM_COMPS, SUBR_NAME )
            IF (NUM_COMPS == 6) THEN                       ! AGRID_R is a physical grid
               X0_R = RGRID(GRID_ID_ROW_NUM_R,1)
               Y0_R = RGRID(GRID_ID_ROW_NUM_R,2)
               Z0_R = RGRID(GRID_ID_ROW_NUM_R,3)
               ECORD_R = GRID(GRID_ID_ROW_NUM_R,3)
               IF (ECORD_R /= 0) THEN
                  DO I=1,NCORD                             ! We know that if ECORD_R > 0, ICORD_R is defined (chk in CORD_PROC)
                     IF (ECORD_R == CORD(I,2)) THEN
                        ICORD_R = I
                        EXIT
                     ENDIF
                  ENDDO
                  CALL GEN_T0L ( GRID_ID_ROW_NUM_R, ICORD_R, THETAD, PHID, T0G_R )
                  DO I=1,3
                     DO J=1,3
                        TTR_R(I  ,J  ) = T0G_R(I,J)
                        TTR_R(I+3,J+3) = T0G_R(I,J)
                     ENDDO
                  ENDDO
               ELSE
                  DO I=1,6
                     TTR_R(I,I) = ONE
                  ENDDO
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1001) AGRID_R
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1001) AGRID_R
               ENDIF
               ECORD_R = 0
               X0_R = ZERO
               Y0_R = ZERO
               Z0_R = ZERO
               DO I=1,6
                  TTR_R(I,I) = ONE
               ENDDO
            ENDIF
         ENDIF

      ELSE IF (REF_PT_TXT == 'BASIC ORIGIN')  THEN

         X0_R = ZERO
         Y0_R = ZERO
         Z0_R = ZERO

         DO I=1,6
            TTR_R(I,I) = ONE
         ENDDO

      ELSE IF (REF_PT_TXT == 'CG') THEN

         X0_R = MODEL_XCG
         Y0_R = MODEL_YCG
         Z0_R = MODEL_ZCG

         DO I=1,6
            TTR_R(I,I) = ONE
         ENDDO

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1409) SUBR_NAME, REF_PT_TXT
         WRITE(F06,1409) SUBR_NAME, REF_PT_TXT
         CALL OUTA_HERE ( 'Y' )

      ENDIF

! Get the 6x6 RB matrices for each grid and transform to global

      DO K=1,NGRID

         AGRID_K = GRID_ID(INV_GRID_SEQ(K))
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(K), NUM_COMPS, SUBR_NAME )

         IF (NUM_COMPS == 6) THEN                          ! Only process physical grids. Let rows of RBGLOBAL = 0 otherwise

            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_K, GRID_ID_ROW_NUM_K )

            DO I=1,6
               DO  J=1,6
                  TTR_K(I,J) = ZERO
               ENDDO
            ENDDO

            X0_K = RGRID(GRID_ID_ROW_NUM_K,1)
            Y0_K = RGRID(GRID_ID_ROW_NUM_K,2)
            Z0_K = RGRID(GRID_ID_ROW_NUM_K,3)

            DX0  = X0_K - X0_R
            DY0  = Y0_K - Y0_R
            DZ0  = Z0_K - Z0_R

            DO I=1,6
               DO J=1,6
                  RB_GRID_BASIC(I,J) = ZERO
               ENDDO
               RB_GRID_BASIC(I,I) = ONE
            ENDDO
            RB_GRID_BASIC(1,5) =  DZ0
            RB_GRID_BASIC(1,6) = -DY0
            RB_GRID_BASIC(2,4) = -DZ0
            RB_GRID_BASIC(2,6) =  DX0
            RB_GRID_BASIC(3,4) =  DY0
            RB_GRID_BASIC(3,5) = -DX0


            ECORD_K = GRID(GRID_ID_ROW_NUM_K,3)
            IF (ECORD_K /= 0) THEN
               DO I=1,NCORD
                  IF (ECORD_K == CORD(I,2)) THEN
                     ICORD_K = I
                     EXIT
                  ENDIF
               ENDDO
               CALL GEN_T0L ( GRID_ID_ROW_NUM_K, ICORD_K, THETAD, PHID, T0G_K )
               DO I=1,3
                  DO J=1,3
                     TTR_K(I  ,J  ) = T0G_K(I,J)
                     TTR_K(I+3,J+3) = T0G_K(I,J)
                  ENDDO
               ENDDO
            ELSE
               DO I=1,6
                  TTR_K(I,I) = ONE
               ENDDO
            ENDIF

            IF (ECORD_R == 0) THEN
               DO I=1,6
                  DO J=1,6
                     DUM1(I,J) = RB_GRID_BASIC(I,J)
                  ENDDO
               ENDDO
            ELSE
               CALL MATMULT_FFF   ( RB_GRID_BASIC, TTR_R, 6, 6, 6, DUM1 )
            ENDIF
            IF (ECORD_K == 0) THEN
               DO I=1,6
                  DO J=1,6
                     RB_GRID_GLOBL(I,J) = DUM1(I,J)
                  ENDDO
               ENDDO
            ELSE
               CALL MATMULT_FFF_T ( TTR_K, DUM1, 6, 6, 6, RB_GRID_GLOBL  )
            ENDIF

            IF ((DEBUG(11) == 1) .OR. (DEBUG(11) == 3)) THEN
               WRITE(F06,111) AGRID_K
               DO I=1,6
                  WRITE(F06,112) (RB_GRID_GLOBL(I,J),J=1,6)
               ENDDO
               WRITE(F06,*)
            ENDIF

            AGRID = GRID_ID(INV_GRID_SEQ(K))
!xx         CALL CALC_TDOF_ROW_NUM ( AGRID, ROW_NUM_START, 'N' )
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            G_SET_DOF_NUM = TDOF(ROW_NUM_START,G_SET_COL)
            DO I=1,6
               DO J=1,6
                  RBGLOBAL_GSET(G_SET_DOF_NUM+I-1,J) = RB_GRID_GLOBL(I,J)
               ENDDO
            ENDDO

         ENDIF

      ENDDO

      IF ((DEBUG(11) == 2) .OR. (DEBUG(11) == 3)) THEN
         WRITE(F06,121)
         IF      (REF_PT_TXT == 'CG') THEN
            WRITE(F06,122)
         ELSE IF (REF_PT_TXT == 'BASIC ORIGIN'  ) THEN
            WRITE(F06,123)
         ELSE IF ((REF_PT_TXT == 'EQCHK REF GRID') .OR. (REF_PT_TXT == 'GRID')) THEN
            IF (AGRID_R == 0) THEN
               WRITE(F06,123)
            ELSE
               WRITE(F06,124) AGRID_R, ECORD_R
            ENDIF
         ENDIF
         WRITE(F06,125)
         L = 0
         COMP(1) = 'T1'
         COMP(2) = 'T2'
         COMP(3) = 'T3'
         COMP(4) = 'R1'
         COMP(5) = 'R2'
         COMP(6) = 'R3'
         DO I=1,NGRID
            AGRID_K = GRID_ID(INV_GRID_SEQ(I))
            AGRID = GRID_ID(INV_GRID_SEQ(I))
            CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               L = L + 1
               IF (J == 1) THEN
                  WRITE(F06,141) AGRID_K, COMP(J), (RBGLOBAL_GSET(L,K),K=1,6)
               ELSE
                  WRITE(F06,142)          COMP(J), (RBGLOBAL_GSET(L,K),K=1,6)
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
  111 FORMAT(' RIGID BODY DISPL MATRIX IN GLOBAL COORDS FOR GRID ',I8)

  112 FORMAT(1X,6(1ES15.6))

  121 FORMAT('                                RBGLOBAL RIGID BODY MATRIX - G-SET RIGID BODY DISPLS IN GLOBAL COORDS')

  122 FORMAT('                                COLUMNS ARE DUE TO UNIT DISPLACEMENTS OF THE MODEL CG IN BASIC COORDS',/)

  123 FORMAT('                           COLUMNS ARE DUE TO UNIT DISPLACEMENTS AT THE BASIC SYSTEM ORIGIN IN BASIC COORDS',/)

  124 FORMAT('                        COLUMNS ARE DUE TO UNIT DISPLCEMENTS AT GRID ',I8,' IN GLOBAL COORD SYSTEM ',I8,/)

  125 FORMAT(13X,'GRID COMP       T1             T2             T3             R1             R2             R3')

  141 FORMAT(9X,I8,2X,A2,6(1ES15.6))

  142 FORMAT(19X,A2,6(1ES15.6))

 1001 FORMAT(' *WARNING    : REFERENCE GRID = ',I8,' FOR CALCULATING RIGID BODY RBGLOBAL MATRIX MUST BE A PHYSICAL GRID NOT A',    &
                           ' SCALAR POINT. BASIC ORIGIN WILL BE USED INSTEAD')

 1408 FORMAT(' *ERROR  1408: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' GRID NUMBER ',I8,' DOES NOT EXIST IN ARRAY ',A)

 1409 FORMAT(' *ERROR  1409: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT INDICATOR OF REF PT = ',A,' FOR CALC OF RIGID BODY DISPL VECTORS')

29381 format(' In RB_DISP_MATRIX_PROC: Act grid = ',i8,' is at row ',i8,' in RGRID')


! **********************************************************************************************************************************

      END SUBROUTINE RB_DISP_MATRIX_PROC


   END MODULE GRID_POINT_WEIGHT


   ! Compatibility entry point for USER_DEFINED_ELEMENTS. The module procedure remains
   ! available to existing GRID_POINT_WEIGHT callers while this external symbol avoids
   ! a module dependency cycle through EMG_MOD.
      SUBROUTINE RB_DISP_MATRIX_PROC ( REF_PT_TXT, REF_PT )

      USE GRID_POINT_WEIGHT, ONLY : RB_DISP_MATRIX_PROC_MOD => RB_DISP_MATRIX_PROC

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: REF_PT_TXT
      INTEGER, INTENT(IN) :: REF_PT

      CALL RB_DISP_MATRIX_PROC_MOD ( REF_PT_TXT, REF_PT )

      END SUBROUTINE RB_DISP_MATRIX_PROC
