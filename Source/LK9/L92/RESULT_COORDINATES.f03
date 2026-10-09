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

   MODULE RESULT_COORDINATES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELMDIS, ELMDIS_PLY, TRANSFORM_NODE_FORCES, TRANSFORM_SHELL_STR, STR_TENSOR_TRANSFORM

   CONTAINS

      SUBROUTINE ELMDIS

! Get displs for one element, one subcase from list of all displ's (in UG_COL). Transform them to local elem coords.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, INT_SC_NUM, meldof, MELGP, NCORD, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  AGRID, CAN_ELEM_TYPE_OFFSET, GRID, CORD, BGRID, ELGP, ELDOF, GRID_ID, OFFSET, OFFDIS,     &
                                         SCNUM, TE, TYPE, UEB, UEG, UEL, UGG
      USE COL_VECS, ONLY              :  UG_COL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EID, AGRID

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE MATERIAL_MATRIX_TRANSFER, ONLY:  MATGET, MATPUT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMDIS'

      INTEGER(LONG)                   :: G_SET_COL         ! Col number in TDOF where the G-set DOF's exist
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number
      INTEGER(LONG)                   :: GLOBAL_CID        ! Global coord. sys. ID for a grid (BGRID(i)) of the element
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I1,I2             ! Calculated indices in arrays
      INTEGER(LONG)                   :: ICORD             ! Internal coord. system corresponding to GLOBAL_CID
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG), PARAMETER        :: NCOLA     = 3     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NCOLB     = 1     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NROWA     = 3     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NROW      = 3     ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG), PARAMETER        :: NCOL      = 1     ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG)                   :: NUM_COMPS         ! Either 6 or 1 depending on whether grid is a physical grid or a SPOINT
      INTEGER(LONG)                   :: PROW              ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG), PARAMETER        :: PCOL      = 1     ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW          ! Row no. in array TDOF to find GDOF DOF number


      REAL(DOUBLE)                    :: DXI               ! An offset distance in direction 1
      REAL(DOUBLE)                    :: DYI               ! An offset distance in direction 2
      REAL(DOUBLE)                    :: DZI               ! An offset distance in direction 3
      REAL(DOUBLE)                    :: T0G(3,3)          ! Coord transformation matrix - basic to global
      REAL(DOUBLE)                    :: DUM1(3),DUM2(3)   ! Dummy arrays needed in transforming from global to basic coords
      REAL(DOUBLE)                    :: THETAD,PHID       ! Returns from subr GEN_T0L (not used here)


! **********************************************************************************************************************************
! Initialize

      DO I=1,MELDOF
         UEL(I) = ZERO
         UEB(I) = ZERO
         UEG(I) = ZERO
         UGG(I) = ZERO
      ENDDO

! Get displs for the grid points for this elem. AGRID(i) contains the actual grid point no's and BGRID(i) contains
! the internal grid point no's for this elem. These are generated in subroutine ELMDAT called by EMG
! called by OFP3 prior to this call. The grid point displ's in global coord's are in array UG_COL created in LINK9.
! This elems grid point displ's in global coord's, UGG, are gotten from UG_COL.

      I2 = 0
      DO I=1,ELGP
!        CALL CALC_TDOF_ROW_NUM ( AGRID(I), ROW_NUM_START, 'N' )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(I), IGRID )
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
            TDOF_ROW = ROW_NUM_START + J - 1               ! TDOF has rows in grid point numerical order
            GDOF = TDOF(TDOF_ROW,G_SET_COL)                ! This is the G-set DOF for BGRID(I), component J
            I2 = I2 + 1
            UGG(I2) = UG_COL(GDOF)
         ENDDO
      ENDDO

      IF (DEBUG(56) > 0) THEN
         WRITE(F06,5000)
         WRITE(F06,5001) TRIM(TYPE), EID, SCNUM(INT_SC_NUM)
         WRITE(F06,*)
         WRITE(F06,5101)
         WRITE(F06,5002) ' UGG FOR G.P. ', AGRID(1), ':       ', (UGG(I),I=1, 6)
         WRITE(F06,5002) ' UGG FOR G.P. ', AGRID(2), ':       ', (UGG(I),I=7,12)
         WRITE(F06,*)
      ENDIF
                                                           ! For all but ELAS, USERIN there is a 3 step process to get UEL from UGG
      IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  ')) THEN

!        ---------------------------------------------------------------------------------------------------------------------------
         I2 = 0                                            ! (1) Transform global displs at grids to global displs at elem nodes
         DO I=1,ELGP                                          ! First assume no offset at this node
            CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               I2 = I2 + 1
               UEG(I2) = UGG(I2)
            ENDDO
            IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
               IF (TYPE /= 'BUSH    ') THEN                !     BUSH local axes are et the grids, not at the BUSH location
                  IF (OFFSET(I) == 'Y') THEN               !     Elem is offset at this node so transform UGG to UEG
                     DXI = OFFDIS(I,1)
                     DYI = OFFDIS(I,2)
                     DZI = OFFDIS(I,3)
                     PROW = 6*(I-1) + 1
                     UEG(PROW)   = UGG(PROW)                     + DZI*UGG(PROW+4) - DYI*UGG(PROW+5)
                     UEG(PROW+1) = UGG(PROW+1) - DZI*UGG(PROW+3)                   + DXI*UGG(PROW+5)
                     UEG(PROW+2) = UGG(PROW+2) + DYI*UGG(PROW+3) - DXI*UGG(PROW+4)
                     UEG(PROW+3) = UGG(PROW+3)
                     UEG(PROW+4) = UGG(PROW+4)
                     UEG(PROW+5) = UGG(PROW+5)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         IF (DEBUG(56) > 0) THEN
            WRITE(F06,5102)
            WRITE(F06,5002) ' UEG for G.P. ', AGRID(1), ':       ', (UEG(I),I=1, 6)
            WRITE(F06,5002) ' UEG for G.P. ', AGRID(2), ':       ', (UEG(I),I=7,12)
            WRITE(F06,*)
         ENDIF

!        ---------------------------------------------------------------------------------------------------------------------------
         I2 = 0                                            ! (2) Transform from global (UEG) to basic (UEB) coords at elem nodes
         DO I=1,ELGP
            GLOBAL_CID = GRID(BGRID(I),3)                  !     GLOBAL_CID local coord sys exists. It was checked in CORD_PROC
            IF (GLOBAL_CID /= 0) THEN                      !     If global is not basic, do coord transformation
               DO J=1,NCORD
                  IF (CORD(J,2) == GLOBAL_CID) THEN
                     ICORD = J                             !     ICORD is the internal coord. sys. ID corresponding to GLOBAL_CID
                     EXIT
                  ENDIF
               ENDDO
               CALL GEN_T0L ( BGRID(I), ICORD, THETAD, PHID, T0G )
               DO J=1,2
                  PROW = 6*(I-1) + 1 + 3*(J-1)
                  CALL MATGET ( UEG,  6*MELGP, 1, PROW, PCOL, NROW, NCOL, DUM1 )
                  CALL MATMULT_FFF ( T0G, DUM1, NROWA, NCOLA, NCOLB, DUM2 )
                  CALL MATPUT ( DUM2, 6*MELGP, 1, PROW, PCOL, NROW, NCOL, UEB )
               ENDDO
               I2 = I2 + NUM_COMPS
            ELSE                                           ! If global is basic, get UEB terms directly from UEG
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               DO J=1,NUM_COMPS
                  I2 = I2 + 1
                  UEB(I2) = UEG(I2)
               ENDDO
            ENDIF
         ENDDO

         IF (DEBUG(56) > 0) THEN
            WRITE(F06,5103)
            WRITE(F06,5002) ' UEB for G.P. ', AGRID(1), ':       ', (UEB(I),I=1, 6)
            WRITE(F06,5002) ' UEB for G.P. ', AGRID(2), ':       ', (UEB(I),I=7,12)
            WRITE(F06,*)
         ENDIF

!        ---------------------------------------------------------------------------------------------------------------------------
         DO I=1,2*ELGP                                     ! (3) Transform from basic (UEB) to local (UEL) elem coords at elem nodes
            I1 = 3*I - 2
            UEL(I1)   = TE(1,1)*UEB(I1)+ TE(1,2)*UEB(I1+1)+ TE(1,3)*UEB(I1+2)
            UEL(I1+1) = TE(2,1)*UEB(I1)+ TE(2,2)*UEB(I1+1)+ TE(2,3)*UEB(I1+2)
            UEL(I1+2) = TE(3,1)*UEB(I1)+ TE(3,2)*UEB(I1+1)+ TE(3,3)*UEB(I1+2)
         ENDDO

         IF (DEBUG(56) > 0) THEN
            WRITE(F06,5104)
            WRITE(F06,5002) ' UEL for G.P. ', AGRID(1), ':       ', (UEL(I),I=1, 6)
            WRITE(F06,5002) ' UEL for G.P. ', AGRID(2), ':       ', (UEL(I),I=7,12)
            WRITE(F06,*)
            WRITE(F06,5000)
         ENDIF

      ELSE                                                 ! ELAS, USERIN have stiff, etc, in global coords so UEL = UGG

         DO I=1,ELDOF
            UEG(I) = UGG(I)
            UEB(I) = UGG(I)
            UEL(I) = UGG(I)
         ENDDO

      ENDIF

      IF (DEBUG(109) > 0) THEN
         CALL DEBUG_ELMDIS
      ENDIF



      RETURN

! **********************************************************************************************************************************

 5000 FORMAT('=============================================================================================================',/)

 5001 FORMAT(23X, 'S U B R O U T I N E   ELMDIS   F O R   ',A,I8,', S/C ',I8)

 5002 FORMAT(A,I3,A,12ES14.6)

 5101 FORMAT(' Displacements at grids in global coords',/,' ---------------------------------------')

 5102 FORMAT(' Displacements at elem nodes in global coords',/,' --------------------------------------------')

 5103 FORMAT(' Displacements at elem nodes in basic coords',/,' -------------------------------------------')

 5104 FORMAT(' Displacements at elem nodes in local elem coords',/,' ------------------------------------------------')


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_ELMDIS

      IMPLICIT NONE

! **********************************************************************************************************************************

      I2 = 0
      DO I=1,ELGP
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
            TDOF_ROW = ROW_NUM_START + J - 1
            GDOF = TDOF(TDOF_ROW,G_SET_COL)
            I2 = I2 + 1
         ENDDO
      ENDDO

      I2 = 0
      DO I=1,ELGP
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            I2 = I2 + 1
         ENDDO
      ENDDO

      I2 = 0
      DO I=1,ELGP
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            I2 = I2 + 1
         ENDDO
      ENDDO

      I2 = 0
      DO I=1,ELGP
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         DO J=1,NUM_COMPS
            I2 = I2 + 1
         ENDDO
      ENDDO

! **********************************************************************************************************************************








! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_ELMDIS

      END SUBROUTINE ELMDIS


      SUBROUTINE ELMDIS_PLY

! Get displs for one ply, or lamina, of a multi-ply element given the UG displs for the overall element (i.e. at the mid plane of
! the laminate). Result goes back into UEL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  CONV_DEG_RAD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ELGP, ELDOF, UEL, ZPLY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMDIS_PLY'

      INTEGER(LONG)                   :: I,j               ! DO loop index


      REAL(DOUBLE)                    :: DUM(6*ELGP)       ! Intermediate variable in the calculation of UEL for the ply



! **********************************************************************************************************************************
      DO I=1,ELGP
         DUM(6*(I-1)+1) = UEL(6*(I-1)+1) + ZPLY*UEL(6*(I-1)+5)
         DUM(6*(I-1)+2) = UEL(6*(I-1)+2) - ZPLY*UEL(6*(I-1)+4)
         DUM(6*(I-1)+3) = UEL(6*(I-1)+3)
         DUM(6*(I-1)+4) = UEL(6*(I-1)+4)
         DUM(6*(I-1)+5) = UEL(6*(I-1)+5)
         DUM(6*(I-1)+6) = UEL(6*(I-1)+6)
      ENDDO

      DO I=1,ELDOF
         UEL(I) = DUM(I)
      ENDDO




      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE ELMDIS_PLY


      SUBROUTINE TRANSFORM_NODE_FORCES ( COORD_SYS )

! Converts node forces for all elements from local to global or local to basic coords. The local to basic transformation is done
! every time this subr is called since that transformation must be done if either basic global is the final system anyway.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MELGP, NCORD
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  CAN_ELEM_TYPE_OFFSET, GRID, CORD, BGRID, ELDOF, ELGP, OFFDIS, OFFSET, PEB, PEG, PEL, TE,  &
                                         TYPE

      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE MATERIAL_MATRIX_TRANSFER, ONLY:  MATGET, MATPUT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TRANSFORM_NODE_FORCES'
      CHARACTER(LEN=*), INTENT(IN)    :: COORD_SYS         ! 'B" for basic, 'G' for global

      INTEGER(LONG)                   :: GLOBAL_CID        ! Global coord. sys. ID for a grid (BGRID(i)) of the element
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Calculated displ component no's for ELAS elems
      INTEGER(LONG)                   :: ICORD             ! Internal coord. system corresponding to GLOBAL_CID
      INTEGER(LONG)                   :: NROWS             ! DO loop limits when calculating elem nodal or engr forces
      INTEGER(LONG), PARAMETER        :: NCOLA     = 3     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NCOLB     = 1     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NROWA     = 3     ! An input to subr MATMULT_FFF called herein
      INTEGER(LONG), PARAMETER        :: NROW      = 3     ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG), PARAMETER        :: NCOL      = 1     ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG)                   :: PROW              ! An input to subr MATPUT, MATGET called herein
      INTEGER(LONG), PARAMETER        :: PCOL      = 1     ! An input to subr MATPUT, MATGET called herein


      REAL(DOUBLE)                    :: DXI               ! An offset distance in direction 1
      REAL(DOUBLE)                    :: DYI               ! An offset distance in direction 2
      REAL(DOUBLE)                    :: DZI               ! An offset distance in direction 3
      REAL(DOUBLE)                    :: T0G(3,3)          ! Coord transformation matrix - basic to global
      REAL(DOUBLE)                    :: DUM1(3)           ! Dummy arrays needed in transforming from global to basic coords
      REAL(DOUBLE)                    :: DUM2(3)           ! Dummy arrays needed in transforming from global to basic coords
      REAL(DOUBLE)                    :: DUM3(ELDOF)       ! Dummy arrays needed in transforming from global to basic coords
      REAL(DOUBLE)                    :: THETAD,PHID       ! Returns from subr GEN_T0L (not used here)



! **********************************************************************************************************************************
      NROWS = ELDOF

      DO I=1,NROWS
         PEB(I) = ZERO
         PEG(I) = ZERO
      ENDDO

      IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  ')) THEN
                                                           ! (1) Transform from local to basic. Use TE transpose to get PEB from PEL
         DO I=1,2*ELGP
            I1 = 3*I - 2
            PEB(I1)   = TE(1,1)*PEL(I1)+ TE(2,1)*PEL(I1+1)+ TE(3,1)*PEL(I1+2)
            PEB(I1+1) = TE(1,2)*PEL(I1)+ TE(2,2)*PEL(I1+1)+ TE(3,2)*PEL(I1+2)
            PEB(I1+2) = TE(1,3)*PEL(I1)+ TE(2,3)*PEL(I1+1)+ TE(3,3)*PEL(I1+2)
         ENDDO

         IF (COORD_SYS == 'G') THEN                        ! (2) Transform from basic to global if output is to be in global

            DO I=1,ELGP
               GLOBAL_CID = GRID(BGRID(I),3)               !     GLOBAL_CID local coord system exists. It was checked in CORD_PROC.
               IF (GLOBAL_CID /= 0) THEN                   !     If global is not basic, do coord transformation
                  DO J=1,NCORD
                     IF (CORD(J,2) == GLOBAL_CID) THEN
                        ICORD = J                          !     ICORD is the internal coord. sys. ID corresponding to GLOBAL_CID
                        EXIT
                     ENDIF
                  ENDDO
                  CALL GEN_T0L ( BGRID(I), ICORD, THETAD, PHID, T0G )
                  DO J=1,2
                     PROW = 6*(I-1) + 1 + 3*(J-1)
                     CALL MATGET ( PEB,  6*MELGP, 1, PROW, PCOL, NROW, NCOL, DUM1 )
                     CALL MATMULT_FFF_T ( T0G, DUM1, NROWA, NCOLA, NCOLB, DUM2 )
                     CALL MATPUT ( DUM2, 6*MELGP, 1, PROW, PCOL, NROW, NCOL, PEG )
                  ENDDO
               ELSE                                        !     If global is basic, get PEB terms directly from PEG
                  PROW = 6*(I-1) + 1
                  DO J=1,6
                     PEG(PROW+J-1) = PEB(PROW+J-1)
                  ENDDO
               ENDIF
            ENDDO

         ENDIF

      ELSE

         DO I=1,ELDOF
            PEG(I) = PEL(I)
            PEB(I) = PEL(I)
         ENDDO

      ENDIF

! Transform element loads at element ends to grids for BAR and ROD. Only need to do this if there are any offsets.
! Still use PEG to denote element forces (but now at grids, not elem ends)

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'ROD     ')) THEN
         DO I=1,ELGP
            PROW = 6*(I-1) + 1
            DO J=1,6
               DUM3(PROW+J-1) = PEG(PROW+J-1)
            ENDDO
            IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN
               IF (OFFSET(I) == 'Y') THEN                     ! Elem is offset at this node so transform PEG (using DUM3)
                  DXI = OFFDIS(I,1)
                  DYI = OFFDIS(I,2)
                  DZI = OFFDIS(I,3)
                  PROW = 6*(I-1) + 1
                  DUM3(PROW)   = PEG(PROW)
                  DUM3(PROW+1) = PEG(PROW+1)
                  DUM3(PROW+2) = PEG(PROW+2)
                  DUM3(PROW+3) = PEG(PROW+3)                   - DZI*PEG(PROW+1) + DYI*PEG(PROW+2)
                  DUM3(PROW+4) = PEG(PROW+4) + DZI*PEG(PROW)                     - DXI*PEG(PROW+2)
                  DUM3(PROW+5) = PEG(PROW+5) - DYI*PEG(PROW)   + DXI*PEG(PROW+1)
               ENDIF
            ENDIF
         ENDDO

         DO I=1,ELDOF                                      ! PEG is now at grids, not elem ends
            PEG(I) = DUM3(I)
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE TRANSFORM_NODE_FORCES


      SUBROUTINE TRANSFORM_SHELL_STR ( T, STR_VEC, SHR_FAC )

! Transforms a shell stress or strain vector with separate membrane, bending, and transverse shear terms.
! Use SHR_FAC = ONE for stress
! Use SHR_FAC = TWO for engineering strain to make it convert to tensor strain for transforming

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      REAL(DOUBLE),  INTENT(IN)       :: T(3,3)
      REAL(DOUBLE),  INTENT(INOUT)    :: STR_VEC(9)
      REAL(DOUBLE),  INTENT(IN)       :: SHR_FAC
      REAL(DOUBLE)                    :: STR_TENSOR(3,3)
      REAL(DOUBLE)                    :: DUM33(3,3)


! **********************************************************************************************************************************

                                                           ! Membrane and transverse shear
      STR_TENSOR(1,1) = STR_VEC(1)           ; STR_TENSOR(1,2) = STR_VEC(3) / SHR_FAC ;   STR_TENSOR(1,3) = STR_VEC(7) / SHR_FAC
      STR_TENSOR(2,1) = STR_VEC(3) / SHR_FAC ; STR_TENSOR(2,2) = STR_VEC(2)           ;   STR_TENSOR(2,3) = STR_VEC(8) / SHR_FAC
      STR_TENSOR(3,1) = STR_VEC(7) / SHR_FAC ; STR_TENSOR(3,2) = STR_VEC(8) / SHR_FAC ;   STR_TENSOR(3,3) = ZERO

      CALL MATMULT_FFF (STR_TENSOR, TRANSPOSE(T), 3, 3, 3, DUM33 )
      CALL MATMULT_FFF (T, DUM33, 3, 3, 3, STR_TENSOR )

      STR_VEC(1) = STR_TENSOR(1,1)
      STR_VEC(2) = STR_TENSOR(2,2)
      STR_VEC(3) = STR_TENSOR(1,2) * SHR_FAC
      STR_VEC(7) = STR_TENSOR(1,3) * SHR_FAC
      STR_VEC(8) = STR_TENSOR(2,3) * SHR_FAC

                                                           ! Bending
      STR_TENSOR(1,1) = STR_VEC(4)           ;   STR_TENSOR(1,2) = STR_VEC(6) / SHR_FAC ;   STR_TENSOR(1,3) = ZERO
      STR_TENSOR(2,1) = STR_VEC(6) / SHR_FAC ;   STR_TENSOR(2,2) = STR_VEC(5)           ;   STR_TENSOR(2,3) = ZERO
      STR_TENSOR(3,1) = ZERO                 ;   STR_TENSOR(3,2) = ZERO                 ;   STR_TENSOR(3,3) = ZERO

      CALL MATMULT_FFF (STR_TENSOR, TRANSPOSE(T), 3, 3, 3, DUM33 )
      CALL MATMULT_FFF (T, DUM33, 3, 3, 3, STR_TENSOR )

      STR_VEC(4) = STR_TENSOR(1,1)
      STR_VEC(5) = STR_TENSOR(2,2)
      STR_VEC(6) = STR_TENSOR(1,2) * SHR_FAC


      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE TRANSFORM_SHELL_STR


      SUBROUTINE STR_TENSOR_TRANSFORM ( STRESS_TENSOR, STRESS_CORD_SYS )

! Transforms an input stress or strain tensor in input local elem coord sys to an output coord sys whose actual coord sys ID is
! STRESS_CORD_SYS

! The 3x3 transform matrix Tse which transforms a unit vector in elem coords (Ue) to a unit vector in stress output coords (Us) is:

!                                                  Us = Tse*Ue                                                         (1)

! Tse can be written as the product of two transformation matrices Ts0 (transform vector in basic to output s coord sys) and
! the transpose of Te0 (which transforms a vector in basic coords to local elem coords).
! That is:

!                                                 Tse = Ts0*T0e = Ts0*Te0'                                             (2)
!  and its transpose is Tes:
!                                                 Tes = Te0*Ts0'                                                       (3)


! where Te0 is output in subr ELMGMi as array TE.  With (2), eqn (1) becomes:

!                                                  Us = (Ts0*Te0')*Ue                                                  (4)

! The 2nd order stress tensor in element local coords is a 3x3 matrix of the 9 (6 independent) stresses:

!                                                       | Sxx  Sxy  Sxz |
!                                                 See = | Syx  Syy  Syz |                                              (5)
!                                                 ~     | Szx  Szy  Szz |

! where Sij = Sji and, for example Sxx is the normal stress in the x direction, Sxy is the shear stress in the xy plane, etc.
! Matrix Tse is used to transform this 2nd order stress tensor from element coords to output coords (Soo) via the equation:

!                                                 Soo = Tse*See*Tse' = Tse'*See*Tes                                    (6)
!                                                 ~         ~              ~

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NCORD
      USE IOUNT1, ONLY                :  ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CORD, RCORD, TE

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'STR_TENSOR_TRANSFORM'

      INTEGER(LONG), INTENT(IN)       :: STRESS_CORD_SYS   ! Actual coord system ID for stress/strain/engr force output
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ICORD             ! Internal coord system ID for STRESS_CORD_SYS
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE), INTENT(INOUT)     :: STRESS_TENSOR(3,3)! 2D stress tensor (eqn 6 above)
      REAL(DOUBLE)                    :: DUM33(3,3)        ! Intermediate array used in calc outputs
      REAL(DOUBLE)                    :: TS0(3,3)          ! Transform matrix from basic coords to stress output coords
      REAL(DOUBLE)                    :: T0S(3,3)          ! TS0'
      REAL(DOUBLE)                    :: TES(3,3)          ! Transform matrix from local elem coords to stress output coords
      REAL(DOUBLE)                    :: TSE(3,3)          ! TES'



! **********************************************************************************************************************************
! Get transformation matrix T0S from stress coord sys to basic if it exists

      IF (STRESS_CORD_SYS == 0) THEN                       ! STRESS_CORD_SYS was basic so set T0S to the identity matrix

        DO I=1,3
          DO J=1,3
            T0S(I,J) = ZERO
          ENDDO
          T0S(I,I) = ONE
        ENDDO

      ELSE                                                 ! Need to transform ES mat'l matrix to basic coords

        ICORD = 0
        DO I=1,NCORD                                       ! Get the internal coord system ID for STRESS_CORD_SYS
          IF (STRESS_CORD_SYS == CORD(I,2)) THEN
            ICORD = I
            EXIT
          ENDIF
        ENDDO

        IF (ICORD > 0) THEN                                ! STRESS_CORD_SYS was found so do transformation

          DO I=1,3                                         ! Get TOS from RCORD array
            DO J=1,3
              K = 3 + 3*(I-1) + J
              TS0(I,J) = RCORD(ICORD,K)
              T0S(J,I) = TS0(I,J)
            ENDDO
          ENDDO

        ELSE                                               ! STRESS_CORD_SYS was not found leave tensor in local elem system

          WRITE(F06,9101) STRESS_CORD_SYS

        ENDIF

      ENDIF

! **********************************************************************************************************************************


      CALL MATMULT_FFF ( TE , T0S  , 3, 3, 3, TES )        ! Calc TES then transform to get TSE
      DO I=1,3
        DO J=1,3
          TSE(I,J) = TES(J,I)
        ENDDO
      ENDDO
                                                           ! Transform input stress tensor from e coords to o coords
      CALL MATMULT_FFF (STRESS_TENSOR, TES, 3, 3, 3, DUM33 )
      CALL MATMULT_FFF (TSE, DUM33, 3, 3, 3, STRESS_TENSOR )




      RETURN

! **********************************************************************************************************************************
 9101 FORMAT(' *INFORMATION: COORD SYSTEM ',I8,' FOR STRESS TRANSFORMATION IS UNDEFINED. LOCAL ELEMENT COORD SYSTEM WILL BE USED')

! **********************************************************************************************************************************

      END SUBROUTINE STR_TENSOR_TRANSFORM


   END MODULE RESULT_COORDINATES
