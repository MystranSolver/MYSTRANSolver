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

   MODULE ELEMENT_TRANSFORMATIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELEM_TRANSFORM_LBG

   CONTAINS

      SUBROUTINE ELEM_TRANSFORM_LBG ( WHICH, ZE, QE )

! Transforms one element stiff, mass, thermal load or pressure load matrix from local to basic to global coords at each
! grid including the effects of element offsets. The element matrix is input in array ZE or QE in local element
! coords. The output is array ZE or QE containing the element stiff or mass matrix in global coords at each grid.  Note that plate
! element offsets were processed in subr EMG since those offsets are in local element coords, while the BAR and BEAM offsets are
! handled here after transforming their matrices from local-basic-global, since BAR and BEAM offsets are specified (in the input
! data) in global coordinates and the BUSH in a unique system

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MELDOF, NCORD, NGRID, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  AGRID, CORD, ELDOF, GRID, GRID_ID, KEG, TE_IDENT, TYPE
      USE MODEL_STUF, ONLY            :  ELGP

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE ELEMENT_DIAGNOSTICS, ONLY   :  ELMOFF, ELMTLB
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE MATERIAL_MATRIX_TRANSFER, ONLY:  MATGET, MATPUT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELEM_TRANSFORM_LBG'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr ELMTLB (to tell it what to transform)
      CHARACTER(LEN=*), INTENT(IN)    :: WHICH             ! 'K' for stiffness, 'M' for mass or 'PTE' for thermal load matrix

      INTEGER(LONG)                   :: ACIDJ             ! Actual global coord sys ID for elem grid GRID_ID_ROW_NUM_J
      INTEGER(LONG)                   :: ACIDK             ! Actual global coord sys ID for elem grid GRID_ID_ROW_NUM_K
      INTEGER(LONG)                   :: BEG_ROW_GET       ! An input to subr MATGET/MATPUT (what row to start get/put rows)
      INTEGER(LONG)                   :: BEG_COL_GET       ! An input to subr MATGET/MATPUT (what col to start get/put rows)
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_J ! Row number in array GRID_ID where AGRID(J) is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_K ! Row number in array GRID_ID where AGRID(K) is found
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: ICID              ! Internal coord sys ID
      INTEGER(LONG), PARAMETER        :: NCOLA     = 3     ! An input to subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOLB             ! An input to subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOL_GET          ! An input to subr MATGET/MATPUT (no. cols to get/put)
      INTEGER(LONG)                   :: NCOL_IN           ! Number of cols in matrix being transformed
      INTEGER(LONG), PARAMETER        :: NROW_GET  = 3     ! An input to subr MATGET/MATPUT (no. rows to get/put)
      INTEGER(LONG), PARAMETER        :: NROWA     = 3     ! An input to subr MATMULT_FFF/MATMULT_FFF_T, called herein


      REAL(DOUBLE) , INTENT(INOUT)    :: QE(MELDOF,NSUB)   ! PTE or PPE if WHICH = 'PTE' or 'PPE'
      REAL(DOUBLE) , INTENT(INOUT)    :: ZE(MELDOF,MELDOF) ! Either the mass or stiff matrix of the element
      REAL(DOUBLE)                    :: DUM11(3,3)        ! An intermadiate matrix in a matrix multiply operation
      REAL(DOUBLE)                    :: DUM12(3,3)        ! An intermadiate matrix in a matrix multiply operation
      REAL(DOUBLE)                    :: DUM21(3,NSUB)     ! An intermadiate matrix in a matrix multiply operation
      REAL(DOUBLE)                    :: DUM22(3,NSUB)     ! An intermadiate matrix in a matrix multiply operation
      REAL(DOUBLE)                    :: THETAD,PHID       ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: TJ(3,3)           ! Coord transform matrix from basic to global for an internal
      REAL(DOUBLE)                    :: TK(3,3)           ! Coord transform matrix from basic to global for an internal grid



! **********************************************************************************************************************************

      IF ((TYPE(1:4) == 'ELAS') .OR. (TYPE == 'USERIN  ')) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1407) SUBR_NAME, ' ELASi, USERIN '
         WRITE(F06,1407) SUBR_NAME, ' ELASi, USERIN '
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Set up option flag for ENG

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      IF (WHICH == 'ME') THEN
         OPT(1) = 'Y'
      ELSE IF (WHICH == 'PTE') THEN
         OPT(2) = 'Y'
      ELSE IF (WHICH == 'KE') THEN
         OPT(4) = 'Y'
      ELSE IF (WHICH == 'PPE') THEN
         OPT(5) = 'Y'
      ELSE IF (WHICH == 'KED') THEN
         OPT(6) = 'Y'
      ELSE
         WRITE(ERR,1402) SUBR_NAME, WHICH
         WRITE(F06,1402) SUBR_NAME, WHICH
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
! Transform from local to basic coords (TE_IDENT = 'Y': TE is ident matrix). Note that ELAS elem is already in global coords

      IF (TE_IDENT /= 'Y') THEN
         CALL ELMTLB ( OPT )
      ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
! Transform from basic to global coords

! Double loop over j_do1:, k_do1: will transform KE, KED or ME from basic to global coords for all but ELAS element

ke_me:IF ((WHICH == 'KE') .OR. (WHICH == 'KED') .OR. (WHICH == 'ME')) THEN

         NCOL_IN = MELDOF
         NCOL_GET = 3
         NCOLB    = 3
j_do1:   DO J=1,ELGP
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(J), GRID_ID_ROW_NUM_J )
            ACIDJ = GRID(GRID_ID_ROW_NUM_J,3)
            IF (ACIDJ /= 0) THEN
k_cord1:       DO K=1,NCORD
                  IF (ACIDJ == CORD(K,2)) THEN
                     ICID = K
                     EXIT k_cord1
                  ENDIF
               ENDDO k_cord1
               CALL GEN_T0L ( GRID_ID_ROW_NUM_J, ICID, THETAD, PHID, TJ )
            ELSE
               DO K=1,3
                  DO L=1,3
                     TJ(K,L) = ZERO
                  ENDDO
                  TJ(K,K) = ONE
               ENDDO
            ENDIF

k_do1:      DO K=J,ELGP

               IF (K == J) THEN                            ! TK transf matrix is same as TJ
                  ACIDK = ACIDJ
                  DO L=1,3
                     DO M=1,3
                        TK(L,M) = TJ(L,M)
                     ENDDO
                  ENDDO
               ELSE
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(K), GRID_ID_ROW_NUM_K )
                  ACIDK = GRID(GRID_ID_ROW_NUM_K,3)
                  IF (ACIDK /= 0) THEN
l_cord1:             DO L=1,NCORD
                     IF (ACIDK == CORD(L,2)) THEN
                           ICID = L
                           EXIT l_cord1
                        ENDIF
                     ENDDO l_cord1
                     CALL GEN_T0L ( GRID_ID_ROW_NUM_K, ICID, THETAD, PHID, TK )
                  ELSE
                     DO L=1,3
                        DO M=1,3
                           TK(L,M) = ZERO
                        ENDDO
                        TK(L,L) = ONE
                     ENDDO
                  ENDIF
               ENDIF

               IF ((ACIDJ /= 0) .OR. (ACIDK /= 0)) THEN ! Do the coord transformation using TJ, TK coord transformation matrices
                  DO L=1,2
                     BEG_ROW_GET = 6*(J-1) + 1 + 3*(L-1)
                     DO M=1,2
                        BEG_COL_GET = 6*(K-1) + 1 + 3*(M-1)
                        CALL MATGET ( ZE, MELDOF, NCOL_IN, BEG_ROW_GET, BEG_COL_GET, NROW_GET, NCOL_GET, DUM11 )
                        CALL MATMULT_FFF   ( DUM11, TK, 3    , 3    , 3    , DUM12 )
                        CALL MATMULT_FFF_T ( TJ, DUM12, 3    , 3    , 3    , DUM11 )
                        CALL MATPUT ( DUM11, MELDOF, NCOL_IN, BEG_ROW_GET, BEG_COL_GET, NROW_GET, NCOL_GET, ZE )
                     ENDDO
                  ENDDO

               ENDIF

            ENDDO k_do1

         ENDDO j_do1

         DO I=2,ELDOF                                      ! Set lower triangular partition equal to upper partition
            DO J=1,I-1
               ZE(I,J) = ZE(J,I)
            ENDDO
         ENDDO

      ENDIF ke_me


pte_1:IF ((WHICH == 'PTE') .OR. (WHICH == 'PPE')) THEN


         IF (WHICH == 'PTE') THEN
            NCOL_IN  = NTSUB
            NCOL_GET = NTSUB
            NCOLB    = NTSUB
         ELSE
            NCOL_IN  = NSUB
            NCOL_GET = NSUB
            NCOLB    = NSUB
         ENDIF

j_do2:   DO J=1,ELGP

            BEG_ROW_GET = 6*(J-1) + 1
            BEG_COL_GET = 1
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(J), GRID_ID_ROW_NUM_J )
            ACIDJ = GRID(GRID_ID_ROW_NUM_J,3)
            IF (ACIDJ /= 0) THEN
k_cord2:       DO K=1,NCORD
                  IF (ACIDJ == CORD(K,2)) THEN
                     ICID = K
                     EXIT k_cord2
                  ENDIF
               ENDDO k_cord2
               CALL GEN_T0L ( GRID_ID_ROW_NUM_J, ICID, THETAD, PHID, TJ )
               DO K=1,2
                  CALL MATGET ( QE   , MELDOF, NCOL_IN, BEG_ROW_GET, BEG_COL_GET, NROW_GET, NCOL_GET, DUM21 )
                  CALL MATMULT_FFF_T ( TJ, DUM21, 3    , 3    , NCOLB, DUM22 )
                  CALL MATPUT ( DUM22, MELDOF, NCOL_IN, BEG_ROW_GET, BEG_COL_GET, NROW_GET, NCOL_GET, QE )
                  BEG_ROW_GET = BEG_ROW_GET + 3
               ENDDO
            ELSE
               DO K=1,3
                  DO L=1,3
                     TJ(K,L) = ZERO
                  ENDDO
                  TJ(K,K) = ONE
               ENDDO
            ENDIF

         ENDDO j_do2

      ENDIF pte_1

! Transform the matrix from global at elem ends to global at grids for BAR, BEAM, BUSH

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'BUSH    ')) THEN
         CALL ELMOFF ( OPT, 'N' )

         IF (WHICH == 'KE') THEN                           ! Set the KE which was just calc'd in global coords with offsets to KEG
            DO I=1,ELDOF
               DO J=1,ELDOF
                  KEG(I,J) = ZE(I,J)
               ENDDO
            ENDDO
         ENDIF


      ENDIF

      IF (TYPE == 'BUSH    ') THEN
         CALL GET_KE_OFFSET                                ! Now transform the global KE (with offsets) back to local
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1402 FORMAT(' *ERROR  1402: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' OPTION "WHICH" MUST BE EITHER KE OR KED (FOR ELEMENT STIFFNESS) OR ME (FOR ELEMENT MASS) MATRIX',     &
                           ' COORDINATE TRANSFORMATION BUT VALUE WAS ',A)

 1407 FORMAT(' *ERROR  1407: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THIS SUBR IS NOT PROGRAMMED FOR ELEM TPYES: ',A)



! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_KE_OFFSET

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  MELDOF
      USE MODEL_STUF, ONLY            :  ELGP, KEO_BUSH, TE, TE_GA_GB


      IMPLICIT NONE


      INTEGER(LONG)                   :: BEG_COL           ! Beginning col of matrix to get partition from
      INTEGER(LONG)                   :: BEG_ROW           ! Beginning row of matrix to get partition from
      INTEGER(LONG)                   :: II,JJ             ! DO loop indices
      INTEGER(LONG)                   :: NCOL              ! No. cols to get/put for subrs MATGET/MATPUT, called herein
      INTEGER(LONG), PARAMETER        :: NCOLA     = 3     ! No. cols in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOLB             ! No. cols in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG), PARAMETER        :: NROW      = 3     ! No. rows to get/put for subrs MATGET/MATPUT, called herein
      INTEGER(LONG), PARAMETER        :: NROWA     = 3     ! No. rows in a matrix for subr MATMULT_FFF/MATMULT_FFF_T, called herein

      REAL(DOUBLE)                    :: DUM11(3,3)        ! An intermediate result
      REAL(DOUBLE)                    :: DUM12(3,3)        ! An intermediate result
      REAL(DOUBLE)                    :: TET(3,3)          ! Transpose of TE
      REAL(DOUBLE)                    :: TET_GA_GB(3,3)    ! Transpose of TE

! **********************************************************************************************************************************
      DO II=1,3
         DO JJ = 1,3
            TET(II,JJ) = TE(II,JJ)
            TET_GA_GB(II,JJ) = TE_GA_GB(II,JJ)
         ENDDO
      ENDDO

      NCOL  = 3
      NCOLB = 3
      DO II=1,2*ELGP
         BEG_ROW = 3*II - 2
         DO JJ=II,2*ELGP
            BEG_COL = 3*JJ - 2
            CALL MATGET ( ZE, MELDOF, MELDOF, BEG_ROW, BEG_COL, NROW, NCOL, DUM11 )
            CALL MATMULT_FFF   ( DUM11, TET_GA_GB, NROWA, NCOLA, NCOLB, DUM12 )
            CALL MATMULT_FFF_T ( TE_GA_GB, DUM12, NROWA, NCOLA, NCOLB, DUM11 )
            CALL MATPUT ( DUM11, SIZE(KEO_BUSH, 2), SIZE(KEO_BUSH, 1), BEG_ROW, BEG_COL, NROW, NCOL, KEO_BUSH )
         ENDDO
      ENDDO
      DO II=1,ELDOF                                        ! Set lower portion of KE using symmetry.
         DO JJ=1,II-1
            KEO_BUSH(II,JJ) = KEO_BUSH(JJ,II)
         ENDDO
      ENDDO




! **********************************************************************************************************************************

      END SUBROUTINE GET_KE_OFFSET


      END SUBROUTINE ELEM_TRANSFORM_LBG

   END MODULE ELEMENT_TRANSFORMATIONS
