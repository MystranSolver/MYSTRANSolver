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

   MODULE RIGID_ELEMENT_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: RIGID_ELEM_PROC

   CONTAINS

      SUBROUTINE RIGID_ELEM_PROC

! Processes RBAR, RBE1, RBE2 rigid elements to get terms for the RMG constraint matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1F, LINK1F, L1F_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NRECARD
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RIGID_ELEM_PROC'
      CHARACTER( 8*BYTE)              :: RTYPE

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: IERR  = 0         ! Count of read errors when rigid elem data file is read
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Read a record from L1F and find out which rigid element it is

      IERR   = 0
      REC_NO = 0
      DO I=1,NRECARD

         READ(L1F,IOSTAT=IOCHK) RTYPE
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                                 ! Error reading RTYPE from rigid elem file. Can't continue
         ENDIF

         IF      (RTYPE == 'RBAR    ') THEN

            WRITE(ERR,*) '*ERROR: Code for RBAR processing not written yet. RBAR element ignored'
            WRITE(F06,*) '*ERROR: Code for RBAR processing not written yet. RBAR element ignored'
            WRITE(SC1,*) '*ERROR: Code for RBAR processing not written yet. RBAR element ignored'
            IERR      = IERR + 1
            FATAL_ERR = FATAL_ERR + 1

         ELSE IF (RTYPE == 'RBE1    ') THEN

            WRITE(ERR,*) '*ERROR: Code for RBE1 processing not written yet. RBE1 element ignored'
            WRITE(F06,*) '*ERROR: Code for RBE1 processing not written yet. RBE1 element ignored'
            WRITE(SC1,*) '*ERROR: Code for RBE1 processing not written yet. RBE1 element ignored'
            IERR      = IERR + 1
            FATAL_ERR = FATAL_ERR + 1
         ELSE IF (RTYPE == 'RBE2    ') THEN

            CALL RBE2_PROC    ( RTYPE, REC_NO, IERR )

         ELSE IF (RTYPE == 'RBE3    ') THEN

            CALL RBE3_PROC    ( RTYPE, REC_NO, IERR )

         ELSE IF (RTYPE == 'RSPLINE ') THEN

            CALL RSPLINE_PROC ( RTYPE, REC_NO, IERR )

         ENDIF

      ENDDO

      IF (IERR > 0) THEN
         WRITE(ERR,9996) SUBR_NAME,IERR
         WRITE(F06,9996) SUBR_NAME,IERR
         CALL OUTA_HERE ( 'Y' )                                    ! Errors reading rigid element data file, so quit
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

      END SUBROUTINE RIGID_ELEM_PROC


      SUBROUTINE RBE2_PROC ( RTYPE, REC_NO, IERR )

! Processes a single RBE2 rigid element, per call, to get terms for the RMG constraint matrix. When the Bulk data was read, the
! RBE2 input data was written to file LINK1F. In this subr, file LINK1F is read and RBE2 terms for array RMG are calculated and
! written to file LINK1J.  Later, in subr SPARSE_RMG, LINK1J will be read to create the sparse array RMG (of all rigid element and
! MPC coefficients) which will be used in LINK2 to reduce the G-set mass, stiffness and load matrices to the N-set.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F, LINK1F, L1F_MSG, L1J
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  GRID, RGRID, GRID_ID, CORD
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE DOF_SET_CONSTRUCTION, ONLY  :  RDOF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RBE2_PROC'
      CHARACTER( 8*BYTE), INTENT(IN)  :: RTYPE             ! The type of rigid element being processed (RBE2)
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF (= 1 if a digit 1-6 is in DDOF)

      INTEGER(LONG), INTENT(INOUT)    :: IERR              ! Count of errors in RIGID_ELEM_PROC
      INTEGER(LONG), INTENT(INOUT)    :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: AGRID_D           ! Dep   grid ID (actual) read from a record of file LINK1F
      INTEGER(LONG)                   :: AGRID_I           ! Indep grid ID (actual) read from a record of file LINK1F
      INTEGER(LONG)                   :: DDOF              ! Dep DOF's for AGRID_D
      INTEGER(LONG)                   :: ECORD_D           ! Global coord ID (actual) for grid AGRID_D
      INTEGER(LONG)                   :: ECORD_I           ! Global coord ID (actual) for grid AGRID_I
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D ! Row number in array GRID_ID where AGRID_D is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_I ! Row number in array GRID_ID where AGRID_I is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no., in TDOF array, of the G-set DOF list
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ICORD_D           ! Internal coord ID corresponding to ECORD_D
      INTEGER(LONG)                   :: ICORD_I           ! Internal coord ID corresponding to ECORD_I
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: IROW              ! A row number in matrix RDI_GLOBAL
      INTEGER(LONG)                   :: JERR              ! Local error count
      INTEGER(LONG)                   :: M_SET_COL_NUM     ! Col no., in TDOF array, of the M-set DOF list
      INTEGER(LONG)                   :: NUM_COMPS_D       ! 6 if AGRID_D is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: NUM_COMPS_I       ! 6 if AGRID_I is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REID              ! RBE2 elem ID read from file LINK1F
      INTEGER(LONG)                   :: RMG_COL_NUM       ! Col no. of a term in array RMG
      INTEGER(LONG)                   :: RMG_ROW_NUM       ! Row no. of a term in array RMG
      INTEGER(LONG)                   :: ROW_NUM           ! A row number in array TDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid


      REAL(DOUBLE)                    :: DELTA_0(3,3)      ! 3 x 3 matrix of diffs in coords bet dep & indep grids in basic coords
      REAL(DOUBLE)                    :: DUM1(3,3)         ! Intermediate result in obtaining RDI_GLOBAL
      REAL(DOUBLE)                    :: DUM2(3,3)         ! Intermediate result in obtaining RDI_GLOBAL
      REAL(DOUBLE)                    :: DUM3(3,3)         ! Intermediate result in obtaining RDI_GLOBAL
      REAL(DOUBLE)                    :: PHID, THETAD      ! Angles output from subr GEN_T0L, called herein but not needed here
      REAL(DOUBLE)                    :: RDI_GLOBAL(6,6)   ! RDI matrix (see explanation below)
      REAL(DOUBLE)                    :: T0G_D(3,3)        ! Transforms a dep   DOF vector in basic coords to global coords
      REAL(DOUBLE)                    :: T0G_I(3,3)        ! Transforms a indep DOF vector in basic coords to global coords



! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

      JERR = 0

      READ(L1F,IOSTAT=IOCHK) REID, AGRID_D, DDOF, AGRID_I

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )
      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I, GRID_ID_ROW_NUM_I )

      REC_NO = REC_NO + 1
      IF (IOCHK == 0) THEN
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_D, NUM_COMPS_D, SUBR_NAME )
         IF (NUM_COMPS_D /= 6) THEN
            IERR = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RBE2', REID, NUM_COMPS_D
            WRITE(F06,1951) 'RBE2', REID, NUM_COMPS_D
         ENDIF
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_I, NUM_COMPS_I, SUBR_NAME )
         IF (NUM_COMPS_I /= 6) THEN
            IERR = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RBE2', REID, NUM_COMPS_I
            WRITE(F06,1951) 'RBE2', REID, NUM_COMPS_I
         ENDIF
      ELSE
         CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

! Return if local error count > 0

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! We know that the indep and dep grids (AGRID_I and AGRID_D) exist. This was checked in subr DOF_PROC.
! Get the basic-to-global trensformation matrices for AGRID_D and AGRID_I

      ECORD_D = GRID(GRID_ID_ROW_NUM_D,3)
      IF (ECORD_D /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_D == CORD(I,2)) THEN
               ICORD_D = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_D, ICORD_D, THETAD, PHID, T0G_D )
      ELSE
         DO I=1,3
            DO J=1,3
               T0G_D(I,J) = ZERO
            ENDDO
            T0G_D(I,I) = ONE
         ENDDO
      ENDIF

      ECORD_I = GRID(GRID_ID_ROW_NUM_I,3)
      IF (ECORD_I /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_I == CORD(I,2)) THEN
               ICORD_I = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_I, ICORD_I, THETAD, PHID, T0G_I )
      ELSE
         DO I=1,3
            DO J=1,3
               T0G_I(I,J) = ZERO
            ENDDO
            T0G_I(I,I) = ONE
         ENDDO
      ENDIF

! Generate DELTA_0

      DO I=1,3
         DO J=1,3
            DELTA_0(I,J) = ZERO
         ENDDO
      ENDDO
      DELTA_0(1,2) =  (RGRID(GRID_ID_ROW_NUM_D,3) - RGRID(GRID_ID_ROW_NUM_I,3))
      DELTA_0(1,3) = -(RGRID(GRID_ID_ROW_NUM_D,2) - RGRID(GRID_ID_ROW_NUM_I,2))
      DELTA_0(2,1) = -DELTA_0(1,2)
      DELTA_0(2,3) =  (RGRID(GRID_ID_ROW_NUM_D,1) - RGRID(GRID_ID_ROW_NUM_I,1))
      DELTA_0(3,1) = -DELTA_0(1,3)
      DELTA_0(3,2) = -DELTA_0(2,3)

! The global constraint equations are
!                                           | UN |
!                     RMG x UG = [RMN | RMM]|    | = 0
!                                           | UM |

! where UM are all of the dependent DOF's and UN are all indep. DOF's

! However, we first have to cast them, for each rigid element, as:

!                                           | UD |
!                     RMG x UG = [ I  | RDI]|    | = 0
!                                           | UI |

! where UD are all dependent (or M-set) DOF's, but UI may contain M-set, as well as N-set, DOF's.

! Generate the 6 x 6 matrix of RDI coeffs for this rigid elem. The 6 x 6 consists of three nonzero 3 x 3 matrices:
!   1) The upper left  3 x 3 partition is (T0G_D)(transpose) x (T0G_I)
!   2) The upper right 3 x 3 partition is (T0G_D)(transpose) x (DELTA_0) x (T0G_I)
!   3) The lower left  3 x 3 partition is null
!   4) The lower right 3 x 3 partition is the same as the upper left 3 x 3 partition

! Multiply (T0G_D)(trans) x (T0G_I) to get DUM1 (upper left and lower right 3 x 3 partitions of the 6 x 6 RDI matrix
! of coeffs for this rigid elem)

      CALL MATMULT_FFF_T ( T0G_D, T0G_I, 3, 3, 3, DUM1 )

! Now get DUM3 (upper right 3 x 3 partition of the 6 x 6 RDI matrix of coefficients for this rigid element)

      CALL MATMULT_FFF   ( DELTA_0, T0G_I, 3, 3, 3, DUM2 ) ! Multiply DELTA_0 x (T0G_I)
      CALL MATMULT_FFF_T ( T0G_D  , DUM2 , 3, 3, 3, DUM3 ) ! Multiply (T0G_D)(transpose) x DUM2

! RDI_GLOBAL is the 6 x 6 matrix of RDI coeffs for this element

      DO I=1,3
         DO J=1,3
            RDI_GLOBAL(I  ,J  ) = -DUM1(I,J)
            RDI_GLOBAL(I  ,J+3) = -DUM3(I,J)
            RDI_GLOBAL(I+3,J  ) =  ZERO
            RDI_GLOBAL(I+3,J+3) = -DUM1(I,J)
         ENDDO
      ENDDO

! Now put all coeffs (I as well as RDI) into the RMG file (L1J). Due to the coord transformation (made above) from
! basic to global, we have to assume that all 6 columns in RDI_GLOBAL are not null. However, we only need the rows
! defined by the M-set DOF's in DDOF

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )
      CALL RDOF ( DDOF, CDOF )

      DO I=1,6

         IF (CDOF(I) == '1') THEN                          ! The I-th component is in DDOF so write this row to RMG

            IROW = I
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            ROW_NUM = ROW_NUM_START + I - 1
            RMG_ROW_NUM = TDOF(ROW_NUM,M_SET_COL_NUM)
            RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)

            IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
               WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,ONE
            ELSE
               IF (RMG_ROW_NUM  == 0) THEN
                  WRITE(ERR,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  IERR = IERR + 1
                  JERR = JERR + 1
               ENDIF
               IF (RMG_COL_NUM == 0) THEN
                  WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  IERR = IERR + 1
                  JERR = JERR + 1
               ENDIF
            ENDIF

            DO J=1,6                                       ! Now write RDI coefficients

               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I, IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               ROW_NUM = ROW_NUM_START + J - 1
               RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)

               IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
                  WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,RDI_GLOBAL(IROW,J)
               ELSE
                  IF (RMG_COL_NUM == 0) THEN
                     WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_I,J
                     WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_I,J
                     FATAL_ERR = FATAL_ERR + 1
                     IERR = IERR + 1
                     JERR = JERR + 1
                  ENDIF
               ENDIF

            ENDDO

         ENDIF

      ENDDO

! Return if JERR > 0

      IF (JERR > 0) THEN
         RETURN
      ENDIF

! If DEBUG(14) = 1, print results

      IF (DEBUG(14) == 1) THEN
         WRITE(F06,3001) RTYPE,REID,AGRID_D,ECORD_D,AGRID_I,ECORD_I
         WRITE(F06,*)
         WRITE(F06,3002)
         DO I=1,3
            WRITE(F06,3003) (DELTA_0(I,J),J=1,3)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3004)
         DO I=1,6
            WRITE(F06,3005) (RDI_GLOBAL(I,J),J=1,6)
         ENDDO
         WRITE(F06,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1509 FORMAT(' *ERROR  1509: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                          &
                    ,/,14X,' IS NOT A M-SET DOF IN TABLE TDOFI')

 1510 FORMAT(' *ERROR  1510: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', INDEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                        &
                    ,/,14X,' IS NOT A G-SET DOF IN TABLE TDOFI')

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')

 3001 FORMAT(/,' OUTPUT FROM SUBROUTINE RBE2_PROC FOR ',A8,' RIGID ELEMENT NUMBER ',I8                                             &
            ,/,' DEPENDENT   GRID IS NUMBER ',I8,' WITH GLOBAL COORD SYSTEM NUMBER ',I8                                            &
            ,/,' INDEPENDENT GRID IS NUMBER ',I8,' WITH GLOBAL COORD SYSTEM NUMBER ',I8)

 3002 FORMAT(1X,' 3 x 3 MATRIX DELTA_0 (DIFFERENCES OF COORDS OF DEP/INDEP GRIDS IN BASIC COORDS)')

 3003 FORMAT(1X,3(1ES15.6))

 3004 FORMAT(1X,' 6 x 6 MATRIX OF RDI MATRIX COEFFS IN GLOBAL COORDS FOR THIS RIGID ELEMENT')

 3005 FORMAT(1X,6(1ES15.6))

 3011 FORMAT(/,' In RIGID_ELEM_PROC: AGRID_D, GRID_ID_ROW_NUM_D, M_ROW, M_SET_COL_NUM, RMG_ROW_NUM                     = ',5(1X,I8))

 3012 FORMAT(/,' In RIGID_ELEM_PROC: AGRID_I, GRID_ID_ROW_NUM_I, N_ROW, N_SET_COL_NUM, RMG_COL_NUM, IROW, J RDI_GLOBAL = ',7I8,    &
                 1ES15.6)

! **********************************************************************************************************************************

      END SUBROUTINE RBE2_PROC



      SUBROUTINE RBE3_PROC ( RTYPE, REC_NO, IERR )

! Form the RMG constraint equations for one RBE3 element. Each active scalar independent DOF contributes one residual to a weighted
! rigid-body least-squares fit. If H is the kinematic vector for that DOF, its contribution is
!
!                  A6 = A6 + WEIGHT*H*TRANSPOSE(H)       and       B6(:,COL) = -WEIGHT*H .
!
! The six reference-grid components satisfy A6*Q + B6*U = 0. Components omitted by REFC are eliminated with a rank-revealing LAPACK
! solve before the retained equations are written to LINK1J. See Appendix E of the MYSTRAN User's Reference Manual.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1F, LINK1F, L1F_MSG, L1J
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MRBE3, NCORD, NGRID, NTERM_RMG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  CORD, GRID_ID, GRID, RGRID
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE DOF_SET_CONSTRUCTION, ONLY  :  RDOF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: NUM_RIGID_DOF = 6
      INTEGER(LONG), PARAMETER        :: NUM_VECTOR_DOF = 3
      INTEGER(LONG), PARAMETER        :: ALL_DOF(NUM_RIGID_DOF) = [1_LONG, 2_LONG, 3_LONG, 4_LONG, 5_LONG, 6_LONG]

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RBE3_PROC'
      CHARACTER( 8*BYTE), INTENT(IN)  :: RTYPE

      INTEGER(LONG), INTENT(INOUT)    :: IERR
      INTEGER(LONG), INTENT(INOUT)    :: REC_NO

      CHARACTER( 1*BYTE)              :: INDEP_DOF_ACTIVE(NUM_RIGID_DOF)
      CHARACTER( 1*BYTE)              :: REF_DOF_ACTIVE(NUM_RIGID_DOF)

      INTEGER(LONG)                   :: AGRID_D
      INTEGER(LONG)                   :: AGRID_I(MRBE3)
      INTEGER(LONG)                   :: B6_COL_RMG(NUM_RIGID_DOF*MRBE3)
      INTEGER(LONG)                   :: COMPS_D
      INTEGER(LONG)                   :: COMPS_I(MRBE3)
      INTEGER(LONG)                   :: DISCARDED_IDX(NUM_RIGID_DOF)
      INTEGER(LONG)                   :: G_SET_COL_NUM
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D
      INTEGER(LONG)                   :: IRBE3
      INTEGER(LONG)                   :: ITERM_RMG
      INTEGER(LONG)                   :: JERR
      INTEGER(LONG)                   :: M_SET_COL_NUM
      INTEGER(LONG)                   :: NB6COLS
      INTEGER(LONG)                   :: NUM_DISCARDED
      INTEGER(LONG)                   :: NUM_RETAINED
      INTEGER(LONG)                   :: REID
      INTEGER(LONG)                   :: RETAINED_IDX(NUM_RIGID_DOF)

      REAL(DOUBLE)                    :: A6(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                    :: A_REDUCED(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                    :: B6(NUM_RIGID_DOF,NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                    :: B_REDUCED(NUM_RIGID_DOF,NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: REFERENCE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                    :: T0D(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE)                    :: WEIGHT(MRBE3)
      REAL(DOUBLE)                    :: WEIGHT_SUM_ON_FILE

      INTERFACE
         SUBROUTINE DGELSY ( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
            IMPORT LONG, DOUBLE
            INTEGER(LONG), INTENT(IN)    :: M, N, NRHS, LDA, LDB, LWORK
            INTEGER(LONG), INTENT(INOUT) :: JPVT(*)
            INTEGER(LONG), INTENT(OUT)   :: RANK, INFO
            REAL(DOUBLE), INTENT(INOUT)  :: A(LDA,*), B(LDB,*), WORK(*)
            REAL(DOUBLE), INTENT(IN)     :: RCOND
         END SUBROUTINE DGELSY
      END INTERFACE

! **********************************************************************************************************************************

      EPS1 = EPSIL(1)
      JERR = 0

      CALL READ_RBE3_INPUT
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )

      CALL GET_GRID_TRANSFORM ( GRID_ID_ROW_NUM_D, T0D )
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF
      REFERENCE_POSITION = RGRID(GRID_ID_ROW_NUM_D,1:NUM_VECTOR_DOF)

      CALL ASSEMBLE_LEAST_SQUARES_SYSTEM
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL SELECT_REFERENCE_COMPONENTS
      CALL REDUCE_TO_REFERENCE_COMPONENTS
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL WRITE_RMG_TERMS ( ITERM_RMG )
      NTERM_RMG = NTERM_RMG + ITERM_RMG

      IF (JERR /= 0) FATAL_ERR = FATAL_ERR + 1

      RETURN

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE READ_RBE3_INPUT

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: NUM_COMPS
      INTEGER(LONG)                   :: OUNT(2)
      INTEGER(LONG)                   :: REID_LOCAL

      OUNT = [ERR, F06]
      AGRID_I = 0
      COMPS_I = 0
      WEIGHT = ZERO

! The element type record was read by RIGID_ELEM_PROC. Read the header and only use its values after IOSTAT has been checked.

      READ(L1F,IOSTAT=IOCHK) REID_LOCAL, AGRID_D, COMPS_D, IRBE3, WEIGHT_SUM_ON_FILE
      REC_NO = REC_NO + 1
      IF (IOCHK /= 0) THEN
         CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF
      REID = REID_LOCAL

      IF ((IRBE3 < 1) .OR. (IRBE3 > MRBE3)) THEN
         WRITE(ERR,1952) REID, IRBE3, MRBE3
         WRITE(F06,1952) REID, IRBE3, MRBE3
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )
      CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_D, NUM_COMPS, SUBR_NAME )
      IF (NUM_COMPS /= NUM_RIGID_DOF) THEN
         WRITE(ERR,1951) 'RBE3', REID, AGRID_D
         WRITE(F06,1951) 'RBE3', REID, AGRID_D
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

      DO I = 1, IRBE3
         READ(L1F,IOSTAT=IOCHK) AGRID_I(I), COMPS_I(I), WEIGHT(I)
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF
      ENDDO

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')
 1952 FORMAT(' *ERROR  1952: RBE3 ',I8,' HAS ',I8,' INDEPENDENT GRID RECORDS; THE SUPPORTED RANGE IS 1 THROUGH ',I8)

      END SUBROUTINE READ_RBE3_INPUT

! ##################################################################################################################################

      SUBROUTINE GET_GRID_TRANSFORM ( GRID_ROW, TRANSFORM )

      INTEGER(LONG), INTENT(IN)      :: GRID_ROW
      REAL(DOUBLE), INTENT(OUT)      :: TRANSFORM(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      INTEGER(LONG)                  :: COORD_ID
      INTEGER(LONG)                  :: COORD_ROW
      REAL(DOUBLE)                   :: PHI
      REAL(DOUBLE)                   :: THETA

      COORD_ID = GRID(GRID_ROW,3)
      IF (COORD_ID == 0) THEN
         TRANSFORM = IDENTITY_3()
         RETURN
      ENDIF

      COORD_ROW = FINDLOC(CORD(1:NCORD,2), COORD_ID, DIM=1)
      IF (COORD_ROW == 0) THEN
         WRITE(ERR,1953) REID, GRID_ID(GRID_ROW), COORD_ID
         WRITE(F06,1953) REID, GRID_ID(GRID_ROW), COORD_ID
         IERR = IERR + 1
         JERR = JERR + 1
         TRANSFORM = ZERO
         RETURN
      ENDIF

      CALL GEN_T0L ( GRID_ROW, COORD_ROW, THETA, PHI, TRANSFORM )

 1953 FORMAT(' *ERROR  1953: RBE3 ',I8,' USES GRID ',I8,' WITH UNRESOLVED DISPLACEMENT COORDINATE SYSTEM ',I8)

      END SUBROUTINE GET_GRID_TRANSFORM

! ##################################################################################################################################

      SUBROUTINE ASSEMBLE_LEAST_SQUARES_SYSTEM

      INTEGER(LONG)                  :: COMPONENT
      INTEGER(LONG)                  :: GRID_INDEX
      INTEGER(LONG)                  :: GRID_ROW
      INTEGER(LONG)                  :: NUM_COMPS
      INTEGER(LONG)                  :: RMG_COLUMN
      INTEGER(LONG)                  :: ROW_START
      REAL(DOUBLE)                   :: H(NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: RELATIVE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: T0I(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: TDI(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      A6 = ZERO
      B6 = ZERO
      B6_COL_RMG = 0
      NB6COLS = 0

      DO GRID_INDEX = 1, IRBE3
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I(GRID_INDEX), GRID_ROW )
         CALL GET_GRID_NUM_COMPS ( GRID_ROW, NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS /= NUM_RIGID_DOF) THEN
            WRITE(ERR,1951) 'RBE3', REID, AGRID_I(GRID_INDEX)
            WRITE(F06,1951) 'RBE3', REID, AGRID_I(GRID_INDEX)
            IERR = IERR + 1
            JERR = JERR + 1
            CYCLE
         ENDIF

         CALL GET_GRID_TRANSFORM ( GRID_ROW, T0I )
         IF (JERR /= 0) CYCLE

         RELATIVE_POSITION = MATMUL(TRANSPOSE(T0D), RGRID(GRID_ROW,1:NUM_VECTOR_DOF) - REFERENCE_POSITION)
         TDI = MATMUL(TRANSPOSE(T0D), T0I)
         ROW_START = TDOF_ROW_START(GRID_ROW)

         CALL RDOF ( COMPS_I(GRID_INDEX), INDEP_DOF_ACTIVE )
         DO COMPONENT = 1, NUM_RIGID_DOF
            IF (INDEP_DOF_ACTIVE(COMPONENT) /= '1') CYCLE

            RMG_COLUMN = TDOF(ROW_START + COMPONENT - 1,G_SET_COL_NUM)
            IF (RMG_COLUMN <= 0) THEN
               WRITE(ERR,1513) SUBR_NAME, AGRID_I(GRID_INDEX), COMPONENT, RMG_COLUMN
               WRITE(F06,1513) SUBR_NAME, AGRID_I(GRID_INDEX), COMPONENT, RMG_COLUMN
               IERR = IERR + 1
               JERR = JERR + 1
               CYCLE
            ENDIF

            H = KINEMATIC_VECTOR(COMPONENT, TDI, RELATIVE_POSITION)

            NB6COLS = NB6COLS + 1
            B6_COL_RMG(NB6COLS) = RMG_COLUMN
            B6(:,NB6COLS) = -WEIGHT(GRID_INDEX)*H
            A6 = A6 + WEIGHT(GRID_INDEX)*OUTER_PRODUCT_6(H)
         ENDDO
      ENDDO

 1513 FORMAT(' *ERROR  1513: PROGRAMMING ERROR IN SUBROUTINE ',A,/,14X,'G-SET COLUMN FOR GRID ',I8,', COMPONENT ',I2,              &
                    ' MUST BE > 0 BUT IS ',I8)
 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')

      END SUBROUTINE ASSEMBLE_LEAST_SQUARES_SYSTEM

! ##################################################################################################################################

      SUBROUTINE SELECT_REFERENCE_COMPONENTS

      LOGICAL                        :: IS_RETAINED(NUM_RIGID_DOF)

      CALL RDOF ( COMPS_D, REF_DOF_ACTIVE )
      IS_RETAINED = (REF_DOF_ACTIVE == '1')
      NUM_RETAINED = COUNT(IS_RETAINED)
      NUM_DISCARDED = NUM_RIGID_DOF - NUM_RETAINED

      RETAINED_IDX = 0
      DISCARDED_IDX = 0
      RETAINED_IDX(1:NUM_RETAINED) = PACK(ALL_DOF, IS_RETAINED)
      IF (NUM_DISCARDED > 0) DISCARDED_IDX(1:NUM_DISCARDED) = PACK(ALL_DOF, .NOT. IS_RETAINED)

      END SUBROUTINE SELECT_REFERENCE_COMPONENTS

! ##################################################################################################################################

      SUBROUTINE REDUCE_TO_REFERENCE_COMPONENTS

      INTEGER(LONG)                  :: INFO
      INTEGER(LONG)                  :: NUM_RIGHT_HAND_SIDES
      INTEGER(LONG)                  :: RANK
      REAL(DOUBLE)                   :: A_DD(5,5)
      REAL(DOUBLE)                   :: A_DD_ORIGINAL(5,5)
      REAL(DOUBLE)                   :: A_RD(NUM_RIGID_DOF,5)
      REAL(DOUBLE)                   :: RANK_CHECK_A(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: RANK_CHECK_B(NUM_RIGID_DOF,1)
      REAL(DOUBLE), PARAMETER        :: RANK_RCOND = 1.0E-10_DOUBLE ! Relative rank tolerance for the retained REFC system
      REAL(DOUBLE)                   :: RESIDUAL_SCALE
      REAL(DOUBLE)                   :: RHS(5,NUM_RIGID_DOF+NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                   :: RHS_ORIGINAL(5,NUM_RIGID_DOF+NUM_RIGID_DOF*MRBE3)

      A_REDUCED = ZERO
      B_REDUCED = ZERO

      IF (NUM_RETAINED == 0) THEN
         WRITE(ERR,1954) REID, COMPS_D
         WRITE(F06,1954) REID, COMPS_D
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF

      IF (NUM_DISCARDED == 0) THEN
         A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED) = A6(RETAINED_IDX(1:NUM_RETAINED),RETAINED_IDX(1:NUM_RETAINED))
         B_REDUCED(1:NUM_RETAINED,1:NB6COLS) = B6(RETAINED_IDX(1:NUM_RETAINED),1:NB6COLS)
      ELSE
         NUM_RIGHT_HAND_SIDES = NUM_RETAINED + NB6COLS
         A_DD = ZERO
         RHS = ZERO
         A_RD = ZERO

         A_DD(1:NUM_DISCARDED,1:NUM_DISCARDED) = &
            A6(DISCARDED_IDX(1:NUM_DISCARDED),DISCARDED_IDX(1:NUM_DISCARDED))
         A_RD(1:NUM_RETAINED,1:NUM_DISCARDED) = &
            A6(RETAINED_IDX(1:NUM_RETAINED),DISCARDED_IDX(1:NUM_DISCARDED))
         RHS(1:NUM_DISCARDED,1:NUM_RETAINED) = &
            A6(DISCARDED_IDX(1:NUM_DISCARDED),RETAINED_IDX(1:NUM_RETAINED))
         RHS(1:NUM_DISCARDED,NUM_RETAINED+1:NUM_RIGHT_HAND_SIDES) = &
            B6(DISCARDED_IDX(1:NUM_DISCARDED),1:NB6COLS)

         A_DD_ORIGINAL = A_DD
         RHS_ORIGINAL = RHS
         CALL SOLVE_MINIMUM_NORM ( A_DD, RHS, NUM_DISCARDED, NUM_RIGHT_HAND_SIDES, EPS1, RANK, INFO )
         IF (INFO /= 0) THEN
            WRITE(ERR,1955) REID, INFO
            WRITE(F06,1955) REID, INFO
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF

! A rank-deficient discarded block is permitted only when all elimination right-hand sides are in its range. DGELSY then supplies
! the unique minimum-norm solution needed for the generalized Schur complement.

         RESIDUAL_SCALE = MAX(ONE, MAXVAL(ABS(RHS_ORIGINAL(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES))))
         IF (MAXIMUM_RESIDUAL(A_DD_ORIGINAL(1:NUM_DISCARDED,1:NUM_DISCARDED), &
                              RHS(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES), &
                              RHS_ORIGINAL(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES)) > SQRT(EPS1)*RESIDUAL_SCALE) THEN
            WRITE(ERR,1956) REID, RANK, NUM_DISCARDED
            WRITE(F06,1956) REID, RANK, NUM_DISCARDED
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF

         A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED) = SCHUR_REDUCTION( &
            A6(RETAINED_IDX(1:NUM_RETAINED),RETAINED_IDX(1:NUM_RETAINED)), &
            A_RD(1:NUM_RETAINED,1:NUM_DISCARDED), RHS(1:NUM_DISCARDED,1:NUM_RETAINED))
         B_REDUCED(1:NUM_RETAINED,1:NB6COLS) = SCHUR_REDUCTION( &
            B6(RETAINED_IDX(1:NUM_RETAINED),1:NB6COLS), A_RD(1:NUM_RETAINED,1:NUM_DISCARDED), &
            RHS(1:NUM_DISCARDED,NUM_RETAINED+1:NUM_RIGHT_HAND_SIDES))
      ENDIF

! A singular retained block cannot define every requested dependent DOF. Diagnose it here instead of inventing unit pivot terms.

      RANK_CHECK_A = ZERO
      RANK_CHECK_B = ZERO
      RANK_CHECK_A(1:NUM_RETAINED,1:NUM_RETAINED) = A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED)
! The rank is judged with a relative tolerance of RANK_RCOND, not machine precision: a reduced system that is singular in exact
! arithmetic (e.g. two independent grids with only 123 whose line misses the reference grid: rotation about that line is free)
! comes out with sigma_min/sigma_max of 1.0E-16 to 1.0E-12 from round-off alone, so a machine-precision test accepted or
! rejected identical elements depending on their orientation.

      CALL SOLVE_MINIMUM_NORM ( RANK_CHECK_A, RANK_CHECK_B, NUM_RETAINED, 1_LONG, RANK_RCOND, RANK, INFO )
      IF ((INFO /= 0) .OR. (RANK < NUM_RETAINED)) THEN
         WRITE(ERR,1957) REID, NUM_RETAINED, RANK
         WRITE(F06,1957) REID, NUM_RETAINED, RANK
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

 1954 FORMAT(' *ERROR  1954: RBE3 ',I8,' HAS NO VALID REFC COMPONENTS IN VALUE ',I8)
 1955 FORMAT(' *ERROR  1955: LAPACK DGELSY FAILED WHILE REDUCING RBE3 ',I8,'; INFO = ',I8)
 1956 FORMAT(' *ERROR  1956: RBE3 ',I8,' HAS AN INCONSISTENT RANK-DEFICIENT DISCARDED-COMPONENT SYSTEM (RANK ',I2,' OF ',I2,')')
 1957 FORMAT(' *ERROR  1957: RBE3 ',I8,' CANNOT DETERMINE ALL ',I2,' REQUESTED REFC COMPONENTS; REDUCED RANK IS ',I2)

      END SUBROUTINE REDUCE_TO_REFERENCE_COMPONENTS

! ##################################################################################################################################

      SUBROUTINE SOLVE_MINIMUM_NORM ( MATRIX, RHS, N, NRHS, RCOND, RANK, INFO )

      REAL(DOUBLE), INTENT(INOUT)    :: MATRIX(:,:)
      REAL(DOUBLE), INTENT(INOUT)    :: RHS(:,:)
      INTEGER(LONG), INTENT(IN)      :: N
      INTEGER(LONG), INTENT(IN)      :: NRHS
      REAL(DOUBLE), INTENT(IN)       :: RCOND             ! Relative tolerance for the rank (DGELSY RCOND)
      INTEGER(LONG), INTENT(OUT)     :: RANK
      INTEGER(LONG), INTENT(OUT)     :: INFO

      INTEGER(LONG)                  :: JPVT(NUM_RIGID_DOF)
      INTEGER(LONG)                  :: LWORK
      REAL(DOUBLE), ALLOCATABLE      :: WORK(:)
      REAL(DOUBLE)                   :: WORK_QUERY(1)

      JPVT = 0
      CALL DGELSY ( N, N, NRHS, MATRIX, SIZE(MATRIX,1,KIND=LONG), RHS, SIZE(RHS,1,KIND=LONG), JPVT, RCOND, RANK, &
                    WORK_QUERY, -1_LONG, INFO )
      IF (INFO /= 0) RETURN

      LWORK = MAX(1_LONG, INT(WORK_QUERY(1),KIND=LONG))
      ALLOCATE(WORK(LWORK))
      JPVT = 0
      CALL DGELSY ( N, N, NRHS, MATRIX, SIZE(MATRIX,1,KIND=LONG), RHS, SIZE(RHS,1,KIND=LONG), JPVT, RCOND, RANK, WORK, LWORK, INFO )
      DEALLOCATE(WORK)

      END SUBROUTINE SOLVE_MINIMUM_NORM

! ##################################################################################################################################

      SUBROUTINE WRITE_RMG_TERMS ( TERM_COUNT )

      INTEGER(LONG), INTENT(OUT)     :: TERM_COUNT

      INTEGER(LONG)                  :: DEP_COMPONENT
      INTEGER(LONG)                  :: I
      INTEGER(LONG)                  :: J
      INTEGER(LONG)                  :: REF_COL(NUM_RIGID_DOF)
      INTEGER(LONG)                  :: REF_ROW_START
      INTEGER(LONG)                  :: RMG_ROW
      REAL(DOUBLE)                   :: COEFFICIENT_SCALE
      REAL(DOUBLE)                   :: WRITE_TOLERANCE

      TERM_COUNT = 0
      REF_COL = 0
      REF_ROW_START = TDOF_ROW_START(GRID_ID_ROW_NUM_D)

      DO I = 1, NUM_RETAINED
         DEP_COMPONENT = RETAINED_IDX(I)
         REF_COL(I) = TDOF(REF_ROW_START + DEP_COMPONENT - 1,G_SET_COL_NUM)
      ENDDO

      DO I = 1, NUM_RETAINED
         DEP_COMPONENT = RETAINED_IDX(I)
         RMG_ROW = TDOF(REF_ROW_START + DEP_COMPONENT - 1,M_SET_COL_NUM)
         IF ((RMG_ROW <= 0) .OR. (REF_COL(I) <= 0)) THEN
            IF (RMG_ROW <= 0) THEN
               WRITE(ERR,1509) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
               WRITE(F06,1509) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
            ENDIF
            IF (REF_COL(I) <= 0) THEN
               WRITE(ERR,1510) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
               WRITE(F06,1510) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
            ENDIF
            IERR = IERR + 1
            JERR = JERR + 1
            CYCLE
         ENDIF

         COEFFICIENT_SCALE = MAX(MAXVAL(ABS(A_REDUCED(I,1:NUM_RETAINED))), MAXVAL(ABS(B_REDUCED(I,1:NB6COLS))))
         WRITE_TOLERANCE = 10.0_DOUBLE*EPS1*COEFFICIENT_SCALE

         DO J = 1, NUM_RETAINED
            IF (ABS(A_REDUCED(I,J)) <= WRITE_TOLERANCE) CYCLE
            WRITE(L1J) RMG_ROW, REF_COL(J), A_REDUCED(I,J)
            TERM_COUNT = TERM_COUNT + 1
         ENDDO

         DO J = 1, NB6COLS
            IF (ABS(B_REDUCED(I,J)) <= WRITE_TOLERANCE) CYCLE
            WRITE(L1J) RMG_ROW, B6_COL_RMG(J), B_REDUCED(I,J)
            TERM_COUNT = TERM_COUNT + 1
         ENDDO
      ENDDO

 1509 FORMAT(' *ERROR  1509: PROGRAMMING ERROR IN SUBROUTINE ',A,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,                           &
                    ', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2,/,14X,' IS NOT AN M-SET DOF IN TABLE TDOFI')
 1510 FORMAT(' *ERROR  1510: PROGRAMMING ERROR IN SUBROUTINE ',A,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,                           &
                    ', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2,/,14X,' IS NOT A G-SET DOF IN TABLE TDOFI')

      END SUBROUTINE WRITE_RMG_TERMS

! ##################################################################################################################################

      PURE FUNCTION KINEMATIC_VECTOR ( COMPONENT, TDI, RELATIVE_POSITION ) RESULT ( H )

! Return the reference-motion coefficients for one scalar independent DOF. For translation, H = [D; R cross D]. For rotation,
! H = [0; D]. D is the independent component direction expressed in the reference grid's displacement coordinate system.

      INTEGER(LONG), INTENT(IN)      :: COMPONENT
      REAL(DOUBLE), INTENT(IN)       :: TDI(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE), INTENT(IN)       :: RELATIVE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: DIRECTION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: H(NUM_RIGID_DOF)

      H = ZERO
      IF (COMPONENT <= NUM_VECTOR_DOF) THEN
         DIRECTION = TDI(:,COMPONENT)
         H(1:NUM_VECTOR_DOF) = DIRECTION
         H(4:NUM_RIGID_DOF) = CROSS_PRODUCT_3(RELATIVE_POSITION, DIRECTION)
      ELSE
         DIRECTION = TDI(:,COMPONENT - NUM_VECTOR_DOF)
         H(4:NUM_RIGID_DOF) = DIRECTION
      ENDIF

      END FUNCTION KINEMATIC_VECTOR

! ##################################################################################################################################

      PURE FUNCTION SCHUR_REDUCTION ( RETAINED, COUPLING, ELIMINATED_SOLUTION ) RESULT ( REDUCED )

      REAL(DOUBLE), INTENT(IN)       :: RETAINED(:,:)
      REAL(DOUBLE), INTENT(IN)       :: COUPLING(:,:)
      REAL(DOUBLE), INTENT(IN)       :: ELIMINATED_SOLUTION(:,:)
      REAL(DOUBLE)                   :: REDUCED(SIZE(RETAINED,1),SIZE(RETAINED,2))

      REDUCED = RETAINED - MATMUL(COUPLING, ELIMINATED_SOLUTION)

      END FUNCTION SCHUR_REDUCTION

! ##################################################################################################################################

      PURE FUNCTION MAXIMUM_RESIDUAL ( MATRIX, SOLUTION, RHS ) RESULT ( MAX_RESIDUAL )

      REAL(DOUBLE), INTENT(IN)       :: MATRIX(:,:)
      REAL(DOUBLE), INTENT(IN)       :: SOLUTION(:,:)
      REAL(DOUBLE), INTENT(IN)       :: RHS(:,:)
      REAL(DOUBLE)                   :: MAX_RESIDUAL
      REAL(DOUBLE)                   :: RESIDUAL(SIZE(RHS,1),SIZE(RHS,2))

      RESIDUAL = MATMUL(MATRIX, SOLUTION) - RHS
      MAX_RESIDUAL = MAXVAL(ABS(RESIDUAL))

      END FUNCTION MAXIMUM_RESIDUAL

! ##################################################################################################################################

      PURE FUNCTION CROSS_PRODUCT_3 ( LEFT, RIGHT ) RESULT ( PRODUCT )

      REAL(DOUBLE), INTENT(IN)       :: LEFT(NUM_VECTOR_DOF)
      REAL(DOUBLE), INTENT(IN)       :: RIGHT(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: PRODUCT(NUM_VECTOR_DOF)

      PRODUCT = [ LEFT(2)*RIGHT(3) - LEFT(3)*RIGHT(2), &
                  LEFT(3)*RIGHT(1) - LEFT(1)*RIGHT(3), &
                  LEFT(1)*RIGHT(2) - LEFT(2)*RIGHT(1) ]

      END FUNCTION CROSS_PRODUCT_3

! ##################################################################################################################################

      PURE FUNCTION IDENTITY_3 () RESULT ( IDENTITY )

      REAL(DOUBLE)                   :: IDENTITY(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      IDENTITY = ZERO
      IDENTITY(1,1) = ONE
      IDENTITY(2,2) = ONE
      IDENTITY(3,3) = ONE

      END FUNCTION IDENTITY_3

! ##################################################################################################################################

      PURE FUNCTION OUTER_PRODUCT_6 ( VECTOR ) RESULT ( PRODUCT )

      REAL(DOUBLE), INTENT(IN)       :: VECTOR(NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: PRODUCT(NUM_RIGID_DOF,NUM_RIGID_DOF)

      PRODUCT = SPREAD(VECTOR,DIM=2,NCOPIES=NUM_RIGID_DOF)*SPREAD(VECTOR,DIM=1,NCOPIES=NUM_RIGID_DOF)

      END FUNCTION OUTER_PRODUCT_6

! **********************************************************************************************************************************

      END SUBROUTINE RBE3_PROC


      SUBROUTINE RSPLINE_PROC ( RTYPE, REC_NO, IERR )

! Processes a single RSPLINE rigid element, per call, to get terms for the RMG constraint matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F, L1F_MSG, LINK1F, L1J
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MRSPLINE, NCORD, NGRID, NTERM_RMG
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  CORD, GRID, RGRID, GRID_ID, CORD
      USE PARAMS, ONLY                :  EPSIL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, GEN_T0L
      USE DOF_SET_CONSTRUCTION, ONLY  :  RDOF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SORTING, ONLY               :  CALC_VEC_SORT_ORDER
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RSPLINE_PROC'
      CHARACTER( 8*BYTE), INTENT(IN)  :: RTYPE             ! The type of rigid element being processed (RSPLINE)
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF (= 1 if a digit 1-6 is in DDOF)

      INTEGER(LONG), INTENT(INOUT)    :: IERR              ! Count of errors in RIGID_ELEM_PROC
      INTEGER(LONG), INTENT(INOUT)    :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: AGRID_D           ! Dep   grid ID (actual) read for the RSPLINE
      INTEGER(LONG)                   :: AGRID_I1          ! Indep grid ID (actual) read for the RSPLINE
      INTEGER(LONG)                   :: AGRID_I2          ! Indep grid ID (actual) read for the RSPLINE
      INTEGER(LONG)                   :: COMPS_D           ! Components for dependent grids, AGRID_D(i)
      INTEGER(LONG)                   :: ECORD_D           ! Global coord ID (actual) for   dep grid AGRID_D
      INTEGER(LONG)                   :: ECORD_I1          ! Global coord ID (actual) for indep grid AGRID_I1
      INTEGER(LONG)                   :: ECORD_I2          ! Global coord ID (actual) for indep grid AGRID_I2
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D ! Row number in array GRID_ID where AGRID_D is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_I1! Row number in array GRID_ID where AGRID_I is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_I2! Row number in array GRID_ID where AGRID_I is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no., in TDOF array, of the G-set DOF list
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: ICORD             ! Internal coord ID
      INTEGER(LONG)                   :: INDEX             ! Index of records read from file LINK1F
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: IROW              ! A row number in matrix RDI_GLOBAL
      INTEGER(LONG)                   :: JERR              ! Local error count
      INTEGER(LONG)                   :: M_SET_COL_NUM     ! Col no., in TDOF array, of the M-set DOF list
      INTEGER(LONG)                   :: NUM_COMPS_D       ! 6 if AGRID_D  is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: NUM_COMPS_I1      ! 6 if AGRID_I1 is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: NUM_COMPS_I2      ! 6 if AGRID_I2 is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REID              ! RBE2 elem ID read from file LINK1F
      INTEGER(LONG)                   :: RMG_COL_NUM       ! Col no. of a term in array RMG
      INTEGER(LONG)                   :: RMG_ROW_NUM       ! Row no. of a term in array RMG
      INTEGER(LONG)                   :: ROW_NUM           ! A row number in array TDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TOTAL_NUM         ! Total number of records read for a single rigid element


      REAL(DOUBLE)                    :: DL_RAT            ! D/L ratio from the B.D. RSPLINE entry

      REAL(DOUBLE)                    :: FR(6,12)          ! Matrix of spline coefficients relating the 6 components of displ at the
!                                                            dependent grid to the 6 comps of displ at the 2 independent grids

      REAL(DOUBLE)                    :: FR11(3,3)         ! Partition of FR in 1st 3 rows and 1st 3 cols
      REAL(DOUBLE)                    :: FR12(3,3)         ! Partition of FR in 1st 3 rows and 2nd 3 cols
      REAL(DOUBLE)                    :: FR13(3,3)         ! Partition of FR in 1st 3 rows and 3rd 3 cols
      REAL(DOUBLE)                    :: FR14(3,3)         ! Partition of FR in 1st 3 rows and 4th 3 cols

      REAL(DOUBLE)                    :: FR21(3,3)         ! Partition of FR in 2nd 3 rows and 1st 3 cols
      REAL(DOUBLE)                    :: FR22(3,3)         ! Partition of FR in 2nd 3 rows and 2nd 3 cols
      REAL(DOUBLE)                    :: FR23(3,3)         ! Partition of FR in 2nd 3 rows and 3rd 3 cols
      REAL(DOUBLE)                    :: FR24(3,3)         ! Partition of FR in 2nd 3 rows and 4th 3 cols

      REAL(DOUBLE)                    :: L12               ! Length of RSPLINE between grids AGRID_I1 and AGRID_I2
      REAL(DOUBLE)                    :: L1D               ! Length of RSPLINE between grids AGRID_I1 and AGRID_D
      REAL(DOUBLE)                    :: PHID, THETAD      ! Angles output from subr GEN_T0L, called herein but not needed here
      REAL(DOUBLE)                    :: T0G_D(3,3)        ! Transforms a dep   DOF vector in basic coords to global coords
      REAL(DOUBLE)                    :: T0G_I1(3,3)       ! Transforms a indep DOF vector in basic coords to global coords
      REAL(DOUBLE)                    :: T0G_I2(3,3)       ! Transforms a indep DOF vector in basic coords to global coords

      REAL(DOUBLE)                    :: TRSPLINE(3,3)     ! Coord transformation matrix for RSPLINE that transforms a vector in
!                                                            basic coords to a vector in RSPLINE coords

      REAL(DOUBLE)                    :: V01(3)            ! Vector in basic coords from basic origin to AGRID_I1
      REAL(DOUBLE)                    :: V02(3)            ! Vector in basic coords from basic origin to AGRID_I2
      REAL(DOUBLE)                    :: V0D(3)            ! Vector in basic coords from basic origin to AGRID_D
      REAL(DOUBLE)                    :: V012(3)           ! Vector in basic coords from AGRID_I1 to AGRID_I2
      REAL(DOUBLE)                    :: V01D(3)           ! Vector in basic coords from AGRID_I1 to AGRID_D
      REAL(DOUBLE)                    :: ZETA              ! Nondimensional distance from AGRID_I1 to AGRID_D



! **********************************************************************************************************************************
! File LINK1F contains data from the logical RSPLINE cards in the input B.D. deck. For each logical RSPLINE card, LINK1F has:

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

      JERR = 0

! Make sure that the grids are all 6 components (i.e. no SPOINTS)

      READ(L1F,IOSTAT=IOCHK) REID, INDEX, TOTAL_NUM, AGRID_I1, AGRID_I2, AGRID_D, COMPS_D, DL_RAT

      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 1' )


! Get the internal grid ID's for AGRID I1, AGRID_I2 and AGRID_D. We know they exist

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I1, GRID_ID_ROW_NUM_I1 )
      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I2, GRID_ID_ROW_NUM_I2 )
      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D , GRID_ID_ROW_NUM_D )

      REC_NO = REC_NO + 1
      IF (IOCHK == 0) THEN
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_I1, NUM_COMPS_I1, SUBR_NAME )
         IF (NUM_COMPS_I1 /= 6) THEN
            IERR  = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RSPLINE', REID, NUM_COMPS_I1
            WRITE(F06,1951) 'RSPLINE', REID, NUM_COMPS_I1
         ENDIF
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_I2, NUM_COMPS_I2, SUBR_NAME )
         IF (NUM_COMPS_I2 /= 6) THEN
            IERR  = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RSPLINE', REID, NUM_COMPS_I2
            WRITE(F06,1951) 'RSPLINE', REID, NUM_COMPS_I2
         ENDIF
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_D, NUM_COMPS_D, SUBR_NAME )
         IF (NUM_COMPS_D /= 6) THEN
            IERR  = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RSPLINE', REID, NUM_COMPS_D
            WRITE(F06,1951) 'RSPLINE', REID, NUM_COMPS_D
         ENDIF
      ELSE
         CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

! Return if local error count > 0

      IF (JERR > 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF


! Get the distances from AGRID_I1 to AGRID_I2 (total length of the RSPLINE) and AGRID_I1 to AGRID_D

      DO I=1,3
         V01(I)  = RGRID(GRID_ID_ROW_NUM_I1,I)
         V02(I)  = RGRID(GRID_ID_ROW_NUM_I2,I)
         V0D(I)  = RGRID(GRID_ID_ROW_NUM_D ,I)
         V012(I) = V02(I) - V01(I)
         V01D(I) = V0D(I) - V01(I)
      ENDDO

      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 2' )

! Get the element-to-basic trensformation matrices for AGRID_D and AGRID_I1, AGRID_I2

      CALL RSPLINE_GEOM ( REID, AGRID_I1, AGRID_I2, V012, L12, TRSPLINE )
      IF (IERR > 0) THEN
         RETURN
      ENDIF
      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 3' )

      L1D  = DSQRT( V01D(1)*V01D(1) + V01D(2)*V01D(2) + V01D(3)*V01D(3) )
      ZETA = L1D/L12
      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 4' )

! Get the spline functions in element coords

      CALL RSPLINE_FUNCTIONS ( ZETA, L12, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24 )
      CALL ASSEMBLE_FR ( FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24, FR  )
      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 5' )

! Transform FR from element coords (x axis along line between the 2 indep grids) to basic coords

      CALL TRANSFORM_FR_E0 ( TRSPLINE, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24 )
      CALL ASSEMBLE_FR ( FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24, FR  )
      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 6' )

! Get the basic to global transformation matrices for the 2 indep grids and the dep grid

      ECORD_I1 = GRID(GRID_ID_ROW_NUM_I1,3)                ! Indep grid, AGRID_I1, basic to global transformation is T0G_I1
      IF (ECORD_I1 /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_I1 == CORD(I,2)) THEN
               ICORD = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_I1, ICORD, THETAD, PHID, T0G_I1 )
      ELSE
         DO I=1,3
            DO J=1,3
               T0G_I1(I,J) = ZERO
            ENDDO
            T0G_I1(I,I) = ONE
         ENDDO
      ENDIF

      ECORD_I2 = GRID(GRID_ID_ROW_NUM_I2,3)                ! Indep grid, AGRID_I2, basic to global transformation is T0G_I2
      IF (ECORD_I2 /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_I2 == CORD(I,2)) THEN
               ICORD = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_I2, ICORD, THETAD, PHID, T0G_I2 )
      ELSE
         DO I=1,3
            DO J=1,3
               T0G_I2(I,J) = ZERO
            ENDDO
            T0G_I2(I,I) = ONE
         ENDDO
      ENDIF

      ECORD_D = GRID(GRID_ID_ROW_NUM_D,3)                  ! Dep grid, AGRID_D, basic to global transformation is T0G_D
      IF (ECORD_D /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_D == CORD(I,2)) THEN
               ICORD = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_D, ICORD, THETAD, PHID, T0G_D )
      ELSE
         DO I=1,3
            DO J=1,3
               T0G_D(I,J) = ZERO
            ENDDO
            T0G_D(I,I) = ONE
         ENDDO
      ENDIF

      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 7' )

! Transform FR from basic coords to global coords at the 2 indep and 1 dep grid

      CALL TRANSFORM_FR_0G ( T0G_I1, T0G_I2, T0G_D, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24 )
      CALL ASSEMBLE_FR ( FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24, FR  )
      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( ' 8' )

! Write coefficients for the RMG constraint matrix to file L1J

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )
      CALL RDOF ( COMPS_D, CDOF )

      DO I=1,6

         IF (CDOF(I) == '1') THEN                          ! The I-th component is in DDOF so write this row to RMG

            IROW = I
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            ROW_NUM = ROW_NUM_START + I - 1
            RMG_ROW_NUM = TDOF(ROW_NUM,M_SET_COL_NUM)
            RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)
            IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
               WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,ONE
               NTERM_RMG = NTERM_RMG + 1
            ELSE IF (RMG_ROW_NUM  == 0) THEN
               WRITE(ERR,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
               WRITE(F06,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ELSE IF (RMG_COL_NUM == 0) THEN
               WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
               WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            DO J=1,6                                       ! Now write FR coefficients for AGRID_I1
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I1, IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               ROW_NUM = ROW_NUM_START + J - 1
               RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)
               IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
                  WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,-FR(IROW,J)
                  NTERM_RMG = NTERM_RMG + 1
               ELSE IF (RMG_ROW_NUM  == 0) THEN
                  WRITE(ERR,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ELSE IF (RMG_COL_NUM == 0) THEN
                  WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_I1,J
                  WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_I1,J
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDDO

            DO J=1,6                                       ! Now write FR coefficients for AGRID_I2
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I2, IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
               ROW_NUM = ROW_NUM_START + J - 1
               RMG_COL_NUM = TDOF(ROW_NUM,G_SET_COL_NUM)
               IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM > 0)) THEN
                  WRITE(L1J) RMG_ROW_NUM,RMG_COL_NUM,-FR(IROW,J+6)
               ELSE IF (RMG_ROW_NUM  == 0) THEN
                  WRITE(ERR,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ELSE IF (RMG_COL_NUM == 0) THEN
                  WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_I2,J
                  WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_I2,J
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDDO

         ENDIF

      ENDDO

! Write lines for end of debug output

      IF (DEBUG(111) > 0) CALL DEB_RSPLINE_PROC ( '99' )



      RETURN

! **********************************************************************************************************************************





 9999 FORMAT(' PROCESSING TERMINATED DUE TO ABOVE ERRORS')

 1509 FORMAT(' *ERROR  1509: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                          &
                    ,/,14X,' IS NOT A M-SET DOF IN TABLE TDOFI')

 1510 FORMAT(' *ERROR  1510: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', INDEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                        &
                    ,/,14X,' IS NOT A G-SET DOF IN TABLE TDOFI')

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')

 3001 FORMAT(/,' OUTPUT FROM SUBROUTINE RSPLINE_PROC FOR ',A8,' RIGID ELEMENT NUMBER ',I8                                          &
            ,/,' DEPENDENT   GRID IS NUMBER ',I8,' WITH GLOBAL COORD SYSTEM NUMBER ',I8                                            &
            ,/,' INDEPENDENT GRID IS NUMBER ',I8,' WITH GLOBAL COORD SYSTEM NUMBER ',I8)

 3002 FORMAT(1X,' 3 x 3 MATRIX DELTA_0 (DIFFERENCES OF COORDS OF DEP/INDEP GRIDS IN BASIC COORDS)')

 3003 FORMAT(1X,3(1ES15.6))

 3004 FORMAT(1X,' 6 x 6 MATRIX OF RDI MATRIX COEFFS IN GLOBAL COORDS FOR THIS RIGID ELEMENT')

 3005 FORMAT(1X,6(1ES15.6))

 3011 FORMAT(/,' In RIGID_ELEM_PROC: AGRID_D, GRID_ID_ROW_NUM_D, M_ROW, M_SET_COL_NUM, RMG_ROW_NUM                     = ',5(1X,I8))

 3012 FORMAT(/,' In RIGID_ELEM_PROC: AGRID_I, GRID_ID_ROW_NUM_I, N_ROW, N_SET_COL_NUM, RMG_COL_NUM, IROW, J RDI_GLOBAL = ',7I8,    &
                 1ES15.6)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE RSPLINE_GEOM ( REID, AGRID1, AGRID2, V012, LENGTH, TRSPLINE )

! Calculate an RSPLINE coord transformation matriox that transforms a vector in basic coords to one in element coords (x axis along
! the line between the 2 RSPLINE indep grids)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RSPLINE_GEOM'
      CHARACTER( 5*BYTE)              :: SORT_ORDER        ! Order in which the VX(i) have been sorted in subr CALC_VEC_SORT_ORDER

      INTEGER(LONG), INTENT(IN)       :: AGRID1            ! Indep grid number 1 on the RSPLINE
      INTEGER(LONG), INTENT(IN)       :: AGRID2            ! Indep grid number 2 on the RSPLINE
      INTEGER(LONG), INTENT(IN)       :: REID              ! Element ID
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I3_IN(3)          ! Integer array used in sorting VX.
      INTEGER(LONG)                   :: I3_OUT(3)         ! Integer array in sort order of VX_SORT. If VX is sorted sp that
!                                                             comp 2 is smallest then comp 3 then comp 1 then I3_OUT is 2, 3, 1


      REAL(DOUBLE) , INTENT(IN)       :: V012(3)           ! Vector in basic coords from 1st to 2nd indep grids on the RSPLINE

      REAL(DOUBLE) , INTENT(OUT)      :: TRSPLINE(3,3)     ! Coord transformation matrix for RSPLINE that transforms a vector in
!                                                            basic coords to a vector in RSPLINE coords
      REAL(DOUBLE) , INTENT(OUT)      :: LENGTH            ! Length of RSPLINE between the 2 independent grids
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: MAGY              ! Magnitude of vector VY
      REAL(DOUBLE)                    :: MAGZ              ! Magnitude of vector VZ
      REAL(DOUBLE)                    :: VX(3)             ! A vector in the elem x dir
      REAL(DOUBLE)                    :: VY(3)             ! A vector in the elem y dir
      REAL(DOUBLE)                    :: VZ(3)             ! A vector in the elem z dir



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize

      DO I=1,3
        DO J=1,3
           TRSPLINE(I,J) = ZERO
        ENDDO
      ENDDO

      DO I=1,3
         VX(I) = ZERO
      ENDDO

! Calculate a vector between ends of the element in basic coords

      VX(1) = V012(1)
      VX(2) = V012(2)
      VX(3) = V012(3)

! Length of element between ends is:

      LENGTH = DSQRT( VX(1)*VX(1) + VX(2)*VX(2) + VX(3)*VX(3) )

      IF (LENGTH < EPS1) THEN
         IERR = IERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1504) RTYPE, REID, AGRID1, AGRID2, LENGTH
         WRITE(F06,1504) RTYPE, REID, AGRID1, AGRID2, LENGTH
         RETURN
      ENDIF

! Unit vector in element X direction except for BUSH element with CID >= 0 (i.e. when BUSH does not have a V vector)

      DO I=1,3
         TRSPLINE(1,I) = VX(I)/LENGTH
      ENDDO

      DO I=1,3
         I3_IN(I)   = I
         I3_OUT(I)  = I3_IN(I)
      ENDDO
      CALL CALC_VEC_SORT_ORDER ( VX, SORT_ORDER, I3_OUT)! Use this rather than SORT_INT1_REAL1 - didn't work for vec 10., 0., 0.
      IF (SORT_ORDER == '     ') THEN                   ! Subr CALC_VEC_SORT_ORDER did not find a sort order
         FATAL_ERR = FATAL_ERR + 1
         IERR = IERR + 1
         WRITE(ERR,1944) SUBR_NAME, 'RSPLINE', REID
         WRITE(F06,1944) SUBR_NAME, 'RSPLINE', REID
         RETURN
      ENDIF
                                                        ! See notes on "Some Basic Vector Operations In IDL" on web site:
                                                        ! http://fermi.jhuapl.edu/s1r/idl/s1rlib/vectors/v_basic.html
      VY(I3_OUT(1)) =  ZERO                             !  (a) Component of VY in direction of min VX is set to zero
      VY(I3_OUT(2)) =  VX(I3_OUT(3))                    !  (b) Other 2 VY(i) are corresponding VX(i) switched with one x(-1)
      VY(I3_OUT(3)) = -VX(I3_OUT(2))
      MAGY  = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))

      IF (DABS(MAGY) < EPS1) THEN
         FATAL_ERR = FATAL_ERR + 1
         IERR = IERR + 1
         WRITE(ERR,1938) SUBR_NAME, 'Y', 'RSPLINE', REID, (VY(I),I=1,3)
         WRITE(F06,1938) SUBR_NAME, 'Y', 'RSPLINE', REID, (VY(I),I=1,3)
         RETURN
      ENDIF

      DO I=1,3
         TRSPLINE(2,I) = VY(I)/MAGY
      ENDDO

      CALL CROSS ( VX, VY, VZ )

      MAGZ  = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))

      IF (DABS(MAGZ) < EPS1) THEN
         FATAL_ERR = FATAL_ERR + 1
         IERR = IERR + 1
         WRITE(ERR,1938) SUBR_NAME, 'Z', 'RSPLINE', REID, (VZ(I),I=1,3)
         WRITE(F06,1938) SUBR_NAME, 'Z', 'RSPLINE', REID, (VZ(I),I=1,3)
         RETURN
      ENDIF

      DO I=1,3
         TRSPLINE(3,I) = VZ(I)/MAGZ
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1504 FORMAT(' *ERROR  1504: ',A, I8,' HAS LENGTH BETWEEN INDEPENDENT GRIDS ',I8,' AND ',I8,' = ',1ES10.2,' MUST BE > 0')

 1938 FORMAT(' *ERROR  1938: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT CALCULATE PROPER VECTOR IN ELEMENT LOCAL COORDINATE DIRECTION ',A,' FOR ',A,' ELEMENT ',I8,'.' &
                    ,/,14X,' THE VECTOR COMPONENTS CALCULATED WERE ',3(1ES14.6))

 1944 FORMAT(' *ERROR  1944: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE VX VECTOR FOR ',A,' ELEMENT ',I8,' WAS LEFT UNSORTED. IT MUST BE SORTED TO DETERMINE VY, VZ')

! **********************************************************************************************************************************

      END SUBROUTINE RSPLINE_GEOM

! ##################################################################################################################################

      SUBROUTINE RSPLINE_FUNCTIONS ( Z, L12, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24  )

! Calculate thespline functions for an RSPLINE element in element coords (x axis along the line between the 2 indep grids)

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ONE, TWO, THREE, FOUR, SIX
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RSPLINE_FUNCTIONS'

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE) , INTENT(IN)       :: L12               ! Length of RSPLINE between the 2 independent grids
      REAL(DOUBLE) , INTENT(IN)       :: Z                 ! Nondim distance to the RSPLINE dependent grid from the 1st indep grid

      REAL(DOUBLE) , INTENT(OUT)      :: FR11(3,3)         ! Partition of FR in 1st 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR12(3,3)         ! Partition of FR in 1st 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR13(3,3)         ! Partition of FR in 1st 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR14(3,3)         ! Partition of FR in 1st 3 rows and 4th 3 cols

      REAL(DOUBLE) , INTENT(OUT)      :: FR21(3,3)         ! Partition of FR in 2nd 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR22(3,3)         ! Partition of FR in 2nd 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR23(3,3)         ! Partition of FR in 2nd 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(OUT)      :: FR24(3,3)         ! Partition of FR in 2nd 3 rows and 4th 3 cols
      REAL(DOUBLE)                    :: Z_2               ! Z squared
      REAL(DOUBLE)                    :: Z_3               ! Z cubed



! **********************************************************************************************************************************
! Initialize

      DO I=1,3
         DO J=1,3
            FR11(I,J) = ZERO   ;   FR12(I,J) = ZERO   ;   FR13(I,J) = ZERO   ;   FR14(I,J) = ZERO
            FR21(I,J) = ZERO   ;   FR22(I,J) = ZERO   ;   FR23(I,J) = ZERO   ;   FR24(I,J) = ZERO
         ENDDO
      ENDDO

      Z_2 = Z*Z
      Z_3 = Z*Z_2

! 1st 3 rows of FR

      FR11(1,1) = ONE - Z
      FR11(2,2) = ONE - THREE*Z_2 + TWO*Z_3
      FR11(3,3) = FR11(2,2)

      FR12(2,3) =  (Z - TWO*Z_2 + Z_3)*L12
      FR12(3,2) = -FR12(2,3)

      FR13(1,1) = Z
      FR13(2,2) = THREE*Z_2 - TWO*Z_3
      FR13(3,3) = FR13(2,2)

      FR14(2,3) = (Z_3 - Z_2)*L12
      FR14(3,2) = -FR14(2,3)

! 2nd 3 rows of FR

      FR21(2,3) = -(Z_2 - Z)*SIX/L12
      FR21(3,2) = -FR21(2,3)

      FR22(1,1) = ONE - Z
      FR22(2,2) = ONE - FOUR*Z +THREE*Z_2
      FR22(3,3) = FR22(2,2)

      FR23(2,3) = (Z_2 - Z)*SIX/L12
      FR23(3,2) = -FR23(2,3)

      FR24(1,1) = Z
      FR24(2,2) = THREE*Z_2 - TWO*Z
      FR24(3,3) = FR24(2,2)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE RSPLINE_FUNCTIONS

! ##################################################################################################################################

      SUBROUTINE ASSEMBLE_FR ( FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24, FR  )

! Assemble the 8 3x3 FRij matrices into the 6x12 matrix FR

      IMPLICIT NONE

      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE) , INTENT(OUT)      :: FR(6,12)          ! Matrix of spline coefficients relating the 6 components of displ at the
!                                                            dependent grid to the 6 comps of displ at the 2 independent grids

      REAL(DOUBLE) , INTENT(IN)       :: FR11(3,3)         ! Partition of FR in 1st 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR12(3,3)         ! Partition of FR in 1st 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR13(3,3)         ! Partition of FR in 1st 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR14(3,3)         ! Partition of FR in 1st 3 rows and 4th 3 cols

      REAL(DOUBLE) , INTENT(IN)       :: FR21(3,3)         ! Partition of FR in 2nd 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR22(3,3)         ! Partition of FR in 2nd 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR23(3,3)         ! Partition of FR in 2nd 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(IN)       :: FR24(3,3)         ! Partition of FR in 2nd 3 rows and 4th 3 cols

! **********************************************************************************************************************************
! Put 3x3 matrices Fij into 6x12 matrix FR

      DO I=1,3
         DO J=1,3
            FR(I  ,J  ) = FR11(I,J)
            FR(I  ,J+3) = FR12(I,J)
            FR(I  ,J+6) = FR13(I,J)
            FR(I  ,J+9) = FR14(I,J)
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,3
            FR(I+3,J  ) = FR21(I,J)
            FR(I+3,J+3) = FR22(I,J)
            FR(I+3,J+6) = FR23(I,J)
            FR(I+3,J+9) = FR24(I,J)
         ENDDO
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE ASSEMBLE_FR

! ##################################################################################################################################

      SUBROUTINE TRANSFORM_FR_E0 ( TR, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24 )

! Transform FR matrix for an RSPLINE from element coords (x axis along lin e between the 2 indep grids) to basic coords

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: TR(3,3)           ! Coord transformation matrix for RSPLINE that transforms a vector in
!                                                            basic coords to a vector in RSPLINE coords

      REAL(DOUBLE) , INTENT(INOUT)    :: FR11(3,3)         ! Partition of FR in 1st 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR12(3,3)         ! Partition of FR in 1st 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR13(3,3)         ! Partition of FR in 1st 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR14(3,3)         ! Partition of FR in 1st 3 rows and 4th 3 cols

      REAL(DOUBLE) , INTENT(INOUT)    :: FR21(3,3)         ! Partition of FR in 2nd 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR22(3,3)         ! Partition of FR in 2nd 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR23(3,3)         ! Partition of FR in 2nd 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR24(3,3)         ! Partition of FR in 2nd 3 rows and 4th 3 cols

      REAL(DOUBLE)                    :: DUM1(3,3)         ! Intermediate matrix

! **********************************************************************************************************************************
      CALL MATMULT_FFF   ( FR11, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR11 )

      CALL MATMULT_FFF   ( FR12, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR12 )

      CALL MATMULT_FFF   ( FR13, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR13 )

      CALL MATMULT_FFF   ( FR14, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR14 )

      CALL MATMULT_FFF   ( FR21, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR21 )

      CALL MATMULT_FFF   ( FR22, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR22 )

      CALL MATMULT_FFF   ( FR23, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR23 )

      CALL MATMULT_FFF   ( FR24, TR, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR, DUM1, 3, 3, 3, FR24 )

! **********************************************************************************************************************************

      END SUBROUTINE TRANSFORM_FR_E0

! ##################################################################################################################################

      SUBROUTINE TRANSFORM_FR_0G ( TR_I1, TR_I2, TR_D, FR11, FR12, FR13, FR14, FR21, FR22, FR23, FR24 )

! Transform FR matrix for an RSPLINE from basic to global coords

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: TR_I1(3,3)        ! Matrix that transforms a vector in global to basic for 1st indep grid
      REAL(DOUBLE) , INTENT(IN)       :: TR_I2(3,3)        ! Matrix that transforms a vector in global to basic for 2nd indep grid
      REAL(DOUBLE) , INTENT(IN)       :: TR_D(3,3)         ! Matrix that transforms a vector in global to basic for dep grid

      REAL(DOUBLE) , INTENT(INOUT)    :: FR11(3,3)         ! Partition of FR in 1st 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR12(3,3)         ! Partition of FR in 1st 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR13(3,3)         ! Partition of FR in 1st 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR14(3,3)         ! Partition of FR in 1st 3 rows and 4th 3 cols

      REAL(DOUBLE) , INTENT(INOUT)    :: FR21(3,3)         ! Partition of FR in 2nd 3 rows and 1st 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR22(3,3)         ! Partition of FR in 2nd 3 rows and 2nd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR23(3,3)         ! Partition of FR in 2nd 3 rows and 3rd 3 cols
      REAL(DOUBLE) , INTENT(INOUT)    :: FR24(3,3)         ! Partition of FR in 2nd 3 rows and 4th 3 cols

      REAL(DOUBLE)                    :: DUM1(3,3)         ! Intermediate matrix

! **********************************************************************************************************************************
      CALL MATMULT_FFF   ( FR11, TR_I1, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR11 )

      CALL MATMULT_FFF   ( FR12, TR_I1, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR12 )

      CALL MATMULT_FFF   ( FR13, TR_I1, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR13 )

      CALL MATMULT_FFF   ( FR14, TR_I1, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR14 )

      CALL MATMULT_FFF   ( FR21, TR_I2, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR21 )

      CALL MATMULT_FFF   ( FR22, TR_I2, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR22 )

      CALL MATMULT_FFF   ( FR23, TR_I2, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR23 )

      CALL MATMULT_FFF   ( FR24, TR_I2, 3, 3, 3, DUM1 )
      CALL MATMULT_FFF_T ( TR_D, DUM1 , 3, 3, 3, FR24 )

! **********************************************************************************************************************************

      END SUBROUTINE TRANSFORM_FR_0G

! ##################################################################################################################################

      SUBROUTINE DEB_RSPLINE_PROC ( WHAT )

! Write debug info for RSPLINE rigid element code

      IMPLICIT NONE

      CHARACTER( 2*BYTE)              :: WHAT              ! Indicator of what to write out

      INTEGER(LONG)                   :: II,JJ             ! DO loop indices

! **********************************************************************************************************************************
      IF      (WHAT == ' 1' ) THEN
         WRITE(F06,1001)
         WRITE(F06,1002) RTYPE, REID, AGRID_D, AGRID_I1, AGRID_I2
         WRITE(F06,*)

      ELSE IF (WHAT == ' 2' ) THEN
         WRITE(F06,2001) 'indep', AGRID_I1, (V01(II) ,II=1,3)
         WRITE(F06,2001) 'indep', AGRID_I2, (V02(II) ,II=1,3)
         WRITE(F06,2001) '  dep', AGRID_D , (V0D(II) ,II=1,3)
         WRITE(F06,*)
         WRITE(F06,2002) AGRID_I1,'indep', AGRID_I2,(V012(II),II=1,3)
         WRITE(F06,2002) AGRID_I1,'  dep', AGRID_D ,(V01D(II),II=1,3)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 3' ) THEN
         WRITE(F06,3001)
         DO II=1,3
            WRITE(F06,3002) (TRSPLINE(II,JJ),JJ=1,3)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 4' ) THEN
         WRITE(F06,4001) L12, ZETA
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 5' ) THEN
         WRITE(F06,5001)
         WRITE(F06,5002)
         DO II=1,3
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         DO II=4,6
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 6' ) THEN
         WRITE(F06,5001)
         WRITE(F06,6002)
         DO II=1,3
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         DO II=4,6
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 7' ) THEN
         WRITE(F06,7001) AGRID_I1, ECORD_I1, AGRID_I2, ECORD_I2, AGRID_D, ECORD_D
         DO II=1,3
            WRITE(F06,7002) (T0G_I1(II,JJ),JJ=1,3), (T0G_I2(II,JJ),JJ=1,3), (T0G_D(II,JJ),JJ=1,3)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == ' 8' ) THEN
         WRITE(F06,5001)
         WRITE(F06,8002)
         DO II=1,3
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         DO II=4,6
            WRITE(F06,5003) (FR(II,JJ),JJ=1,3), (FR(II,JJ),JJ=4,6), (FR(II,JJ),JJ=7,9), (FR(II,JJ),JJ=10,12)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHAT == '99' ) THEN
         WRITE(f06,9001)

      ENDIF

! **********************************************************************************************************************************
 1001 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::::START DEBUG(111) OUTPUT FROM SUBROUTINE RSPLINE_PROC::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1002 FORMAT(14X,A,I9,' with dependent grid ',I8,' connects between independent grids ',I9,' and ',I8)

 2001 FORMAT(1X,'Vector in basic coords from basic system origin to ',A,' grid ',I8,':',3(1ES14.6))

 2002 FORMAT(1X,'Vector in basic coords from indep grid ',I8,' to ',A,' grid ',I8,':',3(1ES14.6))

 3001 FORMAT(19X,'Coord transformation matrix TRSPLINE (transforms a vector to element coords from basic coords):',/)

 3002 FORMAT(45X,3(1ES14.6))

 4001 FORMAT(1X,'Length of the RSPLINE (distance between the 2 independent grids)      :',1ES14.6,/,                               &
             1x,'Nondimensional distance from the 1st indep grid to the dep grid       :',1ES14.6)

 5001 FORMAT(38X,'                       Matrix FR of RSPLINE coefficients for the dependent grid',/,                              &
             38x,'Ud = FR*Ui where Ud are the 6 DOF''s at the dep grid and Ui are 12 DOF''s; 6 at each of the 2 indep grids')

 5002 FORMAT(56X,'(coord system is elem coords with x axis between the 2 indep grids)',/)

 5003 FORMAT(3(1ES14.6),'  |', 3(1ES14.6),'  |', 3(1ES14.6),'  |', 3(1ES14.6))

 6002 FORMAT(56X,'                      (coord system is basic)',/)

 7001 FORMAT(21X,'Coord transformation matrices that transform a vector to basic coords from global coords for:',/,                &
             21X,'--------------------------------------------------------------------------------------------' ,/,                &
             3X,'Indep grid',I8,', global CID',I8,7X,'Indep grid',I8,', global CID',I8,8X,'Dep grid',I8,', global CID',I8,/,       &
             3X,'    (transformation matrix T0G_I1)',15X,'(transformation matrix T0G_I2)',15X,'(transformation matrix T0G_D)',/)

 7002 FORMAT(3(3(1ES14.6),3X))

 8002 FORMAT(56X,'                      (coord system is global)',/)

 9001 FORMAT(' :::::::::::::::::::::::::::::::::::::::END DEBUG(111) OUTPUT FROM SUBROUTINE RSPLINE_PROC:::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEB_RSPLINE_PROC

      END SUBROUTINE RSPLINE_PROC


   END MODULE RIGID_ELEMENT_PROCESSING
