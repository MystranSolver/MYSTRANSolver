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

   MODULE USER_DEFINED_ELEMENTS

   USE DOF_LOOKUP_UTILS, ONLY :  TDOF_COL_NUM

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: KUSER1, USERIN

   CONTAINS

      SUBROUTINE KUSER1 ( OPT, WRITE_WARN )

! Calc's matrices for user supplied subroutine elements

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'KUSER1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not





! **********************************************************************************************************************************
      WRITE(ERR,1934) SUBR_NAME,TYPE
      WRITE(F06,1934) SUBR_NAME,TYPE
      FATAL_ERR = FATAL_ERR + 1
      CALL OUTA_HERE ( 'Y' )



      RETURN

! **********************************************************************************************************************************
 1934 FORMAT(' *ERROR  1934: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CODE NOT WRITTEN FOR ELEMENT TYPE = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE KUSER1


      SUBROUTINE USERIN ( INT_ELEM_ID, OPT, EMG_CALLING_SUBR, WRITE_WARN )

! Reads in matrices for a USERIN element from file IN4FIL specified in Exec Control (with statement IN4 i = in4file file name
! The matrices read are:

!   (1) IN4_MAT mass      matrix which is expanded to full ELDOF DOF size as elem mass ME
!   (2) IN4_MAT stiffness matrix which is expanded to full ELDOF DOF size as elem stiffness KE
!   (3) IN4_MAT load      matrix which is expanded to full ELDOF DOF size as elem stiffness PPE
!   (4) IN4_MAT RBM       Read or generate USERIN_RBM0

! These matrices are read as global DOF matrices. The global coord systems have to be defined on the CUSERIN/PUSERIN entries and
! they are defined relative to the basic coord system of that substructure. This way, there are no coord transformations needed
! before adding the element matrices to the overall global matrices for the integrated structure (i.e. OA structure made up of
! all substructures).

! The mass matrices are converted to same weight or mass units as other  elems by dividing here by WTMASS (if user has no
! WTMASS B.D. entry, PARAM WTMASS is 1.0). After the grid point weight generator is run (subr GPWG), the mass matrix will be
! converted back, same as all other mass, by multiplying by WTMASS. NOTE: WTMASS was checked to be sure it is > 0 when WTMASS B.D.
! entry was read

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN4, IN4_MSG, IN4FIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEDAT0_CUSERIN, MELDOF, NDOFG, NGRID, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI, TDOF_ROW_START
      USE PARAMS, ONLY                :  grdpnt, WTMASS
      USE INPUTT4_MATRICES, ONLY      :  IN4_COL_MAP, IN4_MAT
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, TR6_0
      USE MODEL_STUF, ONLY            :  AGRID, EDAT, EID, ELDOF, ELGP, EPNT, GRID_ID, INTL_PID, KE, ME, PPE, PUSERIN, TYPE,       &
                                         NUM_EMG_FATAL_ERRS
      USE MODEL_STUF, ONLY            :  USERIN_ACT_GRIDS, USERIN_ACT_COMPS, USERIN_CID0, USERIN_IN4_INDEX, USERIN_RBM0,           &
                                         USERIN_NUM_BDY_DOF, USERIN_NUM_ACT_GRDS, USERIN_NUM_SPOINTS,                              &
                                         USERIN_MASS_MAT_NAME, USERIN_LOAD_MAT_NAME, USERIN_RBM0_MAT_NAME, USERIN_STIF_MAT_NAME


      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_SET_CONSTRUCTION, ONLY  :  RDOF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE IN4_FILE_LIFECYCLE, ONLY    :  ALLOCATE_IN4_FILES, DEALLOCATE_IN4_FILES
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE OUTPUT4_FILE_IO, ONLY       :  READ_IN4_FULL_MAT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      EXTERNAL :: RB_DISP_MATRIX_PROC

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'USERIN'
      CHARACTER(LEN=*) , INTENT(IN)   :: EMG_CALLING_SUBR  ! Subr that called EMG which, in turn, called this subr
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! Contains 1 in each of the 6 pos'ns corresponding to a DOF from CGRID(I)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID for which
      INTEGER(LONG)                   :: ELDOF_INDEX       !
      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                            row number in array EDAT where data begins for this element.
      INTEGER(LONG)                   :: G_SET_COL         ! Col no. in array TDOF where the  G-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: I,J,L             ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Error count returned from subr READ_IN4_FULL_MAT
      INTEGER(LONG)                   :: IGRID             ! Internal grid number
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: NRC               ! Row/col size of stiff, mass matrices to be read from IN4 file
      INTEGER(LONG)                   :: IROW              ! Row number where TDOFI data begins for GRID_NUM
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr FILE_OPEN

                                                           ! Array that has USERIN grid num in col 1 and comp number in remaining 7
                                                           ! cols (1 col has all comps, others each indiv comp) for USERIN bdy DOF's
      INTEGER(LONG)                   :: USERIN_CID0_ICID  ! Internal coordinate system ID for USERIN_CID0


      REAL(DOUBLE)                    :: DX                ! X offset of USERIN elem CG from overall model basic sys origin
      REAL(DOUBLE)                    :: DY                ! Y offset of USERIN elem CG from overall model basic sys origin
      REAL(DOUBLE)                    :: DZ                ! Z offset of USERIN elem CG from overall model basic sys origin

                                                           ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM(USERIN_NUM_BDY_DOF,6)

      REAL(DOUBLE)                    :: RBM66(6,6)        ! 6x6 RB mass matrix read from IN4 file (if USERIN_RBM0_MAT_NAME exists)

                                                           ! Upper left RxR partition from MXX
      REAL(DOUBLE)                    :: MRRcb_FULL(USERIN_NUM_BDY_DOF,USERIN_NUM_BDY_DOF)

      REAL(DOUBLE)                    :: R66(6,6)          ! Rotates    basic coord sys of USERIN elem to basic sys of overall model
      REAL(DOUBLE)                    :: T66(6,6)          ! Translates basic coord sys of USERIN elem to basic sys of overall model

      REAL(DOUBLE)                    :: TN(3,3)           ! Transformation matrix from RCORD that transforms a vec from coord sys
!                                                            USERIN_CID0 to overall model basic coord system. Used to transform
!                                                            input RBM66 from substructure basic coords to OA model basic coords

                                                           ! Matrix to transform MRRcb to 6x6 RB mass rel to GRDPNT
      REAL(DOUBLE)                    :: TB6(USERIN_NUM_BDY_DOF,6)



! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

      EPNTK = EPNT(INT_ELEM_ID)

! Size of USERIN mass and stiff matrices on IN4FIL file are NRC x NRC

      NRC = USERIN_NUM_BDY_DOF + USERIN_NUM_SPOINTS

! Get rows from RBGLOBAL_GSET that are for the boundary DOF's and put them into matrix TB6

      DO I=1,USERIN_NUM_BDY_DOF                            ! First, init TB6
         DO J=1,6
            TB6(I,J) = ZERO
         ENDDO
      ENDDO

! In the grid point weight generator we need the initial 6x6 rigid body mass for this elem relative to the GPWG ref point (GRDPNT).
! Thus, TB6 needs to be a transformation matrix rel to GRDPNT (which is basic origin if GRDPNT = 0)

      CALL RB_DISP_MATRIX_PROC ( 'GRID', GRDPNT )
      CALL TDOF_COL_NUM ( 'G ',  G_SET_COL )
      K = 0
      DO I=1,USERIN_NUM_ACT_GRDS
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, USERIN_ACT_GRIDS(I), IGRID )
         IROW = TDOF( TDOF_ROW_START(IGRID), G_SET_COL )
         CALL RDOF ( USERIN_ACT_COMPS(I), CDOF )
         DO J=1,6
            IF (CDOF(J) /= '0') THEN
               K = K + 1
               DO L=1,6
                  TB6(K,L) = RBGLOBAL_GSET(IROW+J-1,L)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! Allocate memory for arrays IN4_MAT (matrix read from INPUTT4 file) and IN4_COL_MAP

      CALL ALLOCATE_IN4_FILES ( 'IN4_MAT'    , NRC, NRC, SUBR_NAME//' at 1st occurance' )
      CALL ALLOCATE_IN4_FILES ( 'IN4_COL_MAP', NRC, 1  , SUBR_NAME//' at 1st occurance' )

! Generate IN4_COL_MAP (used to expand cols of IN4_MAT to this elem's ELDOF size)

      DO I=1,NRC                                           ! Init IN4_COL_MAP
         IN4_COL_MAP(I) = 0
      ENDDO

      K = 0
      DO I=1,USERIN_NUM_ACT_GRDS                           ! Get IN4_COL_MAP for elem DOF's due to physical grid dislp components
         CALL RDOF ( USERIN_ACT_COMPS(I), CDOF )
         DO J=1,6
            ELDOF_INDEX = 6*(I-1) + J
            IF (CDOF(J) /= '0') THEN
                K = K + 1
                IN4_COL_MAP(K) = ELDOF_INDEX
            ENDIF
         ENDDO
      ENDDO
      DO I=1,USERIN_NUM_SPOINTS                           ! Get IN4_COL_MAP for elem DOF's due to SPOINT's
         K = K + 1
         ELDOF_INDEX = ELDOF_INDEX + 1
         IN4_COL_MAP(K) = ELDOF_INDEX
      ENDDO

! Get coord transformation matrix TN that transforms a vec from coord system USERIN_CID0 (dirs of basic sys in run that generates
! the USERIN matrices) to overall model basic system.

      CALL GET_TN_TRANSFORM_MAT ( USERIN_CID0_ICID, TN )

      IF (DEBUG(180) > 0) CALL DEB_USERIN ( 11 )
      IF (DEBUG(180) > 1) CALL DEB_USERIN ( 12 )

! Open the IN4 file containing the mass and stiff matrices for this elem

      CALL FILE_OPEN ( IN4, IN4FIL(USERIN_IN4_INDEX), OUNT, 'OLD', IN4_MSG, 'NEITHER', 'UNFORMATTED', 'READ', 'REWIND',            &
                       'Y', 'N' )

! If OPT(1) is 'Y', get the elem mass matrix from the IN4 file and expand it from boundary DOF size to ELDOF size (6 comps/grid).

      IF ((OPT(1) == 'Y') .AND. (USERIN_MASS_MAT_NAME(1:) /= ' ')) THEN

         CALL READ_IN4_FULL_MAT ( TYPE, EID, USERIN_MASS_MAT_NAME, NRC, NRC, IN4, IN4FIL(USERIN_IN4_INDEX), IN4_MAT,               &
                                  IERR, SUBR_NAME//' (for USERIN mass matrix)' )
         IF (IERR > 0) THEN
            WRITE(ERR,9999) SUBR_NAME, IERR
            WRITE(F06,9999) SUBR_NAME, IERR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DO I=1,NRC                                        ! Expand ME mass matrix to ELDOF size and divide by WTMASS to get
            DO J=1,NRC                                     ! same units as other elems (since input IN4_MAT is mass units)
               ME(IN4_COL_MAP(I),IN4_COL_MAP(J)) = IN4_MAT(I,J)/WTMASS
            ENDDO
         ENDDO

         IF (DEBUG(180) > 0) CALL DEB_USERIN ( 13 )

! If there is a RB mass matrix on the IN4 file read it as the 6x6 RB mass matrix for the elem. Otherwise calc the 6x6
! RB mass matrix from the boundary partition of the mass matrix (NASTRAN matrix MRRGN) and TB6, calc'd above.

         IF (USERIN_RBM0_MAT_NAME(1:) /= ' ') THEN
            CALL READ_IN4_FULL_MAT ( TYPE, EID, USERIN_RBM0_MAT_NAME, 6, 6, IN4, IN4FIL(USERIN_IN4_INDEX), RBM66,                  &
                                     IERR, SUBR_NAME//' (for USERIN_RBM0 matrix)' )
            DO I=1,6                                       ! Divide by WTMASS (default is 1.0) here to get USERIN elem mass same
               DO J=1,6                                    ! as other elems (i.e, either weight or mass). Before the elem mass is
                  RBM66(I,J) = RBM66(I,J)/WTMASS           ! used anywhere but in subr GPWG, it will be converted back
               ENDDO
            ENDDO

            IF (DEBUG(180) > 0) CALL DEB_USERIN ( 14 )

            IF (IERR > 0) THEN
               WRITE(ERR,9999) SUBR_NAME, IERR
               WRITE(F06,9999) SUBR_NAME, IERR
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            IF (USERIN_CID0 /= 0) THEN                     ! If USERIN basic not same as overall model basic, transform RBM66
               CALL TRANSFORM_USERIN_RBM0 ( USERIN_CID0_ICID, TN, RBM66, DX, DY, DZ, R66, T66 )
               IF (DEBUG(180) > 0) CALL DEB_USERIN ( 15 )
            ELSE
               DO I=1,6
                  DO J=1,6
                     USERIN_RBM0(I,J) = RBM66(I,J)
                  ENDDO
               ENDDO
            ENDIF

            IF (DEBUG(180) > 0) CALL DEB_USERIN ( 16 )

         ELSE                                              ! No RB mass matrix was input so we have to generate one
            DO I=1,USERIN_NUM_BDY_DOF                      ! First, partition MRRcb from input MXX mass matrix
               DO J=1,USERIN_NUM_BDY_DOF
                  MRRcb_FULL(I,J) = IN4_MAT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform MRRcb to basic coords
            CALL MATMULT_FFF   ( MRRcb_FULL  , TB6, USERIN_NUM_BDY_DOF, USERIN_NUM_BDY_DOF, 6, DUM )
            CALL MATMULT_FFF_T ( TB6, DUM         , USERIN_NUM_BDY_DOF, 6                 , 6, USERIN_RBM0 )
            DO I=1,6                                       ! Divide by WTMASS to get units same as input mass
               DO J=1,6
                  USERIN_RBM0(I,J) = USERIN_RBM0(I,J)/WTMASS
               ENDDO
            ENDDO

            IF (DEBUG(180) > 0) CALL DEB_USERIN ( 17 )

         ENDIF

      ENDIF

! If OPT(4) is 'Y', get the elem stiff matrix from the IN4 file and expand it from boundary DOF size to ELDOF size (6 comps/grid)

      IF ((OPT(4) == 'Y') .AND. (USERIN_STIF_MAT_NAME(1:) /= ' ')) THEN

         CALL READ_IN4_FULL_MAT ( TYPE, EID, USERIN_STIF_MAT_NAME, NRC, NRC, IN4, IN4FIL(USERIN_IN4_INDEX), IN4_MAT,               &
                                  IERR, SUBR_NAME//' (for USERIN stiffness matrix)' )
         IF (IERR > 0) THEN
            WRITE(ERR,9999) SUBR_NAME, IERR
            WRITE(F06,9999) SUBR_NAME, IERR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DO I=1,NRC                                        ! Expand KE stiffness matrix to ELDOF size
            DO J=1,NRC
               KE(IN4_COL_MAP(I),IN4_COL_MAP(J)) = IN4_MAT(I,J)
            ENDDO
         ENDDO

         IF (DEBUG(180) > 0) CALL DEB_USERIN ( 18 )

      ENDIF

! If OPT(4) is 'Y', get the elem load matrix from the IN4 file and expand it from boundary DOF size to ELDOF size (6 comps/grid)
! (we use OPT(4) since there is no OPT for a load matrix like this one and the only time we will need it is when stiff is needed)

      IF ((OPT(5) == 'Y') .AND. (USERIN_LOAD_MAT_NAME(1:) /= ' ')) THEN

         CALL DEALLOCATE_IN4_FILES ( 'IN4_MAT' )
         CALL ALLOCATE_IN4_FILES ( 'IN4_MAT', NRC, NSUB, SUBR_NAME//' at 2nd occurance' )

         CALL READ_IN4_FULL_MAT ( TYPE, EID, USERIN_LOAD_MAT_NAME, NRC, NSUB, IN4, IN4FIL(USERIN_IN4_INDEX), IN4_MAT,              &
                                  IERR, SUBR_NAME//' (for USERIN load matrix)' )
         IF (IERR > 0) THEN
            WRITE(ERR,9999) SUBR_NAME, IERR
            WRITE(F06,9999) SUBR_NAME, IERR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DO I=1,NRC                                        ! Expand PPE load matrix to ELDOF size
            DO J=1,NSUB
               PPE(IN4_COL_MAP(I),J) = IN4_MAT(I,J)
            ENDDO
         ENDDO

         IF (DEBUG(180) > 0) CALL DEB_USERIN ( 19 )

         CALL DEALLOCATE_IN4_FILES ( 'IN4_MAT' )

      ENDIF

! Deallocate USERIN arrays used internally here.

      CALL DEALLOCATE_IN4_FILES ( 'IN4_COL_MAP' )
      CALL DEALLOCATE_IN4_FILES ( 'IN4_MAT' )

      CALL FILE_CLOSE ( IN4, IN4FIL(USERIN_IN4_INDEX), 'KEEP' )

      IF (DEBUG(180) > 0) CALL DEB_USERIN ( 99 )



      RETURN

! **********************************************************************************************************************************
 1952 FORMAT(' *ERROR  1952: PROGRAMMING ERROR IN SUBROUTINE ',A,' CALLED BY SUBR EMG WHICH WAS CALLED BY SUBR ',A                 &
                    ,/,14X,' INCORRECT OPT = ',5A1,' INPUT. ONLY OPT(1) AND OPT(4) ARE PROGRAMMED')

 9999 FORMAT(' PROCESSING STOPPED IN SUBR ',A,' DUE TO ABOVE ',I8,' ERROR(S)')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_TN_TRANSFORM_MAT ( USERIN_CID0_ICID, TN )

! Searches array CORD to find row where actual coord sys USERIN_CID0 exists and sets the internal coord sys ID (USERIN_CID0_ICID)
! equal to that row number in CORD. Then, get the coord transformation matrix, TN, from that row in array RCORD

      USE PENTIUM_II_KIND
      USE SCONTR, ONLY                :  FATAL_ERR, NCORD
      USE IOUNT1, ONLY                :  ERR, F06
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  CORD, NUM_EMG_FATAL_ERRS, RCORD, USERIN_CID0

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: CORD_FND          ! Indicator of whether we found a coord sys in array CORD

      INTEGER(LONG), INTENT(OUT)      :: USERIN_CID0_ICID  ! Internal coordinate system ID for USERIN_CID0
      INTEGER(LONG)                   :: J                 ! DO loop index

      REAL(DOUBLE) , INTENT(OUT)      :: TN(3,3)           ! Transformation matrix from RCORD

! **********************************************************************************************************************************

      USERIN_CID0_ICID = 0
      IF (USERIN_CID0 /= 0) THEN
         CORD_FND = 'N'
j_do12:     DO J=1,NCORD
            IF (USERIN_CID0 == CORD(J,2)) THEN
               CORD_FND = 'Y'
               USERIN_CID0_ICID = J
               EXIT j_do12
            ENDIF
         ENDDO j_do12
         IF (CORD_FND == 'N') THEN                         ! Coord sys ID on USERIN elem undefined
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1822) 'COORD SYSTEM ', USERIN_CID0, TYPE, EID
            WRITE(F06,1822) 'COORD SYSTEM ', USERIN_CID0, TYPE, EID
         ENDIF
      ENDIF
                                                           ! TN transforms vec from coord sys USERIN_CID0 to overall model basic sys
      IF (USERIN_CID0_ICID /= 0) THEN
         TN(1,1) = RCORD(USERIN_CID0_ICID, 4)   ;   TN(1,2) = RCORD(USERIN_CID0_ICID, 5)   ;   TN(1,3) = RCORD(USERIN_CID0_ICID, 6)
         TN(2,1) = RCORD(USERIN_CID0_ICID, 7)   ;   TN(2,2) = RCORD(USERIN_CID0_ICID, 8)   ;   TN(2,3) = RCORD(USERIN_CID0_ICID, 9)
         TN(3,1) = RCORD(USERIN_CID0_ICID,10)   ;   TN(3,2) = RCORD(USERIN_CID0_ICID,11)   ;   TN(3,3) = RCORD(USERIN_CID0_ICID,12)
      ELSE                                                 ! USERIN_CID0 is basic system of overall model
         TN(1,1) = ONE                          ;   TN(1,2) = ZERO                         ;   TN(1,3) = ZERO
         TN(2,1) = ZERO                         ;   TN(2,2) = ONE                          ;   TN(2,3) = ZERO
         TN(3,1) = ZERO                         ;   TN(3,2) = ZERO                         ;   TN(3,3) = ONE
      ENDIF

 ! *********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

! **********************************************************************************************************************************

      END SUBROUTINE GET_TN_TRANSFORM_MAT

! ##################################################################################################################################

      SUBROUTINE TRANSFORM_USERIN_RBM0 ( USERIN_CID0_ICID, TN, RBM66, DX, DY, DZ, R66, T66 )

! Transforms IN4 input matrix RBM66 from USERIN elem basic coords to overall model basic coords (incl translating from the
! USERIN basic coord origin to the overall model basic origin, The transformation is accomplished in 2 steps:
!  (1) Rotate RBM66 from USERIN basic coord system to be parallel with overall model basic coord system.
!  (2) Translate from USERIN basic coord sys origin to overall model basic coord sys origin. Final result is USERIN_RBM0

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      use iount1, only                :  err, f06
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  RCORD, USERIN_RBM0

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: USERIN_CID0_ICID  ! Internal coordinate system ID for USERIN_CID0
      INTEGER(LONG)                   :: I,J             ! DO loop indices

      REAL(DOUBLE), INTENT(IN)        :: TN(3,3)           ! Transformation matrix from RCORD
      REAL(DOUBLE), INTENT(IN)        :: RBM66(6,6)        ! 6x6 RB mass matrix read from IN4 file (if USERIN_RBM0_MAT_NAME exists)
                                                           ! Upper left RxR partition from MXX
      REAL(DOUBLE), INTENT(OUT)       :: DX                ! X offset of USERIN elem CG from overall model basic sys origin
      REAL(DOUBLE), INTENT(OUT)       :: DY                ! Y offset of USERIN elem CG from overall model basic sys origin
      REAL(DOUBLE), INTENT(OUT)       :: DZ                ! Z offset of USERIN elem CG from overall model basic sys origin
      REAL(DOUBLE), INTENT(OUT)       :: R66(6,6)          ! Rotates    basic coord sys of USERIN elem to basic sys of overall model
      REAL(DOUBLE), INTENT(OUT)       :: T66(6,6)          ! Translates basic coord sys of USERIN elem to basic sys of overall model
      REAL(DOUBLE)                    :: DUM1(6,6)         ! Intermediate matrix in a calc

! **********************************************************************************************************************************
! (1) Rotate RBM66 from USERIN basic coord system to be parallel with overall model basic coord system. Call result temp USERIN_RBM0
!     Rotational transformation matrix R66 has 2 copies of the transpose of 3x3 matrix TN on its diagonal
!     TN transforms a vec from coord system USERIN_CID0 (dirs of basic sys in run that generates the USERIN matrices)
!     to overall model basic system.

      DO I=1,6                                         ! Initialize 6x6 transformation matrix for USERIN RB mass
         DO J=1,6
            R66(I,J) = ZERO
         ENDDO
      ENDDO
      DO I=1,3                                         ! Rows/cols 1-3 of R66 are the transpose of the 3 rows/cols from TN
         DO J=1,3
            R66(I,J) = TN(J,I)
         ENDDO
      ENDDO
      DO I=4,6                                         ! Rows/cols 4-6 of R66 are also the transpose of the 3 rows/cols from TN
         DO J=4,6
            R66(I,J) = TN(J-3,I-3)
         ENDDO
      ENDDO

! RBM66, that was read in via IN4FIL, is rotated to overall sys basic coords as temp USERIN_RBM0 = T66'*RBM66*T66:

      CALL MATMULT_FFF   ( RBM66, R66 , 6, 6, 6, DUM1 )
      CALL MATMULT_FFF_T ( R66      , DUM1, 6, 6, 6, USERIN_RBM0 )

! **********************************************************************************************************************************
! (2) Translate temp USERIN_RBM0 from USERIN basic coord sys origin to overall model basic coord sys origin.
!     This is done with coord transformation matrix T66 which is of the form:

!                                            | I   D |
!                                      T66 = |       |
!                                            | O   I |

! where I is a 3x3 identity matrix, O a 3x3 null matrix and D is the 3x3 matrix that has offsets of USERIN coord sys ICID (basic
! coord sys of USERIN elem) relative to basic coord sys origin of the overall model:

!                                            |  0    DZ -DY |
!                                        D = | -DZ   0   DX | , DX = RCORD(ICID,1), DY = RCORD(ICID,2), DZ = RCORD(ICID,3)
!                                            |  DY  -DX  0  |

      DX = RCORD(USERIN_CID0_ICID,1)                       ! Get DX, DY, DZ offsets of USERIN internal coord sys ICID (USERIN basic
      DY = RCORD(USERIN_CID0_ICID,2)                       ! coord sys definition) relative to overall model basic coord sys
      DZ = RCORD(USERIN_CID0_ICID,3)

      DO I=1,6                                             ! Init T66
         DO J=1,6
            T66(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,3                                             ! Upper left 3x3, and lower right 3x3, matrices I are identity matrices
         T66(I,I)     = ONE
         T66(I+3,I+3) = ONE
      ENDDO
                                                           ! Upper right 3x3 matrix in T66 is D
      T66(1,4) =  ZERO   ;   T66(1,5) =  DZ     ;   T66(1,6) = -DY
      T66(2,4) = -DZ     ;   T66(2,5) =  ZERO   ;   T66(2,6) =  DX
      T66(3,4) =  DY     ;   T66(3,5) = -DX     ;   T66(3,6) =  ZERO
                                                           ! Transform temp USERIN_RBM0 to final USERIN_RBM0
      CALL MATMULT_FFF   ( USERIN_RBM0, T66, 6, 6, 6, DUM1 )
      CALL MATMULT_FFF_T ( T66, DUM1, 6, 6, 6, USERIN_RBM0 )

      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

! **********************************************************************************************************************************

      END SUBROUTINE TRANSFORM_USERIN_RBM0

! ##################################################################################################################################

      SUBROUTINE DEB_USERIN ( WHAT )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT              ! What section below to write out

! **********************************************************************************************************************************
      IF      (WHAT == 11) THEN

         WRITE(F06,98720)
         WRITE(F06,'(49X,A,I8)') 'Data for CUSERIN   element ', EID
         WRITE(F06,'(49X,A   )') '***********************************'
         WRITE(F06,*) ' General data:'
         WRITE(F06,*) ' ------------'

         WRITE(F06,98721) ELGP, ELDOF, OPT, EPNTK, PUSERIN(INTL_PID,1), INTL_PID, USERIN_IN4_INDEX, USERIN_NUM_ACT_GRDS,           &
                          USERIN_NUM_SPOINTS, USERIN_CID0, NRC, USERIN_STIF_MAT_NAME, USERIN_MASS_MAT_NAME, USERIN_LOAD_MAT_NAME

         IF (USERIN_RBM0_MAT_NAME(1:) /= ' ') THEN
            WRITE(F06,98722) USERIN_RBM0_MAT_NAME
         ELSE
            WRITE(F06,98723)
         ENDIF

         WRITE(F06,*) ' GRID''s and SPOINT''s that this element connects to:'
         WRITE(F06,*) ' -------------------------------------------------'
         WRITE(F06,98731) (AGRID(I),I=1,ELGP)
         WRITE(F06,*)

         WRITE(F06,*) ' Matrix TB6: 6 rows from RBGLOBAL_GSET that are for the boundary DOF''s:'
         WRITE(F06,*) ' ---------------------------------------------------------------------'
         WRITE(F06,98741) (J,J=1,6)
         DO I=1,USERIN_NUM_BDY_DOF
            WRITE(F06,98742) I, (TB6(I,J),J=1,6)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 12) THEN

         WRITE(F06,*) ' Merging vector IN4_COL_MAP (used to expand cols of INPUTT4 matrices for this element to ELDOF size):'
         WRITE(F06,*) ' ---------------------------------------------------------------------------------------------------'
         DO I=1,NRC
            WRITE(F06,98711) I, IN4_COL_MAP(I)
         ENDDO
         WRITE(F06,*)

         WRITE(F06,*) ' Matrix TN transforms a vector from coord sys USERIN_CID0 = ',userin_cid0,' to overall model basic sys:'
         WRITE(F06,*) ' -----------------------------------------------------------------------------------------------------'
         WRITE(F06,98741) (J,J=1,3)
         DO J=1,3
            WRITE(F06,98742) J,(TN(J,K),K=1,3)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 13) THEN

         WRITE(F06,*) ' Element mass matrix ME in full ELDOF size:'
         WRITE(F06,*) ' -----------------------------------------'
         WRITE(F06,98741) (j,j=1,eldof)
         DO I=1,USERIN_NUM_ACT_GRDS
            DO J=1,6
               K = 6*(I-1) + J
               WRITE(F06,98742) K,(ME(K,L),L=1,6*USERIN_NUM_ACT_GRDS)
            ENDDO
            WRITE(F06,*)
         ENDDO
         DO J=6*USERIN_NUM_ACT_GRDS+1,ELDOF
            WRITE(F06,98742) J,(ME(J,K),K=1,ELDOF)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 14) THEN

         WRITE(F06,*) ' Rigid body mass matrix read from IN4 file:'
         WRITE(F06,*) ' -----------------------------------------'
         WRITE(F06,98741) (j,j=1,6)
         DO J=1,6
            WRITE(F06,98742) J,(RBM66(J,K),K=1,6)
            IF (J == 3) WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 15) THEN

         WRITE(F06,*) ' Matrix R66: Rotates    basic coord sys of USERIN elem to basic sys of overall model:'
         WRITE(F06,*) ' -----------------------------------------------------------------------------------'
         WRITE(F06,98741) (J,J=1,6)
         DO I=1,6
            WRITE(F06,98742) I,(R66(I,J),J=1,6)
            IF (I == 3) WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

         WRITE(F06,'(a,1es14.6,a)') '  DX = ',dx,' (X offset of USERIN elem basic origin expressed in overall model basic coords)'
         WRITE(F06,'(a,1es14.6,a)') '  DY = ',dy,' (Y offset of USERIN elem basic origin expressed in overall model basic coords)'
         WRITE(F06,'(a,1es14.6,a)') '  DZ = ',dz,' (Z offset of USERIN elem basic origin expressed in overall model basic coords)'
         WRITE(F06,*)

         WRITE(F06,*) ' Matrix T66: Translates basic coord sys of USERIN elem to basic sys of overall model:'
         WRITE(F06,*) ' -----------------------------------------------------------------------------------'
         WRITE(F06,98741) (J,J=1,6)
         DO I=1,6
            WRITE(F06,98742) I,(T66(I,J),J=1,6)
            IF (I == 3) WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 16) THEN

         WRITE(F06,*) ' USERIN_RBM0: 6x6 RB mass matrix divided by WTMASS (from RB matrix input and transformed):'
         WRITE(F06,*) ' ----------------------------------------------------------------------------------------'
         WRITE(F06,98741) (j,j=1,6)
         DO J=1,6
            WRITE(F06,98742) J,(USERIN_RBM0(J,K),K=1,6)
            IF (J == 3) WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 17) THEN

         WRITE(F06,*) ' USERIN_RBM0: 6x6 RB mass matrix divided by WTMASS (generated internally):'
         WRITE(F06,*) ' ------------------------------------------------------------------------'
         WRITE(F06,98741) (j,j=1,6)
         DO J=1,6
            WRITE(F06,98742) J,(USERIN_RBM0(J,K),K=1,6)
            IF (J == 3) WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 18) THEN

         WRITE(F06,*) ' Element stiffness matrix KE in full ELDOF size:'
         WRITE(F06,*) ' ----------------------------------------------'
         WRITE(F06,98741) (j,j=1,eldof)
         DO I=1,USERIN_NUM_ACT_GRDS
            DO J=1,6
               K = 6*(I-1) + J
               WRITE(F06,98742) K,(KE(K,L),L=1,6*USERIN_NUM_ACT_GRDS)
            ENDDO
            WRITE(F06,*)
         ENDDO
         DO J=6*USERIN_NUM_ACT_GRDS+1,ELDOF
            WRITE(F06,98742) J,(KE(J,K),K=1,ELDOF)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 19) THEN

         WRITE(F06,*) ' Element load matrix PPE in full ELDOF size:'
         WRITE(F06,*) ' ------------------------------------------'
         WRITE(F06,98741) (j,j=1,eldof)
         K = 0
         DO I=1,USERIN_NUM_ACT_GRDS
            DO J=1,NSUB
               K = K + 1
               WRITE(F06,98742) K,(KE(K,L),L=1,NSUB)
            ENDDO
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 99) THEN

         WRITE(F06,98799)

      ENDIF

      RETURN

! ********************************************************************************************************************************
86429 FORMAT(6(1ES14.6))

98706 FORMAT(' EID, GRID, COMP, RBGLOBAL_GSET(1-6) = ',3I8,6(1ES14.6))

98707 FORMAT(6(1ES14.6))

98711 FORMAT('         I, IN4_COL_MAP(I)  = ',3I8)

98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::::::::START DEBUG(180) OUTPUT FROM SUBROUTINE USERIN::::::::::::::::::::::::::',&
             ':::::::::::::::::',/)

98721 FORMAT('      ( 1) Number of elem grids plus scalar points = ',i8,/,                                                         &
             '      ( 2) Number of element DOF''s                 = ',i8,/,                                                        &
             '      ( 3) OPT(1-6)                                = ',2x,6a,/,                                                      &
             '      ( 4) Data begins in EDAT at EPNTK            = ',i8,/,                                                         &
             '      ( 5) Actual property ID                      = ',i8,/,                                                         &
             '      ( 6) Internal property ID                    = ',i8,/,                                                         &
             '      ( 7) IN4 index number for matrices           = ',i8,/,                                                         &
             '      ( 8) USERIN_NUM_ACT_GRDS (number grids)      = ',i8,/,                                                         &
             '      ( 9) USERIN_NUM_SPOINTS  (number SPOINTs)    = ',i8,/,                                                         &
             '      (10) USERIN_CID0 (basic coord sys ID)        = ',i8,/,                                                         &
             '      (11) Row/col size of matrices to read        = ',i8,' (size of compacted IN4 matrices, not elem ELDOF)',/,     &
             '      (12) Name of stiffness matrix                = "',a,'"',/,                                                     &
             '      (13) Name of mass      matrix                = "',a,'"',/,                                                     &
             '      (14) Name of load      matrix                = "',a,'"')

98722 FORMAT('      (15) Name of RB mass   matrix (USERIN_RBM0)  = "',a,'"',/)

98723 FORMAT('      (15) Name of RB mass   matrix (USERIN_RBM0)  = (none input so one will be generated internally)',/)

98731 FORMAT(10I10)

98741 FORMAT(4X,32767(13X,I3))

98742 FORMAT(I10,32767(2X,1ES14.6))

98799 FORMAT(' :::::::::::::::::::::::::::::::::::::::::::END DEBUG(180) OUTPUT FROM SUBROUTINE USERIN:::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEB_USERIN

      END SUBROUTINE USERIN

   END MODULE USER_DEFINED_ELEMENTS
