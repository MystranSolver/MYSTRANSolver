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

   MODULE RFORCE_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: RFORCE_PROC

   CONTAINS

      SUBROUTINE RFORCE_PROC

! RFORCE load processor. Forces on grids for an RFORCE are:

!           Fi = -Mi*[W x (W x (Ri - Ra)) - A x (Ri - Ra)]

! where x means a vector cross product and:

!           Mi = 6x6 mass matrix at grid i
!           W  = angular velocity of the model about grid i (variable ANG_VEL)
!           Ri = radius from basic coord system origin to grid i (variable RI)
!           Ra = radius from basic coord system origin to the reference point for angular velocity/accel (variable RA)
!           A  = angular acceleration (variable ANG_ACC)

! The input B.D. RFORCE load data is transferred to system force data in the SYS_LOAD array for DOF's in the G set.

! File LINK1U was written when RFORCE B.D entries were read, with one record for each such card.
! There are NRFORCE total number of records written to file LINK1U with each record containing:

!               SETID         = Load set ID
!               ACID_L        = Local coord sys ID that RFORCE load is given in
!               RFORCE_GRID   = ID of grid that rotational (components 4, 5, 6) RFORCE velocity/accels are about
!               SCALEF_AV     = Scale factor for angular velocity in revolutions per unit time.
!               SCALEF_AA     = Scale factor for angular accel in revolutions per unit time squared.
!               VEC(1-3)      = 3 components of the vector for the velocity and/or accel

! The process in creating array SYS_LOAD from this information is as follows:

!  (1) For each record (1 to NRFORCE) in file LINK1U:

!      (a) Read a record

!      (b) Transform coords from local (on RFORCE card) to basic. Transformation to global is done later (see 2a(iii))

!      (d) Store records in global coords.

!  (2) For each subcase (1 to NSUB):

!      (a) Generate LSID, RSID tables of load set ID's/scale factors for load sets for this subcase:

!          (  i) LLOADC is the max number of pairs of scale factors/load set ID's over all LOAD Bulk Data cards
!                including the pair defined by the set ID and overall scale factor on the LOAD Bulk Data card.
!                LSID and RSID are dimensioned 1 to LLOADC and:
!                   LSID(1) is always the load set ID requested in Case Control for this subcase.
!                   RSID(1) is always 1.0

!          ( ii) If the load set requested in Case Control is not a set ID from a LOAD Bulk data card then the
!                load must be on a separate RFORCE card in which case LSID(1) and RSID(1) are all that is needed
!                to define the load set contribution due to the RFORCE card for this subcase.

!          (iii) If the load set requested in Case Control is a set ID from a LOAD Bulk data card then the ramainder
!                (K = 2,LLOADC) of entries into LSID and RSID will be the pairs of load set ID's/scale factors from
!                that LOAD Bulk Data card (with RSID also multiplied by the overall scale factor on the LOAD Bulk data
!                card. The load set ID's are in array LOAD_SIDS(i,j) created when LOAD Bulk Data cards were read.
!                The scale factors are in array LOAD_FACS(i,j) also created when LOAD Bulk Data cards were read.
!                Note, there may not be as many as LLOADC pairs of set ID's/scale factors on a given LOAD Bulk Data
!                card since LLOADC is the max, from all LOAD Bulk Data cards, of pairs.
!                Thus, the entries in LSID from the last entry (for a given LOAD card) to LLOADC will be zero (LSID
!                was initialized to zero). This fact is used in a DO loop to EXIT when LSID(K) = 0

!      (b) For each store record (1 to NRFORCE)

!          ( ii) Scan LSID and RSID to get the scale factor (SCALE) for the ACCEL_RB components in SETID, if this
!                RFORCE's set ID is in LSID. When found, reset RFORCE vector components to ACCEL_RB = SCALE*ACCEL_RB.
!                At this point ACCEL_RB has the correct magnitudes and is in basic coordinates.

!          (iii) For a grid point, determine if a transformation of ACCEL_RB to global is needed and transform it if so.

!          ( iv) Calculate accel at a grid (ACCEL_I) based on rigid body motion of ACCEL_RB

!          (  v) Get the 6 x 6 mass matrix for a grid point times ACCEL_I to get RFORCE forces at this grid

!          ( vi) Load the RFORCE forces into the SYS_LOAD (systems load) array

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1U, LINK1U, L1U_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LLOADC, NCORD, NRFORCE, NGRID, NLOAD, NSUB, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, PI
      USE PARAMS, ONLY                :  SUPWARN
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  CORD, GRID, GRID_ID, LOAD_FACS, LOAD_SIDS, RCORD, RGRID, SYS_LOAD, SUBLOD

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  READERR
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, GEN_T0L
      USE LINK1_LOAD_SUPPORT, ONLY    :  GET_GRID_6X6_MASS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RFORCE_PROC'
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for
      CHARACTER( 1*BYTE)              :: GRID_MGG_FND      ! Indicator on whether a mass matrix was found in MGG for a given grid
      CHARACTER( 8*BYTE)              :: NAME              ! Name for output error purposes

      INTEGER(LONG)                   :: ACID_L            ! Actual local  coord sys ID on FORCE or MOMENT card
      INTEGER(LONG)                   :: ACID_G            ! Actual global coord sys ID for an actual grid
      INTEGER(LONG)                   :: CID_ERR   = 0     ! Count of coord systems undefined
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of grids undefined
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number for a grid
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: ISUB              ! Subcase index
      INTEGER(LONG)                   :: IRFORCE           ! RFORCE index
      INTEGER(LONG)                   :: ICID              ! Internal coordinate system ID for ACID_L or ACID_G
      INTEGER(LONG)                   :: IERRT             ! Total number of errors found
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: LSID(LLOADC+1)    ! Array of load SID's, for RFORCE cards, needed for one S/C
      INTEGER(LONG)                   :: NCOLA             ! No. cols in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOLB             ! No. cols in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NROWA             ! No. rows in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: READ_ERR  = 0     ! Cum. count of errors as we read, and check cards from file LINK1K
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: RFORCE_GRD        ! ID of grid that rotational RFOECE vel/accels are about
      INTEGER(LONG)                   :: RFORCE_GRD_ROW_NUM! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: ROW_NUM           ! Row no. in array TDOF corresponding to GDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: SETID             ! Load set ID read from record in file LINK1U
      INTEGER(LONG)                   :: SETIDS(NRFORCE)
      INTEGER(LONG)                   :: RFORCE_GRDS(NRFORCE)

      REAL(DOUBLE)                    :: ACCEL_I(6)        ! 6 components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_I_T1(3)     ! 3 transl components of accel due to RFORCE at a grid in basic  coords
      REAL(DOUBLE)                    :: ACCEL_I_T2(3)     ! 3 transl components of accel due to RFORCE at a grid in global coords
      REAL(DOUBLE)                    :: ACCEL_I_R1(3)     ! 3 rotat  components of accel due to RFORCE at a grid in basic  coords
      REAL(DOUBLE)                    :: ACCEL_I_R2(3)     ! 3 rotat  components of accel due to RFORCE at a grid in global coords
      REAL(DOUBLE)                    :: ANG_ACC(3)        ! Angular acceleration in units of radian per unit time squared.
      REAL(DOUBLE)                    :: ANG_VEL(3)        ! Angular velocity in units of radian per unit time.
      REAL(DOUBLE)                    :: DRI(3)            ! Components of the vector formed by RI - RA
      REAL(DOUBLE)                    :: FORCE_I(6)        ! 6 forces at a grid due to the RFORCE loading
      REAL(DOUBLE)                    :: GRID_MGG(6,6)     ! 6 X 6 mass matrix for one grid point
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: RSID(LLOADC+1)    ! Array of RFORCE magnitudes (for LSID set ID's) needed for one S/C
      REAL(DOUBLE)                    :: RA(3)             ! Vector components, in basic coords, of GID
      REAL(DOUBLE)                    :: RI(3)             ! Vector components, in basic coords, of grid i
      REAL(DOUBLE)                    :: SCALE             ! Scale factor for a load (on a LOAD Bulk Data entry)
      REAL(DOUBLE)                    :: T12(3,3)          ! Coord transformation matrix
      REAL(DOUBLE)                    :: VEC_LOCAL(3)      ! 3 components of RFORCE vector at RFORCE_GRD in local coords, ACID_L
      REAL(DOUBLE)                    :: VEC_BASIC(3)      ! 3 components of RFORCE vector at RFORCE_GRD in basic coords, ACID_0
      REAL(DOUBLE)                    :: DUM1(3)           ! Intermediate vector in cross product
      REAL(DOUBLE)                    :: DUM2(3)           ! Intermediate vector in cross product
      REAL(DOUBLE)                    :: DUM3(3)           ! Intermediate vector in cross product
      REAL(DOUBLE)                    :: SCALEF_AVS(NRFORCE)
      REAL(DOUBLE)                    :: SCALEF_AAS(NRFORCE)
      REAL(DOUBLE)                    :: VECS_BASIC(3,NRFORCE)
      


! **********************************************************************************************************************************
      NAME = 'RFORCE  '

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! **********************************************************************************************************************************
! (1) Read record from L1U and transform coords for RFORCE VEC from local sys (on RFORCE bulk data card) to basic.

i_do1:DO IRFORCE=1,NRFORCE
                                                           ! Read a record from L1U
         READ(L1U,IOSTAT=IOCHK) SETIDS(IRFORCE),                                                                                             &
                                ACID_L,                                                                                            &
                                RFORCE_GRDS(IRFORCE),                                                                              &
                                SCALEF_AVS(IRFORCE),                                                                               &
                                SCALEF_AAS(IRFORCE),                                                                               &
                                (VEC_LOCAL(J),J=1,3)

         IF (IOCHK /= 0) THEN
            REC_NO = IRFORCE
            CALL READERR ( IOCHK, LINK1U, L1U_MSG, REC_NO, OUNT )
            READ_ERR = READ_ERR + 1                        ! Increment READ_ERR and go back to read another RFORCE card
            CYCLE i_do1
         ENDIF
                                                           ! The local system that RFORCE is defined in is ACID_L.
         IF (ACID_L /= 0) THEN                             ! ACID_L is not basic, so find it and transform coords
            FOUND = 'N'
j_do12:     DO J=1,NCORD
               IF (CORD(J,2) == ACID_L) THEN
                  FOUND = 'Y'
                  ICID = J                                 ! ICID is the internal coord ID corresponding to ACID_L
                  EXIT j_do12
               ENDIF
            ENDDO j_do12
            IF (FOUND == 'N') THEN
               WRITE(ERR,1822) 'COORD SYSTEM ', ACID_L, NAME, SETIDS(IRFORCE)
               WRITE(F06,1822) 'COORD SYSTEM ', ACID_L, NAME, SETIDS(IRFORCE)
               CID_ERR = CID_ERR + 1                       ! Increment READ_ERR and go back to read another RFORCE card
               FATAL_ERR = FATAL_ERR + 1
               CYCLE i_do1
            ENDIF

            DO J=1,3                                       ! Get coord transf matrix (don't need GEN_T0L since ICID is rect.)
               DO K=1,3
                  T12(J,K) = RCORD(ICID, 3 + 3*(J-1) + K)
               ENDDO
            ENDDO
                                                           ! Transform coordinates
            CALL MATMULT_FFF ( T12, VEC_LOCAL, 3, 3, 1, VEC_BASIC )
         ELSE                                              ! No transformation needed since ACID_L is basic
            VEC_BASIC = VEC_LOCAL
         ENDIF

                                                           ! Store data. RFORCE vec now in basic coords
         VECS_BASIC(:,IRFORCE) = VEC_BASIC

      ENDDO i_do1

      IF ((CID_ERR > 0) .OR. (READ_ERR > 0)) THEN
         IERRT = CID_ERR + READ_ERR
         WRITE(ERR,1599) SUBR_NAME,IERRT
         WRITE(F06,1599) SUBR_NAME,IERRT
         CALL OUTA_HERE ( 'Y' )                                    ! Errors from reading RFORCE data, so quit
      ENDIF

! **********************************************************************************************************************************
! Now process RFORCE loads into SYS_LOAD

                                                          ! Initialize LSID, RSID arrays
      LSID(1:LLOADC) = 0
      RSID(1:LLOADC) = ZERO

      IERRT = 0
i_do2:DO ISUB = 1,NSUB                                     ! Loop through the S/C's

         IF (SUBLOD(ISUB,1) == 0) THEN                     ! If no load for this S/C, CYCLE
            CYCLE i_do2
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(ISUB,1)                          ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the RFORCE card in file LINK1U that matches SUBLOD(I,1)
         DO J = 1,NLOAD                                    ! Then, the actual mag. will come from RSID(1) & the ACCEL_RB components

            IF (LSID(1) == LOAD_SIDS(J,1)) THEN            ! The load requested in CC for this S/C is the j-th LOAD BD card
k_do_211:      DO K = 2,LLOADC                             ! Get load sets defined on this LOAD BD card and put into LSID (if any)
                  IF (LOAD_SIDS(J,K) /= 0) THEN
                     NSID = K
                     LSID(K) = LOAD_SIDS(J,K)
                     RSID(K) = LOAD_FACS(J,1)*LOAD_FACS(J,K)
                  ELSE
                     CYCLE k_do_211                        ! If a LSID field left blank on LOAD BD card, CYCLE
                  ENDIF
               ENDDO k_do_211
            ENDIF

         ENDDO

j_do_22: DO IRFORCE = 1,NRFORCE                            ! Process RFORCE card info that is now in basic coords

            SETID = SETIDS(IRFORCE)
            RFORCE_GRD = RFORCE_GRDS(IRFORCE)

                                                           ! Find the location of the axis from the grid point ID.
            RA = ZERO
            IF (RFORCE_GRD > 0) THEN
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, RFORCE_GRD, RFORCE_GRD_ROW_NUM )
               IF (RFORCE_GRD_ROW_NUM == -1) THEN
                  WRITE(ERR,1822) 'GRID ', RFORCE_GRD, NAME, SETID
                  WRITE(F06,1822) 'GRID ', RFORCE_GRD, NAME, SETID
                  GID_ERR = GID_ERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE j_do_22                            ! Don't apply this faulty RFORCE. Proceed to the next one.
               ELSE
                  RA = RGRID(RFORCE_GRD_ROW_NUM,1:3)
               ENDIF
            ENDIF

            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
            DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                                                           ! Ang accel and vel of model
                  ANG_ACC = 2*PI*SCALE*SCALEF_AAS(IRFORCE)*VECS_BASIC(:,IRFORCE)
                  ANG_VEL = 2*PI*SCALE*SCALEF_AVS(IRFORCE)*VECS_BASIC(:,IRFORCE)
                  EXIT
               ENDIF
            ENDDO

            IF (FOUND /= 'Y') THEN                         ! This RFORCE set ID isn't called for in this S/C, so CYCLE on RFORCE's
               CYCLE j_do_22
            ENDIF

            DO IGRID = 1,NGRID

               RI  = RGRID(IGRID,1:3)
               DRI = RI - RA

               CALL CROSS ( ANG_VEL, DRI , DUM1 )
               CALL CROSS ( ANG_VEL, DUM1, DUM2 )
               CALL CROSS ( ANG_ACC, DRI , DUM3 )
               
               ! DUM3 is the component of linear accleration due to angular acceleration.
               ! However, in Nastran, the applied force due to angular acceleration is
               ! defined to be in the same direction as the acceleration so we negate
               ! it here for compatibility with F=-MA later.

               ! DUM2 is centripetal acceleration which is already consistent with F=-MA.

               ACCEL_I_T1 = DUM2 - DUM3
               ACCEL_I_R1 = -ANG_ACC

               ACID_G = GRID(IGRID,3)                      ! The global coord sys for this grid is ACID_G
               IF (ACID_G /= 0) THEN                       ! ACID_G is not basic so transform coords to global
                  DO L=1,NCORD
                     IF (CORD(L,2) == ACID_G) THEN
                        ICID = L                           ! ICID is the internal coord sys ID corresponding to ACID_G
                        EXIT
                     ENDIF
                  ENDDO
                                                           ! T12 is coord transf matrix that would transf a global vector to basic
!                                                            and its transpose will transform a basic coord sys vector to global
                  CALL GEN_T0L ( IGRID, ICID, THETAD, PHID, T12 )

                  CALL MATMULT_FFF_T ( T12, ACCEL_I_T1, 3, 3, 1, ACCEL_I_T2 )
                  CALL MATMULT_FFF_T ( T12, ACCEL_I_R1, 3, 3, 1, ACCEL_I_R2 )
                  ACCEL_I(1:3) = ACCEL_I_T2                ! ACCEL_I (grid accel) is now in global coords if it was not before
                  ACCEL_I(4:6) = ACCEL_I_R2
               ELSE
                  ACCEL_I(1:3) = ACCEL_I_T1                ! ACCEL_I (grid accel) is now in global coords if it was not before
                  ACCEL_I(4:6) = ACCEL_I_R1
               ENDIF

               IF (GRID(IGRID,6) == 1) THEN                ! Scalar point so do not generate grav load on it. Give warn, not fatal
                  IERRT = IERRT + 1
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,9901) GRID_ID(IGRID)
                  IF (SUPWARN == 'N') THEN
                     WRITE(ERR,9901) GRID_ID(IGRID)
                  ENDIF
               ENDIF
               IF (IERRT == 0) THEN
                  CALL GET_GRID_6X6_MASS ( GRID_ID(IGRID), IGRID, GRID_MGG_FND, GRID_MGG )
               ENDIF

               FORCE_I = ZERO

               IF (GRID_MGG_FND == 'Y') THEN
                  !F = -M A. Negative sign because thes are inertia forces are in 
                  !the accelerating reference frame. Eg. a centripetal acceleration 
                  !causes a centrifugal inertia force.
                  CALL MATMULT_FFF ( GRID_MGG, -ACCEL_I, 6, 6, 1, FORCE_I )
               ENDIF

               ROW_NUM_START = TDOF_ROW_START(IGRID)
               DO L = 1,6
                  CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
                  ROW_NUM = ROW_NUM_START + L - 1
                  GDOF = TDOF(ROW_NUM,G_SET_COL_NUM)
                  SYS_LOAD(GDOF,ISUB) = SYS_LOAD(GDOF,ISUB) + FORCE_I(L)
               ENDDO

            ENDDO

         ENDDO j_do_22

         IF (GID_ERR > 0) THEN
            WRITE(ERR,1599) SUBR_NAME,GID_ERR
            WRITE(F06,1599) SUBR_NAME,GID_ERR
            CALL OUTA_HERE ( 'Y' )                         ! Errors from reading RFORCE data, so quit
         ENDIF

      ENDDO i_do2

      RETURN

! **********************************************************************************************************************************

 1599 FORMAT(/,' PROCESSING TERMINATED IN SUBROUTINE ',A,' DUE TO ABOVE LISTED ',I8,' ERRORS')

 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 9901 FORMAT(' *WARNING    : NO RFORCE LOAD IS BEING PROCESSED FOR GRID ',I8,' SINCE IT IS A SCALAR POINT (SPOINT)')

! **********************************************************************************************************************************

      END SUBROUTINE RFORCE_PROC

   END MODULE RFORCE_PROCESSING
