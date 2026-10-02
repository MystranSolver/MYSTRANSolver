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

   MODULE GRID_COORDINATE_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CORD_PROC, GRID_PROC, SEQ_PROC

   CONTAINS

      SUBROUTINE CORD_PROC

! This subroutine calculates coordinate system transformation matrices defined on CORD2R, CORD2C, or CORD2S bulk data
! cards. The final  transformations in this subroutine are between a defined coord systems' principal axes and the
! basic coordinate system.  The principal axes of a coordinate system are the three orthogonal rectangular axes that
! is the base of the coordinate system.  For example, a cylindrical coord system has an implied rectangular system
! from which the radius and angle are measured. The transformations generated in this subroutine are not the final ones
! needed to get a coord transformation for a specific grid point in a cylindrical or spherical coord system since those
! axes depend on the radius and angles in those systems. These transformations are only between the basic system and the
! principal axes of the  coord system.  The final transformations for specific grid points in cylindrical and spherical
! coordinate systems is generated from the transformations in this subroutine and is done in subroutine GEN_T0L.FOR.

! The transformation matrices are stored, temporarily, in a 3-D array TN(I,J,K) where K is the internal coord sys
! number (NCORD of them). Thus, TN can be looked at as NCORD 2-D arrays with each 2-D array being a 3x3 coord
! transformation matrix.

!               TN: Transforms a vector to its reference system, RID, from its master system CID

! At the conclusion of this subroutine, the transformation matrices are stored in array RGRID
! along with the basic coordinates of the origin of the coord system. This data is needed in GEN_T0L.FOR

! The subroutine is divided into 8 Phases with the description of each phase given in that section below. The final result is arrays
! CORD(I) (modified so that all reference coord systems are basic), RCORD(I) (which has 12 cols: 1-3 are basic coords of origin,
! 4-6 are ist row of TN, 7-9 are 2nd row and 10-12 are 3rd row. TN is the 3x3 coord transformation matrix which transforms a vector
! from CID to basic.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE80, PI, CONV_DEG_RAD
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MRCORD, NCORD, NCORD1, NCORD2, NGRID, FATAL_ERR
      USE PARAMS, ONLY                :  EPSIL, PRTCORD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CORD, GRID, RCORD, RGRID, TN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE SORTING, ONLY               :  SORT_INT1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'CORD_PROC'
      CHARACTER( 1*BYTE)              :: ALL_RIDS_0  = 'Y' ! Set to 'Y' when all coord systems' RID is basic
      CHARACTER( 1*BYTE)              :: CIRC_ERR    = 'N' ! Set to 'Y' if a coord sys has a circular reference back to itself
      CHARACTER( 6*BYTE)              :: CORD_NAME(NCORD)  ! Name of coord system (CORD1R, etc)
      CHARACTER( 1*BYTE)              :: FOUND_GA    = 'N' ! Set to 'Y' if grid A for CORD1R coord system is found in array GRID
      CHARACTER( 1*BYTE)              :: FOUND_GB    = 'N' ! Set to 'Y' if grid B for CORD1R coord system is found in array GRID
      CHARACTER( 1*BYTE)              :: FOUND_GC    = 'N' ! Set to 'Y' if grid C for CORD1R coord system is found in array GRID
      CHARACTER( 1*BYTE)              :: FOUND_RID   = 'N' ! Set to 'Y' if the RID for a coord system is a system in array CORD
      CHARACTER( 1*BYTE)              :: TRANS_DONE(NCORD) ! If 'Y' then the transformation to basic for a coord sys has been done

      INTEGER(LONG)                   :: CASCADE_PROC_ARRAY1(3,3)
                                                           ! Array which has info about how many coord systems have to be cascaded
!                                                            in order to get a transformation from a coord sys (actual number CID)
!                                                            to basic. There are 3 cols to allow holding info for the 3 ref pts in
!                                                            a CORD1R. Row 1 has the actual coordinate system number
!                                                            Row 2 has the col number in RID_ARRAY to use when cascading
!                                                            coord systems to get the coord transformation from CID to basic. Row 3
!                                                            has the number of coord systems that need to be cascaded to get the
!                                                            transformation from CID to basic.
!

      INTEGER(LONG)                   :: CASCADE_PROC_ARRAY(3,NCORD)
                                                           ! Same info from CASCADE_PROC_ARRAY1 but values for all coord systems.
!                                                            For  CORD1R, the 3 cols in CASCADE_PROC_ARRAY1 are reduced to have the
!                                                            col which has the shortest path from basic up to the coord system

      INTEGER(LONG)                   :: CID, CIDJ, CIDK   ! Actual coord system numbers
      INTEGER(LONG)                   :: CIDA, CIDB, CIDC  ! Actual coord system numbers
      INTEGER(LONG)                   :: CID_RID0          ! Actual coord sys number for a system that has basic as its' ref system
      INTEGER(LONG)                   :: CORD_TYPE         ! Integer value for coord system type
      INTEGER(LONG)                   :: GA(NCORD)         ! Grid nos put into array CORD for pt A when B.D. was read for CORD1R
      INTEGER(LONG)                   :: GB(NCORD)         ! Grid nos put into array CORD for pt B when B.D. was read for CORD1R
      INTEGER(LONG)                   :: GC(NCORD)         ! Grid nos put into array CORD for pt C when B.D. was read for CORD1R
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: I1                ! A computed index
      INTEGER(LONG)                   :: ICID              ! The internal coord sys ID for an actual coord sys number
      INTEGER(LONG)                   :: ICIDA,ICIDB,ICIDC ! The internal coord sys ID for an actual coord sys number
      INTEGER(LONG)                   :: ICID_RID0         ! The internal coord sys ID for actual coord sys number CID_RID0
      INTEGER(LONG)                   :: INT1,INT2         ! Intermediate variables
      INTEGER(LONG)                   :: IERROR    = 0     ! Count of coord sys errors as we process coord systems
      INTEGER(LONG)                   :: NUM_RIDS          ! Number of ref coord systems specified on one CORDij (e.g. 3 for CORD1R)
      INTEGER(LONG)                   :: NUM_LEFT_AT_BEG   ! Num of CORD1 systems left to process at some point in Phase 7
      INTEGER(LONG)                   :: NUM_LEFT_AT_END   ! Num of CORD1 systems left to process at some point in Phase 7
      INTEGER(LONG)                   :: RID               ! The reference coord sys from a CORD array entry

      INTEGER(LONG)                   :: RID_ARRAY(NCORD+1,3*NCORD1+NCORD2)
                                                           ! An array that shows the chain of coord sys refs for all coord systems
!                                                            There has to be more cols than number of coord systems since CORD1R can
!                                                            take up 3 cols (one for each ref pt on the CORD1R). The top of the
!                                                            chain is the coord system being traced and the references follow. For
!                                                            CORD1R each of the 3 cols will have the same coord sys value in row 1.

      INTEGER(LONG)                   :: RID_ARRAY_COL     ! Col number in RID_ARRAY


      REAL(DOUBLE)                    :: EMTN(3,3)         ! A coord transf matrix from some coord system to basic
      REAL(DOUBLE)                    :: EPS1              ! A small number
      REAL(DOUBLE)                    :: IVEC(3)           ! A unit vector in the x direction of a coord system
      REAL(DOUBLE)                    :: JVEC(3)           ! A unit vector in the y direction of a coord system
      REAL(DOUBLE)                    :: KVEC(3)           ! A unit vector in the z direction of a coord system
      REAL(DOUBLE)                    :: MAGVI             ! Magnitude of VI
      REAL(DOUBLE)                    :: MAGVJ             ! Magnitude of VJ
      REAL(DOUBLE)                    :: MAGVK             ! Magnitude of VK
      REAL(DOUBLE)                    :: PHI               ! Elevation angle in a sph coord system
      REAL(DOUBLE)                    :: RADIUS            ! Radius in a cyl or sph coord system
      REAL(DOUBLE)                    :: RGA(NCORD,3)      ! Coords of origin (pt A) in a coord system definition
      REAL(DOUBLE)                    :: RGB(NCORD,3)      ! Coords of a pt (B) on the z axis in a coord system definition
      REAL(DOUBLE)                    :: RGC(NCORD,3)      ! Coords of a pt (C) in the x-z plane in a coord system definition
      REAL(DOUBLE)                    :: RO(3,NCORD)       ! Array of basic coords of origin of coord sys rel to origin of basic sys
      REAL(DOUBLE)                    :: ROJ(3)            ! A column from array RO
      REAL(DOUBLE)                    :: RP(3,NCORD)       ! Array of basic coords of origin of coord sys rel to their ref sys
      REAL(DOUBLE)                    :: RPJ(3)            ! A column from array RP
      REAL(DOUBLE)                    :: T0A(3,3)          ! A coord transformation matrix in an intermediate calc
      REAL(DOUBLE)                    :: T0B(3,3)          ! A coord transformation matrix in an intermediate calc
      REAL(DOUBLE)                    :: T0C(3,3)          ! A coord transformation matrix in an intermediate calc
      REAL(DOUBLE)                    :: THETA             ! Azimuth angle in a cyl or sph coord system
      REAL(DOUBLE)                    :: V0A(3)            ! Vector from the basic system origin to pt A on a CORD1 entry
      REAL(DOUBLE)                    :: V0B(3)            ! Vector from the basic system origin to pt B on a CORD1 entry
      REAL(DOUBLE)                    :: V0C(3)            ! Vector from the basic system origin to pt C on a CORD1 entry
      REAL(DOUBLE)                    :: V1(3),V2(3),V3(3) ! Vectors in an intermediate calc
      REAL(DOUBLE)                    :: V_0_CIDA(3)       ! Vector from the basic system origin to origin of coord sys A
      REAL(DOUBLE)                    :: V_0_CIDB(3)       ! Vector from the basic system origin to origin of coord sys B
      REAL(DOUBLE)                    :: V_0_CIDC(3)       ! Vector from the basic system origin to origin of coord sys C
      REAL(DOUBLE)                    :: V_CIDA_CID_PTA(3) ! Vector from the orogin of ref sys CIDA to point A of coord sys CID
      REAL(DOUBLE)                    :: V_CIDB_CID_PTB(3) ! Vector from the origin of ref sys CIDB to point B of coord sys CID
      REAL(DOUBLE)                    :: V_CIDC_CID_PTC(3) ! Vector from the origin of ref sys CIDC to point C of coord sys CID
      REAL(DOUBLE)                    :: VI(3)             ! Vector in the x-z plane of a coord system
      REAL(DOUBLE)                    :: VJ(3)             ! Vector in the y direction in a coord system
      REAL(DOUBLE)                    :: VK(3)             ! Vector in the z direction in a coord system

      INTRINSIC                       :: DCOS, DSIN, DSQRT



! **********************************************************************************************************************************
! Initialize

      EPS1    = DABS(EPSIL(1))

      DO I=1,NCORD
         TRANS_DONE = 'N'
         IF (CORD(I,1) == 11) CORD_NAME(I) = 'CORD1R'
         IF (CORD(I,1) == 12) CORD_NAME(I) = 'CORD1C'
         IF (CORD(I,1) == 13) CORD_NAME(I) = 'CORD1S'
         IF (CORD(I,1) == 21) CORD_NAME(I) = 'CORD2R'
         IF (CORD(I,1) == 22) CORD_NAME(I) = 'CORD2C'
         IF (CORD(I,1) == 23) CORD_NAME(I) = 'CORD2S'
      ENDDO

      DO I=1,NCORD
         DO J=1,3                                          ! Initialize
            RGA(I,J) = ZERO
            RGB(I,J) = ZERO
            RGC(I,J) = ZERO
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Phase 1: For CORD1R: (1) Change grid numbers (that were put into array CORD when B.D. was read) to the reference coordinate system
!                          number for that grid.
!                      (2) Put coords of the reference grids on the CORD1R into array RCORD.

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '11' )

      IERROR = 0

      CALL CORDCHK ( IERROR )                              ! CORDCHK makes sure that all cord system ID's are unique

doi11:DO I=1,NCORD

         CORD_TYPE = CORD(I,1)
         IF ((CORD_TYPE == 11) .OR. (CORD_TYPE == 12) .OR. (CORD_TYPE == 13)) THEN

            FOUND_GA = 'N'                                 ! Chg GA,GB,GC vals in array CORD with the ref coord sys for those grids
            FOUND_GB = 'N'
            FOUND_GC = 'N'

            GA(I) = CORD(I,3)
            GB(I) = CORD(I,4)
            GC(I) = CORD(I,5)

doj11:      DO J=1,NGRID
               IF ((FOUND_GA == 'Y') .AND. (FOUND_GB == 'Y') .AND. (FOUND_GC == 'Y')) EXIT doj11
               IF      (GRID(J,1) == GA(I)) THEN
                  FOUND_GA = 'Y'
                  CORD(I,3) = GRID(J,2)
                  DO K=1,3
                     RCORD(I,K) = RGRID(J,K)
                     RGA(I,K)   = RGRID(J,K)
                  ENDDO
                  CYCLE doj11
               ELSE IF (GRID(J,1) == GB(I)) THEN
                  FOUND_GB = 'Y'
                  CORD(I,4) = GRID(J,2)
                  DO K=1,3
                     RCORD(I,K+3) = RGRID(J,K)
                     RGB(I,K)     = RGRID(J,K)
                  ENDDO
                  CYCLE doj11
               ELSE IF (GRID(J,1) == GC(I)) THEN
                  FOUND_GC = 'Y'
                  CORD(I,5) = GRID(J,2)
                  DO K=1,3
                     RCORD(I,K+6) = RGRID(J,K)
                     RGC(I,K)     = RGRID(J,K)
                  ENDDO
                  CYCLE doj11
               ENDIF
            ENDDO doj11

            IF (FOUND_GA == 'N') THEN
               IERROR = IERROR + 1
               WRITE(ERR,1315) GA, CORD(I,2)
               WRITE(F06,1315) GA, CORD(I,2)
            ENDIF

            IF (FOUND_GB == 'N') THEN
               IERROR = IERROR + 1
               WRITE(ERR,1315) GB, CORD(I,2)
               WRITE(F06,1315) GB, CORD(I,2)
            ENDIF

            IF (FOUND_GC == 'N') THEN
               IERROR = IERROR + 1
               WRITE(ERR,1315) GC, CORD(I,2)
               WRITE(F06,1315) GC, CORD(I,2)
            ENDIF

            IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '12' )

         ENDIF

      ENDDO doi11

! Check IERROR and quit if > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Phase 2: Change coords on cylindrical and spherical systems to coords in terms of the defining rectangular axes of that system

      DO I = 1,NCORD

         CORD_TYPE = CORD(I,1)

         CID  = CORD(I,2)                                  ! CID is the actual coord system no. for internal no. I.
         RID  = CORD(I,3)                                  ! RID is the ref coord sys for coord sys I

         IF (RID /= 0) THEN                                ! If RID is not basic, determine the coord type for RID.
            DO J=1,NCORD
               CIDJ = CORD(J,2)
               IF(CIDJ == RID) THEN
                  CORD_TYPE = CORD(J,1)
                  EXIT
               ENDIF
            ENDDO

            IF      ((CORD_TYPE == 12) .OR. (CORD_TYPE == 22)) THEN
               DO J=1,7,3                                  ! If RID is cyl, replace R and THETA in RCORD w/ X, Y coords that are
                  RADIUS = RCORD(I,J)                      ! the defining rectangular axes of the cyl system with the X-Y plane
                  THETA  = RCORD(I,J+1)*CONV_DEG_RAD       ! at THETA = 0.
                  RCORD(I,J)   = RADIUS*DCOS(THETA)
                  RCORD(I,J+1) = RADIUS*DSIN(THETA)
               ENDDO
            ELSE IF ((CORD_TYPE == 13) .OR. (CORD_TYPE == 23)) THEN
               DO J=1,7,3                                  ! If RID is cyl, replace R, PHI, THETA in RCORD w/ X, Y, Z coords that
                  RADIUS = RCORD(I,J)                      ! are the defining rectangular axes of the cyl system with the X-Y plane
                  THETA  = RCORD(I,J+1)*CONV_DEG_RAD       ! at THETA = 0.
                  PHI    = RCORD(I,J+2)*CONV_DEG_RAD
                  RCORD(I,J)   = RADIUS*DSIN(THETA)*DCOS(PHI)
                  RCORD(I,J+1) = RADIUS*DSIN(THETA)*DSIN(PHI)
                  RCORD(I,J+2) = RADIUS*DCOS(THETA)
               ENDDO
            ENDIF
         ENDIF

      ENDDO

      IF (PRTCORD == 2) CALL PARAM_PRTCORD_OUTPUT ( '21' )

! **********************************************************************************************************************************
! Phase 3: Check coordinate system data for logic.
!          (1) There must be a final reference to the basic system (0)
!          (2) There can be no circular references (i.e. a coord system can reference any number of other coord sys but cannot
!              reference itself in that chain of refs).
!          (3) Create RID_ARRAY which shows the chain of references of each coord system

! A complete check is made and processing is not stopped until all errors of the above kind are found. Bulk Data parameter
! PRTCORD is used to  print out information from these checks as well as the coord transformation matrices later.

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '31' )

      IERROR = 0                                           ! Initialize IERROR for Phase 2

      IF (IERROR > 0) THEN
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Begin loop over all coord systems

      RID_ARRAY_COL = 0
doi21:DO I=1,NCORD

         CID = CORD(I,2)                                   ! CID is the coordinate system number from the coord card

         CORD_TYPE = CORD(I,1)
         IF  ((CORD_TYPE == 11) .OR. (CORD_TYPE == 12) .OR. (CORD_TYPE == 13)) THEN
            NUM_RIDS = 3                                   ! CORD1R has 3 ref coord systems
         ELSE
            NUM_RIDS = 1                                   ! CORD2R,C,S has 1 ref coord system
         ENDIF

dom21:   DO M=1,NUM_RIDS

            RID_ARRAY_COL = RID_ARRAY_COL + 1

            RID = CORD(I,2+M)                              ! RID is the reference coord system number for CID

doj21:      DO J=1,NCORD+1                                 ! Initially set all of RID_ARRAY to some negative number since we know,
               RID_ARRAY(J,RID_ARRAY_COL) = -99            ! after this CID has been traced, it cannot have any negative coord ID's
            ENDDO doj21                                    ! but it must have 0 (basic) as the last coord sys in the chain

            I1  = 1                                        ! We need I1 to increment this index later
            RID_ARRAY(I1,RID_ARRAY_COL) = CID
            I1  = 2
            RID_ARRAY(I1,RID_ARRAY_COL) = RID

            IF (PRTCORD == 2) CALL PARAM_PRTCORD_OUTPUT ( '32' )

            IF (RID /= 0) THEN                             ! If RID is not the basic system, make sure it is defined

               FOUND_RID = 'N'
doj22:         DO J=1,NCORD
                  IF (CORD(J,2) == RID) THEN
                     FOUND_RID = 'Y'
                     EXIT doj22
                  ENDIF
               ENDDO doj22

               IF (FOUND_RID == 'N') THEN                  ! Could not find reference coord sys (RID) for coord sys CID
                  WRITE(ERR,1300) RID,CID
                  WRITE(F06,1300) RID,CID
                  IERROR = IERROR + 1
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF

! Build chain of references for CID. Initially, RID was read from field 3 of the CORD card, above.
! Now go thru all CORD cards until we find one that has RID as the defined system in field 2. When we find this CORD
! card, we reset RID to be the value from field 3 of that card and continue to cascade down until we come to a CORD
! card with 0 (basic) in field 3. When we find this system, exit from the loops. Two loops on NCORD are needed to do
! the check. The inner (K) loop steps thru the CORD cards looking for the one that defines the RID system. When it is
! found, we continue in that loop looking for the next system (if necessary). However, that system may have already
! been passed in the K loop, so the outer loop starts the process over if necessary.

               CIRC_ERR = 'N'
doj23:         DO J=1,NCORD
                  IF (J /= I) THEN
dok21:               DO K=1,NCORD
                        IF (K /= I) THEN
                           IF (CORD(K,2) == RID) THEN
                              CIDK = CORD(K,2)
                              RID  = CORD(K,3)
                              I1 = I1 + 1
                              IF (I1 > NCORD+1) THEN
                                 WRITE(ERR,1305) SUBR_NAME
                                 WRITE(F06,1305) SUBR_NAME
                                 FATAL_ERR = FATAL_ERR + 1
                                 CALL OUTA_HERE ( 'Y' )
                              ENDIF
                              RID_ARRAY(I1,RID_ARRAY_COL) = RID

                              IF (PRTCORD == 2) CALL PARAM_PRTCORD_OUTPUT ( '33' )

                              IF (RID == 0) EXIT dok21

dol21:                        DO L=1,I1-1                  ! Check for circular reference
                                 IF(RID_ARRAY(I1,RID_ARRAY_COL) == RID_ARRAY(L,RID_ARRAY_COL)) THEN
                                    WRITE(ERR,1302) RID_ARRAY(I1,RID_ARRAY_COL)
                                    WRITE(F06,1302) RID_ARRAY(I1,RID_ARRAY_COL)
                                    FATAL_ERR = FATAL_ERR + 1
                                    IERROR = IERROR + 1
                                    CIRC_ERR = 'Y'
                                    EXIT dol21
                                 ENDIF
                              ENDDO dol21

                           ENDIF
                        ENDIF
                        IF ((RID ==  0 ) .OR. (CIRC_ERR == 'Y')) EXIT dok21
                     ENDDO dok21
                     IF ((RID ==  0 ) .OR. (CIRC_ERR == 'Y')) EXIT doj23
                  ENDIF
               ENDDO doj23

               IF (RID_ARRAY(I1,RID_ARRAY_COL) == 0) THEN  ! Check if last sys in RID_ARRAY is basic. If it is, CYCLE back to the I
                  CYCLE dom21                              ! loop to begin with a new coord system.
               ELSE                                        ! Otherwise write error and CYCLE
                  WRITE(ERR,1303) CID
                  WRITE(F06,1303) CID
                  IERROR = IERROR + 1                      ! Coord sys CID does not have a final reference to basic
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE dom21
               ENDIF

            ENDIF

         ENDDO dom21

      ENDDO doi21

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '34' )

! Check IERROR and quit if > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Phase 4: Solve for array CASCADE_PROC_ARRAY. It will be used to do the cascading of coord references to get the transformation
! from CID to basic

      DO J=1,NCORD
         CASCADE_PROC_ARRAY(1,J) = 0
         CASCADE_PROC_ARRAY(2,J) = 0
         CASCADE_PROC_ARRAY(3,J) = NCORD + 1
      ENDDO

      RID_ARRAY_COL = 1
doi31:DO I=1,NCORD

         CORD_TYPE = CORD(I,1)
         IF  ((CORD_TYPE == 11) .OR. (CORD_TYPE == 12) .OR. (CORD_TYPE == 13)) THEN
            NUM_RIDS = 3                                   ! CORD1R,C,S has 3 ref coord systems
         ELSE
            NUM_RIDS = 1                                   ! CORD2R,C,S have 1 ref coord system
         ENDIF

         DO J=1,3
            CASCADE_PROC_ARRAY1(1,J) = 0
            CASCADE_PROC_ARRAY1(2,J) = 0
            CASCADE_PROC_ARRAY1(3,J) = NCORD + 1
         ENDDO

doj31:   DO J=1,NUM_RIDS
doK31:      DO K=NCORD+1,1,-1
               IF (RID_ARRAY(K,RID_ARRAY_COL) == 0) THEN
                  IF (K <= CASCADE_PROC_ARRAY1(3,J)) THEN

                     CASCADE_PROC_ARRAY1(1,J) = RID_ARRAY(1,RID_ARRAY_COL)
                     CASCADE_PROC_ARRAY1(2,J) = RID_ARRAY_COL
                     CASCADE_PROC_ARRAY1(3,J) = K - 1

                     CASCADE_PROC_ARRAY(1,I)  = RID_ARRAY(1,RID_ARRAY_COL)
                     CASCADE_PROC_ARRAY(2,I)  = RID_ARRAY_COL
                     CASCADE_PROC_ARRAY(3,I)  = K - 1

                     CYCLE
                  ENDIF
               ENDIF
            ENDDO dok31
            RID_ARRAY_COL = RID_ARRAY_COL + 1
         ENDDO doj31

         IF (NUM_RIDS == 3) THEN
            DO J=2,3
               IF(CASCADE_PROC_ARRAY1(3,J) < CASCADE_PROC_ARRAY1(3,J-1)) THEN
                  CASCADE_PROC_ARRAY(1,I) = CASCADE_PROC_ARRAY1(1,J)  ! Coord system ID (CID)
                  CASCADE_PROC_ARRAY(2,I) = CASCADE_PROC_ARRAY1(2,J)  ! Row in RID_ARRAY where to find the coord references for CID
                  CASCADE_PROC_ARRAY(3,I) = CASCADE_PROC_ARRAY1(3,J)  !
               ENDIF
            ENDDO
         ENDIF

      ENDDO doi31

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '41' )

! **********************************************************************************************************************************
! Phase 5: Get the unit vectors in each of the CORD2 coord systems in terms of their reference system. Put the results into array
! RCORD where the transformation matrices go. These are NOT the final transformations; they are the transformations from a CID
! to its reference system and ONLY for CORD2 systems (CORD1 systams will be handled later).

      IF (PRTCORD == 2) CALL PARAM_PRTCORD_OUTPUT ( '51' )

      IERROR = 0
doi32:DO I=1,NCORD

         CORD_TYPE = CORD(I,1)                             ! Process only CORD2C,R,S first
         IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN

            RID_ARRAY_COL = CASCADE_PROC_ARRAY(2,I)
            IF (RID_ARRAY_COL <= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1321) SUBR_NAME, 'RID_ARRAY_COL in the Phase 5 calcs', RID_ARRAY_COL
               WRITE(F06,1321) SUBR_NAME, 'RID_ARRAY_COL in the Phase 5 calcs', RID_ARRAY_COL
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            CID = RID_ARRAY(1,RID_ARRAY_COL)               ! Coord sys for which we are getting the transformation to its RID
            RID = RID_ARRAY(2,RID_ARRAY_COL)               ! The  reference system, RID, for CID

            DO J=1,NCORD                                   ! Find the internal coord sys ID for CID
               IF (CORD(J,2) == CID) THEN
                  ICID = J
                  EXIT
               ENDIF
            ENDDO

            DO J=1,3
               VK(J) = RCORD(ICID,3+J) - RCORD(ICID,J)     ! Vector in z dir of CID expressed in it's ref system
               VI(J) = RCORD(ICID,6+J) - RCORD(ICID,J)     ! Vector in x-z plane of CID expressed in it's ref system
            ENDDO
            MAGVK = DSQRT( VK(1)*VK(1) + VK(2)*VK(2) + VK(3)*VK(3) )
            MAGVI = DSQRT( VI(1)*VI(1) + VI(2)*VI(2) + VI(3)*VI(3) )

            IF (MAGVK < EPS1) THEN                         ! If MAGVK = 0 then points A and B on the coord card are too close.
               WRITE(ERR,1306) CORD(ICID,2)
               WRITE(F06,1306) CORD(ICID,2)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

            IF (MAGVI < EPS1) THEN                         ! If MAGVI = 0 then points A and C on the coord card are too close.
               WRITE(ERR,1307) CORD(ICID,2)
               WRITE(F06,1307) CORD(ICID,2)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

            DO J=1,3                                       ! Unit vector in z direction is KVEC
               KVEC(J) = VK(J)/MAGVK
               RCORD(ICID,9+J) = KVEC(J)
            ENDDO

            IF ((MAGVK >= EPS1) .AND. (MAGVI >= EPS1)) THEN! Calc unit vector in y dir if vectors VK, VI are not null
               CALL CROSS ( KVEC, VI, VJ )
               MAGVJ = DSQRT(VJ(1)*VJ(1)+VJ(2)*VJ(2)+VJ(3)*VJ(3))

               IF (MAGVJ < EPS1) THEN                      ! If MAGVJ = 0 then vec from pt A to pt B is || to vec from pts A to C
                  WRITE(ERR,1308) CORD(ICID,2)
                  WRITE(F06,1308) CORD(ICID,2)
                  IERROR = IERROR + 1
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDIF

            IF (IERROR > 0) CYCLE

            DO J = 1,3                                     ! Unit vector in y direction is JVEC
               JVEC(J) = VJ(J)/MAGVJ
               RCORD(ICID,6+J) = JVEC(J)
            ENDDO

            CALL CROSS ( JVEC, KVEC, IVEC )                ! Unit vector in x direction is IVEC


            DO J = 1,3                                     ! Put unit vecs into array TN: Transforms a vector to RID from CID
               TN(J,1,ICID) = IVEC(J)                      ! -- Col 1 of TN is IVEC
               TN(J,2,ICID) = JVEC(J)                      ! -- Col 2 of TN is JVEC
               TN(J,3,ICID) = KVEC(J)                      ! -- Col 3 of TN is KVEC
            ENDDO

            DO J=1,3
               RCORD(ICID,3+J) = TN(1,J,ICID)              ! Row 1 of TN goes into RCORD cols  4- 6
               RCORD(ICID,6+J) = TN(2,J,ICID)              ! Row 1 of TN goes into RCORD cols  7- 9
               RCORD(ICID,9+J) = TN(3,J,ICID)              ! Row 1 of TN goes into RCORD cols 10-12
            ENDDO

            IF (PRTCORD == 2) CALL PARAM_PRTCORD_OUTPUT ( '52' )

         ENDIF

      ENDDO doi32

! Check IERROR and quit if > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Phase 6:  Calculate transformations to basic (zero ) system for all CORD2 coordinate systems.

! The TN transformation matrices generated in Phase 5 give the transformation from CID system to RID system for one coord card.
! In Phase 6 we cascade these to get the transformation from basic to CID system for each CORD2 coord system. That is, consider the
! following example:

!             CID             RID             TN from Phase 2
!     CORD2R   1               0                   TN01
!     CORD2R   2               1                   TN12
!     CORD2C   3               2                   TN23
!     CORD2S   4               1                   TN14

! Where, e.g., TN12 is the transformation matrix from coord sys 2 to coord sys 1 (princ. axes) found in Phase 5.
! The transformation  from 2 to 0 is TN01 x TN12.  The transformations generated here are:

!       TN02 = TN01 x TN12

! and stores the result back into array TN and resets the RID to 0 for CID 2

!       TN03 = TN02 x TN23
!       TN04 = TN01 x TN14

! and TN01 was already determined in Phase 5.

! The procedure is to work backwards.  Begin looking for a coord system that has 0 as a reference (or RID = 0).  The coord
! transformation to basic for that system is already in TN.  If that coord system is referenced on another coord
! card, then we multiply those two matrices together to get the transformation to basic for the second system, and
! so on. This transformation matrix is stored over the original TN found in Phase 5.  Once we have gotten the
! transformation to basic for a CORD2 coord system,  we reset RID (i.e. CORD(I,3)) to zero to indicate that its' TN is
! now referenced to basic.

! The basic coords of the origin of each coord system are also calculated. The calc begins by taking the point A coords from the
! system that has 0 as reference and adding the amounts from each coord sys A point after transforming it to basic.

      IERROR = 0
main: DO                                                   ! Until all RID's are 0
         ALL_RIDS_0 = 'Y'                                  ! Find out if there are coord systems whose RID is not, or has not been
         DO J=1,NCORD                                      ! transformed to, basic
            CORD_TYPE = CORD(J,1)
            IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN
               IF (CORD(J,3) /= 0) THEN
                  ALL_RIDS_0 = 'N'
                  EXIT
               ELSE
                  TRANS_DONE(J) = 'Y'                      ! This coord sys has basic as its RID so set TRANS_DONE
               ENDIF
            ENDIF
         ENDDO

         IF (ALL_RIDS_0 == 'N') THEN                       ! From above search we found at least 1 coord sys whose RID /= 0
                                                           ! The following loops (J, K) search through the coord system list to find
                                                           ! a sys that has RID = 0. When found, it sets this sys CID to be CID_RID0
jloop1:     DO J =1,NCORD                                  ! The inner loop then searches for a sys that uses this CID as its RID.
                                                           ! When found, it exits the outer loop, does the coord transformation, and
               CORD_TYPE = CORD(J,1)                       ! resets the RID on this system to 0
               IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN

                  RID = CORD(J,3)

                  CID_RID0  = 0
                  ICID_RID0 = 0
                  IF (RID == 0) THEN                       ! We have found a sys that has RID = 0. Set CID_RID0 to that systems CID
                     CID       = CORD(J,2)
                     CID_RID0  = CID
                     ICID_RID0 = J

kloop2:              DO K = 1,NCORD                        ! Go through list looking for the sys that has, as its' RID, CID_RID0
                        CORD_TYPE = CORD(K,1)
                        IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN
                           IF (K == J) CYCLE kloop2        ! Skip the sys found in the J loop above since it has RID = 0 already
                           CID = CORD(K,2)
                           RID = CORD(K,3)
                           IF (RID == CID_RID0) THEN       ! We have found one coord sys that has its' RID equal to CID_RID0
                              EXIT jloop1                  ! so exit this loop and get the transformation to basic for system CID
                           ELSE
                           ENDIF
                        ENDIF
                     ENDDO kloop2

                  ELSE

                     CYCLE

                  ENDIF

               ENDIF

            ENDDO jloop1

            IF (ICID_RID0 == 0) THEN
               WRITE(ERR,1316) SUBR_NAME, ICID_RID0
               WRITE(F06,1316) SUBR_NAME, ICID_RID0
               IERROR = IERROR + 1
               CYCLE main
            ENDIF

            RO(1,ICID_RID0) = RCORD(ICID_RID0,1)
            RO(2,ICID_RID0) = RCORD(ICID_RID0,2)
            RO(3,ICID_RID0) = RCORD(ICID_RID0,3)

jloop2:     DO J = 1,NCORD                                 ! Loop on J since there may be more than 1 sys that references CID_RID0

               CORD_TYPE = CORD(J,1)                       ! resets the RID on this system to 0
               IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN

                  IF (J == ICID_RID0) CYCLE jloop2

                  RID = CORD(J,3)

                  IF (RID == CID_RID0) THEN                ! This is a coord sys that has CID_RID0 as its RID, so we can mult
                                                           ! matrices to get its' transformation to basic

                     CALL MATMULT_FFF ( TN(1,1,ICID_RID0), TN(1,1,J), 3, 3, 3, EMTN )

                     RO(1,J) = RCORD(J,1)
                     RO(2,J) = RCORD(J,2)
                     RO(3,J) = RCORD(J,3)

                     DO K=1,3
                        ROJ(K) = RO(K,J)
                     ENDDO
                     CALL MATMULT_FFF ( TN(1,1,ICID_RID0), ROJ, 3, 3, 1, RPJ )
                     DO K=1,3
                        RP(K,J) = RPJ(K)
                     ENDDO

                     DO K = 1,3
                        RCORD(J,K) = RO(K,ICID_RID0) + RP(K,J)! Coords of origin of system J in basic coords
                        DO L = 1,3
                           TN(K,L,J) = EMTN(K,L)
                        ENDDO
                     ENDDO

                     CORD(J,3) = 0                         ! Reset RID on this coord sys to 0 since TN matrix now is trans to basic
                     TRANS_DONE(J) = 'Y'                   ! We have completed the transformation to basic for this coord sys

                  ENDIF

               ENDIF

            ENDDO jloop2

         ELSE

            EXIT main

         ENDIF

      ENDDO main

! Check IERROR and quit if > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Now put transformation for the CORD2 systems matrices in RCORD. First 3 words in RCORD are now the basic coords of the origin of
! the coord system. Next 9 words give the 3 rows (row 1 followed by row 2 then row 3) of the transformation matrix which will
! transform a vector in CID (of coord principal axes) to a vector in basic coord system.

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '61' )

      DO I=1,NCORD
         CORD_TYPE = CORD(I,1)
         IF ((CORD_TYPE == 21) .OR. (CORD_TYPE == 22) .OR. (CORD_TYPE == 23)) THEN
            DO J=1,3
               DO K=1,3
                  L = 3 + 3*(J-1) + K
                  RCORD(I,L) = TN(J,K,I)
               ENDDO
            ENDDO
            IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '62' )
         ENDIF
      ENDDO

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '63' )

! **********************************************************************************************************************************
! Phase 7 Process CORD1 entries to get transformations to basic and put results into array RCORD

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '71' )

main7:DO I=1,NCORD

         CORD_TYPE = CORD(I,1)                             ! Only process CORD1 entries. CORD2's were handled above
         IF ((CORD_TYPE == 11) .OR. (CORD_TYPE == 12) .OR. (CORD_TYPE == 13)) THEN

            NUM_LEFT_AT_BEG = 0                            ! Find out how many CORD1 systems still need to be processed
            DO K=1,NCORD
               IF (TRANS_DONE(K) == 'N') THEN
                  NUM_LEFT_AT_BEG = NUM_LEFT_AT_BEG + 1
               ENDIF
            ENDDO
            IF (NUM_LEFT_AT_BEG == 0) EXIT main7

            IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '72' )

            IERROR   = 0
big_loop:   DO J=1,NCORD                                   ! Find a CORD1 with all RID's already processed and calc RCORD, TN for it
               IF (TRANS_DONE(J) == 'Y') THEN
                  CYCLE big_loop
               ENDIF

               CID  = CORD(J,2)                            ! Coord and reference ID's on this CORD1
               CIDA = CORD(J,3)
               CIDB = CORD(J,4)
               CIDC = CORD(J,5)

               ICIDA = -99                                 ! Get the internal coord ID's for the ref systems on this CORD1
               ICIDB = -99
               ICIDC = -99
               DO K=1,NCORD

                  IF (CORD(K,2) == CID ) THEN
                     ICID  = K
                  ENDIF

                  IF (CORD(J,3) /= 0) THEN
                     IF (CORD(K,2) == CIDA) ICIDA = K
                  ELSE
                     ICIDA = 0
                  ENDIF

                  IF (CORD(J,4) /= 0) THEN
                     IF (CORD(K,2) == CIDB) ICIDB = K
                  ELSE
                     ICIDB = 0
                  ENDIF

                  IF (CORD(J,5) /= 0) THEN
                     IF (CORD(K,2) == CIDC) ICIDC = K
                  ELSE
                     ICIDC = 0
                  ENDIF

               ENDDO

               IF      (ICIDA == -99) THEN                 ! Set error and cycle if we couldn't find the internal ID's for CIDA,B,C
                  IERROR = IERROR + 1
                  INT1 = CIDA   ;   INT2 = ICIDA
               ELSE IF (ICIDB == -99) THEN
                  IERROR = IERROR + 1
                  INT1 = CIDB   ;   INT2 = ICIDB
               ELSE IF (ICIDC== -99) THEN
                  IERROR = IERROR + 1
                  INT1 = CIDC   ;   INT2 = ICIDC
               ENDIF
               IF ((ICIDA == -99) .OR. (ICIDB == -99) .OR. (ICIDC == -99)) THEN
                  IERROR = IERROR + 1
                  WRITE(ERR,1317) CORD_NAME(J), INT1, INT2
                  WRITE(F06,1317) CORD_NAME(J), INT1, INT2
                  CYCLE main7
               ENDIF


               IF (ICIDA /= 0) THEN
                  IF (TRANS_DONE(ICIDA) == 'N') CYCLE big_loop
               ENDIF
               IF (ICIDB /= 0) THEN
                  IF (TRANS_DONE(ICIDB) == 'N') CYCLE big_loop
               ENDIF
               IF (ICIDC /= 0) THEN
                  IF (TRANS_DONE(ICIDC) == 'N') CYCLE big_loop
               ENDIF

               IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '73' )

               DO K=1,3
                  DO L=1,3
                     IF (ICIDA /= 0) THEN
                        TN(K,L,ICIDA) = RCORD(ICIDA,3*(K-1)+L+3) ! This TN matrix will transform a vector to basic (0) from CIDA
                        T0A(K,L) = TN(K,L,ICIDA)
                     ENDIF
                     IF (ICIDB /= 0) THEN
                        TN(K,L,ICIDB) = RCORD(ICIDB,3*(K-1)+L+3) ! This TN matrix will transform a vector to basic (0) from CIDB
                        T0B(K,L) = TN(K,L,ICIDB)
                     ENDIF
                     IF (ICIDC /= 0) THEN
                        TN(K,L,ICIDC) = RCORD(ICIDC,3*(K-1)+L+3) ! This TN matrix will transform a vector to basic (0) from CIDC
                        T0C(K,L) = TN(K,L,ICIDC)
                     ENDIF
                  ENDDO
               ENDDO

               DO K=1,3                                    ! Vectors from basic system origin to origins of reference systems A,B,C

                  IF (ICIDA /= 0) THEN
                     V_0_CIDA(K) = RCORD(ICIDA,K)
                  ELSE
                     V_0_CIDA(K) = ZERO
                  ENDIF

                  IF (ICIDB /= 0) THEN
                     V_0_CIDB(K) = RCORD(ICIDB,K)
                  ELSE
                     V_0_CIDB(K) = ZERO
                  ENDIF

                  IF (ICIDC /= 0) THEN
                     V_0_CIDC(K) = RCORD(ICIDC,K)
                  ELSE
                     V_0_CIDC(K) = ZERO
                  ENDIF

               ENDDO

               DO K=1,3                                    ! Vectors in sys 0 from origin of sys CIDA,B,C to pts A,B,C of sys CID
                  V1(K) = RGA(J,K)
                  V2(K) = RGB(J,K)
                  V3(K) = RGC(J,K)
               ENDDO

               IF (ICIDA /= 0) THEN
                  CALL MATMULT_FFF ( T0A, V1, 3, 3, 1, V_CIDA_CID_PTA )
               ELSE
                  DO K=1,3
                     V_CIDA_CID_PTA(K) = V1(K)
                  ENDDO
               ENDIF

               IF (ICIDB /= 0) THEN
                  CALL MATMULT_FFF ( T0B, V2, 3, 3, 1, V_CIDB_CID_PTB )
               ELSE
                  DO K=1,3
                     V_CIDB_CID_PTB(K) = V2(K)
                  ENDDO
               ENDIF

               IF (ICIDC /= 0) THEN
                  CALL MATMULT_FFF ( T0C, V3, 3, 3, 1, V_CIDC_CID_PTC )
               ELSE
                  DO K=1,3
                     V_CIDC_CID_PTC(K) = V3(K)
                  ENDDO
               ENDIF

! Coords of points A, B, C of CID in their reference systems, CIDA, CIDB, CIDC

               DO K=1,3                                    ! Vectors in sys 0 from sys 0 origin to pts A,B,C on coord sys CID

                  V0A(K) = V_0_CIDA(K) + V_CIDA_CID_PTA(K)
                  V0B(K) = V_0_CIDB(K) + V_CIDB_CID_PTB(K)
                  V0C(K) = V_0_CIDC(K) + V_CIDC_CID_PTC(K)

               ENDDO

               DO K=1,3                                    ! Vectors i sys 0 from A->B and A->C for sys CID
                  VK(K)  = V0B(K) - V0A(K)
                  VI(K)  = V0C(K) - V0A(K)
               ENDDO

               MAGVK = DSQRT( VK(1)*VK(1) + VK(2)*VK(2) + VK(3)*VK(3) )
               MAGVI = DSQRT( VI(1)*VI(1) + VI(2)*VI(2) + VI(3)*VI(3) )

               IF (MAGVK < EPS1) THEN
                  WRITE(ERR,1306) CORD(J,2)
                  WRITE(F06,1306) CORD(J,2)
                  IERROR = IERROR + 1
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF

               IF (MAGVI < EPS1) THEN
                  WRITE(ERR,1307) CORD(J,2)
                  WRITE(F06,1307) CORD(J,2)
                  IERROR = IERROR + 1
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF

               DO K = 1,3                                  ! Unit vector in basic coord z direction for z axis of CID is KVEC
                  KVEC(K) = VK(K)/MAGVK
               ENDDO

               IF ((MAGVK >= EPS1) .AND. (MAGVI >= EPS1)) THEN! Calc unit vector in y dir if vectors VK, VI are not null
                  CALL CROSS ( KVEC, VI, VJ )
                  MAGVJ = DSQRT(VJ(1)*VJ(1)+VJ(2)*VJ(2)+VJ(3)*VJ(3))

                  IF (MAGVJ < EPS1) THEN
                     WRITE(ERR,1308) CORD(J,2)
                     WRITE(F06,1308) CORD(J,2)
                     IERROR = IERROR + 1
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF

               IF (IERROR > 0) CYCLE main7

               DO K = 1,3                                  ! Unit vector in basic coord y direction for z axis of CID is JVEC
                  JVEC(K) = VJ(K)/MAGVJ
               ENDDO
                                                           ! Unit vector in basic coord x direction for z axis of CID is IVEC
               CALL CROSS ( JVEC, KVEC, IVEC )

               IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '74' )

               DO K=1,3                                    ! Put unit vecs into arrays TN and RCORD
                  TN(K,1,ICID)    = IVEC(K)
                  TN(K,2,ICID)    = JVEC(K)
                  TN(K,3,ICID)    = KVEC(K)
               ENDDO

               DO K=1,3
                  RCORD(ICID,K) = V0A(K)                   ! Origin of CID is V0A
                  DO L=1,3
                     M = 3 + 3*(K-1) + L
                     RCORD(ICID,M) = TN(K,L,ICID)
                  ENDDO
               ENDDO

               IF (IERROR == 0) THEN
                  TRANS_DONE(ICID) = 'Y'
                  IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '75' )
              ENDIF

            ENDDO big_loop

            NUM_LEFT_AT_END = 0                            ! Find out how many CORD1 systems still need to be processed
            DO K=1,NCORD
               IF (TRANS_DONE(K) == 'N') THEN
                  NUM_LEFT_AT_END = NUM_LEFT_AT_END + 1
               ENDIF
            ENDDO
            IF (NUM_LEFT_AT_END < NUM_LEFT_AT_BEG) THEN
               CONTINUE
            ELSE
               IERROR = IERROR + 1
               WRITE(ERR,1318)
               DO K=1,NCORD
                  IF (TRANS_DONE(K) == 'N') THEN
                     WRITE(ERR,7001) CORD_NAME(K), CORD(K,2)
                  ENDIF
               ENDDO
               WRITE(F06,1318)
               DO K=1,NCORD
                  IF (TRANS_DONE(K) == 'N') THEN
                     WRITE(F06,7001) CORD_NAME(K), CORD(K,2)
                  ENDIF
               ENDDO
               CALL OUTA_HERE ( 'Y' )
            ENDIF

         ENDIF

      ENDDO main7

! **********************************************************************************************************************************
! Phase 8: Check that all systems have been transformed to basic and rewrite RID's in array CORD to reflect this

      IF (IERROR > 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1304)
         WRITE(F06,1304)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IERROR = 0                                           ! Check to make sure all systems have been transformed to basic
      DO I=1,NCORD
         IF (TRANS_DONE(I) /= 'Y') THEN
            IERROR = IERROR + 1
         ENDIF
      ENDDO
      IF (IERROR == 0) THEN                                ! --- No error so reset all RID's that have not already been reset to 0
         DO I=1,NCORD
            DO J=3,5
               IF (CORD(I,J) /= 0) THEN
                  DO K=1,NCORD
                     IF (CORD(K,2) == CORD(I,J)) ICID = K
                  ENDDO
               ENDIF
               IF (TRANS_DONE(ICID) == 'Y') THEN
                  CORD(I,J) = 0
               ENDIF
            ENDDO
         ENDDO
      ELSE                                                 ! --- Something wrong - all systems have not been transformed to basic
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1319) SUBR_NAME
         DO K=1,NCORD
            IF (TRANS_DONE(K) == 'N') THEN
               WRITE(ERR,7001) CORD_NAME(K), CORD(K,2)
            ENDIF
         ENDDO
         WRITE(F06,1319) SUBR_NAME
         DO K=1,NCORD
            IF (TRANS_DONE(K) == 'N') THEN
               WRITE(F06,7001) CORD_NAME(K), CORD(K,2)
            ENDIF
         ENDDO
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF ((PRTCORD == 1) .OR. (PRTCORD == 2)) CALL PARAM_PRTCORD_OUTPUT ( '81' )

! Check IERROR and quit if > 0



      RETURN

! **********************************************************************************************************************************
 1300 FORMAT(' *ERROR  1300: REFERENCE COORDINATE SYSTEM NUMBER ',I8,' ON  COORDINATE SYSTEM NUMBER ',I8,' IS UNDEFINED')

 1302 FORMAT(' *ERROR  1302: COORDINATE SYSTEM NUMBER ',I8,' HAS A CIRCULAR REFERENCE')

 1303 FORMAT(' *ERROR  1303: COORDINATE SYSTEM NUMBER ',I8,' DOES NOT HAVE A FINAL REFERENCE TO THE BASIC SYSTEM')

 1304 FORMAT(' PROCESSING TERMINATED DUE TO PRIOR COORDINATE SYSTEM ERRORS')

 1305 FORMAT(' *ERROR  1305: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DIMENSION OF ARRAY RID_ARRAY TOO SMALL')

 1306 FORMAT(' *ERROR  1306: POINTS A AND B ON COORDINATE SYSTEM NUMBER ',I8,' ARE TOO CLOSE')

 1307 FORMAT(' *ERROR  1307: POINTS A AND C ON COORDINATE SYSTEM NUMBER ',I8,' ARE TOO CLOSE')

 1308 FORMAT(' *ERROR  1308: VECTOR FROM POINT A TO POINT B, AND VECTOR FROM POINT A TO POINT C ARE PARALLEL ON COORD SYSTEM ',I8)

 1315 FORMAT(' *ERROR  1315: GRID ',I8,' UNDEFINED ON COORDINATE SYSTEM NUMBER ',I8)

 1316 FORMAT(' *ERROR  1316: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VARIABLE ICID_RIDO = ',I8,'. MUST BE A POSITIVE INTEGER')

 1317 FORMAT(' *ERROR  1317: CANNOT FIND COORD TRANSFORMATION MATRIX FOR ',A,I8,' (INTERNAL COORD SYS ID IS = ',I8,')')

 1318 FORMAT(' *ERROR  1318: CANNOT COMPLETE THE TRANSFORMATION TO BASIC FOR THE COORD SYSTEMS LISTED BELOW:')

 1319 FORMAT(' *ERROR  1319: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TRANSFORMATION TO BASIC HAS NOT BEEN COMPLETED FOR THE FOLLOWING COORD SYSTEMS:')

 1321 FORMAT(' *ERROR  1403: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                     ,/,14X,' VARIABLE ',A,' HAS AN INVALID VALUE = ',I8)

 7001 FORMAT('               ',A,2X,I8)


! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CORDCHK ( IERROR )

! Checks array CORD to make sure that there are not more than 1 coord systems with the same ID. It does this by
! sorting the coord system ID's and then checking the sorted dummy array for uniqueness

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NCORD, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CORD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CORDCHK'

      INTEGER(LONG), INTENT(OUT)      :: IERROR            ! Count of the number of duplicate coord system ID's
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: DUMCORD(NCORD)    ! Dummy array of coord system ID's sorted




! **********************************************************************************************************************************
! Create DUMCORD to be an array of the coordinate system ID's

      DO I=1,NCORD
         DUMCORD(I) = CORD(I,2)
      ENDDO

! Sort cord system ID's into numerically increasing order using the shell sort method

      IF (NCORD > 1) THEN
         CALL SORT_INT1 ( SUBR_NAME, 'CORD', NCORD, DUMCORD )
      ENDIF

! Check for duplicate numbers

      IERROR = 0
      DO I=1,NCORD-1
         IF (DUMCORD(I) == DUMCORD(I+1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1309) DUMCORD(I+1)
            WRITE(F06,1309) DUMCORD(I+1)
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1309 FORMAT(' *ERROR  1309: COORDINATE SYSTEM NUMBER ',I8,' IS A DUPLICATE.')

! **********************************************************************************************************************************

      END SUBROUTINE CORDCHK

! ##################################################################################################################################

      SUBROUTINE PARAM_PRTCORD_OUTPUT ( WHICH )

      CHARACTER( 2*BYTE)              :: WHICH             ! Decides what to print out for this call to this subr

                                                           ! Char array with RID_ARRAY values
      CHARACTER( 9*BYTE)              :: CRID_ARRAY(NCORD+1,3*NCORD1+NCORD2)
      CHARACTER( 9*BYTE)              :: CRID_ARRAY_BLANK
      CHARACTER( 6*BYTE)              :: NAME1
      CHARACTER(24*BYTE)              :: NAME2

      INTEGER(LONG)                   :: I,J,K             ! Local DO loop indices
      INTEGER(LONG)                   :: II,JJ             ! Local DO loop indices
      INTEGER(LONG)                   :: IA,IB,IC          ! Array indices

      I = 1
      J = 1
      K = 1

! **********************************************************************************************************************************
      IF      (WHICH == '11') THEN

         WRITE(F06,1101)
         WRITE(F06,1102)
         WRITE(F06,1103) ' 1'
         WRITE(F06,1104)

      ELSE IF (WHICH == '12') THEN
         WRITE(F06,1201)
            NAME1 = '      '
         DO II=1,NCORD
            IF (CORD(II,1) == 11) NAME1 = 'CORD1R'
            IF (CORD(II,1) == 12) NAME1 = 'CORD1C'
            IF (CORD(II,1) == 13) NAME1 = 'CORD1S'
            IF (CORD(II,1) == 21) NAME1 = 'CORD2R'
            IF (CORD(II,1) == 22) NAME1 = 'CORD2C'
            IF (CORD(II,1) == 23) NAME1 = 'CORD2S'
            IF       ((CORD(II,1) == 11) .OR. (CORD(II,1) == 12) .OR. (CORD(II,1) == 13)) THEN
               WRITE(F06,1202) NAME1, (CORD(II,JJ),JJ=2,5)
             ELSE IF ((CORD(II,1) == 21) .OR. (CORD(II,1) == 22) .OR. (CORD(II,1) == 23)) THEN
               WRITE(F06,1202) NAME1, (CORD(II,JJ),JJ=2,3)
            ENDIF
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,1203)
         DO II=1,NCORD
            NAME1(1:) = ' '   ;   NAME2(1:) = ' '
            IF (CORD(II,1) == 11) THEN   ;   NAME1 = 'CORD1R'   ;   NAME2 = ' from coords on GRID'   ;   ENDIF
            IF (CORD(II,1) == 12) THEN   ;   NAME1 = 'CORD1C'   ;   NAME2 = ' from coords on GRID'   ;   ENDIF
            IF (CORD(II,1) == 13) THEN   ;   NAME1 = 'CORD1S'   ;   NAME2 = ' from coords on GRID'   ;   ENDIF
            IF (CORD(II,1) == 21) THEN   ;   NAME1 = 'CORD2R'   ;   NAME2 = ' from coords on CORD'   ;   ENDIF
            IF (CORD(II,1) == 22) THEN   ;   NAME1 = 'CORD2C'   ;   NAME2 = ' from coords on CORD'   ;   ENDIF
            IF (CORD(II,1) == 23) THEN   ;   NAME1 = 'CORD2S'   ;   NAME2 = ' from coords on CORD'   ;   ENDIF
            WRITE(F06,1204) NAME1, CORD(II,2),'A',(RCORD(II,JJ),JJ=1,3), NAME2, GA(I)
            WRITE(F06,1205)                   'B',(RCORD(II,JJ),JJ=4,6), NAME2, GB(I)
            WRITE(F06,1205)                   'C',(RCORD(II,JJ),JJ=7,9), NAME2, GC(I)
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '21') THEN

         WRITE(F06,1102)
         WRITE(F06,1103) ' 2'
         WRITE(F06,2100)

         WRITE(F06,2101)
         DO II=1,NCORD
            NAME1 = '      '
            IF (CORD(II,1) == 11) NAME1 = 'CORD1R'
            IF (CORD(II,1) == 12) NAME1 = 'CORD1C'
            IF (CORD(II,1) == 13) NAME1 = 'CORD1S'
            IF (CORD(II,1) == 21) NAME1 = 'CORD2R'
            IF (CORD(II,1) == 22) NAME1 = 'CORD2C'
            IF (CORD(II,1) == 23) NAME1 = 'CORD2S'
            WRITE(F06,2102) NAME1, CORD(II,2),'A',(RCORD(II,JJ),JJ=1,3), CORD(II,3)
            WRITE(F06,2103)                   'B',(RCORD(II,JJ),JJ=4,6), CORD(II,3)
            WRITE(F06,2103)                   'C',(RCORD(II,JJ),JJ=7,9), CORD(II,3)
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '31') THEN

         WRITE(F06,1102)
         WRITE(F06,1103) ' 3'
         WRITE(F06,3101)
         WRITE(F06,3102)

      ELSE IF (WHICH == '32') THEN

         WRITE(F06,*)
         WRITE(F06,*)
         WRITE(F06,3201) CID
         WRITE(F06,3202) CID,RID

      ELSE IF (WHICH == '33') THEN

         WRITE(F06,3202) CIDK,RID

      ELSE IF (WHICH == '34') THEN

         WRITE(F06,*)
         WRITE(F06,*)
         WRITE(F06,3401)
         DO II=1,NCORD+1
            DO JJ=1,3*NCORD1+NCORD2
               IF      (RID_ARRAY(II,JJ) /= -99) THEN
                  WRITE(CRID_ARRAY(II,JJ),'(I9)') RID_ARRAY(II,JJ)
               ELSE
                  WRITE(CRID_ARRAY(II,JJ),'(A9)') '        '
               ENDIF
            ENDDO
         ENDDO

         DO II = 1,NCORD+1
            CRID_ARRAY_BLANK = 'Y'
            DO JJ=1,3*NCORD1+NCORD2
               IF (CRID_ARRAY(II,JJ)(1:) /= ' ') CRID_ARRAY_BLANK = 'N'
            ENDDO
            IF (CRID_ARRAY_BLANK == 'N') THEN
               IF (II == 1) THEN
                  WRITE(F06,3404)     (CRID_ARRAY(II,JJ),JJ=1,3*NCORD1+NCORD2)
               ELSE
                  WRITE(F06,3405) II-1, (CRID_ARRAY(II,JJ),JJ=1,3*NCORD1+NCORD2)
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '41') THEN

         WRITE(F06,1102)
         WRITE(F06,1103) ' 4'
         WRITE(F06,4101)
         WRITE(F06,4102)
         WRITE(F06,4110) (II,II=1,NCORD)
         WRITE(F06,4111) (CORD_NAME(II),II=1,NCORD)
         WRITE(F06,4112) (' --------',II=1,NCORD)
         WRITE(F06,4104) (CASCADE_PROC_ARRAY(1,II),II=1,NCORD)
         WRITE(F06,4105) (CASCADE_PROC_ARRAY(2,II),II=1,NCORD)
         WRITE(F06,4106) (CASCADE_PROC_ARRAY(3,II),II=1,NCORD)
         WRITE(F06,*)
         WRITE(F06,4107)
         WRITE(F06,4108)
         WRITE(F06,4109)
         WRITE(F06,*)

      ELSE IF (WHICH == '51') THEN
         WRITE(F06,1102)
         WRITE(F06,1103) ' 5'
         WRITE(F06,5101)
         WRITE(F06,5102) 'relative to it''s reference system. (only for CORD2 systems at this point)'

      ELSE IF (WHICH == '52') THEN
         WRITE(F06,5202) CID, 'it''s reference coordiate system', RID, (RCORD(I,JJ),JJ=1,3), CID, RID,                             &
                                   (RCORD(I,JJ),JJ=4,6), RID, CID, (RCORD(I,JJ),JJ=7,9), RID, CID, (RCORD(I,JJ),JJ=10,12), RID, CID

      ELSE IF (WHICH == '61') THEN
         WRITE(F06,1102)
         WRITE(F06,1103) ' 6'
         WRITE(F06,6100)
         WRITE(F06,5102) 'relative to the basic system. (only for CORD2 systems at this point)'

      ELSE IF (WHICH == '62') THEN
         WRITE(F06,5202) CORD(I,2), 'the basic system', CORD(I,3), (RCORD(I,JJ),JJ=1,3),                                           &
                                                        CORD(I,2), CORD(I,3), (RCORD(I,JJ),JJ=4,6),                                &
                                                        CORD(I,3), CORD(I,2), (RCORD(I,JJ),JJ=7,9),                                &
                                                        CORD(I,3), CORD(I,2), (RCORD(I,JJ),JJ=10,12), CORD(I,3), CORD(I,2)

      ELSE IF (WHICH == '63') THEN
         WRITE(F06,*)
         WRITE(F06,*) '  Test to see which coordinate transformation matrices have been cascaded back to basic (0) coord system'
         WRITE(F06,*) '  ----'
         WRITE(F06,6301) (CORD_NAME(II),II=1,NCORD)
         WRITE(F06,6302) (CORD(II,2),II=1,NCORD)
         WRITE(F06,6303) (TRANS_DONE(II),II=1,NCORD)
         WRITE(F06,*)
         WRITE(F06,*)

      ELSE IF (WHICH == '71') THEN
         WRITE(F06,1102)
         WRITE(F06,1103) ' 7'
         WRITE(F06,7101)

      ELSE IF (WHICH == '72') THEN
         WRITE(F06,7201) NUM_LEFT_AT_BEG
         DO II=1,NCORD
            IF (TRANS_DONE(II) == 'N') THEN
               IA = -99
               IB = -99
               IC = -99
               DO JJ=1,NCORD
                  IF (CORD(JJ,2) == CORD(II,3)) IA = JJ
                  IF (CORD(JJ,2) == CORD(II,4)) IB = JJ
                  IF (CORD(JJ,2) == CORD(II,5)) IC = JJ
               ENDDO
               IF ((IA >= 0) .AND. (IB >= 0) .AND. (IC >= 0)) THEN
                  WRITE(F06,7202) CORD_NAME(II), CORD(II,2), CORD(II,3), TRANS_DONE(IA), CORD(II,4), TRANS_DONE(IB),               &
                                                            CORD(II,5), TRANS_DONE(IC)
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '73') THEN
         WRITE(F06,7301) CORD_NAME(J), CID, GA(J), GB(J), GC(J)
         WRITE(F06,*)
         WRITE(F06,7302) CID, GA(J), CIDA, (RGA(J,K),K=1,3)
         WRITE(F06,7302) CID, GB(J), CIDB, (RGB(J,K),K=1,3)
         WRITE(F06,7302) CID, GC(J), CIDC, (RGC(J,K),K=1,3)
         WRITE(F06,*)

      ELSE IF (WHICH == '74') THEN
         WRITE(F06,*) '  Vectors used in the transformation of CORD1 systems:'
         WRITE(F06,*) '  -------'
         WRITE(F06,7402) 'V_0_CIDA', CIDA, (V_0_CIDA(JJ),JJ=1,3)
         WRITE(F06,7402) 'V_0_CIDB', CIDB, (V_0_CIDB(JJ),JJ=1,3)
         WRITE(F06,7402) 'V_0_CIDC', CIDC, (V_0_CIDC(JJ),JJ=1,3)
         WRITE(F06,*)

         WRITE(F06,7403) 'V_CIDB_CID_PTA', CIDA, 'A', CID, (V_CIDA_CID_PTA(JJ),JJ=1,3)
         WRITE(F06,7403) 'V_CIDB_CID_PTB', CIDB, 'B', CID, (V_CIDB_CID_PTB(JJ),JJ=1,3)
         WRITE(F06,7403) 'V_CIDB_CID_PTC', CIDC, 'C', CID, (V_CIDC_CID_PTC(JJ),JJ=1,3)
         WRITE(F06,*)

         WRITE(F06,7404) '     V0A            :in sys 0 from origin of sys        0 to pt A on coord sys ', CID, (V0A(JJ),JJ=1,3)
         WRITE(F06,7404) '     V0B            :in sys 0 from origin of sys        0 to pt B on coord sys ', CID, (V0B(JJ),JJ=1,3)
         WRITE(F06,7404) '     V0C            :in sys 0 from origin of sys        0 to pt C on coord sys ', CID, (V0C(JJ),JJ=1,3)
         WRITE(F06,*)

         WRITE(F06,7405) '     VK             :in sys 0 from point A to point B on coord sys             ', CID, (VK(JJ) ,JJ=1,3)
         WRITE(F06,7405) '     VI             :in sys 0 from point A to point C on coord sys             ', CID, (VI(JJ) ,JJ=1,3)
         WRITE(F06,*)

         WRITE(F06,7406) '     IVEC           :in sys 0 in x direction for coord sys                     ', CID, (IVEC(JJ) ,JJ=1,3)
         WRITE(F06,7406) '     JVEC           :in sys 0 in y direction for coord sys                     ', CID, (JVEC(JJ) ,JJ=1,3)
         WRITE(F06,7406) '     KVEC           :in sys 0 in z direction for coord sys                     ', CID, (KVEC(JJ) ,JJ=1,3)
         WRITE(F06,*)

      ELSE IF (WHICH == '75') THEN
         WRITE(F06,5202) CORD(J,2), 'the basic system',0, (RCORD(J,JJ),JJ= 1, 3), CORD(J,2), 0,                                    &
                                                          (RCORD(J,JJ),JJ= 4, 6), 0, CORD(J,2),                                    &
                                                          (RCORD(J,JJ),JJ= 7, 9), 0, CORD(J,2),                                    &
                                                          (RCORD(J,JJ),JJ=10,12), 0, CORD(J,2)
         WRITE(F06,7501) CORD_NAME(J),CID

      ELSE IF (WHICH == '81') THEN
         WRITE(F06,1102)
         WRITE(F06,1103) ' 8'
         WRITE(F06,8101)

         WRITE(F06,*) 'Test to see which coordinate transformation matrices have been cascaded back to basic (0) coord system'
         WRITE(F06,*) '----'
         WRITE(F06,6301) (CORD_NAME(II),II=1,NCORD)
         WRITE(F06,6302) (CORD(II,2),II=1,NCORD)
         WRITE(F06,6303) (TRANS_DONE(II),II=1,NCORD)
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,8104)
         WRITE(F06,8105)
         DO II=1,NCORD
            CORD_TYPE = CORD(II,1)
            IF ((CORD_TYPE == 11) .OR. (CORD_TYPE == 11) .OR. (CORD_TYPE == 11)) THEN
               WRITE(F06,8106) CORD_NAME(II), (CORD(II,JJ),JJ=2,5)
            ELSE
               WRITE(F06,8106) CORD_NAME(II), (CORD(II,JJ),JJ=2,3)
            ENDIF
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*)

         WRITE(F06,8107)
         WRITE(F06,8108)
         DO II=1,NCORD
            WRITE(F06,8109) ' CID ', CORD(II,2),' |',(RCORD(II,JJ),JJ= 1, 3),' |', (RCORD(II,JJ),JJ= 4, 6),' |'                    &
                                                    ,(RCORD(II,JJ),JJ= 7, 9),' |', (RCORD(II,JJ),JJ=10,12),' |'
         ENDDO
         WRITE(F06,*)

 8107 FORMAT(42X,'A R R A Y   R C O R D   W I T H   O R I G I N S   A N D   3 X 3   T R A N F O R M A T I O N   M A T R I C E S',/,&
             70X,'T H A T   T R A N S F O R M   C I D   T O   B A S I C',/)

 8108 FORMAT(14X,'|         Coords of origin of CID        |  Row 1 of coord transformation matrix  |',                            &
                  '  Row 2 of coord transformation matrix  |  Row 3 of coord transformation matrix')

 8109 FORMAT(A,I8,A,3ES13.6,A,3ES13.6,A,3ES13.6,A,3ES13.6,A)
         WRITE(F06,8901)

      ENDIF

! **********************************************************************************************************************************
 1101 FORMAT(' ___________________________________________________________________________________________________________________'&
            ,'________________'                                                                                                ,//,&
             ' ::::::::::::::::::::::::::::::::::::::START PARAM PRTCORD OUTPUT FROM SUBROUTINE CORD_PROC::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1102 FORMAT(1X,'****************************************************************************************************************',&
                '*******************')

 1103 FORMAT(58X,'P H A S E   ',A,/,58X,'--------------',/)

 1104 FORMAT(' Phase 1: For CORD1R: (1) Change grid numbers (that were put into array CORD when B.D. was read) to the reference',  &
                                      ' coordinate system'                                                                      ,/,&
             '                          number for that grid.'                                                                  ,/,&
             '                      (2) Put coords of the reference grids on the CORD1R into array RCORD.',//)

 1201 FORMAT(' Array CORD after replacing the 3 grid numbers for CORD1 entries with the CORD1,2 reference system for that grid' ,/,&
             ' ---------- 1st col is the coord type (CORD1R,C,S or CORD2R,C,S).'                                                ,/,&
             '            2nd col is the coordinate system number (CID).'                                                       ,/,&
             '            3rd col and on give the ref coord system numbers (RID) for the CID. CORD2 types have 1 RID, CORD1',      &
             ' types have 3 RID''s',//,42X,'Coord type      CID      RID''s ->',/)

 1202 FORMAT('                                            ',A,2X,4I9)

 1203 FORMAT(' Array RCORD (cols 1-3) with coords of the 3 reference points from CORD2 entries and coords of the 3 grids on',      &
             ' CORD1 entries',/,&
             ' -----------',/)

 1204 FORMAT('     For ',A,' system ',I8,' the location of pt '   ,A,' is:',3(1ES14.6),A,I8)

 1205 FORMAT('                                                   ',A,' is:',3(1ES14.6),A,I8)

 2100 FORMAT(' Phase 2: Change coords on cylindrical and spherical systems to coords in terms of the defining rectangular axes',   &
             ' of that system',//)

 2101 FORMAT(' Array RCORD (cols 1-3) with coords of the 3 reference points in terms of the defining rectangular axes of the',     &
             ' reference system',/,' -----------',/)

 2102 FORMAT('     For ',A,' system ',I8,' the location of pt '   ,A,' is:',3(1ES14.6),' in reference system ',I8)

 2103 FORMAT('                                                   ',A,' is:',3(1ES14.6),' in reference system ',I8)

 3101 FORMAT(' Phase 3: Check coordinate system data for logic:'                                                                ,/,&
             '          (1) There must be a final reference to the basic system (0)'                                            ,/,&
             '          (2) There can be no circular references (i.e. a coord system can reference any number of other coord sys', &
                          ' but cannot'                                                                                         ,/,&
             '              reference itself in that chain of refs).'                                                           ,/,&
             '          (3) Create RID_ARRAY which shows the chain of references of each coord system',//)

 3102 FORMAT(' Trace of coordinate chain of references to basic system:'                                                        ,/,&
             ' --------------------------------------------------------')

 3201 FORMAT('   References for coordinate system number ',I8                                                                    ,/&
             '   ------------------------------------------------')

 3202 FORMAT('     Coordinate system ',I8,' references coordinate system ',I8)

 3401 FORMAT('   RID_ARRAY: The cols below show the trace of a coord sys (number CID) through it''s reference systems',            &
             ' (RID''s in the cols below CID)',/,                                                                                  &
             '   ---------  to the basic system (0). ICID is the internal coord system ID for actual coord system CID.',           &
             ' Each CORD1C,R,S system will',/,                                                                                     &
             '              have 3 cols (one for each reference point on the CORD1C,R,S BDF entry).',          &
             ' Only 1 col for each CID will be used in',/,                                                                         &
             '              cascading the RID''s back to CID to get the coord trnasformations from basic to CID.',/)

 3404 FORMAT('     CID     ->',32767A9)

 3405 FORMAT('     RID',I4,' ->',32767A9)

 4101 FORMAT(' Phase 4: Solve for array CASCADE_PROC_ARRAY. It will be used to do the cascading of coord references to get the',   &
                     ' transformation from'                                                                                     ,/,&
             '         CID to basic',//)

 4102 FORMAT('   CASCADE_PROC_ARRAY array: (info from above RID_ARRAY used in the cascading of RID''s from CID to basic(0))',/,    &
             '   -------------------')

 4104 FORMAT('     Row 1:',32767I9)

 4105 FORMAT('     Row 2:',32767I9)

 4106 FORMAT('     Row 3:',32767I9)

 4107 FORMAT('     Row 1 of CASCADE_PROC_ARRAY is the CID coord sys number from RID_ARRAY')

 4108 FORMAT('     Row 2 of CASCADE_PROC_ARRAY is the col number in RID_ARRAY to use when solving for the coord transformation',   &
             ' from basic to CID')

 4109 FORMAT('     Row 3 of CASCADE_PROC_ARRAY is the num of coord system refs above the basic sys (0) in RID_ARRAY',              &
             ' (cascade from basic up to CID)')

 4110 FORMAT('     ICID  ',32767I9)

 4111 FORMAT('     Type  ',32767(3X,A6))

 4112 FORMAT('           ',32767A9)

 5101 FORMAT(' Phase 5: Get the unit vectors in each of the CORD2 coord systems in terms of their reference system. Put the',      &
                      ' results into array'                                                                                     ,/,&
             '          RCORD where the transformation matrices go. These are NOT the final transformations; they are the',        &
                      ' transformations from a'                                                                                 ,/,&
             '          CID to its reference system and ONLY for CORD2 systems (CORD1 systams will be handled later).',//)

 5102 FORMAT(' RCORD array with data for a coord system ',A,/,                                                                     &
             ' -----------',/)

 5202 FORMAT ('   RCORD array for coordinate system CID = ',I9,' with values relative to ',A,                                      &
              ' (RID) = ',I9   ,/,&
              '   -----------'                                                                                                  ,/,&
              '     Cols  1- 3: ',3(1ES13.6),': Origin of coord system ',I8,' as measured in coord system ',I8,/                ,/,&
              '     Cols  4- 6: ',3(1ES13.6),': Row 1 of matrix that transforms a vector to ',I9,' from ',I9                    ,/,&
              '     Cols  7- 9: ',3(1ES13.6),': Row 2 of matrix that transforms a vector to ',I9,' from ',I9                    ,/,&
              '     Cols 10-12: ',3(1ES13.6),': Row 3 of matrix that transforms a vector to ',I9,' from ',I9,/)

 6100 FORMAT(' Phase 6:  Calculate transformations to basic (zero ) system for all CORD2 coordinate systems.',//)

 6101 FORMAT(' Coord transformations between principal directions of CORD2 basic coord systems and all other coord systems:'    ,/,&
             ' -----------------------------------------------------------------------------------------------------------',/)

 6102 FORMAT('   The following matrix will transform a vector to coord system ',I8,' from a vector in coord system ',I8          ,/&
             '   ------------------------------------------------------------------------------------------------------------')

 6103 FORMAT(3(1ES20.8))

 6104 FORMAT(' Basic coords of the origin of CORD2 coord system number ',I8,' are : ',3(1ES16.8))

 6301 FORMAT('   Coordinate system type         :',32767(3X,A6))

 6302 FORMAT('   Coordinate system ID           :',32767I9)

 6303 FORMAT('   Is transformation to basic done?',32767A9)

 7201 FORMAT(' There are',I8,' CORD1 entries remaining to be processed at this time. The table below shows which systems need to', &
             ' be processed,',/,' the reference systems they use, and whether the reference systems have been transformed to basic'&
            ,' (Y or N).',/,' In order to process a CORD1 system all 3 of the ref systems myst have been transformed to basic',/)

 7101 FORMAT(' Phase 7: Process CORD1 entries to get transformations to basic and put results into array RCORD',//)

 7202 FORMAT(3X,A,2X,I8,' has 3 ref systems (Y/N means is it transformed to basic?):',3(I8,' (',A,')'))

 7301 FORMAT(' Working on ',A,2X,I9,' with 3 reference grids: ',3I9,/,' ----------------------------')

 7302 FORMAT('   For coord sys ',I9,' the coords of grid ',I9,' in ref system ',I9,' are:     ',3(1ES14.6))

 7402 FORMAT(5X,A,7X,':in sys 0 from origin of sys        0 to origin of  sys    ', I8,' =',3(1ES14.6))

 7403 FORMAT(5X,A,1X,':in sys 0 from origin of sys ',I8,' to point ',A,' of sys    '   ,I8,' =',3(1ES14.6))

 7404 FORMAT(A,I8,' =', 3(1ES14.6))

 7405 FORMAT(A,I8,' =', 3(1ES14.6))

 7406 FORMAT(A,I8,' =', 3(1ES14.6))

 7501 FORMAT('   Transformation to basic is completed for ',A,2X,I8                                                             ,/,&
             '   =========================================================',//)
 8101 FORMAT(' Phase 8: Check that all systems have been transformed to basic and rewrite RID''s in array CORD to reflect this',//)

 8104 FORMAT(' Arrays CORD and RCORD at end of subr CORD_PROC',/,' ----------------------------------------------',/)

 8105 FORMAT(66X,'A R R A Y   C O R D',/,61X,'(all RID''s should be zero now)',//,60X,'Coord type      CID      RID''s ->',/)

 8106 FORMAT(62X,A,2X,5I9)

 8901 FORMAT(' :::::::::::::::::::::::::::::::::::::::END PARAM PRTCORD OUTPUT FROM SUBROUTINE CORD_PROC:::::::::::::::::::::::::' &
             ,':::::::::::::::::'                                                                                               ,/,&
             ' ___________________________________________________________________________________________________________________'&
            ,'________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE PARAM_PRTCORD_OUTPUT

      END SUBROUTINE CORD_PROC


      SUBROUTINE GRID_PROC

! Performs 7 functions:
!   1) Generate GRID_ID and GRID_SEQ and sort them so that GRID_ID is in numerical order and then GRID_SEQ(I)
!      is the position, in the stack of GRID's in Bulk Data, where GRID_ID(I) exists. GRID_SEQ may change order
!      in subroutine SEQ_PROC depending on the grid point sequencing scheme the user asks for.
!   2) Sort arrays GRID and RGRID so that they are in grid point numerical order and check for duplicate grid ID's.
!   3) Reset coord sys and perm SPC data on GRID cards based on values read from a GRDSET card (if in the data deck)
!   4) Call CORD_PROC.FOR to calc the transformations to basic from each of the coord systems (their principal axes)
!   5) Transform the grid coordinates to the basic coord system
!   6) Write grid data to filename.L1B
!   7) Write some grid data to output file if requested

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  CONV_DEG_RAD
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1B, OP2, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, FATAL_ERR, MCORD, MRCORD, MGRID, MRGRID, NCORD, NGRID
      USE PARAMS, ONLY                :  PRTBASIC
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID, RGRID, GRID_ID, GRID_SEQ, CORD, RCORD, TN

      USE SORTING, ONLY               :  SORT_INT2, SORT_GRID_RGRID
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_GRID_COORDS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GRID_PROC'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes

      INTEGER(LONG)                   :: CP                ! The actual coord sys that a grid is located in.
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: JCORD             ! Internal coord sys ID
      INTEGER(LONG)                   :: JFLD              ! Used in error message to indicate a coord sys ID undefined


      REAL(DOUBLE)                    :: ANG1              ! An angle in a cyl or sph coord sys from the RGRID array
      REAL(DOUBLE)                    :: ANG2              ! An angle in a cyl or sph coord sys from the RGRID array
      REAL(DOUBLE)                    :: RADIUS            ! A radius coord of a grid point in a cyl or sph coord sys (CP)
      REAL(DOUBLE)                    :: RGRIDI_CP(3)      ! Coords of one grid point in the principal rectangular axes of coord
!                                                            sys CP rel to the origin of coord sys CP
      REAL(DOUBLE)                    :: RGRIDI_0(3)       ! Coords of one grid point in basic coords rel to origin of coord sys CP
!                                                            sys CP rel to the origin of coord sys CP

      INTRINSIC                       :: DCOS, DSIN



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

! Part 1:
! -------

! Generate GRID_ID and GRID_SEQ data.
      WRITE(SC1,12345,ADVANCE='NO') '    Initialize arrays GRID and GRID_SEQ                                     ', CR13
      DO I=1,NGRID
         GRID_ID(I)  = GRID(I,1)
         GRID_SEQ(I) = I                                   ! This is the initial GRID_SEQ array (order of GRID input)
      ENDDO                                                ! It may change in subr SEQ_PROC

! Sort these so that the actual G.P. numbers in GRID_ID are in numerical order. ! At this point GRID_SEQ(I) is seq num for internal
! grid I but this will change (in subr SEQ_PROC) unless the final sequence order is to be the order of the grids as read in from the
! input data deck.

      WRITE(SC1,12345,ADVANCE='NO') '    Sort arrays GRID and GRID_SEQ so GRID is in numerical order', CR13
      IF (NGRID > 1) THEN
         CALL SORT_INT2 ( SUBR_NAME, 'GRID_ID, GRID_SEQ', NGRID, GRID_ID, GRID_SEQ )
      ENDIF

! **********************************************************************************************************************************
! Part 2:
! -------

! Sort arrays GRID and RGRID so that they are in grid point numerical order.

      WRITE(SC1,12345,ADVANCE='NO') '    Sort arrays GRID and RGRID so GRID is in numerical order                ', CR13
      CALL SORT_GRID_RGRID ( SUBR_NAME, 'GRID, RGRID', NGRID, GRID, RGRID )

! Check for duplicate GRID ID's and quit if there are any

      WRITE(SC1,12345,ADVANCE='NO') '    Check for duplicate GRID IDs                                            ', CR13
      IERROR = 0
      DO I=1,NGRID-1
         IF (GRID(I+1,1) == GRID(I,1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1344) GRID(I+1,1)
            WRITE(F06,1344) GRID(I+1,1)
         ENDIF
      ENDDO

! Check: If GRID(I,1) /= GRID_ID(I), then coding error

      WRITE(SC1,12345,ADVANCE='NO') '    Check GRID_ID array                                                     ', CR13
      DO I=1,NGRID
         IF (GRID(I,1) == GRID_ID(I)) THEN
            CYCLE
         ELSE
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1346) SUBR_NAME,GRID(I,1),I,GRID_ID(I)
            WRITE(F06,1346) SUBR_NAME,GRID(I,1),I,GRID_ID(I)
         ENDIF
      ENDDO

      IF (IERROR > 0) THEN
         WRITE(ERR,1343) IERROR
         WRITE(F06,1343) IERROR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Part 3:
! -------

! Check to make sure that all coord systems referenced in field 3 (input coord sys) of the grid card except the 0 system exist.

      WRITE(SC1,12345,ADVANCE='NO') '    Check that coord systems in field 3 referenced do exist                 ', CR13
      IERROR = 0
i_do1:DO I=1,NGRID
         CP = GRID(I,2)
         IF (CP /= 0) THEN
            DO J=1,NCORD
               IF (CORD(J,2) == CP) THEN
                  CYCLE i_do1
               ENDIF
            ENDDO
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            JFLD   = 3
            WRITE(ERR,910) GRID(I,2),JFLD,GRID(I,1),JFLD
            WRITE(F06,910) GRID(I,2),JFLD,GRID(I,1),JFLD
         ENDIF
      ENDDO i_do1

! Also check that any coord sys in field 7 (global) are defined. They won't be used at this time, but we want to know they exist.

      WRITE(SC1,12345,ADVANCE='NO') '    Check that coord systems in field 7 referenced do exist                 ', CR13
i_do2:DO I=1,NGRID
         CP = GRID(I,3)
         IF (CP /= 0) THEN
            DO J=1,NCORD
               IF (CORD(J,2) == CP) THEN
                  CYCLE i_do2
               ENDIF
            ENDDO
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            JFLD   = 7
            WRITE(ERR,910) GRID(I,3),JFLD,GRID(I,1),JFLD
            WRITE(F06,910) GRID(I,3),JFLD,GRID(I,1),JFLD
         ENDIF
      ENDDO i_do2
      IF (IERROR > 0) THEN
         WRITE(ERR,1343) IERROR
         WRITE(F06,1343) IERROR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Part 4:
! -------

! Calculate coordinate system transformation matrices if there are any CORD cards

      IF (NCORD /= 0) THEN

         WRITE(SC1,12345,ADVANCE='NO') '    Calc coord sys transformation matrices                                  ', CR13
         CALL CORD_PROC

! **********************************************************************************************************************************
! Part 5:
! -------

! Transform grid coordinates to basic system if the coords were not input in basic (i.e. if GRID(I,2) not 0).

         WRITE(SC1,12345,ADVANCE='NO') '    Transform grid coords to basic system                                   ', CR13
         JCORD = 0
grid_do: DO I=1,NGRID
            CP = GRID(I,2)
            IF (CP == 0) CYCLE grid_do
            DO J = 1,NCORD
               IF (CORD(J,2) == CP) THEN
                  JCORD = J
                  EXIT
               ENDIF
            ENDDO
                                                           ! JCORD=0, means pgm error
            IF (JCORD == 0) THEN                           ! (we should not have gotten here if a coord sys is undefined)
               WRITE(ERR,1345) SUBR_NAME
               WRITE(F06,1345) SUBR_NAME
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
                                                           ! Cylindrical.  The Z in RGRID(I,3) is OK as is.
            IF      ((CORD(JCORD,1) == 12) .OR. (CORD(JCORD,1) == 22)) THEN
               RADIUS     = RGRID(I,1)
               ANG1       = RGRID(I,2)*CONV_DEG_RAD
               RGRID(I,1) = RADIUS*DCOS(ANG1)
               RGRID(I,2) = RADIUS*DSIN(ANG1)
                                                           ! Spherical
            ELSE IF ((CORD(JCORD,1) == 13) .OR. (CORD(JCORD,1) == 23)) THEN
               RADIUS     = RGRID(I,1)
               ANG1       = RGRID(I,2)*CONV_DEG_RAD
               ANG2       = RGRID(I,3)*CONV_DEG_RAD
               RGRID(I,1) = RADIUS*DSIN(ANG1)*DCOS(ANG2)
               RGRID(I,2) = RADIUS*DSIN(ANG1)*DSIN(ANG2)
               RGRID(I,3) = RADIUS*DCOS(ANG1)
            ENDIF

            DO J=1,3                                       ! RGRID(I,1-3) has rectangular coords, in the princ axes of JCORD,
               RGRIDI_CP(J) = RGRID(I,J)                   ! relative to the origin of that coord system
            ENDDO

! The basic coords of G.P. I are now obtained by multiplying the coord system transformation matrix TN (see CORD_PROC.FOR)
! times the rectangular coords in RGRIDI_CP (to get RGRID_0) and adding the basic coords of the origin of coord system
! JCORD (otherwise the calculated basic coords of G.P. I would only be relative to the JCORD origin).  The coords in
! basic of JCORD are in RCORD array (calculated in CORDP.FOR).

            CALL MATMULT_FFF ( TN(1,1,JCORD), RGRIDI_CP, 3, 3, 1, RGRIDI_0 )

            DO J=1,3                                       ! Change RGRID(I,1-3) to have the basic coords of G.P. I
               RGRID(I,J) = RGRIDI_0(J) + RCORD(JCORD,J)
            ENDDO

            GRID(I,2) = CORD(JCORD,3)                      ! CORD(JCORD,3) was set = 0 in subr CORD_PROC after all transformations
!                                                            so set GRID(I,2) to 0 (basic) also
         ENDDO grid_do

      ENDIF

! **********************************************************************************************************************************
! Part 6:

! Store the grid point data in 'L1B'

      WRITE(SC1,12345,ADVANCE='NO') '    Write grid data to file                                                 ', CR13
      DATA_SET_NAME = 'GRID, RGRID'
      WRITE(L1B) DATA_SET_NAME
!      WRITE(OP2) DATA_SET_NAME
      WRITE(L1B) NGRID
!      WRITE(OP2) NGRID
      DO I=1,NGRID
         DO J=1,MGRID
            WRITE(L1B) GRID(I,J)
!            WRITE(OP2) GRID(I,J)
         ENDDO
         DO J=1,MRGRID
            WRITE(L1B) RGRID(I,J)
!            WRITE(OP2) RGRID(I,J)
         ENDDO
      ENDDO

! Write coord data to L1B. Need coord type, CID, 12 words of RCORD. The RID for all systems is now 0

      WRITE(SC1,12345,ADVANCE='NO') '    Write coord data to file                                                ', CR13
      DATA_SET_NAME = 'COORDINATE SYSTEM DATA'
      WRITE(L1B) DATA_SET_NAME
      WRITE(L1B) NCORD
      DO I=1,NCORD
         DO J=1,MCORD
            WRITE(L1B) CORD(I,J)
         ENDDO
         DO J=1,MRCORD
            WRITE(L1B) RCORD(I,J)
         ENDDO
      ENDDO

      WRITE(SC1,*) CR13

! **********************************************************************************************************************************
! Part 7:

! Print out grid data in basic coords

      IF (PRTBASIC == 1) THEN
         CALL WRITE_GRID_COORDS
      ENDIF



      RETURN

! **********************************************************************************************************************************
  910 FORMAT(' *ERROR   910: COORD SYSTEM',I8,' IN FIELD ',I2,' ON GRID ',I8,' (OR FROM GRDSET ENTRY FIELD ',I2,') IS UNDEFINED')

 1343 FORMAT(' PROCESSING TERMINATED DUE TO ',I8,' ERRORS')

 1344 FORMAT(' *ERROR  1344: GRID NUMBER ',I8,' IS A DUPLICATE GRID NUMBER')

 1345 FORMAT(' *ERROR  1345: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' JCORD MUST NOT BE ZERO')

 1346 FORMAT(' *ERROR  1346: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAYS GRID AND GRID_ID SHOULD HAVE THE SAME GRID IDs, BUT THEY DO NOT.'                              &
                    ,/,14X,' GRID NO. ',I8,' AT ROW ',I8,' IS NOT THE SAME AS GRID NO ',I8,' AT THE SAME ROW ')

12345 FORMAT(A, A)

! **********************************************************************************************************************************

      END SUBROUTINE GRID_PROC


      SUBROUTINE SEQ_PROC

! Generates the grid point sequence order.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     SEQ,     L1B
      USE IOUNT1, ONLY                :  WRT_ERR, SEQFIL
      USE IOUNT1, ONLY                :  WRT_ERR, SEQSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, FATAL_ERR, NGRID, NSEQ, PROG_NAME, WARN_ERR
      USE PARAMS, ONLY                :  EPSIL, GRIDSEQ
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  GRID_ID, GRID_SEQ, INV_GRID_SEQ, SEQ1, SEQ2
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE SORTING, ONLY               :  SORT_INT1, SORT_INT1_REAL1, SORT_INT2, SORT_INT2_REAL1, SORT_REAL1_INT1
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, READERR, STMERR
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKCARD, MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  LEFT_ADJ_BDFLD
      USE GRID_COORDINATES, ONLY      :  BD_SEQGP
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_INTEGER_VEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SEQ_PROC'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: TMP_GRID_ID(NGRID)! Set to array GRID_ID for aid in sorting GRID_SEQ
      INTEGER(LONG)                   :: TMP_GRD_SEQ(NGRID)! Set to array GRID_SEQ so we can sort it and get array INV_GRID_SEQ
!                                                            without disturbing GRID_SEQ sequence


      REAL(DOUBLE)                    :: R_GSEQ(NGRID)     ! Real sequence numbers (since SEQGP cards can have real no's). In the
!                                                            end, the sequence array that will be used is integer array GRID_SEQ

      INTRINSIC                       :: DBLE



! **********************************************************************************************************************************
! Coming in to this subr, GRID_SEQ is in the order of the grids as read in the input data deck.

! Generate initial R_GSEQ based on the GRID_SEQ value. R_GSEQ(I) is  the (real) sequence number for Grid Point GRID_ID(I).
! If there are no SEQGP sequencing cards, then this will be the final grid point sequence order (as a real number).
! It will be converted to an array of consecutive integers in GRID_SEQ later.

! NOTE: when Bandit sequencing is requested we use Bandit to generate a set of SEQGP cards (if it runs correctly and finds
!       resequencing is necessary) and then proceed with that complete set of SEQGP cards. Need to do it this way because we have to
!       get INV_GRID_SEQ (later) and to write the seq arrays to L1B.

      IF       (GRIDSEQ(1:6) == 'BANDIT') THEN             ! Call subr AUTO_SEQ_PROC to generate SEQ1, SEQ2 from SEQGP card images
         CALL AUTO_SEQ_PROC
         IF (NSEQ == NGRID) THEN                           ! Bandit did reseq grids. Set R_GSEQ to the SEQ2 from Bandit SEQGP cards
            DO I=1,NGRID
!              R_GSEQ(I) = DBLE(SEQ2(I))                   ! Shouldn't need this, SEQ2 is REAL(DOUBLE)
               R_GSEQ(I) = SEQ2(I)
            ENDDO
         ELSE                                              ! AUTO_SEQ_PROC didn't reseq all grids so reset GRIDSEQ = 'GRID'
            WRITE(ERR,101) NGRID,NSEQ,PROG_NAME
            IF (SUPINFO == 'N') THEN
               WRITE(F06,101) NGRID,NSEQ,PROG_NAME
            ENDIF
            GRIDSEQ = 'GRID    '                           ! Need this to cover case where AUTO_SEQ_PROC returned without completing
         ENDIF
      ENDIF

      IF (GRIDSEQ(1:4) == 'GRID'  ) THEN                   ! Sequence grids in numerical order, but include SEQGP entries (later)
         DO I=1,NGRID
            R_GSEQ(I) = DBLE(I)
         ENDDO
      ELSE IF (GRIDSEQ(1:5) == 'INPUT')  THEN              ! Sequence grids in input order, but include SEQGP entries (later)
         DO I=1,NGRID
           R_GSEQ(I) = DBLE(GRID_SEQ(I))
         ENDDO
      ENDIF

! Check to make sure that all grid points on SEQGP cards are defined

      IERROR = 0
      DO I = 1,NSEQ
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SEQ1(I), IGRID )
         IF (IGRID == -1) THEN
            WRITE(ERR,1361) 'GRID', SEQ1(I), 'SEQGP BULK DATA ENTRY'
            WRITE(F06,1361) 'GRID', SEQ1(I), 'SEQGP BULK DATA ENTRY'
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO
      IF (IERROR > 0 ) THEN
         WRITE(ERR,9999) SUBR_NAME, IERROR
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9999) SUBR_NAME, IERROR
         ENDIF
         CALL OUTA_HERE ( 'Y' )                            ! Some grid ID's on SEQGP cards not defined
      ENDIF

! Reset GRID_SEQ to consecutive integers. R_GSEQ will be the grid sequence number for the time being. At the end,
! the integer sequence numbers will be back in GRID_SEQ

      DO I=1,NGRID
         GRID_SEQ(I) = I
      ENDDO

! If there are SEQGP cards in the data deck, or if auto sequencing has produced them, then arrays SEQ1, SEQ2 exist.
! Sort SEQ1 (grid numbers) with SEQ2 (sequence numbers) so that SEQ1 is in grid point numerical order (like GRID_ID) and then
! check to make sure that no sequence numbers are duplicated.

      IF (NSEQ > 1) THEN

         CALL SORT_INT1_REAL1 ( SUBR_NAME, 'SEQ1, SEQ2', NSEQ, SEQ1, SEQ2 )

         DO I = 1,NSEQ-1
            IF (SEQ1(I) == SEQ1(I+1)) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1399) SEQ1(I+1),SEQ2(I+1)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,1399) SEQ1(I+1),SEQ2(I+1)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

! Change sequence numbers in R_GSEQ to values in SEQ2 (from SEQGP Bulk Data cards). Note: SEQ1 is in grid numerical order
! so R_GSEQ(1) will be for grid GRID_ID(1) and so on.

      IF (NSEQ == NGRID) THEN
         DO I=1,NGRID
            R_GSEQ(I) = SEQ2(I)
         ENDDO
      ELSE
         DO I=1,NGRID
            DO J=1,NSEQ
               IF (SEQ1(J) == GRID_ID(I)) THEN
                  R_GSEQ(I) = SEQ2(J)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! Set TMP_GRID_ID to GRID_ID and then sort TMP_GRID_ID with R_GSEQ so that R_GSEQ is in numerically increasing order.

      DO I=1,NGRID
         TMP_GRID_ID(I) = GRID_ID(I)
      ENDDO


      CALL SORT_REAL1_INT1 ( SUBR_NAME, 'R_GSEQ, TMP_GRID_ID', NGRID, R_GSEQ, TMP_GRID_ID )

! Check to make sure that there are no redundant R_GSEQ sequence numbers


      IERROR = 0
      DO I=1,NGRID-1
         IF ((DABS(R_GSEQ(I+1)) - DABS(R_GSEQ(I))) < EPSIL(1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1360) TMP_GRID_ID(I), TMP_GRID_ID(I+1), R_GSEQ(I)
            WRITE(F06,1360) TMP_GRID_ID(I), TMP_GRID_ID(I+1), R_GSEQ(I)
         ENDIF
      ENDDO

      IF (IERROR > 0) THEN
         WRITE(ERR,9999) SUBR_NAME, IERROR
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9999) SUBR_NAME, IERROR
         ENDIF
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Sort TMP_GRID_ID with R_GSEQ and with GRID_SEQ so that TMP_GRID_ID is in numerically increasing order.

      CALL SORT_INT2_REAL1 ( SUBR_NAME, 'TMP_GRID_ID, GRID_SEQ, R_GSEQ', NGRID, TMP_GRID_ID, GRID_SEQ, R_GSEQ )

! Generate array INV_GRID_SEQ

      DO I=1,NGRID
         TMP_GRD_SEQ(I)  = GRID_SEQ(I)
         INV_GRID_SEQ(I) = I
      ENDDO

      CALL SORT_INT2 ( SUBR_NAME, 'TMP_GRD_SEQ, INV_GRID_SEQ', NGRID, TMP_GRD_SEQ, INV_GRID_SEQ )

! Print table showing R_GSEQ, GRID_SEQ and INV_GRID_SEQ

      IF (DEBUG(13) == 1) THEN
         WRITE(F06,111)
         DO I=1,NGRID
            WRITE(F06,112) GRID_ID(I), I, R_GSEQ(I), GRID_SEQ(I), INV_GRID_SEQ(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,113)
         WRITE(F06,*)
      ENDIF

! **********************************************************************************************************************************
! Now write sequence and GRID_SEQ data to file LINK1B.

      DATA_SET_NAME = 'GRID_SEQ, INV_GRID_SEQ'
      WRITE(L1B) DATA_SET_NAME
      WRITE(L1B) NGRID
      DO I=1,NGRID
         WRITE(L1B) GRID_SEQ(I), INV_GRID_SEQ(I)
      ENDDO

      DATA_SET_NAME = 'SEQ1, SEQ2'
      WRITE(L1B) DATA_SET_NAME
      WRITE(L1B) NSEQ
      DO I=1,NSEQ
         WRITE(L1B) SEQ1(I),SEQ2(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: SUBR AUTO_SEQ_PROC DID NOT SEQUENCE ALL OF THE ',I8,' GRIDS. ONLY ',I8,' GRIDS WERE SEQUENCED.'       &
                  ,/,15X,A,' WILL DEFAULT TO A SEQUENCE THAT IS IN GRID NUMERICAL ORDER',/)

  111 FORMAT(56X,'GRID SEQUENCE DATA',//16X,'GRID ID                  I                   R_GSEQ(I)               GRID_SEQ(I)   ', &
                 '     INV_GRID_SEQ(I)',/,12X,'(Actual grid ID)    (Internal grid ID)     (Grid seq - real num)',                  &
                 '  (Grid seq - integer num)  (* - see below)',/)

  112 FORMAT(3X,2I20,1ES28.6,2I20)

  113 FORMAT(15X,'* INV_GRID_SEQ(I) = internal grid ID that is sequenced I-th')

 1399 FORMAT(' *WARNING    : REDUNDANT VALUE IN G.P. SEQUENCE ARRAY.'                                                              &
                    ,/,14X,' GRID POINT ',I8,' WILL USE SEQUENCE NUMBER ',1ES13.6)

 1360 FORMAT(' *ERROR  1360: SEQUENCE NUMBERS FOR GRIDS ',I8,' AND ',I8,' ARE BOTH ',1ES13.6,'. SEQUENCE NUMBERS MUST BE UNIQUE')


 1361 FORMAT(' *ERROR  1361: UNDEFINED ',A,I8,' ON ',A)

 9999 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE AUTO_SEQ_PROC

! Reads SEQGP card images from bandit output file (filename.SEQ) using subr BD_SEQGP which creates SEQ1, SEQ2 arrays from SEQGP info

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, sc1
      USE SCONTR, ONLY                :  BANDIT_ERR, BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, LSEQ, NGRID, NSEQ,          &
                                         PROG_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE PARAMS, ONLY                :  GRIDSEQ, SEQPRT, SEQQUIT, SUPINFO

      IMPLICIT NONE

      LOGICAL                         :: LEXIST              ! T/F depending on whether file bandit.f07 exists

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'AUTO_SEQ_PROC'
      CHARACTER( 1*BYTE)              :: IS_GRID_SEQD(NGRID) ! 'Y'/'N' indicator of whether a grid was sequenced in file SEQ
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD                ! Card image read from SEQ file (which was generated by Bandit)
      CHARACTER(85*BYTE)              :: VEC_DESCR           ! Char name of array which has the grids that Bandit did not sequence

      INTEGER(LONG)                   :: BANDIT_MAX_SEQ_NUM  ! Max sequence number from SEQ2 that is generated by Bandit sequencing
      INTEGER(LONG)                   :: BANDIT_NSEQ         ! Number of grids sequenced by Bandit
      INTEGER(LONG)                   :: GRID_NOT_SEQD(NGRID)! Array of grid numbers not sequenced by Bandit
      INTEGER(LONG)                   :: GRID_ROW_NUM        ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: IERR0       = 0     ! Local error count
      INTEGER(LONG)                   :: IERR1       = 0     ! Local error count
      INTEGER(LONG)                   :: INT_SEQ2            ! Integer value for real SEQ2 sequence number
      INTEGER(LONG)                   :: IOCHK       = 0     ! IOSTAT error number when opening/reading a file
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD
      INTEGER(LONG)                   :: K                   ! Counter
      INTEGER(LONG)                   :: NUM_GRIDS_NOT_SEQD=0! Number of grids not sequenced by Bandit
      INTEGER(LONG)                   :: NUM_LEFT            ! Number used in writing 4 pairs of SEQ1,2 to SEQ file
      INTEGER(LONG)                   :: NUM_SEQ_FILE_LINES  ! Number of lines in the SEQ file (after STIME)
      INTEGER(LONG)                   :: OUNT(2)             ! File units to write messages to
      INTEGER(LONG)                   :: REC_NO      = 0     ! Indicator of record number when error encountered reading file
      INTEGER(LONG)                   :: XTIME               ! Time stamp read from a file


      INTRINSIC                       :: DBLE, INT



! **********************************************************************************************************************************
      IF (BANDIT_ERR /= 0) THEN
         WRITE(ERR,8888) BANDIT_ERR
         WARN_ERR = WARN_ERR + 1
         WRITE(F06,8888) BANDIT_ERR
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8888) BANDIT_ERR
         ENDIF
         CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, SEQSTAT )
         RETURN
      ENDIF

      NUM_SEQ_FILE_LINES = 0

      OUNT(1) = ERR
      OUNT(2) = F06

      INQUIRE ( FILE=SEQFIL, EXIST=LEXIST )
exist:IF (LEXIST) THEN                                     ! SEQFIL does exist, so open and read data

         IERR0 = 0                                         ! Open SEQ file, read and check STIME
         OPEN (SEQ,FILE=SEQFIL,STATUS='OLD',IOSTAT=IOCHK)
         IF (IOCHK /= 0) THEN
            CALL OPNERR ( IOCHK, SEQFIL, OUNT )
            IERR0 = IERR0 +1
            IF (IOCHK < 0) THEN                            ! File cannot be opened
               WRITE(ERR,9991) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9991) SEQFIL, GRIDSEQ
                  WRITE(F06,9998)
               ENDIF
            ELSE                                           ! ERR reading file
               WRITE(ERR,9999) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               WRITE(F06,9999) SEQFIL, GRIDSEQ
               WRITE(F06,9998)
            ENDIF
            CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, 'KEEP' )
            RETURN
         ELSE                                              ! No OPEN error, so read and check STIME
            READ(SEQ,'(1X,I11)',IOSTAT=IOCHK) XTIME
            IF (IOCHK /= 0) THEN
               REC_NO = 1
               CALL READERR ( IOCHK, SEQFIL, 'STIME', REC_NO, OUNT )
               IERR0 = IERR0 + 1
            ELSE                                           ! No error reading XTIME, so check XTIME = STIME
               IF (XTIME /= STIME) THEN
                  CALL STMERR ( XTIME, SEQFIL, OUNT )
                  IERR0 = IERR0 + 1
               ENDIF
            ENDIF
         ENDIF

 err0:   IF (IERR0 == 0) THEN                              ! No problem opening SEQ file and reading STIME, so process SEQGP entries
            CALL DEALLOCATE_MODEL_STUF ( 'SEQ1,2' )
            LSEQ = NGRID
            NSEQ = 0
            CALL ALLOCATE_MODEL_STUF ( 'SEQ1,2', SUBR_NAME )

            REC_NO = 0
            IERR1  = 0
i_do1:      DO                                             ! Loop reading SEQGP records until end of file
               READ (SEQ,'(A)',IOSTAT=IOCHK) CARD
               REC_NO = REC_NO + 1
               IF      (IOCHK <  0) THEN                   ! EOF/EOR so exit
                  EXIT
               ELSE IF (IOCHK >  0) THEN                   ! Error reading a SEQGP card
                  CALL READERR ( IOCHK, SEQFIL, 'SEQGP cards', REC_NO, OUNT )
                  IERR1 = IERR1 + 1
               ELSE                                        ! READ was OK so process record
                  IF (CARD(1:5) == 'SEQGP   ') THEN        ! If SEQGP image, call BD_SEQGP to read it and to add to SEQ1,2 arrays
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                     CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
                     DO J=1,10
                        CALL LEFT_ADJ_BDFLD ( JCARD(J) )
                     ENDDO
                     CALL MKCARD ( JCARD, CARD )
                     CALL BD_SEQGP ( CARD )
                  ELSE
                     CYCLE i_do1
                  ENDIF
               ENDIF
            ENDDO i_do1
            BANDIT_NSEQ = NSEQ

            IF (IERR1 > 0) THEN
               WRITE(ERR,*)
               WRITE(ERR,9992) SEQFIL, GRIDSEQ
               WRITE(ERR,9998)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,*)
                  WRITE(F06,9992) SEQFIL, GRIDSEQ
                  WRITE(F06,9998)
               ENDIF
               CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, 'KEEP' )
            ENDIF

            IF (NUM_SEQ_FILE_LINES > 0) THEN
               WRITE(ERR,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
               WRITE(ERR,*)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
                  WRITE(F06,*)
               ENDIF
            ELSE
               WRITE(ERR,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
               WRITE(ERR,9994)
               WRITE(ERR,9998)
               WRITE(ERR,*)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9993) SEQFIL, GRIDSEQ, NUM_SEQ_FILE_LINES
                  WRITE(F06,9994)
                  WRITE(F06,9998)
                  WRITE(F06,*)
               ENDIF
            ENDIF

! Finished reading the SEQGP entries in the Bandit SEQ file. Now find out if all grids were sequenced. If not, sequence them
! after the last Bandit sequence number, add them to the SEQ file, and tell user.

            NUM_GRIDS_NOT_SEQD = NGRID - NSEQ

            IF (NUM_GRIDS_NOT_SEQD == 0) THEN                ! Bandit sequenced all grids

               WRITE(ERR,10000) NGRID, GRIDSEQ
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,10000) NGRID, GRIDSEQ
               ENDIF

            ELSE                                             ! Some grids were not sequenced by Bandit

               WRITE(ERR,9996) SEQFIL, GRIDSEQ, NGRID, NSEQ
               WRITE(F06,9996) SEQFIL, GRIDSEQ, NGRID, NSEQ

               DO I=1,NGRID                                  ! Array IS_GRID_SEQD will show which grids were not seq'd by Bandit
                  IS_GRID_SEQD(I) = 'N'
               ENDDO
               DO I=1,NSEQ
                  CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SEQ1(I), GRID_ROW_NUM )
                  IF (GRID_ROW_NUM > 0) THEN
                     IS_GRID_SEQD(GRID_ROW_NUM) = 'Y'
                  ENDIF
               ENDDO

               K = 0                                       ! Array GRID_NOT_SEQD has grid no's for grids not seq'd by Bandit
               DO I=1,NGRID
                  IF (IS_GRID_SEQD(I) == 'N') THEN
                     K = K + 1
                     GRID_NOT_SEQD(K) = GRID_ID(I)
                  ENDIF
               ENDDO
                                                           ! Write list of grids not sequenced by Bandit
               CALL SORT_INT1 ( SUBR_NAME, 'GRID_NOT_SEQD', NUM_GRIDS_NOT_SEQD, GRID_NOT_SEQD )
               VEC_DESCR(1:) = 'GRIDS WHICH BANDIT DID NOT SEQUENCE (PERHAPS THEY HAD NO ELEMENTS CONNECTED TO THEM)'
               IF (NUM_GRIDS_NOT_SEQD > 10) THEN
                  CALL WRITE_INTEGER_VEC ( VEC_DESCR, GRID_NOT_SEQD, NUM_GRIDS_NOT_SEQD )
               ELSE
                  WRITE(F06,10003) VEC_DESCR
                  DO I=1,NUM_GRIDS_NOT_SEQD
                     WRITE(F06,10004) GRID_NOT_SEQD(I)
                  ENDDO
               ENDIF
               WRITE(F06,*)
                                                           ! Sequence the grids not sequenced by Bandit (in grid numerical order)
               BANDIT_MAX_SEQ_NUM = 0                      !  (a) Get max sequence number generated by Bandit
               DO I=1,NSEQ
                  INT_SEQ2 = INT(SEQ2(I))
                  IF (INT_SEQ2 > BANDIT_MAX_SEQ_NUM) THEN
                     BANDIT_MAX_SEQ_NUM = INT_SEQ2
                  ENDIF
               ENDDO
               IF ((NSEQ + NUM_GRIDS_NOT_SEQD) > LSEQ) THEN!  (b) Check to make sure we won't exceed size of allocated arrays SEQ1,2
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1301) SUBR_NAME,NSEQ,NUM_GRIDS_NOT_SEQD,LSEQ
                  WRITE(F06,1301) SUBR_NAME,NSEQ,NUM_GRIDS_NOT_SEQD,LSEQ
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               DO I=1,NUM_GRIDS_NOT_SEQD                   !  (c) Sequence the grids starting with previous NSEQ (from Bandit)
                  NSEQ = NSEQ + 1
                  SEQ1(NSEQ) = GRID_NOT_SEQD(I)            !      NOTE: GRID_NOT_SEQD was sorted above so it is in grid num order
                  SEQ2(NSEQ) = DBLE(BANDIT_MAX_SEQ_NUM + I)
               ENDDO

               NUM_LEFT = NUM_GRIDS_NOT_SEQD               ! Write data for grids not seq'd by Bandit to the SEQ file
               BACKSPACE (SEQ)
               WRITE(SEQ,300) PROG_NAME, NUM_LEFT
               NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 2
               DO I=1,NUM_GRIDS_NOT_SEQD,4
                  IF (NUM_LEFT >= 4) THEN
                     WRITE(SEQ,301) (SEQ1(K),INT(SEQ2(K)),K=I+BANDIT_NSEQ,I+BANDIT_NSEQ+3)
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                     NUM_LEFT = NUM_LEFT - 4
                  ELSE
                     WRITE(SEQ,301) (SEQ1(K),INT(SEQ2(K)),K=I+BANDIT_NSEQ,I+BANDIT_NSEQ+NUM_LEFT-1)
                     NUM_SEQ_FILE_LINES = NUM_SEQ_FILE_LINES + 1
                  ENDIF
               ENDDO

               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,9998)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,9998)
               ENDIF

            ENDIF

! Write SEQ cards to F06, if user requested it

            IF (SEQPRT == 'Y') THEN                        ! See if user wants to write SEQGP card images to F06

               IERR1 = 0
               REWIND (SEQ)
               READ(SEQ,'(1X,I11)',IOSTAT=IOCHK) XTIME
               IF (IOCHK /= 0) THEN
                  REC_NO = 1
                  CALL READERR ( IOCHK, SEQFIL, 'STIME', REC_NO, OUNT )
                  IERR1 = IERR1 + 1
               ELSE                                        ! No error reading XTIME, so check XTIME = STIME
                  IF (XTIME /= STIME) THEN
                     CALL STMERR ( XTIME, SEQFIL, OUNT )
                     IERR1 = IERR1 + 1
                  ENDIF
               ENDIF

               IF (IERR1 == 0) THEN
                  REC_NO = 0
                  WRITE(F06,9997) SEQFIL
i_do2:            DO I=1,NUM_SEQ_FILE_LINES
                     READ (SEQ,'(A)',IOSTAT=IOCHK) CARD
                     REC_NO = REC_NO + 1
                     IF (IOCHK /= 0) THEN
                        WRITE(F06,10001) I
                     ELSE
                        WRITE(F06,10002) I, CARD
                     ENDIF
                  ENDDO i_do2
                  WRITE(F06,*)
               ENDIF

            ENDIF

            CALL FILE_CLOSE ( SEQ, SEQFIL, SEQSTAT )

         ELSE

            CALL AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, SEQSTAT )

         ENDIF err0

      ELSE                                                 ! bandit.f07 does not exist. Write message and quit if SEQQUIT = Y

         WRITE(ERR,9990) SEQFIL, GRIDSEQ
         WRITE(ERR,9994)
         WRITE(ERR,9998)

         IF (SUPINFO == 'N') THEN
            WRITE(F06,9990) SEQFIL, GRIDSEQ
            WRITE(F06,9994)
            WRITE(F06,9998)
         ENDIF

         RETURN

      ENDIF exist



      RETURN

! **********************************************************************************************************************************
  300 FORMAT('$ The following SEQGP entries were generated by ',A,' to sequence',/,                                                &
             '$ the ',I8,' grids not sequenced by Bandit.')

  301 FORMAT('SEQGP   ',8I8)

 1301 FORMAT(' *ERROR  1301: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' (NSEQ + NUM_GRIDS_NOT_SEQD) = ',I8,' MUST BE < LSEQ = ',I8,' SO ALLOCATED MEM FOR ARRAYS SEQ1,2',     &
                           ' WON''T BE EXCEEDED')


 8888 FORMAT(' *WARNING    : BANDIT DID NOT RUN SUCCESSFULLY. THIS IS A NON-CRITICAL OPTIMIZATION. IT QUIT WITH ERROR = ',I8       &
                    ,/,14X,' CHECK FILE BANDIT.OUT IN THE DIRECTORY WHERE MYSTRAN.EXE RESIDES')

 9990 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), DOES NOT EXIST',/)

 9991 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), CANNOT BE OPENED.',/)

 9999 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), CANNOT BE READ.',/)

 9992 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), HAS HAD ABOVE LISTED READ ERROR(S).',/)

 9993 FORMAT(' *INFORMATION: FILE ',A                                                                                              &
                    ,/,14X,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,'), HAS BEEN READ. THERE WERE ',I8,' SEQGP ENTRY IMAGES READ.')

 9994 FORMAT(14X,' IT MAY BE THAT BANDIT FOUND THAT NO RESEQUENCING WAS NEEDED OR DUE TO ERROR IN RUNNING BANDIT.',/)

 9996 FORMAT(' *WARNING    : FILE ',A                                                                                              &
                    ,/,14x,' CONTAINING THE BULK DATA SEQGP ENTRY IMAGES (NEEDED FOR AUTO GRID POINT SEQUENCING REQUESTED BY'      &
                    ,/,14X,' THE USER VIA PARAM GRIDSEQ ',A,') DOES NOT HAVE ALL GRIDS DEFINED ON SEQGP ENTRIES.'                  &
                    ,/,14X,' THE NUMBER OF GRIDS = ',I8,' BUT THE NUMBER OF GRIDS SEQUENCED = ',I8,/)

 9997 FORMAT(' *INFORMATION: THE FOLLOWING IS A LISTING OF THE SEQGP RECORDS READ FROM FILE ',A,':',/)

 9998 FORMAT(14X,' MAKE SURE BANDIT HAS RUN SUCCESSFULLY (CHECK FILE bandit.out IN THE DIRECTORY WHERE MYSTRAN.EXE RESIDES).'/)

10000 FORMAT(' *INFORMATION: ALL BANDIT SEQGP ENTRY IMAGES HAVE BEEN USED TO SUCCESSFULLY RESEQUENCE ',I8,' GRIDS BASED ON PARAM', &
                           ' GRIDSEQ ',A,/)

10001 FORMAT(15X,'Card image ',I7,': could not be written')

10002 FORMAT(15X,'Card image ',I7,': ',A)

10003 FORMAT(24X,A,/)

10004 FORMAT(62X,I8)

! **********************************************************************************************************************************

      END SUBROUTINE AUTO_SEQ_PROC

! ##################################################################################################################################

      SUBROUTINE AUTO_SEQ_PROC_WRAPUP ( SUBR_NAME, CLOSE_STAT )

! Writes message for subr AUTO_SEQ_PROC

      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  SEQQUIT

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: SUBR_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: CLOSE_STAT        ! Status for closing SEQFIL

! **********************************************************************************************************************************
      CALL FILE_CLOSE ( SEQ, SEQFIL, CLOSE_STAT )

      IF (SEQQUIT == 'Y') THEN
         WRITE(ERR,8881) SUBR_NAME, SEQQUIT
         WRITE(F06,8881) SUBR_NAME, SEQQUIT
         CALL OUTA_HERE ( 'Y' )
      ELSE
         WRITE(ERR,8886) SEQQUIT
         WRITE(F06,8886) SEQQUIT
      ENDIF

      RETURN

! **********************************************************************************************************************************
 8881 FORMAT(14X,' PROCESSING TERMINATED IN SUBR ',A,' BASED ON ABOVE PROBLEM AND FIELD 4 OF PARAM GRIDSEQ = ',A,'.'               &
          ,/,14X,' (THIS WAS EITHER ENTERED ON A BULK DATA PARAM GRIDSEQ ENTRY OR IS THE DEFAULT. SEE MYSTRAN DOCUMENTATION)',/)

 8886 FORMAT(14X,' SINCE FIELD 4 OF PARAM GRIDSEQ = ',A,', MYSTRAN WILL DEFAULT TO GRID SEQUENCING BASED ON GRID NUMERICAL ORDER.' &
          ,/,14X,' (THIS WAS EITHER ENTERED ON A BULK DATA PARAM GRIDSEQ ENTRY OR IS THE DEFAULT. SEE MYSTRAN DOCUMENTATION)',/)

! **********************************************************************************************************************************

      END SUBROUTINE AUTO_SEQ_PROC_WRAPUP

      END SUBROUTINE SEQ_PROC

   END MODULE GRID_COORDINATE_PROCESSING
