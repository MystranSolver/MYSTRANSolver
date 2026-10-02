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

   MODULE GRAV_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GRAV_PROC

   CONTAINS

      SUBROUTINE GRAV_PROC

! Gravity load processor

! Transforms the input B.D. gravity load data to system force data in the SYS_LOAD array for dof's in the G set.
! File LINK1P was written when GRAV B.D cards were read, with one record for each such card.
! There are NGRAV total number of records written to file LINK1P with each record containing:

!               SETID         = Load set ID
!               ACID_L        = Local coord sys ID that gravity load is given in
!               GRAV_GRID     = ID of grid that rotational (components 4, 5, 6) grav accels are about
!               ACCEL_RB(1-6) = 6 components of gravity load

! The process in creating array SYS_LOAD from this information is as follows:

!  (1) For each record (1 to NGRAV) in file LINK1P:

!      (a) Read a record

!      (b) Transform coords from local (on GRAV card) to basic. Transformation to global is done later (see 2a(iii))

!      (d) Write records, similar to those in LINK1P, but in global coords, to a scratch file (SCRATCH-99-).

!  (2) For each subcase (1 to NSUB):

!      (a) Generate LSID, RSID tables of load set ID's/scale factors for load sets for this subcase:

!          (  i) LLOADC is the max number of pairs of scale factors/load set ID's over all LOAD Bulk Data cards
!                including the pair defined by the set ID and overall scale factor on the LOAD Bulk Data card.
!                LSID and RSID are dimensioned 1 to LLOADC and:
!                   LSID(1) is always the load set ID requested in Case Control for this subcase.
!                   RSID(1) is always 1.0

!          ( ii) If the load set requested in Case Control is not a set ID from a LOAD Bulk data card then the
!                load must be on a separate GRAV card in which case LSID(1) and RSID(1) are all that is needed
!                to define the load set contribution due to the GRAV card for this subcase.

!          (iii) If the load set requested in Case Control is a set ID from a LOAD Bulk data card then the ramainder
!                (K = 2,LLOADC) of entries into LSID and RSID will be the pairs of load set ID's/scale factors from
!                that LOAD Bulk Data card (with RSID also multiplied by the overall scale factor on the LOAD Bulk data
!                card. The load set ID's are in array LOAD_SIDS(i,j) created when LOAD Bulk Data cards were read.
!                The scale factors are in array LOAD_FACS(i,j) also created when LOAD Bulk Data cards were read.
!                Note, there may not be as many as LLOADC pairs of set ID's/scale factors on a given LOAD Bulk Data
!                card since LLOADC is the max, from all LOAD Bulk Data cards, of pairs.
!                Thus, the entries in LSID from the last entry (for a given LOAD card) to LLOADC will be zero (LSID
!                was initialized to zero). This fact is used in a DO loop to EXIT when LSID(K) = 0

!      (b) For each record in SCRATCH-991 (1 to NGRAV)

!          (  i) Read a record

!          ( ii) Scan LSID and RSID to get the scale factor (SCALE) for the ACCEL_RB components in SETID, if this
!                GRAV's set ID is in LSID. When found, reset grav vector components to ACCEL_RB = SCALE*ACCEL_RB.
!                At this point ACCEL_RB has the correct magnitudes and is in basic coordinates.

!          (iii) For a grid point, determine if a transformation of ACCEL_RB to global is needed and transform it if so.

!          ( iv) Calculate accel at a grid (ACCEL_I) based on rigid body motion of ACCEL_RB

!          (  v) Get the 6 x 6 mass matrix for a grid point times ACCEL_I to get grav forces at this grid

!          ( vi) Load the grav forces into the SYS_LOAD (systems load) array

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, ERR, F06, SCR, L1P, LINK1P, L1P_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LLOADC, NCORD, NGRAV, NGRID, NLOAD, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  CORD, GRID, GRID_ID, LOAD_FACS, LOAD_SIDS, RCORD, RGRID, SYS_LOAD, SUBLOD

      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, READERR
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE LINK1_LOAD_SUPPORT, ONLY    :  GET_GRID_6X6_MASS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GRAV_PROC'
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for
      CHARACTER( 1*BYTE)              :: GRID_MGG_FND      ! Indicator on whether a mass matrix was found in MGG for a given grid
      CHARACTER( 8*BYTE)              :: NAME              ! Name for output error purposes
      CHARACTER(24*BYTE)              :: MESSAG            ! File description.
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG)                   :: ACID      = 0     ! A coord system ID (actual)
      INTEGER(LONG), PARAMETER        :: ACID_0    = 0     ! Basic coord system
      INTEGER(LONG)                   :: ACID_L            ! Actual local  coord sys ID on FORCE or MOMENT card
      INTEGER(LONG)                   :: ACID_G            ! Actual global coord sys ID for an actual grid
      INTEGER(LONG)                   :: GRAV_GRID         ! ID of grid that rotational grav accels are about
      INTEGER(LONG)                   :: CID_ERR   = 0     ! Count of coord systems undefined
      INTEGER(LONG)                   :: GID_ERR   = 0     ! Count of grids undefined
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number for a grid
      INTEGER(LONG)                   :: GRAV_GRID_ROW_NUM ! Row number in array GRID_ID where an actual grid ID is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: ICID              ! Internal coordinate system ID for ACID_L or ACID_G
      INTEGER(LONG)                   :: IERRT             ! Total number of errors found
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: SETID             ! Load set ID read from record in file LINK1P
      INTEGER(LONG)                   :: LSID(LLOADC+1)     ! Array of load SID's, for GRAV cards, needed for one S/C
      INTEGER(LONG)                   :: NCOLA             ! No. cols in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NCOLB             ! No. cols in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NROWA             ! No. rows in a matrix. For subr MATMULT_FFF/MATMULT_FFF_T, called herein
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: READ_ERR  = 0     ! Cum. count of errors as we read, and check cards from file LINK1K
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM           ! Row no. in array TDOF corresponding to GDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! Row no. in array TDOF where data begins for a grid


      REAL(DOUBLE)                    :: ACCEL_I(6)        ! 6 components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_I_T1(3)     ! 3 translational components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_I_T2(3)     ! 3 translational components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_I_R1(3)     ! 3 rotational grav components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_I_R2(3)     ! 3 rotational grav components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_RB(6)       ! 6 components of rigid body accel due to GRAV
      REAL(DOUBLE)                    :: ACCEL_RB_T1(3)    ! 3 translational components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_RB_T2(3)    ! 3 translational components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_RB_R1(3)    ! 3 rotational grav components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: ACCEL_RB_R2(3)    ! 3 rotational grav components of accel due to gravity at a grid
      REAL(DOUBLE)                    :: DXI               ! X (basic) distance between a grid and GRAV_GRID
      REAL(DOUBLE)                    :: DYI               ! Y (basic) distance between a grid and GRAV_GRID
      REAL(DOUBLE)                    :: DZI               ! Z (basic) distance between a grid and GRAV_GRID
      REAL(DOUBLE)                    :: FORCE_I(6)        ! 6 forces at a grid due to the grav loading
      REAL(DOUBLE)                    :: GRID_MGG(6,6)     ! 6 X 6 mass matrix for one grid point
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: RSID(LLOADC+1)     ! Array of grav magnitudes (for LSID set ID's) needed for one S/C
      REAL(DOUBLE)                    :: SCALE             ! Scale factor for a load
      REAL(DOUBLE)                    :: T12(3,3)          ! Coord transformation matrix
      REAL(DOUBLE)                    :: XGID              ! Basic X coord of GRAV_GRID
      REAL(DOUBLE)                    :: YGID              ! Basic Y coord of GRAV_GRID
      REAL(DOUBLE)                    :: ZGID              ! Basic Z coord of GRAV_GRID

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      NAME = 'GRAV    '

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Open a scratch file that will be used to rewrite data from LINK1P after the coords have been transformed to global.

      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',POSITION='REWIND',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
         CALL OUTA_HERE ( 'Y' )                                    ! Error opening scratch file, so quit
      ENDIF
      REWIND (SCR(1))

! **********************************************************************************************************************************
! (1) Read record from L1P and transform coords for GRAV from local sys (on GRAV bulk data card) to basic and write results to
!     SCRATCH-991.

i_do1:DO I=1,NGRAV
                                                           ! Read a record from L1P
         READ(L1P,IOSTAT=IOCHK) SETID, ACID_L, GRAV_GRID, (ACCEL_RB(J),J=1,6)

         IF (IOCHK /= 0) THEN
            REC_NO = I
            CALL READERR ( IOCHK, LINK1P, L1P_MSG, REC_NO, OUNT )
            READ_ERR = READ_ERR + 1                        ! Increment READ_ERR and go back to read another grav card
            CYCLE i_do1
         ENDIF

         IF (GRAV_GRID > 0) THEN
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRAV_GRID, GRAV_GRID_ROW_NUM )
            IF (GRAV_GRID_ROW_NUM == -1) THEN
               WRITE(ERR,1822) 'GRID ', GRAV_GRID, NAME, SETID
               WRITE(F06,1822) 'GRID ', GRAV_GRID, NAME, SETID
               GID_ERR = GID_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDIF

         DO J=1,3
            ACCEL_RB_T1(J) = ACCEL_RB(J)
            ACCEL_RB_R1(J) = ACCEL_RB(J+3)
         ENDDO
                                                           ! The local system that gravity is defined in is ACID_L.
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
               WRITE(ERR,1822) 'COORD SYSTEM ', ACID_L, NAME, SETID
               WRITE(F06,1822) 'COORD SYSTEM ', ACID_L, NAME, SETID
               CID_ERR = CID_ERR + 1                       ! Increment READ_ERR and go back to read another grav card
               FATAL_ERR = FATAL_ERR + 1
               CYCLE i_do1
            ENDIF

            DO J=1,3                                       ! Get coord transf matrix (don't need GEN_T0L since ICID is rect.)
               DO K=1,3
                  L = 3 + 3*(J-1) + K
                  T12(J,K) = RCORD(ICID,L)
               ENDDO
            ENDDO

            NROWA  = 3                                     ! Transform coordinates
            NCOLA  = 3
            NCOLB  = 1
            CALL MATMULT_FFF ( T12, ACCEL_RB_T1, NROWA, NCOLA, NCOLB, ACCEL_RB_T2 )
            CALL MATMULT_FFF ( T12, ACCEL_RB_R1, NROWA, NCOLA, NCOLB, ACCEL_RB_R2 )
         ELSE                                              ! No transformation needed since ACID_L is basic
            DO J=1,3
               ACCEL_RB_T2(J) = ACCEL_RB_T1(J)
               ACCEL_RB_R2(J) = ACCEL_RB_R1(J)
            ENDDO
         ENDIF

         DO J=1,3                                          ! Reset ACCEL_RB(I) to values in basic coords.
            ACCEL_RB(J)   = ACCEL_RB_T2(J)
            ACCEL_RB(J+3) = ACCEL_RB_R2(J)
         ENDDO
                                                           ! Write data to scratch file. GRAV now in basic coords
         WRITE(SCR(1)) SETID,ACID_0,GRAV_GRID,(ACCEL_RB(J),J=1,6)

      ENDDO i_do1

! Quit if CID_ERR, GID_ERR or READ_ERR > 0

      IF ((CID_ERR > 0) .OR. (GID_ERR > 0) .OR. (READ_ERR > 0)) THEN
         IERRT = CID_ERR + GID_ERR + READ_ERR
         WRITE(ERR,1599) SUBR_NAME,IERRT
         WRITE(F06,1599) SUBR_NAME,IERRT
         CALL OUTA_HERE ( 'Y' )                            ! Errors from reading grav data, so quit
      ENDIF

      REWIND (SCR(1))

! **********************************************************************************************************************************
! Now process grav loads into SYS_LOAD

      DO I=1,LLOADC                                         ! Initialize LSID, RSID arrays
         LSID(I) = 0
         RSID(I) = ZERO
      ENDDO

      MESSAG = 'SCRATCH FILE IN FORCEP  '

      WRITE(SC1, * )
      IERRT = 0
i_do2:DO I = 1,NSUB                                        ! Loop through the S/C's

         IF (SUBLOD(I,1) == 0) THEN                        ! If no load for this S/C, CYCLE
            CYCLE i_do2
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(I,1)                             ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the GRAV card in file LINK1P that matches SUBLOD(I,1)
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

j_do_22: DO J = 1,NGRAV                                    ! Process GRAV card info that is now in basic coords
                                                           ! (2-b-  i) Read a GRAV record (in basic coords) from scratch unit
            READ(SCR(1),IOSTAT=IOCHK) SETID,ACID,GRAV_GRID,(ACCEL_RB(K),K=1,6)
            IF (IOCHK /= 0) THEN
               REC_NO = J
               CALL READERR ( IOCHK, SCRFIL, MESSAG, REC_NO, OUNT )
               CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
               CALL OUTA_HERE ( 'Y' )                      ! Error reading scratch file, so quit
            ENDIF

            IF (ACID /= 0) THEN
               WRITE(ERR,1515) SUBR_NAME,ACID
               WRITE(F06,1515) SUBR_NAME,ACID
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error, so quit
            ENDIF

            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
k_do221:    DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                  DO L=1,3
                     ACCEL_RB(L) = SCALE*ACCEL_RB(L)
                  ENDDO
                  EXIT k_do221
               ENDIF
            ENDDO k_do221

            IF (FOUND /= 'Y') THEN                         ! This GRAV set ID isn't called for in this S/C, so CYCLE on GRAV load
               CYCLE j_do_22
            ENDIF

            DO K=1,3
               ACCEL_RB_T1(K) = ACCEL_RB(K)
               ACCEL_RB_R1(K) = ACCEL_RB(K+3)
            ENDDO

k_d0222:    DO K = 1,NGRID
               WRITE(SC1,12345,ADVANCE='NO') K, NGRID, I, CR13

               ACID_G = GRID(K,3)                          ! The global coord sys for this grid is ACID_G

               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(K), IGRID )
               IF (IGRID == -1) THEN
                  IERRT     = IERRT + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', GRID_ID(K), NAME, SETID
                  WRITE(F06,1822) 'GRID ', GRID_ID(K), NAME, SETID
               ENDIF
               IF (GRID(IGRID,6) == 1) THEN                ! Scalar point so do not generate grav load on it. Give warn, not fatal
                  IERRT = IERRT + 1
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,9901) GRID_ID(K)
                  IF (SUPWARN == 'N') THEN
                     WRITE(F06,9901) GRID_ID(K)
                  ENDIF
                  CYCLE k_d0222
               ENDIF

               IF (GRAV_GRID > 0) THEN
                  XGID = RGRID(GRAV_GRID_ROW_NUM,1)
                  YGID = RGRID(GRAV_GRID_ROW_NUM,2)
                  ZGID = RGRID(GRAV_GRID_ROW_NUM,3)
               ELSE
                  XGID = ZERO
                  YGID = ZERO
                  ZGID = ZERO
               ENDIF

               DXI = RGRID(K,1) - XGID
               DYI = RGRID(K,2) - YGID
               DZI = RGRID(K,3) - ZGID


               ACCEL_I(1) = ACCEL_RB(1) + DZI*ACCEL_RB(5) - DYI*ACCEL_RB(6)
               ACCEL_I(2) = ACCEL_RB(2) - DZI*ACCEL_RB(4) + DXI*ACCEL_RB(6)
               ACCEL_I(3) = ACCEL_RB(3) + DYI*ACCEL_RB(4) - DXI*ACCEL_RB(5)
               ACCEL_I(4) = ACCEL_RB(4)
               ACCEL_I(5) = ACCEL_RB(5)
               ACCEL_I(6) = ACCEL_RB(6)
               DO L=1,3
                  ACCEL_I_T1(L) = ACCEL_I(L)
                  ACCEL_I_R1(L) = ACCEL_I(L+3)
               ENDDO

               IF (ACID_G /= 0) THEN                       ! ACID_G is not basic so transform coords of ACCEL_RB vector to global
l_do_2211:        DO L=1,NCORD
                     IF (CORD(L,2) == ACID_G) THEN
                        ICID = L                           ! ICID is the internal coord sys ID corresponding to ACID_G
                        EXIT l_do_2211
                     ENDIF
                  ENDDO l_do_2211
                                                           ! T12 is coord transf matrix that would transf a global vector to basic
!                                                            and its transpose will transform a basic coord sys vector to global
                  CALL GEN_T0L ( K, ICID, THETAD, PHID, T12 )

                  NROWA  = 3
                  NCOLA  = 3
                  NCOLB  = 1
                  CALL MATMULT_FFF_T ( T12, ACCEL_I_T1, NROWA, NCOLA, NCOLB, ACCEL_I_T2 )
                  CALL MATMULT_FFF_T ( T12, ACCEL_I_R1, NROWA, NCOLA, NCOLB, ACCEL_I_R2 )
                  DO L=1,3
                     ACCEL_I(L)   = ACCEL_I_T2(L)          ! ACCEL_I (grid accel) is now in global coords if it was not before
                     ACCEL_I(L+3) = ACCEL_I_R2(L)
                  ENDDO
               ENDIF
                                                           ! Get internal ID for the grid for which we want 6x6 mass
               CALL GET_GRID_6X6_MASS ( GRID_ID(K), IGRID, GRID_MGG_FND, GRID_MGG )

               DO L=1,6
                  FORCE_I(L) = ZERO
               ENDDO

               IF (GRID_MGG_FND == 'Y') THEN
                  NROWA  = 6
                  NCOLA  = 6
                  NCOLB  = 1
                  CALL MATMULT_FFF ( GRID_MGG, ACCEL_I, NROWA, NCOLA, NCOLB, FORCE_I )
               ENDIF

!xx            CALL CALC_TDOF_ROW_NUM ( GRID(K,1), ROW_NUM_START, 'N' )
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID(K,1), IGRID )
               ROW_NUM_START = TDOF_ROW_START(IGRID)
l_do_2214:     DO L = 1,6
                  CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
                  ROW_NUM = ROW_NUM_START + L - 1
                  GDOF = TDOF(ROW_NUM,G_SET_COL_NUM)
                  SYS_LOAD(GDOF,I) = SYS_LOAD(GDOF,I) + FORCE_I(L)
               ENDDO l_do_2214

            ENDDO k_d0222

         ENDDO j_do_22
         REWIND (SCR(1))                                   ! Need to read all of the GRAV records again for the next S/C

      ENDDO i_do2

      WRITE(SC1,*) CR13

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )



      RETURN

! **********************************************************************************************************************************
 1515 FORMAT(' *ERROR  1515: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' COORD SYSTEM ID READ FROM FILE SCRATCH-991 MUST BE BASIC (0) BUT IS ',I8)

 1599 FORMAT(/,' PROCESSING TERMINATED IN SUBROUTINE ',A,' DUE TO ABOVE LISTED ',I8,' ERRORS')

 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 9901 FORMAT(' *WARNING    : NO GRAV LOAD IS BEING PROCESSED FOR GRID ',I8,' SINCE IT IS A SCALAR POINT (SPOINT)')

12345 FORMAT(5X,'Proc grid ',I8,' of ',I8,', subcase ',I8, A)



! **********************************************************************************************************************************

      END SUBROUTINE GRAV_PROC

   END MODULE GRAV_PROCESSING
