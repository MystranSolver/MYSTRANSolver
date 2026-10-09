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

   MODULE FORCE_MOM_PROCESSING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: FORCE_MOM_PROC

   CONTAINS

      SUBROUTINE FORCE_MOM_PROC

! Force/moment processor

! Transforms the input B.D. force and moment data to system force data in the SYS_LOAD array for dof's in the G set.
! File LINK1I was written when FORCE or MOMENT B.D entries were read, with one record for each such card.
! There are NFORCE total no. records written to file LINK1I with each record containing:

!               SETID       = Load set ID on the FORCE or MOMENT card
!               AGRID       = GRID where load is located
!               ACID_L      = Local coord sys ID that load is given in
!               FORMON(1-3) = 3 components of force or moment in coord system ACID_L
!               NAME        = 'FORCE' or 'MOMENT' to describe whether B.D card was a FORCE or MOMENT card

! The process in creating array SYS_LOAD from this information is as follows:

!  (1) For each record (1 to NFORCE) in file LINK1I:

!      (a) Read a record: (SETID, AGRID, ACID_L, FORMON(1-3) and NAME)

!      (b) Transform coords from local to basic

!      (c) Transform from basic to global (as defined on the GRID card for AGRID).

!      (d) Write NFORCE records, similar to those in LINK1I, but now in global coords, to a scratch file (SCRATCH99).

!  (2) For each subcase (1 to NSUB):

!      (a) Generate LSID, RSID tables of load set ID's/scale factors for load sets for this subcase:

!          (  i) LLOADC is the max number of pairs of scale factors/load set ID's over all LOAD Bulk Data cards
!                including the pair defined by the set ID and overall scale factor on the LOAD Bulk Data card.
!                LSID and RSID are dimensioned 1 to LLOADC and:
!                   LSID(1) is always the load set ID requested in Case Control for this subcase.
!                   RSID(1) is always 1.0

!          ( ii) If the load set requested in Case Control is not a set ID from a LOAD Bulk data card then the
!                load must be on a separate FORCE/MOMENT in which case LSID(1) and RSID(1) are all that is needed
!                to define the load set contribution due to the FORCE/MOMENT cards for this subcase.

!          (iii) If the load set requested in Case Control is a set ID from a LOAD Bulk data card then the ramainder
!                (K = 2,LLOADC) of entries into LSID and RSID will be the pairs of load set ID's/scale factors from
!                that LOAD Bulk Data card (with RSID also multiplied by the overall scale factor on the LOAD Bulk data
!                card. The load set ID's are in array LOAD_SIDS(i,j) created when LOAD Bulk Data cards were read.
!                The scale factors are in array LOAD_FACS(i,j) also created when LOAD Bulk Data cards were read.
!                Note, there may not be as many as LLOADC pairs of set ID's/scale factors on a given LOAD Bulk Data
!                card since LLOADC is the max, from all LOAD Bulk Data cards, of pairs.
!                Thus, the entries in LSID from the last entry (for a given LOAD card) to LLOADC will be zero (LSID
!                was initialized to zero). This fact is used in a DO loop to EXIT when LSID(K) = 0

!      (b) For each record in SCRATCH-991 (1 to NFORCE)

!          (  i) Read a record from file: (SETID, AGRID, ACID_G, FORMON(i) (in coord sys ACID_ G), NAME)

!          ( ii) Scan LSID and RSID to get the scale factor for the FORMON components in SETID, if this
!                FORCE/MOMENT's set ID is in LSID

!          (iii) Load these force/moment values into the system load array, SYS_LOAD


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06, SCR, L1I, LINK1I, L1I_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LLOADC, NCORD, NFORCE, NGRID, NLOAD, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  LOAD_SIDS, LOAD_FACS, SYS_LOAD, SUBLOD, GRID, GRID_ID, CORD, RGRID

      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, READERR
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FORCE_MOM_PROC'
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for
      CHARACTER( 1*BYTE)              :: CORD_FND          ! = 'Y' if coord sys ID on FORCE/MOMENT defined, 'N' otherwise
      CHARACTER( 1*BYTE)              :: GRID_FND          ! = 'Y' if grid ID on FORCE/MOMENT defined, 'N' otherwise
      CHARACTER(24*BYTE)              :: MESSAG            ! File description.
      CHARACTER( 8*BYTE)              :: NAME              ! Name to indicate whether we are processing a FORCE or a MOMENT
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG)                   :: ACID_L            ! Actual local  coord sys ID on FORCE or MOMENT card
                                                           ! (-1 FORCE1/MOMENT1, -2 FORCE2/MOMENT2)
      INTEGER(LONG)                   :: GD(4)             ! G1-G4 of FORCE1/2, MOMENT1/2
      INTEGER(LONG)                   :: GROW(4)           ! Their rows in GRID_ID
      REAL(DOUBLE)                    :: DIR(3), A(3), B(3)! Direction of FORCE1/2, MOMENT1/2 in basic, and its two vectors
      REAL(DOUBLE)                    :: DLEN, SIZE_AB     ! Length of DIR, and of the vectors that make it
      INTEGER(LONG)                   :: ACID_G            ! Actual global coord sys ID for AGRID
      INTEGER(LONG)                   :: AGRID             ! Actual grid number from FORCE or MOMENT card
      INTEGER(LONG)                   :: COMP1, COMP2      ! DOF components (1-6)
      INTEGER(LONG)                   :: GDOF              ! G-set DOF no. for actual grid AGRID
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: ICID              ! Internal coordinate system ID for ACID_L or ACID_G
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error count
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: K1                ! Counter
      INTEGER(LONG)                   :: SETID             ! Load set ID read from record in file LINK1I
      INTEGER(LONG)                   :: LSID(LLOADC+1)     ! Array of load SID's, from a LOAD Bulk Data card, for one S/C
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: READ_ERR  = 0     ! Count of errors reading records from FORCE/MOMENT file
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM           ! Row no. in array TDOF corresponding to GDOF
      INTEGER(LONG)                   :: ROW_NUM_START     ! Row no. in array TDOF where data begins for AGRID


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: F1(3), F2(3)      ! 3 force or moment components in intermediate calcs
      REAL(DOUBLE)                    :: FORCEI            ! The force value that goes into SYS_LOAD for a grid/subcase
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: FORMON(3)         ! Array of 3 force/moment mag's read from file LINK1I
      REAL(DOUBLE)                    :: RSID(LLOADC+1)    ! Array of load magnitudes (for LSID set ID's) needed for one S/C
      REAL(DOUBLE)                    :: SCALE             ! Scale factor for a load
      REAL(DOUBLE)                    :: T12(3,3)          ! Coord transformation matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      NAME = '        '                                    ! will be read as FORCE or MOMENT from file L1I

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Open a scratch file that will be used to rewrite data from L1I after the coords have been transformed to
! global. This file is only needed in this subr and is closed and deleted herein.

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
! Successively read records from LINK1I and transform coords:

!   (1) From local (as defined on the FORCE/MOMENT card) to basic
!   (2) From basic to global (global for the grid point that the load is at)

! Then write out a record similar to what was read from LINK1I to a scratch file that has the same info as
! records from LINK1I but with data in global coords


i_do1:DO I=1,NFORCE
                                                           ! (1-a) Read a record from file LINK1I
         READ(L1I,IOSTAT=IOCHK) SETID,AGRID,ACID_L,(GD(J),J=1,4),(FORMON(J),J=1,3),NAME
         IF (IOCHK /= 0) THEN
            REC_NO = I
            CALL READERR ( IOCHK, LINK1I, L1I_MSG, REC_NO, OUNT )
            READ_ERR = READ_ERR + 1
            CYCLE i_do1
         ENDIF

         DO J=1,3
            F1(J) = FORMON(J)
         ENDDO

! From actual grid pt number (AGRID) on the FORCE/MOMENT card, get row number in array GRID_ID where AGRID exists

         GRID_FND = 'Y'
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
         IF (GRID_ID_ROW_NUM == -1) THEN
            GRID_FND = 'N'                                 ! Grid ID on FORCE/MOMENT undefined
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1822) 'GRID ', AGRID, NAME, SETID
            WRITE(F06,1822) 'GRID ', AGRID, NAME, SETID
         ENDIF

         CORD_FND = 'N'
         ICID      = -1
         IF (ACID_L < 0) THEN                              ! FORCE1/2, MOMENT1/2: magnitude F1(1) along a direction in basic
            DO J=1,-2*ACID_L
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GD(J), GROW(J) )
               IF (GROW(J) == -1) THEN
                  GRID_FND = 'N'
                  IERROR = IERROR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'GRID ', GD(J), NAME, SETID
                  WRITE(F06,1822) 'GRID ', GD(J), NAME, SETID
               ENDIF
            ENDDO
            IF (GRID_FND == 'N') CYCLE i_do1
            A = RGRID(GROW(2),1:3) - RGRID(GROW(1),1:3)    ! G1 to G2
            SIZE_AB = SQRT(DOT_PRODUCT(A,A))
            IF (ACID_L == -1) THEN
               DIR = A
            ELSE                                           ! (G2 - G1) x (G4 - G3)
               B = RGRID(GROW(4),1:3) - RGRID(GROW(3),1:3)
               SIZE_AB = SIZE_AB*SQRT(DOT_PRODUCT(B,B))
               DIR(1) = A(2)*B(3) - A(3)*B(2)
               DIR(2) = A(3)*B(1) - A(1)*B(3)
               DIR(3) = A(1)*B(2) - A(2)*B(1)
            ENDIF
            DLEN = SQRT(DOT_PRODUCT(DIR,DIR))
            IF ((SIZE_AB == ZERO) .OR. (DLEN <= EPS1*SIZE_AB)) THEN
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1826) NAME, SETID, AGRID, (GD(J),J=1,-2*ACID_L)
               WRITE(F06,1826) NAME, SETID, AGRID, (GD(J),J=1,-2*ACID_L)
               CYCLE i_do1
            ENDIF
            DO J=1,3
               F2(J) = F1(1)*DIR(J)/DLEN
            ENDDO

         ELSE IF (ACID_L /= 0) THEN                        ! Get local coord sys for this FORCE/MOMENT card
j_do12:     DO J=1,NCORD
               IF (ACID_L == CORD(J,2)) THEN
                  CORD_FND = 'Y'
                  ICID = J
                  EXIT j_do12
               ENDIF
            ENDDO j_do12
            IF (CORD_FND == 'N') THEN                      ! Coord sys ID on FORCE/MOMENT undefined
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1822) 'COORD SYSTEM ', ACID_L, NAME, SETID
               WRITE(F06,1822) 'COORD SYSTEM ', ACID_L, NAME, SETID
            ENDIF

            IF (( CORD_FND == 'N') .OR. (GRID_FND == 'N')) THEN
               CYCLE i_do1                                 ! Can't continue (GRID_ID_ROW_NUM or ACID_L not found),
!                                                            so CYCLE & read next record
            ENDIF
                                                           ! Get coord transformation matrix - local to basic
            CALL GEN_T0L ( GRID_ID_ROW_NUM, ICID, THETAD, PHID, T12 )

            CALL MATMULT_FFF ( T12, F1, 3, 3, 1, F2 )      ! Transform coords from local (forces F1) to basic (forces F2)

         ELSE                                              ! Local coord system is basic, so equate F2, F1
            DO J=1,3
               F2(J) = F1(J)
            ENDDO
         ENDIF

                                                           ! (1-c) Transform from basic to global. Note that F2 are forces in basic
         IF (IERROR == 0) THEN

            ACID_G = GRID(GRID_ID_ROW_NUM,3)               ! Get internal coord sys number for this grid's global coord sys
            ICID = -1
            IF (ACID_G /= 0) THEN                          ! We know ACID_G exists; that was verified in subr GRID_PROC
j_do14:        DO J=1,NCORD
                  IF (ACID_G == CORD(J,2)) THEN
                     ICID = J
                     EXIT j_do14
                  ENDIF
               ENDDO j_do14
                                                           ! Get coord transformation matrix - basic to global
               CALL GEN_T0L ( GRID_ID_ROW_NUM, ICID, THETAD, PHID, T12 )

               CALL MATMULT_FFF_T (T12, F2, 3, 3, 1, F1)   ! Transform coords from basic to global for F2
               DO J=1,3
                  F2(J) = F1(J)                            ! F2 has the force/moment values in global coords
               ENDDO
            ENDIF

            DO J=1,3                                       ! Reset FORMON to F2, which now has forces in global coords
               FORMON(J) = F2(J)
            ENDDO
                                                           ! (1-d) Write results to scratch file. FORMON comps arein global coords
            WRITE(SCR(1)) SETID,AGRID,ACID_G,(FORMON(J),J=1,3),NAME

         ENDIF

      ENDDO i_do1

      IF (READ_ERR > 0) THEN
         WRITE(ERR,9998) READ_ERR,LINK1I
         WRITE(ERR,9998) READ_ERR,LINK1I
         CALL OUTA_HERE ( 'Y' )                            ! Quit due to errors reading FORCE/MOMENT file
      ENDIF

      IF (IERROR > 0) THEN
         WRITE(ERR,9996) SUBR_NAME,IERROR
         WRITE(ERR,9996) SUBR_NAME,IERROR
         CALL OUTA_HERE ( 'Y' )                            ! Quit due to undefined grid and coord sys ID's
      ENDIF

! **********************************************************************************************************************************
! Now process forces in global coords into SYS_LOAD for each subcase.

      DO I=1,LLOADC                                         ! Initialize LSID, RSID arrays
         LSID(I) = 0
         RSID(I) = ZERO
      ENDDO

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: FORCE_MOM_PROC '

i_do2:DO I=1,NSUB                                          ! Loop through the S/C's

         IF (SUBLOD(I,1) == 0) THEN                        ! If no load for this S/C, CYCLE
            CYCLE i_do2
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(I,1)                             ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the FORCE or MOMENT card in file LINK1I that matches SUBLOD(I,1)
         DO J = 1,NLOAD                                    ! Then, the actual mag. will come from RSID(1) & the FORMON components
            IF (SUBLOD(I,1) == LOAD_SIDS(J,1)) THEN
k_do21:        DO K = 2,LLOADC
                  IF (LOAD_SIDS(J,K) == 0) THEN
                     EXIT k_do21
                  ELSE
                     NSID = K                              ! Note: NSID will not get larger than LLOADC
                     RSID(K) = LOAD_FACS(J,1)*LOAD_FACS(J,K)
                     LSID(K) = LOAD_SIDS(J,K)
                   ENDIF
                ENDDO k_do21
            ENDIF
         ENDDO

j_do22:  DO J=1,NFORCE                                     ! Process FORCE / MOMENT card info
                                                           ! (2-b-  i) Read a record from scratch - forces are in global coords
            READ(SCR(1),IOSTAT=IOCHK) SETID,AGRID,ACID_G,(FORMON(K),K=1,3),NAME
            IF (IOCHK /= 0) THEN
               REC_NO = J
               CALL READERR ( IOCHK, SCRFIL, MESSAG, REC_NO, OUNT )
               CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
               CALL OUTA_HERE ( 'Y' )                              ! Error reading scratch file, so quit
            ENDIF

            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
k_do221:    DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                  EXIT k_do221
               ENDIF
            ENDDO k_do221

            IF (FOUND == 'N') THEN                         ! Cycle back on J loop and read another force/moment card
               CYCLE j_do22
            ENDIF

            IF      (NAME == 'FORCE   ') THEN              ! Set component range (for loop below) based on card type
               COMP1 = 1
               COMP2 = 3
            ELSE IF (NAME == 'MOMENT  ') THEN
               COMP1 = 4
               COMP2 = 6
            ELSE
               WRITE(ERR,1516) SUBR_NAME,NAME
               WRITE(F06,1516) SUBR_NAME,NAME
               CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (not FORCE or MOMENT), so quit
            ENDIF
                                                           ! Get GRID_ID_ROW_NUM, we checked it's existence earlier
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )

            IF ((DABS(FORMON(1)) < EPS1) .AND. (DABS(FORMON(2)) < EPS1) .AND. (DABS(FORMON(3)) < EPS1)) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,1513) NAME,SETID
               IF (SUPWARN == 'N') THEN                     ! Issue warning if all force components zero
                  WRITE(F06,1513) NAME,SETID
               ENDIF
            ENDIF

!xx         CALL CALC_TDOF_ROW_NUM ( AGRID, ROW_NUM_START, 'N' )
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            K1 = 0                                         ! (2-b-iii). Put forces and moments into SYS_LOAD array
k_do222:    DO K = COMP1,COMP2
               K1 = K1+1
               IF ((K1 < 1) .OR. (K1 > 3)) THEN
                  WRITE(ERR,1514) SUBR_NAME
                  WRITE(F06,1514) SUBR_NAME
                  FATAL_ERR = FATAL_ERR + 1
                  CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
                  CALL OUTA_HERE ( 'Y' )                   ! Coding error (dim on array FORMON out of bounds), so quit
               ELSE
                  FORCEI = SCALE*FORMON(K1)
                  IF (DABS(FORCEI) < EPS1) THEN
                     CYCLE k_do222
                  ELSE
                     CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
                     ROW_NUM = ROW_NUM_START + K -1
                     GDOF = TDOF(ROW_NUM,G_SET_COL_NUM)
                     SYS_LOAD(GDOF,I) = SYS_LOAD(GDOF,I) + FORCEI
                  ENDIF
               ENDIF
            ENDDO k_do222

         ENDDO j_do22

         REWIND (SCR(1))                                   ! Need to read all of the FORCE/MOMENT records again for the next S/C

      ENDDO i_do2

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )



      RETURN

! **********************************************************************************************************************************






 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1826 FORMAT(' *ERROR  1826: THE DIRECTION OF ',A,' SET ',I8,' AT GRID ',I8,' HAS ZERO LENGTH (GRIDS ',4(I8,1X),')')

 1513 FORMAT(' *WARNING    : ',A8,1X,I8,' HAS ALL ZERO COMPONENTS')

 1514 FORMAT(' *ERROR  1514: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' K1 = ',I8,' BUT MUST BE 1 <= K1 <= 3')

 1516 FORMAT(' *ERROR  1516: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VARIABLE "NAME" SHOULD BE "FORCE" OR "MOMENT" BUT IS: ',A8)

 9998 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS READING FILE:',/,A)

 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE FORCE_MOM_PROC

   END MODULE FORCE_MOM_PROCESSING
