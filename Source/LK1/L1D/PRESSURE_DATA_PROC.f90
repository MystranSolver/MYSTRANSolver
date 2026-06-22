!--- cbeam add --- begin!
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
 
      SUBROUTINE PRESSURE_DATA_PROC
 
! Processes element pressure data. The element overall pressure and/or element G.P. pressure data was written to file L1Q when the
! element pressure Bulk Data entries were read. In this subr, records are read from L1Q and the data processed into arrays PPNT and
! PDATA (arrays used in the element generation routines)
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1Q
      USE IOUNT1, ONLY                :  WRT_ERR,                             LINK1Q
      USE IOUNT1, ONLY                :  WRT_ERR,                             L1Q_MSG
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, DATA_NAM_LEN, FATAL_ERR, JCARD_LEN, LPDAT, LLOADC,            &
                                         MPDAT_PLOAD1, MPDAT_PLOAD2, MPDAT_PLOAD4, MPLOAD4_3D_DATA, NELE, NLOAD, NPCARD,           &
                                         NPLOAD4_3D, NPDAT, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE MODEL_STUF, ONLY            :  LOAD_SIDS, LOAD_FACS, SUBLOD, PDATA, PPNT, PLOAD4_3D_DATA, PTYPE
      USE MODEL_STUF, ONLY            :  ESORT1
 
      USE PRESSURE_DATA_PROC_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PRESSURE_DATA_PROC'
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! B.D elem pressure card image read from file LINK1Q
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes
      CHARACTER( 1*BYTE)              :: EFLAG             ! Flag to decide whether a situation in subr EPPUT is an error or not
      CHARACTER( 1*BYTE)              :: FOUND             ! Indicator on whether we found something we were looking for
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters in CARD
      CHARACTER( 8*BYTE)              :: NAME              ! Card name (PLOAD1,2 or 4)
      CHARACTER( 8*BYTE)              :: PLOAD1_SCALE      ! PLOAD1 scale keyword
      CHARACTER( 8*BYTE)              :: PLOAD1_TYPE       ! PLOAD1 load direction keyword
      CHARACTER( 6*BYTE)              :: PLATE_OR_SOLID    ! 'PLATE' or 'SOLID' element designation from PLOAD4 entries
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Variable to test whether "THRU" option was used on B.D. PLOAD2 card
      CHARACTER( 8*BYTE)              :: THRU              ! ='Y' if THRU option used on TEMPRB, TEMPP1 continuation card
 
      INTEGER(LONG)                   :: EID               ! Actual element ID
      INTEGER(LONG)                   :: EID1,EID2         ! The 2 actual elem ID's in "EID1 THRU EID2" on elem press B.D. card 
      INTEGER(LONG)                   :: EL_PRES_ERR       ! Count of error messages when elements have redundant pressures
      INTEGER(LONG)                   :: EL_REDUNDANT_PRES ! Count of warning messages when elements have redundant pressures
      INTEGER(LONG)                   :: G1,G34            ! Grid numbers from a PLOAD4 entry for solid elements
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IELEM             ! Internal element number for the actual element ID, EID
      INTEGER(LONG)                   :: IERROR    = 0     ! Cum. count of errors as we read, and check cards from file LINK1Q
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: IPLOAD4_3D        ! Count of 3D elements that have PLOAD4 face pressure definition
      INTEGER(LONG)                   :: IPPNT             ! Index in array PPNT (a pointer to where press data for EID starts)
      INTEGER(LONG)                   :: LSID(LLOADC+1)    ! Array of load SID's, from a LOAD Bulk Data card, for one S/C 
      INTEGER(LONG)                   :: NFIELD            ! No. fields on PLOAD1, PLOAD2 B.D. cards that have pressure data
      INTEGER(LONG)                   :: NSID              ! Count on no. of pairs of entries on a LOAD B.D. card (<= LLOADC) 
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN  
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SETID             ! Pressure load set ID read from an elem pressure B.D. card
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file

 
      REAL(DOUBLE)                    :: SCALE             ! Scale factor from a LOAD Bulk Data card
      REAL(DOUBLE)                    :: RPDAT             ! Real pressure value read from file LINK1Q
      REAL(DOUBLE)                    :: RPDAT1            ! Real pressure value read from file LINK1Q
      REAL(DOUBLE)                    :: RPDAT2            ! 2nd real pressure value read from file LINK1Q
      REAL(DOUBLE)                    :: RSID(LLOADC+1)    ! Array of load magnitudes (for LSID set ID's) needed for one S/C
      REAL(DOUBLE)                    :: X1                ! PLOAD1 start fraction/location
      REAL(DOUBLE)                    :: X2                ! PLOAD1 end fraction/location
      REAL(DOUBLE)                    :: X1_STORE          ! Stored PLOAD1 start value (fraction or encoded length)
      REAL(DOUBLE)                    :: X2_STORE          ! Stored PLOAD1 end value (fraction or encoded length)
      LOGICAL                         :: PLOAD1_IS_LENGTH  ! True for LE/LEPR x-location input
      LOGICAL                         :: PLOAD1_PROJECTED  ! True for FRPR/LEPR projected load interpretation


! **********************************************************************************************************************************
      EL_REDUNDANT_PRES  =  0                              ! Initialize warn, err indicators for redundant elem press definition
      EL_PRES_ERR        =  0
      IPLOAD4_3D         =  0
      PLATE_OR_SOLID(1:) = ' '

      OUNT(1) = ERR                                        ! Set units for writing errors (for subr READERR)
      OUNT(2) = F06

      DO I=1,LLOADC                                         ! Initialize LSID, RSID arrays
         LSID(I) = 0
         RSID(I) = ZERO
      ENDDO

isubc:DO I=1,NSUB                                          ! Loop through the S/C's

         IF (SUBLOD(I,1) == 0) THEN                        ! If no load for this S/C, CYCLE
            CYCLE isubc
         ENDIF
                                                           ! (2-a) Generate LSID/RSID tables for this S/C.
         NSID    = 1                                       ! There is always 1 pair (more if there are LOAD B.D cards).
         LSID(1) = SUBLOD(I,1)                             ! Note: If there are no LOAD B.D. cards, LSID(1) and RSID(1) will be
         RSID(1) = ONE                                     ! for the elem pressure data in file LINK1Q that matches SUBLOD(I,1)
j_do1:   DO J = 1,NLOAD                                    ! Then, the actual mag. will come from RSID(1)
            IF (SUBLOD(I,1) == LOAD_SIDS(J,1)) THEN
k_do1:         DO K = 2,LLOADC
                  IF (LOAD_SIDS(J,K) == 0) THEN
                     EXIT k_do1  
                  ELSE
                     NSID = K                              ! Note: NSID will not get larger than LLOADC
                     RSID(K) = LOAD_FACS(J,1)*LOAD_FACS(J,K)
                     LSID(K) = LOAD_SIDS(J,K)
                   ENDIF
                ENDDO k_do1   
            ENDIF
         ENDDO j_do1   

pcards:  DO J=1,NPCARD                                     ! Process elem pressure card info
                                                           ! (2-b-  i) Read a record from scratch - forces are in global coords
            READ(L1Q,IOSTAT=IOCHK) CARD                    ! Read element pressure CARD from LINK1Q
            IF (IOCHK /= 0) THEN
               REC_NO = J + 1
               CALL READERR ( IOCHK, LINK1Q, L1Q_MSG, REC_NO, OUNT )
               IERROR = IERROR + 1
               CYCLE isubc
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )        ! Make 10 fields of 8 chars from CARD
 
            NAME = CARD(1:8)                               ! Set NFIELD based on parent CARD type
            IF      ((NAME(1:7) == 'PLOAD1 ') .OR. (NAME(1:7) == 'PLOAD1*')) THEN
               NFIELD = MPDAT_PLOAD1
            ELSE IF ((NAME(1:7) == 'PLOAD2 ') .OR. (NAME(1:7) == 'PLOAD2*')) THEN
               NFIELD = MPDAT_PLOAD2
            ELSE IF ((NAME(1:7) == 'PLOAD4 ') .OR. (NAME(1:7) == 'PLOAD4*')) THEN
               NFIELD = MPDAT_PLOAD4
            ELSE 
               CYCLE pcards                                 ! Ignore record if not for PLOAD1 or PLOAD2
            ENDIF
            IPPNT = NPDAT + 1                              ! Set index for pointer array, PPNT 
            READ(JCARD(2),'(I8)') SETID                    ! Get pressure load SID
 
            FOUND = 'N'                                    ! (2-b- ii). Scan through LSID to find set that matches SETID read.
k_do2:      DO K = 1,NSID                                  ! There is a match; we made sure all requested loads were in B.D. deck
               IF (SETID == LSID(K)) THEN                  ! We start with K = 1 to cover the case of no LOAD B.D cards
                  SCALE = RSID(K)
                  FOUND = 'Y'
                  EXIT k_do2 
               ENDIF
            ENDDO k_do2  

            IF (FOUND == 'N') THEN                         ! Cycle back on J loop and read another elem pressure card   
               CYCLE pcards   
            ENDIF

! Put pressure data into PDATA
            IF (((NAME(1:6) /= 'PLOAD1') .AND. ((NPDAT + NFIELD) > LPDAT))) THEN
               WRITE(ERR,1523) SUBR_NAME,LPDAT
               WRITE(F06,1523) SUBR_NAME,LPDAT
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (dim of array PDATA too small), so quit
            ENDIF
 
                                                           ! Set the field number where we expect to find the 1st pressure value
            IF      ((NAME(1:7) == 'PLOAD1 ') .OR. (NAME(1:7) == 'PLOAD1*')) THEN
!*** ADDED bt CODEX -- 2026-04-20 -- FOR BEAM DSB ***
               CALL LEFT_ADJ_BDFLD ( JCARD(4) )
               CALL LEFT_ADJ_BDFLD ( JCARD(5) )
               PLOAD1_TYPE  = JCARD(4)(1:8)
               PLOAD1_SCALE = JCARD(5)(1:8)
               IF (PLOAD1_TYPE(1:2) == 'MX') PLOAD1_TYPE = 'MXE'
               IF (PLOAD1_TYPE(1:2) == 'MY') PLOAD1_TYPE = 'MYE'
               IF (PLOAD1_TYPE(1:2) == 'MZ') PLOAD1_TYPE = 'MZE'
               IF (PLOAD1_SCALE == ' ') PLOAD1_SCALE = 'FR'
               PLOAD1_IS_LENGTH = ((PLOAD1_SCALE == 'LE') .OR. (PLOAD1_SCALE == 'LEPR'))
               PLOAD1_PROJECTED = ((PLOAD1_SCALE == 'FRPR') .OR. (PLOAD1_SCALE == 'LEPR'))

               READ(JCARD(3),'(I8)') EID
               CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, EID, IELEM )
               IF (IELEM == -1) THEN
                  WRITE(ERR,1530) EID, SETID
                  WRITE(F06,1530) EID, SETID
                  FATAL_ERR = FATAL_ERR + 1
                  EL_PRES_ERR = EL_PRES_ERR + 1
                  CYCLE pcards
               ENDIF

               READ(JCARD(6),'(F16.0)') X1
               READ(JCARD(7),'(F16.0)') RPDAT
               IF ((JCARD(8)(1:) == ' ') .AND. (JCARD(9)(1:) == ' ')) THEN
                  X2    = X1
                  RPDAT2 = RPDAT
               ELSE
                  READ(JCARD(8),'(F16.0)') X2
                  IF (JCARD(9)(1:) == ' ') THEN
                     RPDAT2 = RPDAT
                  ELSE
                     READ(JCARD(9),'(F16.0)') RPDAT2
                  ENDIF
               ENDIF

               IF ((PLOAD1_SCALE /= 'FR'  ) .AND. (PLOAD1_SCALE /= 'LE'  ) .AND.                                         &
                   (PLOAD1_SCALE /= 'FRPR') .AND. (PLOAD1_SCALE /= 'LEPR')) THEN
                  WRITE(ERR,1524) TRIM(PLOAD1_SCALE), EID, SETID
                  WRITE(F06,1524) TRIM(PLOAD1_SCALE), EID, SETID
                  FATAL_ERR = FATAL_ERR + 1
                  EL_PRES_ERR = EL_PRES_ERR + 1
                  CYCLE pcards
               ENDIF
               IF (.NOT.PLOAD1_IS_LENGTH) THEN
                  IF ((X1 < -1.0D-12) .OR. (X1 > ONE + 1.0D-12) .OR. (X2 < -1.0D-12) .OR. (X2 > ONE + 1.0D-12)) THEN
                     WRITE(ERR,1525) EID, SETID, X1, X2
                     WRITE(F06,1525) EID, SETID, X1, X2
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
               ELSE
                  IF ((X1 < -1.0D-12) .OR. (X2 < -1.0D-12)) THEN
                     WRITE(ERR,1525) EID, SETID, X1, X2
                     WRITE(F06,1525) EID, SETID, X1, X2
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
               ENDIF
               IF (X2 < X1 - 1.0D-12) THEN
                  WRITE(ERR,1526) EID, SETID, X1, X2
                  WRITE(F06,1526) EID, SETID, X1, X2
                  FATAL_ERR = FATAL_ERR + 1
                  EL_PRES_ERR = EL_PRES_ERR + 1
                  CYCLE pcards
               ENDIF

               X1_STORE = X1
               X2_STORE = X2
               IF (PLOAD1_IS_LENGTH) THEN
                  X1_STORE = -(X1 + TWO)
                  X2_STORE = -(X2 + TWO)
               ENDIF
               CALL PLOAD1_PUT ( IELEM, I, IPPNT, EL_REDUNDANT_PRES, EL_PRES_ERR )
               IF (IPPNT == 0) CYCLE pcards

               IF ((PLOAD1_TYPE == 'FY') .OR. (PLOAD1_TYPE == 'FYE') .OR. (PLOAD1_TYPE == 'Y')) THEN
                  IF (PDATA(IPPNT+2) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'FYE/FY/Y'
                     WRITE(F06,1528) EID, SETID, 'FYE/FY/Y'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT  ) = SCALE*RPDAT
                  PDATA(IPPNT+1) = SCALE*RPDAT2
                  PDATA(IPPNT+2) = X1_STORE
                  PDATA(IPPNT+3) = X2_STORE
                  IF ((PLOAD1_TYPE == 'FY') .OR. (PLOAD1_TYPE == 'Y')) THEN
                     PDATA(IPPNT+24) = ONE
                     IF (PLOAD1_PROJECTED) PDATA(IPPNT+24) = TWO
                  ENDIF
               ELSE IF ((PLOAD1_TYPE == 'FZ') .OR. (PLOAD1_TYPE == 'FZE') .OR. (PLOAD1_TYPE == 'Z')) THEN
                  IF (PDATA(IPPNT+6) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'FZE/FZ/Z'
                     WRITE(F06,1528) EID, SETID, 'FZE/FZ/Z'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT+4) = SCALE*RPDAT
                  PDATA(IPPNT+5) = SCALE*RPDAT2
                  PDATA(IPPNT+6) = X1_STORE
                  PDATA(IPPNT+7) = X2_STORE
                  IF ((PLOAD1_TYPE == 'FZ') .OR. (PLOAD1_TYPE == 'Z')) THEN
                     PDATA(IPPNT+25) = ONE
                     IF (PLOAD1_PROJECTED) PDATA(IPPNT+25) = TWO
                  ENDIF
               ELSE IF (PLOAD1_TYPE == 'FXE') THEN
                  IF (PDATA(IPPNT+10) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'FXE'
                     WRITE(F06,1528) EID, SETID, 'FXE'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT+8 ) = SCALE*RPDAT
                  PDATA(IPPNT+9 ) = SCALE*RPDAT2
                  PDATA(IPPNT+10) = X1_STORE
                  PDATA(IPPNT+11) = X2_STORE
               ELSE IF (PLOAD1_TYPE == 'MXE') THEN
                  IF (PDATA(IPPNT+14) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'MXE'
                     WRITE(F06,1528) EID, SETID, 'MXE'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT+12) = SCALE*RPDAT
                  PDATA(IPPNT+13) = SCALE*RPDAT2
                  PDATA(IPPNT+14) = X1_STORE
                  PDATA(IPPNT+15) = X2_STORE
               ELSE IF (PLOAD1_TYPE == 'MYE') THEN
                  IF (PDATA(IPPNT+18) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'MYE'
                     WRITE(F06,1528) EID, SETID, 'MYE'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT+16) = SCALE*RPDAT
                  PDATA(IPPNT+17) = SCALE*RPDAT2
                  PDATA(IPPNT+18) = X1_STORE
                  PDATA(IPPNT+19) = X2_STORE
               ELSE IF (PLOAD1_TYPE == 'MZE') THEN
                  IF (PDATA(IPPNT+22) >= ZERO) THEN
                     WRITE(ERR,1528) EID, SETID, 'MZE'
                     WRITE(F06,1528) EID, SETID, 'MZE'
                     FATAL_ERR = FATAL_ERR + 1
                     EL_PRES_ERR = EL_PRES_ERR + 1
                     CYCLE pcards
                  ENDIF
                  PDATA(IPPNT+20) = SCALE*RPDAT
                  PDATA(IPPNT+21) = SCALE*RPDAT2
                  PDATA(IPPNT+22) = X1_STORE
                  PDATA(IPPNT+23) = X2_STORE
               ELSE
                  WRITE(ERR,1527) TRIM(PLOAD1_TYPE), EID, SETID
                  WRITE(F06,1527) TRIM(PLOAD1_TYPE), EID, SETID
                  FATAL_ERR = FATAL_ERR + 1
                  EL_PRES_ERR = EL_PRES_ERR + 1
                  CYCLE pcards
               ENDIF
            ELSE IF ((NAME(1:7) == 'PLOAD2 ') .OR. (NAME(1:7) == 'PLOAD2*')) THEN
               IPPNT = NPDAT + 1                           ! Set index for pointer array, PPNT
               NPDAT = NPDAT + 1
               READ(JCARD(3),'(F16.0)') RPDAT
               PDATA(NPDAT) = SCALE*RPDAT               
            ELSE IF ((NAME(1:7) == 'PLOAD4 ') .OR. (NAME(1:7) == 'PLOAD4*')) THEN
               IPPNT = NPDAT + 1                           ! Set index for pointer array, PPNT
               NPDAT = NPDAT + 1
               READ(JCARD(K+3),'(F16.0)') RPDAT1
               PDATA(NPDAT) = SCALE*RPDAT1
               DO K = 5,7
                  NPDAT = NPDAT + 1
                  IF (JCARD(K)(1:) /= ' ') THEN
                     READ(JCARD(K),'(F16.0)') RPDAT
                     PDATA(NPDAT) = SCALE*RPDAT
                  ELSE
                     PDATA(NPDAT) = SCALE*RPDAT1
                  ENDIF
               ENDDO
            ENDIF

! Process EID's. First check for the 2 options on specifying elem data. For PLOAD2, either all data are EID's or THRU option is used
!                For PLOAD4, either "THRU" is used or there is only 1 EID
 
            IF      (NAME(1:6) == 'PLOAD1') THEN
               CYCLE pcards

            ELSE IF (NAME(1:6) == 'PLOAD2') THEN
               TOKEN = JCARD(5)(1:8)                       ! Only send the 1st 8 chars of this JCARD. It has been left justified
            ELSE IF (NAME(1:6) == 'PLOAD4') THEN
               TOKEN = JCARD(8)(1:8)
            ENDIF        
            CALL TOKCHK ( TOKEN, TOKTYP )

            THRU = 'N'
            IF (TOKTYP == 'THRU    ') THEN
               THRU = 'Y'
            ENDIF 
  
            IF (NAME(1:6) == 'PLOAD2') THEN

               IF (THRU == 'N') THEN
                  EFLAG = 'Y'
k_do4:            DO K=4,9
                     IF (JCARD(K) == '        ') THEN      ! See if there is any data in other JCARD fields
                        CYCLE k_do4
                     ELSE 
                        READ(JCARD(K),'(I8)') EID
                        CALL EPPUT ( SETID, EID, I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                     ENDIF
                  ENDDO k_do4
               ELSE 
                  READ(JCARD(4),'(I8)') EID1
                  READ(JCARD(6),'(I8)') EID2
                  EFLAG = 'N'
k_do5:            DO K=EID1,EID2
                     CALL EPPUT ( SETID, K,   I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  ENDDO k_do5
               ENDIF
 
            ELSE IF (NAME(1:6) == 'PLOAD4') THEN

               IF (THRU == 'N') THEN
                  READ(JCARD(3),'(I8)') EID
                  CALL EPPUT ( SETID, EID, I, IPPNT, NAME, 'Y', IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  IF (IELEM == -1) CYCLE pcards
                  IF ((JCARD(8)(1:) /= ' ') .AND. (JCARD(9)(1:) /= ' ')) THEN
                     READ(JCARD(8),'(I8)') G1
                     READ(JCARD(9),'(I9)') G34
                     PLATE_OR_SOLID = 'SOLID'
                     IPLOAD4_3D = IPLOAD4_3D + 1
                     PLOAD4_3D_DATA(IPLOAD4_3D,1) = EID
                     PLOAD4_3D_DATA(IPLOAD4_3D,2) = IELEM
                     PLOAD4_3D_DATA(IPLOAD4_3D,3) = I      ! I is the internal subcase number (loop ID is: isubc)
                     PLOAD4_3D_DATA(IPLOAD4_3D,4) = G1
                     PLOAD4_3D_DATA(IPLOAD4_3D,5) = G34
                  ELSE
                     PLATE_OR_SOLID = 'PLATE'
                  ENDIF
               ELSE
                  READ(JCARD(3),'(I8)') EID1
                  READ(JCARD(9),'(I8)') EID2
                  EFLAG = 'N'
k_do6:            DO K=EID1,EID2
                     CALL EPPUT ( SETID, K,   I, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
                  ENDDO k_do6
               ENDIF

            ENDIF

         ENDDO pcards 
 
         IF (EL_REDUNDANT_PRES > 0) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9999)
            IF (SUPWARN == 'Y') THEN
               WRITE(F06,9999)
            ENDIF
         ENDIF
 
         REWIND (L1Q)                                       ! Need to read all of the elem pressure records again for the next S/C
         READ(L1Q,IOSTAT=IOCHK) XTIME
         IF (IOCHK /= 0) THEN
            REC_NO = 1
            CALL READERR ( IOCHK, LINK1Q, L1Q_MSG, REC_NO, OUNT )
            CALL OUTA_HERE ( 'Y' )                         ! Cannot read STIME from temperature data file, so quit
         ENDIF
 
      ENDDO isubc
 
      IF ((IERROR /= 0) .OR. (EL_PRES_ERR > 0)) THEN
         WRITE(ERR,9998) IERROR+EL_PRES_ERR
         WRITE(F06,9998) IERROR+EL_PRES_ERR
         CALL OUTA_HERE ( 'Y' )                            ! Quit due to errors reading elem press file
      ENDIF
 
! **********************************************************************************************************************************
! Now finally write processed pressure data to L1Q (destroying the PLOAD card data on the records)
 
! First close and delete L1Q file
 
      CALL FILE_CLOSE ( L1Q, LINK1Q, 'DELETE' )
 
! Open L1Q for write:
 
      CALL FILE_OPEN ( L1Q, LINK1Q, OUNT, 'REPLACE', L1Q_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
 
      DATA_SET_NAME = 'PPNT'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NELE
      WRITE(L1Q) NSUB
      DO I=1,NELE
         DO J=1,NSUB
            WRITE(L1Q) PPNT(I,J)
         ENDDO
      ENDDO
 
      DATA_SET_NAME = 'PDATA'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NPDAT
      DO I=1,NPDAT
         WRITE(L1Q) PDATA(I)
      ENDDO
 
      DATA_SET_NAME = 'PTYPE'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NELE
      DO I=1,NELE
         WRITE(L1Q) PTYPE(I)
      ENDDO
 
      DATA_SET_NAME = 'PLOAD4_3D_DATA'
      WRITE(L1Q) DATA_SET_NAME
      WRITE(L1Q) NPLOAD4_3D
      DO I=1,NPLOAD4_3D
         WRITE(L1Q) (PLOAD4_3D_DATA(I,J),J=1,MPLOAD4_3D_DATA)
      ENDDO
 


      RETURN

! **********************************************************************************************************************************




 9998 FORMAT('PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS')

 9999 FORMAT(' *WARNING    : CHECK ERR OUTPUT FILE FOR WARNING MESSAGES REGARDING REDUNDANT ELEM PRESSURE DEFINITION')

 1523 FORMAT(' *ERROR  1523: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MUCH ELEMENT PRESSURE DATA. MAX IS LPDAT = ',I8)

 1524 FORMAT(' *ERROR      : PLOAD1 SCALE "',A,'" ON ELEMENT ',I8,' IN LOAD SET ',I8,' IS NOT SUPPORTED. USE FR, LE, FRPR OR LEPR')

 1525 FORMAT(' *ERROR      : PLOAD1 ON ELEMENT ',I8,' IN LOAD SET ',I8,                                                   &
                    ' HAS X1/X2 OUT OF RANGE [0,1]. INPUT WAS ',1ES12.4,', ',1ES12.4)

 1526 FORMAT(' *ERROR      : PLOAD1 ON ELEMENT ',I8,' IN LOAD SET ',I8,                                                   &
                    ' HAS X2 < X1. INPUT WAS X1=',1ES12.4,' X2=',1ES12.4)

 1527 FORMAT(' *ERROR      : PLOAD1 TYPE "',A,'" ON ELEMENT ',I8,' IN LOAD SET ',I8,                                    &
                    ' IS NOT SUPPORTED. USE FYE/FZE/FXE/MXE/MYE/MZE (legacy FY/FZ and Y/Z also accepted)')

 1528 FORMAT(' *ERROR      : PLOAD1 ON ELEMENT ',I8,' IN LOAD SET ',I8,' HAS MORE THAN ONE ',A,' entry in one subcase. ',         &
                    'CURRENT IMPLEMENTATION ALLOWS ONE ENTRY PER COMPONENT/PER ELEMENT/PER SUBCASE')

 1530 FORMAT(' *ERROR  1520: ELEMENT ',I8,' ON PLOADi ',I8,' DOES NOT EXIST OR IS OF WRONG TYPE FOR THE PRESSURE CARD:')

! ##################################################################################################################################
 
      CONTAINS

      SUBROUTINE PLOAD1_PUT ( IELEM, JSUB, IPPNT, EL_REDUNDANT_PRES, EL_PRES_ERR )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, WARN_ERR
      USE MODEL_STUF, ONLY            :  ESORT1, PDATA, PPNT, PTYPE

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: IELEM
      INTEGER(LONG), INTENT(IN)       :: JSUB
      INTEGER(LONG), INTENT(OUT)      :: IPPNT
      INTEGER(LONG), INTENT(INOUT)    :: EL_REDUNDANT_PRES
      INTEGER(LONG), INTENT(INOUT)    :: EL_PRES_ERR

      IF (PPNT(IELEM,JSUB) == 0) THEN
         IPPNT = NPDAT + 1
         IF ((NPDAT + MPDAT_PLOAD1) > LPDAT) THEN
            WRITE(ERR,1531) SUBR_NAME,LPDAT
            WRITE(F06,1531) SUBR_NAME,LPDAT
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         NPDAT = NPDAT + MPDAT_PLOAD1
         DO K=0,MPDAT_PLOAD1-1
            PDATA(IPPNT+K) = ZERO
         ENDDO
         PDATA(IPPNT+ 2) = -ONE
         PDATA(IPPNT+ 6) = -ONE
         PDATA(IPPNT+10) = -ONE
         PDATA(IPPNT+14) = -ONE
         PDATA(IPPNT+18) = -ONE
         PDATA(IPPNT+22) = -ONE
         PPNT(IELEM,JSUB) = IPPNT
         PTYPE(IELEM)     = '2'
      ELSE
         IPPNT = PPNT(IELEM,JSUB)
         IF (PTYPE(IELEM) /= '2') THEN
            WRITE(ERR,1532) ESORT1(IELEM), JSUB, PTYPE(IELEM)
            WRITE(F06,1532) ESORT1(IELEM), JSUB, PTYPE(IELEM)
            FATAL_ERR = FATAL_ERR + 1
            EL_PRES_ERR = EL_PRES_ERR + 1
            IPPNT = 0
            RETURN
         ENDIF
      ENDIF

 1531 FORMAT(' *ERROR  1523: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                             &
                    /,14X,' TOO MUCH ELEMENT PRESSURE DATA. MAX IS LPDAT = ',I8)

 1532 FORMAT(' *ERROR      : ELEMENT ',I8,' IN INTERNAL SUBCASE ',I8,                                                        &
                    ' HAS INCOMPATIBLE MIXED PLOAD TYPES. EXISTING PTYPE = ',A1)

      END SUBROUTINE PLOAD1_PUT
 
! Change log (changes following completion of Version 1.02 on 05/01/03)

! 08/13/04: (1) Add EL_PRES_ERR to arg list (like EL_REDUNDANT_PRES it gives actual error count as opposed to warning count). 
!               Then remove CALL OUTA_HERE and replace with updating EL_PRES_ERR and RETURN
!           (2) Change EL_REDUNDANT_PRES to INOUT arg (was only OUT so total amount of warn errors was not correct)
!           (3) Change arg CARD to SETID and change *ERROR  1320

! ##################################################################################################################################

      SUBROUTINE EPPUT ( SETID, EID, JSUB, IPPNT, NAME, EFLAG, IELEM, EL_REDUNDANT_PRES, EL_PRES_ERR )
 
! Element pressure routine - generates the PPNT(i,J) array and PTYPE(i) array
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NELE, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ESORT1, ETYPE, SUBLOD, PPNT, PTYPE
 
      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EPPUT'
      CHARACTER(LEN=*), INTENT(IN)    :: EFLAG             ! Flag to decide whether a situation in subr EPPUT is an error or not
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! 'PLOAD1', 'PLOAD2' or 'PLOAD4'
 
      INTEGER(LONG), INTENT(IN)       :: EID               ! Actual elem ID
      INTEGER(LONG), INTENT(IN)       :: IPPNT             ! Index in array PPNT (a pointer to where press data for EID starts)
      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Internal subcase number
      INTEGER(LONG), INTENT(IN)       :: SETID             ! Actual elem ID
      INTEGER(LONG), INTENT(OUT)      :: IELEM    ! Internal elem ID for actual elem ID EID
      INTEGER(LONG), INTENT(INOUT)    :: EL_REDUNDANT_PRES      ! Count of warning messages when elements have redundant pressures
      INTEGER(LONG), INTENT(INOUT)    :: EL_PRES_ERR       ! Count of errors where elem ID is wrong (*ERROR  1320)




! **********************************************************************************************************************************
! Initialize outputs

      IELEM = 0

! First convert EID to interal value
 
      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, EID, IELEM )

! If EID not found then: error if EFLAG = Y, or return if EFLAG /= Y 
 
      IF (IELEM == -1) THEN
         IF (EFLAG == 'Y') THEN
            WRITE(ERR,1520) EID, SETID
            WRITE(F06,1520) EID, SETID
            FATAL_ERR = FATAL_ERR + 1
            EL_PRES_ERR = EL_PRES_ERR + 1
         ENDIF
         RETURN
      ENDIF
 
! No error, so put pointer into PPNT(i,j) and give warning if the elem has had a pressure defined previously
 
      IF (PPNT(IELEM,JSUB) /= 0) THEN
         EL_REDUNDANT_PRES = EL_REDUNDANT_PRES + 1
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1521) JSUB,ESORT1(IELEM)
         IF (SUPWARN == 'N') THEN 
            WRITE(F06,1521) JSUB,ESORT1(IELEM)
         ENDIF 
      ENDIF  

      PPNT(IELEM,JSUB) = IPPNT

! Set PTYPE for this element

      IF      (NAME(1:6) == 'PLOAD1') THEN
         PTYPE(IELEM) = '2'
      ELSE IF (NAME(1:6) == 'PLOAD2') THEN
         PTYPE(IELEM) = '1'
      ELSE IF (NAME(1:6) == 'PLOAD4') THEN
         IF      ((ETYPE(IELEM)(1:4) == 'TRIA') .OR. (ETYPE(IELEM)(1:5) == 'TETRA') .OR. (ETYPE(IELEM)(1:5) == 'PENTA')) THEN 
            PTYPE(IELEM) = '3'
         ELSE IF ((ETYPE(IELEM)(1:4) == 'QUAD') .OR. (ETYPE(IELEM)(1:4) == 'HEXA')) THEN
            PTYPE(IELEM) = '4'
         ENDIF
      ENDIF
  

      RETURN

! **********************************************************************************************************************************
 1520 FORMAT(' *ERROR  1520: ELEMENT ',I8,' ON PLOADi ',I8,' DOES NOT EXIST OR IS OF WRONG TYPE FOR THE PRESSURE CARD:')

 1521 FORMAT(' *WARNING    : FOR INTERNAL SUBCASE NUMBER ',I8,' ELEMENT ',I8,' HAS PRESSURE DEFINED MORE THAN ONCE.',              &
                           ' LAST VALUE IN INPUT DECK WILL BE USED')

 1528 FORMAT(' *ERROR      : ELEMENT ',I8,' IN INTERNAL SUBCASE ',I8,' HAS INCOMPATIBLE MIXED PLOAD TYPES. EXISTING PTYPE = ',A1)

 1522 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE PRESSURE DATA ERRORS')

! **********************************************************************************************************************************
 
      END SUBROUTINE EPPUT
 
      END SUBROUTINE PRESSURE_DATA_PROC





