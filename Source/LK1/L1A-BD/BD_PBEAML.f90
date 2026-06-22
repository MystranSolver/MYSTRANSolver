! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE BD_PBEAML ( CARD, LARGE_FLD_INP )

! Processes a minimal phase-1 PBEAML Bulk Data card.
! Current support is intentionally narrow:
!   - selected section types with direct property mapping into internal PBEAM
!   - small/free-field continuation chain only
! The parsed section is mapped into the existing internal PBEAM storage so the
! current CBEAM station runtime can be reused without a separate formulation path.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE PARAMS, ONLY                :  EPSIL, SUPINFO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBEAM, NPBEAM,             &
                                          MPBEAM_STATIONS, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBEAM, PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS, RPBEAM

      USE BD_PBEAML_USE_IFs

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: MAXTOK = 512
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PBEAML'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN=LEN(CARD))        :: CARD_WORK
      CHARACTER(LEN=LEN(CARD))        :: RAW_LINE
      CHARACTER(LEN=LEN(CARD))        :: PARSE_LINE
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN=JCARD_LEN)        :: TOKENS(MAXTOK)
      CHARACTER(LEN=JCARD_LEN)        :: WORDS(32)
      CHARACTER(LEN=JCARD_LEN)        :: ID
      CHARACTER(LEN=JCARD_LEN)        :: SEC_TYPE
      CHARACTER(LEN=JCARD_LEN)        :: SOFLAG
      CHARACTER(LEN=JCARD_LEN)        :: TOKEN_WORK

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: ICOM
      INTEGER(LONG)                   :: NDIM_SEC
      INTEGER(LONG)                   :: IFIRST
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: ITOK
      INTEGER(LONG)                   :: MATERIAL_ID
      INTEGER(LONG)                   :: NWORDS
      INTEGER(LONG)                   :: NTOK
      INTEGER(LONG)                   :: PROPERTY_ID
      INTEGER(LONG)                   :: STATION_COUNT

      REAL(DOUBLE)                    :: DIMS_A(10)
      REAL(DOUBLE)                    :: DIMS_B(10)
      REAL(DOUBLE)                    :: DIMS_CUR(10)
      REAL(DOUBLE)                    :: NSM_A
      REAL(DOUBLE)                    :: NSM_B
      REAL(DOUBLE)                    :: NSM_CUR
      REAL(DOUBLE)                    :: STATION_XL

! **********************************************************************************************************************************
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      ID = JCARD(2)

      IF (LARGE_FLD_INP == 'Y') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1301) 'large-field PBEAML'
         WRITE(F06,1301) 'large-field PBEAML'
         RETURN
      ENDIF

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )
      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )

      SEC_TYPE = ADJUSTL(JCARD(5))
      IF (SEC_TYPE(1:1) == ' ') SEC_TYPE = ADJUSTL(JCARD(4))
      IF (SEC_TYPE(1:1) == ' ') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1302) ID
         WRITE(F06,1302) ID
         RETURN
      ENDIF
      CALL TO_UPPER ( SEC_TYPE )
      NDIM_SEC = GET_SECTION_NDIMS ( SEC_TYPE )
      IF (NDIM_SEC <= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1303) ID, SEC_TYPE
         WRITE(F06,1303) ID, SEC_TYPE
         RETURN
      ENDIF

      NPBEAM = NPBEAM + 1
      PBEAM_NSTATIONS(NPBEAM) = 1
      PBEAM_XL(NPBEAM,1) = ZERO
      PBEAM(NPBEAM,1) = PROPERTY_ID
      PBEAM(NPBEAM,2) = MATERIAL_ID
      PBEAM(NPBEAM,3) = 1

      DO I=1,NPBEAM-1
         IF (PROPERTY_ID == PBEAM(I,1)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1145) 'PBEAML', PROPERTY_ID
            WRITE(F06,1145) 'PBEAML', PROPERTY_ID
            RETURN
         ENDIF
      ENDDO

      DO I=1,MAXTOK
         TOKENS(I)(1:) = ' '
      ENDDO
      NTOK = 0
      CARD_WORK = CARD

collect_tokens: DO
         READ(IN1,'(A)',IOSTAT=IOCHK) RAW_LINE
         IF (IOCHK < 0) EXIT collect_tokens
         IF (IOCHK > 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1305) ID
            WRITE(F06,1305) ID
            RETURN
         ENDIF
         IFIRST = 0
         DO I=1,LEN(RAW_LINE)
            IF (RAW_LINE(I:I) /= ' ') THEN
               IFIRST = I
               EXIT
            ENDIF
         ENDDO
         IF (IFIRST == 0) CYCLE collect_tokens
         IF (RAW_LINE(IFIRST:IFIRST) == '$') CYCLE collect_tokens
         IF (RAW_LINE(IFIRST:IFIRST) /= '+') THEN
            BACKSPACE(IN1)
            EXIT collect_tokens
         ENDIF

         PARSE_LINE = RAW_LINE
         ICOM = INDEX(PARSE_LINE, '$')
         IF (ICOM > 0) PARSE_LINE(ICOM:) = ' '
         PARSE_LINE(IFIRST:IFIRST) = ' '

         CALL PARSE_CHAR_STRING ( PARSE_LINE, LEN_TRIM(PARSE_LINE), 32, JCARD_LEN, NWORDS, WORDS, IERR )
         DO I=1,NWORDS
            TOKEN_WORK = ADJUSTL(WORDS(I))
            IF (TOKEN_WORK(1:1) == ' ') CYCLE
            IF (NTOK >= MAXTOK) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1304) ID, MAXTOK
               WRITE(F06,1304) ID, MAXTOK
               RETURN
            ENDIF
            NTOK = NTOK + 1
            TOKENS(NTOK) = TOKEN_WORK
         ENDDO
      ENDDO collect_tokens

      IF (NTOK < NDIM_SEC) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1305) ID
         WRITE(F06,1305) ID
         RETURN
      ENDIF

      DIMS_A = ZERO
      DIMS_B = ZERO
      DIMS_CUR = ZERO

      DO I=1,NDIM_SEC
         READ(TOKENS(I),*,ERR=900) DIMS_A(I)
      ENDDO
      ITOK = NDIM_SEC + 1
      NSM_A = ZERO
      IF (ITOK <= NTOK) THEN
         IF (LEN_TRIM(TOKENS(ITOK)) > 0) THEN
            IF (.NOT. IS_SO_TOKEN(TOKENS(ITOK))) THEN
               READ(TOKENS(ITOK),*,ERR=900) NSM_A
               ITOK = ITOK + 1
            ENDIF
         ENDIF
      ENDIF

      CALL LOAD_SECTION_A ( SEC_TYPE, NDIM_SEC, DIMS_A, NSM_A )
      DIMS_B = DIMS_A
      NSM_B  = NSM_A
      STATION_COUNT = 1

! --- cbeam_pbeaml_constant begin --- !
! Support the standard constant-section PBEAML form where the only continuation
! fields are the section dimensions (and optional NSM) with no SO/XL station
! data. In that case the section is uniform from end A to end B.
      IF (ITOK > NTOK) THEN
         STATION_COUNT = 2
         PBEAM_NSTATIONS(NPBEAM) = STATION_COUNT
         PBEAM_XL(NPBEAM,STATION_COUNT) = 1.0D0
         CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )
         RETURN
      ENDIF
! --- cbeam_pbeaml_constant end --- !

station_parse: DO WHILE (ITOK <= NTOK)
         SOFLAG = TOKENS(ITOK)
         IF (.NOT. IS_SO_TOKEN(SOFLAG)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1306) ID, TOKENS(ITOK)
            WRITE(F06,1306) ID, TOKENS(ITOK)
            RETURN
         ENDIF
         ITOK = ITOK + 1
         IF (ITOK + NDIM_SEC > NTOK + 1) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1307) ID
            WRITE(F06,1307) ID
            RETURN
         ENDIF

         READ(TOKENS(ITOK),*,ERR=900) STATION_XL
         ITOK = ITOK + 1
         DO I=1,NDIM_SEC
            READ(TOKENS(ITOK),*,ERR=900) DIMS_CUR(I)
            ITOK = ITOK + 1
         ENDDO
         NSM_CUR = ZERO
         IF (ITOK <= NTOK) THEN
            IF (.NOT. IS_SO_TOKEN(TOKENS(ITOK))) THEN
               READ(TOKENS(ITOK),*,ERR=900) NSM_CUR
               ITOK = ITOK + 1
            ENDIF
         ENDIF

         IF (STATION_COUNT < MPBEAM_STATIONS) THEN
            STATION_COUNT = STATION_COUNT + 1
            PBEAM_XL(NPBEAM,STATION_COUNT) = STATION_XL
            PBEAM_NSTATIONS(NPBEAM) = STATION_COUNT
            CALL STORE_SECTION_PROPS ( STATION_COUNT, SEC_TYPE, NDIM_SEC, DIMS_CUR, NSM_CUR )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,1197) 'PBEAML', ID, MPBEAM_STATIONS
            WRITE(F06,1197) 'PBEAML', ID, MPBEAM_STATIONS
         ENDIF

         DIMS_B = DIMS_CUR
         NSM_B  = NSM_CUR
      ENDDO station_parse

      IF (PBEAM_NSTATIONS(NPBEAM) <= 1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1308) ID
         WRITE(F06,1308) ID
         RETURN
      ENDIF

      CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )

      RETURN

  900 CONTINUE
      FATAL_ERR = FATAL_ERR + 1
      WRITE(ERR,1309) ID
      WRITE(F06,1309) ID
      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)
 1197 FORMAT(' *WARNING 1197: ',A,' ID = ',A,' HAS MORE THAN ',I8,' STORED STATIONS. EXTRA x/L VALUES WILL BE IGNORED')
 1301 FORMAT(' *ERROR  1301: ',A,' IS NOT SUPPORTED IN THIS PHASE-1 PBEAML IMPLEMENTATION')
 1302 FORMAT(' *ERROR  1302: PBEAML ID = ',A,' IS MISSING A SECTION TYPE')
 1303 FORMAT(' *ERROR  1303: PBEAML ID = ',A,' HAS UNSUPPORTED SECTION TYPE "',A,'" IN THIS PHASE-1 IMPLEMENTATION')
 1304 FORMAT(' *ERROR  1304: PBEAML ID = ',A,' EXCEEDED TOKEN BUFFER OF ',I8,' CONTINUATION TOKENS')
 1305 FORMAT(' *ERROR  1305: PBEAML ID = ',A,' DOES NOT CONTAIN ENOUGH DIMENSION DATA')
 1306 FORMAT(' *ERROR  1306: PBEAML ID = ',A,' EXPECTED SO TOKEN YES/YESA/NO BUT FOUND "',A,'"')
 1307 FORMAT(' *ERROR  1307: PBEAML ID = ',A,' HAS AN INCOMPLETE STATION DEFINITION')
 1308 FORMAT(' *ERROR  1308: PBEAML ID = ',A,' DID NOT DEFINE ANY STATION BEYOND END A')
 1309 FORMAT(' *ERROR  1309: PBEAML ID = ',A,' HAS NONNUMERIC TYPE=I DIMENSION DATA')
 1310 FORMAT(' *ERROR  1310: PBEAML ID = ',A,' CURRENTLY REQUIRES FREE-FIELD COMMA INPUT IN THIS PHASE-1 IMPLEMENTATION')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      LOGICAL FUNCTION IS_SO_TOKEN ( TOKEN )

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      CHARACTER(LEN=JCARD_LEN)     :: TOKEN_UP

      TOKEN_UP = TOKEN
      CALL TO_UPPER ( TOKEN_UP )
      IF ((TOKEN_UP(1:3) == 'YES') .OR. (TOKEN_UP(1:4) == 'YESA') .OR. (TOKEN_UP(1:2) == 'NO')) THEN
         IS_SO_TOKEN = .TRUE.
      ELSE
         IS_SO_TOKEN = .FALSE.
      ENDIF

      END FUNCTION IS_SO_TOKEN

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION GET_SECTION_NDIMS ( SEC )

      CHARACTER(LEN=*), INTENT(IN) :: SEC

      GET_SECTION_NDIMS = -1
      IF (SEC(1:5) == 'TUBE2') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'ROD ') THEN
         GET_SECTION_NDIMS = 1
      ELSE IF (SEC(1:4) == 'TUBE') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'BAR ') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'BOX ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'H   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'CHAN') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'T   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'L   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'I   ') THEN
         GET_SECTION_NDIMS = 6
      ENDIF

      END FUNCTION GET_SECTION_NDIMS

! ##################################################################################################################################

      SUBROUTINE LOAD_SECTION_A ( SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)             :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)             :: STRE(8)
      REAL(DOUBLE)             :: K1, K2, YC, ZC, YS, ZS, IWARP
      INTEGER(LONG)            :: IERR_LOC

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      RPBEAM(NPBEAM, 1) = AREA
      RPBEAM(NPBEAM, 2) = I1
      RPBEAM(NPBEAM, 3) = I2
      RPBEAM(NPBEAM, 4) = I12
      RPBEAM(NPBEAM, 5) = JTOR
      RPBEAM(NPBEAM, 6) = NSM
      DO I=1,8
         RPBEAM(NPBEAM,6+I) = STRE(I)
      ENDDO
      RPBEAM(NPBEAM,30) = K1
      RPBEAM(NPBEAM,31) = K2
      RPBEAM(NPBEAM,36) = IWARP
      RPBEAM(NPBEAM,38) = YS
      RPBEAM(NPBEAM,39) = ZS
      RPBEAM(NPBEAM,42) = YC
      RPBEAM(NPBEAM,43) = ZC
      CALL STORE_SECTION_PROPS ( 1, SEC_TYPE_IN, NDIM_IN, DIMS, NSM )
      CALL CHECK_BAR_MOIs ( 'PBEAML', ID, RPBEAM(NPBEAM,2), RPBEAM(NPBEAM,3), RPBEAM(NPBEAM,4), IERR_LOC )
      IF (IERR_LOC /= 0) FATAL_ERR = FATAL_ERR + 1

      END SUBROUTINE LOAD_SECTION_A

! ##################################################################################################################################

      SUBROUTINE LOAD_SECTION_B ( SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)             :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)             :: STRE(8)
      REAL(DOUBLE)             :: K1, K2, YC, ZC, YS, ZS, IWARP
      INTEGER(LONG)            :: IERR_LOC

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      RPBEAM(NPBEAM,15) = PBEAM_XL(NPBEAM,PBEAM_NSTATIONS(NPBEAM))
      RPBEAM(NPBEAM,16) = AREA
      RPBEAM(NPBEAM,17) = I1
      RPBEAM(NPBEAM,18) = I2
      RPBEAM(NPBEAM,19) = I12
      RPBEAM(NPBEAM,20) = JTOR
      RPBEAM(NPBEAM,21) = NSM
      DO I=1,8
         RPBEAM(NPBEAM,21+I) = STRE(I)
      ENDDO
      RPBEAM(NPBEAM,30) = K1
      RPBEAM(NPBEAM,31) = K2
      RPBEAM(NPBEAM,37) = IWARP
      RPBEAM(NPBEAM,40) = YS
      RPBEAM(NPBEAM,41) = ZS
      RPBEAM(NPBEAM,44) = YC
      RPBEAM(NPBEAM,45) = ZC
      CALL STORE_SECTION_PROPS ( PBEAM_NSTATIONS(NPBEAM), SEC_TYPE_IN, NDIM_IN, DIMS, NSM )
      CALL CHECK_BAR_MOIs ( 'PBEAML', ID, RPBEAM(NPBEAM,17), RPBEAM(NPBEAM,18), RPBEAM(NPBEAM,19), IERR_LOC )
      IF (IERR_LOC /= 0) FATAL_ERR = FATAL_ERR + 1

      END SUBROUTINE LOAD_SECTION_B

! ##################################################################################################################################

      SUBROUTINE STORE_SECTION_PROPS ( ISTA_IN, SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      INTEGER(LONG), INTENT(IN)    :: ISTA_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)                 :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)                 :: K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      IF ((ISTA_IN < 1) .OR. (ISTA_IN > MPBEAM_STATIONS)) RETURN

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      PBEAM_RPROPS(NPBEAM,ISTA_IN,1) = AREA
      PBEAM_RPROPS(NPBEAM,ISTA_IN,2) = I1
      PBEAM_RPROPS(NPBEAM,ISTA_IN,3) = I2
      PBEAM_RPROPS(NPBEAM,ISTA_IN,4) = I12
      PBEAM_RPROPS(NPBEAM,ISTA_IN,5) = JTOR
      PBEAM_RPROPS(NPBEAM,ISTA_IN,6) = NSM

      CALL WRITE_PBEAML_CONVERTED_PBEAM_DEBUG ( 'S', PBEAM_XL(NPBEAM,ISTA_IN), SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12,   &
                                                JTOR, NSM, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      END SUBROUTINE STORE_SECTION_PROPS

! ##################################################################################################################################

      SUBROUTINE WRITE_PBEAML_CONVERTED_PBEAM_DEBUG ( WHICH, XL, SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, NSM, K1, K2, &
                                                      YC, ZC, YS, ZS, IWARP, STRE )

      CHARACTER(LEN=*), INTENT(IN) :: WHICH
      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: XL, DIMS(10), AREA, I1, I2, I12, JTOR, NSM, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)
      INTEGER(LONG)                :: IDIM
      CHARACTER(LEN=256)           :: LINE

      WRITE(F06,'(A)') ' '
      WRITE(F06,'(A)') '*** PBEAML CONVERTED PBEAM DEBUG *********************************************'
      WRITE(F06,'(A,I0)') '  Property ID      : ', PROPERTY_ID
      WRITE(F06,'(A,I0)') '  Material ID      : ', PBEAM(NPBEAM,2)
      WRITE(F06,'(A,A)')  '  PBEAML section   : ', TRIM(SEC_TYPE_IN)
      IF (WHICH == 'A') THEN
         WRITE(F06,'(A)') '  Snapshot         : End A / station 1'
      ELSE IF (WHICH == 'B') THEN
         CALL WRITE_LABEL_VALUE ( '  Snapshot         : End B at x/L = ', XL )
      ELSE
         CALL WRITE_LABEL_VALUE ( '  Snapshot         : Intermediate station at x/L = ', XL )
      ENDIF
      WRITE(F06,'(A)') '  Input dimensions :'
      DO IDIM=1,NDIM_IN
         WRITE(LINE,'(A,I0,A,A)') '    DIM(', IDIM, ') = ', TRIM(FMT_REAL_SHORT(DIMS(IDIM)))
         WRITE(F06,'(A)') TRIM(LINE)
      ENDDO
      WRITE(F06,'(A)') '  Converted PBEAM fields :'
      CALL WRITE_LABEL_VALUE ( '    Area           = ', AREA )
      CALL WRITE_LABEL_VALUE ( '    I1             = ', I1 )
      CALL WRITE_LABEL_VALUE ( '    I2             = ', I2 )
      CALL WRITE_LABEL_VALUE ( '    I12            = ', I12 )
      CALL WRITE_LABEL_VALUE ( '    J              = ', JTOR )
      CALL WRITE_LABEL_VALUE ( '    NSM            = ', NSM )
      CALL WRITE_LABEL_VALUE ( '    K1 shear       = ', K1 )
      CALL WRITE_LABEL_VALUE ( '    K2 shear       = ', K2 )
      CALL WRITE_LABEL_VALUE ( '    CW             = ', IWARP )
      CALL WRITE_LABEL_VALUE ( '    Neutral axis Y = ', YC )
      CALL WRITE_LABEL_VALUE ( '    Neutral axis Z = ', ZC )
      CALL WRITE_LABEL_VALUE ( '    Shear center Y = ', YS )
      CALL WRITE_LABEL_VALUE ( '    Shear center Z = ', ZS )
      WRITE(F06,'(A)') '  Stress recovery points :'
      CALL WRITE_POINT_VALUE ( '    C = ', STRE(1), STRE(2) )
      CALL WRITE_POINT_VALUE ( '    D = ', STRE(3), STRE(4) )
      CALL WRITE_POINT_VALUE ( '    E = ', STRE(5), STRE(6) )
      CALL WRITE_POINT_VALUE ( '    F = ', STRE(7), STRE(8) )
      WRITE(F06,'(A)') '  NX-style converted PBEAM card image:'
      WRITE(F06,'(A)') '     THE USER SUPPLIED PBEAML BULK DATA ENTRY IS REPRESENTED INTERNALLY AS:'
      WRITE(LINE,'(A,I0,1X,I0)') '  PBEAM      ', PROPERTY_ID, PBEAM(NPBEAM,2)
      CALL WRITE_VALUE_LIST_6 ( TRIM(LINE), AREA, I1, I2, I12, JTOR, NSM )
      CALL WRITE_VALUE_LIST_8 ( '              ', STRE(1), STRE(2), STRE(3), STRE(4), STRE(5), STRE(6), STRE(7), STRE(8) )
      CALL WRITE_VALUE_LIST_8 ( '              ', K1, K2, 0.0D0, 0.0D0, 0.0D0, 0.0D0, IWARP, IWARP )
      CALL WRITE_VALUE_LIST_8 ( '              ', YS, ZS, YS, ZS, YC, ZC, YC, ZC )
      WRITE(F06,'(A)') '***************************************************************************'

      END SUBROUTINE WRITE_PBEAML_CONVERTED_PBEAM_DEBUG

! ##################################################################################################################################

      CHARACTER(LEN=24) FUNCTION FMT_REAL_SHORT ( VALUE )

      REAL(DOUBLE), INTENT(IN) :: VALUE
      CHARACTER(LEN=32)        :: BUFFER
      CHARACTER(LEN=24)        :: MANT
      CHARACTER(LEN=8)         :: EXPSTR
      INTEGER(LONG)            :: EPOS, IEND

      IF (DABS(VALUE) <= 1.0D-12) THEN
         FMT_REAL_SHORT = '0.0'
         RETURN
      ENDIF

      WRITE(BUFFER,'(ES16.8E2)') VALUE
      BUFFER = ADJUSTL(BUFFER)
      EPOS = INDEX(BUFFER,'E')
      IF (EPOS <= 0) THEN
         FMT_REAL_SHORT = TRIM(BUFFER)
         RETURN
      ENDIF

      MANT = BUFFER(1:EPOS-1)
      EXPSTR = BUFFER(EPOS:)
      IEND = LEN_TRIM(MANT)

      DO WHILE (IEND > 1)
         IF (MANT(IEND:IEND) /= '0') EXIT
         IEND = IEND - 1
      ENDDO

      IF (MANT(IEND:IEND) == '.') THEN
         MANT(IEND+1:IEND+1) = '0'
         IEND = IEND + 1
      ENDIF
      MANT = MANT(1:IEND)

      IF ((TRIM(EXPSTR) == 'E+00') .OR. (TRIM(EXPSTR) == 'E-00')) THEN
         FMT_REAL_SHORT = TRIM(MANT)
      ELSE
         FMT_REAL_SHORT = TRIM(MANT)//TRIM(EXPSTR)
      ENDIF

      END FUNCTION FMT_REAL_SHORT

! ##################################################################################################################################

      SUBROUTINE WRITE_LABEL_VALUE ( LABEL, VALUE )

      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      REAL(DOUBLE), INTENT(IN)     :: VALUE

      WRITE(F06,'(A,1X,A)') TRIM(LABEL), TRIM(FMT_REAL_SHORT(VALUE))

      END SUBROUTINE WRITE_LABEL_VALUE

! ##################################################################################################################################

      SUBROUTINE WRITE_POINT_VALUE ( LABEL, VALUE1, VALUE2 )

      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      REAL(DOUBLE), INTENT(IN)     :: VALUE1, VALUE2
      CHARACTER(LEN=128)           :: LINE

      LINE = LABEL//'('//TRIM(FMT_REAL_SHORT(VALUE1))//', '//TRIM(FMT_REAL_SHORT(VALUE2))//')'
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_POINT_VALUE

! ##################################################################################################################################

      SUBROUTINE WRITE_VALUE_LIST_6 ( PREFIX, V1, V2, V3, V4, V5, V6 )

      CHARACTER(LEN=*), INTENT(IN) :: PREFIX
      REAL(DOUBLE), INTENT(IN)     :: V1, V2, V3, V4, V5, V6
      CHARACTER(LEN=256)           :: LINE

      LINE = PREFIX//' '//TRIM(FMT_REAL_SHORT(V1))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V2))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V3))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V4))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V5))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V6))
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_VALUE_LIST_6

! ##################################################################################################################################

      SUBROUTINE WRITE_VALUE_LIST_8 ( PREFIX, V1, V2, V3, V4, V5, V6, V7, V8 )

      CHARACTER(LEN=*), INTENT(IN) :: PREFIX
      REAL(DOUBLE), INTENT(IN)     :: V1, V2, V3, V4, V5, V6, V7, V8
      CHARACTER(LEN=256)           :: LINE

      LINE = PREFIX//' '//TRIM(FMT_REAL_SHORT(V1))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V2))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V3))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V4))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V5))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V6))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V7))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V8))
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_VALUE_LIST_8

! ##################################################################################################################################

      SUBROUTINE CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      CHARACTER(LEN=*), INTENT(IN)  :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)     :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)      :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT)     :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      CHARACTER(LEN=JCARD_LEN) :: SEC

      SEC = SEC_TYPE_IN
      CALL TO_UPPER ( SEC )

      AREA  = ZERO ; I1    = ZERO ; I2    = ZERO ; I12  = ZERO ; JTOR = ZERO
      K1    = 1.D0 ; K2    = 1.D0 ; YC    = ZERO ; ZC   = ZERO ; YS   = ZERO ; ZS = ZERO ; IWARP = ZERO
      STRE  = ZERO

      IF (SEC(1:5) == 'TUBE2') THEN
         CALL CALC_TUBE_SECTION ( DIMS, .TRUE., AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'I   ') THEN
         CALL CALC_I_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      ELSE IF (SEC(1:4) == 'ROD ') THEN
         CALL CALC_ROD_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'TUBE') THEN
         CALL CALC_TUBE_SECTION ( DIMS, .FALSE., AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'BAR ') THEN
         CALL CALC_BAR_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'BOX ') THEN
         CALL CALC_BOX_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'H   ') THEN
         CALL CALC_H_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      ELSE IF (SEC(1:4) == 'CHAN') THEN
         CALL CALC_CHAN_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE )
      ELSE IF (SEC(1:4) == 'T   ') THEN
         CALL CALC_T_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      ELSE IF (SEC(1:4) == 'L   ') THEN
         CALL CALC_L_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      ELSE
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      END SUBROUTINE CALC_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_I_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      REAL(DOUBLE) :: H, B1, B2, TW, T1, T2
      REAL(DOUBLE) :: A1, A2, A3
      REAL(DOUBLE) :: Z1, Z2, Z3, ZBAR, HWEB

      H  = DIMS(1)
      B1 = DIMS(2)
      B2 = DIMS(3)
      TW = DIMS(4)
      T1 = DIMS(5)
      T2 = DIMS(6)

      HWEB = H - T1 - T2
      IF (HWEB < ZERO) HWEB = ZERO

      A1 = B1*T1
      A2 = HWEB*TW
      A3 = B2*T2
      AREA = A1 + A2 + A3

      Z1 = 0.5D0*T2
      Z2 = T2 + 0.5D0*HWEB
      Z3 = H - 0.5D0*T1
      IF (AREA > ZERO) THEN
         ZBAR = (A1*Z3 + A2*Z2 + A3*Z1)/AREA
      ELSE
         ZBAR = ZERO
      ENDIF

      YC = ZERO
      ZC = ZBAR - 0.5D0*H
      YS = ZERO
      ZS = ZERO
      IWARP = ZERO

      I1 = (B1*T1**3)/12.D0 + A1*(Z3-ZBAR)**2                                                  &
         + (TW*HWEB**3)/12.D0 + A2*(Z2-ZBAR)**2                                                &
         + (B2*T2**3)/12.D0 + A3*(Z1-ZBAR)**2

      I2 = (T1*B1**3)/12.D0 + (HWEB*TW**3)/12.D0 + (T2*B2**3)/12.D0
      I12 = ZERO
      JTOR = (B1*T1**3 + B2*T2**3 + HWEB*TW**3)/3.D0
      K1 = 0.5D0
      K2 = 0.5D0

      STRE(1) =  0.5D0*B1
      STRE(2) =  0.5D0*H - ZBAR
      STRE(3) = -0.5D0*B1
      STRE(4) =  0.5D0*H - ZBAR
      STRE(5) = -0.5D0*B2
      STRE(6) = -0.5D0*H - ZBAR
      STRE(7) =  0.5D0*B2
      STRE(8) = -0.5D0*H - ZBAR

      END SUBROUTINE CALC_I_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_H_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

! Temporary phase-1 support: treat H as an I-section rotated 90 degrees with equal flange thicknesses.
! NX reverse-engineering for H is still thin, so keep this mapping conservative and explicit.
      REAL(DOUBLE) :: IDIMS(10)

      IDIMS = ZERO
      IDIMS(1) = DIMS(1)
      IDIMS(2) = DIMS(3)
      IDIMS(3) = DIMS(3)
      IDIMS(4) = DIMS(4)
      IDIMS(5) = DIMS(2)
      IDIMS(6) = DIMS(2)

      CALL CALC_I_SECTION ( IDIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      END SUBROUTINE CALC_H_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_ROD_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: R, PI
      PI = 4.D0*DATAN(1.D0)
      R = DIMS(1)
      AREA = PI*R**2
      I1 = PI*R**4/4.D0
      I2 = I1
      I12 = ZERO
      JTOR = PI*R**4/2.D0
      K1 = 0.9D0 ; K2 = 0.9D0
      STRE = (/ R, ZERO, ZERO, R, -R, ZERO, ZERO, -R /)
      END SUBROUTINE CALC_ROD_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_TUBE_SECTION ( DIMS, USE_THICKNESS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      LOGICAL, INTENT(IN)       :: USE_THICKNESS
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: RO, RI, T, PI
      PI = 4.D0*DATAN(1.D0)
      IF (USE_THICKNESS) THEN
         RO = 0.5D0*DIMS(1)
         T  = DIMS(2)
         RI = RO - T
      ELSE
         RO = DIMS(1)
         RI = DIMS(2)
      ENDIF
      IF (RI < ZERO) RI = ZERO
      AREA = PI*(RO**2 - RI**2)
      I1 = PI*(RO**4 - RI**4)/4.D0
      I2 = I1
      I12 = ZERO
      JTOR = PI*(RO**4 - RI**4)/2.D0
      K1 = 0.5D0 ; K2 = 0.5D0
      STRE = (/ RO, ZERO, ZERO, RO, -RO, ZERO, ZERO, -RO /)
      END SUBROUTINE CALC_TUBE_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_BAR_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: H, B, A, BT, RATIO
      B = DIMS(1)
      H = DIMS(2)
      AREA = B*H
      I1 = B*H**3/12.D0
      I2 = H*B**3/12.D0
      I12 = ZERO
      A = DMAX1(H,B)/2.D0
      BT = DMIN1(H,B)/2.D0
      IF (DABS(A-BT) < 1.D-10) THEN
         JTOR = 2.25D0*A**4
      ELSE
         RATIO = BT/A
         JTOR = A*BT**3*(16.D0/3.D0 - 3.36D0*RATIO*(1.D0 - RATIO**4/12.D0))
      ENDIF
      K1 = 5.D0/6.D0 ; K2 = 5.D0/6.D0
      STRE = (/ 0.5D0*B, 0.5D0*H, -0.5D0*B, 0.5D0*H, -0.5D0*B, -0.5D0*H, 0.5D0*B, -0.5D0*H /)
      END SUBROUTINE CALC_BAR_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_BOX_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: H, B, T1, T2, HI, BI, AM, PERIM, TAVG
      H = DIMS(1); B = DIMS(2); T1 = DIMS(3); T2 = DIMS(4)
      AREA = 2.D0*(B*T1 + (H - 2.D0*T1)*T2)
      HI = H - 2.D0*T1
      BI = B - 2.D0*T2
      IF (HI < ZERO) HI = ZERO
      IF (BI < ZERO) BI = ZERO
      I1 = B*H**3/12.D0 - BI*HI**3/12.D0
      I2 = H*B**3/12.D0 - HI*BI**3/12.D0
      I12 = ZERO
      AM = (B - T2)*(H - T1)
      PERIM = 2.D0*((B - T2) + (H - T1))
      TAVG = ZERO
      IF ((B/T2 + H/T1) > ZERO) TAVG = PERIM/(2.D0*(B/T2 + H/T1))
      IF ((PERIM > ZERO) .AND. (TAVG > ZERO)) THEN
         JTOR = 4.D0*AM**2/(PERIM/TAVG)
      ELSE
         JTOR = ZERO
      ENDIF
      K1 = 0.5D0 ; K2 = 0.5D0
      STRE = (/ 0.5D0*B, 0.5D0*H, -0.5D0*B, 0.5D0*H, -0.5D0*B, -0.5D0*H, 0.5D0*B, -0.5D0*H /)
      END SUBROUTINE CALC_BOX_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_CHAN_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE(8)
      REAL(DOUBLE) :: B, H, TW, TF, AWEB, AFLG, AWEB_EFF, AFLG_EFF, YWEB, YFLG, ZBAR, E
      B  = DIMS(1)
      H  = DIMS(2)
      TW = DIMS(3)
      TF = DIMS(4)
      AWEB = TW*H
      AFLG = 2.D0*B*TF
      AREA = AWEB + AFLG
      YWEB = H/2.D0
      YFLG = H - TF/2.D0
      ZBAR = ZERO
      IF (AREA > ZERO) ZBAR = (AWEB*YWEB + AFLG*YFLG)/AREA
      YC = ZERO
      ZC = ZBAR - H/2.D0
      I1 = TW*H**3/12.D0 + AWEB*(YWEB - ZBAR)**2 + 2.D0*(B*TF**3/12.D0 + B*TF*(YFLG - ZBAR)**2)
      I2 = H*TW**3/12.D0 + 2.D0*TF*B**3/12.D0
      I12 = ZERO
      JTOR = (2.D0*B*TF**3 + H*TW**3)/3.D0
      AWEB_EFF = TW*(H - 2.D0*TF)
      IF (AWEB_EFF < ZERO) AWEB_EFF = ZERO
      AFLG_EFF = ZERO
      IF (B > ZERO) AFLG_EFF = AFLG*(1.D0 - TF/(H + B))
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = AWEB_EFF/AREA
         K2 = AFLG_EFF/AREA
      ENDIF
      E = ZERO
      IF (I2 > ZERO) E = B**2*TF*H**2*TW/(4.D0*I2)
      YS = E
      ZS = ZERO
      STRE = (/ B, 0.5D0*H, ZERO, 0.5D0*H, ZERO, -0.5D0*H, B, -0.5D0*H /)
      END SUBROUTINE CALC_CHAN_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_T_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE(8)
      REAL(DOUBLE) :: B, H, TF, TW, AWEB, AFLG, YWEB, YFLG, ZBAR
      B  = DIMS(1)
      H  = DIMS(2)
      TF = DIMS(3)
      TW = DIMS(4)
      AWEB = TW*(H - TF)
      AFLG = B*TF
      AREA = AWEB + AFLG
      YWEB = 0.5D0*(H - TF)
      YFLG = H - 0.5D0*TF
      ZBAR = ZERO
      IF (AREA > ZERO) ZBAR = (AWEB*YWEB + AFLG*YFLG)/AREA
      YC = ZERO
      ZC = ZBAR - H/2.D0
      I1 = TW*(H - TF)**3/12.D0 + AWEB*(YWEB - ZBAR)**2 + B*TF**3/12.D0 + AFLG*(YFLG - ZBAR)**2
      I2 = (H - TF)*TW**3/12.D0 + TF*B**3/12.D0
      I12 = ZERO
      JTOR = (B*TF**3 + (H - TF)*TW**3)/3.D0
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = AWEB/AREA
         K2 = AFLG/AREA
      ENDIF
      STRE = (/ 0.5D0*TW, ZERO, 0.5D0*TW, 0.5D0*H, -B + 0.5D0*TW, ZERO, 0.5D0*TW, -0.5D0*H /)
      END SUBROUTINE CALC_T_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_L_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE(8)
      REAL(DOUBLE) :: B, H, TF, TW, A1, A2, Y1, Z1, Y2, Z2, YBAR, ZBAR
      B  = DIMS(1)
      H  = DIMS(2)
      TF = DIMS(3)
      TW = DIMS(4)
      A1 = (H - TF)*TW
      A2 = (B - TW)*TF
      AREA = A1 + A2
      Y1 = 0.5D0*(H - TF) ; Z1 = 0.5D0*TW
      Y2 = 0.5D0*TF ; Z2 = TW + 0.5D0*(B - TW)
      YBAR = ZERO ; ZBAR = ZERO
      IF (AREA > ZERO) THEN
         YBAR = (A1*Y1 + A2*Y2)/AREA
         ZBAR = (A1*Z1 + A2*Z2)/AREA
      ENDIF
      YC = YBAR - H/2.D0
      ZC = ZBAR - B/2.D0
      I1 = TW*(H - TF)**3/12.D0 + A1*(Y1 - YBAR)**2 + (B - TW)*TF**3/12.D0 + A2*(Y2 - YBAR)**2
      I2 = (H - TF)*TW**3/12.D0 + A1*(Z1 - ZBAR)**2 + TF*(B - TW)**3/12.D0 + A2*(Z2 - ZBAR)**2
      I12 = A1*(Y1 - YBAR)*(Z1 - ZBAR) + A2*(Y2 - YBAR)*(Z2 - ZBAR)
      JTOR = ((H - TF)*TW**3 + (B - TW)*TF**3)/3.D0
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = A1/AREA
         K2 = A2/AREA
      ENDIF
      STRE = (/ B - ZBAR, 0.5D0*TF, -ZBAR, H - YBAR, -ZBAR, -YBAR, B - ZBAR, -YBAR /)
      END SUBROUTINE CALC_L_SECTION

! ##################################################################################################################################

      SUBROUTINE TO_UPPER ( STR )

      CHARACTER(LEN=*), INTENT(INOUT) :: STR
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: ICODE

      DO K=1,LEN(STR)
         ICODE = IACHAR(STR(K:K))
         IF ((ICODE >= IACHAR('a')) .AND. (ICODE <= IACHAR('z'))) THEN
            STR(K:K) = ACHAR(ICODE - 32)
         ENDIF
      ENDDO

      END SUBROUTINE TO_UPPER

! ##################################################################################################################################

      END SUBROUTINE BD_PBEAML
