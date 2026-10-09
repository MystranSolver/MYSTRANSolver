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

   MODULE MATERIAL_CARDS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_MAT1, BD_MAT2, BD_MAT8, BD_MAT9

   CONTAINS

      SUBROUTINE BD_MAT1 ( CARD, LARGE_FLD_INP )

! Processes MAT1 Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMATL, MRMATLC, NMATL, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL, SUPINFO, SUPWARN
      USE MODEL_STUF, ONLY            :  MATL, RMATL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MAT1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=LEN(CARD))        :: CARDP             ! Parent card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_E           ! The field that contains E (Young's modulus)
      CHARACTER(LEN(JCARD))           :: JCARD_G           ! The field that contains G (shear modulus)
      CHARACTER(LEN(JCARD))           :: JCARD_NU          ! The field that contains NU (Poisson's ratio)

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATL_ID   = 0     ! The ID for this MAT1 (field 2)


      REAL(DOUBLE)                    :: R8INP              ! A real input value read



! **********************************************************************************************************************************
! MAT1 Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Material ID     MATL(nmatl,1)
!           Material type   MATL(nmatl,2) (e.g. 1 indicates MAT1 entry)
!    3      Elas. Mod.- E  RMATL(nmatl,1)
!    4      Shear.Mod.- G  RMATL(nmatl,2)
!    5      Poiss.Rat.-NU  RMATL(nmatl,3)
!    6      Mass Den.- RHO RMATL(nmatl,4)
!    7      Exp. Coef.- A  RMATL(nmatl,5)
!    8      Ref. Temp.- T  RMATL(nmatl,6)
!    9      Damp.Coef.-GE  RMATL(nmatl,7)
! on optional second card:
!    2      Ten. Lim. -ST  RMATL(nmatl,8)
!    3      Com.  Lim.-SC  RMATL(nmatl,9)
!    4      Shr. Lim. -SS  RMATL(nmatl,10)

      CARDP = CARD

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for too many MATL entries

      NMATL = NMATL+1
      IF (NMATL > LMATL) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMATL
         WRITE(F06,1163) SUBR_NAME,JCARD(1),LMATL
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

      CALL I4FLD ( JCARD(2), JF(2), MATL_ID )              ! Check for duplicate ID number
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NMATL-1
            IF (MATL_ID == MATL(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),MATL_ID
               WRITE(F06,1145) JCARD(1),MATL_ID
               EXIT
            ENDIF
         ENDDO
         MATL(NMATL,1) = MATL_ID
         MATL(NMATL,2) = 1                                 ! Type is 1 for MAT1 card
      ENDIF

      DO J = 1,7                                           ! Read 7 fields of real data on the parent card
         R8INP = ZERO
         CALL R8FLD ( JCARD(J+2), JF(J+2), R8INP )
         IF (IERRFL(J+2) == 'N') THEN
            RMATL(NMATL,J) = R8INP
         ENDIF
      ENDDO

      JCARD_E  = JCARD(3)
      JCARD_G  = JCARD(4)
      JCARD_NU = JCARD(5)

! Check on reasonable E, G, NU and calculate values for fields that were left blank

      IF ((IERRFL(3) == 'N') .AND. (IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) THEN
         CALL MAT1_VALUE_CHECK
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      DO J=8,MRMATLC                                       ! Null optional data on 2nd, 3rd cards:
         RMATL(NMATL,J) = ZERO
      ENDDO

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=8,10                                         ! Read optional data
            CALL R8FLD ( JCARD(J-6), JF(J-6), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )! Issue warning if fields 5, 6, 7, 8, 9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)


! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MAT1_VALUE_CHECK

      INTEGER(LONG)                   :: MAT1_WERR(5) = (/0, 0, 0, 0, 0/)
                                                            ! Warning error indicators

      REAL(DOUBLE)                    :: FACTOR             ! E/(2*(1+NU)*G)

! **********************************************************************************************************************************
! Check for reasonable input NU value if NU field was not blank

      IF (JCARD_NU(1:) /= ' ') THEN
         IF ((RMATL(NMATL,3) < ZERO) .OR. (RMATL(NMATL,3) > HALF)) THEN
            MAT1_WERR(1) = 1
         ENDIF
      ENDIF

! If E, G, NU fields are all input, make sure 1-|E/(2*(1+NU)*G)| is greater than .01

      IF ((JCARD_E(1:) /= ' ') .AND. (JCARD_G(1:) /= ' ') .AND. (JCARD_NU(1:) /= ' ')) THEN
         FACTOR = RMATL(NMATL,1)/(TWO*(ONE + RMATL(NMATL,3))*RMATL(NMATL,2))
         IF (DABS(ONE - FACTOR) >= .01D0) THEN
            MAT1_WERR(2) = 1
         ENDIF
      ENDIF

! Calc E, G or NU if all were not input. Note E = G = 0 is a fatal error

      IF (JCARD_E(1:) /= ' ') THEN                                           ! E was input

         IF (JCARD_G(1:) /= ' ') THEN                                        !        G was input
            IF (JCARD_NU(1:) /= ' ') THEN                                    !               NU was input
               CONTINUE
            ELSE                                                             !               NU was not input. Calc NU
               RMATL(NMATL,3) = RMATL(NMATL,1)/(TWO*RMATL(NMATL,2)) - ONE
               WRITE(F06,1382) MATL(NMATL,1),RMATL(NMATL,3)
               IF ((RMATL(NMATL,3) < ZERO) .OR. (RMATL(NMATL,3) > HALF)) THEN
                  MAT1_WERR(3) = 1
               ENDIF
            ENDIF
         ELSE                                                                !        G was not input
            IF (JCARD_NU(1:) /= ' ') THEN                                    !               NU was input. Calc G
               RMATL(NMATL,2) = RMATL(NMATL,1)/(TWO*(ONE + RMATL(NMATL,3)))
               WRITE(ERR,1383) MATL(NMATL,1),RMATL(NMATL,2)
               WRITE(F06,1383) MATL(NMATL,1),RMATL(NMATL,2)
            ELSE                                                             !               NU was not input. Set G, NU = 0
               RMATL(NMATL,2) = ZERO
               RMATL(NMATL,3) = ZERO
               WRITE(ERR,1384) MATL(NMATL,1),RMATL(NMATL,2), RMATL(NMATL,3)
               WRITE(F06,1384) MATL(NMATL,1),RMATL(NMATL,2), RMATL(NMATL,3)
            ENDIF
         ENDIF

      ELSE                                                                   ! E was not input

         IF (JCARD_G(1:) /= ' ') THEN                                        !        G was input
            IF (JCARD_NU(1:) /= ' ') THEN                                    !               NU was input. Calc E
               RMATL(NMATL,1) = TWO*RMATL(NMATL,2)*(ONE + RMATL(NMATL,3))
               WRITE(ERR,1385) MATL(NMATL,1),RMATL(NMATL,1)
               WRITE(F06,1385) MATL(NMATL,1),RMATL(NMATL,1)
            ELSE                                                             !               NU was not input. Set E, NU = 0
               RMATL(NMATL,1) = ZERO
               RMATL(NMATL,3) = ZERO
               WRITE(ERR,1386) MATL(NMATL,1),RMATL(NMATL,1),RMATL(NMATL,3)
               WRITE(F06,1386) MATL(NMATL,1),RMATL(NMATL,1),RMATL(NMATL,3)
            ENDIF
         ELSE                                                                !        G was not input (error: E or G must be input)
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1119) MATL(NMATL,1)
            WRITE(F06,1119) MATL(NMATL,1)
         ENDIF

      ENDIF

! Check MAT1_WARN ERROR and write messages if needed

      IF ((MAT1_WERR(1) /=0) .OR. (MAT1_WERR(2) /=0) .OR. (MAT1_WERR(3) /=0)) THEN
         WRITE(ERR,101) CARDP
         IF (ECHO == 'NONE  ') THEN
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) CARDP
            ENDIF
         ENDIF
      ENDIF

      IF (MAT1_WERR(1) /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1391) RMATL(NMATL,3),MATL(NMATL,1)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1391) RMATL(NMATL,3),MATL(NMATL,1)
         ENDIF
      ENDIF

      IF (MAT1_WERR(2) /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1393)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1393)
         ENDIF
      ENDIF

      IF (MAT1_WERR(3) /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1392) RMATL(NMATL,3),MATL(NMATL,1)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1392) RMATL(NMATL,3),MATL(NMATL,1)
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1119 FORMAT(' *ERROR  1119: CANNOT HAVE E AND G ZERO ON MATERIAL ENTRY ',I8)

 1382 FORMAT(' *INFORMATION: MAT1 ENTRY ',I8,' HAD FIELD FOR NU BLANK. MYSTRAN CALCULATED NU = ',1ES13.6)

 1383 FORMAT(' *INFORMATION: MAT1 ENTRY ',I8,' HAD FIELD FOR G  BLANK. MYSTRAN CALCULATED G  = ',1ES13.6)

 1384 FORMAT(' *INFORMATION: MAT1 ENTRY ',I8,' HAD FIELD FOR G AND NU  BLANK. MYSTRAN SET G         = ',1ES13.6,' NU = ',1ES13.6)

 1385 FORMAT(' *INFORMATION: MAT1 ENTRY ',I8,' HAD FIELD FOR E  BLANK. MYSTRAN CALCULATED E  = ',1ES13.6)

 1386 FORMAT(' *INFORMATION: MAT1 ENTRY ',I8,' HAD FIELD FOR E AND NU BLANK. MYSTRAN SET E         = ',1ES13.6,' NU = ',1ES13.6)

 1391 FORMAT(' *WARNING    : UNREASONABLE VALUE OF NU = ',1ES10.3,' INPUT ON MAT1 ENTRY ID = ',I8)

 1392 FORMAT(' *WARNING    : UNREASONABLE VALUE OF NU= ',1ES10.3,' CALCULATED FROM E AND G ON MAT1 ENTRY ID = ',I8)

 1393 FORMAT(' *WARNING    : VALUES INPUT FOR E, G AND NU ARE SUCH THAT |1 - E/(2*(1+NU)*G)| >= .01 (MAT2 IS RECOMMENDED IN THESE',&
                           ' CASES)')

! **********************************************************************************************************************************
      END SUBROUTINE MAT1_VALUE_CHECK

      END SUBROUTINE BD_MAT1


      SUBROUTINE BD_MAT2 ( CARD, LARGE_FLD_INP )

! Processes MAT2 Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMATL, MRMATLC, NMATL, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE MODEL_STUF, ONLY            :  MATL, RMATL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MAT2'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATL_ID   = 0     ! The ID for this MAT2 (field 2)


      REAL(DOUBLE)                    :: R8INP             ! A real input value read



! **********************************************************************************************************************************
! MAT2 Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Material ID     MATL(nmatl, 1)
!           Material type   MATL(nmatl, 2) (e.g. 2 indicates MAT2 entry)
!    3      G11            RMATL(nmatl, 1)
!    4      G12            RMATL(nmatl, 2)
!    5      G13            RMATL(nmatl, 3)
!    6      G22            RMATL(nmatl, 4)
!    7      G23            RMATL(nmatl, 5)
!    8      G33            RMATL(nmatl, 6)
!    9      RHO            RMATL(nmatl, 7)
! on optional second card:
!    2      A1             RMATL(nmatl, 8)
!    3      A2             RMATL(nmatl, 9)
!    4      A3             RMATL(nmatl,10)
!    5      TREF           RMATL(nmatl,11)
!    6      GE             RMATL(nmatl,12)
!    7      ST             RMATL(nmatl,13)
!    8      SC             RMATL(nmatl,14)
!    9      SS             RMATL(nmatl,15)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for too many MATL entries

      NMATL = NMATL+1
      IF (NMATL > LMATL) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMATL
         WRITE(F06,1163) SUBR_NAME,JCARD(1),LMATL
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

      CALL I4FLD ( JCARD(2), JF(2), MATL_ID )              ! Check for duplicate ID number
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NMATL-1
            IF (MATL_ID == MATL(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),MATL_ID
               WRITE(F06,1145) JCARD(1),MATL_ID
               EXIT
            ENDIF
         ENDDO
         MATL(NMATL,1) = MATL_ID
         MATL(NMATL,2) = 2                                 ! Type is 2 for MAT2 card
      ENDIF

      DO J = 1,7                                           ! Read 7 fields of real data on the parent card
         R8INP = ZERO
         CALL R8FLD ( JCARD(J+2), JF(J+2), R8INP )
         IF (IERRFL(J+2) == 'N') THEN
            RMATL(NMATL,J) = R8INP
         ENDIF
      ENDDO

! Null optional data:

      DO J=8,15
         RMATL(NMATL,J) = ZERO
      ENDDO

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=8,15                                         ! Read optional data
            CALL R8FLD ( JCARD(J-6), JF(J-6), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF




      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)


! **********************************************************************************************************************************

      END SUBROUTINE BD_MAT2


      SUBROUTINE BD_MAT8 ( CARD, LARGE_FLD_INP )

! Processes MAT8 Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMATL, MRMATLC, NMATL
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  MATL, RMATL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MAT8'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: SETID             ! The set ID in field 2

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATL_ID   = 0     ! The ID for this MAT8 (field 2)


      REAL(DOUBLE)                    :: E1                 ! Modulus in longitudinal direction
      REAL(DOUBLE)                    :: E2                 ! Modulus in lateral direction
      REAL(DOUBLE)                    :: NU12               ! Poissons ratio



! **********************************************************************************************************************************
! MAT8 Bulk Data Card routine

!   FIELD   ITEM                                                              ARRAY ELEMENT
!   -----   ------------                                                      -------------
!    2      Material ID                                                         MATL(nmatl,1)
!           Material type                                                       MATL(nmatl,2) (e.g. 8 indicates MAT8 entry)
!    3      E1 elastic modulus                                                 RMATL(nmatl,1)
!    4      E2 elastic modulus                                                 RMATL(nmatl,2)
!    5      NU12 Poisson's ratio                                               RMATL(nmatl,3)
!    6      G12 in plane shear modulus                                         RMATL(nmatl,4)
!    7      G1Z Transverse shear modulus (1-Z plane)                           RMATL(nmatl,5)
!    8      G2Z Transverse shear modulus (2-Z plane)                           RMATL(nmatl,6)
!    9      RHO mass density                                                   RMATL(nmatl,7)
! on optional second card:
!    2      A1 thermal expansion coeff (1 direction)                           RMATL(nmatl,8)
!    3      A2 thermal expansion coeff (2 direction)                           RMATL(nmatl,9)
!    4      TREF temperaure reference                                          RMATL(nmatl,10)
!    5      Xt Longitudinal dir tension allow stress or strain                 RMATL(nmatl,11)
!    6      Xc Longitudinal dir compr allow stress or strain                   RMATL(nmatl,12)
!    7      Yt Lateral dir tension allow stress or strain                      RMATL(nmatl,13)
!    8      Yc Lateral dir compr allow stress or strain                        RMATL(nmatl,14)
!    9      S in-plane shear allowable stress or strain                        RMATL(nmatl,15)
! on optional third card:
!    2      GE structural damping coeff                                        RMATL(nmatl,16)
!    3      F12 Interaction term for failure theory of Tsai-Wu                 RMATL(nmatl,17)
!    4      STRN Indicates whether Xt, Xc, Yt, Yc are stress or strain allows  RMATL(nmatl,18)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for too many MATL entries

      NMATL = NMATL+1
      IF (NMATL > LMATL) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMATL
         WRITE(F06,1163) SUBR_NAME,JCARD(1),LMATL
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

      SETID = JCARD(2)
      CALL I4FLD ( JCARD(2), JF(2), MATL_ID )              ! Check for duplicate ID number
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NMATL-1
            IF (MATL_ID == MATL(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),MATL_ID
               WRITE(F06,1145) JCARD(1),MATL_ID
               EXIT
            ENDIF
         ENDDO
         MATL(NMATL,1) = MATL_ID
         MATL(NMATL,2) = 8                                 ! Type is 8 for MAT8 card
      ENDIF

      DO J = 1,7
         CALL R8FLD ( JCARD(J+2), JF(J+2), RMATL(NMATL,J) )
      ENDDO

      E1 = ZERO
      IF (IERRFL(3) == 'N') THEN
         E1 = RMATL(NMATL,1)
      ENDIF

      E2 = ZERO
      IF (IERRFL(4) == 'N') THEN
         E2 = RMATL(NMATL,2)
      ENDIF

      NU12 = ZERO
      IF (IERRFL(5) == 'N') THEN
         NU12 = RMATL(NMATL,3)
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Null optional data:

      DO J = 8,18
         RMATL(NMATL,J) = ZERO
      ENDDO

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=8,15
            CALL R8FLD ( JCARD(J-6), JF(J-6), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF

! Optional third card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=16,18
            CALL R8FLD ( JCARD(J-14), JF(J-14), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )! Issue warning if fields 5, 6, 7, 8, 9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF

! Check MAT8 values

      CALL MAT8_VALUE_CHECK



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MAT8_VALUE_CHECK

      USE CONSTANTS_1, ONLY           :  ONE

      CHARACTER( 2*BYTE)              :: MODULUS            ! Character to print out for error

      INTEGER(LONG)                   :: JERR    = 0        ! Local erro count

      REAL(DOUBLE)                    :: EPS1               ! A small number to compare real zero
      REAL(DOUBLE)                    :: DEN                ! (1 - NU12*NU21)
      REAL(DOUBLE)                    :: NU21               ! Poissons ratio = N12*(E2/E1)

      INTRINSIC                       :: DABS

! **********************************************************************************************************************************
      EPS1 = DABS(EPSIL(1))

! Make sure E1 and E2 > 0.

      IF (DABS(E1) < EPS1) THEN
         JERR = JERR + 1
         MODULUS = 'E1'
         WRITE(ERR,1115) SETID,MODULUS,E1,EPS1
         WRITE(F06,1115) SETID,MODULUS,E1,EPS1
      ENDIF

      IF (DABS(E2) < EPS1) THEN
         JERR = JERR + 1
         MODULUS = 'E2'
         WRITE(ERR,1115) SETID,MODULUS,E2,EPS1
         WRITE(F06,1115) SETID,MODULUS,E2,EPS1
      ENDIF

! Make sure (1 - NU12*NU21) > 0.

      IF (DABS(E2) > ZERO) THEN
         NU21 = NU12*(E1/E2)
         DEN = ONE - NU12*NU21
         IF (DABS(DEN) < EPS1) THEN
            JERR = JERR + 1
            WRITE(ERR,1116) SETID,DEN,EPS1
            WRITE(F06,1116) SETID,DEN,EPS1
         ENDIF
      ENDIF

! If JERR > 0 set fatal error flag

      IF (JERR > 0) THEN
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! **********************************************************************************************************************************
 1115 FORMAT(' *ERROR  1115: MAT8 MATERIAL BULK DATA ENTRY ',A,' HAS MODULUS ',A,' EQUAL TO ',1ES11.2,'. MUST BE > ',1ES11.2)

 1116 FORMAT(' *ERROR  1116: MAT8 MATERIAL BULK DATA ENTRY ',A,' HAS POISSON RATIO NU12 SUCH THAT 1 - NU12*NU21 = ',1ES11.2,       &
                          '. MUST BE > ',1ES11.2)

! **********************************************************************************************************************************
      END SUBROUTINE MAT8_VALUE_CHECK

      END SUBROUTINE BD_MAT8


      SUBROUTINE BD_MAT9 ( CARD, LARGE_FLD_INP )

! Processes MAT9 Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMATL, MRMATLC, NMATL, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE MODEL_STUF, ONLY            :  MATL, RMATL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MAT9'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN(JCARD))           :: NAME              ! JCARD(1) from parent entry

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATL_ID   = 0     ! The ID for this MAT9 (field 2)


      REAL(DOUBLE)                    :: R8INP             ! A real input value read



! **********************************************************************************************************************************
! MAT9 Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Material ID     MATL(nmatl, 1)
!           Material type   MATL(nmatl, 2) (e.g. 2 indicates MAT9 entry)
!    3      G11            RMATL(nmatl, 1)
!    4      G12            RMATL(nmatl, 2)
!    5      G13            RMATL(nmatl, 3)
!    6      G14            RMATL(nmatl, 4)
!    7      G15            RMATL(nmatl, 5)
!    8      G16            RMATL(nmatl, 6)
!    9      G22            RMATL(nmatl, 7)
! on first continuation card:
!    2      G23            RMATL(nmatl, 8)
!    3      G24            RMATL(nmatl, 9)
!    4      G25            RMATL(nmatl,10)
!    5      G26            RMATL(nmatl,11)
!    6      G33            RMATL(nmatl,12)
!    7      G34            RMATL(nmatl,13)
!    8      G35            RMATL(nmatl,14)
!    9      G36            RMATL(nmatl,15)
! on second continuation card:
!    2      G44            RMATL(nmatl,16)
!    3      G45            RMATL(nmatl,17)
!    4      G46            RMATL(nmatl,18)
!    5      G55            RMATL(nmatl,19)
!    6      G56            RMATL(nmatl,20)
!    7      G66            RMATL(nmatl,21)
!    8      RHO            RMATL(nmatl,22)
!    9      A1             RMATL(nmatl,23)
! on optional third continuation card:
!    2      A2             RMATL(nmatl,24)
!    3      A3             RMATL(nmatl,25)
!    4      A4             RMATL(nmatl,26)
!    5      A6             RMATL(nmatl,27)
!    6      A6             RMATL(nmatl,28)
!    7      TREF           RMATL(nmatl,29)
!    8      GE             RMATL(nmatl,30)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NAME = JCARD(1)
      ID   = JCARD(2)

! Check for too many MATL entries

      NMATL = NMATL+1
      IF (NMATL > LMATL) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMATL
         WRITE(F06,1163) SUBR_NAME,JCARD(1),LMATL
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

      CALL I4FLD ( JCARD(2), JF(2), MATL_ID )              ! Check for duplicate ID number
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NMATL-1
            IF (MATL_ID == MATL(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),MATL_ID
               WRITE(F06,1145) JCARD(1),MATL_ID
               EXIT
            ENDIF
         ENDDO
         MATL(NMATL,1) = MATL_ID
         MATL(NMATL,2) = 9                                 ! Type is 9 for MAT9 card
      ENDIF

      DO J = 1,7                                           ! Read 7 fields of real data on the parent card
         R8INP = ZERO
         CALL R8FLD ( JCARD(J+2), JF(J+2), R8INP )
         IF (IERRFL(J+2) == 'N') THEN
            RMATL(NMATL,J) = R8INP
         ENDIF
      ENDDO

! Mandatory 1st continuation card card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=8,15
            CALL R8FLD ( JCARD(J-6), JF(J-6), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF

! Mandatory 2nd continuation card card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=16,23
            CALL R8FLD ( JCARD(J-14), JF(J-14), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF

! Null optional data:

      DO J=24,MRMATLC
         RMATL(NMATL,J) = ZERO
      ENDDO

! Optional third continuation card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         DO J=24,30
            CALL R8FLD ( JCARD(J-22), JF(J-22), RMATL(NMATL,J) )
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )! Issue warning if field 9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF




      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)


! **********************************************************************************************************************************

      END SUBROUTINE BD_MAT9

   END MODULE MATERIAL_CARDS
