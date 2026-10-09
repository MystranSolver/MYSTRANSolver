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

   MODULE QUADRILATERAL_SHAPE_FUNCTIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SHP2DQ

   CONTAINS

      SUBROUTINE SHP2DQ ( IGAUS, JGAUS, NUM_NODES, CALLING_SUBR, IORD_MSG, IORZZZ, SSI, SSJ, WRT_BUG_THIS_TIME, PSH, DPSHG )

! Generates shape functions for 2D elements.

! The node numbering and axes convention are shown below with XI and ETA ranging from -1 to +1

!                         ETA
!                          |
!                          |
!                          |
!              4 . . . . . 7 . . . . . 3
!              .           |           .
!              .           |           .
!              .           |           .
!              8           ------------6----> XI
!              .                       .
!              .                       .
!              .                       .
!              1 . . . . . 5 . . . . . 2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, ERR, F06, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_SHPJ_BIT, MEFE, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE MODEL_STUF, ONLY            :  EID, EMG_IFE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SHP2DQ'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr (used for debug output)
      CHARACTER(LEN=*), INTENT(IN)    :: IORD_MSG          ! Character name of the integration order (used for debug output)
      CHARACTER(17*BYTE)              :: NAME(5)           ! Used for output annotation
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gauss point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gauss point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: IORZZZ            ! Integration order (used for debug output)
      INTEGER(LONG), INTENT(IN)       :: NUM_NODES         ! Number of element nodes
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: NODES_4     = 4   ! Number of nodes for one type of element
      INTEGER(LONG)                   :: NODES_8     = 8   ! Number of nodes for one type of element


      REAL(DOUBLE) , INTENT(IN)       :: SSI               ! Gauss point location component
      REAL(DOUBLE) , INTENT(IN)       :: SSJ               ! Gauss point location component
      REAL(DOUBLE) , INTENT(OUT)      :: PSH(NUM_NODES)    ! Shape functions for all grid points for this Gauss point
      REAL(DOUBLE) , INTENT(OUT)      :: DPSHG(2,NUM_NODES)! Derivatives of PSH with respect to xi and eta.
      REAL(DOUBLE)                    :: XI(16), ET(16)    ! Elem node locations in isoparametric xi, eta coords
      REAL(DOUBLE)                    :: A1,A2,A3          ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: B1,B2,B3          ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: XI2               ! Squares of xi coords
      REAL(DOUBLE)                    :: ET2               ! Squares of eta coords


! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NUM_NODES
         PSH(I) = ZERO
      ENDDO

      DO I=1,2
         DO J=1,NUM_NODES
            DPSHG(I,J) = ZERO
         ENDDO
      ENDDO

! Generate shape functions for 2D element with 4 grid points

      IF (NUM_NODES == NODES_4) THEN

         XI(1) = -ONE
         XI(2) =  ONE
         XI(3) =  ONE
         XI(4) = -ONE

         ET(1) = -ONE
         ET(2) = -ONE
         ET(3) =  ONE
         ET(4) =  ONE

         DO I=1,4
            PSH(I) = (ONE + SSI*XI(I))*(ONE + SSJ*ET(I))/FOUR
            DPSHG(1,I) = XI(I)*(ONE + SSJ*ET(I))/FOUR
            DPSHG(2,I) = ET(I)*(ONE + SSI*XI(I))/FOUR
         ENDDO

! Generate shape functions and derivatives for 8 node serendipity elem

      ELSE IF (NUM_NODES == NODES_8) THEN

         XI(1) = -ONE
         XI(2) =  ONE
         XI(3) =  ONE
         XI(4) = -ONE
         XI(5) =  ZERO
         XI(6) =  ONE
         XI(7) =  ZERO
         XI(8) = -ONE

         ET(1) = -ONE
         ET(2) = -ONE
         ET(3) =  ONE
         ET(4) =  ONE
         ET(5) = -ONE
         ET(6) =  ZERO
         ET(7) =  ONE
         ET(8) =  ZERO

         DO I=1,8
            XI2 = XI(I)*XI(I)
            ET2 = ET(I)*ET(I)
            A1 = ONE + SSI*XI(I)
            B1 = ONE + SSJ*ET(I)
            A2 = ONE - SSI*SSI
            B2 = ONE - SSJ*SSJ
            A3 = ONE - XI2
            B3 = ONE - ET2
            PSH(I)     = ((A1 - A2)*B1 - B2*A1)*XI2*ET2/FOUR + A2*B1*A3*ET2/TWO + B2*A1*B3*XI2/TWO

            DPSHG(1,I) = ((TWO*SSI + XI(I))*B1 - XI(I)*B2)*XI2*ET2/FOUR - SSI*B1*A3*ET2 + XI(I)*B2*B3*XI2/TWO

            DPSHG(2,I) = ((TWO*SSJ + ET(I))*A1 - ET(I)*A2)*XI2*ET2/FOUR - SSJ*A1*B3*XI2 + ET(I)*A2*A3*ET2/TWO
         ENDDO

! Error: NUM_NODES is not 4 or 8

      ELSE

         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(F06,1932) SUBR_NAME, NUM_NODES, TYPE, EID, NODES_4, NODES_8
            WRITE(ERR,1932) SUBR_NAME, NUM_NODES, TYPE, EID, NODES_4, NODES_8
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1932
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = NUM_NODES
            ENDIF
         ENDIF
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ENDIF

! **********************************************************************************************************************************
! ELDATA output:

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(7) > 0)) THEN

         NAME(1) = 'Nodes  1 thru  4:'
         NAME(2) = 'Nodes  5 thru  8:'
         NAME(3) = 'Nodes  9 thru 12:'
         NAME(4) = 'Nodes 13 thru 16:'
         NAME(5) = 'Nodes 17 thru 20:'

         WRITE(BUG,1101) ELDT_BUG_SHPJ_BIT, TYPE, EID, CALLING_SUBR

         IF      (IORZZZ > 0) THEN
            WRITE(BUG,1102) IGAUS, JGAUS, IORD_MSG,  IORZZZ
         ELSE IF (IORZZZ < 0) THEN
            WRITE(BUG,1103) IGAUS, JGAUS, IORD_MSG, -IORZZZ
         ENDIF
         WRITE(BUG,1104) SSI,SSJ
         WRITE(BUG,1112)

         J = 0
         WRITE(BUG,1105)
         DO I=1,NUM_NODES,4
            J = J + 1
            WRITE(BUG,1106) NAME(J), PSH(I), PSH(I+1), PSH(I+2), PSH(I+3)
         ENDDO
         WRITE(BUG,*)
         WRITE(BUG,*)

         J = 0
         WRITE(BUG,1115)
         DO I=1,NUM_NODES,4
            J = J + 1
            WRITE(BUG,1106) NAME(J), DPSHG(1,I), DPSHG(1,I+1), DPSHG(1,I+2), DPSHG(1,I+3)
         ENDDO
         WRITE(BUG,*)
         WRITE(BUG,*)

         J = 0
         WRITE(BUG,1125)
         DO I=1,NUM_NODES,4
            J = J + 1
            WRITE(BUG,1106) NAME(J), DPSHG(2,I), DPSHG(2,I+1), DPSHG(2,I+2), DPSHG(2,I+3)
         ENDDO
         WRITE(BUG,*)
         WRITE(BUG,*)

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1932 FORMAT(' *ERROR  1932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' SUBROUTINE CALLED WITH ILLEGAL NUMBER OF NODES =',I8,' FOR ',A,' ELEMENT',I8,'. MUST BE ',I3,' OR ',I3)

 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,', Calling subroutine ',A,/,                            &
             ' ========================================================================================',/)

 1102 FORMAT(51X,' Gauss point: I = ',I3,', J = ',I3,/,31X,'with integration order ',A,I2,/,                                       &
             35X,'and using Gaussian integration to determine strain-displ matrix'//)

 1103 FORMAT(51X,' Gauss point: I = ',I3,', J = ',I3,/,31X,'with integration order ',A,I2,/,                                       &
             33X,'and using Jacobian weighted average to determine strain-displ matrix'//)

 1112 format(51x,'Outputs from subroutine SHP2DQ',/,51x,'------------------------------',/)

 1104 FORMAT(53X,'Coordinates of Gauss point',/,57x,'XI              ET',/,48X,2(1ES16.6),//)

 1105 FORMAT(58X,'Shape functions')

 1106 FORMAT(22X,A,1X,4(1ES16.6))

 1115 FORMAT(42X,'Derivatives of shape functions with respect to XI')

 1125 FORMAT(42X,'Derivatives of shape functions with respect to ET')

! **********************************************************************************************************************************

      END SUBROUTINE SHP2DQ

   END MODULE QUADRILATERAL_SHAPE_FUNCTIONS
