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


      SUBROUTINE ELMDAT2 ( INT_ELEM_ID, OPT, WRITE_WARN )

! Generates small arrays of elem data, for use by subroutine EMG, one elem at a time for all elems. Arrays generated are:

!   DT (1 elem temperatures) and PRESS (1 element pressure load)
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LPDAT, MPRESS, MDT, MTDAT_TEMPRB, NSUB, NTSUB 
      USE TIMDAT, ONLY                :  TSEC
! --- cbeam_add begin --- !
      USE SCONTR, ONLY                :  MPDAT_PLOAD1
! --- cbeam_add end --- !
      USE CONSTANTS_1, ONLY           :  ZERO, QUARTER, THIRD
      USE MODEL_STUF, ONLY            :  BGRID, DT, ELGP, ETYPE, GTEMP, PDATA, PPNT, PTYPE, PRESS, TDATA, TPNT, TYPE

      USE ELMDAT2_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMDAT2'
      CHARACTER(LEN=*), INTENT(IN)    :: OPT(6)            ! Array of EMG option indicators
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID for which
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IPPN              ! A pointer into array PPNT
      INTEGER(LONG)                   :: ITPN              ! A pointer into array TPNT


! **********************************************************************************************************************************
      TYPE = ETYPE(INT_ELEM_ID)

! ELGP is the number of G.P.'s for this elem. Call GET_ELGP to find out how many grids there are for elem type TYPE

      CALL GET_ELGP ( INT_ELEM_ID )

! **********************************************************************************************************************************
! Generate DT temp array for this elem if OPT(2) or OPT(3) says to. Need to calc DT for OPT(3) since SEi stress arrays use DT
! The temp data may be on G.P. temp cards or elem. cards.
! The G.P. temperature data is stored in GTEMP and the element temperature data in TDATA (with TPNT pointer array). These arrays
! have been generated in an earlier subr (TEMPERATURE_DATA_PROC)

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN


         DT= ZERO


         IF (NTSUB > 0) THEN

            IF ((TYPE(1:4) /= 'ELAS') .OR. (TYPE /= 'USERIN  ') .OR. (TYPE /= 'BUSH    ')) THEN
                                                           ! Num of entries on element temp. card is 6 if BAR/BEAM and 2 otherwise
               DO J=1,NTSUB                                ! Calculate DT in 1 of 2 ways, depending on where element temps located
                  ITPN = TPNT(INT_ELEM_ID,J)
                  IF (ITPN /= 0) THEN                      ! Temperatures are on elem temp cards

                     IF ((TYPE == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'USER1   ')) THEN
                        DO I=1,MTDAT_TEMPRB
                           DT(I,J) = TDATA(ITPN-1+I)
                        ENDDO 
                     ELSE 
                        DO I=1,ELGP
                           DT(I,J) = TDATA(ITPN)
                        ENDDO 
                        DT(ELGP+1,J) = TDATA(ITPN+1)
                     ENDIF

                  ELSE                                     ! Temperatures are on grid cards

                     DO I=1,ELGP     
                        DT(I,J) = GTEMP(BGRID(I),J)  
                     ENDDO 

                  ENDIF

               ENDDO 

            ENDIF

         ENDIF

      ENDIF
 
! **********************************************************************************************************************************
! Generate PRESS pressure array for this element if OPT(5) says to. The pressure data was on elem pressure B.D. cards and
! has been transferred to arrays PPNT, PDATA by an earlier subr (PRESSURE_DATA_PROC) 

      IF (OPT(5) == 'Y') THEN

         IF (LPDAT > 0) THEN
!---  cbeam add --- begin!
            DO I=1,MPRESS
               DO J=1,NSUB
                  PRESS(I,J) = ZERO
               ENDDO
            ENDDO 

            IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BART    ') .OR. (TYPE == 'BEAM    ')) THEN

               IF (PTYPE(INT_ELEM_ID) == '2') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
!*** ADDED -- 2026-04-20 -- FOR BEAM DSB ***
                        DO J=1,MPDAT_PLOAD1
                           PRESS(J,I) = PDATA(IPPN+J-1)
                        ENDDO
                     ENDIF
                  ENDDO

               ENDIF
!---  cbeam add --- end!
            ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN

               IF      (PTYPE(INT_ELEM_ID) == '1') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
                        PRESS(3,I) = PDATA(IPPN)
                     ENDIF
                  ENDDO

               ELSE IF (PTYPE(INT_ELEM_ID) == '3') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
                        PRESS(3,I) = THIRD*( PDATA(IPPN) + PDATA(IPPN+1) + PDATA(IPPN+2) )
                     ENDIF
                  ENDDO

               ELSE IF (PTYPE(INT_ELEM_ID) == '4') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
                        PRESS(3,I) = QUARTER*( PDATA(IPPN) + PDATA(IPPN+1) + PDATA(IPPN+2) + PDATA(IPPN+3) )
                     ENDIF
                  ENDDO

               ENDIF

            ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

               IF (PTYPE(INT_ELEM_ID) == '3') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
                        PRESS(3,I) = THIRD*( PDATA(IPPN) + PDATA(IPPN+1) + PDATA(IPPN+2) )
                     ENDIF
                  ENDDO

               ELSE IF (PTYPE(INT_ELEM_ID) == '4') THEN

                  DO I=1,NSUB
                     IPPN = PPNT(INT_ELEM_ID,I)
                     IF (IPPN /= 0) THEN
                        PRESS(3,I) = QUARTER*( PDATA(IPPN) + PDATA(IPPN+1) + PDATA(IPPN+2) + PDATA(IPPN+3) )
                     ENDIF
                  ENDDO

               ENDIF

            ENDIF

         ENDIF

      ENDIF 



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ELMDAT2


