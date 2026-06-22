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
 
      SUBROUTINE CALC_ELEM_NODE_FORCES
 
! Calculates elem nodal forces in local elem coord system for one elem and one subcase for all element types.
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, INT_SC_NUM, JTSUB, NCORD, NGRID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, CORD, EID, ELAS_COMP, ELDOF, ELGP, GRID, KE, KEG, KEO_BUSH,                &
                                         PEB, PE_GA_GB, PEG, PEL, PTE, RCORD, SCNUM, SUBLOD, TYPE, UEB, UEG, UEL, TE, TE_GA_GB
 
      USE CALC_ELEM_NODE_FORCES_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_ELEM_NODE_FORCES'
      CHARACTER( 1*BYTE)              :: TEMP_OPT(6)       ! Array of EMG option indicators
 
      INTEGER(LONG)                   :: ACID_G            ! Actual coordinate system ID
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: ICID              ! Internal coord sys no. corresponding to an actual coord sys no. 
      INTEGER(LONG)                   :: I1,I2             ! Calculated displ component no's for ELAS elems
      INTEGER(LONG)                   :: NCOLS             ! Number of rows in element stiffness matrix
      INTEGER(LONG)                   :: NROWS             ! Number of cols in element stiffness matrix
      INTEGER(LONG)                   :: NUM_COMPS_GRID_1  ! No. displ components for 1st grid on ELAS elems
 
      REAL(DOUBLE)                    :: DUM1(3),DUM2(3)   ! Intermediate variables
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: T0G(3,3)          ! Matrix to transform offsets from global to basic  coords 
      REAL(DOUBLE)                    :: TET(3,3)          ! Transpose of TE
      REAL(DOUBLE)                    :: TET_GA_GB(3,3)    ! Transpose of TE
      REAL(DOUBLE)                    :: TR(12,12)         ! Matrix with 4 TE matrices on the diagonal


! **********************************************************************************************************************************
      NROWS = ELDOF
      NCOLS = ELDOF
 
      DO I=1,NCOLS
         PEL(I) = ZERO
      ENDDO 
 
! **********************************************************************************************************************************
! Calc forces for one element. The ELAS and ROD1 elem have very sparse stiffness matrices, so an explicit form is used
! for them. All other element forces are calculated by multiplication of complete stiffness matrix with the displ's.
 
      IF (TYPE(1:4) == 'ELAS') THEN                        ! Calculate forces for ELAS1-4 elems
 
         I1 = ELAS_COMP(1)
         CALL GET_GRID_NUM_COMPS ( AGRID(1), NUM_COMPS_GRID_1, SUBR_NAME ) !AGRID or BGRID ?
         I2 = NUM_COMPS_GRID_1 + ELAS_COMP(2)
         PEL(I1) = KE(I1,I1)*UEL(I1) + KE(I1,I2)*UEL(I2)   ! Note: KE is global and local for the ELAS elems
         PEL(I2) = KE(I2,I1)*UEL(I1) + KE(I2,I2)*UEL(I2)
 
      ELSE IF (TYPE == 'ROD     ') THEN                    ! Calculate forces for ROD1 elem
 
         IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
            PEL(1) = -PTE(1,JTSUB)
            PEL(7) = -PTE(7,JTSUB)
         ENDIF

         PEL( 1) = PEL(1) + KE( 1, 1)*UEL( 1) + KE( 1, 7)*UEL( 7)
         PEL( 4) =          KE( 4, 4)*UEL( 4) + KE( 4,10)*UEL(10)
         PEL( 7) = PEL(7) + KE( 7, 1)*UEL( 1) + KE( 7, 7)*UEL( 7)
         PEL(10) =          KE(10, 4)*UEL( 4) + KE(10,10)*UEL(10)
 
      ELSE IF (TYPE == 'USERIN  ') THEN

         WRITE(F06,9991) TYPE

      ELSE IF (TYPE == 'BUSH    ') THEN

         CALL ELEM_TRANSFORM_LBG ( 'KE', KE, PTE )
         
         DO I=1,6
            TEMP_OPT(I) ='N'
         ENDDO
         TEMP_OPT(4) = 'Y'
         CALL ELMOFF ( TEMP_OPT, 'Y' )         

         DO I=1,NROWS
            DO J=1,NCOLS
               PEG(I) = PEG(I) + KEG(I,J)*UEG(J)
            ENDDO
         ENDDO

         DO I=1,2
            ACID_G = GRID(BGRID(I),3)                      ! Get global coord sys for this grid
            IF (ACID_G /= 0) THEN                          ! Global is not basic so need to transform offset from basic to global
               ICID = 0
               DO J=1,NCORD
                  IF (ACID_G == CORD(J,2)) THEN
                     ICID = J
                     EXIT
                  ENDIF
               ENDDO   
               CALL GEN_T0L ( BGRID(I), ICID, THETAD, PHID, T0G )
               IF (I == 1) THEN
                  DO J=1,3
                     PEB(J)   = T0G(J,1)*PEG(J)   + T0G(J,2)*PEG(J)   + T0G(J,3)*PEG(J) 
                  ENDDO   
                  DO J=1,3
                     PEB(J+3) = T0G(J,1)*PEG(J+3) + T0G(J,2)*PEG(J+3) + T0G(J,3)*PEG(J+3) 
                  ENDDO
               ELSE   
                  DO J=1,3
                     PEB(J+6) = T0G(J,1)*PEG(J+6) + T0G(J,2)*PEG(J+6) + T0G(J,3)*PEG(J+6) 
                  ENDDO   
                  DO J=1,3
                     PEB(J+9) = T0G(J,1)*PEG(J+9) + T0G(J,2)*PEG(J+9) + T0G(J,3)*PEG(J+9) 
                  ENDDO
               ENDIF   
            ELSE                                           ! Global was basic so no transformation of coords needed
               DO J=1,12
                  PEB(J) = PEG(J)
               ENDDO
            ENDIF
         ENDDO

         DO I=1,3
            DO J=1,3
               TET(I,J) = TE(J,I)
            ENDDO
         ENDDO

         DO I=1,12
            DO J=1,12
               TR(I,J) = ZERO
            ENDDO
         ENDDO

         DO I=1,3                                          ! TR is a 12x12 matrix with 4 TET matrices on its diagonal
            DO J=1,3
               TR(I  ,J  ) = TET(I,J)
               TR(I+3,J+3) = TET(I,J)
               TR(I+6,J+6) = TET(I,J)
               TR(I+9,J+9) = TET(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF_T ( TR, PEB, 12, 12, 1, PEL )
         DO I=1,3
            DO J=1,3
               TET_GA_GB(I,J) = TE_GA_GB(J,I)
            ENDDO
         ENDDO

         DO I=1,12
            DO J=1,12
               TR(I,J) = ZERO
            ENDDO
         ENDDO

         DO I=1,3                                          ! TR is a 12x12 matrix with 4 TET matrices on its diagonal
            DO J=1,3
               TR(I  ,J  ) = TET_GA_GB(I,J)
               TR(I+3,J+3) = TET_GA_GB(I,J)
               TR(I+6,J+6) = TET_GA_GB(I,J)
               TR(I+9,J+9) = TET_GA_GB(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF_T ( TR, PEB, 12, 12, 1, PE_GA_GB )

      ELSE                                                 ! Calculate forces for any other type of elem

         DO I=1,NROWS
 
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
               PEL(I) = -PTE(I,JTSUB)
            ENDIF
 
            DO J=1,NCOLS
               PEL(I) = PEL(I) + KE(I,J)*UEL(J)
            ENDDO
         ENDDO
 
      ENDIF


      IF (DEBUG(56) > 0) THEN                              ! Print UEL, PEL to f06 for debug    
         WRITE(F06,5000)
         WRITE(F06,5004) TRIM(TYPE), EID, SCNUM(INT_SC_NUM)
         J = 0 ; K = 0
         DO I=1,ELDOF
            J = J+1
            IF (J == 1) THEN
               K = K+1
               WRITE(F06,5006) '  GRID, COMP, UEL, PEL = ', AGRID(K), J, UEL(I), PEL(I)
            ELSE
               WRITE(F06,5007) '        COMP, UEL, PEL = ', J, UEL(I), PEL(I)
            ENDIF
            IF (J == 6) THEN
               WRITE(F06,*)
               J = 0
            ENDIF
         ENDDO
         WRITE(F06,5000)
      ENDIF
   

      RETURN

! **********************************************************************************************************************************

 5000 FORMAT('=============================================================================================================',/)

 5004 FORMAT(16X, 'S U B R O U T I N E   CALC_ELEM_NODE_FORCES   F O R   ',A,I8,', S/C',I8,/,                                      &
             
             28X,'(displacements and node forces in element coordinates)',/)

 5006 FORMAT(21X,A,I8,1X,I3,2ES15.6)

 5007 FORMAT(21X,A,9X,I3,2ES15.6)

 9991 FORMAT(' *INFORMATION: ELEMENT NODE FORCE CALCULATION NOT PROGRAMMED FOR ',A,' ELEMENTS')

! **********************************************************************************************************************************

      END SUBROUTINE CALC_ELEM_NODE_FORCES

