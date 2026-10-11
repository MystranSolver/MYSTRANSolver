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

   MODULE ELEMENT_RECOVERY_SUPPORT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELEM_STRE_STRN_ARRAYS, CALC_ELEM_NODE_FORCES, GET_COMP_SHELL_ALLOWS, POLYNOM_FIT_STRE_STRN, SHELL_ENGR_FORCE_OGEL

   CONTAINS

      SUBROUTINE ELEM_STRE_STRN_ARRAYS ( STR_PT_NUM )

! Calculates element stress and strain arrays (arrays STRESS and STRAIN). Stresses are calculated for all engineering elements
! (1-D, 2-D, 3-D elements). Strains are calculated for the BUSH element, all 2-D and 3-D elements. The default method for
! calculating stresses (for 2D and 3D elements) is to first calculate strains by multiplying the element strain-displ matrices
! (BEi) times the elem displ's (in local elem coords) and then calc stresses from those strains using material props.

! For the TRIA3 and QUAD4 there is a DEBUG option to calculate stresses directly by multiplying the stress-displ matrices (SEi)
! times the elem displ's

! The arrays STRESS and STRAIN are calculated at the mid plane of 2-D elements and have to be processed later to get element
! specific outputs (e.g. stresses at the top and bottom of the 2-D element). The same is true for some of the 1-D elements (e.g.
! the BAR element stresses at the 4 points on the cross-section have to be processed from the STRESS array generated here)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, INT_SC_NUM, JTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, one, four
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, BE3, DT, EM, EB, ES, ET, ELDOF, PEL, PHI_SQ, STRAIN, STRESS, SUBLOD,    &
                                         TREF, TYPE, UEL, UEB, SE1, SE2, SE3, STE1, STE2, STE3, ELGP, ISOLID
      USE DEBUG_PARAMETERS
      USE PARAMS, ONLY                :  STR_CID, QUAD4TYP
      USE MITC_STUF, ONLY             :  THERM_KAPPA, MSTRPT


      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF
      USE RESULT_COORDINATES, ONLY    :  STR_TENSOR_TRANSFORM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELEM_STRE_STRN_ARRAYS'

      INTEGER(LONG), INTENT(IN)       :: STR_PT_NUM        ! Which point (3rd index in SEi matrices) this call is for
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter

      INTEGER(LONG)                   :: STR_CID_SOLID

      REAL(DOUBLE)                    :: ALPT(6)           ! Col of ALPVEC times temperatures
      REAL(DOUBLE)                    :: ALPTM(3)          ! Col of ALPVEC times temperatures
      REAL(DOUBLE)                    :: ALPTB(3)          ! Col of ALPVEC times temperatures
      REAL(DOUBLE)                    :: ALPTT(3)          ! Col oF ALPVEC times temperatures
      REAL(DOUBLE)                    :: DUM31(3)          ! Array used in an intermediate calc
      REAL(DOUBLE)                    :: DUM32(3)          ! Array used in an intermediate calc
      REAL(DOUBLE)                    :: DUM33(3)          ! Array used in an intermediate calc
      REAL(DOUBLE)                    :: ET3(3,3)          ! Material matrix ET expanded TO 3x3
      REAL(DOUBLE)                    :: STRAIN1(3)        ! 1st 3 rows of array STRAIN
      REAL(DOUBLE)                    :: STRAIN2(3)        ! 2nd 3 rows of array STRAIN
      REAL(DOUBLE)                    :: STRAIN3(3)        ! 3rd 3 rows of array STRAIN
      REAL(DOUBLE)                    :: STRESS1(3)        ! 1st 3 rows of array STRESS
      REAL(DOUBLE)                    :: STRESS2(3)        ! 2nd 3 rows of array STRESS
      REAL(DOUBLE)                    :: STRESS3(3)        ! 3rd 3 rows of array STRESS
      REAL(DOUBLE)                    :: STRESS_THERM(6)   ! Part of array STRESS
      REAL(DOUBLE)                    :: STRESS1_THERM(3)  ! Part of array STRESS1
      REAL(DOUBLE)                    :: STRESS2_THERM(3)  ! Part of array STRESS2
      REAL(DOUBLE)                    :: STRESS3_THERM(3)  ! Part of array STRESS3
      REAL(DOUBLE)                    :: STRESS_MECH(6)    ! Part of array STRESS
      REAL(DOUBLE)                    :: STRESS1_MECH(3)   ! Part of array STRESS1
      REAL(DOUBLE)                    :: STRESS2_MECH(3)   ! Part of array STRESS2
      REAL(DOUBLE)                    :: STRESS3_MECH(3)   ! Part of array STRESS3
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature
      REAL(DOUBLE)                    :: STR_TENSOR(3,3)   ! 2D stress or strain tensor



! **********************************************************************************************************************************
! Initialize

      DO I=1,9
         STRAIN(I) = ZERO
         STRESS(I) = ZERO
      ENDDO

! **********************************************************************************************************************************
! Calc stresses for 1D elements
! Warning: For ELAS1/2/3/4, the STRESS() array contains force, not stress. This is because we must calculate stress from force, not
! the other way around.

      IF ((TYPE(1:3) == 'BAR') .OR. (TYPE(1:4) == 'BUSH') .OR. (TYPE(1:4) == 'ELAS') .OR. (TYPE(1:3) == 'ROD') .OR.                &
          (TYPE(1:5) == 'USER1')) THEN

         DO I=1,3
            STRESS(I) = ZERO
            DO J=1,ELDOF
               STRESS(I) = STRESS(I) + SE1(I,J,STR_PT_NUM)*UEL(J)
               if (dabs(uel(j)) > 1.e-15) then
               endif
            ENDDO
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
               STRESS(I) = STRESS(I) - STE1(I,JTSUB,STR_PT_NUM)
            ENDIF
         ENDDO

         IF ((TYPE(1:3) == 'BAR') .OR. (TYPE(1:4) == 'BUSH')) THEN
            K = 0
            DO I=4,6
               STRESS(I) = ZERO
               K = K + 1
               IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
                  STRESS(I) = -STE2(K,JTSUB,STR_PT_NUM)
               ENDIF
               DO J=1,ELDOF
                  STRESS(I) = STRESS(I) + SE2(K,J,STR_PT_NUM)*UEL(J)
               ENDDO
            ENDDO
         ENDIF


      ENDIF

! **********************************************************************************************************************************
! Calc strains for 1D BUSH element

      IF (TYPE(1:4) == 'BUSH') THEN

         DO I=1,3
            STRAIN(I) = ZERO
            DO J=1,ELDOF
               STRAIN(I) = STRAIN(I) + BE1(I,J,STR_PT_NUM)*UEL(J)
            ENDDO
         ENDDO

         K = 0
         DO I=4,6
            STRAIN(I) = ZERO
            K = K + 1
            DO J=1,ELDOF
               STRAIN(I) = STRAIN(I) + BE2(K,J,STR_PT_NUM)*UEL(J)
            ENDDO
         ENDDO


! **********************************************************************************************************************************
! Calc strains, then stresses for 2D elements

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8') .OR.                                 &
               (TYPE(1:5) == 'SHEAR') .OR. (TYPE(1:5) == 'USER1')) THEN

         DO I=1,3
            STRAIN(I) = ZERO
            STRAIN(I+3) = ZERO
            DO J=1,ELDOF
               STRAIN(I)   = STRAIN(I)   + BE1(I,J,STR_PT_NUM)*UEL(J)
               STRAIN(I+3) = STRAIN(I+3) + BE2(I,J,STR_PT_NUM)*UEL(J)
            ENDDO
         ENDDO

         DO I=1,2
            STRAIN(I+6) = ZERO
            DO J=1,ELDOF
               STRAIN(I+6) = STRAIN(I+6) + BE3(I,J,STR_PT_NUM)*UEL(J)
            ENDDO
         ENDDO

                                                           ! Calc stresses from strains

         STRESS1(:)       = ZERO
         STRESS2(:)       = ZERO
         STRESS3(:)       = ZERO

         STRESS1_MECH(:)  = ZERO
         STRESS2_MECH(:)  = ZERO
         STRESS3_MECH(:)  = ZERO

         STRESS1_THERM(:) = ZERO
         STRESS2_THERM(:) = ZERO
         STRESS3_THERM(:) = ZERO

         DO I=1,3
            STRAIN1(I) = STRAIN(I)
            STRAIN2(I) = STRAIN(I+3)
            STRAIN3(I) = STRAIN(I+6)
         ENDDO

         IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
           TBAR = 0
           DO J=1,ELGP
             TBAR = TBAR + DT(J,JTSUB)
           ENDDO
           TBAR = TBAR / ELGP
           ALPTM(1) = ALPVEC(1  ,1)*(TBAR - TREF(1))
           ALPTM(2) = ALPVEC(2  ,1)*(TBAR - TREF(1))
           ALPTM(3) = ALPVEC(4  ,1)*(TBAR - TREF(1))
           DO I=1,3
             ALPTB(I) = ALPVEC(I  ,2)*DT(5,JTSUB)
             ALPTT(I) = ALPVEC(I+3,3)*(TBAR - TREF(1))
           ENDDO

                                                           ! MITC4/MITC4+ on a curved reference surface also has a free thermal
                                                           ! curvature from uniform mid-surface scaling, the same term subr MITC4
                                                           ! puts in the thermal load vector. It has to come out of the recovered
                                                           ! curvature here as well, or a free uniformly heated curved shell -
                                                           ! which deforms correctly precisely because the load vector carries the
                                                           ! term - reports a bending stress of [EB]{k_th} that is not there. The
                                                           ! stress point ordering matches the BEi third index that subr MITC4
                                                           ! filled, and THERM_KAPPA is zero for a flat element whose director is
                                                           ! normal to it, so flat plates are unchanged.
           IF ((TYPE(1:5) == 'QUAD4') .AND. ((QUAD4TYP == 'MITC4 ') .OR. (QUAD4TYP == 'MITC4+'))) THEN
             IF ((STR_PT_NUM >= 1) .AND. (STR_PT_NUM <= MSTRPT)) THEN
               DO I=1,3
                 ALPTB(I) = ALPTB(I) + THERM_KAPPA(I,STR_PT_NUM)*(TBAR - TREF(1))
               ENDDO
             ENDIF
           ENDIF
         ELSE
            ALPTM(:) = ZERO
            ALPTB(:) = ZERO
            ALPTT(:) = ZERO
         ENDIF


         ET3(:,:) = ZERO
         DO I=1,2
            DO J=1,2
               ET3(I,J) = ET(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( EM , STRAIN1, 3, 3, 1, DUM31 )
         CALL MATMULT_FFF ( EB , STRAIN2, 3, 3, 1, DUM32 )
         CALL MATMULT_FFF ( ET3, STRAIN3, 3, 3, 1, DUM33 )

         STRESS1_MECH =        DUM31
         STRESS2_MECH =        DUM32
         STRESS3_MECH = PHI_SQ*DUM33                       ! Need PHI_SQ on transv shear stress since this calc is from strains and
                                                           ! BE3, not SE3. If DEBUG(176) > 0 then stresses are calc'd from the SE3
                                                           ! below and SE3 has PHI_SQ incorporated in subrs QPLT1, QPLT3, TPLT2.


         IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
            CALL MATMULT_FFF ( EM , ALPTM  , 3, 3, 1, STRESS1_THERM )
            CALL MATMULT_FFF ( EB , ALPTB  , 3, 3, 1, STRESS2_THERM )
            CALL MATMULT_FFF ( ET3, ALPTT  , 3, 3, 1, STRESS3_THERM )
         ENDIF

         CALL MATADD_FFF  ( STRESS1_MECH, STRESS1_THERM, 3, 1, ONE, -ONE, 0, STRESS1 )
         CALL MATADD_FFF  ( STRESS2_MECH, STRESS2_THERM, 3, 1, ONE, -ONE, 0, STRESS2 )
         CALL MATADD_FFF  ( STRESS3_MECH, STRESS3_THERM, 3, 1, ONE, -ONE, 0, STRESS3 )

         DO I=1,3
            STRESS(I)   = STRESS1(I)
            STRESS(I+3) = STRESS2(I)
            STRESS(I+6) = STRESS3(I)
         ENDDO




         IF (DEBUG(176) > 0) THEN                          ! If DEBUG(176) > 0, calc stresses using SEi instead of above from STRAIN
            DO I=1,3                                       ! NOTE: PHI_SQ is incorporated into SE3
               STRESS(I  ) = ZERO
               STRESS(I+3) = ZERO
               DO J=1,ELDOF
                  STRESS(I)   = STRESS(I)   + SE1(I,J,STR_PT_NUM)*UEL(J)
                  STRESS(I+3) = STRESS(I+3) + SE2(I,J,STR_PT_NUM)*UEL(J)
               ENDDO
               IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
                  STRESS(I)   = STRESS(I)   - STE1(I,JTSUB,STR_PT_NUM)
                  STRESS(I+3) = STRESS(I+3) - STE2(I,JTSUB,STR_PT_NUM)
               ENDIF
            ENDDO
            DO I=1,2
               STRESS(I+6) = ZERO
               DO J=1,ELDOF
                  STRESS(I+6) = STRESS(I+6) + SE3(I,J,STR_PT_NUM)*UEL(J)
               ENDDO
               IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
                  STRESS(I+6) = STRESS(I+6) - STE3(I,JTSUB,STR_PT_NUM)
               ENDIF
            ENDDO
         ENDIF





! **********************************************************************************************************************************
! Calc strains, then stresses for 3D elements

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

         DO I=1,6
            STRESS_MECH(I)  = ZERO
            STRESS_THERM(I) = ZERO
         ENDDO

         DO I=1,3
            DO J=1,ELDOF
               STRAIN(I)   = STRAIN(I  ) + BE1(I,J,STR_PT_NUM)*UEL(J)
               STRAIN(I+3) = STRAIN(I+3) + BE2(I,J,STR_PT_NUM)*UEL(J)
            ENDDO
         ENDDO

         DO I=1,6
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
               TBAR = (DT(1,JTSUB) + DT(2,JTSUB) + DT(3,JTSUB) + DT(4,JTSUB))/FOUR
               ALPT(I) = ALPVEC(I,1)*(TBAR - TREF(1))
            ELSE
               ALPT(I) = ZERO
            ENDIF
         ENDDO

         CALL MATMULT_FFF ( ES, STRAIN, 6, 6, 1, STRESS_MECH )

         IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
            CALL MATMULT_FFF ( ES, ALPT, 6, 6, 1, STRESS_THERM )
         ENDIF

         CALL MATADD_FFF  ( STRESS_MECH, STRESS_THERM, 6, 1, ONE, -ONE, 0, STRESS )


      ENDIF

! **********************************************************************************************************************************
! Transform coord for STRESS/STRAIN arrays, if requested (and if for 2D or 3D elements)

      IF (STR_CID /= -1) THEN                         ! User req diff stress/strain/engr force output coord sys than elem local

! Warning. For STR_CID >= 0, the titles of stress and strain in .f06 still say
! L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M even when it's transformed here.
! STR_CID == -2 says M A T E R I A L   C O O R D I N A T E   S Y S T E M for solids.

         IF      ((TYPE (1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'TRIA3')) THEN

            IF (STR_CID /= -2) THEN
! Shells don't work because STR_TENSOR_TRANSFORM should be between setting STR_TENSOR and setting stress
! and shear strain may need factor of 2 before transforming.
               WRITE(ERR,9303)
               WRITE(F06,9303)
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

! Some of this code could be replaced with calls to TRANSFORM_SHELL_STR
                                                              ! Transform 2D membrane and transverse shear stresses
               STR_TENSOR(1,1) = STRESS(1)   ;   STR_TENSOR(1,2) = STRESS(3)   ;   STR_TENSOR(1,3) = STRESS(7)
               STR_TENSOR(2,1) = STRESS(3)   ;   STR_TENSOR(2,2) = STRESS(2)   ;   STR_TENSOR(2,3) = STRESS(8)
               STR_TENSOR(3,1) = STRESS(7)   ;   STR_TENSOR(3,2) = STRESS(8)   ;   STR_TENSOR(3,3) = ZERO

               STRESS(1) = STR_TENSOR(1,1)
               STRESS(2) = STR_TENSOR(2,2)
               STRESS(3) = STR_TENSOR(1,2)
               STRESS(7) = STR_TENSOR(1,3)
               STRESS(8) = STR_TENSOR(2,3)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID )
                                                              ! Transform 2D bending stresses
               STR_TENSOR(1,1) = STRESS(4)   ;   STR_TENSOR(1,2) = STRESS(6)   ;   STR_TENSOR(1,3) = ZERO
               STR_TENSOR(2,1) = STRESS(6)   ;   STR_TENSOR(2,2) = STRESS(5)   ;   STR_TENSOR(2,3) = ZERO
               STR_TENSOR(3,1) = ZERO        ;   STR_TENSOR(3,2) = ZERO        ;   STR_TENSOR(3,3) = ZERO

               STRESS(4) = STR_TENSOR(1,1)
               STRESS(5) = STR_TENSOR(2,2)
               STRESS(6) = STR_TENSOR(1,2)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID )
                                                              ! Transform 2D membrane and transverse shear strains
               STR_TENSOR(1,1) = STRAIN(1)   ;   STR_TENSOR(1,2) = STRAIN(3)   ;   STR_TENSOR(1,3) = STRAIN(7)
               STR_TENSOR(2,1) = STRAIN(3)   ;   STR_TENSOR(2,2) = STRAIN(2)   ;   STR_TENSOR(2,3) = STRAIN(8)
               STR_TENSOR(3,1) = STRAIN(7)   ;   STR_TENSOR(3,2) = STRAIN(8)   ;   STR_TENSOR(3,3) = ZERO

               STRAIN(1) = STR_TENSOR(1,1)
               STRAIN(2) = STR_TENSOR(2,2)
               STRAIN(3) = STR_TENSOR(1,2)
               STRAIN(7) = STR_TENSOR(1,3)
               STRAIN(8) = STR_TENSOR(2,3)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID )
                                                              ! Transform 2D bending strains
               STR_TENSOR(1,1) = STRAIN(4)   ;   STR_TENSOR(1,2) = STRAIN(6)   ;   STR_TENSOR(1,3) = ZERO
               STR_TENSOR(2,1) = STRAIN(6)   ;   STR_TENSOR(2,2) = STRAIN(5)   ;   STR_TENSOR(2,3) = ZERO
               STR_TENSOR(3,1) = ZERO        ;   STR_TENSOR(3,2) = ZERO        ;   STR_TENSOR(3,3) = ZERO

               STRAIN(4) = STR_TENSOR(1,1)
               STRAIN(5) = STR_TENSOR(2,2)
               STRAIN(6) = STR_TENSOR(1,2)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID )

            ENDIF

         ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

            IF (STR_CID == -2) THEN
               STR_CID_SOLID = ISOLID(3)                   ! CORDM from PSOLID card.
            ELSE
               STR_CID_SOLID = STR_CID
            ENDIF

            IF(STR_CID_SOLID == -1) then
                                                           ! STR_CID_SOLID is -1 if ISOLID(3) is -1 which means
                                                           ! material coordinates = element coordinates so don't
                                                           ! transform it.

            ELSE IF(STR_CID_SOLID >= 0) THEN

                                                           ! Transform 3D stresses
               STR_TENSOR(1,1) = STRESS(1)   ;   STR_TENSOR(1,2) = STRESS(4)   ;   STR_TENSOR(1,3) = STRESS(6)
               STR_TENSOR(2,1) = STRESS(4)   ;   STR_TENSOR(2,2) = STRESS(2)   ;   STR_TENSOR(2,3) = STRESS(5)
               STR_TENSOR(3,1) = STRESS(6)   ;   STR_TENSOR(3,2) = STRESS(5)   ;   STR_TENSOR(3,3) = STRESS(3)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID_SOLID )

               STRESS(1) = STR_TENSOR(1,1)
               STRESS(2) = STR_TENSOR(2,2)
               STRESS(3) = STR_TENSOR(3,3)
               STRESS(4) = STR_TENSOR(1,2)
               STRESS(5) = STR_TENSOR(2,3)
               STRESS(6) = STR_TENSOR(1,3)

                                                           ! Transform 3D strains
               STR_TENSOR(1,1) = STRAIN(1)   ;   STR_TENSOR(1,2) = STRAIN(4)/2 ;   STR_TENSOR(1,3) = STRAIN(6)/2
               STR_TENSOR(2,1) = STRAIN(4)/2 ;   STR_TENSOR(2,2) = STRAIN(2)   ;   STR_TENSOR(2,3) = STRAIN(5)/2
               STR_TENSOR(3,1) = STRAIN(6)/2 ;   STR_TENSOR(3,2) = STRAIN(5)/2 ;   STR_TENSOR(3,3) = STRAIN(3)

               CALL STR_TENSOR_TRANSFORM ( STR_TENSOR, STR_CID_SOLID )

               STRAIN(1) = STR_TENSOR(1,1)
               STRAIN(2) = STR_TENSOR(2,2)
               STRAIN(3) = STR_TENSOR(3,3)
               STRAIN(4) = STR_TENSOR(1,2)*2
               STRAIN(5) = STR_TENSOR(2,3)*2
               STRAIN(6) = STR_TENSOR(1,3)*2

            ELSE

               ! We don't know what to do for CORDM <= -2 because it's not defined in Mystran.
               WRITE(ERR,9304)
               WRITE(F06,9304)
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

         ELSE

                                                           ! Other element types are not an error for STR_CID = -2 or -1,
                                                           ! both of which mean leave them as they are.
            IF (STR_CID /= -2) THEN
               WRITE(ERR,9203) TYPE
               WRITE(F06,9203) TYPE
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

 9303 FORMAT(' *ERROR  9303: PARAM,STR_CID not implemented for QUAD and TRIA elements.' )

 9304 FORMAT(' *ERROR  9304: PSOLID field 4, CORDM <= -2 not allowed for solid elements.' )


! ##################################################################################################################################

      END SUBROUTINE ELEM_STRE_STRN_ARRAYS


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
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG

      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T
      USE VECTOR_GEOMETRY, ONLY       :  GEN_T0L
      USE ELEMENT_DIAGNOSTICS, ONLY   :  ELMOFF
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
         CALL GET_GRID_NUM_COMPS ( BGRID(1), NUM_COMPS_GRID_1, SUBR_NAME )
         I2 = NUM_COMPS_GRID_1 + ELAS_COMP(2)
         IF (ELGP == 1) THEN                               ! A grounded spring: one point
            PEL(I1) = KE(I1,I1)*UEL(I1)
         ELSE
            PEL(I1) = KE(I1,I1)*UEL(I1) + KE(I1,I2)*UEL(I2)! Note: KE is global and local for the ELAS elems
            PEL(I2) = KE(I2,I1)*UEL(I1) + KE(I2,I2)*UEL(I2)
         ENDIF

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


      SUBROUTINE GET_COMP_SHELL_ALLOWS ( STRE_ALLOWABLES, STRN_ALLOWABLES )

! Gets allowable stresses and strains for a composite element ply. Arrays ULT_STRE, ULT_STRN were formulated from user supplied
! data on the MATi Bulk Data entries in material processing subrs called by subr EMG.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE TIMDAT, ONLY                :  TSEC
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE MACHINE_PARAMS, ONLY        :  MACH_LARGE_NUM
      USE MODEL_STUF, ONLY            :  ULT_STRE, ULT_STRN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_COMP_SHELL_ALLOWS'



      REAL(DOUBLE), INTENT(OUT)       :: STRE_ALLOWABLES(9)! Stress allowables for the material
      REAL(DOUBLE), INTENT(OUT)       :: STRN_ALLOWABLES(9)! Strain allowables for the material



! **********************************************************************************************************************************
      STRE_ALLOWABLES(1) = ULT_STRE(1,1)                !   Axis   1 tension     stress allowable
      STRE_ALLOWABLES(2) = ULT_STRE(2,1)                !   Axis   1 compression stress allowable
      STRE_ALLOWABLES(3) = ULT_STRE(1,1)                !   Axis   2 tension     stress allowable
      STRE_ALLOWABLES(4) = ULT_STRE(2,1)                !   Axis   2 compression stress allowable
      STRE_ALLOWABLES(5) = MACH_LARGE_NUM               !   Axis   3 tension     stress allowable
      STRE_ALLOWABLES(6) = MACH_LARGE_NUM               !   Axis   3 compression stress allowable
      STRE_ALLOWABLES(7) = ULT_STRE(9,3)                !   Plane 23 shear       stress allowable (from transv shear matl props)
      STRE_ALLOWABLES(8) = ULT_STRE(8,3)                !   Plane 13 shear       stress allowable (from transv shear matl props)
      STRE_ALLOWABLES(9) = ULT_STRE(7,1)                !   Plane 12 shear       stress allowable

      STRN_ALLOWABLES(1) = ULT_STRN(1,1)                !   Axis   1 tension     strain allowable
      STRN_ALLOWABLES(2) = ULT_STRN(2,1)                !   Axis   1 compression strain allowable
      STRN_ALLOWABLES(3) = ULT_STRN(1,1)                !   Axis   2 tension     strain allowable
      STRN_ALLOWABLES(4) = ULT_STRN(2,1)                !   Axis   2 compression strain allowable
      STRN_ALLOWABLES(5) = MACH_LARGE_NUM               !   Axis   3 tension     strain allowable
      STRN_ALLOWABLES(6) = MACH_LARGE_NUM               !   Axis   3 compression strain allowable
      STRN_ALLOWABLES(7) = ULT_STRN(9,3)                !   Plane 23 shear       strain allowable (from transv shear matl props)
      STRN_ALLOWABLES(8) = ULT_STRN(8,3)                !   Plane 13 shear       strain allowable (from transv shear matl props)
      STRN_ALLOWABLES(9) = ULT_STRN(7,1)                !   Plane 12 shear       strain allowable



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_COMP_SHELL_ALLOWS



      SUBROUTINE POLYNOM_FIT_STRE_STRN ( STR_IN, NROW, NCOL, STR_OUT, STR_OUT_PCT_ERR, STR_OUT_ERR_INDEX, PCT_ERR_MAX )

! Extrapolates stress/strain values at the points at which the stress/strain matrices were calculated (i.e. Gauss points or MIN4T
! sub-tria centroids) to element corner nodes. Uses polynomial surface fit algorithm in module LSQ_MYSTRAN with the polynomial
! fit returned from subr SURFACE_FIT, called herein.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MAX_STRESS_POINTS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, THREE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EID, ELGP, TYPE, XEL
      USE PARAMS, ONLY                :  Q4SURFIT, QUAD4TYP

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE VECTOR_GEOMETRY, ONLY       :  PARAM_CORDS_ACT_CORDS, SURFACE_FIT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'POLYNOM_FIT_STRE_STRN'
      CHARACTER(100*BYTE)             :: MESSAGE             ! Message to send to subr SURFACE_FIT for debug output purposes

      INTEGER(LONG), INTENT(IN)       :: NCOL                ! Number of cols in arrays STR_IN, STR_OUT (1 more than the num of elem
!                                                              nodes since 1st col in arrays is due to elem center stress/strain
!                                                              which is not processed by this subr)
      INTEGER(LONG), INTENT(IN)       :: NROW                ! Number of rows in arrays STR_IN, STR_OUT

                                                             ! Stress/strain index number (1-9) where % error in fit is max
      INTEGER(LONG), INTENT(OUT)      :: STR_OUT_ERR_INDEX(MAX_STRESS_POINTS)
      INTEGER(LONG), PARAMETER        :: IORD = 2            ! Gaussian integration order
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: OUNT(2)             ! Output units for SURFACE_FIT
      INTEGER(LONG)                   :: SF_IERR             ! Output error indicator from subr SURFACE_FIT


      REAL(DOUBLE), INTENT(IN)        :: STR_IN(NROW,NCOL)   ! Input stress/strain vals. NROW are num of diff stress/strain vals and
!                                                              NCOL are number of points to use in the poly fit for one value

      REAL(DOUBLE), INTENT(OUT)       :: STR_OUT(NROW,NCOL)  ! Output stress/strain vals. NROW are num of diff stress/strain vals
!                                                              and NCOL are number of points to use in the poly fit for one value

      REAL(DOUBLE), INTENT(OUT)       :: PCT_ERR_MAX         ! Max value from array PCT_ERR (max poly fit err at any input data pt)

                                                             ! 2D array of POLY_PCT_ERR for all stress points
      REAL(DOUBLE)                    :: PCT_ERR(NROW,MAX_STRESS_POINTS)

                                                             ! % diff bet. fitted and actual input data for 1 stress or strain point
!                                                              Only NROW x NCOL-1 vals determined (or needed)
      REAL(DOUBLE)                    :: POLY_PCT_ERR(MAX_STRESS_POINTS)

                                                             ! 1st val is zero for center stress (not fitted) and remainder are for
!                                                              the points fitted (from row 2 to row NCOL; or NCOL-1 values)
      REAL(DOUBLE), INTENT(OUT)       :: STR_OUT_PCT_ERR(MAX_STRESS_POINTS)

      REAL(DOUBLE)                    :: PCT_ERR_1STR        ! % error from 1 of the NROW stress/strain values
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)! Gauss weights
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)! Gauss point coords
      REAL(DOUBLE)                    :: XI(NCOL-1)          ! X coords of the input  data points
      REAL(DOUBLE)                    :: YI(NCOL-1)          ! Y coords of the input  data points
      REAL(DOUBLE)                    :: WI(NCOL-1)          ! Values of the function to fit at the input data points
      REAL(DOUBLE)                    :: XO(NCOL-1)          ! X coords of the output data points
      REAL(DOUBLE)                    :: YO(NCOL-1)          ! Y coords of the output data points
      REAL(DOUBLE)                    :: XEA(NCOL-1,3)       ! Actual local element coords corresponding to XEP
      REAL(DOUBLE)                    :: XEP(NCOL-1,3)       ! Parametric coords of NCOL points
      REAL(DOUBLE)                    :: WO(NCOL-1)          ! Values of the function to fit at the output data points



! **********************************************************************************************************************************
      MESSAGE(1:) = ' '

      OUNT(1) = ERR
      OUNT(2) = F06

      DO J=1,MAX_STRESS_POINTS
         POLY_PCT_ERR(J) = ZERO
         STR_OUT_PCT_ERR(J) = ZERO
         DO I=1,NROW
            PCT_ERR(I,J) = ZERO
         ENDDO
      ENDDO

      DO J=1,NCOL
         STR_OUT_ERR_INDEX(J) = 0
         DO I=1,NROW
            STR_OUT(I,J) = STR_IN(I,J)
         ENDDO
      ENDDO

! Calc actual coords of the points for which the BEi, SEi matrices were calculated

      IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8')) THEN

         IF (NCOL /= 5) THEN                               ! Number of stress/strain points = number of corner points+1
            WRITE(ERR,9202) SUBR_NAME, TYPE, NCOL, 4+1
            WRITE(F06,9202) SUBR_NAME, TYPE, NCOL, 4+1
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
                                                           ! For the MIN4/MITC4+ QUAD4 and QUAD8. XEP parametric coords
                                                           ! are the Gauss point coords
         IF ((QUAD4TYP == 'MIN4  ') .OR.                                                                                           &
             (QUAD4TYP == 'MITC4 ') .OR.                                                                                           &
             (QUAD4TYP == 'MITC4+') .OR.                                                                                           &
             (TYPE(1:6) == 'QUAD4K') .OR.                                                                                          &
             (TYPE(1:5) == 'QUAD8')) THEN
            CALL ORDER_GAUSS ( IORD, SSS, HHH )
            XEP(1,1) =  SSS(1)      ;   XEP(1,2) =  SSS(1)      ;   XEP(1,3) = ZERO
            XEP(2,1) =  SSS(1)      ;   XEP(2,2) =  SSS(2)      ;   XEP(2,3) = ZERO
            XEP(3,1) =  SSS(2)      ;   XEP(3,2) =  SSS(1)      ;   XEP(3,3) = ZERO
            XEP(4,1) =  SSS(2)      ;   XEP(4,2) =  SSS(2)      ;   XEP(4,3) = ZERO
         ELSE                                              ! For the MIN4T QUAD4, XEP parametric coords are centroids of the 4 trias
            XEP(1,1) =  ZERO        ;   XEP(1,2) = -TWO/THREE   ;   XEP(1,3) = ZERO
            XEP(2,1) =  TWO/THREE   ;   XEP(2,2) =  ZERO        ;   XEP(2,3) = ZERO
            XEP(3,1) =  ZERO        ;   XEP(3,2) =  TWO/THREE   ;   XEP(3,3) = ZERO
            XEP(4,1) = -TWO/THREE   ;   XEP(4,2) =  ZERO        ;   XEP(4,3) = ZERO
         ENDIF

         CALL PARAM_CORDS_ACT_CORDS ( 4, IORD, XEP, XEA )  ! Get actual local elem coords corresponding to the parametric coords
         DO I=1,NCOL-1
            XI(I) = XEA(I,1)                               ! XEA are act coords of the XEP.        XI is input to subr SURFACE_FIT
            YI(I) = XEA(I,2)
            XO(I) = XEL(I,1)                               ! XEL are act coords of the elem nodes. XO is input to subr SURFACE_FIT
            YO(I) = XEL(I,2)
         ENDDO

! Now do surface fit on the data (STR_IN) to get stress/strain at the elem corners (STR_OUT)

         PCT_ERR_MAX = ZERO
         DO I=1,NROW
            DO J=2,NCOL
               WI(J-1) = STR_IN(I,J)                       ! The input stress/strain are the values at the Gauss or tria centroids
               STR_OUT_PCT_ERR(J) = ZERO
            ENDDO
            WRITE(MESSAGE( 1:16),'(A )') ' for stress row '
            WRITE(MESSAGE(17:18),'(I2)') I
            WRITE(MESSAGE(19:22),'(A )') ' of '
            WRITE(MESSAGE(23:24),'(I2)') NROW
            WRITE(MESSAGE(25:41),'(A )') ' for MIN4T QUAD4 '
            WRITE(MESSAGE(42:49),'(I8)') EID
            CALL SURFACE_FIT ( NCOL-1, Q4SURFIT, XI, YI, WI, XO, YO, WO, DEBUG(175), MESSAGE, OUNT, POLY_PCT_ERR, PCT_ERR_1STR,    &
                               SF_IERR )
            IF (DABS(PCT_ERR_1STR) > DABS(PCT_ERR_MAX)) THEN
               PCT_ERR_MAX = PCT_ERR_1STR
            ENDIF
            DO J=2,NCOL
               STR_OUT(I,J) = WO(J-1)                      ! These are the stress/strain values extrapolated to the element corners
               PCT_ERR(I,J) = POLY_PCT_ERR(J-1)
            ENDDO
         ENDDO
                                                           ! Sort over rows of PCT_ERR to get largest abs val for all stress/strain
         DO J=2,NCOL                                       ! values for each stress point
            DO I=1,NROW
               IF (DABS(PCT_ERR(I,J)) > DABS(STR_OUT_PCT_ERR(J))) THEN
                  STR_OUT_PCT_ERR(J)   = PCT_ERR(I,J)
                  STR_OUT_ERR_INDEX(J) = I
               ENDIF
            ENDDO
         ENDDO

      ELSE

         RETURN

      ENDIF



      RETURN

! **********************************************************************************************************************************
            WRITE(F06,9202) SUBR_NAME, TYPE, NCOL, ELGP

 9202 FORMAT(' *ERROR  9202: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' FOR ELEMENT TYPE = "',A,'" THE INPUT VALUE OF ARGUMENT NCOL WAS = ',I8,' BUT MUST BE ELGP+1 = ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE POLYNOM_FIT_STRE_STRN



      SUBROUTINE SHELL_ENGR_FORCE_OGEL ( NUM1 )

! Calculates element engineering forces for plate elements from array STRESS (generated in subr ELEM_STRE_STRN_ARRAYS) using FCONV
! (conversion factor from stress to engr force)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  FCONV, STRESS
      USE LINK9_STUFF, ONLY           :  MAXREQ, MAXREQ, OGEL

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SHELL_ENGR_FORCE_OGEL'

      INTEGER(LONG), INTENT(INOUT)    :: NUM1              ! Cum rows written to OGEL prior to running this subr
      INTEGER(LONG)                   :: I                 ! DO loop indices




! **********************************************************************************************************************************
      NUM1 = NUM1 + 1
      IF (NUM1 > MAXREQ) THEN
         WRITE(ERR,9200) SUBR_NAME,MAXREQ
         WRITE(F06,9200) SUBR_NAME,MAXREQ
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (dim of array OGEL too small), so quit
      ENDIF
      DO I=1,3
         OGEL(NUM1,I) = FCONV(1)*STRESS(I)
      ENDDO
      DO I=4,6
         OGEL(NUM1,I) = FCONV(2)*STRESS(I)
      ENDDO
      DO I=7,9
         OGEL(NUM1,I) = FCONV(3)*STRESS(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

! **********************************************************************************************************************************

      END SUBROUTINE SHELL_ENGR_FORCE_OGEL

   END MODULE ELEMENT_RECOVERY_SUPPORT
