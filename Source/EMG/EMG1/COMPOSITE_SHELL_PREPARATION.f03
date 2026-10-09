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

   MODULE COMPOSITE_SHELL_PREPARATION

   USE ELEMENT_LOOKUPS, ONLY :  GET_MATANGLE_FROM_CID

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SHELL_ABD_MATRICES, IS_ELEM_PCOMP_PROPS, GET_PCOMP_SECT_PROPS, ROT_COMP_ELEM_AXES, GET_ELEM_NUM_PLIES

   CONTAINS

      SUBROUTINE SHELL_ABD_MATRICES ( INT_ELEM_ID, WRITE_WARN )

! Generates shell force resultant vs strain matrices to be used in generating ME, KE, PTE, PPE, SEi for shell elements (plate
! elements whose properties are specified on a Bulk Data PSHELL or PCOMP entry). This subroutine is run under two circumstances:

!   1) To get integrated effect of all plies in a shell element (QUAD4 or TRIA3) used in development of the overall stiffness and
!      mass matrices for the element (NOTE: an element using PSHELL properties is considered as a 1 "ply" element)

!   2) To get the individual matrices for a single ply of the element used for stress/strain calcs for that ply

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, ERR, F06, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEMATC, MRMATLC, MPCOMP_PLIES, MPCOMP0, MRPCOMP_PLIES, MRPCOMP0, &
                                         WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, THIRD, HALF, THREE, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL, IORQ1M, QUAD4TYP, PCOMPEQ, PCMPTSTM, SHRFXFAC, SUPWARN, TSTM_DEF

      USE MODEL_STUF, ONLY            :  ALPVEC, EB, EBM, EM, ET, EDAT, EID, EMAT, EPNT, EPROP, ETYPE, FAILURE_THEORY, FCONV,      &
                                         INTL_MID, INTL_PID, MASS_PER_UNIT_AREA, MATL, MEPROP, MTRL_TYPE,                          &
                                         NUM_EMG_FATAL_ERRS, NUM_PLIES, PLY_NUM, PCOMP, PCOMP_LAM, PCOMP_PROPS, RPCOMP, PSHEL,     &
                                         RPSHEL, RHO, RMATL, SHELL_A, SHELL_B, SHELL_D, SHELL_T, SHELL_AALP, SHELL_BALP,           &
                                         SHELL_DALP, SHELL_TALP, SHELL_T_MOD, THETA_PLY, TPLY, TYPE, ULT_STRE, ULT_STRN, ZPLY, ZS


      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE MATERIAL_PROPERTIES, ONLY   :  MATERIAL_PROPS_2D

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SHELL_ABD_MATRICES'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

! Variables common to homogeneous and composite shell elements

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID for which
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices


      REAL(DOUBLE)                    :: DET_SHELL_T        ! Determinant of SHELL_T
      REAL(DOUBLE)                    :: EPS1               ! Small number with which to comapre zero
      REAL(DOUBLE)                    :: NSM                ! Nonstructural mass

! Variables for homogeneous shell elements

      REAL(DOUBLE)                    :: IB                 ! Bending moment of inertia
      REAL(DOUBLE)                    :: TM                 ! Membrane thickness
      REAL(DOUBLE)                    :: TS                 ! Shear thickness

! Variables for composite elements

      INTEGER(LONG)                   :: FT                 ! Failure theory (1=HILL, 2=HOFF, 3=TSAI, 4=STRN)
      INTEGER(LONG)                   :: JPLY               ! = either PLY_NUM or K in the DO loop over K=1,NUM_PLIES_TO_PROC
      INTEGER(LONG)                   :: MTRL_ACT_ID(4)     ! Material ID from MATi B.D. entries, in MATL(INTL_MID(I),1)
      INTEGER(LONG)                   :: NUM_PLIES_TO_PROC  ! = 1 if we are processing only 1 ply or NUM_PLIES if processing all
      INTEGER(LONG)                   :: PLY_PCOMP_INDEX    ! Index in array  PCOMP where data for ply K begins
      INTEGER(LONG)                   :: PLY_RPCOMP_INDEX   ! Index in array RPCOMP where data for ply K begins
      INTEGER(LONG)                   :: SOUTK              ! Stress or strain output request (1=YES or 0=NO) for ply K

      REAL(DOUBLE)                    :: ALPB(3)            ! The 3 rows of ALPVEC for mem/bend  strains
      REAL(DOUBLE)                    :: ALPD(3)            ! The 3 rows of ALPVEC for bending   strains
      REAL(DOUBLE)                    :: ALPM(3)            ! The 3 rows of ALPVEC for membrane  strains
      REAL(DOUBLE)                    :: ALPT(3)            ! The 3 rows of ALPVEC for trans shr strains
      REAL(DOUBLE)                    :: BALP(3)            ! Intermediate matrix
      REAL(DOUBLE)                    :: DALP(3)            ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM(3)             ! Intermediate matrix
      REAL(DOUBLE)                    :: AALP(3)            ! Intermediate matrix
      REAL(DOUBLE)                    :: TALP(3)            ! Intermediate matrix
      REAL(DOUBLE)                    :: PCOMP_TM           ! Membrane thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE)                    :: PCOMP_IB           ! Bending MOI of PCOMP for equivalent PSHELL
      REAL(DOUBLE)                    :: PCOMP_TS           ! Transverse shear thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE)                    :: PLY_A(3,3)         ! Transformed material matrix A for a ply
      REAL(DOUBLE)                    :: PLY_B(3,3)         ! Transformed material matrix B for a ply
      REAL(DOUBLE)                    :: PLY_D(3,3)         ! Transformed material matrix D for a ply
      REAL(DOUBLE)                    :: PLY_T(2,2)         ! Transformed material matrix T for a ply
      REAL(DOUBLE)                    :: SB                 ! Allowable interlaminar shear stress. Required if FT is specified
      REAL(DOUBLE)                    :: TREFK              ! Ref temperature for ply K
      REAL(DOUBLE)                    :: Z0                 ! Coord from ref plane to bottom surface of element
      REAL(DOUBLE)                    :: ZBK,ZTK            ! Coord from ref plane to bot and top of ply K
      REAL(DOUBLE)                    :: ZBK2,ZTK2          ! ZBK^2, ZTK^2
      REAL(DOUBLE)                    :: ZBK3,ZTK3          ! ZBK^3, ZTK^3



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      TYPE  = ETYPE(INT_ELEM_ID)

      IF ((TYPE(1:5) /= 'TRIA3') .AND. (TYPE(1:5) /= 'QUAD4') .AND. (TYPE(1:5) /= 'QUAD8') .AND. (TYPE(1:5) /= 'SHEAR')) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1946) TYPE, SUBR_NAME
         WRITE(F06,1946) TYPE, SUBR_NAME
         RETURN
      ENDIF

      CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )

      SHELL_T_MOD = 'N'                                    ! Reset this to N every time this subr is called

! ---------------------------------------------------------------------------------------------------------------------------------
pcom0:IF (PCOMP_PROPS == 'N') THEN                         ! Element is not a composite - uses PSHELL
!                                                            NOTE: MATERIAL_PROPS_2D is called in procedure EMG for these
         IF (TYPE == 'SHEAR   ') THEN
            TM       = EPROP(1)
            IB       = ZERO
            TS       = TM
            NSM      = EPROP(2)
            ZS(1)    = ZERO
            ZS(2)    = ZERO
            FCONV(1) = TM
            FCONV(2) = ZERO
            FCONV(3) = ZERO
         ELSE
            TM                =  EPROP(1)
            IB                =  EPROP(2)*TM*TM*TM/TWELVE
            TS                =  EPROP(3)*TM
            NSM               =  EPROP(4)
            ZS(1)             =  EPROP(5)
            ZS(2)             =  EPROP(6)
            FCONV(1)          =  TM
            FCONV(2)          = -IB                        ! Note neq sign on FCONV(2): due to sign convention on positive bending
            FCONV(3)          =  TS
         ENDIF

         IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA3')) THEN
            MASS_PER_UNIT_AREA = (RHO(1)*TM + NSM)
         ENDIF

         DO I=1,3
            DO J=1,3
               SHELL_A(I,J) = ZERO
               SHELL_D(I,J) = ZERO
               SHELL_B(I,J) = ZERO
            ENDDO
         ENDDO
         IF (TYPE == 'SHEAR   ') THEN
            SHELL_A(3,3) = TM*EM(3,3)                      ! For SHEAR, all but SHELL_A(3,3) are 0.
         ELSE
            DO I=1,3
               DO J=1,3
                  SHELL_A(I,J) = TM*EM(I,J)                ! Units are force/distance
                  SHELL_D(I,J) = IB*EB(I,J)                ! Units are force*distance^2
                  SHELL_B(I,J) = TM*TM*EBM(I,J)            ! Units are force
               ENDDO
            ENDDO
         ENDIF

         DO I=1,2
            DO J=1,2
               SHELL_T(I,J) = TS*ET(I,J)                   ! Units are force/distance
            ENDDO
         ENDDO

         DO I=1,3
            ALPB(I) = ALPVEC(I,2)
         ENDDO

         ALPM(1) = ALPVEC(1,1)                             ! xx
         ALPM(2) = ALPVEC(2,1)                             ! yy
         ALPM(3) = ALPVEC(4,1)                             ! xy

         CALL MATMULT_FFF ( EB, ALPB, 3, 3, 1, DUM )
         DO I=1,3
            SHELL_DALP(I) = IB*DUM(I)
         ENDDO

         CALL MATMULT_FFF ( EM, ALPM, 3, 3, 1, DUM )
         DO I=1,3
            SHELL_AALP(I) = TM*DUM(I)
         ENDDO

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE pcom0                                           ! This element is a composite with properties defined on PCOMP
!                                                            NOTE: subr MATERIAL_PROPS_2D called here for these
         DO I=1,MEPROP
            EPROP(I) = ZERO
         ENDDO

         DO I=1,3
            ALPB(I) = ZERO
            ALPD(I) = ZERO
            ALPM(I) = ZERO
            ALPT(I) = ZERO
            AALP(I) = ZERO
            BALP(I) = ZERO
            DALP(I) = ZERO
            TALP(I) = ZERO
         ENDDO

         EPROP(4) =  RPCOMP(INTL_PID,1)    ;    Z0    = EPROP(4)
         EPROP(5) =  RPCOMP(INTL_PID,2)    ;    NSM   = EPROP(5)
         EPROP(6) =  RPCOMP(INTL_PID,3)    ;    SB    = EPROP(6)
         EPROP(7) =  RPCOMP(INTL_PID,4)    ;    TREFK = EPROP(7)
         EPROP(8) =  Z0                    ;    ZS(1) = EPROP(8)
         EPROP(9) = -Z0                    ;    ZS(2) = EPROP(9)

         FT = PCOMP(INTL_PID,3)
         IF      (FT == 0) THEN
            FAILURE_THEORY = 'NONE'
         ELSE IF (FT == 1) THEN
            FAILURE_THEORY = 'HILL'
         ELSE IF (FT == 2) THEN
            FAILURE_THEORY = 'HOFF'
         ELSE IF (FT == 3) THEN
            FAILURE_THEORY = 'TSAI'
         ELSE IF (FT == 4) THEN
            FAILURE_THEORY = 'STRE'
         ELSE IF (FT == 5) THEN
            FAILURE_THEORY = 'STRN'
         ENDIF

         IF (PCOMP(INTL_PID,4) == 1) THEN                  ! Check whether elem is sym or nonsym layuo
            PCOMP_LAM = 'SYM'
         ELSE
            PCOMP_LAM = 'NON'                              ! If nonsym layup, make sure int order = 2 (BIG_BB, BIG_BM for QUAD)
            IF ((TYPE(1:6) == 'QUAD4 ') .AND. (QUAD4TYP == 'MIN4T ')) THEN
               IF (IORQ1M /= 2) THEN
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR, 1948) 'IORQ1M', IORQ1M
                  IF (SUPWARN == 'N') THEN
                     WRITE(F06, 1948) 'IORQ1M', IORQ1M
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         DO I=1,MEMATC
            INTL_MID(I) = 0
         ENDDO

         CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )

         IF (WRT_BUG(1) > 0) THEN
            CALL BUG_SHELL_ABD_MATRICES ( 0, 20 )
         ENDIF

         DO I=1,3
            DO J=1,3
               SHELL_A(I,J) = ZERO
               SHELL_B(I,J) = ZERO
               SHELL_D(I,J) = ZERO
            ENDDO
         ENDDO

         DO I=1,2
            DO J=1,2
               SHELL_T(I,J) = ZERO
            ENDDO
         ENDDO

         DO I=1,3
            SHELL_AALP(I) = ZERO
            SHELL_BALP(I) = ZERO
            SHELL_DALP(I) = ZERO
         ENDDO

         DO I=1,2
            SHELL_TALP(I) = ZERO
         ENDDO

         MASS_PER_UNIT_AREA = ZERO
         MASS_PER_UNIT_AREA = NSM
         ZBK = Z0
         PCOMP_TM = ZERO
         PCOMP_IB = ZERO
         PCOMP_TS = ZERO

         IF (PLY_NUM == 0) THEN                            ! PLY_NUM = 0 means we want integrated effect of all plies
            NUM_PLIES_TO_PROC = NUM_PLIES
         ELSE                                              ! PLY_NUM > 0 means process only ply number PLY_NUM
            NUM_PLIES_TO_PROC = 1
         ENDIF

ply_do:  DO K=1,NUM_PLIES_TO_PROC

            IF (PLY_NUM == 0) THEN
               JPLY = K
            ELSE
               JPLY = PLY_NUM
            ENDIF
                                                           ! Indices in PCOMP, RPCOMP arrays where ply K data begins
            PLY_PCOMP_INDEX  = MPCOMP0  + MPCOMP_PLIES*(JPLY - 1)  + 1
            PLY_RPCOMP_INDEX = MRPCOMP0 + MRPCOMP_PLIES*(JPLY - 1) + 1

            INTL_MID(1) = PCOMP(INTL_PID,PLY_PCOMP_INDEX)  ! Shell A, D, T use membrane material props

            INTL_MID(2) = INTL_MID(1)                      ! Need to have INTL_MID(2) nonzero so that subr QDEL1,TREL1
!                                                            will call subr to calc bending stiffness. Also SEi depends on EB

            INTL_MID(3) = INTL_MID(1)                      ! Shell T uses transverse shear material props unless material matrix
!                                                            calcuated below has zero transverse shear modulus in which case
!                                                            INTL_MID(3) will be reset to 0

            INTL_MID(4) = INTL_MID(1)                      ! Bending/membrane coupling

            SOUTK     =  PCOMP(INTL_PID,PLY_PCOMP_INDEX+1)
            TPLY      = RPCOMP(INTL_PID,PLY_RPCOMP_INDEX)  ! Ply thickness
            THETA_PLY = RPCOMP(INTL_PID,PLY_RPCOMP_INDEX+1)! Ply angle from elem material axis to ply longitudinal axis
            ZPLY      = RPCOMP(INTL_PID,PLY_RPCOMP_INDEX+2)! Coord of mid plane of ply relative to mid plane of elem

            IF      ((TYPE == 'QDMEM   ') .OR. (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD4   ') .OR.                                 &
                     (TYPE == 'TRMEM   ') .OR. (TYPE == 'TRIA3K  ') .OR. (TYPE == 'TRIA3   ')) THEN
               MASS_PER_UNIT_AREA = MASS_PER_UNIT_AREA + (RHO(1)*TPLY)
            ELSE IF ((TYPE == 'QDPLT1  ') .OR. (TYPE == 'QDPLT2  ') .OR.                                                           &
                     (TYPE == 'TRPLT1  ') .OR. (TYPE == 'TRPLT2  ')) THEN
               MASS_PER_UNIT_AREA = MASS_PER_UNIT_AREA
            ENDIF

            IF (INTL_MID(1) /= 0) THEN
               MTRL_ACT_ID(1) = MATL(INTL_MID(1),1)
               MTRL_TYPE(1)   = MATL(INTL_MID(1),2)
            ENDIF

            IF (INTL_MID(2) /= 0) THEN
               MTRL_ACT_ID(2) = MATL(INTL_MID(2),1)
               MTRL_TYPE(2)   = MATL(INTL_MID(2),2)
            ENDIF

            IF (INTL_MID(3) /= 0) THEN
               MTRL_ACT_ID(3) = MATL(INTL_MID(1),1)
               MTRL_TYPE(3)   = MATL(INTL_MID(3),2)
            ENDIF

            DO I=1,MRMATLC
               DO J=1,MEMATC
                  IF (INTL_MID(J) /= 0) THEN
                     EMAT(I,J) = RMATL(INTL_MID(J),I)
                  ENDIF
               ENDDO
            ENDDO

            IF (PCOMP(INTL_PID,2) == 0) THEN               ! Props were defined on a PCOMP entry so we need to modify TREF to be
               DO J=1,MEMATC                               ! from array PCOMP rather than from RMATL. OK if defined on PCOMP1
                  IF      (MTRL_TYPE(J) == 1) THEN
                     EMAT( 6,J) = RPCOMP(INTL_PID,4)
                  ELSE IF (MTRL_TYPE(J) == 2) THEN
                     EMAT(11,J) = RPCOMP(INTL_PID,4)
                  ELSE IF (MTRL_TYPE(J) == 8) THEN
                     EMAT(10,J) = RPCOMP(INTL_PID,4)
                  ENDIF
               ENDDO
            ENDIF

            EMAT(MRMATLC+1,3) = SB                         ! SB is the allowable interlaminar shear stress (on PCOMP B.D. entry)
            EMAT(MRMATLC+2,3) = SB
            CALL MATERIAL_PROPS_2D ( WRITE_WARN )
                                                           ! Reset INTL_MID(3) if either transverse shear modulii are zero
            IF ((DABS(ET(1,1)) < EPS1) .OR. (DABS(ET(2,2)) < EPS1)) THEN
               INTL_MID(3) = 0
            ENDIF


            IF (WRT_BUG(1) > 0) THEN
               CALL BUG_SHELL_ABD_MATRICES ( K, 21 )
            ENDIF

            CALL ROT_COMP_ELEM_AXES ( INT_ELEM_ID, K, THETA_PLY, '1-2' )

            ZTK = ZBK + TPLY    ;    ZTK2 = ZTK*ZTK    ;    ZTK3 = ZTK2*ZTK    ;    ZBK2 = ZBK*ZBK    ;    ZBK3 = ZBK2*ZBK

            PCOMP_TM = PCOMP_TM +       (ZTK  - ZBK )
            PCOMP_IB = PCOMP_IB + THIRD*(ZTK3 - ZBK3)
                                                           ! PCMPTSTM is a factor (< 1) for shear to total plate thickness for PCOMP
            PCOMP_TS = PCOMP_TS +       (ZTK  - ZBK )*PCMPTSTM

            IF (TYPE == 'SHEAR   ') THEN                   ! For SHEAR elem there is only SHELL_A and only its 3,3 term is nonzreo
               DO I=1,3
                  DO J=1,3
                     SHELL_A(I,J) = ZERO
                     SHELL_B(I,J) = ZERO
                     SHELL_D(I,J) = ZERO
                  ENDDO
               ENDDO
               SHELL_A(3,3) = SHELL_A(3,3) + PLY_A(3,3)
            ELSE
               DO I=1,3
                  DO J=1,3
                     IF (PLY_NUM == 0) THEN                ! Use the following when all plies are to be integrated
                        PLY_A(I,J)   =       (ZTK  - ZBK )*EM(I,J)
                        PLY_B(I,J)   =  HALF*(ZTK2 - ZBK2)*EM(I,J)
                        PLY_D(I,J)   = THIRD*(ZTK3 - ZBK3)*EM(I,J)
                     ELSE                                  ! Use the following when only individual plies are evaluated separately
                        PLY_A(I,J)   = TPLY*EM(I,J)
                        PLY_B(I,J)   = ZERO
                        PLY_D(I,J)   = (TPLY*TPLY*TPLY)*EM(I,J)/TWELVE
                     ENDIF
                     SHELL_A(I,J) = SHELL_A(I,J) + PLY_A(I,J)
                     SHELL_B(I,J) = SHELL_B(I,J) + PLY_B(I,J)
                     SHELL_D(I,J) = SHELL_D(I,J) + PLY_D(I,J)
                  ENDDO
               ENDDO
            ENDIF

            DO I=1,2
               DO J=1,2
                  PLY_T(I,J)   = (ZTK - ZBK )*ET(I,J)*PCMPTSTM
                  SHELL_T(I,J) = SHELL_T(I,J) + PLY_T(I,J)
               ENDDO
            ENDDO

            ALPM(1) = ALPVEC(1,1)
            ALPM(2) = ALPVEC(2,1)
            ALPM(3) = ALPVEC(4,1)

            DO I=1,3
               ALPD(I) = ALPVEC(I,2)
               ALPB(I) = ALPVEC(I,4)
            ENDDO

            DO I=1,2
               ALPT(I) = ALPVEC(I,3)
            ENDDO

            CALL MATMULT_FFF ( PLY_A, ALPM, 3, 3, 1, AALP )
            CALL MATMULT_FFF ( PLY_B, ALPB, 3, 3, 1, BALP )
            CALL MATMULT_FFF ( PLY_D, ALPD, 3, 3, 1, DALP )
            CALL MATMULT_FFF ( PLY_T, ALPT, 2, 2, 1, TALP )
            DO I=1,3
               SHELL_AALP(I) = SHELL_AALP(I) + AALP(I)
               SHELL_BALP(I) = SHELL_BALP(I) + BALP(I)
               SHELL_DALP(I) = SHELL_DALP(I) + DALP(I)
            ENDDO
            DO I=1,2
               SHELL_TALP(I) = SHELL_TALP(I) + TALP(I)
            ENDDO

            IF (WRT_BUG(1) > 0) THEN
               CALL BUG_SHELL_ABD_MATRICES ( K, 22 )
            ENDIF

            ZBK = ZTK

         ENDDO ply_do

! Now put PCOMP_TM, PCOMP_IB and PCOMP_TS into array EPROP, but put them in similar to the EPROP for PSHELL properties

         EPROP(1) = PCOMP_TM                               ! PCOMP_TM should be > 0 based on checks in subr BD_PCOMP
         EPROP(2) = TWELVE*PCOMP_IB/(PCOMP_TM*PCOMP_TM*PCOMP_TM)
         EPROP(3) = PCOMP_TS/PCOMP_TM

         IF (WRT_BUG(1) > 0) THEN
            CALL BUG_SHELL_ABD_MATRICES ( 0, 23 )
         ENDIF

         FCONV(1)  =  PCOMP_TM
         FCONV(2)  = -PCOMP_IB                             ! Note neq sign on FCONV(2): due to sign convention on positive bending
         FCONV(3)  =  PCOMP_TS
                                                           ! Write equiv PSHELL, MAT2, if not already written for this PCOMP
         IF ((PCOMPEQ > 0) .AND. (PCOMP(INTL_PID,6) == 0)) THEN
            IF ((PCOMP_TM > EPS1) .AND. (PCOMP_IB > EPS1) .AND. (PCOMP_TS > EPS1)) THEN
               CALL WRITE_PCOMP_EQUIV ( PCOMP_TM, PCOMP_IB, PCOMP_TS )
            ELSE
               WRITE(ERR,9800) PCOMP(INTL_PID,1), PCOMP_TM, PCOMP_IB, PCOMP_TS
               WRITE(F06,9800) PCOMP(INTL_PID,1), PCOMP_TM, PCOMP_IB, PCOMP_TS
            ENDIF
         ENDIF

      ENDIF pcom0

! Reset SHELL_T if singular so that we can calc finite KS shear stiffness in subrs QPLT2, TPLT2.

      DET_SHELL_T = SHELL_T(1,1)*SHELL_T(2,2) - SHELL_T(1,2)*SHELL_T(2,1)
      IF (DABS(DET_SHELL_T) < EPS1) THEN
         SHELL_T(1,1) = HALF*SHRFXFAC*( SHELL_A(1,1) + SHELL_A(2,2) )
         SHELL_T(2,2) = SHELL_T(1,1)
         SHELL_T(1,2) = ZERO
         SHELL_T(2,1) = ZERO
         WRITE(ERR,9801) EID, SHELL_T(1,1), SHELL_T(2,2)
!        WRITE(F06,9801) EID, SHELL_T(1,1), SHELL_T(2,2)
         SHELL_T_MOD = 'Y'
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1946 FORMAT(' *ERROR  1946: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE ',A,' IS NOT ONE OF THE ELEMENT TYPES PROCESSED IN SUBR ',A)

 1948 FORMAT(' *ERROR  1948: ',A,I8,' MUST HAVE GAUSSIAN INTEGRATION ORDERS = 2 FOR PCOMP ELEMENTS WITH SYM LAYUP.'                &
                      ,/,14X,' HOWEVER, BULK DATA PARAM ',A,' = ',I3)
 9800 FORMAT(' *INFORMATION: CANNOT OUTPUT EQUIVALENT PSHELL AND MAT2 ENTRIES FOR PCOMP ',I8,' SINCE ONE OR MORE OF THE EQUIV',    &
                           ' PROPERTIES IS ZERO:'                                                                                  &
                    ,/,14X,'    (1) EQUIVALENT MEMBRANE THICKNESS = ',1ES10.3                                                      &
                    ,/,14X,'    (2) EQUIVALENT BENDING INERTIA    = ',1ES10.3                                                      &
                    ,/,14X,'    (3) EQUIVALENT SHEAR THICKNESS    = ',1ES10.3,/)

 9801 FORMAT(' *INFORMATION: ARRAY SHELL_T USED IN FORMULATING TRANSVERSE SHEAR STIFFNESS FOR ELEMENT ', I8,' HAS BEEN RESET DUE', &
             ' TO SINGULARITY.',/, 14X,' THIS RESET COULD HAVE BEEN AVOIDED BY SPECIFYING TRANSVERSE SHEAR MATERIAL PROPERTIES',   &
             ' FOR THIS ELEMENT.',/,14X,' THE RESET VALUES FOR THE DIAGONALS OF THE 2 BY 2 SHELL_T',     &
             ' MATRIX FOR THIS ELEMENT ARE:', 2(1ES14.6),/,14X,' IN AN ATTEMPT TO SIMULATE ZERO TRANSVERSE SHEAR FLEXIBILITY'/)



49832 format(' In SHELL_ABD_MATRICES: TM, IB, TS, FCONV(1), FCONV(2) = ',5(1es14.6))













91304 FORMAT('   SHELL_D row',I2,' = ',3(1ES14.6),'  bending              ')



! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE BUG_SHELL_ABD_MATRICES ( JPLY, WHAT )

      USE MODEL_STUF, ONLY            :  T1M, T1T

      IMPLICIT NONE

      CHARACTER( 8*BYTE)              :: LAM = '        '  ! Field 9 of parent entry
      CHARACTER(4*BYTE)               :: MTRL_NAME         ! Material name (MAT1, etc)

      INTEGER(LONG), INTENT(IN)       :: JPLY              ! Integer ply number
      INTEGER(LONG)                   :: IIROW,IICOL       !
      INTEGER(LONG)                   :: J1,J2             ! Counters
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in 1 PCOMP entry incl sym plies not explicitly defined
      INTEGER(LONG), INTENT(IN)       :: WHAT              ! Which block of code to write

! **********************************************************************************************************************************
      IF      (WHAT == 13) THEN

         WRITE(BUG,91301)

         DO I=1,3                                         ! Write final result, SHELL_A
            WRITE(BUG,91302) I, (SHELL_A(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,3                                         ! Write final result, SHELL_B
            WRITE(BUG,91303) I, (SHELL_B(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,3                                         ! Write final result, SHELL_D
            WRITE(BUG,91304) I, (SHELL_D(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,2                                         ! Write final result, SHELL_T
            WRITE(BUG,91305) I, (SHELL_T(I,J),J=1,2)
         ENDDO
         WRITE(BUG,*)

      ELSE IF (WHAT == 20) THEN

         PCOMP_PLIES = PCOMP(INTL_PID,5)
         WRITE(BUG,92001) EID, PCOMP(INTL_PID,1), PCOMP_PLIES, PCOMP_PLIES, PCOMP(INTL_PID,1)

         WRITE(BUG,92003) 'ZO        = ', RPCOMP(INTL_PID,1),'  z coord of bottom of ply 1'
         WRITE(BUG,92003) 'NSM       = ', RPCOMP(INTL_PID,2),'  nonstructural mass'
         WRITE(BUG,92003) 'SB        = ', RPCOMP(INTL_PID,3),'  interlaminar allowable shear stress'

         IF      (FAILURE_THEORY == 'NONE') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  no failure theory specified'
         ELSE IF (FAILURE_THEORY == 'HILL') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  use HILL failure theory'
         ELSE IF (FAILURE_THEORY == 'HOFF') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  use HOFF failure theory'
         ELSE IF (FAILURE_THEORY == 'TSAI') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  use TSAI failure theory'
         ELSE IF (FAILURE_THEORY == 'STRN') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  use STRN (max strain) failure theory'
         ELSE IF (FAILURE_THEORY == 'STRE') THEN
            WRITE(BUG,92002) 'FT        = ', FT, '  use STRE (max stress) failure theory'
         ENDIF

         WRITE(BUG,92003) 'TREF      = ', RPCOMP(INTL_PID,4),'  reference temperature'
         WRITE(BUG,92003) 'GE        = ', RPCOMP(INTL_PID,5),'  damping coefficient'

         IF      (PCOMP(INTL_PID,4) == 0) THEN
            LAM(1:3) = 'NON'
            WRITE(BUG,92002) 'LAM       = ', PCOMP(INTL_PID,4),'  laminate is not a symmetric layup'
         ELSE IF (PCOMP(INTL_PID,4) == 1) THEN
            LAM(1:3) = 'SYM'
            WRITE(BUG,92002) 'LAM       = ', PCOMP(INTL_PID,4),'  laminate is a symmetric layup'
         ENDIF

         WRITE(BUG,*)
         WRITE(BUG,*)

         IF (LAM(1:3) == 'NON') THEN

            WRITE(BUG,92004)
            WRITE(BUG,92006)
            DO J=1,PCOMP_PLIES
               IICOL = MPCOMP0  + MPCOMP_PLIES*(J - 1)  + 1
               IIROW = PCOMP(INTL_PID,IICOL)
               J1 = MPCOMP0 + MPCOMP_PLIES*(J-1) + 1
               J2 = MRPCOMP0 + MRPCOMP_PLIES*(J-1) + 1
               IF (PCOMP(INTL_PID,J1+1) == 0) THEN
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (NO )', RPCOMP(INTL_PID,J2+2)
               ELSE
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (YES)', RPCOMP(INTL_PID,J2+2)
               ENDIF
            ENDDO
            WRITE(BUG,92007) RPCOMP(INTL_PID,6)
            WRITE(BUG,*)

         ELSE

            WRITE(BUG,92005)
            WRITE(BUG,92006)
            DO J=1,PCOMP_PLIES/2
               IICOL = MPCOMP0  + MPCOMP_PLIES*(J - 1)  + 1
               IIROW = PCOMP(INTL_PID,IICOL)
               J1 = MPCOMP0 + MPCOMP_PLIES*(J-1) + 1
               J2 = MRPCOMP0 + MRPCOMP_PLIES*(J-1) + 1
               IF (PCOMP(INTL_PID,J1+1) == 0) THEN
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (NO )', RPCOMP(INTL_PID,J2+2)
               ELSE
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (YES)', RPCOMP(INTL_PID,J2+2)
               ENDIF
            ENDDO
            WRITE(BUG,*)

            DO J=PCOMP_PLIES/2 + 1,PCOMP_PLIES
               IICOL = MPCOMP0  + MPCOMP_PLIES*(J - 1)  + 1
               IIROW = PCOMP(INTL_PID,IICOL)
               J1 = MPCOMP0 + MPCOMP_PLIES*(J-1) + 1
               J2 = MRPCOMP0 + MRPCOMP_PLIES*(J-1) + 1
               IF (PCOMP(INTL_PID,J1+1) == 0) THEN
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (NO )', RPCOMP(INTL_PID,J2+2)
               ELSE
                  WRITE(BUG,92008) J, MATL(IIROW,1), RPCOMP(INTL_PID,J2), RPCOMP(INTL_PID,J2+1), '0 (YES)', RPCOMP(INTL_PID,J2+2)
               ENDIF
            ENDDO
            WRITE(BUG,92007) RPCOMP(INTL_PID,6)
            WRITE(BUG,*)

         ENDIF

      ELSE IF (WHAT == 21) THEN

         MTRL_NAME = '****'
         IF      (MTRL_TYPE(1) == 1) THEN
            MTRL_NAME = 'MAT1'
         ELSE IF (MTRL_TYPE(1) == 2) THEN
            MTRL_NAME = 'MAT2'
         ELSE IF (MTRL_TYPE(1) == 8) THEN
            MTRL_NAME = 'MAT8'
         ENDIF

         WRITE(BUG,92101) JPLY, TPLY, THETA_PLY, MTRL_NAME, MTRL_ACT_ID(1), INTL_MID(1)
         WRITE(BUG,92102)
         DO I=1,2                                          ! Write EM, ET ORIG material matrices
            WRITE(BUG,92103) I, (EM(I,J),J=1,3), I, (ET(I,J),J=1,2)
         ENDDO
         DO I=3,3
            WRITE(BUG,92104) I, (EM(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,92105)                                  ! Write stress/strain allowables
         WRITE(BUG,92106) ULT_STRE(1,1), ULT_STRE(2,1), ULT_STRE(7,1), ULT_STRE(8,3), ULT_STRE(9,3)
         WRITE(BUG,92107) ULT_STRE(3,1), ULT_STRE(4,1), ULT_STRE(7,1)
         WRITE(BUG,*)
         WRITE(BUG,92108) ULT_STRN(1,1), ULT_STRN(2,1), ULT_STRN(7,1), ULT_STRN(8,3), ULT_STRN(9,3)
         WRITE(BUG,92109) ULT_STRN(3,1), ULT_STRN(4,1), ULT_STRN(7,1)
         WRITE(BUG,*)

      ELSE IF (WHAT == 22) THEN

         WRITE(BUG,92201)
         DO I=1,2                                          ! Write T1M and T1T
            WRITE(BUG,92202) I, (T1M(I,J),J=1,3), I, (T1T(I,J),J=1,2)
         ENDDO
         DO I=3,3
            WRITE(BUG,92203) I, (T1M(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,92204)
         DO I=1,2                                          ! Write transformed material matrices, EM, ET
            WRITE(BUG,92205) I, (EM(I,J),J=1,3), I, (ET(I,J),J=1,2)
         ENDDO
         DO I=3,3
            WRITE(BUG,92206) I, (EM(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,92207) JPLY, ZBK, ZTK
         DO I=1,3                                         ! Write ply A matrix
            WRITE(BUG,92208) I, (PLY_A(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,3                                          ! Write ply B matrix
            WRITE(BUG,92209) I, (PLY_B(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,3                                          ! Write ply D matrix
            WRITE(BUG,92210) I, (PLY_D(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         DO I=1,2                                          ! Write ply T matrix
            WRITE(BUG,92211) I, (PLY_T(I,J),J=1,2)
         ENDDO
         WRITE(BUG,*)

      ELSE IF (WHAT == 23) THEN

         PCOMP_PLIES = PCOMP(INTL_PID,5)
         WRITE(BUG,92301) PCOMP(INTL_PID,1), PCOMP_PLIES

         DO I=1,3                                         ! Write final result, SHELL_A
            WRITE(BUG,92302) I, (SHELL_A(I,J),J=1,3)
         ENDDO
         WRITE(BUG,99991)

         DO I=1,3                                         ! Write final result, SHELL_B
            WRITE(BUG,92303) I, (SHELL_B(I,J),J=1,3)
         ENDDO
         WRITE(BUG,99991)

         DO I=1,3                                         ! Write final result, SHELL_D
            WRITE(BUG,92304) I, (SHELL_D(I,J),J=1,3)
         ENDDO
         WRITE(BUG,99991)

         DO I=1,2                                         ! Write final result, SHELL_T
            WRITE(BUG,92305) I, (SHELL_T(I,J),J=1,2)
         ENDDO
         WRITE(BUG,99991)
         WRITE(BUG,99992)

      ENDIF

! **********************************************************************************************************************************
91301 FORMAT(/,'   Shell A, B, D and T matrices',/,                                                                                &
               '   ----------------------------')
91302 FORMAT('   SHELL_A row',I2,' = ',3(1ES14.6),'  membrane             ')

91303 FORMAT('   SHELL_B row',I2,' = ',3(1ES14.6),'  bend/mem coupling    ')

91304 FORMAT('   SHELL_D row',I2,' = ',3(1ES14.6),'  bending              ')

91305 FORMAT('   SHELL_T row',I2,' = ',2(1ES14.6),'                transverse shear     ')

92001 FORMAT('   Element number',I8,' uses PCOMP ',I8,' with ',I3,' plies:',//,                                                    &
             27X,'D A T A   A P P L I C A B L E   T O   A L L  ',I3,'   P L I E S   F O R   P C O M P  ',I8,/)

92002 FORMAT(39X,A,I3,11X,A)

92003 FORMAT(39X,A,1ES14.6,A)

92004 FORMAT('                                                   I N D I V I D U A L   P L Y   D A T A',/)

92005 FORMAT('                                                   I N D I V I D U A L   P L Y   D A T A',/,                         &
             '                                    (* indicates symmetric ply implicitly defined by the PCOMP entry)',/)

92006 format(29x,'Ply no.         Actual      Thickness     Matl angle      SOUTi           ZI      ',/,                           &
             29x,'               Matl  ID       TPLY          THETAi                  (Coord of mid)',/,                           &
             29x,'-------      ------------  ------------  ------------  ------------  ------------')

92007 FORMAT(56x,'------------',/,54X,1ES14.6,'  total thickness')

92008 FORMAT(30X,I3,' ',9X,I8,3X,2(1ES14.6),4X,A7,3X,1ES14.6)

92101 FORMAT(//,5X,'  Ply number ',I2,' with thickness TPLY =',1ES13.6,' and THETA = ',0PF9.3,' deg uses ',A,1X,I8,                &
                     ' (internal matl ID',I7,')',/,                                                                                &
                5X,'  =======================================================================================================',    &
                '================',/)

92102 FORMAT(43X,'Membrane and transverse shear material matrices in input ply coords')

92103 FORMAT('             EM row',I2,' = ',3(1ES14.6),'                      ET row',I2,' = ',2(1ES14.6))

92104 FORMAT('             EM row',I2,' = ',3(1ES14.6))

92105 FORMAT(33x,'Membrane and transverse shear material stress and strain allowables in input ply coords'                      ,/,&
             29x,'Tension     Compression In-plane shear                                    Transv 13     Transv 23')

92106 FORMAT('        Stress (long) = ',3(1ES14.6),7X,'        Stress (transv) = ',2(1ES14.6))

92107 FORMAT('        Stress (lat ) = ',3(1ES14.6))

92108 FORMAT('        Strain (long) = ',3(1ES14.6),7X,'        Strain (transv) = ',2(1ES14.6))

92109 FORMAT('        Strain (lat ) = ',3(1ES14.6))

92201 FORMAT(45X,'Membrane and transverse shear material transformation matrices')

92202 FORMAT('            TEM row',I2,' = ',3(1ES14.6),'                     TET row',I2,' = ',2(1ES14.6))

92203 FORMAT('            TEM row',I2,' = ',3(1ES14.6))

92204 FORMAT(35X,'Membrane and transverse shear material matrices transformed to elem material coords')

92205 FORMAT('             EM row',I2,' = ',3(1ES14.6),'                      ET row',I2,' = ',2(1ES14.6))

92206 FORMAT('             EM row',I2,' = ',3(1ES14.6))

92207 FORMAT(52X,'A, B, D, and T matrices for ply number',I3,/,44X,'(coords of bot/top of ply = ',2(1ES14.6),')',/)

92208 FORMAT(41x,'Ply A row',I2,' = ',3(1ES14.6),'  membrane')

92209 FORMAT(41X,'Ply B row',I2,' = ',3(1ES14.6),'  bend/mem coupling')

92210 FORMAT(41X,'Ply D row',I2,' = ',3(1ES14.6),'  bending')

92211 FORMAT(41X,'Ply T row',I2,' = ',2(1ES14.6),'                transverse shear')

92301 FORMAT(//,36X,'+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +'                     ,/,&
                36X,'+                                                                                   +'                     ,/,&
                36X,'+  C O M P O S I T E   A, B, D, T   M A T R I C E S   F O R   P C O M P  ',I8,'   +'                       ,/,&
                36X,'+                             (sum over all ',I3,' plies)                              +'                  ,/,&
                36X,'+                                                                                   +')

92302 FORMAT(36X,'+   SHELL_A row',I2,' = ',3(1ES14.6),'  membrane            +')

92303 FORMAT(36X,'+   SHELL_B row',I2,' = ',3(1ES14.6),'  bend/mem coupling   +')

92304 FORMAT(36X,'+   SHELL_D row',I2,' = ',3(1ES14.6),'  bending             +')

92305 FORMAT(36X,'+   SHELL_T row',I2,' = ',2(1ES14.6),'                transverse shear    +')

99991 FORMAT(36X,'+                                                                                   +')

99992 FORMAT(36X,'+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +'/)

! **********************************************************************************************************************************

      END SUBROUTINE BUG_SHELL_ABD_MATRICES

      END SUBROUTINE SHELL_ABD_MATRICES


      SUBROUTINE IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )

! Given a shell (TRIA3 or QUAD4) element's internal ID, determine if its properties are defined on a Bulk Data PCOMP entry

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  DEDAT_T3_SHELL_KEY, DEDAT_Q4_SHELL_KEY, DEDAT_Q8_SHELL_KEY
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, ETYPE, PCOMP_PROPS, TYPE

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID for which
      INTEGER(LONG)                   :: EPNTK              ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                             row number in array EDAT where data begins for this element.
! **********************************************************************************************************************************
      EPNTK = EPNT(INT_ELEM_ID)
      TYPE  = ETYPE(INT_ELEM_ID)

      PCOMP_PROPS = 'N'
      IF      (TYPE(1:5) == 'TRIA3') THEN
         IF (EDAT(EPNTK+DEDAT_T3_SHELL_KEY) == 2) THEN
            PCOMP_PROPS = 'Y'
         ENDIF
      ELSE IF (TYPE(1:5) == 'QUAD4') THEN
         IF (EDAT(EPNTK+DEDAT_Q4_SHELL_KEY) == 2) THEN
            PCOMP_PROPS = 'Y'
         ENDIF
      ELSE IF (TYPE(1:5) == 'QUAD8') THEN
         IF (EDAT(EPNTK+DEDAT_Q8_SHELL_KEY) == 2) THEN
            PCOMP_PROPS = 'Y'
         ENDIF
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE IS_ELEM_PCOMP_PROPS


      SUBROUTINE GET_PCOMP_SECT_PROPS ( PCOMP_TM, PCOMP_IB, PCOMP_TS )

! Calculates section properties, PCOMP_TM, PCOMP_IB, PCOMP_TS, for shell elements that have PCOMP properties.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MPCOMP_PLIES, MPCOMP0, MRPCOMP_PLIES, MRPCOMP0
      USE MODEL_STUF, ONLY            :  EPROP, INTL_PID, NUM_PLIES, RPCOMP, TPLY
      USE PARAMS, ONLY                :  PCMPTSTM
      USE CONSTANTS_1, ONLY           :  ZERO, THIRD
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_PCOMP_SECT_PROPS'

      INTEGER(LONG)                   :: K                  ! DO loop index
      INTEGER(LONG)                   :: PLY_RPCOMP_INDEX   ! Index in array RPCOMP where data for ply K begins


      REAL(DOUBLE), INTENT(OUT)       :: PCOMP_TM           ! Membrane thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE), INTENT(OUT)       :: PCOMP_IB           ! Bending MOI of PCOMP for equivalent PSHELL
      REAL(DOUBLE), INTENT(OUT)       :: PCOMP_TS           ! Transverse shear thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE)                    :: ZBK,ZTK            ! Coord from ref plane to bot and top of ply K
      REAL(DOUBLE)                    :: ZBK2,ZTK2          ! ZBK^2, ZTK^2
      REAL(DOUBLE)                    :: ZBK3,ZTK3          ! ZBK^3, ZTK^3


! **********************************************************************************************************************************
      ZBK      = EPROP(4)
      PCOMP_TM = ZERO
      PCOMP_IB = ZERO
      PCOMP_TS = ZERO

      DO K=1,NUM_PLIES
                                                        ! Indices in PCOMP, RPCOMP arrays where ply K data begins
         PLY_RPCOMP_INDEX = MRPCOMP0 + MRPCOMP_PLIES*(K - 1) + 1

         TPLY = RPCOMP(INTL_PID,PLY_RPCOMP_INDEX)       ! Ply thickness

         ZTK = ZBK + TPLY    ;    ZTK2 = ZTK*ZTK    ;    ZTK3 = ZTK2*ZTK    ;    ZBK2 = ZBK*ZBK    ;    ZBK3 = ZBK2*ZBK

         PCOMP_TM = PCOMP_TM + (ZTK  - ZBK )
         PCOMP_IB = PCOMP_IB + (ZTK3 - ZBK3)*THIRD
         PCOMP_TS = PCOMP_TS + (ZTK  - ZBK )*PCMPTSTM   ! PCMPTSTM is a factor (< 1) for shear to total plate thickness for PCOMP

         ZBK = ZTK

      ENDDO




      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_PCOMP_SECT_PROPS


      SUBROUTINE ROT_COMP_ELEM_AXES ( INT_ELEM_ID, IPLY, THETA, DIRECTION )

! Rotates axes of a ply material and CTE matrices in ply coords to coords along and perpendicular to element material axes, or
! vice versa, depending on input arg DIRECTION

! Ref 1: "Stress Tensor coord transform.doc" my derivation documented in a WORD file located in \MYSTRAN\Documentation (this would
! be T1P in the code below)

! Ref 2: "Practical Analysis of Composites" by J.N. Reddy and A. Miravete, CRC Press, 1995 section 3.3. Eqn 7 in Ref 2 is matrix T1
! whereas my matrix is T1P (in the code below). Note that Ref 2 stress tensor definition and mine have a different order:

!      -------------------------------------------------------------------------
!     | Ref 2 stress tensor for matrix T1 | MYSTRAN stress tensor for matrix T1P|
!     |-----------------------------------|-------------------------------------|
!     |         sig-1 = sig-xx            |            sig-1 = sig-xx           |
!     |         sig-2 = sig-yy            |            sig-2 = sig-yy           |
!     |         sig-3 = sig-zz            |            sig-3 = sig-zz           |
!     |         sig-4 = sig-yz            |            sig-4 = sig-xy           |
!     |         sig-5 = sig-zx            |            sig-5 = sig-yz           |
!     |         sig-6 = sig-xy            |            sig-6 = sig-zx           |
!      -------------------------------------------------------------------------

! Thus, rows and cols 4 and 6 are transposed between Ref 2 matrix T1 and the code here for matrix T1P

! Terms from T1 transform the stresses (in vector, not tensor format) and material matrix and terms from T1' transform strains (in
! vector, not tensor format)

! Note that we can get these terms from the T1 matrix in Ref 2 (rather than the procedure outlined in subr ROT_AXES_MATL_TO_LOC)
! since the transformation here involves only one coord transformation, not 2 in subr ROT_AXES_MATL_TO_LOC (which rotates from
! material to local via basic system since we can't go directly from material to local in that case)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEMATC, DEDAT_Q4_MATANG_KEY, DEDAT_T3_MATANG_KEY
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  CONV_DEG_RAD, ZERO, HALF, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  ALPVEC, EB, EM, ET, EBM, INTL_MID, MTRL_TYPE, STRESS, STRAIN, T1P, T1M, T1T, T2P, T2M, &
                                         T2T, QUAD_DELTA, THETAM, TYPE, EDAT, MATANGLE, EPNT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ROT_COMP_ELEM_AXES'
      CHARACTER(LEN=*), INTENT(IN)    :: DIRECTION         ! =1-2, rotate from ply to elem mat'l axes (when gen ABD matrices)
!                                                            =2-1, rotate stress from elem mat'l to ply axes (when recov stresses)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG), INTENT(IN)       :: IPLY              ! Ply number
      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE), INTENT(IN)        :: THETA             ! Orient angle of long dir of ply i wrt matl axis for the composite elem
      REAL(DOUBLE)                    :: ALP3(3,MEMATC)    ! The 3 rows of ALPVEC for membrane strains
      REAL(DOUBLE)                    :: C,S,SC            ! COS, SIN, SIN*COS of RADIANS_ROT
      REAL(DOUBLE)                    :: C2,S2             ! C*C, S*S
      REAL(DOUBLE)                    :: RADIANS_ROT       ! Angle (radians) to rotate from elem material axes to ply axes
      REAL(DOUBLE)                    :: DUM1(3)           ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(3)           ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM3(3)           ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM33(3,3)        ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM22(2,2)        ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM3M(3,MEMATC)   ! Intermediate matrix
      REAL(DOUBLE)                    :: STRESS1(3)        ! Rows 1-3 of STRESS
      REAL(DOUBLE)                    :: STRESS2(3)        ! Rows 4-6 of STRESS
      REAL(DOUBLE)                    :: STRESS3(3)        ! Rows 7-9 of STRESS
      REAL(DOUBLE)                    :: STRAIN1(3)        ! Rows 1-3 of STRAIN
      REAL(DOUBLE)                    :: STRAIN2(3)        ! Rows 4-6 of STRAIN
      REAL(DOUBLE)                    :: STRAIN3(3)        ! Rows 7-9 of STRAIN
      REAL(DOUBLE)                    :: T1Mt(3,3)         ! Transpose of T1M
      REAL(DOUBLE)                    :: T1Tt(2,2)         ! Transpose of T1T
      REAL(DOUBLE)                    :: T1T3(3,3)         ! T1T expanded to 3x3
      REAL(DOUBLE)                    :: T2T3(3,3)         ! T2T expanded to 3x3
      REAL(DOUBLE)                    :: MATL_AXES_ROTATE  ! Angle in radians to rotate material axes to coincide with local elem x axis
      INTEGER(LONG)                   :: INT41,INT42       ! An integer used in getting MATANGLE
      CHARACTER( 2*BYTE)              :: LOC               ! Location where THETAM is calculated (for DEBUG output purposes)
      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID.




! **********************************************************************************************************************************
! Calc T1P matrix from eqn 3.3-7 in Ref 1. (with order 1,2,3,4,5,6 changed to 1,2,3,6,4,5 to account for the fact that Ref (1) has
! the 6 position for xy stress whereas it is the 4th position here)

      RADIANS_ROT = CONV_DEG_RAD*THETA                     ! THETA is angle (deg) from elem matl axis to ply K longitudinal axis


! **********************************************************************************************************************************
! Add material axes angle to ply angle.


      EPNTK = EPNT(INT_ELEM_ID)

      !----- Copy and pasted from EMG.f03 -----
      THETAM = ZERO

      IF      (TYPE(1:5) == 'QUAD4') THEN
         INT41 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY)     ! Key to say whether matl angle is an actual angle or a coord ID
         INT42 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY+1)   ! Key to say (if INT41 is neg) whether coord ID is basic or otherwise
         IF      (INT41 >  0) THEN                   ! Angle is defined in array MATANGLE at row INT41
            LOC = '#1'
            THETAM = CONV_DEG_RAD*MATANGLE( INT41 )
         ELSE IF (INT41 <  0) THEN                   ! Angle is defined by a coord sys ID whose value is -INT41
            LOC = '#2'
            CALL GET_MATANGLE_FROM_CID ( -INT41 )
         ELSE IF (INT41 ==  0) THEN                  ! Angle is either specified as defined by basic coord sys or angle is 0.
            IF      (INT42 == 1) THEN
               LOC = '#3'
               CALL GET_MATANGLE_FROM_CID ( 0 )
            ELSE IF (INT42 == 0) THEN
               LOC = '#4'
               THETAM = ZERO
            ENDIF
         ENDIF

      ELSE IF (TYPE(1:5) == 'TRIA3') THEN
         INT41 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY)     ! Key to say whether matl angle is an actual angle or a coord ID
         INT42 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY+1)   ! Key to say (if INT41 is neg) whether coord ID is basic or otherwise
         IF      (INT41 >  0) THEN                   ! Angle is defined in array MATANGLE at row INT41
            LOC = '#1'
            THETAM = CONV_DEG_RAD*MATANGLE( INT41 )
         ELSE IF (INT41 <  0) THEN                   ! Angle is defined by a coord sys ID whose value is -INT41
            LOC = '#2'
            CALL GET_MATANGLE_FROM_CID ( -INT41 )
         ELSE IF (INT41 ==  0) THEN                  ! Angle is either specified as defined by basic coord sys or angle is 0.
            IF      (INT42 == 1) THEN
               LOC = '#3'
               CALL GET_MATANGLE_FROM_CID ( 0 )
            ELSE IF (INT42 == 0) THEN
               LOC = '#4'
               THETAM = ZERO
            ENDIF
         ENDIF

      ENDIF
      !----- Copy and pasted from EMG.f03 -----

      IF      (TYPE(1:5) == 'QUAD4') THEN

         MATL_AXES_ROTATE = THETAM - QUAD_DELTA         ! We need angle which would rotate local elem x axis to mat'l axis,
!                                                            NOT the opposite (since we want a transform matrix that will take
!                                                            a vector in elem coords and convert it into a vector in mat'l coords
      ELSE IF (TYPE(1:5) == 'TRIA3') THEN

         MATL_AXES_ROTATE = THETAM

      ENDIF

      RADIANS_ROT = RADIANS_ROT + MATL_AXES_ROTATE

! **********************************************************************************************************************************




      C  = DCOS(RADIANS_ROT)
      S  = DSIN(RADIANS_ROT)
      SC = S*C
      C2 = C*C
      S2 = S*S

! Transformation matrix T1P is 6x6 for transforming 3D solid material matrices from ply coords to element coords. Notice that this
! is the same as matrix T1 in Ref 2 except that the ordering of rows and cols is different since my definition of the stress tensor
! has a different ordering than Ref 2

      T1P(1,1) =  C2    ;  T1P(1,2) =  S2    ;  T1P(1,3) =  ZERO  ;  T1P(1,4) = -TWO*SC   ;  T1P(1,5) =  ZERO  ;  T1P(1,6) =  ZERO

      T1P(2,1) =  S2    ;  T1P(2,2) =  C2    ;  T1P(2,3) =  ZERO  ;  T1P(2,4) =  TWO*S*C  ;  T1P(2,5) =  ZERO  ;  T1P(2,6) =  ZERO

      T1P(3,1) =  ZERO  ;  T1P(3,2) =  ZERO  ;  T1P(3,3) =  ONE   ;  T1P(3,4) =  ZERO     ;  T1P(3,5) =  ZERO  ;  T1P(3,6) =  ZERO

      T1P(4,1) =  SC    ;  T1P(4,2) = -SC    ;  T1P(4,3) =  ZERO  ;  T1P(4,4) =  C2 - S2  ;  T1P(4,5) =  ZERO  ;  T1P(4,6) =  ZERO

      T1P(5,1) =  ZERO  ;  T1P(5,2) =  ZERO  ;  T1P(5,3) =  ZERO  ;  T1P(5,4) =  ZERO     ;  T1P(5,5) =  C     ;  T1P(5,6) =  S

      T1P(6,1) =  ZERO  ;  T1P(6,2) =  ZERO  ;  T1P(6,3) =  ZERO  ;  T1P(6,4) =  ZERO     ;  T1P(6,5) = -S     ;  T1P(6,6) =  C

! T1M and T1T are portions of T1P for membrane and transverse shear (for 2D shell elements)

      T1M(1,1) = T1P(1,1)   ;   T1M(1,2) = T1P(1,2)   ;   T1M(1,3) = T1P(1,4)
      T1M(2,1) = T1P(2,1)   ;   T1M(2,2) = T1P(2,2)   ;   T1M(2,3) = T1P(2,4)
      T1M(3,1) = T1P(4,1)   ;   T1M(3,2) = T1P(4,2)   ;   T1M(3,3) = T1P(4,4)

      T1T(1,1) = T1P(5,5)   ;   T1T(1,2) = T1P(6,5)
      T1T(2,1) = T1P(5,6)   ;   T1T(2,2) = T1P(6,6)

! Transformation matrix T2P is the inverse of T1P.

      T2P(1,1) =  C2    ;  T2P(1,2) =  S2    ;  T2P(1,3) =  ZERO  ;  T2P(1,4) =  TWO*SC   ;  T2P(1,5) =  ZERO  ;  T2P(1,6) =  ZERO

      T2P(2,1) =  S2    ;  T2P(2,2) =  C2    ;  T2P(2,3) =  ZERO  ;  T2P(2,4) = -TWO*S*C  ;  T2P(2,5) =  ZERO  ;  T2P(2,6) =  ZERO

      T2P(3,1) =  ZERO  ;  T2P(3,2) =  ZERO  ;  T2P(3,3) =  ONE   ;  T2P(3,4) =  ZERO     ;  T2P(3,5) =  ZERO  ;  T2P(3,6) =  ZERO

      T2P(4,1) = -SC    ;  T2P(4,2) =  SC    ;  T2P(4,3) =  ZERO  ;  T2P(4,4) =  C2 - S2  ;  T2P(4,5) =  ZERO  ;  T2P(4,6) =  ZERO

      T2P(5,1) =  ZERO  ;  T2P(5,2) =  ZERO  ;  T2P(5,3) =  ZERO  ;  T2P(5,4) =  ZERO     ;  T2P(5,5) =  C     ;  T2P(5,6) = -S

      T2P(6,1) =  ZERO  ;  T2P(6,2) =  ZERO  ;  T2P(6,3) =  ZERO  ;  T2P(6,4) =  ZERO     ;  T2P(6,5) =  S     ;  T2P(6,6) =  C

!  T2M and T2T are portions of T2P for membrane and transverse shear (for 2D shell elements)

      T2M(1,1) = T2P(1,1)   ;   T2M(1,2) = T2P(1,2)   ;   T2M(1,3) = T2P(1,4)
      T2M(2,1) = T2P(2,1)   ;   T2M(2,2) = T2P(2,2)   ;   T2M(2,3) = T2P(2,4)
      T2M(3,1) = T2P(4,1)   ;   T2M(3,2) = T2P(4,2)   ;   T2M(3,3) = T2P(4,4)

      T2T(1,1) = T2P(5,5)   ;   T2T(1,2) = T2P(5,6)
      T2T(2,1) = T2P(6,5)   ;   T2T(2,2) = T2P(6,6)

! Do transformation

      IF      (DIRECTION == '1-2') THEN                    ! Transform material and CTE matrices from ply to elem material axes
                                                           ! (1) Membrane matl matrix
         DO I=1,3
            DO J=1,3
               T1Mt(I,J) = T1M(J,I)
            ENDDO
         ENDDO

         DO I=1,2
            DO J=1,2
               T1Tt(I,J) = T1T(J,I)
            ENDDO
         ENDDO

         DO J=1,MEMATC                                     ! Transform ALPVEC rows for ply membrane from ply to elem axes
            ALP3(1,J) = ALPVEC(1,J)
            ALP3(2,J) = ALPVEC(2,J)
            ALP3(3,J) = ALPVEC(4,J)
         ENDDO

         CALL MATMULT_FFF ( T1M   , EM  , 3, 3, 3, DUM33 ) ! (1) Transform EM  membrane matl matrix from ply to elem coords
         CALL MATMULT_FFF ( DUM33 , T1Mt, 3, 3, 3, EM)

         CALL MATMULT_FFF ( T1M   , EB  , 3, 3, 3, DUM33 ) ! (2) Transform EB  bending matl matrix from ply to elem coords
         CALL MATMULT_FFF ( DUM33 , T1Mt, 3, 3, 3, EB)

         CALL MATMULT_FFF ( T1T   , ET  , 2, 2, 2, DUM22 ) ! (3) Transform TM  transverse shear matl matrix from ply to elem coords
         CALL MATMULT_FFF ( DUM22 , T1Tt, 2, 2, 2, ET)

         CALL MATMULT_FFF ( T1M   , EBM , 3, 3, 3, DUM33 ) ! (4) Transform EBM bending/membrane matl matrix from ply to elem coords
         CALL MATMULT_FFF ( DUM33 , T1Mt, 3, 3, 3, EBM)

         CALL MATMULT_FFF_T ( T2M, ALP3, 3, 3, MEMATC, DUM3M )

         DO J=1,MEMATC
            ALPVEC(1,J) = DUM3M(1,J)
            ALPVEC(2,J) = DUM3M(2,J)
            ALPVEC(3,J) = ZERO
            ALPVEC(4,J) = DUM3M(3,J)
            ALPVEC(5,J) = ZERO
            ALPVEC(6,J) = ZERO
         ENDDO


      ELSE IF (DIRECTION == '2-1') THEN                    ! Transform stresses in elem material axes to ply axes

         DO I=1,3                                          ! Initialize T2T3, T2T3
            DO J=1,3
               T1T3(I,J) = ZERO
               T2T3(I,J) = ZERO
            ENDDO
         ENDDO

         DO I=1,2                                          ! Load all 2x2 of T1T into upper 2x2 of T1T3
            DO J=1,2
!zzzz          T1T3  = T1T(I,J)
               T1T3(I,J)  = T1T(I,J)
            ENDDO
         ENDDO
         T1T3(3,3)  = ZERO

         DO I=1,2                                          ! Load all 2x2 of T2T into upper 2x2 of T2T3
            DO J=1,2
!zzzz          T2T3  = T2T(I,J)
               T2T3(I,J)  = T2T(I,J)
            ENDDO
         ENDDO
         T2T3(3,3)  = ZERO


         DO I=1,3                                          ! Transform axes on stress (in vector form) from elem to ply axes
            DUM1(I) = STRESS(I)
            DUM2(I) = STRESS(I+3)
            DUM3(I) = STRESS(I+6)
         ENDDO

         CALL MATMULT_FFF ( T2M , DUM1, 3, 3, 1, STRESS1 )
         CALL MATMULT_FFF ( T2M , DUM2, 3, 3, 1, STRESS2 )
         CALL MATMULT_FFF ( T2T3, DUM3, 3, 3, 1, STRESS3 )

         DO I=1,3
            STRESS(I)   = STRESS1(I)
            STRESS(I+3) = STRESS2(I)
            STRESS(I+6) = STRESS3(I)
         ENDDO

         DO I=1,3                                          ! Transform axes on strain (in vector form) from elem to ply axes
            DUM1(I) = STRAIN(I)
            DUM2(I) = STRAIN(I+3)
            DUM3(I) = STRAIN(I+6)
         ENDDO

         CALL MATMULT_FFF_T ( T1M , DUM1, 3, 3, 1, STRAIN1 )
         CALL MATMULT_FFF_T ( T1M , DUM2, 3, 3, 1, STRAIN2 )
         CALL MATMULT_FFF_T ( T1T3, DUM3, 3, 3, 1, STRAIN3 )

         DO I=1,3
            STRAIN(I)   = STRAIN1(I)
            STRAIN(I+3) = STRAIN2(I)
            STRAIN(I+6) = STRAIN3(I)
         ENDDO


      ELSE

         FATAL_ERR = FATAL_ERR
         WRITE(ERR,1962) SUBR_NAME, DIRECTION
         WRITE(F06,1962) SUBR_NAME, DIRECTION
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1962 FORMAT(' *ERROR  1962: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INPUT ARGUMENT "DIRECTION" MUST BE = "1-2" OR "2-1" BUT VALUE IS ',A)







! **********************************************************************************************************************************

      END SUBROUTINE ROT_COMP_ELEM_AXES



      SUBROUTINE SOLVE_SHELL_ALP ( SHELL_ALP_ERR )

! Solves for SHELL_ALP which is an effective CTE matrix for the PSHELL equivalent for a PCOMP (akin to ALPVEC's 1st 3 rows for a
! homogeneous shell element)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEMATC
      USE PARAMS, ONLY                :  EPSIL

      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  SHELL_ALP, SHELL_A, SHELL_B, SHELL_D, SHELL_T, SHELL_AALP, SHELL_BALP, SHELL_DALP,        &
                                         SHELL_TALP

      USE FULL_MATRIX_ALGEBRA, ONLY   :  INVERT_FF_MAT, MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_SHELL_ALP'

      INTEGER(LONG), INTENT(OUT)      :: SHELL_ALP_ERR(4)  ! Error indicator if SHELL_ALP was calculated
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: INFO              ! Error indicator from subr INVERT_FF_MAT

      REAL(DOUBLE)                    :: A(3,3)            ! One of the matrices SHELL_A, SHELL_D, SHELL_T, SHELL_B
      REAL(DOUBLE)                    :: AI(3,3)           ! Inverse of A
      REAL(DOUBLE)                    :: B(3)              ! Solution for a col of SHELL_ALP
      REAL(DOUBLE)                    :: C(3)              ! One of the matrices SHELL_AALP, SHELL_DALP, SHELL_TALP, SHELL_BALP
      REAL(DOUBLE)                    :: DETA              ! Determinant of A = SHELL_A
      REAL(DOUBLE)                    :: EPS1              ! Small number with which to comapre zero

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      DO I=1,3
         B(I) = ZERO
         DO J=1,MEMATC
            SHELL_ALP(I,J) = ZERO
         ENDDO
      ENDDO

      DO K=1,4                                             ! Solve for cols of SHELL_ALP. Put in array ALP_COL(I) for the time being

         IF      (K == 1) THEN

            DO I=1,3
               C(I) = SHELL_AALP(I)
               DO J=1,3
                  A(I,J) = SHELL_A(I,J)
               ENDDO
            ENDDO

         ELSE IF (K == 2) THEN

            DO I=1,3
               C(I) = SHELL_DALP(I)
               DO J=1,3
                  A(I,J) = SHELL_D(I,J)
               ENDDO
            ENDDO

         ELSE IF (K == 3) THEN

            DO I=1,3                                       ! Make A for SHELL_T into a 3x3 matrix with SHELL_T in the upper 2x2 and
               C(I) = ZERO                                 ! a 1 in the 3x3 slot. In this manner we can call the inversion routine
               DO J=1,3                                    ! in a single loop
                  A(I,J) = ZERO
               ENDDO
            ENDDO
            DO I=1,2
               C(I) = SHELL_TALP(I)
               DO J=1,2
                  A(I,J) = SHELL_T(I,J)
               ENDDO
            ENDDO
            A(3,3) = ONE

         ELSE IF (K == 4) THEN

            DO I=1,3
               C(I) = SHELL_BALP(I)
               DO J=1,3
                  A(I,J) = SHELL_B(I,J)
               ENDDO
            ENDDO

         ENDIF

         SHELL_ALP_ERR(K) = 0

         DETA = (A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(3,1)*A(2,3)*A(1,2))                                               &
              - (A(1,3)*A(2,2)*A(3,1) + A(1,2)*A(2,1)*A(3,3) + A(1,1)*A(3,2)*A(2,3))

         IF (DABS(DETA) > EPS1) THEN                       ! Det A = SHELL_A is > 0 so call subr to invert A = SHELL_A

            DO I=1,3
               DO J=1,3
                  AI(I,J) = A(I,J)
               ENDDO
            ENDDO

            CALL INVERT_FF_MAT ( SUBR_NAME, 'NAME(K)', AI, 3, INFO )

            IF (INFO == 0) THEN                            ! Inversion was successful so calc SHELL_ALP

               CALL MATMULT_FFF ( AI, C, 3, 3, 1, B )

               IF      (K == 1) THEN

                  DO I=1,3
                     SHELL_ALP(I,1) = B(I)
                  ENDDO

               ELSE IF (K == 2) THEN

                  DO I=1,3
                     SHELL_ALP(I,2) = B(I)
                  ENDDO

               ELSE IF (K == 3) THEN

                  DO I=1,2
                     SHELL_ALP(I+4,3) = B(I)
                  ENDDO

               ELSE IF (K == 4) THEN

                  DO I=1,3
                     SHELL_ALP(I,4) = B(I)
                  ENDDO

               ENDIF

            ELSE

               SHELL_ALP_ERR(K) = 2                        ! Inversion was unsuccessful so set error flag = 2

            ENDIF

         ELSE                                              ! Det A = 0 so set error flag = 1

            SHELL_ALP_ERR(K) = 1

         ENDIF

      ENDDO


      RETURN

! **********************************************************************************************************************************



! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_SHELL_ALP


      SUBROUTINE WRITE_PCOMP_EQUIV ( PCOMP_TM, PCOMP_IB, PCOMP_TS )

! Write equiv PSHELL and MAT2's for a PCOMP used, if requested, based on user Bulk Data PARAM PCOMPEQ

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  TWELVE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEMATC, MID1_PCOMP_EQ, MID2_PCOMP_EQ, MID3_PCOMP_EQ,                        &
                                         MID4_PCOMP_EQ, MID1_PCOMP_EQ, MID2_PCOMP_EQ, MID3_PCOMP_EQ, MID4_PCOMP_EQ
      USE PARAMS, ONLY                :  EPSIL, PCOMPEQ, SUPINFO
      USE MODEL_STUF, ONLY            :  INTL_PID, PCOMP, RHO, SHELL_ALP, SHELL_A, SHELL_B, SHELL_D, SHELL_T, SHELL_T_MOD,         &
                                         TREF, ZS

      USE TIMDAT, ONLY                :  TSEC

      USE TEXT_FIELD_UTILS, ONLY      :  REAL_DATA_TO_C8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_PCOMP_EQUIV'
      CHARACTER( 8*BYTE)              :: C8FLD_GIJ(3,3,4)   ! Char representation of MAT2 Gij entries
      CHARACTER( 8*BYTE)              :: C8FLD_ALP(6,MEMATC)! Char representation of MAT2 Ai (CTE) entries
      CHARACTER( 8*BYTE)              :: C8FLD_RHO(4)       ! Char representation of MAT2 RHO entry
      CHARACTER( 8*BYTE)              :: C8FLD_TREF(4)      ! Char representation of MAT2 TREF entry
      CHARACTER( 8*BYTE)              :: C8FLD_TM           ! Char representation of MAT2 TM entry
      CHARACTER( 8*BYTE)              :: C8FLD_ZS(2)        ! Char representation of MAT2 ZS entry
      CHARACTER(18*BYTE)              :: NAME1(4)           ! Name of a SHELL matrix (A, D, T, or B)
      CHARACTER( 7*BYTE)              :: NAME2(4)           ! Name of a SHELL matrix (A, D, T, or B)
      CHARACTER( 1*BYTE)              :: FINITE_MAT_PROPS(4)! Indicator of whether any SHELL matrix had zero diag terms

      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: ICONT              ! Continuation mnemonic for PSHELL equivalent B.D. entry
      INTEGER(LONG)                   :: IERR(4)            ! Error indicator if SHELL_ALP was calculated

      REAL(DOUBLE), INTENT(IN)        :: PCOMP_TM           ! Membrane thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE), INTENT(IN)        :: PCOMP_IB           ! Bending MOI of PCOMP for equivalent PSHELL
      REAL(DOUBLE), INTENT(IN)        :: PCOMP_TS           ! Transverse shear thickness of PCOMP for equivalent PSHELL
      REAL(DOUBLE)                    :: EPS1               ! Small number
      REAL(DOUBLE)                    :: PCOMP_IBP          ! 12*IB/TM^3
      REAL(DOUBLE)                    :: PCOMP_TSTM         ! TM/TS
      REAL(DOUBLE)                    :: PCOMP_BEN_MAT2(3,3)! MAT2 material matrix for bending for equivalent PSHELL
      REAL(DOUBLE)                    :: PCOMP_MBC_MAT2(3,3)! MAT2 material matrix for mem/ben coupling for equivalent PSHELL
      REAL(DOUBLE)                    :: PCOMP_MEM_MAT2(3,3)! MAT2 material matrix for membrane for equivalent PSHELL
      REAL(DOUBLE)                    :: PCOMP_TSH_MAT2(2,2)! MAT2 material matrix for transverse shear for equivalent PSHELL

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Determine if any SHELL_A, D, T, B  matrices are null

      FINITE_MAT_PROPS(1) = 'Y'
      IF ((DABS(SHELL_A(1,1)) < EPS1) .OR. (DABS(SHELL_A(1,1)) < EPS1) .OR. (DABS(SHELL_A(1,1)) < EPS1) .OR.                       &
          (DABS(SHELL_A(1,1)) < EPS1)) THEN
         FINITE_MAT_PROPS(1) = 'N'
      ENDIF

      FINITE_MAT_PROPS(2) = 'Y'
      IF ((DABS(SHELL_D(1,1)) < EPS1) .OR. (DABS(SHELL_D(1,1)) < EPS1) .OR. (DABS(SHELL_D(1,1)) < EPS1) .OR.                       &
          (DABS(SHELL_D(1,1)) < EPS1)) THEN
         FINITE_MAT_PROPS(2) = 'N'
      ENDIF

      FINITE_MAT_PROPS(3) = 'Y'
      IF ((DABS(SHELL_T(1,1)) < EPS1) .OR. (DABS(SHELL_T(1,1)) < EPS1) .OR. (DABS(SHELL_T(1,1)) < EPS1)) THEN
         FINITE_MAT_PROPS(3) = 'N'
      ENDIF

      FINITE_MAT_PROPS(4) = 'Y'
      IF ((DABS(SHELL_B(1,1)) < EPS1) .OR. (DABS(SHELL_B(1,1)) < EPS1) .OR. (DABS(SHELL_B(1,1)) < EPS1) .OR.                       &
          (DABS(SHELL_B(1,1)) < EPS1)) THEN
         FINITE_MAT_PROPS(4) = 'N'
      ENDIF

! Write message if any CTE's were not able to be calculated, but only if that SHELL matrix had props to be output

      NAME1(1) = 'MEMBRANE'  ;   NAME1(2) = 'BENDING'   ;   NAME1(3) = 'TRANSVERSE SHEAR'   ;   NAME1(4) = 'MEMB/BEND COUPLING'
      NAME2(1) = 'SHELL_A'   ;   NAME2(2) = 'SHELL_D'   ;   NAME2(3) = 'SHELL_T'            ;   NAME2(4) = 'SHELL_B'

      CALL SOLVE_SHELL_ALP ( IERR )                         ! Solve for equiv CTE's for this PCOMP

      WRITE(F06,9900)
      DO I=1,4
         IF      (IERR(I) == 1) THEN
            IF (FINITE_MAT_PROPS(I) == 'Y') THEN
               WRITE(ERR,9801) NAME1(I), PCOMP(INTL_PID,1), NAME2(I)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9801) NAME1(I), PCOMP(INTL_PID,1), NAME2(I)
               ENDIF
            ENDIF
         ELSE IF (IERR(I) == 2) THEN
            IF (FINITE_MAT_PROPS(I) == 'Y') THEN
               WRITE(ERR,9802) NAME1(I), PCOMP(INTL_PID,1), NAME2(I)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,9802) NAME1(I), PCOMP(INTL_PID,1), NAME2(I)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      PCOMP_IBP  = TWELVE*PCOMP_IB/(PCOMP_TM*PCOMP_TM*PCOMP_TM)
      PCOMP_TSTM = PCOMP_TS/PCOMP_TM
      WRITE(F06,9901) PCOMP(INTL_PID,1), PCOMP_TM, PCOMP_IBP, PCOMP_TSTM
      IF (SUPINFO == 'N') THEN
         WRITE(F06,9901) PCOMP(INTL_PID,1), PCOMP_TM, PCOMP_IBP, PCOMP_TSTM
      ENDIF

      DO I=1,3
         DO J =1,3
            PCOMP_MEM_MAT2(I,J) = SHELL_A(I,J)/PCOMP_TM
            PCOMP_MBC_MAT2(I,J) = SHELL_B(I,J)/-(PCOMP_TM * PCOMP_TM)
            PCOMP_BEN_MAT2(I,J) = SHELL_D(I,J)/PCOMP_IB
         ENDDO
      ENDDO

      DO I=1,2
         DO J =1,2
            PCOMP_TSH_MAT2(I,J) = SHELL_T(I,J)/PCOMP_TS
         ENDDO
      ENDDO

      MID1_PCOMP_EQ = MID1_PCOMP_EQ + PCOMP(INTL_PID,1)
      MID2_PCOMP_EQ = MID2_PCOMP_EQ + PCOMP(INTL_PID,1)
      MID3_PCOMP_EQ = MID3_PCOMP_EQ + PCOMP(INTL_PID,1)
      MID4_PCOMP_EQ = MID4_PCOMP_EQ + PCOMP(INTL_PID,1)

      ICONT = MID1_PCOMP_EQ

      IF (PCOMPEQ > 1) THEN

         IF ((IERR(1) == 0) .OR. (IERR(2) == 0) .OR. (IERR(3) == 0) .OR. (IERR(4) == 0)) THEN
            WRITE(F06,9902)
         ENDIF

         IF (FINITE_MAT_PROPS(1) == 'Y') THEN
            IF (IERR(1) == 0) THEN
               WRITE(F06,9903) MID1_PCOMP_EQ, PCOMP_MEM_MAT2(1,1), PCOMP_MEM_MAT2(1,2), PCOMP_MEM_MAT2(1,3),                       &
                                              PCOMP_MEM_MAT2(2,2), PCOMP_MEM_MAT2(2,3), PCOMP_MEM_MAT2(3,3),                       &
                                              SHELL_ALP(1,1), SHELL_ALP(2,1), SHELL_ALP(3,1)
            ELSE
               WRITE(F06,9903) MID1_PCOMP_EQ, PCOMP_MEM_MAT2(1,1), PCOMP_MEM_MAT2(1,2), PCOMP_MEM_MAT2(1,3),                       &
                                              PCOMP_MEM_MAT2(2,2), PCOMP_MEM_MAT2(2,3), PCOMP_MEM_MAT2(3,3)
            ENDIF
         ENDIF

         IF (FINITE_MAT_PROPS(2) == 'Y') THEN
            IF (IERR(2) == 0) THEN
               WRITE(F06,9904) MID2_PCOMP_EQ, PCOMP_BEN_MAT2(1,1), PCOMP_BEN_MAT2(1,2), PCOMP_BEN_MAT2(1,3),                       &
                                              PCOMP_BEN_MAT2(2,2), PCOMP_BEN_MAT2(2,3), PCOMP_BEN_MAT2(3,3),                       &
                                              SHELL_ALP(1,2), SHELL_ALP(2,2), SHELL_ALP(3,2)
            ELSE
               WRITE(F06,9904) MID2_PCOMP_EQ, PCOMP_BEN_MAT2(1,1), PCOMP_BEN_MAT2(1,2), PCOMP_BEN_MAT2(1,3),                       &
                                              PCOMP_BEN_MAT2(2,2), PCOMP_BEN_MAT2(2,3), PCOMP_BEN_MAT2(3,3)
            ENDIF
         ENDIF

         IF (FINITE_MAT_PROPS(3) == 'Y') THEN
            IF (IERR(3) == 0) THEN
               WRITE(F06,9905) MID3_PCOMP_EQ, PCOMP_TSH_MAT2(1,1), PCOMP_TSH_MAT2(1,2), PCOMP_TSH_MAT2(2,2),                       &
                                              SHELL_ALP(5,3), SHELL_ALP(6,3)
            ELSE
               WRITE(F06,9905) MID3_PCOMP_EQ, PCOMP_TSH_MAT2(1,1), PCOMP_TSH_MAT2(1,2), PCOMP_TSH_MAT2(2,2)
            ENDIF
         ENDIF

         IF (FINITE_MAT_PROPS(4) == 'Y') THEN
            IF (IERR(4) == 0) THEN
               WRITE(F06,9906) MID4_PCOMP_EQ, PCOMP_MBC_MAT2(1,1), PCOMP_MBC_MAT2(1,2), PCOMP_MBC_MAT2(1,3),                       &
                                              PCOMP_MBC_MAT2(2,2), PCOMP_MBC_MAT2(2,3), PCOMP_MBC_MAT2(3,3),                       &
                                              SHELL_ALP(1,4), SHELL_ALP(2,4), SHELL_ALP(3,4)
            ELSE
               WRITE(F06,9906) MID4_PCOMP_EQ, PCOMP_MBC_MAT2(1,1), PCOMP_MBC_MAT2(1,2), PCOMP_MBC_MAT2(1,3),                       &
                                              PCOMP_MBC_MAT2(2,2), PCOMP_MBC_MAT2(2,3), PCOMP_MBC_MAT2(3,3)
            ENDIF
         ENDIF

         WRITE(F06,*)

      ENDIF

      CALL GET_CHAR8_OUTPUTS ( C8FLD_GIJ, C8FLD_ALP, C8FLD_RHO, C8FLD_TREF, C8FLD_TM, C8FLD_ZS )

! Write PSHELL

      IF (IERR(4) == 0) THEN
         WRITE(F06,9912) PCOMP(INTL_PID,1), MID1_PCOMP_EQ, C8FLD_TM, MID2_PCOMP_EQ, MID3_PCOMP_EQ,ICONT,                           &
                         ICONT, C8FLD_ZS(2), C8FLD_ZS(2), MID4_PCOMP_EQ
      ELSE
         WRITE(F06,9913) PCOMP(INTL_PID,1), MID1_PCOMP_EQ, C8FLD_TM, MID2_PCOMP_EQ, MID3_PCOMP_EQ,ICONT,                           &
                            ICONT, C8FLD_ZS(1), C8FLD_ZS(2)
      ENDIF

! Write MAT2 for membrane

      IF (FINITE_MAT_PROPS(1) == 'Y') THEN
         IF (IERR(1) == 0) THEN

            WRITE(F06,9922) MID1_PCOMP_EQ, C8FLD_GIJ(1,1,1), C8FLD_GIJ(1,2,1), C8FLD_GIJ(1,3,1),                                   &
                                           C8FLD_GIJ(2,2,1), C8FLD_GIJ(2,3,1), C8FLD_GIJ(3,3,1), C8FLD_RHO(1),ICONT+1,             &
                                           ICONT+1, C8FLD_ALP(1,1), C8FLD_ALP(2,1), C8FLD_ALP(3,1), C8FLD_TREF(1)
         ELSE
            WRITE(F06,9923) MID1_PCOMP_EQ, C8FLD_GIJ(1,1,1), C8FLD_GIJ(1,2,1), C8FLD_GIJ(1,3,1),                                   &
                                           C8FLD_GIJ(2,2,1), C8FLD_GIJ(2,3,1), C8FLD_GIJ(3,3,1), C8FLD_RHO(1),ICONT+1,             &
                                           ICONT+1, C8FLD_TREF(1)
         ENDIF
      ENDIF

! Write MAT2 for bending

      IF (FINITE_MAT_PROPS(2) == 'Y') THEN
         IF (IERR(2) == 0) THEN
            WRITE(F06,9922) MID2_PCOMP_EQ, C8FLD_GIJ(1,1,2), C8FLD_GIJ(1,2,2), C8FLD_GIJ(1,3,2),                                   &
                                           C8FLD_GIJ(2,2,2), C8FLD_GIJ(2,3,2), C8FLD_GIJ(3,3,2), C8FLD_RHO(2),ICONT+2,             &
                                           ICONT+2, C8FLD_ALP(1,2), C8FLD_ALP(2,2), C8FLD_ALP(3,2), C8FLD_TREF(1)
         ELSE

            WRITE(F06,9923) MID2_PCOMP_EQ, C8FLD_GIJ(1,1,2), C8FLD_GIJ(1,2,2), C8FLD_GIJ(1,3,2),                                   &
                                           C8FLD_GIJ(2,2,2), C8FLD_GIJ(2,3,2), C8FLD_GIJ(3,3,2), C8FLD_RHO(2),ICONT+2,             &
                                           ICONT+2, C8FLD_TREF(1)
         ENDIF
      ENDIF

! Write MAT2 for transverse shear

      IF (FINITE_MAT_PROPS(3) == 'Y') THEN
         IF (IERR(3) == 0) THEN
            WRITE(F06,9932) MID3_PCOMP_EQ, C8FLD_GIJ(1,1,3), C8FLD_GIJ(1,2,3), '0.0000+0'      ,                                   &
                                           C8FLD_GIJ(2,2,3), '0.0000+0'      , '0.0000+0', C8FLD_RHO(3), ICONT+3,                  &
                                           ICONT+3, C8FLD_ALP(5,3), C8FLD_ALP(6,3), C8FLD_TREF(1)
         ELSE

            WRITE(F06,9933) MID3_PCOMP_EQ, C8FLD_GIJ(1,1,3), C8FLD_GIJ(1,2,3), '0.0000+0'      ,                                   &
                                           C8FLD_GIJ(2,2,3), '0.0000+0'      , '0.0000+0', C8FLD_RHO(3), ICONT+3,                  &
                                           ICONT+3, C8FLD_TREF(1)
         ENDIF
      ENDIF

! Write MAT2 for membrane/bending coupling

      IF (FINITE_MAT_PROPS(4) == 'Y') THEN
         IF (IERR(4) == 0) THEN
            WRITE(F06,9922) MID4_PCOMP_EQ, C8FLD_GIJ(1,1,4), C8FLD_GIJ(1,2,4), C8FLD_GIJ(1,3,4),                                   &
                                           C8FLD_GIJ(2,2,4), C8FLD_GIJ(2,3,4), C8FLD_GIJ(3,3,4), C8FLD_RHO(4),ICONT+4,             &
                                           ICONT+4, C8FLD_ALP(1,4), C8FLD_ALP(2,4), C8FLD_ALP(3,4), C8FLD_TREF(1)
         ELSE

            WRITE(F06,9923) MID4_PCOMP_EQ, C8FLD_GIJ(1,1,4), C8FLD_GIJ(1,2,4), C8FLD_GIJ(1,3,4),                                   &
                                           C8FLD_GIJ(2,2,4), C8FLD_GIJ(2,3,4), C8FLD_GIJ(3,3,4), C8FLD_RHO(4),ICONT+4,             &
                                           ICONT+4, C8FLD_TREF(1)
         ENDIF
      ENDIF

! Write message if we have changed SHELL_T to reflect approx zero transverse shear flex

      IF (SHELL_T_MOD == 'Y') THEN
         WRITE(ERR,9995) MID3_PCOMP_EQ
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9995) MID3_PCOMP_EQ
         ENDIF
      ENDIF

      WRITE(F06,9900)

      PCOMP(INTL_PID,6) = 1                                 ! Lets future calls to this subr know that PSHELL, MAT2 were written

! **********************************************************************************************************************************
 9801 FORMAT(' *INFORMATION: Cannot calculate equiv CTE''s for ',A,' for PCOMP ',I8,' since the det of matrix ',A,' is zero',/)

 9802 FORMAT(' *INFORMATION: Cannot calculate equiv CTE''s for ',A,' for PCOMP ',I8,' since matrix ',A,' cannot be inverted',/)

 9900 FORMAT('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++', &
             '+++++++++++++++++++++++++++++++++++++++++++++++++')

 9901 FORMAT(' *INFORMATION: Equivalent PSHELL amd MAT2 entries for PCOMP ',I8,' with:'                                            &
                    ,/,14X,' membrane thickness TM = ',1ES13.6,',   12*IB/(TM^3) = ',1ES13.6,',   TS/TM = ',1ES13.6,':',/)

 9902 FORMAT('        Type         MATL ID       G11            G12            G13            G22            G23            G33',  &
             '             A1             A2             A3',/                                                                     &
             ' ------------------ --------  -------------  -------------  -------------  -------------  -------------',            &
             '  -------------  -------------  -------------  -------------')

 9903 FORMAT('      membrane      ',I8,9(1ES15.6))

 9904 FORMAT('      bending       ',I8,9(1ES15.6))

 9905 FORMAT('  transverse shear  ',I8,2(1ES15.6),15X,1ES15.6,30X,2(1ES15.6))

 9906 FORMAT('  mem/bend coupling ',I8,9(1ES15.6))

 9912 FORMAT('PSHELL   ,',I9,',',I9,',',A9,',',I9,',      1.0,',I9,',','      1.0,         , +',I7,/,' +',I7,2(',',A9),I9,/)

 9913 FORMAT('PSHELL   ,',I9,',',I9,',',A9,',',I9,',      1.0,',I9,',','      1.0,         , +',I7,/,' +',I7,2(',',A9),/)

 9922 FORMAT('MAT2     ,',I9,7(',',A9),', +',I7,/,' +',I7,3(',',A9),',',A9,/)

 9923 FORMAT('MAT2     ,',I9,7(',',A9),', +',I7,/,' +',I7,3(',',8X),',',A9,/)

 9932 FORMAT('MAT2     ,',I9,7(',',A9),', +',I7,/,' +',I7,2(',',A9),',',9X,',',A9,/)

 9933 FORMAT('MAT2     ,',I9,7(',',A9),', +',I7,/,' +',I7,2(',',8X),',',9X,',',A9,/)

 9995 FORMAT(' *INFORMATION: The transverse shear modulii on the above MAT2 ',I8,' entry are the MYSTRAN calculated values to',    &
                           ' replace the zero input values by the user that',/,14X,' were meant to simulate zero transverse shear',&
                           ' flexibility')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_CHAR8_OUTPUTS ( C8FLD_GIJ, C8FLD_ALP, C8FLD_RHO, C8FLD_TREF, C8FLD_TM, C8FLD_ZS )

      IMPLICIT NONE

      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_GIJ(3,3,4)    ! Char representation of MAT2 Gij entries
      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_ALP(6,MEMATC) ! Char representation of MAT2 CTE entries
      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_RHO(4)        ! Char representation of MAT2 RHO entry
      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_TREF(4)       ! Char representation of MAT2 TREF entry
      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_TM            ! Char representation of MAT2 TM entry
      CHARACTER(8*BYTE), INTENT(OUT)  :: C8FLD_ZS(2)         ! Char representation of MAT2 ZS entry

      INTEGER(LONG)                   :: II,JJ,KK            ! DO loop indices

! **********************************************************************************************************************************
      DO II=1,3
         DO JJ=1,3
            DO KK=1,4
               C8FLD_GIJ(II,JJ,KK)(1:8) = ' '
            ENDDO
         ENDDO
      ENDDO

      DO II=1,6
         DO JJ=1,MEMATC
            C8FLD_ALP(II,JJ)(1:8) = ' '
         ENDDO
      ENDDO

! Put GIJ values in character format

      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(1,1) ,C8FLD_GIJ(1,1,1) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(1,2) ,C8FLD_GIJ(1,2,1) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(1,3) ,C8FLD_GIJ(1,3,1) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(2,2) ,C8FLD_GIJ(2,2,1) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(2,3) ,C8FLD_GIJ(2,3,1) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MEM_MAT2(3,3) ,C8FLD_GIJ(3,3,1) )

      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(1,1) ,C8FLD_GIJ(1,1,2) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(1,2) ,C8FLD_GIJ(1,2,2) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(1,3) ,C8FLD_GIJ(1,3,2) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(2,2) ,C8FLD_GIJ(2,2,2) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(2,3) ,C8FLD_GIJ(2,3,2) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_BEN_MAT2(3,3) ,C8FLD_GIJ(3,3,2) )

      CALL REAL_DATA_TO_C8FLD ( PCOMP_TSH_MAT2(1,1) ,C8FLD_GIJ(1,1,3) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_TSH_MAT2(1,2) ,C8FLD_GIJ(1,2,3) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_TSH_MAT2(2,2) ,C8FLD_GIJ(2,2,3) )

      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(1,1) ,C8FLD_GIJ(1,1,4) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(1,2) ,C8FLD_GIJ(1,2,4) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(1,3) ,C8FLD_GIJ(1,3,4) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(2,2) ,C8FLD_GIJ(2,2,4) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(2,3) ,C8FLD_GIJ(2,3,4) )
      CALL REAL_DATA_TO_C8FLD ( PCOMP_MBC_MAT2(3,3) ,C8FLD_GIJ(3,3,4) )

! Put CTE's in character format

      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(1,1), C8FLD_ALP(1,1) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(2,1), C8FLD_ALP(2,1) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(3,1), C8FLD_ALP(3,1) )

      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(1,2), C8FLD_ALP(1,2) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(2,2), C8FLD_ALP(2,2) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(3,2), C8FLD_ALP(3,2) )

      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(5,3), C8FLD_ALP(5,3) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(6,3), C8FLD_ALP(6,3) )

      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(1,4), C8FLD_ALP(1,4) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(2,4), C8FLD_ALP(2,4) )
      CALL REAL_DATA_TO_C8FLD ( SHELL_ALP(3,4), C8FLD_ALP(3,4) )

! Put RHO's in character format

      CALL REAL_DATA_TO_C8FLD ( RHO(1), C8FLD_RHO(1) )
      CALL REAL_DATA_TO_C8FLD ( RHO(2), C8FLD_RHO(2) )
      CALL REAL_DATA_TO_C8FLD ( RHO(3), C8FLD_RHO(3) )
      CALL REAL_DATA_TO_C8FLD ( RHO(4), C8FLD_RHO(4) )

! Put TREF's in character format

      CALL REAL_DATA_TO_C8FLD ( TREF(1), C8FLD_TREF(1) )
      CALL REAL_DATA_TO_C8FLD ( TREF(2), C8FLD_TREF(2) )
      CALL REAL_DATA_TO_C8FLD ( TREF(3), C8FLD_TREF(3) )
      CALL REAL_DATA_TO_C8FLD ( TREF(4), C8FLD_TREF(4) )

! Put TM in character format

      CALL REAL_DATA_TO_C8FLD ( PCOMP_TM, C8FLD_TM )

! Put ZS's in character format

      CALL REAL_DATA_TO_C8FLD ( ZS(1), C8FLD_ZS(1) )
      CALL REAL_DATA_TO_C8FLD ( ZS(2), C8FLD_ZS(2) )

! **********************************************************************************************************************************

      END SUBROUTINE GET_CHAR8_OUTPUTS

      END SUBROUTINE WRITE_PCOMP_EQUIV


      SUBROUTINE GET_ELEM_NUM_PLIES ( INT_ELEM_ID )

! Gets shell element number of plies (1 unless elem uses PCOMP props) given the element's internal ID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DEDAT_Q4_SHELL_KEY, DEDAT_T3_SHELL_KEY, NPCOMP, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, EID, EPNT, ETYPE, INTL_PID, NUM_PLIES, PCOMP, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_ELEM_NUM_PLIES'

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID for which
      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                            row number in array EDAT where data begins for this element.
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I1                ! Index into EDAT




! **********************************************************************************************************************************
      NUM_PLIES = 1

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN

         EPNTK    = EPNT(INT_ELEM_ID)
         TYPE     = ETYPE(INT_ELEM_ID)                     ! NOTE: Must keep (this subr not always called when TYPE     is known)
         EID      = EDAT(EPNTK)                            ! NOTE: Must keep (this subr not always called when EID      is known)
         INTL_PID = EDAT(EPNTK+1)                          ! NOTE: Must keep (this subr not always called when INTL_PID is known)

         IF      ((TYPE == 'TRMEM   ') .OR. (TYPE == 'TRPLT1  ') .OR. (TYPE == 'TRPLT2  ') .OR.                                    &
                  (TYPE == 'TRIA3K  ') .OR. (TYPE == 'TRIA3   ')) THEN
             I1 = DEDAT_T3_SHELL_KEY
         ELSE IF ((TYPE == 'QDMEM   ') .OR. (TYPE == 'QDPLT1  ') .OR. (TYPE == 'QDPLT2  ') .OR.                                    &
                  (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD4   ')) THEN
             I1 = DEDAT_Q4_SHELL_KEY
         ENDIF

         IF (EDAT(EPNTK+I1) == 2) THEN
            DO I=1,NPCOMP
               NUM_PLIES = PCOMP(INTL_PID,5)
            ENDDO
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_ELEM_NUM_PLIES

   END MODULE COMPOSITE_SHELL_PREPARATION
