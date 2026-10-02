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

   MODULE RIGID_BODY_STORAGE_LIFECYCLE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ALLOCATE_CB_ELM_OTM, ALLOCATE_CB_GRD_OTM, ALLOCATE_RBGLOBAL, DEALLOCATE_RBGLOBAL

   CONTAINS

      SUBROUTINE ALLOCATE_CB_ELM_OTM ( NAME_IN )

! Calculates how many rows/cols are going to be needed for element related OTM's (Output Transformation Matrices) for Craig-Bampton
! model generation runs and allocates memory to the arrays

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR,                                                                  &
                                         ELOUT_ELFE_BIT, ELOUT_ELFN_BIT, ELOUT_STRE_BIT, ELOUT_STRN_BIT,                           &
                                         IBIT, NELE, NUM_CB_DOFS,                                                                  &
                                         NROWS_OTM_ELFE, NROWS_OTM_ELFN, NROWS_OTM_STRE, NROWS_OTM_STRN,                           &
                                         NROWS_TXT_ELFE, NROWS_TXT_ELFN, NROWS_TXT_STRE, NROWS_TXT_STRN, TOT_MB_MEM_ALLOC


      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE MODEL_STUF, ONLY            :  ELGP, ELOUT, ELMTYP, ETYPE, METYPE, NELGP, TYPE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFE, OTM_ELFN, OTM_STRE, OTM_STRN, TXT_ELFE, TXT_ELFN, TXT_STRE, TXT_STRN
      USE PARAMS, ONLY                :  OTMSKIP

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_CB_ELM_OTM'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Array name of the matrix to be allocated
      CHARACTER(LEN(NAME_IN))         :: NAME              ! Name for output error purposes

      INTEGER(LONG)                   :: ELOUT_ELFE        ! If > 0, there are ELFORCE(ENGR) requests for some elems
      INTEGER(LONG)                   :: ELOUT_ELFN        ! If > 0, there are ELFORCE(NODE) requests for some elems
      INTEGER(LONG)                   :: ELOUT_STRE        ! If > 0, there are STRESS   requests for some elems
      INTEGER(LONG)                   :: ELOUT_STRN        ! If > 0, there are STRAIN   requests for some elems
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NCOLS             ! Number of cols in OTM matrix
      INTEGER(LONG)                   :: NROWS_MAT         ! Number of rows in OTM matrix
      INTEGER(LONG)                   :: NROWS_TXT         ! Number of rows in TXT mmatrix


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: IAND



! **********************************************************************************************************************************
      JERR = 0

! ----------------------------------------------------------------------------------------------------------------------------------
      IF (NAME_IN == 'OTM_ELFE') THEN                      ! Determine size of OTM_ELFE and allocate it and TXT_ELFE

         NROWS_MAT = 0
         NROWS_TXT = 0
         DO I=1,METYPE

             DO J = 1,NELE

               TYPE  = ETYPE(J)

               IF      (ETYPE(J)(1:4) == 'ELAS') THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,1),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NROWS_MAT = NROWS_MAT + 1
                        NROWS_TXT = NROWS_TXT + 1 + OTMSKIP
                     ENDIF
                  ENDIF

               ELSE IF (ETYPE(J)(1:4) == 'BUSH') THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,1),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NROWS_MAT = NROWS_MAT + 6
                        NROWS_TXT = NROWS_TXT + 6 + OTMSKIP
                     ENDIF
                  ENDIF

               ELSE IF (ETYPE(J)(1:3) == 'ROD') THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,1),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NROWS_MAT = NROWS_MAT + 2
                        NROWS_TXT = NROWS_TXT + 2 + OTMSKIP
                     ENDIF
                  ENDIF

               ELSE IF (ETYPE(J)(1:3) == 'BAR') THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,1),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NROWS_MAT = NROWS_MAT + 8
                        NROWS_TXT = NROWS_TXT + 8 + OTMSKIP
                     ENDIF
                  ENDIF

               ELSE IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'SHEAR') .OR.             &
                       (ETYPE(J)(1:6) == 'USERIN')) THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,1),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NROWS_MAT = NROWS_MAT + 8
                        NROWS_TXT = NROWS_TXT + 8 + OTMSKIP
                     ENDIF
                   ENDIF

               ELSE IF((ETYPE(J)(1:4) == 'HEXA'  ) .OR. (ETYPE(J)(1:5) == 'PENTA' ) .OR. (ETYPE(J)(1:5) == 'TETRA' ) .OR.          &
                       (ETYPE(J)(1:5) == 'USER1' ) .OR. (ETYPE(J)(1:6) == 'PLOTEL')) THEN
                  CONTINUE                                 ! No element forces for these elements

               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,963) SUBR_NAME, TYPE
                  WRITE(F06,963) SUBR_NAME, TYPE
                  CALL OUTA_HERE ( 'Y' )

               ENDIF

            ENDDO

         ENDDO

         NROWS_OTM_ELFE = NROWS_MAT
         NROWS_TXT_ELFE = NROWS_TXT
         NCOLS          = NUM_CB_DOFS

         NAME = 'OTM_ELFE'
         IF (ALLOCATED(OTM_ELFE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_ELFE(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_ELFE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_ELFE'
         IF (ALLOCATED(TXT_ELFE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_ELFE(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_ELFE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_ELFN') THEN                 ! Determine size of OTM_ELFN and allocate it and TXT_ELFN

         NROWS_MAT = 0
         NROWS_TXT = 0
         DO I=1,METYPE
            DO J = 1,NELE
               TYPE  = ETYPE(J)
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_ELFN   = IAND(ELOUT(J,1), IBIT(ELOUT_ELFN_BIT))
                  IF (ELOUT_ELFN > 0) THEN
                     NROWS_MAT = NROWS_MAT + 6*NELGP(I)
                     NROWS_TXT = NROWS_TXT + 6*NELGP(I) + OTMSKIP
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         NROWS_OTM_ELFN = NROWS_MAT
         NROWS_TXT_ELFN = NROWS_TXT
         NCOLS          = NUM_CB_DOFS

         NAME = 'OTM_ELFN'
         IF (ALLOCATED(OTM_ELFN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_ELFN(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_ELFN(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_ELFN'
         IF (ALLOCATED(TXT_ELFN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_ELFN(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_ELFN(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_STRE') THEN                 ! Determine size of OTM_STRE and allocate it and TXT_STRE

         NROWS_MAT = 0
         NROWS_TXT = 0
         DO I=1,METYPE
            DO J = 1,NELE
               TYPE  = ETYPE(J)
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_STRE = IAND(ELOUT(J,1),IBIT(ELOUT_STRE_BIT))
                  IF (ELOUT_STRE > 0) THEN
                     IF (TYPE(1:4) == 'ELAS') THEN
                        NROWS_MAT = NROWS_MAT + 1
                        NROWS_TXT = NROWS_TXT + 1 + OTMSKIP
                     ELSE IF (TYPE == 'ROD     ') THEN
                        NROWS_MAT = NROWS_MAT + 4
                        NROWS_TXT = NROWS_TXT + 4 + OTMSKIP
                     ELSE IF (TYPE == 'BAR     ') THEN
                        NROWS_MAT = NROWS_MAT + 18
                        NROWS_TXT = NROWS_TXT + 18 + OTMSKIP
                     ELSE IF (TYPE(1:5) == 'QUAD4') THEN
                        IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
                           NROWS_MAT = NROWS_MAT + 100     ! 10 stresses each for -Z1, Z1 mult by 5 (CENTER + 4 grids)
                           NROWS_TXT = NROWS_TXT + 100 + OTMSKIP
                        ELSE
                           NROWS_MAT = NROWS_MAT + 20
                           NROWS_TXT = NROWS_TXT + 20 + OTMSKIP
                        ENDIF
                     ELSE IF (TYPE(1:5) == 'TRIA3') THEN
                        NROWS_MAT = NROWS_MAT + 20
                        NROWS_TXT = NROWS_TXT + 20 + OTMSKIP
                     ELSE IF (TYPE == 'USERIN  ') THEN
                        NROWS_MAT = NROWS_MAT + 9
                        NROWS_TXT = NROWS_TXT + 9 + OTMSKIP
                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
                        NROWS_MAT = NROWS_MAT + 6
                        NROWS_TXT = NROWS_TXT + 6 + OTMSKIP
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         NROWS_OTM_STRE = NROWS_MAT
         NROWS_TXT_STRE = NROWS_TXT
         NCOLS          = NUM_CB_DOFS

         NAME = 'OTM_STRE'
         IF (ALLOCATED(OTM_STRE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_STRE(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_STRE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_STRE'
         IF (ALLOCATED(TXT_STRE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_STRE(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_STRE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_STRN') THEN                 ! Determine size of OTM_STRN and allocate it and TXT_STRN

         NROWS_MAT = 0
         NROWS_TXT = 0
         DO I=1,METYPE
            DO J = 1,NELE
               TYPE  = ETYPE(J)
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_STRN = IAND(ELOUT(J,1),IBIT(ELOUT_STRN_BIT))
                  IF (ELOUT_STRN > 0) THEN
                     IF      (TYPE(1:5) == 'QUAD4') THEN
                        IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
                           NROWS_MAT = NROWS_MAT + 100     ! 10 strains each for -Z1, Z1 mult by 5 (CENTER + 4 grids)
                           NROWS_TXT = NROWS_TXT + 100 + OTMSKIP
                        ELSE
                           NROWS_MAT = NROWS_MAT + 20
                           NROWS_TXT = NROWS_TXT + 20 + OTMSKIP
                        ENDIF
                     ELSE IF (TYPE(1:5) == 'TRIA3') THEN
                        NROWS_MAT = NROWS_MAT + 20
                        NROWS_TXT = NROWS_TXT + 20 + OTMSKIP
                     ELSE IF (TYPE == 'USERIN  ') THEN
                        NROWS_MAT = NROWS_MAT + 9
                        NROWS_TXT = NROWS_TXT + 9 + OTMSKIP
                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
                        NROWS_MAT = NROWS_MAT + 6
                        NROWS_TXT = NROWS_TXT + 6 + OTMSKIP
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         NROWS_OTM_STRN = NROWS_MAT
         NROWS_TXT_STRN = NROWS_TXT
         NCOLS          = NUM_CB_DOFS

         NAME = 'OTM_STRN'
         IF (ALLOCATED(OTM_STRN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_STRN(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_STRN(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_STRN'
         IF (ALLOCATED(TXT_STRN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_STRN(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_STRN(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  963 FORMAT(' *ERROR   946: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NO CODE FOR ELEMENT TYPE "',A,'"')

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_CB_ELM_OTM



      SUBROUTINE ALLOCATE_CB_GRD_OTM ( NAME_IN )

! Calculates how many rows/cols are going to be needed for grid related OTM's (Output Transformation Matrices) for Craig-Bampton
! model generation runs and allocates memory to the arrays

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC,                                                &
                                         GROUT_ACCE_BIT, GROUT_DISP_BIT, GROUT_SPCF_BIT, GROUT_MPCF_BIT,                           &
                                         IBIT, NDOFR, NGRID, NUM_CB_DOFS, NVEC,                                                    &
                                         NROWS_OTM_ACCE, NROWS_OTM_DISP, NROWS_OTM_MPCF, NROWS_OTM_SPCF,                           &
                                         NROWS_TXT_ACCE, NROWS_TXT_DISP, NROWS_TXT_MPCF, NROWS_TXT_SPCF
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE PARAMS, ONLY                :  OTMSKIP
      USE MODEL_STUF, ONLY            :  GRID, GROUT
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ACCE, OTM_DISP, OTM_MPCF, OTM_SPCF, TXT_ACCE, TXT_DISP, TXT_MPCF, TXT_SPCF

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_COL_VEC'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Array name of the matrix to be allocated
      CHARACTER(LEN(NAME_IN))         :: NAME              ! Name for output error purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IB                ! If > 0, there are displ or applied load output requests
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NCOLS             ! Number of cols in OTM matrix
      INTEGER(LONG)                   :: NROWS_MAT         ! Number of rows in OTM matrix
      INTEGER(LONG)                   :: NROWS_TXT         ! Number of rows in TXT mmatrix


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: IAND



! **********************************************************************************************************************************
      JERR = 0

      IF (NAME_IN == 'OTM_ACCE') THEN

         NROWS_MAT = 0                                     ! Determine how many grids have accel output requested
         NROWS_TXT = 0
         DO I=1,NGRID
            IB = IAND(GROUT(I,1),IBIT(GROUT_ACCE_BIT))
            IF (IB > 0) THEN
               NROWS_MAT = NROWS_MAT + GRID(I,6)
               NROWS_TXT = NROWS_TXT + GRID(I,6) + OTMSKIP
            ENDIF
         ENDDO
         NROWS_OTM_ACCE = NROWS_MAT
         NROWS_TXT_ACCE = NROWS_TXT
         NCOLS = NUM_CB_DOFS

         NAME = 'OTM_ACCE'
         IF (ALLOCATED(OTM_ACCE)) THEN                     ! If not already allocated, allocate OTM_ACCE
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_ACCE(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_ACCE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_ACCE'
         IF (ALLOCATED(TXT_ACCE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_ACCE(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_ACCE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_DISP') THEN

         NROWS_MAT = 0                                     ! Determine how many grids have accel output requested
         NROWS_TXT = 0
         DO I=1,NGRID
            IB = IAND(GROUT(I,1),IBIT(GROUT_DISP_BIT))
            IF (IB > 0) THEN
               NROWS_MAT = NROWS_MAT + GRID(I,6)
               NROWS_TXT = NROWS_TXT + GRID(I,6) + OTMSKIP
            ENDIF
         ENDDO
         NROWS_OTM_DISP = NROWS_MAT
         NROWS_TXT_DISP = NROWS_TXT
         NCOLS = NUM_CB_DOFS

         NAME = 'OTM_DISP'
         IF (ALLOCATED(OTM_DISP)) THEN                     ! If not already allocated, allocateDISP _OTM
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_DISP(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_DISP(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_DISP'
         IF (ALLOCATED(TXT_DISP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_DISP(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_DISP(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_MPCF') THEN

         NROWS_MAT = 0                                     ! Determine how many grids have accel output requested
         NROWS_TXT = 0
         DO I=1,NGRID
            IB = IAND(GROUT(I,1),IBIT(GROUT_MPCF_BIT))
            IF (IB > 0) THEN
               NROWS_MAT = NROWS_MAT + GRID(I,6)
               NROWS_TXT = NROWS_TXT + GRID(I,6) + OTMSKIP
            ENDIF
         ENDDO
         NROWS_OTM_MPCF = NROWS_MAT
         NROWS_TXT_MPCF = NROWS_TXT
         NCOLS = NUM_CB_DOFS

         MB_ALLOCATED = REAL(DOUBLE)*NROWS_MAT*NCOLS/ONEPP6

         NAME = 'OTM_MPCF'
         IF (ALLOCATED(OTM_MPCF)) THEN                     ! If not already allocated, allocate OTM_MPCF
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_MPCF(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_MPCF(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_MPCF'
         IF (ALLOCATED(TXT_MPCF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_MPCF(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_MPCF(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (NAME_IN == 'OTM_SPCF') THEN

         NROWS_MAT = 0                                     ! Determine how many grids have accel output requested
         NROWS_TXT = 0
         DO I=1,NGRID
            IB = IAND(GROUT(I,1),IBIT(GROUT_SPCF_BIT))
            IF (IB > 0) THEN
               NROWS_MAT = NROWS_MAT + GRID(I,6)
               NROWS_TXT = NROWS_TXT + GRID(I,6) + OTMSKIP
            ENDIF
         ENDDO
         NROWS_OTM_SPCF = NROWS_MAT
         NROWS_TXT_SPCF = NROWS_TXT
         NCOLS = NUM_CB_DOFS

         MB_ALLOCATED = REAL(DOUBLE)*NROWS_MAT*NCOLS/ONEPP6

         NAME = 'OTM_SPCF'
         IF (ALLOCATED(OTM_SPCF)) THEN                     ! If not already allocated, allocate OTM_SPCF
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OTM_SPCF(NROWS_MAT,NCOLS),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS_MAT)*REAL(NCOLS)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_MAT
                  DO J=1,NCOLS
                     OTM_SPCF(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TXT_SPCF'
         IF (ALLOCATED(TXT_SPCF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TXT_SPCF(NROWS_TXT),STAT=IERR)
            IF (IERR == 0) THEN
               MB_ALLOCATED = REAL(NROWS_TXT)/ONEPP6
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NROWS_TXT
                  TXT_SPCF(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME, SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_CB_GRD_OTM



      SUBROUTINE ALLOCATE_RBGLOBAL ( SET, CALLING_SUBR )

! Allocate arrays for rigid body displ matrices (except RBM0 used in Craig-Bampton analyses). RBGLOBAL matrices are used in
! stiffness matrix equilibrium checks. The TR6 matrices are used in transforming some Craig-Bampton matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFG, NDOFN, NDOFF, NDOFA, NDOFL, NDOFR, BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, RBGLOBAL_NSET, RBGLOBAL_FSET, RBGLOBAL_ASET, RBGLOBAL_LSET,                &
                                         TR6_CG, TR6_MEFM, TR6_0

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_RBGLOBAL'
      CHARACTER(LEN=*), INTENT(IN)    :: SET               ! Set name of the displ matrix
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(14*BYTE)              :: NAME              ! Specific array name used for output error message

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Number of rows in array
      INTEGER(LONG), PARAMETER        :: NCOLS       = 6   ! Number of cols in array


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      MB_ALLOCATED = ZERO
      JERR = 0

      IF      (SET == 'G ') THEN                           ! Allocate array for G-set rigid body disp matrix

         NROWS = NDOFG
         NAME = 'RBGLOBAL_GSET'
         IF (ALLOCATED(RBGLOBAL_GSET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RBGLOBAL_GSET(NDOFG,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFG
                  DO J=1,6
                     RBGLOBAL_GSET(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'N ') THEN                           ! Allocate array for N-set rigid body disp matrix

         NROWS = NDOFN
         NAME = 'RBGLOBAL_NSET'
         IF (ALLOCATED(RBGLOBAL_NSET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RBGLOBAL_NSET(NDOFN,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFN
                  DO J=1,6
                     RBGLOBAL_NSET(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'F ') THEN                           ! Allocate array for F-set rigid body disp matrix

         NROWS = NDOFF
         NAME = 'RBGLOBAL_FSET'
         IF (ALLOCATED(RBGLOBAL_FSET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RBGLOBAL_FSET(NDOFF,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFF
                  DO J=1,6
                     RBGLOBAL_FSET(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'A ') THEN                           ! Allocate array for A-set rigid body disp matrix

         NROWS = NDOFA
         NAME = 'RBGLOBAL_ASET'
         IF (ALLOCATED(RBGLOBAL_ASET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RBGLOBAL_ASET(NDOFA,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFA
                  DO J=1,6
                     RBGLOBAL_ASET(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'L ') THEN                           ! Allocate array for L-set rigid body disp matrix

         NROWS = NDOFL
         NAME = 'RBGLOBAL_LSET'
         IF (ALLOCATED(RBGLOBAL_LSET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RBGLOBAL_LSET(NDOFL,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFL
                  DO J=1,6
                     RBGLOBAL_LSET(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'R ') THEN                           ! Allocate array for R-set rigid body disp matrix

         NROWS = NDOFR
         NAME = 'TR6_CG'
         IF (ALLOCATED(TR6_CG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TR6_CG(NDOFR,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFR
                  DO J=1,6
                     TR6_CG(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NROWS = NDOFR
         NAME = 'TR6_MEFM'
         IF (ALLOCATED(TR6_MEFM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TR6_MEFM(NDOFR,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFR
                  DO J=1,6
                     TR6_MEFM(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NROWS = NDOFR
         NAME = 'TR6_0'
         IF (ALLOCATED(TR6_0)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TR6_0(NDOFR,6),STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,NDOFR
                  DO J=1,6
                     TR6_0(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,921) SUBR_NAME, 'ALLOCATED', SET
         WRITE(F06,921) SUBR_NAME, 'ALLOCATED', SET
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS)*REAL(NCOLS)/ONEPP6
      CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  921 FORMAT(' *ERROR   921: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DISPL SET FOR MATRIX TO BE ',A,' IS INCORRECT. SET DESIGNATION WAS "',A,'"')

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_RBGLOBAL


      SUBROUTINE DEALLOCATE_CB_ELM_OTM ( NAME )

! Deallocates memory from the elem related OTM (Output Transformation Matrices) used in Craig-Bampton model generation

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE OUTPUT4_MATRICES

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_CB_ELM_OTM'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF      (NAME == 'OTM_ELFE') THEN                    ! Dellocate array for OTM_ELFE

         IF (ALLOCATED(OTM_ELFE)) THEN
            DEALLOCATE (OTM_ELFE,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'OTM_ELFN') THEN                    ! Dellocate array for OTM_ELFN

         IF (ALLOCATED(OTM_ELFN)) THEN
            DEALLOCATE (OTM_ELFN,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'OTM_STRE') THEN                    ! Dellocate array for OTM_STRE

         IF (ALLOCATED(OTM_STRE)) THEN
            DEALLOCATE (OTM_STRE,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_CB_ELM_OTM


      SUBROUTINE DEALLOCATE_CB_GRD_OTM ( NAME )

! Deallocates memory from the grid related OTM (Output Transformation Matrices) used in Craig-Bampton model generation

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ACCE, OTM_DISP, OTM_MPCF, OTM_SPCF

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_CB_GRD_OTM'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF      (NAME == 'OTM_ACCE') THEN                    ! Dellocate array for OTM_ACCE

         IF (ALLOCATED(OTM_ACCE)) THEN
            DEALLOCATE (OTM_ACCE,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'OTM_DISP') THEN                    ! Dellocate array for OTM_DISP

         IF (ALLOCATED(OTM_DISP)) THEN
            DEALLOCATE (OTM_DISP,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'OTM_MPCF') THEN                    ! Dellocate array for OTM_MPCF

         IF (ALLOCATED(OTM_MPCF)) THEN
            DEALLOCATE (OTM_MPCF,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME == 'OTM_SPCF') THEN                    ! Dellocate array for OTM_SPCF

         IF (ALLOCATED(OTM_SPCF)) THEN
            DEALLOCATE (OTM_SPCF,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_CB_GRD_OTM


      SUBROUTINE DEALLOCATE_RBGLOBAL ( SET )

! Deallocate arrays for rigid body displ matrices (except RBM0 used in Craig-Bampton analyses). RBGLOBAL matrices are used in
! stiffness matrix equilibrium checks. The TR6 matrices are used in transforming some Craig-Bampton matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, RBGLOBAL_NSET, RBGLOBAL_FSET, RBGLOBAL_ASET, RBGLOBAL_LSET,                &
                                         TR6_CG, TR6_MEFM, TR6_0

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_RBGLOBAL'
      CHARACTER(LEN=*), INTENT(IN)    :: SET               ! Set name of the displ matrix
      CHARACTER(13*BYTE)              :: NAME              ! Specific array name used for output error message

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF      (SET == 'G ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'RBGLOBAL_GSET'
         IF (ALLOCATED(RBGLOBAL_GSET)) THEN
            DEALLOCATE (RBGLOBAL_GSET,STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'N ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'RBGLOBAL_NSET'
         IF (ALLOCATED(RBGLOBAL_NSET)) THEN
            DEALLOCATE (RBGLOBAL_NSET,STAT=IERR)
               NAME = 'RBGLOBAL_NSET'
            IF (IERR /= 0) THEN
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'F ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'RBGLOBAL_FSET'
         IF (ALLOCATED(RBGLOBAL_FSET)) THEN
            DEALLOCATE (RBGLOBAL_FSET,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'RBGLOBAL_FSET'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (SET == 'A ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'RBGLOBAL_ASET'
         IF (ALLOCATED(RBGLOBAL_ASET)) THEN
            DEALLOCATE (RBGLOBAL_ASET,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'RBGLOBAL_ASET'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

       ELSE IF (SET == 'L ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'RBGLOBAL_LSET'
         IF (ALLOCATED(RBGLOBAL_LSET)) THEN
            DEALLOCATE (RBGLOBAL_LSET,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'RBGLOBAL_LSET'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

       ELSE IF (SET == 'R ') THEN                           ! Dellocate array for RBGLOBAL_GSET

         NAME = 'TR6_CG'
         IF (ALLOCATED(TR6_CG)) THEN
            DEALLOCATE (TR6_CG,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'TR6_CG'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TR6_MEFM'
         IF (ALLOCATED(TR6_MEFM)) THEN
            DEALLOCATE (TR6_MEFM,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'TR6_MEFM'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TR6_0'
         IF (ALLOCATED(TR6_0)) THEN
            DEALLOCATE (TR6_0,STAT=IERR)
            IF (IERR /= 0) THEN
               NAME = 'TR6_0'
               WRITE(ERR,992) NAME,SUBR_NAME
               WRITE(F06,992) NAME,SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,901) SUBR_NAME, 'DEALLOCATED', SET
         WRITE(F06,901) SUBR_NAME, 'DEALLOCATED', SET
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + JERR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  901 FORMAT(' *ERROR   901: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DISPL SET FOR MATRIX TO BE ',A,' IS INCORRECT. SET DESIGNATION WAS "',A,'"')

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_RBGLOBAL

   END MODULE RIGID_BODY_STORAGE_LIFECYCLE
