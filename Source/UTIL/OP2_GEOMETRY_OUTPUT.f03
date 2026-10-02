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

   MODULE OP2_GEOMETRY_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_OP2_HEADER, END_OP2_TABLE, END_OP2_TABLES, WRITE_TABLE_HEADER, WRITE_ITABLE

   CONTAINS
!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM1_GRID(ITABLE, NGRID_ACTUAL, GRID_INDEX)
      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  ERR, OP2
      USE SCONTR, ONLY                :  NGRID
      USE MODEL_STUF, ONLY            :  GRID_ID, GRID, RGRID
      IMPLICIT NONE
      INTEGER(LONG),INTENT(INOUT)                  :: ITABLE
      INTEGER(LONG),INTENT(IN)                     :: NGRID_ACTUAL
      INTEGER(LONG), INTENT(IN), DIMENSION(NGRID)  :: GRID_INDEX      ! oversized, but mehhh
      INTEGER(LONG)                                :: I, J
      INTEGER(LONG)                                :: NVALUES
      INTEGER(LONG),PARAMETER                      :: NUM_WIDE = 8
      INTEGER(LONG),PARAMETER                      :: SEID = 0

      WRITE(ERR,1) "WRITE_OP2_GEOM1_GRID"
      WRITE(ERR,2) "NGRID_ACTUAL=",NGRID_ACTUAL
      IF(NGRID_ACTUAL > 0) THEN
 1    FORMAT("****DEBUG:   ", A)
 2    FORMAT("****DEBUG:   ", A, i4)
 3    FORMAT("****DEBUG:   GRID_ID=",i4, "; ID=",i4, "; cp=", i4, "; x=", f8.3, "; y=", f8.3, "; z=", f8.3, &
             "; cd=", i4, "; PS=", i4,"; NDOF=",i4)

        ! if there are GRIDs
        NVALUES = NUM_WIDE * NGRID_ACTUAL

        DO J=1,NGRID_ACTUAL
          ! nid, cp, x, y, z, cd, ps, seid
          I = GRID_INDEX(J)
          WRITE(ERR,3) GRID_ID(I), GRID(I,1), GRID(I,2), REAL(RGRID(I,1),4), REAL(RGRID(I,2),4), REAL(RGRID(I,3),4), &
                       GRID(I,3), GRID(I,4), GRID(I,6)
        ENDDO
        ! oh dear...
        !
        ! the core of what we want is...
        !WRITE(OP2) (GRID_ID(I), GRID(I,2),                 &
        !            RGRID(I,1), RGRID(I,2), RGRID(I,3),    &
        !            GRID(I,3), GRID(I,4), SEID, J=1,NGRID_ACTUAL)
        !
        ! We then sub in I = GRID_INDEX(J), which is required because
        ! mystran stores the GRIDs and SPOINTs in a single array.
        !
        ! Then we wrap any doubles in REAL(*,4) to cast it as a float.
        !
        ! Finally we tack on the magic code (4501,45,1).
        ! Thanks to that code, we have to add 3 the nvalues...
        !
        WRITE(OP2) NVALUES + 3
        !           nid, cp, x, y, z, cd, ps, seid
        WRITE(OP2) 4501, 45, 1,                                                                                     &
                    (GRID_ID(GRID_INDEX(J)),         GRID(GRID_INDEX(J),2),                                         &
                    REAL(RGRID(GRID_INDEX(J),1),4), REAL(RGRID(GRID_INDEX(J),2),4), REAL(RGRID(GRID_INDEX(J),3),4), &
                    GRID(GRID_INDEX(J),3), GRID(GRID_INDEX(J),4), SEID, J=1,NGRID_ACTUAL)

        CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM1_GRID

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM1_CORD(ITABLE)
      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  NCORD
      USE MODEL_STUF, ONLY            : CORD
      IMPLICIT NONE
      INTEGER(LONG),INTENT(INOUT)      :: ITABLE
      INTEGER(LONG)                    :: I ! loop index
      INTEGER(LONG)                    :: NCORD1R = 0, NCORD1C = 0, NCORD1S = 0
      INTEGER(LONG)                    :: NCORD2R = 0, NCORD2C = 0, NCORD2S = 0
      INTEGER(LONG), DIMENSION(NCORD)  :: CORD1R_INDEX, CORD1C_INDEX, CORD1S_INDEX, CORD2R_INDEX, CORD2C_INDEX, CORD2S_INDEX      ! oversized, but mehhh

      ! Each row of CORD, for a CORD1R,C,S has:
      !    Col. 1: Coord type
      !            11 for CORD1R
      !            12 for CORD1C
      !            13 for CORD1S
      !    Col. 2: CID = Coord system ID
      !    Col. 3: RID = Reference coord system ID for Grid A
      !    Col. 4: RID = Reference coord system ID for Grid B
      !    Col. 5: RID = Reference coord system ID for Grid C
      !
      ! Each row of CORD, for a CORD2R,C,S has the same as for CORD1 above but only needs 1 reference system:
      !    Col. 1: Coord type
      !            21 for CORD2R
      !            22 for CORD2C
      !            23 for CORD2S
      !    Col. 2: CID = Coord system ID
      !    Col. 3: RID = Reference coord system ID for this CORD2R,C,S
      !
      ! Each row of RCORD has:
      !   After Bulk Data has been read:
      !     Col.  1: X coord of Pt A in RID coord system
      !     Col.  2: Y coord of Pt A in RID coord system
      !     Col.  3: Z coord of Pt A in RID coord system
      !     Col.  4: X coord of Pt B in RID coord system
      !     Col.  5: Y coord of Pt B in RID coord system
      !     Col.  6: Z coord of Pt B in RID coord system
      !     Col.  7: X coord of Pt C in RID coord system
      !     Col.  8: Y coord of Pt C in RID coord system
      !     Col.  9: Z coord of Pt C in RID coord system
      !     Col. 10: 0.
      !     Col. 11: 0.
      !     Col. 12: 0.
      !
      !   After subr CORD_PROC has run:
      !     Col.  1: Location of CID origin in basic coord sys X dir
      !     Col.  2: Location of CID origin in basic coord sys Y dir
      !     Col.  3: Location of CID origin in basic coord sys Z dir
      !     Col.  4: TN(1,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
      !     Col.  5: TN(1,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 1st row of TN in cols  4- 6 of RCORD
      !     Col.  6: TN(1,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
      !     Col.  7: TN(2,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
      !     Col.  8: TN(2,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 2nd row of TN in cols  7- 9 of RCORD
      !     Col.  9: TN(2,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
      !     Col. 10: TN(3,1,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]
      !     Col. 11: TN(3,2,M) for M = row no. in RCORD (internal coord. sys. ID = M)  |- 3rd row of TN in cols 10-12 of RCORD
      !     Col. 12: TN(3,3,M) for M = row no. in RCORD (internal coord. sys. ID = M)  ]


      ! Each row of CORD, for a CORD1R,C,S has:
      !    Col. 1: Coord type
      !            11 for CORD1R
      !            12 for CORD1C
      !            13 for CORD1S
      !    Col. 2: CID = Coord system ID
      !    Col. 3: RID = Reference coord system ID for Grid A
      !    Col. 4: RID = Reference coord system ID for Grid B
      !    Col. 5: RID = Reference coord system ID for Grid C
      !
      ! Each row of CORD, for a CORD2R,C,S has the same as for CORD1 above but only needs 1 reference system:
      !    Col. 1: Coord type
      !            21 for CORD2R
      !            22 for CORD2C
      !            23 for CORD2S
      !    Col. 2: CID = Coord system ID
      !    Col. 3: RID = Reference coord system ID for this CORD2R,C,S
      !        ! identify the GRIDs vs. the SPOINTs

      DO I=1,NCORD
        IF (CORD(I,1) == 11) THEN
          NCORD1R = NCORD1R + 1
          CORD1R_INDEX(NCORD1R) = I
        ELSE IF (CORD(I,1) == 12) THEN
          NCORD1C = NCORD1C + 1
          CORD1C_INDEX(NCORD1C) = I
        ELSE IF (CORD(I,1) == 13) THEN
          NCORD1S = NCORD1S + 1
          CORD1S_INDEX(NCORD1S) = I

        ELSE IF (CORD(I,1) == 21) THEN
          NCORD2R = NCORD2R + 1
          CORD2R_INDEX(NCORD2R) = I
        ELSE IF (CORD(I,1) == 22) THEN
          NCORD2C = NCORD2C + 1
          CORD2C_INDEX(NCORD2C) = I
        ELSE IF (CORD(I,1) == 23) THEN
          NCORD2S = NCORD2S + 1
          CORD2S_INDEX(NCORD2S) = I
        ENDIF
      ENDDO


      ! CORD1C : (1701, 17, 6)
      ! CORD1R : (1801, 18, 5)
      ! CORD1S : (1901, 19, 7)
      !IF(NCORD1R > 0) THEN
      !    CALL WRITE_OP2_CORD1(ITABLE, NCORD1R, CORD1R_INDEX, 1, 1801, 18, 5)
      !ENDIF
      !IF(NCORD1C > 0) THEN
      !    CALL WRITE_OP2_CORD1(ITABLE, NCORD1C, CORD1C_INDEX, 2, 1701, 17, 6)
      !ENDIF
      !IF(NCORD1S > 0) THEN
      !    CALL WRITE_OP2_CORD1(ITABLE, NCORD1S, CORD1S_INDEX, 3, 1901, 19, 7)
      !ENDIF

      ! CORD2C : (2001, 20, 9)
      ! CORD2R : (2101, 21, 8)
      ! CORD2S : (2201, 22, 10)
      IF(NCORD2R > 0) THEN
          CALL WRITE_OP2_CORD2(ITABLE, NCORD2R, CORD2R_INDEX, 1, 2101, 21, 8)
      ENDIF
      IF(NCORD2C > 0) THEN
          CALL WRITE_OP2_CORD2(ITABLE, NCORD2C, CORD2C_INDEX, 2, 2001, 20, 9)
      ENDIF
      IF(NCORD2S > 0) THEN
          CALL WRITE_OP2_CORD2(ITABLE, NCORD2S, CORD2S_INDEX, 3, 2201, 22, 10)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM1_CORD

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_CORD1(ITABLE, NCORD1, CORD1_INDEX, RCS_INT, CODEA, CODEB, CODEC)
      !! TODO: this is probably very wrong...
      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  OP2
      USE SCONTR, ONLY                :  NCORD
      USE MODEL_STUF, ONLY            :  CORD, RCORD
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT)                 :: ITABLE                  ! the subtable counter
      INTEGER(LONG), INTENT(IN)                    :: NCORD1                  ! the number of CORD1x cards
      INTEGER(LONG), DIMENSION(NCORD), INTENT(IN)  :: CORD1_INDEX             ! oversized, but mehhh
      INTEGER(LONG), INTENT(IN)                    :: RCS_INT                 ! R, C, or S
      INTEGER(LONG), INTENT(IN)                    :: CODEA, CODEB, CODEC     ! the CORD2x codes
      INTEGER(LONG)                                :: I,J                     ! counter

      INTEGER(LONG)                                :: NVALUES                 ! helper flag
      INTEGER(LONG), PARAMETER                     :: NUM_WIDE = 6            ! helper flag
      INTEGER(LONG), PARAMETER                     :: COORD_INT = 1           ! CORD1x

      ! if there are NCORD1Rs
      NVALUES = NUM_WIDE * NCORD1
!          Col. 1: Coord type
!                  11 for CORD1R
!                  12 for CORD1C
!                  13 for CORD1S
!          Col. 2: CID = Coord system ID
!          Col. 3: RID = Reference coord system ID for Grid A
!          Col. 4: RID = Reference coord system ID for Grid B
!          Col. 5: RID = Reference coord system ID for Grid C

      ! [cid, coord_rcs_int, coord_int, g1, g2, g3]
      ! TODO: is this g1 or g1.rid ???
      WRITE(OP2) CODEA, CODEB, CODEC,                                                      &
                  (CORD(CORD1_INDEX(J),2), RCS_INT, COORD_INT,                             &
                  CORD(CORD1_INDEX(J),3), CORD(CORD1_INDEX(J),3), CORD(CORD1_INDEX(J),3),  &
                  J=1,NCORD1)
      CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      END SUBROUTINE WRITE_OP2_CORD1


!===================================================================================================================================
      SUBROUTINE WRITE_OP2_CORD2(ITABLE, NCORD2, CORD2_INDEX, RCS_INT, CODEA, CODEB, CODEC)
      ! needs a docstring
      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  OP2
      USE SCONTR, ONLY                :  NCORD
      USE MODEL_STUF, ONLY            :  CORD, RCORD
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT)                 :: ITABLE                 ! the subtable counter
      INTEGER(LONG), INTENT(IN)                    :: NCORD2                 ! the number of CORD2x cards
      INTEGER(LONG), DIMENSION(NCORD), INTENT(IN)  :: CORD2_INDEX            ! oversized, but mehhh
      INTEGER(LONG), INTENT(IN)                    :: RCS_INT                ! R, C, or S
      INTEGER(LONG), INTENT(IN)                    :: CODEA, CODEB, CODEC    ! the CORD2x codes
      INTEGER(LONG)                                :: I                      ! counter
      INTEGER(LONG)                                :: NVALUES                ! helper flag
      INTEGER(LONG), PARAMETER                     :: NUM_WIDE = 13          ! helper flag
      INTEGER(LONG), PARAMETER                     :: COORD_INT = 2          ! CORD2x

      ! if there are CORD2Rs
      NVALUES = NUM_WIDE * NCORD2
      WRITE(OP2) NVALUES + 3
      ! [cid, coord_rcs_int, coord_int, rid, e1, e2, e3]
      !
      ! Col. 2: CID = Coord system ID
      ! Col. 3: RID = Reference coord system ID for this CORD2R,C,S
      ! Col.  1: X coord of Pt A in RID coord system
      ! Col.  2: Y coord of Pt A in RID coord system
      ! Col.  3: Z coord of Pt A in RID coord system
      ! Col.  4: X coord of Pt B in RID coord system
      ! Col.  5: Y coord of Pt B in RID coord system
      ! Col.  6: Z coord of Pt B in RID coord system
      ! Col.  7: X coord of Pt C in RID coord system
      ! Col.  8: Y coord of Pt C in RID coord system
      ! Col.  9: Z coord of Pt C in RID coord system

      !! TODO: we want this to be written before CORD_MAP or whatever is called...
      WRITE(OP2) CODEA, CODEB, CODEC,                                                                                &
                  (CORD(CORD2_INDEX(I),2), RCS_INT, COORD_INT, CORD(CORD2_INDEX(I),3),                               &
                  REAL(RCORD(CORD2_INDEX(I),1),4), REAL(RCORD(CORD2_INDEX(I),2),4), REAL(RCORD(CORD2_INDEX(I),3),4), &
                  REAL(RCORD(CORD2_INDEX(I),4),4), REAL(RCORD(CORD2_INDEX(I),5),4), REAL(RCORD(CORD2_INDEX(I),6),4), &
                  REAL(RCORD(CORD2_INDEX(I),7),4), REAL(RCORD(CORD2_INDEX(I),8),4), REAL(RCORD(CORD2_INDEX(I),9),4), &
                  I=1,NCORD2)
      CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      END SUBROUTINE WRITE_OP2_CORD2

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2()
      ! Writes the GEOM2 / elements table
      ! for some reason, we need to write the SPOINTs here too
      !
      USE PENTIUM_II_KIND, ONLY       :  LONG, BYTE
      USE IOUNT1, ONLY                :  ERR, OP2

      ! GEOM2 - elements
      USE SCONTR, ONLY : NCTETRA4, NCTETRA10, NCPENTA6, NCPENTA15, NCHEXA8, NCHEXA20
      USE SCONTR, ONLY : NCQUAD4, NCQUAD4K, NCSHEAR, NCTRIA3, NCTRIA3K
      USE SCONTR, ONLY : NCROD, NCBAR, NCBEAM, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4
      USE SCONTR, ONLY : NCMASS, NCONM2

      USE SCONTR, ONLY     : NELE, NEDAT
      USE MODEL_STUF, ONLY : ETYPE, EOFF, EPNT

      IMPLICIT NONE
      INTEGER(LONG)                    :: ITABLE                        ! the subtable counter
      INTEGER(LONG)                    :: I,J                           ! counter
      CHARACTER(LEN=8*BYTE)            :: TABLE_NAME = "GEOM2"          ! the current op2 table name
      LOGICAL                          :: IS_GEOM2                      ! do we need to write the table

      INTEGER(LONG)  :: NCTETRA, NCPENTA, NCHEXA
      INTEGER(LONG), DIMENSION(NCELAS1, 6)  :: CELAS1
      INTEGER(LONG), DIMENSION(NCELAS2, 6)  :: CELAS2
      INTEGER(LONG), DIMENSION(NCELAS3, 4)  :: CELAS3
      INTEGER(LONG), DIMENSION(NCELAS4, 4)  :: CELAS4
      INTEGER(LONG), DIMENSION(NCROD, 4)    :: CROD
      INTEGER(LONG), DIMENSION(NCTRIA3, 4)  :: CTRIA3
      INTEGER(LONG), DIMENSION(NCQUAD4, 4)  :: CQUAD4
      INTEGER(LONG), DIMENSION(NCSHEAR, 4)  :: CSHEAR
!      INTEGER(LONG), DIMENSION(:, :), allocatable :: CTETRA
!      INTEGER(LONG), DIMENSION(:, :), allocatable :: CPENTA
!      INTEGER(LONG), DIMENSION(:, :), allocatable  :: CHEXA
      INTEGER(LONG)                         :: NCROD_ACTUAL = 0
      INTEGER(LONG)                         :: NCONROD_ACTUAL = 0
      INTEGER(LONG), DIMENSION(NCROD)       :: CROD_INDEX, CONROD_INDEX

      !CHARACTER( 8*BYTE), ALLOCATABLE :: ETYPE(:)    ! NELE  x 1 array of elem types
      !CHARACTER( 1*BYTE), ALLOCATABLE :: EOFF(:)     ! NELE  x 1 array of 'Y' for elem offsets or 'N' if not

      !INTEGER(LONG)     , ALLOCATABLE :: EDAT(:)     ! NEDAT x 1 array of elem connection data
      !INTEGER(LONG)     , ALLOCATABLE :: EPNT(:)     ! NELE  x 1 array of pointers to EDAT where data begins for an elem

      !DO I=1,NEDAT
      !ENDDO

      NCTETRA = NCTETRA4 + NCTETRA10
      NCPENTA = NCPENTA6 + NCPENTA15
      NCHEXA = NCHEXA8 + NCHEXA20

      IF((NCELAS1 > 0) .OR. (NCELAS2 > 0) .OR. (NCELAS3 > 0) .OR. (NCELAS4 > 0) .OR.        &
         (NCROD > 0) .OR. (NCBAR > 0) .OR. (NCBEAM > 0) .OR. (NCBUSH > 0) .OR.              &
         (NCSHEAR > 0) .OR. (NCTRIA3 > 0) .OR. (NCQUAD4 > 0) .OR.                           &
         (NCTETRA > 0) .OR. (NCPENTA > 0) .OR. (NCHEXA > 0)) THEN
        IS_GEOM2 = .TRUE.
      ENDIF


      IS_GEOM2 = .FALSE.
      IF (IS_GEOM2) THEN
!        ALLOCATE ( CTETRA(NCTETRA,12) )
!        ALLOCATE ( CPENTA(NCPENTA,17) )
!        ALLOCATE ( CHEXA(NCHEXA,22) )

        CALL GET_GEOM2(CELAS1, CELAS2, CELAS3, CELAS4, &
                       CROD,                           &
                       CTRIA3, CQUAD4, CSHEAR,         &
!                       CTETRA, CPENTA, CHEXA,          &
                       NCTETRA, NCPENTA, NCHEXA,       &
                       NCROD_ACTUAL, NCONROD_ACTUAL,   &
                       CROD_INDEX, CONROD_INDEX)

        CALL WRITE_OP2_GEOM_HEADER(TABLE_NAME, ITABLE)

        CALL WRITE_OP2_GEOM2_CROD(ITABLE, CROD,                 &
                                  NCROD_ACTUAL, NCONROD_ACTUAL, &
                                  CROD_INDEX, CONROD_INDEX)
        !CALL WRITE_OP2_GEOM2_CELAS1(ITABLE, CELAS1)
        !CALL WRITE_OP2_GEOM2_CELAS2(ITABLE, CELAS2)
        CALL WRITE_OP2_GEOM2_CELAS3(ITABLE, CELAS3)
        CALL WRITE_OP2_GEOM2_CELAS4(ITABLE, CELAS4)
        !CALL WRITE_OP2_GEOM2_CSHEAR(ITABLE, CSHEAR)
        !CALL WRITE_OP2_GEOM2_CTRIA3(ITABLE, CTRIA3)
        !CALL WRITE_OP2_GEOM2_CQUAD4(ITABLE, CQUAD4)
        !CALL WRITE_OP2_GEOM2_CTETRA(ITABLE, CTETRA, NCTETRA)

!       CORD2S_INDEX(1000000)   ! intentional crash
        CALL END_OP2_GEOM_TABLE(ITABLE)
!        DEALLOCATE (CTETRA)
!        DEALLOCATE (CPENTA)
!        DEALLOCATE (CHEXA)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2

!===================================================================================================================================
      SUBROUTINE GET_GEOM2(CELAS1, CELAS2, CELAS3, CELAS4, &
                     CROD,                                 &
                     CTRIA3, CQUAD4, CSHEAR,               &
!                     CTETRA, CPENTA, CHEXA,                &
                     NCTETRA, NCPENTA, NCHEXA,             &
                     NCROD_ACTUAL, NCONROD_ACTUAL,         &
                     CROD_INDEX, CONROD_INDEX)
      ! breaks EDAT into simpler to use arrays
      USE PENTIUM_II_KIND, ONLY       :  LONG, BYTE
      USE IOUNT1, ONLY                :  ERR, OP2
      USE SCONTR, ONLY     : NCTETRA4, NCTETRA10, NCPENTA6, NCPENTA15, NCHEXA8, NCHEXA20
      USE SCONTR, ONLY     : NCQUAD4, NCQUAD4K, NCSHEAR, NCTRIA3, NCTRIA3K
      USE SCONTR, ONLY     : NCROD, NCBAR, NCBEAM, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4
      USE SCONTR, ONLY     : NCMASS, NCONM2
      USE SCONTR, ONLY     : NELE, NEDAT
      USE SCONTR, ONLY : NPROD
      USE MODEL_STUF, ONLY : ETYPE, EOFF, EPNT, EDAT, PROD, RPROD

      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN)  :: NCTETRA ! number of elements
      INTEGER(LONG), INTENT(IN)  :: NCPENTA ! number of elements
      INTEGER(LONG), INTENT(IN)  :: NCHEXA  ! number of elements
      INTEGER(LONG) :: I, J                 ! loop counters

      INTEGER(LONG) :: ICELAS1 = 1         ! loop counter
      INTEGER(LONG) :: ICELAS2 = 1         ! loop counter
      INTEGER(LONG) :: ICELAS3 = 1         ! loop counter
      INTEGER(LONG) :: ICELAS4 = 1         ! loop counter

      INTEGER(LONG) :: ICROD = 1           ! loop counter
      INTEGER(LONG) :: ICBAR = 1           ! loop counter
      INTEGER(LONG) :: ICBEAM = 1          ! loop counter

      INTEGER(LONG) :: ICQUAD4 = 1         ! loop counter
      INTEGER(LONG) :: ICTRIA3 = 1         ! loop counter
      INTEGER(LONG) :: ICSHEAR = 1         ! loop counter

      INTEGER(LONG) :: ICTETRA = 1         ! loop counter
      INTEGER(LONG) :: ICPENTA = 1         ! loop counter
      INTEGER(LONG) :: ICHEXA = 1          ! loop counter
      INTEGER(LONG) :: EPNTK
      INTEGER(LONG) :: NCROD_ACTUAL, NCONROD_ACTUAL

      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCELAS1, 6)  :: CELAS1
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCELAS2, 6)  :: CELAS2
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCELAS3, 4)  :: CELAS3
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCELAS4, 4)  :: CELAS4
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCROD, 4)    :: CROD

      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCTRIA3, 4)  :: CTRIA3
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCQUAD4, 4)  :: CQUAD4
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCSHEAR, 4)  :: CSHEAR
      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCROD)       :: CROD_INDEX, CONROD_INDEX
!      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCTETRA, 12) :: CTETRA
!      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCPENTA, 17) :: CPENTA
!      INTEGER(LONG), INTENT(INOUT), DIMENSION(NCHEXA, 22)  :: CHEXA
      INTEGER(LONG) :: IPROD

 1    FORMAT(A)
      WRITE(ERR,1) "STARTING GET_GEOM2"
      FLUSH(ERR)
2     FORMAT("PROD: I=",i4, "; PID=", i8, "; MID=", i8, "; A=",f8.4,"; J=",f8.4)

 3    FORMAT(i8, "; ETYPE=", A, "; EOFF=",A, "; EPNT=", i8)
 4    FORMAT("  4 int fields: ", A, 4(" ", i4))
 6    FORMAT("  6 int fields: ", A, 6(" ", i4))
 8    FORMAT("  8 int fields: ", A, 8(" ", i4))
10    FORMAT(" 10 int fields: ", A, 10(" ", i4))
12    FORMAT(" 12 int fields: ", A, 12(" ", i4))
100   FORMAT(" ROD: ", A, " ", 3(A, i8, " "))
101   FORMAT(" ROD: ", A, i8)

      WRITE(1, 101) "NPROD=", NPROD
      DO I=1,NPROD
      !A,J,C,NMS
        WRITE(ERR,2) I, PROD(I,1), PROD(I,2), RPROD(I,1), RPROD(I,2)
        FLUSH(ERR)
      ENDDO
      WRITE(ERR,1) "FINISHED WRITING PROD"
      WRITE(ERR, 101) "NELE=", NELE

      DO I=1,NELE
        WRITE(ERR,3) I, ETYPE(I), EOFF(I), EPNT(I)
        EPNTK = EPNT(I)

        !EID  = EDAT(EPNTK)
        !PID  = EDAT(EPNTK+1)
        IF      (ETYPE(I)(1:5) .EQ. 'ELAS1') THEN
          WRITE(ERR,6) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,6)
          FLUSH(ERR)
          ! 6 fields
          DO J=1,6
            CELAS1(ICELAS1, J) = EDAT(EPNTK-1+J)
          ENDDO
          ICELAS1 = ICELAS1 + 1
        ELSE IF (ETYPE(I)(1:5) .EQ. 'ELAS2') THEN
          ! 6 fields
          ! TODO: ge, s
          DO J=1,6
            CELAS2(ICELAS2, J) = EDAT(EPNTK-1+J)
          ENDDO
          !WRITE(ERR,6) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,6)
          FLUSH(ERR)
          ICELAS2 = ICELAS2 + 1

        ELSE IF (ETYPE(I)(1:5) .EQ. 'ELAS3') THEN
          ! 4 fields
          DO J=1,4
            CELAS3(ICELAS3, J) = EDAT(EPNTK-1+J)
          ENDDO
          !WRITE(ERR,4) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,4)
          FLUSH(ERR)
          ICELAS3 = ICELAS3 + 1

        ELSE IF (ETYPE(I)(1:5) .EQ. 'ELAS4') THEN
          ! 4 fields
          DO J=1,4
            CELAS4(ICELAS4, J) = EDAT(EPNTK-1+J)
          ENDDO
          !WRITE(ERR,4) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,4)
          FLUSH(ERR)
          ICELAS4 = ICELAS4 + 1

        ELSE IF (ETYPE(I)(1:3) .EQ. 'ROD') THEN
          !WRITE(ERR,4) "CROD?", (EDAT(EPNTK-1+J), J=1,4)
          ! 4 fields
          CROD(ICROD, J) = EDAT(EPNTK)   ! eid
          CROD(ICROD, J) = EDAT(EPNTK+1) ! pid
          CROD(ICROD, J) = EDAT(EPNTK+2) ! ga
          CROD(ICROD, J) = EDAT(EPNTK+3) ! gb
          IPROD = EDAT(EPNTK+1) ! CROD(ICROD,2)

          ! PROD: I=   1; PID=      20; MID=       1; A=  1.0000; J=  0.0000
          ! PROD: I=   2; PID=     -23; MID=       1; A=  2.0000; J=  0.1000
          IF (PROD(IPROD,1) > 0) THEN
            ! pid > 0
            ! ncrod_actual starts at 0
            ! noncrod_actual starts at 0
            ! icrod starts at 1
            WRITE(ERR,4) "CROD", (EDAT(EPNTK-1+J), J=1,4)
            FLUSH(ERR)
            NCROD_ACTUAL = NCROD_ACTUAL + 1
            WRITE(ERR,100) "CROD", "ICROD=",ICROD,"NCROD_ACTUAL=",NCROD_ACTUAL
            FLUSH(ERR)
            CROD_INDEX(NCROD_ACTUAL) = ICROD
          ELSE
            WRITE(ERR,4) "CONROD", (EDAT(EPNTK-1+J), J=1,4)
            FLUSH(ERR)
            ! eid, pid, ga, gb
            NCONROD_ACTUAL = NCONROD_ACTUAL + 1
            WRITE(ERR,100) "CONROD", "ICROD=",ICROD,"NCONROD_ACTUAL=",NCONROD_ACTUAL
            FLUSH(ERR)
            CONROD_INDEX(NCONROD_ACTUAL) = ICROD
          ENDIF
          ICROD = ICROD + 1


        ELSE IF (ETYPE(I)(1:5) .EQ. 'SHEAR') THEN
          ! 6 fields
          ICSHEAR = ICSHEAR + 1
          DO J=1,6
            CSHEAR(ICSHEAR, J) = EDAT(EPNTK-1+J)
          ENDDO
          !WRITE(ERR,6) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,6)

!        ELSE IF (ETYPE(I)(1:5) .EQ. 'TETRA4') THEN
!          ! 6 fields
!          DO J=1,6
!            CTETRA(ICTETRA, J) = EDAT(EPNTK-1+J)
!          ENDDO
!          DO J=7,12
!            CTETRA(ICTETRA, J) = 0
!          ENDDO
!          WRITE(ERR,6) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,12)
!          ICTETRA = ICTETRA + 1
!        ELSE IF (ETYPE(I)(1:5) .EQ. 'TETRA10') THEN
!          ! 12 fields
!          DO J=1,12
!            CTETRA(ICTETRA, J) = EDAT(EPNTK-1+J)
!          ENDDO
!          WRITE(ERR,12) ETYPE(I)(1:5), (EDAT(EPNTK-1+J), J=1,12)
!          ICTETRA = ICTETRA + 1
        ENDIF
        !EPNT(I)
      ENDDO
      WRITE(ERR,1) "FINISHED GET_GEOM2"
      FLUSH(ERR)
      END SUBROUTINE GET_GEOM2

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CROD(ITABLE, CROD,                 &
                                      NCROD_ACTUAL, NCONROD_ACTUAL, &
                                      CROD_INDEX, CONROD_INDEX)
      ! Writes the CROD and CONROD subtables
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  ERR,OP2
      USE SCONTR, ONLY : NCROD
      USE MODEL_STUF, ONLY : PROD, RPROD
      IMPLICIT NONE
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE                        ! the op2 subtable counter
      INTEGER(LONG), INTENT(IN), DIMENSION(NCROD, 4) :: CROD                          ! stores eid, ipid, n1, n2 from EDAT
      INTEGER(LONG)                                  :: NUM_WIDE                      ! helper flag
      INTEGER(LONG)                                  :: NVALUES                       ! helper flag
      INTEGER(LONG)                                  :: I,J,IPID,MID                  ! counters
      INTEGER(LONG)                                  :: NCROD_ACTUAL, NCONROD_ACTUAL  ! counters
      INTEGER(LONG), INTENT(IN), DIMENSION(NCROD)    :: CROD_INDEX, CONROD_INDEX      ! flags where the CROD and CONRODs are in CROD
      LOGICAL                                        :: WRITE_OP2 = .FALSE.
      LOGICAL                                        :: WRITE_ERR = .TRUE.
      FLUSH(ERR)
 1    FORMAT(A)
 2    FORMAT(" **WRITE_OP2_GEOM2_CROD: ", 4(i8, " "))
      WRITE(ERR,1) "STARTING WRITE_OP2_GEOM_CROD"
      FLUSH(ERR)
      DO I=1,NCROD
        WRITE(ERR,2) CROD(I,1), CROD(I,2), CROD(I,3), CROD(I,4)
        FLUSH(ERR)
      ENDDO
      IF (NCROD_ACTUAL > 0) THEN
        NUM_WIDE = 4
        NVALUES = NUM_WIDE * NCROD_ACTUAL

        ! ROD  4 words: (all read by call to subr ELEPO from subr BD_ROD1)
        !   1) Elem ID
        !   2) Prop ID index
        !   3) Grid A
        !   4) Grid B
        IF (WRITE_ERR) THEN
 4        FORMAT(" **WRITE_OP2_GEOM2_CROD-B4: ", 5(i8, " "))
          FLUSH(ERR)
          DO I=1,NCROD_ACTUAL
            J = CROD_INDEX(I)
            WRITE(ERR,4) J, CROD(J,1), CROD(J,2), CROD(J,3), CROD(J,4)
            FLUSH(ERR)
          ENDDO

 5        FORMAT(" **WRITE_OP2_GEOM2_CROD-B5: ", 5(i8, " "))
          DO I=1,NCROD_ACTUAL
            J = CROD_INDEX(I)
            IPID = CROD(J,2)
            WRITE(ERR,5) J, CROD(J,1), PROD(IPID,1), CROD(J,3), CROD(J,4)
            FLUSH(ERR)
          ENDDO
        ENDIF

        IF(WRITE_OP2) THEN
          WRITE(OP2) NVALUES + 3
          ! eid, pid, n1, n2
          WRITE(OP2) 3001, 30, 48, (                   &
                            CROD(CROD_INDEX(I),1),     &
                       PROD(CROD(CROD_INDEX(I),2),1),  &
                            CROD(CROD_INDEX(I),3),     &
                            CROD(CROD_INDEX(I),4), I=1,NCROD_ACTUAL)
          CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
          FLUSH(OP2)
        ENDIF
      ENDIF

      IF (NCONROD_ACTUAL > 0) THEN
        NUM_WIDE = 8
        NVALUES = NUM_WIDE * NCONROD_ACTUAL

 50       FORMAT(" **WRITE_OP2_GEOM2_CONROD: ", 4(i4, " "), 4(f8.4, " "))
        !IF (WRITE_ERR) THEN
        !  DO I=1,NCONROD_ACTUAL
        !    ! eid, ipid
        !    J = CONROD_INDEX(I)
        !    IPID = CROD(J, 2)
        !    MID = -PROD(IPID, 2)
        !    ! (eid, n1, n2, mid, a, j, c, nsm) = out
        !    WRITE(ERR,50) CROD(J,1), MID, CROD(J,3), CROD(J,4),                        &
        !                  RPROD(IPID,1), RPROD(IPID,2), RPROD(IPID,3), RPROD(IPID,4)
        !    FLUSH(ERR)
        !  ENDDO
        !ENDIF
        !WRITE(OP2) NVALUES + 3
        !
        !! (eid, n1, n2, mid, a, j, c, nsm) = out
        !CALL GET_PCONROD_INDEXS(CONROD_INDEX, NCONROD_ACTUAL)
        !! mid, a, j, c, nsm
        !WRITE(OP2) 1601, 16, 47, (CROD(CONROD_INDEX(I),1), CROD(CONROD_INDEX(I),3), CROD(CONROD_INDEX(I),4), &
        !                        PROD(PCONROD_INDEX(I),1),                                                  &
        !                        RPROD(PCONROD_INDEX(I),1), PROD(PCONROD_INDEX(I),2),                       &
        !                        RPROD(PCONROD_INDEX(I),3), PROD(PCONROD_INDEX(I),4),                       &
        !                        I=1,NCONROD_ACTUAL)
        !CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CROD


!===================================================================================================================================
      !SUBROUTINE GET_PCONROD_INDEXS(CONROD_INDEX, NCONROD_ACTUAL)
      !USE PENTIUM_II_KIND, ONLY         :  LONG
      !INTEGER(LONG), INTENT(INOUT), DIMENSION(NCONROD_ACTUAL) :: PCONROD_INDEX
      !INTEGER(LONG)                                           :: NCONROD_ACTUAL
      !END SUBROUTINE GET_PCONROD_INDEXS

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CELAS1(ITABLE, CELAS1)
      ! writes the CELAS1 subtable
      ! shockingly straightforward
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      USE SCONTR, ONLY : NCELAS1
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCELAS1, 6)           :: CELAS1
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NVALUES                      ! helper flag
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 6                 ! helper flag
      INTEGER(LONG)                                  :: I,J                          ! counters
      IF (NCELAS1 > 0) THEN
        NVALUES = NUM_WIDE * NCELAS1
        WRITE(OP2) NVALUES + 3

        ! ELAS1  6 words: (all read by call to subr ELEPO from subr BD_CELAS1)
        !   1) Elem ID
        !   2) Prop ID
        !   3) Grid A
        !   4) Grid B
        !   5) Components at Grid A
        !   6) Components at Grid B
        !(eid, pid, g1, g2, c1, c2) = out
        WRITE(OP2) 601, 6, 73, ((CELAS1(I,J), J=1,6), I=1,NCELAS1)
        CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CELAS1

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CELAS2(ITABLE, CELAS2)
      ! writes the CELAS2 subtable
      ! TODO needs K, GE, S
      USE PENTIUM_II_KIND, ONLY         :  LONG, DOUBLE
      USE IOUNT1, ONLY                  :  ERR, OP2
      USE SCONTR, ONLY : NCELAS2
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCELAS2, 6)           :: CELAS2
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NVALUES                      ! helper flag
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 8                 ! helper flag
      INTEGER(LONG)                                  :: I,J,IPID                     ! counters
      REAL(DOUBLE)                                   :: K = 1000.0
      REAL(DOUBLE), PARAMETER                        :: GE = 0.0
      REAL(DOUBLE), PARAMETER                        :: S = 0.0
      LOGICAL                                        :: WRITE_OP2 = .FALSE.
      LOGICAL                                        :: WRITE_ERR = .TRUE.
      IF (NCELAS2 > 0) THEN
        NVALUES = NUM_WIDE * NCELAS2

        ! ELAS2  6 words: (all read by call to subr ELEPO from subr BD_CELAS2)
        !   1) Elem ID
        !   2) Prop ID which is set to -EID since real props are on the CELAS2 entry
        !   3) Grid A
        !   4) Grid B
        !   5) Components at Grid A
        !   6) Components at Grid B
        IF (WRITE_ERR) THEN
 6        FORMAT(" **WRITE_OP2_GEOM2_CELAS2: ", i4, " ", f8.1, " ", 4(i4, " "))
          DO I=1,NCELAS2
            IPID = CELAS2(I,2)
            !K = RPELAS(IPID,1)
            !GE = RPELAS(IPID,1)
            !S = RPELAS(IPID,1)
            WRITE(ERR,6) CELAS2(I,1), K, CELAS2(I,3),                    &
                         CELAS2(I,4), CELAS2(I,5), CELAS2(I,6)
            FLUSH(ERR)
          ENDDO
        ENDIF

        ! [eid, k, g1, g2, c1, c2], ge, s
        IF (WRITE_OP2) THEN
          WRITE(OP2) NVALUES + 3
          WRITE(OP2) 701, 7, 74, (CELAS2(I,1), CELAS2(I,2), CELAS2(I,3), &
                                  CELAS2(I,4), CELAS2(I,5), CELAS2(I,6), &
                                  REAL(GE,4), REAL(S,4), I=1,NCELAS2)
          CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
        ENDIF
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CELAS2

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CELAS3(ITABLE, CELAS3)
      ! writes the CELAS3 subtable
      ! shockingly straightforward
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      USE SCONTR, ONLY : NCELAS3
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCELAS3, 4)           :: CELAS3
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NVALUES                      ! helper flag
      INTEGER(LONG)                                  :: I,J                          ! counters
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 4                 ! helper flag
      IF (NCELAS3 > 0) THEN
        NVALUES = NUM_WIDE * NCELAS3
        WRITE(OP2) NVALUES + 3

        ! ELAS3  4 words: (all read by call to subr ELEPO from subr BD_CELAS3)
        !   1) Elem ID
        !   2) Prop ID
        !   3) Scalar point A
        !   4) Scalar point B
        ! (eid, pid, s1, s2) = out
        WRITE(OP2) 801, 8, 75, ((CELAS3(I,J), J=1,4), I=1,NCELAS3)
        CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CELAS3

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CELAS4(ITABLE, CELAS4)
      ! writes the CELAS4 subtable
      ! TODO needs K, not pid
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      USE SCONTR, ONLY : NCELAS4
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCELAS4, 4)           :: CELAS4
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NVALUES                      ! helper flag
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 4                 ! helper flag
      INTEGER(LONG)                                  :: I,J                          ! counters
      IF (NCELAS4 > 0) THEN
        NVALUES = NUM_WIDE * NCELAS4
        WRITE(OP2) NVALUES + 3

        ! ELAS4  4 words: (all read by call to subr ELEPO from subr BD_CELAS4)
        !   1) Elem ID
        !   2) Prop ID which is set to -EID since real props are on the CELAS2 entry
        !   3) Scalar point A
        !   4) Scalar point B
        ! (eid, k, s1, s2) = out
        WRITE(OP2) 901, 9, 76, (CELAS4(I,1), CELAS4(I,2), CELAS4(I,3), CELAS4(I,4), I=1,NCELAS4)
        CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CELAS4

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CSHEAR(ITABLE, CSHEAR)
      ! writes the CSHEAR subtable
      ! shockingly straightforward
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      USE SCONTR, ONLY : NCSHEAR
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCSHEAR, 6)           :: CSHEAR
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NVALUES                      ! helper flag
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 6                 ! helper flag
      INTEGER(LONG)                                  :: I,J                          ! counters
      IF (NCSHEAR > 0) THEN
        NVALUES = NUM_WIDE * NCSHEAR
        WRITE(OP2) NVALUES + 3
        ! eid, pid, n1, n2, n3, n4
        WRITE(OP2) 3101, 31, 61, ((CSHEAR(I,J), J=1,6), I=1,NCSHEAR)
        CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CSHEAR!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM2_CTETRA(ITABLE, CTETRA, NCTETRA)
      ! writes the CTETRA subtable
      ! CTETRA4 / CTETRA10 are the same kind of table despite being stored differently in MYSTRAn
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NCTETRA, 12)          :: CTETRA
      INTEGER(LONG),INTENT(INOUT)                    :: ITABLE
      INTEGER(LONG)                                  :: NCTETRA, NVALUES             ! helper flag
      INTEGER(LONG), PARAMETER                       :: NUM_WIDE = 12                ! helper flag
      INTEGER(LONG)                                  :: I,J                          ! counters
      IF (NCTETRA > 0) THEN
        NVALUES = NUM_WIDE * NCTETRA
        WRITE(OP2) NVALUES + 3

        ! TETRA4   6 words: (read by 1 or more calls to subr ELEPO from subr BD_TETRA)
        !   1) Elem ID
        !   2) Prop ID
        !   3) etc, Grids 1-4

        ! TETRA10  2 words: (read by 1 or more calls to subr ELEPO from subr BD_TETRA)
        !   1) Elem ID
        !   2) Prop ID
        !   3) etc, Grids 1-10
        WRITE(OP2) 5508, 55, 217, ((CTETRA(I,J), J=1,12), I=1,NCTETRA)
        !WRITE(OP2) 5508, 55, 217, &
        !           ((TETRA4(I,J),  J=1,6), 0, 0, 0, 0, 0, 0, I=1,NCTETRA4), &
        !           ((TETRA10(I,J), J=1,12), I=1,NCTETRA10)
          CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM2_CTETRA
!===================================================================================================================================

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_HEADER(POST)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  OP2

      integer, intent(in)  :: POST
      character(len=28)    :: TAPE_CODE
      character(len=8)     :: NASTRAN_VERSION
      IF (POST /= -2) THEN
        !_write_markers(op2, op2_ascii, [3, 0, 7])
        WRITE(OP2) 3
        WRITE(OP2) 3, 24, 2021-2000  ! date
        !WRITE(OP2) 0
        WRITE(OP2) 7
        TAPE_CODE = 'NASTRAN FORT TAPE ID CODE - '
!       write(OP2) 1
!       write(OP2) 7
        WRITE(OP2) TAPE_CODE
!       fop2.write(pack(endian + b'7i 28s i', *[4, 1, 4,
!                                               4, 7, 4,
!                                               28, tape_code, 28]))

!       nastran_version = b'NX8.5   ' if obj.is_nx else b'XXXXXXXX'
        NASTRAN_VERSION = 'XXXXXXXX'
        WRITE(OP2) 2
        WRITE(OP2) NASTRAN_VERSION
!       write(OP2) 'XXXXXXXX'
!       fop2.write(pack(endian + b'4i 8s i', *[4, 2, 4,
!                                              #4, 2, 4,
!                                              #4, 1, 4,
!                                              #4, 8, 4,
!                                              8, nastran_version, 8]))
        WRITE(OP2) -1
        WRITE(OP2) 0
!       fop2.write(pack(endian + b'6i', *[4, -1, 4,
!                                         4, 0, 4,]))
      ELSE
!       write_markers(fop2_ascii, [2, 4])
        WRITE(OP2) 2
        WRITE(OP2) 4
      ENDIF
      END SUBROUTINE WRITE_OP2_HEADER

!===================================================================================================================================
      SUBROUTINE  END_OP2_TABLE(ITABLE)
      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  ERR, OP2
      IMPLICIT NONE
      INTEGER(LONG) :: ITABLE
      WRITE(ERR,9114) ITABLE

 9114 FORMAT(" *DEBUG:       END_OP2_TABLE; ITABLE=", I8)
      WRITE(OP2) ITABLE
      WRITE(OP2) 1
      WRITE(OP2) 0
      WRITE(OP2) 0
      END SUBROUTINE END_OP2_TABLE

!===================================================================================================================================
      SUBROUTINE  END_OP2_TABLES()
      USE IOUNT1, ONLY                :  ERR, OP2, OP2FIL
      IMPLICIT NONE
      LOGICAL                         :: FILE_OPND
 9115 FORMAT(" *DEBUG:       END_OP2_TABLES", A)
      WRITE(ERR,9115) " "
      INQUIRE ( FILE=OP2FIL, OPENED=FILE_OPND )
      IF (FILE_OPND) THEN
        WRITE(OP2) 0
      ENDIF
      END SUBROUTINE END_OP2_TABLES

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM()
      USE PENTIUM_II_KIND, ONLY       :  LONG, BYTE
      USE IOUNT1, ONLY                :  ERR, OP2

      ! GEOM1 - GRID/COORDs
      USE SCONTR, ONLY                :  NGRID
      USE MODEL_STUF, ONLY            :  GRID_ID, GRID

      ! GEOM2 - elements
      !USE MODEL_STUF, ONLY : ELAS1, ELAS2, ELAS3, ELAS4, ROD, TETRA4, TETRA10
      USE SCONTR, ONLY : NCTETRA4, NCTETRA10, NCPENTA6, NCPENTA15, NCHEXA8, NCHEXA20
      USE SCONTR, ONLY : NCQUAD4, NCQUAD4K, NCSHEAR, NCTRIA3, NCTRIA3K
      USE SCONTR, ONLY : NCROD, NCBAR, NCBEAM, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4
      USE SCONTR, ONLY : NCMASS, NCONM2

      USE SCONTR, ONLY     : NELE, NEDAT
      USE MODEL_STUF, ONLY : ETYPE, EOFF, EPNT

      IMPLICIT NONE
      INTEGER(LONG)                    :: ITABLE                        ! the subtable counter
      INTEGER(LONG)                    :: I,J                           ! counter
      CHARACTER(LEN=8*BYTE)            :: TABLE_NAME                    ! the current op2 table name
      INTEGER(LONG)                    :: NGRID_ACTUAL = 0              ! number of GRID;    corrects for oversizing
      INTEGER(LONG)                    :: NSPOINT_ACTUAL = 0            ! number of SPOINTs; corrects for oversizing
      INTEGER(LONG), DIMENSION(NGRID)  :: GRID_INDEX, SPOINT_INDEX      ! oversized, but mehhh
      LOGICAL                          :: IS_GEOM1 = .FALSE.            ! do we need to write the table
      LOGICAL                          :: IS_GEOM3 = .FALSE.            ! do we need to write the table
      LOGICAL                          :: IS_GEOM4 = .FALSE.            ! do we need to write the table

      !---------------------------------------------------------------------------------------------
      IF(NGRID > 0) THEN
        IS_GEOM1 = .TRUE.
      ENDIF

! Each row of GRID is for one grid point and contains:
!              (1) Grid point number    in col 1
!              (2) Input coord system   in col 2
!              (3) Global coord system  in col 3
!              (4) Permanent SPC's      in col 4
!              (5) Line break indicator in col 5 (put this many line breaks in F06 after this grid number)
!              (6) Num of comps         in col 6 (1 for SPOINT, 6 for actual grid)
!
!            The array is sorted in the following order:
!              (1) After the B.D. deck is read it is in GRID input order
!              (2) After subr GRID_PROC it is in grid point numerical order
!
! Each row of RGRID is for one grid point and contains:
!          After Bulk Data has been read: the 3 coords of the grid in the coord sys defined in col 2 of array GRID for this G.P.
!          After subr GRID_PROC has run : the 3 coords of the grid in the basic (0) coord sys of the model
 1    FORMAT("****DEBUG:   WRITE_OP2 GEOM ngrid",i4)
      IF (IS_GEOM1) THEN
        WRITE(ERR,1) NGRID

        TABLE_NAME = "GEOM1"
        CALL WRITE_OP2_GEOM_HEADER(TABLE_NAME, ITABLE)

        ! identify the GRIDs vs. the SPOINTs
        DO I=1,NGRID
          IF (GRID(I,6) == 1) THEN
            NSPOINT_ACTUAL = NSPOINT_ACTUAL + 1
            SPOINT_INDEX(NSPOINT_ACTUAL) = I
          ELSE
            NGRID_ACTUAL = NGRID_ACTUAL + 1
            GRID_INDEX(NGRID_ACTUAL) = I
          ENDIF
        ENDDO

        CALL WRITE_OP2_GEOM1_GRID(ITABLE, NGRID_ACTUAL, GRID_INDEX)
        CALL WRITE_OP2_GEOM1_CORD(ITABLE)

        !DO I=1,NGRID
          ! is GRID_ID(I) the same as GRID(I,1)???
          !
          !(nid, cp, x1, x2, x3, cd, ps, seid)
          !                        nid        cp         x                   y                   z
        !  WRITE(ERR,2) GRID_ID(I), GRID(I,1), GRID(I,2), REAL(RGRID(I,1),4), REAL(RGRID(I,2),4), REAL(RGRID(I,3),4), &
          !            cd         ps         ndof
        !               GRID(I,3), GRID(I,4), GRID(I,6)
        !ENDDO
        CALL END_OP2_GEOM_TABLE(ITABLE)
      ENDIF

      CALL WRITE_OP2_GEOM2()
      CALL WRITE_OP2_GEOM_EPT()
      !CALL WRITE_OP2_GEOM_MPT()

      END SUBROUTINE WRITE_OP2_GEOM

!===================================================================================================================================
      SUBROUTINE GET_PROD_INDEX(PID, IPROD)
      USE PENTIUM_II_KIND, ONLY   :  LONG
      USE SCONTR, ONLY            :  NPROD
      USE MODEL_STUF, ONLY        :  PROD
      USE IOUNT1, ONLY            :  ERR

      INTEGER(LONG), INTENT(IN)    :: PID
      INTEGER(LONG), INTENT(INOUT) :: IPROD
      INTEGER(LONG)                :: I
      WRITE(ERR,1) PID,IPROD
      DO I=1,NPROD
        IF (PROD(I,1) == PID) THEN
          IPROD = I
          WRITE(ERR,2) PID
          EXIT  ! break statement
        ENDIF
      ENDDO
      WRITE(ERR,3) PID,IPROD,PROD(IPROD,1)
 1    FORMAT("GET_PROD_INDEX START:  PID=",i4,"; IPROD=",i4)
 2    FORMAT("GET_PROD_INDEX MID:  *********FOUND PID=",i4)
 3    FORMAT("GET_PROD_INDEX END:  PID=",i4,"; IPROD=",i4,"; PROD(I,1)=",i4)
      END SUBROUTINE GET_PROD_INDEX

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM_EPT()
      ! writes the element property table (EPT)
      USE PENTIUM_II_KIND, ONLY       :  LONG, BYTE
      USE IOUNT1, ONLY                :  ERR, OP2
      USE SCONTR, ONLY : NPELAS, NPROD, NPBUSH, NPCOMP, NPSHEAR, NPSOLID ! NPSHELL,
      USE MODEL_STUF, ONLY : PELAS, PROD, PBAR, PBEAM, PSHEAR, PCOMP, PSOLID ! PSHELL,
      USE MODEL_STUF, ONLY : RPELAS, RPROD, RPROD, RPBAR, RPBEAM, RPSHEAR, RPCOMP ! RPSHELL, RPSOLID

      IMPLICIT NONE
      !INTEGER(LONG), INTENT(IN)   :: NCTETRA ! number of elements
      INTEGER(LONG)                :: I, J                        ! loop counter
      INTEGER(LONG)                :: ITABLE                      ! the subtable id
      INTEGER(LONG)                :: NUM_WIDE                    ! the width of a single card
      INTEGER(LONG)                :: NVALUES                     ! the number of words in the block
      CHARACTER(8*BYTE), PARAMETER :: TABLE_NAME = "EPT"          ! the table_name we're writing
      INTEGER(LONG)                :: PID, NPROD_ACTUAL = 0
      INTEGER(LONG),DIMENSION(NPROD) :: PROD_INDEX

      LOGICAL :: IS_EPT


      !  PELAS  = Array of integer data from PELAS Bulk Data entries. Each row is for one PELAS entry read in B.D. and contains:
      !             ( 1) Col  1: Property ID
      !
      !  RPELAS = Array of real data from PELAS Bulk Data entries. Each row is for one PELAS entry read in B.D. and contains:
      !             ( 1) Col  1: Spring rate, K
      !             ( 2) Col  2: Damping coefficient, GE
      !             ( 3) Col  3: Stress recovery coefficient, S


      !  PROD   = Array of integer data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
      !             ( 1) Col  1: Property ID
      !             ( 2) Col  2: Material ID
      !
      !  RPROD  = Array of real data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
      !             ( 1) Col  1: Cross sectional area, A
      !             ( 2) Col  2: Torsion constant, J
      !             ( 3) Col  3: Torsion stress recov. coeff, C
      !             ( 4) Col  4: Non structural mass, NSM

      IS_EPT = ( (NPROD > 0) )

      IF (IS_EPT) THEN
        CALL WRITE_OP2_GEOM_HEADER(TABLE_NAME, ITABLE)

        IF (NPROD > 0) THEN
          NUM_WIDE = 6
          ! (pid, mid, a, j, c, nsm)
          !  PROD= Array of integer data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
          !    (1) Col  1: Property ID
          !    (2) Col  2: Material ID
          !  RPROD = Array of real data from PROD  Bulk Data entries. Each row is for one PROD  entry read in B.D. and contains:
          !    (1) Col  1: Cross sectional area, A
          !    (2) Col  2: Torsion constant, J
          !    (3) Col  3: Torsion stress recov. coeff, C
          !    (4) Col  4: Non structural mass, NSM
          DO I=1,NPROD
            PID = PROD(I,1)
            IF (PID /= 0) THEN
              NPROD_ACTUAL = NPROD_ACTUAL + 1
              PROD_INDEX(I) = NPROD_ACTUAL
            ENDIF
          ENDDO
          NVALUES = NUM_WIDE * NPROD
          WRITE(OP2) NVALUES + 3
          WRITE(OP2) 902, 9, 29, (PROD(PROD_INDEX(I),1),PROD(PROD_INDEX(I),2),               &
                                 (REAL(RPROD(PROD_INDEX(I),J),4), J=1,4), I=1,NPROD_ACTUAL)
          CALL WRITE_OP2_SUBTABLE_INCREMENT(ITABLE)
        ENDIF
        CALL END_OP2_GEOM_TABLE(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM_EPT

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM_MPT()
      ! writes the material property table (MPT)
      USE PENTIUM_II_KIND, ONLY       :  LONG, BYTE, DOUBLE
      USE IOUNT1, ONLY                :  ERR, OP2
      USE SCONTR, ONLY : NMATL
      USE MODEL_STUF, ONLY : MATL, RMATL
      IMPLICIT NONE
      INTEGER(LONG), DIMENSION(NMATL) :: MAT1_INDEX, MAT2_INDEX, MAT8_INDEX, MAT9_INDEX
      REAL(DOUBLE), DIMENSION(NMATL, 10) :: RMAT1, RMAT2, RMAT8, RMAT9

!  MATL   = Array of integer data from Material Bulk Data entries. Each row is for one Material entry read in B.D. and contains:
!             ( 1) Col  1: Material ID
!             ( 2) Col  2: Material type: 1, in MATL array for MAT1 entry (isotropic)
!                                         2, in MATL array for MAT2 entry (anisotropic, 2D elements)
!                                         8, in MATL array for MAT8 entry (orthotropic)
!                                         9, in MATL array for MAT9 entry (anisotropic, 3D elements)

!  RMATL  = Array of real data from material Bulk Data entries. Each row is for one material entry read in B.D. and contains:
!
!           MAT1 (isotropic):
!           ----------------
!                ( 1) Col  1: E     Young's modulus, E
!                ( 2) Col  2: G     Shear modulus, G
!                ( 3) Col  3: NU    Poisson's ratio, NU
!                ( 4) Col  4: RHO   Mass density, RHO
!                ( 5) Col  5: ALPHA Thermal expansion coeff, ALPHA
!                ( 6) Col  6: TREF  Reference temperature, TREF
!                ( 7) Col  7: GE    Damping coefficient, GE
!                ( 8) Col  8: ST    Tension limit, ST
!                ( 9) Col  9: SC    Compression limit, SC
!                (10) Col 10: SS    Shear limit, SS

!           MAT2 (2D anisotropic):
!           ---------------------
!                ( 1) Col  1: G11   11 term in stress/strain matrix
!                ( 2) Col  2: G12   12 term in stress/strain matrix
!                ( 3) Col  3: G13   13 term in stress/strain matrix
!                ( 4) Col  4: G22   22 term in stress/strain matrix
!                ( 5) Col  5: G23   23 term in stress/strain matrix
!                ( 6) Col  6: G33   33 term in stress/strain matrix
!                ( 7) Col  7: RHO   Mass density
!                ( 8) Col  8: A1    Thermal expansion coeff in material x direction
!                ( 9) Col  9: A2    Thermal expansion coeff in material y direction
!                (10) Col 10: A3    Thermal expansion coeff in shear
!                (11) Col 11: TREF  Reference temperature for thermal expansion
!                (12) Col 12: GE    Damping coeff
!                (13) Col 13: ST    Tension limit, ST
!                (14) Col 14: SC    Compression limit, SC
!                (15) Col 15: SS    Shear limit, SS

!           MAT8 (orthotropic):
!           ------------------
!                ( 1) Col  1: E1    Modulus of elasticity in the longitudinal direction
!                ( 2) Col  2: E2    Modulus of elasticity in the lateral direction
!                ( 3) Col  3: NU12  Poisson's ratio
!                ( 4) Col  4: G12   In-plane shear modulus
!                ( 5) Col  5: G1Z   Transverse shear modulus in 1-Z plane
!                ( 6) Col  6: G2Z   Transverse shear modulus in 2-Z plane
!                ( 7) Col  7: RHO   Mass density
!                ( 8) Col  8: A1    Thermal expansion coeff in the longitudinal direction
!                ( 9) Col  9: A2    Thermal expansion coeff in the lateral direction
!                (10) Col 10: TREF  Reference temperature for thermal expansion
!                (11) Col 11: Xt    Allowable stress or strain in tension in the longitudinal directio
!                (12) Col 12: Xc    Allowable stress or strain in compression in the longitudinal direction
!                (13) Col 13: Yt    Allowable stress or strain in tension in the lateral direction
!                (14) Col 14: Yc    Allowable stress or strain in compression in the lateral direction
!                (15) Col 15: S     Allowable stress or strain for in-plane shear
!                (16) Col 16: GE    Damping coeff
!                (17) Col 17: F12   Interaction term in the tensor polynomial theory of failure
!                (18) Col 18: STRN  For max strain failure theory only. Indicates whether allowables are stress or strain allowables

      INTEGER(LONG)                :: I, J                        ! loop counter
      INTEGER(LONG)                :: ITABLE                      ! the subtable id
      INTEGER(LONG)                :: NUM_WIDE                    ! the width of a single card
      INTEGER(LONG)                :: NVALUES                     ! the number of words in the block
      CHARACTER(8*BYTE), PARAMETER :: TABLE_NAME = "MPT"          ! the table_name we're writing
      INTEGER(LONG)                :: IMAT1 = 0, IMAT2 = 0, IMAT8 = 0, IMAT9 = 0

      LOGICAL :: IS_MPT

      IS_MPT = ( (NMATL > 0) )
      IS_MPT = .FALSE.
      IF (IS_MPT) THEN
        DO I=1,NMATL
          IF      (MATL(I,2) == 1) THEN
            IMAT1 = IMAT1 + 1
            MAT1_INDEX(IMAT1) = I
          ELSE IF (MATL(I,2) == 2) THEN
            IMAT2 = IMAT2 + 1
            MAT2_INDEX(IMAT2) = I
          ELSE IF (MATL(I,2) == 8) THEN
            IMAT8 = IMAT8 + 1
            MAT8_INDEX(IMAT8) = I
          ELSE IF (MATL(I,2) == 9) THEN
            IMAT9 = IMAT9 + 1
            MAT9_INDEX(IMAT9) = I
          ENDIF
        ENDDO
        !  MATL = Array of integer data from Material Bulk Data entries. Each row is for one Material entry read in B.D. and contains:
        ! (1) Col  1: Material ID
        ! (2) Col  2: Material type: 1, in MATL array for MAT1 entry (isotropic)
        !                            2, in MATL array for MAT2 entry (anisotropic, 2D elements)
        !                            8, in MATL array for MAT8 entry (orthotropic)
        !                            9, in MATL array for MAT9 entry (anisotropic, 3D elements)

        ! MAT1_INDEX, MAT2_INDEX, MAT8_INDEX, MAT9_INDEX
        CALL WRITE_OP2_GEOM_HEADER(TABLE_NAME, ITABLE)
        CALL END_OP2_GEOM_TABLE(ITABLE)
      ENDIF
      END SUBROUTINE WRITE_OP2_GEOM_MPT

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_GEOM_HEADER(TABLE_NAME, ITABLE)
      USE PENTIUM_II_KIND, ONLY         :  LONG, BYTE
      USE IOUNT1, ONLY                  :  ERR,OP2
      IMPLICIT NONE
      CHARACTER(LEN=8*BYTE),INTENT(IN)  :: TABLE_NAME
      INTEGER(LONG),INTENT(INOUT)       :: ITABLE
 1    FORMAT("writing table_name ",A)
      WRITE(ERR,1) TABLE_NAME

      CALL WRITE_TABLE_HEADER(TABLE_NAME)
      !WRITE(OP2) 7
      !WRITE(OP2) 1, 2, 3, 4, 5, 6, 7

      !WRITE(OP2) -2
      !WRITE(OP2) 1
      !WRITE(OP2) 0

      !WRITE(OP2) 2
      !WRITE(OP2) 1, 2

      WRITE(OP2) -3
      WRITE(OP2) 1
      WRITE(OP2) 0
      ITABLE = -3
      END SUBROUTINE WRITE_OP2_GEOM_HEADER

!===================================================================================================================================
      SUBROUTINE WRITE_OP2_SUBTABLE_INCREMENT( ITABLE )
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      IMPLICIT NONE
      INTEGER(LONG),INTENT(INOUT)       :: ITABLE

      ITABLE = ITABLE - 1
      WRITE(OP2) ITABLE
      WRITE(OP2) 1
      WRITE(OP2) 0
      END SUBROUTINE WRITE_OP2_SUBTABLE_INCREMENT

!===================================================================================================================================
      SUBROUTINE END_OP2_GEOM_TABLE( ITABLE )
      USE PENTIUM_II_KIND, ONLY         :  LONG
      USE IOUNT1, ONLY                  :  OP2
      IMPLICIT NONE
      INTEGER(LONG),INTENT(INOUT)       :: ITABLE

      !WRITE(OP2) 3
      !WRITE(OP2) 1, 2, 3
      !ITABLE = ITABLE - 1
      !
      !WRITE(OP2) ITABLE
      !WRITE(OP2) 1
      !WRITE(OP2) 0
      !
      !WRITE(OP2) 0
      !ITABLE = ITABLE - 1

      WRITE(OP2) 0
      END SUBROUTINE END_OP2_GEOM_TABLE

!==================================================================================================
      SUBROUTINE WRITE_TABLE_HEADER(TABLE_NAME)
      ! the hard part is getting the date
      !https://docs.oracle.com/cd/E19957-01/805-4942/6j4m3r8t2/index.html
      !
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG
      USE IOUNT1, ONLY           :  ERR, OP2
      CHARACTER(LEN=8*BYTE), INTENT(IN) :: TABLE_NAME ! The table name
      !INTEGER(LONG) :: ITABLE ! the subtable id
      INTEGER(LONG), DIMENSION(8) :: DATE_TIME
      INTEGER(LONG)               :: MONTH, DAY, YEAR
      CHARACTER(LEN=10*BYTE) B(3)

      CALL DATE_AND_TIME(B(1), B(2), B(3), DATE_TIME)
      YEAR  = DATE_TIME(1)
      MONTH = DATE_TIME(2)
      DAY   = DATE_TIME(3)

      WRITE(ERR,9110) TABLE_NAME

!      table0 = [
!        4, 2, 4,
!        8, table_name.encode('ascii'), 8,
!        #4, 0, 4,
!      ]
      WRITE(OP2) 2
      WRITE(OP2) TABLE_NAME
      !WRITE(OP2) 0

      !data_a = [4, -1, 4,]
      !data_c = [4, 7, 4,]
      WRITE(OP2) -1
      WRITE(OP2) 7
!      table1 = [
!        28,
!        102, 0, 0, 0, 512, 0, 0,
!        28,
!      ]
!      WRITE(OP2) 102, 0, 0, 8, 0,   0, 0
      WRITE(OP2) 102, 0, 0, 0, 512, 0, 0
!      data = [
!        4, -2, 4,
!        4, 1, 4,
!        4, 0, 4,
!      ]
      WRITE(OP2) -2
      WRITE(OP2) 1
      WRITE(OP2) 0

!      table2 = [
!        4, 7, 4,
!        28,  # 4i -> 13i
!        # todays date 3/6/2014, 0, 1  ( year=year-2000)
!        month, day, dyear, 0, 1,
!        28,
!      ]
      WRITE(OP2) 7
      WRITE(OP2) 0, 1, MONTH, DAY, YEAR - 2000, 0, 1

      !ITABLE = -3
      !CALL WRITE_ITABLE(ITABLE)
9110 FORMAT(" *DEBUG:       WRITE WRITE_TABLE_HEADER; TABLE_NAME=", A)
      END SUBROUTINE WRITE_TABLE_HEADER

! ##################################################################################################################################
      SUBROUTINE WRITE_ITABLE(ITABLE)
      USE PENTIUM_II_KIND, ONLY  :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY           :  ERR, OP2
      IMPLICIT NONE
      INTEGER(LONG), INTENT(IN) :: ITABLE   ! The subtable id
!      INTEGER(LONG), INTENT(IN) :: NTOTAL   ! the width of the block

      WRITE(OP2) ITABLE
      WRITE(OP2) 1
      WRITE(OP2) 0
!     WRITE(ERR,*) " *INFORMATION: "
!      WRITE(ERR,9114) " *DEBUG:       WRITE ITABLE; ITABLE=", ITABLE
      WRITE(ERR,9114) ITABLE

9114 FORMAT(" *DEBUG:       WRITE ITABLE; ITABLE=", I8)
      END SUBROUTINE WRITE_ITABLE


   END MODULE OP2_GEOMETRY_OUTPUT
