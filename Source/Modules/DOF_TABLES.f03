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

      MODULE DOF_TABLES

! Arrays used to describe the MYSTRAN DOF set. These are used extensively throught the program.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  TSET_CHR_LEN

      IMPLICIT NONE

      SAVE


      CHARACTER(LEN(TSET_CHR_LEN)), ALLOCATABLE                                                                                    &
                                      :: TSET(:,:)

      CHARACTER(LEN(TSET_CHR_LEN)), ALLOCATABLE                                                                                    &
                                      :: USET(:,:)

      CHARACTER(2*BYTE)               :: USETSTR_TABLE(16,2)

!                         1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
      DATA USETSTR_TABLE/'G ','M ','N ','SA','SB','SG','SZ','SE','S ','F ','O ','A ','R ','L ','U1','U2',                          &
                         '0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 ','0 '/

      INTEGER(LONG), ALLOCATABLE      :: TDOF(:,:)         ! NDOFG x MTDOF integer array of DOF numbers for every G.P.and comp in
!                                                            all of the displ sets sorted in grid point numerical order

      INTEGER(LONG), ALLOCATABLE      :: TDOFI(:,:)        ! TDOF sorted in G-set DOF order

      INTEGER(LONG), ALLOCATABLE      :: TDOF_ROW_START(:) ! Row number in TDOF where data begins for a GRID or SPOINT

      END MODULE DOF_TABLES
