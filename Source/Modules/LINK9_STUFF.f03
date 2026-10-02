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

      MODULE LINK9_STUFF

! Grid point and element solution variables for data recovery LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      CHARACTER( 1*BYTE), ALLOCATABLE :: MSPRNT(:,:)           ! Flags for whether to print margins of safety for ROD, BAR
      CHARACTER( 4*BYTE), ALLOCATABLE :: FTNAME(:)             ! Stress failure index name output with stresses/strains

      INTEGER(LONG)                   :: MAXREQ                ! Max number of rows needed for array OGEL

      INTEGER(LONG)     , ALLOCATABLE :: GID_OUT_ARRAY(:,:)    ! Array of integer grid no's for some output in LINK9

      INTEGER(LONG)     , ALLOCATABLE :: EID_OUT_ARRAY(:,:)    ! Array of elem no's (col 1) and num of plies (col 2) for that elem
!                                                                that are printed with certain outputs IN LINK9

      INTEGER(LONG)     , ALLOCATABLE :: POLY_FIT_ERR_INDEX(:)! Index num for POLY_FIT_ERR (i.e. which of the 1 through 9 stress
!                                                                or strain values has the largest error in polynomial fit


      REAL(DOUBLE)      , ALLOCATABLE :: OGEL(:,:)             ! Master array for holding outputs in LINK9 until they are printed

      REAL(DOUBLE)      , ALLOCATABLE :: POLY_FIT_ERR(:)       ! Array of polynom fit errors for elems that extrapolate stress or
!                                                                strain values from one set of output points to another

      END MODULE LINK9_STUFF
