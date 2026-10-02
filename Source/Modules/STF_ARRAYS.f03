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

      MODULE STF_ARRAYS

! Initial, linked list, representation of the G-set stiffness matrix. This is the form in which the individual element stiffness
! matrices are assembled into the system stiffness matrix in subr ESP. Subr SPARSE_KGG taskes this form and creates the sparse
! CRS (compressed row storage) form (I_KGG, J_KGG, KGG) for use in LINK2 and later.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE DERIVED_DATA_TYPES, ONLY    :  INT2_REAL1

      IMPLICIT NONE

      SAVE

      TYPE(INT2_REAL1), ALLOCATABLE   :: STF3(:)

      INTEGER(LONG), ALLOCATABLE      :: STFCOL(:)         ! See explanation below
      INTEGER(LONG), ALLOCATABLE      :: STFPNT(:)         ! See explanation below
      INTEGER(LONG), ALLOCATABLE      :: STFKEY(:)         ! See explanation below

      REAL(DOUBLE) , ALLOCATABLE      :: STF(:)            ! See explanation below


!  STF(I)    = 1D real array of the nonzero terms, above the diag, in the G-set stiffness matrix. Size is NTERM_KGG.

!  STFKEY(I) = 1D integer array of pointers to where, in STF, the 1st stiffness term for row I exists.
!              STFKEY(I) = 0 if row I  is null. Size is NDOFG.

!  STFPNT(I) = 1D integer array of pointers to where, in STF, the 2nd and remaining stiffness terms for row I exist.
!              The 2nd term in row I is at STFPNT(STFKEY(I))          in STF.
!              The 3rd term IN ROW i is at STFPNT(STFPNT(STFKEY(I))), etc.
!              When a zero value of STFPNT is reached in this scheme it means that there are no more terms in row I.
!              Size is NTERM_KGG.

!  STFCOL(I) = 1D integer array of the column numbers of the terms in STF(I). Size is NTERM_KGG.

! An example of this is given below.

! Consider the 6 x 6 stiffness matrix:

!                        1   2   3   4   5   6
!
!                      | 5.  0.  0.  0. -3. -2.| 1
!                      |     9.  0. -5. -4.  0.| 2
!                  K = |         6. -6.  0.  0.| 3
!                      |            11.  0.  0.| 4
!                      |  SYMMETRIC      7.  0.| 5
!                      |                     3.| 6

! The 1D integer arrays STFKEY, STFCOL, STFPNT and the 1D real array STF are generated in subroutine ESP based on
! elem stiffness  matrices (which combine to produce the above global K). The arrays  are produced by calculating
! each element stiffness matrix (in their element number numerical order - order of ESORT2 after subr
! ELESORT). The terms from each of these elem stiffness matrices are written to STF as they are generated (values are
! accumulated at a common DOF when more than 1 element contributes to the stiffness at that DOF).
! For the above global stiffness matrix, the arrays generated in ESP are:


!             I   STFKEY(I)   STFPNT(I)   STFCOL(I)  STF(I)
!
!             1           8           2           3      6.
!----------------------------------------------------------
!             2           4           0           4     -6.
!----------------------------------------------------------
!             3           1           0           4     11.
!----------------------------------------------------------
!             4           3           5           4     -5.
!----------------------------------------------------------
!             5           7           6           2      9.
!----------------------------------------------------------
!             6          11           0           5     -4.
!----------------------------------------------------------
!             7                       0           5      7.
!----------------------------------------------------------
!             8                       9           5     -3.
!----------------------------------------------------------
!             9                      10           1      5.
!----------------------------------------------------------
!            10                       0           6     -2.
!----------------------------------------------------------
!            11                       0           6      3.

! In this context it can be seen that STFPNT can be viewed as an array of pointers, one array for each row. The
! STFPNT values for any one row give the value of the stiffness term (STF) and the column it goes in (STFCOL)

! The way that the above arrays are interpreted to get the stiff matrix K (using IS as a intermediate variable) is:

! Row 1 (I=1).
! ------------
!   Subroutine ESP generates STF, STFKEY, STFCOL, STFPNT:
!
!   IS = STFKEY( I) =  8, then STF( 8) = -3. goes in col STFCOL( 8) = 5

!   IS = STFPNT(IS) =  9, then STF( 9) =  5. goes in col STFCOL( 9) = 1

!   IS = STFPNT(IS) = 10, then STF(10) = -2. goes in col STFCOL(10) = 6

!   IS = STFPNT(IS) =  0, No more terms in row I=1

!   Therefore, row 1 has -3. in col 5, 5. in col 1 and -2. in col 6 which is what is shown in the K matrix above.

!   When subroutine SPARSE_kgg is run the 1st row of tha G-set stiffness, prior to sorting contains the terms above:

!   KGG value = -3.   5.  -2.    (stiffness values in row 1)
!   KGG col   =  5    1    6     (column numbers of terms in row 1)

!   When SPARSE_KGG has sorted the columns, prior to writing the row to file the 1st row is:

!   KGG value =  5.  -3.  -2.    (stiffness terms in row 1)
!   KGG col   =  1    5    6     (column numbers of terms in row 1)


! Repeat the above for rows 2 - 6 to get all terms above the diag of K from the STFKEY, STFCOL, STFPNT, STF arrays.

! By applying the values in STFKEY for each of rows 1 through 5 we can write the above as:

!                       Row      STFPNT      STFCOL    STF

!                         1           9           5     -3.
!                                    10           1      5.
!                                     0           6     -2.
!
!                         2           5           4     -5.
!                                     6           2      9.
!                                     0           5     -4.
!
!                         3           2           3      6.
!                                     0           4     -6.
!
!                         4           0           4     11.
!
!                         5           0           5      7.
!
!                         6           0           6      3.


      END MODULE STF_ARRAYS
