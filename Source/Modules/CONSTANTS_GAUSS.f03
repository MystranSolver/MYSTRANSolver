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

      MODULE CONSTANTS_GAUSS

! Real parameters used in subroutine ORDER for Gaussian integration

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE

      IMPLICIT NONE

      SAVE

      REAL(DOUBLE), PARAMETER         :: SSV(29) = (/ 0.577350269189626D0, & ! SSV( 1)
                                                      0.774596669241483D0, & ! SSV( 2)
                                                      0.000000000000000D0, & ! SSV( 3)
                                                      0.861136311594053D0, & ! SSV( 4)
                                                      0.339981043584856D0, & ! SSV( 5)
                                                      0.906179845938664D0, & ! SSV( 6)
                                                      0.538469310105683D0, & ! SSV( 7)
                                                      0.000000000000000D0, & ! SSV( 8)
                                                      0.932469514203152D0, & ! SSV( 9)
                                                      0.661209386466265D0, & ! SSV(10)
                                                      0.238619186083197D0, & ! SSV(11)
                                                      0.949107912342759D0, & ! SSV(12)
                                                      0.741531185599394D0, & ! SSV(13)
                                                      0.405845151377397D0, & ! SSV(14)
                                                      0.000000000000000D0, & ! SSV(15)
                                                      0.960289856497536D0, & ! SSV(16)
                                                      0.796666477413627D0, & ! SSV(17)
                                                      0.525532409916329D0, & ! SSV(18)
                                                      0.183434642495650D0, & ! SSV(19)
                                                      0.968160239507626D0, & ! SSV(20)
                                                      0.836031107326636D0, & ! SSV(21)
                                                      0.613371432700590D0, & ! SSV(22)
                                                      0.324253423403809D0, & ! SSV(23)
                                                      0.000000000000000D0, & ! SSV(24)
                                                      0.973906528517172D0, & ! SSV(25)
                                                      0.865063366688985D0, & ! SSV(26)
                                                      0.679409568299024D0, & ! SSV(27)
                                                      0.433395394129247D0, & ! SSV(28)
                                                      0.148874338981631D0 /) ! SSV(29)

      REAL(DOUBLE), PARAMETER         :: HHV(29) = (/ 1.000000000000000D0, & ! HHV( 1)
                                                      0.555555555555556D0, & ! HHV( 2)
                                                      0.888888888888889D0, & ! HHV( 3)
                                                      0.347854845137454D0, & ! HHV( 4)
                                                      0.652145154862546D0, & ! HHV( 5)
                                                      0.236926885056189D0, & ! HHV( 6)
                                                      0.478628670499366D0, & ! HHV( 7)
                                                      0.568888888888889D0, & ! HHV( 8)
                                                      0.171324492379170D0, & ! HHV( 9)
                                                      0.360761573048139D0, & ! HHV(10)
                                                      0.467913934572691D0, & ! HHV(11)
                                                      0.129484966168870D0, & ! HHV(12)
                                                      0.279705391489277D0, & ! HHV(13)
                                                      0.381830050505119D0, & ! HHV(14)
                                                      0.417959183673469D0, & ! HHV(15)
                                                      0.101228536290376D0, & ! HHV(16)
                                                      0.222381034453374D0, & ! HHV(17)
                                                      0.313706645877887D0, & ! HHV(18)
                                                      0.362683783378362D0, & ! HHV(19)
                                                      0.081274388361574D0, & ! HHV(20)
                                                      0.180648160694857D0, & ! HHV(21)
                                                      0.260610696402935D0, & ! HHV(22)
                                                      0.312347077040003D0, & ! HHV(23)
                                                      0.330239355001260D0, & ! HHV(24)
                                                      0.066671344308688D0, & ! HHV(25)
                                                      0.149451349150581D0, & ! HHV(26)
                                                      0.219086362515982D0, & ! HHV(27)
                                                      0.269266719309996D0, & ! HHV(28)
                                                      0.295524224714753D0 /) ! HHV(29)

      END MODULE CONSTANTS_GAUSS
