! ##################################################################################################################################

      MODULE BANDIT_MODULE

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                   :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                   :  HOUR, MINUTE, SEC,
     &                                      SFRAC, TSEC
      USE PARAMS, ONLY                   :  DELBAN
      USE BANDIT_FILE_SUPPORT, ONLY      :  BANDIT_FILES

! Notes:
! ------

! (1) The correspondence of array IPARAM and the $ directive entries described in Gordon's documentation is as follows:

!          IPARAM( 1) corresponds to $LOOP       (not documented: $LOOP is used for program testing - loop on several input decks)
!          IPARAM( 2) corresponds to $TABLE
!          IPARAM( 3) corresponds to $METHOD
!          IPARAM( 4) corresponds to $MPC
!          IPARAM( 5) corresponds to $SEQUENCE
!          IPARAM( 6) corresponds to $CRITERION
!          IPARAM( 7) corresponds to $SPRING
!          IPARAM( 8) corresponds to $NASTRAN    (not documented - $NASTRAN YES card causes message to be printed to output file )
!          IPARAM( 9) corresponds to $ELEMENTS
!          IPARAM(10) corresponds to $PRINT
!          IPARAM(11) corresponds to $FRONTAL

! (2) The statement IPARAM(I) = J has the following meaning for various J:

!          J =  1 means IPARAM(I) is the I-th directive (above) with value SEQGP
!          J =  2 means IPARAM(I) is the I-th directive (above) with value ALL
!          J =  3 means IPARAM(I) is the I-th directive (above) with value NO
!          J =  4 means IPARAM(I) is the I-th directive (above) with value YES
!          J =  5 means IPARAM(I) is the I-th directive (above) with value MIN
!          J =  6 means IPARAM(I) is the I-th directive (above) with value MAX
!          J =  7 means IPARAM(I) is the I-th directive (above) with value CM
!          J =  8 means IPARAM(I) is the I-th directive (above) with value GPS
!          J =  9 means IPARAM(I) is the I-th directive (above) with value BOTH
!          J = 10 means IPARAM(I) is the I-th directive (above) with value BAND
!          J = 11 means IPARAM(I) is the I-th directive (above) with value PROFILE
!          J = 12 means IPARAM(I) is the I-th directive (above) with value RMS
!          J = 13 means IPARAM(I) is the I-th directive (above) with value WAVEFRONT

! For the most part, all changes to Gordon's code (except cosmetic ones) are delineated as shown in the example below where Program
! BANDIT was changed to a subr and where file unit in1 from module iount1 is added for use.

! Example of how mods are delineated:
! B////////////////////////////////////////////////////////////////////B
! E////////////////////////////////////////////////////////////////////E


! Mods to Gordon's original bandit.f file (for use as a set of subrs called by MYSTRAN):
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! Aug 03 - Oct 03 Mods:
! *********************

!    General changes:
!    ----------------
!    ( 1) Add ! ##############... comment lines as separators between each subroutine
!    ( 2) Modify subroutine END statements to append subroutine name to the END statement

!    Main:
!    -----
!    ( 1) Change Program BANDIT to subroutine bandit and add statement: "use iount1, only:  in1"
!    ( 2) Add arg IER to subr BANDIT so that MYSTRAN can check for error and give message.
!    ( 3) Do not open IOU5, instead set iou5 = in1 (MYSTRAN input file)
!    ( 4) Change the statement: 135 IOU5=5, to 135 iou5 = in1 to reflect that IOU5 is the MYSTRAN input file.
!    ( 5) Change STOP 13 and STOP 5 to writing a message and returning to MYSTRAN
!    ( 6) Change final STOP to RETURN and change END to end subroutine bandit
!    ( 7) Change open status from 'UNKNOWN' to 'REPLACE' for units IOU6,7,8,9,10,14,16,17. This is done to make sure that MYSTRAN
!         does not attempt to use an old file.

!    BLOCK DATA:
!    -----------
!    ( 1) Remove definition of unit IOU5 (since it is now MYSTRAN IN1
!    ( 2) for all other file numbers, add 100 (e.g. IOU6 goes from unit 6 to unit 106, etc.)

!    Subr DOLLAR:
!    ------------
!    ( 1) chge default for IPARAM(6), the variable for $CRITERION, from RMS to BAND

!    Subr ELTYPE:
!    ------------
!    ( 1) Comment out the code that sets fatal error and return if SEQGP cards are in the input deck. MYSTRAN checks
!         this (before Bandit runs) and ignores SEQGP cards if the user has requested Bandit resequencing.

!    Subr FINISH:
!    ------------
!    ( 1) Add: use iount1, only: seq
!    ( 2) Add several lines regarding MYSTRAN file SEQ (used for writing SEQGP card images as well on Bandit file unit IOU6)

!    Subr SEQGP:
!    -----------
!    ( 1) Add: use iount1, only: seq
!    ( 2) Add WRITE(IOU7,*) just before ENDFILE IOU7 to get 1 blank line written to IOU7 (in case there are no SEQGP card images
!         written to IOU7)
!    ( 3) Add several lines regarding MYSTRAN file SEQ so that the SEQGP card images are written to that file as well as Bandit IOU6


! 11/30/03 Mods:
! **************

!     Subr TIMER:
!     -----------
!     ( 1) Comment out 2 lines in subroutine TIMER (variables tarray, iwall, time not used and Layhey complains about type specifier)
!     ( 2) Change "call second(t)" to "call cpu_time(t)". Calling second was causing the -NaN problem when BANDIT was imbedded in
!          MYSTRAN since there is no intrinsic function in Layhey called  second. Not sure what is being used since there is no
!          intrinsic  function in Layhey documentation called second. There is a subr called second in module ARPACK_UTIL that I
!          modified to call cpu_time and, if we put the statement: use arpack_util above then we get no problem with -NaN. Module
!          ARPACK_UTIL calls cpu_time to get t and this is a valid procedure in Layhey. Thus, switch to call cpu_time here.


! Mods to Gordon's second bandit.f file (after he added RBE's):
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! 01/08/04 Mods:
! **************

!     General changes:
!     ----------------
!     ( 1) Change this collection of subroutines to a module and change its name to BANDIT_MODULE. This required changing several
!          scalar variables to appear as if they are arrays to avoid Lahey giving fatal errs regarding shapes of arrays not matching

!     Subr bandit:
!     ------------
!     ( 1) Earlier changes regarding STOP 13 and STOP 5 no longer apply. Gordon has removed these stops and, apparently, is handling
!          these situations in subr FINISH

!     BLOCK DATA:
!     -----------
!     ( 1) Remove this code from the module (it can't be in a module) and create a new procedure called BLOCK DATA BANDIT_BLOCK_DATA

!     Subr ELTYPE:
!     ------------
!     ( 1) Add integer array NCON_array(1) to be equiv to scalar NCON and call subr READIT with NCON_array(1) instead of NCON

!     Subr FINISH:
!     ------------
!     ( 1) Cosmetic changes (not delineated with ! B////////...) to make subr more readable

!     Subr GIBSTK:
!     ------------
!     ( 1) Remove SORT2 from INTEGER declar. It is an integer function that is now part of this module and can't be redefined here

!     Subr INSERT:
!     ------------
!     ( 1) Add integer array NCARD_array(1) to be equiv to scalar NCARD and call subr READIT with NCARD_array(1) instead of NCARD

!     Subr REED:
!     ----------
!     ( 1) Add integer array EID_array(1) to be equiv to scalar EID and call subr READIT with EID_array(1) instead of EID

!     Subr RIGID:
!     -----------
!     ( 1) Add integer array IG_array(1) to be equiv to scalar IG and call subr READIT with IG_array(1) instead of IG (2 locations)
!     ( 2) Call SCAT with integer array IG_array(1) instead of scalar IG

!     Subr TAXI:
!     ----------
!     ( 1) Add integer array NAXIC_array(1) to be equiv to scalar NAXIC and call subr READIT with NAXIC_array(1) instead of NAXIC

! 01/19/04 Mods:
! **************

!     General changes:
!     ----------------
!     ( 1) Add type declarations for all variables not alre3ady typed so that IMPLICIT NONE switch can be used in the compiler
!          directive.

! 02/01/04 Mods:
! **************

!     Subr BANDIT:
!     ------------
!     ( 1) Add MYSTRAN_NGRID to arg list (INTENT(IN)) in SUBROUTINE BANDIT and change NGRID = KOR/30  to MAX(MYSTRAN_NGRID,KORE/30)
!     ( 2) Just before setting NGRID as stated above, comment out the line KDIM=250 and, just after the line NGRID = MYSTRAN_NGRID
!          set KDIM = MAX(250,NGRID/10)

! 02/21/04 Mods:
! **************

!     Subr RIGID:
!     ------------
!     ( 1) Change dimension of LG(2) to LK(*) per 2/19/04 email from Gordon Everstine

! 03/02/04 Mods:
! **************

!     Subr BANDIT:
!     ------------
!     ( 1) Change to KDIM = MAX(250,MYSTRAN_NGRID) . It apparently is not large enough when there is a large "Nodal degree limit"
!          (as when RBE's connect to many grids). Later in Bandit, KDIM seems to be set to as large as NGRID so use max of the
!          250 that Gordon had here and NGRID.

!     ( 2) Change the 02/01/04 mod to: change original NGRID=KORE/30 to NGRID = MAX(MYSTRAN_NGRID,10).
!          The 10 is a guess to cover small problems since NGRID = MYSTRAN_NGRID doesn't seem to work for small problems and since
!          Bandit uses NGRID/10 later.
!          I had problem with chassis_grav. Bandit wouldn't run because Nodal degree limit of 21 was exceeded. Apparently, since
!          NGRID = MAX(MYSTRAN_NGRID,KORE/30) used KORE/30 = 10,000,000/30 (since KORE/30 > MYSTRAN_NGRID) was using much of the
!          memory for NGRID and then there was not enough for the Nodal degree limit. The Nodal degree limit can be set in Bandit
!          with the $GRID N directive. This sets NGRID to N, so when I input $GRID 7298 (actual number of grids in chassis_grav)
!          Bandit ran successfully.

! 08/01/06 Mods:
! **************

!     Subr BANDIT:
!     ------------

!     ( 1) Add var NEW_BW to be returned from subr SUMUP to be the new bandwidth found
!     ( 2) Change call to SUMUP to have NEW_BW as an arg
!     ( 3) Add var NEW_BW to subr SUMUP ans set NEW_BW = NBW (from COMMON/D/ which is printed out there as the new bandwidth)

! 10/05/06 Mods:
! **************

!     Subr BANDIT_MODULE:
!     ------------------

!     Add call to BANDIT_FILES (new subr to make sure no files are open in another program and all old files are deleted)

! 07/23/08 Mods:
! **************

!     (1) Add DEN to arg list in CALL SUMUP
!     (2) Add DEN to arg list in subr SUMUP (returned to calling program)
!     (3) Add DEN as return arg to BANDIT subr

! 01/26/09 Mods:
! **************

!     Mod to version 6.00c
!     Change FORMAT #10 in subr SEQGP to write SEQGP card images in large field format. Need this when MYSTRAN SEQGP calls
!     subr BD_SEQGP. This was causing an error when GRID numbers in large field format were used

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################
! B////////////////////////////////////////////////////////////////////B
      SUBROUTINE BANDIT ( MYSTRAN_NGRID, NEW_BW, DEN, IER )
      USE IOUNT1, ONLY :  IN1
! E////////////////////////////////////////////////////////////////////E
c
c     Gordon C. Everstine, Gaithersburg, MD, geversti@comcast.net
c     1/8/04
c
c     This is the PC version.  The only machine or compiler dependencies
c     are in Subroutine Timer related to getting the CPU time.  This
c     subroutine must be checked before compiling.
c
c     Initial Bandit version: December 1969
c
c  Recent revisions:
c     4/16/93  Read MPCAX cards.
c     11/8/93  Minor changes needed to suit MS-FORTRAN.
c     4/22/94  Write to file the components and grids (unit 15) if
c              CM is executed.
c     12/5/94  Move connection table output from unit 11 to 16.
c              Make unit 11 a scratch file.
c              Delete selected files at end of job.
c              Add elapsed time report for HP workstations.
c     1/13/95  Add log file, unit 17; increase default KDIM to 250.
c     5/22/96  Speed up element sort in Subroutine FRONT.
c              Change file names from TAPE* to bandit.f*.
c    11/18/96  Change author's mailing address.
c     2/25/98  Add timing call for PC (subroutine timer).
c     6/20/03  Change author's email address.
c    10/27/03  Enable component table with $print max (unit 15).
c    11/20/03  Add additional elastic elements.
c      1/8/04  Add additional rigid elements and some cosmetic changes.
c              Remove extra plus signs, a holdover from BCD/EBCDIC days.
c
C  This is a NASTRAN preprocessor which reads a NASTRAN deck and
C  resequences the grid points to reduce the matrix bandwidth, profile,
C  or wavefront.  Two algorithms are provided,
C  Gibbs-Poole-Stockmeyer (GPS) and Reverse Cuthill-McKee (CM).
C  In addition, the elements can be resequenced for frontal solvers
C  using the Razzaque approach to determine an element sequence,
C  given a grid point sequence determined by GPS or CM.
C
C  The input data deck to BANDIT consists of a standard NASTRAN deck
C  plus one or more $ cards to tell BANDIT what to do.  Output includes
C  a set of SEQGP card images on bandit.f07 and the complete NASTRAN deck
C  (including SEQGP cards) on bandit.f08.  The new element sequence is on
C  bandit.f14.  See "BANDIT User's Guide" by G.C. Everstine.
C
C                    Program Installation
C                    --------------------
C  Set variable MEM (the amount of working storage in real, single-
C  precision words) to the desired value.
c
! B////////////////////////////////////////////////////////////////////B
      INTEGER, PARAMETER :: MEM = 10 000 000
! E////////////////////////////////////////////////////////////////////E
c
C  To get timing information, add in Subroutine TIMER a call to the
C  appropriate CPU clock routine.
c
C  Integer-packing is not available in this version, but may be
C  enabled by acquiring from the author the packing routines and
C  replacing some code wherever the string "CPACK" appears.
c  With large memories, integer-packing is rarely needed today and
c  inhibits vectorization.
C
C---------------------------------------------------------------------
C
C  I/O FILES - - -
C
C  (SEE BLOCK DATA FOR DEFINITION OF FILES IN COMMON BLOCK /IOUNIT/)
C
C  Unit  File name       Format       I/O     Opn STAT                                  Use
C  ----  -----------   -----------   ------   --------  ---------------------------------------------------------------------------
C    5   INFILE          FORMATTED   input    OLD        MYSTRAN input file
C    6   bandit.out'     FORMATTED   output   REPLACE    Bandit printed output file
C    7   bandit.f07'     FORMATTED   output   REPLACE    SEQGP card images
C    8   bandit.f08'     FORMATTED   output   REPLACE    MYSTRAN input file with SEQGP card images
C    9   bandit.f09'     FORMATTED   output   REPLACE    References in NASNUM, SPRING, and FINISH
C   10   bandit.ins'     FORMATTED   input    UNKNOWN    Insert file (card images)
C   11   bandit.f11'     FORMATTED   output   REPLACE    Scratch file - references in DOLLAR, NASNUM, and NEWIN
C   12   bandit.f12'   UNFORMATTED   output   REPLACE    Scratch file - references in MPC, RESET, RIGID, and TIGER
C   13   bandit.f13    UNFORMATTED   output   SCRATCH    Scratch file - store elems for generating elem ordering for frontal solver
C   14   bandit.f14'     FORMATTED   output   REPLACE    Output new element list for frontal solution
C   15   bandit.f15'     FORMATTED   output   REPLACE    Component and grid list (CM must be run to get this list)
C   16   bandit.f16'     FORMATTED   output   REPLACE    Connection table output ($TABLE YES)
C   17   bandit.log'     FORMATTED   output   REPLACE    Some run-time messages to benefit interactive running
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IADD   ,IB     ,IBYTE  ,IDIM   ,IER    ,IFIR   ,
     &         IFL    ,IIG    ,IGNORE ,IGDEG  ,IH     ,INP    ,IOP    ,
     &         IPARAM ,IPASS  ,ISTA   ,ISTART ,II1    ,II3    ,IWALL1 ,
     &         IWALL2

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  K1     ,K2     ,K3     ,K4     ,K5     ,K6     ,K7     ,
     &         K8     ,K9     ,KDIM4  ,KDIM   ,KMOD   ,KNEW   ,KORE   ,
     &         KORIG

      INTEGER  LINES

      INTEGER  MA     ,MB     ,MAXDEG ,MAXGRD ,MINDEG ,MDIM   ,ME     ,
     &         MM

      INTEGER  NAXIC  ,NBITIN ,NBYTE  ,NCM    ,NEDGE  ,NEL    ,NELEM  ,
     &         NEQ    ,NEQR   ,NGRID  ,NLINK  ,NN     ,NW     ,NTYPE  ,
     &         NUM    ,NZERO

      REAL     DUMG   ,DUMW   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      INTEGER   NEW_BW
      REAL      DEN
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'BANDIT'
      INTEGER(LONG), INTENT(IN)       :: MYSTRAN_NGRID
! E////////////////////////////////////////////////////////////////////E

      COMMON /S/ NN,MM,IH,IB,LINES,NEDGE,IADD,MINDEG,NAXIC
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /BITS/ NBITIN,KORE,IFL,NGRID,IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /D/ KORIG,KNEW,IOP,INP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /W/ DUMW(6)
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),
     -              NELEM(160),MDIM
      INTEGER vype,TYPE,WYPE
      COMMON /GRA/ DUMG(3)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      integer KOM(MEM)

! B////////////////////////////////////////////////////////////////////B

! E////////////////////////////////////////////////////////////////////E
C
      IPARAM(1)=0
C
C     OPEN FILES (EXCEPT MAYBE 5 AND 6)
C
! B////////////////////////////////////////////////////////////////////B
      IOU5 = IN1              ! IOU5 IS MYSTRAN INPUT FILE. ALREADY OPEN
! E////////////////////////////////////////////////////////////////////E
      CALL BANDIT_FILES ( IOU6 , IOU7 , IOU8 , IOU9 , IOU11, IOU12,
     &                    IOU13, IOU14, IOU15, IOU16, IOU17 )
c
      OPEN(IOU6,FILE='bandit.out',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU7,FILE='bandit.f07',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU8,FILE='bandit.f08',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU9,FILE='bandit.f09',FORM='FORMATTED',STATUS='replace')

c     OPEN(IOU11,FORM='FORMATTED',STATUS='SCRATCH')
c     OPEN(IOU12,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(IOU16,FILE='bandit.f16',FORM='FORMATTED',STATUS='replace')
c     OPEN(IOU17,FILE='bandit.log',FORM='FORMATTED',STATUS='replace')
c     See below for opening of scratch files
c     (placed there to avoid unexplained problem when looping on
c     multiple data sets)
C
C     TOP OF LOOP ON MULTIPLE DATA DECKS (used for code testing).
C     USE '$LOOP YES' IN FIRST DATA DECK TO ENABLE LOOPING AND
C     '$LOOP NO' IN LAST DATA DECK TO DISABLE LOOPING.
C
25    CALL TIMER(TA,IWALL1,0,IOU6)
      IF(IPARAM(1).EQ.4) Write(iou6,'(80(1H#))')
      WRITE(IOU6,30)
30    FORMAT('Bandit 1/8/04, G.C. Everstine, Gaithersburg, Maryland,',
     -       ' geversti@comcast.net')
c
c     Open scratch files as named files.
c
      OPEN(IOU11,file='bandit.f11',FORM='FORMATTED',STATUS='replace')
      OPEN(IOU12,file='bandit.f12',FORM='UNFORMATTED',STATUS='replace')
C
C     INITIALIZE SOME VARIABLES.
C
      KORE=MEM
      IFL=KORE
      IER=0
C     NAXIC=NUMBER OF HARMONICS ON AXIC CARD (DEFAULT SET HERE)
      NAXIC=99999
C
C     SET DEFAULT NUMBER OF PRINTED LINES PER PAGE.
C
      LINES=55
C     USE 55 FOR STANDARD 11 IN. PAPER, 40 FOR 8.5 IN. PAPER.
C
C     SET DEFAULT FOR KDIM AND MAX NUMBER OF GRID POINTS.
C
! B////////////////////////////////////////////////////////////////////B
!     KDIM=250
      KDIM=MAX(250,MYSTRAN_NGRID)
! E////////////////////////////////////////////////////////////////////E
      KDIM4=4*KDIM
      KORE=KORE-KDIM4
! B////////////////////////////////////////////////////////////////////B
!     NGRID=KORE/30
      NGRID = MAX(MYSTRAN_NGRID,10)
! E////////////////////////////////////////////////////////////////////E
      KORE=KORE+KDIM4
C
C     READ EXECUTIVE AND CASE CONTROL DECKS.
C
      CALL CASE(IER)
      IF(IER.GT.0) GO TO 135
C
C     PRINT ELEMENT LIBRARY IF REQUESTED.
C
c     ICHAR=MA(3)
      IF(IPARAM(9).EQ.4) then
         WRITE(IOU6,40) (I,vype(i),TYPE(I),WYPE(I),ME(I),I=4,NTYPE)
40       FORMAT(/,'Element Library:'/(I6,3X,A1,2A4,I7))
      end if
C
C     INITIALIZE ELEMENT COUNTERS.
C
      DO I=1,NTYPE
         NELEM(I)=0
      end do
C
      IF(IPARAM(5).EQ.4) GO TO 80
C
C     COPY BULK DATA TO UNIT 8 IF RESEQUENCING IS NOT REQUESTED.
C
      CALL NOSEQ(KOM)
C
      GO TO 120
C
C     COMPUTE MAXGRD AND MAXDEG.
C
80    CONTINUE
      KDIM4=4*KDIM
      KORE=KORE-KDIM4
      CALL GRID(IER,IOU6)
      IF(IER.GT.0) GO TO 135
C
C     PRINT CORE ALLOCATION INFORMATION.
C
      CALL COREKO
C
C     PARTITION BLANK COMMON.
C
      II1=MAXGRD/NW
      II3=2*MAXGRD
      K1=1
      K2=K1+KDIM4
      K3=K2+MAXGRD+1
      K4=K3+2*II3+1
      K5=K4+MAXGRD
      K6=K5+MAXGRD
      K7=K6+MAXGRD
      K8=K7+MAXDEG+1
      K9=K8+II1*MAXDEG
C     (K9 DEFINES BEGINNING OF NEXT ARRAY IF THERE WERE ONE.)
C     IG, THE BIG ARRAY IN NASNUM, IS LOCATED LAST IN OPEN CORE TO KEEP
C     CALLING ADDRESSES SMALLER, A BENEFIT ON UNIVAC WHEN USING OVER 65K.
C
C     READ BULK DATA, SET UP CONNECTION TABLE, RESEQUENCE NODES, AND
C            GENERATE SEQGP CARDS.
C
      CALL NASNUM(KOM(K8),II1,KOM(K3),II3,KOM(K4),KOM(K5),KOM(K2),
     -            KOM(K6),KOM(K7),KOM(K1),KDIM4,IER)
      IF(IER.GT.0) GO TO 135
C
C     PRINT OUT ELEMENT COUNTERS.
C
      WRITE(IOU6,90)
90    FORMAT(/,'Element Counts for Data Deck:')
      DO I=4,NTYPE
         IF(NELEM(I).gt.0)
     -        WRITE(IOU6,100) vype(i),TYPE(I),WYPE(I),NELEM(I)
100      FORMAT(5X,A1,2A4,I8)
      end do
C
C     PRINT BANDIT SUMMARY.
C
! B////////////////////////////////////////////////////////////////////B
      CALL SUMUP ( NEW_BW, DEN)
! E////////////////////////////////////////////////////////////////////E
C
120   CONTINUE
      ENDFILE IOU8
      REWIND IOU8
      IF(IPASS.GT.0) WRITE(IOU6,130) IPASS
130   FORMAT(/,'Number of calls to the pack/unpack routines',I15)
C
C     WRAP UP JOB
C
c     Reset iou5 in case there is looping on multiple data sets.
! B////////////////////////////////////////////////////////////////////B
!!135 IOU5=5
  135 IOU5 = IN1
! E////////////////////////////////////////////////////////////////////E
c     delete scratch files
      close(iou11,status='delete')
      close(iou12,status='delete')
      CALL TIMER(TB,IWALL2,0,IOU6)
      TB=TB-TA
      IWALL2=IWALL2-IWALL1
c     IF(IWALL2.GT.0) WRITE(IOU6,138) IWALL2
138   FORMAT(/,'Elapsed time =',I6,' seconds.')
      WRITE(IOU6,140) TB
140   FORMAT(/,'End of BANDIT Job.  Total CP time =',F12.3,' seconds.')
      Write(iou6,*)
C     IF LOOPING ON MULTIPLE DECKS IS REQUESTED, GO BACK TO BEGINNING.
      IF(IPARAM(1).EQ.4) GO TO 25
      if(iparam(7).eq.3) close(iou9,status='delete')
      if(iparam(2).eq.3) close(iou16,status='delete')
c     close(iou17,status='delete')
c     if(ier.gt.0) stop 13
C     $NASTRAN YES CARD RETURNS CONDITION CODE 5 FOR USE BY IBM.
c     (IBM JCL was written to check Bandit's exit code; now obsolete.
c     Stop 13 above served same purpose.)
c     IF(IPARAM(8).EQ.4) STOP 5
! B////////////////////////////////////////////////////////////////////B
      IF(DELBAN.EQ.1) THEN
          close(iou6,status="delete")
          close(iou7,status="delete")
          close(iou8,status="delete")

      ELSE
          close(iou6,status="keep")
          close(iou7,status="keep")
          close(iou8,status="keep")
      END IF


 9000 continue

      RETURN
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      END SUBROUTINE BANDIT
! E////////////////////////////////////////////////////////////////////E

! ##################################################################################################################################
      SUBROUTINE BRIGIT(INV,II3,INT,ILD,IER)
C
C     THIS ROUTINE SORTS THE ORIGINAL GRID NUMBERS AND OUTPUTS THE LIST
C         IN INT, WHERE INT(I)=THE ITH ORIGINAL GRID NUMBER.
C     ALSO OUTPUT IS ILD, WHERE ILD(I)=SORTED INTERNAL NUMBER
C         CORRESPONDING TO THE UNSORTED BANDIT INTERNAL LABEL I.
C
C     INPUT - INV
C     OUTPUT - ILD,INT
C
      INTEGER II3, INV(2,II3),INT(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER I      ,IER    ,IS     ,J      ,KFAC   ,KFM   ,
     &        KMOD   ,L      ,MAXDEG ,MAXGRD ,MINI   ,NN

      REAL     DUMS
! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
C
C     PERFORM A ROUGH SORT OF THE ORIGINAL GRID NUMBERS.
C
      L=0
      KFAC=-1
   20 KFAC=KFAC+1
      MINI=99999999
      KFM=KFAC*KMOD
      DO 50 I=1,KMOD
      IF(INV(1,I).GT.KFM) MINI=MIN(MINI,INV(1,I))
   50 CONTINUE
      KFAC=(MINI-1)/KMOD
      DO 80 I=1,KMOD
      IS=INV(1,I)
      IF(IS.LE.(KFAC*KMOD).OR.IS.GT.(KFAC+1)*KMOD)GO TO 80
      L=L+1
      INT(L)=INV(1,I)
   80 CONTINUE
      IF(L.LT.NN)GO TO 20
C
C     COMPLETE THE SORTING OF THE ORIGINAL GRID NUMBERS.
C
      CALL SORT(INT,NN)
C
C     DETERMINE CORRESPONDENCE (ILD) BETWEEN NORIG AND INT ARRAYS.
C
      DO 40 I=1,NN
      J=INT(I)
      L=INTERN(J,INV,II3,IER)
      IF(IER.GT.0) RETURN
C     J IS DEFINED ABOVE TO ALLOW EXECUTION ON UNIVAC IN OVER 65K.
   40 ILD(L)=I
      RETURN
      END SUBROUTINE BRIGIT

! ##################################################################################################################################
      SUBROUTINE CASE(IER)
C
C     READ EXECUTIVE AND CASE CONTROL DECKS.
C
      INTEGER KA(80)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IER    ,IPARAM ,KT     ,L      ,MA     ,MB     ,
     &         NUM

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      KT=0
      REWIND IOU8
C      WRITE(IOU6,10)
C   10 FORMAT(/,'Echo of Data Through BEGIN BULK Card:')
C
C     INITIALIZE PARAMETERS.
C
      CALL DOLLAR(KA,1,KT,IER)
      IF(IER.GT.0) RETURN
C
C     READ CARD.
C
   20 READ(IOU5,30,END=100) KA
   30 FORMAT(80A1)
      KT=KT+1
      L=LENCAS(KA,80)
C      WRITE(IOU6,60) KT,(KA(I),I=1,L)
C60    FORMAT(I8,'- ',80A1)
C      WRITE(IOU8,30) (KA(I),I=1,L)
C
C     IF COLUMN 1 IS $, PROCESS $-CARD.
C
      IF(KA(1).EQ.MB(1)) then
         CALL DOLLAR(KA,2,KT,IER)
         IF(IER.GT.0) RETURN
         GO TO 20
      end if
C
C     IF THE CARD IS NOT BEGIN BULK, GO BACK AND READ ANOTHER CARD.
C
      IF(NBULK(KA).EQ.0) GO TO 20
C
C     CHECK FOR ILLEGAL PARAMETERS.
C
   80 CALL DOLLAR(KA,3,KT,IER)
      IF(IER.GT.0) RETURN
      RETURN
C
C     END-OF-FILE ENCOUNTERED
C
100   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE CASE

! ##################################################################################################################################
      SUBROUTINE COREKO
C
C     PRINT SUMMARY OF CORE ALLOCATION.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IBYTE  ,IFL    ,KDIM   ,KMOD   ,KORE   ,LL     ,MAXDEG ,
     &         MAXGRD ,NBITIN ,NBYTE  ,NW

      REAL     DUM
! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ NBITIN,KORE,IFL,DUM(2),NW,NBYTE,IBYTE,KDIM
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      LL=IFL-KORE
      WRITE(IOU6,10) KORE,MAXGRD,MAXDEG,KDIM
10    FORMAT(/,'Working Storage:'/
     - 5X,'Length of Working Storage',I10,' words'/
     - 5X,'Grid Point Limit         ',I10/
     - 5X,'Nodal Degree Limit       ',I10/
     - 5X,'$DIMENSION Value         ',I10)
      IF(NW.GT.1) WRITE(IOU6,20) NW
20    FORMAT(5X,'Packing Density',I20,' integers/word')
      RETURN
      END SUBROUTINE COREKO

! ##################################################################################################################################
      SUBROUTINE CUTHIL(NT,NUM,NOM,IO,IP,IG,II1,IC,IDEG,IDIS,IW,
     -  NEW,ICC,ILD,IPP,JUMP,NODESL,KORDIM,IER)
C
C THIS IS THE EXECUTIVE FOR THE CUTHILL-MCKEE GRID POINT RENUMBERING
C      STRATEGY.   THE PRINCIPAL INPUTS ARE THE CONNECTIVITY MATRIX IG
C      AND THE NUMBER OF GRID POINTS (NODES) NN.
C
C INPUT -- NT,NUM,NOM,IO,IP,IG,II1,NN,MAXGRD,ILD,NBITIN,ISTART,ISTA
C OUTPUT -- NEW,ILD,MM,IH0,IHE,KORIG,KNEW,NCM
C SCRATCH -- IC,IDEG,IDIS,IW,ICC,IPP
C SET FOLLOWING DIMENSIONS IN CALLING PROGRAM --
C   IG(II1,M),IC(L),IDEG(L),IDIS(L),IW(L),NEW(L),ICC(L),ILD(L),IP(M)
C        L=MAXGRD EXCEEDS NUMBER OF GRID POINTS
C        II1=MAXGRD/(PACKING DENSITY IN INTEGERS/WORD)
C        M EXCEEDS MAX NODAL DEGREE
C
C NT=MAX NUMBER OF STARTING NODES TO BE CONSIDERED
C NUM AND NOM GIVE THE FRACTION OF THE RANGE FROM MIN DEGREE TO MAX
C        DEGREE TO CONSIDER FOR STARTING NODES
C IO=STRATEGY OPTION - - -
C          1=RMS WAVEFRONT
C          2=BANDWIDTH
C          3=PROFILE
C          4=WAVEFRONT (MAX)
C IP=PRINTING OPTION (0 FOR NO PRINTED OUTPUT FROM CUTHILL)
C IG(I,J) CONTAINS THE GRID POINT LABEL FOR THE JTH NODE ADJACENT TO
C      NODE I  (THE CONNECTIVITY MATRIX).   THE CONNECTION OF A NODE
C      TO ITSELF IS NOT LISTED.
C II1=ROW DIMENSION OF IG
C NN=NUMBER OF GRID POINTS (NODES)
C MM=COLUMN DIMENSION OF IG ON INPUT, MAX NODAL DEGREE ON OUTPUT
C MAXGRD=EFFECTIVE ROW DIMENSION OF IG (NEGLECTING INTEGER PACKING)
C LINE=NUMBER OF PRINTED LINES PER PAGE
C NBITIN=NUMBER OF BITS PER INTEGER (FOR PACKING)
C ISTART(I)=ITH STARTING NODE SELECTED BY USER
C ISTA=NUMBER OF STARTING NODES IN ISTART
C NEW(I)=OLD LABEL FOR GRID POINT NOW LABELLED I
C ILD(I)=NEW LABEL FOR GRID POINT ORIGINALLY LABELLED I
C        ILD AND NEW ARE INVERSES
C ILD MUST BE INPUT TO CUTHILL TO INDICATE AN INITIAL SEQUENCE.
C           NORMALLY, ON INPUT, SET ILD(I)=I FOR ALL I.
C JUMP=1 IF RESEQUENCING ATTEMPTS RESULT IN NO IMPROVEMENT
C     =0 OTHERWISE.
C IH0=ORIG PROFILE
C IHE=NEW PROFILE
C KORIG=ORIG BANDWIDTH
C KNEW=NEW BW
C NCM=NUMBER OF COMPONENTS
C     NODESL IS SCRATCH SPACE.
C IN CALLING PROGRAM, TRY  CALL CUTHIL(80,1,2,2,1, . . . )
C
C    THE FOLLOWING SUBROUTINES WERE WRITTEN BY E. CUTHILL AND
C    J. MCKEE OF NSRDC - - -
C        DEGREE,DIAM,IDIST,KOMPNT,MAXBND,MAXDGR,MINDEG,RELABL,CUTHILL
C     CUTHILL WAS MODIFIED BY G. C. EVERSTINE, DTRC.
C
      INTEGER II1, KORDIM
      DIMENSION IG(II1,1),IC(*),IDEG(*),IDIS(*),IW(*)
      INTEGER NEW(*),ICC(*),ILD(*),IPP(*),NODESL(KORDIM),SUMW
      REAL IM1,IM2
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &        IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &        IOU19  ,IOU20

      INTEGER  I      ,IAJDIM ,IB     ,IC     ,IDEG   ,IDEM   ,IDEM1  ,
     &         IDIM   ,IDIS   ,IER    ,IG     ,IH     ,IH0    ,IHE    ,
     &         ij     ,IO     ,IP     ,IPARAM ,IS     ,ISTA   ,
     &         ISTART ,IW     ,IWALL

      INTEGER  J      ,JMAX   ,JUMP

      INTEGER  K      ,K2     ,KMOD   ,KNEW ,KORIG

      INTEGER  L

      INTEGER  M      ,MA     ,MAD    ,MAXD   ,MAXDEG ,MAXGRD ,MAXLEV ,
     &         MAXW   ,MAXW0  ,MAXW1  ,MEDIAN ,MI     ,MM     ,MODD

      INTEGER  NBITIN ,NC     ,NCM    ,NL     ,NN     ,NNODE  ,NOM    ,
     &                 NT     ,NUM

      REAL     AVERW  ,BRMS   ,BRMS0  ,BRMS1  ,CRIT1  ,CRIT2  ,
     &         DUMBB  ,DUMD   ,DUML   ,DUMO   ,DUMS   ,RMS    ,
     &         RMS0   ,RMS1   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,IH,IB,DUMS(5)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /D/ KORIG,KNEW,IH0,IHE,NCM,DUMD(5)
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /B/ IPARAM(20)
      COMMON /DOL/ ISTART(100),DUMO(100)
      COMMON /DOLL/ IDIM,ISTA,DUML(3)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      CALL TIMER(TA,IWALL,0,IOU6)
C     SET UP SCRATCH SPACE NODESL.
      IDEM=KORDIM/4
      K2=IDEM+1
      IAJDIM=3 *IDEM
      JUMP=0
C     DETERMINE THE DEGREE OF EACH NODE.
      CALL DEGREE(IG,II1,IDEG)
C     DETERMINE THE NUMBER OF COMPONENTS, NCM.
      NCM=KOMPNT(IG,II1,IC,IDEG,IW,ICC)
C     DETERMINE THE MAXIMUM DEGREE OF ANY NODE.
      MAXD=MAXDGR(0,IC,IDEG)
      MM=MAXD
C INITIALIZE NEW ARRAY FROM THE ILD ARRAY.
C   ILD MUST BE INPUT TO CUTHILL.
      DO 30 I=1,NN
      K=ILD(I)
   30 NEW(K)=I
C   COMPUTE ORIGINAL BANDWIDTH AND PROFILE.
C   COMPUTE ORIGINAL WAVEFRONT AND ACTIVE COLUMN DATA.
      CALL WAVEY(IG,II1,ILD,NEW,0,IC,IW,IS,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
      MAXW0=MAXW
      RMS0=RMS
      BRMS0=BRMS
      KORIG=IS
C    IH0=ORIGINAL PROFILE, IS=ORIGINAL BW.
      IH0=IH
C COMPUTE NODAL DEGREE STATISTICS.
      CALL DIST(IDEG,IPP,IP,MEDIAN,MODD)
C     IF REQUESTED, PRINT INTERNAL NUMBER CONNECTION TABLE.
      IF(IP.ne.0) CALL STABLE(IG,II1,IC,IDEG,ILD,IPP)
      WRITE(IOU6,29)
   29 FORMAT(/,'Before Resequencing:')
      WRITE(IOU6,51) IS,IH,MAXW,AVERW,RMS,BRMS
   51 FORMAT(5X,'Bandwidth',I18/5X,'Profile',I20/
     -       5X,'Max Wavefront',I14/5X,'Avg Wavefront',F14.3/
     -       5X,'RMS Wavefront',F14.3/5X,'RMS Bandwidth',F14.3)
      IF(ip.ne.0.and.ISTA.gt.0) then
         WRITE(IOU6,701)
701      FORMAT('Starting nodes supplied by user:')
         WRITE(IOU6,100) (ISTART(I),I=1,ISTA)
      end if
C INITIALIZE ILD AND NEW ARRAYS.
      DO I=1,NN
         NEW(I)=0
         ILD(I)=0
      end do
C
C     GENERATE NUMBERING SCHEME FOR EACH COMPONENT, NC.
C
      DO 500 NC=1,NCM
C     DETERMINE THE RANGE OF DEGREES (MI  TO  MAD) OF NODES OF INTEREST.
      MI=MINDEG(NC,IC,IDEG)
      MAD=MI
      IF(NOM.EQ.0) GO TO 87
   90 MA=MAXDGR(NC,IC,IDEG)
      MAD=MI+((MA-MI)*NUM)/NOM
C     MAKE SURE MAD DOES NOT EXCEED MEDIAN.
      MAD=MIN(MAD,MEDIAN-1)
      MAD=MAX(MAD,MI)
C     DETERMINE BANDWIDTH OR SUM CRITERION FOR EACH NODE MEETING
C        SPECIFIED CONDITION.
   87 CONTINUE
      CALL DIAM(NC,MAD,NL,NODESL,IDEM,MAXLEV,IG,II1,IC,IDEG,IDIS,
     -  IW,ICC)
      IF(IP.EQ.0) GO TO 67
      WRITE(IOU6,39) NC,MAD
   39 FORMAT(/,'Component',I5,', Max Degree Used',I6)
C     COMPUTE THE NUMBER OF NODES IN THIS COMPONENT AND WRITE OUT.
      NNODE=ICC(NC+1)-ICC(NC)
      WRITE(IOU6,41) NNODE
   41 FORMAT('Number of nodes in this component',I8)
      WRITE(IOU6,59) MAXLEV
   59 FORMAT('Starting nodes for minmax number of nodes per level',I6)
      if(NL.gt.0) WRITE(IOU6,100) (NODESL(J),J=1,NL)
  100 FORMAT(10I7)
   67 CONTINUE
      IF(ISTA.LE.0) GO TO 760
      IDEM1=IDEM-1
      M=0
      DO 750 I=1,ISTA
         J=ISTART(I)
         IF(IC(J).NE.NC) GO TO 750
         M=M+1
         DO K=1,IDEM1
            L=IDEM+1-K
            NODESL(L)=NODESL(L-1)
         end do
         NODESL(1)=J
  750 CONTINUE
      NL=MIN(NL+M,IDEM)
      CALL FIXIT(NODESL,NL)
  760 CONTINUE
      IF(IP.EQ.0) GO TO 63
      IF(ISTA.LE.0) GO TO 63
      WRITE(IOU6,730)
  730 FORMAT('Merged list of starting nodes supplied by user',
     -       ' and by BANDIT -')
      WRITE(IOU6,100) (NODESL(I),I=1,NL)
   63 CONTINUE
      JMAX=MIN(NT,NL)
      JMAX=MAX(JMAX,1)
      IM1=1.E8
      IM2=IM1
C
C  CHECK SEQUENCE FOR EACH STARTING NODE SELECTED.
C
      ij=1
      DO 400 J=1,JMAX
      CALL RELABL(1,NODESL(J),IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,
     -            ILD,NODESL(K2),IAJDIM,IER)
      IF(IER.GT.0) RETURN
C  COMPUTE NEW BANDWIDTH,PROFILE,WAVEFRONT DATA.
      CALL WAVEY(IG,II1,ILD,NEW,NC,IC,IW,IB,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
C     IB=BANDWIDTH, IH=PROFILE.
      IF(IP.EQ.0) GO TO 70
      WRITE(IOU6,69) NODESL(J),IB,IH,MAXW,RMS
   69 FORMAT('Starting Node',I6,', Band',I6,', Profile',I8,
     - ', Max W',I6,', RMS W',F9.3 )
   70 CONTINUE
      GO TO (205,210,215,220), IO
  205 CRIT1=RMS
      CRIT2=IH
      GO TO 71
  210 CRIT1=IB
      CRIT2=IH
      GO TO 71
  215 CRIT1=IH
      CRIT2=IB
      GO TO 71
  220 CRIT1=MAXW
      CRIT2=RMS
      GO TO 71
   71 CONTINUE
      IF(IM1-CRIT1) 400,350,300
  300 IM1=CRIT1
      IM2=CRIT2
      IJ=J
      GO TO 400
  350 IF(IM2.LE.CRIT2) GO TO 400
      IM2=CRIT2
      IJ=J
400   CONTINUE
C
C   RECOMPUTE SEQUENCE FOR STARTING NODE WHICH IS BEST FOR CRITERION
C       SELECTED.
C
      CALL RELABL(1,NODESL(IJ),IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,
     -            ILD,NODESL(K2),IAJDIM,IER)
      IF(IER.GT.0) RETURN
C
  500 CONTINUE
C
C
C  DETERMINE NODES OF ZERO DEGREE AND STACK LAST.
      CALL STACK(IDEG,NEW,ILD,IW)
C   COMPUTE BANDWIDTH, PROFILE AND WAVEFRONT DATA.
      CALL WAVEY(IG,II1,ILD,NEW,0,IC,IW,IB,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
C
      WRITE(IOU6,705)
  705 FORMAT(/,'After resequencing by Reverse Cuthill-McKee (CM):')
      WRITE(IOU6,51) IB,IH,MAXW,AVERW,RMS,BRMS
C
C   CHECK CM LABELING AGAINST ORIGINAL LABELING TO SEE IF BETTER.
C
      GO TO (255,260,265,270), IO
  255 IM1=RMS0
      IM2=IH0
      CRIT1=RMS
      CRIT2=IH
      GO TO 711
  260 IM1=IS
      IM2=IH0
      CRIT1=IB
      CRIT2=IH
      GO TO 711
C IB = BANDWIDTH,  IH= PROFILE.
  265 IM1=IH0
      IM2=IS
      CRIT1=IH
      CRIT2=IB
      GO TO 711
  270 IM1=MAXW0
      IM2=RMS0
      CRIT1=MAXW
      CRIT2=RMS
      GO TO 711
  711 CONTINUE
      IF(CRIT1-IM1) 715,742,744
  742 IF(CRIT2.LT.IM2) GO TO 715
C   IF NO IMPROVEMENT RETURN TO ORIGINAL SEQUENCE.
  744 IB=IS
      IH=IH0
      MAXW=MAXW0
      RMS=RMS0
      BRMS=BRMS0
      JUMP=1
      DO 713 I=1,NN
      ILD(I)=I
  713 NEW(I)=I
C
  715 CONTINUE
C   SET FINAL VALUES OF B , P , RMS , W .
      KNEW = IB
      IHE=IH
      MAXW1=MAXW
      RMS1=RMS
      BRMS1=BRMS
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
c     IF(TB.GT.1.E-5) WRITE(IOU6,610) TB
      WRITE(IOU6,610) TB
  610 FORMAT(5X,'CP time',F20.3)
      RETURN
      END SUBROUTINE CUTHIL

! ##################################################################################################################################
      SUBROUTINE DEGREE(IG,II1,IDEG)
C
C     SET UP THE IDEG ARRAY CONTAINING THE DEGREE OF EACH NODE STORED
C     IN THE IG ARRAY.
C     IDEG(I) = DEGREE OF NODE I
C
      INTEGER II1
      DIMENSION IG(II1,1),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IDEG   ,IG    ,J      ,KMOD   ,MAXDEG ,
     &         MAXGRD ,MM     ,NN     ,NBITIN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      DO 100 I=1,NN
      IDEG(I)=0
      DO 80 J=1,MM
      IF(IG(I,J)) 100,100,50
CPACK IF(IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)) 100,100,50
   50 IDEG(I)=IDEG(I)+1
   80 CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE DEGREE

! ##################################################################################################################################
      SUBROUTINE DGREE(NDSTK,NR,NDEG,IOLD,IBW1,IPF1)
C
C     DGREE COMPUTES THE DEGREE OF EACH NODE IN NDSTK AND STORES
C     IT IN THE ARRAY NDEG.  THE BANDWIDTH AND PROFILE FOR THE ORIGINAL
C     OR INPUT RENUMBERING OF THE GRAPH IS COMPUTED ALSO.
C
C     COMPUTE MAXIMUM DEGREE MM AND STORE IN IDEG.
C
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IBW1   ,IDEG   ,IDIF   ,IDPTH  ,IOLD   ,IPF1   ,
     &         IRW    ,ITST   ,J      ,KMOD

      INTEGER  MAXDEG ,MAXGRD ,MM     ,N      ,NBITIN ,NDEG   ,NDSTK  ,
     &         NN     ,NR

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /GRA/ N,IDPTH,IDEG
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /S/ NN,MM,DUMS(7)
! B/2020-12-12/////////////////////////////////////////////////////////B
Cxxxx DIMENSION NDSTK(NR,1),NDEG(*),IOLD(*) gfortran error. Need > 1
      DIMENSION NDSTK(NR,*),NDEG(*),IOLD(*)
! E////////////////////////////////////////////////////////////////////E
      IBW1=0
      IPF1=0
      IDEG=MM
      MM=0
      DO 100 I=1,N
        NDEG(I)=0
        IRW=0
        DO 80 J=1,IDEG
          ITST=NDSTK(I,J)
CPACK     ITST=IUNPK(NDSTK,MAXGRD*(J-1)+I,NBITIN)
          IF(ITST.EQ.0) GO TO 90
   50     NDEG(I)=NDEG(I)+1
          IDIF=IOLD(I)-IOLD(ITST)
          IF(IRW.LT.IDIF) IRW=IDIF
          MM=MAX(MM,J)
   80   CONTINUE
   90   IPF1=IPF1+IRW
        IF(IRW.GT.IBW1) IBW1=IRW
  100 CONTINUE
      IDEG=MM
C     INCLUDE DIAGONAL TERMS IN BANDWIDTH AND PROFILE
      IBW1=IBW1+1
      IPF1=IPF1+N
      RETURN
      END SUBROUTINE DGREE

! ##################################################################################################################################
      SUBROUTINE DIAM(NC,MAXDG,NL,NODESL,IDEM,MAXLEV,
     -                IG,II1,IC,IDEG,IDIS,IW,ICC)
C
C     DETERMINE NL STARTING POINTS AND STORE IN NODESL.
C
      INTEGER II1
      DIMENSION IG(II1,*),IDIS(*),IW(*),ICC(*),IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IC     ,ICC    ,IDEG   ,IDEM   ,IDIS   ,IG     ,
     &         IW     ,KMOD   ,MAXDEG ,MAXDG  ,MAXGRD ,MAXLEV ,
     &         MD     ,ML     ,NBITIN ,NC     ,NL     ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER NODESL(*)
      NL=0
      MAXLEV=600000
      DO 100 I=1,NN
      IF((NC-IC(I)).NE.0) GO TO 100
40    IF((MAXDG-IDEG(I)).LT.0) GO TO 100
  105 MD=IDIST(I,ML,MAXLEV,IG,II1,IC,IDEG,IDIS,IW,ICC)
      IF(MD) 115,115,56
   56 IF(ML-MAXLEV)58,64,100
   58 MAXLEV=ML
      NL=1
      NODESL(1)=I
      GO TO 100
   64 IF(NL.GE.IDEM) GO TO 100
      NL=NL+1
      NODESL(NL)=I
  100 CONTINUE
  110 RETURN
  115 ML=1
      NODESL(1)=I
      MAXLEV=0
      RETURN
      END SUBROUTINE DIAM

! ##################################################################################################################################
      SUBROUTINE DIST(IDEG,HIST,IP,MEDIAN,MODD)
C
C     COMPUTE AND PRINT THE DISTRIBUTION OF NODAL DEGREES WITH MEDIAN
C        AND MODE.
C
C     IDEG(I) = DEGREE OF NODE I
C     HIST(I) = NUMBER OF NODES OF DEGREE I
C     IP      = PRINT OPTION PARAMETER (IF 0, NO PRINTED OUTPUT)
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IP     ,ISUM   ,K      ,MAXI   ,MEDIAN ,MM     ,
     &         MM1    ,MODD   ,NN     ,NN2

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER IDEG(*),HIST(*)
C
C     COMPUTE HISTOGRAM.
C
      MM1=MM+1
      DO 10 I=1,MM1
   10 HIST(I)=0
      DO 20 I=1,NN
      K=IDEG(I)+1
   20 HIST(K)=HIST(K)+1
C
C     COMPUTE MODE (MODD).
C
      MODD=0
      MAXI=0
      DO 25 I=1,MM1
      K=HIST(I)
      IF(K.LE.MAXI) GO TO 25
      MAXI=K
      MODD=I-1
   25 CONTINUE
      IF(IP.EQ.0) GO TO 60
C
C     PRINT HISTOGRAM.
C
      WRITE(IOU6,30)
   30 FORMAT(/,'Distribution of Nodal Degrees:'//5X,
     -       'Degree  Number  Cum. Total')
      ISUM=0
      DO 40 I=1,MM1
      ISUM=ISUM+HIST(I)
      K=I-1
   40 WRITE(IOU6,50) K,HIST(I),ISUM
   50 FORMAT(3X,2I8,I12)
C
C     COMPUTE CUMULATIVE DISTRIBUTION.
C
   60 DO 70 I=2,MM1
   70 HIST(I)=HIST(I)+HIST(I-1)
C     COMPUTE MEDIAN.
      NN2=NN/2
      DO 80 I=1,MM1
      IF(HIST(I).GT.NN2) GO TO 90
   80 CONTINUE
   90 MEDIAN=I-1
      IF(IP.NE.0) WRITE(IOU6,100) MEDIAN,MODD
  100 FORMAT(/5X,'Median',I6/5X,'  Mode',I6)
      RETURN
      END SUBROUTINE DIST

! ##################################################################################################################################
      SUBROUTINE DOLLAR(KA,JUMP,KT,IER)
C
C     INTERPRET DOLLAR SIGN CARD.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &        IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &        IOU19  ,IOU20

      INTEGER  I      ,IADD   ,IDIM   ,IER    ,IFIR   ,IFL    ,IFLD   ,
     &         IGDEG  ,IGNORE ,IIG    ,IPARAM ,ISTA   ,ISTART ,ITYPE  ,
     &         J      ,JUMP   ,K      ,KDIM   ,KMOD   ,KORE   ,KT     ,
     &         LINES                                                  ,
     &         MA     ,MAXDEG ,MAXGRD ,MAXI   ,MB     ,MDIM   ,ME     ,
     &         NBITIN ,NCON   ,NEDGE  ,NELEM  ,NGRID  ,NIP    ,NTYPE  ,
     &         NUM

      REAL     DMY    ,DUM    ,DUMS   ,DUMY

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /S/ DMY(4),LINES,NEDGE,IADD,DUMS,DUMY
      COMMON /BITS/ NBITIN,KORE,IFL,NGRID,DUM(4),KDIM
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),
     -              NELEM(160),MDIM
      INTEGER vype,TYPE,WYPE,KA(80),IP(35)
      DATA MAXI/35/
C
      GO TO (10,30,80), JUMP
C
C     INITIALIZE PARAMETERS.
C
C     IF LOOP PARAMETER IS ALREADY SET ($LOOP YES), INDICATING THAT
C     MULTIPLE DATA SETS ARE TO BE PROCESSED IN ONE EXECUTION, DO NOT
C     RESET IPARAM(1).
C
10    J=1
      IF(IPARAM(1).EQ.4) J=2
      DO 20 I=J,20
20    IPARAM(I)=0
      IDIM=100
      IIG=0
      ISTA=0
      IGDEG=0
      IADD=0
      RETURN
C
C* PARAMETERS - - - -
C
C        KEYWORD 1           KEYWORD 2
C        1 - LOOP            1 - SEQGP
C        2 - TABLE           2 - ALL
C        3 - METHOD          3 - NO
C        4 - MPC             4 - YES
C        5 - SEQUENCE        5 - MIN
C        6 - CRITERION       6 - MAX
C        7 - SPRING          7 - CM
C        8 - NASTRAN         8 - GPS
C        9 - ELEMENTS        9 - BOTH CM AND GPS
C       10 - PRINT          10 - BANDWIDTH
C       11 - FRONTAL        11 - PROFILE
C                           12 - RMS WAVEFRONT
C                           13 - WAVEFRONT (MAX)
C
C     REMAINING PARAMETERS UP TO 20 IN LEFT COLUMN ARE CURRENTLY UNUSED
C
C*   OTHER KEYWORDS - -   IGNORE,GRID,START,DEGREE,INSERT,LINES,
C             DIMENSION,ADD,APPEND
C
C
C     LOOK FOR FIRST KEYWORD.
C
   30 ITYPE=0
C     CHECK COLUMN 2 FOR BLANK.
      IF(KA(2).EQ.MB(2)) RETURN
C     SEE BLOCK DATA FOR ALPHABET KEY.
!    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
!    A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
      IF(KA(2).EQ.MA(12).AND.KA(3).EQ.MA(15)) ITYPE=1      ! LO $LOOP
      IF(KA(2).EQ.MA(20).AND.KA(3).EQ.MA( 1)) ITYPE=2      ! TA $TABLE
      IF(KA(2).EQ.MA(13).AND.KA(3).EQ.MA( 5)) ITYPE=3      ! ME $METHOD
      IF(KA(2).EQ.MA(13).AND.KA(3).EQ.MA(16)) ITYPE=4      ! MP $MPC
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA( 5)) ITYPE=5      ! SE $SEQUENCE
      IF(KA(2).EQ.MA( 3).AND.KA(3).EQ.MA(18)) ITYPE=6      ! CR $RITERION
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA(16)) ITYPE=7      ! SP $SPRING
      IF(KA(2).EQ.MA(14).AND.KA(3).EQ.MA( 1)) ITYPE=8      ! NA $
      IF(KA(2).EQ.MA( 5).AND.KA(3).EQ.MA(12)) ITYPE=9      ! EL $ELEMENTS
      IF(KA(2).EQ.MA(16).AND.KA(3).EQ.MA(18)) ITYPE=10     ! PR $PRINT
      IF(KA(2).EQ.MA( 6).AND.KA(3).EQ.MA(18)) ITYPE=11     ! FR $FRONTAL
      IF(ITYPE.GT.0) GO TO 35
      IF(KA(2).EQ.MA( 7).AND.KA(3).EQ.MA(18)) GO TO 120    ! GR $GRID
      IF(KA(2).EQ.MA( 9).AND.KA(3).EQ.MA( 7)) GO TO 100    ! IG $IGNORE
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA(20)) GO TO 130    ! ST $START
      IF(KA(2).EQ.MA( 4).AND.KA(3).EQ.MA( 5)) GO TO 150    ! DE $DEGREE
      IF(KA(2).EQ.MA( 9).AND.KA(3).EQ.MA(14)) GO TO 160    ! IN $INSERT
      IF(KA(2).EQ.MA(12).AND.KA(3).EQ.MA( 9)) GO TO 180    ! LI $LINES
      IF(KA(2).EQ.MA( 4).AND.KA(3).EQ.MA( 9)) GO TO 190    ! DI $DIMENSION
      IF(KA(2).EQ.MA( 1).AND.KA(3).EQ.MA( 4)) GO TO 230    ! AD $ADD
      IF(KA(2).EQ.MA( 1).AND.KA(3).EQ.MA(16)) GO TO 240    ! AP $APPEND
      RETURN
C
C     LOOK FOR SECOND KEYWORD.
C
   35 DO I=4,70
         IF(KA(I).EQ.MB(2)) GO TO 50
      end do
      RETURN
   50 K=I+1
      DO I=K,71
         IF(KA(I).NE.MB(2)) GO TO 70
      end do
      RETURN
C
70    J=0
      IF(KA(I).EQ.MA(19).AND.KA(I+1).EQ.MA( 5)) J=1        ! SE
      IF(KA(I).EQ.MA( 1).AND.KA(I+1).EQ.MA(12)) J=2        ! AL
      IF(KA(I).EQ.MA(14).AND.KA(I+1).EQ.MA(15)) J=3        ! NO
      IF(KA(I).EQ.MA(25).AND.KA(I+1).EQ.MA( 5)) J=4        ! YE
      IF(KA(I).EQ.MA(13).AND.KA(I+1).EQ.MA( 9)) J=5        ! MI
      IF(KA(I).EQ.MA(13).AND.KA(I+1).EQ.MA( 1)) J=6        ! MA
      IF(KA(I).EQ.MA( 3).AND.KA(I+1).EQ.MA(13)) J=7        ! CM
      IF(KA(I).EQ.MA( 7).AND.KA(I+1).EQ.MA(16)) J=8        ! GP
      IF(KA(I).EQ.MA( 2).AND.KA(I+1).EQ.MA(15)) J=9        ! BO
      IF(KA(I).EQ.MA( 2).AND.KA(I+1).EQ.MA( 1)) J=10       ! BA
      IF(KA(I).EQ.MA(16).AND.KA(I+1).EQ.MA(18)) J=11       ! PR
      IF(KA(I).EQ.MA(18).AND.KA(I+1).EQ.MA(13)) J=12       ! RM
      IF(KA(I).EQ.MA(23).AND.KA(I+1).EQ.MA( 1)) J=13       ! WA
      IF(J.GT.0) IPARAM(ITYPE)=J
      RETURN
C
C     CHECK FOR ILLEGAL PARAMETERS AND SET TO DEFAULTS.
C
80    IF(IPARAM(1).NE.4)                    IPARAM(1)=3
      IF(IPARAM(2).NE.4)                    IPARAM(2)=3
      IF(IPARAM(3).NE.7.AND.IPARAM(3).NE.9) IPARAM(3)=8
      IF(IPARAM(4).NE.4)                    IPARAM(4)=3
      IF(IPARAM(5).NE.3)                    IPARAM(5)=4
! B////////////////////////////////////////////////////////////////////B
! Change default for $CRITERION from RMS to BAND


      IF(IPARAM( 6).NE.11.AND.IPARAM( 6).NE.12.AND.IPARAM(6).NE.13)
     -                                        IPARAM( 6)=10
! E////////////////////////////////////////////////////////////////////E
      IF(IPARAM(7).NE.4)                    IPARAM(7)=3
      IF(IPARAM(8).NE.4)                    IPARAM(8)=3
      IF(IPARAM(9).NE.4)                    IPARAM(9)=3
      IF(IPARAM(10).NE.6)                   IPARAM(10)=5
      IF(IPARAM(11).NE.4)                   IPARAM(11)=3
C     INVOKE CONNECTION TABLE OPTION IF SPRINGS ARE REQUESTED.
      IF(IPARAM(7).EQ.4) IPARAM(2)=4
      KDIM=MIN(KDIM,NGRID)
C     INVOKE RMS CRITERION FOR FRONTAL JOB.
      IF(IPARAM(11).EQ.4) IPARAM(6)=12
c     Force CM to be run if $print max is enabled (to get component list)
      if(iparam(10).eq.6.and.iparam(3).eq.8) iparam(3)=9
C     TO PREVENT USE OF A PARTICULAR $ CARD (SUCH AS $SPRING), ENFORCE
C          THE NO OPTION HERE.
C
      RETURN
C
C     $IGNORE G1,G2, ...     (NODES TO IGNORE)
C
  100 CALL READIT(KA(4),MAXI,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      I=IIG
      IIG=IIG+NIP
      IF(IIG.LE.IDIM) GO TO 105
      WRITE(IOU6,102)
  102 FORMAT(/,'Fatal Error.  Too many points on $-card.')
      call finish(6,IER)
      RETURN
  105 DO 110 J=1,NIP
      K=I+J
  110 IGNORE(K)=IP(J)
      RETURN
C
C     $GRID N    (UPPER BOUND ON NUMBER OF GRID POINTS)
C
  120 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.EQ.0) RETURN
      NGRID=MAX(IP(1),10)
      KDIM=MAX(KDIM,NGRID/10)
      KDIM=MIN(KDIM,NGRID)
      GO TO 195
C
C     $START G1,G2, ...     (USER-SELECTED STARTING NODES FOR CM METHOD)
C
  130 CALL READIT(KA(4),MAXI,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      I=ISTA
      ISTA=ISTA+NIP
      IF(ISTA.LE.IDIM) GO TO 135
      WRITE(IOU6,102)
      call finish(6,IER)
      RETURN
  135 DO 140 J=1,NIP
      K=I+J
  140 ISTART(K)=IP(J)
      RETURN
C
C     $DEGREE N     (TO IGNORE NODES OF DEGREE EXCEEDING N)
C
  150 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.GT.0) IGDEG=IP(1)
      RETURN
C
C     $INSERT       OR      $INSERT N
C
  160 CALL INSERT(KA,KT)
      RETURN
C
C     $LINES N     (NO. OF LINES PER PAGE)
C
  180 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      LINES=MAX(IP(1),10)
      RETURN
C
C     $DIMENSION N    (DIMENSION OF SOME SCRATCH ARRAYS)
C
  190 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      KDIM=MAX(IP(1),10)
195   KDIM=MIN(KDIM,KORE/4)
      RETURN
C
C     $ADD N     (TO ADD N TO NEW SEQUENCE NUMBERS)
C
  230 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.GT.0) IADD=IP(1)
      RETURN
C
C     $APPEND     CNAME     NCON     IFLD
C
C     USER-DEFINED CONNECTION CARD, WHERE
C
C     CNAME=NAME OF CONNECTION CARD (E.G., CBAR) LEFT-ADJUSTED STARTING
C          IN COLUMN 9
C     NCON=NUMBER OF CONNECTIONS ON CARD (I.E., NODES IN ELEMENT)
C     IFLD=NASTRAN FIELD NUMBER ON PARENT CARD IN WHICH FIRST
C          CONNECTION APPEARS
C     NCON AND IFLD MAY APPEAR ANYWHERE IN COLUMNS 17 - 32 SEPARATED BY
C          ONE OR MORE BLANKS.
C
C     REMARKS  -
C          1. NO LONG-FIELD CONNECTION CARDS MAY BE DEFINED.
C          2. CONNECTIONS MUST BE LISTED CONSECUTIVELY ON PARENT AND
C             CONTINUATION CARDS, IF ANY.
C          3. EACH TIME A $APPEND CARD APPEARS, A NEW ELEMENT IS
C             DEFINED.
C          4. IF CNAME MATCHES AN EXISTING ELEMENT NAME, THE OLD
C             ELEMENT IS REPLACED BY THE NEW ONE.
C          5. THE NUMBER OF ELEMENT CONNECTIONS IS FIXED RATHER
C             THAN VARIABLE.
C          6. IF "$FRONTAL YES" IS SPECIFIED, FIELD 2 MUST CONTAIN THE
C             ELEMENT NUMBER (SINCE ELEMENT NUMBERS ARE NEEDED FOR
C             RESEQUENCING).
C
  240 CALL READIT(KA(17),2,16,IP,NIP)
      IF(NIP.LT.2) GO TO 260
      NCON=IP(1)
      IFLD=IP(2)
      IF(NCON.LT.1) GO TO 260
      IF(IFLD.LT.2) GO TO 260
      IF(IFLD.GT.9) GO TO 260
      REWIND IOU11
      WRITE(IOU11,245) (KA(I),I=9,16)
  245 FORMAT(8A1)
      rewind iou11
      READ(IOU11,250) (KA(I),I=1,3)
  250 FORMAT(A1,A4,A3)
      REWIND IOU11
      IF(KA(1).NE.MA(3)) GO TO 260
      IF(KA(2).EQ.MB(2)) GO TO 260
C     CHECK IF CNAME MATCHES NAME ALREADY IN LIST.
      DO 251 I=19,NTYPE
         IF(I.GE.70.AND.I.LE.78) GO TO 251
         IF(I.EQ.93) GO TO 251
         IF(I.EQ.96) GO TO 251
         IF(I.EQ.97) GO TO 251
         IF(KA(2).EQ.TYPE(I).AND.KA(3).EQ.WYPE(I)) GO TO 252
251   CONTINUE
      NTYPE=NTYPE+1
      I=NTYPE
      GO TO 254
252   WRITE(IOU6,253)
253   FORMAT(10X,'This element replaces another of the same name.')
254   IF(NTYPE.GT.MDIM) GO TO 260
      vype(i)=ka(1)
      TYPE(I)=KA(2)
      WYPE(I)=KA(3)
      ME(I)=10*NCON+IFLD
C     TURN ON SWITCH TO PRINT OUT ELEMENT LIBRARY.
      IPARAM(9)=4
C     SET KDIM SO ENOUGH SPACE IS AVAILABLE.
      KDIM=MAX(KDIM, NCON/4+1)
      IF(KDIM.GT.(KORE/4)) GO TO 260
      I=MDIM-NTYPE
      WRITE(IOU6,255) I
  255 FORMAT(10X,'Space exists to define',I4,' more elements.')
      RETURN
C     ERROR FOUND ON $APPEND CARD.
  260 CONTINUE
      WRITE(IOU6,265)
  265 FORMAT(/,'Fatal Error.  Illegal data on $APPEND card.')
      call finish(6,IER)
      RETURN
C
      END SUBROUTINE DOLLAR

! ##################################################################################################################################
      SUBROUTINE ELTYPE(KA,ITYPE,NCON,IFLD,LOOP,MAXI,LEN,LESSOK,IER)
C
C     DETERMINE BULK DATA CARD ELEMENT TYPE.
C
C     KA     = CONTENTS OF BULK DATA CARD (A1,A4,A3,64A1,A1,A4,A3)
C              (INPUT)
C     ITYPE  = ELEMENT TYPE NUMBER (OUTPUT)
C     NCON   = NUMBER OF CONNECTIONS FOR ELEMENT (OUTPUT)
C              (FOR SOME ELEMENTS  (E.G., CELAS1), THE CONNECTIONS ARE
C              NOT LISTED CONSECUTIVELY, IN WHICH CASE NCON BECOMES THE
C              NUMBER OF FIELDS.  LATER, THE MIDDLE FIELD IS BLANKED
C              OUT.)
C     IFLD   = FIELD NUMBER OF FIRST CONNECTION (OUTPUT)
C     LOOP   = 2 IF TWO ELEMENTS CAN BE DEFINED ON ONE CARD,
C              OTHERWISE 1. (OUTPUT)
C              THE ONLY CARD TYPES FOR WHICH TWO ELEMENTS CAN BE DEFINED
c              ON ONE CARD ARE CDAMP3, CDAMP4, CELAS3, CELAS4, CMASS3,
c              CMASS4, CDAMP4*,CELAS4*, CMASS4*, CRIGDR.
c              CROD(10), CTUBE(11), CVISC(12) were previously included
c              in this list and removed.
C     MAXI   = MAXIMUM NUMBER OF GRID POINTS ON A CONNECTION CARD
C              (=kdim4 from main) (INPUT)
C     LEN    = 1 FOR SHORT FIELD DATA CARDS,
C              2 FOR LONG  FIELD DATA CARDS (OUTPUT)
C     LESSOK = .TRUE. IF IT IS OK FOR THE ELEMENT TO HAVE FEWER THAN
C              THE MAXIMUM NUMBER OF NODES PRESENT ON CONNECTION CARD
C              (E.G., CELAS1 OR CPENTA), OTHERWISE, .FALSE. (OUTPUT)
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IER    ,IFLD   ,IPARAM ,ITYPE  ,J      ,L      ,
     &         LEN    ,LOOP   ,MA     ,MAXI   ,MB     ,MDIM   ,ME     ,
     &         NCON   ,NELEM  ,NTYPE  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),
     -              NELEM(160),MDIM
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER vype,TYPE,WYPE,KA(70),EQGP
      LOGICAL LESSOK
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NCON we will use NCON_array
! instead. Needed so Lahey doesn't complain about shape of NCON being
! different than IP (array) in subr READIT
      INTEGER NCON_array(1)
! E////////////////////////////////////////////////////////////////////E
C
C     SEE BLOCK DATA FOR A LISTING OF THE BANDIT ELEMENT LIBRARY.
C     IT CAN BE LISTED AT EXECUTION TIME WITH  $ELEMENT YES  CARD.
C
! B////////////////////////////////////////////////////////////////////B
!xx   DATA EQGP/'EQGP'/
      DATA EQGP/4HEQGP/
! E////////////////////////////////////////////////////////////////////E
      LOOP=1
      LEN=1
      ITYPE=0
      LESSOK=.FALSE.
      IF(KA(1).NE.MA(3)) GO TO 30
C
C     LOOK FOR CONNECTION CARD.
C
      DO I=4,NTYPE
         IF(KA(2).EQ.TYPE(I).AND.KA(3).EQ.WYPE(I)) GO TO 20
      end do
      RETURN
C
   20 ITYPE=I
      NCON=ME(I)/10
      IFLD=ME(I)-10*NCON
      IF(NCON.eq.0) then
         WRITE(IOU6,24) TYPE(I),TYPE(I)
24       FORMAT(/,'Fatal Error.  A',A4,' card does not precede first C',
     -          A4,' card.')
         call finish(6,IER)
         RETURN
      end if
C
C     SET SPECIAL PARAMETERS LOOP, LEN, AND LESSOK.
C
c     IF(I.GE.10.AND.I.LE.12) LOOP=2
      IF(I.GE.13.AND.I.LE.18) LOOP=2
      IF(I.GE.70.AND.I.LE.72) LOOP=2
      IF(I.EQ.93) LOOP=2
      IF(I.GE.70.AND.I.LE.78) LEN =2
      IF(I.GE.04.AND.I.LE.09) LESSOK=.TRUE.
      IF(I.GE.13.AND.I.LE.18) LESSOK=.TRUE.
      IF(I.GE.70.AND.I.LE.75) LESSOK=.TRUE.
      IF(I.EQ.96.OR.I.EQ.97) LESSOK=.TRUE.
      IF(I.GE.121.AND.I.LE.127) LESSOK=.TRUE.
      RETURN
C
C     LOOK FOR ENDDATA, MPC, MPC*, MPCAX, OR SEQGP.
C
   30 IF(KA(1).EQ.MA( 5).AND.KA(2).EQ.TYPE(1))   ITYPE=1              ! enddata
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(2))   ITYPE=2              ! mpc
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(3))   ITYPE=3              ! mpc*
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(112)) ITYPE=112            ! mpcax
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(129)) ITYPE=129            ! rbar
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(131)) ITYPE=131            ! rbe1
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(133)) ITYPE=133            ! rbe2
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(135)) ITYPE=135            ! rrod
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(137)) ITYPE=137            ! rtrplt
      IF(KA(1).EQ.MA(19).AND.KA(2).EQ.EQGP) GO TO 35                  ! seqgp
      IF(KA(1).EQ.MA( 1)) GO TO 50                                    ! a
      RETURN
C
C     SEQGP CARD FOUND.   IF RESEQUENCING REQUESTED, ABORT JOB.
C
   35 IF(IPARAM(5).EQ.3) RETURN
! B////////////////////////////////////////////////////////////////////B
! Comment out this fatal error, MYSTRAN disallows SEQGP card images when
! it reads the input deck prior to Bandit being called
! E////////////////////////////////////////////////////////////////////E
      RETURN
C
C     LOOK FOR ADUMI CARD AND EXTRACT NCON, THE NUMBER OF GRID POINTS.
C
   50 DO 60 I=58,67
      IF(KA(2).EQ.TYPE(I)) GO TO 70
   60 CONTINUE
      RETURN
! B////////////////////////////////////////////////////////////////////B
70    CALL READIT(KA(4),1,8,NCON_array(1),J)
      NCON = NCON_array(1)
! E////////////////////////////////////////////////////////////////////E
      L = ME(I) - (ME(I)/10) * 10
      ME(I) = 10 * NCON + L
C     THE DEFAULT NUMBER OF NODES FOR CDUMMY IS GIVEN IN BLOCK DATA.
C          IT CAN BE CHANGED WITH AN ADUMMY CARD.
C     THE DEFAULT NUMBER OF NODES FOR CDUMI IS 0.
      IF(NCON.LE.MAXI) RETURN
      WRITE(IOU6,80) (KA(I),I=1,3),MAXI
80    FORMAT(/,'Fatal Error.  Number of grid connections on ',A1,A4,A3,
     - 'exceeds limit of',I7)
      WRITE(IOU6,85)
85    FORMAT('Use $DIM N card, where 4*N exceeds the maximum number'/
     -       'of grid points for any one element.')
      call finish(6,IER)
      RETURN
      END SUBROUTINE ELTYPE

! ##################################################################################################################################
      SUBROUTINE FINISH(KUMP,IER)
C
! B////////////////////////////////////////////////////////////////////B
      USE IOUNT1, ONLY   :  SEQ
! E////////////////////////////////////////////////////////////////////B
C     TERMINATE JOB AFTER FATAL ERROR.
C     THE SUBSEQUENT EXECUTION OF NASTRAN IS PREVENTED BY ERASING UNIT 8.
C
c
c     Kump   Reason
c      1     maxgrd or maxdeg exceeded
c      2     $DIM too small for rigid element
c      3     $DIM exceeded in GPS
c      4     EOF encountered
c      5     $DIM exceeded in CM
c      6     all other reasons
c
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IER    ,IPARAM ,ITIME  ,KUMP   ,NN     ,KDIM   ,
     &         KA     ,NCARD

      REAL     DUM    ,DUMO   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /BITS/ DUM(8),KDIM
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /DOL/ KA(20),DUMO(180)
      GO TO (10,60,30,40,30,20), KUMP
C
! *********************************************************************
C KUMP = 4: EOF encountered
! --------
C
40    GO TO 20!WRITE(IOU6,42)
C42    FORMAT(/,'Fatal Error.  End-of-file encountered.')
      GO TO 20
C
! *********************************************************************
C KUMP = 3 or 5: Abort job when scratch dimension exceeded in GPS or CM
! -------------
C
   30 CONTINUE
      KDIM=2*KDIM
      KDIM=MIN(KDIM,NN)
      WRITE(IOU6,32) NN,KDIM
   32 FORMAT(/,'Fatal Error.  Scratch array dimension exceeded.',
     -       '  Resubmit job with'/
     -       15X,'$GRID',I8,'  and'/15X,'$DIM ',I8)
      IF(KUMP.EQ.5) GO TO 20
C
C     RECOVER SEQGP CARDS FROM IOU9 IF STOP OCCURRED IN GPS AFTER
C     FINISHING CM.
C
      IF(IPARAM(3).NE.9) GO TO 20
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      REWIND IOU9
C     NCARD=NUMBER OF SEQGP CARDS.
      NCARD=(NN-1)/4 + 1
C     COPY SEQGP CARDS TO UNIT 7.
      DO 37 I=1,NCARD
      READ(IOU9,36) KA
   36 FORMAT(20A4)
C     KA IS SCRATCH SPACE HERE.
37    WRITE(IOU7,36) KA
! B////////////////////////////////////////////////////////////////////B
      WRITE(SEQ,36) KA
! E////////////////////////////////////////////////////////////////////E
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      WRITE(IOU6,38)
   38 FORMAT(/,'SEQGP cards generated by CM have been recovered',
     -       ' and placed on bandit.f07')
      GO TO 20
C
! *********************************************************************
C KUMP = 1: Quit since MAXGRD or MAXDEG too small
! --------
C
   10 WRITE(IOU6,15)
   15 FORMAT('Use $GRID N card and/or increase memory in code.')
      go to 20
c
! *********************************************************************
c KUMP = 6: Quit since $DIM too small for rigid element
! --------
c
60    KDIM=2*KDIM
      WRITE(IOU6,62) KDIM
62    FORMAT('Resubmit job with'/5X,'$DIM',I8)
      go to 20
C
! *********************************************************************
! KUMP = 2: Quit for all other reasons than above
! --------

C   20 WRITE(IOU8,'(A)') 'Fatal Error.  bandit.f08 deleted.'
   20 if(iparam(1).eq.3) close(iou8,status='delete')
C     WRITE(IOU6,'(A)') '0End of BANDIT Job.'
C     STOP 13
      IER=KUMP
      RETURN
      END SUBROUTINE FINISH

! ##################################################################################################################################
      SUBROUTINE FIXIT(LIST,NL)
C
C     THIS ROUTINE COMPRESSES OUT ZEROS AND MULTIPLE ENTRIES IN A LIST
C     ORIGINALLY OF LENGTH NL.  A CORRECTED LENGTH NL IS RETURNED TO
C     THE CALLING PROGRAM.
C
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,I1     ,J      ,NL     ,NL1

! E////////////////////////////////////////////////////////////////////E
C
C     DELETE ZEROS.
C
      CALL ZERO(LIST,NL)
C
C     DELETE DUPLICATE ENTRIES.
C
      IF(NL.LE.1) RETURN
      NL1=NL-1
      DO 20 I=1,NL1
      I1=I+1
      DO 10 J=I1,NL
      IF(LIST(I).NE.LIST(J)) GO TO 10
      LIST(I)=0
      GO TO 20
10    CONTINUE
20    CONTINUE
C
C     DELETE ZEROS AGAIN.
C
      CALL ZERO(LIST,NL)
      RETURN
      END SUBROUTINE FIXIT

! ##################################################################################################################################
      SUBROUTINE FLIP(LIST,N,INV,II3,IER)
C
C     CONVERT $-ARRAY LIST OF LENGTH N FROM ORIGINAL TO INTERNAL NUMBERS
C
      INTEGER II3, INV(2,II3),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IER    ,N

! E////////////////////////////////////////////////////////////////////E
C     CHECK FOR DUPLICATE AND ZERO ENTRIES AND REDUCE N IF NECESSARY.
      CALL FIXIT(LIST,N)
      IF(N.LE.0) RETURN
      DO 10 I=1,N
   10 LIST(I)=INTERN(LIST(I),INV,II3,IER)
      IF(IER.GT.0) RETURN
      RETURN
      END SUBROUTINE FLIP

! ##################################################################################################################################
      SUBROUTINE FNDIAM(SND1,SND2,NDSTK,NR,NDEG,LVL,LVLS1,LVLS2,
     -  IWK,IDFLT,NDLST,IDIM,IER)
C
C  FNDIAM IS THE CONTROL PROCEDURE FOR FINDING THE PSEUDO-DIAMETER OF
C  NDSTK AS WELL AS THE LEVEL STRUCTURE FROM EACH END
C
C  SND1-        ON INPUT THIS IS THE NODE NUMBER OF THE FIRST
C               ATTEMPT AT FINDING A DIAMETER.  ON OUTPUT IT
C               CONTAINS THE ACTUAL NUMBER USED.
C  SND2-        ON OUTPUT CONTAINS OTHER END OF DIAMETER
C  LVLS1-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND1 AS ROOT
C  LVLS2-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND2 AS ROOT
C  IDFLT-       FLAG USED IN PICKING FINAL LEVEL STRUCTURE, SET
C               =1 IF WIDTH OF LVLS1 @ WIDTH OF LVLS2, OTHERWISE =2
C  LVL,IWK-     WORKING STORAGE
C
      INTEGER FLAG,NR,SND,SND1,SND2
      COMMON /GRA/ N,IDPTH,DUMG
      DIMENSION NDSTK(NR,*),NDEG(*),LVL(*),LVLS1(*),LVLS2(*),IWK(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IDFLT  ,IDIM   ,IDPTH  ,IER    ,IWK    ,LVL    ,
     &         LVLBOT ,LVLN   ,LVLS1  ,LVLS2  ,LVLWTH ,MAXLW  ,mtw1   ,
     &         MTW2   ,N      ,NDEG   ,NDSTK  ,NDXL   ,ndxn

      REAL    DUMG

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER NDLST(IDIM)
C     DIMENSION OF NDLST IS THE MAX NUMBER OF NODES IN LAST LEVEL.
      mtw1=0
      ndxn=0
      FLAG=0
      MTW2=N
      SND=SND1
C  ZERO LVL TO INDICATE ALL NODES ARE AVAILABLE TO TREE
   20 DO 25 I=1,N
        LVL(I)=0
   25 CONTINUE
      LVLN=1
C  DROP A TREE FROM SND
      CALL TREE(SND,NDSTK,NR,LVL,IWK,NDEG,LVLWTH,LVLBOT,LVLN,MAXLW,MTW2)
      IF(FLAG.GE.1) GO TO 110
      FLAG=1
   70 IDPTH=LVLN-1
      MTW1=MAXLW
C  COPY LEVEL STRUCTURE INTO LVLS1
      DO 75 I=1,N
        LVLS1(I)=LVL(I)
   75 CONTINUE
      NDXN=1
      NDXL=0
      MTW2=N
C  SORT LAST LEVEL BY DEGREE  AND STORE IN NDLST
      CALL SORTDG(NDLST,IWK(LVLBOT),NDXL,LVLWTH,NDEG)
      IF(NDXL.LE.IDIM) GO TO 100
C     DIMENSION EXCEEDED; STOP JOB.
      CALL FINISH(3,IER)
      RETURN
  100 CONTINUE
      SND=NDLST(1)
      GO TO 20
  110 IF(IDPTH.GE.LVLN-1) GO TO 120
C  START AGAIN WITH NEW STARTING NODE
      SND1=SND
      GO TO 70
  120 IF(MAXLW.GE.MTW2) GO TO 130
      MTW2=MAXLW
      SND2=SND
C  STORE NARROWEST REVERSE LEVEL STRUCTURE IN LVLS2
      DO 125 I=1,N
        LVLS2(I)=LVL(I)
  125 CONTINUE
  130 IF(NDXN.EQ.NDXL) GO TO 140
C  TRY NEXT NODE IN NDLST
      NDXN=NDXN+1
      SND=NDLST(NDXN)
      GO TO 20
  140 IDFLT=1
      IF(MTW2.LE.MTW1) IDFLT=2
      RETURN
      END SUBROUTINE FNDIAM

! ##################################################################################################################################
      SUBROUTINE FRONT(KG,ILD,NN,EL,MEM,IER)
C
C     GENERATE ELEMENT SEQUENCE FOR FRONTAL SOLVERS.
C
C     G.C. EVERSTINE, NSWCCD 204, 12/3/90 (revised 5/22/96)
C
C     KG      = ARRAY FOR STORING THE CONNECTIONS FOR AN ELEMENT
C               (SCRATCH)
C     ILD(I)  = NEW LABEL FOR NODE WITH ORIGINAL INTERNAL LABEL I
C               (INPUT)
C     NN      = NUMBER OF GRID POINTS IN MODEL (INPUT)
C     EL(1,I) = ORIGINAL ELEMENT ID FOR ELEMENT I
C     EL(2,I) = SMALLEST NODE NUMBER IN NEW SEQUENCE IN ELEMENT I
C               (EL IS A SCRATCH ARRAY)
C     MEM     = MEMORY AVAILABLE FOR EL ARRAY (INPUT)
C
      INTEGER KG(*),ILD(*),EL(2,*),EID,EMOD
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IEL    ,IER    ,J      ,J1     ,K      ,KFLAG  ,
     &         KMIN   ,L      ,LOC    ,MEM    ,mflag  ,NBW    ,NCM    ,
     &         NEED   ,NEL    ,NEL1   ,NEQ    ,NEQR   ,NLINK  ,NN     ,
     &         NP     ,NPT    ,NZERO  ,OBW    ,OP

! E////////////////////////////////////////////////////////////////////E
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
      REWIND IOU13
C
C     CHECK MEMORY
C
c     I=MAX(NEL,NN)
      I=MAX(NEL,NN,mem/4)
      NEED=4*I
      IF(MEM.LT.NEED) then
         WRITE(IOU6,5) NEED,MEM
5        FORMAT(/,'Fatal Error.  Insufficient memory available',
     -          ' for element sort for frontal solver.'/
     -          12X,I10,' words needed,',I10,' words available.')
         CALL FINISH(1,IER)
         RETURN
      end if
C
C     INITIALIZE THE SCATTER SORT ARRAY
C
      EMOD=2*I-IFIX(2.3715*SQRT(FLOAT(I)))
      DO I=1,EMOD
         EL(1,I)=0
         EL(2,I)=0
      end do
C
C     READ ELEMENT CONNECTIONS IN ORIGINAL INTERNAL SORT.
C     FIND LOWEST NUMBERED NODE IN NEW SEQUENCE.
C     LOOK FOR DUPLICATE ELEMENT NUMBERS.
C     FILL UP SCATTER SORT ARRAY (EL) SO AS TO PERFORM ROUGH SORT OF
C        ELEMENTS IN ORDER OF LOWEST NUMBERED NODE IN NEW SEQUENCE.
C        (FOR MOST STRUCTURES, THE ELEMENTS WILL BE IN PROPER SORT
C        AFTER THIS STEP.)
C
      DO 50 IEL=1,NEL
         READ(IOU13) EID,NPT,(KG(I),I=1,NPT)
         KMIN=NN
         DO I=1,NPT
            KMIN=MIN(ILD(KG(I)),KMIN)
         end do
         IF(KMIN.EQ.0) GO TO 180
c        Add more space for tet meshes, which have lots of elements
c        for each grid point.  Ideally, the memory available would
c        be at least 40 times the number of grids.
         kmin=10*kmin
         LOC=KMIN-1
30       LOC=MOD(LOC,EMOD)+1
         IF(EL(1,LOC).EQ.0) GO TO 40
         IF(EL(1,LOC).EQ.EID) GO TO 160
         GO TO 30
40       EL(1,LOC)=EID
         EL(2,LOC)=KMIN
50    CONTINUE
C
C     SQUEEZE OUT THE ZEROS IN THE SORT ARRAY.
C
      K=0
      DO 60 I=1,EMOD
         IF(EL(1,I).EQ.0) GO TO 60
         K=K+1
         EL(1,K)=EL(1,I)
         EL(2,K)=EL(2,I)
60    CONTINUE
      IF(K.NE.NEL) then
         WRITE(IOU6,'(A)') ' Logic error in FRONT'
         call finish(6,IER)
         RETURN
      end if
C
C     COMPLETE SORT OF ELEMENTS.  BY NOW, ELEMENTS ARE PROBABLY
C     ALREADY IN SORT.  IF NOT, THIS BUBBLE SORT WILL FINISH IT.
C
      IF(NEL.LE.1) RETURN
      NEL1=NEL-1
      DO 100 I=1,NEL1
         K=NEL-I
         KFLAG=0
         DO 90 J=1,K
            J1=J+1
            IF(EL(2,J).LE.EL(2,J1)) GO TO 90
            KFLAG=1
            L=EL(1,J)
            EL(1,J)=EL(1,J1)
            EL(1,J1)=L
            L=EL(2,J)
            EL(2,J)=EL(2,J1)
            EL(2,J1)=L
90       CONTINUE
         IF(KFLAG.EQ.0) GO TO 110
100   CONTINUE
C
C     WRITE OUT ELEMENT ORDER ON OUTPUT FILE AND UNIT 14.
c        mflag = 0    write to units 6 and 14
c                1    write to unit 14 only
C
110   mflag=1
      if(mflag.eq.0) then
         WRITE(IOU6,120)
120      FORMAT(/,'Element ordering for frontal solution (bandit.f14)')
         DO I=1,NEL,5
            K=MIN(I+4,NEL)
            WRITE(IOU6,140) I,(EL(1,J),J=I,K)
140         FORMAT(I10,' =',5I10)
         end do
      else
         Write(iou6,145)
145      format('See bandit.f14 for element ordering for',
     -          ' frontal solution.')
      end if
      OPEN(IOU14,FILE='bandit.f14',FORM='FORMATTED',STATUS='replace')
      REWIND IOU14
      WRITE(IOU14,150) (EL(1,J),J=1,NEL)
150   FORMAT(5I10)
      close(iou14)
C
      RETURN
C
c     Fatal error messages.
c
160   WRITE(IOU6,170) EID
170   FORMAT(/,'Fatal Error.  Duplicate element ID',i10)
      call finish(6,IER)
      RETURN
c
180   WRITE(IOU6,190) EID
190   FORMAT(/,'Fatal Error.  Element',I9,
     -       ' has a zero grid point connection.')
      call finish(6,IER)
      RETURN
C
      END SUBROUTINE FRONT

! ##################################################################################################################################
      SUBROUTINE GIBSTK(NDSTK,NR,IOLD,RENUM,NDEG,LVL,LVLS1,LVLS2,CCSTOR,
     - IBW2,IPF2,JUMP,ICRIT,NHIGH,NLOW,NACUM,SIZE,STPT,IDIM,IER)
C
C  GIBBSTOCK USES GRAPH THEORETICAL METHODS TO PRODUCE A PERMUTATION
C  OF AN INPUT ARRAY WHICH REDUCES ITS BANDWIDTH.
C  PROGRAMMED BY H.L.CRANE JR.  COLLEGE OF WILLIAM + MARY   3/74.
C
C      MODIFIED BY G. C. EVERSTINE, DTRC.
C
C  THE FOLLOWING INPUT PARAMETERS ARE REQUIRED--NDSTK,NR,N,IDEG,IOLD
C
C  THESE INTEGER ARRAYS MUST BE DIMENSIONED IN THE CALLING PROGRAM--
C  NDSTK(NR,D1),RENUM(D2+1),NDEG(D2),IOLD(D2),LVL(D2),LVLS1(D2),
C  LVLS2(D2),CCSTOR(D2)   WHERE D1 .GE. MAX DEGREE OF ANY NODE AND
C  D2 AND NR ARE .GE. THE TOTAL NUMBER OF NODES IN THE GRAPH.
C
C  EXPLANATION OF PARAMETERS--
C  NDSTK-       ADJACENCY ARRAY REPRESENTING GRAPH TO BE PROCESSED
C               NDSTK(I,J)=NODE NUMBER OF JTH CONNECTION TO NODE
C               NUMBER I.  A CONNECTION OF A NODE TO ITSELF IS NOT
C               LISTED.  EXTRA POSITIONS MUST HAVE ZERO FILL.
C  NR-          ROW DIMENSION ASSIGNED NDSTK IN CALLING PROGRAM
C  IOLD(I)-     RENUMBERING OF ITH NODE BEFORE GIBBSTOCK PROCESSING
C               IF NO RENUMBERING EXISTS THEN ILD(1)=1,ILD(2)=2, ETC.
C  N-           NUMBER OF NODES IN GRAPH BEING PROCESSED
C  IDEG-        MAX DEGREE OF ANY NODE IN GRAPH BEING PROCESSED
C   JUMP IS SET TO 0 IF EITHER CRITERION IS REDUCED.
C     ICRIT=RESEQUENCING CRITERION
C          1=RMS WAVEFRONT
C          2=BANDWIDTH
C          3=PROFILE
C          4=WAVEFRONT (MAX)
C  ON OUTPUT THESE VARIABLES CONTAIN THE FOLLOWING INFORMATION--
C  RENUM(I)-    THE NEW NUMBER FOR THE ITH NODE
C  NDEG(I)-     THE DEGREE OF THE ITH NODE
C  IDPTH-       NUMBER OF LEVELS IN GIBBSTOCK LEVEL STRUCTURE
C  IBW2-        THE BANDWITH AFTER RENUMBERING
C  IPF2-        THE PROFILE AFTER RENUMBERING
C  THE FOLLOWING ONLY HAVE MEANING IF THE GRAPH WAS ALL ONE COMPONENT
C  LVL(I)-      INDEX INTO LVLS1 TO THE FIRST NODE IN LEVEL I
C               LVL(I+1)-LVL(I)= NUMBER OF NODES IN ITH LEVEL
C  LVLS1-       LEVEL STRUCTURE CHOSEN BY GIBBSTOCK
C  LVLS2(I)-    THE LEVEL ASSIGNED TO NODE I BY GIBBSTOCK
C
C    THE FOLLOWING SUBROUTINES WERE WRITTEN BY N. GIBBS, W. POOLE,
C    P. STOCKMEYER, AND H. CRANE OF THE COLLEGE OF WILLIAM AND MARY - -
C     DGREE,FNDIAM,GIBSTK,NUMBER,PIKLVL,RSETUP,SORTDG,SORT2,TREE.
C
! B////////////////////////////////////////////////////////////////////B
! Can't define SORT2 here. It is an integer function that is part of
! this BANDIT_MODULE module. In the version of BANDIT_SUBRS that is not
! a module but a collection of subroutines in a file, then it must be
! defined here (or in an EXTERNAL command)
      INTEGER STNODE,RVNODE,RENUM,XC,STNUM,SBNUM
! E////////////////////////////////////////////////////////////////////E
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IBW1   ,IBW2   ,ICRIT  ,IDEG   ,IDFLT  ,IDIM   ,
     &         IDPTH  ,IER    ,IPARAM ,IPF1   ,IPF2   ,ISDIR  ,IWALL  ,
     &         JUMP   ,
     &         LOWDG  ,LROOT  ,LVL    ,LVLBOT ,LVLN   ,LVLS1  ,LVLS2  ,
     &         LVLWTH ,
     &         MAXB   ,MAXB0  ,MAXLW  ,MAXW0  ,MAXW1  ,MAXWA  ,MAXWB  ,
     &         MM     ,
     &         N      ,NBW    ,NCM    ,NDEG   ,NDSTK  ,NFLG   ,NN     ,
     &         NP     ,NR     ,NUM    ,NZERO

      REAL     AVERW0 ,AVERWB ,BRMS0  ,BRMS1  ,BRMSA  ,BRMSB  ,
     &         CRIT1  ,CRIT2  ,DUMD   ,DUMS   ,RMS0   ,RMS1   ,
     &         RMSA   ,RMSB   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,DUMD(4)
C     OLD AND NEW MAX AND RMS WAVEFRONT FOR ENTIRE PROBLEM,
C          NOT JUST GIBSTK.
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /B/ IPARAM(20)
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /GRA/ N,IDPTH,IDEG
      INTEGER NHIGH(IDIM),NLOW(IDIM),NACUM(IDIM),SIZE(*),STPT(*)
C     SIZE AND STPT HAVE DIMENSION IDIM/2 AND SHOULD BE CONTIGUOUS IN
C                  CORE WITH SIZE FIRST.
C     XC=NUMBER OF SUB-COMPONENTS RESULTING AFTER REMOVING DIAMETER
C        FROM ONE COMPONENT OF ORIGINAL GRAPH.
      DIMENSION NDSTK(NR,*),LVL(*),LVLS1(*),LVLS2(*),RENUM(*),NDEG(*)
      INTEGER CCSTOR(*),IOLD(*),OBW,OP,XCMAX,SUMW0,SUMWB
      REAL IM1,IM2
      CALL TIMER(TA,IWALL,0,IOU6)
      XCMAX=IDIM/2
      NCM=0
      N=NN
      IBW2=0
      IPF2=0
C  SET RENUM(I)=0 FOR ALL I TO INDICATE NODE I IS UNNUMBERED
      DO 10 I=1,N
        RENUM(I)=0
   10 CONTINUE
C   COMPUTE DEGREE OF EACH NODE AND ORIGINAL B AND P.
      CALL DGREE(NDSTK,NR,NDEG,IOLD,IBW1,IPF1)
C   COMPUTE ORIGINAL ACTIVE COLUMN DATA IF NOT ALREADY SUPPLIED BY
C        CUTHILL.
      IF(IPARAM(3)-9) 20,15,20
   15 MAXWA=MAXW1
      RMSA=RMS1
      BRMSA=BRMS1
      GO TO 25
   20 CONTINUE
      CALL WAVEY(NDSTK,NR,IOLD,LVL,0,LVLS2,LVLS1,MAXB0,MAXW0,AVERW0,
     - SUMW0,RMS0,BRMS0)
      MAXWA=MAXW0
      RMSA=RMS0
      BRMSA=BRMS0
      WRITE(IOU6,29)
   29 FORMAT(/,'Before resequencing:')
      WRITE(IOU6,51) MAXB0,SUMW0,MAXW0,AVERW0,RMS0,BRMS0
   51 FORMAT(5X,'Bandwidth',I18/5X,'Profile',I20/
     -       5X,'Max Wavefront',I14/5X,'Avg Wavefront',F14.3/
     -       5X,'RMS Wavefront',F14.3/5X,'RMS Bandwidth',F14.3)
   25 CONTINUE
C  SBNUM= LOW END OF AVAILABLE NUMBERS FOR RENUMBERING
C  STNUM= HIGH END OF AVAILABLE NUMBERS FOR RENUMBERING
      SBNUM=1
      STNUM=N
C  NUMBER THE NODES OF DEGREE ZERO
      DO 40 I=1,N
        IF(NDEG(I).GT.0) GO TO 40
        RENUM(I)=STNUM
        STNUM=STNUM-1
   40 CONTINUE
C   NODES OF ZERO DEGREE APPEAR LAST IN NEW SEQUENCE.
      NZERO=N-STNUM
      NCM=NZERO
C  FIND AN UNNUMBERED NODE OF MIN DEGREE TO START ON
   50 LOWDG=IDEG+1
      NCM=NCM + 1
      NFLG=1
      ISDIR=1
      DO 70 I=1,N
        IF(NDEG(I).GE.LOWDG) GO TO 70
        IF(RENUM(I).GT.0) GO TO 70
        LOWDG=NDEG(I)
        STNODE=I
   70 CONTINUE
C  FIND PSEUDO-DIAMETER AND ASSOCIATED LEVEL STRUCTURES.
C  STNODE AND RVNODE ARE THE ENDS OF THE DIAM AND LVLS1 AND LVLS2
C  ARE THE RESPECTIVE LEVEL STRUCTURES.
      CALL FNDIAM(STNODE,RVNODE,NDSTK,NR,NDEG,LVL,LVLS1,LVLS2,CCSTOR,
     - IDFLT,SIZE,IDIM,IER)
      IF(IER.GT.0) RETURN
      IF(NDEG(STNODE).LE.NDEG(RVNODE)) GO TO 75
C  NFLG INDICATES THE END TO BEGIN NUMBERING ON
      NFLG=-1
      STNODE=RVNODE
   75 CALL RSETUP(LVL,LVLS1,LVLS2,NACUM,IDIM,IER)
      IF(IER.GT.0) RETURN
C  FIND ALL THE CONNECTED COMPONENTS  (XC COUNTS THEM)
      XC=0
      LROOT=1
      LVLN=1
      DO  80 I=1,N
        IF(LVL(I).NE.0) GO TO 80
        XC=XC+1
        IF(XC.LE.XCMAX) GO TO 85
C     DIMENSION EXCEEDED  . . .  STOP JOB.
      CALL FINISH(3,IER)
      RETURN
   85   CONTINUE
        STPT(XC)=LROOT
        CALL TREE(I,NDSTK,NR,LVL,CCSTOR,NDEG,LVLWTH,LVLBOT,LVLN,MAXLW,N)
        SIZE(XC)=LVLBOT+LVLWTH-LROOT
        LROOT=LVLBOT+LVLWTH
        LVLN=LROOT
   80 CONTINUE
      IF(SORT2(XC,SIZE,STPT).EQ.0) GO TO 90
      CALL PIKLVL(LVLS1,LVLS2,CCSTOR,IDFLT,ISDIR,XC,NHIGH,NLOW,
     -  NACUM,SIZE,STPT)
C  ON RETURN FROM PIKLVL, ISDIR INDICATES THE DIRECTION THE LARGEST
C  COMPONENT FELL.  ISDIR IS MODIFIED NOW TO INDICATE THE NUMBERING
C  DIRECTION.  NUM IS SET TO THE PROPER VALUE FOR THIS DIRECTION.
   90 ISDIR=ISDIR*NFLG
      NUM=SBNUM
      IF(ISDIR.LT.0) NUM=STNUM
C
      CALL NUMBR(STNODE,NUM,NDSTK,LVLS2,NDEG,RENUM,LVLS1,LVL,NR,NFLG,
     -  IBW2,IPF2,CCSTOR,ISDIR,NHIGH,NLOW,NACUM,SIZE,IDIM,IER)
      IF(IER.GT.0) RETURN
C
C  UPDATE STNUM OR SBNUM AFTER NUMBERING
      IF(ISDIR.LT.0) STNUM=NUM
      IF(ISDIR.GT.0) SBNUM=NUM
      IF(SBNUM.LE.STNUM) GO TO 50
C
C  COMPUTE THE NEW BANDWIDTH, PROFILE, AND WAVEFRONT.
C
      CALL WAVEY(NDSTK,NR,RENUM,LVL,0,LVLS2,LVLS1,MAXB,MAXWB,AVERWB,
     1   SUMWB,RMSB,BRMSB)
C
      IBW2=MAXB
      IPF2=SUMWB
      WRITE(IOU6,705)
  705 FORMAT(/,'After resequencing by Gibbs-Poole-Stockmeyer (GPS):')
      WRITE(IOU6,51) MAXB,SUMWB,MAXWB,AVERWB,RMSB,BRMSB
C
C     CHECK NEW NUMBERING AGAINST OLD NUMBERING.
C
      GO TO (130,135,140,145), ICRIT
  130 IM1=RMSA
      IM2=IPF1
      CRIT1=RMSB
      CRIT2=IPF2
      GO TO 92
  135 IM1=IBW1
      IM2=IPF1
      CRIT1=IBW2
      CRIT2=IPF2
      GO TO 92
  140 IM1=IPF1
      IM2=IBW1
      CRIT1=IPF2
      CRIT2=IBW2
      GO TO 92
  145 IM1=MAXWA
      IM2=RMSA
      CRIT1=MAXWB
      CRIT2=RMSB
      GO TO 92
   92 CONTINUE
      IF(CRIT1-IM1) 110,94,97
   94 IF(CRIT2.LT.IM2) GO TO 110
C
   97 CONTINUE
C  IF ORIGINAL NUMBERING IS BETTER THAN NEW ONE, SET UP TO RETURN IT
      DO 100 I=1,N
        RENUM(I)=IOLD(I)
  100 CONTINUE
      IBW2=IBW1
      IPF2=IPF1
      MAXWB=MAXWA
      RMSB=RMSA
      BRMSB=BRMSA
      GO TO 112
C
C    EQUATE CORRESPONDING GPS AND BANDIT VARIABLES.
C
  110 CONTINUE
      JUMP=0
  112 CONTINUE
      IF(IPARAM(3).NE.8) GO TO 115
      OBW=IBW1
      OP=IPF1
  115 NBW=IBW2
      NP=IPF2
      MAXW1=MAXWB
      RMS1=RMSB
      BRMS1=BRMSB
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
c     IF(TB.GT.1.E-5) WRITE(IOU6,610) TB
      WRITE(IOU6,610) TB
  610 FORMAT(5X,'CP time',F20.3)
      RETURN
      END SUBROUTINE GIBSTK

! ##################################################################################################################################
      SUBROUTINE GRID(IER,IOU6)
C
C     PARTITION OPEN CORE AND COMPUTE PROBLEM SIZE LIMITS.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IBYTE  ,IER    ,IOU6   ,IPASS  ,KDIM   ,KMOD   ,KOR    ,
     &         MAXDEG ,MAXGRD ,NBITIN ,NBYTE  ,NGRID  ,NW

      REAL     DUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ NBITIN,KOR,DUM,NGRID,IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /A/ MAXGRD,MAXDEG,KMOD
C
      NW=1
CPACK NW=2
C     NW = INTEGER PACKING DENSITY (NUMBER OF INTEGERS PER WORD)
C     NBITIN IS USED ON CDC.  NBYTE AND IBYTE ARE USED ON PACKED CRAY.
      NBITIN=30
      NBYTE=8/NW
      IBYTE=9-NBYTE
      MAXGRD=NGRID+NW-1
      MAXGRD=MAXGRD-MOD(MAXGRD,NW)
      MAXDEG=(KOR-8*MAXGRD-3)/(MAXGRD/NW+1)
      MAXDEG=MIN(MAXDEG,MAXGRD-1)
      IF(MAXDEG.GT.0) RETURN
C     MAXDEG NEGATIVE  - - -   FATAL ERROR.
      CALL COREKO
      WRITE(IOU6,30)
   30 FORMAT(/,'Fatal Error.  Insufficient memory.')
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE GRID

! ##################################################################################################################################
      FUNCTION IDIST(NS,ML,MAXLEV,IG,II1,IC,IDEG,IDIS,IW,ICC)
C
C     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DISTANCE OF ANY NODE
C          IN COMPONENT IC(NS) FROM THE NODE NS.
C     THE DISTANCE OF EACH NODE IN THIS COMPONENT IS STORED IN THE
C          ARRAY IDIS.
C     THE MAXIMUM NUMBER OF NODES AT THE SAME DISTANCE FROM NS IS
C          STORED IN ML.
C
C     INPUT:  IG,IC,IDEG,ICC,NS,MAXLEV
C     OUTPUT: IDIS,IW,ML
C
      INTEGER II1
      DIMENSION IG(II1,*),IC(*),IDEG(*),IDIS(*),IW(*),ICC(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,ICN    ,IDEG   ,IDIS   ,
     &         IDIST  ,IG     ,II     ,IW     ,
     &         K      ,KI     ,KMOD   ,KO     ,
     &         L      ,LL     ,
     &         MAXDEG ,MAXGRD ,MAXLEV ,ML     ,
     &         N      ,NBITIN ,NN     ,NNC    ,NS

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      ICN=IC(NS)
      NNC=ICC(ICN+1)-ICC(ICN)
      DO 50 I=1,NN
      IF(IC(I)-IC(NS)) 50,40,50
40    IDIS(I)=0
   50 CONTINUE
      LL=1
      L=0
      KI=0
      KO=1
      ML=0
      IW(1)=NS
      IDIS(NS)=-1
  130 KI=KI+1
      IF(KI-LL)135,132,135
  132 L=L+1
      LL=KO+1
      K=KO-KI+1
      IF(K-ML) 135,135,133
  133 ML=K
      IF(ML-MAXLEV) 135,135,220
135    II=IW(KI)
      N=IDEG(II)
      IF(N)140,215,140
  140 DO 200 I=1,N
      IA=IG(II,I)
CPACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IDIS(IA))200,150,200
150   IDIS(IA)=L
      KO=KO+1
      IW(KO)=IA
  200 CONTINUE
      IF(KO-NNC)130,205,205
  205 IDIST=L
      IDIS(NS)=0
      K=KO-LL+1
      IF(K-ML) 206,206,207
  207 ML=K
  206 CONTINUE
      RETURN
  215 L=0
      GO TO 205
  220 IDIST=1
      RETURN
      END FUNCTION IDIST

! ##################################################################################################################################
      SUBROUTINE IGNOR(IG,II1,INV,II3,LIST,N,IER)
C
C     SET UP LIST OF POINTS TO IGNORE IN LIST ARRAY OF LENGTH N.
C
C UPON ENTRY, LIST ALREADY CONTAINS SET OF MPC DEPENDENT NODES, IF ANY.
      INTEGER II1, II3
      DIMENSION IG(II1,*),INV(2,II3),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDIM   ,IER    ,IFIR   ,IG     ,IGDEG  ,IGNORE ,
     &         IIG    ,INV    ,ISTA   ,ISTART ,
     &         J      ,KMOD   ,LIST   ,
     &         MAXDEG ,MAXGRD ,MM     ,N      ,NBITIN ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
C
C     ADD IGNORE POINTS TO LIST.
C
      IF(IIG.LE.0) GO TO 20
      CALL FLIP(IGNORE,IIG,INV,II3,IER)
      IF(IER.GT.0) RETURN
      IF(IIG.LE.0) GO TO 20
      DO 10 I=1,IIG
      J=IGNORE(I)
   10 LIST(J)=J
C
C     ADD POINTS OF DEGREE GT IGDEG TO LIST.
C
   20 IF(IGDEG.LE.0) GO TO 40
      IF(IGDEG.GE.MM) GO TO 40
      DO 30 I=1,NN
      IF(IG(I,IGDEG+1).GT.0) LIST(I)=I
CPACK IF(IUNPK(IG,MAXGRD*IGDEG+I,NBITIN).GT.0) LIST(I)=I
   30 CONTINUE
C
C     COMPRESS OUT ZEROES FROM LIST.
C
   40 N=NN
      CALL ZERO(LIST,N)
      RETURN
      END SUBROUTINE IGNOR

! ##################################################################################################################################
      SUBROUTINE INSERT(KA,KT)
C
C     READ AND INSERT CONTENTS OF INSERT FILE (bandit.ins).
C
C     $INSERT CARD MUST APPEAR BEFORE BEGIN BULK TO BE RECOGNIZED BY
C     BANDIT.  FORMAT IS  $INSERT NCARD , WHERE THE PARAMETER NCARD IS
C     THE NUMBER OF CARDS TO BE INSERTED.  IF NCARD IS MISSING OR ZERO,
C     THE ENTIRE FILE 10 WILL BE INSERTED.  ONLY THE FIRST $INSERT CARD
C     IN DECK WILL BE HONORED.
C
C     KT = RUNNING COUNTER ON TOTAL NUMBER OF CASE CONTROL CARDS
C
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NCARD we will use NCARD_array
! instead. Needed so Lahey doesn't complain about shape of NCARD being
! different than IP (array) in subr READIT
      INTEGER NCARD_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(20)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,II     ,IN     ,KT     ,L      ,NCARD

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DATA IN/0/
      IF(IN.EQ.1) RETURN
      IN=1
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(4),1,69,NCARD_array(1),I)
      NCARD = NCARD_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(NCARD.LE.0) NCARD=999999
C
      OPEN(IOU10,FILE='bandit.ins',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND IOU10
      DO 40 II=1,NCARD
      READ(IOU10,20,END=50) KA
20    FORMAT(20A4)
      KT=KT+1
      L=LENCAS(KA,20)
      WRITE(IOU6,30) KT,(KA(I),I=1,L)
30    FORMAT(I8,'= ',20A4)
      WRITE(IOU8,20) (KA(I),I=1,L)
40    CONTINUE
50    RETURN
      END SUBROUTINE INSERT

! ##################################################################################################################################
      FUNCTION INTERN(IGRID,INV,II3,IER)
C
C     THIS FUNCTION HAS AS ITS VALUE THE INTERNAL NODE LABEL ASSIGNED
C     BY BANDIT TO ORIGINAL GRID POINT IGRID.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IER    ,IGRID  ,II3    ,INTERN ,KMOD   ,LOC    ,MAXDEG ,
     &         MAXGRD

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER INV(2,II3)
      intern=0
      IF(IGRID.LE.0) GO TO 20
      LOC=IGRID-1
   10 LOC=MOD(LOC,KMOD)+1
      IF(INV(1,LOC).EQ.0) GO TO 20
      IF(INV(1,LOC).NE.IGRID) GO TO 10
      INTERN=INV(2,LOC)
      RETURN
C
C     ABORT JOB DUE TO REFERENCE TO NON-EXISTENT GRID POINT.
C
   20 WRITE(IOU6,30) IGRID
   30 FORMAT(/,'Fatal Error.  Grid point',I9,' on $ card not found.')
      call finish(6,IER)
      RETURN
      END FUNCTION INTERN

! ##################################################################################################################################
      FUNCTION KOMPNT(IG,II1,IC,IDEG,IW,ICC)
C
C     THIS FUNCTION HAS AS ITS VALUE THE NUMBER OF COMPONENTS STORED
C     IN THE CONNECTION ARRAY IG.
C     ALSO, IC AND ICC ARE SET UP.
C     IC(I)=COMPONENT INDEX FOR NODE I
C     ICC(I)=THE STARTING POSITION TO BE USED FOR LABELS IN COMPONENT I
C     THUS, ICC(I+1)-ICC(I)= THE NUMBER OF NODES IN COMPONENT I
C
      INTEGER II1
      DIMENSION IG(II1,*),IC(*),IDEG(*),IW(*),ICC(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,IDEG   ,IG     ,II     ,
     &         IS     ,IW     ,
     &         KI     ,KMOD   ,KO     ,KOMPNT ,
     &         MAXDEG ,MAXGRD ,MM     ,
     &         N      ,NBITIN ,NC     ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      kompnt=0
      DO I=1,NN
         ICC(I)=0
         IC(I)=0
      end do
      NC=0
      ICC(1)=1
  105 DO 110 I=1,NN
      IF(IC(I)) 110,120,110
  110 KOMPNT=NC
      RETURN
  120 NC=NC+1
      KI=0
      KO=1
      IW(1)=I
      IC(I)=NC
      IF(NC-1)130,125,125
125   IS=ICC(NC)+1
      ICC(NC+1)=IS
  130 KI=KI+1
      II=IW(KI)
      N=IDEG(II)
      IF(N)140,105,140
  140 DO 200 I=1,N
      IA=IG(II,I)
CPACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IC(IA)) 200,150,200
150   IC(IA)=NC
      KO=KO+1
      IW(KO)=IA
      IS=ICC(NC+1)+1
      ICC(NC+1)=IS
  200 CONTINUE
      IF(KO-KI)105,105,130
      END FUNCTION KOMPNT

! ##################################################################################################################################
      SUBROUTINE LEFT(KA,ITYPE)
C
C     LEFT-ADJUST BULK DATA CARD MNEUMONIC AND RETURN IN KA THE NEW
C     CARD IN FORMAT OF A1,A4,A3,64A1,A1,A4,A3.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,ITYPE  ,MA     ,MB     ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER KA(70)
      ITYPE=0
C
C     RETURN IF MNEUMONIC IS ALREADY LEFT-ADJUSTED.
C              (SEE IF COL. 1 IS BLANK)
C
      IF(KA(1).NE.MB(2)) RETURN
C
      BACKSPACE IOU5
      READ(IOU5,10) (KA(I),I=1,70)
   10 FORMAT(70A1)
C
C     RETURN IF CARD IS ENDDATA CARD.
C
      ITYPE=NDATA(KA)
      IF(ITYPE.EQ.1) RETURN
      BACKSPACE IOU5
C
C     COUNT NUMBER OF LEADING BLANKS.
C
      DO I=2,8
         IF(KA(I).NE.MB(2)) GO TO 40
      end do
      READ(IOU5,30) (KA(I),I=1,70)
   30 FORMAT(A1,A4,A3,65A1,A4,A3)
      RETURN
C
C     VARIABLE FORMATS FOR LEFT-ADJUSTING.
C
   40 KA(2)=MB(2)
      KA(3)=MB(2)
      I=I-1
      GO TO (50,60,70,80,90,100,110), I
C     I=NUMBER OF LEADING BLANKS
C
C     1 BLANK
C
   50 READ(IOU5,55) (KA(I),I=1,70)
   55 FORMAT(1X,A1,A4,A2,65A1,A4,A3)
      RETURN
C
C     2 BLANKS
C
   60 READ(IOU5,65) (KA(I),I=1,70)
   65 FORMAT(2X,A1,A4,66A1,A4,A3)
      RETURN
C
C     3 BLANKS
C
   70 READ(IOU5,75) KA(1),KA(2),(KA(I),I=4,70)
   75 FORMAT(3X,A1,A4,65A1,A4,A3)
      RETURN
C
C     4 BLANKS
C
   80 READ(IOU5,85) KA(1),KA(2),(KA(I),I=4,70)
   85 FORMAT(4X,A1,A3,65A1,A4,A3)
      RETURN
C
C     5 BLANKS
C
   90 READ(IOU5,95) KA(1),KA(2),(KA(I),I=4,70)
   95 FORMAT(5X,A1,A2,65A1,A4,A3)
      RETURN
C
C     6 BLANKS
C
  100 READ(IOU5,105) KA(1),KA(2),(KA(I),I=4,70)
  105 FORMAT(6X,67A1,A4,A3)
      RETURN
C
C     7 BLANKS
C
  110 READ(IOU5,115) KA(1),(KA(I),I=4,70)
  115 FORMAT(7X,66A1,A4,A3)
      RETURN
C
      END SUBROUTINE LEFT

! ##################################################################################################################################
      FUNCTION LENCAS(KA,N)
C
C     DETERMINE THE LENGTH OF AN EXECUTIVE/CASE CONTROL LINE.
C     N = DIMENSION OF ARRAY KA
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,LENCAS ,N

! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(N),BLANK
! B////////////////////////////////////////////////////////////////////B
!xx   DATA BLANK/' '/
      DATA BLANK/1H /
! E////////////////////////////////////////////////////////////////////E
      DO 10 I=N,1,-1
      IF(KA(I).NE.BLANK) GO TO 20
10    CONTINUE
      LENCAS=1
      RETURN
20    LENCAS=I
      RETURN
      END FUNCTION LENCAS

! ##################################################################################################################################
      FUNCTION MAXDGR(NC,IC,IDEG)
C
C     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DEGREE OF ANY NODE OF
C     COMPONENT NC IF NC.GT.0
C     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
C
      INTEGER IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,M      ,MAXDGR ,NC     ,NN

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      M=0
      DO 100 I=1,NN
      IF(NC)40,50,40
40    IF(IC(I)-NC) 100,50,100
50    IF(IDEG(I)-M) 100,100,60
60    M=IDEG(I)
  100 CONTINUE
      MAXDGR=M
      RETURN
      END FUNCTION MAXDGR

! ##################################################################################################################################
      FUNCTION MINDEG(NC,IC,IDEG)
C
C     THIS FUNCTION HAS AS ITS VALUE THE MINIMUM DEGREE OF ANY NODE OF
C     COMPONENT NC IF NC.GT.0
C     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
C
      INTEGER IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,M      ,MINDEG ,NC     ,NN

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      M=600000
      DO 100 I=1,NN
      IF(NC)40,50,40
40    IF(IC(I)-NC) 100,50,100
50    IF(M-IDEG(I)) 100,100,60
60    M=IDEG(I)
  100 CONTINUE
      MINDEG=M
      RETURN
      END FUNCTION MINDEG

! ##################################################################################################################################
      SUBROUTINE MORRIS(LIST,NL,IG,II1)
C
C     THIS ROUTINE DELETES ALL REFERENCE IN THE CONNECTION TABLE IG
C     TO THOSE POINTS IN A LIST OF LENGTH NL.
C
C     NEDGE = NUMBER OF UNIQUE EDGES.
C
C     REVISED 12/4/91 (TO AVOID CRAY COMPILER BUG)
C
      INTEGER II1
      DIMENSION IG(II1,*),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IER    ,IG     ,IJ     ,
     &         J      ,K      ,KMOD   ,L      ,LIST   ,
     &         MAXDEG ,MAXGRD ,MM     ,MM1    ,
     &         N      ,NBITIN ,NEDGE  ,NL     ,NN

      REAL     DUM    ,DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,DUMS(3)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
      IF(NL.LE.0) RETURN
      MM1=MM-1
      DO 60 IJ=1,NL
      I=LIST(IJ)
      DO 50 J=1,MM
      L=IG(I,J)
CPACK L=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(L.EQ.0) GO TO 60
      NEDGE=NEDGE-1
      DO 10 K=1,MM
      IF(IG(L,K).EQ.I) GO TO 15
CPACK IF(IUNPK(IG,MAXGRD*(K-1)+L,NBITIN).EQ.I) GO TO 15
10    CONTINUE
      WRITE(IOU6,'(A)') ' Logic error in MORRIS'
      call finish(6,IER)
15    IF(K.GE.MM) GO TO 40
      DO 30 N=K,MM1
      IG(L,N)=IG(L,N+1)
CPACK IS=IUNPK(IG,MAXGRD*N+L,NBITIN)
CPACK CALL PACK(IG,MAXGRD*(N-1)+L,NBITIN,IS)
30    CONTINUE
40    CONTINUE
      IG(L,MM)=0
      IG(I,J)=0
CPACK CALL PACK(IG,MAXGRD*MM1+L,NBITIN,0)
CPACK CALL PACK(IG,MAXGRD*(J-1)+I,NBITIN,0)
   50 CONTINUE
   60 CONTINUE
      RETURN
      END SUBROUTINE MORRIS

! ##################################################################################################################################
      SUBROUTINE MPC(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IER)
C
C     EXTRACT GRID POINTS FROM MPC EQUATION AND STORE IN KG.
C
C     MAXI=MAXIMUM NUMBER OF GRID POINTS ALLOWED PER ELEMENT.
C     MODIFIED 4/16/93 TO ADD MPCAX CARDS (SHORT FIELD ONLY).
C     ASSUMES FIELD 6 OF FIRST LOGICAL MPCAX CARD IS NOT BLANK.
C
      INTEGER MAXI
      INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IDUM   ,IER    ,II3    ,IPARAM ,ITYPE  ,
     &         J      ,K      ,L      ,
     &         MA     ,MB     ,NIP    ,NOUT   ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      NOUT=IOU12
C
C     L=1 OR 2 FOR SHORT OR LONG FIELDS, RESPECTIVELY
      L=ITYPE-1
      IF(ITYPE.EQ.112) L=1
C     I=COUNTER ON THE GRID POINT COMING UP (I.E., THE NEXT ONE).
      I=1
C
   20 J=4+8*L
      IF(ITYPE.EQ.112.AND.I.EQ.1) J=36
      IF(ITYPE.EQ.112.AND.I.GT.1) J=4
      IF(I.GT.MAXI) GO TO 100
      CALL READIT(KA(J),1,8*L,KG(I),NIP)
      I=I+NIP
      GO TO (30,40), L
   30 J=36
      IF(ITYPE.EQ.112.AND.I.EQ.2) GO TO 65
      GO TO 60
   40 ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,50,END=120) KA
   50 FORMAT(A1,A4,A3,65A1,A4,A3)
C
C     LEFT-ADJUST FIELD 1.
C
      CALL LEFT(KA,IDUM)
      IF(KA(1).NE.MB(4)) GO TO 90
      IF(KA(2).NE.ENDP(1)) GO TO 90
      IF(KA(3).NE.ENDP(2)) GO TO 90
      WRITE(IOU8,50) KA
      J=4
   60 IF(I.GT.MAXI) GO TO 100
      CALL READIT(KA(J),1,8*L,KG(I),NIP)
      I=I+NIP
65    ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,50,END=120) KA
C
C     LEFT-ADJUST FIELD 1.
C
      CALL LEFT(KA,IDUM)
      DO K=3,4
         IF(KA(1).EQ.MB(K)) GO TO 80
      end do
      GO TO 90
   80 IF(KA(2).NE.ENDP(1)) GO TO 90
      IF(KA(3).NE.ENDP(2)) GO TO 90
      WRITE(IOU8,50) KA
      L=1
      IF(K.EQ.5) L=2
      GO TO 20
C
C     END OF LOOP STARTING AT STATEMENT 20.
C
   90 BACKSPACE IOU5
      I=I-1
C
C     CONVERT ORIGINAL GRID NUMBERS TO INTERNAL LABELS.
C
      CALL SCAT(KG,I,INV,II3,NORIG,IER)
      IF(IER.GT.0) RETURN
C
C     DELETE DUPLICATE ENTRIES IN LIST.
C
      DO J=2,I
         IF(KG(J).EQ.KG(1)) KG(J)=0
      end do
      CALL FIXIT(KG,I)
C
C     WRITE OUT LIST OF NODES.
C
      WRITE(NOUT) I,(KG(J),J=1,I)
      RETURN
C
C     FATAL ERROR MESSAGE IF MPC EQUATION HAS MORE THAN MAXI TERMS.
C
  100 WRITE(IOU6,110) MAXI
  110 FORMAT(/,'Fatal Error.  MPC equation has more than',I8,' terms.')
      WRITE(IOU6,115)
  115 FORMAT('Use $DIM N card, where 4*N exceeds the maximum number',
     - ' of terms in any one MPC equation.')
      call finish(6,IER)
      RETURN
C
C     END-OF-FILE ENCOUNTERED
C
120   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE MPC

! ##################################################################################################################################
      SUBROUTINE NASNUM(IG,II1,INV,II3,INT,ICC,ILD,NORIG,IP,KOR,
     -                  KORDIM,IER)
C
C     READ BULK DATA, SET UP CONNECTION TABLE, RESEQUENCE NODES,
C     AND GENERATE SEQGP CARDS.
C
      DIMENSION IG(*),INV(*),INT(*),ICC(*),ILD(*),NORIG(*),IP(*)
      INTEGER KORDIM ,KOR(KORDIM)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IADD   ,IBW2   ,IBYTE  ,ICC    ,ICRIT  ,IDIM   ,
     &         IER    ,IFIR   ,IG     ,IGDEG  ,IGNORE ,II1    ,II3    ,
     &         IIG    ,ILD    ,INT    ,INV    ,IP     ,IPARAM ,IPASS  ,
     &         IPF2   ,IPR    ,ISTA   ,ISTART ,IWALL  ,
     &         J      ,JUMP   ,
     &         K      ,K1     ,K2     ,K3     ,K4     ,K5     ,KDIM   ,
     &         KMOD   ,KORE   ,L      ,
     &         MAXDEG ,MAXGRD ,MEM    ,MINDEG ,MM     ,
     &         N      ,NBW    ,NBYTE  ,NCM    ,NEDGE  ,NEL    ,NEQ    ,
     &         NEQR   ,NLINK  ,NN     ,NORIG  ,NP     ,NW     ,NZERO  ,
     &         OBW    ,OP

      REAL     DUM    ,DUMM   ,DUMS   ,DUMY   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ DUM,KORE,DUMM(2),IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /B/ IPARAM(20)
      COMMON /S/ NN,MM,DUMS(3),NEDGE,IADD,MINDEG,DUMY
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
      CALL TIMER(TA,IWALL,0,IOU6)
C
C     COPY IOU5 TO IOU11 AND REDEFINE IOU5.
C
      CALL NEWIN(IG,IER)
      IF(IER.GT.0) RETURN
C
C     ZERO OUT THE WORKING STORAGE.  THE ARRAY NAME USED HERE MUST
C        BE THE FIRST ONE.
C
      DO 10 I=1,KORE
10    ILD(I)=0
C
      NN=0
      MM=0
      NEDGE=0
      IPASS=0
      MINDEG=500000
      KMOD=2*MAXGRD-IFIX(2.3715*SQRT(FLOAT(MAXGRD)))
      REWIND IOU9
C
C     READ BULK DATA DECK AND SET UP CONNECTION TABLE IG.
C
      CALL REED(IG,II1,INV,II3,NORIG,KOR,KORDIM,INT,IER)
      IF(IER.GT.0) RETURN
      REWIND IOU11
C
      IF(NN.GT.0) GO TO 16
C
      WRITE(IOU6,15)
   15 FORMAT(/,'Fatal Error.  No grid points found on',
     -       ' connection cards.')
      call finish(6,IER)
      RETURN
C
   16 IF(MM.GT.0) GO TO 18
      WRITE(IOU6,17)
   17 FORMAT(/,'Fatal Error.  No connections found.')
      call finish(6,IER)
      RETURN
C
   18 I=NEQ+IIG+IGDEG
      IF(I.LE.0) GO TO 19
C
C     MODIFY IG TO ACCOUNT FOR MPC EQUATIONS.
C
      CALL TIGER(IG,II1,ICC,NORIG,KOR,KORDIM,IER)
      IF(IER.GT.0) RETURN
C
C     SET UP LIST OF POINTS TO IGNORE IN ICC ARRAY.
C
      CALL IGNOR(IG,II1,INV,II3,ICC,N,IER)
      IF(IER.GT.0) RETURN
C
C     DELETE POINTS LISTED IN ICC FROM IG.
C
      CALL MORRIS(ICC,N,IG,II1)
C
C     SORT ORIGINAL GRID NUMBERS AND OUTPUT LIST IN INT.   DETERMINE
C     CORRESPONDENCE ILD BETWEEN UNSORTED AND SORTED INTERNAL LABELS.
C
   19 CALL BRIGIT(INV,II3,INT,ILD,IER)
      IF(IER.GT.0) RETURN
C
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
      WRITE(IOU6,195) TB
195   FORMAT('CP time to set up connection table',F12.3,' seconds')
      WRITE(IOU6,196) NN
196   FORMAT('Number of grid points appearing on connection cards',I9/
     -       'Grid cards are not used.')
      WRITE(IOU6,198) MM
  198 FORMAT('Maximum nodal degree before any nodes are ignored',I8)
C
C     WRITE OUT CONNECTION TABLE IF REQUESTED.
C
      IF(IPARAM(2).EQ.4) CALL PUNCON(IG,II1,IP,ILD)
C
      IF(IPARAM(10).EQ.5) GO TO 20
C
C     PRINT TABLES IF MAXIMUM PRINTING REQUESTED.
C
c     CALL TABLE1(INV,II3,INT,IER)
c     The following four lines replace the call to Table1.
      WRITE(IOU6,23)
23    FORMAT(/4(5X,'Grid   Grid BANDIT')/4(7X,'ID  Count  Label'))
      Write(iou6,24) (int(i),i,intern(int(i),inv,ii3,ier),i=1,nn)
24    FORMAT(4(I9,2I7))
c
      IF(IER.GT.0) RETURN
C
      CALL TABLE2(IG,II1,NORIG,IP)
20    continue
C
C     CONVERT USER-SELECTED STARTING NODES, IF ANY, TO INTERNAL LABELS.
C
      IF(ISTA.GT.0) CALL FLIP(ISTART,ISTA,INV,II3,IER)
      IF(IER.GT.0) RETURN
C
C     SAVE ORIGINAL ORDERING (ILD) IF NECESSARY.
C
      CALL RESET(1,ILD,INT,JUMP)
C
      IPR=0
      IF(IPARAM(10).EQ.6) IPR=1
C     INITIALIZE JUMP.
      JUMP=1
C     JUMP=1, NO IMPROVEMENT OF CRITERION SELECTED
C         =0, IMPROVEMENT
C
C     CHOOSE CRITERION FOR RESEQUENCING
C
      ICRIT=2
      IF(IPARAM(6).EQ.11) ICRIT=3
      IF(IPARAM(6).EQ.12) ICRIT=1
      IF(IPARAM(6).EQ.13) ICRIT=4
C
      I=MAXGRD+2
      J=I+MAXGRD
      K=J+MAXGRD
      IF(IPARAM(3).EQ.8) GO TO 25
C
C     RESEQUENCE NODES WITH CUTHILL-MCKEE ALGORITHM.
C
      CALL CUTHIL(80,1,2,ICRIT,IPR,IG,II1,INV,INV(I),INV(J),INV(K),
     -  INT,ICC,ILD,IP,JUMP,KOR,KORDIM,IER)
c
c     Write to a file the component and original grid id for each grid
c     (to allow sorting to get list of grids in each component).
c     On exit from subroutine cuthil, inv contains array 'ic'.
c     Method CM must be included to get this list.
c
      if(iparam(10).eq.6) then
         open(iou15,file='bandit.f15',form='formatted',status='replace')
         Write(iou15,'(a)') 'Component   Grid'
         Write(iou15,'(i5,i11)') (inv(i),norig(i),i=1,nn)
         Write(iou6,22)
22       format(/,'See bandit.f15 for component list.')
         close(iou15)
      end if
      IF(IER.GT.0) RETURN
C
      IF(IPARAM(3).EQ.7) GO TO 28
C
   25 CONTINUE
C
C     RESET ILD, IF NECESSARY, AND COPY TO INT.
C
      CALL RESET(2,ILD,INT,JUMP)
C
C     SAVE SEQGP CARDS ON IOU9 AFTER EXECUTING CM.  THEN, IN CASE
C     GPS ABORTS DUE TO EXCEEDING SCRATCH DIMENSION, THE CM RESULTS
C     CAN BE RECOVERED AND WRITTEN TO UNIT 7.  SEE SUBROUTINE FINISH.
C     IADD IS IGNORED HERE.
C
      REWIND IOU9
      WRITE(IOU9,26) (NORIG(L),ILD(L),L=1,NN)
   26 FORMAT('SEQGP   ',8I8)
      REWIND IOU9
C
C
C     RESEQUENCE NODES WITH GPS ALGORITHM.
C
      KDIM=KORDIM/4
      K1=1
      K2=K1+KDIM
      K3=K2+KDIM
      K4=K3+KDIM
      K5=K4+KDIM/2
      CALL GIBSTK(IG,II1,INT,ILD,INV(I),INV,INV(J),INV(K),ICC,IBW2,IPF2,
     -  JUMP,ICRIT,KOR(K1),KOR(K2),KOR(K3),KOR(K4),KOR(K5),KDIM,IER)
      IF(IER.GT.0) RETURN
C
   28 CONTINUE
C
C     GENERATE SEQGP CARDS.
C
      CALL SEQGP(NORIG,ILD,INT,JUMP)
C
C     IF REQUESTED, GENERATE SET OF SCALAR SPRINGS (CELAS3) WITH SAME
C     CONNECTIVITY AS ORIGINAL STRUCTURE ON bandit.f09
C
      IF(IPARAM(7).EQ.4) CALL SPRING(IP)
C
C     GENERATE ELEMENT ORDERING FOR FRONTAL SOLVERS.
C
      IF(IPARAM(11).EQ.4) THEN
C        COMPUTE MEMORY AVAILABLE FOR ELEMENT SORT (INV ARRAY).
         MEM=KORE-(MAXGRD+1)
         CALL FRONT(KOR,ILD,NN,INV,MEM,IER)
         IF(IER.GT.0) RETURN
      END IF
C
      RETURN
      END SUBROUTINE NASNUM

! ##################################################################################################################################
      FUNCTION NBULK(KA)
C
C     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A CARD READ BY 80A1 IS
C     THE BEGIN BULK CARD. IF NOT, 0 IS RETURNED.
C
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,K      ,MA     ,MB     ,NBULK  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      NBULK=0
C
C     LOOK FOR FIRST NON-BLANK.
C
      DO I=1,64
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
C
   20 IF(KA(I  ).NE.MA( 2)) RETURN
      IF(KA(I+1).NE.MA( 5)) RETURN
      IF(KA(I+2).NE.MA( 7)) RETURN
      IF(KA(I+3).NE.MA( 9)) RETURN
      IF(KA(I+4).NE.MA(14)) RETURN
      K=I+5
C
C     LOOK FOR FIRST NON-BLANK AFTER -BEGIN-.
C
      DO I=K,69
         IF(KA(I).NE.MB(2)) GO TO 40
      end do
      RETURN
C
   40 IF(KA(I  ).NE.MA( 2)) RETURN
      IF(KA(I+1).NE.MA(21)) RETURN
      IF(KA(I+2).NE.MA(12)) RETURN
      IF(KA(I+3).NE.MA(11)) RETURN
      NBULK=1
      RETURN
      END FUNCTION NBULK

! ##################################################################################################################################
      FUNCTION NDATA(KA)
C
C     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A CARD READ BY 70A1 IS
C     THE ENDDATA CARD.  OTHERWISE, 0 IS RETURNED.
C
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,MA     ,MB     ,NDATA  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      NDATA=0
C
C     LOOK FOR FIRST NON-BLANK.
C
      DO I=1,64
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
C
C     LOOK FOR -ENDDATA-.
C
   20 IF(KA(I  ).NE.MA( 5)) RETURN
      IF(KA(I+1).NE.MA(14)) RETURN
      IF(KA(I+2).NE.MA( 4)) RETURN
      IF(KA(I+3).NE.MA( 4)) RETURN
      IF(KA(I+4).NE.MA( 1)) RETURN
      IF(KA(I+5).NE.MA(20)) RETURN
      IF(KA(I+6).NE.MA( 1)) RETURN
      NDATA=1
      RETURN
      END FUNCTION NDATA

! ##################################################################################################################################
      SUBROUTINE NEWIN(KA,IER)
C
C     COPY BULK DATA DECK FROM UNIT 5 TO UNIT 11.
C     REDEFINE IOU5 TO BE IOU11.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IER    ,L

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER KA(80)
C     KA IS SCRATCH SPACE.
      L=IOU11
      REWIND L
    5 READ(IOU5,10,END=30) KA
   10 FORMAT(80A1)
      WRITE(L,10) KA
      IF(NDATA(KA).EQ.0) GO TO 5
      REWIND L
C     REDEFINE IOU5
      IOU5=L
      RETURN
C
30    CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE NEWIN

! ##################################################################################################################################
      SUBROUTINE NOSEQ(KA)
C
C     WRITE BULK DATA DECK IF RESEQUENCING NOT REQUESTED.
C
      INTEGER KA(80)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IPARAM

! E////////////////////////////////////////////////////////////////////E
C     KA IS SCRATCH SPACE.
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      WRITE(IOU6,10)
   10 FORMAT(/,'Grid point resequencing not requested.' )
   20 READ(IOU5,30,END=40) KA
   30 FORMAT(80A1)
      WRITE(IOU8,30) KA
      IF(NDATA(KA).EQ.0) GO TO 20
40    RETURN
      END SUBROUTINE NOSEQ

! ##################################################################################################################################
      FUNCTION NTHRU(KA,N)
C
C     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A FIELD READ N A1
C     CONTAINS THE CHARACTER STRING "THRU".  IF NOT, 0 IS RETURNED.
C
C     N=NUMBER OF CARD COLUMNS TO SEARCH STARTING AT KA(1)
C
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,MA     ,MB     ,N      ,NA     ,NTHRU  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      NTHRU=0
C
C     LOOK FOR FIRST NON-BLANK.
C
      NA=N-3
      IF(NA.LT.1) RETURN
      DO I=1,NA
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
C
C     LOOK FOR  -THRU-.
C
   20 IF(KA(I  ).NE.MA(20)) RETURN
      IF(KA(I+1).NE.MA( 8)) RETURN
      IF(KA(I+2).NE.MA(18)) RETURN
      IF(KA(I+3).NE.MA(21)) RETURN
      NTHRU=1
      RETURN
      END FUNCTION NTHRU

! ##################################################################################################################################
      SUBROUTINE NUMBR(SND,NUM,NDSTK,LVLS2,NDEG,RENUM,LVLST,LSTPT,NR,
     -                 NFLG,IBW2,IPF2,IPFA,ISDIR,STKA,STKB,STKC,STKD,
     -                 IDIM,IER)
C
C  NUMBR PRODUCES THE NUMBERING OF THE GRAPH FOR MIN BANDWIDTH
C
C  SND-         ON INPUT THE NODE TO BEGIN NUMBERING ON
C  NUM-         ON INPUT AND OUTPUT, THE NEXT AVAILABLE NUMBER
C  LVLS2-       THE LEVEL STRUCTURE TO BE USED IN NUMBERING
C  RENUM-       THE ARRAY USED TO STORE THE NEW NUMBERING
C  LVLST-       ON OUTPUT CONTAINS LEVEL STRUCTURE
C  LSTPT(I)-    ON OUTPUT, INDEX INTO LVLST TO FIRST NODE IN ITH LVL
C               LSTPT(I+1) - LSTPT(I) = NUMBER OF NODES IN ITH LVL
C  NFLG-        =+1 IF SND IS FORWARD END OF PSEUDO-DIAM
C               =-1 IF SND IS REVERSE END OF PSEUDO-DIAM
C  IBW2-        BANDWIDTH OF NEW NUMBERING COMPUTED BY NUMBER
C  IPF2-        PROFILE OF NEW NUMBERING COMPUTED BY NUMBER
C      IBW2 AND IPF2 HERE DO NOT INCLUDE DIAGONAL TERMS.
C  IPFA-        WORKING STORAGE USED TO COMPUTE PROFILE AND BANDWIDTH
C  ISDIR-       INDICATES STEP DIRECTION USED IN NUMBERING(+1 OR -1)
C
      INTEGER IDIM, SND,XA,XB,XC,XD,CX,END,RENUM,TEST
      COMMON /GRA/ N,IDPTH,IDEG
      INTEGER STKA(IDIM),STKB(IDIM),STKC(IDIM),STKD(IDIM)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER NR, IPFA(*)
      DIMENSION NDSTK(NR,*),LVLS2(*),NDEG(*),RENUM(*),LVLST(*),LSTPT(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IBW2   ,IDEG   ,IDPTH  ,IER    ,INX    ,
     &         IPF2   ,IPRO   ,ISDIR  ,
     &         J      ,
     &         KMOD   ,
     &         LND    ,LST    ,LSTPT  ,LVLN   ,LVLS2  ,LVLST  ,
     &         MAXDEG ,MAXGRD ,MAXI   ,
     &         N      ,NBITIN ,NBW    ,NDEG   ,NDSTK  ,NFLG   ,
     &         NSTPT  ,NUM

      REAL     DUMBB

! E////////////////////////////////////////////////////////////////////E
C  SET UP LVLST AND LSTPT FROM LVLS2
      DO 3 I=1,N
        IPFA(I)=0
    3 CONTINUE
      NSTPT=1
      DO 5 I=1,IDPTH
        LSTPT(I)=NSTPT
        DO 5 J=1,N
          IF(LVLS2(J).NE.I) GO TO 5
          LVLST(NSTPT)=J
          NSTPT=NSTPT+1
    5 CONTINUE
      LSTPT(IDPTH+1)=NSTPT
C  THIS ROUTINE USES FOUR STACKS, A,B,C,AND D, WITH POINTERS
C  XA,XB,XC, AND XD.  CX IS A SPECIAL POINTER INTO STKC WHICH
C  INDICATES THE PARTICULAR NODE BEING PROCESSED.
C  LVLN KEEPS TRACK OF THE LEVEL WE ARE WORKING AT.
C  INITIALLY STKC CONTAINS ONLY THE INITIAL NODE, SND.
      LVLN=0
      IF(NFLG.LT.0) LVLN=IDPTH+1
      XC=1
      STKC(XC)=SND
   10 CX=1
      XD=0
      LVLN=LVLN+NFLG
      LST=LSTPT(LVLN)
      LND=LSTPT(LVLN+1)-1
C  BEGIN PROCESSING NODE STKC(CX)
   20 IPRO=STKC(CX)
      RENUM(IPRO)=NUM
      NUM=NUM+ISDIR
      END=NDEG(IPRO)
      XA=0
      XB=0
C  CHECK ALL ADJACENT NODES
      DO 50 I=1,END
        TEST=NDSTK(IPRO,I)
CPACK   TEST=IUNPK(NDSTK,MAXGRD*(I-1)+IPRO,NBITIN)
26      INX=RENUM(TEST)
C  ONLY NODES NOT NUMBERED OR ALREADY ON A STACK ARE ADDED
        IF(INX.EQ.0) GO TO 30
        IF(INX.LT.0) GO TO 50
C  DO PRELIMINARY BANDWIDTH AND PROFILE CALCULATIONS
        NBW=(RENUM(IPRO)-INX)*ISDIR
        IF(ISDIR.GT.0) INX=RENUM(IPRO)
        IF(IPFA(INX).LT.NBW) IPFA(INX)=NBW
        GO TO 50
   30   RENUM(TEST)=-1
C  PUT NODES ON SAME LEVEL ON STKA, ALL OTHERS ON STKB
        IF(LVLS2(TEST).EQ.LVLS2(IPRO)) GO TO 40
        XB=XB+1
      IF(XB.GT.IDIM) GO TO 150
        STKB(XB)=TEST
        GO TO 50
   40   XA=XA+1
      IF(XA.GT.IDIM) GO TO 150
        STKA(XA)=TEST
   50 CONTINUE
C  SORT STKA AND STKB INTO INCREASING DEGREE AND ADD STKA TO STKC
C  AND STKB TO STKD
      IF(XA.EQ.0) GO TO 55
      IF(XA.EQ.1) GO TO 52
      CALL SORTDG(STKC,STKA,XC,XA,NDEG)
      GO TO 55
   52 XC=XC+1
      IF(XC.GT.IDIM) GO TO 150
      STKC(XC)=STKA(XA)
   55 IF(XB.EQ.0) GO TO 65
      IF(XB.EQ.1) GO TO 62
      CALL SORTDG(STKD,STKB,XD,XB,NDEG)
      GO TO 65
   62 XD=XD+1
      IF(XD.GT.IDIM) GO TO 150
      STKD(XD)=STKB(XB)
C  BE SURE TO PROCESS ALL NODES IN STKC
   65 CX=CX+1
      IF(XC.GE.CX) GO TO 20
C  WHEN STKC IS EXHAUSTED LOOK FOR MIN DEGREE NODE IN SAME LEVEL
C  WHICH HAS NOT BEEN PROCESSED
      MAXI=IDEG+1
      SND=N+1
      DO 70 I=LST,LND
        TEST=LVLST(I)
        IF(RENUM(TEST).NE.0) GO TO 70
        IF(NDEG(TEST).GE.MAXI) GO TO 70
        RENUM(SND)=0
        RENUM(TEST)=-1
        MAXI=NDEG(TEST)
        SND=TEST
   70 CONTINUE
      IF(SND.EQ.N+1) GO TO 75
      XC=XC+1
      IF(XC.GT.IDIM) GO TO 150
      STKC(XC)=SND
      GO TO 20
C  IF STKD IS EMPTY WE ARE DONE, OTHERWISE COPY STKD ONTO STKC
C  AND BEGIN PROCESSING NEW STKC
   75 IF(XD.EQ.0) GO TO 100
      DO 80 I=1,XD
        STKC(I)=STKD(I)
   80 CONTINUE
      XC=XD
      GO TO 10
C  DO FINAL BANDWIDTH AND PROFILE CALCULATIONS
  100 DO 120 I=1,N
        IF(IPFA(I).GT.IBW2) IBW2=IPFA(I)
        IPF2=IPF2+IPFA(I)
  120 CONTINUE
      RETURN
C     DIMENSION EXCEEDED  . . .  STOP JOB.
  150 CALL FINISH(3,IER)
      RETURN
      END SUBROUTINE NUMBR

! ##################################################################################################################################
      SUBROUTINE PIKLVL(LVLS1,LVLS2,CCSTOR,IDFLT,ISDIR,XC,NHIGH,NLOW,
     -                  NACUM,SIZE,STPT)
C
C  PIKLVL CHOOSES THE LEVEL STRUCTURE  USED IN NUMBERING GRAPH
C
C  LVLS1-       ON INPUT CONTAINS FORWARD LEVELING INFO
C  LVLS2-       ON INPUT CONTAINS REVERSE LEVELING INFO
C               ON OUTPUT THE FINAL LEVEL STRUCTURE CHOSEN
C  CCSTOR-      ON INPUT CONTAINS CONNECTED COMPONENT INFO
C  IDFLT-       ON INPUT =1 IF WDTH LVLS1@WDTH LVLS2, =2 OTHERWISE
C  NHIGH        KEEPS TRACK OF LEVEL WIDTHS FOR HIGH NUMBERING
C  NLOW-        KEEPS TRACK OF LEVEL WIDTHS FOR LOW NUMBERING
C  NACUM-       KEEPS TRACK OF LEVEL WIDTHS FOR CHOSEN LEVEL STRUCTURE
C  XC-          NUMBER OF CONNECTED COMPONENTS
C  SIZE(I)-     SIZE OF ITH CONNECTED COMPONENT
C  STPT(I)-     INDEX INTO CCSTORE OF 1ST NODE IN ITH CON COMPT
C  ISDIR-       FLAG WHICH INDICATES WHICH WAY THE LARGEST CONNECTED
C               COMPONENT FELL.  =+1 IF LOW AND -1 IF HIGH
C
      INTEGER XC,END
      COMMON /GRA/ N,IDPTH,DUMG
C     DIMENSION OF NHIGH IS MAXIMUM ALLOWABLE NUMBER OF LEVELS.
C     DIMENSION OF SIZE IS MAXIMUM ALLOWABLE NUMBER OF COMPONENTS.
      INTEGER NHIGH(*),NLOW(*),NACUM(*),SIZE(*),STPT(*)
      INTEGER LVLS1(*),LVLS2(*),CCSTOR(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDFLT  ,IDPTH  ,INODE  ,ISDIR  ,IT     ,
     &         J      ,K      ,LVLNH  ,LVLNL  ,
     &         MAX1   ,MAX2   ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
C  FOR EACH CONNECTED COMPONENT DO
      DO 270 I=1,XC
        J=STPT(I)
        END=SIZE(I)+J-1
C  SET NHIGH AND NLOW EQUAL TO NACUM
        DO 205 K=1,IDPTH
          NHIGH(K)=NACUM(K)
          NLOW(K)=NACUM(K)
  205   CONTINUE
C  UPDATE NHIGH AND NLOW FOR EACH NODE IN CONNECTED COMPONENT
        DO 210 K=J,END
          INODE=CCSTOR(K)
          LVLNH=LVLS1(INODE)
          NHIGH(LVLNH)=NHIGH(LVLNH)+1
          LVLNL=LVLS2(INODE)
          NLOW(LVLNL)=NLOW(LVLNL)+1
  210   CONTINUE
        MAX1=0
        MAX2=0
C  SET MAX1=LARGEST NEW NUMBER IN NHIGH
C  SET MAX2=LARGEST NEW NUMBER IN NLOW
        DO 240 K=1,IDPTH
          IF(2*NACUM(K).EQ.NLOW(K)+NHIGH(K)) GO TO 240
          IF(NHIGH(K).GT.MAX1) MAX1=NHIGH(K)
          IF(NLOW(K).GT.MAX2) MAX2=NLOW(K)
  240   CONTINUE
C  SET IT= NUMBER OF LEVEL STRUCTURE TO BE USED
        IT=1
        IF(MAX1.GT.MAX2) IT=2
        IF(MAX1.EQ.MAX2) IT=IDFLT
        IF(IT.EQ.2) GO TO 265
        IF(I.EQ.1) ISDIR=-1
C  COPY LVLS1 INTO LVLS2 FOR EACH NODE IN CONNECTED COMPONENT
        DO 260 K=J,END
          INODE=CCSTOR(K)
          LVLS2(INODE)=LVLS1(INODE)
  260   CONTINUE
C  UPDATE NACUM TO BE THE SAME AS NHIGH
        DO 262 K=1,IDPTH
          NACUM(K)=NHIGH(K)
  262   CONTINUE
        GO TO 270
C  UPDATE NACUM TO BE THE SAME AS NLOW
  265   DO 267 K=1,IDPTH
          NACUM(K)=NLOW(K)
  267   CONTINUE
  270 CONTINUE
      RETURN
      END SUBROUTINE PIKLVL

! ##################################################################################################################################
      SUBROUTINE PUNCON(IG,II1,IP,ILD)
C
C     WRITE CONNECTION TABLE IG ON PERIPHERAL FILE (IOU16).
C
C     ILD(I) = SORTED INTERNAL LABEL FOR NODE WITH UNSORTED INTERNAL
C              LABEL I
C     NN     = NUMBER OF NODES.
C     M      = MAX NODAL DEGREE.
C     IP     = TEMPORARY STORAGE.
C
C     HEADER CARD GIVES NUMBER OF NODES (NN) AND MAX NODAL DEGREE (M)
C     IN 2I5, THEN ONE CARD PER NODE, 24I5.
C     FIELD 1 IS NODE, OTHER FIELDS ARE CONNECTIONS.
C
      INTEGER II1
      DIMENSION IG(II1,*),IP(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IG     ,ILD    ,IP     ,
     &         J      ,K      ,KMOD   ,
     &         M      ,MAXDEG ,MAXGRD ,
     &         NBITIN ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,M,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
      WRITE(IOU16,30) NN,M
      DO 20 I=1,NN
      DO J=1,M
         IP(J)=0
      end do
      DO J=1,M
         K=IG(I,J)
CPACK    K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
         IF(K.EQ.0) GO TO 15
         IP(J)=ILD(K)
      end do
   15 CONTINUE
      WRITE(IOU16,30) ILD(I),(IP(J),J=1,M)
   30 FORMAT(24I5)
C
C     24I5 IS CHOSEN SINCE 24*5=120 IS AN INTEGER MULTIPLE OF 10 (CDC),
C     4 (IBM), AND 6 (UNIVAC), MAKING IT CONVENIENT FOR I/O BETWEEN
C     DISSIMILAR COMPUTERS.
C     120 CHARACTER LINES CAN ALSO BE LISTED ON A PRINTER IF DESIRED.
C     IF THIS FORMAT IS CHANGED, SIMILAR CHANGES ARE REQUIRED IN
C     SUBROUTINES SEQGP AND SPRING.
C
   20 CONTINUE
      WRITE(IOU6,40) IOU16
   40 FORMAT(/,'Connection table generated on bandit.f',I2)
      RETURN
      END SUBROUTINE PUNCON

! ##################################################################################################################################
      SUBROUTINE READIT(KA,N,MAXI,IP,NIP)
C
C     INTERPRET NUMERIC DATA READ IN 80A1 FORMAT.  ONLY NON-NEGATIVE
C     INTEGERS ARE FOUND.
C
C     KA(I) = ITH CHARACTER (A1 FORMAT) (INPUT)
C     N     = NUMBER OF INTEGERS SOUGHT (INPUT)
C     MAXI  = MAXIMUM VALUE OF I FOR SEARCH (INPUT)
C     IP(J) = JTH INTEGER FOUND (OUTPUT)
C     NIP   = NUMBER OF INTEGERS FOUND (OUTPUT)
C
C     TO CONVERT READIT TO A STAND-ALONE UTILITY NOT DEPENDENT ON
C     COMMON BLOCK /ALPHA/, WHICH SUPPLIES THE ARRAY NUM, REPLACE
C     THE COMMON STATEMENT BELOW WITH THE DATA STATEMENT FOR NUM
C     WHICH IS COMMENTED OUT.
C
C     G.C. EVERSTINE, DTRC 128, SEPT. 1973 (REVISED JULY 1978)
C
      INTEGER KA(*),IP(*),FLAG,BLANK
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,KOL    ,
     &         MA     ,MAXI   ,MB     ,N      ,NIP    ,NUM

! E////////////////////////////////////////////////////////////////////E
C     DATA NUM/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
! B////////////////////////////////////////////////////////////////////B
!xx   DATA BLANK/' '/
      DATA BLANK/1H /
! E////////////////////////////////////////////////////////////////////E
      NIP=0
      KOL=0
      DO 40 I=1,N
      IP(I)=0
      FLAG=0
   10 KOL=KOL+1
      IF(KOL.LE.MAXI) GO TO 15
      IF(FLAG.EQ.0) RETURN
      NIP=NIP+1
      RETURN
   15 CONTINUE
      IF(KA(KOL).EQ.BLANK) GO TO 25
      DO 20 J=1,10
      IF(KA(KOL).EQ.NUM(J)) GO TO 30
   20 CONTINUE
   25 IF(FLAG) 10,10,40
   30 FLAG=1
      IP(I)=10*IP(I)+J-1
      GO TO 10
   40 NIP=NIP+1
      RETURN
      END SUBROUTINE READIT

! ##################################################################################################################################
      SUBROUTINE REED(IG,II1,INV,II3,NORIG,KG,MAXI,IDUM,IER)
C
C     READ BULK DATA DECK AND SET UP CONNECTION TABLE IG.
C
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with EID we will use EID_array
! instead. Needed so Lahey doesn't complain about shape of EID being
! different than IP (array) in subr READIT
      INTEGER EID_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER II1, II3, MAXI
      DIMENSION IG(II1,*),INV(2,II3),NORIG(*)
      INTEGER vype,TYPE,WYPE,KA(70),KG(MAXI),IDUM(*),ENDP(2),EID
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,I1     ,I2     ,IER    ,IFLD   ,IG     ,II     ,
     &         INV    ,IPARAM ,ISAVE  ,ITYPE  ,
     &         J      ,K      ,L      ,LEN    ,LOOP   ,
     &         MA     ,MB     ,MDIM   ,ME     ,
     &         NCARD  ,NCON   ,NEL    ,NELEM  ,NEQ    ,NEQR   ,NLINK  ,
     &         NORIG  ,npt    ,NTYPE  ,NUM

      REAL     DUMD

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),
     -              NELEM(160),MDIM
      COMMON /B/ IPARAM(20)
      COMMON /D/ DUMD(6),NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C     MAXI=MAXIMUM NUMBER OF GRID POINTS ALLOWED PER ELEMENT
c          (=kdim4 from main) (input)
      LOGICAL LESSOK,rigide
      IF(IPARAM(11).EQ.4) THEN
         OPEN(IOU13,FORM='UNFORMATTED',STATUS='SCRATCH')
         REWIND IOU13
      END IF
      NEL=0
      NEQ=0
      NEQR=0
      NLINK=0
      NCARD=0
C     NEL   = NUMBER OF ELEMENTS NOT COUNTING RIGID ELEMENTS
C     NEQ   = OVERALL NUMBER OF MPC EQUATIONS, INCLUDING NEQR
C     NEQR  = NUMBER OF MPC EQUATIONS ARISING FROM RIGID LINKS.
C     NLINK = NUMBER OF RIGID LINKS (ELEMENTS)
C     NCARD = NUMBER OF BULK DATA CARDS
C
C     READ BULK DATA CARD.
C
   20 READ(IOU5,30,END=170) KA
   30 FORMAT(A1,A4,A3,65A1,A4,A3)
      NCARD=NCARD+1
c     IF(MOD(NCARD,100).EQ.0) WRITE(IOU17,35) NCARD
35    FORMAT(I9)
C
C     LEFT-ADJUST MNEUMONIC.
C
      CALL LEFT(KA,ITYPE)
C
C     RETURN IF ENDDATA CARD  (A LEFT-ADJUSTED ENDDATA WILL NOT BE
C        FOUND UNTIL LATER).
C
   55 IF(ITYPE.EQ.1) THEN
         WRITE(IOU6,60) NCARD
60       FORMAT(/,'Number of Bulk Data records read',I14)
         RETURN
      END IF
C
C     DETERMINE CARD TYPE.
C
      CALL ELTYPE(KA,ITYPE,NCON,IFLD,LOOP,MAXI,LEN,LESSOK,IER)
      IF(IER.GT.0) RETURN
c
c     Set rigid element flag.
c
      rigide=.false.
      if(itype.ge.91.and.itype.le.93) rigide=.true.
      if(itype.ge.128.and.itype.le.137) rigide=.true.
C
C     IF NOT ENDDATA CARD, WRITE CARD AND CONTINUE.
C
      IF(ITYPE.EQ.1) THEN
         WRITE(IOU6,60) NCARD
         RETURN
      END IF
C
      WRITE(IOU8,30) KA
C
C     IF NO MATCH FOUND, CHECK FOR AXIC, THEN GO BACK AND READ
C             ANOTHER CARD.
C
      IF(ITYPE.EQ.0) THEN
         CALL TAXI(KA,IOU6)
         GO TO 20
      END IF
C
C     CHECK IF MPC CARDS ARE TO BE PROCESSED (MPC, MPC*, MPCAX)
C
      if(itype.le.3.or.itype.eq.112) then
         IF(IPARAM(4).EQ.3) GO TO 20
         NEQ=NEQ+1
C        NEQ=NUMBER OF MPC EQUATIONS
         CALL MPC(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IER)
         IF(IER.GT.0) RETURN
         GO TO 20
      end if
C
C     CHECK FOR AND INTERPRET selected RIGID ELEMENTs
C
      IF(ITYPE.EQ.91.OR.ITYPE.EQ.92.or.
     -      (itype.ge.130.and.itype.le.133)) then
         CALL RIGID(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IDUM,IER,npt)
84       format('frigid',8i6)
         IF(IER.GT.0) RETURN
c        Since RBE1 is read like a rigid element, but processed like a
c        regular element, keep going.
         if(itype.eq.130.or.itype.eq.131) go to 87
         NLINK=NLINK+1
         NELEM(ITYPE)=NELEM(ITYPE)+1
         GO TO 20
      END IF
87    continue
C
C     BLANK OUT FIELD 5 FOR SCALAR ELEMENTS WITH CONNECTIONS NOT LISTED
c         CONSECUTIVELY.
C
      IF(ITYPE.LE.9.or.(ITYPE.GE.73.AND.ITYPE.LE.75)) then
C        LEN=1 FOR SHORT FIELD CARD, 2 FOR LONG FIELD CARD
         I1=4+24*LEN
         I2=I1-1+8*LEN
         DO I=I1,I2
            KA(I)=MB(2)
         end do
      end if
C
C     PROCESS CONNECTION CARD.
C
C     LOOP=NUMBER OF ELEMENTS DEFINABLE ON ONE CARD (1 OR 2, USUALLY 1)
C
      DO 140 II=1,LOOP
C
C     GET ELEMENT ID (EID) (needed for frontal ordering)
C     I1=SUBSCRIPT IN KA ARRAY CORRESPONDING TO BEGINNING OF EID FIELD
      I1=4
      IF(II.EQ.2.AND.LEN.EQ.1) I1=36
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(I1),1,8*LEN,EID_array(1),J)
      EID = EID_array(1)
! E////////////////////////////////////////////////////////////////////E
C
C     DETERMINE GRID CONNECTIONS AND STORE IN KG.
c     THE ACTUAL NUMBER OF NON-ZERO NODES IS RETURNED IN NPT.
c     Grid connections have already been found for RBE1 by Subroutine RIGID.
c
      if(itype.eq.130.or.itype.eq.131) go to 103
        CALL SETUP(KA,NCON,IFLD,KG,NPT,LEN,LESSOK,IER)
        IF(IER.GT.0) RETURN
        IF(II.EQ.2.AND.NPT.EQ.0) GO TO 140
        IF(NPT.EQ.NCON.OR.LESSOK) GO TO 103
        WRITE(IOU6,102) KA
102     FORMAT('Warning.  Missing connection(s) on'/A1,A4,A3,65A1,A4,A3)
103   IF(NPT.LE.0) GO TO 130
C
      if(rigide) then
         nlink=nlink+1
      else
         nel=nel+1
      end if
      NELEM(ITYPE)=NELEM(ITYPE)+1
C
C     CONVERT KG FROM ORIGINAL TO INTERNAL NUMBERS.
C
      CALL SCAT(KG,npt,INV,II3,NORIG,IER)
      IF(IER.GT.0) RETURN
C
C     WRITE EID AND INTERNAL CONNECTIONS ON UNIT 13.
C     SKIP FOR "$FRONTAL NO" and for rigid elements.
C
      IF(IPARAM(11).EQ.3.or.eid.eq.0.or.rigide) GO TO 108
      WRITE(IOU13) EID,NPT,(KG(I),I=1,NPT)
108   CONTINUE
C
C     FILL CONNECTION TABLE IG.
C
      IF(npt.le.1) then
         npt=2
         KG(2)=0
      end if
      K=npt-1
      DO I=1,K
         L=I+1
         DO J=L,npt
            CALL SETIG(KG(I),KG(J),IG,II1,NORIG,IER)
            IF(IER.GT.0) RETURN
         end do
      end do
  130 IFLD=IFLD+4
      IF(LEN.EQ.1) GO TO 140
      IF(LOOP.EQ.1) GO TO 140
      IF(II.EQ.2) GO TO 140
      ENDP(1) = KA(69)
      ENDP(2) = KA(70)
      READ(IOU5,30,END=170) KA
      CALL LEFT(KA,ISAVE)
      IF(KA(1).NE.MB(4))   GO TO 150
      IF(KA(2).NE.ENDP(1)) GO TO 150
      IF(KA(3).NE.ENDP(2)) GO TO 150
      WRITE(IOU8,30) KA
      IFLD=4
140   CONTINUE
C
C     GO BACK AND READ ANOTHER CARD.
C
      GO TO 20
C
C     WARNING MESSAGE.
C
  150 WRITE(IOU6,160) KA
  160 FORMAT(/,'Warning.  The continuation to the card preceding'/
     - A1,A4,A3,65A1,A4,A3/'does not immediately follow its parent',
     - ' and hence will not be found by BANDIT.')
      ITYPE=ISAVE
      GO TO 55
C
C     END-OF-FILE ENCOUNTERED
C
170   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE REED

! ##################################################################################################################################
      SUBROUTINE RELABL(NS,NODES,IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,
     -                  ILD,IAJ,IDIM,IER)
C
C     GENERATE A RELABELING SCHEME STARTING WITH NS NODES FOR WHICH
C     LABELS HAVE BEEN STORED IN ARRAY NODES.
C     SET UP ILD AND NEW.
C     ILD(OLD)=NEW
C     NEW(NEW)=OLD,   THE INVERSE OF ILD
C
      DIMENSION IG(II1,*),IC(*),IDEG(*),IDIS(*),IW(*),NEW(*),ICC(*)
      INTEGER ILD(*)
      COMMON /S/ NN,DUMS(8)
      INTEGER X
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER IDIM,NODES(*),IAJ(IDIM)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,ICN    ,IDEG   ,
     &         IDIS   ,IER    ,IG     ,II     ,II1    ,IJ     ,IW     ,
     &         IZ     ,
     &         J      ,JJ     ,JT     ,
     &         KI     ,KMOD   ,KO     ,
     &         L      ,LL     ,
     &         MAXDEG ,MAXGRD ,
     &         N      ,NBITIN ,NEW    ,NN     ,NNC    ,NS     ,NT

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      I=NODES(1)
      ICN=IC(I)
      NT=ICC(ICN)-1
      DO 50 I=1,NN
      IF(IC(I)-ICN) 50,40,50
40    IDIS(I)=0
50    CONTINUE
      DO 100 J=1,NS
      JJ=NODES(J)
      IDIS(JJ)=-1
      JT=J+NT
      NEW(JT)=JJ
100   ILD(JJ)=JT
      KI=NT
      KO=NS+NT
      LL=KO
      L=1
      J=KO
      NNC=ICC(ICN+1)-1
130   KI=KI+1
      IF(KI-LL)135,132,135
132   L=L+1
      LL=KO+1
135   II=NEW(KI)
      N=IDEG(II)
      IF(N)140,255,140
140   IJ=0
      DO 200 I=1,N
      IA=IG(II,I)
CPACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IDIS(IA)) 200,150,200
150   IJ=IJ+1
      IF(IJ.GT.IDIM) GO TO 270
      IDIS(IA)=L
      KO=KO+1
      IAJ(IJ)=IA
      IW(IJ)=IDEG(IA)
200   CONTINUE
      IF(IJ-1)250,210,220
210   J=KO
      IZ=IAJ(1)
      NEW(KO)=IZ
      ILD(IZ)=KO
      GO TO 250
220   X=0
221   DO 230 I=2,IJ
      IF(IW(I)-IW(I-1))224,230,230
224   CONTINUE
      X=IW(I)
      IW(I)=IW(I-1)
      IW(I-1)=X
225   X=IAJ(I)
      IAJ(I)=IAJ(I-1)
      IAJ(I-1)=X
230   CONTINUE
      IF(X)235,235,220
235   DO 240 I=1,IJ
      J=J+1
      IZ=IAJ(I)
      NEW(J)=IZ
      ILD(IZ)=J
240   CONTINUE
250   IF(KO-NNC)130,255,255
255   CONTINUE
C
C     REVERSE SEQUENCE FOR THIS COMPONENT (ICN).
C
C     ICC IS AN ARRAY USED FOR IDENTIFYING COMPONENTS IN THE NEW ARRAY.
C     ICC(I) CONTAINS THE INDEX FOR THE NEW ARRAY AT WHICH COMPONENT I
C        STARTS.
C     NEW(I) = OLD LABEL FOR NODE NOW LABELLED I.
C     I = NUMBER OF NODES IN PREVIOUS COMPONENTS.
C     J = NUMBER OF NODES IN LATER COMPONENTS.
      I=ICC(ICN)-1
      J=NN-ICC(ICN+1) +1
      IF(J.GT.NN) J=0
C
C     GET REVERSE CM SEQUENCE.
C
      CALL REVERS(NEW,ILD,NN,I,J)
C
      RETURN
C
C     DIMENSION EXCEEDED.  STOP JOB.
C
270   CALL FINISH(5,IER)
      RETURN
      END SUBROUTINE RELABL

! ##################################################################################################################################
      SUBROUTINE RESET(IGOTO,ILD,INT,JUMP)
C
C     THIS ROUTINE MAKES A COPY OF THE ORIGINAL GRID POINT SEQUENCE IF
C     BOTH METHODS ARE TO BE USED.  ALSO, IF CUTHILL ACHIEVES NO BW
C     REDUCTION, THIS SEQUENCE IS RETRIEVED BEFORE CALLING GIBSTK.
C
      INTEGER ILD(*),INT(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IGOTO  ,IPARAM ,JUMP   ,L      ,NN

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      L=IOU12
      GO TO (10,20), IGOTO
C
C     SAVE ORIGINAL ORDERING (ILD) BEFORE CUTHILL IS CALLED.
C
   10 IF(IPARAM(3).NE. 9) RETURN
      REWIND L
      WRITE(L) (ILD(I),I=1,NN)
      RETURN
C
C     RESET ILD, IF NECESSARY, AND COPY TO INT.
C
   20 IF(IPARAM(3).NE. 9) GO TO 30
      REWIND L
C     IF CM MAKES NO IMPROVEMENT IN SEQUENCE, READ ORIGINAL SEQUENCE
C        BACK INTO ILD.
      IF(JUMP.NE.0) READ (L) (ILD(I),I=1,NN)
      REWIND L
c     ENDFILE L
c     REWIND L
C
   30 DO 40 I=1,NN
   40 INT(I)=ILD(I)
      RETURN
      END SUBROUTINE RESET

! ##################################################################################################################################
      SUBROUTINE REVERS(NEW,ILD,NN,N1,N2)
C
C     REVERSE THE NODAL SEQUENCE, OMITTING THE FIRST N1 AND THE LAST
C     N2 POINTS.
C
C     NEW(I) = OLD LABEL FOR NODE NOW LABELLED I (INPUT AND OUTPUT)
C     ILD(I) = NEW LABEL FOR NODE ORIGINALLY LABELED I (OUTPUT)
C     NN     = NUMBER OF NODES (INPUT)
C     N1     = NUMBER OF POINTS AT BEGINNING OF SEQUENCE TO OMIT FROM
C              REVERSAL (INPUT)
C     N2     = NUMBER OF POINTS AT END OF SEQUENCE TO OMIT FROM
C              REVERSAL (INPUT)
C
      INTEGER NEW(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,K      ,L      ,LL     ,M      ,
     &         N1     ,N2     ,NN

! E////////////////////////////////////////////////////////////////////E
C     J = NUMBER OF INTERCHANGES TO MAKE.
      J = (NN - N1 - N2)/2
      IF(J.LE.0) RETURN
      LL = NN - N2 + 1
C     MAKE INTERCHANGES IN NEW ARRAY.
      DO 10 I = 1 , J
      L=LL-I
      K=NEW(L)
      M=N1+I
      NEW(L) = NEW(M)
   10 NEW(M) = K
C     CORRECT ILD, THE INVERSE OF NEW.
      L=1+N1
      M=NN-N2
      DO 20 I=L,M
      K=NEW(I)
   20 ILD(K)=I
      RETURN
      END SUBROUTINE REVERS

! ##################################################################################################################################
      SUBROUTINE RIGID(KA,ITYPE,KG,MAXI,INV,II3,NORIG,LG,IER,npoint)
c
c     Extract grid points from rigid element card, and, in most cases,
c     generate equivalent MPCs.
c
C     CARD CONTENTS ARE STORED IN KA IN THE FORMAT A1,A4,A3,65A1,A4,A3
C
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with IG we will use IG_array
! instead. Needed so Lahey doesn't complain about shape of IG being
! different than IP (array) in subr READIT
      INTEGER IG_array(1)
! E////////////////////////////////////////////////////////////////////E

! B 02/21/04 //////////////////////////////////////////////////////////B
! When running RBE2-09-CBAR-10-A, I got msg that the subscript of array
! LG was out of range (tried to use LG(3) but LG only dimensioned 2.
! Gordon says to change dim to *
!!    INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2),LG(2)
      INTEGER MAXI
      INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2),LG(*)
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IBEG   ,IDUM   ,IEND   ,IER    ,IG     ,II3    ,
     &         IPARAM ,IT     ,ITYPE  ,
     &         K      ,KMOD   ,
     &         L      ,LOCDGP ,LOCIGP ,
     &         MA     ,MAXDEG ,MAXGRD ,MB     ,
     &         N      ,NDGP   ,NEQ    ,NEQR   ,NF     ,NFSDGP ,NLINK  ,
     &         npoint ,NUM

      REAL     DUMD

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /D/ DUMD(7),NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
C     DETERMINE RIGID ELEMENT TYPE.
C
C     IT = 1  CRIGD1 ON LEVEL 16 (of Nastran)
C          2  CRIGD2 ON LEVEL 16
C          3  CRIGD1 ON LEVEL 17 AND ABOVE
C          4  CRIGD2 ON LEVEL 17 AND ABOVE
C          5  CRIGD1 ON LEVEL 17 AND ABOVE (THRU OPTION)
c          6  CRBE1 or RBE1 (read here but treated like regular element)
c          7  CRBE2 or RBE2
c
c     This routine processes the elements listed above.
c     CRIGD1, CRIGD2, and RBE2 generate equivalent MPCs.
c     RBE1 is read in this routine because of the complicated syntax of
c     the element, but RBE1 is otherwise treated like a regular element.
c     The following rigid elements are processed and also handled like
c     regular elements:
c     CRIGDR(93), CRBAR(128), RBAR(129), CRROD(134), RROD(135),
c     CRTRPLT(136), RTRPLT(137).
C     No other rigid elements are currently recognizable.
C
C     CHECK FOR BLANK FIELD 3 TO DISTINGUISH BETWEEN LEVEL 16 AND
C     LEVEL 17 FORMATS FOR CRIGD1 AND CRIGD2
C     (WHICH ARE DIFFERENT, BELIEVE IT OR NOT)
C
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(12),1,8,IG_array(1),I)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
C
      IT=0
      IF(ITYPE.EQ.91.AND.I.EQ.0) IT=1
      IF(ITYPE.EQ.92.AND.I.EQ.0) IT=2
      IF(ITYPE.EQ.91.AND.I.EQ.1) IT=3
      IF(ITYPE.EQ.92.AND.I.EQ.1) IT=4
      IF(IT.EQ.3.AND.NTHRU(KA(28),8).EQ.1) IT=5
      IF(ITYPE.EQ.130.or.itype.eq.131) IT=6
      IF(ITYPE.EQ.132.or.itype.eq.133) IT=7
      IF(IT.EQ.0) RETURN
C
C     SET PARAMETERS FOR INTERPRETING physical CARD
c     (the first card of the logical card)
C
C     LOCIGP=LOCATION IN KA ARRAY OF INDEPENDENT GRID POINT
C     LOCDGP=LOCATION IN KA ARRAY OF FIRST DEPENDENT GRID POINT ON CARD
C     NDGP  =NUMBER OF DEP. G.P. ON CARD
C     NFSDGP=NUMBER OF FIELDS SEPARATING DEP. G.P. ON CARD
C
      GO TO (40,42,44,46,48,50,52), IT
   40 LOCIGP=20
      LOCDGP=28
      NDGP=5
      NFSDGP=1
      GO TO 75
   42 LOCIGP=20
      LOCDGP=36
      NDGP=2
      NFSDGP=2
      GO TO 75
   44 LOCIGP=12
      LOCDGP=20
      NDGP=6
      NFSDGP=1
      GO TO 75
   46 LOCIGP=12
      LOCDGP=20
      NDGP=3
      NFSDGP=2
      GO TO 75
   48 LOCIGP=12
      LOCDGP=20
      NDGP=2
      NFSDGP=2
      GO TO 75
   50 LOCIGP=12
      LOCDGP=28
      NDGP=2
      NFSDGP=2
      GO TO 75
   52 LOCIGP=12
      LOCDGP=28
      NDGP=5
      NFSDGP=1
      GO TO 75
   75 CONTINUE
C
C     GET INDEPENDENT GRID POINT IG.
C
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(LOCIGP),1,8,IG_array(1),N)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
C
C     IS IG THERE?
      IF(IG.eq.0) then
         WRITE(IOU6,20)
20       FORMAT(/,'Fatal Error.  Independent grid point on rigid',
     -          ' element not found.')
         WRITE(IOU6,295) KA
         call finish(6,IER)
         RETURN
      end if
      kg(1)=ig
      npoint=1
C
C     CONVERT IG TO INTERNAL LABEL.
C
! B////////////////////////////////////////////////////////////////////B
! Modify call to SCAT so that 1st arg is an array (like in subr SCAT)
      IG_array(1) = IG
      CALL SCAT(IG_array(1),1,INV,II3,NORIG,IER)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(IER.GT.0) RETURN
C
C     INTERPRET DEPENDENT GRID POINTS ON PARENT OR CONTINUATION CARD.
C
90    CONTINUE
      DO I=1,NDGP
         K=LOCDGP+(I-1)*8*NFSDGP
         CALL READIT(KA(K),1,8,LG(I),NF) ! This is OK. LG is an array
      end do
C     ELIMINATE ZEROS OR DUPLICATE ENTRIES FROM LIST OF DEP. G.P.
      NF=NDGP
      CALL FIXIT(LG,NF)
C     SET UP FOR NEXT CARD.
      GO TO (110,112,114,116,118,124,126), IT
  110 LOCDGP=4
      NDGP=8
      GO TO 150
  112 LOCDGP=4
      NDGP=4
      GO TO 150
  114 LOCDGP=4
      NDGP=8
      GO TO 150
  116 LOCDGP=4
      NDGP=4
      GO TO 150
C     LIST OUT POINTS INCLUDED ON THRU OPTION.
  118 CONTINUE
      IF(NF.LT.2) GO TO 120
      IBEG=LG(1)
      IEND=LG(2)
      N=IEND-IBEG
      IF(N.LT.1) GO TO 120
      IF(N.GT.(MAXGRD-1)) GO TO 120
      DO I=1,N
         LG(I+1)=IBEG+I
      end do
      NF=N+1
      GO TO 150
  120 WRITE(IOU6,121) KA
  121 FORMAT(/,'Fatal Error.  Illegal data on'/A1,A4,A3,65A1,A4,A3)
      call finish(6,IER)
      RETURN
  124 LOCDGP=12
      NDGP=3
      GO TO 150
  126 LOCDGP=4
      NDGP=8
      GO TO 150
150   CONTINUE
c     add points to KG list containing all points on rigid element
      if((npoint+nf).gt.maxi) then
         Write(iou6,155)
155      format(/,'Fatal error.  Insufficient space for rigid element.')
         call finish(2,ier)
         return
      end if
      do i=1,nf
         kg(npoint+i)=LG(i)
      end do
      npoint=npoint+nf
      IF(NF.LE.0) GO TO 235
C
C     CONVERT new DEPENDENT POINTS TO INTERNAL LABELS
c     (except those rigid elements which are read in this routine but
c     treated like regular elements: RBE1)
C
      if(it.ne.6) then
         CALL SCAT(LG,NF,INV,II3,NORIG,IER)
         IF(IER.GT.0) RETURN
C
C        WRITE OUT EQUIVALENT MPCs ON IOU12
C
         DO L=1,NF
c           KG(1)=LG(L)
c           KG(2)=IG
c           I=2
c           WRITE(IOU12) I,(KG(J),J=1,I)
            Write(iou12) 2,LG(L),IG
C           INCREMENT OVERALL MPC EQUATION COUNTER.
            NEQ=NEQ+1
C           INCREMENT COUNTER OF MPC EQUATIONS ARISING FROM RIGID ELEMENTS.
            NEQR=NEQR+1
         end do
      end if
      IF(IT.EQ.5) RETURN
C
C     READ ANOTHER CARD.
C
C     SAVE OLD FIELD 10 FOR COMPARISON TO NEW FIELD 1.
  235 ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,240,END=310) KA
  240 FORMAT(A1,A4,A3,65A1,A4,A3)
C     LEFT-ADJUST FIELD 1.
      CALL LEFT(KA,IDUM)
C     CHECK IF CONTINUATION CARD.
      DO I=3,4
         IF(KA(1).EQ.MB(I)) GO TO 290
      end do
C     NO CONTINUATION.  END OF LOGICAL CARD FOUND.
  280 BACKSPACE IOU5
C     ELIMINATE ZEROS and DUPLICATE ENTRIES FROM LIST OF grids
      call fixit(kg,npoint)
      RETURN
C
C     CONTINUATION FOUND.
  290 CONTINUE
      IF(KA(2).NE.ENDP(1)) GO TO 280
      IF(KA(3).NE.ENDP(2)) GO TO 280
C     IS MATCHING CONTINUATION A LONG-FIELD CARD . . .
      IF(I.NE.5) GO TO 296
      WRITE(IOU6,293)
  293 FORMAT(/,'Fatal Error.  BANDIT cannot read a long field',
     -       ' continuation to a rigid element card.')
      WRITE(IOU6,295) KA
295   FORMAT(A1,A4,A3,65A1,A4,A3)
      call finish(6,IER)
      RETURN
  296 CONTINUE
C     WRITE CARD OUT ON UNIT 8.
      WRITE(IOU8,240)  KA
C     GO BACK AND INTERPRET THIS CARD.
      GO TO 90
C
C     END-OF-FILE ENCOUNTERED
C
310   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE RIGID

! ##################################################################################################################################
      SUBROUTINE RSETUP(LVL,LVLS1,LVLS2,NACUM,IDIM,IER)
C
C     SETUP COMPUTES THE REVERSE LEVELING INFO FROM LVLS2 AND STORES
C     IT INTO LVLS2.  NACUM(I) IS INITIALIZED TO NODES/ITH LEVEL FOR
C     NODES ON THE PSEUDO-DIAMETER OF THE GRAPH.  LVL IS INITIALIZED
C     TO NONZERO FOR NODES ON THE PSEUDO-DIAM AND NODES IN A DIFFERENT
C     COMPONENT OF THE GRAPH.
C
      COMMON /GRA/ N,IDPTH,DUMG
      INTEGER IDIM
      INTEGER NACUM(IDIM),LVL(*),LVLS1(*),LVLS2(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDPTH  ,IER    ,ITEMP  ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
C     IDIM=NUMBER OF LEVELS IN A GIVEN COMPONENT.
      IF(IDPTH.LE.IDIM)  GO TO 20
C     DIMENSION EXCEEDED  . . .  STOP JOB.
      CALL FINISH(3,IER)
      RETURN
   20 CONTINUE
      DO 30 I=1,IDPTH
        NACUM(I)=0
   30 CONTINUE
      DO 140 I=1,N
        LVL(I)=1
        LVLS2(I)=IDPTH+1-LVLS2(I)
        ITEMP=LVLS2(I)
        IF(ITEMP.GT.IDPTH) GO TO 140
        IF(ITEMP.NE.LVLS1(I)) GO TO 100
        NACUM(ITEMP)=NACUM(ITEMP)+1
        GO TO 140
  100   LVL(I)=0
  140 CONTINUE
      RETURN
      END SUBROUTINE RSETUP

! ##################################################################################################################################
      SUBROUTINE SCAT(KG,NCON,INV,II3,NORIG,IER)
C
C     THIS ROUTINE USES SCATTER SORT TECHNIQUES FOR EACH GRID POINT
C     ENCOUNTERED TO DETERMINE WHETHER OR NOT THE POINT HAS BEEN
C     SEEN BEFORE.  IF NOT, INV, NORIG, AND NN ARE UPDATED.
C
C     INV(1,I) CONTAINS AN ORIGINAL GRID POINT NUMBER
C     INV(2,I) CONTAINS THE INTERNAL NUMBER ASSIGNED TO IT (BEFORE
C              SORTING)
C
      INTEGER II3,INV(2,II3),NORIG(*),KG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IER    ,KMOD   ,LOC    ,
     &         MAXDEG ,MAXGRD ,NCON   ,NN     ,NOLD

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,DUMS(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      IF(NCON.LT.1) RETURN
      DO 50 I=1,NCON
      NOLD=KG(I)
      IF(NOLD.EQ.0) GO TO 50
      LOC=NOLD-1
10    LOC=MOD(LOC,KMOD)+1
20    IF(INV(1,LOC).NE.0) GO TO 30
      INV(1,LOC)=NOLD
      NN=NN+1
      IF(NN.GT.MAXGRD) GO TO 60
      NORIG(NN)=NOLD
      INV(2,LOC)=NN
      GO TO 40
30    IF(INV(1,LOC).NE.NOLD) GO TO 10
40    KG(I)=INV(2,LOC)
50    CONTINUE
      RETURN
60    WRITE(IOU6,70) MAXGRD
70    FORMAT(/,'Fatal Error.  Grid point limit of',I9,' exceeded.')
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE SCAT

! ##################################################################################################################################
      SUBROUTINE SEQGP(NORIG,ILD,NEW,JUMP)
C
! B////////////////////////////////////////////////////////////////////B
      USE IOUNT1, ONLY                :  SEQ
! E////////////////////////////////////////////////////////////////////E
C     WRITE SEQGP BULK DATA CARDS ON 7 and 8.
C
C     NN       = NUMBER OF GRID POINTS
C     NORIG(I) = ORIGINAL GRID POINT CORRESPONDING TO BANDIT INTERNAL
C                LABEL I
C     NEW(I)   = SCRATCH ARRAY
C     ILD(I)   = NEW RESEQUENCED LABEL CORRESPONDING TO BANDIT INTERNAL
C                LABEL I
C     JUMP     = 1 IF NO SEQGP CARDS ARE TO BE GENERATED (E.G., IF NO
C                IMPROVEMENT IN BW HAS RESULTED).
C
      INTEGER NORIG(*),ILD(*),NEW(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IADD   ,ILOOP  ,IPARAM ,ITIME  ,
     &         J      ,JUMP   ,
     &         MAXW   ,MM     ,
     &         NAXIC  ,NBW    ,NLOOP  ,NN     ,NP     ,NRMS

      REAL     DMY    ,DUM    ,DUMD   ,DUMS   ,DUMW   ,RMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DMY(5),IADD,DUMS,NAXIC
      COMMON /B/ IPARAM(20)
      COMMON /D/ OBW,NBW,OP,NP,DUMD(6)
      COMMON /W/ DUM(2),MAXW,RMS,DUMW(2)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER OBW,OP
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      IF(NN.LE.0) RETURN
      WRITE(IOU6,2)
2     FORMAT(/,'Field 10 of first SEQGP card shows new grid',
     -       ' point rms wavefront.')
      IF(IADD.NE.0) WRITE(IOU6,3) IADD
    3 FORMAT(/,'Integer added to new sequence numbers',I8)
C
      NRMS=RMS+0.5
C     NRMS = NEW RMS WAVEFRONT ROUNDED TO NEAREST INTEGER
      IF(JUMP.NE.1) GO TO 9
C
      WRITE(IOU6,7)
    7 FORMAT(/,'No SEQGP cards generated.')
      REWIND IOU7
      ENDFILE IOU7
! B////////////////////////////////////////////////////////////////////B
      CLOSE (SEQ, STATUS='DELETE')
! E////////////////////////////////////////////////////////////////////E
      GO TO 60
C
    9 CONTINUE
      WRITE(IOU6,13)
13    FORMAT(/,'See bandit.f07 for SEQGP cards generated.')
C
C     IF SPRINGS ARE REQUESTED WITH $SPRING YES, THE SEQGP CARDS
C     RESEQUENCE THE CONSECUTIVELY NUMBERED OLD NUMBERS RATHER THAN
C     THE ORIGINAL GRID NUMBERS ASSIGNED BY THE USER.
C
      IF(IPARAM(7).EQ.3) GO TO 150
      WRITE(IOU6,125)
  125 FORMAT('Since scalar springs were requested with $SPRING YES,',
     - ' the SEQGP cards resequence the consecutively numbered old'/
     - 'numbers rather than the user-assigned original grid number.'/)
      REWIND iou16
      READ(iou16,130) I,MM
  130 FORMAT(24I5)
C     DETERMINE REPLACEMENT NORIG ARRAY.
      DO I=1,NN
         READ(iou16,130) NORIG(I),(NEW(J),J=1,MM)
C        NEW IS DUMMY SPACE HERE TO INSURE THAT THE RIGHT NUMBER OF
c        CARDS IS READ.
      end do
      REWIND iou16
  150 CONTINUE
C
C     ADD THE NON-NEGATIVE INTEGER IADD TO THE NEW SEQUENCE NUMBERS.
C
      DO I=1,NN
         ILD(I)=ILD(I)+IADD
      end do
C
      NLOOP=0
      IF(NAXIC.NE.99999) NLOOP=IABS(NAXIC)
      DO 30 ILOOP=0,NLOOP
         IF(NAXIC.NE.99999) THEN
            DO I=1,NN
               NORIG(I)=NORIG(I)+1000000
               ILD(I)=ILD(I)+1000000
            end do
            IF(NAXIC.LT.0.AND.ILOOP.NE.NLOOP) GO TO 30
         END IF
         IF(NN.LE.4) THEN
            WRITE(IOU8,10) (NORIG(I),ILD(I),I=1,NN)
            WRITE(IOU7,10) (NORIG(I),ILD(I),I=1,NN)
! B////////////////////////////////////////////////////////////////////B
            WRITE(SEQ, 10) (NORIG(I),ILD(I),I=1,NN)
! E////////////////////////////////////////////////////////////////////E
         ELSE
            WRITE(IOU8,10) (NORIG(I),ILD(I),I=1,4),NRMS,
     -                     (NORIG(J),ILD(J),J=5,NN)
            WRITE(IOU7,10) (NORIG(I),ILD(I),I=1,4),NRMS,
     -                     (NORIG(J),ILD(J),J=5,NN)
! B////////////////////////////////////////////////////////////////////B
            WRITE(SEQ, 10) (NORIG(I),ILD(I),I=1,4),NRMS,
     -                     (NORIG(J),ILD(J),J=5,NN)
! E////////////////////////////////////////////////////////////////////E
10       FORMAT('SEQGP*          ',8I16,I16/('SEQGP*          ',8I16))
         END IF
30    CONTINUE
C
60    WRITE(IOU8,70)
70    FORMAT('ENDDATA')
      REWIND IOU7
      RETURN
      END SUBROUTINE SEQGP

! ##################################################################################################################################
      SUBROUTINE SETIG(KG1,KG2,IG,II1,NORIG,IER)
C
C     THIS ROUTINE SETS IG(KG1,-)=KG2 AND IG(KG2,-)=KG1 IF THIS
C        CONNECTION HAS NOT ALREADY BEEN SET.
C
C     NEDGE = NUMBER OF UNIQUE EDGES.
      INTEGER II1
      DIMENSION IG(II1,*),NORIG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IER    ,IG     ,IS     ,
     &         K      ,KG1    ,KG2    ,KMOD   ,
     &         L      ,LOOP   ,
     &         M      ,MAXDEG ,MAXGRD ,MM     ,
     &         NBITIN ,NEDGE  ,NN     ,NORIG

      REAL     DUM    ,DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,DUMS(3)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C
      IF(KG1.EQ.0 .OR. KG2.EQ.0 .OR. KG1.EQ.KG2) RETURN
      L=KG1
      K=KG2
      DO 50 LOOP=1,2
      IF(LOOP.EQ.1) GO TO 20
      L=KG2
      K=KG1
   20 M=0
   30 M=M+1
      IF(M.GT.MAXDEG) GO TO 60
      IS=IG(L,M)
CPACK IS=IUNPK(IG,MAXGRD*(M-1)+L,NBITIN)
      IF(IS.EQ.0) GO TO 40
      IF(IS.NE.K) GO TO 30
      GO TO 55
40    CONTINUE
      IG(L,M)=K
CPACK CALL PACK(IG,MAXGRD*(M-1)+L,NBITIN,K)
      MM=MAX(MM,M)
      IF(LOOP.EQ.1) NEDGE = NEDGE + 1
   50 CONTINUE
   55 RETURN
C
! B////////////////////////////////////////////////////////////////////B
   60 WRITE(IOU6,70) NORIG(L),M,MAXDEG
70    FORMAT(/,'Fatal Error.  The number of connections for point',I9/
     - ' = ',I9,' and exceeds the nodal degree limit of',I8)
! E////////////////////////////////////////////////////////////////////E
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE SETIG

! ##################################################################################################################################
      SUBROUTINE SETUP(KA,NCON,IFLD,KG,NPT,LEN,LESSOK,IER)
C
C     SET UP ARRAY KG CONTAINING THE LIST OF CONNECTIONS FOR AN ELEMENT.
C
C     DATA CARD ARRAY KA WAS READ A1,A4,A3,64A1,A1,A4,A3
C
C     NCON   = NUMBER OF CONNECTIONS (INPUT)
C     IFLD   = FIELD NUMBER OF FIRST CONNECTION (INPUT)
C     KG     = LIST OF CONNECTIONS (OUTPUT)
C     NPT    = ACTUAL NUMBER OF NON-ZERO CONNECTIONS FOUND (OUTPUT)
C     LEN    = 1 FOR SHORT FIELD DATA CARD, 2 FOR  LONG FIELD DATA CARD
C              (INPUT)
C     LESSOK = .TRUE. IF THE CONNECTION CARD NEED NOT HAVE ALL GRID
C              POINTS PRESENT (E.G., CELAS1 OR CPENTA), OTHERWISE FALSE
C              (INPUT)
C
      INTEGER KA(70),KG(*),ENDP(2)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,      IER    ,IFLD   ,IPARAM ,ITYPE  ,
     &         J      ,K      ,L      ,LEN    ,LEN8   ,
     &         M      ,M1     ,MA     ,MB     ,
     &         N      ,NCON   ,NIP    ,NPT    ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      LOGICAL LESSOK
C
      DO I=1,NCON
         KG(I)=0
      end do
      N=0
      NPT=0
      LEN8=LEN*8
      J = 4 + 8 * (IFLD- 2 ) * LEN
      L = 14 - IFLD - 4 * LEN
C
C     N=NO. OF CONNECTION FIELDS PROCESSED SO FAR
C     J=SUBSCRIPT IN KA CORRESPONDING TO FIRST CONNECTION FIELD ON CARD
C     L=NO. OF FIELDS FROM FIRST CONNECTION FIELD ON CARD TO END OF CARD
C     K=MAX. NO. OF CONNECTIONS POSSIBLE ON A GIVEN PHYSICAL CARD
C
20    K=MIN(L,NCON-N)
      DO M=1,K
         M1=J+LEN8*(M-1)
         I=NPT+1
         CALL READIT(KA(M1),1,LEN8,KG(I),NIP)
         IF(KG(I).NE.0) NPT=NPT+1
      end do
      N=N+K
      IF(N.GE.NCON) RETURN
C     SAVE ENDPUNCHING (FIELD 10).
      ENDP(1)=KA(69)
      ENDP(2)=KA(70)
C     READ NEXT CARD.
      READ(IOU5,40,END=100) KA
40    FORMAT(A1,A4,A3,65A1,A4,A3)
C     LEFT-ADJUST FIELD 1.
      CALL LEFT(KA,ITYPE)
C     CHECK NEW FIELD 1 AGAINST PREVIOUS FIELD 10 FOR MATCH.
      DO I=3,4
         IF(KA(1).EQ.MB(I)) GO TO 60
      end do
      GO TO 70
60    IF(KA(2).NE.ENDP(1)) GO TO 70
      IF(KA(3).NE.ENDP(2)) GO TO 70
      WRITE(IOU8,40) KA
      LEN=1
      IF(I.EQ.5) LEN=2
      L=8/LEN
      J=4
      GO TO 20
C     END OF LOOP STARTING AT STATEMENT 20.
C
C     ABORT JOB IF CONNECTION CARD CONTINUATIONS DO NOT IMMEDIATELY
C        FOLLOW THEIR PARENTS, UNLESS IT IS FOR AN ELEMENT FOR WHICH
C        IT IS OK IF NOT ALL CONNECTIONS ARE PRESENT (E.G., CELAS1
C        OR CPENTA)
C
70    CONTINUE
      IF(LESSOK) GO TO 90
      WRITE(IOU6,80) KA
80    FORMAT(/,'Fatal Error.  The following data card is out of sort.'/
     - A1,A4,A3,65A1,A4,A3/
     - 'Since BANDIT does not sort the data, continuations must'/
     - 'immediately follow their parents.')
      call finish(6,IER)
      RETURN
C
90    CONTINUE
      BACKSPACE IOU5
      RETURN
C
C     END-OF-FILE ENCOUNTERED
C
100   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE SETUP

! ##################################################################################################################################
      SUBROUTINE SORT(LIST,NL)
C
C     SORT A LIST OF INTEGERS OF LENGTH NL.  THIS ROUTINE OPERATES
C        FASTEST FOR THOSE LISTS NOT BADLY OUT OF SORT.
C
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,J1     ,K      ,KFLAG  ,L      ,
     &         NL     ,NL1

! E////////////////////////////////////////////////////////////////////E
      IF(NL.LE.1) RETURN
      NL1=NL-1
      DO 20 I=1,NL1
      K=NL-I
      KFLAG=0
      DO 10 J=1,K
      J1=J+1
      IF(LIST(J).LE.LIST(J1))  GO TO 10
      KFLAG=1
      L=LIST(J)
      LIST(J)=LIST(J1)
      LIST(J1)=L
10    CONTINUE
      IF(KFLAG.EQ.0) RETURN
20    CONTINUE
      RETURN
      END SUBROUTINE SORT

! ##################################################################################################################################
      SUBROUTINE SORTDG(STK1,STK2,X1,X2,NDEG)
C
C     SORTDG SORTS STK2 BY DEGREE OF THE NODE AND ADDS IT TO THE END
C     OF STK1 IN ORDER OF LOWEST TO HIGHEST DEGREE.  X1 AND X2 ARE THE
C     NUMBER OF NODES IN STK1 AND STK2 RESPECTIVELY.
C
      COMMON /GRA/ N,IDPTH,DUMG
      INTEGER NDEG(*),STK1(*),STK2(*),X1,X2,TEMP
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDPTH  ,IND    ,ISTK2  ,ITEST  ,
     &         J      ,JSTK2  ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
      IND=X2
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) GO TO 40
      DO 30 I=1,IND
        J=I+1
        ISTK2=STK2(I)
        JSTK2=STK2(J)
        IF(NDEG(ISTK2).LE.NDEG(JSTK2)) GO TO 30
        ITEST=1
        TEMP=STK2(I)
        STK2(I)=STK2(J)
        STK2(J)=TEMP
   30 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
   40 DO 50 I=1,X2
        X1=X1+1
        STK1(X1)=STK2(I)
   50 CONTINUE
      RETURN
      END SUBROUTINE SORTDG

! ##################################################################################################################################
      INTEGER FUNCTION SORT2(XC,SIZE,STPT)
C
C     SORT2 SORTS SIZE AND STPT INTO DESENDING ORDER ACCORDING TO
C     VALUES OF SIZE.
C
C     XC=NUMBER OF ENTRIES IN EACH ARRAY
C
      INTEGER SIZE(*),STPT(*),TEMP,XC
! B////////////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IND    ,ITEST  ,J

! E////////////////////////////////////////////////////////////////////E
C     THE DIMENSION OF SIZE AND STPT IS THE MAXIMUM ALLOWABLE NUMBER
C            OF CONNECTED COMPONENTS.
      SORT2=0
      IF(XC.EQ.0) RETURN
      SORT2=1
      IND=XC
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) RETURN
      DO 17 I=1,IND
        J=I+1
        IF(SIZE(I).GE.SIZE(J)) GO TO 17
        ITEST=1
        TEMP=SIZE(I)
        SIZE(I)=SIZE(J)
        SIZE(J)=TEMP
        TEMP=STPT(I)
        STPT(I)=STPT(J)
        STPT(J)=TEMP
   17 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
      RETURN
      END FUNCTION SORT2

! ##################################################################################################################################
      SUBROUTINE SPRING(IP)
C
C     GENERATE SCALAR SPRING ELEMENTS CONSISTENT WITH THE CONNECTIVITY
C     MATRIX IG.  CELAS3 ELEMENTS ARE WRITTEN ON UNIT 9.
C     IP = TEMPORARY STORAGE.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,II     ,J      ,MM     ,NN

      REAL     STIFF

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER IP(*),EID,PID
      PID=101
      STIFF=2.71828
      REWIND IOU9
      REWIND iou16
      EID=0
      READ(iou16,10) NN,MM
      IF((NN*MM).EQ.0) GO TO 55
C
C     GENERATE SPRINGS BETWEEN NODES.
C
      DO 25 II=1,NN
         READ(iou16,10) I,(IP(J),J=1,MM)
10       FORMAT(24I5)
         DO 20 J=1,MM
            IF(IP(J).LE.0) GO TO 25
            IF(IP(J).LT.I) GO TO 20
            EID=EID+1
            WRITE(IOU9,15) EID,PID,I,IP(J)
15          FORMAT('CELAS3  ',4I8)
20       CONTINUE
25    CONTINUE
C
C     INSERT BLANK CARD FOR LEVY'S WAVEFRONT PROGRAM.
C
      WRITE(IOU9,'(1x)')
C
C     GENERATE SPRINGS TO GROUND.
C
      PID=PID+1
      DO I=1,NN
         EID=EID+1
         WRITE(IOU9,35) EID,PID,I
35       FORMAT('CELAS3  ',3I8)
      end do
C
C     CREATE PROPERTY CARD.
C
      PID=PID-1
      WRITE(IOU9,50) PID,STIFF,NN,MM
   50 FORMAT('PELAS   ',I8,F8.5,48X,2I4)
      PID=PID+1
      WRITE(IOU9,50) PID,STIFF
C
   55 CONTINUE
      REWIND IOU9
      REWIND iou16
      WRITE(IOU6,60) EID,IOU9
   60 FORMAT(/I8,' CELAS3 spring elements generated on bandit.f',I2)
      RETURN
      END SUBROUTINE SPRING

! ##################################################################################################################################
      SUBROUTINE STABLE(IG,II1,IC,IDEG,ILD,IP)
C
C     PRINT CONNECTION TABLE IN TERMS OF SORTED INTERNAL LABELS.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IC     ,IDEG   ,IG     ,II1    ,ILD    ,IP     ,
     &         J      ,K      ,KMOD   ,
     &         MAXDEG ,MAXGRD ,MDIFF  ,MM     ,
     &         NBITIN ,NN

      REAL     DUM    ,DUMBB

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,MM,DUM(7)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DIMENSION IG(II1,*),IC(*),IDEG(*),ILD(*),IP(*)
      IF((NN*MM).EQ.0) RETURN
      WRITE(IOU6,10)
   10 FORMAT(/,'Sorted Internal Label Connection Table:'/
     - 8X,'Sort       Max'/' Label Label Comp Diff Degr  Connections')
      DO 40 I=1,NN
      DO 20 J=1,MM
   20 IP(J)=0
      MDIFF=0
      DO 30 J=1,MM
      K=IG(I,J)
CPACK K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(K.EQ.0) GO TO 40
      MDIFF=MAX(MDIFF,IABS(ILD(I)-ILD(K)))
   30 IP(J)=ILD(K)
   40 WRITE(IOU6,50) I,ILD(I),IC(I),MDIFF,IDEG(I),(IP(J),J=1,MM)
   50 FORMAT(2I6,24I5/(27X,21I5))
      RETURN
      END SUBROUTINE STABLE

! ##################################################################################################################################
      SUBROUTINE STACK(IDEG,NEW,ILD,IW)
C
C     STACK POINTS OF ZERO DEGREE AT END OF THE NUMBERING.
C
      INTEGER IDEG(*),NEW(*),ILD(*),IW(*)
C     IW IS SCRATCH STORAGE.
      COMMON /S/ NN,DUMS(8)
      COMMON /D/ DUM(5),KT,DUMD(4)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,K      ,KT     ,L      ,NN     ,NN1

      REAL     DUM    ,DUMD   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      KT=0
      NN1=NN-1
C     LIST POINTS OF ZERO DEGREE AND INCREMENT COUNTER KT.
      DO 10 I=1,NN
      IF(IDEG(I).GT.0) GO TO 10
      KT=KT+1
      IW(KT)=ILD(I)
   10 CONTINUE
      IF(KT.LE.0) GO TO 70
C     SORT LIST OF RENUMBERED NUMBERS TO BE STACKED.
      CALL SORT(IW,KT)
C     STACK POINTS OF ZERO DEGREE AT END OF NEW.
      DO 40 L=1,KT
      I=IW(L)-L+1
      K=NEW(I)
      IF(I.GE.NN) GO TO 30
      DO 20 J=I,NN1
   20 NEW(J)=NEW(J+1)
   30 NEW(NN)=K
   40 CONTINUE
C     CORRECT ILD, THE INVERSE OF NEW.
   70 DO 80 I=1,NN
      K=NEW(I)
   80 ILD(K)=I
      RETURN
      END SUBROUTINE STACK

! ##################################################################################################################################
! B////////////////////////////////////////////////////////////////////B
      SUBROUTINE SUMUP ( NEW_BW, DEN )
! E////////////////////////////////////////////////////////////////////E
C
C     PRINT BANDIT SUMMARY.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  IADD   ,IPARAM ,
     &         MAXW0  ,MAXW1  ,MINDEG ,MM     ,
     &         NBW    ,NCM    ,NEDGE  ,NEL    ,NEQ    ,NEQR   ,NLINK  ,
     &         NN     ,NONZ   ,NP     ,NZERO

! B////////////////////////////////////////////////////////////////////B
      INTEGER  NEW_BW
! E////////////////////////////////////////////////////////////////////E
      REAL     AN     ,ANN    ,AV1    ,AV2    ,
     &         BRMS0  ,BRMS1  ,
     &         DEN    ,DUM    ,DUMY   ,
     &         rawf   ,rbw    ,RMS0   ,RMS1   ,rmwf   ,rp     ,
     &         rrbw   ,rrwf

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,IADD,MINDEG,DUMY
      COMMON /B/ IPARAM(20)
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER OBW,OP
! B////////////////////////////////////////////////////////////////////B
      NEW_BW = NBW         ! New bandwidth printed out in this subr
! E////////////////////////////////////////////////////////////////////E
      ANN=FLOAT(NN)
      AV1=FLOAT(OP)/ANN
      AV2=FLOAT(NP)/ANN
      NONZ=2*NEDGE+NN
      AN=ANN*ANN
      DEN=FLOAT(NONZ)*100./AN
      NEQ=NEQ-NEQR
      rbw=float(nbw)/float(obw)
      rp=float(np)/float(op)
      rmwf=float(maxw1)/float(maxw0)
      rawf=av2/av1
      rrwf=rms1/rms0
      rrbw=brms1/brms0
      WRITE(IOU6,20) OBW,NBW,rbw,OP,NP,rp,MAXW0,MAXW1,rmwf,
     -               AV1,AV2,rawf,RMS0,RMS1,rrwf,BRMS0,BRMS1,rrbw
   20 FORMAT(/,'BANDIT Summary:'/39X,'Before',7X,'After',5x,'Ratio'/
     - 5X,'Bandwidth (B)',15X,2I12,f10.3/
     - 5X,'Profile (P)', 17X,2I12,f10.3/
     - 5X,'Maximum Wavefront (C-MAX)',3X,2I12,f10.3/
     - 5X,'Average Wavefront (C-AVG)',3X,2F12.3,f10.3/
     - 5X,'RMS Wavefront (C-RMS)',7X,2F12.3,f10.3/
     - 5X,'RMS Bandwidth (B-RMS)',7X,2F12.3,f10.3/)
      IF(IPARAM(6).EQ.10) WRITE(IOU6,24)
      IF(IPARAM(6).EQ.11) WRITE(IOU6,26)
      IF(IPARAM(6).EQ.12) WRITE(IOU6,28)
      IF(IPARAM(6).EQ.13) WRITE(IOU6,30)
  30  FORMAT(5X,'Criterion',30X,'Max Wavefront')
  24  FORMAT(5X,'Criterion',34X,'Bandwidth')
  26  FORMAT(5X,'Criterion',36X,'Profile')
  28  FORMAT(5X,'Criterion',30X,'RMS Wavefront')
      IF(IPARAM(3).EQ.7) WRITE(IOU6,32)
      IF(IPARAM(3).EQ.8) WRITE(IOU6,34)
      IF(IPARAM(3).EQ.9) WRITE(IOU6,36)
   36 FORMAT(5X,'Method Selected',27X,'CM and GPS')
   32 FORMAT(5X,'Method Selected',35X,'CM')
   34 FORMAT(5X,'Method Selected',34X,'GPS')
      WRITE(IOU6,60) NN,NEL,NLINK,NCM,MM,MINDEG,NEDGE,DEN,NZERO,NEQR,NEQ
   60 FORMAT(5X,'Number of Grid Points (N)',I27/
     -       5X,'Number of Elements (Non-Rigid)',I22/
     -       5X,'Number of Rigid Elements Processed',I18/
     -       5X,'Number of Components',I32/
     -       5X,'Maximum Nodal Degree',I32/
     -       5X,'Minimum Nodal Degree',I32/
     -       5X,'Number of Unique Edges',I30/
     -       5X,'Matrix Density in Percent',F27.4/
     -       5X,'Number of Points of Zero Degree',I21/
     -       5X,'Number of Rigid Element MPC Equations',I15/
     -       5X,'Number of MPC Equations Processed',I19)
      WRITE(IOU6,80)
80    FORMAT(/,'All BANDIT statistics use grid point, rather than',
     - ' DOF, connectivity and'/
     - 'include matrix diagonal terms.  Statistics such as',
     - ' C-MAX, C-AVG, C-RMS,'/
     - 'and N should each be multiplied by the average number',
     - ' of DOF per grid point'/
     - 'before estimating NASTRAN time and core requirements.')
      RETURN
      END SUBROUTINE SUMUP

! ##################################################################################################################################
      SUBROUTINE TABLE2(IG,II1,NORIG,IP)
C
C     PRINT CONNECTION TABLE IN TERMS OF ORIGINAL GRID NUMBERS.
C
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IG     ,II1    ,IP     ,J      ,K      ,KMOD   ,
     &         MAXDEG ,MAXGRD ,MM     ,NBITIN ,NN     ,NORIG

      REAL     DUM    ,DUMBB

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,MM,DUM(7)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DIMENSION IG(II1,*),IP(*),NORIG(*)
      IF((NN*MM).EQ.0) RETURN
      WRITE(IOU6,10)
   10 FORMAT(/,'Original Grid Number Connection Table:'/
     -       ' Label  Grid ID     Connections')
      DO 30 I=1,NN
      DO J=1,MM
         IP(J)=0
      end do
      DO 20 J=1,MM
      K=IG(I,J)
CPACK K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(K.EQ.0) GO TO 30
   20 IP(J)=NORIG(K)
   30 WRITE(IOU6,40) I,NORIG(I),(IP(J),J=1,MM)
   40 FORMAT(I6,14I9/(15X,13I9))
      RETURN
      END SUBROUTINE TABLE2

! ##################################################################################################################################
      SUBROUTINE TAXI(KA,IOU6)
C
C     CHECK FOR AXIC BULK DATA CARD, AND, IF PRESENT, READ NUMBER OF
C     HARMONICS.  NEGATIVE NUMBER OF HARMONICS IMPLIES SINGLE HARMONIC
C     PROBLEM.  IT IS ASSUMED THAT THE MNEUMONIC IS LEFT-ADJUSTED.
C
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NAXIC we will use NAXIC_array
! instead. Needed so Lahey doesn't complain about shape of NAXIC being
! different than IP (array) in subr READIT
      INTEGER NAXIC_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(*),A,XIC,BLANK,MINUS
! B////////////////////////////////////////////////////////////////////B
!xx   DATA A,XIC,BLANK,MINUS/'A','XIC',' ','-'/
      DATA A,XIC,BLANK,MINUS/1HA,3HXIC,1H ,1H-/
! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ DUM(8),NAXIC
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IOU6   ,NAXIC  ,NIP

      REAL     DUM

! E////////////////////////////////////////////////////////////////////E
C
C     CHECK FOR AXIC
C
      IF(KA(1).NE.A  ) RETURN
      IF(KA(2).NE.XIC) RETURN
C
C     LOOK FOR FIRST NON-BLANK IN FIELD 2 OF AXIC CARD
C
      DO 10 I=4,11
      IF(KA(I).NE.BLANK) GO TO 20
10    CONTINUE
      RETURN
C
C     READ AXIC INTEGER AND WRITE MESSAGE
C
! B////////////////////////////////////////////////////////////////////B
20    CALL READIT(KA(I),1,12-I,NAXIC_array(1),NIP)
      NAXIC = NAXIC_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(KA(I).EQ.MINUS) NAXIC=-NAXIC
      WRITE(IOU6,30) NAXIC
30    FORMAT(/,'AXIC',I8,' card read.')
C
      RETURN
      END SUBROUTINE TAXI

! ##################################################################################################################################
      SUBROUTINE TIGER(IG,II1,LIST,NORIG,KG,MAXI,IER)
C
C     THIS ROUTINE MAKES ADDITIONS TO THE CONNECTION TABLE IG TO REFLECT
C     THE PRESENCE OF MPC*S AND STORES THE DEPENDENT POINTS IN LIST.
C
C     NEQ = NUMBER OF MPC EQUATIONS.
C
      INTEGER II1, MAXI
      DIMENSION IG(II1,*),LIST(*),NORIG(*)
      INTEGER KG(MAXI)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,
     &         IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,
     &         IOU19  ,IOU20

      INTEGER  I      ,IER    ,IG     ,IGRID  ,II     ,IOPT   ,
     &         J      ,KMOD   ,L      ,LIST   ,
     &         MAXDEG ,MAXGRD ,MM     ,
     &         NBITIN ,NEQ    ,NN     ,NORIG  ,NTERM

      REAL     DUM    ,DUMBB  ,DUMD   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /D/ DUM(7),NEQ,DUMD(2)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,
     -                IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
C     IOPT=1 IF DEPENDENT POINTS ARE TO BE ACCUMULATED IN A LIST FOR
C            LATER REMOVAL FROM THE CONNECTION TABLE IG (OTHERWISE 0)
      DATA IOPT/0/
      IF(NEQ.EQ.0)RETURN
      REWIND IOU12
C     GENERATE NEW CONNECTIONS.
      DO 30 II=1,NEQ
C     READ MPC EQUATION.
      READ(IOU12) NTERM,(KG(I),I=1,NTERM)
C     IGRID=DEPENDENT GRID POINT IN AN MPC EQUATION.
      IGRID=KG(1)
      IF(IOPT.EQ.1) LIST(IGRID)=IGRID
      DO 20 I=1,MAXDEG
         L=IG(IGRID,I)
CPACK    L=IUNPK(IG,MAXGRD*(I-1)+IGRID,NBITIN)
C        L=A GRID POINT THAT IGRID IS CONNECTED TO BEFORE THE MPC IS APPLIED
         IF(L.LE.0) GO TO 30
         IF(NTERM.LT.2) GO TO 20
         DO J=2,NTERM
            CALL SETIG(L,KG(J),IG,II1,NORIG,IER)
            IF(IER.GT.0) RETURN
         end do
20    CONTINUE
   30 CONTINUE
      REWIND IOU12
c     ENDFILE IOU12
c     REWIND IOU12
      RETURN
      END SUBROUTINE TIGER

! ##################################################################################################################################
      SUBROUTINE TIMER(T,IWALL,IOPT,IOU6)
C
C     RETURN IN T THE CURRENT VALUE OF THE CP CLOCK IN SECONDS.
C     RETURN IN IWALL THE ELAPSED TIME IN INTEGER SECONDS.
C     IOPT = 1 IF PRINTOUT OF TIME IS DESIRED, OTHERWISE 0 (INPUT)
C
c     real*4 tarray(2)
c     integer*4 iwall,time
c
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOPT   ,IOU6   ,IWALL

      REAL     T

! E////////////////////////////////////////////////////////////////////E
      t=0.
      iwall=0
c-----------------------------------------------------------------------
c     Cray:
c     call second(t)
c-----------------------------------------------------------------------
c     HP workstation:
c     etime on the HP requires the +U77 Fortran compile-line option.
c     The etime function value is the sum of user and system time,
c     where tarray(1)=user time, tarray(2)=system time.
c     To use these calls, activate the next 3 lines and the real*4 and
c     integer*4 lines above.
c     t=etime(tarray)
c     t=tarray(1)
c     iwall=time(0.)
c-----------------------------------------------------------------------
c     PC using MS Fortran Powerstation or Compaq Visual Fortran:
c     t=secnds(0.0)
c-----------------------------------------------------------------------
c     PC using G77 compiler:
c     call second(t)
c-----------------------------------------------------------------------
c     PC using Lahey compiler or G77 compiler:
c     G77's 'cpu_time' is an alias for 'second', so either can be used.
      call cpu_time(t)
c-----------------------------------------------------------------------
      IF(IOPT.GT.0) WRITE(IOU6,10) T
10    FORMAT('CP time =',F12.3,' seconds')
      RETURN
      END SUBROUTINE TIMER

! ##################################################################################################################################
      SUBROUTINE TREE(IROOT,NDSTK,NR,LVL,IWK,NDEG,LVLWTH,LVLBOT,
     -                LVLN,MAXLW,IBORT)
C
C  TREE DROPS A TREE IN NDSTK FROM IROOT
C
C  LVL-         ARRAY INDICATING AVAILABLE NODES IN NDSTK WITH ZERO
C               ENTRIES. TREE ENTERS LEVEL NUMBERS ASSIGNED
C               DURING EXECUTION OF OF THIS PROCEDURE
C  IWK-         ON OUTPUT CONTAINS NODE NUMBERS USED IN TREE
C               ARRANGED BY LEVELS (IWK(LVLN) CONTAINS IROOT
C               AND IWK(LVLBOT+LVLWTH-1) CONTAINS LAST NODE ENTERED)
C  LVLWTH-      ON OUTPUT CONTAINS WIDTH OF LAST LEVEL
C  LVLBOT-      ON OUTPUT CONTAINS INDEX INTO IWK OF FIRST
C               NODE IN LAST LEVEL
C  MAXLW-       ON OUTPUT CONTAINS THE MAXIMUM LEVEL WIDTH
C  LVLN-        ON INPUT THE FIRST AVAILABLE LOCATION IN IWK
C               USUALLY ONE BUT IF IWK IS USED TO STORE PREVIOUS
C               CONNECTED COMPONENTS, LVLN IS NEXT AVAILABLE LOCATION.
C               ON OUTPUT THE TOTAL NUMBER OF LEVELS + 1
C  IBORT-       INPUT PARAM WHICH TRIGGERS EARLY RETURN IF
C               MAXLW BECOMES .GE. IBORT
C
      INTEGER NR
      DIMENSION NDSTK(NR,*),LVL(*),IWK(*),NDEG(*)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IBORT  ,INOW   ,IROOT  ,ITEST  ,ITOP   ,IWK    ,IWKNOW ,
     &         J      ,KMOD   ,
     &         LVL    ,LVLBOT ,LVLN   ,LVLTOP ,LVLWTH ,
     &         MAXDEG ,MAXGRD ,MAXLW  ,
     &         NBITIN ,NDEG   ,NDROW  ,NDSTK

      REAL     DUMBB

! E////////////////////////////////////////////////////////////////////E
      MAXLW=0
      ITOP=LVLN
      INOW=LVLN
      LVLBOT=LVLN
      LVLTOP=LVLN+1
      LVLN=1
      LVL(IROOT)=1
      IWK(ITOP)=IROOT
   30 LVLN=LVLN+1
   35 IWKNOW=IWK(INOW)
      NDROW=NDEG(IWKNOW)
      DO 40 J=1,NDROW
        ITEST=NDSTK(IWKNOW,J)
CPACK   ITEST=IUNPK(NDSTK,MAXGRD*(J-1)+IWKNOW,NBITIN)
        IF(LVL(ITEST).NE.0) GO TO 40
        LVL(ITEST)=LVLN
        ITOP=ITOP+1
        IWK(ITOP)=ITEST
   40 CONTINUE
      INOW=INOW+1
      IF(INOW.LT.LVLTOP) GO TO 35
      LVLWTH=LVLTOP-LVLBOT
      IF(MAXLW.LT.LVLWTH) MAXLW=LVLWTH
      IF(MAXLW.GE.IBORT) RETURN
      IF(ITOP.LT.LVLTOP) RETURN
      LVLBOT=INOW
      LVLTOP=ITOP+1
      GO TO 30
      END SUBROUTINE TREE

! ##################################################################################################################################
      SUBROUTINE WAVEY(IG,II1,ILD,NEW,NC,IC,KACT,MAXB,MAXW,AVERW,SUMW,
     -                 RMS,BRMS)
C
C     COMPUTE WAVEFRONT AND ACTIVE COLUMN DATA - -
C     MAXIMUM WAVEFRONT, AVERAGE WAVEFRONT, SUM OF ROW WAVEFRONTS,
C     SUM OF SQUARES OF ROW WAVEFRONTS, RMS WAVEFRONT, BANDWIDTH,
C     RMS BANDWIDTH, AND MINIMUM NODAL DEGREE.
C     DIAGONAL TERMS ARE INCLUDED.
C
C     IG      = CONNECTION TABLE
C     II1     = ROW DIMENSION OF IG
C     ILD(I)  = NEW LABEL FOR NODE WITH ORIGINAL INTERNAL LABEL I
C     NEW(I)  = INTERNAL LABEL CORRESPONDING TO NEW LABEL I
C               NEW AND ILD ARE INVERSES OF EACH OTHER
C     NC      = COMPONENT ID
C               IF NC.LE.0, USE ALL COMPONENTS.
C     IC(I)   = COMPONENT INDEX FOR ORIGINAL NODE I.
C     KACT(I) = LIST OF ACTIVE COLUMN FLAGS (UPDATED FOR EACH ROW)
C                 = 1 IF COL I IS ACTIVE AT GIVEN ROW
C               KACT IS SCRATCH SPACE(TEMPORARY STORAGE)
C     MAXB    = BANDWIDTH
C     MAXW    = MAXIMUM WAVEFRONT
C     AVERW   = AVERAGE WAVEFRONT
C     SUMW    = SUM OF ROW WAVEFRONTS
C     SUMSQ   = SUM OF SQUARES OF ROW WAVEFRONTS
C     BSUMSQ  = SUM OF SQUARES OF ROW BANDWIDTHS
C     RMS     = RMS WAVEFRONT
C     BRMS    = RMS BANDWIDTH
C     MINDEG  = MINIMUM NODAL DEGREE
C     NN      = NUMBER OF NODES
C     MM      = MAX NODAL DEGREE
C     MAXGRD  = EFFECTIVE ROW DIMENSION OF IG
C     NBITIN  = NUMBER OF BITS PER INTEGER(CDC)
C     INPUT   - IG,II1,ILD,NN,MM,MAXGRD,NBITIN,NC,IC
C     OUTPUT  - NEW,KACT,MAXW,AVERW,SUMW,RMS,MAXB,BRMS,MINDEG
C
      COMMON /S/ NN,MM,DUMS(5),MINDEG,DUMY
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER II1
      DIMENSION IG(II1,*),ILD(*),NEW(*),KACT(*)
      INTEGER IC(*),SUMW
      DOUBLE PRECISION SUMSQ,BSUMSQ
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IB     ,IB1    ,IG     ,ILD    ,IWAVE  ,
     &         J      ,K      ,KACT   ,KMOD   ,KT     ,L      ,
     &         M      ,MAXB   ,MAXDEG ,MAXGRD ,MAXW   ,MINDEG ,MM     ,
     &         NBITIN ,NC     ,NEW    ,NN

      REAL     ANN    ,AVERW  ,BRMS   ,DUMBB  ,DUMS   ,DUMY   ,
     &         RMS

! E////////////////////////////////////////////////////////////////////E
C
C     INITIALIZE WAVEFRONT DATA.
C
      MAXB=0
      MAXW=0
      SUMW=0
      SUMSQ=0.D0
      BSUMSQ=0.D0
      AVERW=0.
      RMS=0.
      BRMS=0.
      MINDEG=MIN(MINDEG,MM)
      IF((NN*MM).LE.0) RETURN
C
C     INITIALIZE NEW, THE INVERSE OF ILD
C
      IF(NC.GT.0) GO TO 8
      DO 5 I=1,NN
      K=ILD(I)
      IF(K.LE.0) GO TO 5
      NEW(K)=I
    5 CONTINUE
    8 CONTINUE
C
C     INITIALIZE ACTIVE COLUMN FLAGS (1 FOR ACTIVE)
C
      DO 10 I=1,NN
   10 KACT(I)=0
C
C     COMPUTE WAVEFRONT DATA.
C
      IWAVE=1
      KT=0
      DO 40 I=1,NN
C     COMPUTE NUMBER OF ACTIVE COLUMNS FOR ROW I
      K=NEW(I)
      IF(NC) 18,18,15
   15 IF(K.LE.0) GO TO 40
      IF(NC-IC(K)) 40,18,40
18    KT=KT+1
      IB=0
      DO 28 J=1,MM
      L=IG(K,J)
CPACK L=IUNPK(IG,MAXGRD*(J-1)+K,NBITIN)
      IF(L.EQ.0) GO TO 30
      M=ILD(L)
      IB=MAX(IB,I-M)
      IF(M.LE.I) GO TO 28
      IF(KACT(M).EQ.1) GO TO 28
      IWAVE=IWAVE+1
      KACT(M)=1
   28 CONTINUE
      GO TO 35
   30 CONTINUE
      MINDEG=MIN(MINDEG,J-1)
   35 CONTINUE
      IB1=IB+1
C     IB1=ROW BANDWIDTH FOR ROW I (DIAGONAL INCLUDED)
      MAXB=MAX(MAXB,IB1)
      IF(KACT(I).EQ.1) IWAVE=IWAVE-1
C   IWAVE=CURRENT NUMBER OF ACTIVE COLUMNS FOR ROW I (DIAGONAL INCLUDED)
      MAXW=MAX(MAXW,IWAVE)
      SUMW=SUMW+IWAVE
      SUMSQ=SUMSQ+FLOAT(IWAVE)*FLOAT(IWAVE)
      BSUMSQ=BSUMSQ+FLOAT(IB1)*FLOAT(IB1)
   40 CONTINUE
C
      ANN=FLOAT(KT)
      AVERW=FLOAT(SUMW)/ANN
      RMS=SQRT(SNGL(SUMSQ)/ANN)
      BRMS=SQRT(SNGL(BSUMSQ)/ANN)
      RETURN
      END SUBROUTINE WAVEY

! ##################################################################################################################################
      SUBROUTINE ZERO(LIST,NL)
C
C     DELETE ZEROS FROM A LIST OF INTEGERS ORIGINALLY OF LENGTH NL.
C     A CORRECTED LENGTH NL IS RETURNED.
C
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,KT     ,NL

! E////////////////////////////////////////////////////////////////////E
      IF(NL.LE.0) RETURN
      KT=0
      DO 10 I=1,NL
      IF(LIST(I).EQ.0) GO TO 10
      KT=KT+1
      LIST(KT)=LIST(I)
   10 CONTINUE
      NL=KT
      RETURN
      END SUBROUTINE ZERO

! ##################################################################################################################################

      END MODULE BANDIT_MODULE
