*##################################################################
*     FUNCTIONS
*##################################################################
c +-------------------------------------------------------------------+
c     PROGRAM  : FERMI
c     TYPE     : function
c     PURPOSE  : calculate the Fermi-Dirac distribution
c     VERSION  : 31-01-2006
c     COMMENTS : 
c +-------------------------------------------------------------------+
      FUNCTION fermi(x, mu1, beta1)
      double precision fermi, x, mu1, beta1, xeff
      xeff=x-mu1
      if(abs(xeff).le.0.5d0)then
         fermi = 1.d0/(1.d0+exp(beta1*(xeff)))      
      else
         if(xeff.gt.0.d0)fermi=0.d0
         if(xeff.lt.0.d0)fermi=1.d0
      endif
      return
      END

c +-------------------------------------------------------------------+
c     PROGRAM  : DFERMI
c     TYPE     : function
c     PURPOSE  : calculate the Fermi-Dirac distribution
c     VERSION  : 31-01-2006
c     COMMENTS : 
c +-------------------------------------------------------------------+
      FUNCTION dfermi(x, mu1, beta1)
      REAL*8 dfermi, x, mu1, beta1,zeta1
      zeta1=exp(-beta1*mu1)
      if(dabs(x).gt.3.d0)then   !take care with zeta in this case
         dfermi=0.d0
      else
         dfermi=-(zeta1*beta1*exp(beta1*x))/(1.d0+zeta1*exp(beta1*x))**2
      endif
      return
      END

c     +-------------------------------------------------------------------+
c     PROGRAM  : DENS
c     TYPE     : function
c     PURPOSE  : calculate the non-interacting density of states
c     VERSION  : 31-01-2006
c     COMMENTS : we are using the hypercubic lattice that has a 
c     gaussian dos
c     +-------------------------------------------------------------------+
      FUNCTION ddens(x,if1)
      integer if1
      REAL*8 ddens,x,t1,D
      complex*16 root
      if(if1.eq.2)then
         t1=1.d0/sqrt(2.d0)
         ddens = (1/(t1*sqrt(2*3.141592653589793238d0)))*exp(-(x**2)/
     &        (2*t1**2))
      elseif(if1.eq.1)then
         D=(1.d0,0.d0)      
         root=dcmplx((1.d0-1.d0*((x/D))**2),0.d0)
         root=sqrt(root)
         ddens=(2.d0/(3.141592653589793238d0*D))*root
      endif
      return
      END


*##################################################################
*     SUBROUTINES
*##################################################################
c########################################################################
c     PROGRAM  : FOUR1
c     TYPE     : subroutine
c     PURPOSE  : Fourier(anti)-transform a given function
c     COMMENTS : the M.frequencies are adapted to this subroutine.
c########################################################################
      SUBROUTINE four1(data,nn,isign)
c     implicit real*8(a-h,o-z)
c     DOUBLE PRECISION data(2*nn)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
         if(j.gt.i)then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
 1       if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            goto 1
         endif
         j=j+m
 11   continue
      mmax=2
 2    if (n.gt.mmax) then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*dsin(0.5d0*theta)**2
         wpi=dsin(theta)
         wr=1.d0
         wi=0.d0
         do 13 m=1,mmax,2
            do 12 i=m,n,istep
               j=i+mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
 12         continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
 13      continue
         mmax=istep
         goto 2
      endif
      return
      END

c########################################################################
c     PROGRAM  : INTERP
c     TYPE     : subroutine
c     PURPOSE  : interpolate a G^qmc(-Lfak:Lfak) to 
c     G^smooth(-L:L) using the cubic-spline algorithm. 
c     COMMENTS : 
c########################################################################
      subroutine interp(gtmp,g,Lfak,Lfak1,L)
      implicit real*8(a-h,o-z)
      integer L,Lfak,Lfak1
      double precision gtmp(-Lfak:Lfak),g(-L:L)
      double precision xa(Lfak1),ya(4,Lfak1)
c************************************************************************
c     interpolate gtmp which has Lfak pts to g which has L pts
c************************************************************************
      do 10 i=1,Lfak1
         xa(i)=float(i-1)/float(Lfak)
         ya(1,i) = gtmp(i-1)
 10   continue
      call CUBSPL(xa,ya,Lfak1,0,0)
c*******assign g(i)********
      do 20 i=1,L
         x=float(i)/float(L)
         g(i) = PPVALU(xa,ya,Lfak,4,x,0)
 20   continue
      g(0)=gtmp(0)
      do 40 i=1,L
         g(-i)=-g(L-i)
 40   continue
      return            
      end

c+--------------------------------------------------------------------
c  Extracted from "A practical guide to splines," 1st edition,       +
c  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   -
c--------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION PPVALU (BREAK, COEF, L, K, X, JDERIV)
      IMPLICIT NONE
CALLS INTERV
CALCULATES VALUE AT X OF JDERIV-TH DERIVATIVE OF PP FCT FROM PP-REPR
C
C******  I N P U T  *********
C  BREAK, COEF, L, K.....FORMS THE PP-REPRESENTATION OF THE FUNCTION F
C        TO BE EVALUATED. SPECIFICALLY, THE J-TH DERIVATIVE OF F IS 
C        GIVEN BY
C
C     (D**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I) +
C                             + H*COEF(K,I)/K-J-I))/(K-J-2) ... )/2)/1
C
C        WITH  H = X - BREAK(I),  AND
C
C       I = MAX( 1 , MAX( J , BREAK(J) .LE. X , 1 .LE. J .LE. L ) ).
C
C  X.....THE POINT AT WHICH TO EVALUATE.
C  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUAT-
C        ED. ASSUMED TO BE ZERO OR POSITIVE.
C
C*****  O U T P U T  *********
C  PPVALU.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF F AT X.
C
C*****  M E T H O D  *********         
C    THE INTERVAL INDEX I, APPROPRIATE FOR X, IS FOUND THROUGHT A 
C  CALL TO INTERV. THE FORMULA ABOVE FOR THE JDERIV-YH DERIVATIVE
C  OF F IS THEN EVALUATED (BY NESTED MULTIPLICATION).
C
      INTEGER JDERIV, K, L, I, M, NDUMMY
      DOUBLE PRECISION BREAK(L), COEF(K,L), X, FMMJDR, H
      PPVALU = 0.D0
      FMMJDR = K - JDERIV
C              DERIVATIVES OF ORDER K OR HIGHER ARE IDENTICALLY ZERO.
      IF (FMMJDR .LE. 0.D0)           GO TO 99
C
C              FIND INDEX I OF LARGEST BREAKPOINT TO THE LEFT OF X.
      CALL INTERV ( BREAK, L, X, I, NDUMMY )
C
C      EVALUATE JDERIV-TH DERIVATIVE OF I-TH POLYNOMIAL PIECE AT X.
      H = X - BREAK(I) 
      DO 10 M=K,JDERIV+1,-1
         PPVALU = (PPVALU/FMMJDR)*H + COEF(M,I)
  10     FMMJDR = FMMJDR - 1.D0
  99                                    RETURN
      END


c+--------------------------------------------------------------------
c  Extracted from "A practical guide to splines," 1st edition,       +
c  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   -
c--------------------------------------------------------------------+  

      SUBROUTINE CUBSPL ( TAU, C, N, IBCBEG, IBCEND)
C     *****************************************************************
C     N = NUMBER OF DATA POINTS. ASSUMED TO BE .GE. 2.
C     (TAU(I), C(1,I), I=1,...,N) = ABSCISSAE AND ORDINATES OF THE
C        DATA POINTS. TAU IS ASSUMED TO BE STRICTLY INCREASING.
C     IBCBEG, IBCEND = BOUNDARY CONDITION INDICATORS, AND
C     C(2,1), C(2,N) = BOUNDARY CONDITION INFORMATION. SPECIFICALLY,
C        IBCBEG = 0  MEANS NO BOUNDARY CONDITION AT TAU(1) IS GIVEN.
C           IN THIS CASE, THE NOT-A-KNOT CONDITION IS USED, I.E. THE
C           JUMP IN THE THIRD DERIVATIVE ACROSS TAU(2) IS FORCED TO
C           ZERO, THUS THE FIRST AND THE SECOND CUBIC POLYNOMIAL PIECES
C           ARE MADE TO COINCIDE.)
C        IBCBEG = 1  MEANS THAT THE SLOPE AT TAU(1) IS MADE TO EQUAL
C           C(2,1), SUPPLIED BY INPUT.
C        IBCBEG = 2  MEANS THAT YHE SECOND DERIVATIVE TAU(1) IS 
C           MADE TO EQUAL C(2,1), SUPPLIED BY INPUT.
C        IBCEND = 0, 1, OR 2 HAS ANALOGOUS MEANING CONCERNING THE
C           BOUNDARY CONDITION AT TAU(N), WITH THE ADDITIONAL INFOR-
C           MATION TAKEN FROM C(2,N).
C     *********************** 	OUTPUT  ******************************
C     C(J,I), J=1,...,4; I=1,...,L (= N-1) = THE POLYNOMIAL COEFFICIENTS 
C        OF THE CUBIC INTERPOLATING SPLINE WITH INTERIOR KNOTS (OR
C        JOINTS) TAU(2),...,TAU(N-1). PRECISELY, IN THE 
C        INTERVAL (TAU(I), TAU(I+1)), THE SPLINE F IS GIVEN BY
C           F(X) = C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)/3.)/2.)
C        WHERE H = X - TAU(I). THE FUNCTION PROGRAM *PPVALU* MAY BE 
C        USED TO EVALUATE F OR ITS DERIVATIVES FROM TAU, C, L = N-1,
C        AND K=4.
      IMPLICIT NONE
      INTEGER IBCBEG, IBCEND, N, I, J, L, M
      DOUBLE PRECISION C(4,N), TAU(N), DIVDF1, DIVDF3, DTAU, G
C******  A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SPLOPES S(I) OF
C     F AT TAU(I), I=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS 
C     ELIMINATION, WITH S(I) ENDING UP IN C(2,I), ALL I.
C     C(3,.) AND C(4,.) ARE USED INITIALLY FOR TEMPORARY STORAGE.
      L = N - 1
COMPUTE FIRST DIFFERENCES OF TAU SEQUENCE AND STORE IN C(3,.). ALSO,
COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN C(4,.).
      DO 10 M = 2, N
         C(3,M) = TAU(M) - TAU(M-1)
   10    C(4,M) = (C(1,M) - C(1,M-1))/C(3,M)
CONSTRUCT FIRST EQUATION FROM THE BOUNDARY CONDITION, OF THE FORM
C             C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
      IF (IBCBEG-1)                     11,15,16
   11 IF (N .GT. 2)                     GO TO 12
C     NO CONDITION AT LEFT END AND N = 2.
      C(4,1) = 1.D0
      C(3,1) = 1.D0
      C(2,1) = 2.D0*C(4,2)
                                        GO TO 25
C     NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
   12 C(4,1) = C(3,3)
      C(3,1) = C(3,2) + C(3,3)
      C(2,1) = 
     *      ((C(3,2)+2.D0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))/C(3,1)
                                        GO TO 19
C     SLOPE PRESCRIBED AT LEFT END.
   15 C(4,1) = 1.D0
      C(3,1) = 0.D0
                                        GO TO 18
C     SECOND DERIVATIVE PRESCRIBED AT LEFT END.
   16 C(4,1) = 2.D0
      C(3,1) = 1.D0
      C(2,1) = 3.D0*C(4,2) - C(3,2)/2.D0*C(2,1)
   18 IF (N .EQ. 2)                     GO TO 25
C  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESP. EQUATIONS AND CAR-
C  RY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE M-TH
C  EQUATION READS    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).
   19 DO 20 M=2,L
         G = -C(3,M+1)/C(4,M-1)
         C(2,M) = G*C(2,M-1) + 3.D0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
   20    C(4,M) = G*C(3,M-1) + 2.D0*(C(3,M) + C(3,M+1))
CONSTRUCT LAST EQUATION FROM THE SECOND BOUNDARY CONDITION, OF THE FORM
C           (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
C     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
C     SUBSTITUTION, SINCE C ARRAY HAPPENS TO BE SET UP JUST RIGHT FOR IT
C     AT THIS POINT.
      IF (IBCEND-1)                     21,30,24
   21 IF (N .EQ. 3 .AND. IBCEND .EQ. 0) GO TO 22
C     NOT-A-KNOT AND N .GE. 3, AND EITHER N .GT. 3 OR ALSO NOT-A-KNOT AT
C     LEFT END POINT.
      G = C(3,N-1) + C(3,N)
      C(2,N) = ((C(3,N)+2.D0*G)*C(4,N)*C(3,N-1)
     *           + C(3,N)**2*(C(1,N-1)-C(1,N-2))/C(3,N-1))/G
      G = -G/C(4,N-1)
      C(4,N) = C(3,N-1)
                                        GO TO 29
C     EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND NOT-A-
C     KNOT AT LEFT END POINT).
   22 C(2,N) = 2.D0*C(4,N)
      C(4,N) = 1.D0
                                        GO TO 28
C     SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
   24 C(2,N) = 3.D0*C(4,N) + C(3,N)/2.D0*C(2,N)
      C(4,N) = 2.D0
                                        GO TO 28
   25 IF (IBCEND-1)                     26,30,24
   26 IF (IBCBEG .GT. 0)                GO TO 22
C     NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
      C(2,N) = C(4,N)
                                        GO TO 30
   28 G = -1.D0/C(4,N-1)
COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
   29 C(4,N) = G*C(3,N-1) + C(4,N)
      C(2,N) = (G*C(2,N-1) + C(2,N))/C(4,N)
CARRY OUT BACK SUBSTITUTION
   30 DO 40 J=L,1,-1
   40    C(2,J) = (C(2,J) - C(3,J)*C(2,J+1))/C(4,J)
C****** GENERATE CUBIC COEFFICIENTS IN EACH INTERVAL, I.E., THE DERIV.S
C  AT ITS LEFT ENDPOINT, FROM VALUE AND SLOPE AT ITS ENDPOINTS.
      DO 50 I=2,N
         DTAU = C(3,I)
         DIVDF1 = (C(1,I) - C(1,I-1))/DTAU
         DIVDF3 = C(2,I-1) + C(2,I) - 2.D0*DIVDF1
         C(3,I-1) = 2.D0*(DIVDF1 - C(2,I-1) - DIVDF3)/DTAU
   50    C(4,I-1) = (DIVDF3/DTAU)*(6.D0/DTAU)
                                          RETURN
      END


c+--------------------------------------------------------------------
c  Extracted from "A practical guide to splines," 1st edition,       +
c  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   -
c--------------------------------------------------------------------+

      SUBROUTINE INTERV ( XT, LXT, X, LEFT, MFLAG )
COMPUTES LEFT = MAX( I , 1 .LE. I .LE. LXT .AND. XT(I) .LE. X ).
C
C****** I N P U T  ***************
C  XT.....A REAL SEQUENCE, OF LENGHT LXT, ASSUMED TO BE NONDECREASING
C  LXT.....NUMBER OF TERMS IN THE SEQUENCE XT.
C  X.....THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE XT IS 
C        TO BE DETERMINED.
C
C****** O U T P U T  *************
C  LEFT, MFLAG.....BOTH INTEGERS, WHOSE VALUE IS
C
C   1     -1      IF               X .LT. XT(1)
C   I      0      IF   XT(I)  .LE. X .LT. XT(I+1)
C  LXT     1      IF  XT(LXT) .LE. X
C
C         IN PARTICULAR, MFLAG = 0 IS THE 'USUAL' CASE. MFLAG .NE. 0
C         INDICATES THAT X LIES OUTSIDE THE HALFOPEN INTERVAL
C         XT(1) .LE. Y .LT. XT(LXT). THE ASYMETRIC TREATMENT OF THE 
C         INTERVAL IS DUE TO THE DECISION TO MAKE ALL PP FUNCTIONS CONT-
C         INUOUS FROM THE RIGHT.
C
C*****  M E T H O D  ********
C  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT
C  IT IS CALLED REPEATEDLY, WITH X TAKEN FROM AN INCREASING OR DECREA-
C  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE
C  GRAPHED. THE FIRST GUESS FOR LEFT IS THEREFORE TAKEN TO BE THE VAL-
C  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE L O C A L VARIA-
C  BLE  ILO . A FIRST CHECK ASCERTAINS THAT  ILO .LT. LXT (THIS IS NEC-
C  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI- 
C  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT = 
C  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.
C     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO
C  WHILE ALSO MOVING  ILO  AND  IHI  IN THE DIRECTION OF  X, UNTIL
C                      XT(ILO) .LE. X .LT. XT(IHI) ,
C  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .
C  LEFT = ILO  IS THE RETURNED.
C
      IMPLICIT NONE
      INTEGER  LEFT, LXT,MFLAG, IHI, ILO, ISTEP, MIDDLE
      DOUBLE PRECISION X, XT(LXT)
      DATA ILO /1/
C     SAVE ILO  (A VALID FORTRAN STATEMENT IN THE NEW 1977 STANDARD)
      IHI = ILO + 1
      IF (IHI .LT. LXT)                 GO TO 20
         IF (X .GE. XT(LXT))            GO TO 110
         IF (LXT .LE. 1)                GO TO 90
         ILO = LXT - 1
         IHI = LXT
C
  20  IF (X .GE. XT(IHI))               GO TO 40
      IF (X .GE. XT(ILO))               GO TO 100
C
C             ****  NOW X .LT. XT(ILO). DECREASE  ILO  TO CAPTURE X .
      ISTEP = 1
  31     IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .LE. 1)                GO TO 35
         IF (X .GE. XT(ILO))            GO TO 50
         ISTEP = ISTEP*2
                                        GO TO 31
  35  ILO = 1
      IF (X .LT. XT(1))                 GO TO 90
                                        GO TO 50
C              **** NOW X .GE. XT(IHI). INCREASE  IHI  TO CAPTURE  X . 
  40  ISTEP = 1
  41     ILO = IHI
         IHI = ILO + ISTEP
         IF (IHI .GE. LXT)              GO TO 45
         IF (X .LT. XT(IHI))            GO TO 50
         ISTEP = ISTEP*2
                                        GO TO 41
  45  IF (X .GE. XT(LXT))               GO TO 110
      IHI = LXT
C
C           **** NOW  XT(ILO) .LE. X .LT. XT(IHI). NARROW THE INTERVAL.
  50  MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO)              GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1 .
      IF (X .LT. XT(MIDDLE))            GO TO 53
         ILO = MIDDLE
                                        GO TO 50
  53     IHI = MIDDLE
                                        GO TO 50
C**** SET OUTPUT AND RETURN.
  90  MFLAG = -1
      LEFT = 1
                                        RETURN
 100  MFLAG = 0
      LEFT = ILO
                                        RETURN
 110  MFLAG = 1
      LEFT = LXT
                                        RETURN
      END 
