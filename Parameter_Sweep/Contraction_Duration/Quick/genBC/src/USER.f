c     MUSC 2 Immediate postop
c     Ethan Kung  keo@ucsd.edu

c     Created by Mahdi Esmaily Moghadam 12-01-2010
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c This subroutine initializes the parameters, you need to read the
c comments and specify them carefuly
c--------------------------------------------------------------------
c This is an example for RCR boundary condition with parameters
c Rd, R, and C which are distal, and proximal resistance and capacitor.

      SUBROUTINE INITIALIZE(nTimeStep)
      USE COM
      IMPLICIT NONE
      INTENT(OUT) nTimeStep

      LOGICAl ierr
      INTEGER i, nTimeStep
      REAL(KIND=8), ALLOCATABLE :: tZeroX(:)
c
c********************************************************************
c For instance if pressure in 3D solver is in cgs and here mmHg
c pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
c here is mL/s qConv=1000=1D3. In the case both solver are using same
c unites you can set these two conversion coefficients to 1D0
      pConv = 1.334D3
      qConv = 1D0

c Only when all the surfaces of you model are coupled with NeumannSrfs
c you may set this to .TRUE.
      pCorr = .FALSE.
      qCorr = .FALSE.

c********************************************************************
c Block for your inputs

!      nbr = 24    !SET NUMBER OF PULMONARY BRANCHES#####
!      ninf = 4    !SET NUMBER OF INFERIOR BRANCHES#####
!      nX=19+nbr+ninf+2 !num of variables (19+nbr+12) ########

c These two value should match to that of defined in solver.inp
      nDirichletSrfs = 1
      nNeumannSrfs   = 1
c Number of unknowns that you need inside your lumped parameter network
      nUnknowns      = 49
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 100

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 20

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      ALLOCATE (tZeroX(nUnknowns), srfToXdPtr(nDirichletSrfs))       !
      ALLOCATE (srfToXPtr(nNeumannSrfs))                             !
      tZeroX = 0D0
c--------------------------------------------------------------------

      INCLUDE "initial_values_final.f"

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      INQUIRE (FILE='InitialData', EXIST=ierr)                       !
      IF (.NOT.ierr) THEN                                            !
c         PRINT *, 'Initializing unknowns in LPM'                     !
         OPEN (1, FILE='InitialData',STATUS='NEW',FORM='UNFORMATTED')!
         WRITE (1) 0D0                                               !
         DO i=1, nUnknowns                                           !
            WRITE (1) tZeroX(i)                                      !
         END DO                                                      !
         CLOSE(1)                                                    !
      END IF                                                         !
c--------------------------------------------------------------------

c Surface to X pointer: this defines which Unknown corresponds to which
c suface in "List of Neumann Surfaces" inside solver.inp
c For example, if you have "List of Neumann Surfaces= 2 8 4 ...."
c and you set "srfToXPtr = (/5,3,9,.../)"
C this means X(5) corresponds to surface 2, X(3) <-> surface 8,
c and X(9) <-> surface 4
c Also, Q(1) corresponds to the first element in srfToXPtr, in this
c example, Q(1)<=>X(5), Q(2)<=>X(3)
      srfToXPtr  = (/48/)
      srfToXdPtr = (/44/)

      END SUBROUTINE INITIALIZE

c####################################################################
c Here you should find the f_i=dx_i/dt, based on the following parameters:
c  current x_i:                   x(i)
c  Current time:                  t
c  Flowrates from 3D code:        Q(i)
c  Pressure from Dirichlet faces: P(i)

      SUBROUTINE FINDF(t, x, f, Q, P)
      USE COM
      IMPLICIT NONE
      INTENT(IN) t, Q
      INTENT(OUT) f

      INTEGER, PARAMETER :: nbr= 24    !SET NUMBER OF PULMONARY BRANCHES#####
      INTEGER, PARAMETER :: ninf= 4    !SET NUMBER OF INFERIOR BRANCHES#####
c      INTEGER, PARAMETER :: nFaces= ninf + nbr + 1 = 29   !SET NUMBER OF INFERIOR BRANCHES#####

      REAL(KIND=8) t, x(nUnknowns), f(nUnknowns), Q(nNeumannSrfs),
     2   P(nDirichletSrfs)

!     These are the dumy variables
      REAL(KIND=8)  Qao, Psa, Psv, Tc, Tm, Tsad, Tr,
     2    Tmr, Tsr, Aa, Av, Ar, dAr, Pith, Piab, dPith, dPiab,
     3    Rabivc, Rthivc, Tve, Pbar, elas, elase, exR,
     4    Aacon, Tmcon, Tcond, Tcon, Qsvc, Qinf, Qpul(nbr),
     5    allC(19), Vol, CdC

      COMPLEX(KIND=8) fdP(20), fP(20), fe(20), fee(30)
      INTEGER k

      INCLUDE "parameters_final.f"

!     Time periods
      Tc   = 60.0/fc
      Tm   = MOD(t,Tc)
      Tcon   = 6D1/fcon
      Tmcon   = MOD(t,Tcon)
      Tsad = Tc - Tsas
      Tcond = Tcon - Tcons

!     Atrium
      IF (Tm .LE. T1) THEN
         Aa = 5D-1*(1D0 - COS(2D0*pi*(Tm - T1 + Tsas)/Tsas))
      ELSEIF (T1+Tsad.LE.Tm .AND. Tm.LE.Tc) THEN
         Aa = 5D-1*(1D0 - COS(2D0*pi*(Tm - T1 - Tsad)/Tsas))
      ELSE
         Aa = 0D0
      END IF

      Psa = Aa*(x(14) - Vsa0)/CCsa + csa*(EXP(dsa*(x(14) - Vsa0)) - 1D0)

! Ventricle.  normalized elastance, these functions peak at 0.3/T

!  "Resting" elastance shape
      fe=(/  (0.223791662, 0.000000000), 
     &       (0.041053540 , -0.409484599), 
     &       (-0.231401116 , -0.026813529), 
     &       (0.017514609 , 0.085190160), 
     &       (0.013159473 , -0.050109864), 
     &       (-0.047293329 , 0.000210483), 
     &       (-0.001339371 , 0.022154522),
     &       (-0.000039917 , -0.005933291), 
     &       (-0.012594069 , 0.002455735), 
     &       (-0.001121362 , 0.009165310),
     &       (0.001600754 , 0.000559328), 
     &       (-0.002965519 , 0.001610140), 
     &       (0.000445087 , 0.003430950), 
     &       (0.000958265 , 0.000409491), 
     &       (0.000004164 , 0.000557564), 
     &       (0.000041503 , 0.000217449), 
     &       (-0.000014608 , 0.000330808), 
     &       (-0.000094191 , 0.000170295), 
     &       (-0.000035561 , 0.000198275), 
     &       (-0.000039917 , 0.000147014) /) 

      elas=0D0   
      DO k = 1, 20, 1
         elas=elas + real(fe(k))*cos(2D0*pi*(k-1)* Tm*.3D0/Tsvs ) 
     2   - aimag(fe(k))*sin(2D0*pi*(k-1)*  Tm*.3D0/Tsvs )  !scale by Tsvs
      END DO
! set elas=0 if Tc/Tsvs is too large when HR is low, and go over the elas waveform cycle      
      IF (Tm/Tsvs .GT. 3.33333D0) THEN  !meaning that Tm*.3D0/Tsvs > 1
         elas = 0D0
      END IF

      elas= a*elas + Eoffset   !shift by Emin
      Psv = elas*(x(16) - Vsv0)   

!     This is the only case that x is directly manipulated, ugly, but the
!     only way!!!
      IF (x(15) .LT. 0D0) x(15) = 0D0 
      IF (x(44) .LT. 0D0) x(44) = 0D0 

!     Respiration
      Tr = 6D1/fr
      Tmr = MOD(t,Tr)
      Tsr = Trr*Tr  !don't need this if prescribing fourier waveform

      IF (Tmr .LT. Tsr) THEN
         Ar  = 5D-1*(1D0 - COS(2D0*pi*Tmr/Tsr))
         dAr = (pi/Tsr)*SIN(2D0*pi*Tmr/Tsr)
      ELSE
         Ar  = 0D0
         dAr = 0D0
      END IF
      Pith  = APith*Ar + P0ith
      dPith = APith*dAr

      Piab  = APiab*Ar + P0iab
      dPiab = APiab*dAr
      

!!!!!! Now connect Pith to the heart pressures !!!
!      Psv = Psv + Pith
!      Psa = Psa + Pith


!     Aortic flow
      IF ( Psv .GT. x(4)) THEN
         Qao = (SQRT(Rmyo*Rmyo + 4D0*Kao*(Psv 
     2      - x(4))) - Rmyo)/2D0/Kao
      ELSE
         Qao = 0D0
      END IF

      !     No more adjusting parameters based on collapsibility
      Rabivc = R0abivc  
      Rthivc = R0thivc 

!     Calculate Fontan junction flows
!      Qpul = (x(20+nbr+ninf+1)-x(20:19+nbr))/Rr(2+ninf:1+ninf+nbr)
      Qpul = (x(20+nbr+ninf)-x(20:19+nbr))/(Rd/10)
!       Qpul = Qsvc + x(20+nbr+ninf+1)
!ivc flow
!      Qinf = (x(20+nbr:19+nbr+ninf)-x(20+nbr+ninf)) / Rr(2:1+ninf)

!svc flow
      Qsvc = (x(1)-x(20+nbr+ninf)) / Rr(1)   

!     The main body of equations
!     upper body
!      f(1)  = ((x(2) - x(1))/Rubv - Qsvc)/Csvc + dPith
      f(1)  = ((x(2) - x(1))/Rubv - Qsvc)/Csvc
      f(2)  = (x(3) - (x(2) - x(1))/Rubv)/Cub
      f(3)  = (x(4) - x(3)*Ruba - x(2))/Luba

!     aorta
!      f(4) = (Qao - x(3) - x(5))/Cao + dPith
      f(4) = (Qao - x(3) - x(5))/Cao

      f(5) = (x(4) - x(5)*Rthao - x(6))/Lthao
!      f(6) = (x(5) - x(7) - (x(6) - x(17))/Rla
!     2  - (x(6) - x(18))/Rka)/Cthao + dPith
      f(6) = (x(5) - x(7) - (x(6) - x(17))/Rla
     2  - (x(6) - x(18))/Rka)/Cthao
      f(7) = (x(6) - x(7)*Rabao - x(8))/Labao
!      f(8) = (x(7) - (x(8) - x(19))/Ria - x(9))/Cabao + dPiab
      f(8) = (x(7) - (x(8) - x(19))/Ria - x(9))/Cabao 

!     legs
      f(9) = (x(8) - x(9)*Rlega - x(10))/Llega
      f(10) = (x(9) - (x(10) - x(11))/Rlegc)/Clega
      f(11) = ((x(10) - x(11))/Rlegc - (abs(x(11) - x(12))
     2  + x(11) - x(12))/Rlegv/2D0)/Clegv

!     inferior vena cava
!      f(12) = ((abs(x(11) - x(12)) + x(11) - x(12))/Rlegv/2D0
!     2  - (x(12) - x(13))/Rabivc)/Cabivc + dPiab
!     f(13) and f(20) is modified for stage 3
!      f(13) = ((x(12) - x(13))/Rabivc + (x(18) - x(13))/Rkv
!     2  + (x(17) - x(13))/Rlv - x(44))/Cthivc + dPith
      f(12) = ((abs(x(11) - x(12)) + x(11) - x(12))/Rlegv/2D0
     2  - (x(12) - x(13))/Rabivc)/Cabivc
!     f(13) and f(20) is modified for stage 3
      f(13) = ((x(12) - x(13))/Rabivc + (x(18) - x(13))/Rkv
     2  + (x(17) - x(13))/Rlv - x(44))/Cthivc

!      ventriclce
      f(14) = SUM((x(20:19+nbr) - Psa )/Rd) - x(15)      
      f(15) = (Psa - Psv - Kav*x(15)*x(15))/Lav
      f(16) = x(15) - Qao

!     liver
!      f(17) = ((x(19) - x(17))/Riv + (x(6) - x(17))/Rla
!     2  - (x(17) - x(13))/Rlv)/Cl + dPiab
      f(17) = ((x(19) - x(17))/Riv + (x(6) - x(17))/Rla
     2  - (x(17) - x(13))/Rlv)/Cl

!     kidney
!      f(18) = ((x(6) - x(18))/Rka - (x(18) - x(13))/Rkv)/Ck + dPiab
      f(18) = ((x(6) - x(18))/Rka - (x(18) - x(13))/Rkv)/Ck

!     intestine
!      f(19) = ((x(8) - x(19))/Ria - (x(19) - x(17))/Riv)/Ci + dPiab
      f(19) = ((x(8) - x(19))/Ria - (x(19) - x(17))/Riv)/Ci

      IF (x(15).EQ.0D0 .AND. Psa.LT.(Psv-Rmyo*Qao)) THEN
         f(15) = 0D0
      END IF

! pulmonary branches  20:43 
!      f(20:19+nbr)= (Qpul(1:nbr)- (x(20:19+nbr)-Psa)/Rd)/Ca  +dPith
      f(20:19+nbr)= (Qpul(1:nbr)- (x(20:19+nbr)-Psa)/Rd)/Ca

! inferior branches 44:47
      f(20+nbr:19+nbr+ninf) = f(13)

! Conjuction pressure 48
      f(20+nbr+ninf) = (Q(1)+Qsvc-SUM(Qpul(1:nbr)))/Ccon

! Upper conduit flow 49
!      f(20+nbr+ninf+1) = (P(1)-x(20+nbr+ninf+1)*x(20+nbr+ninf+1)*Kuc-
!     2   x(20+nbr+ninf))/Luc

!      IF (x(20+nbr+ninf+1) .LT. 0D0 .AND. P(1) .LT. x(20+nbr+ninf)) THEN
!         f(20+nbr+ninf+1) = 0D0
!      END IF
       f(20+nbr+ninf+1) = 0

! Lower conduit flow 44
      f(44) = (x(13)-x(44)*Rinf-P(1)-x(44)*x(44)*Klc)/Llc

      IF (x(44) .LE. 0D0 .AND. (x(13)-x(44)*Rinf).LT. P(1)) THEN
         f(44) = 0D0
      END IF

      allC = (/Csvc, Cub, 0D0, Cao, 0D0, Cthao, 0D0, Cabao, 0D0, 
     2   Clega, Clegv, Cabivc, Cthivc, 1D0, 0D0, 1D0, Cl, Ck, Ci/)
!      Vol = SUM(allC*x(1:32)) - 2D0*Cc*Plvp + SUM(Cpa1*x(33:32+nbr)) +
!     2   SUM(Cpa2*x(33+nbr:32+2*nbr))
      Vol = SUM(allC*x(1:19)) + Ccon*x(48) + SUM(Ca*x(20:19+nbr))
!      CdC = SUM(allC*allC)
      CdC = SUM(allC*allC)+ SUM(Ca*Ca) + Ccon*Ccon

      offset(1:19) = (Vtot - Vol)/CdC*allC
      offset(20:43) = (Vtot - Vol)/CdC*Ca
      offset(48) = (Vtot - Vol)/CdC*Ccon


!     Assign the additional parameters to be printed
      Xprint(1) = t
      Xprint(2) = elas !Qivc 
      Xprint(3) = Psa
      Xprint(4) = Psv
      Xprint(5) = x(14)  !Vsa
      Xprint(6) = Qao
      Xprint(7) = Pith
      Xprint(8) = dPith
      Xprint(9) = f(20+nbr+ninf)
      Xprint(10) = f(20+nbr+ninf+1)
      Xprint(11) = 1
      Xprint(12) = SUM((x(20:19+nbr) -  Psa)/Rd) !61
      Xprint(13) = P(1) !62
      Xprint(14) = Q(1) !63
      Xprint(15) = Qsvc !64
      Xprint(16) = x(13)-x(44)*Rinf !65
      Xprint(17) = f(13)!66
      Xprint(18) = x(13) !67
      Xprint(19) = f(44) !68
      Xprint(20) = Vol !69



      RETURN
      END SUBROUTINE FINDF


