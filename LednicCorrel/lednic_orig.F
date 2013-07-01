          PROGRAM FSI
C---Calculates the correlation function R(p1,p2) of two particles with
C-- 4-momenta p1, p2 defined as the ratio of the double inclusive
C-- momentum distribution to te one which would be observed in the case
C-- of no quantum statistics and final state interaction effects.
C-- R= <|B-S amplitude|^2>, where the modulus squared of the Bethe-
C-- Salpeter amplitude is averaged over the total spin S of the particle
C-- pair and integrated over the space-time coordinates of the
C-- production points.
C-- Integration in the source rest frame (SRF),
C-- Bethe-Salpeter amplitude calculated in the pair rest frame (PRF).
C---References:
C-- R.Lednicky, V.L.Lyuboshitz, Yad. Fiz. 35 (1982) 1316
C-- M.Gmitro, J.Kvasil, R.Lednicky, V.L.Lyuboshitz,
C--                             Czech. J. Phys. B36 (1986) 1281
C---The approximation used to calculate Bethe-Salpeter amplitude is set
C-- with the key NS:
C-- NS=1  Square well potential,
C-- NS=3  Exponential potential (for two neutrons at p1=p2 only),
C-- NS=4  scattered wave approximated by the spherical wave,
C-- NS=5  approximation of large k* (valid for symmetric spatial
C--       distribution only); here only for the case of sph. symmetric
C--       gaussian (LAW=1) and for spin-0 particles (now only for
C--       alpha-alpha: LL=4)
C-- see V.L.Lyuboshitz, YaF 41 (1985) 820
C--     B.K.Jenning et al., PR C33 (1986) 1303
C-- NS=2  same as NS=4 but the approx. of equal emission times in PRF
C--       not required (t=0 approx. used in all other cases).
C---The gaussian distribution is assumed for the emission time if
C-- LAW=1-5 while various distribution laws of the spatial coordinates  
C-- of the emission points may be chosen with the help of the key LAW:  
C-- LAW=1  cylindrically symmetrical gaussian,
C-- LAW=2  sphere volume,
C-- LAW=3  sphere surface,
C-- LAW=4  sphere surface + Lambert law;
C-- LAW=5  gaussian correlated tube:
C-- particles are produced near the projectile trajectory in the
C-- nuclear matter described by the gaussian density with the width ZRZ;
C-- the gaussian spread of the space-time coord. of the product. points
C-- around the trajectory r = u*t (u=1 is projectile velocity) is ZT=ZR,
C-- see M.L.Allaberdin, R.Lednicky, V.L.Lyuboshitz, B.S.Yuldashev,
C--                             Yad. Fiz. 46 (1987) 1785.
C-- LAW=6  cylindrically symmetrical gaussian for spatial coordinates
C--        like LAW=1, but exponetial decay law is assumed: exp(-!t!/ZT)
C-- parameters ZR, ZRZ, ZT in fm;
C-- if LAW=2-4 then ZR= radius,
C-- if LAW=1,5 then ZR measures transv. gauss. width and ZRZ measures
C--            gauss. width along the reaction axis,
C--            ZRZ=ZR in the case of spherical symmetry;
C-- r.m.s. radius = 2*ZR**2+ZRZ**2 if LAW=1,
C--               = (3/5)*ZR**2    if LAW=2-4;
C-- ZT measures the width of the time distribution: the time difference
C-- t in the SRF is distributed as exp[-(t/ZT)**2/4] if LAW=1-4
C--                          or as exp(-!t!/ZT) if LAW=6
C---Switch:
C-- ISW=1  complete calculation,
C-- ISW=2  the strong interaction is switched off,
C-- ISW=3  in addition the quantum statistics is switched off;
C-- to switch off the Coulomb force put ICH=0 and substitute the strong
C-- amplitude parameters by the ones not affected by Coulomb interact.
C-- (e.g., the scattering length may be calculated by SEAR).
C---Kinematical variables:
C-- The foll. 5 variab. charact. the particle pair in the source rest
C-- frame (azimuth. sym. with resp. to the reaction axis z is ass.):
C--   Q [GEV], CVK=cos(beta), FIQ, V, CTHT=cos(theta),
C-- V and theta are the pair velocity and production angle in SRF,
C-- Q= 2*k is the relative momentum of the two particles in PRF,
C-- FIQ= azimuth. angle of the vector Q (FIQ=0 if Q lies in the pair
C--      production plane),
C-- beta= angle between the vectors V and Q;
C-- the variables Q0, QT and QL= Q0/V in SRF are also used:
C-- QT= Q*sin(beta), QL= Q*cos(beta)*GAM (GAM= pair Lorentz factor).
C---Note:
C-- ICH= 0 (1) if the Coulomb interaction is absent (present);
C-- ISPIN= JJ= 1,2,..,MSPIN denote incr. values of the contributing pair
C-- total spin S.
C---The following two particle systems may be chosen with the key LL:
C-- LL=1  neutron-neutron     NS=1-4
C-- LL=2  proton-proton       NS=1,2,4
C-- LL=3  neutron-proton      NS=1,2,4
C-- LL=4  alpha-alpha         NS=1,2,4
C--       the d-wave contribution is added in the large r* and large k*
C--       approximation, valid for symmetric spatial distribution only
C--       and for nonrelativistic particles;
C--       here only the case of sph. symmetric gaussian (LAW=1);
C--       formulae of V.L.Lyuboshitz, YaF 41 (1985) 820
C--                   B.K.Jenning et al., PR C33 (1986) 1303
C--       are modified to the case of charged particles (it is important
C--       to calculate exactly the Coulomb WF of free particles)
C-- LL=5  pi+pi-              NS=1,2,4
C-- LL=6  pi0pi0              NS=2,4
C-- LL=7  pi+pi+ (pi-pi-)     NS=2,4
C-- LL=8  neutron-deuteron    NS=2,4
C-- LL=9  proton-deuteron     NS=2,4.
C-- To calculate the CF of two particles (with masses m1, m2 and
C-- charges Z1, Z2) the following information is required:
C-- AM= AMH(LL)= 2*m1*m2/(m1+m2) in GeV/c^2,
C-- DM= DMH(LL)= (m1-m2)/(m1+m2), required if NS=2;
C-- AC= ABH(LL)= Bohr radius= 2*137.036*0.1973270536/(Z1*Z2*AMH) in fm;
C-- AC > 1.D9 if Z1*Z2= 0, AC < 0 if Z1*Z2 < 0;
C-- MSPIN= MSPINH(LL)= number of contributing values of the total pair
C-- spin S;
C-- FD= FDH(LL,JJ), RD= RDH(LL,JJ)= scattering length and effective
C---radius for each value of contrib. total pair spin S, JJ= 1,..,MSPIN;;
C-- the corresponding square well parameters EB= EBH(LL,JJ), RB=
C-- RBH(LL,JJ) (required if NS=1) may be calculated by SEAR;
C-- if the effective range approximation is not valid (as is the case,
C-- e.g., for two-pion system) a code for calculation of the scattering
C-- amplitude should be supplemented;
C-- RHO= RHOH(LL,JJ), SF= SFH(LL,JJ), SE= SEH(LL) are spin factors;
C-- RHO= the probability that the spins j1 and j2 of the two particles
C-- will combine in a total spin S;
C-- RHO= (2*S+1)/[(2j1+1)*(2j2+1)] for unpolarized particles;
C-- RHO= (1-P1*P2)/4 and (3+P1*P2)/4 correspond to S=0 and 1 in the
C-- case of spin-1/2 particles with polarizations P1 and P2;
C-- the spin factors SE and SF depend on the identity of the particles;
C-- for identical particles:
C-- SF= 2*RHO if S is even, SF= 0 if S is odd (see Eq. 10 of L-L);
C-- SE is the sum of the RHO values corresponding to even S,
C-- in particular, SE= 1/4 for unpolarized spin-1/2 particles;
C-- for nonidentical particles:
C-- SF= RHO (see Eq. (27) in L-L) and SE= 1/2.
C-----------------------------------------------------------------------
C---The following keys and quantities should be set by the user prior
C-- the calculation of the correlation function:
C---Choosing the approximation for Bethe-Salpeter amplitude
C     NS=
C---Choosing the two-particle system
C     LL=
C---Setting ISW
C     ISW=
C---Choosing the law and space-time parameters in fm
C     LAW=
C     ZT=
C     ZR=
C     ZRZ=
C---Setting pair velocity and cos(theta)
C     V=
C     CTHT=
C---Setting azimuthal angle of the vector Q=2*k
C     FIQ=
C---Setting CVK=cos(V*k)
C-- CVK= 0 (1) corresponds to transverse (longit.) configuration
C     CVK=
C---Setting Q GeV/c
C     Q=
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TIM1,TIM2
      COMMON/IQS/IQS
      COMMON/PI/PI
      COMMON/BP/B,P
      COMMON/COULPH/EIDC
      COMMON/SPIN/RHO(10),SE,SO,SF(10)
      COMMON /DDD/DDH(5)
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/WW/WRT,WRL,WT
      COMMON/SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),SBKRB(10),SDKK(10)
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON/IN/XR(200),XT(200),ER(200),ET(200),A1,A2
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMMON/FD/FD(10),FDS(10),RD(10)
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      COMMON/CHL/CHL
      COMMON/C/C(10),AKF(10),AKFS(10),AM,AMS,DM
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      COMPLEX*16 C,CHL,F,EIDC
      COMPLEX*8 Z8,CGAMMA
      EXTERNAL EFH
      COMMON/AA/BBR,CCR,AA,BB,BBT,CC,CCT,EPS
      COMMON/SHH/SH,CHH
      COMMON/AAPI/AAPI(20,2)/AAND/AAND(20,4)/PRAA/PRAA(8,2)
      COMMON/A/A(100)/DAK/DAK,IFUN
      DIMENSION FDH(15,10),RDH(15,10),RBH(15,10),EBH(15,10)
      DIMENSION SEH(15),C0H(15),RHOH(15,10),SFH(15,10)
      DIMENSION S0(2),P0(2),PHI0(2)
      DIMENSION AM1H(15),AM2H(15),C1H(15),C2H(15),MSPINH(15)
C--- mass of the first and second particle
      DATA AM1H/.93956563D0,.93827231D0,.93956563D0,3.7294D0,.13957D0,
     C          .13497D0,.13957D0, .9395731D0, .9382796D0, 6*0.D0/
      DATA AM2H/.93956563D0,.93827231D0,.93827231D0,3.7294D0,.13957D0,
     C          .13497D0,.13957D0, 2*1.87561339D0, 6*0.D0/
C---  charge of the first and second particle
      DATA C1H/0.D0,1.D0,0.D0,2.D0, 1.D0,0.D0,1.D0,0.D0,1.D0, 6*0.D0/
      DATA C2H/0.D0,1.D0,1.D0,2.D0,-1.D0,0.D0,3*1.D0, 6*0.D0/
C---MSPIN vs (LL)
      DATA MSPINH/1,1,2,1,1,1,1,2,2,6*0/
C---Spin factors RHO, SF vs (LL,ISPIN) and SE vs (LL)
      DATA RHOH/3*.25D0, 4*1.D0, 2*.3333D0, 6*0.D0,
     C           0.D0, 0.D0,.75D0,4*0.D0, 2*.6667D0, 126*0.D0/
      DATA SFH/2*.5D0,.25D0,2.D0,1.D0,2*2.D0,2*.3333D0, 6*0.D0,
     C           0.D0, 0.D0,.75D0,4*0.D0, 2*.6667D0, 126*0.D0/
      DATA SEH/2*.25D0,.5D0,1.D0,.5D0,2*1.D0,2*.5D0,6*0.D0/
C-- Probability of the total spin S for unpolarized particles
C-- with spins j1 and j2 is RHO(JJ)=(2*S+1)/[(2*j1+1)*(2*j2+1)],
C-- for spin-1/2 particles with polarizations P1 and P2:
C-- RHO(1)= (1-P1*P2)/4, RHO(2)= (3+P1*P2)/4 for S=0,1, respect.;
C-- SF(JJ), SE, SO, C0=SE-SO (SE+SO=1) are the s-wave spin factors;
C-- for identical particles:
C-- SF(JJ)= 2*RHO(1),0,2*RHO(3),0,..,
C-- SE=RHO(1)+RHO(3)+.., SO=RHO(2)+RHO(4)+..,
C-- SE= 1/4, SO= 3/4, C0= -1/2 for unpolarized spin-1/2 particles;
C-- for nonidentical particles:
C-- SF(JJ)= RHO(JJ), SE=SO=1/2, C0=0.
C---Scattering length FD and effective radius RD in fm vs (LL,ISPIN)
      DATA FDH/17.1D0,7.828D0,23.748D0,
     C 2230.1218D0,.225D0,.081D0,-.063D0,
     C  -.65D0,-2.73D0,6*0.D0,
     C    0.D0, 0.D0,-5.424D0,4*0.D0,-6.35D0,-11.88D0,126*0.D0/
      DATA RDH/2.86D0,2.80D0,2.75D0,
     C 1.12139906D0,-44.36D0,64.0D0,784.9D0,
     C  477.9D0, 2.27D0, 6*0.D0,
     C          .0D0, .0D0,1.786D0,4*0.D0,2.0D0,2.63D0,126*0.D0/
C---Corresponding square well parameters RB (width in fm) and
C-- EB =SQRT(-AM*U) (in GeV/c); U is the well height
      DATA RBH/2.688663D0,   2.781030D0, 2.631925D0, 5.023544D0,
     C         1.124673D0, .3925180D0,.09D0, 2.D0, 4.058058D0, 6*0.D0,
     C        2*.0D0, 2.126224D0, 4*0.D0, 2.D0, 4.132163D0, 126*0.D0/
      DATA EBH/.1085619D0,    .1046383D0,   .1128052D0, .1186010D0,
     C         .7947389D0,2.281208D0, 8.7D0, .4D0, .1561219D0, 6*0.D0,
     C        2*.0D0, .1761675D0, 4*0.D0, .4D0, .1150687D0, 126*0.D0/
      write(*,*) 'Input ZR, ZRZ'
      read(*,*) ZR,ZRZ
      OPEN(6,FILE='mesfsi0pp',STATUS='UNKNOWN') ! message file
      CHL=DCMPLX(.5D0,-.5D0)
      PI=4*DATAN(1.D0)
      SPI=DSQRT(PI)
      SQR2=DSQRT(2.D0)
      W=1/.1973270536D0  !fm to 1/GeV
      PRAA(1,1)=PRAA(1,1)*W
      PRAA(1,2)=PRAA(1,2)*W
C---Choosing the two-particle system (LL= )
CC      LL=2       !!set
      LL=2       !!set
C---Choosing the approximation for Bethe-Salpeter amplitude
cc      NS=4       !!set
      NS=1       !!set
      NSS=NS
C---Choosing the law and space-time parameters in fm
CC    LAW=6      !!set
      LAW=1      !!set
CC    ZT=100.0D0   !!set
      ZT=.01D0   !!set


C      ZR=1.5D0   !!set
C      ZRZ=1.5D0  !!set

C---fm to 1/GeV
      SR=ZR*W
      ST=ZT*W
      SRZ=ZRZ*W
 
      SRZS=SRZ**2
      SRS=SR**2
      STS=ST**2
      SEP=DSQRT(SRS+STS)
      SR1S=SRS+2*SRZS
      SR1=DSQRT(SR1S)
C---setting particle masses and charges
      AM1=AM1H(LL)
      AM2=AM2H(LL)
      C1=C1H(LL)
      C2=C2H(LL)
C---calcul. twice the reduced mass (AM), the relative mass difference
C-- (DM) and the Bohr radius (AC)
      AM=2*AM1*AM2/(AM1+AM2)
      AMS=AM*AM
      DM=(AM1-AM2)/(AM1+AM2)
      AC=1.D10
      C12=C1*C2
      IF(C12.NE.0.D0)AC=2*137.036D0/(C12*AM)
C---Setting IQS
      IQS=0
      IF(C1+AM1.EQ.C2+AM2)IQS=1
C---Setting spin factors
      MSPIN=MSPINH(LL)
      MSP=MSPIN
      DO 91 ISPIN=1,10
      RHO(ISPIN)=RHOH(LL,ISPIN)
  91  SF(ISPIN)=SFH(LL,ISPIN)
      SE=SEH(LL)
      SO=1-SE
      C0=SE-SO
C---Setting ISW and ICH
      ISW=1    !!set
cc      ISW=3    !!set  no strong interaction and no QS
      ICH=0
      IF(AC.LT.1.D9)ICH=1
C---Approx. (r << AC) or exact calc. of hypergeom. functions F and G
C-- IHYPER= 0 OR 1
      IHYPER=0
      IF(AC.LT.1.D9)IHYPER=1
C---Analytical (numerical) integration over azimuth of vector RT
C-- IFI= 0(1)
C---Analyt. integr. (Eq. (18) of L-L) can be used if ICH=0 or IHYPER=0
C-- and the production region is spherically symmetrical and LAW=1.
C-- Note that the accuracy of the approx. F=1 is better than 10% at
C-- R/AC < .05. Thus if <R/AC> << .05 the angular dependence of F may
C-- be neglected and the option IFI=0 may be used simultaneously with
C-- IHYPER=1 to save execution time.
C---Note that at R > AC the f-n F substantially diff. from 1 and forces
C-- the CF to 1 (while F=1 yields CF= Coulomb penetr. factor squared).
      IFI=0
      IF(LAW.NE.1) IFI=1
      IF(ZR.NE.ZRZ) IFI=1
      IF(IHYPER.EQ.1) IFI=1
      IF(LL.NE.4)GOTO 39
C---Setting the value B0=S0(Eg) needed to calculate the l=0 phase shift
C-- for the alpha-alpha system
      AKG=DSQRT(AM*PRAA(3,1)/1000) !momentum [GeV/c] corresponding to
      ETAG=1/AKG/AC                !the ground energy Eg=PRAA(3,1) [MeV]
      RHOG=AKG*PRAA(1,1)
      CALL COUL(1,RHOG,ETAG,S0,P0,PHI0)
      PRAA(2,1)=S0(1)
C---HN= normalization factor
  39  IF(LAW.EQ.1)HN=8*3.14159D0*SRZ*SRS*ST
      IF(LAW.EQ.2)HN=(8/9.D0)*SPI*ST*SRS**3
      IF(LAW.EQ.3)HN=8*SPI*ST*SRS
      IF(LAW.EQ.4)HN=3.14159D0*SPI*ST*SRS
      IF(LAW.EQ.5)HN=8*3.14159D0*SRS*SR*SR1
      IF(LAW.EQ.6)HN=8*SPI*SRZ*SRS*ST
      WRT=.25D0/SRS
      WRL=.25D0/SRZS
      IF(LAW.NE.5.AND.LAW.NE.6)WT=.25D0/STS
      IF(LAW.EQ.5)WRL=.25D0/SR1S
      IF(LAW.EQ.5)WT=WRL*SRZS/SRS
      IF(LAW.EQ.6)WT=1/ST
 4    CONTINUE
C---Setting integration accuracy
      EPS=1.0D-18
C      EPS=1.0D-1
      DDD=.02D0
      DS=1.2D0
      DDH(1)=DDD/DS
      DDH(2)=DDH(1)/DS
      DDH(3)=DDH(2)/DS
C---Integration limits (AA,BB,BBT,BBR) and init. steps (CC,CCT,CCR)
      AA=.0D0
      IF(NS.EQ.2.OR.NS.EQ.4)AA=.1D0
      N=20
      BB=3*(SR+SRZ)
      BBT=4*ST
      IF(LAW.EQ.6)BBT=8*ST
C??      IF(LAW.EQ.5)BBT=BT
      CCT=BBT/N
      CC=BB/N
      BBR=2*SR
      CCR=BBR/N
C---Setting pair velocity and cosine of the production angle
      V=.001D0           !!set
      CTHT=1.D0         !!set
      CTHTS=CTHT*CTHT
      STHTS=1-CTHTS
      STHT=DSQRT(STHTS)
      GAS=1/(1-V**2)
      GAM=DSQRT(GAS)
      GAMV=GAM*V
      RO=DSQRT(SRS+(ST*V)**2*(CTHTS+STHTS*SRZS/SRS))
      PVN=4*SPI*GAM*RO*SRS
      HA=V*SEP/RO        !*
      HA1=DASIN(HA)/HA   !*
C---Setting scatt. length (FD), eff. radius (RD) and, if possible,
C-- also the corresp. square well parameters (EB, RB)
      DO 55 JJ=1,MSP
      ISPIN=JJ
      FD(JJ)=FDH(LL,JJ)*W !in 1/GeV
      RD(JJ)=RDH(LL,JJ)*W
      RB(JJ)=RBH(LL,JJ)*W
      EB(JJ)=EBH(LL,JJ)   !in GeV
      IF(LL.NE.8.AND.LL.NE.9)GOTO25
C---Resets FD and RD for a nucleon-deuteron system (LL=8,9)
      JH=LL-7+2*JJ-2
      FD(JJ)=AAND(1,JH)
      RD(JJ)=AAND(2,JH)-2*AAND(3,JH)/AAND(1,JH)
C---Resets FD and RD for a pion-pion system (LL=5,6,7)
  25  IF(LL.GT.7.OR.LL.LT.5)GOTO 24
      IF(LL.EQ.7)FD(JJ)=AAPI(1,2)/AM
      IF(LL.EQ.5)FD(JJ)=(.6667*AAPI(1,1)+.3333*AAPI(1,2))/AM
      IF(LL.EQ.6)FD(JJ)=(.3333*AAPI(1,1)+.6667*AAPI(1,2))/AM
      AKS=0.D0
      DAKS=1.D-5
      AKSH=AKS+DAKS
      AKH=DSQRT(AKSH)
      GPI1H=GPIPI(AKSH,1)
      GPI2H=GPIPI(AKSH,2)
      H=1/FD(JJ)
      IF(LL.EQ.7)C(JJ)=1/DCMPLX(GPI2H,-AKH)
      IF(LL.EQ.5)
     +C(JJ)=.6667D0/DCMPLX(GPI1H,-AKH)+.3333D0/DCMPLX(GPI2H,-AKH)
      IF(LL.EQ.6)
     +C(JJ)=.3333D0/DCMPLX(GPI1H,-AKH)+.6667D0/DCMPLX(GPI2H,-AKH)
      HH=DREAL(1/C(JJ))
      RD(JJ)=2*(HH-H)/DAKS
 24   CONTINUE
C---Calculation continues for any system (any LL)
      FDS(JJ)=FD(JJ)**2
C---CALL SEAR if you need to calculate the square well parameters
C-- i.e., activate the following two statements
cc      CALL SEAR
cc      IF(JJ.EQ.MSP) STOP
C---Comment out GOTO 55 if you want to print B=R-1 at Q=0
C-- calculated analytically in the spherical wave approximation
C-- for a system without Coulomb interaction; i.e., the
C-- strong FSI contr. to B(Q=0)=BB3 (Eq. (22) of L-L) is calculated,
C-- valid for spherically symmetric gaussian law;
C-- useful to check the integration accuracy.
      GOTO 55
      IF(ZR.NE.ZRZ)GOTO 55
C---PVF=correction to CF due to spherical wave approx. (Eq. (21) of L-L)
      PVF=RD(JJ)*FDS(JJ)/PVN
      BB1=.5D0*FDS(JJ)*HA1/(GAM*SR*RO)
      BB1=SF(JJ)*(BB1-PVF)
      BB2=2*SF(JJ)*FD(JJ)*DLOG((1+HA)/(1-HA))/(GAM*2*HA*SPI*RO)
      BB3=BB1+BB2
      PRINT 54,BB1,BB2,BB3
 54   FORMAT(31X,3E10.3)
 55   CONTINUE
C---Setting azimuthal angle of the vector Q=2*k
      FIQ=0.D0       !!set
      CFIQ=DCOS(FIQ)
C---Setting CVK=cos(V*k)
C-- CVK= 0 (1) corresponds to transverse (longit.) configuration
CC    CVK=0.D0       !!set
      CVK=1.D0       !!set
      SVK=DSQRT(1-CVK**2)
      PRINT 52
      WRITE(3,52)
 52   FORMAT(//1X,3HLAW,2X,2HNS,2X,2HLL,1X,3HIFI,1X,
     C3HIHY,5X,3HSRT,4X,3HSRL,7X,2HST,6X,1HV,4X,5H(V*Z),4X,2HFI,
     C3X,5H(V*K))
      PRINT 152,LAW,NS,LL,IFI,IHYPER,ZR,ZRZ,ZT,V,CTHT,FIQ,CVK
      WRITE(3,152) LAW,NS,LL,IFI,IHYPER,ZR,ZRZ,ZT,V,CTHT,FIQ,CVK
 152  FORMAT(2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,F6.1,1X,F6.1,
     C1X,F8.1,2X,F7.3,1X,F6.2,1X,F6.2,1X,F6.2)
      PRINT 153,(FD(JJ)/W,JJ=1,MSP)
      WRITE(3,153) (FD(JJ)/W,JJ=1,MSP)
 153  FORMAT(2X,17HFD(ISPIN) IN FM =,10F9.3)
      PRINT 154,(RD(JJ)/W,JJ=1,MSP)
      WRITE(3,154) (RD(JJ)/W,JJ=1,MSP)
 154  FORMAT(2X,17HRD(ISPIN) IN FM =,10F9.3)
      PRINT 51
      WRITE(3,51)
C---Setting Q in GeV/c
      DQ=0.01000D0        !!set
CC      DQ=0.00100D0        !!set
cc      DQ=0.000100D0        !!set
CC      Q=-DQ+2*.01860D0    !!set
cc      Q=-DQ+.0510D0    !!set
      Q= 0.0D0            !!set
call datime(idate,itime)
      WRITE(6,*) 'date  ',idate,' time   ',itime
      DO 5 JQT=1,30       !!set
C---Setting initial time
c      CALL CPUTIME(TIM1,IRC1)
call datime(idate,itime)
      tim1=itime
      Q=Q+DQ
      QS=Q*Q
      QLH=Q*CVK
      QT=Q*SVK
      Q0=GAMV*QLH
      QL=Q0/V
      QTS=QT**2
      AKS=.25D0*QS
      AK=.5D0*Q
      IF(LL.EQ.4)E=1000*AKS/AM !energy [MeV] of relative motion
      IF(ICH.EQ.0)GOTO 21
      ACH=0
      ACHR=0
      IF(AK.EQ.0.D0.AND.ICH.EQ.1)GOTO31
      XH=AC*AK
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
      ETA2=ETA*ETA
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 10, 15
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8)
      EIDC=Z8/CABS(Z8)
C---HCP, HPR needed (e.g. in GST) if ICH=1
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
  21  CONTINUE
      B1=0
      B2=0
      BDEL=0
      DO 30 JJ=1,MSP
      ISPIN=JJ
      IF(NS.NE.1)GOTO22
C---Calc. quantities for the square well potential
      BK(JJ)=DSQRT(EB(JJ)**2+AKS)
      X=2*RB(JJ)/AC
      H=BK(JJ)*RB(JJ)
      CALL SEQ(X,H)
      SBKRB(JJ)=H*B
      H=AK*RB(JJ)
      CALL GST(X,H)
      SDK(JJ)=SH
      CDK(JJ)=CHH
      SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK
      IF(ICH.EQ.1)SDK(JJ)=ACH*SDK(JJ)
  22  CONTINUE
C-----------------------------------------------------------------------
C---Calc. the strong s-wave scattering amplitude = C(JJ)
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(NS.NE.1)GOTO 230
      IF(LL.NE.4)GOTO 230 ! SW scat. amplitude used for alfa-alfa only
      AKACH=AK
      IF(ICH.EQ.1)AKACH=AK*ACH
      C(JJ)=1/DCMPLX(G(AK),-AKACH) ! amplitude for the SW-potential
      GOTO 23
 230  IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)GOTO20
      IF(LL.EQ.8.OR.LL.EQ.9)GOTO20
      A1=RD(JJ)*FD(JJ)*AKS
      A2=1+.5D0*A1
      IF(ICH.EQ.1)A2=A2-2*HCP*FD(JJ)/AC
      A2S=A2**2
      AKF(JJ)=AK*FD(JJ)
      IF(ICH.EQ.1)AKF(JJ)=AKF(JJ)*ACH
      C(JJ)=FD(JJ)/DCMPLX(A2,-AKF(JJ))
      GOTO23
 20   CONTINUE
C---Calc. scatt. ampl. C(JJ) and eff. radii RD(JJ) for pipi and Nd
C-- systems (RD dep. on k)
      DAKS=1.D-5
      AKSH=AKS+DAKS
      AKH=DSQRT(AKSH)
      JH=LL-7+2*JJ-2
      IF(LL.GT.7)GPI2H=GND(AKSH,JH)
      IF(LL.LE.7)GPI2H=GPIPI(AKSH,2)
      C(JJ)=1.D0/DCMPLX(GPI2H,-AKH)
      IF(LL.GE.7)GOTO28
      GPI1H=GPIPI(AKSH,1)
      IF(LL.EQ.5)C(JJ)=.6667D0/DCMPLX(GPI1H,-AKH)+.3333D0*C(JJ)
      IF(LL.EQ.6)C(JJ)=.3333D0/DCMPLX(GPI1H,-AKH)+.6667D0*C(JJ)
      C(2)=SQR2*(C(JJ)-DCMPLX(GPI1H,-AKH))/3  !*
      HHTR=DREAL(1/C(2))                      !*
 28   HH=DREAL(1/C(JJ))
      IF(LL.GT.7)GPI2=GND(AKS,JH)
      IF(LL.LE.7)GPI2=GPIPI(AKS,2)
      C(JJ)=1.D0/DCMPLX(GPI2,-AK)
      IF(LL.GE.7)GOTO27
      GPI1=GPIPI(AKS,1)
      IF(LL.EQ.5)C(JJ)=.6667D0/DCMPLX(GPI1,-AK)+.3333D0*C(JJ)
      IF(LL.EQ.6)C(JJ)=.3333D0/DCMPLX(GPI1,-AK)+.6667D0*C(JJ)
      C(2)=SQR2*(C(JJ)-DCMPLX(GPI1,-AK))/3  !*
      HTR=DREAL(1/C(2))                     !*
      RDTR=2*(HHTR-HTR)/DAKS                !*
 27   H=DREAL(1/C(JJ))
      RD(JJ)=2*(HH-H)/DAKS
CC      PRINT 54,RD(JJ)/W,FD(JJ)/W
C---Calc. pi+pi-, pi+pi+ or pd s-wave scatt. amplitude
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(ICH.EQ.0)GOTO 23
      AAK=ACH*AK
      HCP2=2*HCP/AC
      C(JJ)=1/(1/C(JJ)-HCP2+DCMPLX(0.D0,AK-AAK))
      IF(LL.EQ.5) C(2)=1/(1/C(2)-HCP2+DCMPLX(0.D0,AK-AAK)) !*
 23   CONTINUE  !Calculation continues for any system (any LL)
      AKFS(JJ)=SF(JJ)*(DREAL(C(JJ))**2+DIMAG(C(JJ))**2)
      FTRS=0.D0                                     !*
      IF(LL.EQ.5.OR.LL.EQ.6)                        !*
     & FTRS=SF(1)*(DREAL(C(2))**2+DIMAG(C(2))**2)   !*
C---B2= correction to B due to the WF distortion inside the range
C-- of strong interaction (generalization of Eq. (21) of L-L to the
C-- spher. nonsymmetrical gaussian law)
      IF((NS.EQ.2.OR.NS.EQ.4).AND.(LAW.EQ.1.AND.ISW.EQ.1))
     & B2=B2-RD(JJ)*AKFS(JJ)/PVN
     &      -RDTR*FTRS/PVN          !*
     &      +2*SPI*SR*HA1*FTRS/PVN  !*
      IF(LL.NE.4)GOTO 30
CC      BDEL=BDEL+SF(JJ)*(DERD(1,E)+5*DERD(2,E))/AKS/(2*PI*PVN)
      BDEL=BDEL+SF(JJ)*(5*DERD(2,E))/AKS/(2*PI*PVN)
 30   CONTINUE
C***********************************************************************
C---Integrating the f-n (R-1) over the distr. of the emission points
      IF(IFI.NE.0) GO TO 9
C---Analytical integration over the azimuth of vector RT
      AI=EFH(0.D0)
      GO TO 8
C---Numerical integration over the azimuth of vector RT;
C-- AI= Integral EFH(FIH) from 0 to PI;
C-- in some cases it might be necessary to integr. from 0 to 2*PI.
  9   FIMAX=PI
      CALL DSIMPV(.0D0,FIMAX,.314159D0,DDD,EPS,EFH,FIH,AI,A1H,ABS)
      AI=AI/FIMAX
  8   CONTINUE
      B1=AI/HN
C---Multiply by 2 since integr. is over the positive time diff. only
      B1=2*B1
 31   CONTINUE
      B3=B1+B2
C---B4= quantum statistics contrib. to B (B4=0 for nonidentical part.)
C-- calculated if LAW= 1; B4 used if IFI=0.
      B4=0
      IF(LAW.GT.1)GOTO 10
C---Calcul. B4 if LAW=1
      HH=SRZS-SRS
      ATT=SRS+HH*STHTS*CFIQ**2
      ALL=V**2*STS+SRZS-HH*STHTS
      ALT=2*HH*CTHT*STHT*CFIQ
      B4=C0*EXF(-ATT*QTS+ALT*QL*QT-ALL*QL*QL)
 10   CONTINUE
C---B8= B= R - 1
      B8=B3
      IF(IFI.EQ.0) B8=B3+B4
      R8=B8+1
      IF(ICH.EQ.1)R8=ACH*R8
      IF(LL.EQ.4)R8=R8+BDEL
C---Calcul. current time
53    continue
      call datime(idate,itime)
      tim2=itime
c      CALL CPUTIME(TIM2,IRC2)
      TIM=(TIM2-TIM1)/1.E6
      PRINT 50,Q,B1,B2,B3,B4,B8+1,R8,TIM
      WRITE(3,50) Q,B1,B2,B3,B4,B8+1,R8,TIM !Q
      WRITE(33,50) Q/2.,R8 !Q/2
51    FORMAT(/7X,1HQ,
     C8X,2HB1,8X,2HB2,8X,2HB3,8X,2HB4,8X,2HRN,9X,1HR,6X,4HTIME)
50    FORMAT(1X,8E10.3)
 5    CONTINUE
      END
      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AAPI/AAPI(20,2)/AAND/AAND(20,4)/PRAA/PRAA(8,2)
      COMMON/AU/AU(100)/AL/AL(100)
C---Parameters for GPIPI (I,J), J=1(2) corresp. to isospin=0(2)
      DATA AAPI/.2600D00, .2500D00, .3480D00,-.0637D00, 16*0.D0,
     1         -.0280D00,-.0820D00, .2795D00,-.0086D00, 16*0.D0/
C---Parameters for GND (I,J), J=1(2) corresp. to nd(pd), S=1/2,
C---                          J=3(4) corresp. to nd(pd), S=3/2
      DATA AAND/-.3295D1,-.8819D3, .4537D4, .28645D5, 16*0.D0,
     1          -.13837D2, .11505D2, .0D0, .10416D2,  16*0.D0,
     2          -.32180D2, .1014D2,  .0D0, .0D0,      16*0.D0,
     3          -.60213D2, .1333D2,  .0D0,-.70309D2,  16*0.D0/
C---Parameters PRAA(I,J) for DELAA (J), J=1(2) corresp. to l=0(2)
C-- a[fm], B, 3*(gam^2[MeV], E[MeV]); B(l=0) should be recalculated
C-- as B(l=0)=S(Eg=.09212 MeV)
      DATA PRAA/7.00D0,.0D0,.09212D0,.087D0,6.30D0,.616D0,22.8D0,1.45D0,
     1          6.75D0,.00D0,2.560D0,.387D0,8.07D0,.746D0,34.1D0,2.68D0/
C---Upper and lower limits on the parameters in SEAR
      DATA AU/1.D3,1.D2,98*0.D0/,AL/1.D-1,.1D0,98*0.D0/
      END
      SUBROUTINE VZ1
C---CALC. THE FSI CONTR. F1 TO B.
C-- IF IFI=1 THEN F1 INCLUDES ALSO THE EFFECT OF QUANTUM STATISTICS.
C-- IF ISW=1 COMPLETE CALCULATION;
C-- IF ISW=2 THE STRONG INTERACTION IS SWITCHED OFF;
C-- IF ISW=3 IN ADDITION THE QUANTUM STATISTICS IS SWITCHED OFF.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/IQS/IQS
      COMMON/SPIN/RHO(10),SE,SO,SF(10)
      COMMON/ES/E,ES
      COMMON/WW/WRT,WRL,WT
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/AA/BBR,CCR,AA,AY(5)
      COMMON/FFF/FS,FA,F12,F21
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON/FD/FD(10),FDS(10),RD(10)
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      COMMON/C/C(10),AKF(10),AKFS(10),AM,AMS,DM
      COMPLEX*16 F,C,G
       COMPLEX*16 FS,FA,F12,F21
      RHS=RH**2
C---RL AND T ARE THE LONGITUDINAL-SPACE AND TIME DISTANCES
C-- IN THE PAIR REST FRAME
      RL=GAM*(RLH-V*TH)
      T=GAM*(TH-V*RLH)
      RLS=RL**2
      RPS=RTS+RLS
      RP=DSQRT(RPS)
      IF(RP.LE.AA)RETURN
      RZ=RLH*CTHT-RT*STHT*CFI
      RX=RLH*STHT+RT*CTHT*CFI
      RTTS=RTS*(1-CFI*CFI)+RX*RX
      RZS=RZ*RZ
      IF(LAW.EQ.1)EX=EXF(-RZS*WRL-RTTS*WRT-TH*TH*WT)
      IF(LAW.EQ.2)EX=EXF(-TH*TH*WT)*(.6667D0*SRS*SR
     --.5D0*SRS*RH+RH*RHS/24)
      IF(LAW.EQ.3)EX=EXF(-TH*TH*WT)/RH
      IF(LAW.EQ.5)EX=EXF(-(RZS+TH*TH)*WRL-(RZ-TH)**2*WT-RTTS*WRT)
      IF(LAW.EQ.6)EX=EXF(-RZS*WRL-RTTS*WRT-TH*WT)
      IF(LAW.NE.4)GOTO2
C---SPHERE SURFACE + LAMBERT LAW
      Y=.25D0*RHS/SRZS
      EX=EXF(-TH*TH*WT)*FW(Y)/RH
2     CONTINUE
      IF(EX.EQ.0.D0)RETURN
      IF(IFI.EQ.0)GOTO 4
      RHOS=AK*RP               !k*r
      HT=.5D0*QT*RT
      HL=.5D0*Q0*RL/GAMV
      HS=HT*DCOS(FI-FIQ)+HL    !k-vector*r-vector
      CALL FF(RHOS,HS)
C---GOTO 5 IF YOU WANT TO SWITCH OFF THE STRONG INTERACTION
C---GOTO 7 IF YOU WANT IN ADDITION SWITCH OFF QUANTUM STATISTICS
      GOTO(4,5,7)ISW
 4    CONTINUE
      F1=-1.D0
      CALL FIRT
      MSP=MSPIN
      DO 10 JJ=1,MSP
      G=F(JJ)*C(JJ)
      IF(IQS.EQ.1.OR.IFI.EQ.0)GOTO 9
      G=F12+G
      F1=F1+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
      GOTO 10
 9    A=DREAL(F(JJ))**2+DIMAG(F(JJ))**2
      F1=AKFS(JJ)*A+F1
      IF(IFI.EQ.0)GOTO6
      F2=SF(JJ)*DREAL(G*DCONJG(FS))
      IF(JJ.GT.1)GOTO 1
      FH=SE*(DREAL(FS)**2+DIMAG(FS)**2)
      FH=FH+SO*(DREAL(FA)**2+DIMAG(FA)**2)
      F1=.5D0*FH + F1
      GOTO1
 6    F2=2*SF(JJ)*DREAL(G)
      IF(QS.EQ.0.)GOTO1
C---ANALYT. INTEGR. OVER THE AZIMUTH OF VECTOR RT (EQ. (18) OF L-L)
C      PRINT 100,RT
 100  FORMAT(1X,E15.5)
      BES=DBESJ0(.5D0*QT*RT)
      HH=.5D0*Q0*RL/GAMV
      F2=F2*BES*DCOS(HH)
      IF(IQS.EQ.0)F2=F2-2*SF(JJ)*DIMAG(G)*BES*DSIN(HH)
 1    F1=F1+F2
 10   CONTINUE
      GOTO 8
 5    IF(IQS.EQ.0)GOTO 7
      FH=SE*(DREAL(FS)**2+DIMAG(FS)**2)
      FH=FH+SO*(DREAL(FA)**2+DIMAG(FA)**2)
      F1=.5D0*FH-1
      GOTO 8
 7    FH=DREAL(F12)**2+DIMAG(F12)**2
      IF(IQS.EQ.1)FH=.5D0*(DREAL(F21)**2+DIMAG(F21)**2+FH)
      F1=FH-1
 8    CONTINUE
      F1=F1*RHS*EX
      RETURN
      END
      FUNCTION GPIPI(X,J)
C--- GPIPI = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to isospin=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AAPI/AAPI(20,2)
      COMMON/C/HELP(40),AM,AM2,DM
      OM=DSQRT(X+AM2)
      XX=X/AM2
      GPIPI=OM/AAPI(1,J)
      GPIPI=GPIPI*(1+(AAPI(3,J)-AAPI(1,J)**2)*XX+AAPI(4,J)*XX*XX)
      GPIPI=GPIPI/(1+(AAPI(3,J)+AAPI(2,J)/AAPI(1,J))*XX)
      RETURN
      END
      FUNCTION GND(X,J)
C--- GND = k*COTG(DELTA), X=k^2
C--- J=1(2) corresp. to nd(pd), S=1/2,
C--- J=3(4) corresp. to nd(pd), S=3/2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AAND/AAND(20,4)
      XX=X
      GND=1/AAND(1,J)+.5*AAND(2,J)*X
      DO 1 I=4,4
      XX=XX*X
   1  GND=GND+AAND(I,J)*XX
      GND=GND/(1+AAND(3,J)*X)
      RETURN
      END
      FUNCTION RINVS(J,E)
C---Calculates 1/R for J=1(2) corr. to l=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AAPI/AAPI(20,2)/AAND/AAND(20,4)/PRAA/PRAA(8,2)
      RINVS=.0D0
      R=0.D0
      DO 1 K=1,3
      K1=2*K+1
      A=PRAA(K1,J)-E
      IF(A.EQ.0.D0)RETURN
   1  R=R+PRAA(K1+1,J)/A
      IF(R.NE.0.D0)RINVS=1/R
      IF(R.EQ.0.D0)RINVS=1.D20
      RETURN
      END
      FUNCTION DELAA(J,E)
C---calc. phase shift for alpha-alpha at E[MeV]; J=1(2) --> l=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PI/PI
      COMMON/AAPI/AAPI(20,2)/AAND/AAND(20,4)/PRAA/PRAA(8,2)
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      DIMENSION S(2),P(2),PHI(2)
      RHO=AK*PRAA(1,J)
      CALL COUL(J,RHO,ETA,S,P,PHI)
      RI=RINVS(J,E)
      DELAA=-PHI(J)+DATAN2(P(J),PRAA(2,J)-S(J)+RI)
      IF(DELAA.GT.PI)DELAA=DELAA-PI
      IF(DELAA.LE.-PI)DELAA=DELAA+PI
      IF(E.LT.PRAA(5,J).AND.DELAA.LT.0.D0)DELAA=DELAA+PI
      IF((J.EQ.1.AND.E.GE.PRAA(5,J)).AND.DELAA.GE.PI/2)DELAA=DELAA-PI
CC      WRITE(3,100) J,E,DELAA*180/PI,RI
 100  FORMAT(' J=',I2,' E=',E12.4,' DEL=',E12.4,' RINV',E12.4)
      RETURN
      END
      FUNCTION DERD(J,E)          ! dDELAA(J,E)/dk
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AAPI/AAPI(20,2)/AAND/AAND(20,4)/PRAA/PRAA(8,2)
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      COMMON/C/HELP(40),AM,AM2,DM
      DIMENSION S(2),P(2),PHI(2)
      D=1.D-10
      D0=DELAA(J,E)
      AK1=AK+D            !step= 100eV/c --> 1eV in E
      E1=1000*AK1*AK1/AM
      D1=DELAA(J,E1)
      DERD=(D1-D0)/D
      RETURN
      END
      SUBROUTINE COUL(J,RHO,ETA,S,P,PHI)
C---J=1(2) corresponds to the orbital angular mom. l=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 RHO1,ETA1,ZL,F,FP,G,GP,SIG
      DIMENSION F(0:3),FP(0:3),G(0:3),GP(0:3),SIG(0:3)
      DIMENSION S(2),P(2),PHI(2)
      ZL=DCMPLX(0.D0,0.D0)
      ETA1=DCMPLX(ETA,0.D0)
      RHO1=DCMPLX(RHO,0.D0)
      CALL WCLBES(RHO1,ETA1,ZL,2*(J-1),F,G,FP,GP,SIG,-1,1,IFAIL,1)
      DO 1 K=1,J
      II=2*(K-1)
      FR=DREAL(F(II))
      GR=DREAL(G(II))
      A=FR*FR+GR*GR
      P(K)=RHO/A
      S(K)=P(K)*(FR*DREAL(FP(II))+GR*DREAL(GP(II)))
CC      IF(GR.EQ.0.D0)GR=1.D-20
      PHI(K)=DATAN2(FR,GR)
CC      WRITE(3,100) K,ETA,RHO,DREAL(SIG(K))
CC      WRITE(3,101) FR,GR,S(K),P(K),PHI(K)
 100  FORMAT(' K',I3,' ETA=',E12.4,' RHO=',E12.4,' SIG',E12.4)
 101  FORMAT(' F,G,S,P,PHI',I2,5E12.4)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE VZ
C---CALC. THE FSI CONTR. F1 TO B.
C-- IF IFI=1 THEN F1 INCLUDES ALSO THE EFFECT OF QUANTUM STATISTICS.
C-- IF ISW=1 COMPLETE CALCULATION;
C-- IF ISW=2 THE STRONG INTERACTION IS SWITCHED OFF;
C-- IF ISW=3 IN ADDITION THE QUANTUM STATISTICS IS SWITCHED OFF.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SPIN/RHO(10),SE,SO,SF(10)
      COMMON/ES/E,ES
      COMMON/WW/WRT,WRL,WT
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/AA/BBR,CCR,AA,AY(5)
      COMMON/FFF/FS,FA,F12,F21
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON/FD/FD(10),FDS(10),RD(10)
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      COMMON/C/C(10),AKF(10),AKFS(10),AM,AM2,DM
      COMPLEX*16 F,C,G
       COMPLEX*16 FS,FA,F12,F21
      F1=0.D0
      RHS=RH**2
C---RL AND T ARE THE LONGITUDINAL-SPACE AND TIME DISTANCES
C-- IN THE PAIR REST FRAME
      RL=GAM*(RLH-V*TH)
      T=GAM*(TH-V*RLH)
      RLS=RL**2
      RPS=RTS+RLS
      RP=DSQRT(RPS)
      IF(RP.LE.AA)RETURN
      RZ=RLH*CTHT-RT*STHT*CFI
      RX=RLH*STHT+RT*CTHT*CFI
      RTTS=RTS*(1-CFI*CFI)+RX*RX
      RZS=RZ*RZ
      IF(LAW.EQ.1)EX=EXF(-RZS*WRL-RTTS*WRT-TH*TH*WT)
      IF(LAW.EQ.2)EX=EXF(-TH*TH*WT)*(.6667D0*SRS*SR
     --.5D0*SRS*RH+RH*RHS/24)
      IF(LAW.EQ.3)EX=EXF(-TH*TH*WT)/RH
      IF(LAW.EQ.5)EX=EXF(-(RZS+TH*TH)*WRL-(RZ-TH)**2*WT-RTTS*WRT)
      IF(LAW.EQ.6)EX=EXF(-RZS*WRL-RTTS*WRT-TH*WT)
      IF(LAW.NE.4)GOTO2
C---SPHERE SURFACE + LAMBERT LAW
      Y=.25D0*RHS/SRZS
      EX=EXF(-TH*TH*WT)*FW(Y)/RH
2     CONTINUE
      IF(EX.EQ.0.D0)RETURN
      IF(IFI.EQ.0)GOTO 4
      RHOS=AK*RP               !k*r
      HT=.5D0*QT*RT
      HL=.5D0*Q0*RL/GAMV
      HS=HT*DCOS(FI-FIQ)+HL    !k-vector*r-vector
      CALL FF(RHOS,HS)
C---GOTO 5 IF YOU WANT TO SWITCH OFF THE STRONG INTERACTION
C---GOTO 7 IF YOU WANT IN ADDITION SWITCH OFF QUANTUM STATISTICS
      GOTO(4,5,7)ISW
 4    CONTINUE
      CALL FIRT
      MSP=MSPIN
      DO 1 JJ=1,MSP
      A=DREAL(F(JJ))**2+DIMAG(F(JJ))**2
      F1=AKFS(JJ)*A+F1
      IF(IFI.EQ.0)GOTO6
      G=F(JJ)*C(JJ)
      F2=SF(JJ)*DREAL(G*DCONJG(FS))
      IF(JJ.GT.1)GOTO 1
      FH=SE*(DREAL(FS)**2+DIMAG(FS)**2)
      FH=FH+SO*(DREAL(FA)**2+DIMAG(FA)**2)
      F1=.5D0*FH + F1 - 1
      GOTO1
 6    F2=2*SF(JJ)*DREAL(F(JJ)*C(JJ))
      IF(QS.EQ.0.)GOTO1
C---ANALYT. INTEGR. OVER THE AZIMUTH OF VECTOR RT (EQ. (18) OF L-L)
C      PRINT 100,RT
 100  FORMAT(1X,E15.5)
      F2=F2*DBESJ0(.5D0*QT*RT)*DCOS(.5D0*Q0*RL/GAMV)
 1    F1=F1+F2
      GOTO 8
 5    FH=SE*(DREAL(FS)**2+DIMAG(FS)**2)
      FH=FH+SO*(DREAL(FA)**2+DIMAG(FA)**2)
      F1=.5D0*FH-1
      GOTO 8
 7    FH=DREAL(F12)**2+DIMAG(F12)**2
      FH=DREAL(F21)**2+DIMAG(F21)**2+FH
      F1=.5D0*FH-1
 8    CONTINUE
      F1=F1*RHS*EX
      RETURN
      END
      SUBROUTINE FIRT
C---CALC. THE F(JJ)
C-- F(JJ)*C(JJ)= DEVIATION OF THE BETHE-SALPETER AMPL. FROM PLANE WAVE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SHH/SH,CHH
      COMMON/BP/B,P
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/C/C(10),AKF(10),AKFS(10),AM,AM2,DM
      COMMON/SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),SBKRB(10),SDKK(10)
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON/Q/QS,QTS,QT,Q0,V,GAM,GAS,GAMV,AK,AKS,CTHT,STHT,FIQ
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMMON/CHL/CH
      EQUIVALENCE(RSS,RP),(TSS,T)
      COMPLEX*16 F,C,CS,CH,CH1,CH2,CH3
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      MSP=MSPIN
      DO 10 JJ=1,MSP
      IF(JJ.GT.1)GOTO 3
      IF(NSS.EQ.3)GOTO7
      X=2*RSS/AC
      IF(AK.NE.0.D0)GOTO2
      SHK=1.D0
      SH=.0D0
      SHH=SH
      CHH=1/RSS
      GOTO3
  2   H=AK*RSS
      CALL GST(X,H)
      SH=SH/RSS
      CHH=CHH/RSS                             ! ReG(k,r)/r
      SHH=SH
      IF(ICH.EQ.1) SHH=ACH*SH                 ! ImG(k,r)/r
  3   IF(NSS.EQ.2)GOTO1
C---F= ASYMPTOTIC FORMULA (T= 0 APPROX.); NS= 4
  6   F(JJ)=DCMPLX(CHH,SHH)                   ! G(k,r)/r
      IF(NSS.EQ.4) GOTO 10
      IF(NSS.NE.1)GOTO 10
C---F INSIDE THE SQUARE-WELL (T= 0 APPROX.); NS= 1
      IF(RSS.GE.RB(JJ)) GOTO 10
      IF(AK.NE.0.D0.AND.JJ.EQ.1)SHK=SH/AK     ! B(k,r)
      H=BK(JJ)*RSS
      CALL GST(X,H)
      SKR=B*BK(JJ)                           ! kappa*B(kappa,r)
      F(JJ)=DCMPLX(CDK(JJ),SDK(JJ))*SKR      ! G(k,d)*kappa*B(kappa,r)
      CH1=(SDKK(JJ)*SKR-SHK*SBKRB(JJ))/C(JJ) !kappa*d*[B(k,d)*B(kappa,r)
                                             !     -B(k,r)*B(kappa,d)]/f
      F(JJ)=(F(JJ)+CH1)/SBKRB(JJ) ! [G(k,d)/d+B(k,d)/f]*B(kappa,r)/
                                  !  /B(kappa,d) - B(k,r)/f
      GOTO 10
  1   CONTINUE
C---F= ASYMPTOTIC FORMULA (T= 0 NOT REQUIRED); NS= 2
      IF(JJ.GT.1)GOTO 8
      IF(TSS.EQ.0.D0)GOTO6
      TSSA=DABS(TSS)
      IF(DM.NE.0.D0)GOTO 11
      H=AM*.5D0/TSSA
      IF(AK.NE.0.D0)GOTO4
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      S1=DFRSIN(HM)
      C1=DFRCOS(HM)
      C2=C1
      S2=S1
      GOTO5
  4   CONTINUE
      H1=AK*TSSA/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      S1=DFRSIN(HM)
      C1=DFRCOS(HM)
      S2=DFRSIN(HP)
      C2=DFRCOS(HP)
      GOTO 5
  11  CONTINUE
      S1=0.D0
      S2=0.D0
      C1=0.D0
      C2=0.D0
      DO 13 I=1,2
      IF(I.EQ.1)TSSH=TSSA*(1+DM)
      IF(I.EQ.2)TSSH=TSSA*(1-DM)
      H=AM*.5D0/TSSH
      IF(AK.NE.0.D0)GOTO 12
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      S1=S1+DFRSIN(HM)/2
      C1=C1+DFRCOS(HM)/2
      IF(I.EQ.1)GOTO 13
      C2=C1
      S2=S1
      GOTO 13
  12  CONTINUE
      H1=AK*TSSH/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      S1=S1+DFRSIN(HM)/2
      C1=C1+DFRCOS(HM)/2
      S2=S2+DFRSIN(HP)/2
      C2=C2+DFRCOS(HP)/2
  13  CONTINUE
  5   C12=C1+S2
      S12=S1+C2
      A12=S1-C2
      A21=S2-C1
      A2=.5D0*(CHH*(A12+A21)+SH*(A12-A21))+SHH
      A1=.5D0*(CHH*(C12+S12)+SH*(C12-S12))
      F(JJ)=.3989422D0*DCMPLX(A1,A2)
      GOTO 10
C---F(Q=0) FOR TWO-NEUTRON EXPONENTIAL POTENTIAL (T= 0 APPROX.); NS= 3
  7   IF(JJ.GT.1)GOTO 8
      X=2.286D0*EXF(-.1417D0*RSS)
      G1=(1.058215D0*DBESJ0(X)-.12864D0*DBESY0(X))/RSS-.011606D0
      F(JJ)=DCMPLX(G1,.0D0)
      GOTO 10
  8   F(JJ)=F(1)
 10   CONTINUE
      RETURN
      END
      FUNCTION EFH(FIH)
C---EFH(FIH)= INTEGRAL EF(E) FROM -1 TO 1;
C-- E=(R*Z) IS THE COSINE OF THE POLAR ANGLE OF THE VECTOR R=R1-R2.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /DDD/DDD,DH(4)
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/AA/BBR,CCR,AA,BB,BBT,CC,CCT,EPS
      EXTERNAL EF
      FI=FIH
      CFI=DCOS(FI)
      CALL DSIMPU(-1.D0,1.D0,.2D0,DDD,EPS,EF,E,AI,A1H,ABS)
      EFH=AI
      RETURN
      END
      FUNCTION EF(E)
C---EF(E)= INTEGRAL FUN(R) FROM 0 (AA) TO B;
C-- R IS THE DISTANCE BETWEEN THE EMISSION POINTS (MODULUS OF THE
C-- VECTOR R=R1-R2).
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/IN/XR(200),XT(200),ER(200),ET(200),A1,A2
      CALL INF(E)
      EF=A1
      RETURN
      END
      SUBROUTINE INF(E)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON /DDD/D1,DD,DH(3)
      COMMON/IN/XR(200),XT(200),ER(200),ET(200),A1,A2
      COMMON/ES/Z,ES
      COMMON/AA/BBR,CCR,AA,BB,BBT,CC,CCT,EPS
      EXTERNAL FUN
      ES=DSQRT(1-E**2)
      B=BB
      IF(LAW.GT.1.AND.LAW.LT.6)B=BBR
      IF(LAW.NE.4)GOTO1
      B=BBR*ES
      IF(B.LT.AA)B=AA
  1   Z=E
      CALL DSIMPS(AA,B,CCR,DD,EPS,FUN,R,A1,A1H,AB)
      RETURN
      END
      FUNCTION FUN(R)
C---FUN(R)= INTEGRAL FIN(X) FROM 0 (AA) TO INFINITY (BBT);
C-- X=T1-T2 IS THE DIFFERENCE OF THE EMISSION TIMES.
C-- FIN(X)= F1 IS CALCULATED IN SUBROUTINE VZ.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /DDD/D1(2),DD,DH(2)
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      COMMON/ES/E,ES
      COMMON/AA/BBR,CCR,AA,BB,BBT,CC,CCT,EPS
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMPLEX*16 F
      EXTERNAL FIN
      RH=R
C---RLH AND RT ARE THE LONGITUDINAL AND TRANSVERSE (WITH RESPECT TO
C-- THE PAIR VELOCITY) COMPONENTS OF THE SPACE DISTANCE OF THE
C-- EMISSION POINTS IN THE SOURCE REST FRAME.
      RLH=E*R
      RT=ES*R
      RTS=RT**2
      CALL DSIMPT(0.D0,BBT,CCT,DD,EPS,FIN,X,AI,AIH,AB)
      FUN=AI
      RETURN
      END
      FUNCTION FIN(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RR/SRS,STS,SRZS,RL,RT,RTS,RLS,T,F(10),F1,F2
      COMMON/RP/RP,RPS,SEG,RLH,TH,RH
      COMPLEX*16 F
C---TH= DIFFERENCE OF THE EMISS. TIMES IN THE SOURCE REST FRAME
      TH=X
      CALL VZ
      FIN=F1
      RETURN
      END
      FUNCTION EXF(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.LT.-15.D0) GO TO 1
      EXF=DEXP(X)
      RETURN
  1   EXF=.0D0
      RETURN
      END
      SUBROUTINE SEQ(X,H)
C---CALC. FUNCTIONS B, P (EQS. (17) OF G-K-L-L);
C-- NEEDED TO CALC. THE CONFLUENT HYPERGEOMETRIC FUNCTION GST.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BP/B,P
      DIMENSION BH(3),PH(3)
      DATA ERR/1.D-7/
      BH(1)=1.D0
      PH(1)=1.D0
      PH(2)=.0D0
      BH(2)=.5D0*X
      B=1+BH(2)
      P=1.D0
      HS=H*H
      J=0
  2   J=J+1
      BH(3)=(X*BH(2)-HS*BH(1))/((J+1)*(J+2))
      PH(3)=(X*PH(2)-HS*PH(1)-(2*J+1)*X*BH(2))/(J*(J+1))
      B=B+BH(3)
      P=P+PH(3)
      Z=DABS(BH(2))+DABS(BH(3))+DABS(PH(2))+DABS(PH(3))
      IF(Z.LT.ERR)RETURN
      BH(1)=BH(2)
      BH(2)=BH(3)
      PH(1)=PH(2)
      PH(2)=PH(3)
      GOTO 2
      END
      SUBROUTINE SEQA(X,H)
C---CALC. FUNCTIONS CHH=REAL(GST), SH=IMAG(GST)/ACH, B=SH/H
C-- IN THE ASYMPTOTIC REGION H=K*R >> 1.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BP/B,P
      COMMON/SHH/SH,CHH
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/COULPH/EIDC
      COMPLEX*16 EIDC,GST
      ARG=H-ETA*DLOG(2*H)
      GST=DCMPLX(DCOS(ARG),DSIN(ARG))
      GST=ACHR*EIDC*GST
      CHH=DREAL(GST)
      SH=DIMAG(GST)/ACH
      B=SH/H
      RETURN
      END
      SUBROUTINE FF(RHO,H)
C---CALC. F12, F21, FS, FA;
C-- F12= F12* PLANE WAVE, WHERE F12 IS THE CONFL. HYPERGEOM. F-N
C-- (EQ. (15) OF G-K-L-L), F12= 1 AT R << AC,
C-- FS=F12+F21, FA=F12-F21. NOTE THAT
C-- FS (FA)= SYMMETRIC (ANTISYMMETRIC) COMBIN. OF THE PLANE WAVES
C-- AT R << AC.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/FFF/FS,FA,F12,F21
      COMPLEX*16 FAS,ALF,ALF1,Z12,Z21,F12,F21,S12,S21,A,FS,FA
      DATA ERR/1.D-5/
      S12=DCMPLX(1.D0,0.D0)
      S21=S12
      F12=S12
      F21=S12
      C=DCOS(H)
      S=DSIN(H)
      IF(ICH.EQ.0)GOTO2
      IF(IHYPER.EQ.0)GOTO 2
      RHOP=RHO+H
      RHOM=RHO-H
      ALF=DCMPLX(.0D0,-ETA)
      ALF1=ALF-1
      IF(RHOP.LT.15.D0)GOTO5
      F12=FAS(RHOP)
      GOTO6
  5   Z12=DCMPLX(.0D0,RHOP)
      J=0
  3   J=J+1
      A=(ALF1+J)/(J*J)
      S12=S12*A*Z12
      F12=F12+S12
      ZR12=DABS(DREAL(S12))
      ZI12=DABS(DIMAG(S12))
      IF((ZR12+ZI12).GT.ERR)GOTO3
 6    IF(RHOM.LT.15.D0)GOTO8
      F21=FAS(RHOM)
      GOTO2
 8    Z21=DCMPLX(.0D0,RHOM)
      J=0
 4    J=J+1
      A=(ALF1+J)/(J*J)
      S21=S21*A*Z21
      F21=F21+S21
      ZR21=DABS(DREAL(S21))
      ZI21=DABS(DIMAG(S21))
      IF((ZR21+ZI21).GT.ERR)GOTO4
 2    CONTINUE
CC      PRINT 201,2*RHO*ETA,RHO+H,DREAL(F12),DIMAG(F12)
CC      A=FAS(RHO+H)
CC      PRINT 202,DREAL(A),DIMAG(A)
 201  FORMAT(1X,4E15.5)
 202  FORMAT(31X,2E15.5)
      F12=F12*DCMPLX(C,-S)
      F21=F21*DCMPLX(C,S)
      FS=F12+F21
      FA=F12-F21
      RETURN
      END
      FUNCTION ACP(X)
C--- ACP = COULOMB PENETRATION FACTOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.EQ.0D0)GOTO 1
      Y=6.2831853D0/X
      IF(Y.GT.100.D0) GO TO 1
      ACP=Y/(DEXP(Y)-1)
      RETURN
   1  ACP=.0D0
      RETURN
      END
      FUNCTION FAS(RKS)
C---CALC. THE HYPERGEOM. F-N F IN THE ASYMPTOTIC REGIME K*R >> 1.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 FAS,EIDC,Z1
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/ETA/ETA,ETA2,FI,CFI,SFI
      COMMON/COULPH/EIDC
      D1=DLOG(RKS)*ETA
      D2=ETA2/RKS
      Z1=DCMPLX(DCOS(D1),DSIN(D1))/EIDC
      FAS=DCMPLX(1.D0,-D2)*Z1
      FAS=(FAS-DCMPLX(DCOS(RKS),DSIN(RKS))*ETA/RKS/Z1)/ACHR
      RETURN
      END
      FUNCTION HC(XA)
C---HC = H-FUNCTION OF LANDAU-LIFSHITZ: h(x)=Re[psi(1-i/x)]+ln(x)
C-- psi(x) is the digamma function (the logarithmic derivative of
C-- the gamma function)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BN(15)
      DATA BN/.8333333333D-1,.8333333333D-2,.396825396825D-2,
     1        .4166666667D-2,.7575757576D-2,.2109279609D-1,
     2        .8333333333D-1,.4432598039D0 ,.305395433D1,
     3        .2645621212D2, .2814601449D3, .3607510546D4,
     4        .5482758333D5, .9749368235D6, .200526958D8/
      X=DABS(XA)
      IF(X.LT..33D0) GOTO 1
C      IF(X.GE.3.5D0) GO TO 2
      S=.0D0
      N=0
   3  N=N+1
      DS=1.D0/N/((N*X)**2+1)
      S=S+DS
      IF(DS.GT.0.1D-12) GOTO 3
C---Provides 7 digit accuracy
      HC=S-.5772156649D0+DLOG(X)
      RETURN
   2  HC=1.2D0/X**2+DLOG(X)-.5772156649D0
      RETURN
   1  X2=X*X
      XP=X2
      HC=0
      IMA=9
      IF(X.LT.0.1D0)IMA=3
      DO 4 I=1,IMA
      HC=HC+XP*BN(I)
   4  XP=X2*XP
      RETURN
      END
      SUBROUTINE GST(X,H)
C---CALC. THE CONFL. HYPERGEOM. F-N = CHH+i*SH
C-- AND THE COULOMB F-S B, P (CALLS SEQ OR SEQA).
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NSS/SR,NSS,LAW,ISW
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,ISPIN,MSPIN
      COMMON/SHH/SH,CHH
      COMMON/BP/B,P
  1   IF(ICH.EQ.0)GOTO 3
C      IF(NSS.EQ.1.OR.IHYPER.EQ.1)GOTO 2
      IF(IHYPER.EQ.1)GOTO 2
  3   SH=DSIN(H)
      CHH=DCOS(H)
      B=1
      IF(H.NE.0.D0)B=SH/H
      P=CHH
      RETURN
  2   CONTINUE
      IF(H.GT.10.D0)GOTO4
      CALL SEQ(X,H)
      SH=H*B
      CHH=P+B*X*(DLOG(DABS(X))+HPR)
      RETURN
  4   CALL SEQA(X,H)
      RETURN
      END
      FUNCTION FW(Y)
      IMPLICIT REAL*8 (A-H,O-Z)
C---FW = WEIGHT DUE TO THE INFLUENCE OF LAMBERT LAW.
      COMMON/ES/E,ES
      DATA PIH/1.570796327D0/
      EE=E*E
      A=(1-Y)*(1-EE)
      IF(A.EQ.0.D0) GOTO 3
      B=1-Y-(1+Y)*EE
      X=DSQRT(Y)*E
      IF(E.LT.0.D0) X=-X
      C=X*DSQRT(A-X*X)
      AS=DASIN(X/DSQRT(A))
    2 FW=B*(PIH-AS)+C
      RETURN
    3 FW=0.D0
      RETURN
      END
      SUBROUTINE SEAR
C---FINDS THE SQUARE WELL/BARRIER PARAMETERS A(1)=EB*RB^2 AND A(2)=RB
C-- BY SOLVING THE EQS. FUNCT(1)=1/FD, FUNCT(2)=RD
C-- WITH THE HELP OF THE FITTING ROUTINE FUMILI;
C-- FUNCT(1)=1/SCATT. AMPL.(K=0) AND FUNCT(2)=EFF. RADIUS ARE
C-- CALCULATED FROM THE SCATT. AMPL. CORRESP. TO THE SQ.-WELL POTENTIAL.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,JJ,MSPIN
      COMMON/A/A(100)/EXDA/EXDA(1500)/PL/PL(100)/NED/NED(2)
      COMMON/AL/AL(100)/AU/AU(100)
      COMMON/SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),SBKRB(10),SDKK(10)
      COMMON/FD/FD(10),FDS(10),RD(10)
      COMMON/DAK/DAK,IFUN
C---DAK= STEP TO CALCULATE dD/dk2 IN FUNCT
      DAK=.1D-0/DSQRT(DABS(FD(JJ)*RD(JJ)))
      IF(DAK*RD(JJ).GT.1.D0)DAK=.1D-0/RD(JJ)
      EXDA(1)=1/FD(JJ)
      EXDA(2)=.1D0
      EXDA(3)=1.D0
      EXDA(4)=RD(JJ)
      EXDA(5)=.15D0
      EXDA(6)=2.D0
      NED(1)=2
      NED(2)=3
      EPS=.1D-6
      DO 1 J=1,2
      IFUN=J
      IH1=ICH
      IF(J.EQ.1)ICH=0
      A(1)=EB(JJ)*RB(JJ)
      A(2)=RB(JJ)
      PL(1)=1.D0
      PL(2)=A(2)
      CALL FUMILI(S,2,2,1,999,EPS,AKA,ALA,0,MC)
      RB(JJ)=A(2)
      EB(JJ)=A(1)/RB(JJ)
      IH=ICH
C---CALC. FZ= SCATT. AMPL. AT K=0 SWITCHING OFF THE COULOMB INT.
      ICH=0
      FZ=.1973270536D0/G(.0D0)
      DH=RB(JJ)*.1973270536D0
      ICH=IH
      IF(J.EQ.1)PRINT 101,ICH,DH,EB(JJ),FZ
      IF(J.EQ.2)PRINT 102,ICH,DH,EB(JJ),FZ
 101  FORMAT(1X,4HICH=,I2,2X,3HRB=,E14.7,1X,2HFM,2X,
     C 3HEB=,E14.7,1X,5HGEV/C,2X,3HF0=,F9.3,1X,2HFM,2X,
     C 'ANALYT. CALC.')
 102  FORMAT(1X,4HICH=,I2,2X,3HRB=,E14.7,1X,2HFM,2X,
     C 3HEB=,E14.7,1X,5HGEV/C,2X,3HF0=,F9.3,1X,2HFM)
   1  ICH=IH1
      RETURN
      END
      FUNCTION FUNCT(X)
C---CALLED BY FITTING ROUTINE FUMILI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/A/A(100)/DAK/DAK,IFUN
      COMMON/SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),SBKRB(10),SDKK(10)
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,JJ,MSPIN
      RB(JJ)=A(2)
      EB(JJ)=A(1)/RB(JJ)
      IF(IFUN.EQ.1)GOTO 1
      GG=G(0.D0)
      IF(X.EQ.2.D0)GOTO 4
      FUNCT=GG
      RETURN
 4    AK=DAK
      AKS=AK*AK
      G1=G(AK)
      FUNCT=2*(G1-GG+HAC)/AKS
      RETURN
   1  SDK(JJ)=DSIN(A(1))
      CDK(JJ)=DCOS(A(1))
      FF=SDK(JJ)/A(1)-CDK(JJ)
      IF(FF.NE.0.D0)FF=CDK(JJ)/FF
      IF(X.EQ.1.D0)FUNCT=FF/A(2)
      IF(X.EQ.2.D0)FUNCT=A(2)*(1+FF/A(1)**2-FF**2/3)
      RETURN
      END
      FUNCTION G(AK)
C---CALLED BY F-N FUNCT(X) AND SUBR. SEAR
C---NOTE THAT SCATT. AMPL.= 1/CMPLX(G(AK),-AK*ACH)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),SBKRB(10),SDKK(10)
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,JJ,MSPIN
      COMMON/BP/B,P/DERIV/BPR,PPR/SHH/SH,CHH
      HAC=.0D0
      IF(ICH.EQ.0)GOTO 1
      XH=AC*AK
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
      ACH=ACP(XH)
      HAC=2*HCP/AC
   1  AKS=AK**2
      BK(JJ)=DSQRT(AKS+EB(JJ)**2) ! kappa=kp
      X=2*RB(JJ)/AC
      H=BK(JJ)*RB(JJ)             ! kp*d
      CALL GST(X,H)
      BRHO=B                      ! B(kp,d)
      SBKRB(JJ)=SH                ! kp*d*B(kp,d)
      CALL DERIW(X,H)
      BRHOP=BPR                   ! B'(kp,d)= dB(kp,r)/dln(r) at r=d
      H=AK*RB(JJ)
      CALL GST(X,H)
      CDK(JJ)=CHH                 ! ReG(k,d)
      BRHOS=B                     !  B(k,d)
      PRHOS=P                     !  P(k,d)
      SDK(JJ)=SH
      IF(ICH.EQ.0)GOTO 2
      SDK(JJ)=ACH*SH              ! ImG(k,d)
      IF(AK.EQ.0.D0.AND.AC.LT.0.D0)SDK(JJ)=3.14159*X*B
   2  SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK ! d*B(k,d)
      CALL DERIW(X,H)              ! PPR=P'(k,d)= dP(k,r)/dln(r) at r=d
      ZZ=PPR-PRHOS
      IF(ICH.EQ.1)ZZ=ZZ+X*(BRHOS+BPR*(DLOG(DABS(X))+HPR))
C   ZZ= P'(k,d)-P(k,d)+x*{B(k,d)+B'(k,d)*[ln!x!+2*C-1+h(k*ac)]}
      GG=(BRHOP*CDK(JJ)-BRHO*ZZ)/RB(JJ)
C   GG= [B'(kp,d)*ReG(k,d)-B(kp,d)*ZZ]/d
      G=GG/(BRHO*BPR-BRHOP*BRHOS)
C    G= GG/[B(kp,d)*B'(k,d)-B'(kp,d)*B(k,d)]
      RETURN
      END
      SUBROUTINE DERIW(X,H)
C---CALLED BY F-N G(AK)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ICH/HPR,AC,ACH,ACHR,HAC,ICH,IFI,IHYPER,JJ,MSPIN
      COMMON/SHH/SH,CHH/BP/B,P/DERIV/BPR,PPR
      HH=.1D-3
      CALL GST(X,H-HH)
      Q1=P
      B1=B
      CALL GST(X,H+HH)
      HHH=HH+HH
      BPR=H*(B-B1)/HHH
      PPR=H*(P-Q1)/HHH
      IF(ICH.EQ.0)RETURN
      CALL GST(X-HH,H)
      Q1=P
      B1=B
      CALL GST(X+HH,H)
      BPR=BPR+X*(B-B1)/HHH
      PPR=PPR+X*(P-Q1)/HHH
      RETURN
      END
      SUBROUTINE DSIMPS(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL%LIKE 1.E-17<,THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL %LIKE 1.E-10<,THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1(REPS1(0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*DABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END
      SUBROUTINE DSIMPT(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1(REPS1(0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*DABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END
      SUBROUTINE DSIMPU(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1(REPS1(0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
       IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*DABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
 
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END
      SUBROUTINE DSIMPV(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1(REPS1(0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*DABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END
      SUBROUTINE DSIMPW(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND -
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON-S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.D0
    4 X0=X
      IF((X0+4.D0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.D0
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.D0
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*DABS(F(K))
      DI1=(F(1)+4.D0*F(3)+F(5))*2.D0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.D0) 17,14,14
   17 H=2.D0*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.D0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.D0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
         END
      SUBROUTINE FUMILI (S,M,N1,N2,N3,EPS,AKAPPA,ALAMBD,IT,MC)
C-----ENTRY FOR CHISQ MINIMISATION
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510PL/PL(100)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      COMMON /A/A(100)/PL/PL0(100)/AU/AMX(100)/AL/AMN(100)/NED/NED(2)
      COMMON/SIGMA/SIGMA(100)/EXDA/EXDA(1500)
      COMMON/R/R(100)/ERROR/ERROR(500)/Z0/Z0(2485)
C-----10.*MAXIMUM RELATIVE PRECISION ON CDC 6000
      DATA RP/1.D-24/
      INDFLG(3)=0
 1    IF (IT.GE.0) WRITE(6,84)
C-      CALL VFILL(AMX,100,1.0D37)
C-      CALL VCOPYN(AMX,AMN,100)
      NN2=0
      N=M
      FIXFLG=0.D0
      ENDFLG=0.D0
      INDFLG(2)=0
      IFIX1=0
      FI=0.D0
      NN3=0
      DO 2 I=1,N
      R(I)=0.D0
      IF (EPS.GT.0.D0) SIGMA(I)=0.D0
 2    PL(I)=PL0(I)
C-----START NEW ITERATION
 3    NN1=1
      T1=1.D0
C-----REPEAT ITERATION WITH SMALLER STEP
 4    S=0.D0
      N0=0
      DO 7 I=1,N
      G(I)=0.D0
      IF (PL0(I)) 7,7,5
 5    N0=N0+1
      IF (PL(I)) 7,7,6
 6    PL0(I)=PL(I)
 7    CONTINUE
      NN0=N0*(N0+1)/2
      IF (NN0.LT.1) GO TO 9
      DO 8 I=1,NN0
      Z(I)=0.D0
 8    CONTINUE
 9    NA=M
      INDFLG(1)=0
C-----CALCULATE OBJECTIVE FUNCTION
      CALL SGZ (M,S)
      SP=RP*DABS(S)
      IF (NN0.LT.1) GO TO 11
      DO 10 I=1,NN0
      Z0(I)=Z(I)
 10   CONTINUE
 11   IF (NN3) 19,19,12
 12   IF (NN1-N1) 13,13,19
 13   T=2*(S-OLDS-GT)
      IF (INDFLG(1)) 16,14,16
 14   IF (DABS(S-OLDS).LE.SP.AND.-GT.LE.SP) GO TO 19
      IF (0.59D0*T+GT) 19,15,15
 15   T=-GT/T
      IF (T-0.25D0) 16,17,17
 16   T=0.25D0
 17   GT=GT*T
      T1=T1*T
      NN2=0
      DO 18 I=1,N
      IF (PL(I).LE.0.) GO TO 18
      A(I)=A(I)-DA(I)
      PL(I)=PL(I)*T
      DA(I)=DA(I)*T
      A(I)=A(I)+DA(I)
 18   CONTINUE
      NN1=NN1+1
      GO TO 4
C-----REMOVE CONTRIBUTION OF FIXED PARAMETERS FROM Z
   19 IF(INDFLG(1).EQ.0) GO TO 20
      ENDFLG=-4.D0
      GO TO 85
 20   K1=1
      K2=1
      I1=1
      DO 30 I=1,N
      IF (PL0(I)) 30,30,21
 21   IF (PL(I).EQ.0.) PL(I)=PL0(I)
      IF (PL(I)) 23,23,24
 22   PL(I)=0.D0
 23   K1=K1+I1
      GO TO 29
C- 24   IF (A(I).GE.AMX(I).AND.G(I).LT.0.D0) GO TO 22
C-      IF (A(I).LE.AMN(I).AND.G(I).GT.0.D0) GO TO 22
 24   IF (A(I).GE.AMX(I)) GO TO 22
      IF (A(I).LE.AMN(I)) GO TO 22
      DO 28 J=1,I
      IF (PL0(J)) 28,28,25
 25   IF (PL(J)) 27,27,26
 26   Z(K2)=Z0(K1)
      K2=K2+1
 27   K1=K1+1
 28   CONTINUE
 29   I1=I1+1
 30   CONTINUE
C-----INVERT Z
      I1=1
      L=I1
      DO 32 I=1,N
      IF (PL(I)) 32,32,31
 31   R(I)=Z(L)
      I1=I1+1
      L=L+I1
 32   CONTINUE
      N0=I1-1
      CALL MCONV (N0)
      IF (INDFLG(1)) 33,34,33
 33   INDFLG(1)=0
      INDFLG(2)=1
      GO TO 49
 34   CONTINUE
C-----CALCULATE THEORETICAL STEP TO MINIMUM
      I1=1
      DO 41 I=1,N
      DA(I)=0.D0
      IF (PL(I)) 41,41,35
 35   L1=1
      DO 40 L=1,N
      IF (PL(L)) 40,40,36
 36   IF (I1-L1) 37,37,38
 37   K=L1*(L1-1)/2+I1
      GO TO 39
 38   K=I1*(I1-1)/2+L1
 39   DA(I)=DA(I)-G(L)*Z(K)
      L1=L1+1
 40   CONTINUE
      I1=I1+1
 41   CONTINUE
C-----CHECK FOR PARAMETERS ON BOUNDARY
      AFIX=0.D0
      IFIX=0
      I1=1
      L=I1
      DO 47 I=1,N
      IF (PL(I)) 47,47,42
 42   SIGI=DSQRT(DABS(Z(L)))
      R(I)=R(I)*Z(L)
      IF (EPS) 44,44,43
 43   SIGMA(I)=SIGI
 44   IF ((A(I).LT.AMX(I).OR.DA(I).LE.0.).AND.(A(I).GT.AMN(I).OR.DA(I).G
     1E.0.D0)) GO TO 46
      AKAP=DABS(DA(I)/SIGI)
      IF (AKAP-AFIX) 46,46,45
 45   AFIX=AKAP
      IFIX=I
      IFIX1=I
 46   I1=I1+1
      L=L+I1
 47   CONTINUE
      IF (IFIX) 48,50,48
 48   PL(IFIX)=-1.D0
 49   FIXFLG=FIXFLG+1
      FI=0.D0
C-----REPEAT CALCULATION OF THEORETICAL STEP AFTER FIXING EACH PARAMETER
      GO TO 19
C-----CALCULATE STEP CORRECTION FACTOR
 50   ALAMBD=1.D0
      AKAPPA=0.D0
      IMAX=0
      DO 60 I=1,N
      IF (PL(I)) 60,60,51
 51   BM=AMX(I)-A(I)
      ABI=A(I)+PL(I)
      ABM=AMX(I)
      IF (DA(I)) 52,52,53
 52   BM=A(I)-AMN(I)
      ABI=A(I)-PL(I)
      ABM=AMN(I)
 53   BI=PL(I)
      IF (BI-BM) 55,55,54
 54   BI=BM
      ABI=ABM
 55   IF (DABS(DA(I))-BI) 58,58,56
 56   AL=DABS(BI/DA(I))
      IF (ALAMBD-AL) 58,58,57
 57   IMAX=I
      AIMAX=ABI
      ALAMBD=AL
 58   AKAP=DABS(DA(I)/SIGMA(I))
      IF (AKAP-AKAPPA) 60,60,59
 59   AKAPPA=AKAP
 60   CONTINUE
C-----CALCULATE NEW CORRECTED STEP
      GT=0.D0
      AMB=1.D18
      IF (ALAMBD) 62,62,61
 61   AMB=0.25D0/ALAMBD
 62   CONTINUE
      DO 67 I=1,N
      IF (PL(I)) 67,67,63
 63   IF (NN2-N2) 66,66,64
 64   IF (DABS(DA(I)/PL(I))-AMB) 66,65,65
 65   PL(I)=4*PL(I)
      T1=4.D0
 66   DA(I)=DA(I)*ALAMBD
      GT=GT+DA(I)*G(I)
 67   CONTINUE
C-----CHECK IF MINIMUM ATTAINED AND SET EXIT MODE
      IF (-GT.GT.SP.OR.T1.GE.1.D0.OR.ALAMBD.GE.1.D0) GO TO 68
      ENDFLG=-1.D0
 68   IF (ENDFLG) 85,69,69
 69   IF (AKAPPA-DABS(EPS)) 70,75,75
 70   IF (FIXFLG) 72,71,72
 71   ENDFLG=1.D0
      GO TO 85
 72   IF (ENDFLG) 85,77,73
 73   IF (IFIX1) 85,85,76
 74   IF (FI-FIXFLG) 76,76,77
 75   IF (FIXFLG) 74,76,74
 76   FI=FI+1
      ENDFLG=0.D0
 85   IF(ENDFLG.EQ.0.D0.AND.NN3.GE.N3) ENDFLG=-3.D0
      IF(ENDFLG.GT.0.D0.AND.INDFLG(2).GT.0) ENDFLG=-2.D0
      CALL MONITO (S,M,NN3,IT,EPS,GT,AKAPPA,ALAMBD)
      IF (ENDFLG) 83,79,83
C-----CHECK IF FIXING ON BOUND IS CORRECT
 77   ENDFLG=1.D0
      FIXFLG=0.D0
      IFIX1=0
      DO 78 I=1,M
 78   PL(I)=PL0(I)
      INDFLG(2)=0
      GO TO 19
C-----NEXT ITERATION
 79   ENDFLG=0.D0
      DO 80 I=1,N
      A(I)=A(I)+DA(I)
 80   CONTINUE
      IF (IMAX) 82,82,81
 81   A(IMAX)=AIMAX
 82   OLDS=S
      NN2=NN2+1
      NN3=NN3+1
      GO TO 3
 83   MC=ENDFLG
      RETURN
C-----ENTRY FOR MAXIMUM LIKLEHOOD
      ENTRY LIKELM(S,M,N1,N2,N3,EPS,AKAPPA,ALAMBD,IT,MC)
      INDFLG(3)=1
      GO TO 1
C
 84   FORMAT('1',1X,'FUNCTION MINIMISATION BY SUBROUTINE FUMILI/LIKE',
     +'LM'/'0',1X,'IN THE FOLLOWING PRINT-OUT'/
     +     '0',1X,'S = VALUE OF OBJECTIVE FUNCTION,',
     + 'EC = EXPECTED CHANGE IN S DURING NEXT ITERATION'/
     +   '0',1X,'KAPPA = ESTIMATED DISTANCE TO MINIMUM,  LAMBDA =',
     + 'STEP LENGTH MODIFIER'///)
      END
      FUNCTION FRSIN(RX)
 
      DOUBLE PRECISION DFRSIN,DFRCOS,DX,X,V,Y,R,S,H,ALFA,B0,B1,B2
      DOUBLE PRECISION A,B,C1,C2,C0
      DOUBLE PRECISION ZERO,ONE,TWO,FIVE,EIGHT
      CHARACTER NAME*6
 
      DIMENSION A(0:16),B(0:15),C1(0:25),C2(0:28)
 
      DATA ZERO /0D0/, ONE /1D0/, TWO /2D0/, FIVE /5D0/, EIGHT /8D0/
      DATA C0 /1.25331 41373 15500 3D0/
 
      DATA NA,NB,NC1,NC2 /16,15,25,28/
 
      DATA A( 0) / 0.76435 13866 41860 002D0/
      DATA A( 1) /-0.43135 54754 76601 793D0/
      DATA A( 2) / 0.43288 19997 97266 531D0/
      DATA A( 3) /-0.26973 31033 83871 110D0/
      DATA A( 4) / 0.08416 04532 08769 354D0/
      DATA A( 5) /-0.01546 52448 44613 820D0/
      DATA A( 6) / 0.00187 85542 34398 220D0/
      DATA A( 7) /-0.00016 26497 76188 875D0/
      DATA A( 8) / 0.00001 05739 76563 833D0/
      DATA A( 9) /-0.00000 05360 93398 892D0/
      DATA A(10) / 0.00000 00218 16584 549D0/
      DATA A(11) /-0.00000 00007 29016 212D0/
      DATA A(12) / 0.00000 00000 20373 325D0/
      DATA A(13) /-0.00000 00000 00483 440D0/
      DATA A(14) / 0.00000 00000 00009 865D0/
      DATA A(15) /-0.00000 00000 00000 175D0/
      DATA A(16) / 0.00000 00000 00000 003D0/
 
      DATA B( 0) / 0.63041 40431 45705 392D0/
      DATA B( 1) /-0.42344 51140 57053 335D0/
      DATA B( 2) / 0.37617 17264 33436 566D0/
      DATA B( 3) /-0.16249 48915 45095 674D0/
      DATA B( 4) / 0.03822 25577 86330 087D0/
      DATA B( 5) /-0.00564 56347 71321 909D0/
      DATA B( 6) / 0.00057 45495 19768 974D0/
      DATA B( 7) /-0.00004 28707 15321 020D0/
      DATA B( 8) / 0.00000 24512 07499 233D0/
      DATA B( 9) /-0.00000 01109 88418 409D0/
      DATA B(10) / 0.00000 00040 82497 317D0/
      DATA B(11) /-0.00000 00001 24498 302D0/
      DATA B(12) / 0.00000 00000 03200 484D0/
      DATA B(13) /-0.00000 00000 00070 324D0/
      DATA B(14) / 0.00000 00000 00001 336D0/
      DATA B(15) /-0.00000 00000 00000 022D0/
 
      DATA C1( 0) / 0.99056 04793 73497 549D0/
      DATA C1( 1) /-0.01218 35098 31478 997D0/
      DATA C1( 2) /-0.00248 27428 23113 060D0/
      DATA C1( 3) / 0.00026 60949 52647 247D0/
      DATA C1( 4) /-0.00000 10790 68987 406D0/
      DATA C1( 5) /-0.00000 48836 81753 933D0/
      DATA C1( 6) / 0.00000 09990 55266 368D0/
      DATA C1( 7) /-0.00000 00750 92717 372D0/
      DATA C1( 8) /-0.00000 00190 79487 573D0/
      DATA C1( 9) / 0.00000 00090 90797 293D0/
      DATA C1(10) /-0.00000 00019 66236 033D0/
      DATA C1(11) / 0.00000 00001 64772 911D0/
      DATA C1(12) / 0.00000 00000 63079 714D0/
      DATA C1(13) /-0.00000 00000 36432 219D0/
      DATA C1(14) / 0.00000 00000 10536 930D0/
      DATA C1(15) /-0.00000 00000 01716 438D0/
      DATA C1(16) /-0.00000 00000 00107 124D0/
      DATA C1(17) / 0.00000 00000 00204 099D0/
      DATA C1(18) /-0.00000 00000 00090 064D0/
      DATA C1(19) / 0.00000 00000 00025 506D0/
      DATA C1(20) /-0.00000 00000 00004 036D0/
      DATA C1(21) /-0.00000 00000 00000 570D0/
      DATA C1(22) / 0.00000 00000 00000 762D0/
      DATA C1(23) /-0.00000 00000 00000 363D0/
      DATA C1(24) / 0.00000 00000 00000 118D0/
      DATA C1(25) /-0.00000 00000 00000 025D0/
 
      DATA C2( 0) / 0.04655 77987 37516 4561D0/
      DATA C2( 1) / 0.04499 21302 01239 4140D0/
      DATA C2( 2) /-0.00175 42871 39651 4532D0/
      DATA C2( 3) /-0.00014 65340 02581 0678D0/
      DATA C2( 4) / 0.00003 91330 40863 0159D0/
      DATA C2( 5) /-0.00000 34932 28659 7731D0/
      DATA C2( 6) /-0.00000 03153 53003 2345D0/
      DATA C2( 7) / 0.00000 01876 58200 8529D0/
      DATA C2( 8) /-0.00000 00377 55280 4930D0/
      DATA C2( 9) / 0.00000 00026 65516 5010D0/
      DATA C2(10) / 0.00000 00010 88144 8122D0/
      DATA C2(11) /-0.00000 00005 35500 7671D0/
      DATA C2(12) / 0.00000 00001 31576 5447D0/
      DATA C2(13) /-0.00000 00000 15286 0881D0/
      DATA C2(14) /-0.00000 00000 03394 7646D0/
      DATA C2(15) / 0.00000 00000 02702 0267D0/
      DATA C2(16) /-0.00000 00000 00946 3142D0/
      DATA C2(17) / 0.00000 00000 00207 1565D0/
      DATA C2(18) /-0.00000 00000 00012 6931D0/
      DATA C2(19) /-0.00000 00000 00013 9756D0/
      DATA C2(20) / 0.00000 00000 00008 5929D0/
      DATA C2(21) /-0.00000 00000 00003 1070D0/
      DATA C2(22) / 0.00000 00000 00000 7515D0/
      DATA C2(23) /-0.00000 00000 00000 0648D0/
      DATA C2(24) /-0.00000 00000 00000 0522D0/
      DATA C2(25) / 0.00000 00000 00000 0386D0/
      DATA C2(26) /-0.00000 00000 00000 0165D0/
      DATA C2(27) / 0.00000 00000 00000 0050D0/
      DATA C2(28) /-0.00000 00000 00000 0009D0/
 
      X=RX
      NAME='FRSIN '
    9 V=ABS(X)
      IF(V .LT. EIGHT) THEN
       Y=V/EIGHT
       H=TWO*Y**2-ONE
       ALFA=H+H
       B1=ZERO
       B2=ZERO
       DO 4 I = NB,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    4  B1=B0
       H=SQRT(V)*Y*(B0-B2)
      ELSE
       R=ONE/V
       Y=FIVE*R
       H=Y+Y-ONE
       ALFA=H+H
       B1=ZERO
       B2=ZERO
       DO 5 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    5  B1=B0
       S=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 6 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    6  B1=B0
       H=C0-SQRT(R)*(S*COS(V)+(B0-H*B2)*SIN(V))
      END IF
      IF(X .LT. ZERO) H=-H
      IF(NAME .EQ. 'DFRSIN') THEN
      GO TO 10
      ELSE
       FRSIN=SNGL(H+(H-DBLE(SNGL(H))))
      END IF
      RETURN
 
      ENTRY FRCOS(RX)
      X=RX
      NAME='FRCOS '
    8 V=ABS(X)
      IF(V .LT. EIGHT) THEN
       Y=V/EIGHT
       H=TWO*Y**2-ONE
       ALFA=H+H
       B1=ZERO
       B2=ZERO
       DO 1 I = NA,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=SQRT(V)*(B0-H*B2)
      ELSE
       R=ONE/V
       Y=FIVE*R
       H=Y+Y-ONE
       ALFA=H+H
       B1=ZERO
       B2=ZERO
       DO 2 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 3 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    3  B1=B0
       H=C0-SQRT(R)*((B0-H*B2)*COS(V)-S*SIN(V))
      END IF
      IF(X .LT. ZERO) H=-H
      IF(NAME .EQ. 'DFRCOS') THEN
      GO TO 11
      ELSE
       FRCOS=SNGL(H+(H-DBLE(SNGL(H))))
      END IF
      RETURN
 
      ENTRY DFRSIN(DX)
      X=DX
      NAME='DFRSIN'
      GO TO 9
  10  DFRSIN = H
      RETURN
 
      ENTRY DFRCOS(DX)
      X=DX
      NAME='DFRCOS'
      GO TO 8
  11  DFRCOS = H
      END
      REAL FUNCTION BESJ0(RX)
      REAL RX
      CHARACTER*6 ENAME
      DOUBLE PRECISION X,Y,V,H,ALFA,ZERO,ONE,TWO,EIGHT,D
      DOUBLE PRECISION PI1,PI2,C1(0:14),C2(0:9),C3(0:10),B0,B1,B2,P,Q,R
      DOUBLE PRECISION DBESJ0,DX
 
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, EIGHT /8.0D0/
      DATA PI1 /0.79788 45608 0287D0/, PI2 /0.78539 81633 9745D0/
 
      DATA C1( 0) /+0.15772 79714 7489D0/
      DATA C1( 1) /-0.00872 34423 5285D0/
      DATA C1( 2) /+0.26517 86132 0334D0/
      DATA C1( 3) /-0.37009 49938 7265D0/
      DATA C1( 4) /+0.15806 71023 3210D0/
      DATA C1( 5) /-0.03489 37694 1141D0/
      DATA C1( 6) /+0.00481 91800 6947D0/
      DATA C1( 7) /-0.00046 06261 6621D0/
      DATA C1( 8) /+0.00003 24603 2882D0/
      DATA C1( 9) /-0.00000 17619 4691D0/
      DATA C1(10) /+0.00000 00760 8164D0/
      DATA C1(11) /-0.00000 00026 7925D0/
      DATA C1(12) /+0.00000 00000 7849D0/
      DATA C1(13) /-0.00000 00000 0194D0/
      DATA C1(14) /+0.00000 00000 0004D0/
 
      DATA C2( 0) /+0.99946 03493 4752D0/
      DATA C2( 1) /-0.00053 65220 4681D0/
      DATA C2( 2) /+0.00000 30751 8479D0/
      DATA C2( 3) /-0.00000 00517 0595D0/
      DATA C2( 4) /+0.00000 00016 3065D0/
      DATA C2( 5) /-0.00000 00000 7864D0/
      DATA C2( 6) /+0.00000 00000 0517D0/
      DATA C2( 7) /-0.00000 00000 0043D0/
      DATA C2( 8) /+0.00000 00000 0004D0/
      DATA C2( 9) /-0.00000 00000 0001D0/
 
      DATA C3( 0) /-0.01555 58546 05337D0/
      DATA C3( 1) /+0.00006 83851 99426D0/
      DATA C3( 2) /-0.00000 07414 49841D0/
      DATA C3( 3) /+0.00000 00179 72457D0/
      DATA C3( 4) /-0.00000 00007 27192D0/
      DATA C3( 5) /+0.00000 00000 42201D0/
      DATA C3( 6) /-0.00000 00000 03207D0/
      DATA C3( 7) /+0.00000 00000 00301D0/
      DATA C3( 8) /-0.00000 00000 00033D0/
      DATA C3( 9) /+0.00000 00000 00004D0/
      DATA C3(10) /-0.00000 00000 00001D0/
      ROUND(D)  =  SNGL(D+(D-DBLE(SNGL(D))))
 
      X=RX
      ENAME=' BESJ0'
      GOTO 9
      ENTRY DBESJ0(DX)
      ENAME='DBESJ0'
      X=DX
 
    9 V=ABS(X)
      IF(V .LT. EIGHT) THEN
       Y=V/EIGHT
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 1 I = 14,0,-1
       B0=C1(I)-ALFA*B1-B2
       B2=B1
    1  B1=B0
       B1=B0-H*B2
      ELSE
       R=ONE/V
       Y=EIGHT*R
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 2 I = 9,0,-1
       B0=C2(I)-ALFA*B1-B2
       B2=B1
    2  B1=B0
       P=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 3 I = 10,0,-1
       B0=C3(I)-ALFA*B1-B2
       B2=B1
    3  B1=B0
       Q=Y*(B0-H*B2)
       B0=V-PI2
       B1=PI1*SQRT(R)*(P*COS(B0)-Q*SIN(B0))
      ENDIF
      IF(ENAME .EQ. ' BESJ0')  THEN
         BESJ0=ROUND(B1)
      ELSE
         DBESJ0=B1
      ENDIF
      RETURN
      END
      REAL FUNCTION BESY0(RX)
      REAL RX,SX
      CHARACTER*6 ENAME
      LOGICAL MFLAG,RFLAG
      DOUBLE PRECISION X,Y,H,ALFA,ZERO,ONE,TWO,EIGHT,B0,B1,B2,P,Q,R,D
      DOUBLE PRECISION PI1,PI2,PI3,CE,C1(0:14),C2(0:14),C3(0:9),C4(0:10)
      DOUBLE PRECISION DBESY0,DX
 
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, EIGHT /8.0D0/
 
      DATA PI1 /0.79788 45608 0287D0/, PI2 /0.63661 97723 6758D0/
      DATA PI3 /0.78539 81633 9745D0/, CE /0.57721 56649 0153D0/
 
      DATA C1( 0) /+0.15772 79714 7489D0/
      DATA C1( 1) /-0.00872 34423 5285D0/
      DATA C1( 2) /+0.26517 86132 0334D0/
      DATA C1( 3) /-0.37009 49938 7265D0/
      DATA C1( 4) /+0.15806 71023 3210D0/
      DATA C1( 5) /-0.03489 37694 1141D0/
      DATA C1( 6) /+0.00481 91800 6947D0/
      DATA C1( 7) /-0.00046 06261 6621D0/
      DATA C1( 8) /+0.00003 24603 2882D0/
      DATA C1( 9) /-0.00000 17619 4691D0/
      DATA C1(10) /+0.00000 00760 8164D0/
      DATA C1(11) /-0.00000 00026 7925D0/
      DATA C1(12) /+0.00000 00000 7849D0/
      DATA C1(13) /-0.00000 00000 0194D0/
      DATA C1(14) /+0.00000 00000 0004D0/
 
      DATA C2( 0) /-0.02150 51114 4966D0/
      DATA C2( 1) /-0.27511 81330 4352D0/
      DATA C2( 2) /+0.19860 56347 0255D0/
      DATA C2( 3) /+0.23425 27461 0902D0/
      DATA C2( 4) /-0.16563 59817 1365D0/
      DATA C2( 5) /+0.04462 13795 4067D0/
      DATA C2( 6) /-0.00693 22862 9152D0/
      DATA C2( 7) /+0.00071 91174 0375D0/
      DATA C2( 8) /-0.00005 39250 7972D0/
      DATA C2( 9) /+0.00000 30764 9329D0/
      DATA C2(10) /-0.00000 01384 5718D0/
      DATA C2(11) /+0.00000 00050 5105D0/
      DATA C2(12) /-0.00000 00001 5258D0/
      DATA C2(13) /+0.00000 00000 0388D0/
      DATA C2(14) /-0.00000 00000 0008D0/
 
      DATA C3( 0) /+0.99946 03493 4752D0/
      DATA C3( 1) /-0.00053 65220 4681D0/
      DATA C3( 2) /+0.00000 30751 8479D0/
      DATA C3( 3) /-0.00000 00517 0595D0/
      DATA C3( 4) /+0.00000 00016 3065D0/
      DATA C3( 5) /-0.00000 00000 7864D0/
      DATA C3( 6) /+0.00000 00000 0517D0/
      DATA C3( 7) /-0.00000 00000 0043D0/
      DATA C3( 8) /+0.00000 00000 0004D0/
      DATA C3( 9) /-0.00000 00000 0001D0/
 
      DATA C4( 0) /-0.01555 58546 05337D0/
      DATA C4( 1) /+0.00006 83851 99426D0/
      DATA C4( 2) /-0.00000 07414 49841D0/
      DATA C4( 3) /+0.00000 00179 72457D0/
      DATA C4( 4) /-0.00000 00007 27192D0/
      DATA C4( 5) /+0.00000 00000 42201D0/
      DATA C4( 6) /-0.00000 00000 03207D0/
      DATA C4( 7) /+0.00000 00000 00301D0/
      DATA C4( 8) /-0.00000 00000 00033D0/
      DATA C4( 9) /+0.00000 00000 00004D0/
      DATA C4(10) /-0.00000 00000 00001D0/
      ROUND(D)  =  SNGL(D+(D-DBLE(SNGL(D))))
 
      ENAME=' BESY0'
      X=RX
      GOTO 9
      ENTRY DBESY0(DX)
      ENAME='DBESY0'
      X=DX
    9 IF(X .LE. ZERO) THEN
       CALL KERMTR('C312.1',LGFILE,MFLAG,RFLAG)
       IF(MFLAG) THEN
        SX=X
        IF(LGFILE .EQ. 0) THEN
         WRITE(*,100) ENAME,SX
        ELSE
         WRITE(LGFILE,100) ENAME,SX
        ENDIF
       ENDIF
       IF(.NOT.RFLAG) CALL ABEND
       IF(ENAME .EQ. ' BESY0')  THEN
          BESY0=ZERO
       ELSE
          DBESY0=ZERO
       ENDIF
       RETURN
      ENDIF
      IF(X .LT. EIGHT) THEN
       Y=X/EIGHT
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 1 I = 14,0,-1
       B0=C1(I)-ALFA*B1-B2
       B2=B1
    1  B1=B0
       P=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 2 I = 14,0,-1
       B0=C2(I)-ALFA*B1-B2
       B2=B1
    2  B1=B0
       B1=PI2*(CE+LOG(X/TWO))*P+B0-H*B2
      ELSE
       R=ONE/X
       Y=EIGHT*R
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 3 I = 9,0,-1
       B0=C3(I)-ALFA*B1-B2
       B2=B1
    3  B1=B0
       P=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 4 I = 10,0,-1
       B0=C4(I)-ALFA*B1-B2
       B2=B1
    4  B1=B0
       Q=Y*(B0-H*B2)
       B0=X-PI3
       B1=PI1*SQRT(R)*(Q*COS(B0)+P*SIN(B0))
      ENDIF
      IF(ENAME .EQ. ' BESY0')  THEN
         BESY0=ROUND(B1)
      ELSE
         DBESY0=B1
      ENDIF
      RETURN
  100 FORMAT(7X,A6,' ... NON-POSITIVE ARGUMENT X = ',E16.6)
      END
      subroutine abend
      write(*,*) 'SKVIIIIIIIIIK !!!!!'
      return
      end
          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
                    PARAMETER(KOUNTE  =  27)
          CHARACTER*6         ERCODE,   CODE(KOUNTE)
          LOGICAL             MFLAG,    RFLAG
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
          DATA      LOGF      /  0  /
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 100, 100 /
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 100, 100 /
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 100, 100 /
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 100, 100 /
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 100, 100 /
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 100, 100 /
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 100, 100 /
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 100, 100 /
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 100, 100 /
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 100, 100 /
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 100, 100 /
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 100, 100 /
          DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 100, 100 /
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 100, 100 /
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 100, 100 /
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 100, 100 /
          DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 100, 100 /
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 100, 100 /
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 100, 100 /
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 100, 100 /
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 100, 100 /
          DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 100,   0 /
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 100,   0 /
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 100,   0 /
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 100,   0 /
          DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 100, 100 /
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 100, 100 /
          LOGF  =  LGFILE
          IF(ERCODE .EQ. ' ')  THEN
             L  =  0
          ELSE
             DO 10  L = 1, 6
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
  10            CONTINUE
  12         CONTINUE
          ENDIF
          DO 14     I  =  1, KOUNTE
             IF(L .EQ. 0)  GOTO 13
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13         KNTM(I)  =  LIMITM
             KNTR(I)  =  LIMITR
  14         CONTINUE
          RETURN
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
          LOG  =  LOGF
          DO 20     I  =  1, KOUNTE
             IF(ERCODE .EQ. CODE(I))  GOTO 21
  20         CONTINUE
          WRITE(*,1000)  ERCODE
          CALL ABEND
          RETURN
  21      RFLAG  =  KNTR(I) .GE. 1
          IF(RFLAG  .AND.  (KNTR(I) .LT. 100))  KNTR(I)  =  KNTR(I) - 1
          MFLAG  =  KNTM(I) .GE. 1
          IF(MFLAG  .AND.  (KNTM(I) .LT. 100))  KNTM(I)  =  KNTM(I) - 1
          IF(.NOT. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1001)  CODE(I)
             ELSE
                WRITE(LOGF,1001)  CODE(I)
             ENDIF
          ENDIF
          IF(MFLAG .AND. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1002)  CODE(I)
             ELSE
                WRITE(LOGF,1002)  CODE(I)
             ENDIF
          ENDIF
          RETURN
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +           ' ERROR MONITOR. RUN ABORTED.')
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +           'CONDITION ',A6)
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
          END
      SUBROUTINE ARITHM (Y)
C-----ARITHM SETS Y AND CALCULATES DY/DA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      COMMON /A/A(100)/PL/PL0(100)/AU/AMX(100)/AL/AMN(100)/NED/NED(2)
      COMMON/SIGMA/SIGMA(100)/EXDA/EXDA(1500)
C-----10.*MAXIMUM RELATIVE PRECISION ON CDC 6000: 1.E-14
      DATA RP/1.D-24/
      Y=FUNCT(X)
      DO 6 I=1,NA
      DF(I)=0.D0
      IF (PL0(I)) 6,6,1
 1    AI=A(I)
      HI=0.01D0*PL0(I)
      PI=RP*DABS(AI)
      IF (HI.LE.PI) HI=PI
      A(I)=AI+HI
      IF (A(I)-AMX(I)) 5,5,2
 2    A(I)=AI-HI
      HI=-HI
      IF (A(I)-AMN(I)) 3,5,5
 3    A(I)=AMX(I)
      HI=AMX(I)-AI
      IF (AMN(I)-AI+HI) 4,5,5
 4    A(I)=AMN(I)
      HI=AMN(I)-AI
 5    DF(I)=(FUNCT(X)-Y)/HI
      A(I)=AI
 6    CONTINUE
      RETURN
      END
      SUBROUTINE MCONV (N)
C-----MCONV INVERTS THE POSITIVE DEFINITE PACKED SYMMETRIC MATRIX Z
C-----BY THE SQUARE-ROOT METHOD
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510PL/PL(100)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      COMMON/R/R(100)/ERROR/ERROR(500)/Z0/Z0(2485)
C-----MAXIMUM REAL NUMBER AND 10.*MAXIMUM RELATIVE PRECISION ON CDC6000
      DATA AM,RP / 1.0D75, 1.0D-24/
      IF (N.LT.1) RETURN
      APS=DSQRT(AM/N)
      AP=1/(APS*APS)
      IR=0
      DO 11 I=1,N
 1    IR=IR+1
      IF (PL(IR)) 1,1,2
 2    NI=I*(I-1)/2
      II=NI+I
      K=N+1
      IF (Z(II).LE.RP*DABS(R(IR)).OR.Z(II).LE.AP) GO TO 19
      Z(II)=1/DSQRT(Z(II))
      NL=II-1
 3    IF (NL-NI) 5,5,4
 4    Z(NL)=Z(NL)*Z(II)
      IF (DABS(Z(NL)).GE.APS) GO TO 16
      NL=NL-1
      GO TO 3
 5    IF (I-N) 6,12,12
 6    K=K-1
      NK=K*(K-1)/2
      NL=NK
      KK=NK+I
      D=Z(KK)*Z(II)
      C=D*Z(II)
      L=K
 7    LL=NK+L
      LI=NL+I
      Z(LL)=Z(LL)-Z(LI)*C
      L=L-1
      NL=NL-L
      IF (L-I) 9,9,7
 8    LL=NK+L
      LI=NI+L
      Z(LL)=Z(LL)-Z(LI)*D
 9    L=L-1
      IF (L) 10,10,8
 10   Z(KK)=-C
      IF (K-I-1) 11,11,6
 11   CONTINUE
 12   DO 14 I=1,N
      DO 14 K=I,N
      NL=K*(K-1)/2
      KI=NL+I
      D=0.D0
      DO 13 L=K,N
      LI=NL+I
      LK=NL+K
      D=D+Z(LI)*Z(LK)
      NL=NL+L
 13   CONTINUE
      KI=K*(K-1)/2+I
      Z(KI)=D
 14   CONTINUE
 15   RETURN
 16   K=I+NL-II
      IR=0
      DO 18 I=1,K
 17   IR=IR+1
      IF (PL(IR)) 17,17,18
 18   CONTINUE
 19   PL(IR)=-2.D0
      R(IR)=0.D0
      INDFLG(1)=IR
      GO TO 15
      END
      SUBROUTINE MONITO (S,M,NN3,IT,EPS,GT,AKAPPA,ALAMBD)
C-----ITERATION PROGRESS OUTPUT
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510PL/PL(100)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      COMMON /A/A(100)/PL/PL0(100)/AU/AMX(100)/AL/AMN(100)/NED/NED(2)
      COMMON/SIGMA/SIGMA(100)/EXDA/EXDA(1500)
      COMMON/R/R(100)/ERROR/ERROR(500)/Z0/Z0(2485)
      IF (IT) 11,3,1
 1    IF (NN3) 4,4,2
 2    IF (NM) 3,4,4
 3    IF (ENDFLG) 4,12,4
C-----PRINTER CARRIAGE CONTROL
 4    I1=6
      IF (M.GT.6) I1=5
      IF (M.GT.12) I1=4
      IF (M.GT.23) I1=1
C   NON-ANSI CARRIAGE CONTROL SUPPRESSED ON IBM
      IF(I1.GT.1) I1=0
      WRITE(6,19) I1,NN3,S,GT,AKAPPA,ALAMBD
      DO 10 I=1,M
      IF (PL0(I)) 9,9,5
 5    IF (PL (I)) 8,7,6
 6    CONTINUE
      WRITE(6,20) I,A(I),SIGMA(I),R(I)
      GO TO 10
 7    WRITE(6,21) I,A(I),SIGMA(I),R(I)
      GO TO 10
 8    IF (PL(I).GE.-1.D0) GO TO 7
      WRITE(6,22) I,A(I)
      GO TO 10
 9    WRITE(6,23) I,A(I)
 10   CONTINUE
 11   NM=-IT
 12   NM=NM+1
      IF (ENDFLG) 13,14,14
 13   I=-ENDFLG
      GO TO (15,16,17,18), I
 14   RETURN
 15   IF(IT.GE.0) WRITE(6,24)
      GO TO 14
 16   IF(IT.GE.0) WRITE(6,25)
      GO TO 14
 17   IF(IT.GE.0) WRITE(6,26)
      GO TO 14
 18   IF(IT.GE.0) WRITE(6,27)
      GO TO 14
C
 19   FORMAT(I1,'IT.',I3,', S= ',E11.5,', EC = ',E11.5,
     +', KAP= ',E10.5,', LAMB= ',E10.5//1X,2(6X,'PARAMETER'),
     +9X,'STANDARD',8X,'CORRELATION'/9X,'NUMBER',9X,'VALUE',11X,
     +'DEVIATION',9X,'FACTOR'/)
 20   FORMAT(9X,I3,4X,3(5X,E12.5))
 21   FORMAT(9X,I3,4X,3(5X,E12.5),' BOUNDARY')
 22   FORMAT(9X,I3,9X,E12.5,5X,'INFINITE ERROR ESTIMATED')
 23   FORMAT(9X,I3,9X,E12.5,5X,'THIS PARAMETER FIXED')
 24   FORMAT('0MINIMISATION TERMINATED: NO FURTHER DECREASE IN S',
     + ' IS OBTAINABLE'/'1')
 25   FORMAT('0MINIMISATION TERMINATED: INFINITE ERRORS ESTIMATED'/'1')
 26   FORMAT('0MINIMISATION TERMINATED: ITERATION LIMIT REACHED'/'1')
 27   FORMAT('0MINIMISATION TERMINATED: NEGATIVE OR ZERO Y ENCOUNTERED',
     + ' AS LOGARITHMIC ARGUMENT',/'1')
      END
      SUBROUTINE SCAL (M,ER)
C-----SUBROUTINE SCAL CALCULATES ER=VARIANCE OF THEORETICAL FUNCTION Y
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510PL/PL(100)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      ER=0.D0
      N=0
      DO 2 J=1,M
      IF (PL(J)) 2,2,1
 1    N=N+1
      DF(N)=DF(J)
 2    CONTINUE
      IF (N.EQ.0) RETURN
      DO 6 I=1,N
      DO 6 L=1,N
      IF (I-L) 3,3,4
 3    K=L*(L-1)/2+I
      GO TO 5
 4    K=I*(I-1)/2+L
 5    ER=ER+Z(K)*DF(I)*DF(L)
 6    CONTINUE
      RETURN
      END
      SUBROUTINE SGZ (M,S)
C-----SUBROUTINE SGZ SETS UP S(OBJECTIVE FUNCTION),G(GRADIENT OF S) AND
C-----Z(APPROXIMATE COVARIANCE MATRIX)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/D510SI/Z(2485),G(100),DA(100),DF(100),ENDFLG,NA,INDFLG(5)
      COMMON/X/X(10)
      COMMON /A/A(100)/PL/PL0(100)/AU/AMX(100)/AL/AMN(100)/NED/NED(2)
      COMMON/SIGMA/SIGMA(100)/EXDA/EXDA(1500)
      K=NED(1)
      K2=1
      DO 12 L1=1,K
      K1=K2
      NX=NED(2)-2
      IF (INDFLG(3)) 1,2,1
 1    NX=NED(2)
      K1=K1-2
 2    CONTINUE
      DO 3 I=1,NX
      KI=K1+1+I
      X(I)=EXDA(KI)
 3    CONTINUE
      CALL ARITHM (Y)
      IF (INDFLG(3)) 4,6,4
 4    IF (Y) 13,13,5
C-----MAXIMUM LIKLEHOOD
 5    S=S-LOG(Y)
      Y=-Y
      SIG=Y
      GO TO 7
C-----CHI SQUARED
 6    SIG=EXDA(K2+1)
      Y=Y-EXDA(K1)
      S=S+((Y/SIG)**2)/2
 7    CONTINUE
      N=0
      DO 9 J=1,M
      IF (PL0(J)) 9,9,8
 8    N=N+1
      DF(N)=DF(J)/SIG
      G(J)=G(J)+DF(N)*(Y/SIG)
 9    CONTINUE
      L=1
      IF (N.LT.1) GO TO 11
      DO 10 I=1,N
      DO 10 J=1,I
      Z(L)=Z(L)+DF(I)*DF(J)
      L=L+1
 10   CONTINUE
 11   K2=K2+NED(2)
 12   CONTINUE
      RETURN
C-----  -VE OR ZERO Y IN MAXIMUM LIKLEHOOD
 13   INDFLG(1)=1
      S=1.D10
      RETURN
      END
      SUBROUTINE  VCOPYN (A,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C CERN PROGLIB# F121    VCOPYN          .VERSION KERNFOR  3.15  820113
C ORIG. 01/04/73
C
      DIMENSION A(*), X(*)
C
      IF(N.LE.0)  RETURN
      DO 10  I=1,N
   10 X(I)= -A(I)
      RETURN
      END
      SUBROUTINE VFILL (A,N,STUFF)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C CERN PROGLIB# F121    VFILL           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= STUFF
C
      RETURN
      END
