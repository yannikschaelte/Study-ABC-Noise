!***********************************************************************************
! This is HH_master.f95
! Written by: Joshua H. Goldwyn
! On: September 1, 2010
! This is a collection of subroutines and functions that define different versions of the stochastic Hodgkin-Huxley Model
! Called by the executable file: HH_run
!***********************************************************************************



!***********************************************************************************
SUBROUTINE STR2INT(string_in,length,int_out)
!
! Converts string of length "length" to integer number
!
CHARACTER(length), INTENT(IN) :: string_in
REAL temp
INTEGER, INTENT(OUT) :: int_out

read(string_in,*) temp
int_out = INT(temp)

END SUBROUTINE STR2INT
!***********************************************************************************


!***********************************************************************************
SUBROUTINE STR2DBLE(string_in,length,dble_out)
!
! Converts string of length "length" to double number
!
CHARACTER(length), INTENT(IN) :: string_in
DOUBLE PRECISION, INTENT(OUT) :: dble_out

read(string_in,*) dble_out

END SUBROUTINE STR2DBLE
!***********************************************************************************


!***********************************************************************************
DOUBLE PRECISION FUNCTION RANDNORM()

! Random Number Generator for Standard Normal Distribution
! Box-Muller Algorithm from Numerical Recipes in C, pg 217

! Mersenne Twister Random Number Generator
USE mtmod

IMPLICIT NONE

! Declare Local Variables
  DOUBLE PRECISION :: x1,x2,r
 
  ! Generate Uniform random in square, until find something in the unit disc
  r=2.
  DO WHILE (r .GE. 1)
    x1 = 2.0*grnd()-1.
    x2 = 2.0*grnd()-1.
    r = x1*x1 + x2*x2
  END DO
  
  RANDNORM = x1*SQRT(-2.0*LOG(r)/r)

END FUNCTION RANDNORM
!***********************************************************************************


!***********************************************************************************
SUBROUTINE NumberChannel(A,NNa,NK)
! Compute Number of Na and K Channels from area of membrane

! Read in parameter values from modules
USE ParameterModule

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: A ! Membrane Area (mu m^2)

! Declare Output Variables
INTEGER, INTENT(OUT) :: NNa, NK  !Total Number of Channels

NNa = NINT(DNa * A)
NK = NINT(DK * A)

END SUBROUTINE NumberChannel
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION alpham(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience
  
IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
alpham = 0.1 * (25.-V)/ (EXP((25.-V)/10.)-1.)

END FUNCTION alpham
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION betam(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
betam = 4. * EXP(-V/18.)

END FUNCTION betam
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION alphah(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
alphah = 0.07 * exp(-V/20.)

END FUNCTION alphah
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION betah(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
betah = 1. / (exp((30.-V)/10.)+1.)

END FUNCTION betah
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION alphan(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
alphan = 0.01 * (10.-V) / (EXP((10.-V)/10.)-1.)

END FUNCTION alphan
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION betan(V)
! See Chapter 2 in Izhikevich, Dynamical Systems in Neuroscience

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
 
betan = 0.125 * exp(-V/80.)

END FUNCTION betan
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION dm(V,m)
! Deterministic dynamics of gating variable m

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham
DOUBLE PRECISION betam

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: m
 
dm = alpham(V)*(1-m) - betam(V)*m

END FUNCTION dm
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION dh(V,h)
! Deterministic dynamics of gating variable h

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphah
DOUBLE PRECISION betah

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: h
 
dh = alphah(V)*(1-h) - betah(V)*h

END FUNCTION dh
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION dn(V,n)
! Deterministic dynamics of gating variable n

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphan
DOUBLE PRECISION betan

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: n
 
dn = alphan(V)*(1-n) - betan(V)*n

END FUNCTION dn
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION dV(V,M,h,N,I)
! Dynamics of voltage V

! Read in parameter values from modules
USE ParameterModule

IMPLICIT NONE

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: M
DOUBLE PRECISION, INTENT(IN) :: h
DOUBLE PRECISION, DIMENSION(4), INTENT(IN) :: N
DOUBLE PRECISION, INTENT(IN) :: I ! Applied Current
 
dV = (I - gK*N(1)*N(2)*N(3)*N(4)*(V-EK) - gNa*M(1)*M(2)*M(3)*h*(V-ENa) - gL*(V-EL)) / C

END FUNCTION dV
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION mnoise(V,m,NNa)
! Standard Deviation of noise in Subunit-based SDE models (m subunit)
  
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham
DOUBLE PRECISION betam

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: m
INTEGER, INTENT(IN) :: NNa
 
mnoise = SQRT((alpham(V)*(1.-m) + betam(V)*m) / NNa)

END FUNCTION mnoise
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION hnoise(V,h,NNa)
! Standard Deviation of noise in Subunit-based SDE models (h subunit)
  
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphah
DOUBLE PRECISION betah

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: h
INTEGER, INTENT(IN) :: NNa
 
hnoise = SQRT((alphah(V)*(1.-h) + betah(V)*h) / NNa)

END FUNCTION hnoise
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION nnoise(V,n,NK)
! Standard Deviation of noise in Subunit-based SDE models (n subunit)
  
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphan
DOUBLE PRECISION betan

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, INTENT(IN) :: n
INTEGER, INTENT(IN) :: NK
 
nnoise = SQRT((alphan(V)*(1.-n) + betan(V)*n) / NK)

END FUNCTION nnoise
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION mnoise_modify(V,NNa)
! Standard deviation of noise for m subunit for modified subunit-based SDE model

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, alphah
DOUBLE PRECISION betam, betah

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
INTEGER, INTENT(IN) :: NNa
 
! Declare Local Variables
DOUBLE PRECISION muh, mum, varh
DOUBLE PRECISION a,b,c

mum = alpham(V) / (alpham(V) + betam(V))
muh = alphah(V) / (alphah(V) + betah(V))
varh = (alphah(V) * betah(V)) / (DBLE(NNa) * (alphah(V) + betah(V))**2.)  

a = 36.*mum**2*muh**2.
b = 9.*mum**4. *muh**2. + 15.*mum**4. * varh
c = mum**6.*varh - (mum**3.*muh*(1.-mum**3.*muh))/DBLE(NNa)

! Order 1/N^2
mnoise_modify = SQRT(2.*(alpham(V)+betam(V)) * (-b + SQRT(b**2.-4.*a*c)) / (2.*a) )

END FUNCTION mnoise_modify
!***********************************************************************************

!***********************************************************************************
DOUBLE PRECISION FUNCTION nnoise_modify(V,NK)
! Standard deviation of noise for n subunit for modified subunit-based SDE model
 
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphan
DOUBLE PRECISION betan

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
INTEGER, INTENT(IN) :: NK
 
! Declare Local Variables
DOUBLE PRECISION mun
DOUBLE PRECISION a,b,c
mun = alphan(V) / (alphan(V) + betan(V))

! O(1/N^2)
a = 168.*mun**4.
b = 16.*mun**6.
c = -(mun**4. * (1.-mun**4.)) / DBLE(NK)
nnoise_modify = SQRT(2.*(alphan(V)+betan(V)) * (-b + SQRT(b**2.-4.*a*c)) / (2.*a) )

END FUNCTION nnoise_modify
!***********************************************************************************

!***********************************************************************************
SUBROUTINE dXNa(V,X,dX)
! Deterministic part of Na channel for channel-based SDE

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, alphah
DOUBLE PRECISION betam, betah

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, DIMENSION(1:7), INTENT(IN) :: X
 
! Declare Local Variable
DOUBLE PRECISION X0

! Declare Output Variables
DOUBLE PRECISION, DIMENSION(7), INTENT(OUT) :: dX

X0 = 1.-SUM(X)
dX(1) = 3.*alpham(V)*X0 + (-2.*alpham(V)-betam(V)-alphah(V)) *X(1) + (2.*betam(V))*X(2) + ( betah(V))*X(5)
dX(2) = 2.*alpham(V)*X(1)  + (-alpham(V)-2*betam(V)-alphah(V)) * X(2) + 3*betam(V)*X(3) + betah(V)*X(6)
dX(3) = alpham(V)*X(2) + (-3*betam(V)-alphah(V)) * X(3) + betah(V)*X(7)
dX(4) = alphah(V)*X0+( -3*alpham(V)-betah(V))*X(4) + (betam(V))*X(5)
dX(5) = alphah(V)*X(1) + 3*alpham(V)*X(4) + (-2*alpham(V)-betam(V)-betah(V))*X(5)+ 2*betam(V)*X(6)
dX(6) = alphah(V)*X(2)+ 2*alpham(V)*X(5) + (-alpham(V)-2*betam(V)-betah(V))*X(6)+3*betam(V)*X(7)
dX(7) = alphah(V)*X(3) + alpham(V) *X(6) + (-3*betam(V)-betah(V))*X(7)

END SUBROUTINE dXNa
!***********************************************************************************   

!***********************************************************************************
SUBROUTINE dXK(V,X,dX)
! Deterministic part of K channel for channel-based SDE

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alphan
DOUBLE PRECISION betan

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
DOUBLE PRECISION, DIMENSION(1:4), INTENT(IN) :: X

! Declare Local Variable
DOUBLE PRECISION X0
 
! Declare Output Variables
DOUBLE PRECISION, DIMENSION(4), INTENT(OUT) :: dX

X0 = 1.-SUM(X)

dX(1) = 4.*alphan(V)*X0 + (-3.*alphan(V)-betan(V))*X(1) + 2.*betan(V)*X(2)
dX(2) = 3.*alphan(V)*X(1) + (-2*alphan(V)-2*betan(V))*X(2) +  3.*betan(V)*X(3)
dX(3) =  2.*alphan(V)*X(2) + (-alphan(V)-3*betan(V))*X(3)+ 4.*betan(V)*X(4)
dX(4) =  alphan(V)*X(3) -4.*betan(V)*X(4)

END SUBROUTINE dXK
!***********************************************************************************   

!***********************************************************************************
SUBROUTINE Nanoisematrix_eq(V, N, sqrtD)
! Square Root of Diffusion matrix for Na channel in channel-based SDE
! Using equilibrium noise approximation - state variables yij are replaced by the equilibrium mean values

! Read in parameter values from modules
USE ParameterModule
  
IMPLICIT NONE

! Declare LAPACK subroutine
EXTERNAL DSYEV

! Declare Functions
DOUBLE PRECISION alphah, alpham
DOUBLE PRECISION betah, betam

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
INTEGER, INTENT(IN) :: N ! number of Na channels

! Declare Local Variables
DOUBLE PRECISION y10, y20, y30, y01, y11, y21, y31, y00
DOUBLE PRECISION, DIMENSION(7,7) :: D  ! diffusion matrix
INTEGER i,j
INTEGER NN,LDA,LWMAX, INFO, LWORK
DOUBLE PRECISION, DIMENSION(7) :: EIGVAL
DOUBLE PRECISION, DIMENSION(1000) :: WORK
DOUBLE PRECISION, DIMENSION(7,7) :: temp

! Declare Output Variables
DOUBLE PRECISION, DIMENSION(7,7), INTENT(OUT) :: sqrtD ! stochastic part of sde

y10 = 3.*betam(V)**2. * alpham(V) * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y20 = 3.*betam(V) * alpham(V)**2. * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y30 = alpham(V)**3. * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y01 = betam(V)**3. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y11 = 3.*betam(V)**2. * alpham(V) * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y21 = 3.*betam(V) * alpham(V)**2. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
y31 = alpham(V)**3. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )

y00 = 1. - (y10+y20+y30+y01+y11+y21+y31)

!print *, y00, y10, y20, y30, y01, y11, y21, y31
D(1,1) = ((betam(V)+2*alpham(V))*y10 + 2*betam(V)*y20 + 3*alpham(V)*y00 + alphah(V)*y10 + betah(V)*y11) / (2.*DBLE(N))
D(1,2) = -(2*alpham(V)*y10 + 2*betam(V)*y20) / (2.*DBLE(N))
D(1,3) = 0
D(1,4) = 0
D(1,5) = -(alphah(V)*y10 + betah(V)*y11) / (2.*DBLE(N))
D(1,6) = 0
D(1,7) = 0

D(2,1) = D(1,2)
D(2,2) = ((2*betam(V)+alpham(V))*y20 + 3*betam(V)*y30 + 2*alpham(V)*y10 + alphah(V)*y20 + betah(V)*y21) / (2.*DBLE(N))
D(2,3) = -(alpham(V)*y20+3*betam(V)*y30) / (2.*DBLE(N))
D(2,4) = 0
D(2,5) = 0
D(2,6) = -(alphah(V)*y20+betah(V)*y21) / (2.*DBLE(N))
D(2,7) = 0

D(3,1) = D(1,3)
D(3,2) = D(2,3)
D(3,3) = (3*betam(V)*y30 + alpham(V)*y20 + alphah(V)*y30 + betah(V)*y31) / (2.*DBLE(N))
D(3,4) = 0
D(3,5) = 0
D(3,6) = 0
D(3,7) = -(alphah(V)*y30 + betah(V)*y31) / (2.*DBLE(N))

D(4,1) = D(1,4)
D(4,2) = D(2,4)
D(4,3) = D(3,4)
D(4,4) = (3*alpham(V)*y01 + betam(V)*y11 + betah(V)*y01 + alphah(V)*y00) / (2.*DBLE(N))
D(4,5) = -(3*alpham(V)*y01 + betam(V)*y11) / (2.*DBLE(N))
D(4,6) = 0
D(4,7) = 0

D(5,1) = D(1,5)
D(5,2) = D(2,5)
D(5,3) = D(3,5)
D(5,4) = D(4,5)
D(5,5) = ((betam(V)+2*alpham(V))*y11 + 2*betam(V)*y21 + 3*alpham(V)*y01 + betah(V)*y11 + alphah(V)*y10) / (2.*DBLE(N))
D(5,6) = -(2*alpham(V)*y11+2*betam(V)*y21)/(2.*DBLE(N))
D(5,7) = 0

D(6,1) = D(1,6)
D(6,2) = D(2,6)
D(6,3) = D(3,6)
D(6,4) = D(4,6)
D(6,5) = D(5,6)
D(6,6) = ((2*betam(V)+alpham(V))*y21+3*betam(V)*y31+2*alpham(V)*y11+betah(V)*y21+alphah(V)*y20) / (2.*DBLE(N))
D(6,7) = -(alpham(V)*y21+3*betam(V)*y31) / (2.*DBLE(N))

D(7,1) = D(1,7)
D(7,2) = D(2,7)
D(7,3) = D(3,7)
D(7,4) = D(4,7)
D(7,5) = D(5,7)
D(7,6) = D(6,7)
D(7,7) = (3*betam(V)*y31 + alpham(V)*y21 + betah(V)*y31 + alphah(V)*y30) / (2.*DBLE(N));

! Factorize D
NN = 7
LDA = 7
LWMAX = 1000
LWORK = -1 

CALL DSYEV( 'Vectors', 'Upper', NN, D, LDA, EIGVAL, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DSYEV( 'Vectors', 'Upper', NN, D, LDA, EIGVAL, WORK, LWORK, INFO )

! Compute Square Root Matrix
DO i=1,7
  DO j=1,7  
    temp(j,i) = D(j,i)*SQRT(EIGVAL(i))
  END DO
END DO

DO i=1,7
  DO j=1,7
    sqrtD(i,j) = DOT_PRODUCT(temp(i,:),D(j,:))
  END DO
END DO

END SUBROUTINE Nanoisematrix_eq
!***********************************************************************************

!***********************************************************************************
SUBROUTINE Knoisematrix_eq(V, N, sqrtD)
! Square Root of Diffusion matrix for K channel in channel-based SDE
! Using equilibrium noise approximation - state variables xi are replaced by the equilibrium mean values

! Read in parameter values from modules
USE ParameterModule
  
IMPLICIT NONE

! Declare LAPACK subroutine
EXTERNAL DSYEV

! Declare Functions
DOUBLE PRECISION alphan
DOUBLE PRECISION betan

! Declare Input Variables
DOUBLE PRECISION, INTENT(IN) :: V
INTEGER, INTENT(IN) :: N ! number of K channels

! Declare Local Variables
DOUBLE PRECISION, DIMENSION(4,4) :: D  ! diffusion matrix
INTEGER i,j
INTEGER NN,LDA,LWMAX, INFO, LWORK
DOUBLE PRECISION, DIMENSION(4) :: EIGVAL
DOUBLE PRECISION, DIMENSION(1000) :: WORK
DOUBLE PRECISION, DIMENSION(4,4) :: temp
DOUBLE PRECISION, DIMENSION(4) :: X ! state of K channels

! Declare Output Variables
DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: sqrtD ! stochastic part of sde

X(1) =  4.*alphan(V)*betan(V)**3./(alphan(V)+betan(V))**4.
X(2) =  6.*alphan(V)**2.*betan(V)**2./(alphan(V)+betan(V))**4.
X(3) =  4.*alphan(V)**3.*betan(V)/(alphan(V)+betan(V))**4.
X(4) =  alphan(V)**4./(alphan(V)+betan(V))**4.

D(1,1) = ((-alphan(V)+betan(V))*X(1)+ (2*betan(V)-4*alphan(V))*X(2) - &
4*alphan(V)*X(3) -4*alphan(V)*X(4) + 4*alphan(V)) / (2.*DBLE(N))
D(1,2) = -(2*betan(V)*X(2) + 3*alphan(V)*x(1)) / (2.*DBLE(N))
D(1,3) = 0
D(1,4) = 0
D(2,1) = D(1,2)
D(2,2) = (3*alphan(V)*X(1) + (2*alphan(V)+2*betan(V))*X(2) + 3*betan(V)*X(3)) /(2.*DBLE(N))
D(2,3) = -(3*betan(V)*X(3) + 2*alphan(V)*X(2)) / (2*N)
D(2,4) = 0
D(3,1) = 0
D(3,2) = D(2,3)
D(3,3) = (2*alphan(V)*X(2) + (alphan(V)+3*betan(V))*X(3) +4*betan(V)*X(4)) / (2.*DBLE(N))
D(3,4) = -(4*betan(V)*X(4) + alphan(V)*X(3)) / (2.*DBLE(N))
D(4,1) = 0
D(4,2) = 0
D(4,3) = D(3,4)
D(4,4) = (alphan(V)*X(3) + 4*betan(V)*X(4))/ (2.*DBLE(N));

! Factorize D
NN = 4
LDA = 4
LWMAX = 1000
LWORK = -1 

CALL DSYEV( 'Vectors', 'Upper', NN, D, LDA, EIGVAL, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DSYEV( 'Vectors', 'Upper', NN, D, LDA, EIGVAL, WORK, LWORK, INFO )

! Compute Square Root Matrix
DO i=1,4
  DO j=1,4  
    temp(j,i) = D(j,i)*SQRT(EIGVAL(i))
  END DO
END DO

DO i=1,4
  DO j=1,4
    sqrtD(i,j) = DOT_PRODUCT(temp(i,:),D(j,:))
  END DO
END DO

END SUBROUTINE Knoisematrix_eq
!***********************************************************************************


!***********************************************************************************
SUBROUTINE HHode(StimParam,nt,dt,maxnisi,PrintOut,seed)
! Forward Euler Solver for HH ODE

! Read in parameter values from modules
USE ParameterModule
USE mtmod

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION dm, dh, dn
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: maxnisi
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam
INTEGER, INTENT(IN) :: seed, PrintOut

! Declare Local Variables
INTEGER i
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike
DOUBLE PRECISION Iin
DOUBLE PRECISION V,m,h,n, Vold

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF

m = alpham(V) / (alpham(V) + betam(V))
h = alphah(V) / (alphah(V) + betah(V))
n = alphan(V) / (alphan(V) + betan(V))

nspike = 0
timelastspike = -9999.
tsincespike = 0.
stop = 0
i = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, m**3*h, n**4
END IF
 
DO WHILE (stop .EQ. 0)

  i = i+1
  t = DBLE(i)*dt
  tsincespike = tsincespike+dt

  Vold = V
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,[m,m,m],h,[n,n,n,n],Iin)
  END IF

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  m = m + dt*dm(Vold,m) 
  h = h + dt*dh(Vold,h) 
  n = n + dt*dn(Vold,n)

  ! Detect spikes and compute ISI
  IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
    nspike = nspike + 1
  
  IF (PrintOut .EQ. 2) THEN
    ! Write ISI as output
    IF (nspike .GT. 1) THEN
      print *, t - timelastspike
    END IF
  END IF

    timelastspike = t
    tsincespike = 0.
  END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, m**3*h, n**4  
  END IF

  ! Stopping Conditions
  IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
    stop = 1 
  ELSE IF (i .EQ. nt) THEN ! reached last time bin
    stop = 1
  ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
    stop = 1
  END IF

END DO

END SUBROUTINE HHode
!***********************************************************************************


!***********************************************************************************
SUBROUTINE HHsde_subunit_identical(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! Subunit SDE with variables raised to powers (Fox and Lu subunit approximation)

! Read in parameter values from modules
USE ParameterModule
USE mtmod
 
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION dm, dh, dn
DOUBLE PRECISION mnoise, nnoise, hnoise
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: maxnisi
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam
INTEGER, INTENT(IN) :: seed, PrintOut

! Declare Local Variables
INTEGER i
DOUBLE PRECISION V,m,h,n, Vold
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike
DOUBLE PRECISION Iin

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)

IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF
m = alpham(V) / (alpham(V) + betam(V))
h = alphah(V) / (alphah(V) + betah(V))
n = alphan(V) / (alphan(V) + betan(V))

nspike = 0
timelastspike = -9999.
tsincespike =0.
stop = 0
i = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, m**3*h, n**4
END IF


DO WHILE (stop .EQ. 0)
  i = i+1
  t = DBLE(i)*dt
  tsincespike = tsincespike+dt
  Vold = V
 
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,[m,m,m],h,[n,n,n,n],Iin)
  END IF

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  m = MIN(1., MAX(0., m + dt*dm(Vold,m) + SQRT(dt)*mnoise(Vold,m,NNa)*RANDNORM()))
  h = MIN(1., MAX(0., h + dt*dh(Vold,h) + SQRT(dt)*hnoise(Vold,h,NNa)*RANDNORM() ))
  n = MIN(1., MAX(0., n + dt*dn(Vold,n) + SQRT(dt)*nnoise(Vold,n,NK)*RANDNORM()))

  ! Detect spikes and compute ISI
  IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
    nspike = nspike + 1

  IF (PrintOut .EQ. 2) THEN
    ! Write ISI as output
    IF (nspike .GT. 1) THEN
      print *, t - timelastspike
    END IF
  END IF

    timelastspike = t
    tsincespike = 0.
  END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, m**3*h, n**4
  END IF 

  ! Stopping Conditions
  IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
    stop = 1 
  ELSE IF (i .EQ. nt) THEN ! reached last time bin
    stop = 1
  ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
    stop = 1
  END IF

END DO

END SUBROUTINE HHsde_subunit_identical
!***********************************************************************************


!***********************************************************************************
SUBROUTINE HHsde_subunit_independent(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! Subunit SDE using products of independent variables (Shuai and Jung subunit approximation)

! Read in parameter values from modules
USE ParameterModule
USE mtmod
  
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION dm, dh, dn
DOUBLE PRECISION mnoise, nnoise, hnoise
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: maxnisi
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam
INTEGER, INTENT(IN) :: seed, PrintOut

! Declare Local Variables
INTEGER i
DOUBLE PRECISION V, Vold, h
DOUBLE PRECISION, DIMENSION(1:3) :: M
DOUBLE PRECISION, DIMENSION(1:4) :: N
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike
DOUBLE PRECISION Iin

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF
M(:) = alpham(V) / (alpham(V) + betam(V))
h = alphah(V) / (alphah(V) + betah(V))
N(:) = alphan(V) / (alphan(V) + betan(V))

nspike = 0
timelastspike = -9999.
tsincespike = 0.
stop = 0
i = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, M(1)*M(2)*M(3)*h, N(1)*N(2)*N(3)*N(4)
END IF

DO WHILE (stop .EQ. 0)
  i = i+1
  t = DBLE(i)*dt
  tsincespike = tsincespike+dt

  Vold = V
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,M(:),h,N(:),Iin)
  END IF

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  M(1) = MIN(1., MAX(0., M(1) + dt*dm(Vold,M(1)) + SQRT(dt)*mnoise(Vold,M(1),NNa)*RANDNORM() ))
  M(2) = MIN(1., MAX(0., M(2) + dt*dm(Vold,M(2)) + SQRT(dt)*mnoise(Vold,M(2),NNa)*RANDNORM() ))
  M(3) = MIN(1., MAX(0., M(3) + dt*dm(Vold,M(3)) + SQRT(dt)*mnoise(Vold,M(3),NNa)*RANDNORM() ))
  h   = MIN(1., MAX(0., h + dt*dh(Vold,h) + SQRT(dt)*hnoise(Vold,h,NNa)*RANDNORM()  ))
  N(1) = MIN(1., MAX(0., N(1) + dt*dn(Vold,N(1)) + SQRT(dt)*nnoise(Vold,N(1),NK)*RANDNORM() ))
  N(2) = MIN(1., MAX(0., N(2) + dt*dn(Vold,N(2)) + SQRT(dt)*nnoise(Vold,N(2),NK)*RANDNORM() ))
  N(3) = MIN(1., MAX(0., N(3) + dt*dn(Vold,N(3)) + SQRT(dt)*nnoise(Vold,N(3),NK)*RANDNORM() ))
  N(4) = MIN(1., MAX(0., N(4) + dt*dn(Vold,N(4)) + SQRT(dt)*nnoise(Vold,N(4),NK)*RANDNORM() ))

    ! Detect spikes and compute ISI
    IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
      nspike = nspike + 1
  
    IF (PrintOut .EQ. 2) THEN
      ! Write ISI as output
      IF (nspike .GT. 1) THEN
        print *, t - timelastspike
      END IF
    END IF

      timelastspike = t
      tsincespike = 0.
    END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, M(1)*M(2)*M(3)*h, N(1)*N(2)*N(3)*N(4)
  END IF

    ! Stopping Conditions
    IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
      stop = 1 
    ELSE IF (i .EQ. nt) THEN ! reached last time bin
      stop = 1
    ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
      stop = 1
    END IF

END DO

END SUBROUTINE HHsde_subunit_independent
!***********************************************************************************


!***********************************************************************************
SUBROUTINE HHsde_channel(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! Fox and Lu Channel Based SDE - the correct SDE Approximation

! Read in parameter values from modules
USE ParameterModule
USE mtmod

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: seed, maxnisi, PrintOut
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam

! Declare Local Variables
INTEGER i, j
DOUBLE PRECISION, DIMENSION(1:7) :: dX7, r7
DOUBLE PRECISION, DIMENSION(1:4) :: dX4, r4
DOUBLE PRECISION, DIMENSION(1:7,1:7) :: nanoise
DOUBLE PRECISION, DIMENSION(1:4,1:4) :: knoise
DOUBLE PRECISION  Iin
DOUBLE PRECISION  V, Vold
DOUBLE PRECISION, DIMENSION(1:7) :: XNa
DOUBLE PRECISION, DIMENSION(1:4) :: XK
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike

! Initialize Mersenne Twiste Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF

XNa(1) = 3.*betam(V)**2. * alpham(V) * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(2) = 3.*betam(V) * alpham(V)**2. * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(3) = alpham(V)**3. * betah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(4) = betam(V)**3. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(5) = 3.*betam(V)**2. * alpham(V) * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(6) = 3.*betam(V) * alpham(V)**2. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XNa(7) = alpham(V)**3. * alphah(V) / ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  )
XK(1) =  4.*alphan(V)*betan(V)**3./(alphan(V)+betan(V))**4.
XK(2) =  6.*alphan(V)**2.*betan(V)**2./(alphan(V)+betan(V))**4.
XK(3) =  4.*alphan(V)**3.*betan(V)/(alphan(V)+betan(V))**4.
XK(4) =  alphan(V)**4./(alphan(V)+betan(V))**4.

timelastspike = -9999.
tsincespike = 0.
stop = 0
i = 0
nspike = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, XNa(7), XK(4)
END IF

DO WHILE (stop .EQ. 0)

  i = i+1
  t = DBLE(i)*dt
  tsincespike = tsincespike+dt

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  Vold = V
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,[XNa(7),oned,oned],oned,[XK(4),oned,oned,oned],Iin)
  END IF
  
  CALL dXNa(Vold,XNa(:),dX7)
  CALL dXK(Vold,XK(:),dX4)
  CALL Nanoisematrix_eq(V,NNa, nanoise)  ! Equil approx to avoid problems in square root step
  CALL Knoisematrix_eq(V,NK, knoise)    ! Equil approx to avoid problems in square root step
  
  ! note that noise in fox and lu is scaled by 2 (eq 34)
  r7 = SQRT(2.)* (/RANDNORM(), RANDNORM(), RANDNORM(), RANDNORM(),&
                   RANDNORM(), RANDNORM(), RANDNORM()/)
  r4 = SQRT(2.)* (/RANDNORM(), RANDNORM(), RANDNORM(), RANDNORM() /)

  DO j=1,7
    XNa(j) = XNa(j) + dt*dX7(j)  + SQRT(dt) * DOT_PRODUCT(nanoise(j,:), r7)
  END DO
  DO j=1,4 
    XK(j) = XK(j) + dt*dX4(j) + SQRT(dt) * DOT_PRODUCT(knoise(:,j), r4)
  END DO

    ! Detect spikes and compute ISI
    IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
      nspike = nspike + 1
  
    IF (PrintOut .EQ. 2) THEN
      ! Write ISI as output
      IF (nspike .GT. 1) THEN
        print *, t - timelastspike
      END IF
    END IF

      timelastspike = t
      tsincespike =0.
    END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, XNa(7), XK(4)
  END IF

    ! Stopping Conditions
    IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
      stop = 1 
    ELSE IF (i .EQ. nt) THEN ! reached last time bin
      stop = 1
    ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
      stop = 1
    END IF

END DO ! end while loop

END SUBROUTINE HHsde_channel
!***********************************************************************************

!***********************************************************************************
SUBROUTINE HHmc(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! Markov Chain - Gillespie

! Read in parameter values from modules
USE ParameterModule
USE mtmod

IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: seed, maxnisi, PrintOut
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam

! Declare Local Variables
INTEGER i
DOUBLE PRECISION, DIMENSION(3) :: M
DOUBLE PRECISION, DIMENSION(4) :: N
DOUBLE PRECISION r
DOUBLE PRECISION tswitch, tupdate
DOUBLE PRECISION totalrate
INTEGER, DIMENSION(0:3,0:1) :: Nastate
INTEGER, DIMENSION(0:4) :: Kstate
DOUBLE PRECISION, DIMENSION(1:28) :: rate
DOUBLE PRECISION Iin
DOUBLE PRECISION V, Vold
INTEGER NNaopen, NKopen
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF

! Initialize Channel States
Nastate(0,0)=FLOOR(NNa * betam(V)**3. * betah(V)  / &
( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  ) )
Nastate(1,0)=FLOOR(NNa * 3.*betam(V)**2. * alpham(V) * betah(V) / &
( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))))
Nastate(2,0)=FLOOR(NNa * 3.*betam(V) * alpham(V)**2. * betah(V) / &
( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))))
Nastate(3,0)=FLOOR(NNa * alpham(V)**3. * betah(V) / &
( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  ) )
Nastate(0,1)=FLOOR(NNa *  betam(V)**3. * alphah(V) / &
( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))  ) )
Nastate(1,1)=FLOOR(NNa * 3.*betam(V)**2. * alpham(V) * alphah(V) /&
 ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))))
Nastate(2,1)=FLOOR(NNa *  3.*betam(V) * alpham(V)**2. * alphah(V) /&
 ( (alpham(V)+betam(V))**3. * (alphah(V)+betah(V))))
Nastate(3,1) = NNa - ( Nastate(0,0)+Nastate(1,0)+Nastate(2,0)+Nastate(3,0)+Nastate(0,1)+Nastate(1,1)+Nastate(2,1) )

Kstate(0) = FLOOR(NK * betan(V)**4./(alphan(V)+betan(V))**4.)
Kstate(1) = FLOOR(NK * 4.*alphan(V)*betan(V)**3./(alphan(V)+betan(V))**4.)
Kstate(2) = FLOOR(NK * 6.*alphan(V)**2.*betan(V)**2./(alphan(V)+betan(V))**4.)
Kstate(3) = FLOOR(NK * 4.*alphan(V)**3.*betan(V)/(alphan(V)+betan(V))**4.)
Kstate(4) = NK - (Kstate(0) + Kstate(1) + Kstate(2) + Kstate(3) )
 
! Initialize number of open channels
NNaopen = Nastate(3,1)
NKopen = Kstate(4)

timelastspike = -9999.
tsincespike = 0.
stop = 0
tupdate = 0.
i = 0
nspike = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, REAL(NNaopen)/ REAL(NNa), REAL(NKopen)/REAL(NK)
END IF

DO WHILE (stop .EQ. 0)
  i = i+1
  t = DBLE(i)*dt
  tswitch = t
  tsincespike = tsincespike+dt

  ! Update Voltage (Forward Euler)
  M(1) = DBLE(NNaopen) / DBLE(NNa)
  M(2) = 1.
  M(3) = 1.
  N(1) = DBLE(NKopen) / DBLE(NK)
  N(2) = 1.
  N(3) = 1.
  N(4) = 1.

  Vold = V
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,M,oned,N,Iin)
  END IF

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  ! Update Channel States
  DO WHILE (tswitch .LT. (t+dt)) 

    ! Determine which state switches by partitioning total rate into its 28 components
    rate(1) = 3.*alpham(Vold) * DBLE(Nastate(0,0))
    rate(2) = rate(1) + 2.*alpham(Vold) * DBLE(Nastate(1,0))
    rate(3) = rate(2) + 1.*alpham(Vold) * DBLE(Nastate(2,0))
    rate(4) = rate(3) + 3.*betam(Vold) * DBLE(Nastate(3,0))
    rate(5) = rate(4) + 2.*betam(Vold) * DBLE(Nastate(2,0))
    rate(6) = rate(5) + 1.*betam(Vold) * DBLE(Nastate(1,0))
    rate(7) = rate(6) + alphah(Vold) * DBLE(Nastate(0,0))
    rate(8) = rate(7) + alphah(Vold) * DBLE(Nastate(1,0))
    rate(9) = rate(8) + alphah(Vold) *DBLE( Nastate(2,0))
    rate(10) = rate(9) + alphah(Vold) * DBLE(Nastate(3,0))
    rate(11) = rate(10) + betah(Vold) * DBLE(Nastate(0,1))
    rate(12) = rate(11) + betah(Vold) * DBLE(Nastate(1,1))
    rate(13) = rate(12) + betah(Vold) * DBLE(Nastate(2,1))
    rate(14) = rate(13) + betah(Vold) * DBLE(Nastate(3,1))
    rate(15) = rate(14) + 3.*alpham(Vold) * DBLE(Nastate(0,1))
    rate(16) = rate(15) + 2.*alpham(Vold) * DBLE(Nastate(1,1))
    rate(17) = rate(16) + 1.*alpham(Vold) * DBLE(Nastate(2,1))
    rate(18) = rate(17) + 3.*betam(Vold) * DBLE(Nastate(3,1))
    rate(19) = rate(18) + 2.*betam(Vold) * DBLE(Nastate(2,1))
    rate(20) = rate(19) + 1.*betam(Vold) * DBLE(Nastate(1,1))
    rate(21) = rate(20) + 4.*alphan(Vold) * DBLE(Kstate(0))
    rate(22) = rate(21) + 3.*alphan(Vold) * DBLE(Kstate(1))
    rate(23) = rate(22) + 2.*alphan(Vold) * DBLE(Kstate(2))
    rate(24) = rate(23) + 1.*alphan(Vold) * DBLE(Kstate(3))
    rate(25) = rate(24) + 4.*betan(Vold) * DBLE(Kstate(4))
    rate(26) = rate(25) + 3.*betan(Vold) * DBLE(Kstate(3))
    rate(27) = rate(26) + 2.*betan(Vold) * DBLE(Kstate(2))
    rate(28) = rate(27) + 1.*betan(Vold) * DBLE(Kstate(1))

    ! Total Transition Rate
    totalrate = rate(28)

    ! Exponential Waiting Time Distribution
    tupdate = -LOG(grnd()) / totalrate

    ! Time of Next Switching Event (Exp Rand Var)
    tswitch = tswitch + tupdate

    IF (tswitch .LT. (t+dt)) THEN

      ! Scaled Uniform RV to determine which state to switch
      r = totalrate*grnd()

      IF (r .LT. rate(1)) THEN
        Nastate(0,0) = Nastate(0,0)-1
        Nastate(1,0) = Nastate(1,0)+1 
      ELSE IF (r .LT. rate(2)) THEN
       Nastate(1,0) = Nastate(1,0)-1
       Nastate(2,0) = Nastate(2,0)+1 
      ELSE IF (r .LT. rate(3)) THEN
       Nastate(2,0) = Nastate(2,0)-1
       Nastate(3,0) = Nastate(3,0)+1 
      ELSE IF (r .LT. rate(4)) THEN
       Nastate(3,0) = Nastate(3,0)-1
       Nastate(2,0) = Nastate(2,0)+1  
      ELSE IF (r .LT. rate(5)) THEN
       Nastate(2,0) = Nastate(2,0)-1
       Nastate(1,0) = Nastate(1,0)+1
      ELSE IF (r .LT. rate(6)) THEN
       Nastate(1,0) = Nastate(1,0)-1
       Nastate(0,0) = Nastate(0,0)+1
      ELSE IF (r .LT. rate(7)) THEN
       Nastate(0,0) = Nastate(0,0)-1
       Nastate(0,1) = Nastate(0,1)+1
      ELSE IF (r .LT. rate(8)) THEN
       Nastate(1,0) = Nastate(1,0)-1
       Nastate(1,1) = Nastate(1,1)+1
      ELSE IF (r .LT. rate(9)) THEN
       Nastate(2,0) = Nastate(2,0)-1
       Nastate(2,1) = Nastate(2,1)+1
      ELSE IF (r .LT. rate(10)) THEN
       Nastate(3,0) = Nastate(3,0)-1
       Nastate(3,1) = Nastate(3,1)+1
      ELSE IF (r .LT. rate(11)) THEN
       Nastate(0,1) = Nastate(0,1)-1
       Nastate(0,0) = Nastate(0,0)+1
      ELSE IF (r .LT. rate(12)) THEN
       Nastate(1,1) = Nastate(1,1)-1
       Nastate(1,0) = Nastate(1,0)+1
      ELSE IF (r .LT. rate(13)) THEN
       Nastate(2,1) = Nastate(2,1)-1
       Nastate(2,0) = Nastate(2,0)+1
      ELSE IF (r .LT. rate(14)) THEN
       Nastate(3,1) = Nastate(3,1)-1
       Nastate(3,0) = Nastate(3,0)+1
      ELSE IF (r .LT. rate(15)) THEN
       Nastate(0,1) = Nastate(0,1)-1
       Nastate(1,1) = Nastate(1,1)+1
      ELSE IF (r .LT. rate(16)) THEN
       Nastate(1,1) = Nastate(1,1)-1
       Nastate(2,1) = Nastate(2,1)+1
      ELSE IF (r .LT. rate(17)) THEN
       Nastate(2,1) = Nastate(2,1)-1
       Nastate(3,1) = Nastate(3,1)+1
      ELSE IF (r .LT. rate(18)) THEN
       Nastate(3,1) = Nastate(3,1)-1
       Nastate(2,1) = Nastate(2,1)+1
      ELSE IF (r .LT. rate(19)) THEN
       Nastate(2,1) = Nastate(2,1)-1
       Nastate(1,1) = Nastate(1,1)+1
      ELSE IF (r .LT. rate(20)) THEN
       Nastate(1,1) = Nastate(1,1)-1
       Nastate(0,1) = Nastate(0,1)+1
      ELSE IF (r .LT. rate(21)) THEN
       Kstate(0) = Kstate(0)-1
       Kstate(1) = Kstate(1)+1
      ELSE IF (r .LT. rate(22)) THEN
       Kstate(1) = Kstate(1)-1
       Kstate(2) = Kstate(2)+1
      ELSE IF (r .LT. rate(23)) THEN
       Kstate(2) = Kstate(2)-1
       Kstate(3) = Kstate(3)+1
      ELSE IF (r .LT. rate(24)) THEN
       Kstate(3) = Kstate(3)-1
       Kstate(4) = Kstate(4)+1
      ELSE IF (r .LT. rate(25)) THEN
       Kstate(4) = Kstate(4)-1
       Kstate(3) = Kstate(3)+1
      ELSE IF (r .LT. rate(26)) THEN
       Kstate(3) = Kstate(3)-1
       Kstate(2) = Kstate(2)+1
      ELSE IF (r .LT. rate(27)) THEN
       Kstate(2) = Kstate(2)-1
       Kstate(1) = Kstate(1)+1
      ELSE
       Kstate(1) = Kstate(1)-1
       Kstate(0) = Kstate(0)+1
      END IF

    END IF
 
  END DO

  ! Count Number of Open Channels
  NNaopen = Nastate(3,1)
  NKopen = Kstate(4)

  ! Variables for voltage equation
  M(1) = DBLE(NNaopen)/DBLE(NNa)
  N(1) = DBLE(NKopen)/DBLE(NK)


    ! Detect spikes and compute ISI
    IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
      nspike = nspike + 1
  
    IF (PrintOut .EQ. 2) THEN
      ! Write ISI as output
      IF (nspike .GT. 1) THEN
        print *, t - timelastspike
      END IF
    END IF

      timelastspike = t
      tsincespike = 0.
    END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, REAL(NNaopen)/ REAL(NNa), REAL(NKopen)/REAL(NK)
  END IF

    ! Stopping Conditions
    IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
      stop = 1 
    ELSE IF (i .EQ. nt) THEN ! reached last time bin
      stop = 1
    ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
      stop = 1
    END IF

END DO ! end while loop

END SUBROUTINE HHmc
!***********************************************************************************


!***********************************************************************************
SUBROUTINE HHsde_subunit_quasistationary(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! Subunit SDE with variables raised to powers (Fox and Lu subunit approximation)
! Noise terms changed to quasistationary values

! Read in parameter values from modules
USE ParameterModule
USE mtmod
 
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION dm, dh, dn
DOUBLE PRECISION mnoise_modify, nnoise_modify, hnoise
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: maxnisi
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam
INTEGER, INTENT(IN) :: seed, PrintOut

! Declare Local Variables
INTEGER i
DOUBLE PRECISION V,m,h,n, Vold
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike
DOUBLE PRECISION Iin

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V = StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF
m = alpham(V) / (alpham(V) + betam(V))
h = alphah(V) / (alphah(V) + betah(V))
n = alphan(V) / (alphan(V) + betan(V))

nspike = 0
timelastspike = -9999.
tsincespike =0.
stop = 0
i = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, m**3*h, n**4
END IF

DO WHILE (stop .EQ. 0)
    i = i+1
    t = DBLE(i)*dt
    tsincespike = tsincespike+dt

    Vold = V
    IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
      V= StimParam(1)
    ELSE
      V = V + dt*dV(V,[m,m,m],h,[n,n,n,n],Iin)
    END IF

    ! Define Stimulus Based on StimParam
    Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

    m = MIN(1., MAX(0., m + dt*dm(Vold,m) + SQRT(dt)*mnoise_modify(Vold,NNa)*RANDNORM()))
    h = MIN(1., MAX(0., h + dt*dh(Vold,h) + SQRT(dt)*hnoise(Vold,h,NNa)*RANDNORM() ))
    n = MIN(1., MAX(0., n + dt*dn(Vold,n) + SQRT(dt)*nnoise_modify(Vold,NK)*RANDNORM()))

    ! Detect spikes and compute ISI
    IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
      nspike = nspike + 1
  
    IF (PrintOut .EQ. 2) THEN
      ! Write ISI as output
      IF (nspike .GT. 1) THEN
        print *, t - timelastspike
      END IF
   END IF

      timelastspike = t
      tsincespike = 0.
    END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, m**3*h, n**4
  END IF

    ! Stopping Conditions
    IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
      stop = 1 
    ELSE IF (i .EQ. nt) THEN ! reached last time bin
      stop = 1
    ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
      stop = 1
    END IF

END DO

END SUBROUTINE HHsde_subunit_quasistationary
!***********************************************************************************

!***********************************************************************************
SUBROUTINE HHsde_channel_quasistationary(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)
! CURRENT QS METHOD: Forward Euler Solve HH ODE + non-Markov noise added 
! Need to load look up tables for noise terms

! Read in parameter values from modules
USE ParameterModule
USE mtmod
  
IMPLICIT NONE

! Declare Functions
DOUBLE PRECISION alpham, betam, alphah, betah, alphan, betan
DOUBLE PRECISION dV
DOUBLE PRECISION dm, dh, dn
DOUBLE PRECISION RANDNORM

! Declare Input Variables
INTEGER, INTENT(IN) :: NNa, NK
INTEGER, INTENT(IN) :: nt
DOUBLE PRECISION, INTENT(IN) :: dt
INTEGER, INTENT(IN) :: maxnisi
DOUBLE PRECISION, DIMENSION(1:5), INTENT(IN) :: StimParam
INTEGER, INTENT(IN) :: seed, PrintOut

! Declare Local Variables
INTEGER i
DOUBLE PRECISION mode,hode,node
DOUBLE PRECISION r, r1, r2
DOUBLE PRECISION taum, tauh, taun, vark, varna
DOUBLE PRECISION, DIMENSION(1:14001) :: k1, k2, k3, k4
DOUBLE PRECISION, DIMENSION(1:14001) :: na1, na2, na3, na4, na5, na6, na7
DOUBLE PRECISION tempk1, tempk2, tempk3, tempk4
DOUBLE PRECISION tempna1, tempna2, tempna3, tempna4, tempna5, tempna6, tempna7
DOUBLE PRECISION,DIMENSION(1:4,0:14000) :: KNoiseData
DOUBLE PRECISION,DIMENSION(1:7,0:14000) :: NaNoiseData
INTEGER Vindex
DOUBLE PRECISION V,xna,xk, Vold
INTEGER stop
INTEGER nspike
DOUBLE PRECISION t, timelastspike, tsincespike
DOUBLE PRECISION Iin

! Initialize Mersenne Twister Random Number
CALL sgrnd(seed)

! Read in Coefficient values
open(unit=1, action='READ', file='NoiseData_Na.txt')    ! note the data file has to be findable when the command is executed
        read(1,1,end=100) NaNoiseData
1       format(f10.8) 
        close(1)
100 continue
na1 = NaNoiseData(1,:)
na2 = NaNoiseData(2,:)
na3 = NaNoiseData(3,:)
na4 = NaNoiseData(4,:)
na5 = NaNoiseData(5,:)
na6 = NaNoiseData(6,:)
na7 = NaNoiseData(7,:)

open(unit=2, action='READ',file='NoiseData_K.txt')   ! note the data file has to be findable when the command is executed
        read(2,2,end=200) KNoiseData
2       format(f10.8) 
        close(2)
200 continue
k1 = KNoiseData(1,:)
k2 = KNoiseData(2,:)
k3 = KNoiseData(3,:)
k4 = KNoiseData(4,:)

! Initial Values
Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt)
IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
  V= StimParam(1)
ELSE
  V = Erest ! Resting Potential
END IF
mode = alpham(V) / (alpham(V) + betam(V))
hode = alphah(V) / (alphah(V) + betah(V))
node = alphan(V) / (alphan(V) + betan(V))

xna = mode**3.*hode
xk = node**4.

r1=0.
r2=0.

tempna1 = 0.
tempna2 = 0.
tempna3 = 0.
tempna4 = 0.
tempna5 = 0.
tempna6 = 0.
tempna7 = 0.

tempk1 = 0.
tempk2 = 0.
tempk3 = 0.
tempk4 = 0.

nspike = 0
timelastspike = -9999.
tsincespike = 0.
stop = 0
i = 0
t = i*dt

IF (PrintOut .EQ. 1) THEN
  print *, t, V, xna, xk
END IF

DO WHILE (stop .EQ. 0)
  i = i+1
  t = DBLE(i)*dt
  tsincespike = tsincespike+dt

  Vold = V
  IF (StimParam(5) .EQ. 1.) THEN ! Voltage Clamp
    V= StimParam(1)
  ELSE
    V = V + dt*dV(V,[xna,oned,oned,oned],oned,[xk,oned,oned,oned],Iin)
  END IF

  ! Define voltage index that is used to find correct values in look up table of coefficients for the noise processes
  Vindex = NINT( (MIN(120., MAX(-20.,Vold ) )+ 20.) * 100.)

  ! Define Stimulus Based on StimParam
  Iin = StimParam(1) + StimParam(2)*RANDNORM()/SQRT(dt) + StimParam(3)*SIN(2.*PI*StimParam(4)*t/1000.)

  ! Assume gating variables evolve deterministically
  mode = mode + dt*dm(Vold,mode) 
  hode = hode + dt*dh(Vold,hode) 
  node = node + dt*dn(Vold,node)

  ! Noise Processes solved separately
  taum = 1./ (alpham(Vold) + betam(Vold))
  tauh = 1./ (alphah(Vold) + betah(Vold))
  taun = 1./ (alphan(Vold) + betan(Vold))

  varna = (alpham(Vold)/(alpham(Vold)+betam(Vold)))**3. * &
(alphah(Vold)/(alphah(Vold)+betah(Vold))) * &
(1. - (alpham(Vold)/(alpham(Vold)+betam(Vold)))**3.*&
(alphah(Vold)/(alphah(Vold)+betah(Vold)))) / DBLE(NNa)

  vark = ((alphan(Vold)/(alphan(Vold)+betan(Vold)))**4. *&
 (1.-(alphan(Vold)/(alphan(Vold)+betan(Vold)))**4.) / DBLE(NK))

  ! 7ple Gaussian Process Process for Na Channel
  r =  RANDNORM()
  tempna1 = EXP(-1.*dt/taum) * tempna1 + na1(Vindex)*sqrt(dt)* r * SQRT(varna)
  tempna2 = EXP(-2.*dt/taum) * tempna2 + na2(Vindex)*sqrt(dt)* r * SQRT(varna)
  tempna3 = EXP(-3.*dt/taum) * tempna3 + na3(Vindex)*sqrt(dt)* r * SQRT(varna)
  tempna4 = EXP(-1.*dt/tauh) * tempna4 + na4(Vindex)*sqrt(dt)* r * SQRT(varna)
  tempna5 = EXP(-dt * ((taum+tauh)/(taum*tauh)) ) *tempna5 + na5(Vindex)*sqrt(dt)*r*&
 SQRT(varna)
  tempna6 = EXP(-dt * ((taum+2.*tauh)/(taum*tauh)) ) *tempna6 + na6(Vindex)*sqrt(dt)*r*&
 SQRT(varna)
  tempna7 = EXP(-dt * ((taum+3.*tauh)/(taum*tauh)) ) *tempna7 + na7(Vindex)*sqrt(dt)*r*&
 SQRT(varna)
  r1 = tempna1 + tempna2 + tempna3 + tempna4 + tempna5 + tempna6 + tempna7

  ! 4ple Guassian Process For K Channel 
  r = RANDNORM()
  tempk1 =  EXP(-1.*dt/taun) * tempk1 + k1(Vindex)*sqrt(dt) * r *SQRT(vark)
  tempk2 =  EXP(-2.*dt/taun) * tempk2 + k2(Vindex)*sqrt(dt) * r *SQRT(vark)
  tempk3 =  EXP(-3.*dt/taun) * tempk3 + k3(Vindex)*sqrt(dt) * r *SQRT(vark)
  tempk4 =  EXP(-4.*dt/taun) * tempk4 + k4(Vindex)*sqrt(dt) * r *SQRT(vark)
  r2 = (tempk1 + tempk2 + tempk3 + tempk4)

  ! Add fluctuations to the deterministic paths
  xna = MIN(1., MAX(0., mode**3.*hode + r1 ))
  xk  = MIN(1., MAX(0., node**4. + r2 ))

    ! Detect spikes and compute ISI
    IF ((V .GT. spikethreshold) .AND. ((t-timelastspike) .GT. minisi)) THEN
      nspike = nspike + 1
  
    IF (PrintOut .EQ. 2) THEN
      ! Write ISI as output
      IF (nspike .GT. 1) THEN
        print *, t - timelastspike
      END IF
    END IF

      timelastspike = t
      tsincespike =0.
    END IF  

  IF (PrintOut .EQ. 1) THEN
    print *, t, V, xna, xk
  END IF

    ! Stopping Conditions
    IF (nspike .EQ. (maxnisi+1)) THEN  ! reached max number of isi
      stop = 1 
    ELSE IF (i .EQ. nt) THEN ! reached last time bin
      stop = 1
    ELSE IF ( tsincespike .GT. maxisi) THEN  ! haven't had a spike in a really long time
      stop = 1
    END IF

END DO

END SUBROUTINE HHsde_channel_quasistationary
!***********************************************************************************