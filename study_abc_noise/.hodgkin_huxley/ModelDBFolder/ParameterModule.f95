!***********************************************************************************
MODULE ParameterModule

IMPLICIT NONE

SAVE

! HH Parameters (Hodgkin and Huxley, 1952)
DOUBLE PRECISION, PARAMETER :: C = 1 ! muF /cm^2
DOUBLE PRECISION, PARAMETER :: gNa = 120 ! mS/cm^2
DOUBLE PRECISION, PARAMETER :: ENa = 115 !  mV  ! changed from 120 to 115 because of apparent typo in Izhikevich book
DOUBLE PRECISION, PARAMETER :: DNa = 60. ! Channels/mu m^2 Channel Density
DOUBLE PRECISION, PARAMETER :: gK = 36 ! mS/cm^2
DOUBLE PRECISION, PARAMETER :: EK = -12 !  mV
DOUBLE PRECISION, PARAMETER :: DK = 18. ! Channels/mu m^2 Channel Density
DOUBLE PRECISION, PARAMETER :: gL = 0.3 ! mS / cm^2
DOUBLE PRECISION, PARAMETER :: EL = 10.6 ! mV
DOUBLE PRECISION, PARAMETER :: Erest = 0.
DOUBLE PRECISION, PARAMETER :: oned = 1.0
DOUBLE PRECISION, PARAMETER :: spikethreshold = 60. ! Threshold for spike count (mV)
DOUBLE PRECISION, PARAMETER :: minisi = 2.0 ! Minimum ISI for spike count (ms)
DOUBLE PRECISION, PARAMETER :: maxisi = 2000. ! Maximum ISI to kill program (ms)
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462643	! PI (Find a better way to compute this?)      

END MODULE ParameterModule
!***********************************************************************************
