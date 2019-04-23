!***********************************************************************************
PROGRAM HH_run
! Run HH models with channel noise
!
! Command Line Inputs:
! Model: 0 = ODE
!	 1 = Markov Chain
!	 2 = SDE Channel (Fox and Lu, 1994)
!	 3 = SDE Subunit Identical (Fox and Lu, 1994)
!	 4 = SDE Subunit Independent (Shuai and Jung, 2002)
!	 5 = SDE Subunit Quasistationary
!	 6 = Channel Quasistationary
! Membrane Area (mu m^2)
! Maximum number of time steps before terminating simulation
! Increment of time step (ms)
! Maximum number of ISIs before terminating simulation
! Stimulus - DC Input
! Stimulus - White noise scaling
! Stimulus - Amplitude of sinusoidal input
! Stimulus - Frequency of sinusoidal input
! Whether simulation is in voltage clamp (0=no, 1=yes)
! What to data to print: 1 = Voltage and proportions of open channels
!			 2 = Interspike intervals
! Random number seed (integer)

! Read in parameter values from modules
USE ParameterModule
  
IMPLICIT NONE

! Declare Local Variables
INTEGER commandcount, length
INTEGER model
DOUBLE PRECISION area
INTEGER NNa, NK
INTEGER nt
DOUBLE PRECISION dt
DOUBLE PRECISION, DIMENSION(1:5) :: StimParam
INTEGER PrintOut
CHARACTER(12) :: commandin
INTEGER maxnisi
REAL seed

! Read in Command Line Inputs
  commandcount = COMMAND_ARGUMENT_COUNT()

! Which Model To Use
  CALL GET_COMMAND_ARGUMENT(1, commandin,length)
  CALL STR2INT(commandin,length,model)

! Membrane Area
  CALL GET_COMMAND_ARGUMENT(2, commandin,length)
  CALL STR2DBLE(commandin,length,area)
  CALL NumberChannel(area,NNa,NK) 

! Max Number of Time Steps
  CALL GET_COMMAND_ARGUMENT(3, commandin,length)
  CALL STR2INT(commandin,length,nt)

! Time Step
  CALL GET_COMMAND_ARGUMENT(4, commandin,length)
  CALL STR2DBLE(commandin,length,dt)

! Max Number of ISI
  CALL GET_COMMAND_ARGUMENT(5, commandin,length)
  CALL STR2INT(commandin,length,maxnisi)

! DC input
  CALL GET_COMMAND_ARGUMENT(6, commandin,length)
  CALL STR2DBLE(commandin,length,StimParam(1))

! White Noise Variance
  CALL GET_COMMAND_ARGUMENT(7, commandin,length)
  CALL STR2DBLE(commandin,length,StimParam(2))

! Sinusoidal Frequency
  CALL GET_COMMAND_ARGUMENT(8, commandin,length)
  CALL STR2DBLE(commandin,length,StimParam(3))

! Sinusoidal Amplitude
  CALL GET_COMMAND_ARGUMENT(9, commandin,length)
  CALL STR2DBLE(commandin,length,StimParam(4))

! Whether in Voltage Clamp
  CALL GET_COMMAND_ARGUMENT(10, commandin,length)
  CALL STR2DBLE(commandin,length,StimParam(5))

! What data to print
  CALL GET_COMMAND_ARGUMENT(11, commandin,length)
  CALL STR2INT(commandin,length,PrintOut)

! Random Number Seed INTEGER
  CALL GET_COMMAND_ARGUMENT(12, commandin,length)
  CALL STR2INT(commandin,length,seed)


!!!!!!!!!! ODE !!!!!!!!!!
IF (model .EQ. 0) THEN
  CALL HHode(StimParam,nt,dt,maxnisi,PrintOut,seed)

!!!!!!!!!! Markov Chain !!!!!!!!!!
ELSE IF (model .EQ. 1) THEN
  CALL HHmc(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

!!!!!!!!!! Fox and Lu Channel SDE !!!!!!!!!!
ELSE IF (model .EQ. 2) THEN
  CALL HHsde_channel(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

!!!!!!!!!! Fox and Lu Subunit SDE  !!!!!!!!!!
ELSE IF (model .EQ. 3) THEN
  CALL HHsde_subunit_identical(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

!!!!!!!!!! Shuai and Jung !!!!!!!!!!
ELSE IF (model .EQ. 4) THEN
  CALL HHsde_subunit_independent(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

!!!!!!!!!! Quasistationary (subunit) !!!!!!!!!!
ELSE IF (model .EQ. 5) THEN
  CALL HHsde_subunit_quasistationary(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

ELSE IF (model .EQ. 6) THEN
  CALL HHsde_channel_quasistationary(StimParam,NNa,NK,nt,dt,maxnisi,PrintOut,seed)

END IF

END PROGRAM HH_run
!***********************************************************************************
