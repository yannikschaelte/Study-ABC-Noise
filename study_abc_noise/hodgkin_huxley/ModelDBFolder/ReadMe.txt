ReadMe.txt

This is the readme for the model code associated with the paper:

Joshua H. Goldwyn, Nikita S. Imennov, Michael Famulare, Eric
Shea-Brown. (submitted) On stochastic differential equation models for
ion channel noise in Hodgkin-Huxley neurons.

Submitted manuscript available at http://arxiv.org/abs/1009.4172

The following files have been uploaded to ModelDB:

HH_master.f95 -- all subroutines and programs used to solve HH
                 equations
HH_run.95 -- program file for simulating HH equations (uses command
                 line inputs, see below)
Makefile -- make executable
MT19937.f90 -- Mersenne twister code (written by Richard Woloshyn and
                 downloaded from
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html)
NoiseData_K.txt and NoiseData_Na.txt -- data files used for noise
                 terms in colored noise quasistationary model
NonMarkov_K.mw and NonMarkov_Na.mw -- Maple files used to create
                 NoiseData_K.txt and NoiseData_Na.txt, respectively
ParameterModule.f95 -- parameter values used in HH model


Usage:

Create HH_run using make command

Note from the ModelDB Administrator: On my Fedora Core 11 machine I
found the lapack and blas library files were not found until I
replaced LDFLAGS in Makefile to read

LDFLAGS = -L/usr/local/epd/lib -I/usr/local/epd/include /usr/lib/liblapack.so.3 /usr/lib/libblas.so.3

Then the command line can be used:

./HH_run [Model Number] [Membrane Area] [# Time Steps] \
[Time Step Size] [# ISIs] [DC] [Noise] [Sine Amplitude] [Sine Frequency] \
[Voltage Clamp] [Data to Print Out] [Random Number Seed]

output is data (format of which determined by [Data to Print Out]
option, see below)

Where :
Model Number: 	0 = ODE
		1 = Markov Chain
		2 = SDE Channel (Fox and Lu, 1994)
		3 = SDE Subunit Identical (Fox and Lu, 1994)
		4 = SDE Subunit Independent (Shuai and Jung, 2002)
		5 = SDE Subunit Quasistationary
		6 = Channel Quasistationary

Applied Current is of the form [DC] + [Noise]*N(0,1) + [Sine
Amplitude]*sin(2*Pi*[Sine Frequency]*t)

Voltage Clamp: 	0 = No
		1 = Yes

Data to Print Out: 1 = t, V, proportion of open Na channels,
                       proportion of open K channels
		   2 = Interspike intervals

Two specific examples of HH_run are provided in the scripts
example1.sh and example2.sh that implement the following

EXAMPLE 1:

Run all models for a constant input (strength of DC input is 7 (micro
amp / cm^2). Write out first 30 interspike intervals (in ms).

EXAMPLE 2:

Run all models in voltage clamp (voltage clamp to 20 mV) for 100 ms
(1E4 time steps with 0.01ms step size). Write out time, voltage,
proportion of open Na channels and proportion of open K channels
