********************************************************************
* PHITS-based Analytical Radiation Model in the Atmosphere (PARMA) *
*                     Version 4.10 (2021/03/21)                    *
********************************************************************

PARMA4 can calculate the atmospheric cosmic-ray fluxes & dose rates.
Details of the model is described in
T.Sato, Analytical model for estimating terrestrial cosmic-ray fluxes nearly anytime and anywhere in the world; Extension of PARMA/EXPACS, 10(12): e0144679 (2015)
T.Sato, Analytical model for estimating the zenith angle dependence of terrestrial cosmic ray fluxes, PLOS ONE, 11(8): e0160390 (2016)
This program is designed to be implemented in some other systems such as route-dose calculation systems. 
If you would like to calculate cosmic-ray fluxes or doses at certain conditions, 
it is much easier to use EXPACS.
You can download EXPACS from the webiste below
http://phits.jaea.go.jp/expacs

**** Message from the author ********
This program was translated from the Fortran version of PARMA.
I learned c++ for this translation purpose, and experts of c++ may feel
strange to some parts of this program.
Your feedback on programming is very appreciated.

**** Table of contents ******************
  Readme.txt:      This file
  main.cpp:        Main program of PARMA for calculating cosmic-ray fluxes & doses
  main-simple.cpp: Main program for calculating flux at a certain condition (easy to start up)
  main-generation.cpp: Main program for generating energy and direction of cosmic-rays at a certain condition (for Monte Carlo simulation)
  subroutines.cpp: Subroutines of PARMA
  /AngOut:         Sample output files for angular distribution calculation
  /condition:      Sample input files for flux & dose calculation
  /dcc:            Database for fluence to dose (or other quantity) conversion coefficients
  /DoseOut:        Sample output files for dose calculation
  /GeneOut:        Sample output file from main-generator.f90
  /input:          Databases used in PARMA
  /SpecOut:        Sample output files for flux calculation

Q1. How to use?
A1. You have to compile the programs using a c++ compiler.

Q2. How to compile and execute?
A2. To start up, it is easier to use main-simple.cpp as the control routine.
    In that case, you have to compile main-simple.cpp and subroutines.cpp at once.
    The followings are an example of commands that should be typed (gcc case)
    > g++ main-simple.cpp subroutines.cpp
    > a.exe

Q3. How to change the condition?
A3. You can change the condition such as time and location by changing
    the parameters written in main-simple.cpp.

Q4. How to use this program as an event generator?
A4. You have to use main-generator.cpp instead of main-simple.cpp as the control routine.
    You have to specify the energy and angle ranges for generation, 
    the number of particles to be generated, and the radius of the irradiation sphere "radi".
    If you would like to irradiate some objects by cosmic-rays in your simulation,
    you have to put them incide the irradiation sphere.
    The variables "e" is energy, "u,v,w" is the direction vector, 
    and "x,y,z" is the position of the generated particle.
    Note that w=1.0 and -1.0 indicate the heaven and ground directions, respectively.
    The results are written in "GeneOut/generation.out" file.
    For normalization, the your simulation results should be multiplied with 
    the total flux (/cm2/s) written in flux.out and pi*radi*radi.

Q5. How to calculate dose?
A5. You have to use main.cpp instead of main-simple.cpp as the control routine.
    To specify the calculation condition, you have to make your own input file
    in "condition" folder by referring sample input files such as "Tokyo-Smin.inp".
    At the 1st line of the input file, you have to specify isout and istyle parameters.
    isout=0: No output for flux data for each condition, =1: Output Flux data in "SpecOut" folder
    istyle= 0:W-index, cutoff rigidity(GV), atmospheric depth(g/cm^2), g(direct input)
          = 1:year, month, day, latitude(deg), longitude(deg), atmospheric detph(g/cm^2), g(direct input)
          = 2:year, month, day, latitude(deg), longitude(deg), altitude (m), g(direct input) 
          = 3:year, month, day, latitude(deg), longitude(deg), altitude (ft), g(direct input)
          = 4:year, month, day, cutoff rigidity(GV), altitude (m), g(direct input)
    After 2nd lines, you have to specify the conditions

Q6. How to change the surrounding conditions, such as ground or aircraft
A6. You have to change "g" parameter.
    If you want to calculate the cosmic-ray spectra in the ideal atmosphere
    (i.e. without considering the surrounding effect), you should set g=10.0
    If you want to calculate the ground level cosmic-ray spectra, you should
    set 0 =< g =< 1. In this case, "g" means the water fraction in ground.
    Ground level muon correction is also considered in this mode.
    If you do not know the water fraction at your specified location, 
    I recommended to use 0.15 for "g".
    If you want to calculate the cosmic-ray spectra in aircraft, you should
    set -10 < g < 0. In this case, the absolute value of "g" indicates
    the mass of the aircraft in the unit of 100 ton. 
    It should be noted that only neutron spectrum is influenced by 
    the surrounding condition.
    If you want to calculate the angular differential cosmic-ray fluxes without
    the albedo from the Earth (i.e. black hole mode), you should set g=100.0.

Q7. How to calculate other quantity such as count rates of neutron monitor?
A7. You can select the evaluation quantity by changing "dccname" parameter in "main.cpp".
    Four quantities can be calculated, which are the effective dose for isotropic
    irradiation, H*(10), dose rate in air, and count rate of 6 tube neutron monitor.

Q8. What is the conditions to use this program?
A8. For non-commercial use, you should refer the following manuscripts in any published use of this program,
    T.Sato, Analytical model for estimating terrestrial cosmic-ray fluxes nearly anytime and anywhere in the world; Extension of PARMA/EXPACS, 10(12): e0144679 (2015)
    T.Sato, Analytical model for estimating the zenith angle dependence of terrestrial cosmic ray fluxes, PLOS ONE, 11(8): e0160390 (2016)
    The commercial use of this program is NOT allowed with a prior agreement with JAEA

Contact
  If you have any questions or requests on this program, 
  please E-mail to nsed-expacs@jaea.go.jp
