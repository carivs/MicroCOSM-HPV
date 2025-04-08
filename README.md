# MicroCOSM-HPV
This repository includes the code base and documentation for MicroCOSM-HPV v1.0

Microsimulation for the Control of South African Morbidity and Mortality (MicroCOSM) is an agent based model developed to simulate the epidemiology of various diseases and the social processes that drive these diseases in South Africa. The version included in this repository focuses on the simulation of HIV, HPV and cervical cancer in South Africa. 

## Description of files

* Microsimulation.cpp and .h are the main code files for the model.
* Other C++ files contain functions for the random number generator and statistical distributions. 
* Files in the /randomHPV folder numbered from 1-96 contain random numbers that generate the 96 best fitting parameter combinations.
* Files in the /input folder contain epidemiological data and fixed parameters.
* CCPreventionStrategiesJ.json is a control file to update vaccinations and screening intervention parameters.

## Usage

To compile the code and create the executable "Microsimulation":

```gcc
g++ -std=c++14 Microsimulation.cpp StatFunctions.cpp mersenne.cpp -o Microsimulation
```

When running the code, random numbers are read from /randomHPV/xRandomUniformHPV.txt, where x=1 to 96. To run the code for one parameter combination, where x=1: 

```
./Microsimulation 1 
```

The code is set up to run in parallel on a cluster, with each parameter combination running on one core. The files will run the 'status quo' of HPV vaccination and cervical cancer screening in South Africa, from 1985 to 2120, 50 iterations per parameter combination (IterationsPerPC). The output will be text files containing the average (over 50 iterations) age standardised rates (ASR) between 1985 and 2120. One parameter combination and 50 iterations should run for about 10 hours.

## Technical Appendix

Data sources, assumptions, calibration methods, and model fits are described in MicroCOSMHPVappendix.pdf.

## Support

Please email me at carivs@sun.ac.za if you have any questions.
