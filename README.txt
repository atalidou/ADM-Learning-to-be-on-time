Notes/comments regarding the programs 

	(1) ADM - Modelling plasticity.cpp
	(2) ADM - Information transmission.cpp
	(3) ADM - Kramer mean firing rate and variance.cpp

to generate results presented in "Homeostatic coordination and up-regulation of neural activity by activity-dependent myelination" by Talidou et al. 

The attached programs are a routine meant to be ran on any C++ compiler. We recomment using  Bloodshed Dev C++ (https://www.bloodshed.net/ ) which is a Open Source C/C++ IDE for Windows for optimal compatibility.


The program (1) implements a simple Euler-Maruyama numerical integration scheme applied to a system of N (neurons)+N^2(axons) delay differential equations with Poisson noise.

Once compiled and run, the program (1) simulates the dynamics of the network for a duration of T steps (1ms) using a sliding time window approach. Once this is done, statistics are extracted and data is outputed for plotting purposes.
----------
The program (2) implements a simple Euler-Maruyama numerical integration scheme applied to a system of N (neurons)+N^2(axons) delay differential equations with Poisson noise.

Once compiled and run, the program simulates the mutual information of the network for a duration of T steps (1ms) using a sliding time window approach. Once this is done, statistics are extracted and data is outputed for plotting purposes.
----------
The program (3) implements a simple Euler-Maruyama numerical integration scheme applied to a system of N (neurons)+N^2(axons) delay differential equations with additive noise.

Once compiled and run, the program simulates the dynamics of the network for a duration of T steps (1ms). Once this is done, statistics (mean value and variance) are extracted and data is outputed for plotting purposes. 


To compile and run program (2) the data 'LENGTHS100.txt' and 'CVS100.txt' are used.

To generate manuscript figures, the values of the parameters were modified accordingly manually - using values specified in the manuscript.

When required, additional loops were implemented to cycle between different parameter values (e.g., input, spatial scales).
Ouput data was imported in a plotting software (OriginPro 2020 (64-bit) SR1 9.7.0.188 (Government)) to generate and edit figures. Powerpoint was used to collect and assemble panels and anotations. 

Mathematical analysis was performed independently, supported by Maple (MapleSoft, Waterloo 2020). 




