# Computation of analytical sensitivity matrix for the frequency-domain EM data: MATLAB code 
Name of the program: ANALYTIC_SENSITIVITY_FREQ_DOMAIN_EM.m
Title of manuscript: Computation of analytical sensitivity matrix for the frequency-domain EM data: MATLAB code.
Author name: Jide Nosakare Ogunbo
Affiliation: The Federal University of Technology, Akure
Address: PMB 704, Akure, Ondo state, Nigeria.
Email address: nosajide@yahoo.com; nosajide@gmail.com; jnogunbo@futa.edu.ng

                           Code description:
The code computes the analytical sensitivity of 1D frequency-domain EM data.
It uses the basic differentiation rules and logarithmic differentiation to
compute the 1D recursive transverse electric frequency=domain electromagnetic
foward response.
Results are speedily and accurately produced. Compared with the finite difference
approach the analytical approach actually gains speed as the number of parameters
increases. Thus it encourages the use of more parameters in computation rather
than discouraging when huge data set is involved.

The code can thus be used in areas where sensitivity matrix is needed, for example
for data/model resolution analysis, experimental design and particularly nonlinear
inversion 1D frequency EM data. It is simply a harvest of speed gain without loss
in accuracy.

The document "Jide_Ogunbo_Instruction guide on running the MATLAB code.doc" gives a
description on the usage of the code.

The DRIVER.m is a sample code for running the code.
The OUTPUT.txt shows the result of the example in the DRIVER.m
