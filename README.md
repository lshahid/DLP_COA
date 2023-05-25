# DLP_COA

## Description

This repository contains the code for the implementation of a Distributed Lumped Parameter (DLP) model for coarctation of the aorta. The DLP model we used is based on Mirramezani & Shadden [1], and Pewowaruk & Roldán-Alzate [2]. This DLP model is a reduced-order model of blood flow based on simplified physics. The code calculates the resistances based on energy losses due to viscous dissipation, unsteadiness, flow separation, vessel curvature and vessel bifurcations. Please see [1] for more details. The motivation of this DLP model is that it can simulate hemodynamic flow and pressure significantly faster than computational fluid dynamics (CFD).

## Instructions to use this code:
1. Run DLP_script.m
2. Select the centerline file created by Mimics.
3. Select the STL of the segmented aorta.
4. Based on the plot of the radius vs point number, select the starting and ending point number of all segments.
5. Select inlets and outlets.
6. Select flow file for inlet.
7. Select file containing resistance and capacitance values for Windkessel outlet conditions.
8. Select whether you want impedance from pulse wave velocity or not.
9. Select which Windkessel parameter to use for each outlet.
10. Type in the name of the .mat file containing the parameters inputted in steps 2-9. Next time you wwant to run the simulation, you can skip steps 2-9 and just load this .mat file.
11. The DLP simulation will run and will output the pressure and flow rates in a .mat file.

## References

[1] Mirramezani, M., Shadden, S.C. A Distributed Lumped Parameter Model of Blood Flow. Ann Biomed Eng 48, 2870–2886 (2020). https://doi.org/10.1007/s10439-020-02545-6

[2] Pewowaruk, R., Roldán-Alzate, A. A distributed lumped parameter model of blood flow with fluid-structure interaction. Biomech Model Mechanobiol 20, 1659–1674 (2021). https://doi.org/10.1007/s10237-021-01468-y

## Acknowledgements

This code was developed by Ryan Pewowaruk, and implemented for coarctation of the aorta by Labib Shahid and Matthew Culver.
