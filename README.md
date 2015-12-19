# interThermalPhaseChangeFoam
CFD simulation platform for liquid-vapor thermal phase change flows 

Copyright 2015: Alexander S Rattner, Mahdi Nabil, Sanjay S. Adhikari

[Multiscale Thermal Fluids and Energy (MTFE) Laboratory](http://sites.psu.edu/mtfe/)

The Pennsylvania State University

## Overview
interThermalPhaseChangeFoam is an extensible open source volume-of-fluid (VOF) based simulation tool for studying two-phase fluid flows with thermally driven phase phenomena. The solver approach is based on the adiabatic two-phase interFoam code developed by [OpenCFD Ltd.](http://openfoam.com/). Target applications for interThermalPhaseChangeFoam include:

* Film condensation
* Dropwise condensation
* Evaporation
* Nucleate boiling

![Dropwise Condensation Example](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/04/DropwiseCond_Sigma_1E-3sm.gif)

If you use this solver in a project or scholarly work, we ask that you include a citation for [Rattner and Garimella (2014)](http://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1829850). 

## Installation
...

## Tutorial cases
* Horizontal Film Condensation (Stefan Problem)
* Smooth Nusselt Falling Film Condensation
* Wavy Nusselt Falling Film Condensation
* Two-dimensional Nucleate Boiling in a Cavity
* Bubble Condensation

## Algorithm
Flow Simulation Algorithm Summary:
  Initialize simulation data and phase change model 
  WHILE t<t_end DO
    1. Update delta_t for stability
    2. Update fluid and turbulence properties
    3. Update phase change model 
    4. Phase equation sub-cycle
    5. DO PIMPLE
      1. Form u equation
      2. PISO
          1. Obtain and correct face fluxes
          2. Solve p-Poisson equation
          3. Correct u
    6. LOOP
    7. Update h(T)
    8. DO Energy Loop
      1. Solve h equation
      2. Update T(h)
    9. LOOP
  LOOP
  
Two sample tutorial cases, i.e. Horizontal film condensation and Smooth Nusselt falling film condensation are validated versus the available analytical solutions in the literature with less than 2% error. The corresponding MATLAB scripts included in the aforementioned tutorial cases folders are CheckStefan.m and CheckNusselt.m, respectively. 

## Phase Change Models
A number of  phase change models are included with the solver, and are described below:
o	HiLoRelaxed – An improved version of the model of Rattner and Garimella (2014) that determines the phase change heat sources so that interface cells recover the saturation temperature at each time step. This model performs a graph scan over mesh cells, and applies phase change on the two-cell thick interface layer about user-specified threshold values of α1. Different high and low threshold values for condensation and evaporation, respectively, can be specified, which has been found to reduce numerical smearing of the interface. Numerical under-relaxation of the phase change rate is supported, which can improve numerical stability.
o	HiLoRelaxedSplit – A modified version of the above model, which splits the liquid and vapor portions of  , and applies them on the respective sides of the interface (Rattner, 2015). This approach yields better conservation of the two phases, and reduces smearing of the interface during evaporation.
o	HiLoNoPCV – A modified version of HiLoRelaxed that sets   to 0. This model can yield accurate results with reduced CFL time step constraints for cases without significant vapor-phase effects on phase change (e.g., falling film condensation in a quiescent vapor medium).
o	Yang – An implementation of the empirical rate parameter model of Yang et al. (2008).

## Example applications
* Film Boiling over Heat Generating Rod:        https://www.youtube.com/watch?v=HrVNSpSxnY4
* Two-Phase Flow Jet with Heat Transfer:        https://www.youtube.com/watch?v=aR9arq0Sc-0
* Two-Phase Flow between Two Concentric Pipes:  https://www.youtube.com/watch?v=SAP9-ohW3TY
* Water Discharge in a Liquid Container:        https://www.youtube.com/watch?v=jRguJd2EDgg

## Contribute
MTFE welcomes collaboration with other investigators studying phase-change flows. Please [contact us](mailto:Alex.Rattner@psu.edu) if you are interested in expanding the solver or find bugs to correct. We can also provide limited support (on a case-by-case basis) or consulting servies.

## License
Information about the GPL V3.0 license is included in the "GNU License" file in the main CFD-PC directory.

## References
1. Rattner, A.S., 2015. Single-pressure absorption refrigeration systems for low-source-temperature applications. Ph.D. Thesis, Georgia Institute of Technology, Atlanta, GA. [Link](https://smartech.gatech.edu/handle/1853/53912).
2. Rattner, A.S., Garimella, S., 2014. Simple mechanistically consistent formulation for volume-of-fluid based computations of condensing flows. Journal of Heat Transfer 136 (7): 71501-1–9. [DOI: 10.1115/1.4026808](http://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1829850).
3. Yang, Z., Peng, X. F., & Ye, P., 2008. Numerical and experimental investigation of two phase flow during boiling in a coiled tube. International Journal of Heat and Mass Transfer, 51(5-6), 1003–1016. (http://doi.org/10.1016/j.ijheatmasstransfer.2007.05.025)



