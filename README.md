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
The installation process is straightforward. The bash script "Allwmake.sh" (https://github.com/MahdiNabil/CFD-PC/blob/master/interThermalPhaseFoam/Allwmake.sh), i.e. the only file that needs to be executed, is included in the main folder "interThermalPhaseFoam". Line 4 of this script builds "incompressibleTwoPhaseThermalMixture" library which will be dynamically linked to OpenFOAM. The next line will build the solver "interThermalPhaseChangeFoam".
The interested users only need to go through the below steps to be able to use this solver:
* Pull the "CFD-PC" folder from the Github repository in a directory of their choice: 
       * $ git pull https://github.com/MahdiNabil/CFD-PC.git
* Go into the directory of the main folder: 
       * $ cd ../CFD-PC/interThermalPhaseFoam
* Execute the provided bash script to build both the dynamic library and solver:
       * $ chmod +x Allwmake.sh
       * $ ./Allwmake.sh

## Tutorial cases
### Stefan (Horizontal Film Condensation)
This tutorial case demonstrates horizontal film condensation on an isothermal subcooled surface (Stefan problem). In this test case, the dynamic effects are relatively minor. Vapor condenses to form a liquid film on the top surface of an isothermal plate (at Tw) in a pure atmosphere. The analytical solution is readily available.

To run the case, call the `InitScript.sh` script. This will generate the mesh, initialize the fields, and begin the simulation. Results (film thickness) are logged to an output `.dat` file during the run. After completing or ending the run, results can be validated with the provide octave script. To run the validation check, call: `octave CheckStefan.m`. Errors are on the order of 10% early in the simulation due to the relative coarseness of the mesh, but reduce as film thickness grows. The case can be reset using the `cleanup.sh` script. An image of the output is presented below.

![Stefan Problem Example](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/12/Stefan_Snapshot.png)

### NusseltSmooth (Smooth Falling Film Condensation)
This tutorial demonstrates condensation of a smooth 2D falling film (Nusselt problem). A mechanistic analytical solution is available for this case. The vertical wall is isothermal, and subcooled. A small vane is employed to prevent inlet film waviness.

The case can be started using the `InitScript.sh` script. The script uses [swak4Foam](https://openfoamwiki.net/index.php/Contrib/swak4Foam) to initialize the film . By default, this script will start a 4-way parallel simulation. Results can be checked using the provided octave script (run `octave CheckNusselt.m`). Wall heat flux errors are generally < 1% after sufficient startup time. Representative film profile and temperature field results are presented below.

![Nusselt Smooth Problem Example](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/12/NusseltSmooth_Snapshot.png)

* Wavy Laminar Nusselt Falling Film Condensation: Smooth falling films on vertical plates are inherently unstable, and wavy behavior is initiated at finite Reynolds numbers. Waves tend to generate thin film regions with reduced heat-transfer resistance, yielding increased condensation rates. 
* Two-dimensional Nucleate Boiling in a Cavity: This test case shows the process of nucleate boiling in a single cavity due to vapor bubble growth and detachment from the heated bottom surface. 
* Bubble Condensation: This test case represents the phase change (shrinkage) of a vapor bubble condensing as it rises in a column of liquid water (due to Buoyancy force).

## Algorithm
At the beginning, the solver loads the mesh, reads in fields and boundary conditions, and initializes submodels for two-phase fluid properties, turbulence (if selected), and the phase-change model. During this initialization stage, the phase-change model constructs the graph of mesh-cell connectivity used to identify interface cells. The main solver loop is then initiated. First, the time step is
dynamically modified to ensure numerical stability. Next, the two-phase fluid mixture properties and turbulence quantities are updated. The phase-change model is then evaluated. The discretized phase-fraction equation is then solved for a user-defined number of subtime steps (typically 2–3) using the multidimensional universal limiter with explicit solution solver. This solver is included in the OpenFOAM library, and performs conservative solution of hyperbolic convective transport equations with defined bounds (0 and 1 for α1). Once the updated phase field is obtained, the program enters the pressure–velocity loop, in which p and u are corrected in an alternating fashion. The process of correcting the pressure and velocity fields in sequence is known as pressure implicit with splitting of operators (PISO). In the OpenFOAM environment, PISO is repeated for multiple iterations at each time step. This process is referred to as merged PISO- semi-implicit method for pressure-linked equations (SIMPLE), or the pressure-velocity loop (PIMPLE) process, where SIMPLE is an iterative pressure–velocity solution algorithm for steady flows. The PIMPLE continues for a user specified number of iterations.Finally, the thermal energy transport subsection is entered. First, the enthalpy field is reevaluated from the temperature field. Then for a user-defined number of steps (typically 2–3), alternating correction of the energy equation and update
of the temperature field (T(h)) are performed. This process is employed because temperature and enthalpy are coupled, but are not directly proportional for many fluids, and thus cannot be solved together in a single system of equations. The main solver loop iterates until program termination. A summary of the simulation algorithm is presented below:
* Flow Simulation Algorithm Summary:
  * Initialize simulation data and phase change model 
  * WHILE t<t_end DO
  * 1. Update delta_t for stability
  * 2. Update fluid and turbulence properties
  * 3. Update phase change model 
  * 4. Phase equation sub-cycle
  * 5. DO PIMPLE
    * 1. Form u equation
    * 2. PISO
        * 1. Obtain and correct face fluxes
        * 2. Solve p-Poisson equation
        * 3. Correct u
  * 6. LOOP
  * 7. Update h(T)
  * 8. DO Energy Loop
    * 1. Solve h equation
    * 2. Update T(h)
  * 9. LOOP
  * LOOP
  
Two sample tutorial cases, i.e. Horizontal film condensation and Smooth Nusselt falling film condensation are validated versus the available analytical solutions in the literature with less than 2% error. The corresponding MATLAB scripts included in the aforementioned tutorial cases folders are CheckStefan.m and CheckNusselt.m, respectively. 

## Phase Change Models
A number of  phase change models are included with the solver, and are described below:
* **HiLoRelaxed** – An improved version of the model of Rattner and Garimella (2014) that determines the phase change heat sources so that interface cells recover the saturation temperature at each time step. This model performs a graph scan over mesh cells, and applies phase change on the two-cell thick interface layer about user-specified threshold values of α1. Different high and low threshold values for condensation and evaporation, respectively, can be specified, which has been found to reduce numerical smearing of the interface. Numerical under-relaxation of the phase change rate is supported, which can improve numerical stability.
* **HiLoRelaxedSplit** – A modified version of the above model, which splits the liquid and vapor portions of the dilatation rate, and applies them on the respective sides of the interface (Rattner, 2015). This approach yields better conservation of the two phases, and reduces smearing of the interface during evaporation.
* **HiLoNoPCV** – A modified version of HiLoRelaxed that sets the dilatation rate to 0. This model can yield accurate results with reduced CFL time step constraints for cases without significant vapor-phase effects on phase change (e.g., falling film condensation in a quiescent vapor medium).
* **HiLoNoPCVAlpha1** – A modified version of HiLoNoPCV that sets both the dilatation rate and phase fraction source term to 0. This model can also yield accurate results with reduced CFL time step constraints for cases without significant vapor-phase effects on phase change 
* **Yang** – An implementation of the empirical rate parameter model of Yang et al. (2008).

## Example applications
* Progression of dropwise condensation for a moderate surface tension fluid
![Dropwise condensation, high sigma](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/04/DropwiseCond_Sigma_1E-3sm.gif)
* Progression of dropwise condensation for a low surface tension fluid (transition to film condensation)
![Dropwise condensation, low sigma](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/04/DropwiseCond_Sigma_1E-4sm.gif)
* Simulation of nucleate boiling on a structured surface
![Nucleate boiling](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/04/NucBoilArray-Merge.gif)
* Simulation of vapor absorption into solution flowing over rectangular cooled tubes (film colored by concentration)
![Falling film absorption](http://sites.psu.edu/mtfe/wp-content/uploads/sites/23865/2015/04/SaltTubesBenchmarksm.gif)
- For more information, please visit:  http://sites.psu.edu/mtfe/simulation-of-phase-change-flows/

## Contribute
MTFE welcomes collaboration with other investigators studying phase-change flows. Please [contact us](mailto:Alex.Rattner@psu.edu) if you are interested in expanding the solver or find bugs to correct. We can also provide limited support (on a case-by-case basis) or consulting servies.

## License
Information about the GPL V3.0 license is included in the "GNU License" file in the main CFD-PC directory:
https://github.com/MahdiNabil/CFD-PC/blob/master/GNU%20Licence

## References
1. Rattner, A.S., 2015. Single-pressure absorption refrigeration systems for low-source-temperature applications. Ph.D. Thesis, Georgia Institute of Technology, Atlanta, GA. [Link](https://smartech.gatech.edu/handle/1853/53912).
2. Rattner, A.S., Garimella, S., 2014. Simple mechanistically consistent formulation for volume-of-fluid based computations of condensing flows. Journal of Heat Transfer 136 (7): 71501-1–9. [DOI: 10.1115/1.4026808](http://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1829850).
3. Yang, Z., Peng, X. F., & Ye, P., 2008. Numerical and experimental investigation of two phase flow during boiling in a coiled tube. International Journal of Heat and Mass Transfer, 51(5-6), 1003–1016. (http://doi.org/10.1016/j.ijheatmasstransfer.2007.05.025)



