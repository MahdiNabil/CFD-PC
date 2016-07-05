#!/bin/bash
 
#Clean script for interThermalPhaseChangeFoam
wclean incompressibleTwoPhaseThermalMixture
wclean interThermalPhaseChangeFoam
cd Libraries
./Allwclean.sh
