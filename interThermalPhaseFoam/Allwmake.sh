#!/bin/bash
 
#Build script for interThermalPhaseChangeFoam
wmake libso incompressibleTwoPhaseThermalMixture
wmake interThermalPhaseChangeFoam
cd Libraries
./Allwmake.sh
