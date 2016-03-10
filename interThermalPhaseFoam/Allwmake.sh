#!/bin/bash
 
#Build script for interThermalPhaseChangeFoam
wmake libso incompressibleTwoPhaseThermalMixture
wmake interThermalPhaseChangeFoam
wmake libso twoPhaseTransportProperties
