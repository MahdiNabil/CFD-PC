#!/bin/bash

#Cleanup script
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
cleanTimeDirectories

rm -r processor*
rm -r dynamicCode
rm -r 0/*
rm constant/polyMesh/cellZones
rm constant/polyMesh/faces
rm constant/polyMesh/faceZones
rm constant/polyMesh/neighbour
rm constant/polyMesh/owner
rm constant/polyMesh/points
rm constant/polyMesh/pointZones
rm -r constant/polyMesh/sets
rm Bubble_Condensation.dat
