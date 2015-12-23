#!/bin/bash

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

#Cleanup script
rm -r dynamicCode
rm -r 0/*
cleanTimeDirectories
rm constant/polyMesh/blockMeshDict
rm constant/polyMesh/faces
rm constant/polyMesh/neighbour
rm constant/polyMesh/owner
rm constant/polyMesh/points
rm constant/polyMesh/cellZones
rm constant/polyMesh/faceZones
rm constant/polyMesh/pointZones
rm -r constant/polyMesh/sets
rm Nucleate_Boiling.dat
