#!/bin/bash

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

#Cleanup script
rm -r processor*
rm -r dynamicCode
rm -r 0/*
rm constant/polyMesh/blockMeshDict
rm constant/polyMesh/faces
rm constant/polyMesh/neighbour
rm constant/polyMesh/owner
rm constant/polyMesh/points
rm WallHeatFlux.dat

cleanTimeDirectories
