#!/bin/bash
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh
checkMesh
mkdir 0
cp -r A/* 0/ 
setFields
decomposePar -force
mpirun -np 12 interThermalPhaseChangeFoam -parallel >Stefan.log&
