#!/bin/bash
blockMesh
checkMesh
mkdir -p 0
cp -r A/* 0/ 
setFields
decomposePar -force
mpirun -np 6 interThermalPhaseChangeFoam -parallel
