#!/bin/bash
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh
checkMesh
cp -r A/* 0/ 
funkySetFields -time 0 -allowFunctionObjects
decomposePar -force
mpirun -np 4 interThermalPhaseChangeFoam -parallel
