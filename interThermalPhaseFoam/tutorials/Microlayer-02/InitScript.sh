#!/bin/bash
unset FOAM_SIGFPE
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh

mkdir -p 0
cp -r A/* 0/
funkySetFields -time 0 -allowFunctionObjects

decomposePar -force
mpirun -np 6 interThermalPhaseChangeFoam -parallel
