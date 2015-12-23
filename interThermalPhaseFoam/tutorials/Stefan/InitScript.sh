#!/bin/bash
blockMesh
checkMesh
mkdir -p 0
cp -r A/* 0/
interThermalPhaseChangeFoam
