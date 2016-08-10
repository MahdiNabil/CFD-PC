#!/bin/bash
unset FOAM_SIGFPE

blockMesh

cp -r A/* 0/

setFields

interThermalPhaseChangeFoam

