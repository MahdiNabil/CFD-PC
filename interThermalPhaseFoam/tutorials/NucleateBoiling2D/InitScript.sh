#!/bin/bash
unset FOAM_SIGFPE
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh

#Need to check if makeAxialMesh is available
HAS_MAKEAXIALMESH=`command -v makeAxialMesh`
if [ ! $HAS_MAKEAXIALMESH ]
then
	echo "Need to install makeAxialMesh"
	mkdir -p $HOME/Downloads
	CURDIR=`pwd`
	cd $HOME/Downloads
	svn checkout svn://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder_2.0/utilities/mesh/manipulation/MakeAxialMesh
	cd MakeAxialMesh
	./Allwmake
	cd $CURDIR
fi

makeAxialMesh -overwrite
collapseEdges -overwrite
cp -r A/* 0/
mv 0/alpha1.org 0/alpha1
setFields
decomposePar -force
mpirun -np 12 interThermalPhaseChangeFoam -parallel >NucleateBoiling2D.log&
