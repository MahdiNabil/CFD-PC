/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    MeshGraph

Description
    Class for constructing graph of mesh, and finding interface cells

\*---------------------------------------------------------------------------*/

#ifndef __MeshGraph
#define __MeshGraph

//include openfoam stuff
#include "fvCFD.H"

//Stuff from stdlib
#include <vector>

namespace Foam
{


class MeshGraph
{
    public:
        //- Struct for sharing interface cell-face pairs:
        struct CellFacePair
        {
            label   f;        //Index of the face           
            label   c1, c2;   //Indices of the two cells sharing the face
            scalar  v1, v2;   //Values of the cells sharing the face
        };

        //- Construct from existing mesh
        MeshGraph( const fvMesh& Min, const volScalarField& Fin );
        
        //- Default destructor
        ~MeshGraph();
        
        //- Resets the graph values / interface states
        void Reset();
        
        //- Traverses graph and finds interface cells (those containing the
        //  intVal)
        void GetInterfaceCells( labelList& IntCells, const scalar& intVal);
        
        //- Like above, but returns cells on both sides of the interface
        void GetDoubleInterfaceCells
        (
            labelList& IntCells,
            const scalar& intVal
        );
        
        //- Another double layer method, but it actually returns the cells &
        //  data that are on the double layer:
        void GetInterfaceCellFacePairs
        (
            std::vector<CellFacePair>& IntCellFacePairs, 
            const scalar& intVal
        );


    private:
        //Connection to driving scalar field
        const volScalarField& F;

        //Graph node
        struct MeshGraphCell
        {
            //Contents
            scalar val;     //Value of current cell
            vector c;                   //Cell center coordinates
        };
        struct MeshGraphFace            //Face connecting two cell centers
        {
            label  c1, c2;              //Cell indices
            vector n;                   //Face normal (points out of owner cell)
            vector c;                   //Face center
        };

        std::vector <MeshGraphCell> Cells;
        std::vector <MeshGraphFace> Faces;  
};

} // End namespace Foam


#endif
