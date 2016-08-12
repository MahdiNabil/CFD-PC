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
    Implementation: class for constructing graph of mesh, and finding interface
    cells

\*---------------------------------------------------------------------------*/

#include "MeshGraph.H"
#include <set>


// Constructor from mesh
MeshGraph::MeshGraph( const fvMesh& Min, const volScalarField& Fin ) :
F(Fin)
{
    // Start by allocating space in the Node vector
    Cells.resize( Min.nCells() );
    // Now populate nodes with centers
    forAll( Min.cells(), cI )
    {
        // Set position      
        Cells[cI].c = vector( Min.C()[cI] );
    }

    // Now allocate faces
    Faces.resize( Min.faceNeighbour().size() ); // number of connected faces
    // Now populate faces
    forAll( Min.faceNeighbour(), fI )
    {
        // Record connecting cells:
        Faces[fI].c1 = Min.faceOwner()[fI];
        Faces[fI].c2 = Min.faceNeighbour()[fI];
        Faces[fI].n  = (Min.Sf()[fI]/Min.magSf()[fI]);
        Faces[fI].c  = Min.Cf()[fI];
    }

}


// Default destructor
MeshGraph::~MeshGraph()
{
    // Clear out vector of nodes and faces
    Cells.clear();
    Faces.clear();
}


// Updates the cell values
void MeshGraph::Reset()
{
    int NC = Cells.size();
    for ( int i = 0; i < NC; i++ )
    {  Cells[i].val = F[i];  }
}

// Traverses graph and finds interface cells (those containing the intVal)
void MeshGraph::GetInterfaceCells(labelList& IntCells, const scalar& intVal)
{
    // Traverse through graph faces
    int N = Faces.size();

    // Initialize a set so we only deal with unique interface cells:
    std::set<label> IntCellsSet;

    for ( int i = 0; i < N; i++ )
    {
        // Get two connecting cells
        MeshGraphCell& C1 = Cells[ Faces[i].c1 ];
        MeshGraphCell& C2 = Cells[ Faces[i].c2 ];

        // Get values and difference
        scalar Val  = C1.val;
        scalar dVal = ( C2.val - Val ) + SMALL;

        // Check for no change (potential issue when using +SMALL)
        if (dVal == 0.0)
        {
            continue;
        }

        // Check if relative distance is oob
        scalar relDist = (intVal - Val)/dVal;
        if ( (relDist < 0) || (relDist > 1) )
        {
            continue;
        }

        // Check where the interface lies relative to the cell face
        // (on which side)
        vector intPos = C1.c + (C2.c - C1.c)*relDist;
        scalar side = (Faces[i].c - intPos) & Faces[i].n;

        // Mark if current cell contains interface
        if ( side >= 0 )
        {   
            IntCellsSet.insert( Faces[i].c1 );
        }

        // Mark if neighbor cell has interface
        if ( side <= 0 )
        {
            IntCellsSet.insert( Faces[i].c2 );  
        }
    }

    // Now populate the list of interface cells
    for 
    (
        std::set<label>::iterator it = IntCellsSet.begin();
        it != IntCellsSet.end();
        ++it
    )
    {
        IntCells.append( *it );
    }
}

// Like above, but returns cells on both sides of the interface...
void MeshGraph::
GetDoubleInterfaceCells( labelList& IntCells, const scalar& intVal)
{
    // Traverse through graph faces
    int N = Faces.size();

    // Initialize a set so we only deal with unique interface cells:
    std::set<label> IntCellsSet;

    for ( int i = 0; i < N; i++ )
    {
        // Get two connecting cells
        MeshGraphCell& C1 = Cells[ Faces[i].c1 ];
        MeshGraphCell& C2 = Cells[ Faces[i].c2 ];

        // Get values and difference
        scalar Val  = C1.val;
        scalar dVal = ( C2.val - Val ) + SMALL;

        // Check for no change (potential issue when using +SMALL)
        if (dVal == 0.0)
        {
            continue;  
        }

        // Check if relative distance is oob
        scalar relDist = (intVal - Val) / dVal;
        if ( (relDist < 0) || (relDist > 1) )
        {  
            continue;
        }

        // Mark both cells as on the interface
        IntCellsSet.insert( Faces[i].c1 );
        IntCellsSet.insert( Faces[i].c2 );
    }

    // Now populate the list of interface cells
    for
    (
        std::set<label>::iterator it = IntCellsSet.begin();
        it != IntCellsSet.end();
        ++it 
    )
    {
        IntCells.append( *it );  
    }
}


// Another double layer method, but it actually returns the cells & data that
// are on the double layer:
void MeshGraph::GetInterfaceCellFacePairs
(
    std::vector<CellFacePair>& IntCellFacePairs,
    const scalar& intVal
)
{
    // Traverse through graph faces
    int N = Faces.size();

    for ( int i = 0; i < N; i++ )
    {
        // Get two connecting cells
        MeshGraphCell& C1 = Cells[ Faces[i].c1 ];
        MeshGraphCell& C2 = Cells[ Faces[i].c2 ];

        // Get values and difference
        scalar Val  = C1.val;
        scalar dVal = ( C2.val - Val ) + SMALL;

        // Check for no change (potential issue when using +SMALL)
        if (dVal == 0.0)
        {
            continue;
        }

        // Check if relative distance is oob
        scalar relDist = (intVal - Val) / dVal;
        if ( (relDist < 0) || (relDist > 1) )
        {
            continue;
        }

        // OK, now this cell face pair is an interface pair, so add it to the
        // vector:
        // Populate struct - 
        CellFacePair curCFP;
        curCFP.f = i;
        // mark higher value as first cell:
        if ( C1.val >= C2.val)
        {       
            curCFP.c1 = Faces[i].c1;
            curCFP.c2 = Faces[i].c2;
            curCFP.v1 = C1.val;
            curCFP.v2 = C2.val;
        }
        else
        {
            curCFP.c1 = Faces[i].c2;
            curCFP.c2 = Faces[i].c1;
            curCFP.v1 = C2.val;
            curCFP.v2 = C1.val;
        }

        // Append struct:
        IntCellFacePairs.push_back(curCFP);
    }
}

