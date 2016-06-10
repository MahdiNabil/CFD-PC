/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.7.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.1;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//For m4 preprocessing
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(roundcalc, [esyscmd(perl -e 'printf ( int( ($1) + 0.5 ) )')])

//Some dimensions here
//define(R_m, 1250E-6);      //Main domain radius (width/2)
//define(H_m, 1250E-6);      //Main domain height
define(R_m, 1250E-5);      //Main domain radius (width/2)
define(H_m, 1250E-5);      //Main domain height
//Discretization params
define(N_m_x, 50);       //Main-zone x-direction
define(N_m_y, 50);       //Main-zone y-direction

convertToMeters 1;

vertices        
(
	//Bot Face
	(0                    0                  0)       //00
	(R_m                  0                  0)       //01
	(R_m                  H_m                0)       //02
	(0                    H_m                0)       //03
	//Top Face
	(0                    0                  H_m)       //04
	(R_m                  0                  H_m)       //05
	(R_m                  H_m                H_m)       //06
	(0                    H_m                H_m)       //07
);

blocks          
(
	hex ( 0  1  2  3  4  5  6  7) (N_m_x N_m_y N_m_y) 
	simpleGrading 
	(
		(
			(0.3 0.7 1) //10% y-dir, 90% cells, expansion = 1
			(0.7 0.3 14)
		)
		(
			(0.3 0.7 1) //10% y-dir, 90% cells, expansion = 1
			(0.7 0.3 14)
		)
		(
			(0.3 0.7 1) //10% y-dir, 90% cells, expansion = 1
			(0.7 0.3 14)
		)

	)
);

edges           
(
);
/*
boundary         
(
	Base
	{
		type            wall;
		(
			( 0  1  5  4)
		);
	}

	FarStream
	{
		type            patch;
		(
			( 1  2  6  5)
		);
	}

	Top
	{
		type     patch;
		(
			( 2  3  7  6)
		);
	}

	Axis
	{
		type     patch;
		(
			( 0  4  7  3)
		);
	}

	FrontNBack
	{
		type     FrontNBack;
		(
			( 0  3  2  1)
			( 4  5  6  7)
		);
	}

);
*/

patches         
(


	patch FarStream
	(
		( 1  2  6  5)
		( 4  5  6  7)
		( 2  3  7  6)
	)


	symmetryPlane sp1
	(
		( 0  1  5  4)
	)


	symmetryPlane sp2
	(
		( 0  4  7  3)
	)

	wall Base
	(
		( 0  3  2  1)
	)
);



mergePatchPairs 
();

// ************************************************************************* //
