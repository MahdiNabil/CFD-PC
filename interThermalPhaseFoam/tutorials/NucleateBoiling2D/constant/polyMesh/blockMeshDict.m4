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
define(R_c, 25E-6);       //Cavity radius
define(H_c, 50E-6);      //Cavity Height
define(R_m, 300E-6);      //Main domain radius (width/2)
define(H_m, 450E-6);      //Main domain height
//Discretization params
define(N_c_x, 8);        //Cavity x-direction
define(N_c_y, 15);        //Cavity y-direction
define(N_m_x, 90);       //Main-zone x-direction
define(N_m_y, 136);       //Main-zone y-direction

convertToMeters 1;

vertices        
(
	//Bot Face
	(0                    0                  -1E-6)       //00
	(R_c                  0                  -1E-6)       //01
	(0                    H_c                -1E-6)       //02
	(R_c                  H_c                -1E-6)       //03
	(R_m                  H_c                -1E-6)       //04
	(0                    calc(H_c+H_m)      -1E-6)       //05
	(R_c                  calc(H_c+H_m)      -1E-6)       //06
	(R_m                  calc(H_c+H_m)      -1E-6)       //07
	//Top Face
	(0                    0                  1E-6)        //08
	(R_c                  0                  1E-6)        //09
	(0                    H_c                1E-6)        //10
	(R_c                  H_c                1E-6)        //11
	(R_m                  H_c                1E-6)        //12
	(0                    calc(H_c+H_m)      1E-6)        //13
	(R_c                  calc(H_c+H_m)      1E-6)        //14
	(R_m                  calc(H_c+H_m)      1E-6)        //15
);

blocks          
(
	hex ( 0  1  3  2  8  9 11 10) (N_c_x N_c_y 1) simpleGrading (1 1 1)       //00
	hex ( 2  3  6  5 10 11 14 13) (N_c_x N_m_y 1) simpleGrading (1 1 1)      //01
	hex ( 3  4  7  6 11 12 15 14) (roundcalc(N_m_x-N_c_x) N_m_y 1) simpleGrading (1 1 1)     //02
);

edges           
(
);

patches         
(
	wall Cavity
	(
		( 0  1  9  8)
	)

	wall Base
	(
		( 3  4 12 11)
		( 1  3 11  9)
	)

	patch FarStream
	(
		( 4  7 15 12)
	)

	patch Top
	(
		( 5 13 14  6)
		( 6 14 15  7)
	)

	patch Axis
	(
		( 0  8 10  2)
		( 2 10 13  5)
	)

	empty FrontNBack
	(
		( 0  2  3  1)
		( 2  5  6  3)
		( 3  6  7  4)
		( 8 10 11  9)
		(10 13 14 11)
		(11 14 15 12)
	)

);

mergePatchPairs 
();

// ************************************************************************* //
