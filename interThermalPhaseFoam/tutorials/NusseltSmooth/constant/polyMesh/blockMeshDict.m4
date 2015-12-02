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
define(dF, 16.72E-5);       //Film region thickness
define(thV, 5E-5);        //Thickness of the inlet vane
define(hV, 2.5E-4);        //Height of the inlet vane
define(hB, 2.5E-2);         //Bulk section height
define(thB, 2.5E-3);        //Thickness of the bulk section
//Discretization params
define(hX, 3.3333E-6);         //Normal direction resolution
define(hXB, 1.1111E-5);		//Normal direction resolution in the bulk
define(hY, 1.1111E-5);           //Length-direction resolution

convertToMeters 1;

vertices        
(
	//Back face
	(0                    0                  0)           //00
	(dF                   0                  0)           //01
	(calc(dF+thV)         0                  0)           //02
	(calc(dF+thV+thB)     0                  0)           //03
	(0                    hB                 0)           //04
	(dF                   hB                 0)           //05
	(calc(dF+thV)         hB                 0)           //06
	(calc(dF+thV+thB)     hB                 0)           //07
	(0                    calc(hB+hV)        0)           //08
	(dF                   calc(hB+hV)        0)           //09
	//Front face
	(0                    0                  dF)          //10
	(dF                   0                  dF)          //11
	(calc(dF+thV)         0                  dF)          //12
	(calc(dF+thV+thB)     0                  dF)          //13
	(0                    hB                 dF)          //14
	(dF                   hB                 dF)          //15
	(calc(dF+thV)         hB                 dF)          //16
	(calc(dF+thV+thB)     hB                 dF)          //17
	(0                    calc(hB+hV)        dF)          //18
	(dF                   calc(hB+hV)        dF)          //19
	//Extra corner block for U domain
	(calc(dF+thV)         calc(hB+hV)        0)           //20
	(calc(dF+thV+thB)     calc(hB+hV)        0)           //21
	(calc(dF+thV)         calc(hB+hV)        dF)          //22
	(calc(dF+thV+thB)     calc(hB+hV)        dF)          //23

);

blocks          
(
	hex ( 0  1  5  4 10 11 15 14) (roundcalc(dF/hX)	    roundcalc(hB/hY)    1) simpleGrading (1 1 1)       //00
	hex ( 1  2  6  5 11 12 16 15) (roundcalc(thV/hX)    roundcalc(hB/hY)    1) simpleGrading (1 1 1)       //01
	hex ( 2  3  7  6 12 13 17 16) (roundcalc(thB/hXB)   roundcalc(hB/hY)    1) simpleGrading (5 1 1)       //02
	hex ( 4  5  9  8 14 15 19 18) (roundcalc(dF/hX)	    roundcalc(hV/hY)    1) simpleGrading (1 1 1)       //03
	hex ( 6  7 21 20 16 17 23 22) (roundcalc(thB/hXB)   roundcalc(hV/hY)    1) simpleGrading (5 1 1)       //04

);

edges           
(
);

patches         
(
	wall VertWall 
	(
		( 0 10 14  4)
		( 4 14 18  8)
	)

	wall Vane
	(
		( 5 15 16  6)
		( 5  9 19 15)
		( 6 16 22 20)
	)

	patch FreeStream
	(
		( 3  7 17 13)
		( 7 21 23 17)
	)

	patch Inlet
	(
		( 8 18 19  9)
	)

	patch Top
	(
		(20 22 23 21)
	)

	patch Bottom
	(
		( 0  1 11 10)
		( 1  2 12 11)
		( 2  3 13 12)
	)

	empty Sides
	(
		( 0  4  5  1)
		( 1  5  6  2)
		( 2  6  7  3)
		( 4  8  9  5)
		(10 11 15 14)
		(11 12 16 15)
		(12 13 17 16)
		(14 15 19 18)
		( 6 20 21  7)
		(16 17 23 22)
	)	
);

mergePatchPairs 
();

// ************************************************************************* //
