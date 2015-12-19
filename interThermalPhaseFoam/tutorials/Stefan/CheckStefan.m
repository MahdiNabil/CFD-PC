%Testing script for Stefan Problem
%Mahdi Nabil, 2015-12-19

%Some constants
g         = 9.81;        %m/s^2

%Fluid material properties (iso-butane 25 C)
rho_L     = 550.6;                     %kg/m^3
rho_V     = 9.12;                      %kg/m^3
k_L       = 0.0892;                    %W/m-K
mu_L      = 2.74E-7 * rho_L;           %kg/m-s
c_L	  	  = 2446; 					   %J/kgK
a_L		  = k_L/(rho_L*c_L);           %m^2/s
L		  = 329365;                    %J/kg

%Boundary Conditions
Tsat	  = 298;                       % Saturated temperature [K]
Tw		  =	293;                       % Wall temperature [K]

%Read in data from file:
D       		  = load('LiquidAccumulation.dat');
t        		  = D(:,1);                    %s
Vol_Liquid        = D(:,2);                    %m^3
Vol_Total         = D(:,3);                    %m^3

%Compare with analytical predictions
IP_an 	  = sqrt((2.*t).*(a_L).*((0.5+(L./(c_L.*(Tsat-Tw)))).^-1));   %Analytical solution for interface position

%Compare result:
dx		  		      = 0.002; 						%m
dz		 			  = 0.0005;						%m
Area    			  = dx*dz;						%m^2
IP_sim			      = Vol_Liquid/Area;		    %m

disp(sprintf('Intergface position simulation: %g m', IP_sim));
disp(sprintf('Intergface position analytical: %g m', IP_an));
disp(sprintf('Relative error: %g', abs(IP_sim-IP_an)/abs(IP_an) ));
