%Testing script for Stefan Problem
%Mahdi Nabil, 2015-12-19

%Some constants
g         = 9.81;        %m/s^2

%Fluid material properties (iso-butane 25 C)
rho_L     = 550.6;                     %kg/m^3
rho_V     = 9.12;                      %kg/m^3
k_L       = 0.0892;                    %W/m-K
mu_L      = 2.74E-7 * rho_L;           %kg/m-s
c_L  	  = 2446;                      %J/kgK
a_L	  = k_L/(rho_L*c_L);           %m^2/s
L	  = 329365;                    %J/kg

%Boundary Conditions
Tsat      = 298;                       % Saturated temperature [K]
Tw        = 293;                       % Wall temperature [K]

%Read in data from file:
D       	  = load('LiquidAccumulation.dat');
t        	  = D(:,1);                    %s
dt                = D(:,2);                    %s
delta_sim         = D(:,3);                    %m

%Compare with analytical predictions
delta_an = sqrt((2.*t).*(a_L).*((0.5+(L./(c_L.*(Tsat-Tw)))).^-1));   %Analytical solution for interface position

%Compute relative integrated errror
delta_t_int_sim = sum(dt.*delta_sim);
delta_t_int_an  = sum(dt.*delta_an);

%Compare result:
disp(sprintf('Final Interface position simulation: %g m', delta_sim(end)));
disp(sprintf('Final Interface position analytical: %g m', delta_an(end)));
disp(sprintf('Relative integrated error: %g', abs(delta_t_int_sim - delta_t_int_an)/delta_t_int_an ));
