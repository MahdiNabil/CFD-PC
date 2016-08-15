%Testing script for Rising condensing bubble case
%Alex Rattner, 2015-12-16

%Some constants
g         = 9.81;        %m/s^2
DeltaT    = 1;           %K

%Geometry info
theta_wedge = 5;         %deg, Axisymmetric wedge angle
Vol_tot     = 1.5688e-09; %m^3, domain volume

%Fluid material properties (iso-butane 25 C)
rho_L     = 900;                     %kg/m^3
rho_V     = 10;                      %kg/m^3
k_L       = 1;                       %W/m-K
cp_L      = 2000;                    %J/kg-K
mu_L      = 8.73876E-6 * rho_L;        %kg/m-s
mu_G      = 5.0E-5;                  %kg/m-s
Pr_L      = mu_L*cp_L/k_L;           %-, Liquid Prandtl number
U_domain  = 0;                       %m/s


%Read in data from file:
D         = load('Bubble_Condensation.dat');
%Trim the first 0.02 s to block out startup effects
ind       = find(D(:,1) > 0.15, 1, 'first');
D         = D(ind:end,:);
%Get out data entries
t         = D(:,1);                    %s
Q         = D(:,2);                    %W
VF        = D(:,3);                    %-, bubble void fraction
U_bub     = D(:,4) + U_domain;         %m/s

%Bubble properties:
Vol_bub   = VF*Vol_tot*(360/theta_wedge);
D_bub     = ( (6/pi)*Vol_bub ).^(1/3);
Re_bub    = rho_L*abs(U_bub).*D_bub/mu_L;
A_bub     = pi*D_bub.^2;
H_sim     = (360/theta_wedge)*Q ./ (A_bub * DeltaT);

%Analytical model, Ranz and Marshall (1952)
Nu_an     = 2 + 0.6*Re_bub.^(0.5) * Pr_L^0.33;
H_an      = Nu_an*k_L./D_bub;

%Bubble rise velocity
%Model from Hadamard (1911) and Rybczynski (1911)
k1        = mu_G/mu_L;
U_bub_an  = (g*D_bub.^2 *(rho_L-rho_V)/(6*mu_L))*((1+k1)/(2+3*k1));


%Compare with analytical predictions
disp(sprintf('H_avg simulation: %g W/m^2-K', mean(H_sim) ));
disp(sprintf('H_avg analytical: %g W/m^2-K', mean(H_an) ));
disp(sprintf('Average simulation error: %g', mean( abs(H_sim-H_an)./abs(H_an) ) ));

