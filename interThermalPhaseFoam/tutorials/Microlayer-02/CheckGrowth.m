%water nucleate boiling radius growth
%Tsat = 573.15K (300C)
clear;
tstep = 10000001;
t = linspace(0, 20000, tstep);                 % [s] time

k = 0.548;            % [W/m-k] liq. thermal conductivity
lambda = 1405E3;      % [J/kg] latent heat of evap.
alpha = 1.312E-7;      % [m^2/s] liq. thermal diffusivity
nu = 1.21E-7;         % [m^2/s] kinematic viscosity
delT = 0.5;              % [K] superheat of wall
rho_v = 46.2;         % [kg/m^3] density of vapor
rho_l = 712.3;        % [kg/m^3] density of liquid
cp_l = 5746;          % [J/kg*K] specific heat of liquid

%analytical radius of bubble [Strenge 1961]
R = sqrt(3) * (2*k*delT)/(lambda*rho_v*sqrt(pi*alpha)) * sqrt(t);  %[m] 

%analytical radius of bubble Thermopedia [Cooper 1969]
Ja = ((rho_l*cp_l*delT)/(rho_v*lambda));
Pr = nu/alpha;
%R = 2.5 * (Ja/Pr) * sqrt(nu*t); 

% Load imulated data
D = load('Nucleate_Boiling.dat');

%load init and final simulated radius 
%D(:,6) = ((D(:,4).*6)./pi).^(1/3);
initRad = D(1,6);

lastRad = D(57, 6);

simGrowth = lastRad - initRad;

%calibrate analytical start/end time to simulation
for i = 1:columns(R)
	if R(i) >= initRad
    Rstartt = t(i);
    disp(R(i))
		break;
	endif
end
Rendt = t(i + 50); %0.1 seconds later

%analytical start/end radius
expR_init = R(i);
expR_end = R(i+50);
expR_growth = expR_end - expR_init;



disp(sprintf('Simulation Initial Radius %g m', initRad));
disp(sprintf('Simulation Last Radius %g m', lastRad));
disp(sprintf('Simulated growth in 0.1 secs %g m', simGrowth));

disp(sprintf('Analytical start time %g s', Rstartt));
disp(sprintf('Analytical end time %g s', Rendt));

disp(sprintf('Analytical Initial Radius %g m', expR_init));
disp(sprintf('Analytical Last Radius %g m', expR_end));
disp(sprintf('Analytical growth in 0.1 secs %g m', expR_growth));

disp(sprintf('     Relative error: %g percent', 100*(simGrowth-expR_growth)/(expR_growth) ));

plot(t,R,D(:,1)+Rstartt,D(:,6))








