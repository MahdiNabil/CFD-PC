%Testing script for Nusselt falling film wall heat flux results
%Alex Rattner, 2015-12-16

%Some constants
g         = 9.81;        %m/s^2
DeltaT    = 5;           %K

%Fluid material properties (iso-butane 25 C)
rho_L     = 550.6;                     %kg/m^3
rho_V     = 9.12;                      %kg/m^3
k_L       = 0.0892;                    %W/m-K
mu_L      = 2.74E-7 * rho_L;           %kg/m-s


%Read in data from file:
D         = load('WallHeatFlux.dat');
t         = D(:,1);                    %s
dt        = D(:,2);                    %s
q_w_sim   = D(:,3);                    %W/m^2
Re_f_sim  = D(:,4);                    %-

%Compare with analytical predictions
delta_an  = (3*mu_L^2/(4*rho_L*(rho_L-rho_V)*g))^(1/3) * Re_f_sim.^(1/3);   %Smooth falling film result
q_w_an    = DeltaT*k_L./delta_an;

%Compare result:
q_avg_an  = sum(dt.*q_w_an)./sum(dt);
q_avg_sim = sum(dt.*q_w_sim)./sum(dt);
disp(sprintf('Q_avg simulation: %g W/m^2', q_avg_sim));
disp(sprintf('Q_avg analytical: %g W/m^2', q_avg_an));
disp(sprintf('Relative error: %g', abs(q_avg_sim-q_avg_an)/abs(q_avg_an) ));
