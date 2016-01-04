%Testing script for Nusselt falling film wall heat flux results
%Alex Rattner, 2015-12-16

%Some constants
g         = 9.81;        %m/s^2
DeltaT    = 5;           %K

%Fluid material properties (iso-butane 25 C)
rho_L     = 550.6;                   %kg/m^3
rho_V     = 9.12;                    %kg/m^3
k_L       = 0.0892;                  %W/m-K
mu_L      = 2.74E-7 * rho_L;         %kg/m-s


%Read in data from file:
D         = load('WallHeatFlux.dat');
%Trim the first 0.1 s for stability
ind       = find(D(:,1) < 0.5, 1, 'last');
D         = D(ind:end,:);
%Get out data entries
t         = D(:,1);                    %s
dt        = D(:,2);                    %s
q_w_sim   = D(:,3);                    %W/m^2
Re_f_sim  = D(:,4);                    %-

%Compare with correlations
H_anA     = k_L*(mu_L^2 / (g*rho_L*(rho_L-rho_V)))^(-1/3) * 0.82* Re_f_sim.^-0.22; %Edwards et al. (1979)
q_w_anA   = H_anA * DeltaT;
H_anB     = k_L*(mu_L^2 / (g*rho_L*(rho_L-rho_V)))^(-1/3) * 1.76 * Re_f_sim.^(-1/3); % Fujita and Ueda (1978)
q_w_anB   = H_anB * DeltaT;

%Compare result:
q_avg_anA  = sum(dt.*q_w_anA)./sum(dt);
q_avg_anB  = sum(dt.*q_w_anB)./sum(dt);
q_avg_sim = sum(dt.*q_w_sim)./sum(dt);
disp(sprintf('Q_avg simulation: %g W/m^2', q_avg_sim));
disp(sprintf('Q_avg analytical: %g W/m^2 (Edwards et al., 1979)', q_avg_anA));
disp(sprintf('     Relative error: %g', abs(q_avg_sim-q_avg_anA)/abs(q_avg_anA) ));
disp(sprintf('Q_avg analytical: %g W/m^2 (Fujita and Ueda, 1978)', q_avg_anB));
disp(sprintf('     Relative error: %g', abs(q_avg_sim-q_avg_anB)/abs(q_avg_anB) ));
