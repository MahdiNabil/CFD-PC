%Testing script for parasitic current heat transfer validation
%Alex Rattner, 2016-08-09

%Some constants
L         = 10E-6;       %m, Length of domain in third dimension
D_drop    = 50E-6;       %m, Diameter of droplet
W         = 100E-6;      %m, Length of domain
DeltaT    = 1;           %K

%Fluid material properties (iso-butane 25 C)
k         = 0.1;                     %W/m-K

%Read in data from file:
D         = load('DataSummary.dat');
%Trim the first 0.005 s to block out startup effects
ind       = find(D(:,1) > 0.02, 1, 'first');
D         = D(ind:end,:);
%Get out data entries
t         = D(:,1);                    %s
dt        = D(:,2);                    %s
U         = D(:,3);                    %m/s
Q         = D(:,4);                    %W

%Analytical model, shape factor
Q_an      = (1/4) * (2*pi*L/log(1.08*W/D_drop)) * k * DeltaT;

%Sim result:
Q_sim     = sum( dt.*Q )./sum(dt);

%Compare with analytical predictions
disp(sprintf('Average parasitic current velocity: %g m/s', mean(U) ));
disp(sprintf('Q_avg simulation: %g W', Q_sim ));
disp(sprintf('Q_avg analytical: %g W', Q_an ));
disp(sprintf('Relative error:   %g %%', 100*(Q_sim-Q_an)/Q_an ));

