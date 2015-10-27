% Assuming this is a 90 degree contact angle problem
% All the properties need to be checked and modified later 
clc
%% properties of the phases
sigma = 1e-3; % [N/m] surface tension  
Tsat = 373; % [K] saturation temperature
rhol = 999.97; % [kg/m3] density of the liquid phase
rhov = 0.804; % [kg/m3] density of the vapor phase
Hlv = 226e4; % [J/kg] latent heat of vaporization
gamma = 1.32; % [] Cp/Cv  The value for water
MW = 18.0;   %[kg/kmol] Molar mass of water
Ru = 8314;   %J/kmol-K
Rg = Ru/MW;  %[J/k-kg] specific ideal gas constant
kl = 0.6; % [W/m-K] thermal conductivity of the liquid

%% constants and initial values
alpha = 0; % [] initial phase fraction
corr = 0.627 / 0.664; % [] correction factor for the heat flux relation
C1 = 1; % [] constant 1 for heat flux relation
C2 = 1; % [] constant 2 for heat flux relation
C3 = 1.1; % [] factor for rmin
Tw = 369; % [K] temperature of the surface
Tv = 373; % [K] temperature of the vapor phase
deltaT = Tv - Tw; % [K]
rmin = C3 * 2 * sigma * Tsat / deltaT / rhol / Hlv; % [m] minimum possible
% radius of the droplet; Form Rose (1998)

%% Derive integral for SGS heat transfer

R_MAX = 0.1; %m, some arbitrary big size

R_maxs = logspace( log10(rmin), log10(R_MAX), 1000 );
q_sgs = zeros(size(R_maxs));

for j = 1:length(q_sgs)
  rmax = R_maxs(j);
  qint = @(r) r.^(-2/3) .* (deltaT - (2*sigma*Tsat)./(r*rhol*Hlv)) ./ ... 
            ((C1*r/kl) + C2*corr*Tsat/(Hlv^2*rhov)*(gamma+1)/(gamma-1)*((Rg*Tsat)/(2*pi))^0.5) ;
    
  q_sgs(j) = (1/(3*rmax^(1/3))) * quad( qint, rmin, rmax);
end











%% Domain and time
siteD = 1e8; % [] site density from GlickMan (1971)
Lx = 1e-2; % [m]
Ly = Lx; % [m] redundant dimension
Lz = Lx; % [m] redundant dimension
dLx = Lx / sqrt(siteD); % [m]
dLy = dLx; % [m] redundant dimension
dLz = dLx; % [m]redundant dimension
rmax = dLx/2; % [m] maximum radius of a drop before coalescence
dt = 1e-19; % [s] timeStep


%% the setup 
t = 0; % [s] initial time
N = 5000; % number of time steps
time = zeros(1,N);
heatFlux = zeros(1,N);
rs = zeros(1,N);
r = rmin; % [m] radius of the droplets
n = siteD; % [] initial number of drops
% volAlpha = 0; % [m3] volume for alpha Generation
volDrop = 4/6*pi*r^3; % [m3] volume of a drop I'll do some geometry
vol = volDrop*n; % [m3] initial volume
Iter = 1; 
while t < N*dt
    % the multiplications and divisions by n are not optimized. It is just 
    % for clarity of the procedure
    if r < rmax
        qb = (deltaT - (2*sigma*Tsat)/(r*rhol*Hlv)) /... 
            ((C1*r/kl) + C2*corr*Tsat/(Hlv^2*rhov)*(gamma+1)/(gamma-1)*((Rg*Tsat)/(2*pi))^0.5)...
            * dt; % [J] heat through the base of one droplet From Rose (1998)
        qbn = qb * n; % [J] heat through all the drops
        vol = vol + qbn/(rhol*Hlv);
        volDrop = vol / n; % [m3] volume of a drop
        r = (6/4*volDrop/pi)^(1/3); % [m] new radius of the drop
        rs(Iter) = r;
        t = t + dt;
        time(Iter) = t;
        heatFlux(Iter) = qbn/dt;
        Iter = Iter + 1;
    else 
        qb = (deltaT - (2*sigma*Tsat)/(r*rhol*Hlv)) /... 
            ((C1*r/kl) + C2*corr*Tsat/(Hlv^2*rhov)*(gamma+1)/(gamma-1)*((Rg*Tsat)/(2*pi))^0.5)...
            * dt % [J] heat through the base of one droplet From Rose (1998)
        qbn = qb * n % [J] heat through all the drops
        volHold = 4/6*pi*n*rmax^3 % [m3] liquid on surface
        volAlpha = 4/6*pi*n*r^3 - volHold + qbn / (rhol*Hlv) % [m3]
        r = rmax % [m]
        alpha = volAlpha / (Lx*Ly*Lz) % []
        n = ceil(n*(1-alpha)) % []
        r = max((6/4*volHold/pi/n)^(1/3), rmin) % [m]
        rs(Iter) = r;
        t = t + dt
        time(Iter) = t;
        heatFlux(Iter) = qbn/dt;
        Iter = Iter + 1;
    end        

end
plot (time, heatFlux)

%% Second order implicit time stepping approach