function [sim,pred,prey] = params_danio
% Parameter values for danio predator-prey.  Assumes axial symmetry for the
% prey.

%% PREY PARAMETERS 

% Load larval body shape ('bod' structure)
load('larva_body_shape.mat')

% Prey length (m)
prey.len = 3.5e-3;

% Prey density (kg m^-3)
prey.rho = 1000;

% Number of body segments defining morphology
prey.num_segs = 10;

% Position of COM
prey.sCOM = 0.25 * prey.len;

% Body positons of segments
prey.s = linspace(0,prey.len,prey.num_segs)';

% Radius of body (m)
prey.r = interp1(bod.s,bod.r,prey.s);

%width = flipud(width);
% x-Sectional area of segments
prey.x_area = pi.*prey.r.^2;     

% Section volume
prey.vol = [0; diff(prey.s)].*pi.*prey.r.^2;

% Wetted area of the segments
prey.wet_area = [0; diff(prey.s)].*2*pi.*prey.r;

% Volume of body
prey.bod_vol =  trapz(prey.s,prey.x_area);

% Mass
prey.mass = prey.bod_vol * prey.rho;

% Added mass coefficient for a cylinder(?) (dimensionless)
prey.add_mass = 0.6; 

% Thrust (N)
prey.th_max    = 0*-2e-6;

% Timing of peak thrust (after trigger, s)
prey.th_delay  = 15e-3;

% Drag coefficient for body segment, set to cylinder drag coefficent for 
% above-critical Re (Hoerner, 1965)
prey.Cd = 1.2;

% Initial prey position (x (m), y (m), orientaton angle (rad))
prey.pos0= [5e-3 0 pi];

% Initial prey speed (x(m/s), y(m/s), angle rate (rad/s))
prey.vel0 = [0 0 0];

% Sensitivity threshold (distance from predator, m)
prey.thresh = 10e-3;


%% PREDATOR PARAMETERS 

% Max suction speed (m/s)
pred.flw_spd.max = 1;

% Time of max speed (s)
pred.flw_spd.t_max = 40e-3;

% Shape factor for speed
pred.flw_spd.alpha = 2;

% Max gape (m)
pred.gape.max = 0.22e-3;

% Time of max gape (s)
pred.gape.t_max = 30e-3;

% Shape factor
pred.gape.alpha = 2;

% Initial position of predator (m)
pred.pos0 = [0 0 0];

% Approach speed (m/s)
pred.app_spd = 10e-2;


%% SIMULATION PARAMETERS

% Water density (kg m^-3)
sim.rho_water = 1000;

% Kinematic viscosity (Pa s)
sim.mu_water = 0.001;

% Maximum simulation duration (s)
sim.dur = .15; 

% Relative tolerance for the solver (dimensionless)
sim.reltol = 1e-9;

% Number of time values to define predator flow
sim.num_time = 200;

% Number of values along the x-axis to define pred flow
sim.num_x = 200;

% Boundaries in x and y directions for defining pred flow
sim.flow_lim = [-3e-3 15e-3 -7e-3 7e-3];

% Scaling constants
sim.sL = prey.len;
sim.sT = 10^-3;
sim.sM = prey.mass*10^6;
sim.sF = sim.sM .* sim.sL ./ sim.sT^2;

end




