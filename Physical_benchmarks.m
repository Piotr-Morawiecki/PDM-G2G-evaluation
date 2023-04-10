% ---------------------------
%
% Script name: Physical_benchmarks.R
%
% Purpose of script:  Generate physical benchmark hydrographs for different
%                     values of peak and mean precipitation rate, P and P0,
%                     which are used in other scripts to assess PDM and
%                     Grid-to-Grid model.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-04-02
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------
%
% Requires:
%
%   - MODELS/Catchment1D.m          1D physical benchmark model class
%
%   - MODELS/scenario_B_settings.m  benchmark settings
%
% WARNING:  This script takes a significant amount of time to compute.
%           Use the precomputed results when available.
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path

dir = 'DATA/Physical_benchmarks/';   % output directory

if ~exist(dir, 'dir')           % create the directory if they do not exist
  mkdir(dir)
end


%% Set values of the saturation zone size (a0) and precipitation (P)

% In case of the simulations with varying mean precipitation rate (P0), its
% value correspond to different sizes of the intially saturated zone (a0),
% which according to the benchmark physical model can be expressed as
%
%                     a0 = 1 - (K_s*L_z*S_x) / (r0*L_x).

% We generate training data for the following a0 values:
a0_values = [0.2, 0.35, 0.5, 0.65, 0.8];

% In case of the simulations with varying simulated precipitation rate (P),
% its value in each simulation will be a different multiples of the mean
% precipitation rate (P0).

% We generate data for P=2*P0, P=4*P0, P=6*P0, P=8*P0, P=12*P0 and P=16*P0
P_multpliers = [2, 4, 6, 8, 12, 16];


%% Run physical benchmark simulation

% Load settings, for which the physical benchmark simulation is performed
% (here we use settings representing a low productive aquifer from
% https://github.com/Piotr-Morawiecki/benchmark-catchment-model repo)
load('scenario_B_settings.mat', 'settings');

% set simulated time in seconds (here we simulate one week of rainfall)
settings.t = 7 * 24 * 3600;

% Set mesh resolution and number of time steps to use in the simulation
settings.nx = 200;
settings.nt = 2100;
  
% Run simulation for each value of the saturated zone size (a0)
for a0 = a0_values
  settings.r0 = settings.K * settings.Sx * settings.Lz ./ ...
    (settings.Lx * (1-a0));
  filename = [dir, 'hydrograph_a0_', num2str(a0), '.mat'];
  fprintf('Running simulation for a0=%.2f\n', a0);
  generate_hydrograph(settings, filename)
end
  
% Run simulation for each value of P0 multiplier (m)
for m = P_multpliers
  settings.r = settings.r0 * m;
  filename = [dir, 'hydrograph_P_', num2str(m), '.mat'];
  fprintf('Running simulation for P=%d*P0\n', m);
  generate_hydrograph(settings, filename)
end


%% Functions

% Function generate_hydrograph() runs a 1D model simulation of the physical
% benchmark model, and exports its results to a file
%
% INPUT:
%
%   settings      - structure with simulation settings required by
%                   catchment_1D.setParametersDimensional() class as well
%                   as simulation time (t), number of time steps (nt) and
%                   spatial mesh resolution (nx)
%
%   output_name   - name of the file to which the results of the simulation
%                   are exported

function generate_hydrograph(settings, output_name)
  
  % To run a scenario we create Catchment1D object
  catchment_1D = Catchment1D();

  % Then we set parameters; setParametersDimensional is used for predefined
  % physical parameters in SI units
  [catchment_1D, settings] = ...
    catchment_1D.setParametersDimensional(settings);

  % We set uniform time steps and uniform element size
  tspan = linspace(0, settings.t / settings.T0, settings.nt);
  xmesh = linspace(0, 1, settings.nx);

  % We find steady state of the system and use it as the initial condition
  [steady_state, ~, catchment_1D] = catchment_1D.findSteadyState(xmesh);

  % Simulation is run to find H(x,t) for x in xmesh and t in tspan
  sol = catchment_1D.solve(steady_state, xmesh, tspan, false);
  
  % Calculate the total river inflow (Q), as well as its surface (Qf) and
  % groundwater (Qs) flow components
  [Q, Qf, Qs] = catchment_1D.computeFlow(xmesh, sol);
  
  % Multiplying them by Q_scale allows us to obtain flow in dimensional
  % quantities [m^2/s]
  Q  = Q  * settings.Q_scale;
  Qs = Qs * settings.Q_scale;
  Qf = Qf * settings.Q_scale;
  
  % Convert dimensionless time to [s]
  t = tspan * settings.T0 / 3600;
  
  % Convert the precipitation rates to [m^2/s]
  P0 = settings.r0 * settings.Lx;
  P = settings.r * settings.Lx;

  % Save results to .mat file; you can skip running this function and use
  % precalculated results for plotting.
  save(output_name, 'P0', 'P', 't', 'Q', 'Qs', 'Qf')
end