% ---------------------------
%
% Script name: Physical_drying_scenario.R
%
% Purpose of script:  Generate physical benchmark hydrograph for the
%                     drying scenario, in which a catchment is initially
%                     in a steady state with a rainfall P0, but at t=0
%                     the precipitation rate drops to 0.
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
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path

dir = 'DATA/Physical_drying_scenario/';   % output directory

if ~exist(dir, 'dir')           % create the directory if they do not exist
  mkdir(dir)
end

% Load settings, for which the physical benchmark simulation is performed
% (here we use settings representing a low productive aquifer from
% https://github.com/Piotr-Morawiecki/benchmark-catchment-model repo)
load('scenario_B_settings.mat', 'settings')

% Change precipitation rate during the simulation to 0
settings.r = 0;

% Change simulated period to entire week
settings.t = 7 * 24 * 3600;

% Set drainable porosity to a constant (unfortunately the physical model
% formulation with the varying f(x) function cannot be applied in the
% drying scenario)

settings.f = 0.1;

% Set computational mesh size and number of time steps
settings.nx = 200;
settings.nt = 7001;

%% Run a simulation

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

%% Export hydrograph for the first 24 hours

% Calculate the total river inflow (Q), as well as its surface (Qf) and
% groundwater (Qs) flow components
q_data = catchment_1D.computeFlow(xmesh, sol);

% Multiplying them by Q_scale allows us to obtain flow in dimensional
% quantities [m^2/s]
q_data  = q_data  * settings.Q_scale;

% Convert dimensionless time to [s]
t_data = tspan * settings.T0 / 3600;

% Extract flow and time steps for the first 24 hours (only they are used
% in 'G2G_drying_scenario.m' script
q_data = q_data(t_data <= 24);
t_data = t_data(t_data <= 24);

% Export data to a .mat file
save([dir, 'drying_scenario.mat'], 't_data', 'q_data', 'settings')

%% Plot the solution

% Here two graphs are plotted - one with few surface water profiles and
% second with with the corresponding groundwater profiles

% Plot surface water height profiles hs(x)
subplot(1,2,1)
plot(xmesh, sol(1:25:151,:))

ylim([0,Inf])   % Restricting y-axis range to positive values allows us
                % not to plot groundwater table shape (it corresponds to
                % negative values in the sol matrix)

% Set custom axes labels, title and legend
xlabel('x')
ylabel('h_s(x,t)')
title('(a) surface water height')
legend('t=0 h', 't=1 h', 't=2 h', 't=3 h', 't=4 h', 't=5 h', 't=6 h')

% Plot surface water height profiles hs(x)
subplot(1,2,2)
plot(xmesh, 1+sol(1:1000:7001,:))   % we add 1, since it corresponds to
                                    % the maximal groundwater depth
                        
% Set custom axes range, labels, title and legend   
ylim([0,Inf])
xlabel('x')
ylabel('H(x,t)')
title('(b) groundwater height')
legend('t=0 d', 't=1 d', 't=2 d', 't=3 d', 't=4 d', 't=5 d', 't=6 d', ...
  't=7 d', 'Location', 'SouthWest')

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [50, 100, 850, 300])
set(gcf,'color','w');

% Export figure to a .png file
exportgraphics(gcf, 'FIGURES/G2G_drying_solution.png', 'Resolution', 300)
