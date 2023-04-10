% ---------------------------
%
% Script name: G2G_drying_scenario.R
%
% Purpose of script:  Script demonstrates differences in the behaviour of
%                     the Grid-to-Grid model and the physical benchmark
%                     model when applied to a drying scenario (with a
%                     positive initial mean rainfall r0>0, but no rainfall
%                     during the simulation, r=0).
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
%   - MODELS/G2G.m                   Grid-to-Grid Model class
%
%   - MODELS/scenario_B_settings.m   benchmark settings
%
%   - DATA/Physical_drying_scenario/drying_scenario.mat
%     benchmark scenario generated with 'Physical_drying_scenario.m' script
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path

load('scenario_B_settings.mat') % load predefined settings; alternatively
                                % you can define your own setting

dir = 'DATA/G2G_drying_scenario/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Demonstration using method of characteristic

% Here we show how does the characteristic diagram of the surface water
% height differ between Grid model with constant flow speed, and physical
% benchmark model, in which the flow speed  is described with a Manning's
% law, q=h^(5/3}. We will assume that initially h(x,t=0)=1-x, which
% correspond to a steady state with a mean precipitation rate r0=1.

k = 5/3;                  % exponent from the Manning's law
x0 = linspace(0,1,101);   % starting locations of characteristic curves
q0 = 1 - x0;              % initially overland flow is equal to 1-x0
c = 0.5;                  % speed used in Grid model

points = 6:10:101;        % few points will be highlighted to demonstrate
                          % the effect of varying speed
                          
t_values = [0, 0.5, 1, 1.5, 2];   % h(x,t) solution will be plotted
                                  % at these times

nt = length(t_values);    % number of time steps

x_data = cell(nt,2);      % array to save x(t) data for each model

% Calculate location of a given characteristic curve assuming:
for i = 1:nt
  
  % 1) that its speed is constant equal to c
  x_data{i,1}  = x0 - c * t_values(i);
  
  % 2) that its speed depends on the value of q0 as in the Manning's law
  x_data{i,2}  = x0 - q0 .^ (1 - 1 / k) * t_values(i);
end

% Plot the resulting h(x,t) profiles for both models
hf = figure(1);
for i = 1:2
  subplot(1,2,i);
  hold on
  % For each timestep plot profiles h(x,t) corresponding to each time step
  for j = 1:nt
    x = x_data{j,i};
    plot(x, q0)
  end
  
  % Highlight points specified by 'points array'
  for j = 1:nt
    x = x_data{j,i};
    plot(x(points), q0(points), '.k', 'MarkerSize', 15)
    
    % For all profiles except for the last one an arrow will be added to
    % represent the speed with which a given characteristic curve
    % propagates
    if j < nt
      
      % Find the speed according to one of the models
      if i == 1
        v = c * ones(size(q0));   % constant speed
      else
        v = q0 .^ (1 - 1 / k);    % speed given by the Manning's law
      end
      v = 0.5 * v;    % the velocity is multiplied by dt=0.5
      
      % Add arrows to the plot
      for p = points
        if v(p) > 0
          h = annotation('arrow', 'Color', 0.8*[1,1,1]);
          h.Parent = hf.CurrentAxes;
          h.X = [x(p) - 0.02, x(p) - v(p) + 0.02];
          h.Y = [q0(p), q0(p)];
        end
      end
    end
    
    % On the second plot add a legend
    if i == 2
      legend('t = 0', 't = 0.5', 't = 1', 't = 1.5', 't = 2')
    end
  end
  hold off
  
  % Set custom axis limits, labels and title
  xlim([0,1])
  ylim([0,1])
  xlabel('x')
  ylabel('q_s(x,t)')
  if i == 1
    title('Grid model (u=const.)')
  else
    title('Physical benchmark (u~q^{2/5})')
  end
end

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [100, 100, 600, 250])
set(gcf,'color','w');

% Export figure to a .png file
exportgraphics(gcf, 'FIGURES/G2G_drying_scenario_demo.png', ...
  'Resolution', 300)


%% Grid-to-Grid model calibration to a physical benchmark hydrograph

% The rest of the script compares the hydrograph obtained using a physical
% benchmark model and compares it with the fitted Grid-to-Grid model
% hydrograph.

% For training and plotting we use scenario saved in the following file,
% which is generated with 'Physical_dyring_scenario.m' script

load('DATA/Physical_drying_scenario/drying_scenario.mat', ...
  't_data', 'q_data', 'settings')

% Calculate the total precipitation over the hillslope (it is required by
% the current Grid-to-Grid model implementation

settings.P0 = settings.r0 * settings.Lx;

% Calibration takes some time, however a calibrated model is saved in a
% file. If you want to recalibrate the model (e.g. after changing the
% scenario settings), set 'recalibrate' to true.

recalibrate = false;

if (recalibrate)

  % Initialise a GridToGrid model object
  model = GridToGrid();
  
  % Now we have to set an initial guess for the Grid-to-Grid parameters.
  % Providing a good initial guess minimises the risk that fitting
  % algorithm will find a suboptimal fit.

  %       c     cb   r c_max  b   k  beta   nx
  x0 = [1e-1, 1e-10, 0, 1e-4, 1, 0.1,   1, settings.nx];
  
  % Set the guessed parameters
  model = model.setParameters(x0);

  % Calibrate the model to match the physical benchmark hydrograph
  [model, summary] = model.calibrate({settings.P0}, {t_data}, {0}, {q_data});
  
  % Save calibrated model and its summary
  save([dir, 'calibrated_model.mat'], 'model', 'summary')
end
  
% Load the calibrated Grid-to-Grid model
load([dir, 'calibrated_model.mat'], 'model', 'summary')

% After the fitting is complete, find the steady state of the fitted model
% for precipitation rate P0, and use it as the initial condition.
model = model.setInitialCondition('steady state', settings.P0);

% Run Grid-to-Grid simulation of the drying scenario (with P=0)
[solution, hydrograph] = model.simulate(0, t_data(end), settings.nt);

%% Compare simulation results

% Here we will plot both physical benchmark model and Grid-to-Grid model
% highlighting their early- and late-time behaviour

hf = figure(2);

% On the first graph show physical benchmark hydrograph
subplot(1,2,1)
plot(t_data, q_data, 'LineWidth', 2)

% Add to it early- and late-time asymptotics (see derivation in our paper;
% you can find the reference in the README file)
hold on
settings.a0 = 1 - 1 / settings.rho0;    % size of the saturated zone

% we have to convert dimensional time to dimensionless time used to derive
% the asymptotic behaviour
t = t_data * 3600 / settings.T0 * settings.mu ^ (1 / settings.k); 

% Plot the early-time asymptotics (linear)
plot(t_data, settings.Q_scale + settings.rho0 * settings.a0 * ...
  (1 - settings.k * settings.rho0 / (settings.rho0 *settings.a0) ^ ...
  (1 / settings.k) * t) * settings.Q_scale, '--c', 'LineWidth', 1.5);

% Plot the late-time asymptotics (power law)
plot(t_data, settings.Q_scale + (settings.a0 ./ (settings.k * t)) .^ ...
  (settings.k / (settings.k - 1)) * settings.Q_scale, '--g', 'LineWidth', 1.5);
hold off

% Add custom axes limits, labels, graph's title and legend
ylim([5e-6, 2e-5]);
xlabel('time, t [h]')
ylabel('total river inflow, Q [m^2/s]')
title('Physical benchmark model')
legend('numerical solution', 'early-time asymptotics (linear)', ...
  'late-time asymptotics (q~t^{-5/2})')


% On the second graph show physical benchmark hydrograph
subplot(1,2,2);
plot(hydrograph.t, hydrograph.Q, 'LineWidth', 2)

% Add to it early- and late-time asymptotics
hold on
base_flow = hydrograph.Q(end);    % baseflow is estimated based on
                                  % the final value of the flow

% Add early-time asymptotics (linear)
plot(hydrograph.t, base_flow + (hydrograph.Q(1)-base_flow) * ...
  (1 - hydrograph.t / (1 / summary.c)), '--c', 'LineWidth', 1)

% Add line representing the base flow
plot(hydrograph.t, base_flow * ones(size(hydrograph.t)), ...
  '--g', 'LineWidth', 1)
hold off

% Add custom axes limits, labels, graph's title and legend
ylim([5e-6, 2e-5]);
xlabel('time, t [h]')
ylabel('total river inflow, Q [m^2/s]')
legend('fitted numerical solution', 'analytic approximation (linear)', ...
  'base flow')
title('Fitted Grid-to-Grid model')

% Add a label showing the numerical diffusion during the transition from
% linear to constant regime
text(5, 1.1e-5, 'effect of the numerical diffusion')
h = annotation('line');
h.Parent = hf.CurrentAxes;
h.X = [5, 7];
h.Y = [7.9e-6, 1e-5];

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [100, 100, 800, 300])
set(gcf,'color','w');

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/G2G_dry_scenario.png', 'Resolution', 300)