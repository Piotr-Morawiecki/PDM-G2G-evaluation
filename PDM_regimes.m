% ---------------------------
%
% Script name: PDM_regimes.R
%
% Purpose of script:  Plots hydrograph representing three different regimes
%                     of the probability-distributed model (PDM):
%                       (1) P0 < P < d_max,
%                       (2) P0 < d_max < P,
%                       (3) d_max < P0 < P
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
%  - MODELS/PDM.m     Probability Distributed Model class
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
dir = 'DATA/PDM_regimes/';      % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Prepare simulation settings

% Here we create a par structure, which will include all parameters
% required by PDM class and further postprocessing of simulation results.

% Choose PDM model parameters
par.c_max = 1;
par.b = 3;
par.kg = 1;
par.st = 0;
par.kf = 0.1;
par.ks = 1;

par.t_max = 100;  % simulated time
par.nt = 10000;   % number of timesteps

% Check the maximum value of drainage rate

par.d_max = (par.c_max / (par.b + 1) - par.st) / par.kg;
fprintf('d_max=%f\n', par.d_max)

% Set individual values of mean and peak precipitation for each of the
% three simulated regimes (for reference d_max = 0.25

par.P0 = [0.1, 0.1, 0.3]; % mean precipitation (for the initial condition)
par.P = [0.2, 0.4, 0.4];  % peak precipitation (for the simulated rainfall)

% Check number of simulations (default: 3)

n = length(par.P);

%% Run PDM simulation

% The simulation results will be stored in 'hydrographs' cell object
hydrographs = cell(1, n);

% For each simulation:
for i = 1:n
  
  % 1) Initialize PDM class object
  model = PDM();

  % 2) Set its parameters to the ones defined in par structure
  model = model.setParameters(par);

  % 3) Set the initial condition to be a steady state corresponding to the
  %    mean rainfall P0
  model = model.setInitialCondition('steady state', par.P0(i));

  % 4) Run PDM simulation for rainfall P of duration t_max with nt time
  %    steps and save the resulting hydrograph in the 'hydrographs' object
  [~, hydrographs{i}] = model.simulate(par.P(i), par.t_max, par.nt);
end

%% Plotting and exporting simulation results

for i = 1:n   % for each simulation
  
  % Plot a hydrograph and save plotted data to a separate 'data' object
  subplot(1,n,i)
  data = plot_hydrograph(hydrographs{i}, par.P(i), par.P0(i), par.d_max);
  
  % Set graph title and y-axis limit
  title(strcat("Regime ", num2str(i), " (P_0=", num2str(par.P0(i)), ...
    ", P=", num2str(par.P(i)), ")"))
  ylim([0,0.55])
  
  % In case of the first plot add a legend
  if i == 1
    legend('overland flow, Q_f(t)', 'groundwater flow, Q_s(t)', ...
      'total flow, Q(t)', 'precipitation rate, P', ...
      'initial precipitation, P_0', 'max drainage, d_{max}', ...
      'Location', 'NorthWest')
  end

  % The simulation results will be exported with a lower time resolution
  % (to reduce output data size). We save the values of flow components at
  % 100 time points, equally distributed in log space between t=1e-2 and
  % t=log10(par.t_max).
  t = logspace(-2, log10(par.t_max), 100);
  
  % The value of each flow component is found through the interpolation
  data = [t, interp1(data(:,1), data(:,2), t), ...
             interp1(data(:,1), data(:,3), t), ...
             interp1(data(:,1), data(:,4), t)];
           
  % The resulting data set is exported to a .dat file
  writematrix(data, [dir, 'hydrograph_regime_', num2str(i), '.dat']);
end

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [50, 100, 1000, 350])
set(gcf,'color','w');

exportgraphics(gcf, 'FIGURES/PDM_regimes.png', 'Resolution', 300)

%% Functions

% Function plot_hydrograph() plots a hydrograph for specified data
%
% INPUT:
%         hydrograph - hydrograph produced by the PDM.simulate() method
%         P          - simulated precipitation rate
%         P0         - mean precipitation rate (from the initial condition)
%         d_max      - maximum drainage rate
%
% OUTPUT:
%         data - plotted data in an array format

function [data] = plot_hydrograph(hydrograph, P, P0, d_max)

  % Plot a hydrograph including fast, slow and total flow components
  plot(hydrograph.t, hydrograph.Qf, 'LineWidth', 2)
  hold on
  plot(hydrograph.t, hydrograph.Qs, 'LineWidth', 2)
  plot(hydrograph.t, hydrograph.Q_total, 'LineWidth', 2)
  
  % Add lines representing 
  yline(P, '--k', 'LineWidth', 1.5)
  yline(P0, ':k', 'LineWidth', 1.5)
  yline(d_max, '-.k', 'LineWidth', 1.5)
  hold off
  
  % Add axes labels and set log scale on x-axis
  xlabel('time')
  ylabel('flow')
  set(gca, 'XScale', 'log')
  
  % Put all plotted data in an array
  data = [hydrograph.t',hydrograph.Qf',hydrograph.Qs',hydrograph.Q_total'];
end