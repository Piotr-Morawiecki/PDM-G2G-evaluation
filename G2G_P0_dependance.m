% ---------------------------
%
% Script name: G2G_P0_dependance.R
%
% Purpose of script:  Plots Grid-to-Grid simulation results calibrated to
%                     the physical benchmark model results for different
%                     values of the mean precipitation rate, P0, used to
%                     specify the initial condition.
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
%  - MODELS/GridToGrid.m     Probability Distributed Model class
%
%  - 'DATA/Physical_benchmarks/hydrograph_a0_*.mat'
%    benchmark hydrographs generated with a 'Physical_benchmarks.m' script
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
input_dir = 'DATA/Physical_benchmarks/';   % input directory with the
                                              % physical benchmark results

dir = 'DATA/G2G_P0_dependance/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Set training and validation set

% In this script a Grid-to-Grid model (G2G) will be fitted to
% results of the physical benchmark model. For this purpose two data sets
% are created:
%   - training data set, which is used to calibrate G2G parameters, and
%   - validation data set, used to assess the accuracy of the fitted G2G.

% The precipitation rate (P0) in each simulation correspond to a different
% value of the intially saturated zone size (a0), which according to
% the benchmark physical model can be expressed as
%
%                     a0 = 1 - (K_s*L_z*S_x) / (P0*L_x).

% We generate training data for a0=0.35 and a0=0.65
training_data = [0.35, 0.65];

% For validation data we use a0 equal to 0.2, 0.35, 0.5, 0.65 and 0.8
validation_data = [0.2, 0.35, 0.5, 0.65, 0.8];

% IMPORTANT:  If you change these values, make sure to generate relevant
%             physical benchmark results using 'Physical_benchmarks.m'
%             script.


%% Load precomputed benchmark scenarios

% The training data will be saved in four cell arrays:

P0_test_data = cell(1, length(training_data));  % for mean precipitation
P_test_data = cell(1, length(training_data));   % for peak precipitation
t_test_data = cell(1, length(training_data));   % for time values, t
Q_test_data = cell(1, length(training_data));   % for total flow, Q(t)

% For each element of the training data set open the corresponding
% physical benchmark hydrograph

for i = 1:length(training_data)
  filename = [input_dir, 'hydrograph_a0_', ...
              num2str(training_data(i)), '.mat'];
  try
    file_content = load(filename);
  catch
    error(['File ', filename, ' does not exist. ', ...
      'Generate it using Physical_benchmarks.m script.']);
  end
  P0_test_data{i} = file_content.P0;
  P_test_data{i} = file_content.P;
  t_test_data{i} = file_content.t;
  Q_test_data{i} = file_content.Q;
end

% The training data will be saved in six cell arrays:

% - the same four cell arrays as in case of the training data set,
P0_valid_data = cell(1, length(training_data));
P_valid_data = cell(1, length(training_data));
t_valid_data = cell(1, length(training_data));
Q_valid_data = cell(1, length(training_data));

% - and two cell arrays to store information about the flow components.
Q_subsurf_valid_data = cell(1, length(training_data));
Q_surf_valid_data = cell(1, length(training_data));

% For each element of the validation data set open the corresponding
% physical benchmark hydrograph

for i = 1:length(validation_data)
  filename = [input_dir, 'hydrograph_a0_', ...
              num2str(validation_data(i)), '.mat'];
  try
    file_content = load(filename);
  catch
    error(['File ', filename, ' does not exist. ', ...
      'Generate it using Physical_benchmarks.m script.']);
  end
  P0_valid_data{i} = file_content.P0;
  P_valid_data{i} = file_content.P;
  t_valid_data{i} = file_content.t;
  Q_valid_data{i} = file_content.Q;
  Q_subsurf_valid_data{i} = file_content.Qs;
  Q_surf_valid_data{i} = file_content.Qf;
end

%% Guess initial parameters

% We start calibration process by setting (arbitrary) starting values of
% Grid-to-Grid parameters (here we set really good starting values to
% reduce the duration of the calibration process)

%       c    cb     r c_max  b  k beta  nx
x0 = [0.5, 1e-4, 0.03, 0.15, 3, 1,   3, 20];

% The quality of this initial guess can be assessed using
% compare_hydrographs() function.
figure(1)
subplot(1,2,1)
compare_hydrographs(x0, P0_test_data, P_test_data, t_test_data, ...
  Q_test_data, validation_data);

xlim([0,168])       % plot it for entire week (168 hours)
xticks(0:24:168)

subplot(1,2,2)
compare_hydrographs(x0, P0_test_data, P_test_data, t_test_data, ...
  Q_test_data, validation_data);

xlim([0,24])        % plot it for the first day (24 hours)
xticks(0:6:24)

drawnow;

% WARNING: If the quality of the model is poor, the calibration method
% may fail to converge or give an suboptimal solution.


%% Calibrate the model

% Calibration may take a significant amount of computing time.

% The results are however stored in a 'calibration_summary.mat' file, and
% so this section of the script can be skipped. You can do this by setting
% 'recalibrate' to false

recalibrate = false;

if recalibrate
  
  % Initialize a G2G class object
  model = GridToGrid();
  
  % Set the initial Grid-to-Grid model parameters
  model = model.setParameters(x0);

  % Fit G2G model parameters to the available training data
  [~, summary] = model.calibrate(P0_test_data, ...
    t_test_data, P_test_data, Q_test_data);

  % Save calibration results to a .mat file.
  save([dir, 'calibration_summary.mat'], 'summary');

else
  % Import previous calibration results
  load([dir, 'calibration_summary.mat'], 'summary');
end

% Display fitted parameters
disp(summary)


%% Plot validation data

summary.nx = 20;

% Compare the validation data for the whole simulated week with the first
% six hours zoomed
figure(2)
[benchmark_data, g2g_data] = compare_hydrographs(summary, ...
  P0_valid_data, P_valid_data, t_valid_data, Q_valid_data, ...
  validation_data);

% plot it for entire week (168 hours)
xlim([0,168])
xticks(0:24:168)

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/G2G_P0_Validation_set.png', 'Resolution', 300)

% Export the plotted data to .dat files (every 10th row to reduce the
% output data size)
writematrix(g2g_data(1:10:end,:), [dir, 'g2g_r0_dependence.dat']);
writematrix(benchmark_data(1:10:end,:), ...
  [dir, 'benchmark_r0_dependence.dat']);


%% Functions

% Function compare_hydrographs() generates a figure, on which the solutions
% of the G2G simulation with the physical benchmark model for the provided
% scenarios.
%
% INPUT:
%   summary   - structure containing G2G model parameters
%               (c_max, b, kg, st, kf, ks)
%   P0_data   - mean precipitation rate
%   P_data    - simulated precipitation rate
%   t_data    - time steps of the simulation
%   Q_data    - measured total flow in each time step
%   a0_values - values of a0 parameter to include in the legend
%
% OUTPUT:
%   benchmark_data  - data plotted for the physical benchmark model
%   g2g_data        - data plotted for the G2G

function [benchmark_data, g2g_data] = compare_hydrographs(summary, ...
  P0_data, P_data, t_data, Q_data, a0_values)

  % Initialize a GridToGrid class object
  model = GridToGrid();
  
  % Set its parameters to the ones provided in the 'summary' structure
  model = model.setParameters(summary);

  % Set colors in which results of each scenario will be represented
  colors = {'r', 'g', 'b', 'c', ' m'};
  
  % Initialize a cell array to store legend labels
  legend_entries = cell(1, length(P0_data));
  
  % Create an empty plot to which results for each validation scenario will
  % be added
  hold on
  
  % For each scenario, G2G and physical benchmark model predictions
  % will be added
  for i = 1:length(P0_data)
    
    % Get mean precipitation (P0), time (t) and flow (Q) data as predicted
    % by the physical benchmark model
    P0 = P0_data{i};
    t = t_data{i};
    Q = Q_data{i};
    
    % Set initial condition to be a steady state for given mean
    % precipitation P0
    model = model.setInitialCondition('steady state', P0);
    
    % Run a time-dependent G2G simulation for given precipitation P,
    % of the same length as in case of physical benchmark, t(end) with
    % 10000 time steps
    [~, hydrograph] = model.simulate(P_data{i}, t(end), 10000);

    % If this is the first hydrograph initialise two data arrays to store
    % data plotted for Grid-to-Grid and physical benchmark model; otherwise
    % add this data to the already existing arrays.
    if i == 1
      benchmark_data = [t' / 3600, Q];
      g2g_data = [hydrograph.t', hydrograph.Q_total'];
    else
      benchmark_data = [benchmark_data, Q];
      g2g_data = [g2g_data, hydrograph.Q_total'];
    end

    % Plot hydrograph for physical benchmark model (solid line)
    plot(t, Q, colors{i})
    
    % and G2G simulation (dashed line)
    plot(hydrograph.t, hydrograph.Q_total, ['--', colors{i}], ...
      'HandleVisibility', 'off');
    
    % Add corresponding legend entry to be displayed later
    legend_entries{i} = ['a_0=', num2str(a0_values(i),2)];
  end
  hold off
  
  % Set custom axes labels
  xlabel('Time, t [h]')
  ylabel('Total flow, Q')
  
  % Add legend
  legend(legend_entries, 'Location', 'southeast')

  % Set the size of the whole figure and its background color to white
  set(gcf, 'Position',  [100, 100, 600, 350])
  set(gcf, 'color', 'w');
end