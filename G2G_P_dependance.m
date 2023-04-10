% ---------------------------
%
% Script name: G2G_P_dependance.R
%
% Purpose of script:  Plots Grid-to-Grid simulation results calibrated to
%                     the physical benchmark model results for different
%                     values of the simulated precipitation rate, P.
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
%  - MODELS/GridToGrid.m     Grid-to-Grid model class
%
%  - 'DATA/Physical_benchmarks/hydrograph_P_*.mat'
%    benchmark hydrographs generated with a 'Physical_benchmarks.m' script
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
input_dir = 'DATA/Physical_benchmarks/';   % input directory with the
                                              % physical benchmark results
                                
dir = 'DATA/G2G_P_dependance/'; % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Set training and validation set

% In this script a Grid-to-Grid model (G2G) will be fitted to
% results of the physical benchmark model. For this purpose two data sets
% are created:
%   - training data set, which is used to calibrate G2G parameters, and
%   - validation data set, used to assess the accuracy of the fitted G2G.

% The precipitation rate (P) in each simulation will be a different multple
% of mean precipitation rate (P0).

% We generate training data for P=2*P0, P=4*P0, P=6*P0 and P=8*P0
training_data = [2, 4, 6, 8];

% For validation data we use P=4*P0, P=8*P0, P=12*P0 and P=16*P0
validation_data = [4, 8, 12, 16];

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
  filename = [input_dir, 'hydrograph_P_', ...
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
  filename = [input_dir, 'hydrograph_P_', ...
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

%     c    cb    r    c_max  b  k beta  nx
x0 = [0.5, 1e-4, 0.03, 0.15, 3, 1,   3, 20];

% The quality of this initial guess can be assessed using
% compare_hydrographs() function.
figure(1)
subplot(1,2,1)
compare_hydrographs(x0, P0_test_data, P_test_data, t_test_data, ...
  Q_test_data, 24);

subplot(1,2,2)
compare_hydrographs(x0, P0_test_data, P_test_data, t_test_data, ...
  Q_test_data, 24);
xlim([0,24])

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

% Add number of mesh elements to the summary

summary.nx = x0(end);

% Compare the validation data for the whole simulated week with the first
% six hours zoomed

figure(2)
t_zoom = 6;
compare_hydrographs(summary, P0_valid_data, P_valid_data, ...
  t_valid_data, Q_valid_data, t_zoom, dir);
exportgraphics(gcf, 'FIGURES/G2G_P_validation_set.png', 'Resolution', 300)


%% Functions

% Function compare_hydrographs() generates a figure, on which the solutions
% of the G2G simulation with the physical benchmark model for the provided
% scenarios.
%
%   INPUT:
%     summary   - structure containing G2G model parameters
%                 (c_max, b, kg, st, kf, ks)
%     P0_data   - mean precipitation rate
%     P_data    - simulated precipitation rate
%     t_data    - time steps of the simulation
%     Q_data    - measured total flow in each time step
%     t_zoom    - (optional) if specified the plotted data will be
%                 zoomed for the first t_zoom (e.g. 6) hours
%     dir       - (optional) output directory

function compare_hydrographs(summary, P0_data, P_data, t_data, Q_data, ...
                             t_zoom, dir)

  % Initialize a GridToGrid class object
  model = GridToGrid();
  
  % Set its parameters to the ones provided in the 'summary' structure
  model = model.setParameters(summary);
  
  % Check number of validation set elements
  N = length(P0_data);
  
  % For each item in the validation set a separate subplot is generated
  for i = 1:N
    
    % Obtain parameters specifying a given scenario
    P0 = P0_data{i};
    P = P_data{i};
    t = t_data{i};
    Q = Q_data{i};
    
    % Set initial conditions of the G2G model to the steady state for
    % the given mean precipitation rate P0
    model = model.setInitialCondition('steady state', P0);
    
    % Run G2G simulation for the specified time (with 10k time steps)
    [~, hydrograph] = model.simulate(P, t(end), 10000);
    
    % Normalise the flow in the physical benchmark hydrograph and G2G
    % hydrograph using formula (Q - P0) / (P - P0); in this parametrisation
    % the flow grows from 0 at t=0 to 1 for t->ininitiy
    benchmark_data = [t', (Q - P0) / (P - P0)];
    g2g_data = [hydrograph.t', (hydrograph.Q_total' - P0) / (P - P0)];
    
    % Export these hydrographs to .dat files if a directory is provided
    if nargin == 7
      writematrix(benchmark_data(1:10:end,:), ...
        [dir, 'benchmark_validation_', num2str(i), '.dat']);
      writematrix(g2g_data(1:100:end,:), ...
        [dir, 'G2G_validation_', num2str(i), '.dat'])
      writematrix(benchmark_data(benchmark_data(:,1)<=t_zoom,:), ...
        [dir, 'benchmark_validation_', num2str(i), '_zoomed.dat']);
      writematrix(g2g_data(g2g_data(:,1)<=t_zoom,:), ...
        [dir, 'G2G_validation_', num2str(i), '_zoomed.dat'])
    end
    
    % Plot graph in a new subplot
    h = subplot(2, N/2, i);
    
    % If a t_zoom option is used highlight zoomed area in yellow
    if nargin > 5
      patch([0 t_zoom t_zoom 0], [0, 0, 1, 1], 'y', 'FaceAlpha', 0.2, ...
        'edgecolor','none')
    end
    
    % Plot dimensionless hydrograph for physical benchmark model and G2G
    % solution
    hold on
    plot(t, (Q - P0) / (P - P0))
    plot(hydrograph.t, (hydrograph.Q_total - P0) / (P - P0))
    hold off
    
    % Set custom axes range, labels, ticks and subfigure's title
    xlim([0, max(t)])
    ylim([0,1])
    xlabel('t [h]')
    ylabel(strcat('\Delta','Q / (P - P_0)'))
    xticks([0,24,72,120,168])
    title(strcat("P = ", num2str(P/1.8172e-5), " P_0"))
    
    % If a t_zoom option is used add a smaller graph with a zoomed version
    % of the hydrograph
    if nargin > 5
      
      % Set the position of the axes for this smaller graph
      pos = get(h, 'Position');
      pos = pos + [0.1, 0.05, -0.15, -0.15];
      h = axes('Position', pos);
      
      % The plot is generated in exactly the same way as the main one
      patch([0 t_zoom t_zoom 0], [0, 0, 1, 1], 'y', 'FaceAlpha', 0.2, ...
        'edgecolor','none')
      hold on
      plot(t, (Q - P0) / (P - P0))
      plot(hydrograph.t, (hydrograph.Q_total - P0) / (P - P0))
      hold off
      
      % Except for x limits, which are set to [0, t_zoom]
      xlim([0, t_zoom])
    end
  end
  
  % Set the size of the whole figure and its background color to white
  set(gcf, 'Position',  [100, 100, 800, 600])
  set(gcf,'color','w');
end
