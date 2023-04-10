% ---------------------------
%
% Script name: PDM_P_dependance.R
%
% Purpose of script:  Plots PDM simulation results calibrated to the
%                     physical benchmark model results for different
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
%  - MODELS/PDM.m     Probability Distributed Model class
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
                                
dir = 'DATA/PDM_P_dependance/'; % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Set training and validation set

% In this script a probability distributed model (PDM) will be fitted to
% results of the physical benchmark model. For this purpose two data sets
% are created:
%   - training data set, which is used to calibrate PDM parameters, and
%   - validation data set, used to assess the accuracy of the fitted PDM.

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

% We start calibration process by guessing relatively good starting PDM
% parameters

%     c_max   b   kg   st   kf    ks
x0 = [ 0.03,  1,  720,  0,  0.5,  0.001];

% The quality of this initial guess can be assessed using
% compare_hydrographs() function.
figure(1)
compare_hydrographs(x0, P0_valid_data, P_valid_data, t_valid_data, ...
  Q_valid_data)

% WARNING: If the quality of the model is poor, the calibration method
% may fail to converge or give an suboptimal solution.

%% Calibrate the model

% Calibration may take a significant amount of computing time.

% The results are however stored in a 'calibration_summary.mat' file, and
% so this section of the script can be skipped. You can do this by setting
% 'recalibrate' to false

recalibrate = false;

if recalibrate
  
  % Initialize a PDM class object
  model = PDM();

  % Fit PDM model parameters to the available training data
  [~, summary] = model.calibrate(P0_valid_data, ...
    t_valid_data, P_valid_data, Q_valid_data, x0);

  % Save calibration results to a .mat file.
  save([dir, 'calibration_summary.mat'], 'summary');

else
  % Import previous calibration results
  load([dir, 'calibration_summary.mat'], 'summary');
end

% Display fitted parameters
disp(summary)

%% Plot validation data

% Compare the validation data for the whole simulated week with the first
% six hours zoomed

figure(2)
t_zoom = 6;
compare_hydrographs(summary, P0_valid_data, P_valid_data, ...
  t_valid_data, Q_valid_data, t_zoom, dir);
exportgraphics(gcf, 'FIGURES/PDM_P_Validation_set.png', 'Resolution', 300)

%% Compare hydrograph components for P=8*P0 simulation

% Plot comparison of surface, subsurface and total flow components between
% PDM simulation results and reference physical benchmark
figure(3)
[phys_data, pdm_data] = compare_hydrograph_components(summary, ...
  P0_valid_data{2}, P_valid_data{2}, t_valid_data{2}, Q_valid_data{2}, ...
  Q_subsurf_valid_data{2}, Q_surf_valid_data{2});

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/PDM_Validation_components.png', ...
  'Resolution', 300)

% Export data used to generate the figure to a .dat file
writematrix(phys_data(1:20:end,:), [dir, 'benchmark_components.dat']);
writematrix(pdm_data(1:100:end,:), [dir, 'pdm_components.dat']);

%% Calculate timescales

% Calculate the max soil moisture storage (S_max) and
% max drainage rate (d_max) for the fitted PDM
summary.S_max = summary.c_max / (summary.b + 1);
summary.d_max = (summary.S_max - summary.st) / summary.kg;

% Calculate the characteristic timescales for soil moisture store (tc),
% fast store (tf) and flow store (ts) for the fitted PDM
timescales.tc = summary.c_max ./ cell2mat(P_valid_data);
timescales.tf = summary.kf;
timescales.ts = summary.ks / summary.d_max;

% Display the timescales (in seconds)
disp(timescales)


%% Functions

% Function compare_hydrographs() generates a figure, on which the solutions
% of the PDM simulation with the physical benchmark model for the provided
% scenarios.
%
%   INPUT:
%     summary   - structure containing PDM model parameters
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

  % Initialize a PDM class object
  model = PDM();
  
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
    
    % Set initial conditions of the PDM model to the steady state for
    % the given mean precipitation rate P0
    model = model.setInitialCondition('steady state', P0);
    
    % Run PDM simulation for the specified time (with 10k time steps)
    [~, hydrograph] = model.simulate(P, t(end), 10000);
    
    % Normalise the flow in the physical benchmark hydrograph and pdm
    % hydrograph using formula (Q - P0) / (P - P0); in this parametrisation
    % the flow grows from 0 at t=0 to 1 for t->ininitiy
    benchmark_data = [t', (Q-P0) / (P-P0)];
    pdm_data = [hydrograph.t', (hydrograph.Q_total'-P0) / (P-P0)];
    
    % Export these hydrographs to .dat files if a directory is provided
    if nargin == 7
      writematrix(benchmark_data(1:10:end,:), ...
        [dir, 'benchmark_validation_', num2str(i), '.dat']);
      writematrix(pdm_data(1:100:end,:), ...
        [dir, 'pdm_validation_', num2str(i), '.dat'])
      writematrix(benchmark_data(benchmark_data(:,1)<=t_zoom,:), ...
        [dir, 'benchmark_validation_', num2str(i), '_zoomed.dat']);
      writematrix(pdm_data(pdm_data(:,1)<=t_zoom,:), ...
        [dir, 'pdm_validation_', num2str(i), '_zoomed.dat'])
    end
    
    % Plot graph in a new subplot
    h = subplot(2, N/2, i);
    
    % If a t_zoom option is used highlight zoomed area in yellow
    if nargin > 5
      patch([0 t_zoom t_zoom 0], [0, 0, 1, 1], 'y', 'FaceAlpha', 0.2, ...
        'edgecolor','none')
    end
    
    % Plot dimensionless hydrograph for physical benchmark model and PDM
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


% Function compare_hydrograph_components() compares the flow components for
% the hydrographs generated with a PDM simulation and physical benchmark
% model.
%
%   INPUT:
%     summary   - structure containing PDM model parameters
%                 (c_max, b, kg, st, kf, ks)
%     P0        - mean precipitation rate
%     P         - simulated precipitation rate
%     t         - time steps of the simulation
%     Q         - total flow in each time step according to physical model
%     Qs        - slow (subsurface) flow component
%     Qf        - fast (surface) flow component
%
%   OUTPUT:
%     phys_data - plotted data for the physical benchmark model
%     pdm_data  - plotted data for the PDM benchmark model

function [phys_data, pdm_data] = compare_hydrograph_components(...
                                            summary, P0, P, t, Q, Qs, Qf)

  % Initialize a PDM class object
  model = PDM();
  
  % Set its parameters to the ones provided in the 'summary' structure
  model = model.setParameters(summary);
  
  % Set initial conditions of the PDM model to the steady state for
  % the given mean precipitation rate P0
  model = model.setInitialCondition('steady state', P0);

  % Run PDM simulation for the specified time (with 10k time steps)
  [~, hydrograph] = model.simulate(P, t(end), 10000);
  
  % On the first subfigure plot total, slow and fast flow components
  % for the physical benchmark model
  subplot(1, 2, 1)
  plot(t / 24, Q) % we divide by 24 to convert time from [hours] to [days]
  hold on
  plot(t / 24, Qs)
  plot(t / 24, Qf)
  hold off
  
  % Set custom axes labels, legend and title
  xlabel('t [days]')
  ylabel('Q [m/s]')
  legend('total flow, Q', 'subsurface flow, Q_s', ...
    'surface flow, Q_f', 'Location', 'East')
  title('Physical banchmark')
  
  % On the second subfigure plot total, slow and fast flow components
  % for the PDM simulation
  subplot(1, 2, 2)
  plot(hydrograph.t / 24, hydrograph.Q_total)
  hold on
  plot(hydrograph.t / 24, hydrograph.Qs)
  plot(hydrograph.t / 24, hydrograph.Qf)
  hold off
  
  % Set custom axes labels and title 
  xlabel('t [days]')
  ylabel('Q [m/s]')
  title('Fitted PDM')
  
  % Set the size of the whole figure and its background color to white
  set(gcf, 'Position',  [100, 100, 800, 600])
  set(gcf,'color','w');

  % Create arrays for the data, that were used for plotting
  
  phys_data = [t'/24, Q, Qs, Qf];
  
  pdm_data = [hydrograph.t'/24, hydrograph.Q_total', ...
              hydrograph.Qs', hydrograph.Qf'];
end