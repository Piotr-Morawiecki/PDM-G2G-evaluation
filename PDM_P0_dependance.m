% ---------------------------
%
% Script name: PDM_P0_dependance.R
%
% Purpose of script:  Plots PDM simulation results calibrated to the
%                     physical benchmark model results for different
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
%  - MODELS/PDM.m     Probability Distributed Model class
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

dir = 'DATA/PDM_P0_dependance/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Set training and validation set

% In this script a probability distributed model (PDM) will be fitted to
% results of the physical benchmark model. For this purpose two data sets
% are created:
%   - training data set, which is used to calibrate PDM parameters, and
%   - validation data set, used to assess the accuracy of the fitted PDM.

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

% We start calibration process by guessing relatively good starting PDM
% parameters

%     c_max   b   kg   st   kf    ks
x0 = [ 0.03,  1,  720,  0,  1,  0.001];

% The quality of this initial guess can be assessed using
% compare_hydrographs() function.
figure(1)
compare_hydrographs(x0, P0_test_data, P_test_data, t_test_data, ...
  Q_test_data, validation_data);

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
  [~, summary] = model.calibrate(P0_test_data, ...
    t_test_data, P_test_data, Q_test_data, x0);

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
[benchmark_data, pdm_data] = compare_hydrographs(summary, ...
  P0_valid_data, P_valid_data, t_valid_data, Q_valid_data, ...
  validation_data);

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/PDM_P0_Validation_set.png', 'Resolution', 300)

% Export the plotted data to .dat files (every 10th row to reduce the
% output data size)
writematrix(pdm_data(1:10:end,:), [dir, 'pdm_r0_dependence.dat']);
writematrix(benchmark_data(1:10:end,:), ...
  [dir, 'benchmark_r0_dependence.dat']);


%% Plot dependance of a0 on P0 for PDM and physical benchmark model

% Here we plot this relationship for the following PDM parameters:
%    c_max   b   kg   st  kf   ks
x0 = [0.03, 2, 210, 0.002, 1, 0.001];
model = PDM().setParameters(x0);

% and physical benchmark model with the default parameters:
load('scenario_B_settings.mat', 'settings');

% Plot the size of the saturated zone a0(P0) as a function of mean
% precipitation rate P0 for PDM and physical benchmark model
[data, arrow_data] = compare_a0_dependance(model, settings);

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/PDM_a_shape.png', 'Resolution', 300)

% Export the plotted data to a .dat file
writematrix(data, [dir,'a_shape_data.dat']);
writematrix(arrow_data, [dir, 'arrow_data.dat']);


%% Functions

% Function compare_hydrographs() generates a figure, on which the solutions
% of the PDM simulation with the physical benchmark model for the provided
% scenarios.
%
% INPUT:
%   summary   - structure containing PDM model parameters
%               (c_max, b, kg, st, kf, ks)
%   P0_data   - mean precipitation rate
%   P_data    - simulated precipitation rate
%   t_data    - time steps of the simulation
%   Q_data    - measured total flow in each time step
%   a0_values - values of a0 parameter to include in the legend
%
% OUTPUT:
%   benchmark_data  - data plotted for the physical benchmark model
%   pdm_data        - data plotted for the PDM

function [benchmark_data, pdm_data] = compare_hydrographs(summary, ...
  P0_data, P_data, t_data, Q_data, a0_values)

  % Initialize a PDM class object
  model = PDM();
  
  % Set its parameters to the ones provided in the 'summary' structure
  model = model.setParameters(summary);

  % Set colors in which results of each scenario will be represented
  colors = {'r', 'g', 'b', 'c', ' m'};
  
  % Initialize a cell array to store legend labels
  legend_entries = cell(1, length(P0_data));
  
  % Create an empty plot to which results for each validation scenario will
  % be added
  subplot(1,1,1)
  hold on
  
  % For each scenario, PDM and physical benchmark model predictions
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
    
    % Run a time-dependent PDM simulation for given precipitation P,
    % of the same length as in case of physical benchmark, t(end) with
    % 10000 time steps
    [~, hydrograph] = model.simulate(P_data{i}, t(end), 10000);

    % If this is the first hydrograph initialise two data arrays to store
    % data plotted for PDM and physical benchmark model; otherwise add
    % this data to already existing arrays.
    if i == 1
      benchmark_data = [t' / 3600, Q];
      pdm_data = [hydrograph.t', hydrograph.Q_total'];
    else
      benchmark_data = [benchmark_data, Q];
      pdm_data = [pdm_data, hydrograph.Q_total'];
    end

    % Plot hydrograph for physical benchmark model (solid line)
    plot(t, Q, colors{i})
    
    % and PDM simulation (dashed line)
    plot(hydrograph.t, hydrograph.Q_total, ['--', colors{i}], ...
      'HandleVisibility', 'off');
    
    % Add corresponding legend entry to be displayed later
    legend_entries{i} = ['a_0=', num2str(a0_values(i),2)];
  end
  hold off
  
  % Set custom x-axis limit, axes ticks and labels
  xlim([0,24])
  xticks(0:6:24)
  xlabel('Time, t [h]')
  ylabel('Total flow, Q')
  
  % Add legend
  legend(legend_entries, 'Location', 'southeast')

  % Set the size of the whole figure and its background color to white
  set(gcf, 'Position',  [100, 100, 600, 350])
  set(gcf, 'color', 'w');
end


%% compare_a0_dependance

% Function compare_a0_dependance() plots the size of the saturated zone
% a0(P0) as a function of mean precipitation P0 according to the PDM and
% physical benchmark model.
%
% INPUT:
%   pdm         PDM class object with PDM model parameters
%   settings    structure with the physical model parameters
%
% OUTPUT:
%   data        plotted data for both models
%   arrow_data  coordinates of arrows, which represent difference between
%               these two models

function [data, arrow_data] = compare_a0_dependance(pdm, settings)
  % Based on our publication (see README file for reference) size of the
  % saturated zone for PDM model is expressed as
  %
  %              a0(P0) = 1 - (C1 - C2 * P0) ^ (b / (b+1))
  % where:
  C1 = 1 - pdm.par.st / pdm.par.S_max;
  C2 = pdm.par.kg / pdm.par.S_max;
  a_pdm = @(P0) 1 - (C1 - C2 * P0) .^ (pdm.par.b / (pdm.par.b+1));

  % The saturated zone size in the physical benchmark model is given by:
  a_benchmark = @(P0) 1 - settings.K * settings.Sx * settings.Lz ./ P0;

  % We plot a0 for P0 in range from 0 to 5e-5 (with 100 plotted points) 
  P0_values = linspace(0,5e-5,100);
  
  % Put data used for plotting in an array
  data = [P0_values', a_pdm(P0_values)', a_benchmark(P0_values)'];

  % Since a0 can be only between 0 and 1, any value of the functions
  % outside this range is set to -0.01 and 1.01 (i.e. slightly below 0 and
  % 1, so that these data points will be outside the plotted area).
  data(data<0) = 0;
  data(data>1) = 1;

  % Plot a_benchmark(P0) and a_pdm(P0) graph
  hf = figure();
  plot(P0_values, data(:,3), 'LineWidth', 2)
  hold on
  plot(P0_values, data(:,2), 'LineWidth', 2);

  % Find intersection between a_benchmark(P0) and a_pdm(P0) functions
  P0_points = [fzero(@(P0) a_pdm(P0) - a_benchmark(P0), 1e-5), ...
               fzero(@(P0) a_pdm(P0) - a_benchmark(P0), 3.5e-5)];

  % Put black points at the location of these intersections
  plot(P0_points, a_benchmark(P0_points), 'ok', 'MarkerFaceColor', 'k')

  % Additionally arrows are plotted between the graph of a_benchmark(P0)
  % and a_pdm(P0) functions. 20 arrows are plotted between P0=K*Sx*Lz, where
  % a_benchmark(P0)=0, and P0=C1/C2, where a_pdm(P0)=1.  
  %P0_arrows = linspace(settings.K * settings.Sx * settings.Lz, C1/C2, 20);
  P0_arrows = linspace(0, 5e-5, 30);

  % Arrow coordinates are save in 'arrow_data' array for exporting
  arrow_data = [];

  % The arrows are added in a for loop
  for p = P0_arrows

    % The arrows are added only if the separation between the lines is wide
    % enough to fit them
    if abs(a_benchmark(p) - a_pdm(p)) > 0.05

      % Add red arrow between the function graphs
      h = annotation('arrow', 'Color', 0.8*[1,1,1]);
      h.Parent = hf.CurrentAxes;
      h.X = [p, p];
      h.Y = [max(0, a_benchmark(p)), min(1, a_pdm(p))];

      % Add a small gap between the arrow and function graphs
      if a_benchmark(p) > a_pdm(p)
        h.Y = h.Y + 0.015 * [-1 1];
        h.Color = 'red';
      else
        h.Y = h.Y + 0.015 * [1 -1];
        h.Color = 'blue';
      end

      % Add arrow coordinates to the 'arrow_data' array
      arrow_data = [arrow_data, h.X', h.Y'];
    end
  end
  hold off

  % Add legend, labels and axes limits
  legend('Physical benchmark', 'PDM', 'Location', 'SouthEast')
  ylabel('Size of the saturated zone, a(P_0)')
  xlabel('Mean precipitation P_0')
  xlim([0,5e-5])
  ylim([0,1])

  % Set the size of the whole figure and its background color to white
  set(gcf, 'Position',  [100, 100, 800, 400])
  set(gcf,'color','w');
end