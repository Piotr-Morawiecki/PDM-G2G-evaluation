% ---------------------------
%
% Script name: PDM_timescales.R
%
% Purpose of script:  Plots simulation result of PDM model together with
%                     the characteristic timescales characteristing
%                     different phases of the hydrograph formation
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
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
dir = 'DATA/PDM_timescales/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end

%% Prepare simulation settings

% Here we create a par structure, which will include all parameters
% required by PDM class and further postprocessing of simulation results.

% Choose PDM model parameters
par.c_max = 0.002;
par.b = 1;
par.kg = 20;
par.st = 0;
par.kf = 0.1;
par.ks = 0.001;

par.t_max = 100;  % simulated time
par.nt = 10000;   % number of timesteps

par.P0 = 2e-05;   % mean precipitation (for the initial condition)
par.P = 2e-04;    % peak precipitation (for the simulated rainfall)

% Calculate maximum value of soil moisture store and drainage rate
par.S_max = par.c_max / (par.b + 1);
par.d_max = (par.S_max - par.st) / par.kg;

% Calculate characteristic timescale of soil moisture store (tc),
% fast store (tf) and slow store (ts)
par.tc = par.c_max / par.P;
par.tf = par.kf;
par.ts = par.ks / par.d_max;

% Display all model parameters
disp(par)

%% Run the Probability Distributed Model (PDM) simulation

%Initialize PDM class object
model = PDM();

% Set its parameters to the ones defined in par structure
model = model.setParameters(par);

% Set the initial condition to be a steady state corresponding to the
% mean rainfall P0
model = model.setInitialCondition('steady state', par.P0);

% Run PDM simulation for rainfall P of duration t_max with nt time steps
[solution, hydrograph] = model.simulate(par.P, par.t_max, par.nt);

%% Exporting simulation results

% We will plot the simulation data with a time presented on the log axis
% with 100 datapoints equally distributed (in log space) between t=1e-2
% and t=log10(par.t_max).
t = logspace(-2, log10(par.t_max), 100);

% Interpolate the simulation results to find storage values and flow
% components at specified times

solution.S = interp1(hydrograph.t, solution.S, t);    % soil moisture store
solution.Sf = interp1(hydrograph.t, solution.Sf, t);  % fast store
solution.Ss = interp1(hydrograph.t, solution.Ss, t);  % slow store

hydrograph.Qf = interp1(hydrograph.t, hydrograph.Qf, t);        % fast flow
hydrograph.Qs = interp1(hydrograph.t, hydrograph.Qs, t);        % slow slow
hydrograph.Q_total = interp1(hydrograph.t, hydrograph.Q_total, t);  % total

% Export the simulation results to a .dat file
writematrix([t', solution.S', solution.Sf', solution.Ss', ...
  hydrograph.Qf', hydrograph.Qs', hydrograph.Q_total'], ...
  [dir, 'PDM_timescales_plot_data.dat'])

%% Plotting simulation results

% Plot soil moisture store variation with time
subplot(3,2,1)
semilogx(t, solution.S, 'LineWidth', 2)
xline(par.tc, '--')             % add line representing characteristic time
text(par.tc*1.2, 6e-4, 't_c')   % of the soil moisture store
xlabel('t')
ylabel('soil storage, S')

% Plot fast (surface) store variation with time
subplot(3,2,3)
semilogx(t, solution.Sf, 'LineWidth', 2)
xline(par.tf, '--')             % add line representing characteristic time
text(par.tf*1.2, 0.5e-5, 't_f') % of the fast store...
xline(par.tc, '--')             % ... and the soil moisture store
text(par.tc*1.2, 0.5e-5, 't_c')
xlabel('t')
ylabel('fast storage, S_f')

% Plot slow (subsurface) store variation with time
subplot(3,2,5)
semilogx(t, solution.Ss, 'LineWidth', 2)
xline(par.ts, '--')               % add line representing characteristic
text(par.ts*1.2, -0.0107, 't_s')  % time of the slow store
xlabel('t')
ylabel('slow storage, S_s')

% Plot hydrograph presenting variation of fast, slow and total flow in time
subplot(1,2,2)
semilogx(t, hydrograph.Qf, 'LineWidth', 2)
hold on
plot(t, hydrograph.Qs, 'LineWidth', 2)
plot(t, hydrograph.Q_total, 'LineWidth', 2)
plot(t, par.P * ones(size(t)), '--k')

% Add lines representing each of three timescales
xline(par.tf, ':')
text(1.2*par.tf, 0.7e-4, 't_f')
xline(par.ts+par.tc, ':')
text(1.2*par.ts+par.tc+2, 0.7e-4, 't_c+t_s')
xline(par.tf+par.tc, ':')
text(1.2*par.tf+par.tc+2, 0.7e-4, 't_c+t_f')
hold off

% Configure other plot settings and add legend
xlabel('time')
ylabel('flow')
ylim([0,2.2e-4])
legend('overland flow, Q_f(t)', 'groundwater flow, Q_s(t)', 'total flow, Q(t)', ...
  'total precipitation, P', 'Location', 'NorthWest')

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [100, 100, 1200, 500])
set(gcf,'color','w');

% Export the figure to a .png file
exportgraphics(gcf, 'FIGURES/PDM_timescales.png', 'Resolution', 300)