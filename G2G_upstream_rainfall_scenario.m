% ---------------------------
%
% Script name: G2G_upstream_rainfall_scenario.R
%
% Purpose of script:  Script demonstrates differences in behaviour of the
%                     Grid-to-Grid model and physical benchmark model when
%                     applied to a scenario with a rainfall, which occurs
%                     only in the upstream part of catchment.
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
%  - MODELS/G2G.m     Grid-to-Grid Model class
%
% ---------------------------

%% Prepare the workspace

clear                           % clear the workspace

addpath(genpath('MODELS'))      % add MODELS folder and its subdirectories
                                % to the search path
                                
dir = 'DATA/G2G_upstream_rainfall_scenario/';   % output directory
if ~exist(dir, 'dir')
  mkdir(dir)                    % create the directory if it does not exist
end


%% Demonstration of the mechanism behind a shock formation

% In this script we use a characteristic method. We describe water height
% with characteristic curves that start from x0 (between 0 and 1) and with
% initial value h0(x0).

% Here we use a normal distribution curve for h0(x0)
x0 = linspace(-3,3,101);
h0 = exp(-x0.^2);

% In kinematic wave approximation we set constant wave speed, c
c = 1;

% In the figure the following five points are tracked
points = [21,36,51,66,81];

% with their location displayed at the following times
t_values = [0, 1.1, 2.2];

% Find number of time steps
nt = length(t_values);

% Prepare cell structre to record x(t) for all points in these time steps
x_data = cell(nt,2);

% Calculate propagation of characteristic curves
for i = 1:nt
  
  % kinematic wave approximation (K.W.A) with v=c=const.
  x_data{i,1}  = x0 - c * t_values(i);
  
  % solution following the Manning's law v=h^(2/3)
  x_data{i,2}  = x0 - h0 .^ (2/3) * t_values(i);
end

% Now 6 figures are plotted - 2 models with 3 time steps presented
hf = figure();
n = 0;
for i = 1:3
  for j = 1:2
    n = n + 1;
    subplot(3,2,n);
    
    % Choose data representing solution of j-th model in the i-th time step
    x = x_data{i,j};
    
    % Plot solution
    plot(x, h0)
    hold on
    
    % Highlight selected points
    plot(x(points), h0(points), '.', 'MarkerSize', 15)
    
    % Find speed of each point
    if j == 1
      v = c * ones(size(h0));   % kinematic wave with constant speed
    else
      v = h0 .^ (2/3);          % kinemaitc approximation of Manning's law
    end
    
    % Add arrows to represent the speed of each highlighted point
    for k = points
      h = annotation('arrow');
      h.Parent = hf.CurrentAxes;
      h.X = [x(k), x(k) - v(k)];
      h.Y = [h0(k), h0(k)];
    end
    hold off
    
    % Set axes limits
    xlim([-5,3])
    ylim([0,1.2])
  end
end

% Add red line representing shock formed in the 2nd model
hold on
xshock = -2.182;
plot([xshock, xshock], [0.01515, 0.999], '-r')
hold off

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [100, 100, 800, 400])
set(gcf, 'color', 'w');

% Export the figure to .png file
exportgraphics(gcf, 'FIGURES/G2G_shock_formation.png', 'Resolution', 300)

%% Compare the solution for PDM and benchmark physical model

% Here the solution for each model is plotted in 2D and 3D for a simulation
% in which there is not initial flow P0=0 and then rainfall starts to fall
% in the upstream area of the domain, x>0.5.

% Plot 1: 2D solution of the physical model (based on the Manning's law)
subplot(2,2,1)
plot_physical_benchmark_model(2, 20)
view(0,90)

% Plot 2: 2D solution of the Grid-to-Grid model (based on K.W.A.)
subplot(2,2,2)
plot_PDM_solution(20, 20)
view(0,90)

% Plot 3: 3D solution of the physical model (based on the Manning's law)
subplot(2,2,3)
plot_physical_benchmark_model(2, 20)
view(-25,20)

% Plot 4: 3D solution of the Grid-to-Grid model (based on K.W.A.)
subplot(2,2,4)
plot_PDM_solution(8, 20)
view(-25,20)

% Set the size of the whole figure and its background color to white
set(gcf, 'Position',  [100, 100, 800, 600])
set(gcf,'color','w');

% Export the figure to .png file
exportgraphics(gcf, 'FIGURES/G2G_shock_diagram.png', 'Resolution', 300)

%% Functions

% Function plot_physical_benchmark_model plots a solution of the physical
% model with the kinematic wave speed given by the Manning's law.

% INPUT:
%   n1  - number of characteristic curves reaching the end of the domain
%         before intersecting with a shock
%   n2  - number of characteristic curves starting from the area covered
%         with the rainfall (x>0.5)

function plot_physical_benchmark_model(n1, n2)

  k = 5/3;          % exponent from the Manning's law
  
  t_max = 1;        % simulated time
  
  nt = 100;         % number of time steps
                    % (since method of characteristic is analytic reducing
                    % this number does not decrease the accuracy, however
                    % reduces the quality of the output figure)
  
  dt = t_max / nt;  % duration of a single time step

  r = 1;            % precipitation rate
  
  xr = 0.5;         % the precipitation is nonzero only for x>xr
  
  % This is equation for the heigh h(x,t) of the water at given position x
  % and time t (up to the moment when the characteristic line intersects
  % with the shock)
  h = @(x,t) fzero(@(h) x - xr + k * h^(k-1) * (t - h/r), [(k-1)/k*r*t, r*t]);

  % Now we want to find location of the shock x(t) and the flow
  % corresponding to this point by solving ODE derived from the
  % Rankine-Hugoniot conditions:
  
  %                       dx/dt = -h^(k-1)
  
  % Initialise vectors to store values t, x(t), q(x,t) at the shock
  x_shock = zeros(1, nt);
  t_shock = linspace(0, t_max, nt+1);
  q_shock = zeros(1, nt+1);
  
  % The shock starts from x=xr and q=0
  x_shock(1) = xr;
  q_shock(1) = 0;
  
  % At each timestep:
  for i = 1:nt
    
    % 1) find height of the shock, h(x,t)
    h_shock = h(x_shock(i), t_shock(i));
    
    % 2) find the corresponding, q=h^k
    q_shock(i) = h_shock .^ k;
    
    % 3) Apply the Rankine-Hugoniot condition to find new x_shock (using
    %    Euler method)
    x_shock(i+1) = x_shock(i) - h_shock .^ (k - 1) * dt;
  end
  
  % Add the last element to q_shock vector
  q_shock(i+1) = h(x_shock(i+1), t_shock(i+1)) .^ k;

  % Now the characteristic curves will be computed and plotted
   
  r = @(x) x>xr;    % function descrining the precipitation rate r(x)
  
  % The n2 curves start from upstream part of the domain [xr,1]
  x0_values = linspace(xr,1,n2+1);
  x0_values = x0_values(2:end);   % no curve is plotted at the location
                                  % of the shock wave
  
  % For each starting point plot a characteristic curve starting from a
  % given point at time t=0 and ending at the location of the shock (or at
  % t=t_max, whichever is first)
  hold on
  for x0 = x0_values
    
    % Set ODE options; event checks whether the characteristic curve
    % crossed the shock
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, ...
      'Events', @(t,y) event(t, y, t_shock, x_shock));
    
    % Solve characteristic eqautions (Lagrange–Charpit equations)
    [t, y, te, ye] = ode45(@(t,y) [-k*y(2).^(k-1); r(y(1))], ...
      [0, t_max], [x0, 0], opts);
    
    % Plot the resulting characteristic curve
    plot3(y(:,1), t, y(:,2).^k, 'Color', 'c')
    
    % If the curve intersected the shock add another characteristic curve
    % (in yellow) reaching the shock at the same location, but starting
    % from the downstream part of the domain (x<xr)
    if length(te) == 1
      last_ye = ye(1);
      plot3([ye(1), ye(1)], [0, te], [0, 0], 'Color', 'y')
      
      % Also add a red vertical dashed line to represent the discontinuity
      % at this location
      plot3([ye(1), ye(1)], [te, te], [0, ye(2).^k], '--', ...
        'Color', [1,0.8,0.8])
    end
  end
  
  % Add n1 characteristic curves starting from x values to which the shock
  % has not arrived before t_max
  x0_values = linspace(0, last_ye, n1+1);
  for x0 = x0_values(1:end-1)
    plot3([x0, x0], [0, t_max], [0, 0], 'Color', 'y')
  end

  % Add red line to represent the shock front
  plot3(x_shock, t_shock, q_shock, 'Color', 'red');
  plot3(x_shock, t_shock, 0 * q_shock, 'Color', 'red');
  hold off
  
  % Set custom axes limits and labels
  xlim([0,1])
  ylim([0,t_max])
  xlabel('x')
  ylabel('t')
  zlabel('h')
end


% Function plot_PDM_solution plots a solution of the Grid-to-Grid
% model with the kinematic wave approximation with a constant wave speed.

% INPUT:
%   n1  - number of characteristic curves starting from the area not
%         covered with the rainfall (x<0.5)
%   n2  - number of characteristic curves starting from the area
%         covered with the rainfall (x>0.5)

function plot_PDM_solution(n1, n2)

  t_max = 1;        % simulated time
  
  xr = 0.5;         % the precipitation is nonzero only for x>xr
  
  r = @(x,t) x>xr;  % function descrining the precipitation rate r(x)
  
  
  % Choose starting points for the characteristic curves
  x0_values = linspace(xr,1,n2);
  x0_values = x0_values(2:end);
  x0_values = [linspace(0,xr,n1), x0_values];
  
  % Compute and plot characteristic curves for each staring location
  hold on
  for x0 = x0_values
    
    % Set ODE options
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    
    % Solve characteristic eqautions (Lagrange–Charpit equations)
    [t, y] = ode45(@(t,y) [-1; r(y(1),t)], [0, t_max], [x0, 0], opts);
    
    % Depending on whether the curve started from initially saturated zone
    % or not plot in in yellow or cyan
    if x0 < xr
      plot3(y(:,1), t, y(:,2), 'Color', 'y')
    else
      plot3(y(:,1), t, y(:,2), 'Color', 'c')
    end
  end
  hold off
  
  % Set custom axes limits and labels
  xlim([0,1])
  ylim([0,t_max])
  xlabel('x')
  ylabel('t')
  zlabel('h')
end


% Function event checks whether the solution already crossed the shock
% front (its value for given time t is obtained via a linear
% interpolation); see description of required output parameters
% in the 'ode45' description in the MatLAB documentation

function [value,stop,direction] = event(t, y, t_shock, x_shock)
  value = y(1) - interp1(t_shock, x_shock, t);
  stop = 1;
  direction = 0;
end