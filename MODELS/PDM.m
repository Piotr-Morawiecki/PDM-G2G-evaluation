% ---------------------------
%
% Class name: PDM
%
% Purpose of class: Implements a probability distributed model (PDM) from
%                   paper "Calibration of a conceptual rainfall-runoff
%                   model for flood frequency estimation by continuous
%                   simulation" by R. Lamb (1999). For the details of this
%                   implementation follow references to our paper from
%                   the REAMDE file.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-04-02
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------

classdef PDM
  properties
    par   % all model parameters are stored in par structure
    C     % current soil storage capacity
    Sf    % current fast (surface) storage value
    Ss    % current slow (subsurface) storage value
  end
  methods
    
    %% Overview
    
    % The class includes the following functions for users to use
    % (detailed descriptions are provided below):
    %
    %   Preprocessing:    - setParameters
    %                     - setInitialCondition
    %
    %   Simulations:      - simulate
    %                     - calibrate
    
    %% setParameters
    
    % Function setParameters() sets model's parameters to ones given
    % in the par structure
    %
    % INPUT:
    %   par   either vector of six numeric parameters or a structure that
    %         should include:
    %
    %           c_max   maximal soil storage capacity,
    %           b       exponent characterising probability distribution,
    %           kg      time constant of the drainage function,
    %           st      threshold storage,
    %           kf      time constant of the fast store,
    %           ks      decay parameter of the slow store.
    %
    % OUTPUT:
    %   obj   PDM class object with added parameters
    
    function obj = setParameters(obj, par)
      try
        if length(par) > 1
        	% Set parameters if they are given as a vector
          obj.par.c_max = par(1);
          obj.par.b = par(2);
          obj.par.kg = par(3);
          obj.par.st = par(4);
          obj.par.kf = par(5);
          obj.par.ks = par(6);
        else
        	% Otherwise consider par as a structure
          obj.par.c_max = par.c_max;
          obj.par.b = par.b;
          obj.par.kg = par.kg;
          obj.par.st = par.st;
          obj.par.kf = par.kf;
          obj.par.ks = par.ks;
        end
        
        % Calculate maximum soil moisture storage (see our paper)
        obj.par.S_max = obj.par.c_max / (obj.par.b + 1);
        
      catch
        % If the parameters were not provided in the expected form display
        % the following error message:
        error('Not all required parameters are included in par structure.')
      end
    end
    
    
    %% setInitialCondition
    
    % Function setInitialCondition() sets values of each store following a
    % given initial condition.
    %
    % INPUT:
    %
    %   type    string specifying the type of initial condition; two
    %           options are available:
    %           1) 'dry'          - sets value of all stores to 0
    %           2) 'steady state' - sets value of all stores, so that the
    %                               system is in the steady state for
    %                               a given precipitation rate P0
    %   
    %   P0      precipitation rate value required when type='steady state'
    %           is picked
    %
    % OUTPUT:
    %
    %   obj     PDM class object with values of store set to satisfy
    %           a given initial condition
    
    function obj = setInitialCondition(obj, type, P0)
      
      if strcmp(type, 'dry')
        
        % If the initial condition is set to 'dry' set the value of
        % soil moisture, fast and slow store to 0.
        obj.C = 0;
        obj.Sf = 0;
        obj.Ss = 0;
      
      elseif strcmp(type, 'steady state')
        
        % Calculate what the soil moisture storage S should be to
        % give drainage d=(S-st)/kg equal to P0
        s0 = obj.par.st + obj.par.kg * P0;
        
        if s0 < obj.par.S_max
          % If this storage can be achieved (i.e. it does not exceed the
          % maximal value of S) calculate soil moisture storage capacity
          % corresponding to given S
          obj.C = obj.par.c_max * ...
            (1 - (1 - s0 / obj.par.S_max) ^ (1 / (obj.par.b + 1)));
        else
          % Otherwise set S and c to their maximum value
          s0 = obj.par.S_max;
          obj.C = obj.par.c_max;
        end
        
        % Calculate the drainage from the soil moisture store
        drainage = (s0 - obj.par.st) / obj.par.kg;
        
        % Calculate the slow store volume, that gives slow (subsurface)
        % flow qs=exp(Ss/ks) equal to the drainage rate
        obj.Ss = obj.par.ks * log(drainage);
        
        % Calculate the surface runoff
        runoff = P0 - drainage;
        
        % Calculate the slow store volume, that gives fast (surface)
        % flow qf=Sf/kf equal to the runoff rate
        obj.Sf = obj.par.kf * runoff;
      
      else
        % If the type of initial condition does not much any of the
        % settings display the following error message
        error('Unknown initial condition');
      end
    end
    
    
    %% Simulate
    
    % Function simulate() allows to run a time-dependent PDM simulation
    %
    % INPUT:
    % 
    %   P           simulated precipitation rate; either:
    %               - single value if precipitation rate is constant, or
    %               - array of length nt, with value of precipitation
    %                 given separately to each time step
    %   t_max       length of simulation
    %   nt          number of time steps
    %
    % OUTPUT:
    %
    %   solution    structure containing values of all stores in each time
    %               step
    %   hydrograph  structure containing values  of total flow, as well as
    %               its slow and fast components in each time steps
    %   obj         PDM class object with the final store values
    
    function [solution, hydrograph, obj] = simulate(obj, P, t_max, nt)
      
      % Calculate length of each time step
      dt = t_max / nt;
      
      % If precipitation is specified with a single value set precipitation
      % rate to be the same for all time steps
      if length(P) == 1
        P = P * ones(1, nt);
      end
      
      % Initialise structure to save store values (an extra entry is
      % included to store the initial values for each store)
      
      solution.C = zeros(1, nt+1);  % soil moisture storage capacity
      solution.S = zeros(1, nt+1);  % soil moisture storage value
      solution.Sf = zeros(1, nt+1); % fast storage value
      solution.Ss = zeros(1, nt+1); % slow storage value
      
      % Update stores in every time step
      for i = 1:nt+1
        
        % Calculate the soil mosisture value based on the current soil
        % moisture capacity
        if obj.C < obj.par.c_max
          S = obj.par.S_max * (1 - (1 - obj.C / obj.par.c_max) .^ (obj.par.b + 1));
        else
          S = obj.par.S_max;
        end
        
        % Save store values (for i=1 initial value is saved)
        solution.C(i) = obj.C;
        solution.S(i) = S;
        solution.Sf(i) = obj.Sf;
        solution.Ss(i) = obj.Ss;
        
        % If the last step was reached end the simulation
        if i == nt + 1
          break
        end
        
        % Calculate the drainage function
        drainage = max(0, (S - obj.par.st) / obj.par.kg);
        
        % Calculate Pi (source term in ODE for the soil moisture storage)
        Pi = P(i) - drainage;
        
        % Save old value of c (for runoff calculation)
        C_old = obj.C;
        
        % Update value of soil moisture capacity, c
        obj.C = min(obj.C + Pi * dt, obj.par.c_max);
        
        % Calculate the volume of the surface runoff from the water balance
        runoff = Pi * dt - obj.par.S_max * ...
          ((1 - C_old / obj.par.c_max) ^ (obj.par.b + 1) - ...
          (1 - obj.C / obj.par.c_max) ^ (obj.par.b + 1));
        
        % Calculare the rate of runoff generation
        runoff = runoff / dt;
        
        % Calculate slow and fast flow
        flow_s = exp(obj.Ss / obj.par.ks);
        flow_f = obj.Sf / obj.par.kf;
        
        % Update values of the slow and fast storage
        obj.Sf = obj.Sf + (runoff - flow_f) * dt;
        obj.Ss = obj.Ss + (drainage - flow_s) * dt;
      end
      
      % For each time step calculate the corresponding time (t),
      % fast flow (Qf), slow flow (Qs) and total flow (Q_total)
      hydrograph.t = linspace(0, t_max, nt+1);
      hydrograph.Qf = solution.Sf / obj.par.kf;
      hydrograph.Qs = exp(solution.Ss / obj.par.ks);
      hydrograph.Q_total = hydrograph.Qs + hydrograph.Qf;
    end
    
    
    %% calibrate
    
    % Function calibrate() finds value of PDM parameters that minimises
    % root mean square error (RMSE) between the PDM predictions and
    % the provided training data using MatLAB simplex search method.
    % Multiple hydrograph can be provided for the training set.
    % 
    % INPUT:
    %
    %   The training dataset has to be descibed by four structures, in
    %   which each element represents different training hydrograph. The
    %   required structures are:
    %
    %     P0_data   structure including the mean precipitation values
    %               (required to specify the initial condition)
    %     P_data    structure including the simulated precipitation values
    %     t_data    structure including the time for each time step
    %     Q_data    structure including the flow at the provided time steps
    %
    %   x0          (optional) vector of PDM parameters (if not specified, 
    %               the current settings assigned to the given PDM class 
    %               object are used); keep in mind that picking a
    %               relatively good initial condition allows to minimise
    %               the risk of function to find a suboptimal solution.
    %
    % OUTPUT:
    %   obj         PDM class object with the fitted parameters
    %   summary     structure with fitted PDM parameters values
    
    function [obj, summary] = calibrate(obj, P0_data, t_data, P_data, ...
        Q_data, x0)
      
      % If x0 is not provided use the ones assigned to the current PDM
      % object
      if nargin < 6
        x0 = [obj.par.c_max, obj.par.b, obj.par.kg, obj.par.st, ...
              obj.par.kf, obj.par.ks];
      end
      
      % Set maximal value of iterations of the min search method
      max_iterations = 100000;
      options = optimset('MaxIter', max_iterations, ...
        'MaxFunEvals', max_iterations);
      
      % Find PDM parameters minimising the accuracy function defined as
      % a RMSE between PDM predictions and provided training data.
      x = fminsearch(@(x) accuracy(x, obj, P0_data, P_data, t_data, ...
        Q_data), x0, options);
      
      % Construct a structure with fitted values of PDM parameters
      summary = struct('c_max', x(1), 'b', x(2), 'kg', x(3), ...
        'st', x(4), 'kf', x(5), 'ks', x(6));
      
      % Update the parameters of the current PDM object
      obj = obj.setParameters(x);
      
      % Function accuracy() returns a RMSE between the provided training
      % data and PDM predictions, which parameters given by 'x' argument.
      
      function [error] = accuracy(x, obj, P0_data, P_data, t_data, Q_data)
        
        % If at least one parameter is negative return infinite error (this
        % will force fminsearch solver to look for positive solutions only)
        if sum(x<0) > 0
          error = Inf;
          fprintf('fminsearch error: %e\n', error);
          return;
        end
        
        % Set PDM parameters to given x value
        obj = obj.setParameters(x);
        
        % Set error value to 0; in each loop add square errors
        % corresponding to each hydrograph in the training data set
        error = 0;
        for i = 1:length(P0_data)
          
          % Check length and duration of i-th hydrograph
          nt = length(t_data{i})-1;
          t_max = max(t_data{i});
          
          % Set PDM initial condition to a steady state corresponding to
          % precipitation rate P0
          obj = obj.setInitialCondition('steady state', P0_data{i});
          
          % Run a PDM simulation to predict the hydrograph
          [~, hydrograph] = obj.simulate(P_data{i}, t_max, nt);
          
          % Calculate mean square error and add it to the 'error' value
          error = error + mean((hydrograph.Q_total - Q_data{i}') .^ 2);
        end
        
        % Find RMSE from all hydrographs (note that each hydrograph has
        % equal contribution to the error, regardless of the duration of
        % this hydrograph)
        error = sqrt(error / length(P0_data));
        
        % Display the current error
        fprintf('fminsearch error: %e\n', error);
      end
    end
  end
end