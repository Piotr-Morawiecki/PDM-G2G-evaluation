% ---------------------------
%
% Class name: GridToGrid
%
% Purpose of class: Implements a Grid-to-Grid (G2G) model from the paper
%                   "Development of a high resolution grid-based river flow
%                   model for use with regional climate model output" by
%                   V. Bell(2007), applied to a one-dimensional benchmark
%                   catchment model. For the details of this implementation
%                   follow references to our paper from the REAMDE file.
%
% Author: Piotr Morawiecki
%
% Date Created: 2023-04-02
%
% Copyright (c) Piotr Morawiecki, 2023
% Email: pwm27@bath.ac.uk
%
% ---------------------------

classdef GridToGrid
  properties
    par   % all model parameters are stored in par structure
    C     % current soil storage capacity
    Q     % surface flow (array of flow values for each mesh element)
    Qb    % subsurface flow (array of flow values for each mesh element)
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
    %           c       speed of the overland flow
    %           cb      speed of the subsurface flow
    %           r       return flow fraction
    %           c_max   maximal soil storage capacity,
    %           b       exponent characterising probability distribution,
    %           k       storage rate constant
    %           beta    drainage exponent (drainage = k S^beta)
    %           Lx      length of the domain (hillslope)
    %           nx      number of mesh elements
    %
    % OUTPUT:
    %   obj   GridToGrid class object with added parameters
    
    function obj = setParameters(obj, par)
      try
        if length(par) > 1
        	% Set parameters if they are given as a vector
          obj.par.c = par(1);
          obj.par.cb = par(2);
          obj.par.r = par(3);
          obj.par.c_max = par(4);
          obj.par.b = par(5);
          obj.par.k = par(6);
          obj.par.beta = par(7);
          
          % The last three parameters are not changed during the
          % calibration, therefore calibrate() function calls
          % setParameters() without these three arguments. User should
          % always define these parameters when creating a new G2G object.
          if length(par) > 7
            obj.par.nx = par(8);
            obj.par.dx = 1 / obj.par.nx;
          end
        
        else
        	% Otherwise consider par as a structure
          obj.par.c = par.c;
          obj.par.cb = par.cb;
          obj.par.r = par.r;
          obj.par.c_max = par.c_max;
          obj.par.b = par.b;
          obj.par.k = par.k;
          obj.par.beta = par.beta;
          obj.par.nx = par.nx;
          obj.par.dx = 1 / obj.par.nx;
        end
        
        % Calculate maximum soil moisture storage and corresponding
        % maximum drainage rate
        obj.par.S_max = obj.par.c_max / (obj.par.b + 1);
        obj.par.d_max = obj.par.k * obj.par.S_max ^ obj.par.beta;
        
      catch
        % If the parameters were not provided in the expected form display
        % the following error message:
        error('Not all required parameters were included in par structure.')
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
    %   obj     GridToGrid class object with values of store set to satisfy
    %           a given initial condition
    
    function obj = setInitialCondition(obj, type, P0)
      
      if strcmp(type, 'dry')
        % If the initial condition is set to 'dry' set the value of
        % soil moisture, fast and slow stores to 0.
        obj.C = 0;
        obj.Q = zeros(1, obj.par.nx);
        obj.Qb = zeros(1, obj.par.nx);
        
      elseif strcmp(type, 'steady state')
        
        % Calculate what the soil moisture storage S should be to
        % give drainage d=(S-st)/kg equal to P0
        s0 = (P0 / obj.par.k) ^ (1 / obj.par.beta);
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
        drainage = obj.par.k * s0 ^ obj.par.beta;
        
        % Calculate the surface runoff
        runoff = P0 - drainage;
        
        % now calculate the steady state values of Q(x) and Qb(x);
        % start by initializing these arrays
        obj.Q = zeros(1, obj.par.nx);
        obj.Qb = zeros(1, obj.par.nx);
        
        % The algorithm calculated q and qb iteratively for each cell
        % starting from the one located furthers for the river at the
        % catchment boundary - therefore no flow goes into this cell:
        q = 0;
        qb = 0;
        
        % Steady state is found using a reccurence relations:
        %
        %   qb(x+1) = (qb(x) + drainage * dx) / (1 + r * dx / cb),
        %   q(x+1) = q(x) + runoff * dx + r * dx / cb * qb(x).
        %
        
        r = obj.par.r * obj.par.dx / obj.par.cb;
        for i = 1:obj.par.nx
          qb = qb / (1 + r) + drainage * obj.par.dx / (1 + r);
          q = q + runoff * obj.par.dx + r * qb;
          obj.Qb(i) = qb;
          obj.Q(i) = q;
        end
        
      else
        % If the type of initial condition does not much any of the
        % settings display the following error message
        error('Unknown initial condition');
      end
    end
    
    
    %% Simulate
    
    % Function simulate() allows to run a time-dependent Grid-to-Grid
    % simulation
    %
    % INPUT:
    % 
    %   P           simulated precipitation rate; either:
    %               - single value if precipitation rate is constant, or
    %               - array of length nt, with value of precipitation
    %                 given separately to each time step
    %   t_max       length of simulation
    %   dt          number of time steps
    %
    % OUTPUT:
    %
    %   solution    structure containing values of all stores in each time
    %               step
    %   hydrograph  structure containing values  of total flow, as well as
    %               its slow and fast components in each time steps
    %   obj         G2G class object with the final store values
    
    function [solution, hydrograph, obj] = simulate(obj, P, t_max, nt)
      
      % Calculate length of each time step
      dt = t_max / nt;
      
      % Calculated dimenionless speed of the surface and subsurface flow
      theta = obj.par.c / obj.par.dx * dt;
      theta_b = obj.par.cb / obj.par.dx * dt;
      
      % If precipitation is specified with a single value set precipitation
      % rate to be the same for all time steps
      if length(P) == 1
        P = P * ones(1, nt);
      end
      
      % Initialise structure to save store values (an extra entry is
      % included to store the initial values for each store)
      solution.C = zeros(nt+1, 1);
      solution.Q = zeros(nt+1, obj.par.nx);
      solution.Qb = zeros(nt+1, obj.par.nx);
      
      % Update stores in every time step
      for i = 1:nt+1
        
        % Save store values (for i=1 initial value is saved)
        solution.C(i) = obj.C;
        solution.Q(i,:) = obj.Q;
        solution.Qb(i,:) = obj.Qb;
        
        % If the last step was reached end the simulation
        if i == nt + 1
          break
        end
        
        % Calculate the soil mosisture value based on the current soil
        % moisture capacity
        if obj.C < obj.par.c_max
          S = obj.par.S_max * (1 - (1 - obj.C / obj.par.c_max) .^ (obj.par.b + 1));
        else
          S = obj.par.S_max;
        end
        
        % Calculate the drainage function
        drainage = obj.par.k * S .^ obj.par.beta;
        
        % Calculate Pi (source term in ODE for the soil moisture storage)
        Pi = P(i) - drainage;
        
        % Save old value of c (for runoff calculation)
        C_old = obj.C;
        
        % Update value of soil moisture capacity, c
        obj.C = obj.C + Pi * dt;
        
        if obj.C > obj.par.c_max
          obj.C = obj.par.c_max;
        elseif obj.C < 0
          obj.C = 0;
        end
        
        % Calculate the volume of the surface runoff from the water balance
        runoff = Pi * dt - obj.par.S_max * ...
          ((1 - C_old / obj.par.c_max) ^ (obj.par.b + 1) - ...
          (1 - obj.C / obj.par.c_max) ^ (obj.par.b + 1));
        
        % Calculare the rate of runoff generation
        if runoff < 0
          runoff = 0;
        end
        runoff = runoff / dt;
        
        % Calculate the return flow
        R = obj.par.r * obj.Qb * obj.par.dx / obj.par.cb;
        
        % Update surface flow
        Q_upstream = [0, obj.Q(1:end-1)];
        obj.Q = (1 - theta) * obj.Q + theta * ...
          (Q_upstream + runoff * obj.par.dx + R);
        
        % Update subsurface flow
        Q_upstream = [0, obj.Qb(1:end-1)];
        obj.Qb = (1 - theta_b) * obj.Qb + theta_b * ...
          (Q_upstream + drainage * obj.par.dx - R);
      end
      
      % For each time step calculate the corresponding time (t),
      % surface flow at the river (Q), subsurface flow at the river (Qs),
      % and total river inflow (Q_total)
      hydrograph.t = linspace(0, t_max, nt+1);
      hydrograph.Q = solution.Q(:,end)';
      hydrograph.Qb = solution.Qb(:,end)';
      hydrograph.Q_total = hydrograph.Q + hydrograph.Qb;
    end
    
    
    %% calibrate
    
    % Function calibrate() finds value of Grid-to-Grid parameters that
    % minimises root mean square error (RMSE) between the G2G predictions
    % and the provided training data using MatLAB simplex search method.
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
    %   x0          (optional) vector of Grid-to-Grid parameters (if not 
    %               specified, the current settings assigned to the given
    %               GridToGrid class object are used); keep in mind that
    %               picking a relatively good initial condition allows to
    %               minimise the risk of function to find a suboptimal
    %               solution.
    %
    %   weights     (optional) array of weigths that apply to errors
    %               calculated based on each training hydrograph
    %
    %   max_iterations  (optional) max number of iteration of the
    %                   minimisation alogrithm
    %
    % OUTPUT:
    %   obj         GridToGrid class object with the fitted parameters
    %   summary     structure with fitted G2G parameters values
    
    
    function [obj, summary] = calibrate(obj, P0_data, t_data, P_data, ...
        Q_data, x0, weights, max_iterations)
      
      % If x0 is not provided use the ones assigned to the current
      % GridToGrid object
      if nargin < 6
        x0 = [obj.par.c, obj.par.cb, obj.par.r, obj.par.c_max, ...
          obj.par.b, obj.par.k, obj.par.beta];
      end
      
      % If weigths are not assigned, set all weights to 1
      % (i.e. treat all hydrograph in the training set equally)
      if nargin < 7
        weights = ones(1, length(P0_data));
      end
      
      % If max number of iterations is not specified set it to 100k
      if nargin < 8
        max_iterations = 100000;
      end
      
      % Set maximal value of iterations of the min search method
      options = optimset('MaxIter', max_iterations, ...
        'MaxFunEvals', max_iterations);
      
      % Find Grid-to-Grid model parameters minimising the accuracy function
      % defined as a RMSE between the Grid-to-Grid model predictions and
      % the provided training data, with additional constraint that fitted
      % parameters should be between 0 and infinity.
      
      x = fminsearch(@(x) accuracy(x, obj, P0_data, t_data, P_data, ...
        Q_data, weights), x0, options);
      
      %x = fmincon(@(x) accuracy(x, obj, P0_data, t_data, P_data, ...
      %  Q_data, weights), x0,[],[],[],[], zeros(size(x0)), ...
      %  Inf(size(x0)), [], options);
      
      % Construct a structure with fitted values of G2G parameters
      summary = struct('c', x(1), 'cb', x(2), 'r', x(3), ...
        'c_max', x(4), 'b', x(5), 'k', x(6), 'beta', x(7));
      
      % Update the parameters of the current GridToGrid object
      obj = obj.setParameters(x);
      
      % Function accuracy() returns a RMSE between the provided training
      % data and G2G predictions, which parameters given by 'x' argument.
      
      function [error] = accuracy(x, obj, P0_data, t_data, P_data, ...
          Q_data, weights)
        
        % If at least one parameter is negative return infinite error (this
        % will force fminsearch solver to look for positive solutions only)
        if sum(x<0) > 0
          error = Inf;
          fprintf('fminsearch error: %e\n', error);
          return;
        end
        
        % Set G2G parameters to given x value
        obj = obj.setParameters(x);
        
        % Set error value to 0; in each loop add square errors
        % corresponding to each hydrograph in the training data set
        error = 0;
        for i = 1:length(P0_data)
          
          % Check length and duration of i-th hydrograph
          nt = length(t_data{i})-1;
          t_max = max(t_data{i});
          
          % Set G2G initial condition to a steady state corresponding to
          % precipitation rate P0
          obj = obj.setInitialCondition('steady state', P0_data{i});
          
          % Run a G2G simulation to predict the hydrograph
          [~, hydrograph] = obj.simulate(P_data{i}, t_max, nt);
          
          % Calculate mean square error and add it to the 'error' value
          error = error + weights(i) * ...
            mean((hydrograph.Q_total - Q_data{i}') .^ 2);
        end
        
        % Find RMSE from all hydrographs (note that each hydrograph has
        % equal contribution to the error, regardless of the duration of
        % this hydrograph)
        error = sqrt(error / sum(weights));
        
        % Display the current error
        fprintf('fminsearch error: %e\n', error);
      end
    end
  end
end