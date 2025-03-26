%% Cost function used by fmincon defined here - input vector information, return a scalar
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
% Hard coded in the 4parameter ODE - pass this through in the future
% arguments:
% parval:    current value of the parameter that is being fitted
% dtimes:    data times: the times at which the data are measured
% alltime:   the time range over which the entire simulation is being run
% A0:        initial state value of A for entire simulation
% data:      the data that are being fit
% all_pars:  a vector containing values for all of the parameters
% which_par: a vector containing the indices of the parameters that are being fit
function C = cost(parval, dtimes, alltime, A0, data, all_pars, which_par) %Calls Model2 ODE
    
        par = all_pars;

        par(which_par) = parval;
    

    % initial value of A taken from the first data point
        Ainit = A0;
        Vinit = 0; %We always assume no vaccine to start

        y0 = [Ainit Vinit];
    % run the ode, and output solution at times given by the data
    % change the name of the ode file if fitting a different model
    %     [T,Y] = ode45(@NAb_VaccineODE,times,y0,[],par);
    
    % Ensure thatt MaxStep = 0.1 to catch treatments, and that solutions cannot be negative
        myOPTIONS = odeset('MaxStep',0.01,'NonNegative',1);
    
        sol = ode15s(@NAb_VaccineODE4Param,alltime,y0,myOPTIONS,par);

    
    % calculate the norm (square root of sum of squares) of the different
    % between solution and data

    % use all the data
        w = ones(size(data')); %row vector


% %     % more weight on final points
        w(1)     = 3;    %For model 2, give weight to first data point.
        w(2)     = 3; 
        w(end)   = 3;
     
    
        ydiscrete = deval(sol,dtimes,1); % Extract the first component (NAbs) of the solution at specified times
        ysmooth = deval(sol,alltime,1); % Create a smoothed version of the solution
        maxpenalty = max(max(ysmooth)-10,0); %PercentNeutralized/10 should not exceed 10.
        %maxpenalty = max(max(ysmooth)-max(data),0); %PercentNeutralized/10 should not exceed 10.
        %Ensure y and data have same dimensions
        ydiscrete = reshape(ydiscrete,length(ydiscrete),1);
        data = reshape(data,length(data),1);
        C = norm((ydiscrete - data).*w)+10*maxpenalty; %penalty for exceeding maximum NAb level possible


end
