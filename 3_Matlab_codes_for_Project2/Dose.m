%% New Dosing (any number of pulses) function defined here
%
% Project: Modeling Neutralizing Antibody Response to mRNA Covid vaccine
% File: Dose.m
% Called by ResistanceEvolution.m
% Date: 2022-02-17
% Author: LdeP
% Contributors: Lisette de Pillis
% Updates: Feb 16, 2022 LdeP
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
%
% function u = Dose(t_all, treatDays, uON, DelTON);
% t_all: The time points at which the model results are evaluated (can be a single point)
% treatDays:   The time points (days) on which a dose is administered; should be a subset of t_all
% uON:  The strength (amplitude) of the treatment dose (usually just 1)
% How long the treatment is given: we assume 1 [day].
%
%
% Function u can be used as a treatment function.
% u(t) will take on either the value 0 (default) or uON, the amount of
% treatment to deliver when turned on.
%
% Inputs:
%   t_all - current time vector over which to evaluate u(t).
%       t_all is meant to be the time vector for the entire simulation if
%       it is a vector.
%       t_all can be a single value extracted from one point in the
%       simulation.
%   treatDays:   The time points (days) on which a dose is administered; should be a subset of t_all
%   uON - value of u when non-zero
% Output:
%   Vector value u(t) for time vector t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_all:        a vector of time values (can be non-integer)
% treatDays:    vector of days on which a dose is given; (integer vector) should be a subset of t_all
% uON:          strength (amplitude) of treatment (scalar)
function u = Dose(t_all, treatDays, uON)
    tint = floor(t_all); % extract whole day numbers (no fractional days)
    u = zeros(size(t_all)); % Default off value of u(t)
   
    for i = 1:length(treatDays) %LdeP May 10, 2022 Update to make this faster
%         treatIndx = find(tint==treatDays(i)); % Save vector of indexes on which treatment i is administered
%         u(treatIndx) = uON;
        u(tint==treatDays(i))=1;
    end
      
end %end function
%%%%%%%%%%%%%%%%%%%%%%%%%%%
