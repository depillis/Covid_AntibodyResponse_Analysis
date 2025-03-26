%% CONCAVE DOWN SATURATED (CD) ODE function defined here
% NAb_Vaccine_ODE.m : ODE system describing Neutralizing Antibody response
% to vaccine doses over time.
% 4 parameter version
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
% Author: Lisette de Pillis
% Date:   November 24, 2021
% Update: Modifying to use new dosing function Dose.m. Have to find a way
% to convert current input arguments into a vector of days on which
% vaccination is given.
%
%
% Function input:
% t = time points of evaluation.
% y = initial conditions for solution: A = y(1), V=y(2)
% ParamVec = [r1, r2, r3, r4, k1, k2, uON, DelTON, DelTOFF, tStart, tEND]
%
% Function output:
% dydt  Vector holding state derivatives dA/dt and dV/dt
% u(t)  Vector holding treatment u(t)
%
% Calls treatment function Dose.m:
%               u = Dose(t, treatDays, uON);
%
% And defines an ODE where 
% Populations/States:
% A(t):    NAbs:            [100 IU/mL]
% V(t):    Vaccine in body: [mL]

% Forcing function:
% u(t):    Vaccine dose     [mL/day] 

% New parameters
% LdeP Read in model parameters
% LdeP Parameters for NAb-Dose response model
% r1      ; % [1/day] Intrinsic growth response of NAb to Vaccine dose
% r2      ; % [1/day] Mass action parameter r2*(A^2)*V
% r3      ; % [1/day] Intrinsic growth for logistic term
% r4      ; % [mL/(100 IU)] 1/(NAb carrying capacity)
% k1      ; % [1/day] Michaelis Menten scaling for Vaccine clearance
% k2      ; % [IU/mL] Michaelis Menten half-saturation in Vaccine clearance term
% uON     ; % [mL/day] value of injection u when turned on
% treatDays; %[day] vector - vector of days on which vaccine is administered
% tStart  ; % [day]    start of treatment time for u(t)
% tEnd    ; % [day]    end of treatment time for u(t) 
% Parameter vector created:
% ParamVec = [r1, r2, r3, r4, k1, k2, uON, DelTON, DelTOFF; tStart; tEnd; tinit; tfinal];

function [dydt, u] = NAb_VaccineODE(t, y, ParamVec)

    % Extract odel parameters
    r1 = ParamVec(1);
    r2 = ParamVec(2);
    r3 = ParamVec(3); 
    r4 = ParamVec(4); 
    
    k1 = ParamVec(5); % Treatment uptake and elimination parameters
    k2 = ParamVec(6); % Treatment uptake and elimination parameters

    % Extract treatment parameters
    uON = ParamVec(7);
    % LdeP treatDays may be a vector
    treatDays = ParamVec(8:end);
   

    % Extract States by name
    A = y(1); %NAb in 100 IU/mL
    V = y(2); %Vaccine in body in mL

%     % Intialize treatment value
%     u=0; %Vaccine dose ini mL/day

    %LdeP assume we don't need a start and end treatment period, just
    %treatDays vector
    u = Dose(t, treatDays, uON);  


    % Antigen population [100 IU/mL]
     dAdt = r1*V + r2*V.*A + A.*(r3 - r4*A); % original 4-parameter ODE
      
    % Vaccine concentration [mL]
    % Injection amount is u, and vaccine gets cleared with k1*V/(k2+V).
    % 
    p=1; %
    dVdt = p*u - k1*V./(k2+V); % Using uON=1.


    dydt = [dAdt; dVdt];

end

