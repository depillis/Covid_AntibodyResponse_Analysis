%% Error function defined here
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
% January 4, 2022 update: Create new "error bar" ranges formula according
%         to the Aditxt "rule of thumb": 
%         Above 700 IU/mL use 5% error
%         150 IU/mL to 700 IU/mL use 10% error;
%         90 IU/mL to 150 IU/mL use 15% error;
%         Below 90 IU/mL use 20% error.

function E = DataErrFraction(data)
% Input: data - a vector of measured data in fraction
%        DataErrFraction(patient(i).FractionNeutralized)
% Output: E - a vector holding the amount of error in each entry

ysize = length(data); % loop through each point
for i = 1:ysize
    if data(i) < (90/2200) % 90 IU/mL
        E(i) = .2*data(i);
    elseif data(i) < (150/2200)  % 150 IU/mL
        E(i) = 0.15*data(i);
    elseif data(i) < (700/2200)  % 700 IU/mL
        E(i) = 0.10*data(i);
    else
        E(i) = 0.05*data(i);
    end
end

E = E'; % Return a column vector

end
