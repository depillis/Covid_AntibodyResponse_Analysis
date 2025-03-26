function [patientstruct] = InitParams(ptnumber,ptstruct)
%InitParams(ptnumber) initializes the r1,...,r4 and A0 parameters for
%patient ptnumber only.
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
%   
SubjectParams_4ParModelMean = readmatrix('Monolix_estimatedIndividualParameters_2022_2.xlsx','Sheet','Monolix_Mean','Range','A2:I101'); % convalesced and naive, params scaled to [0,1] %NAb
SubjectParams_4ParModelMode = readmatrix('Monolix_estimatedIndividualParameters_2022_2.xlsx','Sheet','Monolix_Mode','Range','A2:I101'); % convalesced and naive, params scaled to [0,1] %NAb

i = ptnumber;

%Subj_params = SubjectParams_4ParModelMean(i,2:6);
Subj_params = SubjectParams_4ParModelMode(i,2:6);
ptstruct.GuessedParams=Subj_params; % Holds r1, r2, r3, r4, A0

if ptstruct.days(1) > 0 % then we need to prepend a 0 day for vaccination
%         % Update for Convalesced - these individuals likely have some
%         % positive antibody levels prior to vaccination. If so, it does not
%         % make sense to set the first day level to zero. Instead, we'll
%         % assign the A0 determined by Monolix
         A0=ptstruct.GuessedParams(end);
         patient(i).days = [0; patient(i).days]; %Start simulation at day 0 (if not provided)
         patient(i).FractionNeutralized = [A0; patient(i).FractionNeutralized]; % put NAb determined by Monolix into day 1 (if NAb not known)        
     
         fN = A0; % For cleaner computation below
    
         % IU_per_mL = 0.00007308*pN.^4 - 0.009817268*pN.^3 + 0.490278932*pN.^2 - 2.881399298*pN;
         fIU = inline('0.00007308*x.^4 - 0.009817268*x.^3 + 0.490278932*x.^2 - 2.881399298*x','x');
         IU_per_mL = fIU(fN*100); % fIU converts PERCENT neutralized to IU per mL   
    
    % Some NAb IU/mL level entries are listed as '<4' in the spreadsheet, 
    % and show up as NaN, so instead of using the spreadsheet-computed NAB in IU/mL,
    % convert NAB Flow % neutralized directly to IU/mL with
    % quartic formula Aditxt uses.
         ptstruct.NAb = [IU_per_mL; ptstruct.NAb] % Store this computed conversion in patient.NAb
    % Get rid of negative values that appear because of the quartic - set those to zero
         ptstruct.NAb((ptstruct.NAb<0))=0;
end

k1 = 10; % After Monolix fit
k2 = 50; % After Monolix fit

patientstruct = ptstruct;
end