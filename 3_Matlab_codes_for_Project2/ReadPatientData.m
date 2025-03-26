% ReadPatientData from spreadsheet
% May 25, 2022 - Updated to read spreadsheets for scaled [0,1] NAb % data.
% October 24, 2022 - Built from ReadPatientData2.m, but modified to read in
% file nonlinearMixedModelDataSet_modified.xlsx which does not have actual
% dates, but only "Days Since First Dose."
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
% Reads in patient data from spreadsheet nonlinearMixedModelDataSet_2022.xlsx
% Returns a "patient" struct.
%
%
%
% "patient" struct holds:
%     ID        - patient ID number
%     days      - days on which NAb levels are measured
%     VaccType  - string with vaccine type that was given
%     NAbraw    - IU/mL NAb from spreadsheet
%     PercentNeutralized - percent neutralization from Aditxt flow assay
%     NAb       - IU/mL NAb computed directly from Aditxt formula
%     NAbSCALED - 100 IU/mL NAb
%     treatDays - vector of days on which vaccination is administered
%     Sex       - sex of the subject - "M", "F", "NA"
%     isConvalesced - convalesced status - "Yes" or "No"
%
%
% In the spreadsheet, the current script assumes that these are the columns:
% A(1): "subject" Subject number (this number must be read in as a string since the
% format is either xxxM or xxxR)
% B(2): "dsfd" Days since first dose. Here, all subjects start on Day 0,
% and that is also the day they receive their first dose.
% C(3): "nabFlow" NAb level in [0,1] scale, or "NA" if the row just has
% information about when a dose was administered but no NAb level recorded.
% D(4): "isConvalesced" Whether subject is convalesced or not. Entry is
% "Yes" or "No"
% E(5): "sex" Whether subject is Male (M), Female (F) or Unknown (NA).
% F(6): "dose2Type" will be Moderna or Pfizer
% G(7): "dose2Type" will be Moderna or Pfizer (for now it matches Dose1)
% H(8): "eventid" This is only used by Monolix to indicate whether the data
% on this row represent a dose event (1) or not (0).
% I(9): "amount" This is only used by Monolix to indicate the amount of
% dose (1 if eventid is 1, and NA if eventid is 0).
% J(10): "a0" This is the initial value of the NAb levels when Dose 1 is
% administered.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Colors for plotting
% Define colors
        pink  = '#FF00FF';
        paleblue = '#4DBEEE';
        rust = '#A2142F';
        purple = '#7E2F8E';
        blue = '#0072BD';
        green = '#3E7115';
        violet = '#CA53D3';
        orange = '#F49309';
        colarray = {blue purple green violet rust paleblue pink purple orange}; %cell array holding color hex codes

% LdeP updated to read convalesced/naive spreadsheet 
PtData = readmatrix('nonlinearMixedModelDataSet_2022.xlsx','Sheet','ConvNaive_PrepForMatlabFormat','Range','A2:J514'); % Spreadsheet from Mark

% Read in table format so we can access the text entries in the spreadsheet
% At some point, maybe figure out how to use just this, but for now,
% readmatrix is easier to use with numerics.
PtTable = readtable('nonlinearMixedModelDataSet_2022.xlsx','Sheet','ConvNaive_PrepForMatlabFormat','Range','A2:J514'); % Spreadsheet from Mark


% Specify columns for this spreadsheet:
% Subject number column number
SN = 1;
% Collection date column number - LdeP We have no collection date data.
CD = nan;
% Convalesced status (convalesced "Yes", naive "No")
isConv = 4;
%Sex of Subject
SX = 5;
% Vaccine Type
VaccType = 6;
% First vaccine dose date column number - LdeP We have no date data
V1 = nan;
% Days since first vaccine dose  column number
DSFD = 2;
% Second vaccine dose date column number - LdeP We have no date data
V2 = nan;% 
% Days since second vaccine dose column number - We don't have a second
% dose column, we have to figure it out from DSFD column combined with
% whether that row has eventid ==1.
DSSD = nan;
% Third vaccine dose date column number - this data set doesn't have Dose 3
V3 = nan;
% Days since third vaccine dose column number - no dose 3 in this dataset
DSTD = nan; 
% Flow NAb fraction neutralized column number - WARNING: this was percent
% earlier, so make sure not to scale as we did before
FN = 3;
% IUmL NAb column number - we don't have this.
IUN = nan;
% Initial Condition - this is NEW for this spreadsheet
IC = 10; % The column number of the IC is 10.
% With Monolix-format spreadsheet, we have an eventid column
EventID = 8; %Where the value in this column==1, a dose was administered.

% Extract patient ID numbers
PtIDVector = string(PtTable{:,SN});
PtID = unique(PtIDVector,'rows'); %Convert cell array to string array, return unique rows

% % Remove NaNs
% % PtID(isnan(PtID))=[];

%LdeP Initialize how to count rows to identify subject j's vaccine type
ptj_VaccType_row = 0; %start at 0

% Partition into individual patient matrices using a struct
for j = 1:length(PtID) % Run through each patient number
    % Initialize treatDays vector:
    treatDays = [];
    Pindex = find(PtIDVector==PtID(j)); % Find all rows associated with Patient j
     
    % Build struct "patient"
    patient(j).ID = PtID(j); % Store patient ID as a string array
    patient(j).days = PtData(Pindex,DSFD);    % Column DSFD is number of days since first vaccine. We assume first vaccine is always on day 0.
    
    %LdeP - BEFORE removing NAN days, so that we count row
    %entries correctly, we'll use the Table form of the data to store the
    %name of the vaccine type.
    %LdeP -  Save vaccine type in patient struct
    %PtTable holds the vaccine type in column "VaccType
    
    %LdeP - this method of finding the vaccine type does not
    %work with re-ordered spreadsheets - instead use one of the rows in
    %Pindex, as long as vaccine type is listed in *each row* of the subject
    %entries.
    ptj_VaccType_row = Pindex(1); %First row in which subject "j" entries are
    
    
    VaccTypeStr = PtTable(ptj_VaccType_row,VaccType); % Read vaccine type from the last row of patient(j)'s entries, and column 7
    patient(j).VaccType = string(VaccTypeStr{1,1}); % This accesses vaccine type as a cell array, and converts to a string

%%% LdeP additional fields
    patient(j).Sex = string(table2array(PtTable(ptj_VaccType_row,SX))); % sex is "M", "F", or "NA"
    patient(j).isConvalesced = string(table2array(PtTable(ptj_VaccType_row,isConv))); % convalesced status is "Yes" or "No"

%%% LdeP Just make sure there are no empty DSD1 (days since dose 1) cells in spreadsheet
   
    
    %The original spreadsheet includes rows in which FractionNeutralized
    %is entered as NaN. In these "NaN" rows we should look for the days
    %since first dose to extract on which days a subject was vaccinated.
    %ALTERNATELY: We can look at column "eventid" and extract the day
    %numbers where that value == 1.
    patient(j).FractionNeutralized = PtData(Pindex,FN); % Fraction neutralized in [0,1] scale
    
    % Extract dose days indexes -- choose one of the two methods below.
    % %LdeP METHOD 1:  Extract dose days from fraction neutralized column
    % ptj_DoseDaysIndx = find(isnan(patient(j).FractionNeutralized)==1); % Index of entries when doses were administered
    %LdeP METHOD 2:  Extract dose days from eventid column
    ptj_DoseDaysIndx = find(PtData(Pindex,EventID)); % Index of entries when doses were administered
    
    DoseDays = patient(j).days(ptj_DoseDaysIndx); % Extract day numbers (dsfd) on which vaccines were given    

   %%% LdeP Just make sure there are no empty NAb entry cells in spreadsheet
    % The "FractionNeutralized" vector will contain "NaN"s on the dose
    % days, so we should remove the rows with NaN before continuing.
    % Remove rows in both the .days and .FractionNeutralized entries of the
    % patient struct.
    patient(j).days(ptj_DoseDaysIndx)=[];
    patient(j).FractionNeutralized(ptj_DoseDaysIndx)=[];
    
    fN = patient(j).FractionNeutralized; % For cleaner computation below
    
    % IU_per_mL = 0.00007308*pN.^4 - 0.009817268*pN.^3 + 0.490278932*pN.^2 - 2.881399298*pN;
    fIU = inline('0.00007308*x.^4 - 0.009817268*x.^3 + 0.490278932*x.^2 - 2.881399298*x','x');
    IU_per_mL = fIU(fN*100); % fIU converts PERCENT neutralized to IU per mL   
    
    % Some NAb IU/mL level entries are listed as '<4' in the spreadsheet, 
    % and show up as NaN, so instead of using the spreadsheet-computed NAB in IU/mL,
    % convert NAB Flow % neutralized directly to IU/mL with
    % quartic formula Aditxt uses.

    patient(j).NAb = IU_per_mL; % Store this computed conversion in patient.NAb
    % Get rid of negative values that appear because of the quartic - set those to zero
    patient(j).NAb((patient(j).NAb<0))=0;
    
    
    patient(j).treatDays = DoseDays; %save in patient struct 
       
end 
