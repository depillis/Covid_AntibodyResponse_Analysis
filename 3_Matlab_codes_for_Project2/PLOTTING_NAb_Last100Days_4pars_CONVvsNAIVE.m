%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING_NAb_Last100Days_CONVvsNAIVE.m  - 4 Parameter Model
%
% Copyright 2025, Lisette de Pillis. All rights reserved.
%
% Author: Lisette de Pillis; depillis@hmc.edu
% Please cite original authorship if any part of this code is used for
% other purposes.
%
% December 3, 2022: plot the last (predicted) 100 days ODE solutions for all subjects
% in one plot using different colors to distinguiseh between convalesced
% and naive.
%%%%%%%%%%%%%%%%%%%%%
%
% November 2, 2022 - Modify "full vs sparse" code to instead compare "convalesced vs naive" so we can see how well the 4-parameter model does on convalesced data.
% December 5, 2022 - in addition, we want to make some "spaghetti plots" to
% compare the predicted additional 90-day trajectories between convalesced
% and not-convalesced.
%
% % If MCMC has been run and a struct file saved, read in that struct file
% % that contatins struct MV_Naive_MCMC holding these fields (but for each subject inside the loop):
% %                                 ID: "94MS"
% %                               days: [4×1 double]
% %                           VaccType: "Pfizer"
% %                             NAbraw: [4×1 double]
% %                 PercentNeutralized: [4×1 double]
% %                                NAb: [4×1 double]
% %                          NAbSCALED: [4×1 double]
% %                          treatDays: [0 28]
% %                      GuessedParams: [0.0210 0.2700 0.0130 0.0340]
% %     PercentOver100SimulationPreFit: [2260×2 double]
% %                          MCMCchain: [5000×4 double]
% %                        MCMCresults: [1×1 struct]
% %                     MCMCprediction: [1×1 struct]
% %                                Sex: "F"
% %                      isConvalesced: "Yes"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ReadPatientData; % Read in NAb data for subjects, creates a "patient" struct
%Save patient struct in "AllPatients" struct
AllPatients = patient;
k1 = 10; % 
k2 = 50; % 

addpath('mcmcstat');
%Plotting colors
% Plot hypothetical immune protection cut-offs
% Define colors
purple = '#7E2F8E';
blue = '#0072BD';
seablue = "#1274D9";
green = '#3E7115';
violet = '#CA53D3';
orange = '#F49309';
yellowdark = "#F4BC0A";
forestgreen = "#46A30F";
redd = '#F4290A';
red  = '#FF0000'; 
pink  = '#FF00FF';
paleblue = '#4DBEEE';
rust = '#A2142F';
colarray = {blue red green redd purple green violet rust paleblue pink blue  orange}; %cell array holding color hex codes

%Define Line Style list
VaccLineStyle={'--',':','-.','-'}; %

% Set size of treatment trigger
uON = 1;
DelTON = 1; % Have treatment on one day at a time;

%% Parameter fitting - select indices of which parameters will be fit
whichpars = 1:4; % Fit only r1, r2, r3, r4 for Model 1
myOPTIONS = odeset('MaxStep',0.01,'NonNegative',1);

lgdfontsize = 22; % Font size of legends
titlefontsize = 24; % Font size of titles
axisfontsize = 22; % Font size of axes

% Function to convert percent on a [0,100] scale to IU/mL
fIU = @(x) 0.00007308*x.^4 - 0.009817268*x.^3 + 0.490278932*x.^2 - 2.881399298*x;

i_hist = 0; %Initiate the histogram counting index
i_histMod = 0;   %Moderna index
i_histPf = 0;    %Pfizer index
i_histConv = 0;  %Convalesced index
i_histNaive = 0; %Not convalesced (naive) index

pnumbers = 1:numel(patient); %All subjects
%pnumbers = 1:5; %Test a subset


% Automate figure numbers - initialize
fignum = 1;

for i = pnumbers
    patienti = AllPatients(i); %holds convalesced status and sex
    % LdeP Oct 31, 2022 Load patient data indiviually
    cd SAVED_RUNS_MCMC/
    PTTitle = sprintf('Patient_%s_MCMC_Fraction_extended_chain_5000',patient(i).ID)
    load(PTTitle);
    close % closes the figure that was just opened
    cd ..
    sprintf('Plotting Subject %s',patient(i).ID)

    % Overwrite the y-label
    ylabel('Fraction Neutralized','FontSize',18,'interpreter','latex');% between 0 and 1



    %%%%%%%%%%
    % LdeP Nov 5, 2022 - 
    % After ReadPatientData, patient(i) struct has fields
    %       ID: "004R"
    %                    days: [18×1 double]
    %                VaccType: "Moderna"
    %                     Sex: [1×1 table]
    %           isConvalesced: [1×1 table]
    %     FractionNeutralized: [18×1 double]
    %                     NAb: [18×1 double]
    %               treatDays: 0
    % But after loading the .mat file that holds the MCMC runs, patient(i)
    % is overwritten with
    %ID: "004R"
    %                         days: [19×1 double]
    %                     VaccType: "Moderna"
    %          FractionNeutralized: [19×1 double]
    %                          NAb: [19×1 double]
    %                    treatDays: 0
    %                GuessedParams: [0.0063 3.5842e-04 8.2520e-04 0.0042 0.9977]
    %     FractionSimulationPreFit: [4250×2 double]
    %                    MCMCchain: [5000×4 double]
    %                  MCMCresults: [1×1 struct]
    %               MCMCprediction: [1×1 struct]
    %                 FitErrorRMSE: 0.0133
    %                 FitCostValue: 1.6559
    % So we need to include the new fields "sex" and "isConvalesced" in
    % patient(i)
    %%%%%%%%%%

    patient(i).Sex = patienti.Sex;
    patient(i).isConvalesced = patienti.isConvalesced;
    
    %Allows for reading in the MODE of monolix-fit parameters - make selection inside InitializeParameters2.m code.
    %Comment out the next two lines if you just want to generate these pics with the "mean" parameter values.
    pt_i = patient(i); %copy part of the struct
    patient(i) = InitParams(i,pt_i); %overwrite GuessedParams using Monolix spreadsheet

    
    PtData  = patient(i).FractionNeutralized;%NAb% scaled to [0,1]
    times   = patient(i).days;

    % Time interval
    tinit = 0;  % Initial time of simulation
    tfinal = patient(i).days(end)+100;% Final time of simulation plus another 100 days (~3 months)
    StepsPerDay=10;
    tspan = linspace(tinit,tfinal,StepsPerDay*tfinal); % Specify output time points for solution

    % Initial conditions for states A and V; These will be the same for M and MS
    A0 = patient(i).GuessedParams(end); %A0 has been fit by Monolix
    V0 = 0;
    Yinit = [A0, V0];

    % Set par vector to hold parameters needed for ODE

    % Initial parameters are already set in
    % patient(i).GuessedParams= [r1, r2, r3, r4, k1, k2];    
    % Set "par" vectors:
    par = [patient(i).GuessedParams(whichpars), k1, k2, uON, patient(i).treatDays'];%LdeP 6 parameters

    r1 = par(1);
    r2 = par(2);
    r3 = par(3);
    r4 = par(4);
    k1 = par(5);
    k2 = par(6);
    
    % Call ODE solver with initial parameter values
    % Note: If we have loaded the patient struct, we may alreay have the
    % field "FractionSimulationPreFit" whose first column is time points,
    % and second column is ODE simulated solutions; But we re-do the
    % solving here since we want to extract the simulated data points,
    % which allows us to determine the RMSE for each fit. Also, the saved
    % ODE solution has been created with the "MEAN" Monolix parameters, and
    % here we are re-doing the fits using the "MODE" Monolix parameters, so
    % the solutiosn will look slightly different. 
    sol = ode15s(@(t,y) NAb_VaccineODE4Param(t,y,par), tspan, Yinit, myOPTIONS);
    
    % Extract solutions for subject
    t = tspan';
    y_A   = deval(sol,t,1)'; % Vector of continuous simulated solution for NAb 
    y_A_V = [y_A deval(sol,t,2)']; % Matrix holding A and V solutions 
    y_A_data = deval(sol,times,1)'; % Vector of simulated NAb solutions just on sample days
    
   %Find error according to cost %April 12 2022
    arg1 = reshape(y_A_data,length(y_A_data),1); %y_A_data holds simulated NAb values in [0,1] fraction on full data collection days
    arg2 = reshape(PtData,length(PtData),1);
    FitError = norm(arg1 - arg2);  % Norm of distance between computed and actual datapoints
    FitErrorRMSE = FitError/length(arg1);
    FitCostValue = cost2(par(whichpars),times,tspan, A0, PtData,par,whichpars); %Cost2 function value distance between data and simulation

    patient(i).FitErrorRMSE = FitErrorRMSE;
    patient(i).FitCostValue = FitCostValue;


    i_hist = i_hist+1; %Increase the index for creating the histogram when computing comparison between M and MS

    % For population level analysis, save all the RMSEs 
    RMSE(i_hist) = FitErrorRMSE;
    % Save RMSE in Moderna or Pfizer vector
    % Update Moderna or Pfizer vector counter
    if patient(i).VaccType=='Moderna'
               i_histMod=i_histMod+1; %Increase Moderna index
               RMSE_Moderna(i_histMod)=FitErrorRMSE;
    else
               i_histPf=i_histPf+1; %Increase Pfizer index
               RMSE_Pfizer(i_histPf)=FitErrorRMSE;
    end
    
    % Update Convalesced or Naive vector counter
    convalescedStatus = patienti.isConvalesced;
    if convalescedStatus =='Yes'
               i_histConv=i_histConv+1; %Increase Convalesced index
               RMSE_isConvalesced(i_histConv)=FitErrorRMSE;
    else
               i_histNaive=i_histNaive+1; %Increase Naive index
               RMSE_isNaive(i_histNaive)=FitErrorRMSE;
    end

% For now, comment out the commands to plot the ODE with the data points 
%%%%%% Begin %%%%%%%
%    fignum = fignum+1;
%    f = figure(fignum); f.Position(3:4) = [1250 600]; % Increase actual window size to capture entire title
%    clf;
%    % To plot continuous simulation result,
%    plot(t,y_A_V(:,1),'k-', 'LineWidth',4,'MarkerSize',3,'DisplayName','Fit to Data Set'); %NAb in [0,1]
%    %plot(patient(i).FractionSimulationPreFit(:,1),patient(i).FractionSimulationPreFit(:,2))
%     
%    % % %    
%     
%     % To plot data points on top of that,
%     hold on;
%     %plot(patient(i).days,patient(i).FractionNeutralized,'b*'); %Not needed
%     %if we are plotting the "errorbars" below
% 
%     %Mark measurement error
%     err = DataErrFraction(PtData); % Function to return error depending on size of vector entry
%     errorbar(times,PtData,err,'bo','MarkerSize',12,'LineWidth',3,'MarkerFaceColor','blue','DisplayName','Data Meas. + Error'); %Error bars
% 
%     axis([0 tspan(end) 0 1.1]); %axis([XMIN XMAX YMIN YMAX])
% 
%     % LdeP Loop to plot treatment days
%     for vind = 1:length(patient(i).treatDays)
%         %colarray{vind} extracts hex code for the "vind" color in colarray cell array
%         VaccNum=sprintf('Vacc%d',vind);
%         % xl=xline(patient(i).treatDays(vind),'LineStyle','-','Color',colarray{vind},'LineWidth',3,'DisplayName',VaccNum);
%         xl=xline(patient(i).treatDays(vind),'Color',colarray{vind},'LineStyle',VaccLineStyle{vind},'LineWidth',5,'DisplayName',VaccNum); %colarray{vind} extracts hex code for the "vind" color in colarray cell array
%     end
%     xlabel('Time in Days','FontSize',axisfontsize,'interpreter','latex');
%     ylabel('Fraction Neutralized','FontSize',axisfontsize,'interpreter','latex');% between 0 and 1
%     TT1 = sprintf('NAb Evolution Using Monolix Parameters \n  r1:%2.6g, r2:%2.6g,  r3:%2.6g, r4:%2.6g, k1:%2.6g, k2:%2.5g\n RMSE: %2.4g',r1, r2, r3, r4, k1,k2,FitErrorRMSE);%LdeP Added May 25, 2022
%%% END %%%

    if convalescedStatus == 'Yes'
        ConvString = 'Convalesced';
    else
        ConvString = 'Naive';
    end

    Scenario = sprintf('Vaccinated %d times - Subject %s, %s \n Vaccine: %s ',length(patient(i).treatDays), patient(i).ID, ConvString, patient(i).VaccType);
    description = [Scenario TT1];
%%% Begin (again)
%     title(description,'FontSize',titlefontsize,'interpreter','latex');
% 
%     %LdeP March 2022 Display Legend
%     lgd = legend;
%         % Add legend properties
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = lgdfontsize;
%     %lgd.Location = 'best';
%     lgd.Location = 'northeast';
%%% END (again)

    %%===================
    %Plot the last 100 days (simulated) 
    f = figure(1);    
    %f.Position(3:4) = [1250 600];
    if convalescedStatus =='Yes'
               solcolor = 'red';
               solinetype = '-';%':' Change line type here
               % Saving the last hundred days NAb data for each convalesced
               % subject in a column
               ConvMatrixLast100Days(:,i)= y_A_V(end-100*StepsPerDay:end,1);
    else
               solcolor = 'blue';
               solinetype = '-';%'-.' Change line type here
               % Saving the last hundred days NAb data for each naive
               % subject in a column
               NaiveMatrixLast100Days(:,i)= y_A_V(end-100*StepsPerDay:end,1);

    end
    plot(t(end-100*StepsPerDay:end),y_A_V(end-100*StepsPerDay:end,1),'Color',solcolor,'LineStyle',solinetype,'LineWidth',5);
    %axis([tspan(end)-100 tspan(end) 0 1.1]); %axis([XMIN XMAX YMIN YMAX])
    hold on; % Hold on until we can plot the last 90 days of the next subject
    %%===================

    %%% For now, we won't save the plots of the ODEs with their data points
    %%% - that can be done with other code. Instead, we want to make
    %%% spagghetti plots. 
    %
    %Save figure

%%% Begin(again)
%     cd FIGURES
%        
%        figure(fignum)
%            PTTitle = sprintf('Patient_%s_MCMC_Fraction_extended_chain_5000',patient(i).ID)
% 
%        filename2 = sprintf('Patient_%s_%sConvalesced', patient(i).ID, convalescedStatus);%NAb sparse vs full ODE plots
%        saveas(gcf,filename2,'fig')
%        saveas(gcf,filename2,'epsc')
%        saveas(gcf,filename2,'jpg')
% 
%     cd ..
%%% END (again)

end

% On the spaghetti plot
axis([0 tspan(end) 0 1.0])
title('Convalesced: Red, Naive: Blue','FontSize',titlefontsize,'interpreter','latex')
xlabel('Last 100 Days','FontSize',axisfontsize,'interpreter','latex'); 
ylabel('Fraction Neutralized','FontSize',axisfontsize,'interpreter','latex')

%%%%%%%%% CONTINUE HERE TO PLOT RMSE BOXCHART %%%%%%%%%%%%%%%%
%RMSE_MtxConvVsNaive = [RMSE_isConvalesced' RMSE_isNaive'];
pause;

fignum = fignum+1;
figure(fignum)
bbfig  = figure(fignum);
bbfig.Position(3:4) = [1250 600];

% Dummy column appended
meanCONV = mean(RMSE_isConvalesced);
meanNAIV = mean(RMSE_isNaive);
medCONV = median(RMSE_isConvalesced);
medNAIV = median(RMSE_isNaive);

% %Plot with dummy mean
% CONVmtx = [RMSE_isConvalesced' meanNAIV*ones(length(RMSE_isConvalesced),1)];
% NAIVmtx = [meanCONV*ones(length(RMSE_isNaive),1) RMSE_isNaive'];
% %Plot with dummy median
CONVmtx = [RMSE_isConvalesced' medNAIV*ones(length(RMSE_isConvalesced),1)];
NAIVmtx = [medCONV*ones(length(RMSE_isNaive),1) RMSE_isNaive'];
clf;
bb=boxchart(CONVmtx,'Notch','on','WhiskerLineStyle','--', 'MarkerSize',10,'MarkerStyle','x','MarkerColor','r');
hold on;
boxchart(NAIVmtx,'Notch','on','WhiskerLineStyle','--', 'MarkerSize',10,'MarkerStyle','x','MarkerColor','r');

%xticklabels(["RMSE FullODE/FullData";"RMSE SparseODE/FullData"]);
%Find the medians and quantiles
medRMSEsCONV = quantile(RMSE_isConvalesced,0.50); % CONV RMSE median
SummaryRMSEsCONV = quantile(RMSE_isConvalesced,[.025 .25 .50 .75 .975]); % a useful summary of x
IQR_RMSEsCONV = iqr(RMSE_isConvalesced);

medRMSEsNAIVE = quantile(RMSE_isNaive,0.50); % NAIVE RMSE median
SummaryRMSEsNAIVE = quantile(RMSE_isNaive,[.025 .25 .50 .75 .975]); % a useful summary of x
IQR_RMSEsNAIVE = iqr(RMSE_isNaive);


set(gca,'XTickLabel',{'RMSE Convalesced','RMSE Naive'},'FontSize',axisfontsize);
%set(gca,'XTickLabel',{fullfulllabel,fullsparselabel},'FontSize',axisfontsize);
% Add text labels
OnPlotFontSize = 18;
Convmed = sprintf('Median: %4.3f', SummaryRMSEsCONV(3))
Naivemed = sprintf('Median: %4.3f', SummaryRMSEsNAIVE(3))
Convquarts = sprintf('Quartiles [%4.3f, %4.3f]',SummaryRMSEsCONV(2),SummaryRMSEsCONV(4));
Naivequarts = sprintf('Quartiles [%4.3f, %4.3f]',SummaryRMSEsNAIVE(2),SummaryRMSEsNAIVE(4));

text(0.1,SummaryRMSEsCONV(3),Convmed,'FontSize',OnPlotFontSize); %Convalesced
text(0.1,SummaryRMSEsCONV(4),Convquarts,'FontSize',OnPlotFontSize);

text(2.3,SummaryRMSEsNAIVE(3),Naivemed,'FontSize',OnPlotFontSize); %Naive
text(2.3,SummaryRMSEsNAIVE(4),Naivequarts,'FontSize',OnPlotFontSize);

ylabel('RMSE');
largestA=max(RMSE_isNaive); largestB=max(RMSE_isConvalesced); largestY=max(largestA,largestB);
ylim([0 1.25*largestY]);
%Overlay with raw data - determine where on the x-axis the data will fall
xdataloc1 = ones(numel(RMSE_isConvalesced),1).*(1+rand(numel(RMSE_isConvalesced),1)+3.5)/5;
xdataloc2 = ones(numel(RMSE_isNaive),1).*(1+rand(numel(RMSE_isNaive),1)+8.5)/5;


%hold on;
scatter(xdataloc1,RMSE_isConvalesced','r','filled')
scatter(xdataloc2,RMSE_isNaive','b','filled')

twolinetitle = sprintf('Comparison of RMSE Distribution: \n Convalesced ODE Fit to Naive ODE Fit')
title(twolinetitle,'Interpreter','Latex')


% filename = 'RMSE_ConvalescedVsNaive_BoxCHARTcomparison';
% cd FIGURES/
%     saveas(gcf,filename,'fig')
%     saveas(gcf,filename,'epsc')
%     saveas(gcf,filename,'jpg')
% cd ..


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



