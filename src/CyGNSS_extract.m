function CyGNSS_extract(configurationPath)
close all
clearvars -except configurationPath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code to extracts CyGNSS observables by also computing reflectivity
% and Trailing Edge and outputting them as trackwise files. Code developed
% by Emanuele Santi (e.santi@ifac.cnr.it) by reapprising previous
% implementations hosted at Università la Sapienza - Rome (IT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ex=exist('configurationPath') ;
if ex ==0
    mode="GUI" ;
    
    [configurationfile configurationPath] = uigetfile('../*.cfg', 'Select input configuration file') ; 
    configurationPath= [ configurationPath configurationfile]  ; 
else
    if ~isfile(configurationPath)
        throw(MException('INPUT:ERROR', "Cannot find configuration file. Please check the command line and try again."))
    end
    mode="input" ;
end

[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax] = ReadConfFile(configurationPath);

switch mode
    case "GUI" 

% ************* Start GUI

    disp('GUI mode')
% 
Answer{1}= char(Taskname) ;         Answer{2}=char(string(initdate)) ;
Answer{3}=char(enddate)  ;          Answer{4}= char(string(savespace)) ;
Answer{5}=char(CyGinpath)  ;         Answer{6}= char(string(CyGoutpath)) ;
Answer{7}= char(string(logpath)) ;
Answer{8}=char(string(LatMin))  ;   Answer{9}= char(string(LatMax)) ;
Answer{10}=char(string(LonMin))  ;   Answer{11}=char(string(LonMax)) ;
% Answer{6}=char(string(LatSouth)); Answer{12}= char(DDM) ;

% Set up prompts and dialog config
prompt = {'Outfileprefix: ', ...
          'First day to extract: ', ...
          'Last day to extract: ', ...
          'Savespace? [Yes/No]: ', ...
          'DataInputRootPath: ', ...
          'DataOutputRootPath: ', ...
          'LogsOutputRootPath: ', ...
          'Southernmost latitude: ', ...
          'Northernmost latitude: ', ...
          'Westernmost longitude: ', ...
          'Easternmost longitude: '};
      
name = 'Extraction og CyGNSS L1b data';
numlines = repmat([1 90], 11, 1);
opts.Resize = 'on';
opts.WindowStyle = 'normal';
opts.Interpreter = 'tex';

defaultanswer={Answer{1},Answer{2},...
                 Answer{3},Answer{4},Answer{5},Answer{6},Answer{7},...
                 Answer{8},Answer{9},Answer{10},Answer{11}} ; 

% Launch input dialog
Answer = inputdlg(prompt, name, numlines, defaultanswer, opts);

project_name= Answer{1};

Taskname= Answer{1};
initdate= Answer{2};
enddate= Answer{3};
savespace= Answer{4};
CyGinpath= Answer{5};
CyGoutpath=Answer{6};
logpath=Answer{7};
LatMin=Answer{8};
LatMax=Answer{9};
LonMin=Answer{10};
LonMax=Answer{11};


init_SM_Day=datetime(Answer{2}) ; 
final_SM_Day=datetime(Answer{3}) ; 

% init_SM_Day=datetime(str2num(Answer{5}), str2num(Answer{6}), str2num(Answer{7})) ;
% final_SM_Day=datetime(str2num(Answer{8}), str2num(Answer{9}), str2num(Answer{10})) ;

%scrivere il configuration
WriteConfig(configurationPath, Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI % ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "input" 
    disp('input mode')

[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax] = ReadConfFile(configurationPath);

% end
end

%%%%%%%%%%%%%%%%%%%%%%%%% DEFINING GENERAL PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CyGinpath='/Users/pgmadonia/PhD/coding/HydroGNSSCalVal/test_data/V3.2/2021/';                            % Input CyGNSS DoY folders
%CyGoutpath='/Users/pgmadonia/PhD/coding/HydroGNSSCalVal/test_data/extracted/';                            % Path to save output matfile containing daily trackwise CYGNSS data (all 8 satellites)
%CyGfigurepath='/Users/pgmadonia/PhD/coding/HydroGNSSCalVal/test_data/extracted/figures/';                 % Path to save output figures for verification
%logpath='/Users/pgmadonia/PhD/coding/HydroGNSSCalVal/test_data/extracted/logs';                          % Error log path 

%%%%%%%%%%%%%%%%%% VERIFYING OUTPUT DIRECTORIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
CyGfigurepath=CyGoutpath ; 
verifydir(CyGoutpath)
verifydir(CyGfigurepath) % ????????  CyGfigurepath NEVER USED !!!!
verifydir(logpath)

%%%%%%%%%%%%%%%%%%%%%% AGB, LCC DEM, and SLOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % auxfile='auxiliary.mat';                                                   

%%%%%%%%%%%%%%%%%%% DEFINING GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Defining the steps of processing %%%%%%%%%%%%%%%%%%%%%%%%
savespace='yes';                                                            % to apply the CYGNSS land flag before saving (it significantly reduces the size of output file and speeds up the processing
aggregate_data = false;                                                     % to aggregate data from different days and save it in a single file

%%%%%%%%%%%%%%%%%%%%%%%%%% Geographic limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LatMin=-90;
LatMax=90;
LonMin=-180;
LonMax=180;

%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdate='02/10/2021';
enddate='03/10/2021';
initdatenum=datenum(initdate,'dd/mm/yyyy');
enddatenum=datenum(enddate,'dd/mm/yyyy');
datelist=initdatenum:enddatenum;

%%%%%%%%%%%%%%%%%%% Parameters for CyGNSS extraction %%%%%%%%%%%%%%%%%%%%%%
CA_chip_delay = 0.2552;                                                    % around 1/4 of CA code chip
delay_vector = 0:CA_chip_delay:16*CA_chip_delay;
Doppler_bins=1:1:17; 
Power_threshold=0.7;
lambda=0.1903;                                                             % 0.19 m --> 19 cm  
nsat=8;                                                                    % CyGNSS constellation
resolution=36;                                                             % Km for ease grid converter

%%%%%%%%%%%%%%%%%%%%% INITIALIZING EMPTY VARIABLES FOR AGGREGATED SINGLE OUTPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if aggregate_data
    daterangechar = [datestr(initdatenum,'yyyymmdd') '-' datestr(enddatenum,'yyyymmdd')]; % date range for aggregated output file
    disp(['% Processing data from ' daterangechar ' and saving in a single output file'])

    agg_SCID=[];                                % CYGNSS sat ID
    agg_SoD=[];                                 % second of the day
    agg_DoY=[];                                 % day of the year
    agg_PRN=[];                                 % PRN --> Prn code = prn -->trasmettitore (1 10 22 etc..)
    agg_SPLAT=[];                               % SP lat on ground
    agg_SPLON=[];                               % SP lon on ground
    agg_THETA=[];                               % incidence angle
    agg_PHI_Initial_sp_az_orbit=[];             % azimuth angle in specular point orbit frame
%    agg_GAIN=[];                                % gain of receiver antenna [dBi]
    agg_EIRP=[];                                % EIRP [W]
    agg_SNR=[];                                 % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    agg_PA=[];                                  % peak power
    agg_NF=[];                                  % noise floor
    agg_RXRANGE=[];                             % Rx range [m]
    agg_TXRANGE=[];                             % Tx range [m]
    agg_NST=[];                                 % overall quality
    agg_LF=[];                                  % land flag
    agg_QC=[];                                  % Quality Flag
    agg_DDM_NBRCS=[];                           % NBRCS
    agg_KURTOSIS=[];                            % Kurtosis
    agg_KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    agg_TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    agg_REFLECTIVITY_LINEAR=[];                 % Reflectivity
    agg_BRCS=[];                                % added by Hamed to save full ddm
else
    disp('% Processing each day separately and saving individual output files');
end

%%%%%%%%%%%%%%%%%%%%% STARTING THE MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% % figure;

for ii=1:length(datelist)     % loop on all the days 
    datechar=datestr(datelist(ii),'yyyymmdd');
    Year=datechar(1:4);
    d=datetime(datechar,'InputFormat','yyyyMMdd');
    doy=day(d,'dayofyear');
    disp(['% now processing day ' num2str(ii) ' out of ' num2str(length(datelist)) ', DoY: ' num2str(doy) ', date: ' datechar]);
    
    %%%%%%%%%%%%%%%%%%%%%% Defining paths for each day %%%%%%%%%%%%%%%%%%%%%%%%%
    DoY_infolderpath     = [CyGinpath    , num2str(doy), '/'];        % Path to DoY folders containing input CyGNSS .nc data


    %%%%%%%%%%%%%%%%%%%%%% CyGNSS data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chkCyGNSSfile=dir([DoY_infolderpath 'cyg0*.ddmi.s' datechar '*.nc']);
    if  ~isempty(chkCyGNSSfile)
        disp('% Extracting CyGNSS data ...')
        [DoY,SoD,SCID,PRN,SPLAT,SPLON,THETA,EIRP,SNR,PHI_Initial_sp_az_orbit, ...
            REFLECTIVITY_LINEAR,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,PA,QC,NF,LF,BRCS]= ...
            extract_CyGNSS(nsat,datechar,doy,DoY_infolderpath,logpath,lambda,Doppler_bins,savespace,delay_vector,Power_threshold);            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aggregate_data
            disp(['% cat variables from day ', datechar ' to aggregated output file']);
            agg_DoY=cat(1,agg_DoY,DoY(:));
            agg_SoD=cat(1,agg_SoD,SoD(:));
            agg_SCID=cat(1,agg_SCID,SCID(:));
            agg_PRN=cat(1,agg_PRN, PRN(:));
            agg_SPLAT=cat(1,agg_SPLAT, SPLAT(:));
            agg_SPLON=cat(1,agg_SPLON, SPLON(:));
            agg_THETA=cat(1,agg_THETA, THETA(:));
            agg_EIRP=cat(1,agg_EIRP, EIRP(:));
            agg_SNR=cat(1,agg_SNR, SNR(:));
            agg_PHI_Initial_sp_az_orbit=cat(1,agg_PHI_Initial_sp_az_orbit, PHI_Initial_sp_az_orbit(:));
            agg_REFLECTIVITY_LINEAR=cat(1,agg_REFLECTIVITY_LINEAR,REFLECTIVITY_LINEAR(:));
            agg_KURTOSIS=cat(1,agg_KURTOSIS, KURTOSIS(:));
            agg_KURTOSIS_DOPP_0=cat(1,agg_KURTOSIS_DOPP_0, KURTOSIS_DOPP_0(:)); 
            agg_TE_WIDTH=cat(1,agg_TE_WIDTH, TE_WIDTH(:)); 
            agg_DDM_NBRCS=cat(1,agg_DDM_NBRCS, DDM_NBRCS(:)); 
            agg_PA=cat(1,agg_PA, PA(:));
            agg_QC=cat(1,agg_QC, QC(:)); 
            agg_NF=cat(1,agg_NF, NF(:));
            agg_LF=cat(1,agg_LF,LF(:));
            agg_BRCS=cat(3, agg_BRCS, BRCS);                                      
            % agg_RXRANGE=cat(1,agg_RXRANGE,RXRANGE); % these variables are extracted in extract_CyGNSS function, but then they are not passed to the function output. Ask Hamed why
            % agg_TXRANGE=cat(1,agg_TXRANGE,TXRANGE);
            % agg_NST=cat(1,agg_NST,NST);
        else
            disp('% saving CyGNSS data');
            save([CyGoutpath datechar '_2.mat'], 'Year', 'DoY', 'SoD', 'SCID', ...  
                'PRN', 'SPLAT', 'SPLON', 'THETA', 'EIRP', 'SNR', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'NF','LF', '-v7.3')
        end
    %%%%%%%%%%%%%%%%%%%%% Displaying Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %          scattermap(real(10.*log10(REFLECTIVITY_LINEAR)),SPLAT,SPLON,datechar,-40,0)
% %          print(gcf,[CyGfigurepath datechar '.png'],'-dpng','-r300')   
    else
        disp('% CyGNSS data files missing for the selected date, output file not saved .... ')
    end      
end
%%%%%%%%%%%%%%%%%%%%%%%%% END OF THE MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if aggregate_data
    disp('% Saving aggregated data in a single output file')
    
    % Renaming variables
    DoY=agg_DoY;
    SoD=agg_SoD;
    SCID=agg_SCID;
    PRN=agg_PRN;
    SPLAT=agg_SPLAT;
    SPLON=agg_SPLON;
    THETA=agg_THETA;
    EIRP=agg_EIRP;
    SNR=agg_SNR;
    PHI_Initial_sp_az_orbit=agg_PHI_Initial_sp_az_orbit;
    REFLECTIVITY_LINEAR=agg_REFLECTIVITY_LINEAR;
    KURTOSIS=agg_KURTOSIS;
    KURTOSIS_DOPP_0=agg_KURTOSIS_DOPP_0; 
    TE_WIDTH=agg_TE_WIDTH; 
    DDM_NBRCS=agg_DDM_NBRCS; 
    PA=agg_PA;
    QC=agg_QC; 
    NF=agg_NF;
    LF=agg_LF;
    BRCS=agg_BRCS;   

    % Saving aggregated data
    save([CyGoutpath 'aggregated_' daterangechar '_2.mat'], 'Year', 'DoY', 'SoD', 'SCID', ...  
                'PRN', 'SPLAT', 'SPLON', 'THETA', 'EIRP', 'SNR', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'NF','LF', '-v7.3');
end
s=duration(0,0,toc);
close all
disp(['total duration is ' char(duration(0,0,toc))])
end