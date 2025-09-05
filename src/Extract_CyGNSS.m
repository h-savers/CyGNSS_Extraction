%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code to extracts CyGNSS observables by also computing reflectivity
% and Trailing Edge and outputting them as trackwise files. Code developed
% by Emanuele Santi (e.santi@ifac.cnr.it) by reapprising previous
% implementations hosted at Università la Sapienza - Rome (IT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
addpath('src/CyGNSS_Extraction/new_version/functions/')

% ************* Start GUI
conf_file = 'D:\Hamed\HydroGNSS_CalVal\CyGNSS_Extraction\conf\Configuration.mat';

% Set up prompts and dialog config
prompt = {'Taskname: ', ...
          'initdate: ', ...
          'enddate: ', ...
          'savespace: ', ...
          'mainpath: ', ...
          'CyGoutpath: ', ...
          'LatMin: ', ...
          'LatMax: ', ...
          'LonMin: ', ...
          'LonMax: ', ...
          'aggregate_data'};

name = 'Soil moisture L2 processor by Sapienza-CRAS';
numlines = repmat([1 30], 11, 1);
opts.Resize = 'on';
opts.WindowStyle = 'normal';
opts.Interpreter = 'tex';

% Try to load previous inputs if they exist
defaultanswer = repmat({''}, 1, 11);
if isfile(conf_file)
    loaded = load(conf_file, 'Answer');
    if isfield(loaded, 'Answer') && numel(loaded.Answer) >= 11
        defaultanswer = loaded.Answer;
    end
end

% Launch input dialog
Answer = inputdlg(prompt, name, numlines, defaultanswer, opts);

% If user clicked OK, save it
if ~isempty(Answer)
    save(conf_file, 'Answer', '-append');
end

project_name= Answer{1};
mainpath= Answer{5};
CyGoutpath = Answer{6};
verifydir(CyGoutpath)

logpath = fullfile(CyGoutpath, 'log');
if ~exist(logpath, 'dir')
    mkdir(logpath);
end
verifydir(logpath)

CyGfigurepath = fullfile(CyGoutpath, 'figure');
if ~exist(CyGfigurepath, 'dir')
    mkdir(CyGfigurepath);
end
verifydir(CyGfigurepath)


save('D:\Hamed\HydroGNSS_CalVal\CyGNSS_Extraction\conf\Configuration.mat', 'Answer', '-append') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI % ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                  

%%%%%%%%%%%%%%%%%%% DEFINING GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
savespace=Answer{4};                                                            % to apply the CYGNSS land flag before saving (it significantly reduces the size of output file and speeds up the processing
aggregate_data = strcmpi(strtrim(Answer{11}), "true");                                                     % to aggregate data from different days and save it in a single file
%%%%%%%%%%%%%%%%%%%%%%%%%% Geographic limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LatMin=Answer{7};
LatMax=Answer{8};
LonMin=Answer{9};
LonMax=Answer{10};

%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdate=Answer{2};
enddate=Answer{3};
initdatenum=datenum(initdate,'yyyy-mm-dd');
enddatenum=datenum(enddate,'yyyy-mm-dd');
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
    agg_SNR_L1=[];                              % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    agg_PA=[];                                  % peak power
    agg_NF=[];                                  % noise floor
    agg_RXRANGE=[];                             % Rx range [m]
    agg_TXRANGE=[];                             % Tx range [m]
    agg_NST=[];                                 % overall quality
    agg_QC=[];                                  % Quality Flag
    agg_DDM_NBRCS=[];                           % NBRCS
    agg_KURTOSIS=[];                            % Kurtosis
    agg_KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    agg_TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    agg_REFLECTIVITY_LINEAR_L1=[];              % Reflectivity
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
    DoY_infolderpath     = [mainpath, '/' , num2str(doy), '/'];        % Path to DoY folders containing input CyGNSS .nc data


    %%%%%%%%%%%%%%%%%%%%%% CyGNSS data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chkCyGNSSfile=dir([DoY_infolderpath 'cyg0*.ddmi.s' datechar '*.nc']);
    if  ~isempty(chkCyGNSSfile)
        disp('% Extracting CyGNSS data ...')
        [DoY,SoD,spacecraft_num,pseudo_random_noise,SPLAT,SPLON,THETA,EIRP,SNR_L1,PHI_Initial_sp_az_orbit, ...
            REFLECTIVITY_LINEAR_L1,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,PA,QC,noise_floor,BRCS]= ...
            extract_CyGNSS(nsat,datechar,doy,DoY_infolderpath,logpath,lambda,Doppler_bins,savespace,delay_vector,Power_threshold);            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aggregate_data
            disp(['% cat variables from day ', datechar ' to aggregated output file']);
            agg_DoY=cat(1,agg_DoY,DoY(:));
            agg_SoD=cat(1,agg_SoD,SoD(:));
            agg_SCID=cat(1,agg_SCID,spacecraft_num(:));
            agg_PRN=cat(1,agg_PRN, pseudo_random_noise(:));
            agg_SPLAT=cat(1,agg_SPLAT, SPLAT(:));
            agg_SPLON=cat(1,agg_SPLON, SPLON(:));
            agg_THETA=cat(1,agg_THETA, THETA(:));
            agg_EIRP=cat(1,agg_EIRP, EIRP(:));
            agg_SNR_L1=cat(1,agg_SNR_L1, SNR_L1(:));
            agg_PHI_Initial_sp_az_orbit=cat(1,agg_PHI_Initial_sp_az_orbit, PHI_Initial_sp_az_orbit(:));
            agg_REFLECTIVITY_LINEAR_L1=cat(1,agg_REFLECTIVITY_LINEAR_L1,REFLECTIVITY_LINEAR_L1(:));
            agg_KURTOSIS=cat(1,agg_KURTOSIS, KURTOSIS(:));
            agg_KURTOSIS_DOPP_0=cat(1,agg_KURTOSIS_DOPP_0, KURTOSIS_DOPP_0(:)); 
            agg_TE_WIDTH=cat(1,agg_TE_WIDTH, TE_WIDTH(:)); 
            agg_DDM_NBRCS=cat(1,agg_DDM_NBRCS, DDM_NBRCS(:)); 
            agg_PA=cat(1,agg_PA, PA(:));
            agg_QC=cat(1,agg_QC, QC(:)); 
            agg_NF=cat(1,agg_NF, noise_floor(:));
            agg_BRCS=cat(3, agg_BRCS, BRCS);                                      
            % agg_RXRANGE=cat(1,agg_RXRANGE,RXRANGE); % these variables are extracted in extract_CyGNSS function, but then they are not passed to the function output. Ask Hamed why
            % agg_TXRANGE=cat(1,agg_TXRANGE,TXRANGE);
            % agg_NST=cat(1,agg_NST,NST);
        else
            disp('% saving CyGNSS data');
            save([CyGoutpath, '/', datechar '.mat'], 'Year', 'DoY', 'SoD', 'spacecraft_num', ...  
                'pseudo_random_noise', 'SPLAT', 'SPLON', 'THETA', 'EIRP', 'SNR_L1', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR_L1', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'noise_floor', '-v7.3')
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
    spacecraft_num=agg_SCID;
    pseudo_random_noise=agg_PRN;
    SPLAT=agg_SPLAT;
    SPLON=agg_SPLON;
    THETA=agg_THETA;
    EIRP=agg_EIRP;
    SNR_L1=agg_SNR_L1;
    PHI_Initial_sp_az_orbit=agg_PHI_Initial_sp_az_orbit;
    REFLECTIVITY_LINEAR_L1=agg_REFLECTIVITY_LINEAR_L1;
    KURTOSIS=agg_KURTOSIS;
    KURTOSIS_DOPP_0=agg_KURTOSIS_DOPP_0; 
    TE_WIDTH=agg_TE_WIDTH; 
    DDM_NBRCS=agg_DDM_NBRCS; 
    PA=agg_PA;
    QC=agg_QC; 
    noise_floor=agg_NF;
    BRCS=agg_BRCS;   

    % Saving aggregated data
    save([CyGoutpath, '/', 'aggregated_' daterangechar '.mat'], 'Year', 'DoY', 'SoD', 'spacecraft_num', ...  
                'pseudo_random_noise', 'SPLAT', 'SPLON', 'THETA', 'EIRP', 'SNR_L1', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR_L1', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'noise_floor', '-v7.3');
end
s=duration(0,0,toc);
close all
disp(['total duration is ' char(duration(0,0,toc))])