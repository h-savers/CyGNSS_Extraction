function CyGNSS_extract(configurationPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code to extracts CyGNSS observables by also computing reflectivity
% and Trailing Edge and outputting them as trackwise files. Code developed
% by Emanuele Santi (e.santi@ifac.cnr.it) by reapprising previous
% implementations hosted at Università la Sapienza - Rome (IT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clearvars -except configurationPath

ex=exist('configurationPath') ;
if ex ==0
    mode="GUI" ;
% check if configuration file exist in conf folder, otherwise ask for it
    if exist('./conf/configuration.cfg')
    configurationPath='./conf/configuration.cfg' ; 
    else
    [configurationfile configurationPath] = uigetfile('../*.cfg', 'Select input configuration file') ; 
    configurationPath= [ configurationPath configurationfile]  ;
    end
%
else
    if ~isfile(configurationPath)
        throw(MException('INPUT:ERROR', "Cannot find configuration file. Please check the command line and try again."))
    end
    mode="input" ;
end
[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
    snr_th, rx_gain_th, inc_angl_th, nsnr_th] = ReadConfFile(configurationPath);
%
switch mode
    case "GUI" 
%
% ************* Start GUI
%
    disp('GUI mode')
% 
Answer{1}= char(Taskname) ;         Answer{2}=char(string(initdate)) ;
Answer{3}=char(enddate)  ;          Answer{4}= char(string(savespace)) ;
Answer{5}=char(CyGinpath)  ;        Answer{6}= char(string(CyGoutpath)) ;
Answer{7}= char(string(logpath)) ;
Answer{8}=char(string(LatMin))  ;   Answer{9}= char(string(LatMax)) ;
Answer{10}=char(string(LonMin))  ;  Answer{11}=char(string(LonMax)) ; 
Answer{12}=char(aggregate_data) ;   Answer{13}=char(out_format) ; 
%
% Set up prompts and dialog config
prompt = {'Outfileprefix: ', ...
          'First day to extract: ', ...
          'Last day to extract: ', ...
          'Savespace? [Yes / No]: ', ...
          'DataInputRootPath: ', ...
          'DataOutputRootPath: ', ...
          'LogsOutputRootPath: ', ...
          'Southernmost latitude: ', ...
          'Northernmost latitude: ', ...
          'Westernmost longitude: ', ...
          'Easternmost longitude: ', ...
          'Multi-day data stack [Yes / No]:',...
          'Output file format [Matlab / netcdf]:'};
%      
name = 'Extraction of CyGNSS L1b data';
numlines = repmat([1 90], 13, 1);
opts.Resize = 'on';
opts.WindowStyle = 'normal';
opts.Interpreter = 'tex';
defaultanswer={Answer{1},Answer{2},...
                 Answer{3},Answer{4},Answer{5},Answer{6},Answer{7},...
                 Answer{8},Answer{9},Answer{10},Answer{11},...
                 Answer{12}, Answer{13}} ; 
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
aggregate_data=Answer{12};
out_format=Answer{13};
% init_SM_Day=datetime(Answer{2}) ; 
% final_SM_Day=datetime(Answer{3}) ; 
% write the new configuration file
WriteConfig(configurationPath, Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
     snr_th, rx_gain_th, inc_angl_th, nsnr_th);
aggregate_data = strcmpi(strtrim(Answer{12}), "Yes");  % numeric switch to aggregate data from different days and save it in a single file
%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initdatenum=datenum(initdate,'yyyy-mm-dd');
% enddatenum=datenum(enddate,'yyyy-mm-dd');
% datelist=initdatenum:enddatenum;


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI % ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "input" 
    disp('input mode')
[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format, ...
    snr_th, sp_rx_gain_th, inc_angl_th, nsnr_th] = ReadConfFile(configurationPath);
% end switch between GUI and input
end




% clear all
% close all
% addpath('src/CyGNSS_Extraction/new_version/functions/')
% 
% % ************* Start GUI
% conf_file = 'D:\Hamed\HydroGNSS_CalVal\CyGNSS_Extraction\conf\Configuration.mat';
% 
% % Set up prompts and dialog config
% prompt = {'Taskname: ', ...
%           'initdate: ', ...
%           'enddate: ', ...
%           'savespace: ', ...
%           'CyGinpath: ', ...
%           'CyGoutpath: ', ...
%           'LatMin: ', ...
%           'LatMax: ', ...
%           'LonMin: ', ...
%           'LonMax: ', ...
%           'aggregate_data'};
% 
% name = 'Soil moisture L2 processor by Sapienza-CRAS';
% numlines = repmat([1 30], 11, 1);
% opts.Resize = 'on';
% opts.WindowStyle = 'normal';
% opts.Interpreter = 'tex';

% % Try to load previous inputs if they exist
% defaultanswer = repmat({''}, 1, 11);
% if isfile(conf_file)
%     loaded = load(conf_file, 'Answer');
%     if isfield(loaded, 'Answer') && numel(loaded.Answer) >= 11
%         defaultanswer = loaded.Answer;
%     end
% end
% 
% % Launch input dialog
% Answer = inputdlg(prompt, name, numlines, defaultanswer, opts);
% 
% % If user clicked OK, save it
% if ~isempty(Answer)
%     save(conf_file, 'Answer', '-append');
% end
% 
% project_name= Answer{1};
% CyGinpath= Answer{5};
% CyGoutpath = Answer{6};
% verifydir(CyGoutpath)
% 
% logpath = fullfile(CyGoutpath, 'log');
% if ~exist(logpath, 'dir')
%     mkdir(logpath);
% end
% verifydir(logpath)
% 
% CyGfigurepath = fullfile(CyGoutpath, 'figure');
% if ~exist(CyGfigurepath, 'dir')
%     mkdir(CyGfigurepath);
% end
% verifydir(CyGfigurepath)
% 
% 
% save('D:\Hamed\HydroGNSS_CalVal\CyGNSS_Extraction\conf\Configuration.mat', 'Answer', '-append') ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GUI % ends here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                  
% 
% %%%%%%%%%%%%%%%%%%% DEFINING GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% savespace=Answer{4};                                                            % to apply the CYGNSS land flag before saving (it significantly reduces the size of output file and speeds up the processing
% aggregate_data = strcmpi(strtrim(Answer{11}), "true");                                                     % to aggregate data from different days and save it in a single file
% %%%%%%%%%%%%%%%%%%%%%%%%%% Geographic limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LatMin=Answer{7};
% LatMax=Answer{8};
% LonMin=Answer{9};
% LonMax=Answer{10};

%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    agg_GAIN=[];                                % gain of receiver antenna [dBi]
    agg_EIRP=[];                                % EIRP [W]
    agg_SNR_L1_L=[];                              % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
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
    agg_REFLECTIVITY_LINEAR_L1_L=[];            % Reflectivity
    agg_BRCS=[];                                % added by Hamed to save full ddm
    agg_REFLECTIVITY_PEAK=[] ;                  % Reflectivity peak from CyGNSS L1b products
    agg_QC_2=[];                                % Second Quality Flag
    agg_COHERENCY_RATIO=[] ;                    % Coherency ratio from CyGNSS L1b products
    agg_DDM_LES=[] ;                            % DDM_LES from CyGNSS L1b products

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
    disp(['% now processing day ' num2str(ii) ' out of ' num2str(length(datelist)) ', DoY: ' num2str(doy) ', date: ' datestr(datelist(ii),'yyyy-mm-dd')]);
    
    %%%%%%%%%%%%%%%%%%%%%% Defining paths for each day %%%%%%%%%%%%%%%%%%%%%%%%%
    DoY_infolderpath     = [CyGinpath, '\' , Year, '\', num2str(doy), '\'];        % Path to DoY folders containing input CyGNSS .nc data


    %%%%%%%%%%%%%%%%%%%%%% CyGNSS data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chkCyGNSSfile=dir([DoY_infolderpath 'cyg0*.ddmi.s' datechar '*.nc']);
    if  ~isempty(chkCyGNSSfile)
        disp('% Extracting CyGNSS data ...')
        [DoY,SoD,spacecraft_num,pseudo_random_noise,SPLAT,SPLON,THETA,GAIN, EIRP,SNR_L1_L,PHI_Initial_sp_az_orbit, ...
            REFLECTIVITY_LINEAR_L1_L,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,PA,QC,noise_floor,BRCS,...
            REFLECTIVITY_PEAK, QC_2 , COHERENCY_RATIO, DDM_LES]= ...
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
            agg_GAIN=cat(1,agg_GAIN, GAIN) ; 
            agg_EIRP=cat(1,agg_EIRP, EIRP(:));
            agg_SNR_L1_L=cat(1,agg_SNR_L1_L, SNR_L1_L(:));
            agg_PHI_Initial_sp_az_orbit=cat(1,agg_PHI_Initial_sp_az_orbit, PHI_Initial_sp_az_orbit(:));
            agg_REFLECTIVITY_LINEAR_L1_L=cat(1,agg_REFLECTIVITY_LINEAR_L1_L,REFLECTIVITY_LINEAR_L1_L(:));
            agg_KURTOSIS=cat(1,agg_KURTOSIS, KURTOSIS(:));
            agg_KURTOSIS_DOPP_0=cat(1,agg_KURTOSIS_DOPP_0, KURTOSIS_DOPP_0(:)); 
            agg_TE_WIDTH=cat(1,agg_TE_WIDTH, TE_WIDTH(:)); 
            agg_DDM_NBRCS=cat(1,agg_DDM_NBRCS, DDM_NBRCS(:)); 
            agg_PA=cat(1,agg_PA, PA(:));
            agg_QC=cat(1,agg_QC, QC(:)); 
            agg_NF=cat(1,agg_NF, noise_floor(:));
            agg_BRCS=cat(3, agg_BRCS, BRCS);   
            agg_REFLECTIVITY_PEAK=cat(1, agg_REFLECTIVITY_PEAK, REFLECTIVITY_PEAK(:)) ; 
            agg_QC_2=cat(1,agg_QC_2, QC_2(:)); 
            agg_COHERENCY_RATIO=cat(1, agg_COHERENCY_RATIO, COHERENCY_RATIO(:)) ;
            agg_DDM_LES=cat(1, agg_DDM_LES,DDM_LES(:)) ; 

            % agg_RXRANGE=cat(1,agg_RXRANGE,RXRANGE); % these variables are extracted in extract_CyGNSS function, but then they are not passed to the function output. Ask Hamed why
            % agg_TXRANGE=cat(1,agg_TXRANGE,TXRANGE);
            % agg_NST=cat(1,agg_NST,NST);
        else
                        % Quality flag 2 - currently using bits:
                        % 17 (low_confidence_gps_eirp_estimate)
                        % 22 (gps_pvt_sp3_error)
                        % 28 (low_quality_gps_ant_knowledge)
                        % 30 (anomalous_sampling_period)
                    oqf1=(bitget(QC,17) | bitget(QC,22) | bitget(QC,28) | bitget(QC,30));
                    % Quality flag 2 - currently using bits:
                        % 12 (overall)
                        % 18 (preliminary_gps_ant_knowledge)
                    oqf2=(bitget(QC_2,12) | bitget(QC_2,18));
                    NOT_TOBE_USED=(oqf1|oqf2) ;                            % Not to be uses sample logical QC index. It is '1' if sample is not recommended

                    NOT_RECOMMENDED= SNR_L1_L > snr_th & ...               % Not recommende logical QC index. It is '1' if sample is suspicious
                    GAIN > rx_gain_th & ...
                    THETA < inc_angl_th ; % & ...
                    % nsnr_dB < nsnr_th;                             
%
       %%%%% save individual day data 
       disp('% saving CyGNSS data daily');
%
            save([CyGoutpath, '/', project_name '_' daterangechar '.mat'], 'Year', 'DoY', 'SoD', 'spacecraft_num', ...  
                'pseudo_random_noise', 'SPLAT', 'SPLON', 'THETA', 'GAIN', 'EIRP', 'SNR_L1_L', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR_L1_L', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'noise_floor',...
                'REFLECTIVITY_PEAK', 'QC_2',  'COHERENCY_RATIO', 'DDM_LES', 'NOT_TOBE_USED', 'NOT_RECOMMENDED', '-v7.3');
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
    
    % Renaming variables
    DoY=agg_DoY;
    SoD=agg_SoD;
    spacecraft_num=agg_SCID;
    pseudo_random_noise=agg_PRN;
    SPLAT=agg_SPLAT;
    SPLON=agg_SPLON;
    THETA=agg_THETA;
    GAIN=agg_GAIN ; 
    EIRP=agg_EIRP;
    SNR_L1_L=agg_SNR_L1_L;
    PHI_Initial_sp_az_orbit=agg_PHI_Initial_sp_az_orbit;
    REFLECTIVITY_LINEAR_L1_L=agg_REFLECTIVITY_LINEAR_L1_L;
    KURTOSIS=agg_KURTOSIS;
    KURTOSIS_DOPP_0=agg_KURTOSIS_DOPP_0; 
    TE_WIDTH=agg_TE_WIDTH; 
    DDM_NBRCS=agg_DDM_NBRCS; 
    PA=agg_PA;
    QC=agg_QC; 
    noise_floor=agg_NF;
    BRCS=agg_BRCS;   
    REFLECTIVITY_PEAK=agg_REFLECTIVITY_PEAK ; 
    QC_2=agg_QC_2; 
    COHERENCY_RATIO=agg_COHERENCY_RATIO ; 
    DDM_LES=agg_DDM_LES ; 
    % apply quality check and filtering
                        % Quality flag 1 - currently using bits:
                        % 17 (low_confidence_gps_eirp_estimate)
                        % 22 (gps_pvt_sp3_error)
                        % 28 (low_quality_gps_ant_knowledge)
                        % 30 (anomalous_sampling_period)
                    oqf1=(bitget(QC,17) | bitget(QC,22) | bitget(QC,28) | bitget(QC,30));
                        % Quality flag 2 - currently using bits:
                        % 12 (overall)
                        % 18 (preliminary_gps_ant_knowledge)
                    oqf2=(bitget(QC_2,12) | bitget(QC_2,18));
                    NOT_TOBE_USED=(oqf1|oqf2) ;                            % Not to be uses sample logical QC index. It is '1' if sample is not recommended

                    NOT_RECOMMENDED=not( SNR_L1_L > snr_th & ...               % Not recommende logical QC index. It is '1' if sample is suspicious
                    GAIN > rx_gain_th & ...
                    THETA < inc_angl_th) ; % & ...
                    % nsnr_dB < nsnr_th;                             
    %
    % Saving aggregated data
    disp('% Saving aggregated data in a single output file')
%
    save([CyGoutpath, '/', project_name '_' daterangechar '.mat'], 'Year', 'DoY', 'SoD', 'spacecraft_num', ...  
                'pseudo_random_noise', 'SPLAT', 'SPLON', 'THETA', 'GAIN', 'EIRP', 'SNR_L1_L', 'PHI_Initial_sp_az_orbit', ...
                'REFLECTIVITY_LINEAR_L1_L', 'KURTOSIS', 'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','PA','QC', 'noise_floor',...
                 'REFLECTIVITY_PEAK', 'QC_2',  'COHERENCY_RATIO', 'DDM_LES', 'NOT_TOBE_USED', 'NOT_RECOMMENDED','-v7.3');
end
s=duration(0,0,toc);
close all
disp(['total duration is ' char(duration(0,0,toc))])