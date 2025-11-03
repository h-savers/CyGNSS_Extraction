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
Answer{7}=char(string(logpath)) ;
Answer{8}=char(string(LatMin))  ;   Answer{9}= char(string(LatMax)) ;
Answer{10}=char(string(LonMin))  ;  Answer{11}=char(string(LonMax)) ; 
Answer{12}=char(aggregate_data) ;   Answer{13}=char(out_format) ; 
%
% Set up prompts and dialog config
prompt = {'Outfileprefix: ', ...
          'First day to extract: ', ...
          'Last day to extract: ', ...
          'Filter out data over ocean? [Yes / No]: ', ...
          'DataInputRootPath: ', ...
          'DataOutputRootPath: ', ...
          'LogsOutputRootPath: ', ...
          'Southernmost latitude [>= -90 deg]: ', ...
          'Northernmost latitude [<= 90 deg]: ', ...
          'Westernmost longitude [>= -180 deg]: ', ...
          'Easternmost longitude [<= 180 deg)]: ', ...
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


% Check if the output file format has been correctly insert
if ~contains(out_format, 'matlab', 'IgnoreCase', true) || ...
   ~contains(out_format, 'netcdf', 'IgnoreCase', true)
    disp('Output file format is not supported, please check the input');
end
% Check if the "Filter out field" has been correclt insert
if strcmpi(savespace,"Yes")
    data_coverage='land';
elseif strcmpi(savespace,"No")
    data_coverage='global';
else
    disp('The Filter out field was filled in incorrectly. Please enter either Yes or No')
end

% init_SM_Day=datetime(Answer{2}) ; 
% final_SM_Day=datetime(Answer{3}) ; 
% write the new configuration file
WriteConfig(configurationPath, Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
     snr_th, rx_gain_th, inc_angl_th, nsnr_th);
aggregate_data = strcmpi(strtrim(Answer{12}), "Yes");  % numeric switch to aggregate data from different days and save it in a single file
LatMin=double(string(LatMin)) ; LatMax=double(string(LatMax)) ; 
LonMin=double(string(LonMin)) ; LonMax=double(string(LonMax)) ; 

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
    agg_UTC_Time=[];
    agg_SoD=[];                                 % second of the day
    agg_DoY=[];                                 % day of the year
    agg_transmittingSpacecraft=[];
    agg_PRN=[];                                 % PRN --> Prn code = prn -->trasmettitore (1 10 22 etc..)
    agg_spAzimuthAngleDegOrbit=[];
    agg_specularPointLat=[];                               % SP lat on ground
    agg_specularPointLon=[];                               % SP lon on ground
    agg_incidenceAngleDeg=[];                               % incidence angle
    agg_PHI_Initial_sp_az_orbit=[];             % azimuth angle in specular point orbit frame
    agg_rxAntennaGain=[];                                % gain of receiver antenna [dBi]
    agg_EIRP=[];                             % EIRP [W]
    agg_SNR_L1_L=[];                            % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    agg_PA_L1_L=[];                             % peak power
    agg_NF=[];                                  % noise floor
    agg_RXRANGE=[];                             % Rx range [m]
    agg_TXRANGE=[];                             % Tx range [m]
    agg_NST=[];                                 % overall quality
    agg_qualityControlFlags=[];                                  % Quality Flag
    agg_DDM_NBRCS=[];                           % NBRCS
    agg_KURTOSIS=[];                            % Kurtosis
    agg_KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    agg_TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    agg_reflectivityLinear_L1_L=[];            % Reflectivity
    agg_BRCS=[];                                % added by Hamed to save full ddm
    agg_REFLECTIVITY_PEAK_L1_L=[] ;             % Reflectivity peak from CyGNSS L1b products
    agg_receivingAntenna=[] ;
    agg_qualityControlFlags_2=[];                                % Second Quality Flag
    agg_coherencyRatio=[] ;                    % Coherency ratio from CyGNSS L1b products
    agg_DDM_LES=[] ;                            % DDM_LES from CyGNSS L1b products
    agg_powerRatio=[] ;                        % Power raio from Mohammad M. Al-Khaldi et al., 2021

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
    DoY_infolderpath     = [CyGinpath, '/' , Year, '/', num2str(doy, '%03.0f'), '/'];        % Path to DoY folders containing input CyGNSS .nc data


    %%%%%%%%%%%%%%%%%%%%%% CyGNSS data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chkCyGNSSfile=dir([DoY_infolderpath 'cyg0*.ddmi.s' datechar '*.nc']);
    if  ~isempty(chkCyGNSSfile)
        disp('% Extracting CyGNSS data ...')
        [UTC_Time,DoY,SoD,receivingSpacecraft,transmittingSpacecraft,pseudoRandomNoise,spAzimuthAngleDegOrbit,specularPointLat,specularPointLon,incidenceAngleDeg,rxAntennaGain, EIRP,SNR_L1_L,PHI_Initial_sp_az_orbit, ...
            reflectivityLinear_L1_L,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,powerAnalogW,qualityControlFlags,noiseFloorCounts,BRCS,...
            REFLECTIVITY_PEAK_L1_L, receivingAntenna, qualityControlFlags_2 , coherencyRatio, DDM_LES, powerRatio]= ...
            extract_CyGNSS(nsat,datechar,doy,DoY_infolderpath,logpath,lambda,Doppler_bins,savespace,delay_vector,Power_threshold);            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aggregate_data
            disp(['% cat variables from day ', datechar ' to aggregated output file']);
            agg_DoY=cat(1,agg_DoY,DoY(:));
            agg_SoD=cat(1,agg_SoD,SoD(:));
            agg_UTC_Time=cat(1,agg_UTC_Time,UTC_Time(:));
            agg_SCID=cat(1,agg_SCID,receivingSpacecraft(:));
            agg_transmittingSpacecraft=cat(1,agg_transmittingSpacecraft,transmittingSpacecraft(:));
            agg_PRN=cat(1,agg_PRN, pseudoRandomNoise(:));
            agg_spAzimuthAngleDegOrbit=cat(1,agg_spAzimuthAngleDegOrbit,spAzimuthAngleDegOrbit(:));
            agg_specularPointLat=cat(1,agg_specularPointLat, specularPointLat(:));
            agg_specularPointLon=cat(1,agg_specularPointLon, specularPointLon(:));
            agg_incidenceAngleDeg=cat(1,agg_incidenceAngleDeg, incidenceAngleDeg(:));
            agg_rxAntennaGain=cat(1,agg_rxAntennaGain, rxAntennaGain) ; 
            agg_EIRP=cat(1,agg_EIRP, EIRP(:));
            agg_SNR_L1_L=cat(1,agg_SNR_L1_L, SNR_L1_L(:));
            agg_PHI_Initial_sp_az_orbit=cat(1,agg_PHI_Initial_sp_az_orbit, PHI_Initial_sp_az_orbit(:));
            agg_reflectivityLinear_L1_L=cat(1,agg_reflectivityLinear_L1_L,reflectivityLinear_L1_L(:));
            agg_KURTOSIS=cat(1,agg_KURTOSIS, KURTOSIS(:));
            agg_KURTOSIS_DOPP_0=cat(1,agg_KURTOSIS_DOPP_0, KURTOSIS_DOPP_0(:)); 
            agg_TE_WIDTH=cat(1,agg_TE_WIDTH, TE_WIDTH(:)); 
            agg_DDM_NBRCS=cat(1,agg_DDM_NBRCS, DDM_NBRCS(:)); 
            agg_PA_L1_L=cat(1,agg_PA_L1_L, powerAnalogW(:));
            agg_qualityControlFlags=cat(1,agg_qualityControlFlags, qualityControlFlags(:)); 
            agg_NF=cat(1,agg_NF, noiseFloorCounts(:));
            agg_BRCS=cat(3, agg_BRCS, BRCS);   
            agg_REFLECTIVITY_PEAK_L1_L=cat(1, agg_REFLECTIVITY_PEAK_L1_L, REFLECTIVITY_PEAK_L1_L(:)) ; 
            agg_receivingAntenna=cat(1,agg_receivingAntenna, receivingAntenna(:)) ;
            agg_qualityControlFlags_2=cat(1,agg_qualityControlFlags_2, qualityControlFlags_2(:)) ; 
            agg_coherencyRatio=cat(1, agg_coherencyRatio, coherencyRatio(:)) ;
            agg_DDM_LES=cat(1, agg_DDM_LES,DDM_LES(:)) ; 
            agg_powerRatio=cat(1, agg_powerRatio,powerRatio(:)) ; 


            % agg_RXRANGE=cat(1,agg_RXRANGE,RXRANGE); % these variables are extracted in extract_CyGNSS function, but then they are not passed to the function output. Ask Hamed why
            % agg_TXRANGE=cat(1,agg_TXRANGE,TXRANGE);
            % agg_NST=cat(1,agg_NST,NST);
        else
                        % Quality flag 2 - currently using bits:
                        % 17 (low_confidence_gps_eirp_estimate)
                        % 22 (gps_pvt_sp3_error)
                        % 28 (low_quality_gps_ant_knowledge)
                        % 30 (anomalous_sampling_period)
                    oqf1=(bitget(qualityControlFlags,17) | bitget(qualityControlFlags,22) | bitget(qualityControlFlags,28) | bitget(qualityControlFlags,30));
                    % Quality flag 2 - currently using bits:
                        % 12 (overall)
                        % 18 (preliminary_gps_ant_knowledge)
                    oqf2=(bitget(qualityControlFlags_2,12) | bitget(qualityControlFlags_2,18));
                    notToBeUsed=(oqf1|oqf2) ;                            % Not to be uses sample logical qualityControlFlags index. It is '1' if sample is not recommended

                    notRecommended= SNR_L1_L > snr_th & ...               % Not recommende logical qualityControlFlags index. It is '1' if sample is suspicious
                    rxAntennaGain > rx_gain_th & ...
                    incidenceAngleDeg < inc_angl_th ; % & ...
                    % nsnr_dB < nsnr_th;                             
%
       %%%%% save individual day data 
  disp('% saving CyGNSS data daily');
  daterangechar=['_' datechar] ; 
%
  disp('% Saving aggregated data in a single output file')
    if LatMin ~= -90 & LatMax  ~= 90 & LonMin  ~= -180 & LonMax ~= 180 
        subgeo=find(specularPointLat >= LatMin & specularPointLat <= LatMax & specularPointLon >= LonMin & specularPointLon >= LonMax ) ; 
        DoY=DoY(subgeo) ; SoD=SoD(subgeo) ; receivingSpacecraft=receivingSpacecraft(subgeo) ; transmittingSpacecraft=transmittingSpacecraft(subgeo);
        pseudoRandomNoise=pseudoRandomNoise(subgeo) ; spAzimuthAngleDegOrbit=spAzimuthAngleDegOrbit(subgeo); specularPointLat=specularPointLat(subgeo) ; specularPointLon=specularPointLon(subgeo) ; incidenceAngleDeg=incidenceAngleDeg(subgeo) ;
        rxAntennaGain=rxAntennaGain(subgeo) ; EIRP=EIRP(subgeo) ; SNR_L1_L=SNR_L1_L(subgeo) ; PHI_Initial_sp_az_orbit=PHI_Initial_sp_az_orbit(subgeo) ;
        reflectivityLinear_L1_L=reflectivityLinear_L1_L(subgeo) ; KURTOSIS=KURTOSIS(subgeo) ; KURTOSIS_DOPP_0=KURTOSIS_DOPP_0(subgeo) ; 
        TE_WIDTH=TE_WIDTH(subgeo) ; DDM_NBRCS=DDM_NBRCS(subgeo) ; powerAnalogW=powerAnalogW(subgeo) ; qualityControlFlags=qualityControlFlags(subgeo) ; noiseFloorCounts=noiseFloorCounts(subgeo) ;
        REFLECTIVITY_PEAK_L1_L=REFLECTIVITY_PEAK_L1_L(subgeo) ; receivingAntenna=receivingAntenna(subgeo) ; qualityControlFlags_2=qualityControlFlags_2(subgeo) ;  coherencyRatio=coherencyRatio(subgeo) ;
        DDM_LES=DDM_LES(subgeo) ; powerRatio=powerRatio(subgeo) ; notToBeUsed=notToBeUsed(subgeo) ; notRecommended=notRecommended(subgeo) ;
    end
            s=duration(0,0,toc);
            save_file(out_format, CyGoutpath, project_name, daterangechar, initdate, ...
             enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
             UTC_Time, DoY, SoD, receivingSpacecraft, transmittingSpacecraft, ...
             pseudoRandomNoise, spAzimuthAngleDegOrbit,specularPointLat, specularPointLon, incidenceAngleDeg, ...
             rxAntennaGain, EIRP, SNR_L1_L, PHI_Initial_sp_az_orbit, reflectivityLinear_L1_L, ...
             KURTOSIS, KURTOSIS_DOPP_0, TE_WIDTH, DDM_NBRCS, powerAnalogW, qualityControlFlags, ...
             noiseFloorCounts, REFLECTIVITY_PEAK_L1_L, receivingAntenna, qualityControlFlags_2, ...
             coherencyRatio, DDM_LES, powerRatio, notToBeUsed, notRecommended);
        end
    %%%%%%%%%%%%%%%%%%%%% Displaying Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %          scattermap(real(10.*log10(REFLECTIVITY_LINEAR)),specularPointLat,specularPointLon,datechar,-40,0)
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
    UTC_Time=agg_UTC_Time;
    receivingSpacecraft=agg_SCID;
    transmittingSpacecraft=agg_transmittingSpacecraft;
    pseudoRandomNoise=agg_PRN;
    spAzimuthAngleDegOrbit=agg_spAzimuthAngleDegOrbit;
    specularPointLat=agg_specularPointLat;
    specularPointLon=agg_specularPointLon;
    incidenceAngleDeg=agg_incidenceAngleDeg;
    rxAntennaGain=agg_rxAntennaGain ; 
    EIRP=agg_EIRP;
    SNR_L1_L=agg_SNR_L1_L;
    PHI_Initial_sp_az_orbit=agg_PHI_Initial_sp_az_orbit;
    reflectivityLinear_L1_L=agg_reflectivityLinear_L1_L;
    KURTOSIS=agg_KURTOSIS;
    KURTOSIS_DOPP_0=agg_KURTOSIS_DOPP_0; 
    TE_WIDTH=agg_TE_WIDTH; 
    DDM_NBRCS=agg_DDM_NBRCS; 
    powerAnalogW=agg_PA_L1_L;
    qualityControlFlags=agg_qualityControlFlags; 
    noiseFloorCounts=agg_NF;
    BRCS=agg_BRCS;   
    REFLECTIVITY_PEAK_L1_L=agg_REFLECTIVITY_PEAK_L1_L ;
    receivingAntenna=agg_receivingAntenna ; 
    qualityControlFlags_2=agg_qualityControlFlags_2; 
    coherencyRatio=agg_coherencyRatio ; 
    DDM_LES=agg_DDM_LES ; 
    powerRatio=agg_powerRatio ; 
    % apply quality check and filtering
                        % Quality flag 1 - currently using bits:
                        % 17 (low_confidence_gps_eirp_estimate)
                        % 22 (gps_pvt_sp3_error)
                        % 28 (low_quality_gps_ant_knowledge)
                        % 30 (anomalous_sampling_period)
                    oqf1=(bitget(qualityControlFlags,17) | bitget(qualityControlFlags,22) | bitget(qualityControlFlags,28) | bitget(qualityControlFlags,30));
                        % Quality flag 2 - currently using bits:
                        % 12 (overall)
                        % 18 (preliminary_gps_ant_knowledge)
                    oqf2=(bitget(qualityControlFlags_2,12) | bitget(qualityControlFlags_2,18));
                    notToBeUsed=(oqf1|oqf2) ;                            % Not to be uses sample logical qualityControlFlags index. It is '1' if sample is not recommended

                    notRecommended=( SNR_L1_L < snr_th | ...               % Not recommende logical qualityControlFlags index. It is '1' if sample is suspicious
                    rxAntennaGain < rx_gain_th | ...
                    incidenceAngleDeg > inc_angl_th) ; % & ...
                    % nsnr_dB < nsnr_th;                             
    %
    s=duration(0,0,toc);
    
    % Saving aggregated data
    save_file(out_format, CyGoutpath, project_name, daterangechar, initdate, ...
    enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
    UTC_Time, DoY, SoD, receivingSpacecraft, transmittingSpacecraft, ...
    pseudoRandomNoise, spAzimuthAngleDegOrbit,specularPointLat, specularPointLon, incidenceAngleDeg, ...
    rxAntennaGain, EIRP, SNR_L1_L, PHI_Initial_sp_az_orbit, reflectivityLinear_L1_L, ...
    KURTOSIS, KURTOSIS_DOPP_0, TE_WIDTH, DDM_NBRCS, powerAnalogW, qualityControlFlags, ...
    noiseFloorCounts, REFLECTIVITY_PEAK_L1_L, receivingAntenna, qualityControlFlags_2, ...
    coherencyRatio, DDM_LES, powerRatio, notToBeUsed, notRecommended);
%
end

close all
    disp(['total duration is ' char(duration(0,0,toc))])