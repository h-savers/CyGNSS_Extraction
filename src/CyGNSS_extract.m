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
[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, calibration_file, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
    snr_th, rx_gain_th, inc_angl_th, nsnr_th, coherency_th] = ReadConfFile(configurationPath);
%
switch mode
    case "GUI" 
%
% ************* Start GUI
%
    disp('GUI mode')
% 
Answer{1}=char(Taskname) ;         Answer{2}=char(string(initdate)) ;
Answer{3}=char(enddate)  ;          Answer{4}=char(string(savespace)) ;
Answer{5}=char(CyGinpath)  ;        Answer{6}=char(string(CyGoutpath)) ;
Answer{7}=char(string(logpath)) ;   Answer{8}=char(string(calibration_file)) ;
Answer{9}=char(string(LatMin))  ;   Answer{10}=char(string(LatMax)) ;
Answer{11}=char(string(LonMin))  ;  Answer{12}=char(string(LonMax)) ; 
Answer{13}=char(aggregate_data) ;   Answer{14}=char(out_format) ; 

%
% Set up prompts and dialog config
prompt = {'Outfileprefix: ', ...
          'First day to extract: ', ...
          'Last day to extract: ', ...
          'Filter out data over ocean? [Yes / No]: ', ...
          'DataInputRootPath: ', ...
          'DataOutputRootPath: ', ...
          'LogsOutputRootPath: ', ...
          'Calibration CSV lookup table full name (''No'' if no calibration):', ...
          'Southernmost latitude [>= -90 deg]: ', ...
          'Northernmost latitude [<= 90 deg]: ', ...
          'Westernmost longitude [>= -180 deg]: ', ...
          'Easternmost longitude [<= 180 deg)]: ', ...
          'Multi-day data stack [Yes / No]:',...
          'Output file format [Matlab / netcdf]: '} ; 
%      
name = 'Extraction of CyGNSS L1b data';
numlines = repmat([1 90], 14, 1);
opts.Resize = 'on';
opts.WindowStyle = 'normal';
opts.Interpreter = 'tex';
defaultanswer={Answer{1},Answer{2},...
                 Answer{3},Answer{4},Answer{5},Answer{6},Answer{7},...
                 Answer{8},Answer{9},Answer{10},Answer{11},...
                 Answer{12}, Answer{13}, Answer{14} }; 
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
calibration_file=Answer{8};
LatMin=Answer{9};
LatMax=Answer{10};
LonMin=Answer{11};
LonMax=Answer{12};
aggregate_data=Answer{13};
out_format=Answer{14};

% ----- Added by Federico P. ----- %

% Check if the output file format has been correctly insert  
if contains(out_format, 'matlab', 'IgnoreCase', true) || ...
   contains(out_format, 'netcdf', 'IgnoreCase', true)
else
    disp('ERROR: output file format is not supported, please check the input');
    return
end
% Check if the "Filter out field" has been correclt insert
if strcmpi(savespace,"Yes")
    data_coverage='land';
elseif strcmpi(savespace,"No")
    data_coverage='global';
else
    disp('ERROR: the Filter out field was filled in incorrectly. Please enter either Yes or No')
end

% ----- End ----- %
if exist(CyGoutpath)==0 ,  disp('The output folder does not exist. Please enter an existing folder'), return, end
% init_SM_Day=datetime(Answer{2}) ; 
% final_SM_Day=datetime(Answer{3}) ; 
% write the new configuration file
WriteConfig(configurationPath, Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, calibration_file, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
     snr_th, rx_gain_th, inc_angl_th, nsnr_th, coherency_th);
aggregate_data = strcmpi(strtrim(Answer{13}), "Yes");  % numeric switch to aggregate data from different days and save it in a single file
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
[Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, calibration_file, LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
    snr_th, rx_gain_th, inc_angl_th, nsnr_th, coherency_th] = ReadConfFile(configurationPath);
aggregate_data = strcmpi(aggregate_data, "Yes");  % numeric switch to aggregate data from different days and save it in a single file
Taskname=char(Taskname) ; savespace=char(savespace); CyGinpath=char(CyGinpath) ; CyGoutpath=char(CyGoutpath) ; logpath=char(logpath) ; calibration_file=char(calibration_file) ; out_format=char(out_format) ; 
project_name=Taskname ; 
if strcmpi(savespace,"Yes")
    data_coverage='land';
elseif strcmpi(savespace,"No")
    data_coverage='global';
else
    disp('The Filter out field was filled in incorrectly. Please enter either Yes or No')
end
% end switch between GUI and input
end


% clear all
% close all
% addpath('src/CyGNSS_Extraction/new_version/functions/')
% 
% % ************* Start GUI
% conf_file = 'D:/Hamed/HydroGNSS_CalVal/CyGNSS_Extraction/conf/Configuration.mat';
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
% save('D:/Hamed/HydroGNSS_CalVal/CyGNSS_Extraction/conf/Configuration.mat', 'Answer', '-append') ;
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
% just_scale_factor: only the tx/rx ranges are needed to compute the scale
% factor, so the DDM / Doppler / wavelength parameters are not required here.

%%%%%%%%%%%%%%%%%%%%% INITIALIZING EMPTY VARIABLES FOR AGGREGATED SINGLE OUTPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if aggregate_data
    daterangechar = [datestr(initdatenum,'yyyymmdd') '-' datestr(enddatenum,'yyyymmdd')]; % date range for aggregated output file
    %disp(['% Processing data from ' daterangechar ' and saving in a single output file'])

    agg_SCID=[];                                % receiving (CYGNSS) sat ID
    agg_timeUTC=[];
    agg_SoD=[];                                 % second of the day
    agg_DoY=[];                                 % day of the year
    agg_transmittingSpacecraft=[];              % transmitting spacecraft
    agg_PRN=[];                                 % PRN code
    agg_SPLAT=[];                               % SP lat on ground
    agg_SPLON=[];                               % SP lon on ground
    agg_scaleFactor=[] ;                        % Scale factor (SS_r) computed from tx/rx ranges

else
    disp('% Processing each day separately and saving individual output files');
end

%%%%%%%%%%%%%%%%%%%%% STARTING THE MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% % figure;

for ii=1:length(datelist)     % loop on all the days 
    datechar=datestr(datelist(ii),'yyyymmdd');
    year=datechar(1:4);
    d=datetime(datechar,'InputFormat','yyyyMMdd');
    doy=day(d,'dayofyear');
    disp(['% now processing day ' num2str(ii) ' out of ' num2str(length(datelist)) ', DoY: ' num2str(doy) ', date: ' datestr(datelist(ii),'yyyy-mm-dd')]);
    
    %%%%%%%%%%%%%%%%%%%%%% Defining paths for each day %%%%%%%%%%%%%%%%%%%%%%%%%
    DoY_infolderpath     = [CyGinpath, '/' , year, '/', num2str(doy, '%03.0f'), '/'];        % Path to DoY folders containing input CyGNSS .nc data


    %%%%%%%%%%%%%%%%%%%%%% CyGNSS data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chkCyGNSSfile=dir([DoY_infolderpath 'cyg0*.ddmi.s' datechar '*.nc']);
    if  ~isempty(chkCyGNSSfile)
        %disp('% Extracting CyGNSS data ...')
        [mission,L1b_product,L1b_product_version,timeUTC, ...
            dayOfYear,secondOfDay,receivingSpacecraft,transmittingSpacecraft,pseudoRandomNoise, ...
            specularPointLat,specularPointLon,scaleFactor]= ...
            extract_CyGNSS(datechar,doy,DoY_infolderpath,logpath,savespace);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aggregate_data
            %disp(['% cat variables from day ', datechar ' to aggregated output file']);
            agg_DoY=cat(1,agg_DoY,dayOfYear(:));
            agg_SoD=cat(1,agg_SoD,secondOfDay(:));
            agg_timeUTC=cat(1,agg_timeUTC,timeUTC(:));
            agg_SCID=cat(1,agg_SCID,receivingSpacecraft(:));
            agg_transmittingSpacecraft=cat(1,agg_transmittingSpacecraft,transmittingSpacecraft(:));
            agg_PRN=cat(1,agg_PRN, pseudoRandomNoise(:));
            agg_SPLAT=cat(1,agg_SPLAT, specularPointLat(:));
            agg_SPLON=cat(1,agg_SPLON, specularPointLon(:));
            agg_scaleFactor=cat(1,agg_scaleFactor, scaleFactor(:)) ;
        else
            %%%%% save individual day data
            disp('% saving CyGNSS data daily');
            daterangechar=['_' datechar] ;
            if LatMin ~= -90 & LatMax  ~= 90 & LonMin  ~= -180 & LonMax ~= 180
                subgeo=find(specularPointLat >= LatMin & specularPointLat <= LatMax & specularPointLon >= LonMin & specularPointLon <= LonMax ) ;
                timeUTC=timeUTC(subgeo) ; receivingSpacecraft=receivingSpacecraft(subgeo) ; transmittingSpacecraft=transmittingSpacecraft(subgeo) ;
                pseudoRandomNoise=pseudoRandomNoise(subgeo) ; specularPointLat=specularPointLat(subgeo) ; specularPointLon=specularPointLon(subgeo) ; scaleFactor=scaleFactor(subgeo) ;
            end
            s=duration(0,0,toc);
            save_file(mission, out_format, L1b_product, L1b_product_version,...
            CyGoutpath, project_name, daterangechar, initdate, ...
            enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
            timeUTC, receivingSpacecraft, transmittingSpacecraft, ...
            pseudoRandomNoise, specularPointLat, specularPointLon, scaleFactor);
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
    timeUTC=agg_timeUTC;
    receivingSpacecraft=agg_SCID;
    transmittingSpacecraft=agg_transmittingSpacecraft;
    pseudoRandomNoise=agg_PRN;
    specularPointLat=agg_SPLAT;
    specularPointLon=agg_SPLON;
    scaleFactor=agg_scaleFactor ;
    %
    % Saving aggregated data
    if LatMin ~= -90 & LatMax  ~= 90 & LonMin  ~= -180 & LonMax ~= 180
        subgeo=find(specularPointLat >= LatMin & specularPointLat <= LatMax & specularPointLon >= LonMin & specularPointLon <= LonMax ) ;
        timeUTC=timeUTC(subgeo) ; receivingSpacecraft=receivingSpacecraft(subgeo) ; transmittingSpacecraft=transmittingSpacecraft(subgeo) ;
        pseudoRandomNoise=pseudoRandomNoise(subgeo) ; specularPointLat=specularPointLat(subgeo) ; specularPointLon=specularPointLon(subgeo) ; scaleFactor=scaleFactor(subgeo) ;
    end
        s=duration(0,0,toc);
        save_file(mission, out_format, L1b_product, L1b_product_version,...
        CyGoutpath, project_name, daterangechar, initdate, ...
        enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
        timeUTC, receivingSpacecraft, transmittingSpacecraft, ...
        pseudoRandomNoise, specularPointLat, specularPointLon, scaleFactor);
       
%
end
close all
disp(['total duration is ' char(duration(0,0,toc))])
