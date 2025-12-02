%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts CyGNSS observables by also computing reflectivity
% and Trailing Edge and outputting them as trackwise files.

% Inputs:
%   nsat: number of satellites (8 for CyGNSS)
%   datechar: date in 'dd-mm-yyyy' format
%   doy: day of the year
%   inpath: input path for CyGNSS data
%   logpath: path for logging errors
%   lambda: wavelength of the signal
%   Doppler_bins: Doppler bins for processing
%   savespace: flag to apply the CYGNSS land flag before saving (reduces file size)
%   delay_vector: vector of delays for processing
%   Power_threshold: threshold for power to compute Trailing Edge


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mission,L1b_product,L1b_product_version,timeUTC, ...
        DoY,SoD,SCID,SV_NUM,PRN,SPLAT,SPLON,THETA,GAIN,EIRP,SNR,PHI_Initial_sp_az_orbit, ...
        REFLECTIVITY_LINEAR,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,PA_PEAK,QC,NF,BRCS, RECEIVING_ANTENNA...
        SP_AZIMUTH_ANGLE, REFLECTIVITY_PEAK, REFLECTIVITY_PEAK_CALIBRATED, QC_2, COHERENCY_RATIO, DDM_LES, PR,PSEUDOSTD,BIT_RATIO,COEFFICIENT_OF_VARIATION]= ...
        extract_CyGNSS(datechar,doy,inpath,logpath,calibration_file,lambda,Doppler_bins,savespace,delay_vector,Power_threshold)

    %%%%%%%%%%%%%%%%%%%% INITIALISING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timeUTC=[];
    Year=[]; 
    SoD=[];                                 % second of the day
    DoY=[];                                 % day of the year
    SCID=[];                                % CYGNSS sat ID
    SV_NUM=[];                              % transmitting Spacecraft number
    PRN=[];                                 % PRN --> Prn code = prn -->trasmettitore (1 10 22 etc..)
    SPLAT=[];                               % SP lat on ground
    SPLON=[];                               % SP lon on ground
    THETA=[];                               % incidence angle
    PHI_Initial_sp_az_orbit=[];             % azimuth angle in specular point orbit frame
    GAIN=[];                                % gain of receiver antenna [dBi]
    EIRP=[];                                % EIRP [W]
    SNR=[];                                 % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    PA=[];                                  % peak power
    PA_PEAK=[];
    NF=[];                                  % noise floor
    RXRANGE=[];                             % Rx range [m]
    TXRANGE=[];                             % Tx range [m]
    NST=[];                                 % overall quality
    QC=[];                                  % Quality Flag
    DDM_NBRCS=[];                           % NBRCS
    KURTOSIS=[];                            % Kurtosis
    KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    REFLECTIVITY_LINEAR=[];                 % Reflectivity
    SP_AZIMUTH_ANGLE=[];                    % Specular point azimuth
    RECEIVING_ANTENNA=[];
    BRCS=[];                                % added by Hamed to save full ddm
    PSEUDOSTD=[] ;                          % pseudo standard deviation of DDM assumed as a bivariate pdf
    BIT_RATIO=[] ;                          % added by Federico - Bit Ratio

%%%%%%%%%%%%%%%%%%%%% added by Mauro
    REFLECTIVITY_PEAK=[] ;                  % Reflectivity in lienar units from L1b product     
    REFLECTIVITY_PEAK_CALIBRATED=[] ;       % Reflectivity peak calibrated
    QC_2=[] ;                               % Second quality control flag from L1b product
    COHERENCY_RATIO=[] ;                    % Coherency ration from L1b produt
    DDM_LES=[]  ;                           % DDM_LES from L1b product
    PR=[] ;                                 % Power ration from Mohammad M. Al-Khaldi et al., 2021
    COEFFICIENT_OF_VARIATION=[] ;

    chkfile=dir([inpath 'cyg0*.nc']);
    length(chkfile);
    sat_index=[];

    for i=1:length(chkfile)        
        sat_index=[sat_index, str2num(chkfile(i).name(5))];
    end

%%%%%%%%%%%%%%%%%%%%%
    for jj=sat_index     % loop on  the 8 satellites   
        chkfile=dir([inpath 'cyg0' num2str(jj) '.ddmi.s' datechar '*.nc']);                    % to avoid end of execution in case file is missing
        if ~isempty(chkfile)    
            infile=chkfile.name;  
            disp(['% reading satellite ' num2str(jj) ' - file ' infile ])
            [mission, L1b_product, L1b_product_version,sp_lat,sp_lon,scid,sv_num,ts,nst_full, ...
            prn,theta,phi_Initial_sp_az_orbit,sp_rx_gain,eirp,snr,nf,rxrange,txrange,ddm_nbrcs,...
            qc,pa,peak,Reflectivity_linear,Kurtosis, Kurtosis_dopp0, brcs, reflectivity_peak, ...
            receivingantenna, sp_azimuth_angle_deg_north, qc_2, coherency_ratio, ddm_les, ...
            raw_counts,bit_ratio,calibration_table,sampling_seconds]=readnc_CyGNSS_v2(inpath,infile,calibration_file,lambda,Doppler_bins,savespace); 
            
            disp('% computing coefficient of variation')
            Nlook=sampling_seconds*10^3 ;
            SNRlinear=10.^snr/10 ;
            KP=sqrt((1+1./SNRlinear).^2+(1./SNRlinear).^2)/sqrt(Nlook) ;

            disp('% computing  Trailing Edge')                     
            TE_width=computeTE(double(pa),delay_vector,Power_threshold);  
            disp('% computing  peak ratio')                       
%             [pr, index] = detect_coherence_v2(pa,snr) ;
            [pr, index] = detect_coherence_v2(single(raw_counts),snr) ;
            disp('% computing pseudo standard deviation of DDM')
            pseudostd=ddm_pseudovariance(pa) ;     %    
            disp('% computing the reflectivity calibration')

            lookup = [scid(:) sv_num(:) receivingantenna(:)]; % create a composite matrix of receiving spacecraft number, spacecraft number and receiving antenna
            decimation_factor = ones(size(lookup,1),1); % create a matrix in a correct size for the decimation factor for the calibration
            [~, mask_idx] = ismember(lookup, calibration_table(:,1:3), 'rows'); % check the corrispondence between the two matrices
            valid = mask_idx > 0; % check the logic indices
            decimation_factor(valid) = calibration_table(mask_idx(valid), 4); % extract the decimation factor from the calibration_table
            valid = (reflectivity_peak ~= -9999); % check if the reflectivity_peak vector has -9999 values, in that case do not calibrate the results
            reflectivity_peak_calibrated = -9999 * ones(size(reflectivity_peak));
            reflectivity_peak_calibrated(valid) = reflectivity_peak(valid) .* decimation_factor(valid);

            dayofyear=zeros(size(sp_lat)) + doy;  % to have the same size as sp_lat
            year=zeros(size(sp_lat)) + str2double(datechar(1:4)); % to have the same size as sp_lat also for the year, since we are interested in the UTC time
        % cat variables
            disp('% cat variables ')
            SCID=cat(1,SCID,scid(:));
            Year=cat(1,Year,year(:));
            SoD=cat(1,SoD,ts(:));
            DoY=cat(1,DoY,dayofyear(:));
            timeUTC=datetime(Year, 1, 1) + days(DoY - 1) + seconds(SoD); % Let's calculate the UTC Time
            SV_NUM=cat(1,SV_NUM,sv_num(:));
            PRN=cat(1,PRN, prn(:));
            SPLAT=cat(1,SPLAT, sp_lat(:));
            SPLON=cat(1,SPLON, sp_lon(:));
            THETA=cat(1,THETA, theta(:));
            PHI_Initial_sp_az_orbit=cat(1,PHI_Initial_sp_az_orbit, phi_Initial_sp_az_orbit(:));
            GAIN=cat(1,GAIN, sp_rx_gain(:));
            EIRP=cat(1,EIRP, eirp(:));
            SNR=cat(1,SNR, snr(:));
            QC=cat(1,QC, qc(:)); 
            PA=cat(3,PA, pa);
            PA_PEAK=cat(1,PA_PEAK, peak);
            NF=cat(1,NF, nf(:));
            DDM_NBRCS=cat(1,DDM_NBRCS, ddm_nbrcs(:)); 
            KURTOSIS=cat(1,KURTOSIS, Kurtosis(:));
            KURTOSIS_DOPP_0=cat(1,KURTOSIS_DOPP_0, Kurtosis_dopp0(:)); 
            TE_WIDTH=cat(1,TE_WIDTH, TE_width(:)); 
            REFLECTIVITY_LINEAR=cat(1,REFLECTIVITY_LINEAR,Reflectivity_linear(:));
            RXRANGE=cat(1,RXRANGE,rxrange);
            TXRANGE=cat(1,TXRANGE,txrange);
            NST=cat(1,NST,nst_full);
            SP_AZIMUTH_ANGLE=cat(1,SP_AZIMUTH_ANGLE,sp_azimuth_angle_deg_north);
            RECEIVING_ANTENNA=cat(1,RECEIVING_ANTENNA,receivingantenna);
            BRCS=cat(3, BRCS, brcs);                                       % added by Hamed to keep full ddm
            REFLECTIVITY_PEAK=cat(1,REFLECTIVITY_PEAK,reflectivity_peak);  
            REFLECTIVITY_PEAK_CALIBRATED=cat(1,REFLECTIVITY_PEAK_CALIBRATED,reflectivity_peak_calibrated);
            QC_2=cat(1,QC_2, qc_2(:));
            COHERENCY_RATIO=cat(1,COHERENCY_RATIO,coherency_ratio);
            DDM_LES=cat(1,DDM_LES,ddm_les);    
            PR=cat(1,PR, pr) ; 
            PSEUDOSTD=cat(1,PSEUDOSTD, pseudostd) ; 
            BIT_RATIO=cat(1,BIT_RATIO, bit_ratio) ;
            COEFFICIENT_OF_VARIATION=cat(1,COEFFICIENT_OF_VARIATION,KP) ;

        else
            diary([logpath 'log_' datestr(now,'dd-mm-yyyy') '.txt'])
            disp(['% WARNING: cyg0' num2str(jj) ' satellite missing for the date ' datechar])
            diary off
        end
    end
end   