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
function [UTC_Time,DoY,SoD,SCID,transmittingSpacecraft,PRN,spAzimuthAngleDegOrbit,specularPointLat,specularPointLon,incidenceAngleDeg,rxAntennaGain, EIRP,SNR,PHI_Initial_sp_az_orbit, ...
        reflectivityLinear_L1_L,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,powerAnalogW,qualityControlFlags,NF,BRCS, ...
        REFLECTIVITY_PEAK, receivingAntenna, qualityControlFlags_2, coherencyRatio, DDM_LES, PR]= ...
        extract_CyGNSS(nsat,datechar,doy,inpath,logpath,lambda,Doppler_bins,savespace,delay_vector,Power_threshold)
    %%%%%%%%%%%%%%%%%%%% INITIALISING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SCID=[];                                % CYGNSS sat ID
    UTC_Time=[];                            % UTC Time
    Year=[];                                % Year
    SoD=[];                                 % second of the day
    DoY=[];                                 % day of the year
    transmittingSpacecraft=[];              % GPS unique vehicle number that trnsmitt the PRN code
    PRN=[];                                 % PRN --> Prn code = prn -->trasmettitore (1 10 22 etc..)
    spAzimuthAngleDegOrbit=[];              %Specular Point Azimuth angle respect to the orcit
    specularPointLat=[];                               % SP lat on ground
    specularPointLon=[];                               % SP lon on ground
    incidenceAngleDeg=[];                               % incidence angle
    PHI_Initial_sp_az_orbit=[];             % azimuth angle in specular point orbit frame
    rxAntennaGain=[];                                % gain of receiver antenna [dBi]
    EIRP=[];                                % EIRP [W]
    SNR=[];                                 % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    powerAnalogW=[];                                  % peak power
    NF=[];                                  % noise floor
    RXRANGE=[];                             % Rx range [m]
    TXRANGE=[];                             % Tx range [m]
    NST=[];                                 % overall quality
    qualityControlFlags=[];                 % Quality Flag
    DDM_NBRCS=[];                           % NBRCS
    KURTOSIS=[];                            % Kurtosis
    KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    reflectivityLinear_L1_L=[];                 % Reflectivity

    BRCS=[];                                % added by Hamed to save full ddm

%%%%%%%%%%%%%%%%%%%%% added by Mauro
    REFLECTIVITY_PEAK=[] ;                  % Reflectivity in lienar units from L1b producxt  
    receivingAntenna=[] ;                   % Receiving Antenna (In CyGNSS it can be  0=none, 1=never used,  2= nadir starboard,  3= nadir port [integer].)
    qualityControlFlags_2=[] ;                               % Second quality control flag from L1b product
    coherencyRatio=[] ;                    % Coherency ration from L1b produt
    DDM_LES=[]  ;                           % DDM_LES from L1b product
    PR=[] ;                                 % Power ration from Mohammad M. Al-Khaldi et al., 2021

%%%%%%%%%%%%%%%%%%%%%
    for jj=1:nsat     % loop on  the 8 satellites   
        chkfile=dir([inpath 'cyg0' num2str(jj) '.ddmi.s' datechar '*.nc']);                    % to avoid end of execution in case file is missing
        if ~isempty(chkfile)    
            infile=chkfile.name;  
            disp(['% reading satellite ' num2str(jj) ' - file ' infile ])
             [sp_lat,sp_lon,theta,scid,strans,ts,nst_full,prn,azimuth_angle,phi_Initial_sp_az_orbit,sp_rx_gain, ...
    eirp,snr,nf,rxrange,txrange,ddm_nbrcs,qc,pa,Reflectivity_linear,Kurtosis,...
    Kurtosis_dopp0, brcs, reflectivity_peak, receivingantenna, qc_2, coherencyratio, ddm_les]=readnc_CyGNSS_v2(inpath,infile,lambda,Doppler_bins,savespace); 
        
            disp('% computing  Trailing Edge')                     
            TE_width=computeTE(pa,delay_vector,Power_threshold);  
            disp('% computing  peak ratio')                       
            [pr, index] = detect_coherence_v2(pa,snr) ;
            dayofyear=zeros(size(sp_lat)) + doy;  % to have the same size as sp_lat
            disp dayofyear
            year=zeros(size(sp_lat)) + str2double(datechar(1:4)); % to have the same size as sp_lat also for the year, since we are interested in the UTC time
        % cat variables
            disp('% cat variables ')
            SCID=cat(1,SCID,scid(:));
            Year=cat(1,Year,year(:));
            SoD=cat(1,SoD,ts(:));
            DoY=cat(1,DoY,dayofyear(:));
            transmittingSpacecraft=cat(1,transmittingSpacecraft,strans(:));
            UTC_Time=datetime(Year, 1, 1) + days(DoY - 1) + seconds(SoD);
            PRN=cat(1,PRN, prn(:));
            spAzimuthAngleDegOrbit=cat(1,spAzimuthAngleDegOrbit, azimuth_angle(:)); 
            specularPointLat=cat(1,specularPointLat, sp_lat(:));
            specularPointLon=cat(1,specularPointLon, sp_lon(:));
            incidenceAngleDeg=cat(1,incidenceAngleDeg, theta(:));
            PHI_Initial_sp_az_orbit=cat(1,PHI_Initial_sp_az_orbit, phi_Initial_sp_az_orbit(:));
            rxAntennaGain=cat(1,rxAntennaGain, sp_rx_gain(:));
            EIRP=cat(1,EIRP, eirp(:));
            SNR=cat(1,SNR, snr(:));
            qualityControlFlags=cat(1,qualityControlFlags, qc(:)); 
            powerAnalogW=cat(1,powerAnalogW, pa(:));
            NF=cat(1,NF, nf(:));
            DDM_NBRCS=cat(1,DDM_NBRCS, ddm_nbrcs(:)); 
            KURTOSIS=cat(1,KURTOSIS, Kurtosis(:));
            KURTOSIS_DOPP_0=cat(1,KURTOSIS_DOPP_0, Kurtosis_dopp0(:)); 
            TE_WIDTH=cat(1,TE_WIDTH, TE_width(:)); 
            reflectivityLinear_L1_L=cat(1,reflectivityLinear_L1_L,Reflectivity_linear(:));
            RXRANGE=cat(1,RXRANGE,rxrange);
            TXRANGE=cat(1,TXRANGE,txrange);
            NST=cat(1,NST,nst_full);

            BRCS=cat(3, BRCS, brcs);                                       % added by Hamed to keep full ddm

            REFLECTIVITY_PEAK=cat(1,REFLECTIVITY_PEAK,reflectivity_peak);   
            receivingAntenna=cat(1,receivingAntenna,receivingantenna);

            qualityControlFlags_2=cat(1,qualityControlFlags_2, qc_2(:));
            coherencyRatio=cat(1,coherencyRatio,coherencyratio);
            DDM_LES=cat(1,DDM_LES,ddm_les);    
            PR=cat(1,PR, pr) ; 
        else
            diary([logpath 'log_' datestr(now,'dd-mm-yyyy') '.txt'])
            disp(['% WARNING: cyg0' num2str(jj) ' satellite missing for the date ' datechar])
            diary off
        end
    end
end   