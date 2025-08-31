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
function [DoY,SoD,SCID,PRN,SPLAT,SPLON,THETA,EIRP,SNR,PHI_Initial_sp_az_orbit, ...
        REFLECTIVITY_LINEAR,KURTOSIS,KURTOSIS_DOPP_0,TE_WIDTH,DDM_NBRCS,PA,QC,NF,LF, BRCS]= ...
        extract_CyGNSS(nsat,datechar,doy,inpath,logpath,lambda,Doppler_bins,savespace,delay_vector,Power_threshold)
    %%%%%%%%%%%%%%%%%%%% INITIALISING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SCID=[];                                % CYGNSS sat ID
    SoD=[];                                 % second of the day
    DoY=[];                                 % day of the year
    PRN=[];                                 % PRN --> Prn code = prn -->trasmettitore (1 10 22 etc..)
    SPLAT=[];                               % SP lat on ground
    SPLON=[];                               % SP lon on ground
    THETA=[];                               % incidence angle
    PHI_Initial_sp_az_orbit=[];             % azimuth angle in specular point orbit frame
    GAIN=[];                                % gain of receiver antenna [dBi]
    EIRP=[];                                % EIRP [W]
    SNR=[];                                 % SNR of reflected signal - NOTE: calculated from the uncalibrated DDM in counts [dB]
    PA=[];                                  % peak power
    NF=[];                                  % noise floor
    RXRANGE=[];                             % Rx range [m]
    TXRANGE=[];                             % Tx range [m]
    NST=[];                                 % overall quality
    LF=[];                                  % land flag
    QC=[];                                  % Quality Flag
    DDM_NBRCS=[];                           % NBRCS
    KURTOSIS=[];                            % Kurtosis
    KURTOSIS_DOPP_0=[];                     % Kurtosis zero-doppler
    TE_WIDTH = [];                          % Trailing Edge (Carreno-Luengo 2020)
    REFLECTIVITY_LINEAR=[];                 % Reflectivity

    BRCS=[];                                % added by Hamed to save full ddm

    for jj=1:nsat     % loop on  the 8 satellites   
        chkfile=dir([inpath 'cyg0' num2str(jj) '.ddmi.s' datechar '*.nc']);                    % to avoid end of execution in case file is missing
        if ~isempty(chkfile)    
            infile=chkfile.name;  
            disp(['% reading satellite ' num2str(jj) ' - file ' infile ])
            [sp_lat,sp_lon,scid,ts,nst_full,prn,theta,phi_Initial_sp_az_orbit,gain, ...
                eirp,snr,nf,rxrange,txrange,ddm_nbrcs,qc,lf,pa,reflectivity_linear,Kurtosis,Kurtosis_dopp0, brcs]=readnc_CyGNSS(inpath,infile,lambda,Doppler_bins,savespace); 
        
            disp('% computing  Trailing Edge') %Kurtosis, Kurtosis zero doppler and
            TE_width=computeTE(pa,delay_vector,Power_threshold);         
            dayofyear=zeros(size(sp_lat)) + doy;  % to have the same size as sp_lat
        % cat variables
            disp('% cat variables ')
            SCID=cat(1,SCID,scid(:));
            SoD=cat(1,SoD,ts(:));
            DoY=cat(1,DoY,dayofyear(:));
            PRN=cat(1,PRN, prn(:));
            SPLAT=cat(1,SPLAT, sp_lat(:));
            SPLON=cat(1,SPLON, sp_lon(:));
            THETA=cat(1,THETA, theta(:));
            PHI_Initial_sp_az_orbit=cat(1,PHI_Initial_sp_az_orbit, phi_Initial_sp_az_orbit(:));
            GAIN=cat(1,GAIN, gain(:));
            EIRP=cat(1,EIRP, eirp(:));
            SNR=cat(1,SNR, snr(:));
            LF=cat(1,LF,lf(:));
            QC=cat(1,QC, qc(:)); 
            PA=cat(1,PA, reflectivity_linear(:));
            NF=cat(1,NF, nf(:));
            DDM_NBRCS=cat(1,DDM_NBRCS, ddm_nbrcs(:)); 
            KURTOSIS=cat(1,KURTOSIS, Kurtosis(:));
            KURTOSIS_DOPP_0=cat(1,KURTOSIS_DOPP_0, Kurtosis_dopp0(:)); 
            TE_WIDTH=cat(1,TE_WIDTH, TE_width(:)); 
            REFLECTIVITY_LINEAR=cat(1,REFLECTIVITY_LINEAR,reflectivity_linear(:));
            RXRANGE=cat(1,RXRANGE,rxrange);
            TXRANGE=cat(1,TXRANGE,txrange);
            NST=cat(1,NST,nst_full);

            BRCS=cat(3, BRCS, brcs);                                       % added by Hamed to keep full ddm
        else
            diary([logpath 'log_' datestr(now,'dd-mm-yyyy') '.txt'])
            disp(['% WARNING: cyg0' num2str(jj) ' satellite missing for the date ' datechar])
            diary off
        end
    end
end   