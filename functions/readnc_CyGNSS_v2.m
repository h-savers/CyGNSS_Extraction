%%%%%%%%%%%%%%%%%%%%%%%%%% function readCyGNSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads the data from the netCDF files and computes 
% Reflectivity, Kurtosis and Kurtosis zero-Doppler.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sp_lat,sp_lon,theta,scid,strans,ts,nst_full,prn,azimuth_angle,phi_Initial_sp_az_orbit,sp_rx_gain, ...
    eirp,snr,nf,rxrange,txrange,ddm_nbrcs,qc,pa,Reflectivity_linear,Kurtosis,...
    Kurtosis_dopp0, brcs, reflectivity_peak, receivingantenna, qc_2, coherencyratio, ddm_les]= ...
    readnc_CyGNSS(inpath,filename,lambda,Doppler_bins,savespace)

     % Reading data
     toread=[inpath,filename];
     % using ncread
     disp('% opening file ')
     ncid = netcdf.open(toread, 'NC_NOWRITE');
% trackNcids = netcdf.inqGrps(ncid);

     disp('% reading lat/lon')
     
     varID=netcdf.inqVarID(ncid, 'sp_lat')  ;
     sp_lat=double(netcdf.getVar(ncid,varID)) ;                            % Latitude   

     varID=netcdf.inqVarID(ncid, 'sp_lon')  ;
     sp_lon=double(netcdf.getVar(ncid,varID)) ;
     sp_lon = rem((sp_lon+180),360)-180;                                   % Longitude
     % sp_lon = rem((ncread(toread,'sp_lon')+180),360)-180;                % Longitude old version
                                      
     disp('% reading satellite ID and time ')

     varID=netcdf.inqVarID(ncid, 'spacecraft_num')  ;
     scid=netcdf.getVar(ncid,varID) ; 
     scid=double(scid).*ones(size(sp_lat));     % spacecraft id
     % scid=double(ncread(toread,'spacecraft_num')).*ones(size(sp_lat));   % spacecraft id ol version
     
     varID=netcdf.inqVarID(ncid, 'sv_num')  ;
     strans=netcdf.getVar(ncid,varID) ; 

     varID=netcdf.inqVarID(ncid, 'prn_code')  ;
     prn=double(netcdf.getVar(ncid,varID)) ;                               % full prn 
     % prn=ncread(toread,'prn_code');                                      % full prn old version 
     
     varID=netcdf.inqVarID(ncid, 'sp_az_orbit')  ;
     azimuth_angle=double(netcdf.getVar(ncid,varID)) ;                     % The mean of the orbit frame azimuth angles of the specular points of the DDMs
     % prn=ncread(toread,'prn_code');                                      % full prn old version 

     varID=netcdf.inqVarID(ncid, 'ddm_timestamp_utc')  ;
     ts=netcdf.getVar(ncid,varID) ; 
     ts=repmat(ts',4,1);                                                   % time in seconds since 2017-03-27 00:00:00.999261529
     % ts=repmat(ncread(toread,'ddm_timestamp_utc')',4,1);                 % time in seconds since 2017-03-27 00:00:00.999261529 old version
            
     disp('% reading observables ')
     varID=netcdf.inqVarID(ncid, 'nst_att_status')  ;
     nst_full=netcdf.getVar(ncid,varID) ; 
     nst_full=double(repmat(nst_full',4,1)) ;                              % attitude flag (not used over land)
     % nst_full=repmat(ncread(toread,'nst_att_status')',4,1);              % attitude flag (not used over land) old version

     varID=netcdf.inqVarID(ncid, 'sp_inc_angle')  ;
     theta=double(netcdf.getVar(ncid,varID)) ;                             % Incidence angle
     % incidenceAngleDeg=ncread(toread,'sp_inc_angle');                                % Incidence angle old version 

     varID=netcdf.inqVarID(ncid, 'sp_az_orbit')  ;
     phi_Initial_sp_az_orbit=double(netcdf.getVar(ncid,varID)) ;           % Azimuth angle
     %phi_Initial_sp_az_orbit=ncread(toread,'sp_az_orbit');                % Azimuth angle old version 
     
     varID=netcdf.inqVarID(ncid, 'sp_rx_gain')  ;
     sp_rx_gain=double(netcdf.getVar(ncid,varID)) ;                        % receiver gain [dB]
     %gain=ncread(toread,'sp_rx_gain');                                    % gain [dB] old version 
     
     varID=netcdf.inqVarID(ncid, 'gps_eirp')  ;
     eirp=double(netcdf.getVar(ncid,varID)) ;                              % eirp [dB]
     %eirp=ncread(toread,'gps_eirp');                                      % gps eirp old version         

     varID=netcdf.inqVarID(ncid, 'ddm_snr')  ;
     snr=double(netcdf.getVar(ncid,varID)) ;                               % snr of reflected signal [dB]
     %snr=ncread(toread,'ddm_snr');                                        % snr of reflected signal [dB] old version        

     varID=netcdf.inqVarID(ncid, 'ddm_noise_floor')  ;
     nf=double(netcdf.getVar(ncid,varID)) ;                                % noise floor from uncalibrated DDM of count
     %nf=ncread(toread,'ddm_noise_floor');                                 % noise floor from uncalibrated DDM of counts old version  

     varID=netcdf.inqVarID(ncid, 'rx_to_sp_range')  ;
     rxrange=double(netcdf.getVar(ncid,varID)) ;                           % rx range
     %rxrange=ncread(toread,'rx_to_sp_range');                             % rx range old version 

     varID=netcdf.inqVarID(ncid, 'tx_to_sp_range')  ;
     txrange=double(netcdf.getVar(ncid,varID)) ;                           % tx range
     %txrange=ncread(toread,'tx_to_sp_range');                             % tx range old version 
     
     varID=netcdf.inqVarID(ncid, 'ddm_nbrcs')  ;
     ddm_nbrcs=double(netcdf.getVar(ncid,varID)) ;                         % Normalized BRCS
     %ddm_nbrcs=ncread(toread,'ddm_nbrcs');                                % Normalized BRCS old version 

     varID=netcdf.inqVarID(ncid, 'brcs')  ;
     brcs=double(netcdf.getVar(ncid,varID)) ;                              % added by Hamed to save full ddm
     % brcs=ncread(toread,'brcs');                                         % added by Hamed to save full ddm old version

     varID=netcdf.inqVarID(ncid, 'reflectivity_peak')  ;
     reflectivity_peak=double(netcdf.getVar(ncid,varID)) ;                 % Peak linear reflectivityc

     varID=netcdf.inqVarID(ncid, 'ddm_ant')  ;
     receivingantenna=double(netcdf.getVar(ncid,varID)) ;                 % Reading receiving Antenna pararmeters

     disp('% reading quality flags ')

     varID=netcdf.inqVarID(ncid, 'quality_flags')  ;
     qc=netcdf.getVar(ncid,varID) ;                                        % quality flag bits
     % qualityControlFlags=ncread(toread,'quality_flags');                                  % quality flag bits old version
     
     varID=netcdf.inqVarID(ncid, 'quality_flags_2')  ;
     qc_2=netcdf.getVar(ncid,varID) ;                                        % quality flag bits number 2

     disp('% reading analog power')

     varID=netcdf.inqVarID(ncid, 'power_analog')  ;
     pa=double(netcdf.getVar(ncid,varID)) ;                                % power analog is the DDM in Watt unit
     % pa=ncread(toread,'power_analog');

     disp('% DDM shape parameters')

     varID=netcdf.inqVarID(ncid, 'coherency_ratio')  ;
     coherencyratio=double(netcdf.getVar(ncid,varID)) ;                   % Estimation of the ratio of received power between the central bins and periphery bins of the raw_counts DDM after the elimination of noise bins [4]. A higher ratio is more indicative of signal coherence.'

     varID=netcdf.inqVarID(ncid, 'ddm_les')  ;
     ddm_les=double(netcdf.getVar(ncid,varID)) ;                           % Leading edge slope

     lf=(bitget(qc,11));                                                   % SP over land flag
     maxpa=squeeze(max(pa(1:size(pa,1),1:size(pa,2),:,:)));
     peak=squeeze(max(maxpa(1:size(pa,2),:,:)));                           % peak of each DDM (tested vs for loop)
     
     % forcing to NaN data out of land
     if strcmp(savespace,'yes') | strcmp(savespace,'Yes')
             pos=find(lf>0); %lf(lf==0)=NaN;
             sp_lat=sp_lat(pos);
             sp_lon=sp_lon(pos);
             scid=scid(pos);
             strans=strans(pos);
             ts=ts(pos);
             nst_full=nst_full(pos);
             prn=prn(pos);
             azimuth_angle=azimuth_angle(pos);
             theta=theta(pos);
             phi_Initial_sp_az_orbit=phi_Initial_sp_az_orbit(pos);
             sp_rx_gain=sp_rx_gain(pos);
             eirp=eirp(pos);
             snr=snr(pos);
             nf=nf(pos);
             rxrange=rxrange(pos);
             txrange=txrange(pos);
             ddm_nbrcs=ddm_nbrcs(pos);
             qc=qc(pos);
             peak=peak(pos);
             pa=pa(:,:,pos);

             brcs=brcs(:, :, pos);                                         % added by Hamed to save full ddm

             reflectivity_peak=reflectivity_peak(pos);
             receivingantenna=receivingantenna(pos);
             qc_2=qc_2(pos);
             coherencyratio=coherencyratio(pos);
             ddm_les=ddm_les(pos);


     else
             sp_lat=sp_lat(:);
             sp_lon=sp_lon(:);
             scid=scid(:);
             strans=strans(:);
             ts=ts(:);
             nst_full=nst_full(:);
             prn=prn(:);
             azimuth_angle=azimuth_angle(:);
             theta=theta(:);
             phi_Initial_sp_az_orbit=phi_Initial_sp_az_orbit(:);
             sp_rx_gain=sp_rx_gain(:);
             eirp=eirp(:);
             snr=snr(:);
             nf=nf(:);
             rxrange=rxrange(:);
             txrange=txrange(:);
             ddm_nbrcs=ddm_nbrcs(:);
             qc=qc(:);
             peak=peak(:);
             pa=reshape(pa(:,:,:,:),size(pa,1),size(pa,2),size(pa,3)*size(pa,4));

             brcs=brcs(:, :, :);                                         % added by Hamed to save full ddm
             reflectivity_peak=reflectivity_peak(:);
             receivingantenna=receivingantenna(:);
             qc_2=qc_2(:);
             coherencyratio=coherencyratio(:);
             ddm_les=ddm_les(:);


    end
     
     % Computing Reflectivity, Kurtosis and Kurtosis zero-Doppler
     
     disp('% computing Reflectivity')
     sp_rx_gain_linear=10.^(sp_rx_gain/10);
     Reflectivity_linear=(((4.*pi).^2).*peak.*((rxrange+txrange).^2))./(eirp.*sp_rx_gain_linear.*lambda.^2);
%      disp('% Done ')    
     
     disp('% computing Kurtosis ')
     pareshaped=reshape(pa(:,:,:),size(pa,1)*size(pa,2),size(pa,3));
     Kurtosis=squeeze(kurtosis(pareshaped));
     
     disp('% computing Kurtosis zero-doppler')
     zero_Dopp = squeeze(pa(6,:,:));
     Doppler_bins=repmat(Doppler_bins',1,size(pa,3));
     zero_Dopp((zero_Dopp <0))=0 ; 
     Dopp_norm = zero_Dopp./sum(zero_Dopp);
     c=sum(Doppler_bins.*Dopp_norm); 
     var=sum(Dopp_norm.*((Doppler_bins-c).^2)) ;
     Kurtosis_dopp0=squeeze(sum(Dopp_norm.*(power(Doppler_bins-c,4))./power(var,2)));

     %      
%      disp('% Done ')  
end