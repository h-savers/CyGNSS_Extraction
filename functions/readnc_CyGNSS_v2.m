%%%%%%%%%%%%%%%%%%%%%%%%%% function readCyGNSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads the data from the netCDF files and computes 
% Reflectivity, Kurtosis and Kurtosis zero-Doppler.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mission, L1b_product, L1b_product_version,sp_lat,sp_lon,scid,sv_num,ts,nst_full, ...
    prn,theta,phi_Initial_sp_az_orbit,sp_rx_gain,eirp,snr,nf,rxrange,txrange,ddm_nbrcs,...
    qc,pa,Reflectivity_linear,Kurtosis, Kurtosis_dopp0, brcs, reflectivity_peak, ...
    receivingantenna, sp_azimuth_angle_deg_north, qc_2, coherency_ratio, ddm_les, ...
    raw_counts]= ...
    readnc_CyGNSS_v2(inpath,filename,lambda,Doppler_bins,savespace)

     % Reading data
     toread=[inpath,filename];
     % using ncread
     disp('% opening file ')
     ncid = netcdf.open(toread, 'NC_NOWRITE');
% trackNcids = netcdf.inqGrps(ncid);
     
     % Read global attributes
     mission=ncreadatt(toread,'/',"project");
     L1b_product=ncreadatt(toread,'/',"title");
     L1b_product_version=ncreadatt(toread,'/',"l1_algorithm_version");

     disp('% reading lat/lon')
     
     varID=netcdf.inqVarID(ncid, 'sp_lat')  ;
     sp_lat= (netcdf.getVar(ncid,varID)) ;                                 % Latitude   

     varID=netcdf.inqVarID(ncid, 'sp_lon')  ;
     sp_lon= (netcdf.getVar(ncid,varID)) ;
     sp_lon = rem((sp_lon+180),360)-180;                                   % Longitude
     % sp_lon = rem((ncread(toread,'sp_lon')+180),360)-180;                % Longitude old version
                                      
     disp('% reading satellite ID and time ')

     varID=netcdf.inqVarID(ncid, 'spacecraft_num')  ;
     scid=netcdf.getVar(ncid,varID) ; 
     scid=int8(ones(size(sp_lat)))*scid ;     % spacecraft id
     % scid=double(ncread(toread,'spacecraft_num')).*ones(size(sp_lat));   % spacecraft id ol version
     
     varID=netcdf.inqVarID(ncid, 'sv_num')  ;
     sv_num= (netcdf.getVar(ncid,varID)) ;                                    % full space vehicle number that trasmitted prn_code 

     varID=netcdf.inqVarID(ncid, 'prn_code')  ;
     prn= (netcdf.getVar(ncid,varID)) ;                                    % full prn 
     % prn=ncread(toread,'prn_code');                                      % full prn old version 
     
     varID=netcdf.inqVarID(ncid, 'ddm_timestamp_utc')  ;
     ts=netcdf.getVar(ncid,varID) ; 
     ts=repmat(ts',4,1);                                                   % time in seconds since 2017-03-27 00:00:00.999261529
     % ts=repmat(ncread(toread,'ddm_timestamp_utc')',4,1);                 % time in seconds since 2017-03-27 00:00:00.999261529 old version
            
     disp('% reading observables ')
     varID=netcdf.inqVarID(ncid, 'nst_att_status')  ;
     nst_full=netcdf.getVar(ncid,varID) ; 
     nst_full= (repmat(nst_full',4,1)) ;                                   % attitude flag (not used over land)
     % nst_full=repmat(ncread(toread,'nst_att_status')',4,1);              % attitude flag (not used over land) old version

     varID=netcdf.inqVarID(ncid, 'sp_inc_angle')  ;
     theta= (netcdf.getVar(ncid,varID)) ;                                  % Incidence angle
     % theta=ncread(toread,'sp_inc_angle');                                % Incidence angle old version 

     varID=netcdf.inqVarID(ncid, 'sp_az_orbit')  ;
     phi_Initial_sp_az_orbit= (netcdf.getVar(ncid,varID)) ;                % Azimuth angle
     %phi_Initial_sp_az_orbit=ncread(toread,'sp_az_orbit');                % Azimuth angle old version 
     
     varID=netcdf.inqVarID(ncid, 'sp_rx_gain')  ;
     sp_rx_gain= (netcdf.getVar(ncid,varID)) ;                             % receiver gain [dB]
     %gain=ncread(toread,'sp_rx_gain');                                    % gain [dB] old version 
     
     varID=netcdf.inqVarID(ncid, 'gps_eirp')  ;
     eirp= (netcdf.getVar(ncid,varID)) ;                                   % eirp [dB]
     %eirp=ncread(toread,'gps_eirp');                                      % gps eirp old version         

     varID=netcdf.inqVarID(ncid, 'ddm_snr')  ;
     snr= (netcdf.getVar(ncid,varID)) ;                                    % snr of reflected signal [dB]
     %snr=ncread(toread,'ddm_snr');                                        % snr of reflected signal [dB] old version        

     varID=netcdf.inqVarID(ncid, 'ddm_noise_floor')  ;
     nf= (netcdf.getVar(ncid,varID)) ;                                     % noise floor from uncalibrated DDM of count
     %nf=ncread(toread,'ddm_noise_floor');                                 % noise floor from uncalibrated DDM of counts old version  

     varID=netcdf.inqVarID(ncid, 'rx_to_sp_range')  ;
     rxrange= (netcdf.getVar(ncid,varID)) ;                                % rx range
     %rxrange=ncread(toread,'rx_to_sp_range');                             % rx range old version 

     varID=netcdf.inqVarID(ncid, 'tx_to_sp_range')  ;
     txrange= (netcdf.getVar(ncid,varID)) ;                                % tx range
     %txrange=ncread(toread,'tx_to_sp_range');                             % tx range old version 
     
     varID=netcdf.inqVarID(ncid, 'ddm_nbrcs')  ;
     ddm_nbrcs= (netcdf.getVar(ncid,varID)) ;                              % Normalized BRCS
     %ddm_nbrcs=ncread(toread,'ddm_nbrcs');                                % Normalized BRCS old version 

     varID=netcdf.inqVarID(ncid, 'brcs')  ;
     brcs= (netcdf.getVar(ncid,varID)) ;                                   % added by Hamed to save full ddm
     % brcs=ncread(toread,'brcs');                                         % added by Hamed to save full ddm old version

     varID=netcdf.inqVarID(ncid, 'reflectivity_peak')  ;
     reflectivity_peak= (netcdf.getVar(ncid,varID)) ;                      % Peak linear reflectivity

     varID=netcdf.inqVarID(ncid, 'ddm_ant')  ;
     receivingantenna=(netcdf.getVar(ncid,varID)) ;                 % Reading receiving Antenna pararmeters


     disp('% reading spacecraft position and velocity to compute the specular point azimuth angle ')

     varID=netcdf.inqVarID(ncid, 'sc_pos_x')  ;
     sc_pos_x=single(netcdf.getVar(ncid,varID)) ;  
     sc_pos_x = repmat(sc_pos_x',4,1);

     varID=netcdf.inqVarID(ncid, 'sc_pos_y')  ;
     sc_pos_y=single(netcdf.getVar(ncid,varID)) ;  
     sc_pos_y = repmat(sc_pos_y',4,1);

     varID=netcdf.inqVarID(ncid, 'sc_pos_z')  ;
     sc_pos_z=single(netcdf.getVar(ncid,varID)) ;  
     sc_pos_z = repmat(sc_pos_z',4,1);

     varID=netcdf.inqVarID(ncid, 'sc_vel_x')  ;
     sc_vel_x=single(netcdf.getVar(ncid,varID)) ;  
     sc_vel_x = repmat(sc_vel_x',4,1);
     
     varID=netcdf.inqVarID(ncid, 'sc_vel_y')  ;
     sc_vel_y=single(netcdf.getVar(ncid,varID)) ;  
     sc_vel_y = repmat(sc_vel_y',4,1);

     varID=netcdf.inqVarID(ncid, 'sc_vel_z')  ;
     sc_vel_z=single(netcdf.getVar(ncid,varID)) ;          
     sc_vel_z = repmat(sc_vel_z',4,1);

     disp('% reading quality flags ')

     varID=netcdf.inqVarID(ncid, 'quality_flags')  ;
     qc=netcdf.getVar(ncid,varID) ;                                        % quality flag bits
     % qc=ncread(toread,'quality_flags');                                  % quality flag bits old version
     
     varID=netcdf.inqVarID(ncid, 'quality_flags_2')  ;
     qc_2=netcdf.getVar(ncid,varID) ;                                      % quality flag bits number 2

     disp('% reading analog power')

     varID=netcdf.inqVarID(ncid, 'power_analog')  ;
     pa= (netcdf.getVar(ncid,varID)) ;                                     % power analog is the DDM in Watt unit
     % pa=ncread(toread,'power_analog');

     varID=netcdf.inqVarID(ncid, 'raw_counts')  ; 
     raw_counts= (netcdf.getVar(ncid,varID)) ;                             % DDM bin raw counts

     disp('% DDM shape parameters')

     varID=netcdf.inqVarID(ncid, 'coherency_ratio')  ;
     coherency_ratio= (netcdf.getVar(ncid,varID)) ;                        % Estimation of the ratio of received power between the central bins and periphery bins of the raw_counts DDM after the elimination of noise bins [4]. A higher ratio is more indicative of signal coherence.'

     varID=netcdf.inqVarID(ncid, 'ddm_les')  ;
     ddm_les= (netcdf.getVar(ncid,varID)) ;                                % Leading edge slope

     lf=(bitget(qc,11));                                                   % SP over land flag
     maxpa=squeeze(max(pa(1:size(pa,1),1:size(pa,2),:,:)));
     peak=squeeze(max(maxpa(1:size(pa,2),:,:)));                           % peak of each DDM (tested vs for loop)
     
     % forcing to NaN data out of land
     if strcmp(savespace,'yes') | strcmp(savespace,'Yes')
             pos=find(lf>0); %lf(lf==0)=NaN;
             sp_lat=sp_lat(pos);
             sp_lon=sp_lon(pos);
             scid=scid(pos);
             sv_num=sv_num(pos);
             ts=ts(pos);
             nst_full=nst_full(pos);
             prn=prn(pos);
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
             raw_counts=raw_counts(:,:,pos) ;
             brcs=brcs(:, :, pos);                                         % added by Hamed to save full ddm
             reflectivity_peak=reflectivity_peak(pos);
             sc_pos_x=sc_pos_x(pos);
             sc_pos_y=sc_pos_y(pos);
             sc_pos_z=sc_pos_z(pos);
             sc_vel_x=sc_vel_x(pos);
             sc_vel_y=sc_vel_y(pos);
             sc_vel_z=sc_vel_z(pos);
             receivingantenna=receivingantenna(pos);
             qc_2=qc_2(pos);
             coherency_ratio=coherency_ratio(pos);
             ddm_les=ddm_les(pos);


     else
             sp_lat=sp_lat(:);
             sp_lon=sp_lon(:);
             scid=scid(:);
             sv_num=sv_num(:);
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
             raw_counts=reshape(raw_counts(:,:,:,:),size(raw_counts,1),size(raw_counts,2),size(raw_counts,3)*size(raw_counts,4));


             brcs=brcs(:, :, :);                                           % added by Hamed to save full ddm
             reflectivity_peak=reflectivity_peak(:);
             sc_pos_x=sc_pos_x(:);
             sc_pos_y=sc_pos_y(:);
             sc_pos_z=sc_pos_z(:);
             sc_vel_x=sc_vel_x(:);
             sc_vel_y=sc_vel_y(:);
             sc_vel_z=sc_vel_z(:);
             receivingantenna=receivingantenna(:);
             qc_2=qc_2(:);
             coherency_ratio=coherency_ratio(:);
             ddm_les=ddm_les(:);


    end
     
     % Computing Reflectivity, Kurtosis and Kurtosis zero-Doppler
     disp('% computing the azimuth angle of the specular point')
     [sp_azimuth_angle_deg_north]=compute_azimuth_angle(sc_pos_x,sc_pos_y,sc_pos_z,sc_vel_x,sc_vel_y,sc_vel_z,phi_Initial_sp_az_orbit);

     disp('% computing Reflectivity')
     sp_rx_gain_linear=10.^(sp_rx_gain/10);
     Reflectivity_linear=(((4.*pi).^2).*peak.*((single(rxrange)+single(txrange)).^2))./(eirp.*sp_rx_gain_linear.*lambda.^2);
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
%    disp('% Done ')  
end