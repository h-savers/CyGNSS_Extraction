%%%%%%%%%%%%%%%%%%%%%%%%%% function readCyGNSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads the data from the netCDF files and computes 
% Reflectivity, Kurtosis and Kurtosis zero-Doppler.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sp_lat,sp_lon,scid,ts,nst_full,prn,theta,phi_Initial_sp_az_orbit,gain, ...
    eirp,snr,nf,rxrange,txrange,ddm_nbrcs,qc,pa,Reflectivity_linear,Kurtosis,Kurtosis_dopp0, brcs]= ...
    readnc_CyGNSS(inpath,filename,lambda,Doppler_bins,savespace)

     % Reading data
     toread=[inpath,filename];
     % using ncread
     disp('% opening file ')
     disp('% reading lat/lon')
     sp_lat = ncread(toread,'sp_lat');                                     % Latitude     
     sp_lon = rem((ncread(toread,'sp_lon')+180),360)-180;                  % Longitude
                                      
     disp('% reading satellite ID and time ')
     scid=double(ncread(toread,'spacecraft_num')).*ones(size(sp_lat));     % spacecraft id
     ts=repmat(ncread(toread,'ddm_timestamp_utc')',4,1);                   % time in seconds since 2017-03-27 00:00:00.999261529
            
     disp('% reading observables ')
     nst_full=repmat(ncread(toread,'nst_att_status')',4,1);                % attitude flag (not used over land)
     prn=ncread(toread,'prn_code');                                        % full prn
     theta=ncread(toread,'sp_inc_angle');                                  % Incidence angle
     phi_Initial_sp_az_orbit=ncread(toread,'sp_az_orbit');                 % Azimuth angle
     gain=ncread(toread,'sp_rx_gain');                                     % gain [dB]
     eirp=ncread(toread,'gps_eirp');                                       % gps eirp        
     snr=ncread(toread,'ddm_snr');                                         % snr of reflected signal [dB]       
     nf=ncread(toread,'ddm_noise_floor');                                  % noise floor from uncalibrated DDM of counts 
     rxrange=ncread(toread,'rx_to_sp_range');                              % rx range
     txrange=ncread(toread,'tx_to_sp_range');                              % tx range
     ddm_nbrcs=ncread(toread,'ddm_nbrcs');                                 % Normalized BRCS
     
     disp('% reading quality flags ')
     qc=ncread(toread,'quality_flags');                                    % quality flag bits
     
     disp('% reading analog power')
     pa=ncread(toread,'power_analog');


     brcs=ncread(toread,'brcs');                                           % added by Hamed to save full ddm


     lf=(bitget(qc,11));                                                     % SP over land flag
     maxpa=squeeze(max(pa(1:size(pa,1),1:size(pa,2),:,:)));
     peak=squeeze(max(maxpa(1:size(pa,2),:,:)));                           % peak of each DDM (tested vs for loop)
     
     % forcing to NaN data out of land
     if strcmp(savespace,'yes')
             pos=find(lf>0); %lf(lf==0)=NaN;
             sp_lat=sp_lat(pos);
             sp_lon=sp_lon(pos);
             scid=scid(pos);
             ts=ts(pos);
             nst_full=nst_full(pos);
             prn=prn(pos);
             theta=theta(pos);
             phi_Initial_sp_az_orbit=phi_Initial_sp_az_orbit(pos);
             gain=gain(pos);
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

     else
             sp_lat=sp_lat(:);
             sp_lon=sp_lon(:);
             scid=scid(:);
             ts=ts(:);
             nst_full=nst_full(:);
             prn=prn(:);
             theta=theta(:);
             phi_Initial_sp_az_orbit=phi_Initial_sp_az_orbit(:);
             gain=gain(:);
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

    end
     
     % Computing Reflectivity, Kurtosis and Kurtosis zero-Doppler
     
     disp('% computing Reflectivity')
     gain_linear=10.^(gain/10);
     Reflectivity_linear=(((4.*pi).^2).*peak.*((rxrange+txrange).^2))./(eirp.*gain_linear.*lambda.^2);
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