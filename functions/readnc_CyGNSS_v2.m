%%%%%%%%%%%%%%%%%%%%%%%%%% function readnc_CyGNSS_v2 %%%%%%%%%%%%%%%%%%%%%%
% just_scale_factor version: reads ONLY the variables needed to compute and
% locate the scale factor (SS_r), to keep the extraction fast and light.
% Heavy DDM arrays (power_analog, raw_counts, brcs) and all derived products
% (reflectivity, kurtosis, trailing edge, coherence, calibration) are skipped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mission, L1b_product, L1b_product_version, sp_lat, sp_lon, scid, sv_num, ts, prn, SS_r] = ...
    readnc_CyGNSS_v2(inpath, filename, savespace)

     % Reading data
     toread=[inpath,filename];
     ncid = netcdf.open(toread, 'NC_NOWRITE');

     % Read global attributes
     mission=ncreadatt(toread,'/',"project");
     L1b_product=ncreadatt(toread,'/',"title");
     L1b_product_version=ncreadatt(toread,'/',"l1_algorithm_version");

     % Lat / Lon
     varID=netcdf.inqVarID(ncid, 'sp_lat')  ;
     sp_lat= (netcdf.getVar(ncid,varID)) ;                                 % Latitude
     varID=netcdf.inqVarID(ncid, 'sp_lon')  ;
     sp_lon= (netcdf.getVar(ncid,varID)) ;
     sp_lon = rem((sp_lon+180),360)-180;                                   % Longitude

     % Spacecraft / PRN / time
     varID=netcdf.inqVarID(ncid, 'spacecraft_num')  ;
     scid=netcdf.getVar(ncid,varID) ;
     scid=int8(ones(size(sp_lat)))*scid ;                                  % receiving spacecraft id

     varID=netcdf.inqVarID(ncid, 'sv_num')  ;
     sv_num= (netcdf.getVar(ncid,varID)) ;                                 % transmitting space vehicle number

     varID=netcdf.inqVarID(ncid, 'prn_code')  ;
     prn= (netcdf.getVar(ncid,varID)) ;                                    % PRN code

     varID=netcdf.inqVarID(ncid, 'ddm_timestamp_utc')  ;
     ts=netcdf.getVar(ncid,varID) ;
     ts=repmat(ts',4,1);                                                   % seconds since reference epoch

     % Ranges needed for the scale factor
     varID=netcdf.inqVarID(ncid, 'rx_to_sp_range')  ;
     rxrange= (netcdf.getVar(ncid,varID)) ;                                % rx range
     varID=netcdf.inqVarID(ncid, 'tx_to_sp_range')  ;
     txrange= (netcdf.getVar(ncid,varID)) ;                                % tx range

     % Quality flags (only for the over-land filter)
     varID=netcdf.inqVarID(ncid, 'quality_flags')  ;
     qc=netcdf.getVar(ncid,varID) ;
     lf=(bitget(qc,11));                                                   % SP over land flag

     netcdf.close(ncid);

     % forcing to land-only when requested (significantly reduces output size)
     if strcmpi(savespace,'yes')
             pos=find(lf>0);
             sp_lat=sp_lat(pos);
             sp_lon=sp_lon(pos);
             scid=scid(pos);
             sv_num=sv_num(pos);
             ts=ts(pos);
             prn=prn(pos);
             rxrange=rxrange(pos);
             txrange=txrange(pos);
     else
             sp_lat=sp_lat(:);
             sp_lon=sp_lon(:);
             scid=scid(:);
             sv_num=sv_num(:);
             ts=ts(:);
             prn=prn(:);
             rxrange=rxrange(:);
             txrange=txrange(:);
     end

     % Computing the scale factor (SS_r) from tx/rx ranges
     Aeff=189.6.*10.^6; %ricavata con la formula inversa 189.6*10^6 m^2
     txrange_d=double(txrange);
     rxrange_d=double(rxrange);
     SS_r=(Aeff.*((txrange_d+rxrange_d).^2))./(4.*pi.*(txrange_d.^2).*(rxrange_d.^2));
end
