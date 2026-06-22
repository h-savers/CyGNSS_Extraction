%%%%%%%%%%%%%%%%%%%%%%%%%% function save_file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just_scale_factor version: saves ONLY the scale factor (scaleFactor) and
% the geolocation/time identifiers needed to match it, as Matlab or NetCDF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_file(mission, out_format, L1b_product, L1b_product_version,...
        CyGoutpath, project_name, daterangechar, initdate, ...
        enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
        timeUTC, receivingSpacecraft, transmittingSpacecraft, ...
        pseudoRandomNoise, specularPointLat, specularPointLon, scaleFactor)

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
timestamp_global_attr = datestr(now, 'yyyy-mm-dd HH:MM:SS');

constellation = repmat("GPS", length(timeUTC), 1);  % column array of strings

if strcmpi(out_format,"Matlab")
      save(fullfile(CyGoutpath, [project_name '_' daterangechar '_' timestamp '.mat']), ...
        'timeUTC', 'receivingSpacecraft', 'transmittingSpacecraft', 'constellation',...
        'pseudoRandomNoise', 'specularPointLat', 'specularPointLon', 'scaleFactor', '-v7.3');

elseif strcmpi(out_format,"netcdf")
    % Create NetCDF file
        netcdf_cyg = netcdf.create([CyGoutpath, '/', project_name, '_' daterangechar, '_' timestamp, '.nc'],'netcdf4');

        % Define the dimension of the NetCDF file
        dimid = netcdf.defDim(netcdf_cyg, 'record', length(timeUTC));

        % Assign global attributes
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'File name', [project_name '_' daterangechar '_' timestamp '.nc']);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'File generation time [yyyy-mm-dd HH:MM:SS]', string(timestamp_global_attr));
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Mission', mission);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b product', L1b_product);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b product version', L1b_product_version);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Initial day', initdate);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Final day', enddate);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Data coverage', data_coverage);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Geographical area', [num2str(LatMin) ' - ' num2str(LatMax) ' / ' num2str(LonMin) ' - ' num2str(LonMax)]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Processing time', char(s));

        % Define variables
        var_timeUTC = netcdf.defVar(netcdf_cyg,'timeUTC','NC_STRING',dimid);
        netcdf.putAtt(netcdf_cyg, var_timeUTC, 'timeUTC', 'Time of observation in UTC');

        var_specularPointLat = netcdf.defVar(netcdf_cyg,'specularPointLat','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLat, 'specularPointLat', 'Latitude of observation [deg]');

        var_specularPointLon = netcdf.defVar(netcdf_cyg,'specularPointLon','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLon, 'specularPointLon', 'Longitude of observation [deg]');

        var_constellation = netcdf.defVar(netcdf_cyg,'constellation','NC_STRING',dimid);
        netcdf.putAtt(netcdf_cyg, var_constellation, 'constellation', 'Name of the constellation [ASCII]');

        var_receivingSpacecraft = netcdf.defVar(netcdf_cyg,'receivingSpacecraft','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_receivingSpacecraft, 'receivingSpacecraft', 'ID of the receiving spacecraft');

        var_transmittingSpacecraft = netcdf.defVar(netcdf_cyg,'transmittingSpacecraft','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_transmittingSpacecraft, 'transmittingSpacecraft', 'ID of the transmitting spacecraft [#]');

        var_pseudoRandomNoise = netcdf.defVar(netcdf_cyg,'pseudoRandomNoise','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_pseudoRandomNoise, 'pseudoRandomNoise', 'PRN code [#]');

        var_scaleFactor = netcdf.defVar(netcdf_cyg,'scaleFactor','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_scaleFactor, 'scaleFactor', 'Scale factor (SS_r) computed from tx/rx ranges as Aeff*(tx_to_sp_range+rx_to_sp_range)^2 / (4*pi*tx_to_sp_range^2*rx_to_sp_range^2), with Aeff=189.6e6 m^2');

        % End definition mode
        netcdf.endDef(netcdf_cyg);

        % Write data to variables
        netcdf.putVar(netcdf_cyg, var_timeUTC, string(timeUTC));
        netcdf.putVar(netcdf_cyg, var_specularPointLat, specularPointLat);
        netcdf.putVar(netcdf_cyg, var_specularPointLon, specularPointLon);
        netcdf.putVar(netcdf_cyg, var_constellation, constellation);
        netcdf.putVar(netcdf_cyg, var_receivingSpacecraft, receivingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_transmittingSpacecraft, transmittingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_pseudoRandomNoise, pseudoRandomNoise);
        netcdf.putVar(netcdf_cyg, var_scaleFactor, scaleFactor);

%% ------- Close netcdf file ------- %%
        netcdf.close(netcdf_cyg) ;
end
end
