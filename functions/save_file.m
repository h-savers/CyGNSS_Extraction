%%%%%%%%%%%%%%%%%%%%%%%%%% function readCyGNSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads the data from the CyGNSS_extract.m function and
% extract all the defined variables as Matlab or NetCDF files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_CyG_data(out_format, CyGoutpath, project_name, daterangechar, initdate, ...
    enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
    UTC_Time, DoY, SoD, receivingSpacecraft, transmittingSpacecraft, ...
    pseudoRandomNoise, spAzimuthAngleDegOrbit,specularPointLat, specularPointLon, incidenceAngleDeg, ...
    rxAntennaGain, EIRP, SNR_L1_L, PHI_Initial_sp_az_orbit, reflectivityLinear_L1_L, ...
    KURTOSIS, KURTOSIS_DOPP_0, TE_WIDTH, DDM_NBRCS, powerAnalogW, qualityControlFlags, ...
    noiseFloorCounts, REFLECTIVITY_PEAK_L1_L, receivingAntenna, qualityControlFlags_2, ...
    coherencyRatio, DDM_LES, powerRatio, notToBeUsed, notRecommended)

disp('% Saving aggregated data in a single output file');

if strcmpi(out_format,"Matlab")
    save(fullfile(CyGoutpath, [project_name '_' daterangechar '.mat']), ...
        'UTC_Time', 'receivingSpacecraft', 'transmittingSpacecraft', ...
        'pseudoRandomNoise', 'spAzimuthAngleDegOrbit', 'specularPointLat', 'specularPointLon', ...
        'incidenceAngleDeg', 'rxAntennaGain', 'EIRP', 'SNR_L1_L', ...
        'PHI_Initial_sp_az_orbit', 'reflectivityLinear_L1_L', 'KURTOSIS', ...
        'KURTOSIS_DOPP_0', 'TE_WIDTH', 'DDM_NBRCS','powerAnalogW', ...
        'qualityControlFlags', 'noiseFloorCounts','REFLECTIVITY_PEAK_L1_L', ...
        'receivingAntenna', 'qualityControlFlags_2', 'coherencyRatio', ...
        'DDM_LES', 'powerRatio', 'notToBeUsed', 'notRecommended','-v7.3');
    
elseif strcmpi(out_format,"netcdf")
    % Create NetCDF file
        netcdf_cyg = netcdf.create([CyGoutpath, '/', project_name, '_' daterangechar, '.nc'],'netcdf4');
        
        % Define the dimension of the NetCDF file
        dimid = netcdf.defDim(netcdf_cyg, 'record', length(DoY));

        % Assign global attributes
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'File name', [CyGoutpath, '/', project_name, '_' daterangechar, '.nc']);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Processing time [hh:mm:ss]', [string(s)]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Mission', 'CYGNSS');
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b product', 'NaN');
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b prodcut version', 'NaN');
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Initial dacy', [initdate]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Final day', [enddate]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Data coverage', [data_coverage]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Geographical area', [num2str(LatMin) ' - ' num2str(LatMax) ' / ' num2str(LonMin) ' - ' num2str(LonMax)]);

        % Define all variables
        var_UTC_time = netcdf.defVar(netcdf_cyg,'UTC_Time','NC_STRING',dimid);
        netcdf.putAtt(netcdf_cyg, var_UTC_time, 'long_name', 'Time of observation in UTC');
        
        var_receivingSpacecraft = netcdf.defVar(netcdf_cyg,'receivingSpacecraft','NC_SHORT',dimid);
        netcdf.putAtt(netcdf_cyg, var_receivingSpacecraft, 'long_name', 'ID of the receiving spacecraft');
        
        var_transmittingSpacecraft = netcdf.defVar(netcdf_cyg,'transmittingSpacecraft','NC_SHORT',dimid);
        netcdf.putAtt(netcdf_cyg, var_transmittingSpacecraft, 'long_name', 'ID of the transmitting spacecraft [#]');
        
        var_pseudoRandomNoise = netcdf.defVar(netcdf_cyg,'pseudoRandomNoise','NC_SHORT',dimid);
        netcdf.putAtt(netcdf_cyg, var_pseudoRandomNoise, 'long_name', 'PRN code [#]');
        
        var_spAzimuthAngleDegOrbit = netcdf.defVar(netcdf_cyg,'spAzimuthAngleDegOrbit','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_spAzimuthAngleDegOrbit, 'long_name', 'Azimuth angle at specular point in the orbit frame [deg]');
        
        var_specularPointLat = netcdf.defVar(netcdf_cyg,'specularPointLat','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLat, 'long_name', 'Latitude of observation [deg]');
        
        var_specularPointLon = netcdf.defVar(netcdf_cyg,'specularPointLon','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLon, 'long_name', 'Longitude of observation [deg]');
        
        var_incidenceAngleDeg = netcdf.defVar(netcdf_cyg,'incidenceAngleDeg','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_incidenceAngleDeg, 'long_name', 'Incidence angle at specular point [deg]');
        
        var_rxAntennaGain = netcdf.defVar(netcdf_cyg,'rxAntennaGain','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_rxAntennaGain, 'long_name', 'Receiver antenna gain toward specular point [dB]');
        
        var_EIRP = netcdf.defVar(netcdf_cyg,'EIRP','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_EIRP, 'long_name', 'Effective Isotropic Radiated Power [Watt]');
        
        var_SNR_L1_L = netcdf.defVar(netcdf_cyg,'SNR_L1_L','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_SNR_L1_L, 'long_name', 'Signal to noise ratio of signal L1 and left polarization [dB]');
        
        var_reflectivityLinear_L1_L = netcdf.defVar(netcdf_cyg,'reflectivityLinear_L1_L','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_reflectivityLinear_L1_L, 'long_name', 'Reflection coefficient of GPS signal L1, left polarization [dB]');
        
        var_qualityControlFlags = netcdf.defVar(netcdf_cyg,'qualityControlFlags','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_qualityControlFlags, 'long_name', 'Quality check flags of the observation [int32]');
        
        var_noiseFloorCounts = netcdf.defVar(netcdf_cyg,'noiseFloorCounts','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_noiseFloorCounts, 'long_name', 'Noise floor of the DDM [Watt]');
        
        var_recevingAntenna = netcdf.defVar(netcdf_cyg,'receivingAntenna','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_recevingAntenna, 'long_name', ['Antenna collecting the signal. In CyGNSS it can be ' ...
            '0=none, 1=never used, 2=nadir starboard, 3=nadir port [integer]']);
        
        var_coherencyRatio = netcdf.defVar(netcdf_cyg,'coherencyRatio','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_coherencyRatio, 'long_name', ['Estimation of the ratio of received power between the central bins ' ...
            'and periphery bins of the raw_counts DDM after the elimination of noise bins']);
        
        var_powerRatio = netcdf.defVar(netcdf_cyg,'powerRatio','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_powerRatio, 'long_name', ['Estimation of the degree of coherency based on the algorithm proposed in ' ...
            'https://doi.org/10.1109/TGRS.2020.3009784 applied to ddm raw counts']);
        
        var_notToBeUsed = netcdf.defVar(netcdf_cyg,'notToBeUsed','NC_BYTE',dimid);
        netcdf.putAtt(netcdf_cyg, var_notToBeUsed, 'long_name', ['Quality flag based on L1b flags indicating reflection not to be used ' ...
            '(1 indicates bad quality)']);
        
        var_notRecommended = netcdf.defVar(netcdf_cyg,'notRecommended','NC_BYTE',dimid);
        netcdf.putAtt(netcdf_cyg, var_notRecommended, 'long_name', ['Quality flag based on L1B variables beyond thresholds defined in ' ...
            'the configuration file. When 1 (i.e., up) reflection is suspicious']);
        
        % Variable NOT used, not in the Excel
        var_PHI_Initial_sp_az_orbit = netcdf.defVar(netcdf_cyg,'PHI_Initial_sp_az_orbit','NC_DOUBLE',dimid);
        var_KURTOSIS = netcdf.defVar(netcdf_cyg,'KURTOSIS','NC_DOUBLE',dimid);
        var_KURTOSIS_DOPP_0 = netcdf.defVar(netcdf_cyg,'KURTOSIS_DOPP_0','NC_DOUBLE',dimid);
        var_TE_WIDTH = netcdf.defVar(netcdf_cyg,'TE_WIDTH','NC_DOUBLE',dimid);
        var_DDM_NBRCS = netcdf.defVar(netcdf_cyg,'DDM_NBRCS','NC_DOUBLE',dimid);
        var_PA_L1_L = netcdf.defVar(netcdf_cyg,'powerAnalogW','NC_DOUBLE',dimid);
        var_REFLECTIVITY_PEAK_L1_L = netcdf.defVar(netcdf_cyg,'REFLECTIVITY_PEAK_L1_L','NC_DOUBLE',dimid);
        var_qualityControlFlags_2 = netcdf.defVar(netcdf_cyg,'qualityControlFlags_2','NC_UINT',dimid);
        var_DDM_LES = netcdf.defVar(netcdf_cyg,'DDM_LES','NC_DOUBLE',dimid);

        % End definition mode
        netcdf.endDef(netcdf_cyg);
    
        % Write data to variables
        UTC_Time_str = string(UTC_Time);
        netcdf.putVar(netcdf_cyg, var_UTC_time, UTC_Time_str);
        netcdf.putVar(netcdf_cyg, var_receivingSpacecraft, receivingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_transmittingSpacecraft, transmittingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_pseudoRandomNoise, pseudoRandomNoise);
        netcdf.putVar(netcdf_cyg, var_spAzimuthAngleDegOrbit, spAzimuthAngleDegOrbit);
        netcdf.putVar(netcdf_cyg, var_specularPointLat, specularPointLat);
        netcdf.putVar(netcdf_cyg, var_specularPointLon, specularPointLon);
        netcdf.putVar(netcdf_cyg, var_incidenceAngleDeg, incidenceAngleDeg);
        netcdf.putVar(netcdf_cyg, var_rxAntennaGain, rxAntennaGain);
        netcdf.putVar(netcdf_cyg, var_EIRP, EIRP);
        netcdf.putVar(netcdf_cyg, var_SNR_L1_L, SNR_L1_L);
        netcdf.putVar(netcdf_cyg, var_reflectivityLinear_L1_L, reflectivityLinear_L1_L);
        netcdf.putVar(netcdf_cyg, var_qualityControlFlags, qualityControlFlags);
        netcdf.putVar(netcdf_cyg, var_noiseFloorCounts, noiseFloorCounts);
        netcdf.putVar(netcdf_cyg, var_recevingAntenna, receivingAntenna);
        netcdf.putVar(netcdf_cyg, var_coherencyRatio, coherencyRatio);
        netcdf.putVar(netcdf_cyg, var_powerRatio, powerRatio);
        netcdf.putVar(netcdf_cyg, var_notToBeUsed, uint8(notToBeUsed));
        netcdf.putVar(netcdf_cyg, var_notRecommended, uint8(notRecommended));
        % netcdf.putVar(netcdf_cyg, var_PHI_Initial_sp_az_orbit, PHI_Initial_sp_az_orbit);
        % netcdf.putVar(netcdf_cyg, var_KURTOSIS, KURTOSIS);
        % netcdf.putVar(netcdf_cyg, var_KURTOSIS_DOPP_0, KURTOSIS_DOPP_0);
        % netcdf.putVar(netcdf_cyg, var_TE_WIDTH, TE_WIDTH);
        % netcdf.putVar(netcdf_cyg, var_DDM_NBRCS, DDM_NBRCS);
        % netcdf.putVar(netcdf_cyg, var_PA_L1_L, powerAnalogW);
        % netcdf.putVar(netcdf_cyg, var_REFLECTIVITY_PEAK_L1_L, REFLECTIVITY_PEAK_L1_L);
        % netcdf.putVar(netcdf_cyg, var_qualityControlFlags_2, qualityControlFlags_2);
        % netcdf.putVar(netcdf_cyg, var_DDM_LES, DDM_LES);
end
