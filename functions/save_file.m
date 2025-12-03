%%%%%%%%%%%%%%%%%%%%%%%%%% function readCyGNSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads the data from the CyGNSS_extract.m function and
% extract all the defined variables as Matlab or NetCDF files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_file(mission, out_format, L1b_product, L1b_product_version,...
        CyGoutpath, project_name, daterangechar, initdate, ...
        enddate, LonMin, LonMax, LatMin, LatMax, data_coverage, s, ...
        timeUTC, receivingSpacecraft, transmittingSpacecraft, ...
        pseudoRandomNoise, spAzimuthAngleDegOrbit,specularPointLat, specularPointLon, incidenceAngleDeg, ...
        rxAntennaGain_L1_L, EIRP_L1, SNR_L1_L, reflectivityLinear_L1_L, ...
        kurtosisDDM, kurtosisDopp0, teWidth, NBRCS_L1_L, powerAnalogW_L1_L, qualityControlFlags_L1_L, ...
        noiseFloorCounts_L1_L, reflectivityPeak_L1_L, reflectivityPeakRecal_L1_L, receivingAntenna, qualityControlFlags_2_L1_L, bitRatio, ...
        spAzimuthAngleDegNorth, coherencyRatio_L1_L, ddmLes, powerRatio_L1_L, notToBeUsed, notRecommended, coefficientOfVariation);

powerRatio_L1_L=single(powerRatio_L1_L) ; 

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
timestamp_global_attr = datestr(now, 'yyyy-mm-dd HH:MM:SS');

constellation = repmat("GPS", length(timeUTC), 1);  % crea un array colonna di stringhe

if strcmpi(out_format,"Matlab")
      save(fullfile(CyGoutpath, [project_name '_' daterangechar '_' timestamp '.mat']), ...
        'timeUTC', 'receivingSpacecraft', 'transmittingSpacecraft', 'constellation',...
        'pseudoRandomNoise', 'spAzimuthAngleDegOrbit', 'specularPointLat', 'specularPointLon', ...
        'incidenceAngleDeg', 'rxAntennaGain_L1_L', 'EIRP_L1', 'SNR_L1_L', ...
        'reflectivityLinear_L1_L', 'kurtosisDDM', ...
        'kurtosisDopp0', 'teWidth', 'NBRCS_L1_L','powerAnalogW_L1_L', ...
        'qualityControlFlags_L1_L', 'noiseFloorCounts_L1_L','reflectivityPeak_L1_L', 'reflectivityPeakRecal_L1_L',...
        'receivingAntenna', 'qualityControlFlags_2_L1_L', 'bitRatio', 'spAzimuthAngleDegNorth', 'coherencyRatio_L1_L', ...
        'ddmLes', 'powerRatio_L1_L', 'notToBeUsed', 'notRecommended','coefficientOfVariation','-v7.3');
    
elseif strcmpi(out_format,"netcdf")
    % Create NetCDF file
        netcdf_cyg = netcdf.create([CyGoutpath, '/', project_name, '_' daterangechar, '_' timestamp, '.nc'],'netcdf4');
        
        % Define the dimension of the NetCDF file
        dimid = netcdf.defDim(netcdf_cyg, 'record', length(timeUTC));
    
        % Assign global attributes
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'File name', [project_name '_' daterangechar '_' timestamp '.nc']);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'File generation time [yyyy-mm-dd HH:MM:SS]', [string(timestamp_global_attr)]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Mission', [mission]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b product', [L1b_product]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'L1b product version', [L1b_product_version]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Initial day', [initdate]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Final day', [enddate]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Data coverage', [data_coverage]);
        netcdf.putAtt(netcdf_cyg, netcdf.getConstant("NC_GLOBAL"), 'Geographical area', [num2str(LatMin) ' - ' num2str(LatMax) ' / ' num2str(LonMin) ' - ' num2str(LonMax)]);

        % Define all variables
        var_timeUTC = netcdf.defVar(netcdf_cyg,'timeUTC','NC_STRING',dimid);
        netcdf.putAtt(netcdf_cyg, var_timeUTC, 'timeUTC', 'Time of observation in UTC');
        
        var_specularPointLat = netcdf.defVar(netcdf_cyg,'specularPointLat','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLat, 'specularPointLat', 'Latitude of observation [deg]');

        var_specularPointLon = netcdf.defVar(netcdf_cyg,'specularPointLon','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_specularPointLon, 'specularPointLon', 'Longitude of observation [deg]');
        
        var_incidenceAngleDeg = netcdf.defVar(netcdf_cyg,'incidenceAngleDeg','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_incidenceAngleDeg, 'incidenceAngleDeg', 'Incidence angle at specular point [deg]');

        var_spAzimuthAngleDegNorth = netcdf.defVar(netcdf_cyg,'spAzimuthAngleDegNorth','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_spAzimuthAngleDegNorth, 'spAzimuthAngleDegNorth', 'Azimuth angle at specular point with respect to North [deg]');

        var_SNR_L1_L = netcdf.defVar(netcdf_cyg,'SNR_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_SNR_L1_L, 'SNR_L1_L', 'Signal to noise ratio of signal L1 and left polarization [dB]');

        var_powerRatio_L1_L = netcdf.defVar(netcdf_cyg,'powerRatio_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_powerRatio_L1_L, 'powerRatio_L1_L', ['Estimation of the degree of coherency based on the algorithm proposed in ' ...
            'https://doi.org/10.1109/TGRS.2020.3009784 applied to ddm raw counts']);

        var_coherencyRatio_L1_L = netcdf.defVar(netcdf_cyg,'coherencyRatio_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_coherencyRatio_L1_L, 'long_name', ['Estimation of the ratio of received power between the central bins ' ...
            'and periphery bins of the raw_counts DDM after the elimination of noise bins']);
        
        var_reflectivityPeak_L1_L = netcdf.defVar(netcdf_cyg,'reflectivityPeak_L1_L','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_reflectivityPeak_L1_L, 'reflectivityPeak_L1_L', 'Reflection coefficient of GPS signal L1, left polarization [dB] read from L1b CyGNSS data');

        var_reflectivityPeakRecal_L1_L = netcdf.defVar(netcdf_cyg,'reflectivityPeakRecal_L1_L','NC_DOUBLE',dimid);
        netcdf.putAtt(netcdf_cyg, var_reflectivityPeakRecal_L1_L, 'reflectivityPeakRecal_L1_L', 'Reflection coefficient of GPS signal L1, left polarization [dB] re-calibrated');

        var_NBRCS_L1_L = netcdf.defVar(netcdf_cyg,'NBRCS_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_NBRCS_L1_L, 'NBRCS_L1_L', 'Normalized bistatic radar cross section of GPS signal L1, left polarization [dB]');

        var_qualityControlFlags_L1_L = netcdf.defVar(netcdf_cyg,'qualityControlFlags_L1_L','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_qualityControlFlags_L1_L, 'qualityControlFlags_L1_L', 'First quality check flags of the observation [int32]');

        var_qualityControlFlags_2_L1_L = netcdf.defVar(netcdf_cyg,'qualityControlFlags_2_L1_L','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_qualityControlFlags_2_L1_L, 'qualityControlFlags_2_L1_L', 'Second quality check flags of the observation [int32]');

        var_notToBeUsed = netcdf.defVar(netcdf_cyg,'notToBeUsed','NC_BYTE',dimid);
        netcdf.putAtt(netcdf_cyg, var_notToBeUsed, 'notToBeUsed', ['Quality flag based on L1b flags indicating reflection not to be used (1 indicates bad quality).']);
        
        var_notRecommended = netcdf.defVar(netcdf_cyg,'notRecommended','NC_BYTE',dimid);
        netcdf.putAtt(netcdf_cyg, var_notRecommended, 'notRecommended', ['Quality flag based on L1B variables that are beyond thresholds defined in the configuration file. When it is 1 (i.e., flag is up) reflection is suspicious.']);

        var_uncertainty = netcdf.defVar(netcdf_cyg,'errorUncertainty','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_uncertainty, 'errorUncertainty', ['Coefficient of variation of observation estimated from SNR and incoherent integration']);

        var_constellation = netcdf.defVar(netcdf_cyg,'constellation','NC_STRING',dimid);
        netcdf.putAtt(netcdf_cyg, var_constellation, 'constellation', 'Name of the constallation [ASCII]');

        var_recevingAntenna = netcdf.defVar(netcdf_cyg,'receivingAntenna','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_recevingAntenna, 'receivingAntenna', ['Antenna collecting the signal']);

        var_receivingSpacecraft = netcdf.defVar(netcdf_cyg,'receivingSpacecraft','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_receivingSpacecraft, 'receivingSpacecraft', 'ID of the receiving spacecraft');

        var_transmittingSpacecraft = netcdf.defVar(netcdf_cyg,'transmittingSpacecraft','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_transmittingSpacecraft, 'transmittingSpacecraft', 'ID of the transmitting spacecraft [#]');

        var_pseudoRandomNoise = netcdf.defVar(netcdf_cyg,'pseudoRandomNoise','NC_INT',dimid);
        netcdf.putAtt(netcdf_cyg, var_pseudoRandomNoise, 'pseudoRandomNoise', 'PRN code [#]');
        
        var_spAzimuthAngleDegOrbit = netcdf.defVar(netcdf_cyg,'spAzimuthAngleDegOrbit','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_spAzimuthAngleDegOrbit, 'spAzimuthAngleDegOrbit', 'Azimuth angle at specular point in the orbit frame [deg]');

        var_EIRP_L1 = netcdf.defVar(netcdf_cyg,'EIRP_L1','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_EIRP_L1, 'EIRP_L1', 'Effective Isotropic Radiated Power [Watt]');

        var_rxAntennaGain_L1_L = netcdf.defVar(netcdf_cyg,'rxAntennaGain_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_rxAntennaGain_L1_L, 'rxAntennaGain_L1_L', 'Receiver antenna gain toward specular point [dB]');

        var_powerAnalogW_L1_L = netcdf.defVar(netcdf_cyg,'powerAnalogW_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_powerAnalogW_L1_L, 'powerAnalogW_L1_L', 'Peak power of the DDM. It the true power that would have been measured by an ideal (analog) power sensor and corrected for quantization effects [Watt]');
 
        var_noiseFloorCounts_L1_L = netcdf.defVar(netcdf_cyg,'noiseFloorCounts_L1_L','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_noiseFloorCounts_L1_L, 'noiseFloorCounts_L1_L', 'Noise floor of the DDM [Watt]');
        
        var_bitRatio = netcdf.defVar(netcdf_cyg,'bitRatio','NC_FLOAT',dimid);
        netcdf.putAtt(netcdf_cyg, var_bitRatio, 'bitRatio', 'Port low/high bit counter ratio defined as (plus_1_cnts + minus_1_cnts) / (plus_3_cnts + minus_3_cnts).');

        %var_reflectivityLinear_L1_L = netcdf.defVar(netcdf_cyg,'reflectivityLinear_L1_L','NC_DOUBLE',dimid);
        %netcdf.putAtt(netcdf_cyg, var_reflectivityLinear_L1_L, 'reflectivityLinear_L1_L', 'Reflection coefficient of GPS signal L1, left polarization [dB]');

        % End definition mode
        netcdf.endDef(netcdf_cyg);
    
        % Write data to variables
        timeUTC_str = string(timeUTC);
        netcdf.putVar(netcdf_cyg, var_timeUTC, timeUTC_str);
        netcdf.putVar(netcdf_cyg, var_specularPointLat, specularPointLat);
        netcdf.putVar(netcdf_cyg, var_specularPointLon, specularPointLon);
        netcdf.putVar(netcdf_cyg, var_incidenceAngleDeg, incidenceAngleDeg);
        netcdf.putVar(netcdf_cyg, var_spAzimuthAngleDegNorth, spAzimuthAngleDegNorth);
        netcdf.putVar(netcdf_cyg, var_SNR_L1_L, SNR_L1_L);
        netcdf.putVar(netcdf_cyg, var_powerRatio_L1_L, powerRatio_L1_L);
        netcdf.putVar(netcdf_cyg, var_coherencyRatio_L1_L, coherencyRatio_L1_L);
        netcdf.putVar(netcdf_cyg, var_reflectivityPeak_L1_L, reflectivityPeak_L1_L);
        netcdf.putVar(netcdf_cyg, var_reflectivityPeakRecal_L1_L, reflectivityPeakRecal_L1_L);
        netcdf.putVar(netcdf_cyg, var_NBRCS_L1_L, NBRCS_L1_L);
        netcdf.putVar(netcdf_cyg, var_qualityControlFlags_L1_L, qualityControlFlags_L1_L);
        netcdf.putVar(netcdf_cyg, var_qualityControlFlags_2_L1_L, qualityControlFlags_2_L1_L);
        netcdf.putVar(netcdf_cyg, var_notToBeUsed, uint8(notToBeUsed));
        netcdf.putVar(netcdf_cyg, var_notRecommended, uint8(notRecommended));
        netcdf.putVar(netcdf_cyg, var_uncertainty, coefficientOfVariation);
        netcdf.putVar(netcdf_cyg, var_constellation, constellation);
        netcdf.putVar(netcdf_cyg, var_recevingAntenna, receivingAntenna);
        netcdf.putVar(netcdf_cyg, var_receivingSpacecraft, receivingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_transmittingSpacecraft, transmittingSpacecraft);
        netcdf.putVar(netcdf_cyg, var_pseudoRandomNoise, pseudoRandomNoise);
        netcdf.putVar(netcdf_cyg, var_spAzimuthAngleDegOrbit, spAzimuthAngleDegOrbit);
        netcdf.putVar(netcdf_cyg, var_EIRP_L1, EIRP_L1);
        netcdf.putVar(netcdf_cyg, var_rxAntennaGain_L1_L, rxAntennaGain_L1_L);
        netcdf.putVar(netcdf_cyg, var_powerAnalogW_L1_L, powerAnalogW_L1_L);
        netcdf.putVar(netcdf_cyg, var_noiseFloorCounts_L1_L, noiseFloorCounts_L1_L);
        netcdf.putVar(netcdf_cyg, var_bitRatio, bitRatio);
        %netcdf.putVar(netcdf_cyg, var_reflectivityLinear_L1_L, reflectivityLinear_L1_L)

%% ------- Close netcdf file ------- %%
        netcdf.close(netcdf_cyg) ; 
end
