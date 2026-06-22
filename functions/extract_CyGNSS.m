%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just_scale_factor version: extracts ONLY the scale factor (SS_r) together
% with the geolocation/time identifiers needed to make it matchable, looping
% over all CyGNSS satellites available for the day. No DDM-derived products.
%
% Inputs:
%   datechar: date in 'yyyymmdd' format
%   doy: day of the year
%   inpath: input path for CyGNSS data
%   logpath: path for logging errors
%   savespace: flag to apply the CYGNSS land flag before saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mission, L1b_product, L1b_product_version, timeUTC, ...
        DoY, SoD, SCID, SV_NUM, PRN, SPLAT, SPLON, SCALE_FACTOR] = ...
        extract_CyGNSS(datechar, doy, inpath, logpath, savespace)

    %%%%%%%%%%%%%%%%%%%% INITIALISING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timeUTC=[];
    Year=[];
    SoD=[];                                 % second of the day
    DoY=[];                                 % day of the year
    SCID=[];                                % receiving (CYGNSS) sat ID
    SV_NUM=[];                              % transmitting spacecraft number
    PRN=[];                                 % PRN code
    SPLAT=[];                               % SP lat on ground
    SPLON=[];                               % SP lon on ground
    SCALE_FACTOR=[];                        % Scale factor (SS_r) computed from tx/rx ranges

    chkfile=dir([inpath 'cyg0*.nc']);
    sat_index=[];
    for i=1:length(chkfile)
        sat_index=[sat_index, str2num(chkfile(i).name(5))];
    end

    for jj=sat_index     % loop on the satellites available for the day
        chkfile=dir([inpath 'cyg0' num2str(jj) '.ddmi.s' datechar '*.nc']);   % to avoid end of execution in case file is missing
        if ~isempty(chkfile)
            infile=chkfile.name;
            disp(['% reading satellite ' num2str(jj) ' - date ' datechar '- file ' infile ])
            [mission, L1b_product, L1b_product_version, sp_lat, sp_lon, scid, sv_num, ts, prn, SS_r] = ...
                readnc_CyGNSS_v2(inpath, infile, savespace);

            dayofyear=zeros(size(sp_lat)) + doy;                  % to have the same size as sp_lat
            year=zeros(size(sp_lat)) + str2double(datechar(1:4)); % year for UTC time

            % cat variables
            SCID=cat(1,SCID,scid(:));
            Year=cat(1,Year,year(:));
            SoD=cat(1,SoD,ts(:));
            DoY=cat(1,DoY,dayofyear(:));
            timeUTC=datetime(Year, 1, 1) + days(DoY - 1) + seconds(SoD); % UTC time
            SV_NUM=cat(1,SV_NUM,sv_num(:));
            PRN=cat(1,PRN, prn(:));
            SPLAT=cat(1,SPLAT, sp_lat(:));
            SPLON=cat(1,SPLON, sp_lon(:));
            SCALE_FACTOR=cat(1,SCALE_FACTOR, SS_r(:));
        else
            diary([logpath 'log_' datestr(now,'dd-mm-yyyy') '.txt'])
            disp(['% WARNING: cyg0' num2str(jj) ' satellite missing for the date ' datechar])
            diary off
        end
    end
end
