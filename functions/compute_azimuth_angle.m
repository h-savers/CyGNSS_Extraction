function [sp_azimuth_angle_deg_north] = compute_azimuth_angle(sc_pos_x,sc_pos_y,sc_pos_z,sc_vel_x,sc_vel_y,sc_vel_z,azimuth_angle)

    wgs84 = wgs84Ellipsoid('meters'); % CYGNSS use the WGS84 ellipsoid
    
    [lat, lon, h] = ecef2geodetic(wgs84, sc_pos_x, sc_pos_y, sc_pos_z);
    phi = deg2rad(lat);    % rad
    lambda = deg2rad(lon); % rad
    
    R11 = -sin(lambda); % Rotation matrix
    R12 = cos(lambda);
    R13 = zeros(size(lambda));
    
    R21 = -sin(phi).*cos(lambda);
    R22 = -sin(phi).*sin(lambda);
    R23 = cos(phi);
    
    R31 = cos(phi).*cos(lambda);
    R32 = cos(phi).*sin(lambda);
    R33 = sin(phi);
    
    V_ECEF = [sc_vel_x'; sc_vel_y'; sc_vel_z']; % Velocity in ECEF

    V_local = zeros(3, length(sc_pos_x)); 
    
    for k = 1:length(sc_pos_x)
        Rk = [R11(k), R12(k), R13(k);
              R21(k), R22(k), R23(k);
              R31(k), R32(k), R33(k)];
        V_local(:,k) = Rk * V_ECEF(:,k);
    end
    
    vN = V_local(2,:);
    vE = V_local(1,:);
    heading_deg = atan2(vE, vN) * 180/pi;
    sp_azimuth_angle_deg_north=heading_deg+azimuth_angle;
    sp_azimuth_angle_deg_north(sp_azimuth_angle_deg_north>360)=sp_azimuth_angle_deg_north(sp_azimuth_angle_deg_north>360)-360;