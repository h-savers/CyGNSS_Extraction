function WriteConfig(configurationPath,Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, calibration_file, ...
    LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
    snr_th, rx_gain_th, inc_angl_th, nsnr_th, coherency_th)
conffileID = fopen(configurationPath, 'W') ; 
% conffileID = fopen(configurationPath) ; 
fprintf(conffileID,'%s',['Taskname=' Taskname] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['initdate=' initdate] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['enddate=' enddate] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['savespace=' savespace] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['CyGinpath=' CyGinpath] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['CyGoutpath=' CyGoutpath] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['logpath=' logpath] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['calibration_file=' calibration_file] ); fprintf(conffileID,'\n') ; 


fprintf(conffileID,['LatMin=' char(string(LatMin))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LatMax=' char(string(LatMax))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LonMin=' char(string(LonMin))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LonMax=' char(string(LonMax))] ); fprintf(conffileID,'\n') ; 

fprintf(conffileID,'%s',['aggregate_data=' aggregate_data] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['out_format=' out_format] ); fprintf(conffileID,'\n') ; 


% fprintf(conffileID,'%s',['Dayinit=' Dayinit] ); fprintf(conffileID,'\n') ; 
% fprintf(conffileID,'%s',['Dayfinal=' Dayfinal] ); fprintf(conffileID,'\n') ; 
% fprintf(conffileID,['DDM=' DDM] ); fprintf(conffileID,'\n') ; 

fprintf(conffileID,['snr_dB_th=' char(string(snr_th))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['rx_gain_dB_th=' char(string(rx_gain_th))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['inc_angl_th=' char(string(inc_angl_th))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['nsnr_th=' char(string(nsnr_th))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['coherency_th=' char(string(coherency_th))] ); fprintf(conffileID,'\n') ; 


fclose(conffileID) ;
end