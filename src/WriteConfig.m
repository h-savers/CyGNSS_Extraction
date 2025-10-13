function WriteConfig(configurationPath,Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax)
conffileID = fopen(configurationPath, 'W') ; 
% conffileID = fopen(configurationPath) ; 
fprintf(conffileID,'%s',['Taskname=' Taskname] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['initdate=' initdate] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['enddate=' enddate] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['savespace=' savespace] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s',['CyGinpath=' CyGinpath] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['CyGoutpath=' CyGoutpath] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,'%s', ['logpath=' logpath] ); fprintf(conffileID,'\n') ; 


fprintf(conffileID,['LatMin=' char(string(LatMin))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LatMax=' char(string(LatMax))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LonMin=' char(string(LonMin))] ); fprintf(conffileID,'\n') ; 
fprintf(conffileID,['LonMax=' char(string(LonMax))] ); fprintf(conffileID,'\n') ; 


% fprintf(conffileID,'%s',['Dayinit=' Dayinit] ); fprintf(conffileID,'\n') ; 
% fprintf(conffileID,'%s',['Dayfinal=' Dayfinal] ); fprintf(conffileID,'\n') ; 
% fprintf(conffileID,['DDM=' DDM] ); fprintf(conffileID,'\n') ; 
fclose(conffileID) ;
end