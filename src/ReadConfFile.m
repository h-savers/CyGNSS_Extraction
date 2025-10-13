function [Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath, LatMin, LatMax, LonMin, LonMax] = ReadConfFile(configurationPath)
%%%%%%%  Read configuration file
%
            lines = string(splitlines(fileread(configurationPath)));
%         
            ConfigRightLine= contains(lines,'Taskname')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            Taskname= extractAfter(lines(ConfigRightLine),startIndex) ;         
%%         
            ConfigRightLine= contains(lines,'initdate')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            initdate= extractAfter(lines(ConfigRightLine),startIndex) ;
%%         
            ConfigRightLine= contains(lines,'enddate')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            enddate= extractAfter(lines(ConfigRightLine),startIndex) ;
%%         
            ConfigRightLine= contains(lines,'savespace')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            savespace= extractAfter(lines(ConfigRightLine),startIndex) ;
%%         
            ConfigRightLine= contains(lines,'CyGinpath')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            CyGinpath= extractAfter(lines(ConfigRightLine),startIndex) ;
%%                  
            ConfigRightLine= contains(lines,'CyGoutpath')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            CyGoutpath= extractAfter(lines(ConfigRightLine),startIndex) ; % max distance between SP and SMAP grid cell in meters
            %%                  
            ConfigRightLine= contains(lines,'LatMin')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            LatMin= extractAfter(lines(ConfigRightLine),startIndex) ; % max distance between SP and SMAP grid cell in meters
            LatMin=double(LatMin) ; 
%%                  
            ConfigRightLine= contains(lines,'LatMax')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            LatMax= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            LatMax=double(LatMax) ; 
            %%                  
            ConfigRightLine= contains(lines,'LonMin')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            LonMin= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            LonMin=double(LonMin) ; 

            %%                  
            ConfigRightLine= contains(lines,'LonMin')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            LonMin= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            LonMin=double(LonMin) ; 

            %%                  
            ConfigRightLine= contains(lines,'LonMax')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            LonMax= extractAfter(lines(ConfigRightLine),startIndex) ; %
            LonMax=double(LonMax) ; 

            % %%                  
            ConfigRightLine= contains(lines,'logpath')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            logpath= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            %%
end