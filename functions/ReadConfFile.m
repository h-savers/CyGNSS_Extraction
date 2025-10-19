function [Taskname, initdate, enddate, savespace, CyGinpath, CyGoutpath, logpath,...
    LatMin, LatMax, LonMin, LonMax, aggregate_data, out_format,...
    snr_th, rx_gain_th, inc_angl_th, nsnr_th] = ReadConfFile(configurationPath)
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
            %%                                   
            ConfigRightLine= contains(lines,'logpath')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            logpath= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            %%                   
            ConfigRightLine= contains(lines,'aggregate_data')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            aggregate_data= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            %%                   
            ConfigRightLine= contains(lines,'out_format')  ;  
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            out_format= extractAfter(lines(ConfigRightLine),startIndex) ; % 
            %%   
            ConfigRightLine= contains(lines,'snr_dB_th')  ;                   % dB, minimum SNR
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            snr_th= extractAfter(lines(ConfigRightLine),startIndex) ; 
            snr_th=double(snr_th) ; 
            

            ConfigRightLine= contains(lines,'rx_gain_dB_th')  ;               % dBi, minimum receiver gain
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            rx_gain_th= extractAfter(lines(ConfigRightLine),startIndex) ; 
            rx_gain_th=double(rx_gain_th) ; 

            ConfigRightLine= contains(lines,'inc_angl_th')  ;               % degrees, maximum incidence angle
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            inc_angl_th= extractAfter(lines(ConfigRightLine),startIndex) ; 
            inc_angl_th=double(inc_angl_th) ;  

            ConfigRightLine= contains(lines,'nsnr_th')  ;                  % dB, maximum NSNR
            ConfigRightLine= find(ConfigRightLine==1)  ;   
            startIndex= regexp(lines(ConfigRightLine),'=') ; 
            nsnr_th= extractAfter(lines(ConfigRightLine),startIndex) ;
            nsnr_th=double(nsnr_th) ; 

end