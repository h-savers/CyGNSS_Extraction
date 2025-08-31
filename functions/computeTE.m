%%%%%%%%%%%%%%%%%%%%%%%%%% function computeTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes trailing edge from the data read by readCyGNSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TE_width=computeTE(pa,delay_vector,Power_threshold)
     SPnumbers=size(pa,3);
     Bins_With_Peak=zeros(1,SPnumbers);
     c=zeros(1,SPnumbers);  
     var=zeros(1,SPnumbers);  
     waveform_peak=zeros(1,SPnumbers);  
     TE_width=zeros(1,SPnumbers); 
     peak_WF=zeros(1,SPnumbers); 
     peak_delay_WF=zeros(1,SPnumbers); 
     WF_threshold_input=zeros(1,SPnumbers); 
     diff_peak_threshold=zeros(1,SPnumbers); 
     min_diff_threshold=zeros(1,SPnumbers); 
     closestIndex=zeros(1,SPnumbers); 
     WF_threshold_value=zeros(1,SPnumbers); 
     ind_thr_value=zeros(1,SPnumbers); 
     ind_peak=zeros(1,SPnumbers); 
     for ss=1:SPnumbers  
             DDM=pa(:,:,ss);                    
             Bins_With_Peak(ss)=sum(DDM(:)==max(DDM(:)));       % count Number of bins with peak value in the whole DDM
             zero_Dopp = DDM(6,:);
             zero_Dopp((zero_Dopp <0))=0 ; 
             %Computing Trailing Edge
             waveform_peak(ss) = max(zero_Dopp);
             if (sum(isnan(zero_Dopp))~=0) %check if DDM zero doppler in the region has an empty value
                 TE_width(ss) = NaN;
             else
                 cs = spline(delay_vector, zero_Dopp);
                 nn = size(zero_Dopp,2)*100 + 1; %index 1:1701
                 resampled_delay = linspace(delay_vector(1),delay_vector(end), nn);
                 resampled_WF = ppval(cs,resampled_delay);
                 [peak_WF(ss),ind_peak(ss)] = max(resampled_WF);
                 peak_delay_WF(ss) = resampled_delay(ind_peak(ss));
                  WF_threshold_input(ss) = peak_WF(ss)*Power_threshold; %power threshold 70% 
                 trail_peak_WF_ones = ones*peak_WF(ss); %vector of the WF_peak value
                 WF_trail = resampled_WF(ind_peak(ss):end); %vector representing values over the peak till the end
                 check_greater_than_threshold = WF_trail >= WF_threshold_input(ss);
                 if (all(check_greater_than_threshold) == 1) %condition if threshold is over the 1700 lag samples
                     TE_width(ss) = NaN;
                 elseif (ind_peak(ss) >= 1690) %check if peak is too shifted to the end of lags
                     TE_width(ss) = NaN;
                 else
                     flip_WF_trail = flip(WF_trail);
                     flip_threshold = flip_WF_trail(flip_WF_trail <= WF_threshold_input(ss)); 
                     WF_trail_second = flip(flip_threshold); %Selected value from the tail of resampled WF to select the nearest value of threshold that is farther from peak 
                     diff_peak_trail = trail_peak_WF_ones - WF_trail_second;
                     diff_peak_threshold(ss) = peak_WF(ss) - WF_threshold_input(ss); %power loss at 70% power reduction
                     min_diff_threshold(ss) = min(abs(diff_peak_trail-diff_peak_threshold(ss)));
                     closestIndex(ss) = find(min_diff_threshold(ss)) + (length(WF_trail) - length(WF_trail_second));
                     WF_threshold_value(ss) = WF_trail(closestIndex(ss)); %closest value to the power threshold 70%, that means "WF_threshold_input" 
                     eligible_values_for_threshold = find(resampled_WF == WF_threshold_value(ss));
                     ind_thr_value(ss) = eligible_values_for_threshold(end); %select the farthest value from peak
                     TE_width(ss) = (resampled_delay(ind_thr_value(ss)) - peak_delay_WF(ss))* 1.5* 10^8 *10^-6; %trailing edge width in metres   
                 end
             end
     end
end