function ddm_std=ddm_pseudo_std(pa)
     [Doppler_bins,Range_bins, numSPtrack] = size(pa) ; 
     pa_clean=pa ; 
     pa_clean(pa_clean <0)= 0 ; 
     ddm_mean=squeeze(sum(pa_clean, [1,2])) ;
     pa_norm = pa_clean./sum(pa_clean, [1,2]);
     max_pa=max(pa_clean, [], [1,2]) ;

     % [ipeak, jpeak]=find(pa(:,:,k)==max(max(pa(:,:,k)))) ; 
     for k=1:numSPtrack, mass=max(max(pa(:,:,k))); [ipeak(k), jpeak(k)]=find(pa(:,:,k)==mass, 1,'first') ; end
     ipeak_mat=reshape(ipeak,1,1,[]).*ones(Doppler_bins, Range_bins, 1);
     jpeak_mat=reshape(jpeak,1,1,[]).*ones(Doppler_bins,Range_bins,1);
    
     i_mat=[1:1:Doppler_bins]; i_mat=repmat(i_mat', 1, Range_bins, numSPtrack) ;
     j_mat=[1:1:Range_bins]; j_mat=repmat(j_mat, Doppler_bins, 1, numSPtrack) ;
     dist_square=(i_mat-ipeak_mat).^2+(j_mat-jpeak_mat).^2 ; 

     ddm_std=squeeze(sum(pa_norm.*dist_square, [1,2]));    
     ddm_std=sqrt(ddm_std) ; 
    
end
