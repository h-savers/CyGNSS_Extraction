%%%%%%%%%%%%%%%%%%%%% plotting the CyGNSS map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scattermap(targetvar,lat,lon,datechar,lowlim,highlim)
  targetvar(~isfinite(targetvar))=NaN;
  scatter(lon,lat,10,targetvar,'filled');
  axis([-180 180 -80 80])
  caxis([lowlim highlim])
  title(datechar)
  xlabel('lon')
  ylabel('lat')
  colorbar
  grid on
end
