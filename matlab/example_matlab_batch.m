% precip is in mm hr-1 for rain, snow, graupel, ice
% if you do a cumsum of the lowest level of the precip array, multiplied by
% the factor (10./3600), it will give total precip in mm.
ARRAY1=[290. 285. 280. 275.]; % cloud base
ARRAY2=[268. 265. 262. 259. 256. 253. 250. 247. 244.]; % cloud top
ARRAY3=[0.1 1. 10. 100. 1000. 10000.]; % number of ice crystals

% make two empty arrays to take the total precip
tot_precip_hm_on=zeros(length(ARRAY1),length(ARRAY2),length(ARRAY3));
tot_precip_hm_off=zeros(length(ARRAY1),length(ARRAY2),length(ARRAY3));

for k=1:length(ARRAY3)
    for j=1:length(ARRAY2)
        for i=1:length(ARRAY1)
            % first read the hm on runs
            nc=netcdf(['/tmp/output_',num2str(i-1),'_',num2str(j-1),'_',num2str(k-1),'_hm_on.nc']);
            dt=nc{'time'}(2)-nc{'time'}(1);
            % this array below could be plotted vs the time array to plot
            % precip vs time:
            precip_vs_time=cumsum(nc{'precip'}(:,1,1).*dt./3600);
            %plot(nc{'time'}(:),precip_vs_time)
            % save the last element in the array
            tot_precip_hm_on(i,j,k)=precip_vs_time(end);
            
            close(nc);
            
            % first read the hm on runs
            nc=netcdf(['/tmp/output_',num2str(i-1),'_',num2str(j-1),'_',num2str(k-1),'_hm_off.nc']);
            % this array below could be plotted vs the time array to plot
            % precip vs time:
            precip_vs_time=cumsum(nc{'precip'}(:,1,1).*dt./3600);
            %plot(nc{'time'}(:),precip_vs_time)
            tot_precip_hm_off(i,j,k)=precip_vs_time(end);
            close(nc);
            
        end
    end
end

figure('name','HM on')
for i=1:6
    subplot(2,3,i)
%    [c,h]=contourf(ARRAY2,ARRAY1,tot_precip_hm_on(:,:,i),[0:5:90]);
    [c,h]=contourf(ARRAY2,ARRAY1,tot_precip_hm_on(:,:,i),[0:1:90]);
    caxis([0 90]);
    clabel(c,h);
end

figure('name','HM off')
for i=1:6
    subplot(2,3,i)
%    [c,h]=contourf(ARRAY2,ARRAY1,tot_precip_hm_off(:,:,i),[0:5:90]);
    [c,h]=contourf(ARRAY2,ARRAY1,tot_precip_hm_off(:,:,i),[0:1:90]);
    caxis([0 90]);
    clabel(c,h);
end