% precip is in mm hr-1 for rain, snow, graupel, ice
% if you do a cumsum of the lowest level of the precip array, multiplied by
% the factor (10./3600), it will give total precip in mm.
ARRAY1=[290. 285. 280. 275.]; % cloud base
ARRAY2=[253.]; % cloud top
ARRAY3=[0.1 1. 10. 100. 1000. 10000.]; % number of ice crystals

% make two empty arrays to take the total precip
tot_precip_hm_on=zeros(length(ARRAY1),length(ARRAY2),length(ARRAY3));
tot_precip_hm_off=zeros(length(ARRAY1),length(ARRAY2),length(ARRAY3));


figure('name','hm plot');
l=1;
for k=1:length(ARRAY3)
    for j=1:length(ARRAY2)
        for i=2:2 %length(ARRAY1)
            subplot(3,2,l);
            % first read the hm on runs
            nc=netcdf(['/tmp/output_',num2str(i-1),'_',num2str(j-1),'_',num2str(k-1),'_hm_on.nc']);
            precip1=nc{'precip'}(:,:,1);
            close(nc);
            
            % first read the hm on runs
            nc=netcdf(['/tmp/output_',num2str(i-1),'_',num2str(j-1),'_',num2str(k-1),'_hm_off.nc']);
            precip2=nc{'precip'}(:,:,1);
            pcolor(nc{'time'}(:)./60,nc{'z'}(1,:)./1000,(precip1-precip2)');shading flat
            l=l+1;
            hold on;
            [c,h]=contour(nc{'time'}(:)./60,nc{'z'}(1,:)./1000,nc{'q'}(:,:,2)',linspace(1e-5,1e-3,10),'k');
            caxis([-0.5 0.5]);
            if(l>6)
                xlabel('time (mins)');
            end
            %if(mod(l+2,4)==0)
                ylabel({'z (km)',['INP: ',num2str(ARRAY3(k)./1000),' L^{-1}']})
            %end
            %if(l<=5)
                text(0.1,0.9,{['Cloud-base temp: '],[num2str(ARRAY1(i)-273),' degC']},'fontsize',7,'units','normalized')
            %end
            close(nc);
        end
    end
end
h=colorbar;
set(h,'position',[0.9304    0.1405    0.0071    0.7833]);
ylabel(h,'difference in precipitation rate: HM - no HM (mm hr^{-1})');
