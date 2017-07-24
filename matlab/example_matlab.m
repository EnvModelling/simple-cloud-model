% open file
nc=netcdf('/tmp/output.nc');


figure;
subplot(231);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,2)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'q_l (kg kg^{-1})');

subplot(232);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,3)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'q_r (kg kg^{-1})');


subplot(233);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,4)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'q_s (kg kg^{-1})');

subplot(234);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,5)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'q_g (kg kg^{-1})');

subplot(235);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,6)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'q_i (kg kg^{-1})');


subplot(236);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'precip'}(:,:,1)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'P (mm hr^{-1})');


figure;
subplot(131);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,7)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'n_i (# kg^{-1})');

subplot(132);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,8)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'n_s (# kg^{-1})');

subplot(133);
pcolor(nc{'time'}(:)./3600,nc{'z'}(:,:)'./1000,nc{'q'}(:,:,9)');shading flat
xlabel('time (hrs)');
ylabel('z (km)')
h=colorbar
ylabel(h,'n_g (# kg^{-1})');



% close file
close(nc);
