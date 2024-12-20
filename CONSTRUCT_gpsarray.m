% file to create fake GNSS stations in Mexico
% AUTHORS:
% Rishav Mallick, JPL, 2024

% provide center lat,lon to sample around (need not be mesh transformation
% coordinates)
lat0 = 18;
lon0 = -98;

ngps = 200;
gps = [lat0 + 1.*randn(ngps,1), lon0 + 2.*randn(ngps,1)];

% read megathrust mesh file
eM = unicycle.greens.nikkhoo15(30e3,0.25);
% read megathrust file in geographic coordinates
rcvll = unicycle.geometry.triangleReceiver('mexico_mt_latlon',eM);

figure(1),clf
rcvll.plotPatch(), hold on
plot(gps(:,2),gps(:,1),'kv')
axis tight equal
view(0,90)

% write to file
writetable(table(gps(:,2),gps(:,1)),'obs_sites.csv','WriteVariableNames',false)