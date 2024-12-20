% BEM solution to interseismic steady creep on a megathrust
% the user must provide 
% (1) a triangular fault mesh in E,N,U & lat,lon coordinates, 
% (2) euler pole (lat,lon,rotation rate in deg/Ma) 
% (3a) locked domain on the megathrust as outlines
% (3b) freely creeping domain on the megathrust as a depth cutoff
% 
% AUTHORS:
% Rishav Mallick, JPL, 2024

clear
import unicycle.*

% read file with locked asperities
outlinefilename = 'fault/mexico_outlines.csv';
stationfilename = 'data/obs_sitesK.csv';

% provide slab2 file for plotting
slabdepthfilename = 'fault/cam_slab2_dep_02.24.18-neglon.grd';

data = readtable(outlinefilename);
% provide reference point to transform from lat,lon to x,y (use same
% reference coordinates as the mesh)
lat0 = 16;
lon0 = -100;

% how deep is your locking? in [m]
deeplockingextent = 80e3;

% euler pole for CO-NA (lat,lon,rotation rate in deg/Ma)
% http://www.geology.wisc.edu/homepages/chuck/public_html/PDF/demets_grl.pdf
latp = 28.1;
lonp = -123.1;
degmyr = 1.402;

% read megathrust mesh file
eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver('fault/mexico_mt_xyz',eM);% in local Cartesian coordinates
% read megathrust file in geographic coordinates
rcvll = unicycle.geometry.triangleReceiver('fault/mexico_mt_latlon',eM);

figure(10),clf
rcv.plotPatch(rcv.xc(:,3)./1e3), hold on
rcv.plotPatch()
axis tight equal, grid on, box on
cb=colorbar;cb.Label.String = 'depth (km)';
view(0,90)
clim([-120 0])
colormap(turbo(12))
xlabel('East (m)'), ylabel('North (m)'), zlabel('depth (m)')
set(gca,'FontSize',15,'LineWidth',1.5)

% compute traction kernels
% Kss,Ksd,Ksn - traction kernels due to strike-slip on fault
% Kds,Kdd,Kdn - traction kernels due to dip-slip on fault
if exist('kernels/tractionkernels.mat','file')~=2
    [Kss,Kds,Ksd,Kdd,Ksn,Kdn] = rcv.tractionKernels(rcv);
else
    load('kernels/tractionkernels.mat')
end
% combine kernels
bigK = [Kss,Kds;Ksd,Kdd];

% compute eigen values to assess stability
if false
    tic
    [Evec,Evals] = eig(bigK);
    toc
    lambda = diag(Evals);
end
%% interseismic simulation
% label locked patches with a logical array
locked = false(rcv.N,1);
for i = 2:5
    % get index for asperity outline
    index = data{:,6} == i;
    [xnode, ynode, ~] = geodetic2enu(data{index,2},data{index,1},0,...
    lat0,lon0,0,wgs84Ellipsoid);
    in = inpolygon(rcv.xc(:,1),rcv.xc(:,2),xnode,ynode);
    locked = locked|in;
end

% combine locked and deep creep patches (deeper than a user-provided depth)
deepcreep = abs(rcv.xc(:,3)) > deeplockingextent;
pinned = locked | deepcreep;

%% calculate long-term stress rate from long-term slip rate
% vs0 = zeros(rcv.N,1);
% vd0 = ones(rcv.N,1);
% compute horizontal plate motion at each point on the mesh
[vE,vN] = pole_velocity(rcvll.xc(:,2),rcvll.xc(:,1),latp,lonp,degmyr);
% rotate from E,N coordinates to strike-slip,dip-slip and normalize
[vs0,vd0] = rotate_velocities(-vE,-vN,rcv);
Vplate = sqrt(vE.^2 + vN.^2);

%% compute long-term stress rate accounting for locked and fully creeping patches
% specify τ_0 = K(x2;x1).(v0-v)(x1) + K(x2;x2).(v0 - 0) +
% K(x2;x3).(v0-v)(x3) 
% where x1 - locked (v=0), x2 - bem, x3 - deep creep at (v=v0)
vs_i = zeros(rcv.N,1);
vd_i = zeros(rcv.N,1);
vs_i(deepcreep) = vs0(deepcreep);
vd_i(deepcreep) = vd0(deepcreep);
vs_i(locked) = 0;
vd_i(locked) = 0;

% compute lopng-term stress rate
bigKmod = [Kss(~pinned,:),Kds(~pinned,:);Ksd(~pinned,:),Kdd(~pinned,:)];
tau0 = -bigKmod*[(vs0-vs_i);(vd0-vd_i)];

%% BEM solution for interseismic slip rate distribution
vsol = zeros(rcv.N*2,1);
% solve matrix BEM equation: equilibrium equations only at BEM nodes
% (~pinned)
sol = -[Kss(~pinned,~pinned),Kds(~pinned,~pinned);...
        Ksd(~pinned,~pinned),Kdd(~pinned,~pinned)] \ tau0;
vsol([~pinned;~pinned]) = sol;
% extract strike-slip and dip-slip component
vs = vsol(1:rcv.N);
vd = vsol(rcv.N+1:end);

% slip rate for locked patches and deep creep patches are imposed
vs(deepcreep) = vs_i(deepcreep);
vd(deepcreep) = vd_i(deepcreep);
vnorm = sqrt(vs.^2 + vd.^2);

%% compute station velocities
data = readtable(stationfilename);
[xnode, ynode, ~] = geodetic2enu(data{:,2},data{:,1},0,...
    lat0,lon0,0,wgs84Ellipsoid);
% compute displacement kernels
[Gs,Gd]=rcv.displacementKernels([xnode,ynode,zeros(length(xnode),1)],3);
% separate into e,n,u components
Ge = [Gs(1:3:end,:),Gd(1:3:end,:)];
Gn = [Gs(2:3:end,:),Gd(2:3:end,:)];
Gz = [Gs(3:3:end,:),Gd(3:3:end,:)];
% compute respective surface velocity field components
ve_bem = Ge*[(vs-vs0);(vd-vd0)];
vn_bem = Gn*[(vs-vs0);(vd-vd0)];
vz_bem = Gz*[(vs-vs0);(vd-vd0)];

%% plot results in Cartesian coordinates
figure(1),clf
subplot(2,1,1)
rcv.plotPatch(vs./Vplate), hold on
% rcv.plotPatch()
rcv.plotSlipVectors(vs./Vplate,vd./Vplate,1e4,'k')
axis tight equal, grid on, box on
cb=colorbar;cb.Label.String = 'v_s';
view(0,90)
clim([-1 1]*0.5)
xlabel('East (m)'), ylabel('North (m)'), zlabel('depth (m)')
set(gca,'FontSize',15,'LineWidth',1.5)

subplot(2,1,2)
rcv.plotPatch(vd./Vplate), hold on
rcv.plotPatch()
axis tight equal, grid on, box on
cb=colorbar;cb.Label.String = 'v_d';
view(0,90)
clim([-0 1])
colormap(turbo(100))
xlabel('East (m)'), ylabel('North (m)'), zlabel('depth (m)')
set(gca,'FontSize',15,'LineWidth',1.5)

figure(2),clf
set(gcf,'Color','w')
rcv.plotPatch(vnorm./Vplate), hold on
rcv.plotPatch()
plot(xnode,ynode,'kv','MarkerFaceColor','w')
quiver([xnode;xnode],[ynode;ynode],[ve_bem;ve_bem.*0],[vn_bem;vz_bem],'k','LineWidth',2)
axis tight equal, grid on, box on
cb=colorbar;cb.Label.String = '(v/v_{pl})';cb.LineWidth=1.5;cb.TickDirection='out';
view(0,90)
clim([0 1])
colormap(parula(10))
xlabel('East (m)'), ylabel('North (m)'), zlabel('depth (m)')
set(gca,'FontSize',20,'LineWidth',1.5,'TickDir','both')

% return

%% plot in geographical coordinates
[lon,lat,slabdepth] = unicycle.export.grdread(slabdepthfilename);

figure(3),clf
scf = 0.04;% scaling factor for quiver plots

set(gcf,'Color','w')
rcvll.plotPatch(vnorm./Vplate), hold on
rcvll.plotPatch(), alpha 0.5
[c,h] = contour(lon,lat,slabdepth,[-5,-10,-20,-30,-50,-80,-100,-150],'-','EdgeColor',[0.8,0.1,0.5],'LineWidth',2);
clabel(c,h,'FontWeight','normal','Fontsize',20,'LabelSpacing',700)

plot(data{:,1},data{:,2},'kv','MarkerFaceColor','w')
quiver(data{:,1},data{:,2},ve_bem.*scf,vn_bem.*scf,0,'r','LineWidth',2)
% quiver([data{:,1};data{:,1}],[data{:,2};data{:,2}],[ve_bem;data{:,3}.*1000],[vn_bem;data{:,4}.*1000],1,'r','LineWidth',2)
% [repmat([1,0,0],length(ve_obs),1);repmat([1,0,1],length(ve_obs),1)]
quiver(data{:,1},data{:,2},data{:,3}.*scf,data{:,4}.*scf,0,'k','LineWidth',2)
axis tight equal, grid off, box on
cb=colorbar;cb.Label.String = '(v/v_{pl})';cb.LineWidth=1.5;cb.TickDirection='out';
view(0,90)
clim([0 1])
colormap(parula(10))
xlim([-102,-93]), ylim([14,20])
xlabel('Longitude'), ylabel('Latitude'), zlabel('depth (m)')
set(gca,'FontSize',20,'LineWidth',1.5,'TickDir','both')
print('Figures/datafits','-djpeg','-r200')
return
%% plot results using mapping toolbox
figure(4),clf
set(gcf,'Color','w')
for i = 1:rcvll.N
    latplot = [rcvll.x(rcvll.vertices(i,:),2);rcvll.x(rcvll.vertices(i,1),2)];
    lonplot = [rcvll.x(rcvll.vertices(i,:),1);rcvll.x(rcvll.vertices(i,1),1)];
    aoi=geopolyshape(latplot,lonplot);
    geoplot(aoi,'FaceAlpha',0.3,ColorData=vnorm(i)), hold on
end
% geoscatter(rcvll.xc(:,2),rcvll.xc(:,1),log10(rcvll.area+1).*2e2,vnorm, 'filled')
geobasemap grayterrain
cb=colorbar;cb.Label.String = '(v/v_{pl})';cb.LineWidth=1.5;cb.TickDirection='out';
clim([0,1]),colormap(parula(10))
set(gca,'FontSize',20,'TickDir','both')
grid off