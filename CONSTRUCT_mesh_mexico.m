% script to read a slab2 grd file + polygon outline and construct a mesh
% AUTHOR:
% Rishav Mallick, 2024, JPL

clear
import unicycle.*

outfile = 'fault/mexico';
slabdepthfilename = 'fault/cam_slab2_dep_02.24.18-neglon.grd';
% outlinefilename = 'slab2-outline-simplified.gmt';
% asperitiesfilename = 'asperities_outline.gmt';
outlinefilename = 'fault/mexico_outlines.csv';

% provide cutoff depth for fault mesh in km
cutoff_depth = -155;

% shift depth upwards uniformally to approximate a flat Earth
zshift = 4;

% select a central point to convert from lat,lon -> xyz
% normally this is a mean value of all the points we are meshing
lat0 = 16;
lon0 = -100;

[lon,lat,slabdepth] = unicycle.export.grdread(slabdepthfilename);
% outlinedata = readtable(outlinefilename,"FileType","text","NumHeaderLines",11);
% asperitiesdata = readtable(asperitiesfilename,"FileType","text","NumHeaderLines",11);

g = create_decsgouts(outlinefilename);
% g_asp = create_decsgouts(asperitiesfilename);
% g_asp(6,:) = 0;
% g_asp(7,:) = 1;

% plot input data
figure(1),clf
pcolor(lon,lat,slabdepth), shading flat, hold on
pdegplot(g,"EdgeLabels","on","FaceLabels","on")
axis tight equal

%% create 2-d mesh using outline
Hmax = 0.5;
Hmin = 0.2;

% construct mesh as a 2-d object
modelmesh = createpde;
% g = decsg(Rdata);
% g2 = [[2;-100;-98;18;18;1;0],[2;-98;-98;18;17;1;0],[2;-98;-100;17;18;1;0]];
% gconcat = [g,g_asp];
geometryFromEdges(modelmesh,g);

meshobj = generateMesh(modelmesh,"Hmax",Hmax,"Hmin",Hmin,"Hgrad",1.2,...
          "Hedge",{[41:42],0.05},"GeometricOrder","linear");

% interpolate from slab2 onto mesh
[longrid,latgrid] = meshgrid(lon,lat);
interpobj = griddedInterpolant(longrid',latgrid',slabdepth',"nearest");
depth_interp = interpobj(meshobj.Nodes(1,:)',meshobj.Nodes(2,:)');

% deal with NaN values in depth
index = find(isnan(depth_interp));
% find triangles that contain the 'NaN' node and calculate a median value
for i = 1:length(index)
    triindex = (index(i) == meshobj.Elements(1,:)') | ...
        (index(i) == meshobj.Elements(2,:)')| (index(i) == meshobj.Elements(3,:)');
    nearestval = meshobj.Elements(:,find(triindex));
    depth_interp(index(i)) = nanmedian(depth_interp(nearestval(:)));
end

%% remove triangles that are deeper than a cutoff depth
index = find((abs(depth_interp) > abs(cutoff_depth)) |...
                (meshobj.Nodes(1,:)' > -93.5) | (meshobj.Nodes(1,:)' < -102.1));
deleterows = [];
for i = 1:length(index)
    triindex = (index(i) == meshobj.Elements(1,:)') | ...
        (index(i) == meshobj.Elements(2,:)')| (index(i) == meshobj.Elements(3,:)');
    deleterows = [deleterows;find(triindex)];
end

%% store triangles as [lon,lat,depth(km)]
p = [meshobj.Nodes' depth_interp + zshift];
t = meshobj.Elements';
t(deleterows,:) = [];
% convert from lat,lon to xyz
% lat0 = mean(p(:,2));
% lon0 = mean(p(:,1));

% convert everything to meters
[xnode, ynode, znode] = geodetic2enu(p(:,2),p(:,1),p(:,3).*1e3,...
    lat0,lon0,0,wgs84Ellipsoid);

if false
    figure(2),clf
    pdeplot(modelmesh), hold on
    scatter(p(:,1),p(:,2),100,p(:,3),'filled')
    plot(outlinedata{:,1},outlinedata{:,2},'k-','LineWidth',1)
    clim([-100 0])
    colorbar
    axis tight equal
    colormap(turbo(20))

    figure(3),clf
    trimesh(t,p(:,1),p(:,2),p(:,3)), hold on
    scatter3(p(:,1),p(:,2),p(:,3),50,p(:,3))
    clim([-100 0])
    axis tight
    view(8,75)
    colormap(turbo(10))
    zlim([-150 0])

    figure(4),clf
    trisurf(t,xnode,ynode,p(:,3).*1e3), hold on
    axis tight equal
    view(0,85)
    colormap(turbo(100))
    colorbar
    clim([cutoff_depth 0].*1e3)
end

figure(12),clf
pdegplot(modelmesh.Geometry,VertexLabels="off",EdgeLabels="on"), hold on
plot(modelmesh.Geometry.Vertices(:,1),modelmesh.Geometry.Vertices(:,2),'.')
set(gca,'Fontsize',20,'Linewidth',1.5), grid on, box on

%% write outputs to meshfile

% output in xyz in meters format for unicycle calculations
writetable(table([(1:length(xnode))',ynode,xnode,-p(:,3).*1e3]),...
    [outfile '_mt_xyz.ned'],'WriteVariableNames',false,"FileType",'text','Delimiter','\t')
writetable(table([(1:length(t(:,1)))',t,zeros(length(t(:,1)),1)])...
    ,[outfile '_mt_xyz.tri'],'WriteVariableNames',false,"FileType",'text','Delimiter','\t')

% output in lat,lon,km coodrinates for plotting
writetable(table([(1:length(xnode))',p(:,2),p(:,1),-p(:,3)]),...
    [outfile '_mt_latlon.ned'],'WriteVariableNames',false,"FileType",'text','Delimiter','\t')
writetable(table([(1:length(t(:,1)))',t,zeros(length(t(:,1)),1)])...
    ,[outfile '_mt_latlon.tri'],'WriteVariableNames',false,"FileType",'text','Delimiter','\t')


%% open using unicycle
eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver([outfile '_mt_xyz'],eM);

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

figure(11),clf
histogram(sqrt(2.*rcv.area./1e6),20)
xlabel('triangle length scale (km)'), ylabel('Frequency')
set(gca,'Fontsize',20,'Linewidth',1.5), grid on, box on


