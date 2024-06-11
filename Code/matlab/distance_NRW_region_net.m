%--------------------------------------------------------------------------
% check the influence of different choices of data points for NRW motorway
% networks  
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% March 1, 2024 
%--------------------------------------------------------------------------
clc
clear all
close all

dir='/Users/working/Documents/program/JupyterNotebook/TrafficDimension/';
nrw=readtable([dir,'regioninfo.csv']);
regionconnect=readmatrix([dir,'region_connect.csv']);
neighbour=readmatrix([dir,'region_neighbour.csv']);
n=size(nrw,1);

% number of connections in the region network
Nc=sum(regionconnect,'all')/2;
% total possible connections
Np=n*(n-1)/2;
% total possible neighbouring connections
Nn=sum(neighbour,'all')/2;
% fraction of empirical region connections
Fcn=Nc/Nn;
% average population density pd_ave
pd_ave=mean(nrw.pop_den);

%% make distance between two neighbouring regions as weight of an edge
regionconnect_weight=zeros(size(regionconnect));
for i=1:n
    for j=1:n
        if i<=j && regionconnect(i,j)==1
            regionconnect_weight(i,j)=geodistkm([nrw.lat(i) nrw.lon(i)],[nrw.lat(j) nrw.lon(j)]);
            regionconnect_weight(j,i)=regionconnect_weight(i,j);
        end
    end
end

G=graph(regionconnect_weight,nrw,'omitselfloops');
G_neig=graph(neighbour,nrw,'omitselfloops');

%% visualize all regions in a map
clear s
shapefile=[dir,'nrw_admin_level8_boundary.shp'];
s = shaperead(shapefile);
x1 = [s.X];
y1 = [s.Y];
mapshow(x1,y1,Color=[0.8 0.8 0.8],LineWidth=1)
hold on
scatter(nrw.lon,nrw.lat,20,'k','filled')
for kk=1:size(G.Edges,1)
    ss=G.Edges.EndNodes(kk,1);
    tt=G.Edges.EndNodes(kk,2);
    plot([nrw.lon(ss) nrw.lon(tt)],[nrw.lat(ss) nrw.lat(tt)],'-o','markersize',3,...
        'Color',[1.00,0.41,0.16],'LineWidth',1,'MarkerFaceColor',[1.00,0.41,0.16],...
        'MarkerEdgeColor',[1.00,0.41,0.16]) %
end
xticklabels([])
yticklabels([])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
box on
axis tight normal


%% calculate network and geodetic distances
k=0;
for i=1:n
    for j=1:n
        if i>j
            % seach shortest network distances
            % with Dijkstra algorithm that requires all edge weights to be nonnegative.
            [~,netd(i,j)] = shortestpath(G,i,j,'Method','positive');
            % calculate geodetic distances
            geod(i,j)=geodistkm([nrw.lat(i) nrw.lon(i)],[nrw.lat(j) nrw.lon(j)]);
            % transform a 2d matrix into a 1d vector
            k=k+1;
            netd_vec(k)=netd(i,j);
            geod_vec(k)=geod(i,j);
        end
    end
end
% remove the distances equal to infinity
netd_vec(isinf(netd_vec)==1)=[];
geod_vec(isinf(netd_vec)==1)=[];
% remove the distances equal to nan
netd_vec(isnan(netd_vec)==1)=[];
geod_vec(isnan(netd_vec)==1)=[];
% remove the distances equal to zero
netd_vec(netd_vec==0)=[];
geod_vec(netd_vec==0)=[];

% the average network distance nd_ave
nd_ave=mean(netd_vec);


%% plot distribution 
% plot distribution without scaling factor
colors={[1 0 0],[0 0.45 0.74]}; 
fa=0.4;
figure(1)
tiledlayout(1,2,"TileSpacing",'compact','Padding','compact')
nexttile
histogram(netd_vec,50,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(geod_vec,50,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlim([0 400])
ylim([0 0.008])
axisformat('\bf distance [km]','\bf pdf','','',1,[],{'$l^\mathrm{(n)}$','$l^\mathrm{(g)}$'})

% calculate scaling factor
clear ed
sf0=[0.01:0.01:10];
de=ceil(max(netd_vec,[],'all')/50);
ed=[0:de:50*de];
for ii=1:length(sf0)
    [N_net,~]=histcounts(netd_vec,ed,'normalization','pdf');
    [N_geo,~]=histcounts(geod_vec*sf0(ii),ed,'normalization','pdf');
    errors(ii)=sum((N_net-N_geo).^2);
    clear N_net N_geo
end
[~,idx_sf] = min(errors,[],'all');
sf=sf0(idx_sf);
H=HellDistForDistrs(netd_vec,geod_vec*sf,ed,de);

% plot distribution with scaling factors
nexttile
histogram(netd_vec,ed,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(geod_vec*sf,ed,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlim([0 400])
yl=ylim;
text(250, 0.4*(yl(1)+yl(2)),['$\alpha=',num2str(sf),'$', newline, '$H=',num2str(round(H,3)),'$'],'interpreter','latex','FontSize',14)
ylim([0 0.008])
%yticks(0:0.002:0.011)
axisformat(' ','','','',1,[],{'$l^\mathrm{(n)}$','$\alpha l^\mathrm{(g)}$'})







%% functions
% calclulate Hellinger distance between two distributions
function H=HellDistForDistrs(newdist_net,newdist_geo_scale,edges,de)
[N1,~]=histcounts(newdist_net,edges,'normalization','pdf');
[N2,~]=histcounts(newdist_geo_scale,edges,'normalization','pdf');
H=HellingerDistance(N1,N2,de);
end

