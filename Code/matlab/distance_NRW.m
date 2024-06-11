%--------------------------------------------------------------------------
% find scaling factors and draw distributions of distances for
% the NRW motorway network and check the influence of different choices of 
% data points for NRW motorway networks
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% March 1, 2024 
%--------------------------------------------------------------------------

clc
clear all
close all

% directory of coordinates at motorway networks
directory='/Users/working/Documents/data/Traffic/map/filted_motorway_coordinates/';
% directory of network distances
directory_nd='/Users/working/Documents/program/MATLAB/TrafficDimension/data_networkdistance/networkdistance/';

files='nordrhein-westfalen_motorway.csv';
tit='Nordrhein-Westfalen';
tit2='nordrhein-westfalen';

%% prepare data
% read coordinates of each motorway network
mw=readtable([directory,files]);
% set distance interval for selected locations on motorways
interval=floor(size(mw,1)/2000);

%% read the short list of motorway sections
clear mw_short
mw_short=readtable([directory,tit,'_mw_less.csv']) ;     

%% load geodetic distance between sections
i=0;
for ii=1:1:size(mw_short,1)
    i=i+1;
    j=0;
    for jj=1:1:size(mw_short,1)
        j=j+1;
        lat1=mw_short.lat(ii);
        lon1=mw_short.lon(ii);
        lat2=mw_short.lat(jj);
        lon2=mw_short.lon(jj);
        distkm_geo(i,j)=geodistkm([lat1 lon1],[lat2 lon2]);
    end
end
distkm_geo(1:size(distkm_geo,2)+1:end)=0;


%% real network distance between sections
distkm_net = readmatrix([directory_nd,tit2,'_networkdistance_motorway.csv']);


%% filter coordinates for used sections
mw_short2=mw_short;
distkm_mat=distkm_net;
distkm_mat(1:size(distkm_mat,2)+1:end)=NaN;
i_del=[];
for i=1:size(distkm_net,1)
    if all(isnan(distkm_mat(i,:)),'all')==1
        i_del=[i_del i];
    end
end
mw_short2(i_del,:)=[];
%writetable(mw_short2,[directory,tit2,'_mw_less2.csv'])

distkm_net(i_del,:)=[];
distkm_net(:,i_del)=[];
distkm_geo(i_del,:)=[];
distkm_geo(:,i_del)=[];

%% size of data points
datainfo=table();
datainfo.location=tit;
datainfo.totalpoints=size(mw,1);
datainfo.usedpoints=size(mw_short,1);
datainfo.finalusedpoints=size(mw_short2,1);

%% ------------------------------------------------------------------------
% distributions of distances without and with scaling factors
% for NRW motorway networks 
% -------------------------------------------------------------------------

k=0;
n=size(distkm_net,1);
for i=1:n
    for j=1:n
        if i>j
            % transform a 2d matrix into a 1d vector
            k=k+1;
            netd_vec(k)=distkm_net(i,j);
            geod_vec(k)=distkm_geo(i,j);
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

% find the best scaling factor for pdf
sf0=[1:0.01:2];
de=ceil(max(distkm_net,[],'all')/50);
edges=[0:de:50*de];
for ii=1:length(sf0)
    [N_net,~]=histcounts(netd_vec,edges,'normalization','pdf');
    [N_geo_scale,~]=histcounts(geod_vec*sf0(ii),edges,'normalization','pdf');
    error(ii)=sum((N_net-N_geo_scale).^2);
    clear newdist_net newdist_geo_scale N_net N_geo_scale
end
[~,idx_sf] = min(error,[],'all');
sf=sf0(idx_sf);
H=HellDistForDistrs(netd_vec,geod_vec*sf,edges,de);

%% plot distribution without scaling factor
% colors={[1 0 0],[0 0.45 0.74]}; 
% fa=0.4;
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
figure(1)
tiledlayout(1,2,"TileSpacing",'compact','Padding','compact')
nexttile
histogram(netd_vec,edges,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(geod_vec,edges,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlim([0 400])
ylim([0 0.01])
axisformat('distance [km]','pdf','','',1,[],{'l^{(n)}','l^{(g)}'})
ax=gca;
ax.YAxis.Exponent=-3;

% plot distribution with scaling factors
nexttile
histogram(netd_vec,edges,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(geod_vec*sf,edges,'Normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlim([0 400])
yl=ylim;
text(250, 0.4*(yl(1)+yl(2)),['\alpha=',num2str(sf), newline, 'H=',num2str(round(H,3))],'interpreter','tex','FontSize',14)
ylim([0 0.01])
axisformat(' ','','','',1,[],{'l^{(n)}','\alphal^{(g)}'})
ax=gca;
ax.YAxis.Exponent=-3;




%% ------------------------------------------------------------------------
% check the influence of choices of data points
% shuffle the indices of sections and randomly select 2000 sections. 
% see how the scaling factors and Hellinger distances change with different
% selections
% -------------------------------------------------------------------------

%% shuffle location and select 2000 data points
% only run for once
% directory='/Users/working/Documents/data/Traffic/map/';
% for k=1:100
%     rng('shuffle')
%     idx_shu = randperm(size(mw,1),2000);
%     mw_short_shuffle{k}=mw(idx_shu,:);
%     writetable(mw_short_shuffle{k},...
%         [directory,'NRW_shuffle/NRW_mw_less_shuffle_',num2str(k),'.csv'],...
%         'Delimiter',',')
% end

%% calculating scaling factors and Hellinger distances
clc
clear all
dir1='/Users/working/Documents/data/Traffic/map/NRW_shuffle/';
dir2='/Users/working/Documents/program/MATLAB/TrafficDimension/data_networkdistance/NRW_shuffle/';
sf0=[1:0.01:2];
for k=1:100
    % read the short list of motorway sections
    mw_short_shuffle=readtable([dir1,'NRW_mw_less_shuffle_',num2str(k),'.csv']);
    % calculate the geodetic distances
    dist_geo_mat=CalGeodeticDistance(mw_short_shuffle);
    % load the network distances
    dist_net_mat=readmatrix([dir2,'NRW_networkdistance_motorway_shuffle_',num2str(k),'.csv']);
    % determine intervals and edges for each histogram 
    de=ceil(max(dist_net_mat,[],'all')/50);
    edges=0:de:50*de;
    % find the best scaling factor
    sf(k)=FindBestScalingFactor(sf0,dist_net_mat,dist_geo_mat,edges);
    % transform matrices to vectors
    kk=0;
     for j=1:size(dist_net_mat,2)-1
        for i=j+1:size(dist_net_mat,1)
            if isnan(dist_net_mat(i,j))~=1
                kk=kk+1;
                dist_net_vec(kk)=dist_net_mat(i,j);
                dist_geo_scale_vec(kk)=dist_geo_mat(i,j)*sf(k);
            end
        end
     end
    % calculate the Hellinger distance
    H(k)=HellDistForDistrs(dist_net_vec,dist_geo_scale_vec,edges,de);

    % clear temporary data to release physical memory
    clear mw_short_shuffle dist_geo_mat dist_net_mat de edges 
    clear dist_net_vec dist_geo_scale_vec
end
%% save data
%save('scalingfactors_NRW_shuffle_2000points.mat','sf','H')
%% load data
clc
clear all
load('scalingfactors_NRW_shuffle_2000points.mat')

%% plot
colors={[0.85 0.33 0.10],[0 0.45 0.74]}; 
tiledlayout(2,1,"TileSpacing","none","Padding","compact")
nexttile
stem(sf,'ko','MarkerFaceColor',colors{1},'BaseValue',1.35,'MarkerSize',8)
axisformat('','\alpha','','',0,[],{})
ylim([1.325 1.38])
yticks([1.34:0.02:1.38])
xticklabels([])
grid minor 
grid on
nexttile
stem(H,'ko','MarkerFaceColor',colors{2},'BaseValue',0.044,'MarkerSize',8)
grid minor 
grid on
axisformat('k-th selection (arbitrary k)','H','','',0,[],{})
ylim([0.03 0.054])
yticks([0.03:0.01:0.05])






%% ------------------------------------------------------------------------
% check the influence of numbers of data points
% shuffle the indices of sections and randomly select 200, 500, ...,10000 sections
% see how the scaling factors and Hellinger distances change with different
% number of sections
% -------------------------------------------------------------------------

% only run for once
% directory='/Users/working/Documents/data/Traffic/map/';
% nums=[200 500:500:10000];
% for k=1:length(nums)
%     rng('shuffle')
%     idx_num = randperm(size(mw,1),nums(k));
%     mw_number_shuffle{k}=mw(idx_num,:);
%     writetable(mw_number_shuffle{k},...
%         [directory,'NRW_number/NRW_mw_less_shuffle_number_',num2str(nums(k)),'.csv'],...
%         'Delimiter',',')
% end

%% calculating scaling factors and Hellinger distances
clc
clear all
dir1='/Users/working/Documents/data/Traffic/map/NRW_number/';
dir2='/Users/working/Documents/program/MATLAB/TrafficDimension/data_networkdistance/NRW_number/';
sf0=[1:0.01:2];
nums=[200 500:500:10000];
for k=1:length(nums)
    % read the short list of motorway sections
    mw_number_shuffle=readtable([dir1,'NRW_mw_less_shuffle_number_',num2str(nums(k)),'.csv']);
    % calculate the geodetic distances
    dist_geo_mat=CalGeodeticDistance(mw_number_shuffle);
    % load the network distances
    dist_net_mat=readmatrix([dir2,'NRW_networkdistance_motorway_shuffle_number_',num2str(nums(k)),'.csv']);
    % determine intervals and edges for each histogram 
    de=ceil(max(dist_net_mat,[],'all')/50);
    edges=0:de:50*de;
    % find the best scaling factor
    sf(k)=FindBestScalingFactor(sf0,dist_net_mat,dist_geo_mat,edges);
    % transform matrices to vectors
    kk=0;
     for j=1:size(dist_net_mat,2)-1
        for i=j+1:size(dist_net_mat,1)
            if isnan(dist_net_mat(i,j))~=1
                kk=kk+1;
                dist_net_vec(kk)=dist_net_mat(i,j);
                dist_geo_scale_vec(kk)=dist_geo_mat(i,j)*sf(k);
            end
        end
     end
    % calculate the Hellinger distance
    H(k)=HellDistForDistrs(dist_net_vec,dist_geo_scale_vec,edges,de);

    % clear temporary data to release physical memory
    clear mw_number_shuffle dist_geo_mat dist_net_mat de edges 
    clear dist_net_vec dist_geo_scale_vec
end

%% save data
%save('scalingfactors_NRW_shuffle_number.mat','sf','H')
%% load data
clc
clear all
nums=[200 500:500:10000];
load('scalingfactors_NRW_shuffle_number.mat')


%% plot
colors={[0.85 0.33 0.10],[0 0.45 0.74]}; 
tiledlayout(2,1,"TileSpacing","none","Padding","compact")
nexttile
stem(nums,sf,'ko','MarkerFaceColor',colors{1},'BaseValue',1.35,'MarkerSize',8)
axisformat('','\alpha','','',0,[],{})
ylim([1.325 1.38])
yticks([1.34:0.02:1.38])
xticklabels([])
grid minor 
grid on
nexttile
stem(nums,H,'ko','MarkerFaceColor',colors{2},'BaseValue',0.044,'MarkerSize',8)
grid minor 
grid on
axisformat('number of selected sections','H','','',0,[],{})
ylim([0.035 0.061])







 %% functions
 % calculate geodetic distance between sections
 function distkm_geo=CalGeodeticDistance(mwlist)
 for i=1:size(mwlist,1)-1
     for j=i+1:1:size(mwlist,1)
         lat1=mwlist.lat(i);
         lon1=mwlist.lon(i);
         lat2=mwlist.lat(j);
         lon2=mwlist.lon(j);
         distkm_geo(i,j)=geodistkm([lat1 lon1],[lat2 lon2]);
         distkm_geo(j,i)=distkm_geo(i,j);
     end
 end
 distkm_geo(1:size(distkm_geo,2)+1:end)=0;
 end

% find the best scaling factor
function sf=FindBestScalingFactor(sf0,dist_net_mat,dist_geo_mat,edges)
for ii=1:length(sf0)
    kk=0;
    for j=1:size(dist_net_mat,2)-1
        for i=j+1:size(dist_net_mat,1)
            if isnan(dist_net_mat(i,j))~=1
                kk=kk+1;
                newdist_net(kk)=dist_net_mat(i,j);
                newdist_geo_scale(kk)=dist_geo_mat(i,j)*sf0(ii);
            end
        end
    end
    [N_net,~]=histcounts(newdist_net,edges,'normalization','pdf');
    [N_geo_scale,~]=histcounts(newdist_geo_scale,edges,'normalization','pdf');
    error(ii)=sum((N_net-N_geo_scale).^2);
    clear newdist_net newdist_geo_scale N_net N_geo_scale
end
[~,idx_sf] = min(error,[],'all');
sf=sf0(idx_sf);
clear idx_sf error
end

% calclulate Hellinger distance between two distributions
function H=HellDistForDistrs(newdist_net,newdist_geo_scale,edges,de)
[N1,~]=histcounts(newdist_net,edges,'normalization','pdf');
[N2,~]=histcounts(newdist_geo_scale,edges,'normalization','pdf');
H=HellingerDistance(N1,N2,de);
end


