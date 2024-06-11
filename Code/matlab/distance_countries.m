%--------------------------------------------------------------------------
% find scaling factors and draw distributions of distances for 9 regions 
% (including China, US, Germany, France, Spain, Great Britain, California, 
% Ontario, North Rhine-Westphalia)
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

files={'china_motorway.csv','france_motorway.csv',...
    'germany_motorway.csv','spain_motorway.csv','us_motorway.csv',...
    'california_motorway.csv','great-britain_motorway.csv',...
    'ontario_motorway.csv'}; 
tit={'China','France','Germany','Spain','US',...
     'California','Great-Britain','Ontario'}; 
tit2={'china','france','germany','spain','us',...
    'california','great-britain','ontario'}; 
%% prepare data
for k=1:length(files)
    % read coordinates of each motorway network
    mw{k}=readtable([directory,files{k}]);  
    % remove the coordinates in Alaska from the USA's mainland
    if strcmp(tit2{k},'us')==1
        idx=find(mw{k}.lon<-128);
        mw{k}(idx,:)=[];
        clear idx
    end
    % set distance interval for selected locations on motorways
    interval(k)=floor(size(mw{k},1)/2000);
end

%% extract less data points and save them as csv files
for k=1:length(files)
    mw_short{k}=table();
    warning('off')
    i=0;
    for ii=1:interval(k):size(mw{k},1)
        i=i+1;
        mw_short{k}.lat(i)=mw{k}.lat(ii);
        mw_short{k}.lon(i)=mw{k}.lon(ii);
    end
    %writetable(mw_short{k},[directory,tit{k},'_mw_less.csv'],'Delimiter',',')  
end

%% read the short list of motorway sections
clear mw_short
for k=1:length(files)
    mw_short{k}=readtable([directory,tit{k},'_mw_less.csv']) ;
end

%% calculate geodetic distance between sections
for k=1:length(files)
    i=0;
    for ii=1:interval(k):size(mw{k},1)
        i=i+1;
        j=0;
        for jj=1:interval(k):size(mw{k},1)
            j=j+1;
            lat1=mw{k}.lat(ii);
            lon1=mw{k}.lon(ii);
            lat2=mw{k}.lat(jj);
            lon2=mw{k}.lon(jj);
            % distance in km based on Haversine formula
            % (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
            distkm_geo{k}(i,j)=geodistkm([lat1 lon1],[lat2 lon2]);   
        end
    end
    distkm_geo{k}(1:size(distkm_geo{k},2)+1:end)=0;
end
%% save geodetic distance
% save('distkm_geo_countries.mat','distkm_geo')

%% load geodetic distance between sections
load("distkm_geo_countries.mat");

%% real network distance between sections
for k=1:length(files)
    distkm_net{k} = readmatrix([directory_nd,tit2{k},'_networkdistance_motorway.csv']);
end

%% filter coordinates for used sections
for k=1:length(files)
    mw_short2{k}=mw_short{k};
    distkm_mat=distkm_net{k};
    distkm_mat(1:size(distkm_mat,2)+1:end)=NaN;
    i_del=[];
    for i=1:size(distkm_net{k},1)
        if all(isnan(distkm_mat(i,:)),'all')==1 
           i_del=[i_del i];
        end
    end
    mw_short2{k}(i_del,:)=[];
    %writetable(mw_short2{k},[directory,tit2{k},'_mw_less2.csv'])
end

%% size of data points
datainfo=table();
for k=1:length(files)
    datainfo.location{k}=tit{k};
    datainfo.totalpoints(k)=size(mw{k},1);
    datainfo.usedpoints(k)=size(mw_short{k},1);   
    datainfo.finalusedpoints(k)=size(mw_short2{k},1);  
end

%% number of distances
for k=1:length(files)
    kk=0;
    for j=1:size(distkm_net{k},2)-1
        for i=j+1:size(distkm_net{k},1)
            if isnan(distkm_net{k}(i,j))~=1
                kk=kk+1;
                newdist_net(kk)=distkm_net{k}(i,j);
                newdist_geo(kk)=distkm_geo{k}(i,j);
            end
        end
    end
    num_netd(k,1)=length(newdist_net);
    num_geod(k,1)=length(newdist_geo);
    clear newdist_net newdist_geo
end

%% find the best scaling factor for pdf
sf0=[1:0.01:2];
for k=1:length(files)
        de(k)=ceil(max(distkm_net{k},[],'all')/50);
        edges=[0:de(k):50*de(k)];
        for ii=1:length(sf0)
            kk=0;
            for j=1:size(distkm_net{k},2)-1
                for i=j+1:size(distkm_net{k},1)
                    if isnan(distkm_net{k}(i,j))~=1
                        kk=kk+1;
                        newdist_net(kk)=distkm_net{k}(i,j);
                        newdist_geo_scale(kk)=distkm_geo{k}(i,j)*sf0(ii);
                    end
                end
            end
            [N_net,~]=histcounts(newdist_net,edges,'normalization','pdf');
            [N_geo_scale,~]=histcounts(newdist_geo_scale,edges,'normalization','pdf');
            error(ii)=sum((N_net-N_geo_scale).^2);
            clear newdist_net newdist_geo_scale N_net N_geo_scale
        end
        [~,idx_sf] = min(error,[],'all');
        sf(k)=sf0(idx_sf);
        edges_all{k}=edges;
        clear idx_sf error edges 
end
%save('scalingfactors_countries.mat','sf')

%% load scaling factors
for k=1:length(files)
        de(k)=ceil(max(distkm_net{k},[],'all')/50);
        edges=[0:de(k):50*de(k)];
        edges_all{k}=edges;
end
load('scalingfactors_countries.mat')

%% plot histogram without (upper row) and with (bottom row) scaling factors
close all
tit3={'China','France','Germany','Spain','USA',...
     'California','Great Britain','Ontario'};
% distribution of distances
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
ym={[0 6.5e-4],[0 2e-3],[0 2.5e-3],[0 2e-3],[0 5e-4],...
    [0 3e-3],[0 3.5e-3],[0 5e-3]};
figure(1)
t=tiledlayout(21,1,'Tilespacing','tight','Padding','compact')
t1=tiledlayout(t,1,5);
t1.Layout.Tile=1;
t1.Layout.TileSpan=[5,1];
t2=tiledlayout(t,1,5);
t2.Layout.Tile=6;
t2.Layout.TileSpan=[5,1];
t3=tiledlayout(t,1,10);
t3.Layout.Tile=12;
t3.Layout.TileSpan=[5,3];
t4=tiledlayout(t,1,10);
t4.Layout.Tile=17;
t4.Layout.TileSpan=[5,3];

for k=1:length(files)
    kk=0;
    for j=1:size(distkm_net{k},2)-1
        for i=j+1:size(distkm_net{k},1)
            if isnan(distkm_net{k}(i,j))~=1
                kk=kk+1;
                newdist_net(kk)=distkm_net{k}(i,j);
                newdist_geo(kk)=distkm_geo{k}(i,j);
                newdist_geo_scale(kk)=distkm_geo{k}(i,j)*sf(k);
            end
        end
    end
    [N1,~]=histcounts(newdist_net,edges_all{k},'normalization','pdf');
    [N2,~]=histcounts(newdist_geo_scale,edges_all{k},'normalization','pdf');
    H(k)=HellingerDistance(N1,N2,de(k)); 
    
    if k<=5
        nexttile(t1)
    else
        nexttile(t3)
    end
    histogram(newdist_net,edges_all{k},'normalization','pdf',...
        'FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(newdist_geo,edges_all{k},'normalization','pdf',...
        'FaceAlpha',fa,'FaceColor',colors{2});
    if k==5 || k==8
       axisformat('','','',tit3{k},1,[],{'l^{(n)}','l^{(g)}'})
    elseif k==6 
       axisformat('','pdf','',tit3{k},0,[],{})
    else   
        axisformat('','','',tit3{k},0,[],{})
    end
    ax=gca;
    %ax.YTickMode='manual';
    xlm{k}=xlim;
    ylim(ym{k});
    a=xlim;
    b=ylim;

    if k<=5
        nexttile(t2)
    else
        nexttile(t4)
    end
    histogram(newdist_net,edges_all{k},'normalization','pdf',...
        'FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(newdist_geo_scale,edges_all{k},'normalization','pdf',...
        'FaceAlpha',fa,'FaceColor',colors{2});
    if k==5 || k==8
       axisformat('','','','',1,[],{'l^{(n)}','\alphal^{(g)}'})
    elseif k==7
       axisformat('distance [km]','','','',0,[],{})
    elseif  k==1 
       axisformat('',' ','','',0,[],{})
    else   
        axisformat('','','','',0,[],{})
    end
    ax=gca;
    %ax.YAxis.Exponent = -3;
    %ax.YTickMode='manual';
    xlim(xlm{k});
    ylim(ym{k});
    if k<=5
        text(a(2)*0.45,b(2)*0.85,['\alpha=',num2str(sf(k)),newline,...
            'H=',num2str(round(H(k),3))],'interpreter','tex','FontSize',12)
    else
        text(a(2)*0.65,b(2)*0.85,['\alpha=',num2str(sf(k)),newline,...
            'H=',num2str(round(H(k),3))],'interpreter','tex','FontSize',12)
    end

    clear newdist_net newdist_geo newdist_geo_scale
end




