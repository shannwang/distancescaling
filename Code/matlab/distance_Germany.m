%--------------------------------------------------------------------------
% find scaling factors and draw distributions of distances for 16 states in
% Germany  
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

files={'nordrhein-westfalen_motorway.csv','bayern_motorway.csv','hessen_motorway.csv',...
    'schleswig-holstein_motorway.csv','baden-wuerttemberg_motorway.csv',...
    'mecklenburg-vorpommern_motorway.csv',...
    'niedersachsen_motorway.csv','brandenburg_motorway.csv',...
    'rheinland-pfalz_motorway.csv',...
    'sachsen-anhalt_motorway.csv','sachsen_motorway.csv','thueringen_motorway.csv',...
    'saarland_motorway.csv','hamburg_motorway.csv','berlin_motorway.csv','bremen_motorway.csv'};
tit={'Nordrhein-Westfalen','Bayern','Hessen','Schleswig-Holstein','Baden-Wuerttemberg',...
    'Mecklenburg-Vorpommern','Niedersachsen','Brandenburg','Rheinland-Pfalz',...
    'Sachsen-Anhalt','Sachsen','Thueringen','Saarland','Hamburg','Berlin','Bremen'};
tit2={'nordrhein-westfalen','bayern','hessen','schleswig-holstein','baden-wuerttemberg',...
    'mecklenburg-vorpommern','niedersachsen','brandenburg','rheinland-pfalz',...
    'sachsen-anhalt','sachsen','thueringen','saarland','hamburg','berlin','bremen'};
tit_en={'North Rhine-Westphalia','Bavaria','Hesse','Schleswig-Holstein','Baden-Württemberg',...
    'Mecklenburg-Western Pomerania','Lower Saxony','Brandenburg','Rhineland-Palatinate',...
    'Saxony-Anhalt','Saxony','Thuringia','Saarland','Hamburg','Berlin','Bremen'};


%% prepare data
for k=1:length(files)
    % read coordinates of each motorway network
    mw{k}=readtable([directory,files{k}]);  
    % set distance interval for selected locations on motorways
    interval(k)=floor(size(mw{k},1)/2000);
end

%% extract less data points and save them as csv files
for k=1:length(files)
    mw_short{k}=table();
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
    mw_short{k}=readtable([directory,tit{k},'_mw_less.csv']);      
end


%% calculate geodetic distance between sections
for k=1:length(files)
    i=0;
    for ii=1:1:size(mw_short{k},1)
        i=i+1;
        j=0;
        for jj=1:1:size(mw_short{k},1)
            j=j+1;
            lat1=mw_short{k}.lat(ii);
            lon1=mw_short{k}.lon(ii);
            lat2=mw_short{k}.lat(jj);
            lon2=mw_short{k}.lon(jj);
            distkm_geo{k}(i,j)=geodistkm([lat1 lon1],[lat2 lon2]);
        end
    end
    distkm_geo{k}(1:size(distkm_geo{k},2)+1:end)=0;
end
% save geodetic distance
% save('distkm_geo_Germany.mat','distkm_geo')

%% load geodetic distance between sections
load("distkm_geo_Germany.mat");

%% real network distance between sections
for k=1:length(files)
    distkm_net{k} = readmatrix([directory_nd,tit2{k},'_networkdistance_motorway.csv']);
end

%% filter coordinates for used sections which will be plotted on a map
for k=1:length(files)
    mw_short2{k}=mw_short{k};
    distkm_mat=distkm_net{k};
    % let diagonal elements zero to be NaN
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
        de=ceil(max(distkm_net{k},[],'all')/50);
        edges=[0:de:50*de];
        for ii=1:length(sf0)
            kk=0;
            % use the data in the low triangle of the matrix 
            for j=1:size(distkm_net{k},2)-1
                for i=j+1:size(distkm_net{k},1)
                    % if it is not a NaN value
                    if isnan(distkm_net{k}(i,j))==0
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
        clear idx_sf error edges de
end
% save scaling factor
% save('scalingfactors_16Germanstates.mat','sf')

%% load scaling factor
load('scalingfactors_16Germanstates.mat')

%% common edges of bins in histograms
for k=1:length(files)
        de(k)=ceil(max(distkm_net{k},[],'all')/50);
        edges=[0:de(k):50*de(k)];
        edges_all{k}=edges;
end

%% plot distribution without scaling factors
%close all
% distribution of distances
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
%colors={[1 0 0],[0 0.45 0.74]};
%fa=0.4;
ym={[0 12e-3],[0 6e-3],[0 11e-3],[0 2e-2],...
    [0 10e-3],[0 10e-3],[0 7e-3],[0 11e-3],...
    [0 10e-3],[0 15e-3],[0 10e-3],[0 15e-3],...
    [0 3e-2],[0 8e-2],[0 7e-2],[0 2e-1]};
figure(1)
tiledlayout(4,4,'Tilespacing','tight','Padding','compact')
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

    nexttile
    histogram(newdist_net,edges_all{k},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(newdist_geo,edges_all{k},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    if k==13
        axisformat(' ',' ','',tit_en{k},0,[],{})
    elseif k==5
       axisformat('','pdf','',tit_en{k},0,[],{})
    elseif k==14
       axisformat('distance [km]','','',tit_en{k},0,[],{})
    elseif k==16
       axisformat(' ','','',tit_en{k},1,[],{'l^{(n)}','l^{(g)}'})
    else   
        axisformat('','','',tit_en{k},0,[],{})
    end
    %title(tit_en{k})
    xlm{k}=xlim;
    ylim(ym{k});
    ax = gca;
    ytk{k}=ax.YTick;
    if k<=12
        ax.YAxis.Exponent = -3;
    else
        ax.YAxis.Exponent = -2;
    end
    ax.YTickMode='manual';
    clear newdist_net newdist_geo
end

%% plot distribution with scaling factors
% distribution of distances
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
figure(2)
tiledlayout(4,4,'Tilespacing','compact','Padding','compact')
for k=1:length(files)
    kk=0;
    for j=1:size(distkm_net{k},2)-1
        for i=j+1:size(distkm_net{k},1)
            if isnan(distkm_net{k}(i,j))~=1
                kk=kk+1;
                newdist_net(kk)=distkm_net{k}(i,j);
                newdist_geo_scale(kk)=distkm_geo{k}(i,j)*sf(k);
            end
        end
    end

    [N1,~]=histcounts(newdist_net,edges_all{k},'normalization','pdf');
    [N2,~]=histcounts(newdist_geo_scale,edges_all{k},'normalization','pdf');
    H(k)=HellingerDistance(N1,N2,de(k));  

    figure(2)
    nexttile
    histogram(newdist_net,edges_all{k},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(newdist_geo_scale,edges_all{k},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});

    if k==13
        axisformat(' ',' ','','',0,[],{});
    elseif k==5
       axisformat('','pdf','','',0,[],{});
    elseif k==14
       axisformat('distance [km]','','','',0,[],{})
    elseif k==15
       axisformat(' ','','','',1,[],{'l^{(n)}','\alphal^{(g)}'});
   else   
        axisformat('','','','',0,[],{});
    end
    xlim(xlm{k});
    ylim(ym{k});
    text(xlm{k}(2)*0.55,ym{k}(2)*0.8,['\alpha=',num2str(sf(k)),newline,'H=',num2str(round(H(k),3)) ],...
        'interpreter','tex','FontSize',12)
    title(tit_en{k})
    ax = gca;
    ax.YTick=ytk{k};
    if k<=12
        ax.YAxis.Exponent = -3;
    else
        ax.YAxis.Exponent = -2;
    end 
    ax.YTickMode='manual';
    clear newdist_net newdist_geo_scale a b
end






