%--------------------------------------------------------------------------
% construct partially and fully random motorway networks in NRW with
% different connection fractions, draw distributions of distances with
% and without scaling factors, and draw distributions related to 
% population densities
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% March 1, 2024 
%--------------------------------------------------------------------------

clc
clear all
close all
dir='/Users/working/Documents/program/JupyterNotebook/TrafficDimension/';
% all region information in NRW
nrw=readtable([dir,'regioninfo.csv']);
% marginal regions
mar_reg=readlines([dir,'region_boundary.csv']);
% regions belonging to ADE 4
reg_ADE4=readlines([dir,'region_name_ADE4.csv']);
% adjcency matrix of empirical region network
regionconnect=readmatrix([dir,'region_connect.csv']);
% adjcency matrix of fully connected network
neighbour=readmatrix([dir,'region_neighbour.csv']);

% number of region connections
Nc=sum(regionconnect,'all')/2;
% total possible neighbouring connections
Nn=sum(neighbour,'all')/2;
% fraction of empirical region connections
Fcn=Nc/Nn;
% numer of nodes
Nd=size(regionconnect,1);


%% empircal graph
for i=1:size(nrw,1)
    for j=1:size(nrw,1)
        geod(i,j)=geodistkm([nrw.lat(i) nrw.lon(i)],[nrw.lat(j) nrw.lon(j)]);
    end
end
logPiPj=(log(nrw.pop_den)*log(nrw.pop_den)');
logPiPj(1:size(logPiPj,2)+1:end)=NaN;
logPiPj_new=tril(logPiPj,-1);
sumlogPiPj=sum(logPiPj_new,'all','omitnan');
% probability of a pair i and j
logpij=logPiPj/sumlogPiPj;

PiPj=(nrw.pop_den*nrw.pop_den');
PiPj(1:size(PiPj,2)+1:end)=NaN;
PiPj_new=tril(PiPj,-1);
sumPiPj=sum(PiPj_new,'all','omitnan');
% probability of a pair i and j
pij=PiPj/sumPiPj;

% weight of a graph is in terms of distances
neig_weight=neighbour.*geod;
regionconnect_weight=regionconnect.*geod;

G_emp=graph(regionconnect_weight,nrw,'omitselfloops');
G_neig=graph(neig_weight,nrw,'omitselfloops');
[s_emp,t_emp] = findedge(G_emp);
[s_neig,t_neig] = findedge(G_neig);

[~,idx_pop,~]= intersect(nrw.name,reg_ADE4,'stable');
[~,idx_mar,~]= intersect(nrw.name,mar_reg,'stable');

% population density matrix
pd_ave=mean(nrw.pop_den);
pop_den_norm=G_emp.Nodes.pop_den/pd_ave;
pop_den_mat=pop_den_norm*pop_den_norm';


%% draw distributions related to population densities
tiledlayout(1,2,"TileSpacing","loose","Padding","compact")
nexttile
omega=pij;
histogram(omega,100,'normalization','pdf')
omega_vec=omega(:);
pd = fitdist(omega_vec,'Exponential');
y = pdf(pd,sort(omega_vec,'ascend'));
hold on
plot(sort(omega_vec,'ascend'),y,'r-','LineWidth',2);
axisformat('\omega_{ij}','pdf','','',1,[],{'empirical','exponential'})
xlim([0 0.0002])

nexttile
eta=logpij;
histogram(eta,100,'normalization','pdf')
eta_vec=eta(:);
pd = fitdist(eta_vec,'Lognormal');
y = pdf(pd,sort(eta_vec,'ascend'));
hold on
hh=plot(sort(eta_vec,'ascend'),y,'r-','LineWidth',2);
axisformat('\eta_{ij}','pdf','','',1,[],{'empirical','log-normal'})



%% construct random networks
clear g g_rand
tic
epsilon=logpij/max(logpij,[],'all','omitnan');
fc=[0.1:0.1:0.6];
for i=1:length(fc)
    rng(2);
    G_full=G_neig;
    g{i}=ConstructNetwork(G_full,nrw,fc(i),Nn,Nd,pij,epsilon);
    g_rand{i}=FullyRandNRWNet(G_full,nrw,fc(i),Nn,Nd);
end
toc
%%
save('graph_fraction.mat','g','g_rand');
%%
epsilon=logpij/max(logpij,[],'all','omitnan');
fc=[0.1:0.1:0.6];
load('graph_fraction.mat');



%% map information
shapefile=[dir,'nrw_admin_level8_boundary.shp'];
s = shaperead(shapefile);
X1 = [s.X];
Y1 = [s.Y];

%% draw networks
close all
tiledlayout(2,3,"TileSpacing","none","Padding","tight")
gg=g_rand; % for fully random network   
%gg=g;       % for partially random network
for i=1:length(fc)
    nexttile
    title(['f=',num2str(fc(i)),''],'interpreter','tex','fontsize',14,'FontWeight','normal')
    DrawNetwork(gg{i},X1,Y1)
end

%% for each network, calculate scaling factor and hellinger distances
tic
clear sf H G Netd_vec Geod_vec edges
for i=1:length(g)
    %G=g_rand{i}; % for fully random network
    G=g{i};       % for partially random network
    [Netd_vec{i},Geod_vec{i}]=CalGeodNetDistance(G,geod);
    [sf(i),H(i),edges{i}]=CalScalingFactor(Netd_vec{i},Geod_vec{i});
    clear G 
end
toc

%% plot distributions of distances
% without scaling factors
xl=[1000 1200 2000 2000 1700 1700]*2.5;
yl=[0.01 0.015 0.01 0.003 0.002 0.002];
% colors={[1 0 0],[0 0.45 0.74]}; 
% fa=0.4;
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
figure(1)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(fc)
    nexttile
    histogram(Netd_vec{ip},edges{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(Geod_vec{ip},edges{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlm=xlim;
    ylm=ylim;
    xl(ip)=xlm(2);
    yl(ip)=ylm(2);
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax = gca;
    if ip<=2
        ax.YAxis.Exponent = -2;
    else
        ax.YAxis.Exponent = -3;
    end
    % ax.YAxis.Exponent = -3;
    if ip==4
        axisformat(' ','pdf','',['f=',num2str(fc(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(fc(ip))],1,[],...
            {'l^{(n)}','l^{(g)}'})
    else
        axisformat('','','',['f=',num2str(fc(ip))],0,[],{})
    end    
end

% with scaling factors
figure(2)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(fc)
    nexttile
    histogram(Netd_vec{ip},edges{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(Geod_vec{ip}*sf(ip),edges{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax = gca;
    if ip<=2
        ax.YAxis.Exponent = -2;
    else
        ax.YAxis.Exponent = -3;
    end
    %ax.YAxis.Exponent = -3;

    if ip==4
        axisformat(' ','pdf','',['f=',num2str(fc(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(fc(ip))],1,[],...
            {'l^{(n)}','\alphal^{(g)}'})
    else
        axisformat('','','',['f=',num2str(fc(ip))],0,[],{})
    end

    text(xl(ip)*0.65,yl(ip)*0.8,['\alpha=',num2str(sf(ip)),...
         newline,'H=',num2str(round(H(ip),3)) ],...
        'interpreter','tex','FontSize',12)
end


%% network and distribution for f=0.2
close all
g_ful=g_rand{2}; % for fully random network   
g_par=g{2};       % for partially random network
i=2;

figure(1)
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
nexttile
DrawNetwork(g_ful,X1,Y1)

figure(2)
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
nexttile
DrawNetwork(g_par,X1,Y1)

%%
clear sf H G Netd_vec Geod_vec edges
G=g_par
[Netd_vec,Geod_vec]=CalGeodNetDistance(G,geod);
[sf,H,edges]=CalScalingFactor(Netd_vec,Geod_vec);
clear G

% plot distributions of distances
% without scaling factors
xl=[1000 1200 2000 2000 1700 1700]*2.5;
yl=[0.01 0.015 0.01 0.003 0.002 0.002];
% colors={[1 0 0],[0 0.45 0.74]}; 
% fa=0.4;
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
figure(1)
tiledlayout(1,2,'Tilespacing','compact','Padding','compact')

nexttile
histogram(Netd_vec,edges,'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(Geod_vec,edges,'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlm=xlim;
ylm=ylim;
xl=xlm(2);
yl=ylm(2);
xlim([0 xl]);
ylim([0 yl]);
ax = gca;
ax.YAxis.Exponent = -2;
axisformat('distance [km]','pdf','','',1,[],{'l^{(n)}','l^{(g)}'})


% with scaling factors
nexttile
histogram(Netd_vec,edges,'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
hold on
histogram(Geod_vec*sf,edges,'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
xlim([0 xl]);
ylim([0 yl]);
ax = gca;
ax.YAxis.Exponent = -2;
axisformat(' ','','','',1,[],{'l^{(n)}','\alphal^{(g)}'})
text(xl*0.65,yl*0.4,['\alpha=',num2str(sf),newline,'H=',num2str(round(H,3)) ],...
    'interpreter','tex','FontSize',14)










%%
function g=ConstructNetwork(G_full,nrw,fc,Nn,Nd,pij,epsilon)
    gmat=zeros(Nd,Nd);
    g=graph(gmat,nrw,'omitselfloops');
    nume=0;
    Nc=round(fc*Nn); 
    while nume<Nc
        nd_pair=randperm(Nd,2);
        if rand(1)*max(pij,[],'all','omitnan')<pij(nd_pair(1),nd_pair(2))
            sp= shortestpath(G_full,nd_pair(1),nd_pair(2));
            for k=1:length(sp)-1
                idx=findedge(G_full,sp(k),sp(k+1));
                wt=G_full.Edges.Weight(idx);
                if numedges(g)==0  
                    g=addedge(g,sp(k),sp(k+1),wt);
                    G_full.Edges.Weight(idx)=G_full.Edges.Weight(idx)*epsilon(nd_pair(1),nd_pair(2));
                else
                    % avoid repeatly add the same edge
                    idx_new=findedge(g,sp(k),sp(k+1));
                    if sum(idx_new)==0 
                        g=addedge(g,sp(k),sp(k+1),wt);
                        G_full.Edges.Weight(idx)=G_full.Edges.Weight(idx)*epsilon(nd_pair(1),nd_pair(2));
                    end   
                end
            end
        end
        nume= numedges(g);
        if numedges(g)==Nc 
           break;
        end
    end
    while nume>Nc
        [s_rm,t_rm]=findedge(g);
        deg_s=degree(g,s_rm);
        idx_rm=find(deg_s==1);
        if isempty(idx_rm)
            deg_t=degree(g,t_rm);
            idx_rm=find(deg_t==1);
        end
        i=randperm(length(idx_rm),1);
        g=rmedge(g,s_rm(idx_rm(i)),t_rm(idx_rm(i)));
        nume= numedges(g);
    end 
end

function g=FullyRandNRWNet(G_full,nrw,fc,Nn,Nd)
gmat=zeros(Nd,Nd);
g=graph(gmat,nrw,'omitselfloops');
Nc=round(fc*Nn);
Ntot=numedges(G_full);
idx=randperm(Ntot,Nc);
ed=G_full.Edges(idx,:);
g=addedge(g,ed);
end


function DrawNetwork(G,X1,Y1)
mapshow(X1,Y1,Color=[0.8 0.8 0.8],LineWidth=1)
hold on
%scatter(nrw.lon,nrw.lat,20,'k','filled')
for kk=1:size(G.Edges,1)
    ss=G.Edges.EndNodes(kk,1);
    tt=G.Edges.EndNodes(kk,2);
    h1=plot([G.Nodes.lon(ss) G.Nodes.lon(tt)],[G.Nodes.lat(ss) G.Nodes.lat(tt)],'-o',...
        'Color',[0 0.45 0.74],'LineWidth',1.5,'MarkerFaceColor',[0 0.45 0.74],...
        'MarkerEdgeColor',[0 0.45 0.74],'MarkerSize',2) ;
    h1.Color(4) = 0.5;
    hold on
end
box on
axis tight normal
axis off
set(gcf,'color','w')
end

function H=HellDistForDistrs(newdist_net,newdist_geo_scale,edges,de)
[N1,~]=histcounts(newdist_net,edges,'normalization','pdf');
[N2,~]=histcounts(newdist_geo_scale,edges,'normalization','pdf');
% calclulate Hellinger distance
H=sqrt(1-sqrt(N1)*sqrt(N2')*de);
end

function [Netd_vec,Geod_vec]=CalGeodNetDistance(G,geod)
% calculate the matrix of network distance
netd = distances(G);
% get the values from upper triangle of the distance matrix
Netd=triu(netd,1);
Geod=triu(geod,1);
idx=find(isinf(Netd)==0 & Netd~=0);
% remove the distances equal to infinity and equal to zero
Netd_vec=Netd(idx);
Geod_vec=Geod(idx);
end

function [sf,H,edges]=CalScalingFactor(Netd_vec,Geod_vec)
sf0=[0.01:0.01:3];
% calculate interval of each bin 
de=ceil(max(Netd_vec,[],'all')/50);
% determine upper and lower edges of all 50 bins
edges=[0:de:50*de];
% calculate pdfs of network and geodetic distance 
% and the difference of two kinds of pdfs 
for ii=1:length(sf0)
    [N_net,~]=histcounts(Netd_vec,edges,'normalization','pdf');
    [N_geo,~]=histcounts(Geod_vec*sf0(ii),edges,'normalization','pdf');
    errors(ii)=sum((N_net-N_geo).^2);
    clear N_net N_geo
end
% find the position of the minimal error
[~,idx_sf] = min(errors,[],'all');
% find the best scaling factor
sf=sf0(idx_sf);
% calclulate Hellinger distance H
H=HellDistForDistrs(Netd_vec,Geod_vec*sf,edges,de);
end




