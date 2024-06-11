%--------------------------------------------------------------------------
% construct random grid networks
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% March 1, 2024 
%--------------------------------------------------------------------------
clc
clear all
close all

rng(3)
m=30; %row
n=30; % column
N=m*n; % number of nodes
p=[0.1:0.1:0.6];
k=0;
for i=1:m
    for j=1:n
        k=k+1;
        posx(k)=j;
        posy(k)=i;
        % map x to lon 5-10, and y to lat 50-55;
        xlon(k)=(10-5).*((j-0.5)/n)+5;
        ylat(k)=(55-50).*((i-0.5)/m)+50;
        z(k)=0.5;
    end
end
NodeTable = table(xlon',ylat',z','VariableNames',{'lon','lat','size'});
G=graph();
G=addnode(G,NodeTable);
for ip=1:length(p)
        g{ip}=GridNetwork(m,n,G,p(ip));
end

%% draw grid network, e.g. 8x8 grid network
figure(1)
tiledlayout(2,3,"TileSpacing","compact","Padding","compact")
for ip=1:length(p)
    nexttile
    [path,d] = shortestpath(g{ip},1,16);
    DrawGridNetwork(g{ip},path)
    if ip<4
        xticklabels([])
    end
    if rem(ip,3)~=1
        yticklabels([])
    end
    title(['f=',num2str(p(ip)),', l_{1,16}=',num2str(round(d,1)),' km'],'Interpreter','tex','FontSize',12)
    axis on
end

%% draw 8x8 grid network with f=0.2
figure(1)
tiledlayout(1,1,"TileSpacing","compact","Padding","loose")
for ip=2
    nexttile
    DrawGridNetwork2(g{ip})
    xlabel(' ','fontsize',16)
    axis on
end

%% number of edges
ne=@(x,y) (x-1)*y+(y-1)*x+(x-1)*(y-1)
ne(m,n)

%% calculate network and geodetic distances for 30 nodes
close all
clear netd geod
for ip=1:length(p)
    k=0;
    for i=1:N-1
        for j=i+1:N 
            [~,nd] = shortestpath(g{ip},i,j);
            if isinf(nd)==0
                k=k+1;
                netd{ip}(k)=nd;
                geod{ip}(k)=geodistkm([g{ip}.Nodes.lat(i),g{ip}.Nodes.lon(i)],...
                [g{ip}.Nodes.lat(j),g{ip}.Nodes.lon(j)]);
            end
        end
    end
end

%% number of distances
for ip=1:length(p)
    num_netd(ip)=length(netd{ip});
    num_geod(ip)=length(geod{ip});
end

%% find the best scaling factor for pdf
sf0=[0.5:0.01:3];
for ip=1:length(p)
        de(ip)=ceil(max(netd{ip},[],'all')/50);
        edges_all{ip}=[0:de(ip):50*de(ip)];
        for ii=1:length(sf0)
            [N_net,~]=histcounts(netd{ip},edges_all{ip},'normalization','pdf');
            [N_geo_scale,~]=histcounts(geod{ip}*sf0(ii),edges_all{ip},'normalization','pdf');
            error(ii)=sum((N_net-N_geo_scale).^2);
            clear N_net N_geo_scale
        end
        [~,idx_sf] = min(error,[],'all');
        sf(ip)=sf0(idx_sf);
        clear idx_sf error 
end


%% plot distributions of distances

% without scaling factors
% xl=[300 300 1200 1500 900 900];
% yl=[0.08 0.040 0.005 0.004 0.004 0.004];
% colors={[1 0 0],[0 0.45 0.74]}; 
% fa=0.4;
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
figure(1)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(p)
    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlm=xlim;
    ylm=ylim;
    xl(ip)=xlm(2);
    yl(ip)=ylm(2);
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax = gca;
    if ip<=2
        ax.YAxis.Exponent = -2;
    % else
    %     ax.YAxis.Exponent = -4;
    end

    if ip==4
        axisformat(' ','pdf','',['f=',num2str(p(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(p(ip))],0,[],{})
    elseif ip==6
        axisformat('','','',['f=',num2str(p(ip))],1,[],...
            {'l^{(n)}','l^{(g)}'})
    else
        axisformat('','','',['f=',num2str(p(ip))],0,[],{})
    end    
end

figure(2)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(p)
    [N1,~]=histcounts(netd{ip},edges_all{ip},'normalization','pdf');
    [N2,~]=histcounts(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf');
    H(ip)=HellingerDistance(N1,N2,de(ip));

    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax = gca;
    if ip<=2
        ax.YAxis.Exponent = -2;
    % else
    %     ax.YAxis.Exponent = -4;
    end

    if ip==4
        axisformat(' ','pdf','',['f=',num2str(p(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(p(ip))],0,[],{})
    elseif ip==6
        axisformat('','','',['f=',num2str(p(ip))],1,[],...
            {'l^{(n)}','\alphal^{(g)}'})
    else
        axisformat('','','',['f=',num2str(p(ip))],0,[],{})
    end

    text(xl(ip)*0.65,yl(ip)*0.8,['\alpha=',num2str(sf(ip)),...
         newline,'H=',num2str(round(H(ip),3))],...
        'interpreter','tex','FontSize',12)
end


%% draw distribution with f=0.2
close all
figure(1)
tiledlayout(2,1,'Tilespacing','compact','Padding','compact')
for ip=2
    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax=gca;
    ax.YAxis.Exponent=-2;
    xtk=xticks;
    axisformat(' ','pdf','','',1,[],{'l^{(n)}','l^{(g)}'})

    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax=gca;
    ax.YAxis.Exponent=-2;
    xticks=xtk;
    text(xl(ip)*0.65,yl(ip)*0.4,['\alpha=',num2str(sf(ip)),...
         newline,'H=',num2str(round(H(ip),3)) ],...
        'interpreter','tex','FontSize',14)
    axisformat('distance [km]','pdf','','',1,[],{'l^{(n)}','\alphal^{(g)}'})

end

%%
function G=GridNetwork(m,n,G,p)
k=0;
% horizontal edges
for i=1:m % row
    for j=1:n-1 % column
        k=k+1;
        % 1. Generate a graph with a rectangular grid topology.
        s(k)=(i-1)*n+j;
        t(k)=(i-1)*n+j+1;
    end
end
% vertical edges
for i=1:m-1
    for j=1:n
        k=k+1;
        s(k)=j+(i-1)*n;
        t(k)=j+i*n;
    end
end
% diagonal edges
for i=1:m
    for j=1:n
        if rem(i,2)==1 && rem(j,2)==1 
            if j==1 && i==1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j+1;
            elseif j==1 && i<m && i>1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j+1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j+1;
            elseif i==1 && j>1 && j<n
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j+1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j-1;
            elseif j==n && rem(n,2)==1 && i==1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j-1;
            elseif j==n && rem(n,2)==1 && i<m && i>1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j-1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j-1;
            elseif i==m && rem(m,2)==1 && j==1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j+1;
            elseif i==m && rem(m,2)==1 && j>1 && j<n
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j+1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j-1;
            elseif i==m && rem(m,2)==1 && j==n && rem(n,2)==1
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j-1;
            elseif i>1 && i<m && j>1 && j<n
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j+1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=i*n+j-1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j+1;
                k=k+1;
                s(k)=(i-1)*n+j;
                t(k)=(i-2)*n+j-1;
            end
        end
    end
end
r=rand(1,length(s));
pc=quantile(r,p);
idx=(r<=pc);
G=addedge(G,s(idx),t(idx));
for i=1:size(G.Edges.EndNodes,1)
    is=G.Edges.EndNodes(i,1);
    it=G.Edges.EndNodes(i,2);
    G.Edges.Weight(i)=geodistkm([G.Nodes.lat(is) G.Nodes.lon(is)],...
        [G.Nodes.lat(it) G.Nodes.lon(it)]);
end
end


function DrawGridNetwork(G,P)
for i=1:numnodes(G)
    plot(G.Nodes.lon(i),G.Nodes.lat(i),'ko','MarkerFaceColor','w','MarkerSize',...
        ceil(3+G.Nodes.size(i)*5))
    hold on
end
xlim([5-0.1 10+0.1])
ylim([50-0.1 55+0.1])
xticks([5:1:10])
xticklabels({'5^{o}E','6^{o}E','7^{o}E',...
    '8^{o}E','9^{o}E','10^{o}E'})
yticks([50:1:55])
yticklabels({'50^{o}N','51^{o}N','52^{o}N',...
    '53^{o}N','54^{o}N','55^{o}N'})
 axisformat('','','','',0,[],{})

for k=1:size(G.Edges,1)
    s=G.Edges.EndNodes(k,1);
    t=G.Edges.EndNodes(k,2);
    h1=plot([G.Nodes.lon(s) G.Nodes.lon(t)],[G.Nodes.lat(s) G.Nodes.lat(t)],'-',...
        'Color',[0 0.4470 0.7410],'LineWidth',2) ;
    h1.Color(4) = 0.6;
    hold on
    plot(G.Nodes.lon(s),G.Nodes.lat(s),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(s)*10))
    plot(G.Nodes.lon(t),G.Nodes.lat(t),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(t)*10))
end

for k=1:length(P)-1
    s=P(k);
    t=P(k+1);
    plot([G.Nodes.lon(s),G.Nodes.lon(t)],[G.Nodes.lat(s),G.Nodes.lat(t)],'r-','LineWidth',1.5)
end
end

function DrawGridNetwork2(G)
for i=1:numnodes(G)
    plot(G.Nodes.lon(i),G.Nodes.lat(i),'ko','MarkerFaceColor','w','MarkerSize',...
        ceil(3+G.Nodes.size(i)*5))
    hold on
end
xlim([5-0.1 10+0.1])
ylim([50-0.1 55+0.1])
xticks([5:1:10])
xticklabels({'5^{o}E','6^{o}E','7^{o}E',...
    '8^{o}E','9^{o}E','10^{o}E'})
yticks([50:1:55])
yticklabels({'50^{o}N','51^{o}N','52^{o}N',...
    '53^{o}N','54^{o}N','55^{o}N'})
 axisformat('','','','',0,[],{})

for k=1:size(G.Edges,1)
    s=G.Edges.EndNodes(k,1);
    t=G.Edges.EndNodes(k,2);
    h1=plot([G.Nodes.lon(s) G.Nodes.lon(t)],[G.Nodes.lat(s) G.Nodes.lat(t)],'-',...
        'Color',[0 0.4470 0.7410],'LineWidth',2) ;
    h1.Color(4) = 0.6;
    hold on
    plot(G.Nodes.lon(s),G.Nodes.lat(s),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(s)*10))
    plot(G.Nodes.lon(t),G.Nodes.lat(t),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(t)*10))
end
end
