%--------------------------------------------------------------------------
% construct random geometric graphs with random nodes and random connections
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% Oct 29, 2024 
%
% Simulation procedure:
% 1. Randomly choose n nodes in a 2d plane.
% 2. Connect two nodes by an edge if its distance is less than or equal to
%    a threshold.
% 3. Calculate geodetic distance between two nodes with an edge
% 4. Search the shortest route between two nodes to obtain the network
%    distance
%--------------------------------------------------------------------------
% an example of network with 100 locations
% set random locations
clc
clear all
close all
rng(6)
n=100;
s=31;
t=45;
r=rand(2,n);
x=r(1,:);
y=r(2,:);
% map x and y to lat 50-55, lon 5-10;
xlon=(10-5).*x+5;
ylat=(55-50).*y+50;
z=0.5*ones(1,n);
NodeTable = table(xlon',ylat',z','VariableNames',{'lon','lat','size'});

% set random edges
diagdist=geodistkm([(55-50)*0+50,(10-5)*0+5],[(55-50)*1+50,(10-5)*1+5]);
pp=[0.001:0.001:1];
for ip=1:length(pp)
    p=pp(ip);
    A=zeros(n,n);
    for i=1:n-1
        for j=i+1:n
            A(i,j)=geodistkm([ylat(i),xlon(i)],[ylat(j),xlon(j)]);
            A(j,i)=A(i,j);
        end
    end
    
    thre=p*diagdist;
    A(A>thre)=0;
    Aset{ip}=A;
    fraction(ip)=sum(sum(A>0))/(n*(n-1));

    clear A p
end

tiledlayout(2,3,"TileSpacing","compact","Padding","compact")
for ip=1:6
    [~,idx_min]=min(abs(fraction-0.1*ip));
    p=fraction(idx_min);
    A=Aset{idx_min};
    G=graph(A,NodeTable,'omitselfloops');
    [P,d] = shortestpath(G,s,t);
    
    nexttile
    DrawNetwork(G,P,s,t);
    if ip<4
        xticklabels([])
    end
    if rem(ip,3)~=1
        yticklabels([])
    end
    title(['f=',num2str(round(p,1)),',  l_{',num2str(s),',',num2str(t),'}=',...
        num2str(round(d,2)),' km'],'Interpreter','tex','FontSize',12)
    clear G p A P d idx_min
end



%% draw network with f=0.2
clc
clear all
close all
rng(6)
n=100;
s=31;
t=45;
r=rand(2,n);
x=r(1,:);
y=r(2,:);
% map x and y to lat 50-55, lon 5-10;
xlon=(10-5).*x+5;
ylat=(55-50).*y+50;
z=0.5*ones(1,n);
NodeTable = table(xlon',ylat',z','VariableNames',{'lon','lat','size'});
% set random edges
diagdist=geodistkm([(55-50)*0+50,(10-5)*0+5],[(55-50)*1+50,(10-5)*1+5]);
pp=[0.001:0.001:1];
for ip=1:length(pp)
    p=pp(ip);
    A=zeros(n,n);
    for i=1:n-1
        for j=i+1:n
            A(i,j)=geodistkm([ylat(i),xlon(i)],[ylat(j),xlon(j)]);
            A(j,i)=A(i,j);
        end
    end
    thre=p*diagdist;
    A(A>thre)=0;
    Aset{ip}=A;
    fraction(ip)=sum(sum(A>0))/(n*(n-1));
    clear A p
end

[~,idx_min]=min(abs(fraction-0.2));
p=fraction(idx_min);
A=Aset{idx_min};
G=graph(A,NodeTable,'omitselfloops');
nexttile
DrawNetwork2(G);
xlabel(' ','fontsize',16)
clear G p q pc ic A P d





%% calculate network and geodetic distances for 100 nodes
clc
clear all
close all
rng(6)
n=100;
r=rand(2,n);
x=r(1,:);
y=r(2,:);
% map x and y to lat 50-55, lon 5-10;
xlon=(10-5).*x+5;
ylat=(55-50).*y+50;
z=0.5*ones(1,n);
NodeTable = table(xlon',ylat',z','VariableNames',{'lon','lat','size'});

% set random edges
diagdist=geodistkm([(55-50)*0+50,(10-5)*0+5],[(55-50)*1+50,(10-5)*1+5]);
pp=[0.001:0.001:1];
for ip=1:length(pp)
    p=pp(ip);
    A=zeros(n,n);
    for i=1:n-1
        for j=i+1:n
            A(i,j)=geodistkm([ylat(i),xlon(i)],[ylat(j),xlon(j)]);
            A(j,i)=A(i,j);
        end
    end
    thre=p*diagdist;
    A(A>thre)=0;
    Aset{ip}=A;
    fraction(ip)=sum(sum(A>0))/(n*(n-1));
    clear A p
end
clear pp

pp=[0.1:0.1:0.6];
for ip=1:length(pp)
    [~,idx_min]=min(abs(fraction-0.1*ip));
    p=fraction(idx_min);
    A=Aset{idx_min};
    g{ip}=graph(A,NodeTable,'omitselfloops');
end

for ip=1:length(pp)
    k=0;
    for i=1:n-1
        for j=i+1:n 
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
for ip=1:length(pp)
    num_netd(ip)=length(netd{ip});
    num_geod(ip)=length(geod{ip});
end


%% find the best scaling factor for pdf
sf0=[1:0.01:3];
for ip=1:length(pp)
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
close all
% colors={[1 0 0],[0 0.45 0.74]}; 
% fa=0.4;
colors={[1 0.41 0.16], [0.30 0.75 0.93] };
fa=0.5;
% without scaling factors
xl=[700 700 700 700 700 700];
yl=[0.004 0.004 0.004 0.004 0.004 0.004];
figure(1)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(pp)
    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    %ax = gca;
    %ax.YAxis.Exponent = -3;
    if ip==4
        axisformat(' ','pdf','',['f=',num2str(pp(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(pp(ip))],1,[],...
            {'l^{(n)}','l^{(g)}'})
    else
        axisformat('','','',['f=',num2str(pp(ip))],0,[],{})
    end    
end

figure(2)
tiledlayout(2,3,'Tilespacing','tight','Padding','compact')
for ip=1:length(pp)
    [N1,~]=histcounts(netd{ip},edges_all{ip},'normalization','pdf');
    [N2,~]=histcounts(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf');
    H(ip)=HellingerDistance(N1,N2,de(ip));

    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    %ax = gca;
    %ax.YAxis.Exponent = -3;
    if ip==4
        axisformat(' ','pdf','',['f=',num2str(pp(ip))],0,[],{})
    elseif ip==5
        axisformat('distance [km]','','',['f=',num2str(pp(ip))],1,[],...
            {'l^{(n)}','\alphal^{(g)}'})
    else
        axisformat('','','',['f=',num2str(pp(ip))],0,[],{})
    end

    text(xl(ip)*0.65,yl(ip)*0.8,['\alpha=',num2str(sf(ip)),...
         newline,'H=',num2str(round(H(ip),3)) ],...
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
    xtk=xticks;
    axisformat(' ','pdf','','',1,[],{'l^{(n)}','l^{(g)}'})

    nexttile
    histogram(netd{ip},edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{1});
    hold on
    histogram(geod{ip}*sf(ip),edges_all{ip},'normalization','pdf','FaceAlpha',fa,'FaceColor',colors{2});
    xlim([0 xl(ip)]);
    ylim([0 yl(ip)]);
    ax=gca;
    xticks=xtk;
    text(xl(ip)*0.65,yl(ip)*0.4,['\alpha=',num2str(sf(ip)),...
         newline,'H=',num2str(round(H(ip),3)) ],...
        'interpreter','tex','FontSize',14)
    axisformat('distance [km]','pdf','','',1,[],{'l^{(n)}','\alphal^{(g)}'})

end



    
  
    

%%
function DrawNetwork(G,P,s0,t0)
for i=1:numnodes(G)
    plot(G.Nodes.lon(i),G.Nodes.lat(i),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(i)*10))
    % if i==s || i==t
    %     text(G.Nodes.lon(i)+0.1,G.Nodes.lat(i),num2str(i),'Interpreter','tex',...
    %         'FontSize',12)
    % end
    hold on
end
xlim([5-0.2 10+0.2])
ylim([50-0.2 55+0.2])
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
        'Color',[0 0.4470 0.7410],'LineWidth',1) ;
    h1.Color(4) = 0.6;
    hold on
end

for k=1:length(P)-1
    s=P(k);
    t=P(k+1);
    plot([G.Nodes.lon(s),G.Nodes.lon(t)],[G.Nodes.lat(s),G.Nodes.lat(t)],'r-','LineWidth',1.5)
end

for i=1:numnodes(G)
    if i==s0 || i==t0
        text(G.Nodes.lon(i)+0.1,G.Nodes.lat(i),num2str(i),'Interpreter','tex',...
            'FontSize',12)
    end
    hold on
end

end

function DrawNetwork2(G)
for i=1:numnodes(G)
    plot(G.Nodes.lon(i),G.Nodes.lat(i),'ko','MarkerFaceColor','b','MarkerSize',...
        ceil(3+G.Nodes.size(i)*10))
    % text(G.Nodes.lon(i)+0.1,G.Nodes.lat(i),num2str(i),'Interpreter','tex',...
    %     'FontSize',12)
    hold on
end
xlim([5-0.2 10+0.2])
ylim([50-0.2 55+0.2])
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
        'Color',[0 0.4470 0.7410],'LineWidth',1) ;
    h1.Color(4) = 0.6;
    hold on
end
end

function H=HellingerDistance(p,q,de)
H=sqrt(1-sqrt(p)*sqrt(q')*de);
end
