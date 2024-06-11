%--------------------------------------------------------------------------
% Distributions of geodetic distances for locations uniformaly distributed
% in a disk, and the empirical distribution of geodetic distances for
% China, France, Germany and North Rhine-Westphalia, fitted with analytic 
% equation (S4) in the Supplementary Information
%
% author: Shanshan Wang
% email: shanshan.wang@uni-due.de
% March 1, 2024 
%--------------------------------------------------------------------------


clc
clear all
close all

rng(10)
% define parameters
n=1000; % number of points
R=100;  % radius of the circle 

% generate uniformly distributed Cartesian coordinates x and y 
rg=rand(n,2);
xrg=rg(:,1);
yrg=rg(:,2);
x=(2*xrg-1)*R;
y=(2*yrg-1)*R;

% circle
theta=0:0.01*pi:2*pi;
X=R.*cos(theta);
Y=R.*sin(theta);

% filter out points
z=sqrt(x.^2+y.^2);
idx=find(z>R);
x(idx)=[];
y(idx)=[];

% plot
close all
figure('Color','w')
tiledlayout(1,2,'TileSpacing','loose','Padding','compact')
nexttile
fill(X,Y,[0 0.45 0.74],'FaceAlpha',0.2)
hold on
plot(x,y,'ko','MarkerFaceColor','k','MarkerSize',5)
grid on
axis([-100 100 -100 100],'equal') 
xlim([-100 100])
ylim([-100 100])
axis off


%% calculate Euclidean distance
clear distkm_geo xmati xmatj ymati ymatj
xmati=repmat(x,1,size(x,1));
xmatj=repmat(x',size(x,1),1);
ymati=repmat(y,1,size(y,1));
ymatj=repmat(y',size(y,1),1);
distkm_geo=sqrt((xmati-xmatj).^2+(ymati-ymatj).^2);

distkm=tril(distkm_geo,-1);
distkm=distkm(:);
distkm=distkm(distkm>0);


%% draw distribution
% plot
close all
figure('Color','w')
tiledlayout(1,2,'TileSpacing','loose','Padding','compact')
nexttile
fill(X,Y,[0 0.45 0.74],'FaceAlpha',0.2)
hold on
plot(x,y,'ko','MarkerFaceColor','k','MarkerSize',5)
axis([-100 100 -100 100],'equal') 
xlim([-100 100])
ylim([-100 100])
axis off

nexttile
histogram(distkm,50,'Normalization','pdf','FaceAlpha',0.4,'FaceColor',[0 0.45 0.74])
axisformat('distance [km]','pdf','','',0,[],{})


%% for gaussian distribution
clc
clear all
close all

p=@(l,sigma) l/(2*sigma.^2).*exp(-l.^2./(4.*sigma.^2));

%% load geodetic distance between sections
clear distkm_geo
tit={'China','France','Germany','North Rhine-Westphalia'}; 
distkm_geo_cou=load("distkm_geo_countries.mat");
distkm_geo_gem=load("distkm_geo_Germany.mat");
distkm_geo=distkm_geo_cou.distkm_geo;
distkm_geo{9}=distkm_geo_gem.distkm_geo{1};

tiledlayout(2,2,"TileSpacing","compact","Padding","compact")
idx_k=[1 2 3 9];
for ik=1:length(idx_k)
    clear newdist_geo y edges x sigma
    k=idx_k(ik);
    kk=0;
    for j=1:size(distkm_geo{k},2)-1
        for i=j+1:size(distkm_geo{k},1)
            if ~isnan(distkm_geo{k}(i,j))
                kk=kk+1;
                newdist_geo(kk)=distkm_geo{k}(i,j);
            end
        end
    end

    nexttile
    [y,edges]=histcounts(newdist_geo,50,'normalization','pdf');
    histogram(newdist_geo,edges,'normalization','pdf','FaceAlpha',0.4,'FaceColor',[0 0.45 0.74]);
    x=(edges(2:end)+edges(1:end-1))/2;
    s=(std(newdist_geo)-40):0.01:(std(newdist_geo)+200);
    clear dy
    for is=1:length(s)
        dy(is)=sum((p(x,s(is))-y).^2);
    end
    [~,idx_min]=min(dy);
    sigma=s(idx_min);
    hold on
    plot(x,p(x,sigma),'r-','LineWidth',2)
    if ik==1
        axisformat('',' ','',[tit{ik}, newline ...
            '\sigma_{fit}=',num2str(round(sigma,2)), newline...
            '\sigma_{emp}=',num2str(round(std(newdist_geo),2))],1,[],{'empirical','analytic'});
        legend('location','southeast')
    elseif ik==2 
        axisformat('','','',[tit{ik}, newline ...
            '\sigma_{fit}=',num2str(round(sigma,2)), newline...
            '\sigma_{emp}=',num2str(round(std(newdist_geo),2))],0,[],{});
    elseif ik==3 
        axisformat('geodetic distance [km]','pdf','',[tit{ik}, newline...
            '\sigma_{fit}=',num2str(round(sigma,2)), newline...
            '\sigma_{emp}=',num2str(round(std(newdist_geo),2))],0,[],{});
    else 
        axisformat(' ','','',[tit{ik},newline ...
            '\sigma_{fit}=',num2str(round(sigma,2)),newline...
            '\sigma_{emp}=',num2str(round(std(newdist_geo),2))],0,[],{});
        ax = gca;
        ax.YAxis.Exponent = -3;
    end
end



%%
function y=Theta(x)
for i=1:length(x)
    if x(i)>=0
        y(i)=1;
    else
        y(i)=0;
    end
end
end