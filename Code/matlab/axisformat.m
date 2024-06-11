function axisformat(xlabelname,ylabelname,zlabelname,titlename,needlegend,legendvector,legendname)
% set fontsize
fs=16; 
if needlegend==1
    leg=legend(legendvector,legendname,'interpreter','tex','fontsize',14,'location','northeast'); % 14
end
ax=gca;
ax.TickLabelInterpreter ='tex';
ax.FontSize=10; 
xlabel(xlabelname,'interpreter','tex','fontsize',fs);
ylabel(ylabelname,'interpreter','tex','fontsize',fs);
zlabel(zlabelname,'interpreter','tex','fontsize',fs);
title(titlename,'interpreter','tex','fontsize',14,'FontWeight','normal') 
set(gcf,'color','w');
end
