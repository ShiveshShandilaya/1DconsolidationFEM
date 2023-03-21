function initial_figure = Plot_Initial(ncoord,eConn,le)

initial_figure=figure;
plot(ncoord(:,1),ncoord(:,2),'o');
title('initial condition');
grid on
set(gca,'ytick',[0:le(1):max(eConn(:,2))]);
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);