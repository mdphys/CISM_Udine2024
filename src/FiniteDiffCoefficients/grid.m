%---- grids 

clear all
close all

M  = 20 ;
L  = 1.2 ;

h = L/M ;

grid1 = (0:M)*h ;
grid2 = (L^2 - ((0:M)*h).^2)/L ;
grid3 = [0;L*rand(18,1);L] ;


subplot(3,1,1)
plot(grid1,0*grid1,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
set(gca,TickLabelInterpreter = 'latex') ;
title('Uniform Grid','interpreter','latex') ;
set(gca,'YTick',[],'fontsize',14) ;
subplot(3,1,2) 
plot(grid2,0*grid2,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
set(gca,TickLabelInterpreter = 'latex') ;
title('Smooth Grid','interpreter','latex')
set(gca,'YTick',[],'fontsize',14) ;
subplot(3,1,3)
plot(grid3,0*grid3,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
set(gca,TickLabelInterpreter = 'latex') ;
title('Non-Smooth Grid','interpreter','latex')
set(gca,'YTick',[],'fontsize',14) ;
xlabel('$x$ (m)','interpreter','latex')





