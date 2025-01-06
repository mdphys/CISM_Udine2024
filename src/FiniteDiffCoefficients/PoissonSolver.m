%------------------------------------------------------------------------
%                   Approximation of Derivatives
%          Task: Solve Poisson's Equation using Grids
%                       Dr M Ducceschi
%                   University of Bologna
%                       11 May 2024
%-------------------------------------------------------------------------


clear all
close all
clc



M  = 20 ;
L  = 1.2 ;

h = L/M ;

grid1   = (1:M-1)'*h ;
grid2   = (L^2 - ((1:M-1)'*h).^2)/L ;
grid2   = grid2(end:-1:1) ;
grid3   = L*rand(M-1,1) ;
grid3   = sort(grid3) ;
gridEx  = (0:500)'*L/500 ;

w1    = -(1-tanh(grid1).^2) ;
w2    = -(1-tanh(grid2).^2) ;
w3    = -(1-tanh(grid3).^2) ;

s1    = -log(cosh(grid1)) + L^(-1)*grid1*log(cosh(L))  ;
s2    = -log(cosh(grid2)) + L^(-1)*grid2*log(cosh(L))  ;
s3    = -log(cosh(grid3)) + L^(-1)*grid3*log(cosh(L))  ;

solEx    = -log(cosh(gridEx)) + L^(-1)*gridEx*log(cosh(L))  ;

D1 = zeros(M-1,M-1) ;
D2 = zeros(M-1,M-1) ;
D3 = zeros(M-1,M-1) ;

for m = 2 : M - 2

    D1m       = grid1(m)-grid1(m-1) ; 
    D1p       = grid1(m+1)-grid1(m) ;
    D1(m,m)   = -2/(D1m*D1p) ;
    D1(m,m+1) = 2/(D1p)/(D1m+D1p) ;
    D1(m,m-1) = 2/(D1m)/(D1m+D1p) ;


    D2m       = grid2(m)-grid2(m-1) ; 
    D2p       = grid2(m+1)-grid2(m) ;
    D2(m,m)   = -2/(D2m*D2p) ;
    D2(m,m+1) = 2/(D2p)/(D2m+D2p) ;
    D2(m,m-1) = 2/(D2m)/(D2m+D2p) ;


    D3m       = grid3(m)-grid3(m-1) ; 
    D3p       = grid3(m+1)-grid3(m) ;
    D3(m,m)   = -2/(D3m*D3p) ;
    D3(m,m+1) = 2/(D3p)/(D3m+D3p) ;
    D3(m,m-1) = 2/(D3m)/(D3m+D3p) ;


end

D1m           = grid1(1) ; 
D1p           = grid1(2)-grid1(1) ;
D1(1,1)       = -2/(D1m*D1p) ;
D1(1,2)       = 2/(D1p*(D1m+D1p)) ;
D1m           = grid1(M-1)-grid1(M-2) ;
D1p           = L-grid1(M-1) ;
D1(end,end)   = -2/(D1m*D1p) ;
D1(end,end-1) = 2/(D1m*(D1m+D1p)) ;


D2m           = grid2(1) ; 
D2p           = grid2(2)-grid2(1) ;
D2(1,1)       = -2/(D2m*D2p) ;
D2(1,2)       = 2/(D2p*(D2m+D2p)) ;
D2m           = grid2(M-1)-grid2(M-2) ;
D2p           = L-grid2(M-1) ;
D2(end,end)   = -2/(D2m*D2p) ;
D2(end,end-1) = 2/(D2m*(D2m+D2p)) ;


D3m           = grid3(1) ; 
D3p           = grid3(2)-grid3(1) ;
D3(1,1)       = -2/(D3m*D3p) ;
D3(1,2)       = 2/(D3p*(D3m+D3p)) ;
D3m           = grid3(M-1)-grid3(M-2) ;
D3p           = L-grid3(M-1) ;
D3(end,end)   = -2/(D3m*D3p) ;
D3(end,end-1) = 2/(D3m*(D3m+D3p)) ;


s1Approx      = D1 \ w1 ;
s2Approx      = D2 \ w2 ;
s3Approx      = D3 \ w3 ;




subplot(3,1,1)
%plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
plot([0;grid1;L],[0;s1;0],'-b') ; hold on;
plot([0;grid1;L],[0;s1Approx;0],'--k') ; 
plot(grid1,0*grid1,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
set(gca,TickLabelInterpreter = 'latex') ;
title('Uniform Grid','interpreter','latex') ;
ylim([-0.05,0.15])
set(gca,'fontsize',12)
legend('Exact', 'Num.')
subplot(3,1,2)
%plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
plot([0;grid2;L],[0;s2;0],'-b') ; hold on;
plot([0;grid2;L],[0;s2Approx;0],'--r') ; 

plot(grid2,0*grid2,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
set(gca,TickLabelInterpreter = 'latex') ;
title('Smooth Grid','interpreter','latex') ;
set(gca,'fontsize',12)
legend('Exact', 'Num.')
ylim([-0.05,0.15])
subplot(3,1,3)
%plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
plot([0;grid3;L],[0;s3;0],'-b') ; hold on ;
plot([0;grid3;L],[0;s3Approx;0],'--g') ; 
plot(grid3,0*grid3,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
set(gca,TickLabelInterpreter = 'latex') ;
title('Random Grid','interpreter','latex') ;
ylim([-0.05,0.15])
set(gca,'fontsize',12)
legend('Exact', 'Num.')



