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

Ntest  = 20 ;
Mvec   = logspace(1,3,Ntest) ;
errvec = zeros(4,Ntest) ;
L      = 1.2 ;

for nT = 1 : Ntest

    M  = round(Mvec(nT)) ;


    h = L/M ;

    grid1   = (1:M-1)'*h ;
    grid2   = (L^2 - ((1:M-1)'*h).^2)/L ;
    grid2   = grid2(end:-1:1) ;
    grid3   = L*rand(M-1,1) ;
    grid3   = sort(grid3) ;

    w1    = -(1-tanh(grid1).^2) ;
    w2    = -(1-tanh(grid2).^2) ;
    w3    = -(1-tanh(grid3).^2) ;

    s1    = -log(cosh(grid1)) + L^(-1)*grid1*log(cosh(L))  ;
    s2    = -log(cosh(grid2)) + L^(-1)*grid2*log(cosh(L))  ;
    s3    = -log(cosh(grid3)) + L^(-1)*grid3*log(cosh(L))  ;


    D1 = zeros(M-1,M-1) ;
    D2 = zeros(M-1,M-1) ;
    D3 = zeros(M-1,M-1) ;

    v                 = ones(M-1,1) ;
    D1fourth          = -2.5*diag(v) + 4/3*diag(v(1:end-1),-1) + 4/3*diag(v(1:end-1),1) -1/12*diag(v(1:end-2),2) -1/12*diag(v(1:end-2),-2) ;
    D1fourth(1,1)     = -2 ; D1fourth(1,2) = 1 ; D1fourth(1,3) = 0 ;
    D1fourth(end,end) = -2 ; D1fourth(end,end-1) = 1 ; D1fourth(end,end-2) = 0 ;
    D1fourth = D1fourth / h^2 ;

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


    s1fourthApprox      = D1fourth \ w1 ;
    s1Approx            = D1 \ w1 ;
    s2Approx            = D2 \ w2 ;
    s3Approx            = D3 \ w3 ;

    errvec(1,nT) = sqrt(((s1Approx-s1)'*(s1Approx-s1))/(M-1)) ;
    errvec(2,nT) = sqrt(((s2Approx-s2)'*(s2Approx-s2))/(M-1)) ;
    errvec(3,nT) = sqrt(((s3Approx-s3)'*(s3Approx-s3))/(M-1)) ;
    errvec(4,nT) = sqrt(((s1fourthApprox-s1)'*(s1fourthApprox-s1))/(M-1)) ;


    if nT == 2
        subplot(3,1,1)
        %plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
        plot([0;grid1;L],[0;s1Approx;0]-[0;s1;0],'k') ; hold on;
        plot(grid1,0*grid1,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
        plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
        set(gca,TickLabelInterpreter = 'latex') ;
        title('Uniform Grid','interpreter','latex') ;
        set(gca,'fontsize',12)
        ylim([-5e-4,5e-4])
        subplot(3,1,2)
        %plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
        plot([0;grid2;L],[0;s2Approx;0]-[0;s2;0],'r') ; hold on;
        plot(grid2,0*grid2,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
        plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
        set(gca,TickLabelInterpreter = 'latex') ;
        title('Smooth Grid','interpreter','latex') ;
        set(gca,'fontsize',12)
        ylim([-5e-4,5e-4])
        subplot(3,1,3)
        %plot(gridEx,solEx,'color',[0.5,0.5,0.5],'linewidth',5); hold on;
        plot([0;grid3;L],[0;s3Approx;0]-[0;s3;0],'g') ; hold on;
        plot(grid3,0*grid3,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','k') ;
        plot([0,L],[0,0],'linestyle','none','marker','o','markerfacecolor','c','markeredgecolor','k')
        set(gca,TickLabelInterpreter = 'latex') ;
        title('Random Grid','interpreter','latex') ;
        set(gca,'fontsize',12)

    end

end


figure
loglog(L./Mvec,errvec(1,:),'k') ; hold on ;
loglog(L./Mvec,errvec(2,:),'r') ;
loglog(L./Mvec,errvec(3,:),'g') ;
loglog(L./Mvec,errvec(4,:),'--k') ;
set(gca,TickLabelInterpreter = 'latex') ;
title('RMSD','interpreter','latex') ;
xlabel('$\Delta_x$','interpreter','latex') ;
ylabel('$\sqrt{\frac{1}{M-1}\sum_{i=1}^{M-1}(u_i - \hat u_i)^2}$','interpreter','latex') ;
legend('Uniform', 'Smooth', 'Random','Uniform (Fourth-Order)') ;
set(gca,'fontsize',14)




