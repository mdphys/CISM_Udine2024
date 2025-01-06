%------------------------------------------------------------------------
%                   Approximation of Derivatives
%          Task: approximate the derivative of exp(sin(x))
%                       Dr M Ducceschi
%                   University of Bologna
%                       13 Dec 2023
%-------------------------------------------------------------------------


clear all
close all
clc

h             = 0.001 ;
x0            = 1.5 ;
x             = -1:h:5.0 ;
f             = exp(sin(x)) ;


%--- case 1: three sample points
samplepoints = [1.05, 1.87, 3.00]' ;
lt           = length(samplepoints) ;

uvec = exp(sin(samplepoints)) ;
VM   = zeros(lt,lt) ;
for m = 1 : lt
    for n = 1 : lt
        VM(m,n) = (samplepoints(m)-x0)^(n-1) ;
    end
end

%-- build polynomial through Vandermonde matrix
alvec = VM \ uvec ;
p = 0 ;
for m = 1 : lt
    p = p + alvec(m)*(x-x0).^(m-1) ;
end

%-- build Lagrange Polynomial
q = 0 ;
for m = 1 : lt
    lmNum = 1 ;
    lmDen = 1 ;
    for n = 1 : lt
        if m ~= n
            lmNum = lmNum.*(x-samplepoints(n)) ; lmDen = lmDen*(samplepoints(m)-samplepoints(n)) ;
        end
    end
    q = q + lmNum/lmDen*uvec(m) ;
end


subplot(1,2,1)
plot(x,f,'LineWidth',1.2,'Color',[0 0 0]); hold on;
plot(x,p,'LineStyle','-',...
    'Color','r','linewidth',1.2) ;
plot(x,q,'LineStyle','--',...
    'Color','g','linewidth',1.2) ;
plot(samplepoints,uvec,'linestyle','none','marker','square','markerfacecolor','k','markersize',8)
plot(x0,exp(sin(x0)),'linestyle','none','marker','square','markerfacecolor','g','markersize',8)

set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
legend('$u(x)$','$p(x)$ Vander','$p(x)$ Lagrange','Interpreter','latex')
ylim([0,3])


%--- case 1: four sample points
samplepoints = [1.20, 1.72, 2.10, 2.47]' ;
lt           = length(samplepoints) ;

uvec = exp(sin(samplepoints)) ;
VM   = zeros(lt,lt) ;
for m = 1 : lt
    for n = 1 : lt
        VM(m,n) = (samplepoints(m)-x0)^(n-1) ;
    end
end

alvec = VM \ uvec ;
p = 0 ;
for m = 1 : lt
    p = p + alvec(m)*(x-x0).^(m-1) ;
end

%-- build Lagrange Polynomial
q = 0 ;
for m = 1 : lt
    lmNum = 1 ;
    lmDen = 1 ;
    for n = 1 : lt
        if m ~= n
            lmNum = lmNum.*(x-samplepoints(n)) ; lmDen = lmDen*(samplepoints(m)-samplepoints(n)) ;
        end
    end
    q = q + lmNum/lmDen*uvec(m) ;
end

subplot(1,2,2)
plot(x,f,'LineWidth',1.2,'Color',[0 0 0]); hold on;
plot(x,p,'LineStyle','-',...
    'Color','r','linewidth',1.2) ;
plot(x,q,'LineStyle','--',...
    'Color','g','linewidth',1.2) ;
plot(samplepoints,uvec,'linestyle','none','marker','square','markerfacecolor','k','markersize',8)
plot(x0,exp(sin(x0)),'linestyle','none','marker','square','markerfacecolor','g','markersize',8)

set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
legend('$u(x)$','$p(x)$ Vander','$p(x)$ Lagrange','Interpreter','latex')
ylim([0,3])



