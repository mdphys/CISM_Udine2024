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


%--- h1
h = 0.75 ;
samplepoints = [0.7*x0,x0+0.5*h,x0+2*h]' ;
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

subplot(1,2,1)
plot(x,f,'LineWidth',1.2,'Color',[0 0 0]); hold on;
plot(x,p,'LineStyle','-',...
    'Color','r','linewidth',1.2) ;
pl = line([x0,x0],[max(f),min(f)],'color','g','linestyle','-') ;

for n = 1 : lt
    fn = samplepoints(n) ;
    pl = line([fn,fn],[max(f),min(f)],'color','k','linestyle','--') ;
    if n > 1
        set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
legend('$u(x)$','$p(x)$','expansion point','Interpreter','latex')
ylim([0,3])


%--- h2
h = 0.75 ;
samplepoints = [x0-0.4*h,x0+0.3*h,x0+0.8*h,x0+1.3*h]' ;
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

subplot(1,2,2)
plot(x,f,'LineWidth',1.2,'Color',[0 0 0]); hold on;
plot(x,p,'LineStyle','-',...
    'Color','r','linewidth',1.2) ;
pl = line([x0,x0],[max(f),min(f)],'color','g','linestyle','-') ;

for n = 1 : lt
    fn = samplepoints(n) ;
    pl = line([fn,fn],[max(f),min(f)],'color','k','linestyle','--') ;
    if n > 1
        set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
legend('$u(x)$','$p(x)$','expansion point','Interpreter','latex')
ylim([0,3])



