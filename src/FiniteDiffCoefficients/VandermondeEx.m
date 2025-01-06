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

h           = 0.001 ;
x0          = 1.5 ;
x           = 0.0:h:3.0 ;
f           = exp(sin(x)) ;
fp          = exp(sin(x0)) + (x-x0)*exp(sin(x0))*cos(x0) ;


%--- h1 
h = 0.2 ;
a = exp(sin(x0)) ;   
b = exp(sin(x0 + h)) ;   
c = exp(sin(x0 + 2*h)) ;  

al0 = a ;
al1 = (2*b)/h - (3*a)/(2*h) - c/(2*h) ;
al2 = a/(2*h^2) - b/h^2 + c/(2*h^2) ;

subplot(1,2,2)
plot(x,f,x,fp,'LineWidth',1.2,'Color',[0 0 0]); hold on; 
p  = al0 + al1*(x-x0) + al2*(x-x0).^2 ;
pp =  al0 + al1*(x-x0) ;
plot(x,p,'LineStyle','--',...
    'Color',[0 0 0]) ; 
plot(x,pp,'LineStyle','--',...
    'Color',[0 0 0]) ; 
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
legend('$u(x)$','$u(x_0)+u^\prime(x_0)(x-x_0)$','$p(x),\, h=0.2$', '$p(x_0) + p^\prime(x_0)(x-x_0),\, h=0.2$','Interpreter','latex')
ylim([1,3])

%--- h2
h = 0.4 ;
y = 0:h:3 ;
a = exp(sin(x0)) ;   
b = exp(sin(x0 + h)) ;   
c = exp(sin(x0 + 2*h)) ;  

al0 = a ;
al1 = (2*b)/h - (3*a)/(2*h) - c/(2*h) ;
al2 = a/(2*h^2) - b/h^2 + c/(2*h^2) ;

subplot(1,2,1)
plot(x,f,x,fp,'LineWidth',1.2,'Color',[0 0 0]); hold on; 
p = al0 + al1*(x-x0) + al2*(x-x0).^2 ;
pp =  al0 + al1*(x-x0) ;
plot(x,p,'LineStyle',':',...
    'Color',[0 0 0]) ; 
plot(x,pp,'LineStyle',':',...
    'Color',[0 0 0]) ; 

set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
xlabel('$x$','Interpreter','latex') ;
ylabel('$u(x)$','Interpreter','latex') ;
ylim([1,3])


legend('$u(x)$','$u(x_0)+u^\prime(x_0)(x-x_0)$','$p(x),\, h=0.2$', '$p(x_0) + p^\prime(x_0)(x-x_0),\, h=0.2$','Interpreter','latex')
ylim([1,3])