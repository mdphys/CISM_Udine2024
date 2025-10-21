%-------------------------------------------------------------------------
%                                 CISM
%        APPROXIMATION OF DERIVATIVES – VANDERMONDE INTERPOLATION
%             CISM Course "Physics of Musical Instruments"
%                          Michele Ducceschi
%                       University of Bologna
%                            13 Dec 2023
%-------------------------------------------------------------------------
%
% Purpose
% -------
% Let  u(x) = exp(sin x).  Around an expansion point  x0, we:
%   1) Choose a small set of sample points {x_j}.
%   2) Build the Vandermonde system in the shifted monomial basis
%         (x − x0)^(k−1).
%   3) Compute the polynomial interpolant  p(x)  matching  u(x_j).
%   4) Extract  p'(x0)  from the coefficient of  (x − x0)  (i.e., a₂).
%
% Comparison
% ----------
% Two distinct node distributions are compared and visualized alongside:
%   • the exact function  u(x)
%   • the interpolating polynomial  p(x)
%   • the expansion point  x0  and the chosen sample nodes
%
% Note
% ----
% In this centered monomial basis,
%       p(x) = a₁ + a₂ (x − x0) + a₃ (x − x0)² + ...
% the derivative approximation at  x0  is simply:
%       p'(x0) = a₂
%-------------------------------------------------------------------------

clear all
close all
clc

%% -------------------- Plotting grid and exact function ------------------
hplot = 1e-3;                    % resolution only for plotting
x0    = 1.5;                     % expansion point
x     = -1:hplot:5.0;            % plotting grid
u     = @(z) exp(sin(z));        % function
f     = u(x);                    % exact curve on the plotting grid
ymin  = 0; ymax = 3;             % y-limits for plots

% Exact derivative for reference (not required but useful to print)
uprime_exact = exp(sin(x0)) * cos(x0);

figure('Color','w')

%% ============================ CASE h1 ==================================
% Three nodes, deliberately asymmetric wrt x0
h  = 0.75;
Xs = [0.7*x0, x0 + 0.5*h, x0 + 2*h]';   % sample points
lt = numel(Xs);

% Vandermonde matrix centered at x0: VM(m,n) = (Xs(m) - x0)^(n-1)
VM   = zeros(lt, lt);
for m = 1:lt
    dx = Xs(m) - x0;
    VM(m,1) = 1;
    for n = 2:lt
        VM(m,n) = dx^(n-1);
    end
end

% Solve for coefficients a so that p(x) = sum a_k (x - x0)^(k-1)
uvec  = u(Xs);
acoef = VM \ uvec;

% Evaluate interpolant
dx_all = x - x0;
p = zeros(size(x));
for k = 1:lt
    p = p + acoef(k) * dx_all.^(k-1);
end

% Derivative estimate at x0 from polynomial coefficient
uprime_est_case1 = acoef(2);

% ----- Plot (left) -----
subplot(1,2,1)
plot(x, f, 'k-', 'LineWidth', 1.2); hold on
plot(x, p, 'r-', 'LineWidth', 1.2)
% vertical line at x0
line([x0 x0], [ymin ymax], 'Color','g','LineStyle','-','DisplayName','expansion point')
% vertical lines at sample points (hide extra entries from legend)
for n = 1:lt
    ln = line([Xs(n) Xs(n)], [ymin ymax], 'Color','k','LineStyle','--');
    if n > 1
        set(get(get(ln,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
legend('$u(x)$','$p(x)$','expansion point','Interpreter','latex','Location','best')
title(sprintf('Case 1: %d nodes,  h=%.2f  (p''(x_0)=%.6f,  u''(x_0)=%.6f)', ...
      lt, h, uprime_est_case1, uprime_exact), 'Interpreter','none')
ylim([ymin, ymax]); 

%% ============================ CASE h2 ==================================
% Four nodes, clustered to the right of x0
h  = 0.75;
Xs = [x0 - 0.4*h, x0 + 0.3*h, x0 + 0.8*h, x0 + 1.3*h]';
lt = numel(Xs);

% Vandermonde (centered at x0)
VM   = zeros(lt, lt);
for m = 1:lt
    dx = Xs(m) - x0;
    VM(m,1) = 1;
    for n = 2:lt
        VM(m,n) = dx^(n-1);
    end
end

uvec  = u(Xs);
acoef = VM \ uvec;

% Evaluate interpolant
p = zeros(size(x));
for k = 1:lt
    p = p + acoef(k) * dx_all.^(k-1);
end

uprime_est_case2 = acoef(2);

% ----- Plot (right) -----
subplot(1,2,2)
plot(x, f, 'k-', 'LineWidth', 1.2); hold on
plot(x, p, 'r-', 'LineWidth', 1.2)
line([x0 x0], [ymin ymax], 'Color','g','LineStyle','-','DisplayName','expansion point')
for n = 1:lt
    ln = line([Xs(n) Xs(n)], [ymin ymax], 'Color','k','LineStyle','--');
    if n > 1
        set(get(get(ln,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
legend('$u(x)$','$p(x)$','expansion point','Interpreter','latex','Location','best')
title(sprintf('Case 2: %d nodes,  h=%.2f  (p''(x_0)=%.6f,  u''(x_0)=%.6f)', ...
      lt, h, uprime_est_case2, uprime_exact), 'Interpreter','none')
ylim([ymin, ymax]); 

%% ---------------------------- Summary ----------------------------------
fprintf('Expansion point x0 = %.6f\n', x0);
fprintf('  Exact   u''(x0) = %.10f\n', uprime_exact);
fprintf('  Case 1  p''(x0) = %.10f\n', uprime_est_case1);
fprintf('  Case 2  p''(x0) = %.10f\n', uprime_est_case2);
