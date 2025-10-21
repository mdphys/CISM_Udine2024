%------------------------------------------------------------------------
%                   Approximation of Derivatives
%          Task: approximate the derivative of exp(sin(x))
%         CISM Course "Physics of Musical Instruments"
%                         Michele Ducceschi
%                      University of Bologna
%                         13 Dec 2023
%-------------------------------------------------------------------------
%
% Goal
% ----
% Given u(x) = exp(sin x), approximate u'(x0) at a point x0 using
% polynomial interpolation built from a small set of sample points.
% We compare:
%   • Vandermonde (monomial) interpolation centered at x0
%   • Lagrange interpolation (equivalent polynomial in a different basis)
%
% Key fact
% --------
% If p(x) = a0 + a1 (x - x0) + a2 (x - x0)^2 + ... interpolates u(x),
% then the derivative estimate at x0 is simply p'(x0) = a1.
%
% We demonstrate two choices of sample points and also visualize the
% interpolants p(x) alongside the true function u(x).
%
%-------------------------------------------------------------------------

clear all
close all
clc

% ------------------------ Global parameters -----------------------------
h   = 1e-3;                    % resolution for plotting
x0  = 1.5;                     % point where the derivative is approximated
x   = -1:h:5.0;                % plotting grid
u   = exp(sin(x));             % exact function values on plotting grid

% Exact derivative at x0 for reference
uprime_exact = exp(sin(x0)) * cos(x0);

% ============================ CASE 1 ====================================
% Three sample points
samplepoints = [1.05, 1.87, 3.00]';
lt           = length(samplepoints);

% Function values at sample points
uvec = exp(sin(samplepoints));

% -------- Vandermonde system around x0 (monomial basis at x0) -----------
% VM(m,n) = (samplepoints(m) - x0)^(n-1), n=1..lt
VM = zeros(lt, lt);
for m = 1:lt
    dx = samplepoints(m) - x0;
    VM(m,1) = 1;
    for n = 2:lt
        VM(m,n) = dx^(n-1);
    end
end

% Solve for coefficients a_k so that p(x) = sum a_k (x - x0)^(k-1)
acoef = VM \ uvec;

% Interpolating polynomial (Vandermonde/monomials)
p = zeros(size(x));
dx_all = x - x0;
for k = 1:lt
    p = p + acoef(k) * dx_all.^(k-1);
end

% Lagrange interpolant q(x) built directly from the nodes
q = zeros(size(x));
for m = 1:lt
    lmNum = ones(size(x));
    lmDen = 1;
    for n = 1:lt
        if m ~= n
            lmNum = lmNum .* (x - samplepoints(n));
            lmDen = lmDen * (samplepoints(m) - samplepoints(n));
        end
    end
    q = q + (lmNum / lmDen) * uvec(m);
end

% Derivative estimates at x0
% (the Lagrange and Vandermonde interpolants are the same polynomial, so
% their derivative at x0 is the same; we use the monomial coefficient a1)
uprime_est_case1 = acoef(2);   % since p'(x0) = a1

% ============================ CASE 2 ====================================
% Four sample points
samplepoints2 = [1.20, 1.72, 2.10, 2.47]';
lt2           = length(samplepoints2);

uvec2 = exp(sin(samplepoints2));

% Vandermonde around x0
VM2 = zeros(lt2, lt2);
for m = 1:lt2
    dx = samplepoints2(m) - x0;
    VM2(m,1) = 1;
    for n = 2:lt2
        VM2(m,n) = dx^(n-1);
    end
end

acoef2 = VM2 \ uvec2;

% Interpolating polynomial (case 2)
p2 = zeros(size(x));
for k = 1:lt2
    p2 = p2 + acoef2(k) * dx_all.^(k-1);
end

% Lagrange interpolant (case 2)
q2 = zeros(size(x));
for m = 1:lt2
    lmNum = ones(size(x));
    lmDen = 1;
    for n = 1:lt2
        if m ~= n
            lmNum = lmNum .* (x - samplepoints2(n));
            lmDen = lmDen * (samplepoints2(m) - samplepoints2(n));
        end
    end
    q2 = q2 + (lmNum / lmDen) * uvec2(m);
end

% Derivative estimate at x0 (case 2)
uprime_est_case2 = acoef2(2);

% ------------------------------ Plots -----------------------------------
figure('Color','w')

% Case 1
subplot(1,2,1)
plot(x, u, 'LineWidth', 1.2, 'Color', [0 0 0]); hold on
plot(x, p,  '-',  'Color','r','LineWidth',1.2)
plot(x, q,  '--', 'Color','g','LineWidth',1.2)
plot(samplepoints, uvec, 'ks', 'MarkerFaceColor','k', 'MarkerSize',8)
plot(x0, exp(sin(x0)), 'ks', 'MarkerFaceColor','g', 'MarkerSize',8)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
legend('$u(x)$','$p(x)$ Vandermonde','$p(x)$ Lagrange','Interpreter','latex','Location','best')
title(sprintf('Case 1: 3 nodes (u''(x_0) exact=%.6f, est=%.6f)', uprime_exact, uprime_est_case1), ...
      'Interpreter','none')
ylim([0,3])

% Case 2
subplot(1,2,2)
plot(x, u,  'LineWidth',1.2,'Color',[0 0 0]); hold on
plot(x, p2, '-',  'Color','r','LineWidth',1.2)
plot(x, q2, '--', 'Color','g','LineWidth',1.2)
plot(samplepoints2, uvec2, 'ks', 'MarkerFaceColor','k', 'MarkerSize',8)
plot(x0, exp(sin(x0)), 'ks', 'MarkerFaceColor','g', 'MarkerSize',8)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
legend('$u(x)$','$p(x)$ Vandermonde','$p(x)$ Lagrange','Interpreter','latex','Location','best')
title(sprintf('Case 2: 4 nodes (u''(x_0) exact=%.6f, est=%.6f)', uprime_exact, uprime_est_case2), ...
      'Interpreter','none')
ylim([0,3])

% --------------------------- Print summary ------------------------------
fprintf('Derivative at x0 = %.6f\n', x0);
fprintf('  Exact:                 %.10f\n', uprime_exact);
fprintf('  Case 1 (3 nodes):      %.10f\n', uprime_est_case1);
fprintf('  Case 2 (4 nodes):      %.10f\n', uprime_est_case2);
