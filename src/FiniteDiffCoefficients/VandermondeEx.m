%-------------------------------------------------------------------------
%                                 CISM
%         APPROXIMATION OF DERIVATIVES – QUADRATIC INTERPOLATION
%             CISM Course "Physics of Musical Instruments"
%                          Michele Ducceschi
%                       University of Bologna
%                            13 Dec 2023
%-------------------------------------------------------------------------
%
% Purpose
% -------
% Let  u(x) = exp(sin x).  Approximate  u'(x0)  using a quadratic
% interpolant constructed through three samples:
%       a = u(x0),   b = u(x0 + h),   c = u(x0 + 2h).
%
% Method
% ------
% Build  p(x) = al0 + al1 (x − x0) + al2 (x − x0)^2
% such that:
%       p(x0) = a,   p(x0 + h) = b,   p(x0 + 2h) = c.
%
% The coefficients are obtained in closed form:
%       al0 = a
%       al1 = (2b)/h − (3a)/(2h) − c/(2h)
%       al2 =  a/(2h²) − b/h² + c/(2h²)
%
% Hence the derivative estimate at x0 is:
%       p'(x0) = al1
%
% Visualization
% -------------
% Compare the quadratic interpolant p(x) with:
%   • the true function u(x)
%   • the first-order Taylor linearization
%         u(x0) + u'(x0)(x − x0)
%
% Results are shown for two step sizes:  h = 0.4  and  h = 0.2.
%-------------------------------------------------------------------------

clear all
close all
clc

%% ---------------------- Setup and exact data ---------------------------
x0   = 1.5;                             % point where we approximate u'(x0)
u    = @(x) exp(sin(x));                % function
up   = @(x) exp(sin(x)).*cos(x);        % exact derivative
u0   = u(x0);                           % u(x0)
up0  = up(x0);                          % u'(x0)

% Grid for visualization
hx_plot = 1e-3;                          % plotting resolution
x       = 0.0:hx_plot:3.0;

f  = u(x);                               % exact function values
fp = u0 + up0*(x - x0);                  % first-order Taylor at x0

%% ---------------------- Case h = 0.4 (left plot) ----------------------
h = 0.4;
a = u0;
b = u(x0 + h);
c = u(x0 + 2*h);

[al0, al1, al2] = quad_coeffs_from_samples(a, b, c, h);
p_h04  = al0 + al1*(x - x0) + al2*(x - x0).^2;   % quadratic interpolant
pp_h04 = al0 + al1*(x - x0);                     % its linear part

% Derivative estimate at x0 and error
uprime_est_h04 = al1;
err_h04        = abs(uprime_est_h04 - up0);

%% ---------------------- Case h = 0.2 (right plot) ---------------------
h = 0.2;
a = u0;
b = u(x0 + h);
c = u(x0 + 2*h);

[al0, al1, al2] = quad_coeffs_from_samples(a, b, c, h);
p_h02  = al0 + al1*(x - x0) + al2*(x - x0).^2;
pp_h02 = al0 + al1*(x - x0);

uprime_est_h02 = al1;
err_h02        = abs(uprime_est_h02 - up0);

%% ----------------------------- Plots -----------------------------------
figure('Color','w')

% --- Left: h = 0.4
subplot(1,2,1)
plot(x, f,  'k-',  'LineWidth',1.2); hold on
plot(x, fp, 'k-',  'LineWidth',1.2, 'Color',[0.3 0.3 0.3])
plot(x, p_h04,  'k:',  'LineWidth',1.2)
plot(x, pp_h04, 'k:',  'LineWidth',1.2, 'Color',[0.4 0.4 0.4])
plot(x0, u0, 'ks', 'MarkerFaceColor','g', 'MarkerSize',8)      % anchor point
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
title(sprintf('Quadratic interpolation, $h=0.4$  (|$p''(x_0)-u''(x_0)$|=%.2e)', err_h04), ...
      'Interpreter','latex')
legend( ...
  '$u(x)$', ...
  '$u(x_0)+u''(x_0)(x-x_0)$', ...
  '$p(x),\,h=0.4$', ...
  '$p(x_0)+p''(x_0)(x-x_0),\,h=0.4$', ...
  'Interpreter','latex', 'Location','best')
ylim([1,3]); 

% --- Right: h = 0.2
subplot(1,2,2)
plot(x, f,  'k-',  'LineWidth',1.2); hold on
plot(x, fp, 'k-',  'LineWidth',1.2, 'Color',[0.3 0.3 0.3])
plot(x, p_h02,  'k--', 'LineWidth',1.2)
plot(x, pp_h02, 'k--', 'LineWidth',1.2, 'Color',[0.4 0.4 0.4])
plot(x0, u0, 'ks', 'MarkerFaceColor','g', 'MarkerSize',8)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex')
ylabel('$u(x)$','Interpreter','latex')
title(sprintf('Quadratic interpolation, $h=0.2$  (|$p''(x_0)-u''(x_0)$|=%.2e)', err_h02), ...
      'Interpreter','latex')
legend( ...
  '$u(x)$', ...
  '$u(x_0)+u''(x_0)(x-x_0)$', ...
  '$p(x),\,h=0.2$', ...
  '$p(x_0)+p''(x_0)(x-x_0),\,h=0.2$', ...
  'Interpreter','latex', 'Location','best')
ylim([1,3]); 

%% --------------------------- Print summary ------------------------------
fprintf('Derivative at x0 = %.6f\n', x0);
fprintf('  Exact u''(x0):         %.10f\n', up0);
fprintf('  Estimate (h=0.4):      %.10f   |error| = %.3e\n', uprime_est_h04, err_h04);
fprintf('  Estimate (h=0.2):      %.10f   |error| = %.3e\n', uprime_est_h02, err_h02);

%% ======================= Local helper function ==========================
function [al0, al1, al2] = quad_coeffs_from_samples(a, b, c, h)
%QUAD_COEFFS_FROM_SAMPLES  Coefficients of quadratic interpolant p(x)
% through the samples a = u(x0), b = u(x0+h), c = u(x0+2h),
% expressed in the monomial basis p(x) = al0 + al1 (x-x0) + al2 (x-x0)^2.
    al0 = a;
    al1 = (2*b)/h - (3*a)/(2*h) - c/(2*h);
    al2 =  a/(2*h^2) - b/h^2    + c/(2*h^2);
end
