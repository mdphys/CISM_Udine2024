%-------------------------------------------------------------------------
%                                 CISM
%                  APPROXIMATION OF DERIVATIVES (POISSON)
%             CISM Course "Physics of Musical Instruments"
%                          Michele Ducceschi
%                       University of Bologna
%                            11 May 2024
%-------------------------------------------------------------------------
%
% Purpose
% -------
% Solve the 1D Poisson equation
%       s''(x) = w(x),   x ∈ (0, L),   with  s(0) = s(L) = 0
% using second-order finite differences on three grid types:
%   1) Uniform
%   2) Smoothly varying
%   3) Random (sorted)
%
% Exact Solution
% --------------
% The manufactured exact solution is:
%       s(x) = -log(cosh(x)) + (x/L)*log(cosh(L))
% which implies:
%       w(x) = s''(x) = -(1 − tanh(x)^2) = −sech(x)^2
%
% Notes
% -----
% - The script assembles the discrete Laplacian for nonuniform grids.
% - Compares numerical and analytical solutions for validation.
%-------------------------------------------------------------------------

clear; close all; clc;

%% -------------------------- Parameters ---------------------------------
M = 20;             % Number of intervals
L = 1.2;            % Domain length
h = L/M;            % Reference spacing for the uniform grid

% High-resolution grid for plotting the exact solution
gridEx = linspace(0, L, 1000).';
solEx  = s_exact(gridEx, L);

%% ------------------------- Define Grids --------------------------------
% Uniform grid (interior nodes only)
grid1 = (1:M-1).' * h;

% Smooth grid (monotone increasing after flipping)
grid2 = (L^2 - ((1:M-1).' * h).^2) / L;
grid2 = flipud(grid2);   % ensure ascending order

% Random grid (sorted interior nodes)
rng(0);                  % for reproducibility
grid3 = sort(L * rand(M-1, 1));

%% --------------------- Right-hand side w(x) ----------------------------
w1 = w_rhs(grid1);
w2 = w_rhs(grid2);
w3 = w_rhs(grid3);

%% -------------------- Discrete operators (A) ---------------------------
A1 = assemble_laplacian_1d_nonuniform(grid1, L);
A2 = assemble_laplacian_1d_nonuniform(grid2, L);
A3 = assemble_laplacian_1d_nonuniform(grid3, L);

%% ------------------------- Solve systems -------------------------------
s1 = A1 \ w1;
s2 = A2 \ w2;
s3 = A3 \ w3;

%% --------------------------- Errors ------------------------------------
e1 = s_exact(grid1, L) - s1;
e2 = s_exact(grid2, L) - s2;
e3 = s_exact(grid3, L) - s3;

L2_1 = sqrt(trapz([0;grid1;L], [0; e1; 0].^2));
L2_2 = sqrt(trapz([0;grid2;L], [0; e2; 0].^2));
L2_3 = sqrt(trapz([0;grid3;L], [0; e3; 0].^2));

fprintf('L2 errors:  uniform=%.3e,  smooth=%.3e,  random=%.3e\n', L2_1, L2_2, L2_3);

%% --------------------------- Plotting ----------------------------------
figure('Color','w');

subplot(3,1,1);
hold on;
plot(gridEx, solEx, '-', 'LineWidth', 1.5);                 % exact curve
plot([0; grid1; L], [0; s1; 0], '--', 'LineWidth', 1.2);    % numerical
plot(grid1, zeros(size(grid1)), 'ks', 'MarkerFaceColor','k'); % nodes
plot([0 L], [0 0], 'o', 'MarkerFaceColor',[0.7 1 1], 'MarkerEdgeColor','k');
title('Uniform Grid','Interpreter','latex');
xlabel('$x$','Interpreter','latex'); ylabel('$s(x)$','Interpreter','latex');
legend('Exact','Numerical','Nodes','BCs','Interpreter','latex','Location','best');
set(gca,'TickLabelInterpreter','latex','FontSize',12); ylim([-0.05, 0.15]); grid on;

subplot(3,1,2);
hold on;
plot(gridEx, solEx, '-', 'LineWidth', 1.5);
plot([0; grid2; L], [0; s2; 0], '--', 'LineWidth', 1.2);
plot(grid2, zeros(size(grid2)), 'ks', 'MarkerFaceColor','k');
plot([0 L], [0 0], 'o', 'MarkerFaceColor',[0.7 1 1], 'MarkerEdgeColor','k');
title('Smooth Grid','Interpreter','latex');
xlabel('$x$','Interpreter','latex'); ylabel('$s(x)$','Interpreter','latex');
legend('Exact','Numerical','Nodes','BCs','Interpreter','latex','Location','best');
set(gca,'TickLabelInterpreter','latex','FontSize',12); ylim([-0.05, 0.15]); grid on;

subplot(3,1,3);
hold on;
plot(gridEx, solEx, '-', 'LineWidth', 1.5);
plot([0; grid3; L], [0; s3; 0], '--', 'LineWidth', 1.2);
plot(grid3, zeros(size(grid3)), 'ks', 'MarkerFaceColor','k');
plot([0 L], [0 0], 'o', 'MarkerFaceColor',[0.7 1 1], 'MarkerEdgeColor','k');
title('Random Grid','Interpreter','latex');
xlabel('$x$','Interpreter','latex'); ylabel('$s(x)$','Interpreter','latex');
legend('Exact','Numerical','Nodes','BCs','Interpreter','latex','Location','best');
set(gca,'TickLabelInterpreter','latex','FontSize',12); ylim([-0.05, 0.15]); grid on;

sgtitle('1D Poisson: Exact vs. Finite Difference on Nonuniform Grids', ...
    'Interpreter','latex','FontSize',16);

%% ======================= Local helper functions ========================
function A = assemble_laplacian_1d_nonuniform(xi, L)
%ASSEMBLE_LAPLACIAN_1D_NONUNIFORM  Second-derivative matrix on (0,L)
% for a *nonuniform* grid with homogeneous Dirichlet BCs at x=0 and x=L.
% INPUT:
%   xi : column vector of interior nodes (ascending), size (M-1) x 1
%   L  : domain length
% OUTPUT:
%   A  : sparse (M-1) x (M-1) matrix approximating d^2/dx^2
%
% The stencil is the standard second-order 3-point formula on nonuniform
% grids (derived from local quadratic interpolation).

    xi = xi(:);                 % ensure column
    M1 = numel(xi);             % number of interior nodes
    A  = spalloc(M1, M1, 3*M1); % allocate tridiagonal

    % Interior nodes: i = 2..M-2 in 1-based indexing
    for m = 2:M1-1
        hm   = xi(m)   - xi(m-1);         % left spacing
        hp   = xi(m+1) - xi(m);           % right spacing
        A(m,m)   = -2/(hm*hp);
        A(m,m-1) =  2/(hm*(hm+hp));
        A(m,m+1) =  2/(hp*(hm+hp));
    end

    % Left boundary row (couples to first interior neighbor and boundary x=0)
    hm = xi(1) - 0;                       % distance to left boundary
    hp = xi(2) - xi(1);                   % right spacing
    A(1,1) = -2/(hm*hp);
    A(1,2) =  2/(hp*(hm+hp));

    % Right boundary row (couples to last interior neighbor and boundary x=L)
    hm = xi(end) - xi(end-1);             % left spacing
    hp = L - xi(end);                     % distance to right boundary
    A(end,end)   = -2/(hm*hp);
    A(end,end-1) =  2/(hm*(hm+hp));
end

function y = w_rhs(x)
%W_RHS  Right-hand side w(x) = s''(x) for manufactured solution
    y = -(1 - tanh(x).^2);   % = -sech^2(x)
end

function s = s_exact(x, L)
%S_EXACT  Exact solution satisfying s(0)=s(L)=0
    s = -log(cosh(x)) + (x./L) * log(cosh(L));
end
