%------------------------------------------------------------------------
%                   Approximation of Derivatives
%          Task: Solve Poisson's Equation using Grids
%                 CISM Course "Physics of Musical Instruments"
%                         Michele Ducceschi
%                      University of Bologna
%                         11 May 2024
%-------------------------------------------------------------------------
%
% This script studies the effect of grid choice on the numerical solution
% of the 1D Poisson problem with homogeneous Dirichlet boundary conditions:
%
%          s''(x) = w(x),    x in (0, L),      s(0) = s(L) = 0.
%
% We use a manufactured (exact) solution
%     s(x) = -log(cosh(x)) + (x/L)*log(cosh(L)),
% which implies        w(x) = s''(x) = -(1 - tanh(x)^2) = -sech(x)^2.
%
% The PDE is discretized with a second–order three-point finite-difference
% stencil on three grids:
%   1) Uniform grid
%   2) Smoothly varying grid (monotone)
%   3) Random grid (sorted)
%
% For the uniform grid we also compare a 4th–order centered scheme.
% For a range of mesh sizes M, we solve the linear systems and collect
% RMS errors against the exact solution on the interior nodes.
%
% Plots:
%   - Pointwise error profiles (for a representative mesh)
%   - Convergence plot: RMS error versus Δx in log–log scale
%
%-------------------------------------------------------------------------

clear all
close all
clc

% -------------------- Experiment parameters ------------------------------
Ntest  = 20;                            % number of mesh refinements
Mvec   = round(logspace(1,3,Ntest));    % M spans ~[10, 1000]
L      = 1.2;                           % domain length
errvec = zeros(4,Ntest);                % RMS errors: [uni2, smooth2, rand2, uni4]

rng(0);                                 % reproducibility for random grid

% -------------------- Sweep over mesh sizes ------------------------------
for nT = 1:Ntest

    M  = Mvec(nT);                      % number of intervals
    h  = L/M;                           % uniform spacing (for ref/uni grid)

    % ---- Construct interior-node grids (ascending order) ----------------
    grid1 = (1:M-1)' * h;                                        % uniform
    grid2 = (L^2 - ((1:M-1)' * h).^2) / L; grid2 = grid2(end:-1:1); % smooth
    grid3 = sort(L * rand(M-1,1));                               % random

    % ---- Exact solution and RHS at interior nodes -----------------------
    s1 = s_exact(grid1, L);
    s2 = s_exact(grid2, L);
    s3 = s_exact(grid3, L);

    w1 = w_rhs(grid1);
    w2 = w_rhs(grid2);
    w3 = w_rhs(grid3);

    % ---- Discrete Laplacians (Dirichlet BCs built-in) -------------------
    % Second-order on nonuniform grids:
    D1 = assemble_laplacian_1d_nonuniform(grid1, L);
    D2 = assemble_laplacian_1d_nonuniform(grid2, L);
    D3 = assemble_laplacian_1d_nonuniform(grid3, L);

    % Fourth-order (only valid/implemented for *uniform* spacing):
    D1fourth = assemble_laplacian_1d_uniform_4th(M-1, h);

    % ---- Solve linear systems A * s_h = w -------------------------------
    s1Approx       = D1       \ w1;
    s2Approx       = D2       \ w2;
    s3Approx       = D3       \ w3;
    s1fourthApprox = D1fourth \ w1;

    % ---- RMS errors over interior nodes --------------------------------
    % (simple discrete RMS; consistent across grids for comparison)
    errvec(1,nT) = sqrt(((s1Approx       - s1)'*(s1Approx       - s1))/(M-1)); % uniform, 2nd
    errvec(2,nT) = sqrt(((s2Approx       - s2)'*(s2Approx       - s2))/(M-1)); % smooth, 2nd
    errvec(3,nT) = sqrt(((s3Approx       - s3)'*(s3Approx       - s3))/(M-1)); % random, 2nd
    errvec(4,nT) = sqrt(((s1fourthApprox - s1)'*(s1fourthApprox - s1))/(M-1)); % uniform, 4th

    % ---- Plot representative pointwise error profiles -------------------
    % Use a moderately fine mesh early in the sweep for readability
    if nT == 2
        figure('Color','w')
        subplot(3,1,1)
        plot([0;grid1;L],[0;s1Approx;0]-[0;s1;0],'k'); hold on
        plot(grid1,0*grid1,'ks','markerfacecolor','k','markeredgecolor','k')
        plot([0,L],[0,0],'o','markerfacecolor','c','markeredgecolor','k')
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12)
        title('Uniform Grid (pointwise error)','Interpreter','latex')
        ylim([-5e-4,5e-4])

        subplot(3,1,2)
        plot([0;grid2;L],[0;s2Approx;0]-[0;s2;0],'r'); hold on
        plot(grid2,0*grid2,'ks','markerfacecolor','k','markeredgecolor','k')
        plot([0,L],[0,0],'o','markerfacecolor','c','markeredgecolor','k')
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12)
        title('Smooth Grid (pointwise error)','Interpreter','latex')
        ylim([-5e-4,5e-4])

        subplot(3,1,3)
        plot([0;grid3;L],[0;s3Approx;0]-[0;s3;0],'g'); hold on
        plot(grid3,0*grid3,'ks','markerfacecolor','k','markeredgecolor','k')
        plot([0,L],[0,0],'o','markerfacecolor','c','markeredgecolor','k')
        set(gca, 'TickLabelInterpreter','latex', 'FontSize',12)
        title('Random Grid (pointwise error)','Interpreter','latex')
        ylim([-5e-4,5e-4])
    end
end

% -------------------- Convergence plot (RMS error) -----------------------
figure('Color','w')
loglog(L./Mvec, errvec(1,:), 'k',  'DisplayName','Uniform (2nd)'); hold on
loglog(L./Mvec, errvec(2,:), 'r',  'DisplayName','Smooth (2nd)')
loglog(L./Mvec, errvec(3,:), 'g',  'DisplayName','Random (2nd)')
loglog(L./Mvec, errvec(4,:), '--k','DisplayName','Uniform (4th)')
set(gca, 'TickLabelInterpreter','latex', 'FontSize',14)
title('RMSD','Interpreter','latex')
xlabel('$\Delta_x$','Interpreter','latex')
ylabel('$\sqrt{\frac{1}{M-1}\sum_{i=1}^{M-1}(u_i - \hat u_i)^2}$','Interpreter','latex')
legend('Location','south east')

% ========================== Local functions ==============================
function A = assemble_laplacian_1d_nonuniform(xi, L)
%ASSEMBLE_LAPLACIAN_1D_NONUNIFORM
% Build the (M-1)x(M-1) second-derivative matrix on (0,L) for *nonuniform*
% interior nodes xi (ascending), enforcing s(0)=s(L)=0 implicitly.
% Second-order accuracy via three-point, nonuniform stencil.

    xi = xi(:);
    M1 = numel(xi);
    A  = spalloc(M1, M1, 3*M1);   % sparse tridiagonal

    % Interior rows
    for m = 2:M1-1
        hm         = xi(m)   - xi(m-1);   % left spacing
        hp         = xi(m+1) - xi(m);     % right spacing
        A(m,m)     = -2/(hm*hp);
        A(m,m-1)   =  2/(hm*(hm+hp));
        A(m,m+1)   =  2/(hp*(hm+hp));
    end

    % Left boundary row (couples to x=0 and first neighbor at xi(2))
    hm = xi(1) - 0;
    hp = xi(2) - xi(1);
    A(1,1) = -2/(hm*hp);
    A(1,2) =  2/(hp*(hm+hp));

    % Right boundary row (couples to x=L and left neighbor at xi(end-1))
    hm = xi(end) - xi(end-1);
    hp = L - xi(end);
    A(end,end)   = -2/(hm*hp);
    A(end,end-1) =  2/(hm*(hm+hp));
end

function A = assemble_laplacian_1d_uniform_4th(Nint, h)
%ASSEMBLE_LAPLACIAN_1D_UNIFORM_4TH
% Fourth-order centered finite-difference Laplacian on a *uniform* grid
% with Dirichlet BCs on (0,L). Returns an Nint x Nint sparse matrix.
%
% Stencil (interior):  (-1/12, 4/3, -2.5, 4/3, -1/12)/h^2
% Ends are adjusted for Dirichlet BCs with one-sided 2nd-order closure.

    e = ones(Nint,1);
    A = spdiags([-1/12*e(1:end-2), 4/3*e(1:end-1), -2.5*e, 4/3*e(1:end-1), -1/12*e(1:end-2)], ...
                [-2, -1, 0, 1, 2], Nint, Nint);

    % Simple 2nd-order boundary closures (as in the original code)
    if Nint >= 2
        A(1, :) = 0;
        A(1,1) = -2; A(1,2) = 1;              % scaled later by 1/h^2
        A(end, :) = 0;
        A(end,end) = -2; A(end,end-1) = 1;    % scaled later by 1/h^2
    end

    A = A / h^2;
end

function y = w_rhs(x)
%W_RHS  Right-hand side w(x) = s''(x) for manufactured solution.
    y = -(1 - tanh(x).^2);   % = -sech^2(x)
end

function s = s_exact(x, L)
%S_EXACT  Exact solution satisfying s(0)=s(L)=0.
    s = -log(cosh(x)) + (x./L) * log(cosh(L));
end
