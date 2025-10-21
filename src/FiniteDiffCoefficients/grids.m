%-------------------------------------------------------------------------
%                                 CISM
%                        GRID GENERATION EXAMPLES
%             CISM Course "Physics of Musical Instruments"
%                          Michele Ducceschi
%                       University of Bologna
%                            11 May 2024
%-------------------------------------------------------------------------
%
% Purpose
% -------
% Illustrate different types of one-dimensional grids used in numerical
% discretizations (finite difference or finite element methods).
%
% Content
% -------
% The script generates and plots three grid types:
%   1) Uniform grid   – evenly spaced nodes
%   2) Smooth grid    – smoothly varying spacing
%   3) Random grid    – irregular node distribution
%
% Notes
% -----
% - These grids are often used to study numerical accuracy and stability
%   in finite-difference schemes.
% - The smooth grid uses a quadratic mapping for node clustering.
%-------------------------------------------------------------------------


clear all
close all
clc

%----------------------------------------------------------------------
% INPUT PARAMETERS
%----------------------------------------------------------------------
M  = 20 ;      % Number of intervals
L  = 1.2 ;     % Domain length

h = L / M ;    % Uniform grid spacing

%----------------------------------------------------------------------
% GRID DEFINITIONS
%----------------------------------------------------------------------
grid1 = (0:M) * h ;                              % Uniform grid
grid2 = (L^2 - ((0:M)*h).^2) / L ;               % Smoothly varying grid
grid3 = [0; L*rand(M-1,1); L] ;                  % Random non-uniform grid

%----------------------------------------------------------------------
% PLOTTING
%----------------------------------------------------------------------
figure('Color','w')

subplot(3,1,1)
plot(grid1,0*grid1,'ks','MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'TickLabelInterpreter','latex','YTick',[],'FontSize',14)
title('Uniform Grid','Interpreter','latex')

subplot(3,1,2)
plot(grid2,0*grid2,'ks','MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'TickLabelInterpreter','latex','YTick',[],'FontSize',14)
title('Smooth Grid','Interpreter','latex')

subplot(3,1,3)
plot(grid3,0*grid3,'ks','MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'TickLabelInterpreter','latex','YTick',[],'FontSize',14)
title('Random Grid','Interpreter','latex')

sgtitle('1D Grid Examples','Interpreter','latex','FontSize',16)
