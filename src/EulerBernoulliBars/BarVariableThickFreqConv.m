%-------------------------------------------------------------------------
%                                 CISM
%      Modes of a Flexural Bar with Variable Cross Section (EB Beam)
%         CISM Course "Physics of Musical Instruments" — Convergence
%                          Michele Ducceschi
%                       University of Bologna
%                            18 Dec 2024
%-------------------------------------------------------------------------
%
% Purpose
% -------
% Study how grid resolution and adaptivity affect the first few natural
% frequencies of an Euler–Bernoulli bar with *variable thickness*.
% We compare three meshing strategies:
%   (1) Uniform coarse mesh      (uses h_max)
%   (2) Uniform fine   mesh      (uses h_min)
%   (3) Adaptive mesh            (h follows bending wave speed)
%
% Geometry
% --------
% Total length: L = 2*xb + xw
% Thickness y(x):
%   - flat (y = ymax) outside the central varying region of width xw
%   - smoothly varying in the central region with exponent p
% Section:
%   A(x) = wb * y(x),    I(x) = (wb * y(x)^3)/12
%
% Discretization
% --------------
% Build a nonuniform second-derivative operator (D2, D2^T) from element
% sizes h_i, assemble the bilaplacian with variable stiffness:
%       K = D2^T * diag(I_interior) * D2
% Lumped mass matrix: M = diag(rho * A)
% Solve the standard eigenproblem:
%       (rho M)^{-1} (E K + eps I) v = (2π f)^2 v
%
% Outputs
% -------
% - freqsTot: first 15 natural frequencies (Hz) for each test
% - Plots: (a) thickness profile; (b) frequency-vs-DOFs trends for
%          the top three plotted modes (highest 3 among the first 15)
%
%-------------------------------------------------------------------------

clear all
close all
clc

% ------------------------- Global study setup ----------------------------
Ntot = 30;                    % number of mesh refinements per strategy
Nkeep = 15;                   % number of eigenfrequencies to store/plot

% Storage across strategies (will overwrite per strategy loop)
freqsTot = zeros(Nkeep, Ntot);
Mvec     = zeros(Ntot, 3);    % number of elements per test & strategy
Mdiff    = zeros(Ntot-1, 3);  % mid-DOF abscissae for difference trends

% ---------------------------- Loop strategies ---------------------------
for nAdaptive = 1:3
    % Strategy flags:
    %  1: uniform (coarse = hmax)
    %  2: uniform (fine   = hmin)
    %  3: adaptive (based on local c_phi)
    freqsTot(:) = 0;  % reset storage for this strategy

    for nTest = 1:Ntot
        % fprintf('Strategy %d — test %d/%d\n', nAdaptive, nTest, Ntot);

        %==================== Problem parameters ==========================
        % Cross-section
        ymax = 4e-2;                  % max thickness [m]
        ymin = 1e-2;                  % min thickness [m]
        wb   = 6e-2;                  % width [m]

        % Thickness profile
        xb = 0.1;                     % flat end length [m]
        xw = 0.5;                     % width of varying region [m]
        p  = 4;                       % even exponent in profile

        % Material
        E   = 16e9;                   % Young's modulus [Pa]
        rho = 1e3;                    % density [kg/m^3]

        % Meshing controls
        if nAdaptive == 1 || nAdaptive == 2
            adaptivemeshOn = false;
        else
            adaptivemeshOn = true;
        end
        ppw   = 5*nTest;              % points per wavelength @ f_ref
        f_ref = 1.0e4;                % reference frequency [Hz]
        eps_stab = 1e-9;              % small regularization for K

        %==================== Derived quantities ==========================
        Delta_y = ymax - ymin;
        xc = xb + 0.5*xw;
        L  = 2*xb + xw;

        % max/min section metrics
        Amax = wb * ymax;
        Imax = (wb * ymax^3) / 12;
        Amin = wb * ymin;
        Imin = (wb * ymin^3) / 12;

        % bending wave speed proxy c_phi ~ (E I / rho A)^(1/4) * sqrt(2π f_ref)
        cphi_max = sqrt( sqrt(E*Imax/(rho*Amax)) * (2*pi*f_ref) );
        cphi_min = sqrt( sqrt(E*Imin/(rho*Amin)) * (2*pi*f_ref) );

        % target element sizes
        hmax = 2*pi*cphi_max / (ppw*(2*pi*f_ref));
        hmin = 2*pi*cphi_min / (ppw*(2*pi*f_ref));

        %==================== Mesh generation =============================
        x_nodes = 0;
        hvec    = [];

        % initialize section for first step
        Acur = Amax; Icur = Imax;
        x = 0;

        while x < L
            if adaptivemeshOn
                cphi = sqrt( sqrt(E*Icur/(rho*Acur)) * (2*pi*f_ref) );
                h    = 2*pi*cphi / (ppw*(2*pi*f_ref));
                % clamp h within [hmin, hmax] for robustness
                h = max(hmin, min(hmax, h));
            else
                if nAdaptive == 1
                    h = hmax; % coarse uniform
                else
                    h = hmin; % fine   uniform
                end
            end

            x_next = x + h;
            if x_next > L
                h     = L - x;
                x_next= L;
            end

            % append
            x_nodes = [x_nodes; x_next]; %#ok<AGROW>
            hvec    = [hvec; h];         %#ok<AGROW>
            x       = x_next;

            % update section for next step (use current x)
            if (x >= xb) && (x < xb + xw)
                y = ymin + Delta_y/(0.5*xw)^p * abs(x - xc)^p;
            else
                y = ymax;
            end
            Icur = (wb*y^3)/12;
            Acur = wb*y;
        end

        % Smooth h for nicer operators (optional; preserves total length)
        if numel(hvec) > 1
            hvec(1:end-1) = 0.5*hvec(1:end-1) + 0.5*hvec(2:end);
            hvec(end) = L - sum(hvec(1:end-1));
        end

        % Rebuild nodal fields
        M      = numel(x_nodes) - 1;    % #elements; #nodes = M+1
        Mvec(nTest, nAdaptive) = M;

        x    = 0;
        grid = 0;
        Avec = Amax; Ivec = Imax; tvec = ymax;

        for m = 1:M
            x = x + hvec(m);
            grid = [grid; x]; %#ok<AGROW>

            if (x >= xb) && (x < xb + xw)
                y = ymin + Delta_y/(0.5*xw)^p * abs(x - xc)^p;
            else
                y = ymax;
            end
            I = (wb*y^3)/12;  A = wb*y;

            Avec = [Avec; A]; %#ok<AGROW>
            Ivec = [Ivec; I]; %#ok<AGROW>
            tvec = [tvec; y]; %#ok<AGROW>
        end

        %==================== Operators & eigenproblem =====================
        % Use interior I for K = D2^T * diag(I_interior) * D2
        I_interior = Ivec(2:end-1);

        % Nonuniform 2nd derivative: D2 (size (M-1) x (M+1))
        D2 = spalloc(M-1, M+1, 3*(M-1));
        for m = 1:M-1
            hm = hvec(m); hp = hvec(m+1);
            D2(m, m)   =  2/(hm*(hm+hp));
            D2(m, m+1) = -2/(hm*hp);
            D2(m, m+2) =  2/(hp*(hm+hp));
        end

        % Adjoint-like D2T (size (M+1) x (M-1)) with boundary closures
        D2T = spalloc(M+1, M-1, 3*(M-1));
        for m = 1:M-1
            hm = hvec(m); hp = hvec(m+1);

            if m == 1
                D2T(1,1) = 1/(hm*hm);
            else
                D2T(m, m) = 2/(hm*(hvec(m-1)+hm));
            end
            D2T(m+1, m) = -2/(hm*hp);

            if m == M-1
                D2T(M+1, M-1) = 1/(hp*hp);
            else
                D2T(m+2, m) = 2/(hp*(hp + hvec(m+2)));
            end
        end

        % Bilaplacian & mass matrices
        K  = D2T * spdiags(I_interior, 0, M-1, M-1) * D2;
        Mm = spdiags(rho * Avec, 0, M+1, M+1);

        % Generalized -> standard eigenproblem: KM v = (2pi f)^2 v
        KM = Mm \ (E * (K + eps_stab*speye(M+1)));

        % Compute eigenpairs (dense here for clarity; sparse eigs is possible)
        [V, Dlam] = eig(full(KM));
        lam       = diag(Dlam);
        [lam, idx] = sort(lam, 'ascend');
        freqs = sqrt(max(lam, 0)) / (2*pi);   % guard against tiny negatives
        freqs = freqs(idx);

        % Store first Nkeep frequencies
        nk = min(Nkeep, numel(freqs));
        freqsTot(1:nk, nTest) = freqs(1:nk);

        % (Optional) quick visualization of node placement for a sample test
        if nTest == 3
            figure(1)
            subplot(2,1,2)
            plot(grid, 0.1*nAdaptive + ones(numel(grid),1), 'k+'); hold on
            set(gca,'TickLabelInterpreter','latex')
            xlabel('$x$ (m)','Interpreter','latex')
            ylabel('strategy id','Interpreter','latex')
            title('Sample node layouts (test 3)','Interpreter','latex')
            drawnow
        end
    end % nTest

    %-------------------- Post-process per strategy -----------------------
    freqsTot = abs(freqsTot);
    % simple mid-DOF abscissae (if e.g., differences are to be plotted)
    Mdiff(:, nAdaptive) = 0.5*Mvec(2:end, nAdaptive) + 0.5*Mvec(1:end-1, nAdaptive);

    % Thickness profile (from the *last* built grid in this loop)
    figure(1)
    subplot(2,1,1)
    plot(grid, tvec, 'LineWidth',1.2); hold on
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex')
    ylabel('$y$ (m)','Interpreter','latex')
    title('Thickness profile','Interpreter','latex')
end % strategies

%========================== Frequency trends ==============================
figure(2)
for n = 1:3
    subplot(3,1,n); hold on; grid on
    if any(Mvec(:,1)>0)
        plot(Mvec(:,1), freqsTot(end-n+1,:), 'g-o', 'DisplayName','Uniform (coarse)')
    end
    if any(Mvec(:,2)>0)
        plot(Mvec(:,2), freqsTot(end-n+1,:), 'b+--', 'DisplayName','Uniform (fine)')
    end
    if any(Mvec(:,3)>0)
        plot(Mvec(:,3), freqsTot(end-n+1,:), 'r-*', 'DisplayName','Adaptive')
    end
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Elements $M$','Interpreter','latex')
    ylabel(sprintf('$f_{%d}$ (Hz)', Nkeep-n+1),'Interpreter','latex')
    title(sprintf('Trend of high-order mode f_{%d}', Nkeep-n+1), 'Interpreter','latex')
    legend('Location','best','Interpreter','latex')
end
