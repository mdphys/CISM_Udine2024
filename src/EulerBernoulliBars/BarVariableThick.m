%-------------------------------------------------------------------------
%                            CISM
%          Modes of A Flexural Bar with Variable Cross Section
%                    Dr Michele Ducceschi
%                   University of Bologna
%                        10 May 2024
%-------------------------------------------------------------------------


clear all
close all
clc
for nBar = 1 : 3
    % This scritp allows computing the modal shapes and frequencies of a bar
    % with a variable cross section. The initial flat length is xb;
    % the slope extends over a length pb; the plateau has length pr;
    % hence, the total length of the bar is
    %
    %                      L = 2*xb + 2*pb + pr ;

    % see the diagram below
    %
    %
    %      |   xb    |pb|     pr    |pb|   xb    |
    %
    %       ---------                  ---------
    %       |          \ ----------- /          |
    %       -------------------------------------
    %
    %
    % The largest and smalles thickness values are ymax and ymin, respectively.
    % A rectangular cross section is assumed, with width wb. The area is
    %
    %                           A = wb * y,
    %
    % and the area moment of inertia is
    %
    %                      I = 1/12 * wb * y^3,
    %
    % where y is the local thickness.



    %-------------------------------------------------------------------------
    %-- custom bar parameters

    % cross section
    ymax           = 4e-2 ;           %-- largest thickness [m]
    if nBar==1
        ymin           = 3e-2 ;           %-- smallest thickness [m]
    elseif nBar == 2
        ymin           = 2e-2 ;
    else
        ymin = 0.9e-2 ;
    end
    wb             = 6e-2 ;           %-- width of cross sec [m]

    % thickness profile (see diagram above)
    xb             = 0.1  ;           %-- [m]
    xw             = 0.5 ;            %-- [m]
    p              = 4 ;              %-- exponent in thickness profile

    % material parameters
    E              = 16e9 ;           %-- young's mod [Pa]
    rho            = 1e3 ;            %-- volume density [kg/m^3]


    adaptivemeshOn = 1 ;

    ppw            = 50 ;              %-- points per wavelength
    f_ref          = 10000 ;           %-- fmax
    %------------------------------------------------------------------------



    %-------------------------------------------------------------------------
    %-- derived parameters (don't touch below!)
    Delta_y        = ymax - ymin ;
    xc             = xb + 0.5*xw ;
    L              = 2*xb + xw ;

    % largest and smallest grid spacings
    Amax           = wb*ymax ;
    Imax           = 1/12 * wb * ymax^3 ;
    cphi_max       = sqrt(sqrt(E*Imax/rho/Amax)*(2*pi*f_ref)) ;
    hmax           = 2*pi*cphi_max/ppw/(2*pi*f_ref)  ;

    Amin           = wb*ymin ;
    Imin           = 1/12 * wb * ymin^3 ;
    cphi_min       = sqrt(sqrt(E*Imin/rho/Amin)*(2*pi*f_ref)) ;
    hmin           = 2*pi*cphi_min/ppw/(2*pi*f_ref)  ;
    %------------------------------------------------------------------------



    %-------------------------------------------------------------------------
    %-- mesh init
    l              = 0 ;
    grid           = 0 ;
    hvec           = [] ;
    hprev          = hmax ;
    I              = Imax ;
    A              = Amax ;

    if adaptivemeshOn == 1

        hprev         = hmax  ;
    else
        if nAdaptive == 1
            h = hmax ;
        else
            h = hmin ;
        end
    end

    while l < 2*L

        if adaptivemeshOn == 1
            cphi      = sqrt(sqrt(E*I/rho/A)*(2*pi*f_ref)) ;
            h         = 2*pi*cphi/ppw/(2*pi*f_ref)  ;
        else
            if nAdaptive == 1
                h = hmax ;
            else
                h = hmin ;
            end
        end

        l     = l+h;
        if l > L
            l = l-h ;
            h = L - grid(end);
            l = l+h ;
            grid  = [grid; l] ;
            hvec  = [hvec;h] ;
            break
        else
            grid  = [grid; l] ;
            hvec  = [hvec;h] ;
        end


        if l >= xb && l < xb + xw
            y = (ymin + Delta_y/(0.5*xw)^p*(abs(l-xc)).^p) ;
        else
            y   = ymax ;
        end

        I     = (wb*y^3)/12 ;
        A     = wb*y ;
        hprev = h ;

    end
    %-------------------------------------------------------------------------
    hvec(1:end-1) = 0.5*hvec(1:end-1) + 0.5*hvec(2:end) ;
    hvec(end) = L - sum(hvec(1:end-1));

    M = length(grid)-1;
    l = 0 ;
    y = ymax ;
    grid = 0 ;
    Avec = Amax;
    Ivec = Imax ;
    tvec = ymax ;

    for m = 1 : M

        l     = l+hvec(m);
        grid = [grid; l] ;
        if l >= xb && l < xb + xw
            y = (ymin + Delta_y/(0.5*xw)^p*(abs(l-xc)).^p) ;
        else
            y   = ymax ;
        end

        I     = (wb*y^3)/12 ;
        A     = wb*y ;
        Avec  = [Avec;A] ;
        Ivec  = [Ivec;I] ;
        tvec  = [tvec;y] ;
    end
    %Avec(2:end-1) = 0.5*Avec(3:end) + 0.5*Avec(1:end-2) ;
    %Ivec(2:end-1) = 0.5*Ivec(3:end) + 0.5*Ivec(1:end-2) ;


    %-------------------------------------------------------------------------
    %-- eigenvalue problem
    Ivec      = Ivec(2:end-1) ;



    D2 = zeros(M-1,M+1) ;
    D2T = zeros(M+1,M-1) ;

    for m = 1 : M-1
        D2(m,m)   = 2/hvec(m)/(hvec(m)+hvec(m+1)) ;
        D2(m,m+1) = -2/hvec(m)/hvec(m+1) ;
        D2(m,m+2) = 2/hvec(m+1)/(hvec(m)+hvec(m+1)) ;
    end

    % hvec = [hvec(1);hvec;hvec(end)];
    for m = 1 : M-1
        if m == 1
            D2T(1) = 1/hvec(1)/hvec(1);
        else
            D2T(m,m)   = 2/hvec(m)/(hvec(m-1)+hvec(m))  ;
        end
        D2T(m+1,m) = -2/hvec(m)/hvec(m+1) ;
        if m == M-1
            D2T(end) = 1/hvec(end)/hvec(end);
        else
            D2T(m+2,m) = 2/hvec(m+1)/(hvec(m+1)+hvec(m+2)) ;
        end
    end





    %D2(1)=1*D2(1) ; D2(end)=1*D2(end);
    biHarm             = D2T*diag(Ivec)*D2 ;

    KM                 = (rho * diag(Avec)) \ (E * (biHarm + 1e-9*eye(M+1)) ) ;
    [V,D]              = eig(KM) ;
    freqs              = sqrt(diag(D)) ;
    [freqs,indfreqs]   = sort(freqs) ;
    V                  = V(:,indfreqs) ;
    %  FDfreqs            = 2/k*asin(0.5*k*freqs) ;

    freqs              = freqs / 2 / pi  ;

    %-------------------------------------------------------------------------


    %-------------------------------------------------------------------------
    % %-- plot stuff
    figure(1) ;
    %  plot(grid, zeros*grid,'linestyle','none','marker','square','markerfacecolor','k','markeredgecolor','none') ; hold on ;
    if nBar == 1
        plot(grid,tvec./1e-2,'k'); hold on ;
    elseif nBar == 2
        plot(grid,tvec./1e-2,'b');
    else
        plot(grid,tvec./1e-2,'r');
    end
    set(gca,'TickLabelInterpreter', 'latex','fontsize',13) ;
    xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (cm)','interpreter','latex')


    % figure
    % plot(freqs); hold on; plot(FDfreqs);

    figure(2)
    for n = 3 : 6
        subplot(2,2,n-2)
        if nBar == 1
            modeshape = V(:,n) ;
            if modeshape(1) < 0
                modeshape = -modeshape ;
            end
            plot(grid,modeshape,'k'); hold on; xlim([0,L])
        elseif nBar == 2
            modeshape = V(:,n) ;
            if modeshape(1) < 0
                modeshape = -modeshape ;
            end
            plot(grid,modeshape,'b'); hold on ; xlim([0,L])
        elseif nBar == 3
            modeshape = V(:,n) ;
            if modeshape(1) < 0
                modeshape = -modeshape ;
            end
            plot(grid,modeshape,'r'); hold on ; xlim([0,L])
        end
    end

    [freqs(3)/freqs(3), freqs(4)/freqs(3), freqs(5)/freqs(3)]
    
    figure(3)
    if nBar == 1
        plot(grid, grid./grid*freqs(3)/freqs(3),'k') ; hold on ;
        plot(grid, grid./grid*freqs(4)/freqs(3),'k') ;
        plot(grid, grid./grid*freqs(5)/freqs(3),'k') ;
    elseif nBar == 2
        plot(grid, grid./grid*freqs(3)/freqs(3),'b') ; hold on ;
        plot(grid, grid./grid*freqs(4)/freqs(3),'b') ;
        plot(grid, grid./grid*freqs(5)/freqs(3),'b') ;
    elseif nBar == 3
        plot(grid, grid./grid*freqs(3)/freqs(3),'r') ; hold on ;
        plot(grid, grid./grid*freqs(4)/freqs(3),'r') ;
        plot(grid, grid./grid*freqs(5)/freqs(3),'r') ;
    end

end



