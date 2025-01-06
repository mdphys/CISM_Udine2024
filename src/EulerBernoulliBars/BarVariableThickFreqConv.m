%-------------------------------------------------------------------------
%                            CISM
%          Modes of A Flexural Bar with Variable Cross Section
%                    Dr Michele Ducceschi
%                   University of Bologna
%                        18 Dec 2024
%-------------------------------------------------------------------------


clear all
close all
clc

Ntot     = 30 ;

freqsTot = zeros(15,Ntot) ;
Mvec = zeros(Ntot,3) ;
Mdiff = zeros(Ntot-1,3) ;

for nAdaptive = 1 : 3

    for nTest = 1 : Ntot

        nTest

        % This scritp allows computing the modal shapes and frequencies of a bar
        % with a variable cross section. The total length of the bar is
        %
        %                      L = 2*xb + l ;

        % see the diagram below
        %
        %
        %      |   xb    -------- xw -------   xb    |
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
        % The thickness varies as y = y_min + \Delta y / xw^p abs(x)^p,
        %                         with -xw/2 \leq x \leq xw/2



        %-------------------------------------------------------------------------
        %-- custom bar parameters

        % cross section
        ymax           = 4e-2 ;           %-- largest thickness [m]
        ymin           = 1e-2 ;           %-- smallest thickness [m]
        wb             = 6e-2 ;           %-- width of cross sec [m]

        % thickness profile (see diagram above)
        xb             = 0.1  ;           %-- [m]
        xw             = 0.5 ;            %-- [m]
        p              = 4 ;              %-- exponent in thickness profile

        % material parameters
        E              = 16e9 ;           %-- young's mod [Pa]
        rho            = 1e3 ;            %-- volume density [kg/m^3]

        % mesh select
        if nAdaptive == 1 || nAdaptive == 2
            adaptivemeshOn = 0 ;              %-- 1 == adaptive mesh
        else
            adaptivemeshOn = 1 ;
        end
        ppw            = 5*nTest ;              %-- points per wavelength
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
        A              = Amax ;
        I              = Imax ;
        hprev          = hmax ;
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
        Mvec(nTest,nAdaptive) = M ;
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

        figure(1)
        if nTest == 3
            subplot(2,1,2)
            plot(grid,0.1*nAdaptive + ones(length(grid),1),'marker','+') ; hold on ;
            drawnow
        end

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
        freqsTot(1:15,nTest) = freqs(1:15) ;

        %-------------------------------------------------------------------------



    end


    freqsTot = abs(freqsTot) ;
    freqsDiff = freqsTot(:,2:end)-freqsTot(:,1:end-1) ;
    Mdiff(:,nAdaptive) = 0.5*Mvec(2:end,nAdaptive) + 0.5*Mvec(1:end-1,nAdaptive) ;


    figure(1)
    subplot(2,1,1)
    plot(grid, tvec) ;


    figure(2)


    for n = 1 : 3

        if nAdaptive == 1
            subplot(3,1,n)
            plot(Mvec(:,1),(freqsTot(end-n,:)),'g','marker','o'); hold on ;
        elseif nAdaptive == 2
            subplot(3,1,n)
            plot(Mvec(:,2),(freqsTot(end-n,:)),'b','marker','+'); hold on ;
        else
            subplot(3,1,n)
            plot(Mvec(:,3),(freqsTot(end-n,:)),'r','marker','*'); hold on ;

        end
    end

end





