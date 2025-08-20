function[ustar] = Solve_ustar(umac,umac_temp,vmac_temp,pstar)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    % ucell = zeros(nx+2,ny+2);
    % ucell(2:nx+1,:) = 0.5 * (umac_temp(1:nx,:) + umac_temp(2:nx+1,:));
    % vcell = zeros(nx+2,ny+2);
    % vcell(:,2:ny+1) = 0.5 * (vmac_temp(:,1:ny) + vmac_temp(:,2:ny+1));
    N = (nx-1) * ny;
    % Create an allocate space for sparse matrix
    A = spalloc(N,N,5*N);
    % Diffusion coefficients.
    De_u = visc * dy / dx;
    Dw_u = visc * dy / dx;
    Dn_u = visc * dx / dy;
    Ds_u = visc * dx / dy;
    % Mass flow rates.
    me_u = rho * dy * 0.5*(umac_temp(3:nx+1, 2:ny+1) + umac_temp(2:nx, 2:ny+1));
    mw_u = rho * dy * 0.5*(umac_temp(1:nx-1, 2:ny+1) + umac_temp(2:nx, 2:ny+1));
    mn_u = rho * dx * 0.5*(vmac_temp(2:nx, 2:ny+1)   + vmac_temp(3:nx+1, 2:ny+1));
    ms_u = rho * dx * 0.5*(vmac_temp(2:nx, 1:ny)     + vmac_temp(3:nx+1, 1:ny));
    % The velocity matrices from AP_u*up = AE_u*ue + AW_u*uw + AN_u*un + AS_u*us
    AE_u = zeros(nx+1,ny+2);
    AW_u = zeros(nx+1,ny+2);
    AN_u = zeros(nx+1,ny+2);
    AS_u = zeros(nx+1,ny+2);
    AWW_u = zeros(nx+1,ny+2);
    ASS_u = zeros(nx+1,ny+2);
    %% Change calculation for each scheme.
    switch string(sch)
        case 'CD'
            % Velocity matrices using CD (central differencing) scheme.
            AE_u(2:nx, 2:ny+1) = De_u - me_u/2;
            AW_u(2:nx, 2:ny+1) = Dw_u + mw_u/2;
            AN_u(2:nx, 2:ny+1) = Dn_u - mn_u/2;
            AS_u(2:nx, 2:ny+1) = Ds_u + ms_u/2;
            AP_u = AE_u + AW_u + AN_u + AS_u + rho*dx*dy/dt;
            % Save 'AP_u' as 'AUp' in 'GlobalsSIMPLE.m'.
            AUp(2:nx,2:ny+1) = AP_u(2:nx,2:ny+1);

            %% Get the AP_u and rhs in the AP_u*ustar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs = -dy*diff(pstar) + rho*dx*dy/dt * umac_temp;
            % (2) E
            % rhs = pressure term + AE_u*uE
            rhs(nx, 2:ny+1) = rhs(nx,2:ny+1) + AE_u(nx,2:ny+1).*umac(nx+1,2:ny+1);
            AE_u(nx, 2:ny+1) = 0;
            % (3) W
            % rhs = pressure term + AE_u*uE + AW_u*uW
            rhs(2, 2:ny+1) = rhs(2,2:ny+1) + AW_u(2,2:ny+1).*umac(1,2:ny+1);
            AW_u(2, 2:ny+1) = 0;
            % (4) N
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN
            rhs(2:nx, ny+1) = rhs(2:nx,ny+1) + AN_u(2:nx,ny+1).*umac(2:nx,ny+2);
            AN_u(2:nx, ny+1) = 0;
            % (5) S
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS
            rhs(2:nx, 2) = rhs(2:nx,2) + AS_u(2:nx,2).*umac(2:nx,1);
            AS_u(2:nx, 2) = 0;
            % (6) Relaxation
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS 
            %       + relaxation
            rhs(2:nx, 2:ny+1) = rhs(2:nx,2:ny+1) + (1 - relaxU)/relaxU...
                                *AP_u(2:nx,2:ny+1).*umac_temp(2:nx,2:ny+1);
            
            idx = 0; % idx = 0;
            stride = nx - 1;
            for j = 2:ny + 1
                for i = 2:nx
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        A(idx, idxS) = -AS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxW) = -AW_u(i, j);
                    end
                    A(idx, idx) = AP_u(i, j) / relaxU;

                    if (i ~= nx)
                        A(idx, idxE) = -AE_u(i, j);
                    end
                    if (j ~= ny+1)
                        A(idx,idxN) = -AN_u(i,j);
                    end
                end
            end
        case 'QUICK'
            % Velocity matrices using QUICK scheme.
            AE_u(2:nx,2:ny+1) = De_u + max(-me_u,0);
            AW_u(2:nx,2:ny+1) = Dw_u + max(mw_u,0);
            AN_u(2:nx,2:ny+1) = Dn_u + max(-mn_u,0);
            AS_u(2:nx,2:ny+1) = Ds_u + max(ms_u,0);
            AP_u = AE_u + AW_u + AN_u + AS_u + rho*dx*dy/dt;
            % Save 'AP_u' as 'AUp' in 'GlobalsSIMPLE.m'.
            AUp(2:nx,2:ny+1) = AP_u(2:nx,2:ny+1);

            %% Get the AP_u and rhs in the AP_u*ustar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs = -dy*diff(pstar) + rho*dx*dy/dt * umac_temp;
            % (2) E
            % rhs = pressure term + AE_u*uE
            rhs(nx, 2:ny+1) = rhs(nx,2:ny+1) + AE_u(nx,2:ny+1).*umac(nx+1,2:ny+1);
            AE_u(nx, 2:ny+1) = 0;
            % (3) W
            % rhs = pressure term + AE_u*uE + AW_u*uW
            rhs(2, 2:ny+1) = rhs(2,2:ny+1) + AW_u(2,2:ny+1).*umac(1,2:ny+1);
            AW_u(2, 2:ny+1) = 0;
            % (4) N
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN
            rhs(2:nx, ny+1) = rhs(2:nx,ny+1) + AN_u(2:nx,ny+1).*umac(2:nx,ny+2);
            AN_u(2:nx, ny+1) = 0;
            % (5) S
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS
            rhs(2:nx, 2) = rhs(2:nx,2) + AS_u(2:nx,2).*umac(2:nx,1);
            AS_u(2:nx, 2) = 0;
            % (6) Flux limiter
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS 
            %       + flux limiter
            re = zeros(1, ny); % zeros(nx, 2:ny+1);
            rw = zeros(1, ny); % zeros(2, 2:ny+1);
            rn = zeros(nx-1, 1); % zeros(2:nx, ny+1);
            rs = zeros(nx-1, 1); % zeros(2:nx, 2);
            for t = 2:ny+1
                re(nx, t)  = (umac(nx, t+1) - umac(nx, t))/(umac(nx, t) - umac(nx, t-1));
                rw(2, t)   = (umac(2, t+1) - umac(2, t))/(umac(2, t) - umac(2, t-1));
                psi_re     = max([0, min([2*re(nx, t), 0.75*re(nx, t)+0.25, 2])]);
                psi_rw     = max([0, min([2*rw(2, t), 0.75*rw(2, t)+0.25, 2])]);
                rhs(nx, t) = rhs(nx, t) - rho*dy*0.5*(umac_temp(nx,t) + umac_temp(nx,t-1))*psi_re/2*(umac(nx, t) - umac(nx, t - 1));
                rhs(2, t)  = rhs(2, t)  + rho*dy*0.5*(umac_temp(2,t) + umac_temp(2,t-1))*psi_rw/2*(umac(2, t) - umac(2, t - 1));
            end
            for s = 2:nx
                rn(s, ny+1) = (umac(s + 1, ny+1) - umac(s, ny+1))/(umac(s, ny+1) - umac(s - 1, ny+1));
                rs(s, 2)    = (umac(s+1, 2) - umac(s, 2))/(umac(s, 2) - umac(s-1, 2));
                psi_rn     = max([0, min([2*rn(s, ny+1), 0.75*rn(s, ny+1)+0.25, 2])]);
                psi_rs     = max([0, min([2*rs(s, 2), 0.75*rs(s, 2)+0.25, 2])]);
                rhs(s,ny+1) = rhs(s,ny+1) - rho*dy*0.5*(umac_temp(s, ny+1) + umac_temp(s-1, ny+1))*psi_rn/2*(umac(s, ny+1) - umac(s-1, ny+1));
                rhs(s,2)    = rhs(s,2)    + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(umac(s, 2) - umac(s-1, 2));
            end
            % (7) Relaxation
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS 
            %       + flux limiter + relaxation
            rhs(2:nx, 2:ny+1) = rhs(2:nx,2:ny+1) + (1 - relaxU)/relaxU...
                                *AP_u(2:nx,2:ny+1).*umac_temp(2:nx,2:ny+1);
            idx = 0;
            stride = nx - 1;
            for j = 2:ny + 1
                for i = 2:nx
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        A(idx, idxS) = -AS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxW) = -AW_u(i, j);
                    end
                    A(idx, idx) = AP_u(i, j) / relaxU;

                    if (i ~= nx)
                        A(idx, idxE) = -AE_u(i, j);
                    end
                    if (j ~= ny+1)
                        A(idx,idxN) = -AN_u(i,j);
                    end
                end
            end
            %{
            A = spalloc(N,N,7*N);
            % Velocity matrices using QUICK scheme.
            AE_u(2:nx, 2:ny+1) = De_u - 3/8*me_u;
            AW_u(2:nx, 2:ny+1) = Dw_u + 3/4*mw_u + 1/8*me_u;
            AN_u(2:nx, 2:ny+1) = Dn_u - 3/8*mn_u;
            AS_u(2:nx, 2:ny+1) = Ds_u + 3/4*ms_u + 1/8*mn_u;
            AWW_u(2:nx, 2:ny+1) = -1/8*mw_u;
            ASS_u(2:nx, 2:ny+1) = -1/8*ms_u;
            AP_u = AE_u + AW_u + AN_u + AS_u + AWW_u + ASS_u + rho*dx*dy/dt;
            % Save 'AP_u' as 'AUp' in 'GlobalsSIMPLE.m'.
            AUp(2:nx, 2:ny+1) = AP_u(2:nx, 2:ny+1);

            %% Get the AP_u and rhs in the AP_u*ustar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs = -dy*diff(pstar) + rho*dx*dy/dt * umac_temp;
            % (2) E
            % rhs = pressure term + AE_u*uE
            rhs(nx, 2:ny+1) = rhs(nx,2:ny+1) + AE_u(nx,2:ny+1).*umac(nx+1,2:ny+1);
            AE_u(nx, 2:ny+1) = 0;
            % (3) W
            % rhs = pressure term + AE_u*uE + AW_u*uW
            rhs(2, 2:ny+1) = rhs(2,2:ny+1) + AW_u(2,2:ny+1).*umac(1,2:ny+1);
            AW_u(2, 2:ny+1) = 0;
            % (4) N
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN
            rhs(2:nx, ny+1) = rhs(2:nx,ny+1) + AN_u(2:nx,ny+1).*umac(2:nx,ny+2);
            AN_u(2:nx, ny+1) = 0;
            % (5) S
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS
            rhs(2:nx, 2) = rhs(2:nx,2) + AS_u(2:nx,2).*umac(2:nx,1);
            AS_u(2:nx, 2) = 0;
            % (6) WW
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS +
            %       AWW_u*uWW
            rhs(1, 2:ny+1) = rhs(1,2:ny+1) + AWW_u(1,2:ny+1).*umac(1,2:ny+1);
            AWW_u(1, 2:ny+1) = 0;
            % (7) SS
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS +
            %       AWW_u*uWW + ASS_u*uSS
            rhs(2:nx, 1) = rhs(2:nx,1) + ASS_u(2:nx,1).*umac(2:nx,1);
            ASS_u(2:nx, 1) = 0;
            % (8) Flux limiter
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS +
            %       AWW_u*uWW + ASS_u*uSS + flux limiter
            re = zeros(1, ny); % zeros(nx, 2:ny+1);
            rw = zeros(1, ny); % zeros(2, 2:ny+1);
            rn = zeros(nx-1, 1); % zeros(2:nx, ny+1);
            rs = zeros(nx-1, 1); % zeros(2:nx, 2);
            for t = 2:ny+1
                re(nx, t)  = (umac(nx, t+1) - umac(nx, t))/(umac(nx, t) - umac(nx, t-1));
                rw(2, t)   = (umac(2, t+1) - umac(2, t))/(umac(2, t) - umac(2, t-1));
                psi_re     = 0.75*re(nx, t) + 0.25;
                psi_rw     = 0.75*rw(2, t) + 0.25;
                rhs(nx, t) = rhs(nx, t) - rho*dy*0.5*(umac_temp(nx,t) + umac_temp(nx,t-1))*psi_re/2*(umac(nx, t) - umac(nx, t - 1));
                rhs(2, t)  = rhs(2, t)  + rho*dy*0.5*(umac_temp(2,t) + umac_temp(2,t-1))*psi_rw/2*(umac(2, t) - umac(2, t - 1));
            end
            for s = 2:nx
                rn(s, ny+1) = (umac(s + 1, ny+1) - umac(s, ny+1))/(umac(s, ny+1) - umac(s - 1, ny+1));
                rs(s, 2)    = (umac(s+1, 2) - umac(s, 2))/(umac(s, 2) - umac(s-1, 2));
                psi_rn      = 0.75*re(s, ny+1) + 0.25;
                psi_rs      = 0.75*rs(s, 2)    + 0.25;
                rhs(s,ny+1) = rhs(s,ny+1) - rho*dy*0.5*(umac_temp(s, ny+1) + umac_temp(s-1, ny+1))*psi_rn/2*(umac(s, ny+1) - umac(s-1, ny+1));
                rhs(s,2)    = rhs(s,2)    + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(umac(s, 2) - umac(s-1, 2));
            end
            % (9) Relaxation
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS +
            %       AWW_u*uWW + ASS_u*uSS + flux limiter + relaxation
            rhs(2:nx, 2:ny+1) = rhs(2:nx,2:ny+1) + (1 - relaxU)/relaxU...
                                *AP_u(2:nx,2:ny+1).*umac_temp(2:nx,2:ny+1);
            
            idx = 0;
            stride = nx - 1;
            for j = 2: 3
                for i = 2: nx
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;

                    if (j ~= 2)
                        A(idx, idxS) = -AS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxW) = -AW_u(i, j);
                    end
                    %{
                    if (j ~= 2)
                        A(idx, idxSS) = -ASS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxWW) = -AWW_u(i, j);
                    end
                    %}
                    A(idx,idx) = AP_u(i,j) / relaxU;

                    if (i ~= nx)
                        A(idx,idxE) = -AE_u(i,j);
                    end
                    if (j ~= ny+1)
                        A(idx,idxN) = -AN_u(i,j);
                    end
                end
            end

            %idx = 1; %idx = 0; 
            %stride = nx - 1;
            for j = 4: ny+1
                for i = 2: nx
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    idxWW = idx - 2;
                    idxSS = idx - stride - 1;
                    if (j ~= 2)
                        A(idx, idxS) = -AS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxW) = -AW_u(i, j);
                    end
                    if (j ~= 2)
                        A(idx, idxSS) = -ASS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxWW) = -AWW_u(i, j);
                    end
                    %{
                    if idx > 2 
                        % if idx = 3, then j = 2, i = 4, but j ~=2, so that
                        % idx = 81, when j = 3, i = 2
                        if (j ~= 2)
                            A(idx, idxSS) = -ASS_u(i, j);
                        end
                        if (i ~= 2)
                            A(idx, idxWW) = -AWW_u(i, j);
                        end
                    end
                    %}
                    A(idx,idx) = AP_u(i,j) / relaxU;

                    if (i ~= nx)
                        A(idx,idxE) = -AE_u(i,j);
                    end
                    if (j ~= ny+1)
                        A(idx,idxN) = -AN_u(i,j);
                    end
                end
            end
            %}
        case 'MUSCL'
            % Velocity matrices using MUSCL scheme.
            AE_u(2:nx,2:ny+1) = De_u + max(-me_u,0);
            AW_u(2:nx,2:ny+1) = Dw_u + max(mw_u,0);
            AN_u(2:nx,2:ny+1) = Dn_u + max(-mn_u,0);
            AS_u(2:nx,2:ny+1) = Ds_u + max(ms_u,0);
            AP_u = AE_u + AW_u + AN_u + AS_u + rho*dx*dy/dt;
            % Save 'AP_u' as 'AUp' in 'GlobalsSIMPLE.m'.
            AUp(2:nx,2:ny+1) = AP_u(2:nx,2:ny+1);

            %% Get the AP_u and rhs in the AP_u*ustar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs = -dy*diff(pstar) + rho*dx*dy/dt * umac_temp;
            % (2) E
            % rhs = pressure term + AE_u*uE
            rhs(nx, 2:ny+1) = rhs(nx,2:ny+1) + AE_u(nx,2:ny+1).*umac(nx+1,2:ny+1);
            AE_u(nx, 2:ny+1) = 0;
            % (3) W
            % rhs = pressure term + AE_u*uE + AW_u*uW
            rhs(2, 2:ny+1) = rhs(2,2:ny+1) + AW_u(2,2:ny+1).*umac(1,2:ny+1);
            AW_u(2, 2:ny+1) = 0;
            % (4) N
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN
            rhs(2:nx, ny+1) = rhs(2:nx,ny+1) + AN_u(2:nx,ny+1).*umac(2:nx,ny+2);
            AN_u(2:nx, ny+1) = 0;
            % (5) S
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS
            rhs(2:nx, 2) = rhs(2:nx,2) + AS_u(2:nx,2).*umac(2:nx,1);
            AS_u(2:nx, 2) = 0;
            % (6) Flux limiter
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS 
            %       + flux limiter
            re = zeros(1, ny); % zeros(nx, 2:ny+1);
            rw = zeros(1, ny); % zeros(2, 2:ny+1);
            rn = zeros(nx-1, 1); % zeros(2:nx, ny+1);
            rs = zeros(nx-1, 1); % zeros(2:nx, 2);
            for t = 2:ny+1
                re(nx, t)  = (umac(nx, t+1) - umac(nx, t))/(umac(nx, t) - umac(nx, t-1));
                rw(2, t)   = (umac(2, t+1) - umac(2, t))/(umac(2, t) - umac(2, t-1));
                psi_re     = max(0, min([2*re(nx, t), (re(nx, t) + 1)/2, 2]));
                psi_rw     = max(0, min([2*rw(2, t), (rw(2, t) + 1)/2, 2]));
                rhs(nx, t) = rhs(nx, t) - rho*dy*0.5*(umac_temp(nx,t) + umac_temp(nx,t-1))*psi_re/2*(umac(nx, t) - umac(nx, t - 1));
                rhs(2, t)  = rhs(2, t)  + rho*dy*0.5*(umac_temp(2,t) + umac_temp(2,t-1))*psi_rw/2*(umac(2, t) - umac(2, t - 1));
            end
            for s = 2:nx
                rn(s, ny+1) = (umac(s + 1, ny+1) - umac(s, ny+1))/(umac(s, ny+1) - umac(s - 1, ny+1));
                rs(s, 2)    = (umac(s+1, 2) - umac(s, 2))/(umac(s, 2) - umac(s-1, 2));
                psi_rn      = max(0, min([2*rn(s, ny+1), (rn(s, ny+1) + 1)/2, 2]));
                psi_rs      = max(0, min([2*rs(s, 2), (rs(s, 2) + 1)/2, 2]));
                rhs(s,ny+1) = rhs(s,ny+1) - rho*dy*0.5*(umac_temp(s, ny+1) + umac_temp(s-1, ny+1))*psi_rn/2*(umac(s, ny+1) - umac(s-1, ny+1));
                rhs(s,2)    = rhs(s,2)    + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(umac(s, 2) - umac(s-1, 2));
            end
            % (7) Relaxation
            % rhs = pressure term + AE_u*uE + AW_u*uW + AN_u*uN + AS_u*uS 
            %       + flux limiter + relaxation
            rhs(2:nx, 2:ny+1) = rhs(2:nx,2:ny+1) + (1 - relaxU)/relaxU...
                                *AP_u(2:nx,2:ny+1).*umac_temp(2:nx,2:ny+1);
            idx = 0;
            stride = nx - 1;
            for j = 2:ny + 1
                for i = 2:nx
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        A(idx, idxS) = -AS_u(i, j);
                    end
                    if (i ~= 2)
                        A(idx, idxW) = -AW_u(i, j);
                    end
                    A(idx, idx) = AP_u(i, j) / relaxU;

                    if (i ~= nx)
                        A(idx, idxE) = -AE_u(i, j);
                    end
                    if (j ~= ny+1)
                        A(idx,idxN) = -AN_u(i,j);
                    end
                end
            end
        otherwise
            msg = 'Error occurred.';
            error(msg);
    end
    
    %% Solve ustar from AP_u*ustar = rhs
    rhs = reshape(rhs(2:nx,2:ny+1),N,1);
    sol = A \ rhs;
    
    % get u*
    ustar = umac;
    ustar(2:nx,2:ny+1) = reshape(sol,nx-1,ny);
    % enforce BC
    ustar = ApplyBC_umac(ustar);
    
    return
end