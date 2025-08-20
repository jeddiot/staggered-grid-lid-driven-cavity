function [ vstar ] = Solve_vstar (vmac,umac_temp,vmac_temp,pstar)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    N = nx * (ny-1);
    % H-operator
    H = spalloc(N,N,5*N);
    % rhs = zeros(N,1);
    rhs = zeros(nx+2,ny+1);
    % Diffusion coefficients.
    De_v = visc * dy / dx;
    Dw_v = visc * dy / dx;
    Dn_v = visc * dx / dy;
    Ds_v = visc * dx / dy;
    % Mass flow rates.
    me_v = rho * dy * 0.5*(umac_temp(2:nx+1, 2:ny) + umac_temp(2:nx+1, 3:ny+1));
    mw_v = rho * dy * 0.5*(umac_temp(1:nx, 2:ny)   + umac_temp(1:nx, 3:ny+1));
    mn_v = rho * dx * 0.5*(vmac_temp(2:nx+1, 2:ny) + vmac_temp(2:nx+1, 3:ny+1));
    ms_v = rho * dx * 0.5*(vmac_temp(2:nx+1, 2:ny) + vmac_temp(2:nx+1, 1:ny-1));
    % The velocity matrices from AP_v*vp = AE_v*ve + AW_v*vw + AN_v*vn + AS_v*vs;
    AE_v = zeros(nx+2,ny+1);
    AW_v = zeros(nx+2,ny+1);
    AN_v = zeros(nx+2,ny+1);
    AS_v = zeros(nx+2,ny+1);
    AWW_v = zeros(nx+2,ny+1);
    ASS_v = zeros(nx+2,ny+1);
    % Change calculation for each scheme.
    switch string(sch)
        case 'CD'
            % Velocity matrices using CD (central differencing) scheme.
            AE_v(2:nx+1, 2:ny) = De_v - me_v/2;
            AW_v(2:nx+1, 2:ny) = Dw_v + mw_v/2;
            AN_v(2:nx+1, 2:ny) = Dn_v - mn_v/2;
            AS_v(2:nx+1, 2:ny) = Ds_v + ms_v/2;
            AP_v = AE_v + AW_v + AN_v + AS_v + rho*dx*dy/dt;
            % Save 'AP_v' as 'AVp' in 'GlobalsSIMPLE.m'.
            AVp(2:nx+1,2:ny) = AP_v(2:nx+1,2:ny);

            %% Get the AP_v and rhs in the AP_v*vstar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs(2:nx+1,2:ny) = -diff(pstar(2:nx+1,2:ny+1)')' * dx + ...
                rho*dx*dy/dt * vmac_temp(2:nx+1,2:ny);
            % (2) E
            % rhs = pressure term + AE_v*vE
            rhs(nx+1,2:ny) = rhs(nx+1,2:ny) + AE_v(nx+1,2:ny) .* vmac(nx+2,2:ny);
            AE_v(nx+1,2:ny) = 0;
            % (3) W
            % rhs = pressure term + AE_v*vE + AW_v*vW
            rhs(2,2:ny) = rhs(2,2:ny) + AW_v(2,2:ny) .* vmac(1,2:ny);
            AW_v(2,2:ny) = 0;
            % (4) N
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN
            rhs(2:nx+1,ny) = rhs(2:nx+1,ny) + AN_v(2:nx+1,ny) .* vmac(2:nx+1,ny+1);
            AN_v(2:nx+1,ny) = 0;
            % (5) S
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            rhs(2:nx+1,2) = rhs(2:nx+1,2) + AS_v(2:nx+1,2) .* vmac(2:nx+1,1);
            AS_v(2:nx+1,2) = 0;
            % (6) Relaxation
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            %       + relaxation
            rhs(2:nx+1,2:ny) = rhs(2:nx+1,2:ny) + ...
                (1-relaxV)/relaxV * AP_v(2:nx+1,2:ny) .* vmac_temp(2:nx+1,2:ny);
            rhs = reshape(rhs(2:nx+1,2:ny), N,1);
            
            idx = 0;
            stride = nx;
            for j = 2:ny
                for i = 2:nx+1
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        H(idx,idxS) = -AS_v(i,j);
                    end
                    if (i ~= 2)
                        H(idx,idxW) = -AW_v(i,j);
                    end
                    H(idx,idx) = AP_v(i,j) / relaxV;
                    if (i ~= nx+1)
                        H(idx,idxE) = -AE_v(i,j);
                    end
                    if (j ~= ny)
                        H(idx,idxN) = -AN_v(i,j);
                    end
                end
            end
        case 'QUICK'
            % Velocity matrices using QUICK scheme.
            AE_v(2:nx+1,2:ny) = De_v + max(-me_v,0);
            AW_v(2:nx+1,2:ny) = Dw_v + max(mw_v,0);
            AN_v(2:nx+1,2:ny) = Dn_v + max(-mn_v,0);
            AS_v(2:nx+1,2:ny) = Ds_v + max(ms_v,0);
            AP_v = AE_v + AW_v + AN_v + AS_v + rho*dx*dy/dt;
            % Save 'AP_v' as 'AVp' in 'GlobalsSIMPLE.m'.
            AVp(2:nx+1,2:ny) = AP_v(2:nx+1,2:ny);

            %% Get the AP_v and rhs in the AP_v*vstar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs(2:nx+1, 2:ny) = -diff(pstar(2:nx+1, 2:ny+1)')' * dx + ...
                                rho*dx*dy/dt * vmac_temp(2:nx+1, 2:ny);
            % (2) E
            % rhs = pressure term + AE_v*vE
            rhs(nx+1, 2:ny) = rhs(nx+1, 2:ny) + AE_v(nx+1, 2:ny) .* vmac(nx+2, 2:ny);
            AE_v(nx+1, 2:ny) = 0;
            % (3) W
            % rhs = pressure term + AE_v*vE + AW_v*vW
            rhs(2, 2:ny) = rhs(2, 2:ny) + AW_v(2, 2:ny) .* vmac(1, 2:ny);
            AW_v(2, 2:ny) = 0;
            % (4) N
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN
            rhs(2:nx+1,ny) = rhs(2:nx+1,ny) + AN_v(2:nx+1,ny) .* vmac(2:nx+1,ny+1);
            AN_v(2:nx+1,ny) = 0;
            % (5) S
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            rhs(2:nx+1,2) = rhs(2:nx+1,2) + AS_v(2:nx+1,2) .* vmac(2:nx+1,1);
            AS_v(2:nx+1,2) = 0;
            % (6) Flux limiter
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            %       + flux limiter
            re = zeros(1, ny-1);
            rw = zeros(1, ny-1);
            rn = zeros(nx, 1);
            rs = zeros(nx, 1);
            for t = 2:ny
                re(nx + 1, t) = (vmac(nx + 1, t+1) - vmac(nx + 1, t))/(vmac(nx + 1, t) - vmac(nx + 1, t-1));
                rw(2, t) = (vmac(2, t+1) - vmac(2, t))/(vmac(2, t) - vmac(2, t-1));
                psi_re = max([0, min([2*re(nx+1, t), 0.75*re(nx+1, t)+0.25, 2])]);
                psi_rw = max([0, min([2*rw(2, t), 0.75*rw(2, t)+0.25, 2])]);
                rhs(nx+1, t) = rhs(nx+1, t) - rho*dy*0.5*(umac_temp(nx+1, t) + umac_temp(nx+1, t-1))*psi_re/2*(vmac(nx + 1, t) - vmac(nx + 1, t - 1));
                rhs(2, t) = rhs(2, t) + rho*dy*0.5*(umac_temp(2, t) + umac_temp(2, t-1))*psi_rw/2*(vmac(2, t) - vmac(2, t - 1));
            end
            for s = 2:nx+1
                rn(s, ny) = (vmac(s + 1, ny) - vmac(s, ny))/(vmac(s, ny) - vmac(s - 1, ny));
                rs(s, 2) = (vmac(s+1, 2) - vmac(s, 2))/(vmac(s, 2) - vmac(s-1, 2));
                psi_rn = max([0, min([2*rn(s, ny), 0.75*rn(s, ny)+0.25, 2])]);
                psi_rs = max([0, min([2*rn(s, 2), 0.75*rn(s, 2)+0.25, 2])]);
                rhs(s,ny) = rhs(s,ny) - rho*dy*0.5*(umac_temp(s, ny) + umac_temp(s-1, ny))*psi_rn/2*(vmac(s, ny) - vmac(s-1, ny));
                rhs(s,2) = rhs(s,2) + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(vmac(s, 2) - vmac(s-1, 2));
            end
            % (7) Relaxation
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            %       + flux limiter + relaxation
            rhs(2:nx+1,2:ny) = rhs(2:nx+1,2:ny) + ...
                               (1-relaxV)/relaxV * AP_v(2:nx+1,2:ny) .* vmac_temp(2:nx+1,2:ny);
            rhs = reshape(rhs(2:nx+1,2:ny), N,1);
            
            idx = 0; %idx = 0;
            stride = nx;
            for j = 2:ny
                for i = 2:nx+1
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        H(idx,idxS) = -AS_v(i,j);
                    end
                    if (i ~= 2)
                        H(idx,idxW) = -AW_v(i,j);
                    end
                    H(idx,idx) = AP_v(i,j) / relaxV;
                    if (i ~= nx+1)
                        H(idx,idxE) = -AE_v(i,j);
                    end
                    if (j ~= ny)
                        H(idx,idxN) = -AN_v(i,j);
                    end
                end
            end
            %{
            H = spalloc(N,N,7*N);
            % Velocity matrices using QUICK scheme.
            AE_v(2:nx+1, 2:ny) = De_v - 3/8*me_v;
            AW_v(2:nx+1, 2:ny) = Dw_v + 3/4*mw_v + 1/8*me_v;
            AN_v(2:nx+1, 2:ny) = Dn_v - 3/8*mn_v;
            AS_v(2:nx+1, 2:ny) = Ds_v + 3/4*ms_v + 1/8*mn_v;
            AWW_v(2:nx+1, 2:ny) = -1/8*mw_v;
            ASS_v(2:nx+1, 2:ny) = -1/8*ms_v;
            AP_v = AE_v + AW_v + AN_v + AS_v + AWW_v + ASS_v + rho*dx*dy/dt;
            % Save 'AP_u' as 'AUp' in 'GlobalsSIMPLE.m'.
            AVp(2:nx+1, 2:ny) = AP_v(2:nx+1, 2:ny);

            %% Get the AP_v and rhs in the AP_v*vstar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs(2:nx+1, 2:ny) = -diff(pstar(2:nx+1, 2:ny+1)')' * dx + ...
                                rho*dx*dy/dt * vmac_temp(2:nx+1, 2:ny);
            % (2) E
            % rhs = pressure term + AE_v*vE
            rhs(nx+1, 2:ny) = rhs(nx+1, 2:ny) + AE_v(nx+1, 2:ny) .* vmac(nx+2, 2:ny);
            AE_v(nx+1, 2:ny) = 0;
            % (3) W
            % rhs = pressure term + AE_v*vE + AW_v*vW
            rhs(2, 2:ny) = rhs(2, 2:ny) + AW_v(2, 2:ny) .* vmac(1, 2:ny);
            AW_v(2, 2:ny) = 0;
            % (4) N
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN
            rhs(2:nx+1,ny) = rhs(2:nx+1,ny) + AN_v(2:nx+1,ny) .* vmac(2:nx+1,ny+1);
            AN_v(2:nx+1,ny) = 0;
            % (5) S
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            rhs(2:nx+1,2) = rhs(2:nx+1,2) + AS_v(2:nx+1,2) .* vmac(2:nx+1,1);
            AS_v(2:nx+1,2) = 0;
            % (6) WW
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS +
            %       AWW_v*vWW
            rhs(1, 2:ny) = rhs(1, 2:ny) + AW_v(1, 2:ny) .* vmac(1, 2:ny);
            AW_v(1, 2:ny) = 0;
            % (6) SS
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS +
            %       AWW_v*vWW + ASS_v*vSS
            rhs(2:nx+1,1) = rhs(2:nx+1,1) + ASS_v(2:nx+1,1).* vmac(2:nx+1,1);
            ASS_v(2:nx+1,1) = 0;

            % (7) Flux limiter
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS +
            %       AWW_v*vWW + ASS_v*vSS+ flux limiter
            re = zeros(1, ny-1);
            rw = zeros(1, ny-1);
            rn = zeros(nx, 1);
            rs = zeros(nx, 1);
            for t = 2:ny
                re(nx + 1, t) = (vmac(nx + 1, t+1) - vmac(nx + 1, t))/(vmac(nx + 1, t) - vmac(nx + 1, t-1));
                rw(2, t) = (vmac(2, t+1) - vmac(2, t))/(vmac(2, t) - vmac(2, t-1));
                psi_re = 0.75*re(nx + 1, t) + 0.25;
                psi_rw = 0.75*re(2, t) + 0.25;
                rhs(nx+1, t) = rhs(nx+1, t) - rho*dy*0.5*(umac_temp(nx+1, t) + umac_temp(nx+1, t-1))*psi_re/2*(vmac(nx + 1, t) - vmac(nx + 1, t - 1));
                rhs(2, t) = rhs(2, t) + rho*dy*0.5*(umac_temp(2, t) + umac_temp(2, t-1))*psi_rw/2*(vmac(2, t) - vmac(2, t - 1));
            end
            for s = 2:nx+1
                rn(s, ny) = (vmac(s + 1, ny) - vmac(s, ny))/(vmac(s, ny) - vmac(s - 1, ny));
                rs(s, 2) = (vmac(s+1, 2) - vmac(s, 2))/(vmac(s, 2) - vmac(s-1, 2));
                psi_rn = 0.75*re(s, ny) + 0.25;
                psi_rs = 0.75*re(s, 2) + 0.25;
                rhs(s,ny) = rhs(s,ny) - rho*dy*0.5*(umac_temp(s, ny) + umac_temp(s-1, ny))*psi_rn/2*(vmac(s, ny) - vmac(s-1, ny));
                rhs(s,2) = rhs(s,2) + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(vmac(s, 2) - vmac(s-1, 2));
            end
            % (8) Relaxation
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS +
            %       AWW_v*vWW + ASS_v*vSS+ flux limiter + relaxation
            rhs(2:nx+1,2:ny) = rhs(2:nx+1,2:ny) + ...
                               (1-relaxV)/relaxV * AP_v(2:nx+1,2:ny) .* vmac_temp(2:nx+1,2:ny);
            rhs = reshape(rhs(2:nx+1,2:ny), N,1);
            
            idx = 0;
            stride = nx;
            for j = 2:3
                for i = 2:nx+1
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        H(idx,idxS) = -AS_v(i,j);
                    end
                    if (i ~= 2)
                        H(idx,idxW) = -AW_v(i,j);
                    end
                    %{
                    if idx > 2
                        if (j ~= 2)
                            H(idx,idxSS) = -ASS_v(i,j);
                        end
                        if (i ~= 2)
                            H(idx,idxWW) = -AWW_v(i,j);
                        end
                    end
                    %}
                    H(idx,idx) = AP_v(i,j) / relaxV;

                    if (i ~= nx+1)
                        H(idx,idxE) = -AE_v(i,j);
                    end
                    if (j ~= ny)
                        H(idx,idxN) = -AN_v(i,j);
                    end
                end
            end
            for j = 4:ny
                for i = 2:nx+1
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    idxWW = idx - 2;
                    idxSS = idx - stride - 1;
                    
                    if (j ~= 2)
                        H(idx,idxS) = -AS_v(i,j);
                    end
                    if (i ~= 2)
                        H(idx,idxW) = -AW_v(i,j);
                    end

                    if idx > 2
                        if (j ~= 2)
                            H(idx,idxSS) = -ASS_v(i,j);
                        end
                        if (i ~= 2)
                            H(idx,idxWW) = -AWW_v(i,j);
                        end
                    end
                    
                    H(idx,idx) = AP_v(i,j) / relaxV;
                    if (i ~= nx+1)
                        H(idx,idxE) = -AE_v(i,j);
                    end
                    if (j ~= ny)
                        H(idx,idxN) = -AN_v(i,j);
                    end
                end
            end
            %}
        case 'MUSCL'
            % Velocity matrices using MUSCL scheme.
            AE_v(2:nx+1,2:ny) = De_v + max(-me_v,0);
            AW_v(2:nx+1,2:ny) = Dw_v + max(mw_v,0);
            AN_v(2:nx+1,2:ny) = Dn_v + max(-mn_v,0);
            AS_v(2:nx+1,2:ny) = Ds_v + max(ms_v,0);
            AP_v = AE_v + AW_v + AN_v + AS_v + rho*dx*dy/dt;
            % Save 'AP_v' as 'AVp' in 'GlobalsSIMPLE.m'.
            AVp(2:nx+1,2:ny) = AP_v(2:nx+1,2:ny);

            %% Get the AP_v and rhs in the AP_v*vstar = rhs
            % (1) Pressure term
            % rhs = pressure term 
            rhs(2:nx+1, 2:ny) = -diff(pstar(2:nx+1, 2:ny+1)')' * dx + ...
                                rho*dx*dy/dt * vmac_temp(2:nx+1, 2:ny);
            % (2) E
            % rhs = pressure term + AE_v*vE
            rhs(nx+1, 2:ny) = rhs(nx+1, 2:ny) + AE_v(nx+1, 2:ny) .* vmac(nx+2, 2:ny);
            AE_v(nx+1, 2:ny) = 0;
            % (3) W
            % rhs = pressure term + AE_v*vE + AW_v*vW
            rhs(2, 2:ny) = rhs(2, 2:ny) + AW_v(2, 2:ny) .* vmac(1, 2:ny);
            AW_v(2, 2:ny) = 0;
            % (4) N
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN
            rhs(2:nx+1,ny) = rhs(2:nx+1,ny) + AN_v(2:nx+1,ny) .* vmac(2:nx+1,ny+1);
            AN_v(2:nx+1,ny) = 0;
            % (5) S
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            rhs(2:nx+1,2) = rhs(2:nx+1,2) + AS_v(2:nx+1,2) .* vmac(2:nx+1,1);
            AS_v(2:nx+1,2) = 0;
            % (6) Flux limiter
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            %       + flux limiter
            re = zeros(1, ny-1);
            rw = zeros(1, ny-1);
            rn = zeros(nx, 1);
            rs = zeros(nx, 1);
            for t = 2:ny
                re(nx + 1, t) = (vmac(nx + 1, t+1) - vmac(nx + 1, t))/(vmac(nx + 1, t) - vmac(nx + 1, t-1));
                rw(2, t) = (vmac(2, t+1) - vmac(2, t))/(vmac(2, t) - vmac(2, t-1));
                psi_re = max(0, min([2*re(nx + 1, t), (re(nx + 1, t) + 1)/2, 2]));
                psi_rw = max(0, min([2*rw(2, t), (rw(2, t) + 1)/2, 2]));
                rhs(nx+1, t) = rhs(nx+1, t) - rho*dy*0.5*(umac_temp(nx+1, t) + umac_temp(nx+1, t-1))*psi_re/2*(vmac(nx + 1, t) - vmac(nx + 1, t - 1));
                rhs(2, t) = rhs(2, t) + rho*dy*0.5*(umac_temp(2, t) + umac_temp(2, t-1))*psi_rw/2*(vmac(2, t) - vmac(2, t - 1));
            end
            for s = 2:nx+1
                rn(s, ny) = (vmac(s + 1, ny) - vmac(s, ny))/(vmac(s, ny) - vmac(s - 1, ny));
                rs(s, 2) = (vmac(s+1, 2) - vmac(s, 2))/(vmac(s, 2) - vmac(s-1, 2));
                psi_rn = max(0, min([2*rn(s, ny), (rn(s, ny) + 1)/2, 2]));
                psi_rs = max(0, min([2*rs(s, 2), (rs(s, 2) + 1)/2, 2]));
                rhs(s,ny) = rhs(s,ny) - rho*dy*0.5*(umac_temp(s, ny) + umac_temp(s-1, ny))*psi_rn/2*(vmac(s, ny) - vmac(s-1, ny));
                rhs(s,2) = rhs(s,2) + rho*dy*0.5*(umac_temp(s, 2) + umac_temp(s-1, 2))*psi_rs/2*(vmac(s, 2) - vmac(s-1, 2));
            end
            % (7) Relaxation
            % rhs = pressure term + AE_v*vE + AW_v*vW + AN_v*vN + AS_v*vS
            %       + flux limiter + relaxation
            rhs(2:nx+1,2:ny) = rhs(2:nx+1,2:ny) + ...
                               (1-relaxV)/relaxV * AP_v(2:nx+1,2:ny) .* vmac_temp(2:nx+1,2:ny);
            rhs = reshape(rhs(2:nx+1,2:ny), N,1);
            
            idx = 0; %idx = 0;
            stride = nx;
            for j = 2:ny
                for i = 2:nx+1
                    idx = idx + 1;
                    idxE = idx + 1;
                    idxW = idx - 1;
                    idxN = idx + stride;
                    idxS = idx - stride;
                    
                    if (j ~= 2)
                        H(idx,idxS) = -AS_v(i,j);
                    end
                    if (i ~= 2)
                        H(idx,idxW) = -AW_v(i,j);
                    end
                    H(idx,idx) = AP_v(i,j) / relaxV;
                    if (i ~= nx+1)
                        H(idx,idxE) = -AE_v(i,j);
                    end
                    if (j ~= ny)
                        H(idx,idxN) = -AN_v(i,j);
                    end
                end
            end
        otherwise
            msg = 'Error occurred.';
            error(msg);
    end
    
    %% Solve ustar from AP_v*vstar = rhs
    sol = H \ rhs;
    
    % restore v*
    vstar = vmac;
    vstar(2:nx+1,2:ny) = reshape(sol,nx,ny-1);
    
    % enforce BC
    vstar = ApplyBC_vmac(vstar);
    return
end