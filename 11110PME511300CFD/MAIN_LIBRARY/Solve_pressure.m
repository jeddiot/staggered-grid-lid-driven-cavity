function[unew,vnew,pnew] = Solve_pressure(ustar,vstar,pstar)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    N = nx * ny;
    Lp = spalloc(N,N,N*5);
    % rhs = zeros(N,1);
    rhs = dy * diff(ustar(1:nx+1,2:ny+1)) + dx * diff(vstar(2:nx+1,1:ny+1)')';
    rhs = reshape(-rho * rhs, N, 1);
    
    stride = nx;
    idx = 0;
    for j = 2:ny+1
        for i = 2:nx+1
            idx = idx + 1;
            ide = idx + 1;
            idw = idx - 1;
            idn = idx + stride;
            ids = idx - stride;
            
            % RHS
            % rhs(idx) = -rho*((ustar(i,j)-ustar(i-1,j))*dy+(vstar(i,j)-vstar(i,j-1))*dx);
            
            ape = 0;
            apw = 0;
            apn = 0;
            aps = 0;
            % assemble
            if (j ~= 2)
                aps = rho * dx^2 / AVp(i,j-1);
                Lp(idx,ids) = -aps;
            end
            if (i ~= 2)
                apw = rho * dy^2 / AUp(i-1,j);
                Lp(idx,idw) = -apw;
            end
            if (i ~= nx+1)
                ape = rho * dy^2 / AUp(i,j);
                Lp(idx,ide) = -ape;
            end
            if (j ~= ny+1)
                apn = rho * dx^2 / AVp(i,j);
                Lp(idx,idn) = -apn;
            end
            Lp(idx,idx) = ape + apw + apn + aps;
        end
    end
    
    % inject ref. pressure
    Lp(1,:) = 0;
    Lp(1,1) = 1;
    rhs(1) = 0;
    
    sol = Lp \ rhs;
    
    % restore pressure-correction
    pdash = zeros(nx+2,ny+2);
    pdash(2:nx+1,2:ny+1) = reshape(sol,nx,ny);
    pdash = ApplyBC_pres(pdash);
    
    % new velocity
    udash = zeros(nx+1,ny+2);
    udash(2:nx,2:ny+1) = -dy * diff(pdash(2:nx+1,2:ny+1)) ./ AUp(2:nx,2:ny+1);
    unew = ustar + udash;
    % unew = ustar;
    % for j = 2:ny+1
        % for i = 2:nx
            % udash = -dy / AUp(i,j) * (pdash(i+1,j)-pdash(i,j));
            % unew(i,j) = unew(i,j) + udash;
        % end
    % end
    unew = ApplyBC_umac(unew);
    
    vdash = zeros(nx+2,ny+1);
    vdash(2:nx+1,2:ny) = -dx * diff(pdash(2:nx+1,2:ny+1)')' ./ AVp(2:nx+1,2:ny);
    vnew = vstar + vdash;
    % vnew = vstar;
    % for j = 2:ny
        % for i = 2:nx+1
            % vdash = -dx / AVp(i,j) * (pdash(i,j+1)-pdash(i,j));
            % vnew(i,j) = vnew(i,j) + vdash;
        % end
    % end
    vnew = ApplyBC_vmac(vnew);
    
    % new pressure
    pnew = pstar + relaxP * pdash;
    pnew = ApplyBC_pres(pnew);
    
    return
end