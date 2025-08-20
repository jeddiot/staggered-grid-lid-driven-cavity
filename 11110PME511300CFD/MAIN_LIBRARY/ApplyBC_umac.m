function[u] = ApplyBC_umac(u)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    %% Set the boundary x-directional velocities.
    % Left boundary.
    u(1,2:ny+1) = 0;
    % Right boundary
    u(nx+1,2:ny+1) = 0;
    % Bottom boundary.
    u(1:nx+1,1) = 0;
    % Top boundary.
    u(1:nx+1,ny+2) = ULid; %2*ULid - u(1:ny+1,nx+1);
    return
end