function[v] = ApplyBC_vmac(v)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    %% Set the boundary v-directional velocities.
    % Left boundary.
    v(1, 1:ny+1) = 0; %-v(2, 1:nx+1);
    % Right boundary
    v(nx+2, 1:ny+1) = 0; %-v(ny+1, 1:nx+1);
    % Bottom boundary.
    v(2:nx+1, 1) = 0;
    % Top boundary.
    v(2:nx+1, ny+1) = 0;
    return
end