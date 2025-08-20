function [ phi ] = ApplyBC_pres(phi)
    %% Import globals from 'GlobalsSIMPLE.m'
    GlobalsSIMPLE;
    %% Set the homogeneous BCs.
    phi(1,2:ny+1) = phi(2,2:ny+1);
    phi(nx+2,2:ny+1) = phi(nx-1,2:ny+1);
    phi(2:nx+1,1) = phi(2:nx+1,2);
    phi(2:nx+1,ny+2) = phi(2:nx+1,ny+1);
    return
end