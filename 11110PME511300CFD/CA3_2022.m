% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
% C                                                                 C
% C           11110PME511300 Computational Fluid Dynamics           C
% C                                                                 C
% C                   Finite Volume Method Program                  C
% C                       In SIMPLE Algorithm                       C
% C                      for Lid-Driven Cavity                      C
% C                                                                 C
% C                         Cheng-Chun Yang                         C
% C                                                                 C
% C                       ALL RIGHTS RESERVED                       C
% C                                                                 C
% C           DEPARTMENT OF POWER MECHANICAL ENGINEERING            C
% C             NATIONAL TSING HUA UNIVERSITY, TAIWAN               C
% C                                                                 C
% C                          Jan, 07, 2023                          C
% C                                                                 C
% CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
% *******************************************************************
%% Clear the previous runs
clear; clc; close all; format shortE;
%% Set the fonts to LaTeX
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex');
%% Add, remove paths
% Add path to subroutine folder.
path_MAINCODE = 'MAIN_LIBRARY'; 
rmpath(path_MAINCODE); addpath(path_MAINCODE);
%% Set up Parameters
% Import globals from 'GlobalsSIMPLE.m'.
GlobalsSIMPLE;
% Max iteration steps.
max_steps = 5000;
% Domain lengths.
Lx = 1; Ly = 1;
% Mesh grids number.
nmesh = [81, 161];
for ncell = nmesh
    % Spacial discretization.
    nx = ncell; ny = ncell; dx = Lx / nx; dy = Ly / ny;
    % Relaxation parameter.
    relax = 0.8; relaxU = relax; relaxV = relax; relaxP = 1 - relax;
    % Top lid velocity.
    ULid = 1;
    % Fluid density
    rho = 1;
    % Reynolds number.
    Reynolds = [100, 1000, 5000];
    scheme = {'CD','QUICK','MUSCL'};
    for Re = Reynolds
        for sch = scheme
            % Fluid viscosity.
            visc = rho*Lx*ULid / Re;
            % Length for each temporal discretization.
            dt = Inf;
            % Declare the storages.
            u_mac = zeros(nx+1,ny+2);
            v_mac = zeros(nx+2,ny+1);
            u_star = zeros(nx+1,ny+2);
            v_star = zeros(nx+2,ny+1);
            p = zeros(nx+2, ny+2);
            pstar = zeros(nx+2, ny+2);
            pdash = zeros(nx+2, ny+2);
            AUp = zeros(nx+1, ny+2);
            AVp = zeros(nx+2, ny+1);
            
            u_mac = ApplyBC_umac(u_mac);
            v_mac = ApplyBC_vmac(v_mac);
            %{
            umac_temp = u_mac;
            vmac_temp = v_mac;
            pold = p;
            %}
            %% Start the iteration.
            disp(['Steady Lid-driven cavity, Re=', num2str(Re), ' using '] + string(sch));
            for istep = 1:max_steps
                umac_temp = u_mac;
                vmac_temp = v_mac;
                pold = pstar;
                % Solve u* and v*.
                % Solve ustar from AP_u*ustar = rhs => ustar = AP_u\rhs
                u_star = Solve_ustar(u_mac,umac_temp,vmac_temp,pstar);
                % Solve vstar from AP_v*vstar = rhs => vstar = AP_v\rhs
                v_star = Solve_vstar(v_mac,umac_temp,vmac_temp,pstar);
                % Solve pressure.
                [u_mac,v_mac,pstar] = Solve_pressure(u_star,v_star,pstar);
                % Solve u' and v'
                if (mod(istep,10) == 0)
                    ucorr = u_mac - umac_temp;
                    vcorr = v_mac - vmac_temp;
                    % pcorr = pstar - pold;
                    tol_abs = 1e-4*ULid;
                    corr = max([norm(ucorr(:),inf), norm(vcorr(:),inf)]);
                    % Display for each 10 steps.
                    disp(['step=', int2str(istep),', corr=', num2str(corr)]);
                    % Convergence criterion.
                    if (corr < tol_abs)
                        break;
                    end
                end
            end
            xcs = linspace(dx/2,Lx-dx/2,nx); ycs = linspace(dy/2,Ly-dy/2,ny);
            %% Obtain solutions
            % x-directional Velocity
            u_sol = 0.5 * (u_mac(1:nx,2:ny+1) + u_mac(2:nx+1,2:ny+1));
            % y-directional Velocity
            v_sol = 0.5 * (v_mac(2:nx+1,1:ny) + v_mac(2:nx+1,2:ny+1));
            % Pressure magnitude
            p_sol = pstar(2:nx+1,2:ny+1);
            % Velocity magnitude
            uMag_sol = sqrt(u_sol'.^2 + v_sol'.^2);
            %% Plot figures
            PostProcess(xcs, ycs, u_sol, v_sol, p_sol, uMag_sol);
        end
    end
end