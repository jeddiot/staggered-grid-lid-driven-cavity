% Mesh grids number.
global nx ny;
% Domain lengths.
global Lx Ly;
% Length for each grid.
global dx dy;
% Length for each temporal discretization.
global dt;
% Fluid viscosity and density.
global visc rho;
% Top lid velocity.
global ULid;
% Reynolds number.
global Re;
% Relaxation parameter.
global relaxU relaxV relaxP;
% 'AVp' is exactly the 'AP_v' from AP_v*vp = AE_v*ve + AW_v*vw + AN_v*vn + AS_v*vs;
% 'AUp' is exactly the 'AP_u' from AP_u*up = AE_u*ue + AW_u*uw + AN_u*un + AS_u*us
global AUp AVp;
% Figure index.
global fig;
% Save the current scheme
global sch;
%{
global umac vmac;
global ustar vstar;
global p pstar pdash;
global uold vold;
global pold;
% pressure matrix
global Aw As Ap An Ae bp;
%}