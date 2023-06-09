%% Briefing
% Author:       Daniel Botezelli

% Bibliography: The Finite Volume Method in Computational Fluid Dynamics:
%               An Advanced Introduction with OpenFOAM® and Matlab, 
%               by F. Moukalled, M. Darwish, L. Mangani. (2015)

% Default:      Channel-Flow at Re = 100. Lid-driven at Re = 1000.

% Date:         Coded and validated in 2023.

%% Pre-Processing
clear; clc; close all; 

% Model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 2;  % 1. Developed Channel-Flow    %%
            % 2. Lid-Driven Cavity         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Refinement
switch model
    case 1
        Nx = 200;
        Ny = 20;
    case 2
        Nx = 150;
        Ny = 150;
end

% Solver
n_sor_m = 500;
n_sor_c = 1000;
max_iteration = 5000;
residual_gap = 20;
tol = 1e-6;
transient = false;

% Matrix
Ac_u = zeros(Ny + 2, Nx + 2);
Ac_v = zeros(Ny + 2, Nx + 2);
Ac_p = zeros(Ny + 2, Nx + 2);
Ae = zeros(Ny + 2, Nx + 2);
Aw = zeros(Ny + 2, Nx + 2);
An = zeros(Ny + 2, Nx + 2);
As = zeros(Ny + 2, Nx + 2);
Bc = zeros(Ny + 2, Nx + 2);
Bc_u = zeros(Ny + 2, Nx + 2);
Bc_v = zeros(Ny + 2, Nx + 2);

% Velocity
u = zeros(Ny + 2, Nx + 2);
v = zeros(Ny + 2, Nx + 2);
uf = zeros(Ny + 2, Nx + 1);
vf = zeros(Ny + 1, Nx + 2);

% Pressure
p = zeros(Ny + 2, Nx + 2);
p_prime = zeros(Ny + 2, Nx + 2);

% Gradient
grad_u_i = zeros(Ny + 2, Nx + 2);
grad_u_j = zeros(Ny + 2, Nx + 2);
grad_v_i = zeros(Ny + 2, Nx + 2);
grad_v_j = zeros(Ny + 2, Nx + 2);
grad_p_i = zeros(Ny + 2, Nx + 2);
grad_p_j = zeros(Ny + 2, Nx + 2);
grad_p_prime_i = zeros(Ny + 2, Nx + 2);
grad_p_prime_j = zeros(Ny + 2, Nx + 2);

% Booleans
channelflow_model = false;
liddriven_model = false;

switch model
    % Channel-flow
    case 1
        U = 1;
        u(2:Ny + 1, 1) = U;
        uf(2:Ny + 1, 1) = U;
        mu = 1e-2;
        rho = 1;
        lambda = 0.7;
        lambda_p = 0.3;
        Lx = 10.0;
        Ly = 1.0;
        dx = Lx/(Nx);
        dy = Ly/(Ny);
        channelflow_model = true;

    % Lid-driven
    case 2 
        U = 1.0;
        u(Ny + 2, 1:Nx + 2) = U;
        uf(Ny + 2, 1:Nx + 1) = U;
        mu = 1e-3;
        rho = 1;
        lambda = 0.7;
        lambda_p = 0.3;
        Lx = 1.0;
        Ly = 1.0;
        dx = Lx/(Nx + 2);
        dy = Ly/(Ny + 2);
        liddriven_model = true;
end

% Meshgrid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

%% SIMPLE algorithm
for i_iteration = 1:max_iteration
    % Gradient
    [grad_u_i, grad_u_j] = gradgauss(Nx, Ny, u, dx, dy, grad_u_i, grad_u_j);
    [grad_v_i, grad_v_j] = gradgauss(Nx, Ny, v, dx, dy, grad_v_i, grad_v_j);
    [grad_p_i, grad_p_j] = gradgauss(Nx, Ny, p, dx, dy, grad_p_i, grad_p_j);

    % Construct momentum eq   
    [Ac_u, Ac_v, Ae, Aw, An, As, Bc_u, Bc_v] = momentum(Nx, Ny, dx, dy, Ac_u, Ac_v, Ae, Aw, An, As, Bc_u, Bc_v, ...
                                                        u, v, uf, vf, grad_p_i, grad_p_j, rho, mu, lambda, U, channelflow_model ,liddriven_model);
    
    % Solve momentum
    u = sor(Nx, Ny, Ac_u, Ae, Aw, An, As, Bc_u, u, n_sor_m);
    v = sor(Nx, Ny, Ac_v, Ae, Aw, An, As, Bc_v, v, n_sor_m);
   
    % Rhie & Chow interpolation
    [u, v, uf, vf] = rhie_and_chow(Nx, Ny, u, v, uf, vf, p, grad_p_i, grad_p_j, ...
                                   Ac_u, Ac_v, dx, dy, channelflow_model);

    % Construct continuity eq
    [Ac_p, Ae, Aw, An, As, Bc, p_prime] = continuity(Nx, Ny, dx, dy, Ac_u, Ac_v, Ac_p, Ae, Aw, An, As, Bc, p_prime, uf, vf, rho, liddriven_model);

    % Solve continuity
    p_prime = sor(Nx, Ny, Ac_p, Ae, Aw, An, As, Bc, p_prime, n_sor_c);

    % % Gradient
    [grad_p_prime_i, grad_p_prime_j] = gradgauss_prime(Nx, Ny, p_prime, dx, dy, grad_p_prime_i, grad_p_prime_j);

    % Pressure correction
    [u, v, uf, vf, p] = pressure_correction(Nx, Ny, u, v, uf, vf, p, grad_p_i, grad_p_j, p_prime, grad_p_prime_i, grad_p_prime_j,...
                                            Ac_u, Ac_v, lambda_p, dx, dy, channelflow_model, rho, U, Ly);
    % Real-time
    if transient == true, tproc_velocity(u, v, X, Y, x, y, Nx, Ny, Lx, channelflow_model); end

    % Residual
    res = residual(Nx, Ny, p_prime, Ac_u, dy);
    if or(mod(i_iteration, residual_gap) == 0, i_iteration == 1), fprintf('%d\t%.4e\n', i_iteration, res); end
    if res < tol, fprintf('%5d\t%.4e\n', i_iteration, res), break; end    
end

%% Post-Processing
if ~transient, postproc_velocity(u, v, X, Y, x, y, Nx, Ny, Lx, channelflow_model, liddriven_model); end
if liddriven_model, postproc_centerline(Nx, Ny, y, u); uy = [u(2:Ny + 1, round(Nx/2)), y']; end
