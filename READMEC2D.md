%% Compressible 2D viscous
clear 
close all
clc

%% discretization
% Discretization in x direction 
Lx = 0.01;
Nx = 128;
hx = Lx / Nx;
alphax = 0.5;
x = linspace(0, Lx, Nx)';
% Discretization in y direction 
Ly = 0.01;
Ny = 128;
hy = Ly / Ny;
alphay = 0.5;
y = linspace(0, Ly, Ny)';
% Temporal discretization 
T = 100;

%% Initial conditions
p0 = 101325;
rho = ones(Nx, Ny);
p = p0 * ones(Nx, Ny);

gamma = 1.4;
mu = 1.58 * 1e-5; % mixture viscosity and scalar diffusivity 

% Conserved variables
U1 = rho;
U2 = rho .* u;
U3 = rho .* v;
U4 = (1 / (gamma - 1)) * p + 0.5 * rho .* (u.^2 + v.^2);
U5 = rho .* z;

%% Temporal loop
tempo = zeros(T, 1);
tt = 0;
tn = 1;

for t = 2:T
    disp(t)

    a = sqrt(gamma * p ./ rho);

    % Calculate maximum eigenvalues
    lambda_max_x = max([max(abs(u - a)), max(abs(u)), max(abs(u + a))]);
    lambda_max_y = max([max(abs(v - a)), max(abs(v)), max(abs(v + a))]);

    % Calculate time step size
    taux = alphax * hx / lambda_max_x;
    tauy = alphay * hy / lambda_max_y;
    tau = min([taux, tauy]);

    tempo(t) = tempo(t - 1) + tau;

    %% Riemann problem solution at interfaces
    % Initialize variables for Riemann problem solutions
    rho_PRx = zeros(Nx, Ny);
    u_PRx = zeros(Nx, Ny);
    p_PRx = zeros(Nx, Ny);
    v_PRx = zeros(Nx, Ny);
    z_PRx = zeros(Nx, Ny);
    rho_PRy = zeros(Nx, Ny);
    v_PRy = zeros(Nx, Ny);
    p_PRy = zeros(Nx, Ny);
    u_PRy = zeros(Nx, Ny);
    z_PRy = zeros(Nx, Ny);

    % Calculate Riemann problem solutions at interfaces
    rho_PRx(2:Nx-1, 2:Ny-1) = PR(p(2:Nx-1, 2:Ny-1), u(2:Nx-1, 2:Ny-1), rho(2:Nx-1, 2:Ny-1), p(3:Nx, 2:Ny-1), u(3:Nx, 2:Ny-1), rho(3:Nx, 2:Ny-1), gamma, gamma, v(2:Nx-1, 2:Ny-1), z(2:Nx-1, 2:Ny-1), v(3:Nx, 2:Ny-1), z(3:Nx, 2:Ny-1), 1e-8);
    rho_PRy(2:Nx-1, 2:Ny-1) = PR(p(2:Nx-1, 2
