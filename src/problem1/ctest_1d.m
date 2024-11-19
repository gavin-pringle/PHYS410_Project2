%% 1.4 - Convergence Testing

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 10;
% Discretization level
level = 6;
% Delta t by Delta x ratio
lambda = 1;

% Derived parameters 
nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
dx = x(2) - x(1);
dt = lambda * dx;
nt = round(tmax / dt) + 1;
t = (0 : nt-1) * dt;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 1;
idpar = zeros(1,3);
idpar(1) = 0.5; % m OR x0
idpar(2) = 1; % delta
idpar(3) = 1; % p 

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 01;
vpar = zeros(1,3);
vpar(1) = 0.8; % xmin 
vpar(2) = 0.9; % xmax 
vpar(3) = -10;  % V_c

% Compute the solution
[x t psi psire psiim psimod prob v] = ...
        sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

surf(x, t, psimod);
xlabel('x'), ylabel('t'), zlabel('|z|')

%plot(prob(:, nx))