%% 2.5 - 2d video of scattering off rectangular barrier

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 0.05;
% Discretization level
level = 8;
% Delta t by Delta x ratio
lambda = 0.05;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 0;
idpar = [2, 3];
mx = idpar(1); my = idpar(2);

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 0;
vpar = zeros(1,5);