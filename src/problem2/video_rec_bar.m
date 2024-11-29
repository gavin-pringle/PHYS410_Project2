%% 2.5 - 2d video of scattering off rectangular barrier

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 0.1;
% Discretization level
level = 7;
% Delta t by Delta x ratio
lambda = 0.05;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 1;
%x0      = idpar(1);      y0 = idpar(2);    
%delta_x = idpar(3); delta_y = idpar(4); 
%p_x     = idpar(5);     p_y = idpar(6);   
idpar = [0.6, 0.2, 0.075, 0.075, 0.0, 30];

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 1;
%x_min = vpar(1);   x_max = vpar(2);    
%y_min = vpar(3);   y_max = vpar(4); 
%Vc    = vpar(5); 
vpar = [0.1, 0.7, 0.7, 0.8, 1e10];

% Compute solution 
[x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);




