%% 2.4 - 2d Convergence Testing

close all;
clear; clc;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Test #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation maximum time 
tmax = 0.05;
% Discretization levels
minlevel = 6;
maxlevel = 7; % CHANGE THIS!!!!
% Delta t by Delta x ratio
lambda = 0.05;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 0;
idpar = [2, 3]; % mx = 2, my = 3

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 0;
vpar = zeros(1,5);

% Perform computation at various levels of discretization, store
% results in cell arrays ...
for l = minlevel : maxlevel
    % Compute the solution
    [x{l} y{l} t{l} psi{l} psire{l} psiim{l} psimod{l} v{l}] = ...
        sch_2d_adi(tmax, l, lambda, idtype, idpar, vtype, vpar);

    [nt{l}, nx{l}, ny{l}] = size(psi{l});

    % Since idtype == 0, compute exact solution
    psixct{l} = zeros(nt{l}, nx{l}, nx{l});
    
end

