%% 1.4 - Convergence Testing

close all;
clear; clc;
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Test #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation maximum time 
tmax = 0.25;
% Discretization levels
minlevel = 6;
maxlevel = 9;
% Delta t by Delta x ratio
lambda = 0.1;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 0;
idpar = [3]; % m = 3

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 0;
vpar = zeros(1,3);

% Perform computation at various levels of discretization, store
% results in cell arrays ...
for l = minlevel : maxlevel
    % Compute the solution
    [x{l} t{l} psi{l} psire{l} psiim{l} psimod{l} prob{l} v{l}] ...
        = sch_1d_cn(tmax, l, lambda, idtype, idpar, vtype, vpar);

    [nt{l}, nx{l}] = size(psi{l});

    % Since idtype == 0, compute exact solution
    psixct{l} = zeros(nt{l}, nx{l});
    for n = 1 : nt{l}
        psixct{l}(n,:) = exp(-1i * idpar(1)^2 * pi^2 * t{l}(n)) ...
                       * sin(idpar(1) * pi * x{l});
    end
end

% Calculating the level-to-level differences, taking every second 
% value of the larger length array
dpsi6 = downsample(psi{7}, 2) - psi{6}; 
dpsi7 = downsample(psi{8}, 2) - psi{7}; 
dpsi8 = downsample(psi{9}, 2) - psi{8}; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Test #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%