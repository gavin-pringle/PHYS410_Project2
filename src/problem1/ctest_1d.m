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

    % Compute t for each array
    dx{l} = x{l}(2) - x{l}(1);
    dt{l} = lambda * dx{l};
    nt{l} = round(tmax / dt{l}) + 1;
    t{l} = (0 : nt{l}-1) * dt{l};
end

% Calculating the level-to-level differences, taking every second 
% value of the larger length array
dpsi6 = downsample(downsample(psi{7}, 2).', 2).' - psi{6}; 
dpsi7 = downsample(downsample(psi{8}, 2).', 2).' - psi{7}; 
dpsi8 = downsample(downsample(psi{9}, 2).', 2).' - psi{8}; 

% Compute l-2 norm of each dpsi, resulting in functions of t
rms_dpsi6 = rms(abs(dpsi6), 2);
rms_dpsi7 = rms(abs(dpsi7), 2);
rms_dpsi8 = rms(abs(dpsi8), 2);

% Plot scaled errors for rho = 4
fig2 = figure;
rho = 4;
hold on 
plot(t{6}, rms_dpsi6, 'LineWidth', 2);
plot(t{7}, rho*rms_dpsi7, 'LineWidth', 2);
plot(t{8}, rho^2*rms_dpsi8, 'LineWidth', 2);
xlabel("Time");
ylabel("Difference between level");
legend('||dΨ^6||', '4 * ||dΨ^7||', '4^2 * ||dΨ^8||', 'Location', 'best');
title("Convergence Test: rho = 4");
ax = gca;
ax.FontSize = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Test #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%