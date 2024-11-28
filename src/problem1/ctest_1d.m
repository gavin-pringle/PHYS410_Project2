%% 1.4 - 1d Convergence Testing

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
    % Compute exact errors and their rms values for later
    Epsi{l} = psixct{l} - psi{l};
    rms_Epsi{l} = rms(abs(Epsi{l}), 2);
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

% Plot scaled errors for different discretization levels
fig1 = figure;
rho = 4;
hold on 
plot(t{6}, rms_dpsi6, 'LineWidth', 2);
plot(t{7}, rho*rms_dpsi7, 'LineWidth', 2);
plot(t{8}, rho^2*rms_dpsi8, 'LineWidth', 2);
xlabel("Time");
ylabel("l-2 norm of difference between level");
legend('||dΨ^6||', '4 * ||dΨ^7||', '4^2 * ||dΨ^8||', 'Location', 'best');
title({"1d Schrodinger equation convergence test - Exact family"
       "l-2 norm of difference between level l solutions"
       "idtype = 0, vtype = 0, tmax = 0.25, lambda = 0.1, 6 <= l <= 9"});
ax = gca;
ax.FontSize = 12;

% Plot scaled exact errors for different discretization levels
fig2 = figure;
rho = 4;
hold on 
plot(t{6}, rms_Epsi{6}, 'LineWidth', 2);
plot(t{7}, rho*rms_Epsi{7}, 'LineWidth', 2);
plot(t{8}, rho^2*rms_Epsi{8}, 'LineWidth', 2);
plot(t{9}, rho^3*rms_Epsi{9}, 'LineWidth', 2);
xlabel("Time");
ylabel("l-2 norm of exact error");
legend('||EΨ^6||', '4 * ||EΨ^7||', '4^2 * ||EΨ^8||', '4^3 * ||EΨ^9||'...
        , 'Location', 'best');
title({"1d Schrodinger equation convergence test - Exact family"
       "l-2 norm of exact error for each level l"
       "idtype = 0, vtype = 0, tmax = 0.25, lambda = 0.1, 6 <= l <= 9"});
ax = gca;
ax.FontSize = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Test #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Simulation maximum time 
tmax = 0.01;
% Discretization levels
minlevel = 6;
maxlevel = 9;
% Delta t by Delta x ratio
lambda = 0.01;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 1;
idpar = [0.50 0.075 0.0]; 

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

% Plot scaled errors for different discretization levels
fig3 = figure;
rho = 4;
hold on 
plot(t{6}, rms_dpsi6, 'LineWidth', 2);
plot(t{7}, rho*rms_dpsi7, 'LineWidth', 2);
plot(t{8}, rho^2*rms_dpsi8, 'LineWidth', 2);
xlabel("Time");
ylabel("l-2 norm of difference between level");
legend('||dΨ^6||', '4 * ||dΨ^7||', '4^2 * ||dΨ^8||', 'Location', 'best');
title({"1d Schrodinger equation convergence test - Boosted Gaussian"
       "l-2 norm of difference between level l solutions"
       "idtype = 1, vtype = 0, tmax = 0.01, lambda = 0.01, 6 <= l <= 9"});
ax = gca;
ax.FontSize = 12;