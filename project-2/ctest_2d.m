%% 2.4 - 2d Convergence Testing

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 0.05;
% Discretization levels
minlevel = 6;
maxlevel = 9;
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

% Perform computation at various levels of discretization, store
% results in cell arrays ...
for l = minlevel : maxlevel
    % Compute the solution
    [x{l} y{l} t{l} psi{l} psire{l} psiim{l} psimod{l} v{l}] = ...
        sch_2d_adi(tmax, l, lambda, idtype, idpar, vtype, vpar)

    [nt{l}, nx{l}, ny{l}] = size(psi{l});

    % Define meshgrid for computing exact psi(x,y,t)
    [X{l}, Y{l}] = meshgrid(x{l}.', y{l}.');

    % Compute exact solution
    psixct{l} = zeros(nt{l}, nx{l}, ny{l});
    for n = 1 : nt{l}
        psixct_n{l} = exp(-1i*(mx^2 + my^2)*pi^2*t{l}(n)) * ...
                      sin(mx*pi*X{l}) .* sin(my*pi*Y{l});
        % Enforce boundary conditions
        psixct_n{l}(1, :)     = 0.0;
        psixct_n{l}(:, 1)     = 0.0;
        psixct_n{l}(nx{l}, :) = 0.0;
        psixct_n{l}(:, ny{l}) = 0.0;
        psixct{l}(n, :, :) = reshape(psixct_n{l}, [1, nx{l}, ny{l}]);
    end

    % Compute exact errors and their rms values for later
    Epsi{l} = psixct{l} - psi{l};
    Epsi_2d{l} = reshape(Epsi{l}, nt{l}, nx{l}*ny{l});
    rms_Epsi{l} = rms(abs(Epsi_2d{l}), 2);

    % Downsample each psi for differencing 
    psi_ds{l} = psi{l}(1:2:end, 1:2:end, 1:2:end);

    % Flatten each 3d array into a 2d array
    psi_2d{l} = reshape(psi{l}, nt{l}, nx{l}*ny{l});
    psi_ds_2d{l} = reshape(psi_ds{l}, (nt{l}-1)/2+1, ...
                          ((nx{l}-1)/2+1)*((ny{l}-1)/2+1));
end

% Calculating the level-to-level differences, taking every second 
% value of the larger length array
dpsi6 = psi_ds_2d{7} - psi_2d{6}; 
dpsi7 = psi_ds_2d{8} - psi_2d{7}; 
dpsi8 = psi_ds_2d{9} - psi_2d{8}; 

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
title({"2d Schrodinger equation convergence test - Exact family"
       "l-2 norm of difference between level l solutions"
       "idtype = 0, vtype = 0, tmax = 0.05, lambda = 0.05, 6 <= l <= 9"});
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
title({"2d Schrodinger equation convergence test - Exact family"
       "l-2 norm of exact error for each level l"
       "idtype = 0, vtype = 0, tmax = 0.05, lambda = 0.05, 6 <= l <= 9"});
ax = gca;
ax.FontSize = 12;