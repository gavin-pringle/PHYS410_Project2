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
    [X{l}, Y{l}] = meshgrid(x{l}, y{l});

    % Compute exact solution
    psixct{l} = zeros(nt{l}, nx{l}, ny{l});
    for n = 1 : nt{l}
        psixct_n{l} = exp(-1i*(mx^2 + my^2)*pi^2*t{l}(n)) .* ...
                      sin(mx*pi*X{l}) .* sin(my*pi*Y{l});
        psixct{l}(n, :, :) = reshape(psixct_n{l}, [1, nx{l}, ny{l}]);
    end

    % Compute exact errors and their rms values for later
    Epsi{l} = psixct{l} - psi{l};
    Epsi_2d{l} = reshape(Epsi{l}, nt{l}, nx{l}*ny{l});
    rms_Epsi{l} = rms(abs(Epsi_2d{l}), 2);
end

