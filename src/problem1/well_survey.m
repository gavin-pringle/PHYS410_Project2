%% 1.5.2 - Well Survey

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 0.10;
% Discretization level
level = 9;
% Delta t by Delta x ratio
lambda = 0.01;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 1;
idpar = [0.40, 0.075, 0.0];

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
vtype = 1;
xmin = 0.6;
xmax = 0.8;
lnV0 = linspace(2, 10, 10); % CHANGE THIS

% Survey range  
x1 = 0.6;
x2 = 0.8;

for idx = 1 : length(lnV0)
    % Get this iteration's vpar
    vpar = [xmin, xmax, -exp(lnV0(idx))];

    % Compute the solution
    [x{idx} t{idx} psi{idx} psire{idx} psiim{idx} psimod{idx} prob{idx} v{idx}] ...
       = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

    % Compute temporal average of probability matrix 
    P_bar{idx} = mean(prob{idx});
    % Normalize the temporal average
    P_bar{idx} = P_bar{idx} / P_bar{idx}(end);

    % Compute indices of x1 and x2 in the array x
    % This only needs to be done once 
    if idx == 1
        [~, x1_loc] = min(abs(x{idx} - x1));
        [~, x2_loc] = min(abs(x{idx} - x2));
    end

    % Compute excess fractional probability and its logarithm
    Fe_bar{idx} = (P_bar{idx}(x1_loc) - P_bar{idx}(x2_loc)) / ...
                    (x{idx}(x1_loc) - x{idx}(x2_loc));
    lnFe_bar{idx} = log(Fe_bar{idx});
end

plot(lnV0, cell2mat(lnFe_bar))