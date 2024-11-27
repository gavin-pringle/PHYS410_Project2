% sch_2d_adi: Solves 2D Schrödinger equation using ADI scheme.
% 
% Inputs:
%
%   tmax:   Maximum integration time
%   level:  Discretization level
%   lambda: dt/dx
%   idtype: Selects initial data type
%   idpar:  Vector of initial data parameters
%   vtype:  Selects potential type
%   vpar:   Vector of potential parameters
%
% Outputs:
%
%   x:      Vector of x coordinates [nx]
%   y:      Vector of y coordinates [ny]
%   t:      Vector of t coordinates [nt]
%   psi:    Array of computed psi values [nt x nx x ny]
%   psire:  Array of computed psi_re values [nt x nx x ny]
%   psiim:  Array of computed psi_im values [nt x nx x ny]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx x ny]
%   v:      Array of potential values [nx x ny]
function [x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)

    % Define mesh and derived parameters
    nx = 2^level + 1;
    ny = nx;
    x = linspace(0.0, 1.0, nx);
    y = linspace(0.0, 1.0, nx);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    t = (0 : nt-1) * dt;

    % Initialize solution, and set initial data
    psi = zeros(nt, nx, ny);
    if idtype == 0
        % Exact family 
        psi(1, :, :) = sin(idpar(1) * pi * x)' * sin(idpar(2) * pi * y);
    elseif idtype == 1
        % Boosted Gaussian
        psi(1, :) = exp(1i * idpar(3) * x) .* ...
                    exp(-((x - idpar(1)) ./ idpar(2)) .^ 2);
    else
        fprintf('sch_2d_adi: Invalid idtype=%d\n', idtype);
        return
    end
    % Set first and last values of initial data to zero
    psi(1, 1) = 0; 
    psi(1, nx) = 0; 

end