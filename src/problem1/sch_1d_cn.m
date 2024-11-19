% sch_1d_cn: Solves 1D Schrödinger equation using O(dt^2,dx^2) 
% Crank-Nicolson implicit scheme.
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
%   t:      Vector of t coordinates [nt]
%   psi:    Array of computed psi values [nt x nx]
%   psire:  Array of computed psi_re values [nt x nx]
%   psiim:  Array of computed psi_im values [nt x nx]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx]
%   prob:   Array of computed running integral values [nt x nx]
%   v:      Array of potential values [nx]
function [x t psi psire psiim psimod prob v] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)

    % Define mesh and derived parameters
    nx = 2^level + 1;
    x = linspace(0.0, 1.0, nx);
    dx = x(2) - x(1);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    t = [0 : nt-1] * dt;

    % Initialize solution, and set initial data
    psi = zeros(nt, nx);
    if idtype == 0
        % Exact family 
        psi(1, :) = sin(idpar(1) * pi * x);
    elseif idtype == 1
        % Boosted Gaussian
        psi(1, :) = exp(i * idpar(3) * x) .* ...
                    exp(-((x - idpar(1)) ./ idpar(2)) .^ 2);
    else
        fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
        return
    end
 
    % Initialize potential, and 
    v = zeros(1,nx);
    if vtype == 0
        % No potential - leave unchanged
    elseif vtype == 1
        % Rectangular barrier or well 
        v(x > vpar(1) & x < vpar(2)) = vpar(3);
    else
        fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
        return
    end

    % Initialize storage for sparse matrix and RHS
    dl = zeros(nx,1);
    d  = zeros(nx,1);
    du = zeros(nx,1);
    f  = zeros(nx,1);

    % Set up tridiagonal system
    dl = 0.5 / dx^2 * ones(nx, 1);
    d  = (i / dt - 1 / dx^2 - 0.5 * v.') .* ones(nx,1);
    du = dl;
    % Fix up boundary cases
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;
    % Define sparse matrix
    A = spdiags([dl d du], -1:1, nx, nx);

    % Compute solution using CN scheme CHANGE THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for n = 1 : nt-1
        % Define RHS of linear system ...
        f(2:nx-1) = u(n, 2:nx-1) / dt + 0.5 * (u(n, 1:nx-2) - 2 * u(n, 2:nx-1) + u(n, 3:nx)) / dx^2;
        f(1) = 0.0;
        f(nx) = 0.0;
        % Solve system, thus updating approximation to next time 
        % step ...
        u(n+1, :) = A \ f;
  
        if trace && ~mod(n,trace)
           fprintf('diff_1d_cn: Step %d of %d\n', n, nt);
        end
     end
end