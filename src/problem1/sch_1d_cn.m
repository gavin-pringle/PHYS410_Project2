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
    t = (0 : nt-1) * dt;

    % Initialize solution, and set initial data
    psi = zeros(nt, nx);
    if idtype == 0
        % Exact family 
        psi(1, :) = sin(idpar(1) * pi * x);
    elseif idtype == 1
        % Boosted Gaussian
        psi(1, :) = exp(1i * idpar(3) * x) .* ...
                    exp(-((x - idpar(1)) ./ idpar(2)) .^ 2);
    else
        fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
        return
    end
    % Set first and last values of initial data to zero
    psi(1, 1) = 0; 
    psi(1, nx) = 0; 

    % Initial storage for prob and calculate for initial time
    prob = zeros(nt, nx);
    for j = 1 : nx
        prob(1, j) = trapz(x(1:j), abs(psi(1, 1:j)).^2);
    end 

    % Initialize potential
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

    % Initialize storage for RHS
    f  = zeros(nx,1);

    % Set up tridiagonal system
    dl = 0.5/dx^2 * ones(nx, 1);
    d  = (1i/dt - 1/dx^2 - 0.5*v.') .* ones(nx,1);
    du = dl;
    % Fix up boundary cases
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;
    % Define sparse matrix
    A = spdiags([dl d du], -1:1, nx, nx);

    % Compute solution using CN scheme
    for n = 1 : nt-1
        % Define RHS of linear system
        f(2:nx-1) = u(n, 2:nx-1) .* (1i/dt + 1/dx^2 + 0.5*v(2:nx-1)) ...
                    + (-0.5/dx^2) * (u(n, 1:nx-2) + u(n, 3:nx));
        f(1) = 0.0;
        f(nx) = 0.0;
        % Solve system, thus updating approximation to next time step
        psi(n+1, :) = A \ f;
        % Set first and last values to zero
        psi(n+1, 1) = 0; 
        psi(n+1, nx) = 0; 

        % Calculate prob each time step
        for j = 1 : nx
            prob(n+1, j) = trapz(x(1:j), abs(psi(n+1, 1:j)).^2);
        end 
    end

    % Compute real, imaginary, and modulus of each entry in psi
    psire  = real(psi);
    psiim  = imag(psi);
    psimod = abs(psi); 
end