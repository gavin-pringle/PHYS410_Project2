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
    u = zeros(nt, nx);
    if idtype == 0
        % Exact family 
        
    elseif idtype == 1
        % Boosted Gaussian

    else
        fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
        return
    end
 
    % Initialize potential
    if vtype == 0
        % No potential
        
    elseif vtype == 1
        % Rectangular barrier or well 
        
    else
        fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
        return
    end
