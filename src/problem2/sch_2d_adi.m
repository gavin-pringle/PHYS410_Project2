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