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
