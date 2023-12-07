%% LPVELIMIV - Eliminate parameters which only have a single grid point
%
%  
%% Syntax
%
%    E = lpvelimiv(M)
%
%% Description
%
% |E = lpvelimiv(M)| eliminates all singleton independent variables from 
% the grid-based LPV object |M|, i.e. the i^th independent variable of |M| 
% is removed if |M.IVData{i}| has length one. 


