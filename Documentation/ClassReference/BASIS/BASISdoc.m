%% BASIS - Basis function object
%
%  
%% Syntax
%
%   b = basis(F,NAME1,PARTIAL1,NAME2,PARTIAL2,...)
%   b = basis(F,DERIVATIVE)
%
%% Description
% 
% |basis| functions are needed for rate-bounded LPV analysis and synthesis.
% These are functions of the independent variables present in the system
% being analyzed.  |basis| functions are specified as scalar |pmat|.  All
% partial derivatives must also be provided by the user.
% 
% |b = basis(F,NAME1,PARTIAL1,NAME2,PARTIAL2,...)| creates a |basis| function
% object with the function |F|, which is a scalar |pmat|. If there are N 
% independent variables in |F|, then 2*N additional arguments are specified, 
% which are pairs of (1) an independent variable name (as |char|), 
% and (2) the corresponding partial derivative (as |pmat|) of |F| with  
% respect to that independent variable.
% 
% |basis(F,DERIVATIVE)| is the same as |basis(F,NAME1,DERIVATIVE)| if |F| has
% only one independent variable.
% 
%% Example

theta = pgrid('theta',0:linspace(0,2*pi,20));
psi = pgrid('psi',linspace(0,2*pi,10));
F = cos(theta)*sin(2*psi);
pTheta = -sin(theta)*sin(2*psi);
pPsi = 2*cos(theta)*cos(2*psi);
B = basis(F,'theta',pTheta,'psi',pPsi)
