function [L,Gamma,Info] = lpvFCsyn(P,nmeas,varargin)
% LPVESTSYN    Parameter dependent full control controller synthesis.
% 
% [F,GAM,INFO] = LPVFCSYN(P,NMEAS) computes a parameter-varying
% full-control gain matrix for the parameter-varying system P. NMEAS 
% specifies the number of available measurements being output from P. F is 
% the full control gain matrix for the plant P, which minimizes the L2 norm. 
% GAM is the minimum L2 norm achived by F.
%
% [F,GAM,INFO] = LPVFCSYN(P,NMEAS,Yb) performs a rate-bounded synthesis.
% Yb is a BASIS object specifying the basis functions to be used in the
% synthesis.
%
% See also: lpvsfsyn, lpvsyn, lpvncfsyn.



% XXX - OPTIONS are not currently being used. 
% [L,GAM,INFO] = LPVESTSYN(P,NMEAS,...,opt) permits the user to pass in 
% See also: lpvsynOptions.

% Input Parsing:
nin = nargin;
narginchk(2, 4)


% Convert to dual state-feedback problem.
[A,B,C,D] = ssdata(P);
Pdual = ss(A',C',B',D');

% Solve dual state-feedback problem.
if nin==2
    [F,Gamma,Info] = lpvsfsyn(Pdual,nmeas); 
else 
    [F,Gamma,Info] = lpvsfsyn(Pdual,nmeas,varargin{:}); 
end

% Convert back to estimation problem
L=F';