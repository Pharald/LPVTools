function [L,Gamma,Info] = lpvestsyn(P,varargin)
% LPVESTSYN    Synthesize a parameter-varying estimator for PSS systems.
% 
% [L,GAM,INFO] = LPVESTSYN(P,'L2') computes a parameter-varying state  
% estimator for the parameter-varying system P. L takes in all the outputs 
% of P and outputs an estimate of the states of P. L is the constant   
% estimation matrix for the plant P, which minimizes the induced L2 norm 
% of the error  in the state-estimate. GAM is the minimum L2 norm achived 
% by L. INFO is a struct with additional data. 
%
% [L,GAM,INFO] = LPVESTSYN(P,'LQG') computes the constant estimation matrix 
% for the plant P, which minimizes the stochastic LPV bound on the state
% estimation error. The stochastic LPV bound on the state estimation error  
% is defined as the expected value of the average instantaneous power of 
% error signal, assuming system inputs are zero mean, white-noise processes  
% with unit intensity.
%
% [L,GAM,INFO] = LPVESTSYN(P,Yb,ALG) performs a rate-bounded synthesis. 
% Yb is a BASIS object specifying the basis functions to be used in the 
% synthesis. ALG can be either 'L2' or 'LQG'. A call without the ALG
% argument is equivalent to [L,GAM,INFO] = LPVESTSYN(P,Yb,'L2')
%
% See also: lpvsynOptions, lpvsfsyn, lpvsyn, lpvncfyn, lpvmixsyn, lpvloopshape.


% XXX - OPTIONS are not currently being used. 
% [L,GAM,INFO] = LPVESTSYN(P,NMEAS,...,opt) permits the user to pass in 
% See also: lpvsynOptions.

% Input Parsing:
nin = nargin;
narginchk(1, 3)

% Convert to dual state-feedback problem.
[A,B,C,D] = ssdata(P);
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

Pdual = ss(A',[eye(nx) C'],B',[zeros(nu,nx) D']);

% Solve dual state-feedback problem.
if nin==2
    [F,Gamma,Info] = lpvsfsyn(Pdual,ny); 
else 
    [F,Gamma,Info] = lpvsfsyn(Pdual,ny,varargin{:}); 
end

% Convert back to estimation problem
L=F';