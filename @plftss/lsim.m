function [Y,T,X] = lsim(varargin)
% LSIM  Pointwise time-domain response of PLFTSS objects to arbitrary inputs
%
% [Y,TOUT,X] = LSIM(SYS,U,T,DOMAIN) plots the time-domain response of the  
% PLFTSS SYS to the input signal described by U and T at each point in the  
% domain specified in the RGRID object DOMAIN. T is a NT-by-1 column vector 
% which describes the time vector for the simulation. U is a NT-by-NU matrix 
% which describes the trajectory of each of the NU inputs to SYS. 
%
% [Y,TOUT,X] = LSIM(SYS,U,T) returns the simulated response data.
% The output TOUT is a NOUT-by-1 vector which descibes the time vector for 
% the simulation results. Y is a NOUT-by-NY matrix which describes the 
% trajectories of the NY outputs of SYS, and X is a NOUT-by-NX matrix which 
% describes the trajectories of the NX states of SYS, at each point in 
% DOMAIN.
%
% [Y,TOUT,X] = LSIM(SYS,U,T,X0,DOMAIN) specifies the initial state vector X0 at time T(1).
% X0 is set to zero when omitted.
%
% LSIM(SYS1,SYS2,...,U,T,X0,DOMAIN) plots the responses of multiple systems 
% on a single plot.
%
% See also: lsim, step, impulse, initial.

DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end
for i = 1:numel(varargin)
    if isa(varargin{i},'plftss')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

if nargout ==0
    lsim(varargin{:});
else
    [Y,T,X] = lsim(varargin{:});
end
    