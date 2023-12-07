function [Y,T,X] =  step(varargin)
% STEP  Pointwise step response for PLFTSS objects
%
% STEP(SYS,DOMAIN) plots the step response for SYS at each point in the 
% domain specified in the RGRID object DOMAIN. 
%
% [Y,T,X] = STEP(SYS) returns the simulated response data.
% The output T is a NT-by-1 PMAT which descibes the time 
% vector for the simulation results. Y is a NT-by-NY PMAT which describes 
% the trajectories of the NY outputs of SYS, and X is a NT-by-NX PMAT 
% which describes the trajectories of the NX states of SYS, at each point 
% in DOMAIN.
%
% [Y,T,X] = STEP(SYS,TFINAL,DOMAIN) simulates the step responses from t=0 
% to the final time t=TFINAL.
%
% [Y,T,X] = STEP(SYS,T,DOMAIN) uses the user-supplied time vector T for the 
% simulations.
%
% STEP(SYS1,SYS2,...,T,DOMAIN) plots the step responses of multiple systems 
% on a single plot.
%
% See also: step, impulse, initial, lsim.

  
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
    step(varargin{:});
else
    [Y,T,X] = step(varargin{:});
end
    


