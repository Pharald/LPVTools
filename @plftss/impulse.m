function [Y,T,X] = impulse(varargin)
% IMPULSE  Pointwise impulse response for PLFTSS objects
%
% IMPULSE(SYS,DOMAIN) plots the impulse response for SYS at each point in  
% the domain specified in the RGRID object DOMAIN.  
%
% [Y,T,X] = IMPULSE(SYS,DOMAIN) returns the simulated response data.
% The output T is a NT-by-1 PMAT which descibes the time 
% vector for the simulation results. Y is a NT-by-NY PMAT which describes 
% the trajectories of the NY outputs of SYS, and X is a NT-by-NX PMAT 
% which describes the trajectories of the NX states of SYS, at each point 
% in DOMAIN.
%
% [Y,T,X] = IMPULSE(SYS,TFINAL,DOMAIN) simulates the impulse responses from 
% t=0 to the final time t=TFINAL.
%
% [Y,T,X] = IMPULSE(SYS,T,DOMAIN) uses the user-supplied time vector T for 
% the simulations.
%
% IMPULSE(SYS1,SYS2,...,T,DOMAIN) plots the responses of multiple systems 
% on a single plot.
%
% See also: impulse, initial, step, lsim.

  
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
    impulse(varargin{:});
else
    [Y,T,X] = impulse(varargin{:});
end
    


