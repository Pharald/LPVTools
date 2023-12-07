function [Y,T,X] = initial(varargin)
% INITIAL  Pointwise initial condition response for PLFTSS objects
%
% [Y,T,X] = INITIAL(SYS,X0,DOMAIN) plots the undriven response of the 
% PLFTSS SYS starting from the initial condition X0. The initial condition 
% response is plotted for each system at each point in the 
% domain specified in the RGRID object DOMAIN. 
%
% [Y,T,X] = INITIAL(SYS,X0) returns the simulated response data.
% The output T is a NT-by-1 PMAT which descibes the time 
% vector for the simulation results. Y is a NT-by-NY PMAT which describes 
% the trajectories of the NY outputs of SYS, and X is a NT-by-NX PMAT 
% which describes the trajectories of the NX states of SYS, at each point 
% in DOMAIN.
%
% [Y,T,X] = INITIAL(SYS,X0,TFINAL,DOMAIN) simulates the initial condition 
% responses from t=0 to the final time t=TFINAL.
%
% [Y,T,X] = INITIAL(SYS,X0,T,DOMAIN) uses the user-supplied time vector T 
% for the simulations.
%
% INITIAL(SYS1,SYS2,...,X0,T,DOMAIN) plots the responses of multiple systems 
% on a single plot.
%
% See also: initial, step, impulse, lsim.

  
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
    initial(varargin{:});
else
    [Y,T,X] = initial(varargin{:});
end
    


