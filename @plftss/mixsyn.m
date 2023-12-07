function [K,CL,GAM,INFO] = mixsyn(varargin)
% MIXSYN Pointwise H-infinity mixed sensitivity synthesis for PLFTSS
%
% [K,CL,GAM,INFO]=mixsyn(G,W1,W2,W3,...,DOMAIN) performs an H-infty mixed
% sensitivity synthesis at each point in the domain specified by the 
% RGRID object DOMAIN. See LTI/MIXSYN for details. 
%
% See also: mixsyn, hinfsyn, ncfsyn, h2syn, loopsyn.

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

[K,CL,GAM,INFO] = mixsyn(varargin{:});