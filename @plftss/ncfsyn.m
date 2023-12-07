function [K,CL,GAM,INFO] = ncfsyn(varargin)
% NCFSYN Pointwise H-infinity normalized coprime factor synthesisfor a PLFTSS.
%
% [K,CL,GAM,INFO]=ncfsyn(G,W1,W2,..,DOMAIN) performs an H-infinity normalized 
% coprime factor controller synthesis at each point in the domain 
% specified by the RGRID object DOMAIN. See LTI/NCFSYN for details.
%
% See also: ncfsyn, hinfsyn, mixsyn, h2syn, loopsyn.

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

[K,CL,GAM,INFO] = ncfsyn(varargin{:});