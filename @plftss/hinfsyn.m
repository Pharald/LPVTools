function [K,CL,GAM,INFO] = hinfsyn(varargin)
% HINFSYN Pointwise H-infinity controller synthesis for PLFTSS
%
% [K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,...,DOMAIN) performs an H-infty
% design at each point in the domain specified by the RGRID object DOMAIN.
% See LTI/HINFSYN for details.
%
% See also: hinfsyn, mixsyn, ncfsyn, h2syn, loopsyn.

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

[K,CL,GAM,INFO] = hinfsyn(varargin{:});