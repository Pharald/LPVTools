function [K,CL,GAM,INFO] = h2syn(varargin)
% H2SYN Pointwise H2 controller synthesis for PLFTSS
%
% [K,CL,GAM,INFO] = h2syn(P,NMEAS,NCON,DOMAIN) performs an H2
% design at each point in the domain specified by the RGRID object DOMAIN.
% See LTI/H2SYN for details.
%
% See also: h2syn, hinfsyn, mixsyn, ncfsyn, loopsyn.

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

[K,CL,GAM,INFO] = h2syn(varargin{:});