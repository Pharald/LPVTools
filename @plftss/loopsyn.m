function [K,CL,GAM,INFO] = loopsyn(varargin)
% LOOPSYN Pointwise loopsyn controller synthesis for PLFTSS
%
% [K,CL,GAM,INFO]=loopsyn(G,Gd,..,DOMAIN) performs a loopsyn controller
% synthesis at each point in the domain specified by the RGRID object
% DOMAIN. See LTI/LOOPSYN for details.
%
% See also: loopsyn, hinfsyn, mixsyn, ncfsyn, h2syn.

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

[K,CL,GAM,INFO] = loopsyn(varargin{:});