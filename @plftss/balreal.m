function [sysb,G,T,Ti] = balreal(varargin)
% BALREAL   Pointwise Gramian-based balancing for PLFTSS objects
%
% [SYSB,G] = BALREAL(SYS,DOMAIN) computes a balanced state-space 
% realization SYSB for the stable portion of the PLFTSS SYS. The balancing 
% is done at each point in the domain specified by the RGRID object DOMAIN.
% G is a PMAT containing the vector of Hankel Singular Values of SYS at 
% each point in the domain.
%
% [SYSB,G,T,Ti] = BALREAL(SYS,OPTIONS,DOMAIN) allows options to be set.
% See HSVDOPTIONS for details on balancing options.  In addition,
% T and Ti are the balancing similarity transformation and its inverse at
% each point in the domain of SYS.
%
% See also: balreal, hsvdOptions, gram, modred.

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

[sysb,G,T,Ti] = balreal(varargin{:});