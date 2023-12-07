function W = gram(varargin)
% GRAM   Pointwise computes Gramians for PLFTSS objects
%
% Wc = GRAM(SYS,'c',DOMAIN) computes the controllability gramian of SYS
% at each point in the domain specified by the RGRID object DOMAIN.
% 
% Wo = GRAM(SYS,'o',DOMAIN) computes its observability gramian.
%
% Rc = GRAM(SYS,'cf',DOMAIN) returns the Cholesky factor of Wc.
%
% Ro = GRAM(SYS,'of',DOMAIN) returns the Cholesky factor of Wo.
%
% See also: gram, balreal, lpvgram.

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

W=gram(varargin{:});