function varargout = lyap(varargin)
% LYAP   Pointwise solve continuous-time Lyapunov equation for PLFTMATs
%
% X = LYAP(A,Q,DOMAIN) solves the Lyapunov matrix equation: 
%   A*X + X*A' + Q = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% X = LYAP(A,B,C,DOMAIN) solves the Sylvester equation: A*X + X*B + C = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% X = LYAP(A,Q,[],E,DOMAIN) solves the generalized Lyapunov equation:
%   A*X*E' + E*X*A' + Q = 0 
% at each point in the domain specified by the RGRID object DOMAIN.
%
% See also: lyap, lyapchol, dlyap, dlyapchol.

DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plftss') || isa(varargin{i},'plftmat')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

varargout = cell(nargout,1);
[varargout{:}]=lyap(varargin{:});