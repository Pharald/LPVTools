function varargout = dlyap(varargin)
% LYAP   Pointwise solve discrete-time Lyapunov equation for PLFTMATs
%
% X = LYAP(A,QDOMAIN) solves the Lyapunov matrix equation: 
%     A*X*A' - X + Q = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% X = LYAP(A,B,C,DOMAIN) solves the Sylvester equation: A*X*B - X + C = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% X = LYAP(A,Q,[],E,DOMAIN) solves the generalized Lyapunov equation:
%   A*X*A' - E*X*E' + Q = 0 
% at each point in the domain specified by the RGRID object DOMAIN.
%
% See also: dlyap, dlyapchol, lyap, lyapchol.

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
[varargout{:}]=dlyap(varargin{:});