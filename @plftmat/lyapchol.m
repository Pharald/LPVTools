function varargout = lyapchol(varargin)
% LYAPCHOL   Pointwise square-root solution of CT Lyapunov equation for PLFTMATs
%
% R = LYAPCHOL(A,B,DOMAIN) computes a Cholesky factorization X = R'*R of 
% the solution X to the continuous-time Lyapunov matrix equation:  
%    A*X + X*A' + B*B' = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% R = LYAPCHOL(A,B,E,DOMAIN) computes a Cholesky factorization X = R'*R of
% the solution X to the generalized Lyapunov equation:
%   A*X*E' + E*X*A' + B*B' = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% See also: lyapchol, lyap, dlyapchol, dlyap.

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
[varargout{:}]=lyapchol(varargin{:});