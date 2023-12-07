function  varargout = dlyapchol(varargin)
% DLYAPCHOL  Pointwise square-root solution of DT Lyapunov equation for PLFTMATs
%
% R = DLYAPCHOL(A,B,DOMAIN) computes a Cholesky factorization X = R'*R of 
% the solution X to the discrete-time Lyapunov matrix equation:  
%    A*X*A'- X + B*B' = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% R = DLYAPCHOL(A,B,E,DOMAIN) computes a Cholesky factorization X = R'*R of
% the solution X to the generalized Lyapunov equation:
%   A*X*A' - E*X*E' + B*B' = 0
% at each point in the domain specified by the RGRID object DOMAIN.
%
% See also: dlyapchol, dlyap, lyapchol, lyap.

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
[varargout{:}]=dlyapchol(varargin{:});