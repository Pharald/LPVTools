function varargout = dlqr(varargin)
%DLQR  Linear-quadratic regulator design for discrete-time PLFT systems.
%
%   [K,S,E] = DLQR(A,B,Q,R,N,DOMAIN) calculates the optimal gain matrix K 
%   at each point in the domain specified by the RGRID object DOMAIN.
%   This function performs a pointwise LQR design for discrete-time models 
%   of the form x[n+1] = Ax[n] + Bu[n]. All matrices can be PLFTMATs. 
%   The quadratic cost function is specified by the matrices Q, R, and N. 
%   K is the optimal LQ gain matrix at each point domain. S is the solution 
%   of the algebraic Riccati equation and E contains the closed-loop 
%   eigenvalues both returned at each point in the combined domain.
%
% See also: dlqr, lqr

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
[varargout{:}]=dqlr(varargin{:});