function [K,S,E] = lqr(varargin)
% LQR   Pointwise Linear Quadratic regulator design for PLFTMAT and PLFTSS.
%
% [K,S,E] = LQR(SYS,Q,R,N,DOMAIN) performs a LQR design at each point in 
% the  domain specified by the RGRID object DOMAIN where SYS is a PLFTSS.
%
% [K,S,E] = LQR(A,B,Q,R,DOMAIN) performs the pointwise LQR design for
% continuous-time models of the form dx/dt = Ax + Bu.  The state matrices
% A and B are PLFTMATs and the design is performed at each point in the
% domain specified by the RGRID object DOMAIN.
%    
% See also: lqr, dlqr

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

[K,S,E]=lqr(varargin{:});