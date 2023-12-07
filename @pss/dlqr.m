function varargout = dlqr(SYS,Q,R,N)
% DLQR   Pointwise Linear Quadratic regulator design for discrete-time PSS.
%
% [K,S,E] = LQR(A,B,Q,R,N) performs a pointwise LQR design for
% discrete-time models of the form x[n+1] = Ax[n] + Bu[n].  The quadratic 
% cost function is specified by the matrices Q, R, and N.  All matrices
% can be PMATs and the design is performed at each point in the combined
% domains of (A,B,Q,R,N).   K is the optimal LQ gain matrix at each point 
% in the combined domain. S is the solution of the algebraic Riccati 
% equation and E contains the closed-loop eigenvalues both returned at 
% each point in the combined domain.
%
% See also: dlqr, lqr.

% Input error checking
narginchk(3, 4)
[A,B,C,D,Ts]=ssdata(SYS);
szB = size(B);
if nargin==3
    N = zeros(szB(1:2));
end

varargout = cell(nargout,1);
if Ts==0
    [varargout{:}]=lqr(A,B,Q,R,N);
else
    [varargout{:}]=dlqr(A,B,Q,R,N);     
end

