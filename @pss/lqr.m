function varargout = lqr(SYS,Q,R,N)
% LQR   Pointwise Linear Quadratic regulator design for PMAT and PSS objects.
%
% [K,S,E] = LQR(SYS,Q,R,N) performs a LQR design at each point in the 
% domain of the PSS object SYS. The quadratic cost function is specified
% by the matrices Q, R, and N.  If Q, R, and N are PMATs then the LQR
% design is done at each point in the combined domains of SYS, Q, R, and
% N. K is the optimal LQ gain matrix at each point in the combined domain.
% S is the solution of the algebraic Riccati equation and E contains the 
% closed-loop eigenvalues both returned at each point in the combined 
% domain.
%
% [K,S,E] = LQR(A,B,Q,R,N) performs the pointwise LQR design for
% continuous-time models of the form dx/dt = Ax + Bu.  The state matrices
% A and B can be PMATs.
%    
% See also: lqr, dlqr.

% Input error checking
narginchk(3, 4)
[A,B,C,D,Ts]=ssdata(SYS);
szB = size(B);
if nargin==3
    N = zeros(szB(1:2));
end

% Array dimensions not allowed - LTI function will error out.
nad = numel(size(SYS))-2;
if nad>0
    error('LQR does not allow PSS with array dimensions.');
end

varargout = cell(nargout,1);
if Ts==0
    [varargout{:}]=lqr(A,B,Q,R,N);
else
    [varargout{:}]=dlqr(A,B,Q,R,N);     
end

