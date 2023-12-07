function [K,S,E] = dlqr(A,B,Q,R,N)
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
nin = nargin;
nout = nargout;
narginchk(4, 5)

% Define inputs on a common domain
szB = size(B); % OK since only 1st 2 size dimensions are used
if nin==4
    N = zeros(szB(1:2));
end
[Aext,Bext,Qext,Rext,Next] = domunion(A,B,Q,R,N);
AData = Aext.DataPrivate;
BData = Bext.DataPrivate;
QData = Qext.DataPrivate;
RData = Rext.DataPrivate;
NData = Next.DataPrivate;

% Call LQR at each point in the combined domain
szA = privatesize(Aext);
K = zeros([szB(2) szB(1) szA(3:end)]);
if nout==1
    for i=1:prod(szA(3:end))
        K(:,:,i) = dlqr( AData(:,:,i), BData(:,:,i), QData(:,:,i), ...
            RData(:,:,i), NData(:,:,i) );
    end
elseif nout==2
    S = zeros([szA(1) szA(1) szA(3:end)]);
    for i=1:prod(szA(3:end))
        [K(:,:,i),S(:,:,i)] = dlqr( AData(:,:,i), BData(:,:,i), ...
            QData(:,:,i), RData(:,:,i), NData(:,:,i) );
    end
    S = pmat( S, Aext.DomainPrivate );
else
    S = zeros([szA(1) szA(1) szA(3:end)]);
    E = zeros([szA(1) 1 szA(3:end)]);
    for i=1:prod(szA(3:end))
        [K(:,:,i),S(:,:,i),E(:,:,i)] = dlqr( AData(:,:,i), ...
            BData(:,:,i), QData(:,:,i), RData(:,:,i), NData(:,:,i) );
    end
    S = pmat( S, Aext.DomainPrivate );
    E = pmat( E, Aext.DomainPrivate );
end
K = pmat( K, Aext.DomainPrivate );    
