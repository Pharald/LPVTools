function [K,S,E] = lqr(A,B,Q,R,N)
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
nin = nargin;
nout = nargout;
narginchk(3, 5)
nargoutchk(0, 3)

% If first argument is an ss then lift and call PSS/LQR
if isa(A,'ss');
    if nin==3
        if nout==1
            K=lqr( pss(A),B,Q);
        elseif nout==2
            [K,S]=lqr( pss(A),B,Q);
        else
            [K,S,E]=lqr( pss(A),B,Q);
        end
    elseif nin==4
        if nout==1
            K=lqr( pss(A),B,Q,R);
        elseif nout==2
            [K,S]=lqr( pss(A),B,Q,R);
        else
            [K,S,E]=lqr( pss(A),B,Q,R);
        end
    end
    return;
end

% Define inputs on a common domain
szB = size(B); %OK sinze SIZE usings only 1st 2 entries
if nin==3
    % PMAT syntax requires four inputs: (A,B,Q,R)
    error(nargchk(4, 5, nin, 'struct'));
elseif nin==4
    N = zeros(szB(1:2));
end
[Aext,Bext,Qext,Rext,Next] = domunion(A,B,Q,R,N);
AData = Aext.DataPrivate;
BData = Bext.DataPrivate;
QData = Qext.DataPrivate;
RData = Rext.DataPrivate;
NData = Next.DataPrivate;

% Array dimensions not allowed - LTI function will error out.
if (numel(size(Aext.Data))-2-Aext.Domain.NumIV)>0
    error('LQR does not allow PMAT with array dimensions.');
end

% Call LQR at each point in the combined domain
szA = privatesize(Aext);
K = zeros([szB(2) szB(1) szA(3:end)]);
if nout==1
    for i=1:prod(szA(3:end))
        K(:,:,i) = lqr( AData(:,:,i),BData(:,:,i), QData(:,:,i), ...
            RData(:,:,i), NData(:,:,i) );
    end
elseif nout==2
    S = zeros([szA(1) szA(1) szA(3:end)]);
    for i=1:prod(szA(3:end))
        [K(:,:,i),S(:,:,i)] = lqr( AData(:,:,i), BData(:,:,i), ...
            QData(:,:,i), RData(:,:,i), NData(:,:,i) );
    end
    S = pmat( S, Aext.DomainPrivate );
else
    S = zeros([szA(1) szA(1) szA(3:end)]);
    E = zeros([szA(1) 1 szA(3:end)]);
    for i=1:prod(szA(3:end))
        [K(:,:,i),S(:,:,i),E(:,:,i)] = lqr( AData(:,:,i), ...
            BData(:,:,i), QData(:,:,i), RData(:,:,i), NData(:,:,i) );
    end
    S = pmat( S, Aext.DomainPrivate );
    E = pmat( E, Aext.DomainPrivate );
end
K = pmat( K, Aext.DomainPrivate );
