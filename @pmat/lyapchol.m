function  R = lyapchol(A,B,E)
% LYAPCHOL   Pointwise square-root solution of CT Lyapunov equation for PMATs
%
% R = LYAPCHOL(A,B) computes a Cholesky factorization X = R'*R of 
% the solution X to the continuous-time Lyapunov matrix equation:  
%    A*X + X*A' + B*B' = 0
% at each point in the combined domain of A and B.  
%
% R = LYAPCHOL(A,B,E) computes a Cholesky factorization X = R'*R of
% the solution X to the generalized Lyapunov equation:
%   A*X*E' + E*X*A' + B*B' = 0
% at each point in the combined domain of A, B and E.  
%
% See also: lyapchol, lyap, dlyapchol, dlyap.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(2, 3)

if nin==2
    % R = LYAPCHOL(A,B) 
    [Aext,Bext] = domunion(A,B);
    AData = Aext.DataPrivate;
    BData = Bext.DataPrivate;
    
    szA = privatesize(Aext);
    R = zeros(szA);    
    for i=1:prod(szA(3:end))
        R(:,:,i) = lyapchol( AData(:,:,i), BData(:,:,i));
    end
    R = pmat( R, Aext.DomainPrivate );    
else
    % R = LYAPCHOL(A,B,E)
    [Aext,Bext,Eext] = domunion(A,B,E);
    AData = Aext.DataPrivate;
    BData = Bext.DataPrivate;
    EData = Eext.DataPrivate;
    
    szA = privatesize(Aext);
    R = zeros(szA);    
    for i=1:prod(szA(3:end))
        R(:,:,i) = lyapchol( AData(:,:,i), BData(:,:,i), EData(:,:,i));
    end
    R = pmat( R, Aext.DomainPrivate );   
end
