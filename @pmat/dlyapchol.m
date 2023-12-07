function  R = dlyapchol(A,B,E)
% DLYAPCHOL   Pointwise square-root solution of DT Lyapunov equation for PMATs
%
% R = DLYAPCHOL(A,B) computes a Cholesky factorization X = R'*R of 
% the solution X to the discrete-time Lyapunov matrix equation:  
%    A*X*A'- X + B*B' = 0
% at each point in the combined domain of A and B.  
%
% R = DLYAPCHOL(A,B,E) computes a Cholesky factorization X = R'*R of
% the solution X to the generalized Lyapunov equation:
%   A*X*A' - E*X*E' + B*B' = 0
% at each point in the combined domain of A, B and E.  
%
% See also: dlyapchol, dlyap, lyapchol, lyap.

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
        R(:,:,i) = dlyapchol( AData(:,:,i), BData(:,:,i));
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
        R(:,:,i) = dlyapchol( AData(:,:,i), BData(:,:,i), EData(:,:,i));
    end
    R = pmat( R, Aext.DomainPrivate );   
end
