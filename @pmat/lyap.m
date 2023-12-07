function  X = lyap(A,Q,C,E)
% LYAP   Pointwise solve continuous-time Lyapunov equation for PMATs
%
% X = LYAP(A,Q) solves the Lyapunov matrix equation: A*X + X*A' + Q = 0
% at each point in the combined domain of A and Q.
%
% X = LYAP(A,B,C) solves the Sylvester equation: A*X + X*B + C = 0
% at each point in the combined domain of A, B and C.
%
% X = LYAP(A,Q,[],E) solves the generalized Lyapunov equation:
%   A*X*E' + E*X*A' + Q = 0 
% at each point in the combined domain of A, Q and E. 
%
% See also: lyap, lyapchol, dlyap, dlyapchol.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(2, 4)

if nin==2
    % X = LYAP(A,Q)
    [Aext,Qext] = domunion(A,Q);
    AData = Aext.DataPrivate;
    QData = Qext.DataPrivate;
    
    szA = privatesize(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = lyap( AData(:,:,i), QData(:,:,i));
    end
    X = pmat( X, Aext.DomainPrivate );    
elseif nin==3
    % X = LYAP(A,B,C)  
    B = Q;
    [Aext,Bext,Cext] = domunion(A,B,C);
    AData = Aext.DataPrivate;
    BData = Bext.DataPrivate;
    CData = Cext.DataPrivate;
    
    szA = privatesize(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = lyap( AData(:,:,i), BData(:,:,i), CData(:,:,i));
    end
    X = pmat( X, Aext.DomainPrivate );           
else
    % X = LYAP(A,Q,[],E)    
    if isa(C,'pmat') 
        if isempty(C)
            C = [];
        else
            C = double(C);
        end
    end
    [Aext,Qext,Eext] = domunion(A,Q,E);
    AData = Aext.DataPrivate;
    QData = Qext.DataPrivate;
    EData = Eext.DataPrivate;
    
    szA = privatesize(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = lyap( AData(:,:,i), QData(:,:,i), C, EData(:,:,i));
    end
    X = pmat( X, Aext.DomainPrivate );        
end

