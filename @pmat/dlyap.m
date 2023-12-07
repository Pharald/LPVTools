function  X = dlyap(A,Q,C,E)
% LYAP   Pointwise solve discrete-time Lyapunov equation for PMATs
%
% X = LYAP(A,Q) solves the Lyapunov matrix equation: A*X*A' - X + Q = 0
% at each point in the combined domain of A and Q.
%
% X = LYAP(A,B,C) solves the Sylvester equation: A*X*B - X + C = 0
% at each point in the combined domain of A, B and C.
%
% X = LYAP(A,Q,[],E) solves the generalized Lyapunov equation:
%   A*X*A' - E*X*E' + Q = 0 
% at each point in the combined domain of A, Q and E. 
%
% See also: dlyap, dlyapchol, lyap, lyapchol.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(2, 4)

if nin==2
    % X = LYAP(A,Q)
    [Aext,Qext] = domunion(A,Q);
    AData = Aext.DataPrivate;
    QData = Qext.DataPrivate;
    
    szA = size(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = dlyap( AData(:,:,i), QData(:,:,i));
    end
    X = pmat( X, Aext.DomainPrivate );    
elseif nin==3
    % X = LYAP(A,B,C)  
    B = Q;
    [Aext,Bext,Cext] = domunion(A,B,C);
    AData = Aext.DataPrivate;
    BData = Bext.DataPrivate;
    CData = Cext.DataPrivate;
    
    szA = size(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = dlyap( AData(:,:,i), BData(:,:,i), CData(:,:,i));
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
    
    szA = size(Aext);
    X = zeros(szA);    
    for i=1:prod(szA(3:end))
        X(:,:,i) = dlyap( AData(:,:,i), QData(:,:,i), C, EData(:,:,i));
    end
    X = pmat( X, Aext.DomainPrivate );        
end

