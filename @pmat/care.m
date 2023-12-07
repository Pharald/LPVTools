function varargout = care(A,B,Q,varargin)
% CARE  Pointwise solve continuous-time algebraic Riccati equations for PMATs
%
% [X,L,G,REPORT] = CARE(A,B,Q,R,S,E) computes the unique stabilizing
% solution X of the continuous-time algebraic Riccati equation                               -1
%        A'XE + E'XA - (E'XB + S) inv(R)  (B'XE + S') + Q = 0 .
% at each point in the combined domains of A, B, Q, R, S and E.
% X is returned as a PMAT that specifies the solution at each point in the 
% combined domain. Similarly, G is a PMAT that specifies the gain matrix, 
% L is a PMAT that contains the eigenvalues,  and REPORT is a PSTRUCT that 
% contains the diagnosis reports. See numerics/care for additional calling 
% options.
%
% See also: care, dare.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(3, inf)

% TODO: GJB 24Jul12  Broken, needs to be fixed.

% Find flags and R/S/E in the variable length input
% Note: CARE allows 'factor'/'nobalance' to appear in any input slot.
idx = cellfun('isclass',varargin,'char');
FLG = varargin(idx);
V = varargin(~idx);
nV = length(V);

% Set defaults for R/S/E
szB = size(B); %OK since only 1st 2 size dimensions used
n = szB(1);
m = szB(2);
if nV<1 || isempty(V{1})
    R = eye(m);
else
    R = V{1};
end
if nV<2 || isempty(V{2})
    S = zeros(n,m);
else
    S = V{2};
end
if nV<3 || isempty(V{3})
    E = eye(n);
else
    E = V{3};
end

% Find common domain for all data
[Acom,Bcom,Qcom,Rcom,Scom,Ecom] = domunion(A,B,Q,R,S,E);
AData = Acom.DataPrivate;
BData = Bcom.DataPrivate;
QData = Qcom.DataPrivate;
RData = Rcom.DataPrivate;
SData = Scom.DataPrivate;
EData = Ecom.DataPrivate;

% Initialize output dimensions based on call to double/care
szA = privatesize(Acom);

temp = cell(nout,1);
[temp{:}] = care( AData(:,:,1), BData(:,:,1), QData(:,:,1), ...
    RData(:,:,1), SData(:,:,1), EData(:,:,1), FLG{:});
for j=1:nout
    s = size(temp{j});
    varargout{j} = zeros( [s(1) s(2) szA(3:end)] );
end

% Solve Riccati equation at each point in the domain
% TODO PJS 4/7/2011: Revisit. Can this double-loop be vectorized?
for i=1:prod(szA(3:end))
    [temp{:}] = care( AData(:,:,i), BData(:,:,i), QData(:,:,i), ...
        RData(:,:,i), SData(:,:,i), EData(:,:,i), FLG{:});
    for j=1:nout
        varargout{j}(:,:,i) = temp{j};
    end
end

% Convert outputs to PMATs
Domain = Acom.DomainPrivate;
for j=1:nout
    varargout{j} = pmat( varargout{j}, Domain );
end




