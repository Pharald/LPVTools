function C = cond(M,P)
% COND   Condition number with respect to inversion for PMAT objects.
%
% COND(M) is the condition number of M at each point in the domain of M. 
%
% COND(M,P) returns the condition number in the P-norm.
%
% See also: cond, rcond, norm.

% TODO: GJB 25Jul12  Should not allow P-norm to change over parm vars.

% Check # of input/output arguments
nin = nargin;
narginchk(1, 2)
if nin==1
    szm = privatesize(M);
    C = zeros([1 1 szm(3:end)]);
    Data = M.DataPrivate;
    for i=1:prod(szm(3:end))
        C(1,1,i) = cond(Data(:,:,i));
    end
    C = pmat(C,M.DomainPrivate);    
elseif strcmp(P,'fro')
    % NOTE PJS 4/2/2011: This function is overloaded to allow the P
    % to be specified as a PMAT except for the case where P='fro'.
    % We cannot handle chars that depend on IVs. Thus P='fro' must
    % be handled separately as char input.  Revisit?
    
    szm = privatesize(M);
    C = zeros([1 1 szm(3:end)]);
    Data = M.DataPrivate;
    for i=1:prod(szm(3:end))
        C(1,1,i) = cond(Data(:,:,i),'fro');
    end
    C = pmat(C,M.DomainPrivate);        
else
    % Define M and P on a common domain 
    [Mext,Pext] = domunion(M,P);
    Data = Mext.DataPrivate;
    PData = Pext.DataPrivate;
    szm = privatesize(Mext);
    
    C = zeros([1 1 szm(3:end)]);
    for i=1:prod(szm(3:end))
        C(1,1,i) = cond(Data(:,:,i),PData(:,:,i));
    end
    C = pmat(C,Mext.DomainPrivate);    
end
