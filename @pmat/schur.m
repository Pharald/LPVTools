function [u,t] = schur(mat,varargin)
% SCHUR   Schur decomposition of a PMAT object.
%
% [U,T] = SCHUR(X) computes the Schur decomposition of X at each 
% point in the domain of X.
%
% T = SCHUR(X) returns just the Schur matrix T.
%
% If X is real then SCHUR(X,'real') and SCHUR(X,'complex') return 
% the real and complex Schur form as documented in the SCHUR help 
% for double matrices.
%
% See also: schur, qz.


% TODO PJS 4/2/2011: Implement functions listed in the "See also".

% Input / output error checking
nin = nargin;
nout = nargout;
narginchk(1, 2)

szm = privatesize(mat);
Data = mat.DataPrivate;
if nout==1
    % T = SCHUR(M) or  T = SCHUR(M,'real') or ...,'complex')
    u = zeros(szm);
    for i=1:prod(szm(3:end))
        u(:,:,i) = schur(Data(:,:,i), varargin{:});
    end
    u = pmat(u,mat.DomainPrivate);
else
    % [U,T] = SCHUR(M) or [U,T] = SCHUR(M,'real') or ...,'complex')
    u = zeros(szm);
    t = zeros(szm);
    for i=1:prod(szm(3:end))
        [u(:,:,i),t(:,:,i)] = schur(Data(:,:,i), varargin{:});
    end
    u = pmat(u,mat.DomainPrivate);
    t = pmat(t,mat.DomainPrivate);
end
