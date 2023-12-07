function varargout = eig(varargin)
% EIG   Eigenvalues and eigenvectors for a PMAT object.
%
% E = EIG(X) returns a PMAT vector containing the eigenvalues of a square
% matrix X at each point in the domain of X.
%
% [V,D] = EIG(X) returns a diagonal PMAT D of eigenvalues and a
% full PMAT V whose columns are the corresponding eigenvectors so
% that X*V = V*D at each point in the domain of X.
%
% E = EIG(A,B) and [V,D]=EIG(A,B) compute the generalized eigenvalues
% and eigenvectors at each point in the combined domains of A and B.
%
% The function also supports the following syntax as documented in the
% EIG help:
%    [V,D] = EIG(X,'nobalance')
%    EIG(A,B,'chol')
%    EIG(A,B,'qz')
%
% See also: eig.

% TODO PJS 4/2/2011: Implement functions listed in the 
% standard "See also: eig, condeig, eigs, ordeig".

% Input / output error checking
nin = nargin;
nout = nargout;
narginchk(1, 3)
nargoutchk(0, 3)

% Check for eig or geneig problem
geflag = false;
A = pmat( varargin{1} );
if nin>=2 && ~ischar(varargin{2})
    geflag = true;
    B = pmat( varargin{2} );
end

if geflag
    % Generalized Eigenval/Eigenvec Problem
    [A,B] = domunion(A,B);
    szA = [privatesize(A) 1];
    niv = A.DomainPrivate.NumIV;
    AData = A.DataPrivate;
    BData = B.DataPrivate;
    if nargout==1
        % E = EIG(A,B) or E = EIG(A,B,'chol') or ...,'qz')
        out = zeros([szA(1) 1 szA(3:end)]);
        for i=1:prod(szA(3:end))
            out(:,1,i) = eig(AData(:,:,i),BData(:,:,i),varargin{3:end});
        end
        varargout{1} = pmat(out,A.DomainPrivate);        
    elseif nargout==2
        % [V,D] = EIG(A,B) or [V,D] = EIG(A,B,'chol') or ...,'qz')
        V = zeros([szA(1) szA(1) szA(3:end)]);
        D = zeros([szA(1) szA(1) szA(3:end)]);
        for i=1:prod(szA(3:end))
            [V(:,:,i),D(:,:,i)] = eig(AData(:,:,i),BData(:,:,i),varargin{3:end});
        end
        varargout{1} = pmat(V,A.DomainPrivate);
        varargout{2} = pmat(D,A.DomainPrivate);        
    end    
else
    % Standard Eigenval/Eigenvec Problem
    szA = [privatesize(A) 1];
    niv = A.DomainPrivate.NumIV;
    AData = A.DataPrivate;
    if nargout==1
        % E = EIG(A) or E = EIG(X,'nobalance')
        out = zeros([szA(1) 1 szA(3:end)]);
        for i=1:prod(szA(3:end))
            out(:,1,i) = eig(AData(:,:,i),varargin{2:end});
        end
        varargout{1} = pmat(out,A.DomainPrivate);        
    elseif nargout==2
        % [V,D] = EIG(A) or [V,D] = EIG(X,'nobalance')
        V = zeros([szA(1) szA(1) szA(3:end)]);
        D = zeros([szA(1) szA(1) szA(3:end)]);
        for i=1:prod(szA(3:end))
            [V(:,:,i),D(:,:,i)] = eig(AData(:,:,i),varargin{2:end});
        end
        varargout{1} = pmat(V,A.DomainPrivate);
        varargout{2} = pmat(D,A.DomainPrivate);        
    end
end


