function [u,s,v] = svd(mat,varargin)
% SVD   Singular value decomposition of a PMAT object.
%
% [U,S,V] = SVD(X) computes the singular value decomposition of X at
% each point in the domain of X.
%
% S = SVD(X) returns a PMAT vector S containing the singular values of X
% at each point in the domain of X.
%
% [U,S,V] = SVD(X,0) and [U,S,V] = SVD(X,'econ') return economy-sized
% decompositions as documented in the SVD help for double matrices.
%
% See also: svd.

% TODO PJS 4/2/2011: Implement functions listed in the 
% standard "See also: svd, svds, gsvd".

% Input / output error checking
nin = nargin;
nout = nargout;
narginchk(1, 2)

szm = privatesize(mat);
min12 = min(szm(1:2));

niv = mat.DomainPrivate.NumIV;
Data = mat.DataPrivate;
if nout==1
    %  S = SVD(X)
    s = zeros([min12 1 szm(3:2+niv)]);
    for i=1:prod(szm(3:end))
        s(:,1,i) = svd(Data(:,:,i));
    end
    u = pmat(s,mat.DomainPrivate);
else
    % [U,S,V] = SVD(X) or [U,S,V] = SVD(X,0) or ...,'econ')
    m = szm(1);
    n = szm(2);
    Zmm = zeros([m m szm(3:2+niv)]);
    Zmn = zeros([m n szm(3:2+niv)]);
    Znn = zeros([n n szm(3:2+niv)]);
    Znm = zeros([n m szm(3:2+niv)]);
        
    % Initialize output matrices based on calling syntax
    if nin==1
        % [U,S,V] = SVD(X)
        u = Zmm; s = Zmn; v = Znn;
    elseif ischar( varargin{2} )
        % [U,S,V] = SVD(X,'econ')
        if m>n
            u = Zmn; s = Znn; v = Znn;
        elseif m==n
            u = Zmm; s = Zmn; v = Znn;                        
        else
            u = Zmm; s = Zmm; v = Znm;            
        end
    else
        % [U,S,V] = SVD(X,0)
        if isa( varargin{2} , 'pmat')
            varargin{2} = double(varargin{2});
        end
        if m>n
            u = Zmn; s = Znn; v = Znn;
        else
            u = Zmm; s = Zmn; v = Znn;            
        end
    end
    
    % Compute SVD at each point in domain of M
    for i=1:prod(szm(3:end))
        [u(:,:,i),s(:,:,i),v(:,:,i)] = svd(Data(:,:,i),varargin{:});
    end
    u = pmat(u,mat.DomainPrivate);
    s = pmat(s,mat.DomainPrivate);
    v = pmat(v,mat.DomainPrivate);
end

