function varargout = svd(varargin)
% SVD   Singular value decomposition of a PLFTMAT object.
%
% [U,S,V] = SVD(X,DOMAIN) computes the singular value decomposition of X 
% at each point in the domain specified by the RGRID object DOMAIN.
%
% S = SVD(X,DOMAIN) returns a PMAT vector S containing the singular values 
% of X at each point in the domain specified by the RGRID object DOMAIN.
%
% [U,S,V] = SVD(X,0,DOMAIN) and [U,S,V] = SVD(X,'econ',DOMAIN) return 
% economy-sized decompositions as documented in the SVD help for double 
% matrices.
%
% See also: svd.

DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plftss') || isa(varargin{i},'plftmat')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

varargout = cell(nargout,1);
[varargout{:}]=svd(varargin{:});