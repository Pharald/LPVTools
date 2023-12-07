function varargout = eig(varargin)
% EIG   Eigenvalues and eigenvectors for a PLFTMAT object.
%
% E = EIG(X,DOMAIN) returns a PMAT vector containing the eigenvalues of 
% a square PLFTMAT X at each point in the domain specified by the RGRID
% object DOMAIN.
%
% [V,D] = EIG(X,DOMAIN) returns a diagonal PMAT D of eigenvalues and a
% full PMAT V whose columns are the corresponding eigenvectors so
% that X*V = V*D at each point in the domain specified by the RGRID
% object DOMAIN.
%
% E = EIG(A,B,DOMAIN) and [V,D]=EIG(A,B,DOMAIN) compute the generalized 
% eigenvaluesand eigenvectors at each point in the domain specified by the 
% RGRID object DOMAIN.
%
% The function also supports the following syntax as documented in the
% EIG help:
%    [V,D] = EIG(X,'nobalance')
%    EIG(A,B,'chol')
%    EIG(A,B,'qz')
%
% See also: eig.

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
[varargout{:}]=eig(varargin{:});