function varargout = schur(varargin)
% SCHUR   Schur decomposition of a PLFTMAT object.
%
% [U,T] = SCHUR(X,DOMAIN) computes the Schur decomposition of X at each 
% point in the domain specified by the RGRID object DOMAIN.
%
% T = SCHUR(X,DOMAIN) returns just the Schur matrix T.
%
% If X is real then SCHUR(X,'real',DOMAIN) and SCHUR(X,'complex',DOMAIN) 
% return the real and complex Schur form as documented the SCHUR help 
% for double matrices.
%
% See also: schur, qz.

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
[varargout{:}]=schur(varargin{:});