function out = blkdiag(varargin)
% BLKDIAG   Block diagonal concatenation of PMATs
%
% M = BLKDIAG(M1,M2, ...) returns the block diagonal concatenation of
% M1, M2,... at each point in the combined domains of M1, M2,...
%
% See also: blkdiag, diag, horzcat, vertcat.

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    out = binop(varargin{1},varargin{2},'blkdiag');
    if nargin>2
        out = blkdiag(out,varargin{3:end});
    end    
end

