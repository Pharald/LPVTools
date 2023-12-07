function out = blkdiag(varargin)
% BLKDIAG   Block diagonal concatenation of PSS objects.
%
% S = BLKDIAG(S1,S2, ...) returns the block diagonal concatenation of
% S1, S2,... at each point in the combined domains of S1, S2,...
%
% See also:  blkdiag, append, series, parallel, feedback.

% Check # of input arguments
% error(nargchk(1, inf, nargin, 'struct'))  %COMPAT FIX 2016a -JT22.06.2016
narginchk(1, inf) 

if nargin==1
    out = varargin{1};
else
    out = binop(varargin{1},varargin{2},'blkdiag');
    if nargin>2
        out = blkdiag(out,varargin{3:end});
    end    
end

