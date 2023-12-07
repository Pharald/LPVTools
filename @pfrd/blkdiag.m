function out = blkdiag(varargin)
% BLKDIAG   Block diagonal concatenation of PFRDs.
%
% S = BLKDIAG(S1,S2, ...) returns the block diagonal concatenation of
% S1, S2,... at each point in the combined domains of S1, S2,...
%
% See also:  blkdiag, append, series, parallel, feedback.

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



% ----- OLD CODE

% function out = blkdiag(varargin)
% 
% if nargin==1
%    out = varargin{1};
% elseif nargin==2
%    if isa(varargin{1},'pss') && isa(varargin{2},'pss')
%       try
%          out = binop(varargin{1},varargin{2},'blkdiag');
%       catch
%          error(lasterr)
%       end
%    else
%       out = blkdiag(pss(varargin{1}),pss(varargin{2}));
%    end
% else
%    out = blkdiag(blkdiag(pss(varargin{1}),pss(varargin{2})),varargin{3:end});
% end

