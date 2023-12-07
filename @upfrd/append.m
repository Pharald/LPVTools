function out = append(varargin)
% APPEND   Block diagonal concatenation of UPFRD
%
% S = APPEND(S1,S2, ...) returns the block diagonal concatenation of
% S1, S2,... at each point in the combined domains of S1, S2,...
% This is the same functionality as BLKDIAG.
%
% See also:  append, blkdiag, series, parallel, feedback.

try
   out = blkdiag(varargin{:});
catch
   error(lasterr)
end
