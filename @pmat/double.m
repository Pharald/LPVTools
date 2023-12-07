function a = double(a)
% DOUBLE   Return Data of a PMAT object.
%
% DOUBLE(M) returns the M.DATA for the PMAT object M.
%
% See also: depromote.
a = pmat(a);
a = a.DataPrivate;
