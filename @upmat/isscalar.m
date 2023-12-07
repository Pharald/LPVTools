function b = isscalar(M)
% ISSCALAR   True if input is a scalar UPMAT object
%
% ISSCALAR(M) returns 1 (true) if the UPMAT M has one row and one column,
% i.e. S=SIZE(M) has S(1:2) equal to [1 1].
%
% See also: isscalar, isvector, isrow, iscolumn.

szM = size(M);
if szM(1)==1 && szM(2)==1
    b = true;
else
    b = false;
end

