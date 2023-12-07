function b = iscolumn(M)
% ISCOLUMN   True if input is a column vector UPMAT object
%
% ISCOLUMN(M) returns 1 (true) if the UPMAT M is a column vector, i.e. 
% S=SIZE(M) has S(2) == 1.
%
% See also: iscolumn,  isscalar, isvector, isrow.

szM = size(M);
if szM(2)==1 
    b = true;
else
    b = false;
end

