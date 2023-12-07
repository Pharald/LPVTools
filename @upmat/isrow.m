function b = isrow(M)
% ISROW   True if input is a row vector UPMAT object
%
% ISROW(M) returns 1 (true) if the UPMAT M is a row vector, i.e. 
% S=SIZE(M) has S(1) == 1.
%
% See also: isrow,  isscalar, isvector, iscolumn.

szM = size(M);
if szM(1)==1 
    b = true;
else
    b = false;
end

