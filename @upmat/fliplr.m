function b = fliplr(a)
% FLIPLR   Flip matrix in left/right direction for UPMAT objects.
%
% FLIPLR(A) preserves the rows of A and flips the columns of A in 
% the left/right direction at each point in the domain of A.
%
% See also: fliplr, flipud.

b = a;
b.DataPrivate  = flipdim(a.DataPrivate,2);
