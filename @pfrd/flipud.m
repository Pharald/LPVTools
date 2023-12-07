function b = flipud(a)
% FLIPUD   Flip output channels in up/down direction for PFRD objects.
%
% FLIPUD(A) preserves the columns of A and flips the rows of A in 
% the up/down direction at each point in the domain of A.
%
% See also: flipud, fliplr.

% TODO PJS 4/1/2011: Implement rot90 and flipdim?

b = a;
b.DataPrivate  = flipdim(a.DataPrivate,1);
