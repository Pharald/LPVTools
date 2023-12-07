function B = isuncertain(A)
%ISUNCERTAIN True if object is uncertain.
%
%  B = ISUNCERTAIN(A) is true if A is an uncertain object and false 
%  otherwise. The uncertain parameter-varying objects are UPMAT, UPFRD, 
%  UPSS, PLFTMAT, and PLFTSS. 
%
% See also: isuncertain.

B = isUncertain(A.Data);
