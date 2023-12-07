function b = isempty(M)
% ISEMPTY True for empty UPFRD objects
%
% ISEMPTY(M) returns 1 (true) if the UPFRD M has no rows or no columns, 
% and 0 (false) otherwise.
%
% See also: isempty.

szM = size(M);
if ~all(szM)
    b = true;
else
    b = false;
end

