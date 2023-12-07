function b = isempty(M)
% ISEMPTY True for empty UPMAT objects
%
% ISEMPTY(M) returns 1 (true) if the UPMAT M has no rows or no columns, 
% and 0 (false) otherwise.
%
% See also: isempty.

szM = size(M);
if ~all(szM)
    b = true;
else
    b = false;
end

