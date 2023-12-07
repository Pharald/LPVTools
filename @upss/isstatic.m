function s = isstatic(sys)
% ISSTATIC  Checks if UPSS is static or dynamic.
%  
% ISSTATIC(SYS) returns TRUE if the UPSS SYS is static and FALSE if
% SYS has dynamics.
%  
% See also: isstatic, isct, isdt.

s = isstatic(sys.DataPrivate);

    
