function s = isdt(sys)
% ISDT  True for discrete-time UPSS.
%  
% ISDT(SYS) returns true if the UPSS SYS is discrete-time
% and false otherwise. Returns true for empty systems and static gains.
%  
% See also: isct, isdt, isstatic.

s = isdt(sys.DataPrivate);
