function s = isdt(sys)
% ISDT  True for discrete-time PSS.
%  
% ISDT(SYS) returns true if the PSS SYS is discrete-time
% and false otherwise. Returns true for empty systems and static gains.
%  
% See also: isdt, isct, isstatic.

s = isdt(sys.DataPrivate);

    
