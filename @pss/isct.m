function s = isct(sys)
% ISCT  True for continuous-time PSS.
%  
% ISCT(SYS) returns true if the PSS SYS is continuous-time, otherwise false.
%  
% See also: isct, isdt, isstatic.

s = isct(sys.DataPrivate);

    
