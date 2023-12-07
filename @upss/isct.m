function s = isct(sys)
% ISCT  True for continuous-time UPSS.
%  
% ISCT(SYS) returns true if the UPSS SYS is continuous-time, otherwise false.
%  
% See also: isct, isdt, isstatic.

s = isct(sys.DataPrivate);

    
