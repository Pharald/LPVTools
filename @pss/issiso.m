function s = issiso(sys)
% ISSISO  True for SISO PSS objects.
%  
% ISSISO(M) returns true if the PSS M is single-input and single-output 
% (SISO), and false otherwise.
%  
% See also: issiso.

s = issiso(sys.DataPrivate);

    
