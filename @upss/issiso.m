function s = issiso(sys)
% ISSISO  True for SISO UPSS objects.
%  
% ISSISO(M) returns true if the UPSS M is single-input and single-output 
% (SISO), and false otherwise.
%  
% See also: issiso.

s = issiso(sys.DataPrivate);

    
