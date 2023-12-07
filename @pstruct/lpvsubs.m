function B = lpvsubs(A,Names,Values)
% LPVSUBS  Evaluate PSTRUCT at points in the domain, and demote results to STRUCT.
%
% B = LPVSUBS(A,NAME,VALUES) evaluates a PSTRUCT A at the domain points
% specified by the NAME/VALUES pair.  NAME is an N-by-1 cell array of
% characters. NAME must contain all of the independent variable names in
% A, namely A.IVName, but may also contain others but they will be ignored. 
% VALUES is an N-by-NPTS double array of the corresponding values.  B is a 
% struct array, with row and column dimensions from A.  The 3rd dimension 
% of B is NPTS.  
%
% See also: lpvinterp, lpvsplit, lpvsample.

% Parse inputs
nin = nargin;
narginchk(3, 3)

Npts = size(Values,2);
for i=1:Npts
    Valuesi = Values(:,i);
    NVPair = [Names Valuesi]';
    try
        Bi = lpvsplit(A,NVPair{:});
    catch
        % XXX PJS Provide more descriptive error message
        error('Name/Value pairs incorrectly specified')
    end
    B(:,:,i) = depromote(lpvelimiv(Bi));
end

