function E = lpvelimiv(M)
% LPVELIMIV  Eliminate singleton independent variables.
%
% E = LPVELIMIV(M) eliminates all singleton independent variables from 
% the PFRD M, i.e. the i^th independent variable of M is removed if
% M.IVData{i} has length one. 
%
% See also: squeeze.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 1)
nargoutchk(0, 1)

% Find singleton dimensions
LIVData = M.Domain.LIVData;
rmidx = find(LIVData==1)';

% Remove singleton dimensions
niv = M.Domain.NumIV;
nad = numel(size(M))-2;
keepidx = setdiff( 1:(niv+nad), rmidx );

% Ordering: Parameter vars, Array dims
permidx = [keepidx, rmidx];
if length(permidx) == 1
    Mdata = M.Data;
else
    Mdata = permute(M.Data,permidx);
end
Dom = lpvelimiv(M.Domain,M.Domain.IVName(rmidx));

% Pack up data
E = pfrd(Mdata,Dom);
