function out = power(A,B)
% POWER   Array power for PMAT objects
%
% POWER(A,B) is the result of A.^B at each point in the combined
% domains of A and B.
%
% See also: power, mpower.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'power');

% % Check input arguments
% error(nargchk(2, 2, nargin, 'struct'))
% 
% % Define A and B on a common domain
% [Aext,Bext] = domunion(A,B);
% 
% % Evaluate atan2 at each point in the combined domain
% out = Aext;
% out.Data = realpow(Aext.Data,Bext.Data);

