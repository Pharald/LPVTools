function out = uminus(A)
% UMINUS   Unary minus for PFRD objects.
%
% -M is the unary minus of M at each point in the domain of M.
%
% See also: uminus, uplus.
    
out = A;
out.DataPrivate = -out.DataPrivate;
