function B = usubs(A,varargin)
% USUBS   Substitue uncertainty for values in a UPMAT
%
% P = USUBS(G,NAME1,VAL1,NAME2,VAL2,...) substitutes the uncertainty in G
% with specific values. See help in InputOutputModel/usubs for details of
% valid function calls.
% 
% See also: usubs, usample, lpvsubs, lpvsample.


B = upmat(usubs(A.DataPrivate,varargin{:}),A.DomainPrivate);