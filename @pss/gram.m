function W = gram(sys,type)
% GRAM   Pointwise computes Gramians for PSS objects
%
% Wc = GRAM(SYS,'c') computes the controllability gramians of the PSS
% SYS at each point in the domain of SYS.
% 
% Wo = GRAM(SYS,'o') computes its observability gramians.
%
% Rc = GRAM(SYS,'cf') returns the Cholesky factor of Wc.
%
% Ro = GRAM(SYS,'of') returns the Cholesky factor of Wo.
%
% See also: gram, balreal, lpvgram.

W = gram(sys.DataPrivate,type);
W = pmat(W,sys.DomainPrivate);
