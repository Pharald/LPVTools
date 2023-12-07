function K = dcgain(S)
% DCGAIN  Pointwise DC gain of a UPSS object
%
% K = DCGAIN(SYS) computes the steady-state gain of the UPSS SYS at each
% point in the domain of SYS. If SYS has NY outputs and NU inputs then
% K is returned as an NY-by-NU PMAT containing the steady-state gains of 
% SYS at each point in the domain of SYS.
%
% See also: dcgain, norm, evalfr, freqresp.


% USS/DCGAIN works on uncertain state space arrays
K = dcgain(S.DataPrivate);
K = pmat(K,S.DomainPrivate);