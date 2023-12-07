function [A,B,C,D,Ts] = ssdata(S)
% SSDATA  Access to state-space data of a PSS 
%
% [A,B,C,D] = SSDATA(SYS) returns the the A, B, C, D matrices of the PSS.  
% The state matrices are returned as PMATs.
% 
% [A,B,C,D,Ts] = SSDATA(SYS) also returns the sampling time Ts.
%
% See also: ssdata.

[a,b,c,d,Ts] = ssdata(S.DataPrivate);
A = pmat(a,S.DomainPrivate);
B = pmat(b,S.DomainPrivate);
C = pmat(c,S.DomainPrivate);
D = pmat(d,S.DomainPrivate);
