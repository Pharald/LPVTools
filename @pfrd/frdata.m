function [Response,Freqs,Ts] = frdata(S)
% FRDATA  Access to frequency response data of a PFRD
% 
% [RESPONSE,FREQ,TS] = FRDATA(SYS) returns the response data, frequency
% samples, and sampling time Ts of the PFRD.  The response data is returned
% as a PMAT.
%
% See also: frdata.

[Response,Freqs,Ts] = frdata(S.DataPrivate);
Response = pmat(Response,S.DomainPrivate);
