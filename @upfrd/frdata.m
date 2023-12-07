function [Response,Freqs,Ts] = frdata(S)
% FRDATA  Access to frequency response data of a UPFRD
% 
% [RESPONSE,FREQ,TS] = FRDATA(SYS) returns the response data, frequency
% samples, and sampling time Ts of the UPFRD.  The response data is 
% returned as a UPMAT.
%
% See also: frdata.

[Response,Freqs,Ts] = frdata(S.DataPrivate);
Response = upmat(Response,S.DomainPrivate);
