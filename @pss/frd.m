function g = frd(sys,omega)
% FRD  Convert a PSS to a PFRD
% 
% G = FRD(SYS,FREQS) converts the parameter dependent system SYS into a
%   parameter-varying frequency response data model G defined on the domain 
%   of SYS. FREQS is a Nf-by-1 DOUBLE of frequency values. G is returned as
%   a PFRD containing the frequency response of SYS at each point in FREQS.
%
% See also: frd, pfrd.

tmp = frd(sys.DataPrivate,omega);
g = pfrd(tmp,sys.DomainPrivate);
