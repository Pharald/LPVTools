function g = frd(sys,omega)
% FRD  Convert a UPSS to a UPFRD
% 
% G = FRD(SYS,FREQS) converts the uncertain parameter dependent system SYS 
%   into a uncertain parameter-varying frequency response data model G, 
%   defined on the domain of SYS. FREQS is a Nf-by-1 DOUBLE of frequency 
%   values. the output G is a UPFRD containing the frequency response of 
%   SYS at each point in FREQS.
%
% See also: ufrd, upfrd.

tmp = frd(sys.DataPrivate,omega);
g = pfrd(tmp,sys.DomainPrivate);
