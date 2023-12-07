function out = frd(Response,Freqs,varargin)
% FRD  Call PFRD constructor
% 
% SYS = FRD(RESPONSE,FREQS) creates a PFRD, a parameter-varying frequency
%   response data model defined on an N-dimensional rectangular grid.
%   FREQS is a Nf-by-1 DOUBLE of frequency values. RESPONSE is a PMAT
%   Ny-by-Nu-by-Nf complex matrix which defines the Ny-by-Nu frequency
%   response of SYS at the Nf frequency values specified in FREQS. The
%   resulting PFRD will be defined on the domain of RESPONSE.
%
% SYS = FRD(RESPONSE,FREQS,TS) creates a discrete time PFRD model with sampling
%   time TS. If TS is not specified the PFRD will be continous-time.
%
% See also: frd, pfrd.

if nargin==2
   out = pfrd(Response,Freqs);
else
   out = pfrd(Response,Freqs,varargin{:});
end
