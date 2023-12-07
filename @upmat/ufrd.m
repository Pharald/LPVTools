function out = ufrd(Response,Freqs,varargin)
% UFRD  Call UPFRD constructor
% 
% SYS = UFRD(RESPONSE,FREQS) creates a UPFRD, a uncertain parameter-varying 
%   frequency response data model defined on an N-dimensional rectangular 
%   grid. FREQS is a Nf-by-1 DOUBLE of frequency values. RESPONSE is a UPMAT
%   Ny-by-Nu-by-Nf complex matrix which defines the Ny-by-Nu frequency
%   response of SYS at the Nf frequency values specified in FREQS. The
%   resulting UPFRD will be defined on the domain of RESPONSE.
%
% SYS = UFRD(RESPONSE,FREQS,TS) creates a discrete time UPFRD model with 
%   sampling time TS. If TS is not specified the UPFRD will be continous-time.
%
% See also: ufrd, upfrd.

if nargin==2
   out = upfrd(Response,Freqs);
else
   out = upfrd(Response,Freqs,varargin{:});
end