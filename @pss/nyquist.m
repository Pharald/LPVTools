function out = nyquist(varargin)
% NYQUIST   Pointwise Nyquist frequency response for PSS objects
%
% NYQUIST(SYS) draws Nyquist plots of SYS at each point in the domain 
% of SYS.
%
% NYQUIST(SYS,{WMIN,WMAX}) draws the Nyquist plot for frequencies 
% between WMIN and WMAX (in radians/second).
%
% NYQUIST(SYS,W) uses the user-supplied vector W of frequencies, in
% radians/second, at which the Nyquist response is to be evaluated.  
%
% NYQUIST(SYS1,SYS2,...,W) graphs the Nyquist response of several systems
% SYS1,SYS2,... on a single plot. 
%
% See also: nichols, bode, bodemag, nyquist, sigma, freqresp.

% TODO PJS 5/1/2011: Revisit the handling of the output arguments. 

nin = nargin;
nout = nargout;
if nout>0
   strerr = ' with a PSS or PFRD cannot be called with output arguments.';
   error([upper(mname) strerr])
end

incell = cell(1,nin);
for i=1:nin
   if isa(varargin{i},'pss')
      incell{i} = varargin{i}.DataPrivate;
   elseif isa(varargin{i},'pfrd')
      incell{i} = varargin{i}.DataPrivate;
   else
      incell{i} = varargin{i};
   end
end

nyquist(incell{:});



