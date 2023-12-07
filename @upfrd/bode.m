function out = bode(varargin)
% BODE  Pointwise Bode frequency response for UPFRD objects
%
% BODE(SYS) draws Bode plots of SYS at each point in the domain of SYS.
%
% BODE(SYS,{WMIN,WMAX}) draws the Bode plot for frequencies between WMIN 
% and WMAX (in radians/second).
%
% BODE(SYS,W) uses the user-supplied vector W of frequencies, in
% radian/second, at which the Bode response is to be evaluated.  
%
% BODE(SYS1,SYS2,...,W) graphs the Bode response of several systems
% SYS1,SYS2,... on a single plot. 
%
% See also: bode, bodemag, nichols, nyquist, sigma, freqresp.

% TODO PJS 5/1/2011: Revisit the handling of the output arguments. 


nin = nargin;
nout = nargout;
if nout>0
   strerr = ' with an UPSS or UPFRD cannot be called with output arguments.';
   error([upper(mname) strerr])
end

incell = cell(1,nin);
for i=1:nin
   if isa(varargin{i},'upss')
      incell{i} = varargin{i}.DataPrivate;
   elseif isa(varargin{i},'upfrd')
      incell{i} = varargin{i}.DataPrivate;
   else
      incell{i} = varargin{i};
   end
end

bode(incell{:});



