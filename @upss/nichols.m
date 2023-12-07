function out = nichols(varargin)
% NICHOLS   Pointwise Nichols frequency response for UPSS objects
%
% NICHOLS(SYS) draws a family of Nichols plots of SYS at each point in the 
% domain of SYS.
%
% NICHOLS(SYS,{WMIN,WMAX}) draws the Nichols plots for frequencies 
% between WMIN and WMAX (in radians/second).
%
% NICHOLS(SYS,W) uses the user-supplied vector W of frequencies, in
% radians/second, at which the Nichols response is to be evaluated.  
%
% NICHOLS(SYS1,SYS2,...,W) graphs the Nichols responses of several systems
% SYS1,SYS2,... on a single plot. 
%
% See also: nichols, bode, bodemag, nyquist, sigma, freqresp.

% TODO PJS 5/1/2011: Revisit the handling of the output arguments. 

nin = nargin;
nout = nargout;
if nout>0
   strerr = ' with a UPSS or UPFRD cannot be called with output arguments.';
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

nichols(incell{:});



