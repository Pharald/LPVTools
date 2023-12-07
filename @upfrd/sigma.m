function out = sigma(varargin)
% SIGMA   Pointwise singular value plot for UPFRD objects
%
% SIGMA(SYS)  produces a singular value plot of the frequency response
% of the UPFRD SYS at each point in the domain of SYS.
%
% SIGMA(SYS,{WMIN,WMAX}) draws the singular value plot for frequencies 
% between WMIN and WMAX (in radians/second).
%
% SIGMA(SYS,W) uses the user-supplied vector W of frequencies, in
% radians/second, at which the frequency response is to be evaluated.  
%
% SIGMA(SYS,W,TYPE) or SIGMA(SYS,[],TYPE) draws the following
% modified singular value plots depending on the value of TYPE:
%          TYPE = 1     -->     Singular value of  inv(SYS)
%          TYPE = 2     -->     Singular value of  I + SYS
%          TYPE = 3     -->     Singular value of  I + inv(SYS) 
%
% SIGMA(SYS1,SYS2,...,W) graphs the singular value response of several 
% systems SYS1,SYS2,... on a single plot. 
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

sigma(incell{:});



