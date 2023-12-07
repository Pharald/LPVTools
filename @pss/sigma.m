function out = sigma(varargin)
% SIGMA   Pointwise singular value plot for PSS objects
%
% SIGMA(SYS)  produces a family of singular value plots of the frequency 
% responses of the PSS SYS at each point in the domain of SYS.
%
% SIGMA(SYS,{WMIN,WMAX}) draws the singular value plots for frequencies 
% between WMIN and WMAX (in radians/second).
%
% SIGMA(SYS,W) uses the user-supplied vector W of frequencies, in
% radians/second, at which the frequency responses is to be evaluated.  
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

sigma(incell{:});



