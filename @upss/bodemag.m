function out = bodemag(varargin)
% BODEMAG  Bode magnitude plot for UPSS objects
%
% BODEMAG(SYS) plots the magnitude of the frequency response of the PSS
% SYS at each point in the domain of SYS.
%
% BODEMAG(SYS,{WMIN,WMAX}) draws the magnitude plot for frequencies 
% between WMIN and WMAX (in radians/second).
%
% BODEMAG(SYS,W) uses the user-supplied vector W of frequencies, in
% radian/second, at which the frequency response is to be evaluated.  
%
% BODEMAG(SYS1,SYS2,...,W) graphs the magnitude response of several 
% systems SYS1,SYS2,... on a single plot. 
%
% See also: bodemag, bode, nichols, nyquist, sigma, freqresp.

% TODO PJS 5/1/2011: Revisit the handling of the output arguments. 

nin = nargin;
nout = nargout;
if nout>0
   strerr = ' with a PSS/PFRD/UPSS cannot be called with output arguments.';
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

bodemag(incell{:});



