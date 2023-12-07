function out = nichols(varargin)
% NICHOLS   Pointwise Nichols frequency response for PLFTSS objects
%
% NICHOLS(SYS,DOMAIN) draws a family of Nichols plots of SYS at each point 
% in the domain specified in the RGRID object DOMAIN.
%
% NICHOLS(SYS,{WMIN,WMAX},DOMAIN) draws the Nichols plots for frequencies 
% between WMIN and WMAX (in radians/second).
%
% NICHOLS(SYS,W,DOMAIN) uses the user-supplied vector W of frequencies, in
% radians/second, at which the Nichols response is to be evaluated.  
%
% NICHOLS(SYS1,SYS2,...,W,DOMAIN) graphs the Nichols responses of several 
% systems SYS1,SYS2,... on a single plot. 
%
% See also: nichols, bode, bodemag, nyquist, sigma, freqresp.

DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end
for i = 1:numel(varargin)
    if isa(varargin{i},'plftss')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

if nargout ==0
    nichols(varargin{:});
else
    out = nichols(varargin{:});
end