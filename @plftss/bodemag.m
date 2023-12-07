function out = bodemag(varargin)
% BODEMAG  Bode magnitude plot for PLFTSS objects
%
% BODEMAG(SYS,DOMAIN) plots the magnitude of the frequency response of the 
% PLFTSS SYS at each point in the domain specified in the RGRID object 
% DOMAIN.
%
% BODEMAG(SYS,{WMIN,WMAX},DOMAIN) draws the magnitude plot for frequencies 
% between WMIN and WMAX (in radians/second).
%
% BODEMAG(SYS,W,DOMAIN) uses the user-supplied vector W of frequencies, in
% radian/second, at which the frequency response is to be evaluated.  
%
% BODEMAG(SYS1,SYS2,...,W,DOMAIN) graphs the magnitude response of several 
% systems SYS1,SYS2,... on a single plot. 
%
% See also: bodemag, bode, nichols, nyquist, sigma, freqresp.

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
    bodemag(varargin{:});
else
    out = bodemag(varargin{:});
end