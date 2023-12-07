function out = nyquist(varargin)
% NYQUIST   Pointwise Nyquist frequency response for PLFTSS objects
%
% NYQUIST(SYS,DOMAIN) draws Nyquist plots of SYS at each point in the domain  
% specified in the RGRID object DOMAIN.
%
% NYQUIST(SYS,{WMIN,WMAX},DOMAIN) draws the Nyquist plot for frequencies 
% between WMIN and WMAX (in radians/second).
%
% NYQUIST(SYS,W,DOMAIN) uses the user-supplied vector W of frequencies, in
% radians/second, at which the Nyquist response is to be evaluated.  
%
% NYQUIST(SYS1,SYS2,...,W,DOMAIN) graphs the Nyquist response of several 
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
    nyquist(varargin{:});
else
    out = nyquist(varargin{:});
end