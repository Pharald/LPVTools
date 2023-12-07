function out = bode(varargin)
% BODE  Pointwise Bode frequency response for PLFTSS objects
%
% BODE(SYS,DOMAIN) draws Bode plots of SYS at each point in the domain  
% specified in the RGRID object DOMAIN.
%
% BODE(SYS,{WMIN,WMAX},DOMAIN) draws the Bode plot for frequencies between  
% WMIN and WMAX (in radians/second).
%
% BODE(SYS,W,DOMAIN) uses the user-supplied vector W of frequencies, in
% radian/second, at which the Bode response is to be evaluated.  
%
% BODE(SYS1,SYS2,...,W,DOMAIN) graphs the Bode response of several systems
% SYS1,SYS2,... on a single plot. 
%
% See also: bode, bodemag, nichols, nyquist, sigma, freqresp.

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
    bode(varargin{:});
else
    out = bode(varargin{:});
end



