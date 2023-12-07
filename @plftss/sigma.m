function out = sigma(varargin)
% SIGMA   Pointwise singular value plot for PLFTSS objects
%
% SIGMA(SYS,DOMAIN)  produces a family of singular value plot of the  
% frequency responses of the PLFTSS SYS at each point in the domain 
% specified in the RGRID object DOMAIN.
%
% SIGMA(SYS,{WMIN,WMAX},DOMAIN) draws the singular value plots for  
% frequencies between WMIN and WMAX (in radians/second).
%
% SIGMA(SYS,W,DOMAIN) uses the user-supplied vector W of frequencies, in
% radians/second, at which the frequency responses is to be evaluated.  
%
% SIGMA(SYS,W,TYPE,DOMAIN) or SIGMA(SYS,[],TYPE,DOMAIN) draws the following
% modified singular value plots depending on the value of TYPE:
%          TYPE = 1     -->     Singular value of  inv(SYS)
%          TYPE = 2     -->     Singular value of  I + SYS
%          TYPE = 3     -->     Singular value of  I + inv(SYS) 
%
% SIGMA(SYS1,SYS2,...,W,DOMAIN) graphs the singular value response of  
% several systems SYS1,SYS2,... on a single plot. 
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
    sigma(varargin{:});
else
    out = sigma(varargin{:});
end