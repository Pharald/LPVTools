function H = semilogx(varargin)
% SEMILOGX   Semi-log scale plot for PMAT objects.
%
% SEMILOGX(X,Y) plots the vector PMAT Y versus the vector PMAT X at each 
% point in the combined domains of X and Y together on a single plot. 
% A logarithmic (base 10) scale is used for the X-axis.
%
% See also: semilogx, plot, semilogy, loglog.

% Call plot engine
if nargout>0
    H = ploteng('semilogx',varargin{:});
else
    ploteng('semilogx',varargin{:});
end
