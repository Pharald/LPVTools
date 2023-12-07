function H = loglog(varargin)
% LOGLOG  Log-log scale plot for PMAT objects.
%
% LOGLOG(X,Y) plots the vector PMAT Y versus the vector PMAT X at each 
% point in the combined domains of X and Y together on a single plot. 
% A logarithmic (base 10) scale is used for the X- and Y-axis.
%
% See also: loglog, plot, semilogx, semilogy.

% Call plot engine
if nargout>0
    H = ploteng('loglog',varargin{:});
else
    ploteng('loglog',varargin{:});
end
