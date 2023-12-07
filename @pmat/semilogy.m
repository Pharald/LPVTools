function H = semilogy(varargin)
% SEMILOGY   Semi-log scale plot for PMAT objects.
%
% SEMILOGY(X,Y) plots the vector PMAT Y versus the vector PMAT X at each 
% point in the combined domains of X and Y together on a single plot. 
% A logarithmic (base 10) scale is used for the Y-axis.
%
% See also: semilogy, plot, semilogx, loglog.

% Call plot engine
if nargout>0
    H = ploteng('semilogy',varargin{:});
else
    ploteng('semilogy',varargin{:});
end
