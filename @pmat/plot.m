function H = plot(varargin)
% PLOT   Linear plot for PMAT objects.
%
% PLOT(X,Y) plots the vector PMAT Y versus the vector PMAT X at each point 
% in the combined domains of X and Y.
%
% PLOT(Y) plots the columns of Y versus their index at each point in the
% domain of Y.
%
% PLOT(X1,Y1,S1,X2,Y2,S2,X3,Y3,S3,...) combines the plots defined by
% the (X,Y,S) triples where the X's and Y's are PMATs or DOUBLEs and the 
% S's are strings specifying the line types, plot symbols and colors.
%
% PLOT(AX,...) plots into the axes with handle AX.
%
% The (X,Y) pairs or (X,Y,S) triples can be followed by  parameter/value
% pairs to specify additional properties of the lines.
%
% See also: plot, semilogx, semilogy, loglog, lpvplot, rcplot.

% Call plot engine
if nargout>0
    H = ploteng('plot',varargin{:});
else
    ploteng('plot',varargin{:});
end

