function varargout = plot3(varargin)
% PLOT3  3-D plot for PMAT objects.
%
% PLOT3(X,Y,Z) produces a three dimensional plot of a line drawn through  
% the coordinates specified by the elements in the vectors X, Y and Z. 
% The vectors X, Y and Z must all have the same length. Plots from each 
% grid point in the combined domains of X, Y, and Z are overlaid on a 
% single figure. See graph3d\plot3 for additional details and call options.
%
% See also: plot3.

% TODO PJS 5/12/2011: Implement SURF, MESH
% (We already have a version of PMAT/SURF but it is not overloaded in
%  pointwise over the domain).


% Check if first input in varargin is an axis handle.
% Note: The syntax ISHGHANDLE(H,'axes') is undocumented.  However, the
% AXESCHECK function is documented and performs a similar task of
% splitting out an AXES object.
if isscalar(varargin{1}) && ishghandle(varargin{1},'axes')
    aflg = 1;
    AH = varargin{1};
    varargin(1)=[];
else
    aflg = 0;
end

% Find locations of character strings and plot data
s = cellfun(@ischar, varargin);
ns = length(s);
H = [];
% PLOT(X,Y,S...,P1,V1,...) syntax

% Find property value/pairs
idx = ns;
go = 1;
while go==1
    if s(idx-1)==1
        idx = idx-2;
    else
        go=0;
    end
end
Data = varargin(1:idx);
if isempty(Data)
    error('Invalid first data argument')
end
PVPairs = varargin(idx+1:end);

% Plot data at each point in the combined domains
nv = length(Data);
while nv>0
    if nv<=2
        error('Data must be a list of (X,Y,Z) triples.');
    elseif nv==3 || ~ischar(Data{3})
        X = Data{1};
        Y = Data{2};
        Z = Data{3};
        SPV = PVPairs;
        Data(1:3) = [];
    else
        X = Data{1};
        Y = Data{2};
        Z = Data{3};
        SPV = [Data(4) PVPairs];
        Data(1:3) = [];
    end
    nv = length(Data);
    
    [Xext,Yext,Zext] = domunion( pmat(X), pmat(Y), pmat(Z) );
    Xdata = Xext.DataPrivate;
    Ydata = Yext.DataPrivate;
    Zdata = Zext.DataPrivate;
    szX = [privatesize(Xext) 1];
    
    if aflg==1
        for i=1:prod(szX(3:end))
            Hi = plot3(AH,Xdata(:,:,i),Ydata(:,:,i),Zdata(:,:,i),SPV{:});
            hold(AH,'on');
            H = [H; Hi];
        end
        hold(AH,'off');
    else
        for i=1:prod(szX(3:end))
            Hi = plot3(Xdata(:,:,i),Ydata(:,:,i),Zdata(:,:,i),SPV{:});
            hold on;
            H = [H; Hi];
        end
        hold off;
    end
    
end

% Store plot handle for output
if nargout>0
    varargout = {H};
end



