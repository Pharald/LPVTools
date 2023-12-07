function H = ploteng(ptype,varargin)
% PLOTENG   Engine for plotting PMAT objects
%
% PLOTENG(PTYPE,X,Y) plots the vector PMAT Y versus the vector PMAT X at
% each point in the combined domains of X and Y. PTYPE is a character
% string specifying the plot type and can be 'plot', 'semilogx',
% 'semilogy', 'loglog', or 'polar'.
%
% PLOT(PTYPE,Y) plots the columns of Y versus their index at each point
% in the domain of Y.
%
% PLOT(PTYPE,X1,Y1,S1,X2,Y2,S2,X3,Y3,S3,...) combines the plots defined by
% the (X,Y,S) triples where the X's and Y's are PMATs or DOUBLEs and the
% S's are strings specifying the line types, plot symbols and colors.
%
% PLOT(PTYPE,AX,...) plots into the axes with handle AX.
%
% The (X,Y) pairs or (X,Y,S) triples can be followed by  parameter/value
% pairs to specify additional properties of the lines.
%
% See also: plot, semilogx, semilogy, loglog.

% TODO PJS 5/10/2011: Implement plotyy? What about graph3d functions?

% Get function handle for the plot type
if ~ismember(ptype,{'plot';'semilogx';'semilogy';'loglog'; 'polar'})
    error(['PTYPE must be ''plot'', ''semilogx'', ''semilogy'', '...
        '''loglog'', or ''polar''.']);
end
phandle = str2func(ptype);

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
if ns==1 || s(2)==1
    % PLOT(Y,S,P1,V1,...) or PLOT(Y,P1,V1,...)
    Y = varargin{1};
    Ydata = Y.DataPrivate;
    szY = [privatesize(Y) 1];
    if aflg==1
        % Get current hold status
        hflg = ishold(AH);
            
        % Plot data
        for i=1:prod(szY(3:end))
            Hi = phandle(AH,Ydata(:,:,i),varargin{2:end});
            hold(AH,'on');
            H = [H; Hi];
        end

        % Set original hold status        
        if ~hflg
            hold(AH,'off');
        end
    else
        % Get current hold status
        hflg = ishold;
        
        % Plot data
        for i=1:prod(szY(3:end))
            Hi = phandle(Ydata(:,:,i),varargin{2:end});
            hold on;
            H = [H; Hi];
        end
        
        % Set original hold status
        if ~hflg
            hold off;
        end
    end
    
else
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

    % Get current hold status
    if aflg==1
        hflg = ishold(AH);
    else
        hflg = ishold;
    end
    
    % Plot data at each point in the combined domains
    nv = length(Data);
    while nv>0
        if nv==1
            error('Data must be a single matrix Y or a list of pairs X,Y');
        elseif nv==2 || ~ischar(Data{3})
            X = Data{1};
            Y = Data{2};
            SPV = PVPairs;
            Data(1:2) = [];
        else
            X = Data{1};
            Y = Data{2};
            SPV = [Data(3) PVPairs];
            Data(1:3) = [];
        end
        nv = length(Data);
        
        [Xext,Yext] = domunion( pmat(X), pmat(Y) );
        Xdata = Xext.DataPrivate;
        Ydata = Yext.DataPrivate;
        szX = [privatesize(Xext) 1];
        
        if aflg==1
            for i=1:prod(szX(3:end))
                Hi = phandle(AH,Xdata(:,:,i),Ydata(:,:,i),SPV{:});
                hold(AH,'on');
                H = [H; Hi];
            end
        else
            for i=1:prod(szX(3:end))
                Hi = phandle(Xdata(:,:,i),Ydata(:,:,i),SPV{:});
                hold on;
                H = [H; Hi];
            end
        end        
    end
    
    % Set original hold status
    if ~hflg
        if aflg==1
            hold(AH,'off');
        else
            hold off;
        end
    end        

end


% % Simple code for PLOT(X,Y,P1,V1,....) syntax
% [Xext,Yext] = domunion(X,Y);
% Xdata = Xext.DataPrivate;
% Ydata = Yext.DataPrivate;
% szX = [size(Xext) 1];
% npts = prod(szX(3:end));
% for i=1:npts
%    plot(Xdata(:,:,i),Ydata(:,:,i),varargin{:}); hold on;
% end
% hold off;
