function [ugrid,varargout] = domunion(varargin)
% DOMUNION  Construct unique rectangular grid
%
% [Ugrid,Aidx,Bidx] = DOMUNION(A,B) constructs a single rectangular
% grid from RGRID objects A and B.  Ugrid is an RGRID object whose
% set of independent variables is the union of the variables in A and B. 
% Aidx and Bidx provide a mapping of variables in A and B into the
% ordering used by Ugrid. Specifically, Ugrid.IVName is equal to
% VA(Aidx) and VB(Bidx) where VA and VB are defined as:
%    VA = [A.IVName; setdiff(Ugrid.IVName,A.IVName)]
%    VB = [B.IVName; setdiff(Ugrid.IVName,B.IVName)]
%
% The syntax [Ugrid,idx1,idx2,...,idxN] = DOMUNION(A1,...,AN) constructs 
% a single rectangular grid from RGRID objects A1, A2, ..., AN. The
% set of independent variables in Ugrid is the union of the variables in 
% A1, ..., AN. idx1, ..., idxN provide a mapping from the independent 
% variables in the input domains into the variables in Ugrid. For example,
% Ugrid.IVName is equal to V1(idx1) where V1 is defined as:
%    V1 = [A1.IVName; setdiff(Ugrid.IVName,A1.IVName)]

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, nin+1)
ArrayName = varargin{1}.ArrayName;

% Find union of all independent variables variables
uIVNameAll = [];
uIVDataAll = [];
uIVRBAll = [];
for i1=1:nin
    uIVNameAll = [uIVNameAll; varargin{i1}.IVName];
    uIVDataAll = [uIVDataAll; varargin{i1}.IVData];
    uIVRBAll = [uIVRBAll; varargin{i1}.IVRateBounds];
end
[uIVName,idx,jdx] = unique(uIVNameAll);
uIVData = uIVDataAll(idx);
uIVRB = uIVRBAll(idx,:);
% Make sure short list has max value of IVData for array dims
adidx = find(strncmp(ArrayName,uIVName,length(ArrayName)));
for ii=1:numel(adidx)
   newidx = find(jdx==adidx(ii));
   tmax = max(cellfun(@numel,uIVDataAll(newidx)));
   uIVData(adidx(ii)) = {(1:tmax)'};
end
   
uNumIV = length(uIVName);
if isempty(uIVName)
    ugrid = rgrid;
    varargout = cell(nin,1);
    return;
end

% Check for compatible IVData and construct indices that map into uIVName
varargout = cell(nin,1);
for i = 1:nin
    % Get i^th variable name/data
    iIVName = varargin{i}.IVName;
    iIVData = varargin{i}.IVData;
    iIVRB   = varargin{i}.IVRateBounds;
      
    % TODO 7/25/2012: This does not handle scalar expansion on the
    % array dimensions.  Need to revisit this.
    
    % Check compatibility of IVData
    % cellfun(@isequal,... is slower than direct call to isequal
    [~,idx] = ismember(iIVName,uIVName);
    adidx = strncmp(ArrayName,iIVName,length(ArrayName));
    ividx = ~adidx;
    inum = cellfun(@numel,iIVData(adidx));
    unum = cellfun(@numel,uIVData(idx(adidx)));
    
    if ~isequal(iIVData(ividx),uIVData(idx(ividx)))
        tmp = cellfun(@isequal,iIVData(ividx),uIVData(idx(ividx)));
        idx = find( ~tmp );
        var = iIVName{ividx(idx(1))};
        error(['Incompatible grid data for variable ' var]);
    end
    if ~all(inum==unum | inum==1 | unum == 1)
         error(['Incompatible array data']);
    end 
    if ~isequal(iIVRB(ividx,:),uIVRB(idx(ividx),:))
        error('Incompatible rate bound data')
    end

    % Start from uIVName(idx) = iIVName.
    % Construct C = [iIVName; setdiff(uIVName,iIVName)]
    %tmp = uIVName;
    %tmp(idx) = [];
    %C = [iIVName; tmp];
    
    % Construct idx2 such that uIVName(idx2) = C;    
    tmp = 1:uNumIV;
    tmp(idx) = [];
    idx2 = [idx; tmp(:)];
    
    % Construct idx3 such that uIVName = C(idx3)
    idx3(idx2) = 1:uNumIV;    
    varargout{i} = idx3;
end
    
% Construct grid with union of independent variables
ugrid = rgrid(uIVName,uIVData,uIVRB);


