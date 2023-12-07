function B = lft2grid(A,varargin)
% Transform a TVREAL, PLFTMAT or PLFTSS object into a grid based LPV object
% PMAT, PSS, UPMAT, or UPSS. The transformation is performed by evaluating
% the PLFT object at a grid of parameter values.
%
% lft2grid(L) evaluates L at 10 values of each independent parameter,
% sampled uniformly from the range of each parameter.
%
% lft2grid(L,N) evaluates L at N values of each independent parameter,
% sampled uniformly from the range of each parameter.
%
% lft2grid(L,DOMAIN) evaluates L at each point containted in DOMAIN. 
% DOMAIN is an rgrid object that must contain the same independent 
% variables as L.


% Parse input to determine syntax case
nin = nargin;

if ~isa(A,'plft')
    error('Variable A needs to be a PLFT object')
end

if isa(A,'tvreal')
    A = plftmat(A);
end
        

if nin==1
    
    N = 10;
    B = lft2grid(A,N);
    return
    
elseif nin==2
    
    if isa(varargin{1} , 'double')
        N = varargin{1};
        TVNames = A.RateBounds(:,1);
        num = length(TVNames);
        TVGrid = cell(num,1);
        
        unc = A.Data.Uncertainty;
        
        for i = 1:num
            IVRange = unc.(TVNames{i}).Range;
            TVGrid{i} = linspace(IVRange(1),IVRange(2),N);
        end
        B = lft2grid(A,TVNames,TVGrid);
        return
        
    elseif isa(varargin{1},'rgrid')
        
        dom = varargin{1};
        TVNames = dom.IVName;
        num = length(TVNames);
        URateBounds = A.RateBounds;
        
        % Error check on consistency of rate bounds
        for i = 1:num
            
            iTVNames = TVNames{i};
            DRateBounds = dom.IVRateBounds(i,:);
            idx = find(strcmp(iTVNames,URateBounds(:,1)));
            
            if ~isempty(idx) && any(DRateBounds~=URateBounds{idx,2})
                error(['Conflicting name or rate bound information for parameter '...
                    iTVNames]);
            end
        end
        
        B = lft2grid(A,TVNames,dom.IVData);
        return
        
    else
        error('error')
    end
    
elseif nin>=3
    
    npairs = (nin-1)/2;
    if nin > 3 && floor(npairs) ~= ceil(npairs)
        error('lft2grid requires name and value pairs')
    end
    
    if isa(varargin{1},'char')
        C = reshape(varargin,[2,npairs]);
        B = lft2grid(A,C(1,:),C(2,:));
        return
    elseif ~iscellstr(varargin{1})
        error('lft2grid requires name and value pairs')
    end
end

TVNames = A.RateBounds(:,1);
NameCell = varargin{1};
GridCell = varargin{2};

% Identify input parameter names that don't exist in A as TVReals.
[~,idx] = setdiff( NameCell,TVNames );
NameCell(idx) = [];
GridCell(idx) = [];

% Make sure the input contains the names ofevery TVReal in A
[~,idxTV,idxNC] = intersect(TVNames,NameCell);
if length(idxTV)~=length(TVNames)
    error(' Input must include every TVReal that is contained in A')
end
RateBounds = cell2mat(A.RateBounds(idxTV,2));
NameCell = NameCell(idxNC);
GridCell = GridCell(idxNC);

%XXX determine how to handle the case when user requests a grid over a
%range of parameter values for which the LFT model TVReal is not defined.

GridCell = GridCell(:);
NameCell = NameCell(:);
if iscell(GridCell) && ndims(GridCell)==2 && size(GridCell,2)==1
    goodd = zeros(size(GridCell,1),1);
    for i=1:size(GridCell,1)
        GridCell{i} = GridCell{i}(:);
        if isa(GridCell{i},'double') && ndims(GridCell{i})==2 && ...
                size(GridCell{i},2)==1 && isreal(GridCell{i}) && ...
                all( diff(GridCell{i}) > 0 )
            goodd(i) = 1;
        end
    end
    if any(goodd==0)        
        error('Contents of GridCell should be single-column, sorted DOUBLE');
    elseif iscell(NameCell) && size(NameCell,1)~= size(GridCell,1)        
        error('Different length of IVName and IVData cell data');
    end
else    
    error('IVData should be a cell array');
end

% identify Array dimensions
sza = size(A);
% if length(sza)==2
%     nad = 2;
%     szAD = [1 1];
% elseif length(sza) == 3
%     nad = 2;
%     szAD = [sza(3),1];
if length(sza)==2
     nad = 0;
     szAD = [];
else
    szAD = sza(3:end);
    nad =length(szAD);
end

num = length(NameCell);
if isempty(GridCell) && isempty(NameCell)
    % Handle empty case. Ex) G = plftss(ss(-1,2,1,0)), lft2grid(G,rgrid)
    B = usubs(A.Data,'',[]);
    RG = rgrid;
else
    tmp = cell(size(GridCell));
    [tmp{:}] = ndgrid( GridCell{:} );
    if ~isempty(szAD)
        for i=1:num
            tmp{i} = repmat(tmp{i},[ones(1,num) szAD]);
            tmp{i} = shiftdim(tmp{i},num);        
        end
    end
    NV = [NameCell tmp]';        
    B = usubs(A.Data,NV{:});
    RG = rgrid(NameCell,GridCell,RateBounds);
end

% Arrange dimensions to be [data dim, IV dims, Array dims]
if isa(B,'double')
    B = permute(B,[1 2 2+nad+(1:num) 2+(1:nad)]);
else
    B = permute(B,[nad+(1:num) (1:nad) num+nad+1 num+nad+2]);
end


if isa(B,'double')    
    B = pmat(B,RG);
elseif isa(B,'ss')
    B = pss(B,RG);
elseif isa(B,'frd')
    B = pfrd(B,RG);
elseif isa(B,'umat')
    B = upmat(B,RG);
elseif isa(B,'uss')
    B = upss(B,RG);    
elseif isa(B,'ufrd')
    B = upfrd(B,RG);
end


%XXX What happens when user supplies a parameter names that doesn't
%exist as a TVReal in the model. [options: ignore or fan out array]

%XXX what happens when user supplies request for grid cell points outside 
% the domain of the TVReals in A.
