function [E,keepidx] = lpvelimiv(R,IVNAME)
% LPVELIMIV  Eliminate independent variables.
%
% [E,KEEPIDX] = LPVELIMIV(R) eliminates all singleton independent 
% variables from the RGRID object R, i.e. the i^th independent variable of 
% R is removed if R.IVData{i} has length one.  KEEPIDX lists the variables  
% in R that are retained in E.
%
% [E,KEEPIDX] = LPVELIMIV(R,IVNAME) eliminates the independent variables 
% specified in IVNAME from the RGRID object R.   IVNAME can be a single 
% char array or a cell-array of CHARs.  KEEPIDX lists the variables in 
% R that are retained in E.
%
% See also: squeeze.

% Check # of input/output arguments
nin = nargin;
narginchk(1, 2)

% Get rgrid data
ivD = R.IVData;
ivN = R.IVName;
ivRB = R.IVRateBounds;
niv = numel(ivD);

% Find locations of array dimensions
%adidx =strncmp(R.ArrayName,R.IVName,length(R.ArrayName));
%idx = find(idx);

% Remove IVs
if isempty(R)
    E = R;
    keepidx = [];    
elseif nin==1
    % Remove singleton independent variables
    keepidx = find(R.LIVData>1);    
%     keepidx = (R.LIVData>1) | adidx;
%     keepidx = find(keepidx);
    E = rgrid( ivN(keepidx), ivD(keepidx) , ivRB(keepidx,:));    
else
    % Remove variables specified in IVNAME
    if isa(IVNAME,'char')
        IVNAME = {IVNAME};
    end
    
    % TODO PJS 4/4/2011: Vectorize for-loop below.
    idx = false(niv,1);
    for i=1:numel(IVNAME)
        idx = idx | strcmp(IVNAME{i},ivN);
    end
    keepidx = 1:niv;
    if numel(find(idx))>0
        keepidx(idx)=[];
        ivD(idx) = [];
        ivN(idx) = [];
        ivRB(idx,:) = [];
        
        % NOTE PJS 4/1/2011: Should we error if IVNAME contains variables
        % that are not contained in input rgrid R?
        
        %else
        %   error('Invalid IV name specification.');
    end
    E = rgrid(ivN,ivD,ivRB);
end

keepidx = keepidx(:)';