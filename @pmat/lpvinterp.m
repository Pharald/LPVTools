function out = lpvinterp(m,varargin)
% LPVINTERP  Interpolate a PMAT
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,...) interpolates a PMAT A on
% the domain specified by the NAME/VALUES pairs. Each VALUE is a vector of 
% values at which to interpolate A along the domain direction of the
% corresponding NAME. If an independent variable of A is not listed
% in the inputs, then all values along this domain direction are retained
% in the output B.
%
% B = LPVINTERP(A,NAME,VALUES) is an alternative syntax. NAME is an N-by-1
% cell array of characters and VALUES is an N-by-1 cell array of values.
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD) includes a final 
% input argument called METHOD, which specifies the interpolation method to 
% be used. METHOD can be: 'nearest', 'linear', 'spline', or 'cubic'. The
% default is 'linear'. 
%
% See also: lpvsubs, lpvsplit, lpvsample.

% TODO PJS 4/4/2011: Revisit. Clean up code. 
% Also, what should be done for values outside of Domain.IVData?
% INTERP1 allows another input argument for to specify an EXTRAPVAL or
% to specify the use of extrapolation.

% Parse inputs
if nargin==3 && isa(varargin{1},'cell') && isa(varargin{2},'cell')
    % B = LPVINTERP(A,NAME,VALUES)
    newIvData = varargin{2};
    cnames = varargin{1};
    imeth = 'linear';
elseif nargin==4 && isa(varargin{1},'cell') && isa(varargin{2},'cell')
    % B = LPVINTERP(A,NAME,VALUES,METHOD)
    newIvData = varargin{2};
    cnames = varargin{1};
    imeth = varargin{3};
elseif floor(nargin/2)==ceil(nargin/2)
    % B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD)
    imeth = varargin{nargin-1};
    cnames = varargin(1:2:end-2);
    newIvData = varargin(2:2:end-1);
else
    % B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,....)
    imeth = 'linear';
    cnames = varargin(1:2:end-1);
    newIvData = varargin(2:2:end);
end

% TODO: GJB 24Jul12
% Need to add error checking to make sure requested data is within
% the bounds of the original data

% Check interpolation method
if ~any(strcmp(imeth,{'linear','cubic','nearest','spline','interpnlinear'}))
    error('Choose an allowable interpolation method')
end

% Get PMAT Data
szm = size(m.Data);
ivn = m.Domain.IVName;
ivd = m.Domain.IVData;
ivr = m.Domain.IVRateBounds;
divd = m.Domain.DIVData;
niv = m.Domain.NumIV;
nEachIV = m.Domain.LIVData;
szout = szm;

% Remove LPV registered ArrayNames before interpolation
idx =strncmp(rgrid.ArrayName,cnames,length(rgrid.ArrayName));
cnames(idx) = [];
newIvData(idx) = [];

% Make all the iv name strings have the same length:
ctmp = char([ivn;cnames(:)]);
civn = ctmp(1:niv,:);
ccnames = ctmp(niv+1:end,:);

% Find locations of NAMEs in the PMAT list of IVs
[~,keepoldiv] = setdiff(civn,ccnames,'rows');
[~,idxorig,idxnew] = intersect(civn,ccnames,'rows');

% Note to self: idxorig and keepoldiv should constitute a complete list of
% the indices of Domain.IVName.

% XXX 2/10/15 - AH - Add a check and generate explicit error when the user  
% inputs IV points outside of the domain.

% Extract Data
if niv==0
    out = m;
else
    % varcell stores the IV data points used for interpolation. 
    % varcell is a cell array ordered to match the IV order in Domain.IVName
    % Each IV has a corresponding cell in varcell, which contains the 
    % points that will be used for interpolation in that IV.
    varcell = cell(1,niv);
    % For IVs that are not to be interpolated, map IVData directly over to varcell.
    varcell(keepoldiv) = ivd(keepoldiv);
    % For IVs that are to be interpolated, map user supplied data over to varcell.
    varcell(idxorig) = newIvData(idxnew);
    
    % szout specifies the size of the data structure post-interpolation.    
    for i=1:length(idxorig);
        szout(idxorig(i)+2) = length(newIvData{idxnew(i)});
    end       
          
    % Define size of array dimensions and total number of AD points
    szAD = szm(3+niv:end);
    ADpts = prod(szAD);   
    % Define indices corresponding to matrix and IV data:
    MIVidx = 1:2+niv;
    
    switch imeth
        case 'linear'
            % XXX For linear interpolation, write our own fast linear
            % interpolation on matrices rather than calling for-loop.
            if niv==1
                VARCELL = varcell;
            else
                [VARCELL{1:niv}] = ndgrid(varcell{:});
            end
            szM = [numel(VARCELL{1}) niv];
            Mat = zeros(szM);
            for k=1:niv
                Mat(:,k) = reshape(VARCELL{k},[szM(1) 1]);
            end                                                
            
            % Loop though ADs. M will be arranged [:,:,IV,ADcolumn]
            M = zeros([szout(MIVidx),ADpts]); 
            template = zeros(szout(MIVidx));
            id = repmat({':'},1,numel(szout(MIVidx)));
            for k = 1:ADpts
                tmp = template;
                for i=1:szM(1)
                    tmp(:,:,i) = ndLinInterp(m.Data(id{:},k),niv,nEachIV,ivd,divd,Mat(i,:));
                end
                M(id{:},k) = tmp;
            end
            % Reshape from [:,:,IV,ADcolumn] to [:,:,IV,AD]
            M = reshape(M,szout);
            out = pmat(M,rgrid(ivn,varcell,ivr));
        otherwise 
            % NOTE PJS 4/7/2011: INTERPNLINEAR is an undocumented METHOD
            % to force a call to interp1 with 'linear' rather than calling
            % our own fast linear interpolation code.
            if isequal(imeth,'interpnlinear')
                imeth = 'linear';
            end
            
            % Will interpolate matrix element by element.
            % Define number of matrix elements
            tnpts = prod(szm([1 2])); 
            % Define data structure for results
            pout = zeros([tnpts szout(3:niv+2) ADpts]);
            % Define size of IV dimensions
            origivl = szm(3:2+niv); 
            
            % Collapse matrix dimensions into a column vector, and array 
            % dimensions into a single dim.            
            datapermute = reshape(m.Data,[tnpts origivl ADpts]); 
            id = repmat({':'},1,niv+1);
            if niv==1
                IVD = ivd;
                VARCELL = varcell;                                
                % Loop through AD and interpolate IVs:
                for k = 1:ADpts 
                    DataAtAD = datapermute(id{:},k);
                    holder = zeros([tnpts szout(3:niv+2)]);
                    for i=1:tnpts
                        % Loop through matrix elements, interpolating each:
                        tmp = interp1(IVD{1},DataAtAD(i,:),VARCELL{1},imeth);
                        holder(i,:) = tmp;
                    end
                    pout(id{:},k) = holder;
                end
            else
                [IVD{1:niv}] = ndgrid(ivd{:});
                [VARCELL{1:niv}] = ndgrid(varcell{:});
                % Loop through AD and interpolate IVs:
                for k = 1:ADpts 
                    DataAtAD = datapermute(id{:},k);
                    holder = zeros([tnpts szout(3:niv+2)]);
                    for i=1:tnpts
                        tmp = interpn(IVD{:},reshape(DataAtAD(i,:),origivl),VARCELL{:},imeth);
                        holder(i,:) = tmp;
                    end
                    pout(id{:},k) = holder;
                end
            end
            % Reshape from [:,:,IV,ADcolumn] to [:,:,IV,AD]
            pout = reshape(pout,szout);
            out = pmat(pout,rgrid(ivn,varcell,ivr));
    end
end





% The function below (ndLinInterp) is in LPVutil
% function R = ndLinInterp(M,nIV,N,IVData,dIVData,pValue)
% % M: nR-by-nC-by-[IVDims] DOUBLE array
% % nIV: scalar, double, number of IVs
% % N: nIV-by-1, length of each IVGrid
% % IVData: nIV-by-1 CELL, IVData{k} has k'th IVGrid, N(k)-by-1
% % dIVData: nIV-by-1 CELL, dIVData{k} has diff(k'th IVGrid), (N(k)-1)-by-1
% % pValue: nIV-by-1 DOUBLE, parameter value at which interpolation should
% % take place.
%
% szM = size(M);
% nVertices = 2^nIV;
% factor = ones(nVertices,1);
% idxcell = cell(1,nIV);
% for k=1:nIV
%    [i,alpha] = LOCALfindslotalpha(N(k),IVData{k},pValue(k),dIVData{k});
%    if ~isnan(i)
%       if i==N(k)
%          i = N(k)-1;
%          alpha = 1;
%       end
%       vecOnes = ones(2^(k-1),1);
%       idxcell{k} = [i i+1];
%       factor = factor.*kron(ones(2^(nIV-k),1),[(1-alpha)*vecOnes;alpha*vecOnes]);
%    end
% end
% R = M(:,:,idxcell{:});
% R = reshape(R,[szM(1)*szM(2) nVertices]);
% R = R*factor;
% R = reshape(R,[szM(1) szM(2)]);
%
%
% function [i,alpha] = LOCALfindslotalpha(N,vec,val,dvec)
% % N integer
% % vec 1-by-N (or N-by-1), sorted
% % val, scalar, vec(1) <= val <= vec(N)
% % dvec = diff(vec)
%
% i = max(find(val>=vec));
% if ~isempty(i)
%    if i<N
%       alpha = (val - vec(i))/dvec(i);
%    elseif val==vec(N)
%       alpha = 0;
%    else
%       i = NaN;
%       alpha = NaN;
%    end
% else
%    i = NaN;
%    alpha = NaN;
% end
