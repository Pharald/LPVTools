% RGRID   Create a rectangular grid object
%   R = RGRID(IVName,IVData,IVRateBounds) creates a rectangular grid object
%   with independent variables and grid data specified by IVName, IVData
%   and IVRateBounds. For an N-dimensional rectangular grid, IVName is
%   an N-by-1 cell array of characters that specify the names of the
%   independent variables. IVData is a N-by-1 cell array of column vectors
%   that specify the grid data along each dimension. IVRateBounds is a
%   N-by-1 double array with two columns, where each row corresponds to a
%   parameter listed in IVNames, and each elements in the first column
%   specifies a lower rate bound and each element in the second column
%   specifies a upper rate bound. Each IVData{i} should be a vector
%   of sorted, real data.  If the RGRID contains only one independent
%   variable then IVName can be specified as a single char, IVData can be
%   specified as a single vector of sorted real data, and IVRateBounds can
%   be specified as a 1-by-2 row vector of real numbers.
% 
%   % EXAMPLE: (CUT/PASTE)
%   % Create an RGRID object with independent variable 'a', grid data 4:10.
%   r1 = rgrid('a',4:10)
% 
%   % Create an RGRID object with independent variable 'a', grid data 4:10
%   % and parameter rate bounds [-1 1].
%   r1 = rgrid('a',4:10,[-1,1])
% 
%   % Create a 2-dimensional RGRID object
%   r2 = rgrid( {'a', 'b'}, {linspace(-2,2,12), 1:5},[-1 1;-4 8] )
%   % Access independent variable name along first dimension
%   r2.IVName{1}
%   % Access independent variable rate bounds for parameter 'a'
%   r2.a.IVRateBounds
%   % Replace independent variable rate bounds for parameter 'b'
%   r2.b.IVRateBounds = [-10 10]
%
% See also: pgrid, pmat, pss, pfrd, upmat, upss, upfrd, pstruct.



% NOTE PJS:
% 1) I added a subsasgn to get around the set/isvalid issue in the
% constructor.  Specifically, the object is temporarily invalid in the
% constructor while we set the properties. We need to do the isvalid
% check only at the end.  However, we want to do an isvalid check
% whenever the user sets any field.  Route any user .-reference through
% subsasgn but let the constructor call the built-in set function.
%
% 2) I created a preliminary () subsref. For example if r is a 2-dim rgrid
% object then r(2,3) returns the (2,3) point in the grid. This makes it
% easy to grab the grid data points. However, the return format gets a bit
% messy in the general case.  If r is an N-dim grid then r(idx1,...,idxN)
% returns an N-by-m1-by-...-by-mN array where mi is the length of idxi.
% Also, I didn't create a () subsasgn because it seems like this should be
% set using a .-reference into IVData. Should either/both be implemented?

% TODO PJS 4/1/2011: Revisit ()-subsref

% NOTE 4/8/2011: Should the structure have the IVs as fields (similar
% to uncertainties.) rather than storing the IVs in cell arrays.

classdef rgrid
    % Class definition for rectangular grid object1
    
    properties (Dependent)
        % NumIV is the number of parameters (not countinng dummy
        % parameters associated with array dimensions) defined in object
        NumIV;
        
        % Parameter Gateway for user
        Parameter;
    end
    
    properties (Hidden)
        % Precompute information for fast computations later.
        
        % DIVData stores the diff of each entry of IVData
        DIVData = cell(0,1);
        
        % LIVData stores the length of each entry of IVData
        LIVData = zeros(0,1)
    end
    
    properties (Constant,Hidden)
        ArrayName = 'MusynLPVArrayName';
    end
    
    properties
        IVName = cell(0,1);
        IVData = cell(0,1);
        IVRateBounds = zeros(0,2);
    end
    
    methods
        % Constructor
        function obj = rgrid(IVName,IVData,IVRateBounds,varargin)
            if nargin==1 && isa(IVName,'rgrid');
                obj = IVName;
            elseif nargin==1 && isa(IVName,'pgrid');
                obj = rgrid(IVName.Name,IVName.GridData,IVName.RateBounds);
            elseif nargin==1 && isvector(IVName) && all(floor(IVName)==ceil(IVName))
                ADStart = 10000;
                ArrayDims = IVName;
                Ndim = numel(IVName);
                IVName = cell(Ndim,1);
                IVData = cell(Ndim,1);
                for i=1:Ndim
                    IVName{i} = [obj.ArrayName int2str(ADStart+i)];
                    IVData{i} = 1:ArrayDims(i);
                end
                IVRateBounds = repmat([-inf,inf],[Ndim,1]);
                obj = rgrid(IVName,IVData,IVRateBounds);
            elseif nargin>=2 && ( isa(IVData,'rgrid') || isa(IVData,'pgrid') )
                % Syntax: c=rgrid(a1,a2,a3,a4) where each ai is an rgrid
                % or pgrid.  This will combine all ai into a single rgrid.
                allgrid{1,1} = IVName;
                allgrid{2,1} = IVData;
                nin = nargin;
                if nin>=3
                    allgrid{3,1} = IVRateBounds;
                    allgrid = [allgrid; varargin(:)];
                end
                
                % Combine all data
                allIVName = cell(0,1);
                allIVData = cell(0,1);
                allIVRateBounds = [];
                for i=1:nin
                    % Lift all pgrids to rgrids
                    allgrid{i} = rgrid( allgrid{i} );
                    
                    allIVName = [allIVName; allgrid{i}.IVName];
                    allIVData = [allIVData; allgrid{i}.IVData];
                    allIVRateBounds = [allIVRateBounds; allgrid{i}.IVRateBounds];
                end
                
                obj = rgrid(allIVName,allIVData,allIVRateBounds);
                
            elseif (nargin==2)
                if isa(IVName,'rgrid')
                    % Call to add dummy IVs for array dims
                    Domain = IVName;
                    ArrayDims = IVData;
                    Ndim = length(ArrayDims);
                    
                    % Number start for Array Dimensions so that alphabetizing
                    % retains order by numbering
                    ADStart = 10000;
                    
                    % Check if Domain has array dimension with dummy name, if
                    % so find number of array dimensions and store in PTR.
                    idx=strncmp(obj.ArrayName,Domain.IVName,length(obj.ArrayName));
                    ptr = ADStart+sum(idx);
                    
                    tmpCell = cell(Ndim,3);
                    for i=1:Ndim
                        tmpCell{i,1} = [obj.ArrayName int2str(ptr+i)];
                        tmpCell{i,2} = 1:ArrayDims(i);
                    end
                    IVName = [Domain.IVName; tmpCell(:,1)];
                    IVData = [Domain.IVData; tmpCell(:,2)];
                    IVRateBounds = [Domain.IVRateBounds; repmat([-inf,inf],[Ndim,1])];
                    obj = rgrid(IVName,IVData,IVRateBounds);
                elseif iscell(IVName)
                    % Set properties of rgrid object
                    obj.IVName = IVName;
                    obj.IVData = IVData;
                    N = length(IVData);
                    obj.IVRateBounds = repmat([-inf,inf],[N,1]);
                elseif ischar(IVName)
                    % Simple syntax for a single var rgrid:
                    %    rgrid(CharString,DoubleVec,1-by-2 Double)
                    % XXX Handle mixed cases where some args are cells
                    % and other args are not.
                    obj.IVName = {IVName};
                    obj.IVData = {IVData};
                    obj.IVRateBounds = [-inf inf];
                end
                
                % Use isvalid to perform all error checking
                obj.IVName = obj.IVName(:);
                IVData = obj.IVData(:);
                for i = 1:length(IVData)
                    IVData{i} = IVData{i}(:);
                end
                obj.IVData = IVData;
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin == 3
                if iscell(IVName)
                    % Set properties of rgrid object
                    obj.IVName = IVName;
                    obj.IVData = IVData;
                    obj.IVRateBounds = IVRateBounds;
                elseif ischar(IVName)
                    % Simple syntax for a single var rgrid:
                    %    rgrid(CharString,DoubleVec,1-by-2 Double)
                    obj.IVName = {IVName};
                    obj.IVData = {IVData};
                    obj.IVRateBounds = IVRateBounds;
                end
                
                % Use isvalid to perform all error checking
                obj.IVName = obj.IVName(:);
                IVData = obj.IVData(:);
                for i = 1:length(IVData)
                    IVData{i} = IVData{i}(:);
                end
                obj.IVData = IVData;
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin>3
                error('Incorrect number of input arguments');
            end
        end
        
        % set: No validity checks
        function obj = set.IVName(obj,Value)
            if iscell(Value)
                % Convert cell array to column vector
                Value = Value(:);
            elseif ischar(Value)
                % Allow single IVName to be specified as a char
                Value = {Value};
            end
            obj.IVName = Value;
        end
        
        % set: No validity checks
        function obj = set.IVData(obj,Value)
            if iscell(Value)
                % Convert cell array and its contents to column vectors
                Value = Value(:);
                for i=1:length(Value)
                    Value{i} = Value{i}(:);
                end
            elseif isa(Value,'double');
                % Allow single IVData to be specified as a double
                Value = {Value(:)};
            end
            obj.IVData = Value;
            
            obj.DIVData = cellfun(@diff,Value,'UniformOutput',false);
            obj.LIVData = cellfun(@length,Value);
        end
        
        % set: No validity checks
        function obj = set.IVRateBounds(obj,Value)
            obj.IVRateBounds = Value;
        end
        
        % Set properties that are publicaly available to the user
        % (Other properties will be getable/setable but only the
        %  user should work through the Parameter gateway )
        function out = properties(obj)
            % PROPERTIES  Display property names for RGRID.
            
            p = {'Parameter'};
            if nargout==0
                sp = '    ';
                Np = length(p);
                for i=1:Np
                    p{i} = [sp p{i}];
                end
                p = [{'Properties for class RGRID:'}; p];
                disp(char(p));
            else
                out = p;
            end
        end
                
        
        function obj = subsasgn(obj,L,RHS)
        % SUBSASGN  Subscripted assignment for RGRID objects.
        %
        % See also: subsasgn, subsref.    
            switch L(1).type
                case '.'
                    % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
                    % The dot-reference here calls the set function
                    try
                        if length(L) == 1
                            tmp = RHS;
                        else
%                             idx = strcmp(L(1).subs,obj.IVName);
%                             if any(idx)
%                                 tmp = rgrid( obj.IVName{idx}, obj.IVData{idx}, ...
%                                     obj.IVRateBounds(idx,:));
%                             else
%                                 tmp = obj.(L(1).subs);
%                             end
%                             
%                             %tmp = obj.(L(1).subs);
                            tmp = subsref(obj,L(1));
                            tmp = subsasgn(tmp,L(2:end),RHS);
                        end
                        if strcmp(L(1).subs,'Parameter') 
                            % Parameter Gateway SUBSASGN
                            fn = fieldnames(tmp);
                            for i=1:numel(fn);
                                fni = fn{i};
                                if ~strncmp( fni, obj.ArrayName, length(obj.ArrayName) )
                                    [~,idx] = ismember(fni,obj.IVName);
                                    if idx>0
                                        GDi = tmp.(fni).GridData;
                                        obj.IVName{idx} = tmp.(fni).Name;
                                        obj.IVData{idx} = GDi;
                                        obj.IVRateBounds(idx,:) = tmp.(fni).RateBounds;
                                        obj.LIVData(idx) = numel(GDi);
                                        obj.DIVData{idx} = diff(GDi);
                                    end
                                end
                            end                            
                        else                        
                            idx = strcmp(L(1).subs,obj.IVName);
                            if any(idx)
                                obj.IVName{idx} = tmp.IVName{1};
                                obj.IVData{idx} = tmp.IVData{1};
                                obj.IVRateBounds(idx,:) = tmp.IVRateBounds;
                            else
                                obj.(L(1).subs) = tmp;
                            end
                        end
                    catch
                        error(lasterr);
                    end
                    
                    % Check validity
                    [pflag,errstr] = isvalid(obj);
                    if pflag==0
                        error(errstr);
                    end
                case '()'
                    error('() SUBSASGN not supported for rgrid objects.')
                case '{}'
                    error('Cell-like {} SUBSASGN not supported for rgrid objects.')
            end
        end
        
        % subsref
        function obj = subsref(obj,L)
        % SUBSREF  Subscripted reference for RGRID objects.
        %
        % See also: subsref, subsasgn.    
            switch L(1).type
                case '.'
                    if strcmp(L(1).subs,'Parameter')
                        % Parameter Gateway: 
                        % R.Parameter returns a struct with fields as pgrids 
                        % XXX Issue if IVName contains 'Parameter'
                        niv = length(obj.IVName);
                        m = struct;
                        for i=1:niv
                            n = obj.IVName{i};
                            if ~strncmp( n, obj.ArrayName, length(obj.ArrayName) )
                                m.(n) = pgrid( n, obj.IVData{i}, obj.IVRateBounds(i,:) );
                            end
                        end
                        obj = m;
                    else
                        %XXX object returns a cell array when call to
                        %obj.parametername.IVData -This is inconvinient
                        idx = strcmp(L(1).subs,obj.IVName);
                        if any(idx)
                            % .-Ref to Parameter: R.p returns as rgrid
                            obj = rgrid( obj.IVName{idx}, obj.IVData{idx}, ...
                                obj.IVRateBounds(idx,:));
                        else
                            % .-Ref to Ojbect property
                            obj = obj.(L(1).subs);
                        end
                    end
                case '()'
                    %  error('() SUBREF not supported for rgrid objects.')
                    niv = length(obj.IVName);
                    if niv>1
                        tmp = cell(niv,1);
                        [tmp{:}] = ndgrid(obj.IVData{:});
                    else
                        tmp = obj.IVData;
                    end
                    tmp = permute(cat(niv+1,tmp{:}),[niv+1 1:niv]);
                    L(1).subs = [':' L(1).subs];
                    obj = subsref(tmp,L(1));
                case '{}'
                    error('Cell-like {} SUBSREF not supported for rgrid objects.')
            end
            if length(L)>1
                obj = subsref(obj,L(2:end));
            end
        end
        
        % get
        function niv = get.NumIV(obj)
            niv = length(obj.IVName);
        end
        
        
        % Get array dimension as a character string
        function [cs,nad] = ad2char(obj)
            % ad2char Return array dimensions as a string
            %   [cs,nad] = ad2char(X) returns a string cs describing the
            %   number of array dimensions in X, and a scalar nad which
            %   describes the number of array dimensions in X.
            sza = size(obj);
            idx =strncmp(obj.ArrayName,obj.IVName,length(obj.ArrayName));
            idx = find(idx);
            nad = length(idx);
            if nad==1 && length(obj.IVData{idx})==1
                nad = 0;
            end
            
            cs = '';
            for i=1:nad
                j = idx(i);
                cs = [cs int2str(sza(j)) 'x'];
            end
            
            if nad==0
                cs = '1x1';
            elseif nad==1
                cs = [cs '1'];
            else
                cs = cs(1:end-1);
            end
        end
        
        % Display
        function s = display(obj)
            niv = length(obj.IVName);
            if niv ==0
                s = 'Empty rgrid';
            else
                s = 'RGRID with the following parameters:';
            end
            
            for i=1:niv
                Namei = obj.IVName{i};
                RBi = obj.IVRateBounds(i,:);
                GDi = obj.IVData{i};
                
                if ~strncmp( Namei, obj.ArrayName, length(obj.ArrayName) )
                    D = [Namei ': Gridded real, '];
                    
                    Npts = numel(GDi);
                    R = [GDi(1) GDi(end)];
                    R = sprintf('[%.3g,%.3g]',R(1),R(2));
                    D = [D int2str(Npts) ' points in ' R ', '];
                    RBi = sprintf('[%.3g,%.3g]',RBi(1),RBi(2));
                    D = [D 'rate bounds ' RBi '.'];
                    
                    s = char(s,['  ' D]);
                end
            end
            
            if nargout==0
                disp(s);
            end
        end
        
        % Size
        function [varargout] = size(obj,arg2)
            % SIZE   Size of an RGRID array
            % 
            %   S = SIZE(R), for an N-dimensional rgrid object R, returns
            %   a 1-by-N vector S with S(i) specifying the number of grid
            %   points in the i^th grid dimension.
            %
            %   S = SIZE(R,DIM) returns the number of grid points in the
            %   dimension specified by the scalar DIM. S=1 if DIM is
            %   greater than the number of grid dimensions of R.
            %
            % See also: size. 
            
            % NOTE 6/14/2012: Note this is a non-standard convention for SIZE. Rename
            % this function (RGRIDSIZE) and rewrite SIZE to have the standard
            % convention with respect to trailing singleton dimensions, etc.
            %niv = length(obj.IVData);
            out = obj.LIVData';
            
            %          if niv==1
            %             out = [obj.LIVData' 1];
            %          else
            %             out = obj.LIVData';
            %          end
            if nargin==1
                arg2 = nan;
            end
            varargout = csize(out,arg2,nargin,nargout);
            
        end
        
        % hasArray
        function out = hasArray(m)
            % hasArray True for rgrid objects with array dimensions.
            %
            %   hasArray(X) returns 1 if X has any array dimensions and
            %   0 otherwise.
            out = any(strncmp(m.ArrayName,m.IVName,length(m.ArrayName)));
        end
        
        % upidx
        % TODO AKP, July 25, 2012, change name to IVADidx
        function [IVidx,ADidx] = upidx(m)
            % UPIDX  Find locations of IV and array dimensions in data.
            % [IVidx,ADidx] = upidx(M) returns the indices of IV
            % dimensions in M as IVidx, and array dimensions of M as ADidx.
            adidx = strncmp(m.ArrayName,m.IVName,length(m.ArrayName));
            IVidx = find(~adidx);
            ADidx = find(adidx);
        end
        
        % isempty
        function out = isempty(obj)
            % ISEMPTY True for empty RGRID.
            %
            % ISEMPTY(M) returns 1 (true) if the RGRID M has no IVs,
            % and 0 (false) otherwise.
            %
            % See also: isempty.
            out = isempty(obj.IVName);
        end
        
        % rmprivateiv
        % TODO AKP, July 25, 2012, change name to removeAD
        function out = rmprivateiv(obj)
            % RMPRIVATEIV  Remove all array dimensions from rgrid object.
            idx = strncmp(obj.ArrayName,obj.IVName,length(obj.ArrayName));
            out = rgrid(obj.IVName(~idx),obj.IVData(~idx));
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if RGRID object is valid.
            
            errstr = [];
            pflag = 1;
            
            % IVName must be an N-by-1 cell array of character strings
            IVN = obj.IVName;
            if iscell(IVN) && ndims(IVN)==2 && size(IVN,2)==1
                goodn = zeros(size(IVN,1),1);
                for i=1:size(IVN,1)
                    if ischar(IVN{i}) && ndims(IVN{i})==2 && size(IVN{i},1)==1
                        goodn(i) = 1;
                    end
                end
                if any(goodn==0)
                    pflag = 0;
                    errstr = 'Contents of IVName should be single-row CHAR';
                end
            else
                pflag = 0;
                errstr = 'IVName should be a cell array';
            end
            
            % IVData must be an N-by-1 cell array of sorted, real vectors
            IVD = obj.IVData;
            if iscell(IVD) && ndims(IVD)==2 && size(IVD,2)==1
                goodd = zeros(size(IVD,1),1);
                for i=1:size(IVD,1)
                    if isa(IVD{i},'double') && ndims(IVD{i})==2 && ...
                            size(IVD{i},2)==1 && isreal(IVD{i}) && ...
                            all( diff(IVD{i}) > 0 )
                        goodd(i) = 1;
                    end
                end
                if any(goodd==0)
                    pflag = 0;
                    errstr = 'Contents of IVData should be single-column, sorted DOUBLE';
                elseif iscell(IVN) && size(IVN,1)~= size(IVD,1)
                    pflag = 0;
                    errstr = 'Different length of IVName and IVData cell data';
                end
            else
                pflag = 0;
                errstr = 'IVData should be a cell array';
            end
            
            % IVRateBounds must be an N-by-2 double array, where
            % IVRateBounds(i,2) > IVRateBounds(i,1) for all i
            IVRB = obj.IVRateBounds;
            if isa(IVRB,'double') && ndims(IVRB)==2 && size(IVRB,2)==2
                goodd = zeros(size(IVD,1),1);
                for i=1:size(IVD,1)
                    if isa(IVD{i},'double') && ndims(IVD{i})==2 && ...
                            size(IVD{i},2)==1 && isreal(IVD{i}) && ...
                            all( diff(IVD{i}) > 0 )
                        goodd(i) = 1;
                    end
                end
                if any( IVRB(:,2) <= IVRB(:,1) )
                    pflag = 0;
                    errstr = 'IVRateBounds data should increase by column';
                elseif iscell(IVN) && size(IVN,1)~= size(IVRB,1)
                    pflag = 0;
                    errstr = 'Different length of IVName and IVRateBounds data';
                end
            else
                pflag = 0;
                errstr = 'IVRateBounds should be a double array with two columns';
            end
            
            % Check for repeated names:
            U = unique(obj.IVName);
            if length(U)~=obj.NumIV
                error('Repeated IV name is not allowed in rgrid')
            end
            
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % XXX PJS: Replace vertcat with rgrid(a1,a2,..) syntax
        %         function out = vertcat(varargin)
        %             if nargin==1
        %                 out = varargin{1};
        %             else
        %                 out = rgrid;
        %                 v1 = varargin{1};
        %                 v2 = varargin{2};
        %                 out.IVName = [v1.IVName; v2.IVName];
        %                 out.IVData = [v1.IVData; v2.IVData];
        %                 out.IVRateBounds = [v1.IVRateBounds; v2.IVRateBounds];
        %                 if nargin>2
        %                     out = vertcat(out,varargin{3:end});
        %                 end
        %             end
        %         end
        
    end % end of methods
end % end of classdef


