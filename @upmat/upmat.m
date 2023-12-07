% UPMAT   Create an uncertain parameter-varying matrix
%
%   M = UPMAT(Data,Domain) creates an uncertain parameter-varying matrix
%   defined on an N-dimensional rectangular grid. Domain is an RGRID object
%   that specifies the N independent variables and the rectangular grid
%   domain.  Data is an (N+2) dimensional double array that specifies the
%   matrix data. Data(:,:,i1,...,iN) is the value of the matrix evaluated
%   at the point Domain(i1,....,iN).
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 matrix defined on a 1-dimensional grid
%   IVData = linspace(-2,2,20);
%   Domain = rgrid('a',IVData);
%   au = ureal('au',-2.3);
%   for i=1:length(IVData)
%       Data(1:2,1:2,i) = [1 au*IVData(i); IVData(i)^2 cos(IVData(i))];
%   end
%   M = upmat(Data,Domain)
%
% See also: rgrid, pgrid, pmat, pss, pfrd, upss, upfrd, pstruct.

% TODO PJS 6/14/2011: If I have an N-by-M matrix X and a N-by-M Domain D,
% then it seems like pmat(X,D) should return a 1-by-1 PMAT defined on
% the N-by-M domain. Currently the user must shift the dimensions of X
% so that it is a 1-by-1-by-N-by-M.

classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?pss,?uss,?ufrd,?umat}) upmat
    % Class definition for uncertain, parameter-varying matrix
    
    properties (Dependent)
        % Public data with array dimensions trailing after indep vars
        Data;
        Domain;
        NominalValue;
        Uncertainty;
        Parameter;
    end
    
    
    properties (Hidden = true)
        DomainPrivate = rgrid;
        DataPrivate = [];
    end
    
    % NOTE PJS 4/8/2011: M.DomainPrivate.Name should return the
    % IVData associated with the IV Name. Do this in SUBSREF/SUBSASGN
    % Domain and Data would be invalid IV Names.
    
    methods
        % Constructor
        function obj = upmat(a1,a2)
            if nargin==1
                if isa(a1,'upmat')
                    obj = a1;
                elseif isa(a1,'pmat')
                    obj.DomainPrivate = a1.DomainPrivate;
                    obj.DataPrivate =  umat(a1.DataPrivate);
                elseif ( isa(a1,'ureal') || isa(a1,'ucomplex') || ...
                        isa(a1,'ucomplexm') ||  isa(a1,'umat')  )
                    obj = upmat(umat(a1),rgrid);
                    %                 TODO - AH 10/10/12 - make decision about promoting double
                    %                 directly to upmat
                    %             elseif isa(a1,'double')
                    %                 obj.DataPrivate = umat(a1);
                else
                    error('Data must be a PMAT or UPMAT')
                end
            elseif nargin>1
                % Set properties of a UPMAT
                niv = a2.NumIV;
                szMD = size(a1);
                szad = szMD(2+niv+1:end);
                obj.DomainPrivate = rgrid(a2,szad);
                obj.DataPrivate =  a1;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin==0
                obj.DataPrivate = umat;
            else
                error('Invalid construction');
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if UPMAT object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = 'Domain must be an rgrid object.';
            elseif ~( isa(obj.DataPrivate,'umat') )
                pflag = 0;
                errstr = 'Data must be a UMAT object.';
            else
                szDomain = size(obj.DomainPrivate);
                szData = size(obj.DataPrivate);
                lszD = length(szDomain);
                lszM = length(szData);
                szDomain = [szDomain ones(1,lszM-2-lszD)];
                szData = [szData ones(1,lszD-lszM+2)];
                if ~( all(szDomain == szData(3:end)) || ...
                        (isempty(obj.DomainPrivate) && ndims(obj.DataPrivate)==2)   )
                    pflag = 0;
                    errstr = 'Dimensions of Domain and Data are incompatible.';
                end
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        
        % get
        function SA = get.Data(obj)
            [IVidx,ADidx] = upidx(obj.DomainPrivate);
            % Ordering: Parameter vars, Array dims
            tmp = [IVidx' ADidx'];
            if numel(tmp)==0 ||  numel(tmp)==1
                SA = obj.DataPrivate;
            else
                SA = permute(obj.DataPrivate,tmp);
            end
        end
        
        % get
        function out = get.Domain(obj)
            [IVidx,~] = upidx(obj.DomainPrivate);
            out = rgrid(obj.DomainPrivate.IVName(IVidx),...
                obj.DomainPrivate.IVData(IVidx),...
                obj.DomainPrivate.IVRateBounds(IVidx,:));
        end
        
        % XXX should user be able to set domain parameters by hand.
        % Set Domain
        function out = set.Domain(obj,Domain)
            InData = obj.Data;
            out = upmat(InData,Domain);
        end
        
        % get
        function out = get.NominalValue(obj)
            Domain = obj.Domain;
            Data = obj.Data.NominalValue;
            out = pmat(Data,Domain);
        end
        
        % get
        function out = get.Uncertainty(obj)
            out = obj.Data.Uncertainty;
        end
        
        
        % hasArray
        function out = hasArray(obj)
            % hasArray True for UPMAT objects with array dimensions.
            %   hasArray(X) returns 1 if X has any array dimensions and
            %   0 otherwise.
            out = hasArray(obj.DomainPrivate);
        end
        
        % ndims
        function out = ndims(obj)
            % NDIMS  Number of array dimensions.
            %
            % N = ndims(M) returns the number of array dimension in M.
            % 
            % See also: ndims.
            [~,ai] = upidx(obj.DomainPrivate);
            out = 2 + numel(ai);
            if out==3
                out = 4;
            end
        end
        
        
        % Display
        function s=display(obj)
            % Make Header
            niv = obj.Domain.NumIV;
            szm = size(obj.Data);
                                    
            % Make header string
            [adcs,nad] = ad2char(obj.DomainPrivate);
            if nad == 0
                s1 = 'UPMAT with ';
            else
                s1 = [adcs ' array of UPMATs with '];
            end
            s1 = [s1 sprintf('%d',szm(1)) ' rows and ' ...
                sprintf('%d',szm(2)) ' columns.'];
            
            % Make display string for parameters
            if niv>0
                dispStr = display(obj.Domain);
                tmp = 'The UPMAT consists of the following blocks:';
                dispStr = char(tmp,dispStr(2:end,:));
            else
                dispStr = '';
            end

            % Make display string for uncertainty
            [~,~,BlkStruct]=lftdata(obj.DataPrivate);
            Nblk = length(BlkStruct);
            for i=1:Nblk
                % Get short text description for each block
                % XXX PJS: Is getDescription a public function?
                Namei = BlkStruct(i).Name;
                Copiesi = BlkStruct(i).Occurrences;
                Blocki = obj.Data.Blocks.(Namei);
                D = getDescription(Blocki,Copiesi);
                dispStr = char(dispStr,['  ' D]);
            end
            
            s = char(s1,dispStr);
            if nargout==0
                disp(s)
                
%                 % TODO PJS 4/4/2011: Revisit. The code below displays the
%                 % matrix data if the the PMAT has only 1 IV point.
%                 if prod( szm(3:end) ) == 1 && ~isempty(obj)
%                     %             disp([inputname(1) ' = ']);
%                     %             disp(obj.Data)
%                     eval([inputname(1) ' = obj.Data'])
%                     
%                 end
            end
        end
        
        
%         function display(obj)
%             % Make Header
%             niv = obj.Domain.NumIV;
%             szm = size(obj.Data);
%             if szm(1)==1
%                 ioStr = '1 Row, ';
%             else
%                 ioStr = [int2str(szm(1)) ' Rows, '];
%             end
%             if szm(2)==1
%                 ioStr = [ioStr '1 Column, '];
%             else
%                 ioStr = [ioStr int2str(szm(2)) ' Columns, '];
%             end
%             if niv==1
%                 ioStr = [ioStr '1 IV'];
%             else
%                 ioStr = [ioStr int2str(niv) ' IVs'];
%             end
%             hdr = [ioStr];
%             
%             % Make header string
%             [adcs,nad] = ad2char(obj.DomainPrivate);
%             if nad==0
%                 disp(['UPMAT: ' hdr]);
%             else
%                 disp([adcs ' array of UPMATs: ' hdr]);
%             end
%             
%             % Make display string
%             if niv>0
%                 dispStr = display(obj.Domain);
%             else
%                 dispStr = '';
%             end
%             disp(dispStr)
%             
%             % TODO PJS 4/4/2011: Revisit. The code below displays the
%             % matrix data if the the PMAT has only 1 IV point.
%             if prod( szm(3:end) ) == 1 && ~isempty(obj)
%                 disp([inputname(1) ' = ']);
%                 disp(obj.Data)
%             end
%             
%             
%             
%             % GJB 26Jan12
%             % Handle Display of Uncertainty
%             zobjuss = obj.DataPrivate;
%             if isuncertain(zobjuss)
%                 warning off
%                 zobjussS = struct(zobjuss(:,:,1));
%                 warning on
%                 BD = getSummary(zobjussS.Data_.Blocks);
%                 nblk = numel(BD);
%                 if nblk==0
%                     disp('No uncertain block')
%                 else
%                     disp('The model uncertainty consists of the following blocks:')
%                     for ct=1:nblk
%                         disp(['  ' BD{ct}]);
%                     end
%                 end
%             end
%         end % end of display
        
    end % end of methods
    
    
    methods (Access = private)
        % privatesize
        function varargout = privatesize(m,arg2)
            out = size(m.DataPrivate);
            if length(out)==3
                out = [out 1];
            end
            
            if nargin==1
                arg2 = nan;
            end
            
            nout = nargout;
            if nout==2
                % Direct call to csize with nout = 2 puts the product of the
                % column and IV dims into the second output. Instead SIZE with
                % two output args should return only the row/col dimensions.
                tmpout = csize(out,arg2,nargin,3);
                varargout = {tmpout{1}, tmpout{2}};
            else
                varargout = csize(out,arg2,nargin,nout);
            end
        end
    end
    
end % end of classdef





