% PSTRUCT   Create a parameter-varying structure
%
%   M = PSTRUCT(Data,Domain) creates a parameter-varying structure defined on
%   an N-dimensional rectangular grid. Domain is an RGRID object that
%   specifies the N independent variables and the rectangular grid domain.
%   Data is an (N+2) dimensional structured array that specifies the
%   data. Data(:,:,i1,...,iN) is the value of the struct array evaluated
%   at the point Domain(i1,....,iN).
%   
%   Note: Use M.FieldName to access the field named 'FieldName' in M.
%         If possible, the content of the field is returned as an
%         object (e.g. pmat, pss), otherwise it is returned as a cell array.
%         
%
% See also: struct, rgrid, pgrid, pmat, pss, pfrd, upmat, upss, upfrd.

classdef pstruct
    % Class definition for parameter-varying structure
    
    properties (Dependent)
        % Public data with array dimensions trailing after indep vars
        Data = [];
        Domain;
        Parameter;
    end
    
    
    properties
        DomainPrivate = rgrid;
        DataPrivate = [];
    end
    
    
    methods
        % Constructor
        function obj = pstruct(DataPrivate,Domain)
            if nargin==1 && isa(DataPrivate,'pstruct');
                obj = DataPrivate;
            elseif nargin>0
                if nargin==1
                    Domain = rgrid;
                end
                
                % Set properties of a pstruct
                niv = Domain.NumIV;
                szMD = size(DataPrivate);
                szad = szMD(2+niv+1:end);
                obj.DomainPrivate = rgrid(Domain,szad);
                obj.DataPrivate =  DataPrivate;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if PSTRUCT object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = 'Domain must be an rgrid object.';
            elseif ~( isa(obj.DataPrivate,'struct') )
                pflag = 0;
                errstr = 'Data must be a struct array.';
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
        function Data = get.Data(obj)
            DomainPrivate = obj.DomainPrivate;
            ArrayName = DomainPrivate.ArrayName;
            idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
            % Ordering: Parameter vars, Array dims
            Data = permute(obj.DataPrivate,[1; 2; 2+find(~idx); 2+find(idx)]);
        end
        
        % get
        function Domain = get.Domain(obj)
            DomainPrivate = obj.DomainPrivate;
            IVNamePrivate = DomainPrivate.IVName;
            IVDataPrivate = DomainPrivate.IVData;
            IVRateBoundsPrivate = DomainPrivate.IVRateBounds;
            
            ArrayName = DomainPrivate.ArrayName;
            idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
            Domain = rgrid( IVNamePrivate(~idx), IVDataPrivate(~idx),...
                            IVRateBoundsPrivate(~idx,:));
        end
        
        
        % hasArray
        function out = hasArray(obj)
            % hasArray True for PSTRUCT objects with array dimensions.
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
          
                                    
            % Make header string 
            niv = obj.Domain.NumIV;
            sza = size(obj.Data);
            lea = length(sza);
            nad = obj.DomainPrivate.NumIV-niv;
            
             adcs = '';
             for i = [1 2 2+niv+1:numel(sza)]
                 adcs = [adcs int2str(sza(i)) 'x'];
             end
            s1 = [adcs(1:end-1) ' array of PSTRUCTs.'];
                        
            % Make display string
            if niv>0
                dispStr = display(obj.Domain);
                tmp = 'The PSTRUCT consists of the following blocks:';
                dispStr = char(tmp,dispStr(2:end,:));
            else
                dispStr = '';
            end
            s = char(s1,dispStr);
            disp(s);
            
            finam = fieldnames(obj.DataPrivate);
            disp(['The PSTRUCT has the following fields: ']);
            for kr = 1:length(finam)
                disp([' ' finam{kr}])
            end            
        end        
                
%         function display(obj)
%             
%             % Make Header
%             niv = obj.Domain.NumIV;
%             if niv==1
%                 ioStr = ['1 IV'];
%             else
%                 ioStr = [int2str(niv) ' IVs'];
%             end
%             hdr = [ioStr];
%             
%             
%             % Make header string    
%             sza = size(obj.Data);
%             lea = length(sza);
%             if lea == niv
%                 adcs = '1x1';
%             else
%                 adcs = '';
%                 for i=1:lea-niv-1
%                     adcs = [adcs int2str(sza(i)) 'x'];
%                 end
%                 adcs = [adcs int2str(sza(i+1))];
%             end             
%             disp([adcs ' array of PSTRUCTs: ' hdr]);
%             
%             
%             % Make IV display string
%             if niv>0
%                 dispStr = display(obj.Domain);
%             else
%                 dispStr = '';
%             end
%             disp(dispStr)
%             finam = fieldnames(obj.DataPrivate);
%             disp(['with fields: ']);
%             for kr = 1:length(finam)
%                 disp([' ' finam{kr}])
%             end
%             
%             %           % TODO - Fix Header
%             %          % Make Header
%             %          niv = obj.Domain.NumIV;
%             %          szm = size(obj.Data);
%             %          if szm(1)==1
%             %             ioStr = '1 Row, ';
%             %          else
%             %             ioStr = [int2str(szm(1)) ' Rows, '];
%             %          end
%             %          if szm(2)==1
%             %             ioStr = [ioStr '1 Column, '];
%             %          else
%             %             ioStr = [ioStr int2str(szm(2)) ' Columns, '];
%             %          end
%             %          if niv==1
%             %             ioStr = [ioStr '1 IV'];
%             %          else
%             %             ioStr = [ioStr int2str(niv) ' IVs'];
%             %          end
%             %          hdr = [ioStr];
%             %
%             %
%             %          % Make header string
%             %          [adcs,nad] = ad2char(obj.DomainPrivate);
%             %          if nad==0
%             %             disp(['PSTRUCT: ' hdr]);
%             %          else
%             %             disp([adcs ' array of PSTRUCTs: ' hdr]);
%             %          end
%             %
%             %          % Make display string
%             %          if niv>0
%             %             dispStr = display(obj.Domain);
%             %          else
%             %             dispStr = '';
%             %          end
%             %          disp(dispStr)
%             
% %             % TODO PJS 4/4/2011: Revisit. The code below displays the
% %             % matrix data if the the PSTRUCT has only 1 IV point.
% %             if prod( szm(3:end) ) == 1 && ~isempty(obj)
% %                 disp([inputname(1) ' = ']);
% %                 disp(obj.Data)
% %             end
%         end
        
        function out = fieldnames(obj)
            % FIELDNAMES   Get PSTRUCT fieldnames 
            %
            % B = FIELDNAMES(M) returns a cell array of strings B, naming 
            % the fields of the PSTRUCT M.
            %
            % See also: lpvgetfield.
            out = fieldnames(obj.DataPrivate);
        end
        
        
    end % end of methods
    
    % Numeric element-by-element utility methods
    methods
        
        
%         function out = ceil(mat)
%             % CEIL  Ceil for PSTRUCT objects.
%             %
%             % CEIL(M) rounds the elements of M to the nearest integer towards +infty
%             % at each point in the domain of M.
%             %
%             % See also: ceil, floor, round, fix.
%             out = mat;
%             out.DataPrivate = ceil(mat.DataPrivate);
%         end
        
        function out = conj(mat)
            % CONJ   Complex conjugate for PSTRUCT objects.
            %
            % CONJ(M) is the complex conjugate of the elements of M at each point in
            % the domain of M.
            %
            % See also: conj, real, imag, i, j.
            out = mat;
            out.DataPrivate = conj(out.DataPrivate);
        end
        
        function b = transpose(a)
            % TRANSPOSE  Non-conjugate transpose for PSTRUCT objects.
            %
            % TRANSPOSE(A) returns the non-conjugate transpose of A at
            % each point in the domain of A.
            %
            % See also: transpose, ctranspose, permute.
            
            out = mat;
            out.DataPrivate = transpose(out.DataPrivate);
        end
        
    end  % element-by-element numerical methods
    
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
                % Direct call to csize with nout = 2 put sthe product of the column
                % and IV dims into the second output. Instead SIZE with two
                % output args should return only the row/col dimensions.
                tmpout = csize(out,arg2,nargin,3);
                varargout = {tmpout{1}, tmpout{2}};
            else
                varargout = csize(out,arg2,nargin,nout);
            end
        end
    end
    
end % end of classdef



