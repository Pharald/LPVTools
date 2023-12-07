% UPSS   Create an uncertain parameter-varying state-space model
%
%   S = UPSS(Data,Domain) creates an uncertain parameter-varying
%   state-space model defined on an N-dimensional rectangular grid. Domain
%   is an RGRID object that  specifies the N independent variables and the
%   rectangular grid domain. Data is an N-dimensional uncertain
%   state-space array that specifies the uncertain state-space data. Note
%   that Data must contain the same uncertainty structure across the
%   array dimensions. Data(:,:,i1,...,iN) is the model evaluated at the 
%   point Domain(i1,....,iN).
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 1-by-1 uncertain, state-space model defined on a
%   % 1-dimensional grid
%   IVData = linspace(2,20,10);
%   Domain = rgrid('a',IVData);
%   theta = ureal('theta', 10,'Range',[8 12]);
%   for i=1:length(IVData)
%       Data(1,1,i) = ss(-IVData(i)*theta,IVData(i),1,0);
%   end
%   US = upss(Data,Domain)
% 
%   % Overlay Bode plots at each independent variable
%   bode(US);
%
% See also: uss, rgrid, pgrid, pmat, pss, pfrd, upmat, upfrd, pstruct.

% NOTE PJS 4/30/2011: Made TF/ZPK/FRD inferior to PSS. Also, modified
% constructor to lift TF/ZPK to a PSS.

% TODO PJS 5/1/2011: Should pfrd be inferior to pss? Or vice versa?

classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?pss,?uss,?ufrd,?umat}) upss
    % Class definition for uncertain, parameter-varying USS
    %
    
    properties (Hidden = true)
        DomainPrivate = rgrid;
        DataPrivate = uss;
    end
    
    properties (Dependent)
        Domain
        Data
        StateName;
        InputName;
        OutputName;
        NominalValue;
        Uncertainty;
        Parameter;
    end
    
    methods
        
        % Constructor
        function obj = upss(Data,Domain,varargin)
            if nargin==1
                if isa(Data,'upss');
                    obj = Data;
                elseif ( isa(Data,'uss') ||  isa(Data,'ultidyn') || ...
                        isa(Data,'udyn') || isa(Data,'ureal') || ...
                        isa(Data,'ucomplex') || isa(Data,'ucomplexm') || ...
                        isa(Data,'tf') || isa(Data,'zpk') || ...
                        isa(Data,'ss') ) || isa(Data,'double')
                    obj = upss(uss(Data),rgrid);
                elseif isa(Data,'upmat')
                    obj = upss(uss(Data.DataPrivate),Data.DomainPrivate);
                elseif isa(Data,'pss')
                    obj = upss(uss(Data.DataPrivate),Data.DomainPrivate);
                elseif isa(Data,'pmat')
                    obj = upss(uss(Data.DataPrivate),Data.DomainPrivate);
                else
                    error('Unsupported conversion to UPSS.');
                end
            elseif nargin==2
                if isa(Data,'ss') || isa(Data,'zpk') || isa(Data,'tf') ||...
                        isa(Data,'double')
                    Data = uss(Data);
                end
                % Set properties of a pss
                niv = Domain.NumIV;
                szP = size(Data);
                szad = szP(2+niv+1:end);
                obj.DomainPrivate = rgrid(Domain,szad);
                obj.DataPrivate = Data;
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin>=4
                A = Data;
                B = Domain;
                C = varargin{1};
                D = varargin{2};
                if nargin==5
                    Ts = varargin{3};
                else
                    Ts = 0;
                end
                [Aext,Bext,Cext,Dext] = domunion(A,B,C,D);
                obj.DataPrivate = ss(Aext.DataPrivate,Bext.DataPrivate,Cext.DataPrivate,Dext.DataPrivate,Ts);
                obj.DomainPrivate = Aext.DomainPrivate;
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if UPSS object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = ['Domain must be an rgrid object.'];
            elseif ~isa(obj.DataPrivate,'uss')
                pflag = 0;
                errstr = ['Data must be a USS array.'];
            else
                szDomain = size(obj.DomainPrivate);
                szData = size(obj.DataPrivate);
                lszD = length(szDomain);
                lszM = length(szData);
                szDomain = [szDomain ones(1,lszM-2-lszD)];
                szData = [szData ones(1,lszD-lszM+2)];
                if ~( all( szDomain == szData(3:end) ) || ...
                        ( isempty(obj.DomainPrivate) && ndims(obj.DataPrivate)==2)   )
                    pflag = 0;
                    errstr = ['Dimensions of Domain and Data are incompatible.'];
                end
            end
            
            %             % Check for varying state dimension
            %             S = obj.Data;
            %             a=ssdata(S,'cell');
            %             [nr,~]=cellfun(@size,a);
            %             if any( nr~=nr(1) )
            %                 error('Varying state dimension is not allowed.');
            %             end
            % Alternate fast method
            try
                S = obj.DataPrivate;
                tmp = S.a;
            catch
                error('Varying state dimension is not allowed.');
            end
            
            % Check for delays
            % NOTE PJS 4/30/2011: What is required to handle delays?
            %             if ~( isempty(S.InternalDelay) && all(S.InputDelay==0) && ...
            %                     all(S.OutputDelay==0) )
            % GJB 24JAN12, eliminated just to get things working at school under
            % R2010a
            %             if hasdelay(S)
            %                 error('Input, Output, and Internal Delays not allowed.');
            %             end
            
            % TODO PJS 4/30/2011: Check for varying state names?
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % StateName
        function out = get.StateName(obj)
            out = obj.DataPrivate.StateName;
        end
        function obj = set.StateName(obj,Val)
            obj.DataPrivate.StateName = Val;
        end
        
        % InputName
        function out = get.InputName(obj)
            out = obj.DataPrivate.InputName;
        end
        function obj = set.InputName(obj,Val)
            obj.DataPrivate.InputName = Val;
        end
        
        % Domain
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
            out = upss(InData,Domain);
        end
        
        % NominalValie
        function out = get.NominalValue(obj)
            out = pss(obj.DataPrivate.NominalValue,obj.DomainPrivate);
        end
        
        % get
        function out = get.Uncertainty(obj)
            out = obj.Data.Uncertainty;
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
        
        % OutputName
        function out = get.OutputName(obj)
            out = obj.DataPrivate.OutputName;
        end
        function obj = set.OutputName(obj,Val)
            obj.DataPrivate.OutputName = Val;
        end
        
        
% TODO AH 9/5/14 - Decide on desired functionality. Listing the
% properties of the underlying USS implies that the user can access
% them. However, this functionality is not implemented yet. 
% See PLFTSS for example of how the properties of the underlying
% USS object can be accessed. In that case its straight forward 
% because the PLFTSS object is built on top of a single USS object.
% How should this work for a grid-based object? Convert numeric
% arrays to UPMATs etc... 
%
%         function out = properties(m)
%             % PROPERTIES  Display property names for UPSS.
%             
%             p = {'Domain';'Data'};
%             p = [properties(m.DataPrivate); p];
%             if nargout==0
%                 sp = '    ';
%                 Np = length(p);
%                 for i=1:Np
%                     p{i} = [sp p{i}];
%                 end
%                 p = [{'Properties for class upss:'}; p];
%                 disp(char(p));
%             else
%                 out = p;
%             end
%         end
        
        
        % Display
        function s=display(obj)
            % Make Header
            niv = obj.Domain.NumIV;
            S = obj.DataPrivate;
            szm = size(S);
            ns = size(S.a,1);
            Ts = S.Ts;
                                    
            % Make header string
            [adcs,nad] = ad2char(obj.DomainPrivate);
            if nad == 0
                s1 = 'UPSS with ';
            else
                s1 = [adcs ' array of UPSSs with '];
            end
            s1 = [s1 int2str(ns) ' States, '];
            s1 = [s1 int2str(szm(1)) ' Outputs, '];
            s1 = [s1 int2str(szm(2)) ' Inputs, '];
            if Ts==0
                s1 = [s1 'Continuous System.'];
            else
                s1 = [s1 'Discrete System, Ts = ' num2str(Ts) '.'];
            end
            
            % Make display string
            if niv>0
                dispStr = display(obj.Domain);
                tmp = 'The UPSS consists of the following blocks:';
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
%             szm = size(obj.DataPrivate);
%             ns = size(obj.DataPrivate.a,1);
%             Ts = obj.DataPrivate.Ts;
%             if ns==1
%                 ioStr = ['1 State, '];
%             else
%                 ioStr = [int2str(ns) ' States, '];
%             end
%             if szm(1)==1
%                 ioStr = [ioStr '1 Output, '];
%             else
%                 ioStr = [ioStr int2str(szm(1)) ' Outputs, '];
%             end
%             if szm(2)==1
%                 ioStr = [ioStr '1 Input, '];
%             else
%                 ioStr = [ioStr int2str(szm(2)) ' Inputs, '];
%             end
%             if Ts==0
%                 ioStr = [ioStr 'Continuous System, '];
%                 %elseif Ts==-1
%                 %    ioStr = [ioStr 'Discrete System, unspecified Ts, '];
%             else
%                 ioStr = [ioStr 'Discrete System, Ts = ' num2str(Ts) ', '];
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
%                 disp(['UPSS: ' hdr]);
%             else
%                 disp([adcs ' array of UPSSs: ' hdr]);
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
%             %             % TODO PJS 4/30/2011: Revisit. The code below displays the
%             %             % SS data if the the PSS has only 1 IV point.
%             %             if prod( szm(3:end) ) == 1 && ~isempty(obj)
%             %                 disp([inputname(1) ' = ']);
%             %                 S
%             %                 %disp(S);
%             %             end
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
%         end
        
        
    end % end of methods
    
    methods (Access = private)
        % privatesize
        function varargout = privatesize(P,arg2)
            out = size(P.DataPrivate);
            if length(out)==3
                out = [out 1];
            end
            if nargin==1
                arg2 = nan;
            end
            nout = nargout;
            if nout==2
                tmpout = csize(out,arg2,nargin,3);
                varargout = {tmpout{1}, tmpout{2}};
            else
                varargout = csize(out,arg2,nargin,nout);
            end
        end
    end
    
    
end % end of classdef





