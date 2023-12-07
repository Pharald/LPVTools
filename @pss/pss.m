% PSS   Create a parameter-varying state-space model
%
%   S = PSS(Data,Domain) creates a parameter-varying state-space
%   model defined on an N-dimensional rectangular grid. Domain is an RGRID
%   object that  specifies the N independent variables and the rectangular
%   grid domain. Data is an N-dimensional state-space array that
%   specifies the state space data. Data(:,:,i1,...,iN) is the model
%   evaluated at the point Domain(i1,....,iN).
% 
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 1-by-1 state-space model defined on a 1-dimensional grid
%   IVData = linspace(2,20,10);
%   Domain = rgrid('a',IVData);
%   for i=1:length(IVData)
%       Data(1,1,i) = ss(-IVData(i),IVData(i),1,0);
%   end
%   S = pss(Data,Domain)
% 
%   % Overlay Bode plots at each independent variable
%   bode(S)
%
% See also: ss, pgrid, rgrid, pmat, pfrd, upmat, upss, upfrd.

% NOTE PJS 4/30/2011: Made TF/ZPK/FRD inferior to PSS. Also, modified
% constructor to lift TF/ZPK to a PSS.

% TODO PJS 5/1/2011: Should pfrd be inferior to pss? Or vice versa?

classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?ureal,?ucomplex,...
        ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) pss
    % Class definition for parameter-varying SS
    
    properties (Hidden = true)
        DomainPrivate = rgrid;
        DataPrivate = ss;
    end
    
    properties (Dependent)
        Domain;
        Data;
        Parameter;
    end
    
    methods
        % Constructor
        function obj = pss(DataPrivate,Domain,varargin)
            if nargin==1
                if isa(DataPrivate,'pss');
                    obj = DataPrivate;
                elseif ( isa(DataPrivate,'ss') || isa(DataPrivate,'tf') ...
                        || isa(DataPrivate,'zpk') ...
                        || isa(DataPrivate,'double') )
                    
                    obj = pss(ss(DataPrivate),rgrid);
                elseif isa(DataPrivate,'pmat')
                    obj = pss(ss(DataPrivate.DataPrivate),DataPrivate.DomainPrivate);
                else
                    error('Unsupported conversion to PSS.');
                end
            elseif nargin==2
                % Set properties of a pss
                niv = Domain.NumIV;
                szAD = size(DataPrivate);
                szad = szAD(2+niv+1:end);
                obj.DomainPrivate = rgrid(Domain,szad);
                obj.DataPrivate =  ss(DataPrivate);
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin>=4
                A = DataPrivate;
                B = Domain;
                C = varargin{1};
                D = varargin{2};
                % User is trying to pass through a pss object to copy its
                % properties into the new objects (i.e. InputName,etc.)         
                for i = 3:length(varargin)
                   if isa(varargin{i},'pss')
                      varargin{i} = varargin{i}.DataPrivate; 
                   end
                end
                [Aext,Bext,Cext,Dext] = domunion(A,B,C,D);
                obj.DataPrivate = ss(Aext.DataPrivate,Bext.DataPrivate,Cext.DataPrivate,Dext.DataPrivate,varargin{3:end});
                obj.DomainPrivate = Aext.DomainPrivate;
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if PSS object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = ['Domain must be an rgrid object.'];
            elseif ~isa(obj.DataPrivate,'ss')
                pflag = 0;
                errstr = ['Data must be a SS array.'];
            else
                % Data has r-by-c-by-Domain-by-ExtraDim
                szDomain = size(obj.DomainPrivate);
                szData = size(obj.DataPrivate);
                lszD = length(szDomain);
                lszS = length(szData);
                szDomain = [szDomain ones(1,lszS-2-lszD)];
                szData = [szData ones(1,lszD-lszS+2)];
                if ~( all(szDomain == szData(3:end)) || ...
                        (isempty(obj.DomainPrivate) && ndims(obj.DataPrivate)==2)   )
                    pflag = 0;
                    errstr = 'Dimensions of Domain and Data are incompatible.';
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
            
%             % Alternate fast method
%             try
%                 S = obj.DataPrivate;
%                 idx = find( isnan(S.d(1,1,:)) );
%                 S(:,:,idx) = [];
%                 tmp = S.a;
%             catch
%                 error('Varying state dimension is not allowed.');
%             end                        
            
            % Check for delays
            % NOTE PJS 4/30/2011: What is required to handle delays?
            %             if ~( isempty(S.InternalDelay) && all(S.InputDelay==0) && ...
            %                     all(S.OutputDelay==0) )
            if hasdelay(S)
                error('Input, Output, and Internal Delays not allowed.');
            end
            
            % TODO PJS 4/30/2011: Check for varying state names?
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end           
                
        % get
        function Data = get.Data(obj)
            DomainPrivate = obj.DomainPrivate;
            if ~isempty(DomainPrivate)
                ArrayName = DomainPrivate.ArrayName;
                idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
                % Ordering: Parameter vars, Array dims
                permidx = [find(~idx); find(idx)];
                if length(permidx) == 1
                    Data = obj.DataPrivate;
                else
                    Data = permute(obj.DataPrivate,permidx);
                end
            else
                Data = obj.DataPrivate;
            end
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
        

        % PJS TODO 1/10/2014: Add PROPERTIES function to make the lower
        % level USS properties accessible, e.g. state matrices? PSS
        % has such a properties function.
        
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
                s1 = 'PSS with ';
            else
                s1 = [adcs ' array of PSSs with '];
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
                tmp = 'The PSS consists of the following blocks:';
                dispStr = char(tmp,dispStr(2:end,:));
            else
                dispStr = '';
            end
            s = char(s1,dispStr);
            if nargout==0
                disp(s)
                
                % TODO PJS 4/4/2011: Revisit. The code below displays the
                % matrix data if the the PMAT has only 1 IV point.
                if prod( szm(3:end) ) == 1 && ~isempty(obj)
                    %             disp([inputname(1) ' = ']);
                    %             disp(obj.Data)
                    eval([inputname(1) ' = obj.Data'])
                    
                end
            end
        end
                
%         function display(obj)
%             % Make Header
%             niv = obj.Domain.NumIV;
%             S = obj.DataPrivate;
%             szm = size(S);
%             ns = size(S.a,1);
%             Ts = S.Ts;
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
%             %          if niv==1
%             %             ioStr = [ioStr '1 IV'];
%             %          else
%             %             ioStr = [ioStr int2str(niv) ' IVs'];
%             %          end
%             hdr = [ioStr];
%             
%             % Make header string
%             [adcs,nad] = ad2char(obj.DomainPrivate);
%             if nad==0
%                 disp(['PSS: ' hdr]);
%             else
%                 disp([adcs ' array of PSS: ' hdr]);
%             end
%             
%             
%             % Make display string
%             if niv>0
%                 dispStr = display(obj.Domain);
%             else
%                 dispStr = '';
%             end
%             disp(dispStr)
%             
%             % TODO PJS 4/30/2011: Revisit. The code below displays the
%             % SS data if the the PSS has only 1 IV point.
%             if prod( size(obj.Domain) ) == 1 && ~isempty(obj)
%                 eval([inputname(1) ' = S'])
%             end
%         end
        
    end % end of methods
    
    
    methods (Access = private)
        % privatesize
        function varargout = privatesize(s,arg2)
            out = size(s.DataPrivate);
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





