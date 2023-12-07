% PLFTMAT   Create a parameter-varying matrix in LFT framework.
%
% M = PLFT(Data,RateBounds) creates a parameter-varying matrix.
% Data is a UMAT. RateBounds is a N-by-2 cell array listing the
% rate bound information for each independent variable in the PLFTMAT.  
% RateBounds{i,1} is the character string name of the i-th independent 
% variable and RateBounds{i,2} is a sorted real vector of form [Low, High] 
% specifying its rate bounds.  RateBounds must only contain names of UREAL 
% objects that exist in Data and this indicates that the UREALs are actually 
% TVREALs representing the independent variables.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 matrix that depends on TVREAL a.
%   a = tvreal('a',[-2 2],[-1 1]);   
%   M = [1, a;a^2, -a]
%
% See also: tvreal, plftss.

% XXX PJS: Errors thrown by UMAT should be rethrown as PLFTMAT.
% XXX PJS: Binary operations simply call the corresponding UMAT operation.
%    We need to be careful that the UMAT binop does not perform an
%    AutoSimplify that is not allowed for time-varying parameters.

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,...
          ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) plftmat < plft
    % Class definition for parameter varying matrix in LFT framework
    
    properties (Dependent)
        Parameter = [];
        Uncertainty = [];
        NominalValue = [];
    end
   
    
    methods
        % Constructor
        function obj = plftmat(Data,RateBounds)
            
            obj.Data = umat;
            
            if nargin==1
                if isa(Data,'plftmat');
                    obj = Data;
                elseif ( isa(Data,'double') || isa(Data,'ureal') ...
                        || isa(Data,'ucomplex') || isa(Data,'ucomplexm') )
                    obj.Data = umat(Data);
                elseif isa(Data,'umat');
                    obj.Data = Data;
                elseif isa(Data,'tvreal');
                    obj.Data = umat(Data.Data);
                    obj.RateBounds = {Data.Data.Name Data.RateBounds};
                elseif ~isstatic(Data)
                    error('Cannot convert dynamic models to a plftmat.');
                else
                    error('Unsupported conversion to PLFTMAT.');
                end
            elseif nargin>1
                obj.Data = Data;
                obj.RateBounds = RateBounds;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            end
            
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if PLFTMAT object is valid.
            errstr = [];
            pflag = 1;
            
            RateBounds = obj.RateBounds;
            Data = obj.Data;
            Unc = fieldnames( Data.Uncertainty );
            
            Np = size(RateBounds,1);
            if length(unique(RateBounds(:,1))) ~= Np
                pflag = 0;
                errstr = 'List of names in RateBounds must be unique.';
            end
            
            for i=1:Np
                Namei = RateBounds{i,1};
                RBi = RateBounds{i,2};
                
                % XXX PJS Should rate bound interval contain 0?
                if ~( isa(RBi,'double') && length(RBi)==2 && diff(RBi)> 0 )
                    pflag = 0;
                    errstr = ['The value of "RateBounds" must be of the form '...
                        '[LOW,HIGH] with LOW < HIGH.'];
                end
                
                idx = find(  strcmp(Namei,Unc) );
                if isempty(idx) || ~isa( Data.Blocks.(Namei), 'ureal')
                    pflag = 0;
                    errstr = ['The parameter ' Namei ' does not exist ' ...
                        'as a UREAL in Data.'];
                end
            end
            
            if ~isa(Data,'umat')
                pflag = 0;
                errstr = 'Data must be a UMAT.';
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % Display
        % XXX PJS: Initial implementation. Needs more thought.
        function s = display(obj)
            szo = size(obj);
            [cs,nad] = ad2char(obj);
            if nad == 0
                s1 = 'PLFTMAT with ';
            else
                s1 = [cs ' array of PLFTMATs with '];
            end
            s1 = [s1 sprintf('%d',szo(1)) ' rows and ' ...
                sprintf('%d',szo(2)) ' columns.'];
            
            [~,~,BlkStruct]=lftdata(obj.Data);
            Nblk = length(BlkStruct);
            if Nblk ==0
                s2 = 'The PLFTMAT has no blocks.';
            else
                s2 = 'The PLFTMAT consists of the following blocks:';
            end
            s = char(s1,s2);
            
            
            RateBounds = obj.RateBounds;
            for i=1:Nblk
                % Get short text description for each block
                % XXX PJS: Is getDescription a public function?
                Namei = BlkStruct(i).Name;
                Copiesi = BlkStruct(i).Occurrences;
                Blocki = obj.Data.Blocks.(Namei);
                D = getDescription(Blocki,Copiesi);
                
                % Fix Description for TVREAL:
                % Remove nominal value and add rate bounds
                idx = find(  strcmp(Namei,RateBounds(:,1)) );
                if ~isempty(idx)
                    D1 = D;
                    RBLow = RateBounds{idx,2}(1);
                    RBHigh = RateBounds{idx,2}(2);
                    cidx = strfind(D,',');
                    
                    D = [Namei ': Time-varying real' D1(cidx(2):cidx(end))];
                    D = [D ' rate bounds = '];
                    D = [D sprintf('[%.3g,%.3g]',RBLow,RBHigh) ];
                    D = [D D1(cidx(end):end)];
                end                
                s = char(s,['  ' D]);
            end
            
            % XXX Add explanatory footer            
            % 'Type "get(N)" to see all properties and "N.Parameter" to 
            %    interact with the time-varying parameters.'

            if nargout==0
                disp(s)
            end
        end

        
        function out = properties(m)
            % PROPERTIES  Display property names for PLFTMAT.
            
            % Enforce hidden properties
            p = fieldnames(m);
            [~,idx]=ismember({'Data'; 'RateBounds'},p);
            p(idx) = [];
            if nargout==0
                sp = '    ';
                Np = length(p);
                for i=1:Np
                    p{i} = [sp p{i}];
                end
                p = [{'Properties for class plftmat:'}; p];
                disp(char(p));
            else
                out = p;
            end
        end        
        
        
        %         function [b,samples] = gridureal(a,varargin)
        %             [b,samples] = gridureal(plftmat(a),varargin{:});
        %         end
        %
        
        function out = ss(a,b,c,d,varargin)
            % Call PLFTSS constructor
            %
            % SYS = ss(A,B,C,D) creates a parameter-varying state-space 
            %     model. SYS will be a PLFTSS.
            %
            % SYS = ss(A) creates a parameter dependent static gain matrix.
            %
            % See also: ss, plftss.
            if nargin==1
                out = plftss(a);
            elseif nargin==4
                out = plftss(a,b,c,d);
            else
                out = plftss(a,b,c,d,varargin{:});
            end
        end
        
        function out = isuncertain(obj)
        %ISUNCERTAIN True if object is uncertain.
        %
        %  B = ISUNCERTAIN(A) is true if A is an uncertain object and false
        %  otherwise. The uncertain parameter-varying objects are UPMAT, 
        %  UPFRD, UPSS, PLFTMAT, and PLFTSS.
        %
        % See also: isuncertain.

            S.type = '.';
            S.subs = 'Uncertainty';
            fn = fieldnames(subsref(obj,S)); 
            out = true;
            if isempty(fn)
                out = false;
            end                
        end
        
    end  % methods
end % end of classdef



