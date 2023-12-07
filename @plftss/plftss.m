% PLFTSS   Create a parameter-varying state-space model in LFT framework.
%
% M = PLFTSS(Data,RateBounds) creates a parameter-varying state-space model.
% Data is a USS. RateBounds is a N-by-2 cell array listing the
% rate bound information for each independent variable in the PLFTSS.  
% RateBounds{i,1} is the character string name of the i-th independent 
% variable and RateBounds{i,2} is a sorted real vector of form [Low, High] 
% specifying its rate bounds.  RateBounds must only contain names of UREAL 
% objects that exist in Data and this indicates that the UREALs are actually 
% TVREALs representing the independent variables.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 1-by-1 state-space model that depends on TVREAL b.
%   b = tvreal('b',[2 20]);
%   S = ss(-b,b,1,0)
%
% See also: tvreal, plftmat.

% XXX PJS: Errors thrown by USS should be rethrown as PLFTSS.
% XXX PJS: Binary operations simply call the corresponding USS operation.
%    We need to be careful that the USS binop does not perform an
%    AutoSimplify that is not allowed for time-varying parameters.

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,...
          ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) plftss < plft
    % Class definition for parameter varying matrix in LFT framework       

    properties (Dependent)        
        Parameter = [];
        Uncertainty = [];
        NominalValue = [];
        
        % Properties from USS (list as dependent to support tab complete)
        % Unsupported: 'SamplingGrid'; 'Blocks'
        a
        b
        c
        d
        e
        StateName
        StateUnit
        InputDelay
        OutputDelay
        Ts
        TimeUnit
        InputName
        InputUnit
        InputGroup
        OutputName
        OutputUnit
        OutputGroup
        Name
        Notes
        UserData
        InternalDelay
    end
    
    
    methods
        % Constructor
        function obj = plftss(Data,RateBounds,varargin)
            
            obj.Data = uss;
            
            if nargin==1
                RateBounds = cell(0,2);
                if isa(Data,'plftss');
                    obj = Data;
                elseif ( isa(Data,'ss') || isa(Data,'tf') ...
                        || isa(Data,'zpk') || isa(Data,'double') )
                    obj = plftss(uss(Data),RateBounds);
                elseif ( isuncertain(Data) && ~isa(Data,'ufrd') )
                    obj = plftss(uss(Data),RateBounds);
                elseif isa(Data,'tvreal')
                    RateBounds = {Data.Name Data.RateBounds};
                    obj = plftss(uss(Data.Data) , RateBounds);
                elseif isa(Data,'plftmat')
                    obj = plftss(uss(Data.Data),Data.RateBounds);
                else
                    error('Unsupported conversion to PSS.');
                end
            elseif nargin==2
                % Set properties of a PLFTSS
                obj.Data = Data;
                obj.RateBounds = RateBounds;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin>=4
                A = plftmat( Data );
                B = plftmat( RateBounds );
                C = plftmat( varargin{1} );
                D = plftmat( varargin{2} );
                
                % Find unique TVREAL list
                RateBounds = [A.RateBounds; B.RateBounds; ...
                    C.RateBounds; D.RateBounds];
                [~,idx1,idx2] = unique(RateBounds(:,1));
                RBUnique = RateBounds(idx1,:);
                
                % Check for coflicting rate bound info
                RBInterval1 = cell2mat( RateBounds(:,2) );
                RBInterval2 = cell2mat( RBUnique(idx2,2) );
                [i,~] = find( RBInterval1~=RBInterval2 );
                if ~isempty(i)
                    Name = RateBounds{i(1),1};
                    error(['The block with name "' Name '"  has two or more' ...
                        ' conflicting rate bound definitions.']);
                end
                
                % Set object data
                obj.Data = ss(A.Data,B.Data,C.Data,D.Data,varargin{3:end});
                obj.RateBounds = RBUnique;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            end
            
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if PLFTSS object is valid.
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
            
            if ~isa(Data,'uss')
                pflag = 0;
                errstr = 'Data must be a USS.';
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % Display
        % XXX PJS: Initial implementation. Needs further consideration.
        function s = display(obj)
            szo = size(obj);
            nx = order(obj);
            [cs,nad] = ad2char(obj);
            if nad == 0
                if isct(obj)
                    s1 = 'Continuous-time PLFTSS with ';
                else
                    s1 = 'Discrete-time PLFTSS with ';
                end
            else
                if isct(obj)
                    s1 = [cs ' array of continuous-time PLFTSSs with '];
                else
                    s1 = [cs ' array of discrete-time PLFTSSs with '];
                end
            end
            s1 = [s1 sprintf('%d',szo(1)) ' outputs, ' ...
                sprintf('%d',szo(2)) ' inputs, ' ...
                sprintf('%d',nx(1)) ' states.'];
            
            [~,~,BlkStruct]=lftdata(obj.Data);
            Nblk = length(BlkStruct);
            if Nblk ==0
                s2 = 'The model has no blocks.';
            else
                s2 = 'The model consists of the following blocks:';
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
                
                % Fix Description for TVREAL
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
            % PROPERTIES  Display property names for PLFTSS.
            
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
                p = [{'Properties for class plftss:'}; p];
                disp(char(p));
            else
                out = p;
            end
        end        
        
        
        % XXX PJS:
        % PLFTSS should have most (all?) of the same methods that exist
        % for a USS. I'll start with simple unary/binary operations.
        
                
        %         function [b,samples] = gridureal(a,varargin)
        %             [b,samples] = gridureal(pmatlft(a),varargin{:});
        %         end   
        
        function out = isuncertain(obj)
        %ISUNCERTAIN True if object is uncertain.
        %
        %  B = ISUNCERTAIN(A) is true if A is an uncertain object and false
        %  otherwise. The uncertain parameter-varying objects are UPMAT, UPFRD,
        %  UPSS, PLFTMAT, and PLFTSS.
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

        function b = issiso(a)
            % ISSISO  True for SISO PLFTSS objects.
            %
            % ISSISO(M) returns true if the PLFTSS M is single-input and 
            % single-output (SISO), and false otherwise.
            %  
            % See also: issiso.            
            b = issiso(a.Data);
        end                                

        function b = ss2ss(a,T)
            % SS2SS  Change of state coordinates for PLFTSS objects
            %
            % SYS = SS2SS(SYS,T) performs the similarity transformation z=Tx on the
            % state vector x of the PLFTSS model.  The transformation T must be a 
            % constant matrix.
            %  
            % see also: ss2ss, canon, balreal.
            b = a;
            b.Data = ss2ss(a.Data,T);
        end                                
        
        function b = repsys(a,varargin)
            % REPSYS   Replicate and tile for PLFTSS objects
            %   
            % B = REPSYS(A,M,N) creates a larger PLFTSS B consisting of an M-by-N tiling 
            % of copies of A. B = REPSYS(A,[M N]) also creates an M-by-N tiling of A.
            %
            % See also: repsys.   
            b = a;
            b.Data = repsys(a.Data,varargin{:});
        end                                
        
        function out = order(m)
            % ORDER  Compute the system order for PLFTSS
            %
            % N = ORDER(SYS) returns the order of the system SYS. For a 
            %      PLFTSS, N corresponds to the number of states in SYS.
            %
            % See also: order.
            out = order(m.Data);
        end
        
        function out = isct(m)
            % ISCT  True for continuous-time PLFTSS.
            %
            % ISCT(SYS) returns true if the PLFTSS SYS is continuous-time, 
            % and false otherwise.
            %
            % See also: isct, isdt, isstatic.
            out = isct(m.Data);
        end
        
        function out = isdt(m)
            % ISDT  True for discrete-time PLFTSS.
            %
            % ISDT(SYS) returns true if the PLFTSS SYS is discrete-time
            % and false otherwise. Returns true for empty systems and 
            % static gains.
            %  
            % See also: isdt, isct, isstatic.
            out = isdt(m.Data);
        end
        
        function out = isstatic(m)
            % ISSTATIC  Checks if PLFTSS is static or dynamic.
            %
            % ISSTATIC(SYS) returns TRUE if the PLFTSS SYS is static and FALSE if
            % SYS has dynamics.
            %
            % See also: isstatic, isct, isdt.
            out = isstatic(m.Data);
        end
               
        
        function [A,B,C,D,Ts] = ssdata(m)
            % SSDATA  Access to state-space data of a PLFTSS 
            %
            % [A,B,C,D] = SSDATA(SYS) returns the the A, B, C, D matrices of the PLFTSS.  
            % The state matrices are returned as PLFTMATs.
            % 
            % [A,B,C,D,Ts] = SSDATA(SYS) also returns the sampling time Ts.
            %
            % See also: ssdata.
            RateBounds = m.RateBounds;
            ssmat = cell(4,1);
            
            ssmat{1} = m.Data.a;
            ssmat{2} = m.Data.b;
            ssmat{3} = m.Data.c;
            ssmat{4} = m.Data.d;
            Ts = m.Data.Ts;
            
            % Remove tvreals that do not appear in ss matrix
            for i=1:4
                mati = ssmat{i};
                RateBoundsi = RateBounds; 
                if isuncertain(mati)
                    if isa(mati,'umat') || isa(mati,'uss') || ...
                            isa(mati,'ufrd')
                        U = fieldnames(mati.Uncertainty);
                    else
                        U = mati.Name;
                    end                    
                else
                    U = [];
                end
                [~,idx]=setdiff(RateBoundsi(:,1),U);
                RateBoundsi(idx,:)=[];
                ssmat{i} = plftmat(umat(mati),RateBoundsi);
            end
            A = ssmat{1};
            B = ssmat{2};
            C = ssmat{3};
            D = ssmat{4};            
        end                

    end  % methods
end % end of classdef



