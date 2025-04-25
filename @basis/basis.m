% BASIS   Create a basis function object
%
%   Basis functions are needed for rate-bounded LPV analysis and synthesis.
%   These are functions of the independent variables present in the system
%   being analyzed.  Basis functions are specified as scalar PMATs.  All
%   partial derivatives must also be provided by the user.
%
%   b = BASIS(F,NAME1,PARTIAL1,NAME2,PARTIAL2,...) creates a basis function
%   object with the function F, which is a scalar PMAT. If there are N 
%   independent variables in F, then 2*N additional arguments are specified, 
%   which are pairs of (1) an independent variable name (as CHAR), 
%   and (2) the corresponding partial derivative (as PMAT) of F with  
%   respect to that independent variable.
%
%   BASIS(F,DERIVATIVE) is the same as BASIS(F,NAME1,DERIVATIVE) if F has
%   only one independent variable.
% 
%   % EXAMPLE: (CUT/PASTE)
%   theta = pgrid('theta',0:linspace(0,2*pi,20));
%   psi = pgrid('psi',linspace(0,2*pi,10));
%   F = cos(theta)*sin(2*psi);
%   pTheta = -sin(theta)*sin(2*psi);
%   pPsi = 2*cos(theta)*cos(2*psi);
%   B = basis(F,'theta',pTheta,'psi',pPsi)


% XXX Address extraneous data in the underlying pmat data which is slowing
% down computation compared to original non-class based structure based 
% basis function representation.

classdef basis
    % Class definition for basis function object

    properties
        BasisFunction = pmat(zeros(0,1));
        Partials = pmat;
        IVName = cell(0,1);
    end
    
    methods
        % Constructor
        function obj = basis(BasisFunction,varargin)
            nvar = numel(varargin)/2;
            if nargin==0                
                return
            elseif nargin==1 && isa(BasisFunction,'basis');
                obj = BasisFunction;
            elseif nargin==2
                if isa(BasisFunction, 'plftmat')
                    % HP 25/04/25: This is a rudimentary first shot. For
                    % lft based basis function it can only be a single
                    % plftmat as basis function and plftmat for its
                    % partial. No other functionalties work at the moment
                    % for lft basis functions
                    obj.BasisFunction = plftmat(BasisFunction);
                    nvar = length(fieldnames(BasisFunction.Parameter));
                    obj.Partials = plftmat(varargin{1});
                else
                    BasisFunction = pmat(BasisFunction);
                    nvar = BasisFunction.Domain.NumIV;
                    Partials = pmat(varargin{1});

                    if nvar == 0 && numel(Partials.Data)==1
                        if Partials.Data == 0
                            obj.BasisFunction = BasisFunction;
                            obj.Partials = pmat(zeros(1,0));
                            % XXX handles calls like: b1 = basis(7,0)
                        else
                            error('Partial associated with a constant BasisFunction must be zero')
                        end

                    elseif nvar ==1
                        [BFext,Pext] = domunion(BasisFunction,Partials);
                        obj.BasisFunction = BFext;
                        obj.Partials = Pext;
                        obj.IVName = BFext.Domain.IVName;
                    else
                        error('BasisFunction must have at most one IV for 2 input argument call')
                    end
                end
                
                
            elseif nargin ==3 && ~ischar(varargin{1})
                [BFext,Pext] = domunion(BasisFunction,varargin{1});
                obj.BasisFunction = BFext;
                obj.Partials = Pext;
                obj.IVName = varargin{2};
            elseif floor(nvar)==ceil(nvar) && nargin>1
                BasisFunction = pmat(BasisFunction);
                BasisFunction = BasisFunction(:);
                nbas = numel(BasisFunction);
                Partials = pmat(zeros(nbas,nvar));
                IVName = cell(nvar,1);
                for i = 1:nvar
                    IVName(i) = varargin(2*i-1);
                    tmp = pmat(varargin{2*i});
                    Partials(:,i) = tmp(:);
                end
                [BFext,Pext] = domunion(BasisFunction,Partials);
                obj.BasisFunction = BFext;
                obj.Partials = Pext;
                obj.IVName = IVName;
            else
                error('Incorrect number of input arguments.');
            end
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end                
        
        % subsasgn
        function obj = subsasgn(obj,L,RHS)
        % SUBSASGN  Subscripted assignment for BASIS objects.
        %
        % See also: subsasgn, subsref.    
            switch L(1).type
                case '.'
                    error('. SUBSASGN not supported for basis objects.')
                case '()'
                    % Define m and RHS on same domain
                    L1 = L(1);
                    if isempty(obj)
                        obj = basis;
                    end
                    z = [obj;RHS];
                    nb = size(obj.BasisFunction,1);                    
                    BF = subsasgn(z.BasisFunction(1:nb),L1,z.BasisFunction(nb+1:end));
                    L1.subs = [L1.subs,':'];
                    P  = subsasgn(z.Partials(1:nb,:),L1,z.Partials(nb+1:end,:));
                    IVN = z.IVName;
                    obj = basis(BF,P,IVN);
                case '{}'
                    error('Cell-like {} SUBSASGN not supported for basis objects.')
            end
            % Check validity
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        
        function obj = subsref(obj,L)
        % SUBSREF  Subscripted reference for BASIS objects.
        %
        % See also: subsref, subsasgn.    
            switch L(1).type
                case '.'
                    obj = obj.(L(1).subs);
                case '()'
                    BF = obj.BasisFunction(L(1).subs{:});
                    P = obj.Partials(L(1).subs{:},:);
                    IVN = obj.IVName;
                    obj = basis(BF,P,IVN);
                case '{}'
                    error('Cell-like {} SUBSREF not supported for basis objects.')
            end
            if length(L)>1
                obj = subsref(obj,L(2:end));
            end
        end
        
        % Display
        function display(obj)
            dispStr = 'Empty basis';
            nbasis = size(obj.BasisFunction,1);
            nvar = size(obj.IVName,1);
            
            % check for plftmat. If it is in plftmat just diplay the
            % plftmat. 
            if isa(obj.BasisFunction,'plftmat')
                obj.BasisFunction

            else


                % Describe number of IVs
                if nvar==1
                    IVStr = ['1 PGRID'];
                else
                    IVStr = [int2str(nvar) ' PGRIDs'];
                end

                if nbasis>0
                    dispStr = ['BASIS: ' int2str(nbasis) ' basis functions '];
                    dispStr = [dispStr 'and ' int2str(nvar) ' partial derivatives'];
                    dispStr = [dispStr ' with respect to ' IVStr];
                end
                disp(dispStr);

                % Display individual IVs
                Dom = obj.BasisFunction.Domain;
                dispStr = display(Dom);
                if nvar>0
                    tmp = 'The BASIS object consists of the following blocks:';
                    for i=2:size(dispStr,1)
                        stri = dispStr(i,:);
                        idx = strfind(stri,', rate bounds');
                        tmp = char(tmp,stri(1:idx-1));
                    end
                    dispStr = tmp;
                    %                 dispStr = char(tmp,dispStr(2:end,:));
                else
                    dispStr = '';
                end
                disp(dispStr)
            end
        end
        
        % Size
        function [varargout] = size(obj,arg2)
            % SIZE   Size of a BASIS object.
            %
            % For an M-by-N BASIS object A with P independent variables, 
            % S = SIZE(A) returns a 1-by-2 row vector. S(1) and S(2) are the 
            % row and column dimensions of A (S(1) = M and S(2) = N). 
            %
            % See also: size, length.    
            
            if nargin==1
                out = size(obj.BasisFunction);
            else
                out = size(obj.BasisFunction,arg2);
            end
            if nargout==2
                if nargin ==2
                    error('Too many output arguments.');
                else
                    varargout = {out(1),out(2)};
                end
            else
                varargout = {[out]};
            end
            
       % OLD CODE
%             if nargin==1
%                 [out1,out2] = size(obj.BasisFunction);
%             else
%                 % XXX-potential bug with b7(2) = b6b
%                 [out1,out2] = size(obj.BasisFunction,arg2);
%             end
%             if nargout==2
%                 varargout = {out1,out2};
%             else
%                 varargout = {[out1,out2]};
%             end
        end

        % length
        function out = length(b)
        % LENGTH   Length of a BASIS object.
        %
        % L = LENGTH(A) returns the maximum of the row and column dimension  
        % of A. This number corresponds to the number of basis functions in A.
        %
        % See also: length, size.    
            szb = size(b);
            out = max(szb(1:2));
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
        % ISVALID Determine if BASIS object is valid.

        % 25/04/25 HP: no check for plftmat yet
            errstr = [];
            pflag = 1;
            
            % IVName must be an N-by-1 cell array of character strings
            IVN = obj.IVName;
            if iscell(IVN) && ndims(IVN)==2 && size(IVN,2)==1
                nvar = size(IVN,1);
                goodn = zeros(nvar,1);
                for i=1:size(IVN,1)
                    if ischar(IVN{i}) && ndims(IVN{i})==2 && size(IVN{i},1)==1
                        goodn(i) = 1;
                    end
                end
                if any(goodn==0)
                    pflag = 0;
                    errstr = 'Contents of IVName should be single-row CHAR';
                end               
                % XXX Check that IVName doesn't contain any repeated IV
                % names               
            else
                pflag = 0;
                errstr = 'IVName should be a cell array';
            end
            

            % BasisFunction must be an nbasis-by-1 PMAT
            BF = obj.BasisFunction;

            if ~isa(BF,'plftmat')
                if isa(BF,'pmat') && size(BF,2)==1
                    nbasis = size(BF,1);
                else
                    pflag = 0;
                    errstr = 'BasisFunction should be an Nbasis-by-1 PMAT';
                end
                % All variables in BF must appear in IVN
                if ~isempty( setdiff(BF.Domain.IVName,IVN) )
                    pflag = 0;
                    errstr = 'Variables in BasisFunction must be a subset of Variables in IVName';
                end
            end
            % Partials must be an nbasis-by-nvar PMAT
            P = obj.Partials;

            if ~isa(P,'plftmat')

                if isa(P,'pmat') && size(P,1)==nbasis && size(P,2)==nvar
                    % All variables in Partial must appear in BF
                    if ~isempty( setdiff(P.Domain.IVName,BF.Domain.IVName) )
                        pflag = 0;
                        errstr = 'Variables in Partials must be a subset of Variables in BasisFunctions';
                    end

                    % Any variables in IVN that don't appear in BF must have
                    % zero partials
                    lidx = ismember(IVN,BF.Domain.IVName);
                    idx = find(~lidx);
                    if ~isempty(idx)
                        tmp = P(:,idx);
                        if any(tmp.Data)
                            pflag = 0;
                            errstr = 'Any variables in IVName that do not appear in BasisFunctions must have zero partials.';
                        end
                    end
                else
                    pflag = 0;
                    errstr = 'Partials should be a Nbasis-by-Nvar PMAT';
                end
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        function out = horzcat(varargin)
            out = vertcat(varargin{:});
        end
            
        function out = vertcat(varargin)
        % VERTCAT  Vertical concatenation of BASIS objects.
        %
        % S = VERTCAT(S1,S2,...) performs a concatenation operation of
        % S = [S1; S2; ...]  at each point in the combined domains 
        % of S1, S2, ...   Note that the concatenation operation will 
        % expand each basis function and its partial derivatives to the 
        % common domain of S1, S2, .... 
        %
        % See also: vertcat.
        
            if nargin==1
                out = varargin{1};
            else
                [v1ext,v2ext] = domunion(varargin{1},varargin{2});
                Partials = [v1ext.Partials;v2ext.Partials];
                IVName = v1ext.IVName;
                BasisFunction = [v1ext.BasisFunction; v2ext.BasisFunction];
                out = basis(BasisFunction,Partials,IVName);
                if nargin>2
                    out = vertcat(out,varargin{3:end});
                end
            end
        end
        
        % XXX - Determine if it makes sense to add more mathematical
        % operations.
        
        function out = uplus(obj)
        % UPLUS   Unary plus for BASIS objects.
        %
        % +M is the unary plus of M at each point in the domain of M.
        %
        % See also: uplus    
            out = obj;
        end
        
        function out = uminus(obj)
        % UMINUS   Unary minus for BASIS objects.
        %
        % -M is the unary minus of M at each point in the domain of M.
        %
        % See also: uminus    
            out = obj;
            out.BasisFunction = -out.BasisFunction;
            out.Partials = -out.Partials;
        end
        
        function out = plus(obj1,obj2)
        % PLUS  Plus for BASIS objects
        %
        % PLUS(A,B) is the result of A+B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: plus.    
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF = obj1ext.BasisFunction+obj2ext.BasisFunction;
            P = obj1ext.Partials+obj2ext.Partials;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = minus(obj1,obj2)
        % MINUS  Minus for BASIS objects
        %
        % MINUS(A,B) is the result of A-B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: minus.    
            out = plus(obj1,-obj2);
        end
        
        function out = times(obj1,obj2) 
        % TIMES  Array multiply for BASIS objects
        %
        % TIMES(A,B) is the result of A,*B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also:  times, mtimes, mldivide, mrdivide, rdivide, ldivide.    
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = BF1.*BF2;
            P = P1.*BF2 + BF1.*P2;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = mtimes(obj1,obj2)
        % MTIMES  Multiply for BASIS objects
        %
        % MTIMES(A,B) is the result of A*B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also:  mtimes, times, mldivide, mrdivide, rdivide, ldivide.   
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = BF1*BF2;
            P = P1*BF2 + BF1*P2;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = power(obj1,m)
        % POWER   Array power for BASIS objects
        %
        % POWER(A,B) is the result of A.^B at each point in the domain of A. 
        % B is restricted to be DOUBLE. The partial derivatives of the 
        % resulting object, with regards to the constituent parameters of A 
        % are derived automatically using the chain rule.
        %
        % See also: power, mpower.  
            if ~isnumeric(m)
                error('The second input must be a DOUBLE')
            end
            BF = obj1.BasisFunction.^m;
            P = m*obj1.BasisFunction.^(m-1).*obj1.Partials
            IVName = obj1.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = mpower(obj1,m)
        % MPOWER   Matrix power for BASIS objects
        %
        % MPOWER(A,B) is the result of A^B at each point in the domain of A. 
        % B is restricted to be DOUBLE. The partial derivatives of the 
        % resulting object, with regards to the constituent parameters of A 
        % are derived automatically using the chain rule.
        %
        % See also: mpower, power.  
            if ~isnumeric(m)
                error('The second input must be a DOUBLE')
            end
            BF = obj1.BasisFunction^m;
            P = m*obj1.BasisFunction^(m-1)*obj1.Partials;
            IVName = obj1.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = ldivide(obj1,obj2)
        % LDIVIDE  Left array division for BASIS objects
        %
        % LDIVIDE(A,B) is the result of A.\B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: ldivide, rdivide, mldivide, mrdivide, times, mtimes.      
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = ldivide(BF1,BF2);
            num = -P1.*BF2 + BF1.*P2;
            den = BF1.^2;
            P = den.\num;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = rdivide(obj1,obj2)
        % RDIVIDE  Right array division for BASIS objects
        %
        % RDIVIDE(A,B) is the result of A./B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: rdivide, ldivide, mrdivide, mldivide, times, mtimes.      
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = rdivide(BF1,BF2);
            num = P1.*BF2 - BF1.*P2;
            den = BF2.^2;
            P = num./den;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = mldivide(obj1,obj2)
        % MLDIVIDE  Left division for BASIS objects
        %
        % MLDIVIDE(A,B) is the result of A\B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: mldivide, mrdivide, ldivide, rdivide, mtimes, times.    
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = mldivide(BF1,BF2);
            num = -P1*BF2 + BF1*P2;
            den = BF1^2;
            P = den\num;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function out = mrdivide(obj1,obj2)
        % MRDIVIDE  Right division for BASIS objects
        %
        % MRDIVIDE(A,B) is the result of A/B at each point in the combined
        % domains of A and B. The partial derivatives of the resulting 
        % object, with regards to the constituent parameters of A and B, 
        % are derived automatically using the chain rule.
        %
        % See also: mrdivide, mldivide, rdivide, ldivide, mtimes, times.    
            [obj1ext,obj2ext] = domunion(obj1,obj2);
            
            BF1 = obj1ext.BasisFunction;
            BF2 = obj2ext.BasisFunction;
            P1 = obj1ext.Partials;
            P2 = obj2ext.Partials;
            BF = mrdivide(BF1,BF2);
            num = P1*BF2 - BF1*P2;
            den = BF2^2;
            P = num/den;
            IVName = obj1ext.IVName;
            out = basis(BF,P,IVName);
        end
        
        function [v1ext,v2ext] = domunion(v1,v2)           
            % DOMUNION   Define BASIS objects on a common domain
            %
            % Let A and B be BASIS objects.  If A depends on independent 
            % variables (X,Y) and B depends on independent variables (X,Z) 
            % then [Aext,Bext]=domunion(A,B) returns BASIS objects Aext and 
            % Bext that have a common domain with independent variables 
            % (X,Y,Z). Partials of Aext with regards Z will be zero, and  
            % partials of Bext with regards Y will be zero.

            % Lift inputs to make sure they are both BASIS objects
            if isnumeric(v1)
                v1 = basis(v1,0);
            elseif isnumeric(v2)
                v2 = basis(v2,0);
            end
            
            [uv,~,I] = unique([v1.IVName; v2.IVName]);
            nuv = numel(uv);
            nvar1 = numel(v1.IVName);
            nvar2 = numel(v2.IVName);
            nbasis1 = size(v1.BasisFunction,1);
            nbasis2 = size(v2.BasisFunction,1);
            
            I1 = I(1:nvar1);
            newpartial1 = pmat(zeros(nbasis1,nuv));
            newpartial1(:,I1) = v1.Partials;
            v1ext = basis(v1.BasisFunction,newpartial1,uv);
            I2 = I(nvar1+1:end);
            newpartial2 = pmat(zeros(nbasis2,nuv));
            newpartial2(:,I2) = v2.Partials;
            v2ext = basis(v2.BasisFunction,newpartial2,uv);
        end
        
        function b = isempty(M)
            % ISEMPTY True for empty BASIS objects
            %
            % ISEMPTY(M) returns 1 (true) if the BASIS M has no rows or no columns,
            % and 0 (false) otherwise.
            %
            % See also: isempty.
            szM = size(M);            
            if ~all(szM)
                b = true;
            else
                b = false;
            end
        end
        
    end % end of methods
    
    
    
end % end of classdef




