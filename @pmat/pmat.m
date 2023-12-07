% PMAT   Create a parameter-varying matrix
%
%   M = PMAT(Data,Domain) creates a parameter-varying matrix defined on
%   an N-dimensional rectangular grid. Domain is an RGRID object that
%   specifies the N independent variables and the rectangular grid domain.
%   Data is an (N+2) dimensional double array that specifies the matrix
%   data. Data(:,:,i1,...,iN) is the value of the matrix evaluated at
%   the point Domain(i1,....,iN).
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 matrix defined on a 1-dimensional grid
%   IVData = linspace(-2,2,20);
%   Domain = rgrid('a',IVData);  
%   for i=1:length(IVData)
%       Data(1:2,1:2,i) = [1 IVData(i); IVData(i)^2 cos(IVData(i))];
%   end
%   M = pmat(Data,Domain)
% 
%   % Plot entries of M versus the independent variable
%   lpvplot(M);
%
% See also: rgrid, pgrid, pss, pfrd, upmat, upss, upfrd, pstruct.

% TODO PJS 6/14/2011: If I have an N-by-M matrix X and a N-by-M Domain D,
% then it seems like pmat(X,D) should return a 1-by-1 PMAT defined on
% the N-by-M domain. Currently the user must shift the dimensions of X
% so that it is a 1-by-1-by-N-by-M.

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd}) pmat
    % Class definition for parameter-varying matrix
    
    properties (Dependent)
        % Public data with array dimensions trailing after indep vars
        Data = [];
        Domain;
        Parameter;
    end
    
    
    properties (Hidden=true)
        DomainPrivate = rgrid;
        DataPrivate = [];
    end
    
    % NOTE PJS 4/8/2011: M.DomainPrivate.Name should return the
    % IVData associated with the IV Name. Do this in SUBSREF/SUBSASGN
    % Domain and DataPrivate would be invalid IV Names.
    
    methods
        % Constructor
        function obj = pmat(DataPrivate,Domain)
            if nargin==1 && isa(DataPrivate,'pgrid');
                Domain = rgrid(DataPrivate);
                DataPrivate = shiftdim(DataPrivate.GridData(:),-2);
                obj = pmat(DataPrivate,Domain);
            elseif nargin==1 && isa(DataPrivate,'pmat');
                obj = DataPrivate;
            elseif nargin>0
                if nargin==1
                    Domain = rgrid;
                end
                
                % Set properties of a pmat
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
            % ISVALID Determine if PMAT object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = 'Domain must be an rgrid object.';
            elseif ~( isa(obj.DataPrivate,'double') || isa(obj.DataPrivate,'logical') )
                pflag = 0;
                errstr = 'Data must be a double or logical array.';
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
        
        % XXX should user be able to set domain parameters by hand.
        % Set Domain
        function out = set.Domain(obj,Domain)
            InData = obj.Data;
            out = pmat(InData,Domain);
        end
        
        % hasArray
        function out = hasArray(obj)
            % hasArray True for PMAT objects with array dimensions.
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
                s1 = 'PMAT with ';
            else
                s1 = [adcs ' array of PMATs with '];
            end
            s1 = [s1 sprintf('%d',szm(1)) ' rows and ' ...
                sprintf('%d',szm(2)) ' columns.'];
            
            % Make display string
            if niv>0
                dispStr = display(obj.Domain);
                tmp = 'The PMAT consists of the following blocks:';
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
        
        
    end % end of methods
    
    % Numeric element-by-element utility methods
    methods
        function out = abs(mat)
            % ABS   Absolute value for PMAT objects.
            %
            % ABS(M) is the absolute value of the elements of M at each point
            % in the domain of M.
            %
            % See also: abs, sign, angle.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = abs(out.DataPrivate);
        end
        
        function out = acos(mat)
            % ACOS   Inverse cosine with result in radians for PMAT objects.
            %
            % ACOS(M) is the arccosine of the elements of M at each point in
            % the domain of M.
            %
            % See also: acos, cos.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = acos(out.DataPrivate);
        end
        
        function out = acot(mat)
            % ACOT   Inverse cotangent with result in radians for PMAT objects.
            %
            % ACOT(M) is the inverse cotangent of the elements of M at each point in
            % the domain of M.
            %
            % See also: acot, cot.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = acot(out.DataPrivate);
        end
        
        function out = acsc(mat)
            % ACSC   Inverse cosecant with result in radians for PMAT objects.
            %
            % ACSC(M) is the inverse cosecant of the elements of M at each
            % point in the domain of M.
            %
            % See also: acsc, csc.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = acsc(out.DataPrivate);
        end
        
        function out = angle(mat)
            % ANGLE   Phase angle for PMAT objects.
            %
            % ANGLE(M) is the phase angle, in radians, of the elements of M at
            % each point in the domain of M.
            %
            % See also: angle, abs.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = angle(out.DataPrivate);
        end
        
        function out = asec(mat)
            % ASEC   Inverse secant with result in radians for PMAT objects.
            %
            % ASEC(M) is the inverse secant of the elements of M at each point
            % in the domain of M.
            %
            % See also: asec, sec.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = asec(out.DataPrivate);
        end
        
        function out = asin(mat)
            % ASIN   Inverse sine with result in radians for PMAT objects.
            %
            % ASIN(M) is the arcsine of the elements of M at each point in the
            % domain of M.
            %
            % See also: asin, sin.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = asin(out.DataPrivate);
        end
        
        function out = atan(mat)
            % ATAN   Inverse tangent with result in radians for PMAT objects.
            %
            % ATAN(M) is the arctangent of the elements of M at each point in
            % the domain of M.
            %
            % See also: atan, atan2, tan.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = atan(out.DataPrivate);
        end
        
        function out = atan2(Y,X)
            % ATAN2   Four quadrant inverse tangent for PMAT objects.
            %
            % ATAN2(Y,X) is the four quadrant arctangent of the real parts 
            % of the elements of X and Y at each point in the common domain 
            % of X and Y.
            %
            % See also: atan2, atan.
            
            % Check # of input arguments
            narginchk(2, 2)
            out = binop(Y,X,'atan2');
        end
        
        function out = ceil(mat)
            % CEIL  Ceil for PMAT objects.
            %
            % CEIL(M) rounds the elements of M to the nearest integer
            % towards +infty at each point in the domain of M.
            %
            % See also: ceil, floor, round, fix.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = ceil(out.DataPrivate);
        end
        
        function out = conj(mat)
            % CONJ   Complex conjugate for PMAT objects.
            %
            % CONJ(M) is the complex conjugate of the elements of M at each
            % point in the domain of M.
            %
            % See also: conj, real, imag, i, j.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = conj(out.DataPrivate);
        end
        
        function out = cos(mat)
            % COS   Cosine of the argument in radians for PMAT objects.
            %
            % COS(M) is the cosine of the elements of M at each point in the
            % domain of M.
            %
            % See also: cos, acos.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = cos(out.DataPrivate);
        end
        
        function out = cot(mat)
            % COT   Cotangent of the argument in radians for PMAT objects.
            %
            % COT(M) is the cotangent of the elements of M at each point in
            % the domain of M.
            %
            % See also: cot, acot.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = cot(out.DataPrivate);
        end
        
        function out = csc(mat)
            % CSC   Cosecant of the argument in radians for PMAT objects.
            %
            % CSC(M) is the cosecant of the elements of M at each point in the
            % domain of M.
            %
            % See also: csc, acsc.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = csc(out.DataPrivate);
        end
        
        function out = exp(mat)
            % EXP   Exponential for PMAT objects.
            %
            % EXP(M) is the exponential of the elements of M at each point in
            % the domain of M.
            %
            % See also: exp, log, log10, expm.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = exp(out.DataPrivate);
        end
        
        function out = fix(mat)
            % FIX  Fix for PMAT objects.
            %
            % FIX(M) rounds the elements of M to the nearest integer towards
            % zero at each point in the domain of M.
            %
            % See also: fix, floor, round, ceil.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = fix(out.DataPrivate);
        end
        
        function out = floor(mat)
            % FLOOR  Floor for PMAT objects.
            %
            % FLOOR(M) rounds the elements of M to the nearest integer
            % towards -infty at each point in the domain of M.
            %
            % See also: floor, round, ceil, fix.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = floor(out.DataPrivate);
        end
        
        function out = imag(mat)
            % IMAG   Imaginary part of a PMAT object.
            %
            % IMAG(M) is the imaginary part of M at each point in the
            % domain of M.
            %
            % See also: imag, real, isreal, conj, angle, abs.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = imag(out.DataPrivate);
        end
        
        function out = log(mat)
            % LOG   Natural logarithm for PMAT objects.
            %
            % LOG(M) is the natural logarithm of the elements of M at each
            % point in the domain of M.
            %
            % See also: log, log10, exp, logm.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = log(out.DataPrivate);
        end
        
        function out = log10(mat)
            % LOG10  Base 10 logarithm for PMAT objects.
            %
            % LOG10(M) is the base 10 logarithm of the elements of M at each
            % point in the domain of M.
            %
            % See also: log10, log, exp, logm.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = log10(out.DataPrivate);
        end
        
        function out = real(mat)
            % REAL   Real part of a PMAT object.
            %
            % REAL(M) is the real part of M at each point in the domain of M.
            %
            % See also: real, isreal, imag, conj, angle, abs.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = real(out.DataPrivate);
        end
        
        function out = round(mat)
            % ROUND  Round for PMAT objects.
            %
            % ROUND(M) rounds the elements of M to the nearest integer at
            % each point in the domain of M.
            %
            % See also: round, floor, ceil, fix.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = round(out.DataPrivate);
        end
        
        function out = sec(mat)
            % SEC   Secant of the argument in radians for PMAT objects.
            %
            % SEC(M) is the secant of the elements of M at each point in the
            % domain of M.
            %
            % See also: sec, asec.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = sec(out.DataPrivate);
            
        end
        
        function out = sign(mat)
            % SIGN  Signum function for PMAT objects.
            %
            % SIGN(M) is the sign of the elements of M at each point in the
            % domain of M.
            %
            % See also: sign, abs.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = sign(out.DataPrivate);
        end
        
        function out = sin(mat)
            % SIN   Sine of the argument in radians for PMAT objects.
            %
            % SIN(M) is the sine of the elements of M at each point in the
            % domain of M.
            %
            % See also: sin, asin.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = sin(out.DataPrivate);
        end
        
        function out = sqrt(mat)
            % SQRT   Square root for PMAT objects.
            %
            % SQRT(M) is the square root of the elements of M at each point
            % in the domain of M.
            %
            % See also: sqrt, sqrtm.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = sqrt(out.DataPrivate);
        end
        
        function out = tan(mat)
            % TAN   Tangent of the argument in radians for PMAT objects.
            %
            % TAN(M) is the tangent of the elements of M at each point in the
            % domain of M.
            %
            % See also: tan, atan, atan2.
            out = pmat(mat);     % PROMOTE TO PMAT
            out.DataPrivate = tan(out.DataPrivate);
        end
        
        function out = and(A,B)
            % &  Logical AND for PMAT objects
            %
            % A & B performs a logical AND of PMATs A and B at each point in
            % the combined domains of A and B.  This returns a PMAT of the
            % same size with elements set to logical 1 where the relation is
            % true and the elements to logical 0 where it is not.
            %
            % See also: and.
            narginchk(2, 2)
            out = binop(A,B,'and');
        end
        
        function out = or(A,B)
            % |  Logical OR for PMAT objects
            %
            % A | B performs a logical OR of PMATs A or B at each point in the
            % combined domains of A and B.  This returns a PMAT of the same
            % size with elements set to logical 1 where the relation is true
            % and elements to logical 0 where it is not.
            %
            % See also: or.
            narginchk(2, 2)
            out = binop(A,B,'or');
        end
        
        function out = xor(A,B)
            % |  Logical Exclusive OR for PMAT objects
            %
            % XOR(A,B) performs a logical exclusive OR of PMATs A or B at each
            % point in the combined domains of A and B.  This returns a PMAT
            % of the same size with elements set to logical 1 where the
            % relation is true and elements to logical 0 where it is not.
            %
            % See also: xor.
            narginchk(2, 2)
            out = binop(A,B,'xor');
        end
        
        function out = uplus(mat)
            % UPLUS   Unary plus for PMAT objects.
            %
            % +M is the unary plus of M at each point in the domain of M.
            %
            % See also: uplus, uminus.
            out = pmat(mat);     % PROMOTE TO PMAT
        end
        
        function out = uminus(m)
            % UMINUS   Unary minus for PMAT objects.
            %
            % -M is the unary minus of M at each point in the domain of M.
            %
            % See also: uminus, uplus.
            out = pmat(m);     % PROMOTE TO PMAT
            out.DataPrivate = -out.DataPrivate;
        end
        
        function out = ne(A,B)
            % ~=  Not equal for PMAT objects
            %
            % A ~= B does element by element comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with elements set to logical 1 where the
            % relation is true and elements set to logical 0 where it is not.
            %
            % See also: ne.
            narginchk(2, 2)
            out = binop(A,B,'ne');
        end
        
        function out = eq(A,B)
            % ==  Equal for PMAT objects
            %
            % A == B does element by element comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with elements set to logical 1 where the
            % relation is true and elements set to logical 0 where it is not.
            %
            % See also: eq.
            narginchk(2, 2)
            out = binop(A,B,'eq');
        end
        
        function out = ge(A,B)
            % >= Greater than or equal for PMAT objects
            %
            % A >= B does element by element comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with elements set to logical 1 where the
            % relation is true and elements set to logical 0 where it is not.
            %
            % See also: ge.
            narginchk(2, 2)
            out = binop(A,B,'ge');
        end
        
        function out = gt(A,B)
            % >  Greater than for PMAT objects
            %
            % A > B does egtment by egtment comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with egtments set to logical 1 where the
            % relation is true and egtments set to logical 0 where it is not.
            %
            % See also: gt.
            narginchk(2, 2)
            out = binop(A,B,'gt');
        end
        
        function out = le(A,B)
            % <= Less than or equal for PMAT objects
            %
            % A <= B does element by element comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with elements set to logical 1 where the
            % relation is true and elements set to logical 0 where it is not.
            %
            % See also: le.
            narginchk(2, 2)
            out = binop(A,B,'le');
        end
        
        function out = lt(A,B)
            % < Less than for PMAT objects
            %
            % A < B does eltment by eltment comparisons between PMATs A and B
            % at each point in the combined domains of A and B.  This returns
            % a PMAT of the same size with eltments set to logical 1 where the
            % relation is true and eltments set to logical 0 where it is not.
            %
            % See also: lt.
            narginchk(2, 2)
            out = binop(A,B,'lt');
        end
        
        function out = not(A)
            % |  Logical NOT for a PMAT object
            %
            % ~A performs a logical NOT of a PMAT A at each point in the
            % domain of A. This returns a PMAT of the same size with elements
            % set to logical 1 where the relation is true and elements to
            % logical 0 where it is not.
            %
            % See also: not.
            narginchk(1, 1)
            out = pmat(A);     % PROMOTE TO PMAT
            out.DataPrivate = ~out.DataPrivate;
        end
        
    end  % element-by-element numerical methods
    
    methods (Access = private)
        % privatesize
        function varargout = privatesize(m,arg2)
            out = size(m.DataPrivate);
            if numel(out)==3
                out = [out 1];
            elseif numel(out) == 2
                out = [out 1 1];
            end
            
            % AH - 7/18/13 - Address the case when Domain has a singleton
            % IV dimension at end that gets cut off in call to SIZE:
            if numel(out)<2+m.Domain.NumIV
                lout = numel(out);
                addon = ones(1,2+m.Domain.NumIV-lout);
                out = [out addon];
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



