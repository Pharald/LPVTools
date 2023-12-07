% PLFT   Create a parameter-varying superclass object in LFT framework.
%
% M = PLFT(Data,RateBounds) creates a parameter-varying superclass object.
% Data is a UMAT or USS. RateBounds is a N-by-2 cell array listing the
% rate bound information for each independent variable in the PLFT.  
% RateBounds{i,1} is the character string name of the i-th independent 
% variable and RateBounds{i,2} is a sorted real vector of form [Low, High] 
% specifying its rate bounds.  RateBounds must only contain names of UREAL 
% objects that exist in Data and this indicates that the UREALs are actually 
% TVREALs representing the independent variables.
%
%  % EXAMPLE: (CUT/PASTE)
%  % Create a UREAL which will represent uncertainty in the model.
%  U = ureal('u',2,'Range',[1 3]);
%  % Create a UREAL which will represent the independent variable in the model.
%  P = ureal('p',5,'Range',[-5 10]);
%  % Create a matrix that depends on both "p" and "u".
%  M = [u, p*u];
%  % Create a PLFT object which has uncertainty "u" and independent
%  variable "p".
%  Ratebounds = {'p',[-2 2]};
%  S = plft(M,Ratebounds)
%
% See also: tvreal, plftmat, plftss.

% XXX PJS: Errors thrown by UMAT should be rethrown as PLFT.
% XXX PJS: Binary operations simply call the corresponding UMAT operation.
%    We need to be careful that the UMAT binop does not perform an
%    AutoSimplify that is not allowed for time-varying parameters.

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd}) plft
    % Class definition for parameter varying superclass in LFT framework
    

    properties (Hidden = true)
    end

    properties 
        % We want Data to be hidden. In addition, we want RateBounds to
        % be hidden for PLFTMAT/PLFTSS and a public property for TVREAL.
        % They are both listed as public because they need to be accessed
        % across the various children classes (PLFTMAT/PLFTSS/TVREAL).
        % I created a custom PROPERTIES method for PLTMAT/PLFTSS/TVREAL to 
        % so that they don't appear in properties list. However, Data
        % and RateBounds are still publically accessible for all children
        % classes (and they appear in FIELDNAMES method). 
        % There should be a better solution.
        Data = [];
        RateBounds = cell(0,2);
    end
    
    methods (Hidden=true)
        % Get array dimension as a character string
        function [cs,nad] = ad2char(obj)
            sza = size(obj);
            ad = sza(3:end);
            nad = length(ad);
            cs = [int2str(ad(:)) repmat('x',[nad 1])];
            cs = reshape(cs',[1 2*nad]);
            if nad==0
                cs = '1x1';
            elseif nad==1
                cs = [cs '1'];
            else
                cs = cs(1:end-1);
            end
        end
    end
    
    methods
        % Constructor
        function obj = plft(Data,RateBounds)
           
            if nargin == 1                
                if isa(Data,'plft')
                    obj.Data = Data.Data;
                    obj.RateBounds = Data.RateBounds;
                else
                    obj.Data = Data;
                end
            elseif nargin ==2
                obj.Data = Data;
                obj.RateBounds = RateBounds;
            end
                
        end                
        
                        
        % XXX PJS:
        % PLFT should have most (all?) of the same methods that exist
        % for a UMAT. I'll start with simple unary/binary operations.
        
        function [varargout] = size(m,varargin)
            % SIZE   Size of a PLFT object (TVREAL, PLFTMAT & PLFTSS).
            %
            % For an M-by-N PLFT A with P independent variables and D array dimensions, 
            % S = SIZE(A) returns a 1-by-(2+D) row vector. S(1) and S(2) are the row 
            % and column dimensions of A. For i>2, S(i) is the i-th array dimension 
            % of A. 
            %
            % See also: size, length.    
            nout = nargout;
            if nout == 0
                size(m.Data,varargin{:})
            else
                varargout = cell(nout,1);               
                [varargout{:}] = size(m.Data,varargin{:});
            end
        end
        
        function [varargout] = iosize(m,varargin)
            % IOSIZE   Size of inputs/outputs for a PLFT
            %
            % [NY,NU] = iosize(M) returns the number inputs NU, and the 
            % number of outputs NY, for the PLFT object M. For static  
            % and PLFTs NY and NU corresponds to the size of the rows
            % and columns of M.
            % 
            % See also: size.
            nout = nargout;
            if nout == 0
                iosize(m.Data,varargin{:})
            else
                varargout = cell(nout,1);               
                [varargout{:}] = iosize(m.Data,varargin{:});
            end
        end
        
        function out = horzcat(varargin)
        % HORZCAT   Horizontal concatenation of PLFT objects
        %
        % S = HORZCAT(S1,S2,...) performs a concatenation operation of
        % S = [S1 , S2, ...].
        %
        % See also: horzcat, vertcat.
        
            if nargin==1
                out = varargin{1};
            else       
                varargin{1} = plft(varargin{1});
                out = binop(varargin{1},varargin{2},'horzcat');
                if nargin>2
                    out = horzcat(out,varargin{3:end});
                end
            end
        end
        
        function out = vertcat(varargin)
        % VERTCAT  Vertical concatenation of PLFT objects.
        %
        % S = VERTCAT(S1,S2,...) performs a concatenation operation of
        % S = [S1; S2; ...].
        %
        % See also: vertcat, horzcat.
            if nargin==1
                out = varargin{1};
            else   
                varargin{1} = plft(varargin{1});
                out = binop(varargin{1},varargin{2},'vertcat');
                if nargin>2
                    out = vertcat(out,varargin{3:end});
                end
            end
        end
        
        function out = uplus(m)
        % UPLUS   Unary plus for PLFT objects.
        %
        % +M is the unary plus of M.
        %
        % See also: uplus , uminus.   
            out = m;
        end
        
        function out = uminus(m)
        % UMINUS   Unary minus for PLFT objects.
        %
        % -M is the unary minus of M.
        %
        % See also: uminus, uplus.    
            if isa(m,'tvreal')
                m = plftmat(m);
            end
            out = m;
            out.Data = -out.Data;
        end
        
        function out = inv(m)
            % INV   Matrix inverse for PLFT objects.
            %
            % INV(M) is the inverse of M at each point in the domain of M.
            %
            % See also: inv, slash, pinv.
            if isa(m,'tvreal')
                m = plftmat(m);
            end
            out = m;
            out.Data = inv( out.Data );
        end
        
        function out = plus(m1,m2)
        % PLUS  Plus for PLFT objects
        %
        % PLUS(A,B) is the result of A+B. For PLFTSS objects A and B, this 
        % is equivalent to connecting A and B in parallel.
        %
        % See also: plus, parallel.    
            out = binop(m1,m2,'plus');
        end
        
        function out = minus(m1,m2)
        % MINUS  Minus for PLFT objects
        %
        % MINUS(A,B) is the result of A-B.
        %
        % See also: minus, plus.    
            out = binop(m1,m2,'minus');
        end
        
%         function out = times(m1,m2)
%         % TIMES  Array multiply for PLFT objects
%         %
%         % TIMES(A,B) is the result of A.*B.
%         %
%         % See also: times, mtimes.    
%             out = binop(m1,m2,'times');
%         end
        
        function out = mtimes(m1,m2)
        % MTIMES  Multiply for PLFT objects
        %
        % MTIMES(A,B) is the result of A*B. For PLFTSS objects A and B, 
        % this is equivalent to connecting A and B in series.
        %
        % See also:  mtimes, series, mldivide, mrdivide, inv.    
            out = binop(m1,m2,'mtimes');
        end
        
%         function out = rdivide(m1,m2)
%             out = binop(m1,m2,'rdivide');
%         end
%         
%         function out = ldivide(m1,m2)
%             out = binop(m1,m2,'ldivide');
%         end
        
        function out = mrdivide(m1,m2)
        % MRDIVIDE  Right division for PLFT objects
        %
        % MRDIVIDE(A,B) is the result of A/B.
        %
        % See also: mrdivide, mldivide, rdivide, ldivide, inv, mtimes.    
            out = binop(m1,m2,'mrdivide');
        end
        
        function out = mldivide(m1,m2)
        % MLDIVIDE  Left division for PLFT objects
        %
        % MLDIVIDE(A,B) is the result of A\B.
        %
        % See also: mldivide, mrdivide, ldivide, rdivide, inv, mtimes.    
            out = binop(m1,m2,'mldivide');
        end
        
%         function out = power(m1,m2)
%             out = binop(m1,m2,'power');
%         end
        
        function out = mpower(m1,m2)
            % MPOWER   Matrix power for PLFT objects
            %
            % MPOWER(A,B) is the result of A^B.
            %
            % See also: mpower, power.
            out = binop(m1,m2,'mpower');
        end
        
        function out = blkdiag(varargin)
            % BLKDIAG   Block diagonal concatenation of PLFT objects.
            %
            % M = BLKDIAG(M1,M2, ...) returns the block diagonal concatenation of
            % M1, M2,... 
            %
            % See also: blkdiag, diag, horzcat, vertcat.
            if nargin==1
                out = varargin{1};
            else
                varargin{1} = plft(varargin{1});
                out = binop(varargin{1},varargin{2},'blkdiag');
                if nargin>2
                    out = blkdiag(out,varargin{3:end});
                end
            end
        end
                
        function out = append(varargin)
            % APPEND   Block diagonal concatenation of PLFT
            %
            % S = APPEND(S1,S2, ... ,SN) returns the block diagonal concatenation of
            % S1, S2,..., SN. This is the same functionality as BLKDIAG.
            %
            % See also:  append, blkdiag, series, parallel, feedback.
            if nargin==1
                out = varargin{1};
            else
                varargin{1} = plft(varargin{1});
                out = binop(varargin{1},varargin{2},'append');
                if nargin>2
                    out = append(out,varargin{3:end});
                end
            end
        end
        
        % XXX PJS: Need to double check that stack can be called
        % recursively similar to horzcat, vertcat, etc.
        function out = stack(dim,varargin)  
            % STACK   Stack PLFT objects into an array
            %
            % M = stack(ARRAYDIM,M1,M2,...) creates an array M with the
            % models M1, M2, ... arranged along array dimensions ARRAYDIM.
            % For example stack(2,M1,M2) will produce a 1x2 array with 
            % M(:,:,1,1) = M1 and M(:,:,1,2) = M2. The models M1, M2, ... 
            % must have identical I/O dimensions.
            %             
            % See also: stack, blkdiag, append, repmat.
            if nargin<3
                out = varargin{1};
            else    
                varargin{1} = plft(varargin{1});
                out = binop(varargin{1},varargin{2},'stack',dim);
                if nargin>3
                    out = stack(dim,out,varargin{3:end});
                end
            end
        end
        
        function out = cat(dim,varargin)
            % CAT Concatenate arrays of PLFT objects.
            % 
            % cat(DIM,A,B) concatenates the PLFT objects A and B along the
            % dimension DIM. cat(1,A,B) is equivalent to vertcat(A,B) or
            % [A;B]. cat(2,A,B) is equivalent to horzcat(A,B) or [A,B].
            %
            % See also: cat.
            switch dim
                case 1
                    out = vertcat(varargin{:});
                case 2
                    out = horzcat(varargin{:});
                otherwise
                    out = stack(dim-2,varargin{:});
            end
        end
        
        %XXX
        %         function [b,samples] = gridureal(a,varargin)
        %             [b,samples] = gridureal(plft(a),varargin{:});
        %         end                                              
        
        function out = feedback(m1,m2,varargin)
            % FEEDBACK  Feedback connection of two PLFT objects.
            %
            % M = FEEDBACK(M1,M2) computes a closed-loop model M for the negative
            % feedback interconnection with M1 in the forward loop and M2 in the
            % feedback path, as shown in the figure below. To apply positive 
            % feedback, use the syntax M = FEEDBACK(M1,M2,+1).
            %
            %     ----->0------>[ M1 ]----------+---->
            %           |-                      |
            %           |                       |
            %           +-------[ M2 ]<----------
            %
            % M = FEEDBACK(M1,M2,FEEDIN,FEEDOUT,SIGN) builds a more general feedback
            % interconnection of PLFT objects.  Refer to the help for 
            % InputOutputModel/FEEDBACK for more details on this general interconnection 
            % syntax.
            %
            % See also: feedback, lft, parallel, series.
            out = binop(m1,m2,'feedback',varargin{:});
        end

        function out = series(m1,m2,varargin)
            % SERIES  Series connection of two PLFT objects.
            %
            % M = SERIES(M1,M2,OUTPUTS1,INPUTS2) connects the PLFT objects M1 and M2 in
            % series. The vectors of indices OUTPUTS1 and INPUTS2 specify which
            % outputs of M1 and which inputs of M2 are connected together. M is the
            % series interconnection at each point in the combined domains of M1 and
            % M2. Refer to the help for InputOutputModel/SERIES for more 
            % details on this interconnection.
            %
            % If OUTPUTS1 and INPUTS2 are omitted, SERIES connects M1 and M2 in
            % cascade and returns M = M2 * M1.
            %
            % See also: series, append, parallel, feedback, lft.
            out = binop(m1,m2,'series',varargin{:});
        end
        
        function out = parallel(m1,m2,varargin)
            % PARALLEL  Parallel connection of two PLFT objects.
            %
            % M = PARALLEL(M1,M2,IN1,IN2,OUT1,OUT2) connects the PLFT objects M1 and M2
            % in parallel. The inputs specified by IN1 and IN2 are connected and the
            % outputs specified by OUT1 and OUT2 are summed. Refer to the help for 
            % InputOutputModel/PARALLEL for more details on this interconnection.
            %
            % If IN1,IN2,OUT1,OUT2 are omitted, PARALLEL forms the standard parallel
            % interconnection of M1 and M2 and returns M = M1 + M2.
            %
            % See also: parallel, append, series, feedback, lft.
            out = binop(m1,m2,'parallel',varargin{:});
        end
        
        function out = lft(m1,m2,varargin)
            % LFT  Generalized feedback interconnection of PLFT objects.
            %
            % M = LFT(M1,M2,NU,NY) forms the feedback interconnection of M1 above M2.
            % NU specifies the number of outputs of M2 that are fed into M1.
            % NY specifies the number of outputs of M1 that are fed into M2.
            %
            % See also: lft, feedback.
            out = binop(m1,m2,'lft',varargin{:});
        end
        
        function b = repmat(a,varargin)
            % REPMAT   Replicate and tile for PLFT objects
            %   
            % B = REPMAT(A,M,N) creates a larger PLFT B consisting of an M-by-N tiling 
            % of copies of A. B = REPMAT(A,[M N]) also creates an M-by-N tiling of A.
            %
            % B = REPMAT(A,M) creates a B that is an M-by-M tiling of A.
            %
            % See also: repmat, repsys.
            b = a;
            b.Data = repmat(a.Data,varargin{:});
        end
        
        function b = flipud(a)
            % FLIPUD   Flip in up/down direction for PLFT objects.
            %
            % FLIPUD(A) preserves the columns of A and flips the rows of A 
            % in the up/down direction.
            %
            % See also: flipud, fliplr.
            b = a;
            b.Data = flipud(a.Data);
        end
        
        function b = fliplr(a)
            % FLIPLR   Flip in left/right direction for PLFT objects.
            %
            % FLIPLR(A) preserves the rows of A and flips the columns of A 
            % in the left/right direction.
            %
            % See also: fliplr, flipud.
            b = a;
            b.Data = fliplr(a.Data);
        end
        
        function b = permute(a,order)
            % PERMUTE  Permute row/column/array dimensions of a PLFT object.
            %
            % B = PERMUTE(A,ORDER) arranges the array dimensions of A to be in the 
            % order specified by ORDER.  It is not possible to permute between the 
            % row/column dimensions of A and its array dimensions.
            %
            % See also: permute, ipermute, transpose, ctranspose.    
            b = a;
            b.Data = permute(a.Data,order);
        end
        
        function b = ipermute(a,order)
            % IPERMUTE  Inverse permute row/column/array dimensions of a PLFT object.
            %
            % A = IPERMUTE(B,ORDER) is the inverse of PERMUTE. The array
            % dimensions of B are rearranges so that B = PERMUTE(A,ORDER).
            % Note that it is not possible to permute between the row/column  
            % dimensions of a PLFT and its array dimensions
            %
            % See also: ipermute, permute, transpose, ctranspose.
            b = a;
            b.Data = ipermute(a.Data,order);
        end
        
        function b = reshape(a,varargin)
            % RESHAPE   Reshape PLFT object arrays
            %
            % B = RESHAPE(A,M,N) reshapes A into an M-by-N PLFT object. The elements of A are 
            % taken columnwise from A.  A must have M*N elements. B = RESHAPE(A,[M N]) 
            % is an alternative syntax           
            %
            % B = RESHAPE(A,M,N,P,..,Q) reshapes A into the M-by-N-by-P-by-...-by-Q
            % array B. M*N*P*...*Q must equal PROD(SIZE(A)).
            %
            % B = RESHAPE(A, ..., [], ... ) leaves one desired dimension unspecified, 
            % and replaces the corresponding entry with []. This dimension is computed 
            % automatically so that the product of the desired dimensions matches PROD(SIZE(A)).
            % Only one occurance of [] can be used.
            %
            % See also: squeeze, shiftdim, colon.
            b = a;
            b.Data = reshape(a.Data,varargin{:});
        end
        
        function out = ctranspose(m)
        % CTRANSPOSE  Complex conjugate transpose for PLFT objects.
        %
        % B=CTRANSPOSE(A) returns the complex conjugate transpose of A.
        %
        % For continuous-time PLFTSS, if A has the transfer function H(s) at a point 
        % in the domain then B is the system with transfer function H(-s).'
        %
        % For discrete-time PLFTSS, if A has the transfer function H(z) at a point 
        % in the domain then B is the system with transfer function H(1/z).'
        %
        % See also: ctranspose, transpose, permute.    
            out = ctranspose(m.Data);
        end
        
        function out = transpose(m)
        % TRANSPOSE  Transpose for PLFT objects.
        %
        % B=TRANSPOSE(A) returns the non-conjugate transpose of A.
        % 
        % If A has the transfer function H(s) then B is the system with 
        % transfer function H(s).'
        %
        % See also: transpose, ctranspose, permute.  
            out = transpose(m.Data);
        end
        
        function b = simplify(a,varargin)
            % SIMPLIFY   Simplify a PLFT object.
            %
            % B,simplify(A) remove redundant copies of uncertain elements 
            % and parameters in A. See  InputOutputModel/simplify for
            % further details and additional call options.
            %
            % See also: UncertainBlock.
            b = a;
            b.Data = simplify(a.Data,varargin{:});
        end
        
        function arg = end(obj,slot,nslots)
            % END   Overloaded END for PLFT objects.
            %
            % END(MAT,DIM,NUMINDICES) returns the index corresponding to the last
            % entry along the dimension DIM in the PLFT.  NUMINDICES is the number of
            % indices used in the indexing expression.
            %
            % See also: end.
            szm = size(obj.Data); % Use SIZE since only working on Array dims
            if slot < nslots
                arg = szm(slot);
            else
                if slot ==1
                    error('NUMINDICES must be  larger than 1 for uncertain objects.');
                elseif slot == 2
                    arg = szm(slot);
                else
                    % Handle mat(3:end) or mat(:,:,:,:,4:end) where remaining
                    % dimensions are aggregated into last index.
                    arg = prod(szm(slot:end));
                end
            end
        end
        
        function b = isempty(M)
            % ISEMPTY True for empty PLFT objects
            %
            % ISEMPTY(M) returns 1 (true) if the PLFT M has no rows or no columns,
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
        
    end  % methods
end % end of classdef



