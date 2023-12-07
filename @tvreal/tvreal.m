% TVREAL   Create a time-varying real parameter
%
% The syntax specifies the name, range, and ratebounds for a time-varying
% real parameter:
%   A = tvreal(Name,Range,RateBounds)
% where
%   1) Name is a character string specifying the name
%   2) Range is a 1-by-2 row vector specifying the lower and upper limits
%      for the TVREAL.
%   3) RateBounds is a 1-by-2 row vector specifying lower and upper
%      bounds on the derivative of the parameter with respect to time.
%      Set RateBounds(1)=-inf and/or  RateBounds(2)=+inf to denote an
%      unbounded rate of change. RateBounds are optional, and a two argument 
%      call: A = tvreal(Name,GridData) will set them to a default value of
%      [-inf,+inf].
%
%   % EXAMPLE 1: (CUT/PASTE)  
%   % Create a TVREAL "a" which has range [-2 2] and ratebounds [-1 1].
%   a = tvreal('a',[-2 2],[-1 1])
% 
%   % EXAMPLE 2: (CUT/PASTE)
%   % Create a TVREAL "b" which has range [2 20] and default 
%   % ratebounds [-inf inf].
%   b = tvreal('b',[2 20])
% 
%   % Example 3: (CUT/PASTE)
%   % Use TVREAL as a building block: Create a 2-by-2 matrix that depends 
%   % on the TVREAL "a".
%   M = [1, a;a^2, -a]
% 
%   % Example 4: (CUT/PASTE)
%   % Use TVREAL as a building block to build a PLFTSS: Create a 1-by-1 
%   % state-space model that depends on TVREALs "a" and "b".
%   S = ss(-a,b,1,0)
% 
% See also: plftmat, plftss.



% XXX PJS: Need to be careful about the AutoSimplify Mode because some of
% the UREAL auto-simplify routines may not apply for time-vayring params.
% Should probably restrict the AutoSimplify choice in the ISVALID.
% The impact of the AutoSimplify issue is noted in several places.

% XXX PJS: Errors thrown by UREAL should be rethrown as TVREAL.

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,...
          ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) tvreal < plft
    % Class definition for time-varying parameter    
    
    properties (Dependent)
        Name = [];
        Range = [];
    end
    
    methods
        % Constructor
        function obj = tvreal(varargin)
            
            obj.Data = ureal;
            obj.RateBounds = {obj.Data.Name,[-inf, inf]};
            
            if nargin==1
                obj = varargin{1};
            elseif nargin==2
                Name = varargin{1};
                Range = varargin{2};
                RateBounds = [-inf, +inf];
            elseif nargin==3
                Name = varargin{1};
                Range = varargin{2};
                RateBounds = varargin{3};
            elseif nargin~=0
                error('Syntax is A = tvreal(Name,Range,RateBounds)');
            end
            if nargin>1
                % Nominal is centered. This gives in an affine normalizing
                % transformation which simplifies RateBounds transform.
                nom = mean(Range);
                obj.Data = ureal(Name,nom,'Range',Range);
                obj.RateBounds = {Name,RateBounds};
            end
            
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if TVREAL object is valid.
            errstr = [];
            pflag = 1;
            RBInterval = obj.RateBounds{2};
            RBName = obj.RateBounds{1};
            % XXX PJS Should rate bound interval contain 0?
            if ~( isa(RBInterval,'double') && length(RBInterval)==2 ...
                    && diff(RBInterval)> 0 )
                pflag = 0;
                errstr = ['The value of "RateBounds" must be of the form '...
                    '[LOW,HIGH] with LOW < HIGH.'];
            end
            
            if ~strcmp(RBName,obj.Data.Name)
                pflag = 0;
                errstr = 'Rate bound name must match ureal name';
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % Display
        function s=display(obj)
            Data = obj.Data;
            
            s = ['Time-varying real parameter "' Data.Name '"'];
            s = [s ' with'];
            R = Data.Range;
            R = sprintf('[%.3g,%.3g]',R(1),R(2));
            s = [s ' range ' R ' and'];
            RB = obj.RateBounds{2};
            RB = sprintf('[%.3g,%.3g]',RB(1),RB(2));
            s = [s ' rate bounds ' RB '.'];
            
            if nargout==0
                disp(s)
            end
        end
          
        function out = properties(m)
            % PROPERTIES  Display property names for TVREAL.
            
            % Enforce hidden properties
            p = fieldnames(m);
            [~,idx]=ismember({'Data'},p);
            p(idx) = [];
            if nargout==0
                sp = '    ';
                Np = length(p);
                for i=1:Np
                    p{i} = [sp p{i}];
                end
                p = [{'Properties for class tvreal:'}; p];
                disp(char(p));
            else
                out = p;
            end
        end
        
        % XXX PJS:
        % TVREAL should have most (all?) of the same methods that exist
        % for a UREAL. I'll start with simple unary/binary operations.
        
        function [nv,ndist] = actual2normalized(blk,av)
            % ACTUAL2NORMALIZED  Transforms actual values to normalized values.
            % 
            % ACTUAL2NORMALIZED is used to compute the normalized value of
            % an uncertainty block. See ureal/actual2normalized for
            % details.                    
            % 
            % See also:  actual2normalized, normalized2actual.
            
            % XXX PJS RateBounds are also affect during normalization
            [nv,ndist] = actual2normalized(blk.Data,av);
        end
        
        function av = normalized2actual(blk,nv)
            % NORMALIZED2ACTUAL  Transforms normalized values to their actual value.
            % 
            % NORMALIZED2ACTUAL is used to compute the actual values of
            % a normalized uncertainty block. See ureal/normalized2actual for
            % details.                    
            % 
            % See also:  normalized2actual, actual2normalized.
            
            % XXX PJS RateBounds are also affect during conversion
            av = normalized2actual(blk.Data,nv);
        end
        
        %         function [b,samples] = gridureal(a,varargin)
        %             [b,samples] = gridureal(pmatlft(a),varargin{:});
        %         end
        
        function sys = ss(A,varargin)
            % Call PLFTSS constructor
            %
            % SYS = ss(A,B,C,D) creates a parameter-varying state-space 
            %     model. SYS will be a PLFTSS.
            %
            % SYS = ss(A) creates a parameter dependent static gain matrix.
            %
            % See also: ss, plftss.
            
            sys = ss(plftmat(A),varargin{:});
        end
        
    end  % methods
    
end % end of classdef



