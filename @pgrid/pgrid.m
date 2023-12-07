% PGRID   Create a time-varying real parameter defined on a grid of points
% 
% A = pgrid(Name,GridData,RateBounds) specifies the name, grid points, and 
% ratebounds for a time-varying real parameter where:
%   1) Name is a character string specifying the name
%   2) GridData is a N-by-1 column vector of sorted values that specify
%      the parameter grid.
%   3) RateBounds is a 1-by-2 row vector specifying lower and upper
%      bounds on the derivative of the parameter with respect to time.
%      Set RateBounds(1)=-inf and/or  RateBounds(2)=+inf to denote an
%      unbounded rate of change. RateBounds are optional, and a two argument 
%      call: A = pgrid(Name,GridData) will set them to a default value of
%      [-inf,+inf].
%
%   % EXAMPLE 1: (CUT/PASTE)  
%   % Create a PGRID "a" which has grid points [1,2,3,4] and ratebounds [-1 1]
%   a = pgrid('a',1:4,[-1 1])
% 
%   % EXAMPLE 2: (CUT/PASTE)
%   % Create a PGRID "b" which has grid points 2:2:20 and default 
%   % ratebounds [-inf inf].
%   b = pgrid('b',2:2:20)
% 
%   % Example 3: (CUT/PASTE)
%   % Use PGRID as a building block to build a PMAT: Create a 2-by-2 matrix 
%   % defined on the 1-dimensional grid specified by the PGRID "a".
%   M = [1, a;a^2, cos(a)]
% 
%   % Example 4: (CUT/PASTE)
%   % Use PGRID as a building block to build a PSS: Create a 1-by-1 
%   % state-space model defined on the 2-dimensional grid specified by the 
%   % PGRIDs "a" and "b".
%   S = ss(-a,b,1,0)
%
% See also rgrid, pmat, pss, pfrd, upmat, upss, upfrd, pstruct.


classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,...
          ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) pgrid < pmat
    % Class definition for time-varying parameter
    
    properties
        Name = 'UNNAMED';
        GridData = 0;
        RateBounds = [-inf inf];
    end
    
    properties (Dependent)
        Range = [];
    end
    
    methods
        % Constructor
        function obj = pgrid(varargin)
            
            if nargin==1
                obj = varargin{1};
            elseif nargin==2 || nargin==3
                obj.Name = varargin{1};
                if isa(obj.GridData,'double')
                    obj.GridData = varargin{2}(:);
                else
                    error('"GridData" must be a double');
                end
                if nargin==3
                    obj.RateBounds = varargin{3};
                end
            elseif nargin~=0
                error('Syntax is A = pgrid(Name,GridData,RateBounds)');
            end
            [pflag,errstr] = isvalid(obj);
            
            if pflag==0
                error(errstr);
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if PGRID object is valid.
            
            errstr = [];
            pflag = 1;
            
            if ~ischar(obj.Name)
                pflag = 0;
                errstr = '"Name" must be a character string.';
            end
            
            GD = obj.GridData;
            if ~( isa(GD,'double') &&  all(diff(GD)>0) )
                pflag = 0;
                errstr = '"GridData" must be a sorted vector of doubles.';
            end
            
            % XXX PJS Should rate bound interval contain 0?
            RB = obj.RateBounds;
            if ~( isa(RB,'double') && length(RB)==2 && diff(RB)> 0 )
                pflag = 0;
                errstr = ['The value of "RateBounds" must be of the form '...
                    '[LOW,HIGH] with LOW < HIGH.'];
            end
            
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % Display
        function s=display(obj)
            GD = obj.GridData;
            RB = obj.RateBounds;
            
            s = ['Gridded real parameter "' obj.Name '"'];
            s = [s ' with '];
            
            Npts = numel(GD);
            R = [GD(1) GD(end)];
            R = sprintf('[%.3g,%.3g]',R(1),R(2));
            s = [s int2str(Npts) ' points in ' R ' and'];
            RB = sprintf('[%.3g,%.3g]',RB(1),RB(2));
            s = [s ' rate bounds ' RB '.'];
            
            if nargout==0
                disp(s)
            end
        end
                        
        function sys = ss(A,varargin)
            % Call PSS constructor
            %
            % SYS = ss(A,B,C,D) creates a parameter-varying state-space model on
            %     the domain of A, B, C and D. SYS will be a PSS.
            %
            % SYS = ss(A)  creates a parameter dependent static gain matrix on the
            %     domain of A.
            %
            % See also: ss, pss.
            sys = ss(pmat(A),varargin{:});
        end

        function out = get.Range(A)
            out = [A.GridData(1) A.GridData(end)];
        end
        
        function [varargout] = size(obj,arg2)
            % SIZE   Size of a PGRID object.
            %
            % S = SIZE(A) returns the row and column dimensions of a PGRID
            % object. S(1) and S(2) are the row and column dimension of S.
            %
            % See also: size, pmat/size.
            
            nout = nargout;
            if nargin==1 && nout<=1
                varargout = {[1 1]};
            elseif nargin==2 && nout<=1
                varargout = {1};
            elseif nargin==1 && nout>1
                varargout = mat2cell(ones(nout,1),ones(nout,1));
            else
                error('Invalid syntax for size of a PGRID.');
            end
            
        end
        
        function out = properties(m)
            % PROPERTIES  Display property names for PGRID.
            
            % Enforce hidden properties
            p = fieldnames(m);
            [~,idx]=ismember({'Data';'Domain';'Parameter'},p);
            p(idx) = [];
            if nargout==0
                sp = '    ';
                Np = length(p);
                for i=1:Np
                    p{i} = [sp p{i}];
                end
                p = [{'Properties for class pgrid:'}; p];
                disp(char(p));
            else
                out = p;
            end
        end
                
        
    end  % methods
    
end % end of classdef



