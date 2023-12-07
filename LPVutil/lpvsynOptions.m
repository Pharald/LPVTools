% LPVSYNOPTIONS  Create a options object for LPV synthesis and analysis.
% 
% opt = lpvsynOptions(Name1,Value1,Name2,Value2,...) creates a options
% object for parameter-varying synthesis and analysis. The LPVSYNOPTIONS
% object is used to specify the parameters of the optimization routines 
% used in the synthesis and analysis functions.
%
% DESCRIPTION
%   Creates an options object for LPVSYN, LPVMIXSYN, LPVNCFSYN, LPVESTSYN,
%   LPVSFSYN, LPVLOOPSHAPESYN, LPVNORM, LPVSTOCHSYN, LPVNORM, LPVWCGAIN.
%
% INPUTS
%   NAME, VALUE pairs: The options property specified by the character string
%      NAME is set to VALUE. The setable options properties are specified
%      below.  The default choice is specified in brackets.
%
%     NAME          |   VALUE     |    Description
%--------------------------------------------------------------------------
%   'Solver'        | ['lmilab']  | Optimization solver to be used.
%   -----------------------------------------------------------------------
%   'SolverOptions' | []          | Options passed directly to the solver.
%   -----------------------------------------------------------------------
%   'SolverInit'    | []          | Initial decision variables for LMI solver.
%   -----------------------------------------------------------------------
%   'Gammalb'       | [1e-6]      | Lower bound on closed-loop induced L2 norm.
%   -----------------------------------------------------------------------
%   'Gammaub'       | [1e6]       | Upper bound on closed-loop induced L2 norm.
%   -----------------------------------------------------------------------
%   'Xlb            | [1e-6]      | X Riccati variable lower bounded by Xlb*I
%   -----------------------------------------------------------------------
%   'Xub            | [1e6]       | X Riccati variable upper bounded by Xub*I
%   -----------------------------------------------------------------------
%   'Ylb            | [1e-6]      | Y Riccati variable lower bounded by Ylb*I
%   -----------------------------------------------------------------------
%   'Yub            | [1e6]       | Y Riccati variable upper bounded by Yub*I
%   -----------------------------------------------------------------------
%   'Method         | ['BackOff'] | String specifying the solution method:
%                   |             |----------------------------------------
%                   | 'MinGamma'  | Minimize the closed-loop induced L2 norm.
%                   |             |----------------------------------------
%                   | 'MaxFeas'   | Maximize the feasibility of the X Riccati 
%                   |             | subject to contraints Gammalb and Gammaub 
%                   |             | on the closed-loop induced L2 norm.
%                   |             |----------------------------------------      
%                   | 'BackOff'   | First solve the 'MinGamma' problem for 
%                   |             | GammaOpt and then solve a second stage
%                   |             | 'MaxFeas' problem with 
%                   |             | Gammaub = BackOffFactor*GammaOpt.
%                   |             | This two-stage solution improves
%                   |             | the numerical conditioning of the 
%                   |             | controller reconstruction.
%                   |             |----------------------------------------  
%                   | 'PoleCon'   | Constrain the closed-loop poles.
%   -----------------------------------------------------------------------              
%   'BackOffFactor' | [1.2]       | Multiplicative factor ( >= 1 ) to back off
%                   |             | the minimum gamma when Method = 'BackOff'.
%   -----------------------------------------------------------------------            
%
% OUTPUT
%   opt: lpvsynOptions object
%
% SYNTAX
%   opt = lpvsynOptions
%     Creates an lpvsynOptions object initialized with default values.
%   opt = lpvsynOptions(Name1,Value1,Name2,Value2,...)
%     Creates an lpvsynOptions object with options specified by the 
%     (Name,value) pairs.
%     
%
% See also lpvsyn.

% 10/26/2010 PJS  Initial Coding


%      -Solver ['lmilab']: Optimization solver to be used.
%          TODO PJS 5/20/2011: Implement Sedumi using LMITRANS
%      -SolverOptions []: Options passed directly to the solver.
%          TODO PJS 5/20/2011: Default LMILab Settings?
%      -SolverInit []:  Initial decision variables for LMI solver.
%      -Gammalb [1e-6]: Lower bound on closed-loop induced L2 norm.
%      -Gammaub [1e6]: Upper bound on closed-loop induced L2 norm.
%      -Xlb [1e-6]: X Riccati variable lower bounded by Xlb*I
%      -Xub [1e6]: X Riccati variable upper bounded by Xub*I
%      -Ylb [1e-6]: Y Riccati variable lower bounded by Ylb*I
%      -Yub [1e6]: Y Riccati variable upper bounded by Yub*I
%      -Method ['BackOff']: String specifying the solution method:
%         'MinGamma': Minimize the closed-loop induced L2 norm.
%         'MaxFeas': Maximize the feasibility of the X Riccati variable
%              subject to contraints Gammalb and Gammaub on the closed-
%              loop induced L2 norm.
%         'BackOff': First solve the 'MinGamma' problem for GammaOpt and
%              then solve a second stage 'MaxFeas' problem with Gammaub =
%              BackOffFactor*GammaOpt. This two-stage solution improves
%              the numerical conditioning of the controller reconstruction.
%         'PoleCon': Constrain the closed-loop poles.
%      -BackOffFactor [1.2]: Multiplicative factor ( >= 1 ) to back off
%         the minimum gamma when Method = 'BackOff'.

classdef lpvsynOptions
    
    properties
        Solver = 'lmilab';
        SolverOptions = [];
        SolverInit = [];
        Gammalb = 1e-6;
        Gammaub = 1e6;
        Xlb = 1e-6;
        Xub = 1e6;
        Ylb = 1e-6;
        Yub = 1e6;
        Method = 'BackOff';
        BackOffFactor = 1.2;
        Jlb = 1e-6;
        Jub = 1e6;
        Llb = 1e-6;
        Lub = 1e6;
        balance = 'on';
    end
    
    methods
        % Constructor
        function opt = lpvsynOptions(varargin)
            % Check # of input args
            nin = nargin;
            if ceil(nin/2)~=floor(nin/2)
                errstr1 = 'lpvsynOptions must have an even number of inputs';
                errstr2 = ' with Name/Value pairs specified together.';
                error([errstr1 errstr2]);
            end
            
            % Set Name/Value pairs:
            % Rely on default error if Name is not a public property
            for i1=1:(nin/2)
                Name = varargin{2*(i1-1)+1};
                Value = varargin{2*i1};
                opt.(Name) = Value;
            end
        end
        
        % Set: Solver
        function opt = set.Solver(opt,value)
            AllowableVal = {'lmilab'; 'sedumi'};
            %'sdpam'; 'csdp'; 'dsdp';'sdpt3';'sdplr'; 'setup'
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.Solver = value;
            else
                errstr1 = 'Solver can be ''lmilab'' or ''sedumi''.';
                errstr2 = [];
                %errstr2 = ' ''csdp'', ''dsdp'', ''sdplr'', or ''sdpt3''.';
                error([errstr1 errstr2]);
            end
        end
        
        % Set: SolverOptions
        function opt = set.SolverOptions(opt,value)
            opt.SolverOptions = value;
        end
        
        % Set: SolverInit
        function opt = set.SolverInit(opt,value)
            opt.SolverInit = value;
        end
        
        % Set: Gammalb
        function opt = set.Gammalb(opt,value)
            if isscalar(value) && isa(value,'double') && value>=0 ...
                               && value<opt.Gammaub
                opt.Gammalb = value;
            else
                error('Gammalb must be a non-negative, scalar number less than Gammaub.');
            end
        end
        
        % Set: Gammaub
        function opt = set.Gammaub(opt,value)
            if isscalar(value) && isa(value,'double') && value>0 ...
                               && value>opt.Gammalb
                opt.Gammaub = value;
            else
                error('Gammaub must be a positive, scalar number greater than Gammalb.');
            end
        end
        
        % Set: Xlb
        function opt = set.Xlb(opt,value)
            if isscalar(value) && isa(value,'double') && value>=0 ...
                               && value<opt.Xub
                opt.Xlb = value;
            else
                error('Xlb must be a non-negative, scalar number less than Xub.');
            end
        end
        
        % Set: Xub
        function opt = set.Xub(opt,value)
            if isscalar(value) && isa(value,'double') && value>0 ...
                               && value>opt.Xlb
                opt.Xub = value;
            else
                error('Xub must be a positive, scalar number greater than Xlb.');
            end
        end
        
        % Set: Ylb
        function opt = set.Ylb(opt,value)
            if isscalar(value) && isa(value,'double') && value>=0 ...
                               && value<opt.Yub
                opt.Ylb = value;
            else
                error('Ylb must be a non-negative, scalar number less than Yub.');
            end
        end
        
        % Set: Yub
        function opt = set.Yub(opt,value)
            if isscalar(value) && isa(value,'double') && value>0 ...
                               && value>opt.Xlb
                opt.Yub = value;
            else
                error('Yub must be a positive, scalar number greater than Ylb.');
            end
        end
        
         % Set: Jlb
        function opt = set.Jlb(opt,value)
            if isscalar(value) && isa(value,'double') && value>=0 ...
                               && value<opt.Jub
                opt.Jlb = value;
            else
                error('Jlb must be a non-negative, scalar number less than Jub.');
            end
        end
        
        % Set: Jub
        function opt = set.Jub(opt,value)
            if isscalar(value) && isa(value,'double') && value>0 ...
                               && value>opt.Jlb
                opt.Jub = value;
            else
                error('Jub must be a positive, scalar number greater than Jlb.');
            end
        end
        
         % Set: Llb
        function opt = set.Llb(opt,value)
            if isscalar(value) && isa(value,'double') && value>=0 ...
                               && value<opt.Lub
                opt.Llb = value;
            else
                error('Llb must be a non-negative, scalar number less than Lub.');
            end
        end
        
        % Set: Lub
        function opt = set.Lub(opt,value)
            if isscalar(value) && isa(value,'double') && value>0 ...
                               && value>opt.Llb
                opt.Lub = value;
            else
                error('Lub must be a positive, scalar number greater than Llb.');
            end
        end       
        
        
        % Set: Method
        function opt = set.Method(opt,value)
            AllowableVal = {'MinGamma'; 'MaxFeas'; 'BackOff'; 'PoleCon'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.Method = value;
            else
                errstr1 = 'Method can be ''MinGamma'', ''MaxFeas'', ';
                errstr2 = '''BackOff'', or ''PoleCon''.';
                error([errstr1 errstr2]);
            end
        end
        
        % Set: BackOffFactor
        function opt = set.BackOffFactor(opt,value)
            if isscalar(value) && isa(value,'double') && value>=1
                opt.BackOffFactor = value;
            else
                error('BackOffFactor must be a scalar >= 1.');
            end
        end

                % Set: BackOffFactor
                function opt = set.balance(opt,value)
                    AllowableVal = {'on'; 'off'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.balance = value;
            else
                error('balance can be ''on'' or ''off''.');
            end
        end
        
    end % methods
end % classdef
