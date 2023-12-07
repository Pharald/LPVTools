function [Y,T,X] = lsim(varargin)
% LSIM  Pointwise time-domain response of PSS objects to arbitrary inputs
%
% LSIM(SYS,U,T) plots the time-domain response of the PSS SYS to 
% the input signal described by U and T.  T is a NT-by-1 column vector 
% which describes the time vector for the simulation. U is a NT-by-NU 
% matrix which describes the trajectory of each of the NU inputs to SYS.
% 
% [Y,TOUT,X] = LSIM(SYS,U,T) returns the simulated response data.
% The output TOUT is a NOUT-by-1 PMAT which descibes the time vector for 
% the simulation results. Y is a NOUT-by-NY PMAT which describes the 
% trajectories of the NY outputs of SYS, and X is a NOUT-by-NX PMAT which 
% describes the trajectories of the NX states of SYS, at each point in the
% domain.
%
% [Y,TOUT,X] = LSIM(SYS,U,T,X0) specifies the initial state vector X0 at time T(1).
% X0 is set to zero when omitted.
%
% LSIM(SYS1,SYS2,...,U,T,X0) plots the responses of multiple systems on a
% single plot.
%
% See also: lsim, step, impulse, initial.


nin = nargin;
nout = nargout;

if nout >0
    % For output case, we will cycle through the domain and collect simulations
    % We will aggregate the simulation results and package the data into PMATs.
    
    % Makse sure system shares a common domain with other inputs.
    varargin{1} = domunion(varargin{:});
    sys = varargin{1};
        
    % Merge domains of function arguments: careful that domunion on multiple
    % arguments will upgrade all to the class with the highest precedence.
    for i = 2:nin            
        if  isa(varargin{i},'pmat') || isnumeric(varargin{i})
            varargin{i} = domunion(varargin{i},pmat(sys.Data.A,sys.Domain));
        end
    end    
    
    % Get Data
    szA = size(sys.Data.A);
    szIVAD = szA(3:end);
    nx = szA(1);
    ny = size(sys,1);
    nu = size(sys,2);    
    Domain = sys.Domain;
    npts = prod(szIVAD);
    
    number_of_simpoints = size(varargin{2},1);
    Y = zeros(number_of_simpoints,ny,npts);
    T = zeros(number_of_simpoints, 1,npts);
    X = zeros(number_of_simpoints,nx,npts);
    % Loop over domain
    for k =1:npts
        % Create cell of same size as varargin
        function_argument_at_k = cell(1,nin);
        
        % Populate function_argument_at_k with the data corresponding to point k
        for vin =1:nin
            if isa(varargin{i},'pss') || isa(varargin{i},'pfrd') || ...
               isa(varargin{i},'pmat') || isa(varargin{i},'upss') || ...
               isa(varargin{i},'upfrd') || isa(varargin{i},'upmat')
                
                function_argument_at_k{vin} = varargin{vin}.Data(:,:,k);                
            else
                function_argument_at_k{vin} = varargin{vin};
            end
        end
        % Run LTI simulation at grid point k
        [Y(:,:,k),T(:,:,k),X(:,:,k)] = lsim(function_argument_at_k{:});
    end
    
    % Reshape output to be [number_of_simpoints,~,IV,AD]
    Y = reshape(Y,[number_of_simpoints ny szIVAD]);
    T = reshape(T,[number_of_simpoints 1  szIVAD]);
    X = reshape(X,[number_of_simpoints nx szIVAD]);
    
    % Pack up data
    Y = pmat(Y,Domain);
    T = pmat(T,Domain);
    X = pmat(X,Domain);
    
else
    incell = cell(1,nin);
    for i=1:nin
        if isa(varargin{i},'pss') || isa(varargin{i},'upss')
            incell{i} = varargin{i}.DataPrivate;
        elseif isa(varargin{i},'pfrd') || isa(varargin{i},'upfrd')
            incell{i} = varargin{i}.DataPrivate;
        else
            incell{i} = varargin{i};
        end
    end
    
    lsim(incell{:});
end




%% OLD CODE
% nin = nargin;
% nout = nargout;
% if nout>0
%    strerr = ' with a PSS cannot be called with output arguments.';
%    error([upper(mname) strerr])
% end
%
% incell = cell(1,nin);
% for i=1:nin
%    if isa(varargin{i},'pss')
%       incell{i} = varargin{i}.DataPrivate;
%    elseif isa(varargin{i},'pfrd')
%       incell{i} = varargin{i}.DataPrivate;
%    else
%       incell{i} = varargin{i};
%    end
% end
%
% lsim(incell{:});

