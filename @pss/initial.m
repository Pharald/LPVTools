function [Y,T,X] = initial(varargin)
% INITIAL  Pointwise initial condition response for PSS objects
%
% INITIAL(SYS,X0) plots the undriven response of the PSS SYS
% starting from the initial condition X0. The initial condition response
% is plotted for each system in the domain of SYS. 
%
% [Y,T,X] = INITIAL(SYS,X0) returns the simulated response data.
% The output T is a NT-by-1 PMAT which descibes the time 
% vector for the simulation results. Y is a NT-by-NY PMAT which describes 
% the trajectories of the NY outputs of SYS, and X is a NT-by-NX PMAT 
% which describes the trajectories of the NX states of SYS, at each point 
% in the domain.
%
% [Y,T,X] = INITIAL(SYS,X0,TFINAL) simulates the initial condition responses  
% from t=0 to the final time t=TFINAL.
%
% [Y,T,X] = INITIAL(SYS,X0,T) uses the user-supplied time vector T for the 
% simulations.
%
% INITIAL(SYS1,SYS2,...,X0,T) plots the responses of multiple systems on a 
% single plot.
%
% See also: initial, step, impulse, lsim.

% TODO Refactor Code


nin = nargin;
nout = nargout;
flg = false;

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
    xinit = varargin{2};
    szA = size(sys.Data.A);
    szIVAD = szA(3:end);
    if isempty(szIVAD)
        szIVAD = 1;
    end
    nx = szA(1);
    ny = size(sys,1);
    nu = size(sys,2);    
    Domain = sys.Domain;
    npts = prod(szIVAD);
    niv = Domain.NumIV;
    
    if nin<3
        % Calling: initial(SYS,X0)
        flg = true;
        Tin = [];        
    else
        Tin = varargin{3};
        if isscalar(Tin)
            % Calling: initial(SYS,X0,TEND)
            flg = true;            
        else
            % Calling: initial(SYS,X0,T)
            
            % Set the number of simulation points.
            szT = size(Tin);
            tnum = max(szT(1:2));
        end
    end                
  
        
    if flg
        % User supplied a TEND. Must force all initial simulation to use
        % same time vector to be able to package the output into PMATs.
               
        Pmax = 0;
        fastest_sys = sys.Data(:,:,1);       
        % Loop through the pointwise systems to find system
        % with fastest pole:
        for p = 1:npts
            P  = pole(sys.Data(:,:,p));
            largestP= max(abs(P));
            if largestP>=Pmax
                fastest_sys = sys.Data(:,:,p);
                fastest_xinit = xinit.Data(:,:,p);
                Pmax = largestP;
            end
        end
                        
        % Simulate LTI initial response of system with fastest pole
        if isempty(Tin)
            [~,Tnew,~] = initial(fastest_sys,fastest_xinit);  
        else            
            TEND = lpvmax(Tin);
            [~,Tnew,~] = initial(fastest_sys,fastest_xinit,TEND.Data);    
        end
        % Replace Tin=TEND with a time vector.
        varargin{3} = Tnew;  
        
        % Set the number of simulation points.
        tnum = numel(Tnew);
    end
                                            
    
    number_of_simpoints = tnum;
    Y = zeros(number_of_simpoints,ny,npts);
    T = zeros(number_of_simpoints, 1,npts);
    X = zeros(number_of_simpoints,nx,npts);
    % Loop over domain
    for k =1:npts
        % Create cell of same size as varargin
        function_argument_at_k = cell(1,nin);
        
        % Populate function_argument_at_k with the data corresponding to point k
        for vin =1:nin
            if isa(varargin{vin},'pss') || isa(varargin{vin},'pfrd') || ...
               isa(varargin{vin},'pmat') || isa(varargin{vin},'upss') || ...
               isa(varargin{vin},'upfrd') || isa(varargin{vin},'upmat')
                
                function_argument_at_k{vin} = varargin{vin}.Data(:,:,k);                
            else
                function_argument_at_k{vin} = varargin{vin};
            end
        end         
        
        % Run LTI simulation at grid point k
        [Y(:,:,k),T(:,:,k),X(:,:,k)] = initial(function_argument_at_k{:});                   
        
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
    
    initial(incell{:});
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
% initial(incell{:});

