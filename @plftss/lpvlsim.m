function [y,t,x,u,ptrajout] = lpvlsim(P,ptraj,uIn,tIn,x0,opt)
% LPVLSIM  Simulate the time response of a PLFTSS along a parameter trajectory
%
% [Y,T,X,U,TRAJ] = LPVLSIM(G,PTRAJ,UIN,TIN) simulates the time-response of 
% the system G, subject to the input signal defined by UIN and TIN, and the 
% parameter tracjetory defined in PTRAJ. G is a PSS with Ny outputs, Nx states,
% Nu inputs, and N independent variables IVName1,...,IVNameN. TIN is a sorted 
% column vector of time values, and UIN is a length(TIN)-by-Nu matrix of 
% corresponding inputs. PTRAJ is a struct which defines the time-variation 
% of the parameters (independent variables) in G. The field PTRAJ.time 
% contains a sorted row vector of time-values. PTRAJ must also have a field 
% for each independend variable in G, such that PTRAJ.IVName1, ... ,PTRAJ.IVNameN 
% each contain a row vector of parameter trajectories corresponding to 
% PTRAJ.time. Y is a length(T)-by-NY matrix whose columns correspond to
% the outputs of G, X is a length(T)-by-Nx matrix whose columns 
% correspond to the state trajectories of G, U is a length(T)-by-Nu matrix 
% whose columns correspond to the inputs of G, and T is a column vector of 
% time values corresponding to Y, X and U. TRAJ contains the corresponding 
% parameter trajectories. 
%
% [Y,T,X,U,TRAJ] = LPVLSIM(G,PTRAJ,UIN,TIN,X0) simulates the time-response 
% of the system G starting from the initial condition X0.
%
% See also: lpvstep, lpvinitial, lpvimpulse, lsim.



% AH - 11/10/13 - Simulation time is now restricted to whichever is shorter
% of tIn or ptraj.time.

% TODO PJS 6/2/2011: Currently uses default linear interpolation to
% evaluate along parameter trajectory. Create options with solver,
% solverOptions, interpolation method)


% Parse Inputs
nin = nargin;
nout = nargout;
narginchk(4, 6)
if nin == 4
    x0 = [];
    opt.solver = [];
    opt.imeth =  [];
    opt.solverOptions = [];
elseif nin==5
    opt.solver = [];
    opt.imeth =  [];
    opt.solverOptions = [];
end

% Deal with array dimensions. Array dimensions only allowed if nargout = 0
szP = size(P);
if nout ==0
    if numel(szP)>2
        % Array dimensions present. Loop through the array
        L.type = '()';
        nAD = prod(szP(3:end));
        for i = 1:nAD
            L.subs = {':',':',i};
            Pi = subsref(P,L);
            lpvlsim(Pi,ptraj,uIn,tIn,x0,opt);
            if nAD>1
                if i==1
                    for j=1:szP(1)
                        % switch to "hold on" in all subplots
                        subplot(szP(1),1,j);
                        hold on;
                    end
                end
                if i==nAD
                    for j=1:szP(1)
                        % switch to "hold off" in all subplots
                        subplot(szP(1),1,j);
                        hold off;
                    end
                end
            end
        end
        return
    end
else
    if numel(szP)>2
        error(['The "lpvlsim" command operates on a non-array dimensional'...
            ' model when used with output arguments.'])
    end
end

% Default values
if isempty(x0)
    ns = size(P.Data.a,1);
    x0 = zeros(ns,1);
end
if isempty(opt)
    opt.solver = [];
    opt.imeth =  [];
    opt.solverOptions = [];
end
if isempty(opt.solver)
    solver = @ode45;
else
    solver = opt.solver;
end
if isempty(opt.solverOptions)
    solverOptions = odeset;
else
    solverOptions = opt.solverOptions;
end

if isempty(opt.imeth)
    IMeth = 'linear';
else
    IMeth = opt.imeth;
end

% Error checking
if isfield(ptraj,'time')
    plist = fieldnames( rmfield(ptraj,'time') );
else
    error('ptraj must contain a field "time"')
end

% Get LFT Data: System should have no array dims at this point.
[nU,ABMat,CDMat,IVData,~,niv,LIVData] = lpvblkinit(P,plist);
[M11ab,M12ab,M21ab,M22ab]=deal(ABMat{:});
[M11cd,M12cd,M21cd,M22cd]=deal(CDMat{:});

% Create a matrix of upper and lower bounds for each IV:
R = zeros(niv,2);
for i=1:niv
    R(i,:) = [IVData{i}(1) IVData{i}(end)];
end


% Map data from u and param time to consistent time grid
% XXX TODO PJS: LOCALpgrid2simgrid is different for PLFTSS and PSS
% implementations. Standardize?
[Values,names,newtime,newu] = LOCALpgrid2simgrid(ptraj,P,uIn,tIn);
Nt = length(newtime);

% Set max allowable timestep
dt = diff(newtime);
if isempty(solverOptions.MaxStep)
    % TODO PJS 6/2/2011 Revisit
    solverOptions.MaxStep = max(dt)/4;   %min(dt)/4;
end

% Check dimensions
if size(newu,2)~= szP(2)
    error('The number of columns of u must equal the number of inputs of P.');
end

% Call ode solver to simulate system
[t,x] = feval(solver,@localLTVsystem,newtime,x0,solverOptions);

%% Calculate output, y

% Loop through the Nt simulation time points and compute outputs
y = zeros(Nt,szP(1));
u = zeros(Nt,szP(2));
p = zeros(Nt,niv);
for i = 1:Nt
    ti = t(i);
    % For i-th point in t, find index of time point in newtime which lies
    % closest to it, but below it. Call this idx. Then compute a fraction 
    % which describes where t(i) lies between newtime(idx) and newtime(idx+1) 
    idx = max(find(ti>=newtime));
    if idx<Nt
        frac = (ti-newtime(idx))/(newtime(idx+1)-newtime(idx));
    elseif idx==Nt
        idx = Nt-1;
        frac = 1;
    else
        frac = NaN;
    end
    
    % Find value of u and p which corresponds to t(i)
    ut = [(1-frac) frac]*newu([idx, idx+1],:);
    pt = [(1-frac) frac]*Values([idx, idx+1],:);
    
    % Project pValue into the domain.
    pt = min( max(pt', R(:,1)), R(:,2) );
    
    
    % Normalize parameter value to [-1,1].
    % TODO PJS 1/24/2014: The normalization assumes centered nominal!
    pt = 2*(pt - R(:,1))./(R(:,2) - R(:,1)) -1;
    Del = pt( LIVData{2} );
    r = size(M11ab,1);
    M21Del = lrscale(M21cd,[],Del(:));
    M11Del = lrscale(M11cd,[],Del(:));
    if all(M11Del==0)
        CD = M22cd + M21Del*M12cd;
    else
        CD = M22cd + M21Del*( (eye(r)-M11Del)\M12cd );
    end
    
    % Denormalize parameter data
    pt = (pt.*( R(:,2) - R(:,1) ) + ( R(:,2) + R(:,1)) )/2;
    
    % Compute output and store data
    u(i,:) = ut';
    p(i,:) = pt';
    y(i,:) = CD*[x(i,:)'; ut(:)];
end

% Store simulated trajectory as a structure
ptrajout = ptraj;
ptrajout.time = newtime;
for i=1:niv
    ptrajout.(names{i}) = p(:,i);
end

%% Plot results if narogut = 0
if nout ==0
    [ny,~] = size(P);
    for i = 1:ny
        subplot(ny,1,i)
        plot(t,y(:,i));
        if i ==1
            title('Input');
        end
        if i == ny
            xlabel('Time [sec]');
        end
        ylabel('Output');
    end
end

%% Nested local "odefun" allows access to @LTVsim workspace
    function xdot = localLTVsystem(t,x)
        % Linear interpolation to evaluate u and p at t
        idx = max(find(t>=newtime));
        if idx<Nt
            frac = (t-newtime(idx))/(newtime(idx+1)-newtime(idx));
        elseif idx==Nt
            idx = Nt-1;
            frac = 1;
        else
            frac = NaN;
        end
        ut = [(1-frac) frac]*newu([idx, idx+1],:);
        pt = [(1-frac) frac]*Values([idx, idx+1],:);
        
        % Project pValue into the domain
        % 11/10/2013: We checked in LOCALpgrid2simgrid that ptraj remains
        % in the domain but (for numerical reasons) the trajectory can
        % slightly exit the domain. Projection needed to avoid interp error
        pt = min( max(pt', R(:,1)), R(:,2) );
        
        % Normalize parameter value to [-1,1].
        % TODO PJS 1/24/2014: The normalization assumes centered nominal!
        r = size(M11ab,1);
        pt = 2*(pt - R(:,1))./(R(:,2) - R(:,1)) -1;
        Del = pt( LIVData{1} );
        M21Del = lrscale(M21ab,[],Del(:));
        M11Del = lrscale(M11ab,[],Del(:));
        if all(M11Del==0)
            AB = M22ab + M21Del*M12ab;
        else
            AB = M22ab + M21Del*( (eye(r)-M11Del)\M12ab );
        end
        xdot = AB*[x; ut(:)];
    end % End nested function


end % End LPVLSIM


function [newtraj,names,newtime,newu] = LOCALpgrid2simgrid(ptraj,P,uIn,tIn)
% ptraj = structure with fields 'time', 'p1', 'p2' etc.
% Dom = rgrid object describing the underlying parameter grid.

% Get time and plist names
time = ptraj.time;
ptraj = rmfield(ptraj,'time');
names = fieldnames(ptraj);

% Check ptraj time vector
if numel(time)<2
    error('Length of ptraj.time vector must be greater than 1.');
end
dt = diff(time);
if any(dt<=0) || ~all(isfinite(time)) || ~isreal(time)
    error(['ptraj.time data must be a vector of real, finite, '...
        'and monotonically increasing values.']);
end
for i=1:numel(names)
    namei = names{i};
    if numel(ptraj.(namei))~=numel(time)
        error(['Length of trajectory ' namei '  does not equal the length ptraj.time.']);
    end
end

% Check (tIn,uIn) time vector
if numel(tIn)<2
    error('Length of tIn vector must be greater than 1.');
end
dt = diff(tIn);
if any(dt<=0) || ~all(isfinite(tIn)) || ~isreal(tIn)
    error(['tIn data must be a vector of real, finite, '...
        'and monotonically increasing values.']);
end
if size(uIn,1)~=numel(tIn)
    error('The number of rows of u must equal the length of the time vector.');
end


% Create single consistent time vector
% AH - 11/10/13 - Limits simulation to the shorter of the maximum times
% specified in tIn and ptraj.time.
shortertime = min(time(end),tIn(end));
newtime = unique([time(:); tIn(:)]);
idxtstop = newtime <= shortertime;
newtime = newtime(idxtstop);
Nt = numel(newtime);

% Interpolate inputs on newtime
Nu = size(uIn,2);
newu = zeros(Nt,Nu);
for i=1:Nu
    newu(:,i) = interp1(tIn,uIn(:,i),newtime);
end

% Interpolate parameters on newtime
Np = numel(names);
newtraj = zeros(Nt,Np);
L.type = '.';
L.subs = 'Blocks';
Blocks = subsref(P,L);
for i=1:Np
    namei = names{i};
    newtraj(:,i) = interp1(time,ptraj.(namei),newtime);
    
    % Check that trajectory lies within domain data
    Domi = Blocks.(namei).Range;
    if min(newtraj(:,i))<min(Domi) || max(newtraj(:,i))>max(Domi)
        error('Parameter trajectory lies outside of parameter domain')
    end
end









% % Error checking
% if isfield(ptraj,'time')
%     time = ptraj.time;
%     ptraj = rmfield(ptraj,'time');
% else
%     error('ptraj must contain a field "time"')
% end
%
% numt = numel(time);
% if numt<2
%     error('Length of ptraj.time vector must be greater than 1.');
% end
% dt = diff(time);
% if any(dt<=0) || ~all(isfinite(time)) || ~isreal(time)
%     error(['ptraj.time data must be a vector of real, finite, '...
%         'and monotonically increasing values.']);
% end
%
% trajnames = fieldnames(ptraj);
% ntraj = numel(trajnames);
%
%
% % Get PMAT Data
% ivn = Dom.IVName;
% niv = Dom.NumIV;
% IVData = Dom.IVData;
%
% [uvn,Itraj,IDom] = intersect(trajnames,ivn);
% if numel(uvn)~=ntraj || numel(uvn)~=niv
%     error('ptraj must contain a field for each of the IVs in Dom')
% end
%
% addtime = [];
% for i = 1:ntraj
%     % Pull out trajectory and domain data associated with i-th parameter.
%     trajnamei = trajnames{i};
%     traji = ptraj.(trajnamei);
%     Domi = IVData{strcmp(trajnamei,ivn)};
%     if min(traji)<min(Domi) || max(traji)>max(Domi)
%         error('Parameter trajectory lies outside of parameter domain')
%     end
%
%     for j = 1:numt-1
%         idx = find(Domi>traji(j) & Domi<traji(j+1));
%         if ~isempty(idx)
%             frac = (Domi(idx)-traji(j))/(traji(j+1)-traji(j));
%             addtime = [addtime;(1-frac)*time(j)+frac*time(j+1)];
%         end
%         idx = find(Domi<traji(j) & Domi>traji(j+1));
%         if ~isempty(idx)
%             frac = (traji(j)-Domi(idx))/(traji(j)-traji(j+1));
%             addtime = [addtime;(1-frac)*time(j)+frac*time(j+1)];
%         end
%     end
% end
%
% newtime = unique([time(:);addtime(:)]);
% tl = numel(newtime);
% trajmat = zeros(tl,niv);
% names = cell(niv,1);
% for i = 1:ntraj
%     names{i} = trajnames{i};
%     traji = ptraj.(names{i});
%     traji = interp1(time,traji,newtime);
%     trajmat(:,i) = traji(:);
% end
%
%
%
end