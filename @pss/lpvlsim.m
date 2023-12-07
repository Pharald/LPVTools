function [y,t,x,u,ptrajout] = lpvlsim(P,ptraj,uIn,tIn,x0,opt)
% LPVLSIM  Simulate the time response of a PSS along a parameter trajectory
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

% Determine if user has input a ptraj that contains function handles, 
% or if he is using a standard call with a structure of vectors.
fn = fieldnames(ptraj);
if any(ismember(fn,'Names')) && any(ismember(fn,'Functions'))
    % User is using the function handle call.
    fhflag = true;    
else
    % User is using the standard call with ptraj a structure of vectors.
    fhflag = false;
end

% Deal with array dimensions. Array dimensions only allowed if nargout = 0
szP = size(P);
if nout ==0
    if numel(szP)>2
        % Array dimensions present. Loops through the array.
        szP = size(P);
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

if fhflag 
    newtime = tIn(:); % Column vector of time values:
    newu = uIn;
    Names = ptraj.Names(:);   % column cell-array of names
    if iscell(ptraj.Functions)
        rhoFH = ptraj.Functions(:); % Cell array of function handles
    else
        rhoFH = {ptraj.Functions};
    end
else
    % Map data from u and param time to consistent time grid
    [Values,Names,newtime,newu] = LOCALpgrid2simgrid(ptraj,P.Domain,uIn,tIn);   
end
 Nt = length(newtime);

% Set max allowable timestep
dt = diff(newtime);
if isempty(solverOptions.MaxStep)
    % TODO PJS 6/2/2011 Revisit
    solverOptions.MaxStep = max(dt)/4;   %min(dt)/4;
end

% Get LFT Data: System should have no array dims at this point.
% XXX: Can we replace some of this code using the LPVBLKINIT function?

% XXX Need to handle array dimensions. The Iv dimensions are alays first in
% underlaying data array, followed by the array dimensions. Need to use
% code similar to in DOMUNION to select off the parameter dependent parts
% of the system from the array dimensions. P.Domain vs P.DomainPrivate.
szP = size(P);
% IVNameP = P.DomainPrivate.IVName;
% IVDataP = P.DomainPrivate.IVData;
% DIVDataP = P.DomainPrivate.DIVData;
% NumIVP = P.DomainPrivate.NumIV;
% LIVDataP = P.DomainPrivate.LIVData;
% XXX Assume that at this point we are looping through array and its a
% singleton dimension that doesnt matter. Look at top of file to see how
% array dimenions are handled through a recursive call [todo]
IVNameP = P.Domain.IVName;
IVDataP = P.Domain.IVData;
DIVDataP = P.Domain.DIVData;
NumIVP = P.Domain.NumIV;
LIVDataP = P.Domain.LIVData;

% Check dimensions
if size(newu,2)~= szP(2)
    error('The number of columns of u must equal the number of inputs of P.');
end
if fhflag ==  false
    if size(Values,2)~= NumIVP 
        error('The number fields in ptraj must equal the number of IVs in P.');
    end
end

% Grab SS data for simulation
[a,b,c,d] = ssdata(P);
ABpmat = [a b];
CDpmat = [c d];

% Check Names
if numel(unique(Names))~= numel(Names)
    error('Specified Names cannot have duplicate entries');
end
[tfidx,loc] = ismember(Names,IVNameP);
% XXX - AKP, March 2015.  User-supplied parameter trajetory should be able 
% to have additional parameters that are not in SYSTEM.  Code below does not
% currently allow this.
if all(tfidx) && numel(Names)==numel(IVNameP)
    IVData = IVDataP(loc);
    DIVData = DIVDataP(loc);
    LIVData = LIVDataP(loc);
    RIVData = zeros(NumIVP,2);
    for i=1:NumIVP
        RIVData(i,:) = [IVData{i}(1) IVData{i}(end)];
    end
    
%     ABMat = permute(ABpmat.DataPrivate,[1 2 2+loc']);
%     CDMat = permute(CDpmat.DataPrivate,[1 2 2+loc']);
% XXX work on DAta instead of DataPrivate - see comment online 105
    ABMat = permute(ABpmat.Data,[1 2 2+loc']);
    CDMat = permute(CDpmat.Data,[1 2 2+loc']);
else
    error('Specified Names must include all IVs in P');
end

% Call ode solver to simulate system
% [t,x] = solver(@localLTVsystem,newtime,x0,solverOptions);
[t,x] = feval(solver,@localLTVsystem,newtime,x0,solverOptions);

% ode45 function will treat the call as a time span call if the time vector
% only has two points. Hence, extract only  the first and last point if Nt=2.
if Nt == 2;
   t = t([1 end]);
   x = x([1 end],:);   
end

%% Calculate output, y
y = zeros(Nt,szP(1));
u = zeros(Nt,szP(2));
p = zeros(Nt,NumIVP);
for i = 1:Nt
%     ti = t(i);
%     idx = max(find(ti>=newtime));
%     if idx<Nt
%         frac = (ti-newtime(idx))/(newtime(idx+1)-newtime(idx));
%     elseif idx==Nt
%         idx = Nt-1;
%         frac = 1;
%     else
%         frac = NaN;
%     end
    
%     ut = [(1-frac) frac]*newu([idx, idx+1],:);
    ti = newtime(i);
    ut = newu(i,:);
    xt = x(i,:)';    
    
    % Compute the value of the parameter at the i-th time step
    if fhflag        
        % XXX - Allowing user to make parameter a function of any state.
        % For quasi-LPV one or more of the parameters is a state of the system.
        pt = zeros(1,numel(Names));
        for n=1:numel(Names)
            pt(1,n) = rhoFH{n}(xt,ut,ti);
        end
    else
%         pt = [(1-frac) frac]*Values([idx, idx+1],:);
        pt = Values(i,:);
    end
    
    % Project pValue into the domain.
    pt = min( max(pt', RIVData(:,1)), RIVData(:,2) );
    
    CD = ndLinInterp(CDMat,NumIVP,LIVData,IVData,DIVData,pt(:));
    
    % Compute output and store data
    u(i,:) = ut';
    p(i,:) = pt';
    y(i,:) = CD*[xt; ut(:)];
end

% Store simulated trajectory as a structure
ptrajout.time = newtime;
for i=1:NumIVP
    ptrajout.(Names{i}) = p(:,i);
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
        % Compute the value of the parameter at the i-th time step
        if fhflag
            % XXX - Allowing user to make parameter a function of any state.
            % For quasi-LPV one or more of the parameters is a state of the system.
            pt = zeros(1,numel(Names));
            for n=1:numel(Names)
                pt(1,n) = rhoFH{n}(x,ut,t);
            end
        else
            pt = [(1-frac) frac]*Values([idx, idx+1],:);
        end        
        
        % Project pValue into the domain
        % 11/10/2013: We checked in LOCALpgrid2simgrid that ptraj remains
        % in the domain but (for numerical reasons) the trajectory can
        % slightly exit the domain. Projection needed to avoid interp error
        pt = min( max(pt', RIVData(:,1)), RIVData(:,2) );
        
        AB = ndLinInterp(ABMat,NumIVP,LIVData,IVData,DIVData,pt(:));
        
        xdot = AB*[x;ut(:)];
    end % End nested function


end % End LPVLSIM


function [newtraj,names,newtime,newu] = LOCALpgrid2simgrid(ptraj,Dom,uIn,tIn)
% ptraj = structure with fields 'time', 'p1', 'p2' etc.
% Dom = rgrid object describing the underlying parameter grid.

% Error checking
if isfield(ptraj,'time')
    time = ptraj.time;
    ptraj = rmfield(ptraj,'time');
else
    error('ptraj must contain a field "time" when not using function handles')
end
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

% Check names
ivn = Dom.IVName;
uvn = intersect(names,ivn);
if numel(uvn)~=numel(names) || numel(uvn)~=numel(ivn)
    error('ptraj must contain a field for each of the IVs in Dom')
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
IVData = Dom.IVData;
Np = numel(names);
newtraj = zeros(Nt,Np);
for i=1:Np
    namei = names{i};
    newtraj(:,i) = interp1(time,ptraj.(namei),newtime);
    
    % Check that trajectory lies within domain data
    Domi = IVData{strcmp(namei,ivn)};
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