function [y,t,x,u,ptrajout] = lpvinitial(P,ptraj,x0,T,opt)
% LPVINITIAL  Parameter dependent step response for PSS objects
% 
% [Y,T,X,U,PTRAJOUT] = LPVINITIAL(SYS,PTRAJ,XO) computes the parameter 
% dependent response of SYS to an initial value X0. SYS is a PSS 
% with Ny outputs, Nx states, Nu inputs, and N independent variables 
% IVName1,...,IVNameN. PTRAJ is a struct which defines the time-variation 
% of the parameters (independent variables) in SYS. The field PTRAJ.time 
% contains a sorted row vector of time-values. PTRAJ must also have a field 
% for each independend variable in SYS, such that 
% PTRAJ.IVName1, ... ,PTRAJ.IVNameN each contain a row vector of parameter 
% trajectories corresponding to PTRAJ.time. The output Y is a 
% length(T)-by-NY-by-Nu matrix such that Y(:,i,j) corresponds to the i-th 
% output of SYS due to a step command in the j-th input channel. Similarly 
% X is a length(T)-by-Nx-by-Nu matrix describing the state trajectories of 
% SYS, U is a length(T)-by-Nu-by-Nu matrix describing the trajectory of the 
% inputs to SYS, and T is a column vector of time values corresponding to 
% Y, X and U. PTRAJOUT contains the corresponding parameter trajectories. 
%
% LPVINITIAL(SYS,PTRAJ,X0) generates plots of the parameter dependent 
% response of SYS to an initial value X0.
% 
% [Y,T,X,U,PTRAJOUT] = LPVINITIAL(SYS,ptraj,X0,TFINAL) simulates the response 
% of SYS to an initial value X0 up to the time TFINAL.
%
% [Y,T,X,U,PTRAJOUT] = LPVINITIAL(SYS,ptraj,X0,T) simulates the response of 
% SYS to an initial value X0 using a user supplied time vector T.
% 
% See also: lpvlsim, lpvstep, lpvimpulse, lsim.

% XXX - OPTIONS call not implemented yet.
% [Y,T,X,U,PTRAJOUT] = lpvstep(SYS,ptraj,...,OPTIONS)

% TODO - AH - 11/14/13 - Need to add functionality that pads time vector
% when it is too course to capture very fast dynamics. Choose time vector
% padding based on magnitude of fastest pole?

% XXX - Fix handling of nargin>=3 case.

narginchk(3,5)

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

opt = [];
if nargin ==3
    if fhflag
        error('Three argument call is not allowed when PTRAJ contains a function handle')
    elseif isfield(ptraj,'time') 
        T = ptraj.time(end);
        if nargout == 0
            lpvinitial(P,ptraj,x0,T);
        else
            [y,t,x,u,ptrajout] = lpvinitial(P,ptraj,x0,T);
        end
        return
    else
        error('ptraj must contain a field "time" when not using function handles')
    end
elseif nargin>3
    if isnumeric(T) && isscalar(T)
        if fhflag
            tIn = linspace(0,T,1000)';
        else
            tIn = [0;T];  
        end
    elseif isnumeric(T) && isvector(T)
        tIn = T(:);
    elseif ~isnumeric(T)
        % XXX address this when opt is implemented.
        opt = T;
        tIn = ptraj.time;
    end
end

uIn = zeros(numel(tIn),size(P,2));

% Deal with array dimensions. Array dimensions only allowed if nargout = 0
if nargout == 0
    lpvlsim(P,ptraj,uIn,tIn,x0,opt);
else
    if ~isequal(P.Domain,squeeze(P.DomainPrivate))
        error(['The "lpvinitial" command operates on a non-array dimensional'...
            ' model when used with output arguments.'])
    else
        [y,t,x,u,ptrajout] = lpvlsim(P,ptraj,uIn,tIn,x0,opt);
    end
end
