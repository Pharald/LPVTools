function [F,Gamma,Info] = lpvsfsyn(P,ncont,Xb,alg,opt)
% LPVSFSYN  Parameter-dependent state feedback controller synthesis for PSS
% 
% [F,GAM,INFO] = LPVSFSYN(P,NCON,'L2') computes a parameter-varying 
% state-feedback controller for the parameter-varying system P. 
% NCON specifies the number of available control inputs in P. F is the 
% state-feedback controller for the plant P, which minimizes the L2 norm 
% from the input of P to its output. GAM is the minimum L2 norm 
% achived by F. INFO is a struct with additional data. 
%
% [F,GAM,INFO] = LPVSFSYN(P,NCON,'LQG') computes a parameter-varying 
% state-feedback controller F, which minimizes the stochastic LPV bound.  
% The stochastic LPV bound is defined as the expected value of the average 
% instantaneous power of the output of P, assuming its inputs are zero mean, 
% white-noise processes with unit intensity.
%
% [F,GAM,INFO] = LPVSFSYN(P,NCON,Xb,Yb,ALG) performs a rate-bounded synthesis.
% Xb and Yb are BASIS objects specifying the basis functions to be used in 
% the synthesis. ALG can be either 'L2' or 'LQG'. A call without the ALG
% argument is equivalent to [F,GAM,INFO] = LPVSFSYN(P,NCON,Xb,Yb,'L2').
%
% See also: lpvsynOptions, lpvestsyn, lpvsyn, lpvncfyn, lpvmixsyn, lpvloopshape.



% XXX - OPTIONS are currently only used if 'L2' is selected. 
% [L,GAM,INFO] = LPVESTSYN(P,NMEAS,...,opt) permits the user to pass in 
% See also: lpvsynOptions.

% TODO PJS 5/17/2011:  Combine rate bounds into a single niv-by-3 cell array?
% RateUB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}
% RateLB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}

% Parse Inputs
% XXX Check later
% XXX OPT is not being used.
nin = nargin;
narginchk(2, 4)
nout = nargout;
if nin==2
    opt = lpvsynOptions;
    Xb = [];
    alg = 'L2';
elseif nin==3
    if isa(Xb,'lpvsynOptions')
        opt = Xb;
        Xb = [];
        alg = 'L2';  
    elseif isa(Xb,'basis')
        alg = 'L2'; 
        opt = lpvsynOptions;
    elseif isa(Xb,'char')
        alg = Xb;
        Xb = [];
        opt = lpvsynOptions;
    else
        error(['The third argment in a call to lpvsfsyn must be '...
            'either a BASIS object, a CHAR, or a lpvsynOptions'])
    end
elseif nin==4
    if isa(alg,'lpvsynOptions')
        opt = alg;
        alg = 'L2';        
    elseif isa(alg,'char')
        opt = lpvsynOptions;
    else
        error(['The fourth argment in a call to lpvsfsyn must be '...
            'either a a CHAR, or a lpvsynOptions'])
    end
end

% Always assume one basis function: constant
if isempty(Xb)
   Xb = basis(1,0); 
end

% XXX Does not work for arrays of LPV systems.
% Check for this and error out if there are array (non-LPV) dimensions
% if hasArray(sys.Domain), error, end

% Map the input data into engine data form:
[Pdata,RateBounds,Fbasis,Fgrad] = basis2data(P,Xb);

% Single balancing transformation
% TODO PJS 8/30/2013: Do this before or after orthogonalization?
% (orthogonalization occurs in the engine)
% Need to undo the coordinate transformation in the state-feedback gain.
% nd = size(P,2) - ncont;
% ne = size(P,1);
% blk = [nd ne; ncont nx];
% P = lpvbalance(P);


if strcmpi(alg,'L2')
    % Run the L2 synthesis engine:
    [F,Gamma,Info] = lpvL2sfsynengine(Pdata,ncont,Fbasis,Fgrad,RateBounds,opt);

    % TODO HP 07/12/2023: Add lpvsynoptions for stochastic synthesis engine
    % The options are only used for L2 synthesis at the moment
elseif strcmpi(alg,'LQG')
    % Run the Stochastic synthesis engine:
    [F,Gamma,Info] = lpvLQGsfsynengine(Pdata,ncont,Fbasis,Fgrad,RateBounds);
end



nx = size(Pdata.A,1);
Dom = P.Domain;
F = pmat(reshape(F,[ncont; nx; Dom.LIVData]'),Dom);




