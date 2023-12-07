function [gamopt,Info] = lpvnorm(P,omega,polelist)
% LPVNORM  Compute bound on norm for PLFTSS systems.
%
% LPVNORM computes a bound on the norm of a PLFTSS system over the set of
% all permissible trajectories of the independend variables which the  
% PLFTSS depends on. 
%
% [GAMMA,INFO] = lpvnorm(P) computes an upper bound GAMMA on 
% the induced L2 norm of the PLFTSS P. If the system P is rate-unbounded 
% the analysis computes a induced L2 norm bound GAMMA that is valid for 
% arbitrarily fast variations of the system parameters. If the system is
% rate-bounded the induced L2 norm bound is valid for any permissible 
% parameter trajectory that doesn't violate those parameter rate-bounds.
% 
% [GAMMA,INFO] = lpvnorm(P,OMEGA) allows the user to specify a custom 
% frequency vector OMEGA for the analysis.
% 
% [GAMMA,INFO] = lpvnorm(P,OMEGA,POLELIST) allows the user to specify basis
% functions for the Integral Quadratic Constraints (IQCs) used in the
% analysis. POLELIST is a 1xN DOUBLE row vector, of negative values. Each 
% value in POLELIST corresponds to a pole of a stable transfer function 
% that is used as a weight on all signals in the IQCs. A default POLELIST, 
% with five pole values, is used when a POLELIST is not supplied by the 
% user. The five pole values are selected automatically from the frequency 
% range of the system dynamics.
%
% See also: norm, lpvwcgain, lpvsyn.

nin = nargin;
narginchk(1,3);

if nin == 1
    omega = [];
    polelist = [];
elseif nin==2
    polelist = [];
end

rb = cell2mat(P.RateBounds(:,2));
rb = min(abs(rb(:)));
% Check if system is uncertain
if isuncertain(P)
    S.type = '.';
    S.subs = 'NominalValue';
    P = subsref(P,S);
end

if isinf(rb) && nin ==1
    % LMI Call  only done if ratebounds are not specified.
    % XXX fix outputs - this is missing lowerbound on lpvnorm (see lpvwcgain)
    [gamopt,Info] = local_lpvnorm(P);    
else
    % IQC Call
    % XXX fix outputs - Info is actually being populated by wcu right now.
    % XXX fix flag setting
    lbflag = false;
    [gamopt,~,Info] = lpvwcgain(P,omega,polelist,lbflag);  
    gamopt = gamopt.UpperBound;
end

end


function [gamopt,Info] = local_lpvnorm(P)

%%
% Initialize Outputs
Xopt = [];
gamopt = [];
Info = [];

%%
% Problem Dimensions
% XXX PJS The code below pulls out all blocks and then checks to see
% if they are all tvreal.  It would be easier to pull out only parameters
% and then check that the resulting M is nominal (not uncertain).
szP = iosize(P);
ny = szP(1);
nu = szP(2);

[M,DELTA,BLKSTRUCT,NORMUNC] = lftdata(P,[],'All');
Nblk = length(BLKSTRUCT);
nri = zeros(Nblk,1);
for i=1:Nblk
    if ~strcmp( BLKSTRUCT(i).Type , 'tvreal' )
        error('Plant may only have tvreal blocks');
    else
        % XXX PJS: This assumes that BLKSTRUCT and NORMUNC have the
        % same ordering for the blocks (with repeats in NORMUNC)
        nri(i) = BLKSTRUCT(i).Occurrences;
    end
end
nr = sum(nri);
nx = order(M);

%% Balance the system
% % XXX Balance I/O channels associated with parameters?
% Min = M;
% M = ssbal(M);

%%
% Unpack Data
[A,B,C,D] = ssdata(M);

%%
% Compute L2 gain bound
% gam is the induced L2 gain bound assuming the parameter range is
% normalized to [-1,1].

% Set up the LMIs
setlmis([]);

X = lmivar(1,[nx 1]);
if nr>0
    J = lmivar(1,[nri(:) ones(Nblk,1)]);
end
[gam2,ndec] = lmivar(1,[1 1]);

% Bounded Real LMI
tmp = [eye(nx) zeros(nx,nr+nu)];
lmiterm([1 1 1 X],[A'; B'],tmp,'s');
tmp = [zeros(nx+nr,nu); eye(nu)];
lmiterm([1 1 1 gam2],-tmp,tmp');
tmp = [C';D']*[zeros(nr,ny); eye(ny)];
lmiterm([1 1 1 0],tmp*tmp');
if nr>0
    tmp = [zeros(nx,nr); eye(nr); zeros(nu,nr)];
    lmiterm([1 1 1 J],-tmp,tmp');
    tmp = [C';D']*[eye(nr); zeros(ny,nr)];
    lmiterm([1 1 1 J],tmp,tmp');
end

% Xlb*I < X < Xub*I
lmiterm([-2 1 1 X],1,1);
cnt = 3;
% if opt.Xlb>0
%     lmiterm([2 1 1 0],opt.Xlb*eye(nx));
% end
% if isfinite(opt.Xub)
%     lmiterm([-cnt 1 1 0],opt.Xub*eye(nx));
%     lmiterm([cnt 1 1 X],1,1);
%     cnt = cnt+1;
% end

% Jlb*I < J < Jub*I
if nr> 0
    lmiterm([-cnt 1 1 J],1,1);
    cnt = cnt+1;
end
% if nr> 0 && opt.Jlb>0
%     lmiterm([cnt 1 1 0],opt.Jlb*eye(nr));
% end
% if nr> 0 && isfinite(opt.Jub)
%     lmiterm([-cnt 1 1 0],opt.Jub*eye(nr));
%     lmiterm([cnt 1 1 J],1,1);
%     cnt = cnt+1;
% end

% Set Objective:
cobj = zeros(ndec,1);
cobj(end) = 1;

% Get LMI Options
% TODO PJS 5/29/2011: Default options for other solvers?
%XXX
opt = lpvsynOptions;
if ~isempty(opt.SolverOptions)
    LMIopt = opt.SolverOptions;
elseif isequal(opt.Solver,'lmilab')
    % Default settings for LMI Lab
    LMIopt = zeros(5,1);
    LMIopt(5) = 1;        % Toggle display
else
    LMIopt = [];
end


% Solve LMI
lmisys = getlmis;
FeasFlag = 1;
if isequal(opt.Solver,'lmilab')
    [copt,xopt] = mincx(lmisys,cobj,LMIopt);
    if isempty(copt)
        FeasFlag = 0;
    end
elseif isequal(opt.Solver,'sedumi')
    % Convert to Sedumi format
    % TODO PJS 5/29/2011: Currently ignores x0
    [F0,Fi,blk] = lmitrans(lmisys);
    K.s = blk;
    [xdual,xopt,info]=sedumi(Fi,-cobj,F0,K,LMIopt);
    copt = cobj'*xopt;
    if info.pinf==1 || info.dinf==1 || info.numerr~=0
        % TODO PJS 5/29/2011: Also check info.feasratio?
        FeasFlag = 0;
    end
else
    % TODO PJS 5/20/2011: Implement other solvers with LMITRANS
    error('Specified solver is currently not available.');
end

% Handle Infeasible LMI Case
% TODO PJS 5/20/2011: What should we return in this case?
if ~FeasFlag
    Xopt = [];
    gamopt = inf;
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'J',[]);
    return;
end

gamopt=sqrt(dec2mat(lmisys,xopt,gam2));
Xopt=dec2mat(lmisys,xopt,X);
if nr>0
    Jopt=dec2mat(lmisys,xopt,J);
else
    Jopt = [];
end

Info = struct('Xopt',Xopt,'xopt',xopt,'copt',copt,'lmisys',lmisys,'Jopt',Jopt);
end
