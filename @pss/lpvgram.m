function W = lpvgram(P,type,Wgt,solflag)
% LPVGRAM   Compute Gramians for PSS objects
%
% Wc = LPVGRAM(SYS,'c') computes the controllability gramian of the PSS SYS.
% The output Wc is a constant DOUBLE matrix, which satisfies the 
% LMI: A*Wc+Wc*A' +B*B' < 0 at each point in the domain of SYS, where A is 
% the state matrix of SYS and B is its input matrix.
%
% Wo = LPVGRAM(SYS,'o') computes the observability gramian of the PSS SYS. 
% The output Wo is a constant DOUBLE matrix, which satisfies the 
% LMI: A'*Wo+Wo*A +C'*C < 0 at each point in the domain of SYS, where A is 
% the state matrix of SYS and C is its output matrix.
%
% W = LPVGRAM(SYS,OPTION,WEIGHT) applies a matrix weighting WEIGHT when 
% solving for the gramian. For a controllabilty gramian the LMI becomes: 
%                 A*WEIGHT*Wc+Wc*WEIGHT*A' +B*B' < 0
% For a observability gramian the LMI becomes: 
%                 A'*WEIGHT*Wo+Wo*WEIGHT*A +C'*C < 0
% If no WEIGHT is specified a default value of eye(size(A)) is used, and
% the resulting gramian is diagonal.
%
% W = LPVGRAM(SYS,...,INVERT) provides an alternative implementation of the
% algorithm which solves for the gramians. If INVERT is True the LMI
% conditions are changed to solve for the inverse of the gramians, which
% can improve the accuracy of the solution for certain systems. The default 
% implementation assumes INVERT=FALSE.
%
% See also: gram, lpvbalreal, lpvbalancmr.  

% XXX Add error checking once input options have been figured out.

% AH 9/10/14 page 169 in Woods thesis reccommends weighting by LTI gramians 
% Should we implement this? Should this be default action?

if isequal(type,'c')
    cflag = true;
elseif isequal(type,'o')
    cflag = false;
else
    error('Gramian type must be set to ''c'' or ''o'' ')
end


% Get state space data
[A,B,C,D] = ssdata(P);
nX = size(A,1);
nU = size(B,2);
nY = size(C,1);

if nargin== 2
    Wgt = eye(nX);
    solflag = true;
elseif nargin == 3    
    if islogical(Wgt)
        solflag = Wgt;
    else       
        solflag = true;
    end
end
Wgt = 0.5*(Wgt+Wgt');

% Create LMI variables
setlmis([]);

if solflag
    % Solve for the gramian
    [Xdec,ndec,Xvar] = lmivar(1,[nX 1]);
else
    % Solve for the inverse of the gramian (eg 7.63 & 7.64 in G Wood thesis
    [Xidec,ndec,Xivar] = lmivar(1,[nX 1]);
end
% Create lmis at each point in the domain
Dom = P.DomainPrivate;
PIVName = P.DomainPrivate.IVName;
npts = prod(Dom.LIVData);
cnt = 1;
for i=1:npts
    % Get parameter value
    if isempty(Dom)
        pvaluec = cell(0,1);
    else
        pvaluec = num2cell(Dom(i));  % single index into RGRID gives value
    end
    
    % Evaluate state matrices at parameter
    % TODO PJS 5/16/2011: Replace with LPVSUBS? LPVSUBS returns data as a
    % double but it does an interpolation rather than a pure extraction.
    AV = double(lpvsplit(A,PIVName,pvaluec));
    if cflag
        BV = double(lpvsplit(B,PIVName,pvaluec));
    else
        CV = double(lpvsplit(C,PIVName,pvaluec));
    end
        
    if solflag
        % Define gramian LMI
        if cflag
            % LMI: AX+XA' +BB' < 0
            lmiterm([cnt 1 1 Xdec],AV,1,'s');
            lmiterm([cnt 1 1 0],BV*BV');
        else
            % LMI: A'X+XA +C'C < 0
            lmiterm([cnt 1 1 Xdec],AV',1,'s');
            lmiterm([cnt 1 1 0],CV'*CV);
        end
    else 
        % Define inverse gramian LMI
        if cflag
            % LMI: XiA+A'Xi +XiBB'Xi < 0
            lmiterm([cnt 1 1 Xidec],1,AV,'s');
            lmiterm([cnt 1 2 Xidec],1,BV);
            lmiterm([cnt 2 2 0],-eye(nU));
        else
            % LMI: XiA'+AXi +XiC'CXi < 0
            lmiterm([cnt 1 1 Xidec],AV,1,'s');
            lmiterm([cnt 1 2 Xidec],1,CV');
            lmiterm([cnt 2 2 0],-eye(nY));
        end
    end        
    cnt = cnt+1;
    
end

% 0 < X 
if solflag
    lmiterm([-cnt 1 1 Xdec],1,1);
else
    lmiterm([-cnt 1 1 Xidec],1,1);
end
cnt = cnt+1;

% Objective: min trace(X) = min c'xvar
if solflag
    cobj = zeros(ndec,1);
    cobj(diag(Xvar)) = diag(Wgt);
    L = tril(Xvar,-1);
    ldx = find(L);
    cobj(L(ldx)) = 2*Wgt(ldx);
else
    cobj = zeros(ndec,1);
    cobj(diag(Xivar)) = diag(Wgt);
    L = tril(Xivar,-1);
    ldx = find(L);
    cobj(L(ldx)) = 2*Wgt(ldx);
    cobj = -cobj;
end
% % Get LMI Options
% XXX - Retain this? For now specify options as empty
LMIopt = [1e-6 0 0 0 1];
% % TODO PJS 5/29/2011: Default options for other solvers?
% if ~isempty(opt.SolverOptions)
%     LMIopt = opt.SolverOptions;
% elseif isequal(opt.Solver,'lmilab')
%     % Default settings for LMI Lab
%     LMIopt = zeros(5,1);
%     %LMIopt(1) = 1/(gmax-gmin); % Tol setting in old code
%     if isequal(Method,'MaxFeas')
%         LMIopt(2) = 40;   % Max # of iters for MaxFeas problem
%     elseif RateBndFlag
%         %LMIopt(2) = 600;  % Setting in old LPVOFSYN1
%         LMIopt(2) = 250;  % Max # of iters for rate bounded syn
%     else
%         LMIopt(2) = 150;  % Max # of iters for non-rate bounded syn
%     end
%     LMIopt(5) = 1;        % Toggle display
% else
%     LMIopt = [];
% end

% Solve LMI
lmisys = getlmis;
[copt,xopt] = mincx(lmisys,cobj,LMIopt,[],-1e10);

% XXX - Allow user to specify solver?
% if isequal(opt.Solver,'lmilab')
%     [copt,xopt] = mincx(lmisys,cobj,LMIopt,x0);
%     if isempty(copt)
%         FeasFlag = 0;
%     end
% elseif isequal(opt.Solver,'sedumi')
%     % Convert to Sedumi format
%     % TODO PJS 5/29/2011: Currently ignores x0
%     [F0,Fi,blk] = lmitrans(lmisys);
%     K.s = blk;
%     [xdual,xopt,info]=sedumi(Fi,-cobj,F0,K,LMIopt);
%     copt = cobj'*xopt;
%     if info.pinf==1 || info.dinf==1 || info.numerr~=0
%         % TODO PJS 5/29/2011: Also check info.feasratio?
%         FeasFlag = 0;
%     end
%     
% else
%     % TODO PJS 5/20/2011: Implement other solvers with LMITRANS
%     error('Specified solver is currently not available.');
% end
if solflag
    W = dec2mat(lmisys,xopt,Xdec);
else
    Wi = dec2mat(lmisys,xopt,Xidec);
    W = inv(Wi);
end
% W = pmat(W,P.DomainPrivate);