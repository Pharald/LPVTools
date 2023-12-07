function [m,delta,blkstruct,normunc] = lftdata(A,list,opt)
% LFTDATA Decomposes the uncertainty and/or parameters from a PLFT
%
% A PLFT object (PLFTSS or PLFTMAT) has a rational dependence on TVREAL
% parameters and (possibly) uncertainty.  The PLFT objects are modeled
% as a fixed matrix/LTI system M in feedback with a normalized, block-
% diagonal matrix GAMMA=diag(DELTA,RHO) of uncertainties DELTA and
% parameters RHO.  LFTDATA extracts either the DELTA and/or RHO from M.
%
% [MRHO,DELTA,BLKSTRUCT,NORMUNC] = lftdata(A) takes a PLFT A and extracts
% MRHO and DELTA such that A=LFT(DELTA,MRHO) where DELTA is the uncertainty
% and MRHO is a (not uncertain) PLFT with parameter dependence but no
% uncertainty. BLKSTRUCT and NORMUNC specify the block structure and
% normalized uncertainty as detailed in the help for uss/lftdata.
%
% [MRHO,DELTA,BLKSTRUCT,NORMUNC] = lftdata(A,LIST) only pulls out the
% uncertainty elements specified in LIST (given as a cell array of
% uncertainty names).
%
% [MDELTA,RHO,BLKSTRUCT,NORMUNC] = lftdata(A,LIST,'Parameters') takes a
% PLFT A and extracts MDELTA and RHO such that A=LFT(RHO,MDELTA) where
% RHO is a (not uncertain) parameter dependent block and MDELTA has no
% parameter dependence but is possibly uncertain. LIST specifies the
% parameters to be extracted. If LIST is empty, all parameters are extracted.
%
% [M,GAMMA,BLKSTRUCT,NORMUNC] = lftdata(A,LIST,'All') takes a PLFT A and
% extracts M and GAMMA such that A=LFT(GAMMA,M) where M is a fixed
% matrix/LTI system and GAMMA is both uncertain and parameter dependent.
% LIST may be empty to specify that all parameters and uncertainties
% should be extracted.
%
% See also lftdata, lft.

% XXX PJS: Do we need to document the 'Uncertainties' option since
% this is the default behavior?
%
% [MRHO,DELTA,BLKSTRUCT,NORMUNC] = lftdata(A,LIST,'Uncertainties') is
% equivalent to the default behavior given by
% [MRHO,DELTA] = lftdata(A,LIST).


% Check # of I/O arguments
nin = nargin;
nout = nargout;
narginchk(1, 3)
nargoutchk(0, 4)

% Set default inputs
if nin==1
    list = cell(0,1);
    opt = 'Uncertainties';
elseif nin==2
    opt = 'Uncertainties';
end

% Get underlying umat data, rate bound info, and param/unc names
AData = A.Data;
RateBounds = A.RateBounds;
Params = RateBounds(:,1);
Unc = setdiff( fieldnames(AData.Uncertainty), Params );

% Determine parameters/uncertainties to extract
if strncmp(opt,'All',1)
    if isempty(list)
        newlist = [Params; Unc];
    else
        newlist = list;
    end
elseif strncmp(opt,'Uncertainties',1)
    if isempty(list)
        newlist = Unc;
    else
        newlist = intersect(list,Unc);
    end
elseif strncmp(opt,'Parameters',1)
    if isempty(list)
        newlist = Params;
    else
        newlist = intersect(list,Params);
    end
else
    error(['Option specified by third input must be ''All'',' ...
        ' ''Parameters'', or ''Uncertainties''.']);
end

% Extract data
[m,delta,blkstruct,normunc] = lftdata(AData,newlist);

% Re-package M as PLFT if it is parameter dependent
% m will be returned as a DOUBLE/SS/UMAT/USS if it is not param dependent.
if isuncertain(m)
    % Determine uncertainties / parameters remaining in m
    if isa(m,'umat') || isa(m,'uss')
        U = fieldnames(m.Uncertainty);
    else
        U = cell(0,1);
    end
    mRB = RateBounds;
    [~,idx]=setdiff(mRB(:,1),U);
    mRB(idx,:)=[];
    
    if ~isempty(mRB)
        % m is parameter dependent
        switch class(m)
            case 'umat'
                m = plftmat(m,mRB);
            case 'uss'
                m = plftss(m,mRB);
        end
    end
end

% Re-package DELTA as PLFT if it is parameter dependent
% DELTA will be returned as a UMAT/USS if it is not param dependent.
dRB = RateBounds;
U = fieldnames(delta.Uncertainty);
[~,idx]=setdiff(dRB(:,1),U);
dRB(idx,:)=[];
if ~isempty(dRB)
    % Delta is parameter dependent
    switch class(delta)
        case 'umat'
            delta = plftmat(delta,dRB);
        case 'uss'
            delta = plftss(delta,dRB);
    end
end


% Update blkstruct to include tvreal info
Nblk = length(blkstruct);
for i=1:Nblk
    blkname = blkstruct(i).Name;
    idx = find( strcmp(blkname,RateBounds(:,1)) );
    if ~isempty(idx)
        blkstruct(i).Type = 'tvreal';
    end
end

% Update normunc to include tvreal info
% Note: Block names in normunc are appended with 'Normalized'.
% XXX PJS: The time-derivative of the normalized parameter depends on
% both the value and rate of the (unormalized) parameter. The code below
% computes the largest rate bounds for the normalized parameter over
% the entire range of values. Should we enforce "centered" tvreals?
Nblk = length(normunc);
for i=1:Nblk
    blkname = normunc{i}.Name(1:end-10);
    idx = find( strcmp(blkname,RateBounds(:,1)) );
    if ~isempty(idx)
        % Get ureal
        U = A.Data.Uncertainty.(blkname);
        URange = U.Range;
        URateBounds = RateBounds{idx,2};
        
        % Get normalizer N:
        %    U = Fu(N,D) where D is normalized and U is unnormalized
        % XXX PJS: Is there an easier way to get the LFT normalizer N?
        N = lftdata(U);
        detN = det(N);
        N11 = N(1,1);
        N12 = N(1,2);
        N21 = N(2,1);
        
        % Rates are related by:
        %   Ddot = (N12*N21)*Udot / ( U*N11-detN )^2
        % Find largest range of rates for normalized parameter D
        % XXX PJS: Double check this. I believe the min/max values of the
        % denominator must occur on the interval edges.
        d1 = min( (URange(1)*N11-detN)^2 , (URange(2)*N11-detN)^2 );
        d2 = max( (URange(1)*N11-detN)^2 , (URange(2)*N11-detN)^2 );
        
        RB11 = (N12*N21)*URateBounds(1)/d1;
        RB12 = (N12*N21)*URateBounds(1)/d2;
        RB21 = (N12*N21)*URateBounds(2)/d1;
        RB22 = (N12*N21)*URateBounds(2)/d2;
        
        RBHigh = max( [RB11,RB12,RB21,RB22] );
        RBLow = min( [RB11,RB12,RB21,RB22] );
        
        normunc{i} = tvreal([blkname 'Normalized'],[-1 1],[RBLow, RBHigh]);
    end
end

