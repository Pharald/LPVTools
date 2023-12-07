function S = IQCcell(c,PerfFlag)

% PerfFlag = false --> IQC for udyn (not performance)
% PerfFlag = true --> IQC for performance 
if nargin==1
    PerfFlag = false;
end

S.IQCfunction = @LOCALcell;
S.IQCparams.cell = c;
S.IQCparams.PerfFlag = PerfFlag; 
S.PsiFlag = false;

function [lmiout,IDcenter,Pi] = LOCALcell(lmisys,ncopies,omega,params)

c = params.cell;
PerfFlag = params.PerfFlag;

nIQCs = numel(c);
if ncopies>1
    error('only 1 copy for these typos of blocks');
end

setlmis(lmisys)
ndec = decnbr(lmisys);
nlmis = lminbr(lmisys);
IDcenter = zeros(nIQCs,1);
Nomega = numel(omega);
for j=1:nIQCs
    if (~PerfFlag) || (PerfFlag && nIQCs>1 && j==1)
        % For Perf IQC, create only 1 dec variable but it will
        % multiply the second perf IQC.
        IDcenter(j) = lmivar( 3, ndec+1);
        ndec = ndec+1;
        
        if (~PerfFlag) 
            % Conic combinations of IQCs require x>=0.
            % Don't need this constraint for perf vars            
            lmiterm([-(nlmis+1) 1 1 IDcenter(j)],1,1);
            nlmis=nlmis+1;
        end
    end
        
    if isa(c{j},'ss') || isa(c{j},'tf') || isa(c{j},'double')
        c{j}= freqresp(ss(c{j}),omega);
    elseif isa(c{j},'function_handle')
        fh = c{j};
        tmp = zeros([size(fh(omega(1))) Nomega]);
        for i = 1:Nomega
            tmp(:,:,i) = fh(omega(i));
        end
        c{j} = tmp;
    else
        error('Only ss, tf, double, function handle');
    end
end
lmiout = getlmis;
Pi = c;
