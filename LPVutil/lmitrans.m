% [F0,Fi,blk] = lmitrans(lmisys)
%
% This function converts from LMILAB to a canonical LMI form.
% The function input, LMISYS, is an LMILAB description of
% the LMI constraint, L(x) < R(x). x is an n dimensional vector
% of decision variables and L(x), R(x) are m by m symmetric
% matrices. The function outputs are the sparse matrix F0 and an
% n by 1 cell array of sparse matrices, FI. All matrices are
% m by m and symmetric.  They specify the same LMI constraint
% in the following canonical form:
%             x1*Fi{1} + .... + xn*Fi{n} < F0
% The matrices inherit a block diagonal structure which is given
% by the output BLK. Specifically, BLK is a vector which specifies
% the block dimensions of each LMI specified in LMISYS.
%
% Input:
%  LMISYS    LMILAB description of the system of LMI constraints
%
% Output:
%  F0        Constant matrix
%  FI        Cell array of matrices
%  BLK       A vector containing the block dimensions of the
%            LMIs given in LMISYS
%
% NOTE: Block structure contained in individual LMIs of LMISYS
%       is not reflected in BLK. Therefore it may be preferable
%       (for cpu time of some LMI solvers) to use 'lmiterm' to
%       specify each block with a separate LMI rather than one LMI
%       with many blocks on the diagonal.
%
% See also:

% Coding:
% 6/28/2001   Pete Seiler
% 5/29/2011   PJS  Store F0/Fi as sparse columns + other speedups

% TODO PJS 5/29/2011--
% This is still in a state of testing. This version
% avoids the use cell arrays to store Fi and is much faster than
% the previous version.  The outputs F0/Fi are stored in Sedumi
% format (i.e blocks are stacked on top of each other and can be
% unpacked from the blk dimension).
%   call:
%     [F0,Fi,blk] = lmitrans(lmisys);
%     K.s = blk;
%     [xdual,xopt,info]=sedumi(Fi,-cobj,F0,K);
%

function [F0,Fi,blk] = lmitrans(lmisys)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Is function called correctly?
if nargin ~=1 || nargout~=3
    error('usage: [F0,Fi,blk] = lmitrans(lmisys)');
end

% Is lmisys a valid LMI description from LMILAB?
if size(lmisys,1)<10 || size(lmisys,2)>1
    error('lmisys is not a valid LMI system');
elseif length(lmisys)~=10+lmisys(1)*lmisys(4)+6*lmisys(3)+...
        lmisys(2)*lmisys(5)+lmisys(7),
    error('lmisys is not a valid LMI system');
elseif any(lmisys(1:8)<0),
    error('lmisys is not a valid LMI system');
elseif any(~lmisys(1:3)),
    error('No matrix variable or term in this LMI system');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Get LMI information from lmisys and resolve
%  any undetermined block / outer factor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lmisys stores LMI information in 5 block pieces:
% 1)header  2)lmiset  3)lmiterm  4)lmivar  5)data
% LOCALlmidata splits out the last 4 blocks of lmisys
[lmiset,lmivar,lmiterm,data]=LOCALlmidata(lmisys);
%ldata=lmisys(7);        % length of data segment
%nlmi=size(lmiset,2);    % number of lmis
nvar=size(lmivar,2);    % number of variables defined using lmivar.m
ndec=lmisys(8);         % number of decision variables

% Get decision variable info
varinfo = cell(nvar,2);
for ii=1:nvar
    % Form entry-wise dependence on decision variables
    V = decinfo(lmisys,ii);
    varinfo{ii,1} =  V;
    
    % Find set of decision variables for this variable
    V = V( V~=0 );
    varinfo{ii,2} = unique( abs(V) );
end

% Resolve any undetermined block and/or outer factor dimensions
% This piece of the code is a modified version of the Gahinet /
% Nemirovskii code in the LMILAB function 'nnsetup.m'
i1=0;
for lmi=lmiset
    i1 = i1+1;
    lmilb = lmi(1);              % LMI number
    insize = lmi(3);             % inner block dimension
    blkdims = lmi(7:6+lmi(6));   % dimensions of LMI blocks
    
    % Fully determine the block sizes and inner size
    % Note: insize>0 only when forced by outer factor
    ind = find(blkdims<=0); % undetermined block sizes
    lneg = length(ind);
    if insize<=0          % inner size undetermined
        blkdims(ind)=ones(length(ind),1);
        insize=sum(blkdims);
        lmiset(3,i1)=insize;
    elseif lneg==1        % one scalar block
        if insize-sum(blkdims(blkdims>0)) <= 0
            error('LMI #%d: the inner and outer factors have incompatible dimensions',lmilb)
        end
        blkdims(ind) = insize-sum(blkdims(blkdims>0));
    elseif lneg,           % several scalar blocks -> ambiguous
        blkdims(ind) = ones(length(ind),1);
        if sum(blkdims) ~= insize
            error(['LMI #%d: ambiguous block dimensioning (multiple scalar terms)\n' ...
                '        please dimension the relevant coefficient matrices'],lmilb);
        end
    elseif sum(blkdims)~=insize
        error('LMI #%d: the inner and outer factors have incompatible dimensions',lmilb);
    end
    if lneg
        lmiset(6+ind,i1) = blkdims(ind);
    end
    if lmi(2)<=0
        lmiset(2,i1)=insize; % set outsize if undetermined
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert from LMILAB from to the canonical form
% specified above, i.e. find Fi and F0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form Fivec, F0vec, and blk
%Fi=cell(ndec,1); F0=sparse([]); blk=[];
Fi=sparse(0,ndec); F0=sparse([]); blk=[];
for lmi=lmiset
    lmilb = lmi(1);            % LMI Number
    outsize = lmi(2);          % Outer Dimensions of LMI
    insize = lmi(3);           % Inner Dimensions of LMI
    blkdims = lmi(7:6+lmi(6)); % Dimensions of LMI blocks
    
    % Grab columns of lmiterm associated with this lmi
    trmcols = lmi(4):lmi(5);
    idx = find( abs(lmiterm(1,trmcols)) == lmilb );
    term = lmiterm(:,trmcols(idx));
    
    % Form outer factors (if they exist) for left/right sides
    lf=1;
    rf=1;
    ofidx = find(term(2,:)==0);
    for i1=1:length(ofidx)
        lstart = term(5,ofidx(i1));
        lend = lstart+term(6,ofidx(i1));
        rA = data(lstart+1);
        cA=data(lstart+2);
        A = reshape( data(lstart+3:lend), [cA rA] )';
        
        if term(1,ofidx(i1)) > 0
            lf = A;  % left side outer factor
        else
            rf = A;  % right side outer factor
        end
    end % end: form outer factors
    
    % Form constant term, F0, for this LMI
    f0 = sparse(outsize,outsize);
    cnstidx = find( term(2,:)~=0 & term(4,:)==0 );
    for i1=1:length(cnstidx)
        lstart = term(5,cnstidx(i1));
        lend = lstart+term(6,cnstidx(i1));
        rA = data(lstart+1);
        cA=data(lstart+2);
        A = reshape( data(lstart+3:lend), [cA rA] )';
        
        row = term(2,cnstidx(i1));
        col = term(3,cnstidx(i1));
        if rA ==1 && cA ==1  % Expand scalar terms
            A = A*eye(blkdims(row));
        end
        
        interm = zeros(insize,insize);
        if row==col
            tempr = sum(blkdims(1:row-1)) + (1:blkdims(row));
            interm( tempr , tempr ) = A;
        else
            tempr = sum(blkdims(1:row-1)) + (1:blkdims(row));
            tempc = sum(blkdims(1:col-1)) + (1:blkdims(col));
            interm( tempr , tempc ) =A;
            interm = interm+interm';
        end
        
        if term(1,cnstidx(i1)) > 0
            f0 = f0 - sparse(lf'*interm*lf);  % left side constant term
        else
            f0 = f0 + sparse(rf'*interm*rf);  % right side constant term
        end
    end % end: form F0 term
    
    % Form Fi terms for this LMI
    fi = zeros(outsize^2,ndec);
    varidx = find( term(2,:)~=0 & term(4,:)~=0 );
    varIDs = unique( abs(term(4,varidx)) );
    for i1=1:length(varIDs)
        % Grab variable data from lmivar
        temp = find( lmivar(1,:) == varIDs(i1) );
        var = lmivar(:,temp);
        
        % Get entry-wise dependence on dec. vars Xdec
        % Also get  unique list of dec vars in Xdec, decvars
        Xdec = varinfo{var(1),1};
        decvars = varinfo{var(1),2};
        
        % Find terms which use this variable
        varidx = find( term(2,:)~=0 & abs(term(4,:))==varIDs(i1) );
        for i2=1:length(varidx)
            lstart = term(5,varidx(i2));
            lend = lstart+term(6,varidx(i2));
            rA = data(lstart+1);
            cA=data(lstart+2);
            temp = lstart+2+rA*cA;
            A = reshape(data(lstart+3:temp),[cA rA])';
            
            rB = data(temp+1);
            cB = data(temp+2);
            B = reshape(data(temp+3:lend-1),[cB rB])';
            bool = data(lend);
            
            if term(4,varidx(i2))<0
                Xdec2 = Xdec';   %Term stored is A*X'*B
            else
                Xdec2 = Xdec;    %Term stored is A*X*B
            end
            
            row = term(2,varidx(i2));
            col = term(3,varidx(i2));
            tempr = sum(blkdims(1:row-1))+ (1:blkdims(row));
            tempc = sum(blkdims(1:col-1))+ (1:blkdims(col));
            for i3 = 1:length(decvars)
                % find inside term
                Xk = (Xdec2==decvars(i3)) - (Xdec2==-decvars(i3));
                
                axb = A*Xk*B;
                sz = size(axb);
                if sz(1)==1 && sz(2)==1 % Expand scalar terms
                    axb = axb*eye(blkdims(row));
                end
                
                % Note: It is faster to leave interm as full and only do
                % the sparse conversion when adding fi data to Fi.
                interm = zeros(insize,insize);
                if row==col
                    if bool
                        interm( tempr , tempr ) = axb;
                    else
                        interm( tempr , tempr ) = axb + axb';
                    end
                else
                    interm( tempr , tempc ) =axb;
                    interm = interm+interm';
                end
                
                % multiply by outer factors and store in ficols
                if ~isequal(lf,1) %lf~=1
                    fiterm = sparse(lf'*interm*lf);
                else
                    fiterm = interm;
                end
                if term(1,varidx(i2)) > 0
                    % left side
                    fi(:,decvars(i3)) = fi(:,decvars(i3)) + fiterm(:);
                else
                    % right side
                    fi(:,decvars(i3)) = fi(:,decvars(i3)) - fiterm(:);
                end
            end
        end
    end % end: form Fi terms
    
    % Add data to list
    F0 = [F0; f0(:)];
    Fi = [Fi; sparse(fi)];
    blk = [blk , outsize];
end % end: form data




%-------------------------------------------------------------
% Separate out LMI data
function [lmiset,lmivar,lmiterm,data]=LOCALlmidata(lmisys)

% Header
header = lmisys(1:10);
nlmi=header(1);             % # of lmis
nvar=header(2);             % # of vars
nterm=header(3);            % # of terms
rs=header(4);               % row sizes of lmiset and lmivar
rv=header(5);               % row sizes of lmiset and lmivar
ls=nlmi*rs;                 % length of lmiset
lv=nvar*rv;                 % length of lmivar
lt=6*nterm;                 % length of lmiterm

% lmiset
ptr=10;
lmiset=reshape(lmisys(ptr+1:ptr+ls),rs,nlmi);

% lmivar
ptr=ptr+ls;
lmivar=reshape(lmisys(ptr+1:ptr+lv),rv,nvar);

% lmiterm
ptr=ptr+lv;
lmiterm=reshape(lmisys(ptr+1:ptr+lt),6,nterm);

% data
ptr=ptr+lt;
data = lmisys(ptr+1:end);

