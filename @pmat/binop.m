function C = binop(A,B,op,varargin)
% BINOP Binary operations for a PMAT
% 
% BINOP is a utility function that handles binary operations (e.g. times,
% ldivide, lft, etc.) on PMATs. The user should not call BINOP directly. 
% Instead each binary operation has a front end function that the user 
% calls (e.g. *,+,lft,horzcat). The frontend function then call calls BINOP 
% to do the required calculations/operations. 

%% PMAT/BINOP
% GJB 26Jan12 Added handling of uncertain objects

% NOTE PJS 4/7/2011: FEVAL is a bit slower than directly performing the 
% operation, eg.
%   if M is a matrix then M*M is about 4x faster than
%   feval('mtimes',M,M) and M*M is about 2x faster than
%   fh(M,M) where fh=@mtimes.
% The main BINOPs are implemented directly (e.g. plus calls double/+).


% Use switchyard to handle lifting if classes are different
if ~isequal(class(A),class(B))
   [A,B]=switchyard(A,B);
   nin = nargin;
   if nin==3
      C=binop(A,B,op);
   else
      C=binop(A,B,op,varargin{:});
   end
   return
end

% Define A and B on a common domain 
% (domunion will lift to pmat if required)
[Aext,Bext] = domunion(A,B);
szA = [privatesize(Aext) 1];
Adata = Aext.DataPrivate;
Bdata = Bext.DataPrivate;

% Initialize output Data (Use double/binop to find correct dims)
% NOTE PJS 4/7/2011: Some binops don't require this initialization,
% e.g. horzcat, and hence minor speed-ups are possible.
funchan = eval(['@' op]);
C1 = funchan(Adata(:,:,1), Bdata(:,:,1) , varargin{:} );
Cdata = zeros( [size(C1) szA(3:end)] );
Cdata(:,:,1) = C1;

% Perform binary operation at each point in the combined domain
switch op
    case 'horzcat'
        % TODO PJS 6/28/2011: [A,[]] will error if A is a nontrivial PMAT.
        Cdata = cat(2,Adata,Bdata);
    case 'vertcat'
        Cdata = cat(1,Adata,Bdata);
    case 'plus';
        % NOTE PJS 4/1/2011: for-loop properly handles scalar expansion,
        % e.g. A+B if A is a 2-by-2 PMAT and B is a 1-by-1 PMAT. 
        % Cdata = Adata+Bdata will error out for this scalar expansion
        % case if the the combined domain has more than 1 point.
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)+Bdata(:,:,i);
        end
        %Cdata = Adata+Bdata;
    case 'minus';
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)-Bdata(:,:,i);
        end
        %Cdata = Adata-Bdata;
    case 'mtimes'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)*Bdata(:,:,i);
        end
    case 'times'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i).*Bdata(:,:,i);
        end
        %Cdata = Adata.*Bdata;
    case 'mpower'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)^Bdata(:,:,i);
        end
    case 'power'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i).^Bdata(:,:,i);
        end
        %Cdata = Adata.^Bdata;
    case 'mldivide' 
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)\Bdata(:,:,i);
        end        
    case 'mrdivide'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)/Bdata(:,:,i);
        end
    case 'ldivide'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i).\Bdata(:,:,i);
        end
        %Cdata = Adata.\Bdata;
    case 'rdivide'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)./Bdata(:,:,i);
        end
        %Cdata = Adata./Bdata;
    case 'blkdiag'
        szB = [privatesize(Bext) 1];
        Z12 = zeros([szA(1),szB(2:end)]);
        Z21 = zeros([szB(1),szA(2:end)]);
        Cdata = [Adata Z12; Z21 Bdata];
    case 'starp'
        % TODO PJS: Is this obsolete if we implement lft binop?
        dim1 = varargin{1};
        dim2 = varargin{2};
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = starp(Adata(:,:,i),Bdata(:,:,i),dim1,dim2);
        end        
    case 'lft'
        % TODO PJS: Call LFTDOUBLE?
        % NOTE PJS: LFT documentation states that it automatically works
        %   on arrays of models. Do we need the for-loop here?
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = lft(Adata(:,:,i),Bdata(:,:,i),varargin{:});
        end        
    case 'eq'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)==Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'ne'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)~=Bdata(:,:,i);
        end
        Cdata = logical(Cdata);        
    case 'lt'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)<Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'gt'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)>Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'le'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)<=Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'ge'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)>=Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'and'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i) & Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'or'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i) | Bdata(:,:,i);
        end
        Cdata = logical(Cdata);
    case 'xor'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = xor(Adata(:,:,i), Bdata(:,:,i));
        end
        Cdata = logical(Cdata);
    otherwise
        % TODO PJS Put this in a try-catch?
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = funchan(Adata(:,:,i), Bdata(:,:,i), varargin{:});
        end
end
C = pmat(Cdata,Aext.DomainPrivate);

% --------  OLD CODE

% function out = binop(A,B,op,dim1,dim2)
%
%
% niva = A.DomainPrivate.NumIV;
% nivb = B.DomainPrivate.NumIV;
% % this should be a more general call (BINOP) to a domain
% % method.  that way, as we introduce new domains,
% % only their binop programs need to be written.
% [aidx,bidx,abrgrid] = domainbin(A.DomainPrivate,B.DomainPrivate);
% Adata = permute(A.DataPrivate,[1 2 2+aidx]);
% Bdata = permute(B.DataPrivate,[1 2 2+bidx]);
%
% switch op
% case 'mtimes'
%     mat = ndtimesr(Adata,Bdata);
% case 'horzcat'
%    % XXX We can't call built-in because they don't do scalar expansion
%    % mat = [Adata Bdata]
%    % mat = cat(1,Adata,Bdata);
%
%     mat = ndhorzr(Adata,Bdata);
% case 'vertcat'
%     mat = ndvertr(Adata,Bdata);
% case 'plus'
%     mat = ndplusr(Adata,Bdata);
% case 'minus'
%     mat = ndplusr(Adata,-Bdata);
% case 'starp'
%    szt = size(Adata);
%    szb = size(Bdata);
%    dimout1 = szt(1) - dim1;
%    dimin1 = szt(2) - dim2;
%    dimout2 = szb(1) - dim2;
%    dimin2 = szb(2) - dim1;
%    top11 = zeros([dimout1 dimin1 szt(3:end)]);
%    top12 = zeros([dimout1 szt(2)-dimin1 szt(3:end)]);
%    top21 = zeros([szt(1)-dimout1 dimin1 szt(3:end)]);
%    top22 = zeros([szt(1)-dimout1 szt(2)-dimin1 szt(3:end)]);
%    bot11 = zeros([dim2 dim1 szb(3:end)]);
%    bot12 = zeros([dim2 szb(2)-dim1 szb(3:end)]);
%    bot21 = zeros([szb(1)-dim2 dim1 szb(3:end)]);
%    bot22 = zeros([szb(1)-dim2 szb(2)-dim1 szb(3:end)]);
%    top11(:) = Adata(1:dimout1,1:dimin1,:);
%    top12(:) = Adata(1:dimout1,dimin1+1:end,:);
%    top21(:) = Adata(dimout1+1:end,1:dimin1,:);
%    top22(:) = Adata(dimout1+1:end,dimin1+1:end,:);
%    bot11(:) = Bdata(1:dim2,1:dim1,:);
%    bot12(:) = Bdata(1:dim2,dim1+1:end,:);
%    bot21(:) = Bdata(dim2+1:end,1:dim1,:);
%    bot22(:) = Bdata(dim2+1:end,dim1+1:end,:);
%    tb = ndtimesr(top22,bot11);
% %   [nrtb nctb] = size(tb);
%    imtb = ndplusr(eye(dim1),-tb);
% %  should check the invertibility of imtb at this point, but
% %  we have had trouble with COND on large, sparse matrices
%    bt = ndtimesr(bot11,top22);
%    %[nrbt,ncbt] = size(bt);
%    imbt = ndplusr(eye(dim2),-bt);
%    x = ndmldivr(imbt,bot12);
%    y = ndmldivr(imtb,top21);
%    upper1 = ndplusr(top11,ndtimesr(ndtimesr(top12,bot11),y));
%    upper2 = ndtimesr(top12,x);
%    lower1 = ndtimesr(bot21,y);
%    lower2 = ndplusr(bot22,ndtimesr(ndtimesr(bot21,top22),x));
%    mat = ndvertr(ndhorzr(upper1,upper2),ndhorzr(lower1,lower2));
% case 'mldivide'
%     mat = ndmldivr(Adata,Bdata);
% case 'mrdivide'
%     mat = ndmrdivr(Adata,Bdata);
% case 'blkdiag'
%     mat = ndbdiagr(Adata,Bdata);
% otherwise
%     error([op ' is not a valid command']);
% end
%
% out = pmat(mat,abrgrid);

%szorigm = size(mat);
%dnr = pvget(abrgrid,'DNR');
%niv = pvget(abrgrid,'niv');
%ivdata = pvget(abrgrid,'IvData');
%ivname = pvget(abrgrid,'IvName');
%tmp = ones(1,2+niv);
%repmatdata = ones(1,2+niv);
%mat = removend(mat);
%% Fix IvData and IvName if reduction occurred...
%szremm = size(mat);
%tmp(1:length(szremm)) = szremm;
%for i=1:niv
%   if dnr(i)==1 & tmp(2+i)==1
%      repmatdata(2+i) = length(ivdata{i});
%      tmp(2+i) = length(ivdata{i});
%   end
%end
%mat = repmat(mat,repmatdata);
%idx = find(tmp(3:2+niv)>1);
%if length(idx)>0
%   % only the actual IV dependence that remains is in the result
%   ridx = [1 2 idx+2];
%   mat = reshape(mat,tmp(ridx));
%   newdom = rgrid(ivname(idx),ivdata(idx),dnr(idx));
%   m = pmat(m.NdData,newdom,m.RCData);
%else

% Fix result into PMAT
% mat = removend(mat);
% szm = size(mat);
% idx = find(szm(3:end)>1);
% if length(idx)>0
%    IVN = abrgrid.IVName;
%    IVD = abrgrid.IVData;
%    newgrid = rgrid(IVN(idx),IVD(idx));
% else
%    newgrid = rgrid;
% end
% % at this point, it has been reduced without regard to DNR.  Any
% % DNR==1 infor in ABRGRID can be used to resurrect IVs that should
% % not have been reduced.
% [newgrid,rdata] = unreduce(newgrid,abrgrid);
% out = depromote(pmat(repmat(reshape(mat,szm([1 2 2+idx])),rdata),...
%    newgrid,rcmanagedata));
