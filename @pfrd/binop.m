function C = binop(A,B,op,varargin)
% BINOP Binary operations for a PFRD
% 
% BINOP is a utility function that handles binary operations (e.g. times,
% ldivide, lft, etc.) on PFRDs. The user should not call BINOP directly. 
% Instead each binary operation has a front end function that the user 
% calls (e.g. *,+,lft,horzcat). The frontend function then call calls BINOP 
% to do the required calculations/operations.

%% PFRD/BINOP

% TODO PJS 4/30/2011: Private function? Utility?

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
Freq = Adata.Frequency;

% Initialize output Data (Use ss/binop to find correct dims)
% TODO PJS Some binops don't require this initialization (e.g. horzcat)
funchan = eval(['@' op]);
C1 = funchan(Adata(:,:,1), Bdata(:,:,1) , varargin{:} );
Cdata = repsys(C1,[1 1 szA(3:end)]);
Cdata(:,:,1) = C1;

% Perform binary operation at each point in the combined domain
switch op
    case 'horzcat'
        Cdata = [Adata,Bdata];
    case 'vertcat'
        Cdata = [Adata; Bdata];
    case 'plus';
        % TODO PJS 5/1/2001: If A is a 2-by-2 PFRD and B is a 1-by-1 PFRD
        % then does Cdata = Adata+Bdata properly handle the scalar 
        % expansion? If not then we need to use a for-loop
        % (Compare with PMAT/BINOP and PSS/BINOP).
        Cdata = Adata+Bdata;
    case 'minus';
        Cdata = Adata-Bdata;
    case 'mtimes'
        Cdata = Adata*Bdata;
    case 'times'
        Cdata = Adata.*Bdata;
    case 'mpower'
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = Adata(:,:,i)^Bdata(:,:,i);
        end
    case 'power'
        Cdata = Adata.^Bdata;
    case 'mldivide'
        Cdata = Adata\Bdata;
    case 'mrdivide'
        Cdata = Adata/Bdata;
    case 'blkdiag'
        Cdata = blkdiag(Adata,Bdata);
    case 'append'
        Cdata = append(Adata,Bdata);
    case 'lft'
        % NOTE PJS: LFT documentation states that it automatically works
        %   on arrays of models. Do we need the for-loop here?
        %for i=2:prod(szA(3:end))
        %    Cdata(:,:,i) = lft(Adata(:,:,i),Bdata(:,:,i),varargin{:});
        %end
        Cdata = lft(Adata,Bdata,varargin{:});
    otherwise
        % TODO PJS Put this in a try-catch?
        for i=2:prod(szA(3:end))
            Cdata(:,:,i) = funchan(Adata(:,:,i), Bdata(:,:,i), varargin{:});
        end
end
C = pfrd(Cdata,Aext.DomainPrivate);

% --------  OLD CODE

% function out = binop(A,B,op,dim1,dim2)
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
%    case 'mtimes'
%       sys = Adata*Bdata;
%    case 'horzcat'
%       sys = [Adata Bdata];
%    case 'vertcat'
%       sys = [Adata; Bdata];
%    case 'plus'
%       sys = Adata+Bdata;
%    case 'minus'
%       sys = Adata-Bdata;
%    case 'lft'
%       sys = lft(Adata,Bdata,dim1,dim2);
%    case 'mldivide'
%       sys = Adata\Bdata;
%    case 'mrdivide'
%       sys = Adata/Bdata;
%    case 'blkdiag'
%       sys = blkdiag(Adata,Bdata);
%    otherwise
%       error([op ' is not a valid command']);
% end
% out = psys(sys,abrgrid);
%
