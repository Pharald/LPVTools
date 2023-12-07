function [m,delta,blkstruct,normunc] = lftdata(A,list)
% LFTDATA   Decomposes a UPSS.
% 
% [M,DELTA] = lftdata(A) decomposes an uncertain parameter varying system
% into two blocks, M and Delta, such that A = LFT(DELTA,M).
% Delta is a normalized, block-diagonal matrix of uncertain elements, while
% M is a PSS.
%
% [M,DELTA] = lftdata(A,LIST) only pulls out the uncertainty elements 
% specified in LIST (given as a cell array of uncertainty names). DELTA 
% contains only the uncertainties listed in LIST, while M is a UPSS which 
% retains the others. 
% 
% [M,DELTA,BLKSTRUCT] = lftdata(A) returns a N-by-1 structure array
% BLKSTRUCT, such that BLKSTRUCT(i) describes the i-th normalized uncertain 
% element.
% 
% [M,DELTA,BLKSTRUCT,NORMUNC] = lftdata(A) returns the cell array NORMUNC
% of normalized uncertain elements.
% 
% See also lftdata, lft.

% Get underlying umat data and rate bound info
narginchk(1,2)
AData = A.Data;
ADomain = A.Domain;

% Extract lftdata from umat
if nargin==1
    [m,delta,blkstruct,normunc] = lftdata(AData);
else
    [m,delta,blkstruct,normunc] = lftdata(AData,list);
end

% Re-package m as an grid based p object
switch class(m)
    case 'double'
        m = pmat(m,ADomain);
    case 'ss'
        m = pss(m,ADomain);
    case 'frd'
        m = pfrd(m,ADomain);
    case 'umat'
        m = upmat(m,ADomain);
    case 'uss'
        m = upss(m,ADomain);
    case 'ufrd'
        m = upfrd(m,ADomain);
    otherwise
        % XXX PJS: Update error to provide more info.
        error('Unsupported operation.');
end

%XXX lftdata on uncertain objects is not applicable to model arrays with
% varying uncertainty descriptions.Hence, Uncertainty can't vary across the
% grid.

% switch class(delta)
%     case 'umat'
%         delta = upmat(delta,ADomain);
%     case 'uss'
%         delta = upss(delta,ADomain);
%     case 'ufrd'
%         delta = upfrd(delta,ADomain);
%     otherwise
%         % XXX PJS: Update error to provide more info.
%         error('Unsupported operation.');
% end

% switch class(normunc)
%     case 'umat'
%         normunc = upmat(normunc,ADomain);
%     case 'uss'
%         normunc = upss(normunc,ADomain);
%     case 'ufrd'
%         normunc = upfrd(normunc,ADomain);
%     otherwise
%         % XXX PJS: Update error to provide more info.
%         error('Unsupported operation.');
% end



