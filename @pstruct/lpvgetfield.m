function out = lpvgetfield(obj,varargin)
% LPVGETFIELD  Get PSTRUCT field contents.
%
% B = LPVGETFIELD(M,FIELDNAME) returns the contents of the field FIELDNAME
% in the PSTRUCT M. LPVGETFIELD(M,FIELDNAME) is equivalent to M.FIELDNAME 
% for STRUCT objects. B is a parameter-varying object corresponding to the 
% type of data in the requested field.
%
% See also: fieldnames.

% Determine class of field content
% Determine if size of field content is consistent across IV/array
% Assemble the
szstruct = size(obj.Data);
numIV = obj.Domain.NumIV;
ptsIV = size(obj.Domain);

tempfirst = getfield(obj.Data(1),varargin{:});
out(:,:,1) = tempfirst;
szval = size(tempfirst);
clval = class(out);
cellout = false;

if prod(szstruct)>1
    
    for k = 2:prod(szstruct)
        
        temp = getfield(obj.DataPrivate(k),varargin{:});
        if ( all(class(temp) == clval) ) && ( all(size(temp)  ==  szval))
            out(:,:,k) = temp;
        else
            cellout = true;
            break
        end
        
    end
    
    if cellout == true || ischar(tempfirst)
        % Return a comma seperated list.
        out = cell(prod(szstruct),1);
        out{1} = tempfirst;
        for k = 2:prod(szstruct)
            temp = getfield(obj.Data(k),varargin{:});
            out{k} = temp;
        end
    else
        
        % Data in Pstruct is ordered [AD(1) AD(2) IVs AD(3:end)]
        % The first two indices correspond to the first two array dimensions.
        % The following indices correspond to the IVs, and the rest of the
        % indices correspond to the remaining array dimensions.
        % Prepare to reorder data as [IV AD]
        ord = [ 2+[1:numIV] 1 2 2+numIV+1:numel(szstruct)];
    
        switch lower(clval)
            case {'double'}
                out = reshape(out,[szval szstruct]);
                out = permute(out,[1 2 2+ord]);
                out = pmat(out,obj.Domain);
            case {'umat'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = upmat(out,obj.Domain);
            case {'frd'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = pfrd(out,obj.Domain);
            case {'ufrd'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = upfrd(out,obj.Domain);
            case {'ss'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = pss(out,obj.Domain);
            case{'uss'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = upss(out,obj.Domain);
            case{'struct'}
                out = reshape(out,szstruct);
                out = permute(out,ord);
                out = pstruct(out,obj.Domain);
            case{'char'}
                error(['Data type not accessable as a parameter varying object.\n'...
                    'Access data through manual input of grid index.'])
        end
    end
end
%
% try
%
%     out =
% catch
%     error(' This part still in works')
% end
%
%
%
% end
%
%
%
%
%
% size(WCG.DataPrivate) = 2     1     3     2
% 1x1 pstruct with i IVs that are 3 and 2 point,
% WCG.DataPrivate(1).LowerBound
% WCG.DataPrivate(2).LowerBound