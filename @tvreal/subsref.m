function m = subsref(m,L)
% SUBSREF  Subscripted reference for TVREAL objects.
%
% See also: subsref, subsasgn.

switch L(1).type
    case '.'
        % Case-insentive match for publically setable/getable fields
        % Handle ambiguity caused by 'Range' and 'RateBounds'.
        N = L(1).subs;
        idx = find( strncmpi(N,{'Name';'Range';'RateBounds';'Data'}, ...
            length(N)) );
                
        if length(idx) == 2
            error(['The specified property name "' N '" is '...
                'ambiguous. Specify more characters for the '...
                'property name.']);
        end
        
        if isempty(idx)
            error(['No property of the class "tvreal" matches '...
                'the string "' N '". Use PROPERTIES to get the '...
                'list of properties for this class.']);            
        else         
            switch idx
                case 1
                    m = m.Data.Name;
                case 2
                    m = m.Data.Range;
                case 3
                    m = m.RateBounds{2};
                case 4
                    m = m.Data;
            end
            
        end
        
    case '()'
        % Perform subsref on data. Currently if s is a ureal then
        % s(1) or s(1,1) simply returns s as a UMAT
        mData = subsref(m.Data,L(1));
        mRateBounds = m.RateBounds;
        
        % Repackage final result as an plft object
        switch class(mData)
            case 'ureal'
                m = tvreal(mData.Name,mData.Range,mRateBounds);
            case 'umat'
                m = plftmat(mData,mRateBounds);
            case 'uss'
                m = plftss(mData,mRateBounds);
            case 'ufrd'
                m = plftfrd(mData,mRateBounds);
            otherwise
                % XXX PJS: Update error to provide more info.
                error('Unsupported ()-reference');
        end
        
    otherwise
        error(['{}-like reference is not supported for ' class(m) 'objects.'])
end

if length(L)>1
    m = subsref(m,L(2:end));
end
