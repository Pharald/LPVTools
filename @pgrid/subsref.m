function m = subsref(m,L)
% SUBSREF  Subscripted reference for PGRID objects.
%
% See also: subsref, subsasgn.

switch L(1).type
    case '.'
        % Case-insentive match for publically setable/getable fields
        % Handle ambiguity caused by 'Range' and 'RateBounds'.
        N = L(1).subs;
        idx = find( strncmpi(N,{'Name';'GridData';'RateBounds';'Range'}, ...
            length(N)) );
                
        if length(idx) == 2
            error(['The specified property name "' N '" is '...
                'ambiguous. Specify more characters for the '...
                'property name.']);
        end
        
        if isempty(idx)
            error(['No property of the class "pgrid" matches '...
                'the string "' N '". Use PROPERTIES to get the '...
                'list of properties for this class.']);            
        else         
            switch idx
                case 1
                    m = m.Name;
                case 2
                    m = m.GridData;
                case 3
                    m = m.RateBounds;
                case 4
                    GD = m.GridData;
                    m = [GD(1) GD(end)];
            end
            
        end
        
    case '()'
        % Currently if s is a ureal then s(1) or s(1,1) simply returns s 
        % as a UMAT. Mimic this behavior.
        if isequal( L(1).subs, {1} ) || isequal( L(1).subs, {1,1} )
            m = pmat( shiftdim(m.GridData,-2), rgrid(m) );
        else
            error ('Subscript is out of range.');
        end        
    otherwise
        error(['{}-like reference is not supported for ' class(m) 'objects.'])
end

if length(L)>1
    m = subsref(m,L(2:end));
end
