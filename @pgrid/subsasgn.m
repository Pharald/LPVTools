function m = subsasgn(m,L,RHS)
% SUBSASGN  Subscripted assignment for PGRID objects.
%
% See also: subsasgn, subsref.

switch L(1).type
    case '.'
        try
            % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
            if length(L) == 1
                tmp = RHS;
%                 if isstruct(m)
%                     % Handle following case:  m.A = B
%                     % where m is a structure with A and/or B as a TVREAL.
%                     % This happens when doing a subsasgn on the Parameter
%                     % property of a PMAT/PSS/PFRD.
%                     error('XXX Do we need this for gridded objects');
%                     m.(L.subs) = RHS;
%                     return;
%                 end
            else
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
            
            % Case-insentive match for publically setable/getable fields
            % Handle ambiguity caused by 'Range' and 'RateBounds'.
            N = L(1).subs;
            idx = find( strncmpi(N,{'Name';'GridData';'RateBounds'; ...
                'Range'}, length(N)) );
            
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
                        m.Name = tmp;
                    case 2
                        if numel(tmp)==numel(m.GridData)
                            m.GridData = tmp;
                        else
                            error(['Number of points must remain the '...
                                'same when replacing the GridData.']);
                        end
                    case 3
                        m.RateBounds = tmp;
                    case 4
                        error('Range is not a setable property of a PGRID.');
                end
                
            end
            
        catch
            error(lasterr);
        end
        
    case '()'
        % ()-reference is not supported for UREAL.
        % By extension it is not supported for PGRID.
        error('Subscripted assignment is not supported for class "pgrid"');
        
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end

