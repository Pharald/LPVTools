function m = subsasgn(m,L,RHS)
% SUBSASGN  Subscripted assignment for TVREAL objects.
%
% See also: subsasgn, subsref.

switch L(1).type
    case '.'
        try
            % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
            if length(L) == 1
                tmp = RHS;
                if isstruct(m)
                    % Handle following case:  m.A = B
                    % where m is a structure with A and/or B as a TVREAL.
                    % This happens when doing a subsasgn on the Parameter
                    % property of a PLFTMAT/PLFTSS.
                    m.(L.subs) = RHS;
                    return;
                end
            else
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
            
            % Case-insentive match for publically setable/getable fields
            % Handle ambiguity caused by 'Range' and 'RateBounds'.
            N = L(1).subs;
            idx = find( strncmpi(N,{'Name';'Range';'RateBounds'}, ...
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
                        m = tvreal(tmp,m.Data.Range,m.RateBounds{2});
                    case 2
                        m = tvreal(m.Data.Name,tmp,m.RateBounds{2});
                    case 3
                        m.RateBounds{2} = tmp;
                end
                
            end
            
        catch
            error(lasterr);
        end
        
    case '()'
        % ()-reference is not supported for UREAL.
        % By extension it is not supported for TVREAL.
        error('Subscripted assignment is not supported for class "tvreal"');
        
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end

