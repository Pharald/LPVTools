function m = subsasgn(m,L,RHS)
% SUBSASGN  Subscripted assignment for PLFTMAT objects.
%
% See also: subsasgn, subsref.

switch L(1).type
    case '.'
        try
            % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
            if length(L) == 1
                tmp = RHS;
            else
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
            
            % Case-insentive match for publically setable/getable fields
            N = L(1).subs;
            idx = find( strncmpi(N,{'Parameter';'Uncertainty';'NominalValue'},...
                length(N)) );
            
            if isempty(idx)
                error(['No property of the class "plftmat" matches '...
                    'the string "' N '". Use PROPERTIES to get the '...
                    'list of properties for this class.']);
            else
                switch idx
                    case 1
                        % Set m.Parameter
                        RB = m.RateBounds;
                        fn = fieldnames(tmp);
                        for i=1:numel(fn);
                            fni = fn{i};
                            [~,idx] = ismember(fni,RB(:,1));
                            if idx>0
                                RB{idx,1} = tmp.(fni).Name;
                                RB{idx,2} = tmp.(fni).RateBounds;
                                tmp.(fni) = tmp.(fni).Data;
                            else
                                tmp = rmfield(tmp,fni);
                            end
                        end
                        mData = usubs(m.Data,tmp);
                        [~,idx]=unique(RB(:,1));
                        
                        m = plftmat(mData,RB(idx,:));
                    case 2
                        % Set m.Uncertainty
                        RB = m.RateBounds;
                        U = setdiff(fieldnames(m.Data.Uncertainty),RB(:,1));
                        fn = fieldnames(tmp);
                        for i=1:numel(fn);
                            fni = fn{i};
                            [~,idx] = ismember(fni,U);
                            if idx==0
                                tmp = rmfield(tmp,fni);
                            end
                        end
                        mData = usubs(m.Data,tmp);
                        m = plftmat(mData,RB);
                    case 3
                        error(['You cannot set the read-only ' ...
                            'property ''NominalValue'' of a plftmat.'])
                end
            end
            
        catch
            error(lasterr);
        end
        
    case '()'
        try
            % Peform all subsasgn but L(1)
            if length(L)==1
                tmp = RHS;
            else
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
        catch
            error(lasterr);
        end
        
        % Perform subsasgn
        % Note: m and tmp can be a tvreal, umat, etc. Use BINOP code
        % to perform the subsasgn on data and repack as a plft object.
        m = binop(m,tmp,'subsasgn',L(1));
        
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end

