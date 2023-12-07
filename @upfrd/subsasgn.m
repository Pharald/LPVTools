function m = subsasgn(m,L,RHS)
% SUBSASGN  Subscripted assignment for UPFRD objects.
%
% See also: subsasgn, subsref.

switch L(1).type
    case '.'
        % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
        try
            if length(L) == 1
                tmp = RHS;
            else
                %tmp = m.(L(1).subs);
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
            if strcmp(L(1).subs,'Parameter')
                % Parameter Gateway SUBSASGN
                fn = fieldnames(tmp);
                D = m.DomainPrivate;
                for i=1:numel(fn);
                    fni = fn{i};
                    if ~strncmp( fni, D.ArrayName, length(D.ArrayName) )
                        [~,idx] = ismember(fni,D.IVName);
                        if idx>0
                            GDi = tmp.(fni).GridData;
                            D.IVName{idx} = tmp.(fni).Name;
                            D.IVData{idx} = GDi;
                            D.IVRateBounds(idx,:) = tmp.(fni).RateBounds;
                            D.LIVData(idx) = numel(GDi);
                            D.DIVData{idx} = diff(GDi);
                        end
                    end
                end
                m.DomainPrivate = D;
            else
                m.(L(1).subs) = tmp;
            end
        catch
            error(lasterr);
        end
    case '()'
        % Define m and RHS on same domain
        r = upfrd(RHS);
        StinkyPeteFlag = true;
        [m,r] = domunion(m,r,StinkyPeteFlag);
        RData = r.Data;
        Data = m.Data;  % Use Data so that array dims are at end
        
        L1 = L(1);
        niv = m.Domain.NumIV;
        if numel(L1.subs)>2
            L1.subs = [L1.subs([1 2]) repmat({':'},1,niv) L1.subs(3:end)];
        end
        Data = subsasgn(Data,L1,RData);
        m = upfrd(Data,m.Domain);
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end
