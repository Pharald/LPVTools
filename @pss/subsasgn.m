function m = subsasgn(m,L,RHS)
%SUBSASGN  Subscripted assignment for PSS objects.
%
% See also: subsasgn, subsref.

switch L(1).type
    case '.'
        % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
        try
            if length(L) == 1
                tmp = RHS;
            else
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
                %---------------
                N = L(1).subs;
                % prop1 = Low level properties that need to be handled.
                prop1 = {'a' 'b' 'c' 'd'};
                % prop2 = low level properties that do not need to be handled.
                prop2  = {'StateName' 'StateUnit' 'Ts' 'TimeUnit' ...
                    'InputName' 'InputUnit' 'InputGroup' ...
                    'OutputName' 'OutputUnit' 'OutputGroup' ...
                    'Name' 'Notes' 'UserData'};
                % prop3 = low level properties that will error out
                prop3  = {'e' 'Scaled' 'InternalDelay' 'InputDelay' ...
                    'OutputDelay'};
                % prop4 = High-level properties of the pss.
                prop4 = {'Domain' 'DomainPrivate' 'Data' 'DataPrivate'};
                
                % Determine matches:
                [flg,listidx] = propmatch(N,prop1,prop2,prop3,prop4);
                % flg = 0 (No match), 1 (ambiguous match), 2 (partial match),
                %       3 (exact match)
                % listidx = If there is a match (flg=2 or flg=3) then listidx(1) and
                %    listidx(2) are the list and list entry where the match occurs.
                %    listidx is returned as empty if there is no match (flg=0) or
                %    ambiguous match (flg=1)
                
                
                if flg == 0
                    error(['No property of the class "pss" matches '...
                        'the string "' N '". Use PROPERTIES to get the '...
                        'list of properties for this class.']);
                elseif flg ==1
                    error(['The specified property name "' N '" is '...
                        'ambiguous. Specify more characters for the '...
                        'property name.']);
                else
                    if listidx(1) == 1
                        % Execute for low level properties that need to be handled
                        [a1,b1,c1,d1] = ssdata(m);
                        switch listidx(2)
                            case 1
                                m = ss(tmp,b1,c1,d1,m);
                            case 2
                                m = ss(a1,tmp,c1,d1,m);
                            case 3
                                m = ss(a1,b1,tmp,d1,m);
                            case 4
                                m = ss(a1,b1,c1,tmp,m);
                        end
                    elseif listidx(1) ==2
                        % Execute for low level properties that don't need handling
                        m.DataPrivate.(L(1).subs) = tmp;
                    elseif listidx(1) ==3
                        error(['The property "' N '" is not settable for "pss" objects.']);
                    else
                        % Execute for high level properties
                        switch listidx(2)
                            case 1 % Handling 'Domain'
                                m = pss(m.Data,tmp);
                            case 2 % Handling 'DomainPrivate'
                                m.DomainPrivate = tmp;
                            case 3 % Handling 'Data'
                                m = pss(tmp,m.Domain);
                            case 4 % Handling 'DataPrivate'
                                m.DataPrivate = tmp;
                        end
                    end
                end
            end
        catch
            error(lasterr);
        end
    case '()'
        % Define m and RHS on same domain
        r = pss(RHS);
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
        m = pss(Data,m.Domain);
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end
