function m = subsref(m,L)
% SUBSREF  Subscripted reference for PSS objects.
%
% See also: subsref, subsasgn.

switch L(1).type
    case '.'
        if strncmp(L(1).subs,'value',3)
            % Default call is m.val(IVN,IVD) where IVN and IVD are cell
            % arrays of names/values.  Convert other calls to this type.
            if isnumeric(L(2).subs{1})
                % m.val(4,5)
                IVN = m.Domain.IVName;
                IVval = L(2).subs;
                L(2).subs = {IVN(:),IVval(:)};
            elseif ischar(L(2).subs{1})
                n = numel(L(2).subs);
                if floor(n/2)==ceil(n/2)
                % m.val('a',4,'b',5)
                    tmp = reshape(L(2).subs,[2 n/2])';
                    L(2).subs = { tmp(:,1), tmp(:,2) };
                else
                % m.val('a',4,'b',5,'linear')
                    tmp = reshape(L(2).subs(1:end-1),[2 (n-1)/2])';
                    L(2).subs = { tmp(:,1), tmp(:,2), L(2).subs{end} };                    
                end                
            end
            
            IVN = m.Domain.IVName;
            IVD = m.Domain.IVData;
            un = L(2).subs{1};
            ud = cell2mat(L(2).subs{2});
            nv = numel(IVN);                                   
            
            allidx = repmat(-1,[nv 1]);
            for i=1:nv
                [tf,idx] = ismember( IVN(i),un  );
                if tf
                    tmp = find(IVD{i} == ud(idx));
                    if ~isempty(tmp)
                        allidx(i) = tmp;
                    end
                end
            end
            
            if all( allidx>0 )
                allidx = num2cell(allidx);                                 
                nad = m.DomainPrivate.NumIV-nv;
                id = repmat({':'},1,nad);                
                % Evaluate desired point
                m2 = m.Data(:,:,allidx{:},id{:});
                % Eliminate singleton IV dims
                m2 = permute(m2,[1+nv:nv+nad 1:nv]);
            else
                m2 = lpvinterp(m,L(2).subs{:});
                m2 = depromote(lpvelimiv(m2)); 
                if isequal(class(m),class(m2))
                    error('The ''value'' method must evaluate all parameters')
                end
            end                                                                                 
            m = m2;
            L(2) = [];
            
        elseif strncmp(L(1).subs,'index',3)
            % Default call is m.index(IVN,IVD) where IVN and IVD are cell
            % arrays of names/indices.  Convert other calls to this type.
            if isnumeric(L(2).subs{1})
                % m.val(4,5)
                IVN = m.Domain.IVName;
                IVidx= L(2).subs;
                L(2).subs = {IVN(:),IVidx(:)};
            elseif ischar(L(2).subs{1})
                n = numel(L(2).subs);
                if floor(n/2)==ceil(n/2)
                % m.val('a',4,'b',5)
                    tmp = reshape(L(2).subs,[2 n/2])';
                    L(2).subs = { tmp(:,1), tmp(:,2) };
                else
                    error('The ''index'' method must take in IV indices.')                    
                end                
            end
            m2 = lpvsplit(m,L(2).subs{:},'index');
            m2 = depromote(lpvelimiv(m2));
            if isequal(class(m),class(m2))
                error('The ''index'' method must evaluate all parameters')
            end                                                                                          
            m = m2;
            L(2) = [];
                
                
        elseif strcmp(L(1).subs,'Parameter')
            % Parameter Gateway for SUBSREF
            D = m.DomainPrivate;
            niv = D.NumIV;
            s = struct;
            for i=1:niv
                n = D.IVName{i};
                if ~strncmp( n, D.ArrayName, length(D.ArrayName) )
                    s.(n) = pgrid( n, D.IVData{i}, D.IVRateBounds(i,:) );
                end
            end
            m=s;
        else            
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
                if listidx(1) ==1
                    % Execute for low level properties that need to be handled
                    m = pmat(m.Data.(L(1).subs),m.Domain);
                elseif listidx(1) ==2
                    % Execute for low level properties that don't need handling
                    m = m.Data.(L(1).subs);
                elseif listidx(1) ==3
                    error(['The property "' N '" is not accessable for "pss" objects.']);
                else
                    % Execute for high level properties
                    switch listidx(2)
                        case 1 % Handling 'Domain'
                            m = m.Domain;
                            %                         DomainPrivate = m.DomainPrivate;
                            %                         IVNamePrivate = DomainPrivate.IVName;
                            %                         IVDataPrivate = DomainPrivate.IVData;
                            %                         IVRateBoundsPrivate = DomainPrivate.IVRateBounds;
                            %
                            %                         ArrayName = DomainPrivate.ArrayName;
                            %                         idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
                            %                         m = rgrid( IVNamePrivate(~idx), IVDataPrivate(~idx),...
                            %                             IVRateBoundsPrivate(~idx,:));
                        case 2 % Handling 'DomainPrivate'
                            m = m.DomainPrivate;
                        case 3 % Handling 'Data'
                            m = m.Data;
                            %                         DomainPrivate = m.DomainPrivate;
                            %                         if ~isempty(DomainPrivate)
                            %                             ArrayName = DomainPrivate.ArrayName;
                            %                             idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
                            %                             % Ordering: Parameter vars, Array dims
                            %                             m = permute(m.DataPrivate,[find(~idx); find(idx)]);
                            %                         else
                            %                             m = m.DataPrivate;
                            %                         end
                        case 4 % Handling 'DataPrivate'
                            m = m.DataPrivate;
                    end
                end
            end
        end
    case '()'
        L1 = L(1);
        if length(L1.subs)==1
            error(['Use two or more subscripts to index ' ...
                'into MIMO models or model arrays, as in the "sys(2,1)" command.']);
        end
        niv = m.Domain.NumIV;
        Data = m.Data;      % Use Data so that array dims are at end
        if ~(numel(L1.subs)==2)
            L1.subs = [L1.subs([1 2]) repmat({':'},1,niv) L1.subs(3:end)];
        end
        Data = subsref(Data,L1);
        m = pss(Data,m.Domain);
        
        if length(L1.subs)==1 && isnumeric(L1.subs) && ...
                size(L1.subs{1},1) == 1 && size(m,1)>1
            m = m';
        end
    otherwise
        error(['{}-like reference is not supported for ' class(m) 'objects.'])
end

if length(L)>1
    m = subsref(m,L(2:end));
end

