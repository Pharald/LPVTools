function m = subsref(m,L)
% SUBSREF  Subscripted reference for PFRD objects.
%
% See also: subsref, subsasgn.

switch L(1).type
    case '.'
        try
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
                m = m.(L(1).subs);
            end
        catch
            error(lasterr);
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
        m = pfrd(Data,m.Domain);
        
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
