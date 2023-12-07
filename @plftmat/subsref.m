function m = subsref(m,L)
% SUBSREF  Subscripted reference for PLFTMAT objects.
%
% See also: subsref, subsasgn.

switch L(1).type
    case '.'
        % Case-insentive match for publically setable/getable fields
        N = L(1).subs;
        idx = find( strncmpi(N,{'Parameter';'Uncertainty';...
            'NominalValue';'RateBounds';'Data'},length(N)) );
        
        if isempty(idx)
            error(['No property of the class "plftmat" matches '...
                'the string "' N '". Use PROPERTIES to get the '...
                'list of properties for this class.']);
        else
            switch idx             
                case 1
                    % Get m.Parameter
                    U = m.Data.Uncertainty;
                    RateBounds = m.RateBounds;
                    m = struct;
                    for i=1:size(RateBounds,1)
                        n = RateBounds{i,1};
                        RB = RateBounds{i,2};
                        m.(n) = tvreal( n, U.(n).Range, RB);
                    end
                case 2
                    % Get m.Uncertainty
                    U = m.Data.Uncertainty;
                    RateBounds = m.RateBounds;
                    m = rmfield(U,RateBounds(:,1));
                case 3
                    % Get m.NominalValue (Only replace true uncerts by 
                    % nominal. Leave tvreals unchanged)
                    U = m.Data.Uncertainty;
                    RateBounds = m.RateBounds;
                    U = rmfield(U,RateBounds(:,1));
                    fn = fieldnames(U);
                    for i=1:numel(fn)                        
                        U.(fn{i}) = U.(fn{i}).NominalValue;
                    end
                    m.Data = usubs(m.Data,U);
                    m.Data = umat(m.Data);
                case 4
                    m = m.RateBounds;
                case 5
                    m = m.Data;
            end
        end
                
    case '()' 
        % Perform subsref on Data
        mData = subsref(m.Data,L(1));
        
        % Remove RB info if mData does not depend on the tvreal
        mRateBounds = m.RateBounds;
        U = fieldnames(mData.Uncertainty);
        [~,idx]=setdiff(mRateBounds(:,1),U);
        mRateBounds(idx,:)=[];
        
        % Repackage final result as an plft object
        switch class(mData)
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
