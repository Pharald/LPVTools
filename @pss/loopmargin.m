function varargout = loopmargin(varargin)
%LOOPMARGIN  Compute stability margins of PSS feedback loops.
%
% [CM,DM,MM] = LOOPMARGIN(L) computes the stability margins of the N-by-N
% parameter-varying system L in the negative feedback loop (use -L for
% positive feedback loops):
%
%         u --->O---->[ L ]----+---> y
%             - |              |
%               +<-------------+
%
% The stability margins are computed at each point in the domain of L.
%
%
% [CMI,DMI,MMI,CMO,DMO,MMO,MMIO] = LOOPMARGIN(P,C) computes the stability
% margins of the feedback loop formed by the N-by-Nu system P and the
% Nu-by-N controller C:
%
%         u --->O--->[ C ]-->[ P ]---+---> y
%             - |                    |
%               +<-------------------+
%
% The outputs are PSTRUCTs describing the results at each point on the
% domain of the original system:
%
%   CM = CMO = N-by-1 PSTRUCT of classical loop-at-a-time gain and phase
%              margins for each feedback channel at the output of P.
%
%   DM = DMO = N-by-1 PSTRUCT of loop-at-a-time disk margins for each
%              feedback channel at the output of P.
%
%   MM = MMO = 1-by-1 PSTRUCT of the multi-loop disk margin at the output
%              of P.
%
%   CMI      = N-by-1 PSTRUCT of classical loop-at-a-time gain and phase
%              margins for each feedback channel at the input of P.
%
%   DMI      = N-by-1 PSTRUCT of loop-at-a-time disk margins for each
%              feedback channel at the input of P.
%
%   MMI      = 1-by-1 PSTRUCT of the multi-loop disk margin at the input
%              of P.
%
%   MMIO     = 1-by-1 PSTRUCT of the multi-loop disk margin for
%              simulteneous/independent variation in all channels of the
%              feedback loop.
%
% To compute a subset of the outputs, an optional string argument can be
% passed to the function, e.g. [MM,CM] = LOOPMARGIN(L,'m,c') will compute
% only the multiloop disk margins and the classical margins, and 
% [DMI,DMO] = LOOPMARGIN(P,C,'di,do') will compute only the disk margins on
% the input and output channels. 
%
% Refer to the documentation for DynamicSystem/loopmargin for details.
%
% See Also: loopmargin, loopsens.


nin = nargin;
nout = nargout;
P = varargin{1};
if isa(varargin{1},'pss')
    if nin > 1 && ~ischar(varargin{2})
        % DOMUNION ensures that IV ordering is identical.
        [P,varargin{2}] = domunion(P,varargin{2});
    end
    Domain = P.Domain;
    varargin{1} = P.Data;
    szO = size(P.Data);
else
    P = pss(P);
end


if nin > 1 && isa(varargin{2},'pss')
    K = varargin{2};
    [P,K] = domunion(P,K);
    szO = size(K.Data);
    Domain = K.Domain;
    varargin{2} = K.Data;
end


% Array dimensions currently not allowed
nad = numel(size(P))-2;
if nad>0
    error('LOOPMARGIN does not allow PSS with array dimensions.');
end

if nin ==1 || (nin==2 && ~ischar(varargin{2}))
    % Call loopmargin with all LTI objects
    % [CMI,DMI,MMI,CMO,DMO,MMO,MMIO] = loopmargin(varargin{:});
    output = cell(1,7);
    [output{:}] = loopmargin(varargin{:});
        
    if nout >= 1
        % CMI
        sz = size(output{1});
        output{1} = pstruct( reshape(output{1},[sz(1) 1 szO(3:end)]) , Domain);
    end
    if nout >= 2
        % DMI
        sz = size(output{2});
        output{2} = pstruct( reshape(output{2},[sz(1) 1 szO(3:end)]) , Domain);
    end
    if nout >= 3
        % MMI
        output{3} = pstruct( reshape(output{3},[1 1 szO(3:end)]) , Domain);
    end
    if nout >= 4
        % CMO
        sz = size(output{4});
        output{4} = pstruct( reshape(output{4},[sz(1) 1 szO(3:end)]) , Domain);
    end
    if nout >= 5
        % DMO
        sz = size(output{5});
        output{5} = pstruct( reshape(output{5},[sz(1) 1 szO(3:end)]) , Domain);
    end
    if nout >= 6
        % MMO
        output{6} = pstruct( reshape(output{6},[1 1 szO(3:end)]) , Domain);
    end
    if nout == 7
        % MMIO
        output{7} = pstruct( reshape(output{7},[1 1 szO(3:end)]) , Domain);
    end
    
else
    % User is specifying a subset of desired outputs using optional
    % string argument
    
    % Split up the input string
    if ischar(varargin{2})
        %loopmargin(P,chars) call
        strings = strsplit(varargin{2},',');
    else
        %loopmargin(P,K,chars) call
        strings = strsplit(varargin{3},',');
    end
    
    % Compute the loopmargin results and package results into output
    output = cell(1,numel(strings));
    [output{:}] = loopmargin(varargin{:});
        
    for i = 1:numel(output)
        
        % loopmargin(P,chars) call
        % Keywords: 'ci','di','mi','co','do','mo','mm'
        
        % loopmargin(P,K,chars) call
        % Keywords: 'c','d','m'
        
        if strcmpi(strings{i}(1),'m')
            output{i} = pstruct( reshape(output{i},[1 1 szO(3:end)]) , Domain);
        else
            sz = size(output{i});
            output{i} = pstruct( reshape(output{i},[sz(1) 1 szO(3:end)]) , Domain);
        end
    end
    
end

varargout = output;