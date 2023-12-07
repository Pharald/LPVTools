function [B,SAMPLES] = usample(A,varargin)
% USAMPLE   Generates random samples of PLFT subclassess.
% 
%   B = USAMPLE(A,N) generates N random samples of the uncertainty in A 
%   and evaluates the PLFT A at each one to generate a PLFT B. The output 
%   B is a [size(A) N] array of PLFTs.
% 
%   B = USAMPLE(A) takes a single sample of the uncertainty and is
%   identical to B = USAMPLE(A,1).
%
%   [B,SAMPLES] = USAMPLE(A,N) returns the N samples of the uncerainty in A.
%   SAMPLES is a N-by-1 struct array, whose fields correspond to the
%   uncertainties in A. 
%
%   [B,SAMPLES] = USAMPLE(A,NAMES,N) only takes samples of the
%   uncertainties specified in NAMES. NAMES can be a string for a single 
%   uncertainty, or a cell array of strings for multiple uncertainties. If 
%   NAMES does not contain the names of all the uncertainties in A, the 
%   output B is a UPMAT containining the remaining uncertainties.
%
%   [B,SAMPLES] = USAMPLE(A,NAME1,N1,NAME2,N2,...) takes N1 samples of the
%   uncertainty NAME1, N2 samples of uncertainty NAME2, etc. The output B 
%   has size [size(A),N1,N2,...].
%   
%   [B,SAMPLES] = USAMPLE(A,...,WMAX) restricts the frequency range of the 
%   uncertain dynamics contained in a sample of A, when A contains ULTIDYN 
%   dynamics. The poles in the samples of the uncertain dynamics will have 
%   magnitude less than WMAX.  
%
%   See also usample, usubs.


% Grab all rate bound information
if isa(A,'plft')    
    ARateBounds = A.RateBounds;
    AData = A.Data;
else
    % XXX first input could be a uvars object, need to handle this.
    % XXX more descriptive error.
    error('Input type not allowed')
end

if nargout == 1
    BData = usample(AData,varargin{:});    
else
    [BData,SAMPLES] = usample(AData,varargin{:});
end


% Remove RB info if corresponding tvreal does not appear in the result
% XXX PJS: CData.Uncertainty will error out if C is not a umat/uss/ufrd.
if isuncertain(BData)
    if isa(BData,'umat') || isa(BData,'uss') || ...
            isa(BData,'ufrd')
        U = fieldnames(BData.Uncertainty);
    else
        U = BData.Name;
    end        
else
    U = [];
end
[~,idx]=setdiff(ARateBounds(:,1),U);
ARateBounds(idx,:)=[];

% Package data in the proper object
switch class(BData)
    case {'double','ss','frd'}
        B = BData;
    case 'umat'
        B = plftmat(BData,ARateBounds);
    case 'uss'
        B = plftss(BData,ARateBounds);
    case 'ufrd'
        B = plftfrd(BData,ARateBounds);
    otherwise
        % XXX PJS: Update error to provide more info.
        error('Unsupported operation.');
end

