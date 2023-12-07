function LTF = loopsens(P,C)
% LOOPSENS  Feedback loop sensitivity functions for PSS objects
%  
% LOOPTRANSFER = LOOPSENS(P,C) constructs the loop transfer functions for
% the multivariable feedback loop consisting of the PSS C in negative 
% feedback with the PSS P, as shown in the diagram below. The loop transfer 
% functions are computed at each point in the combined domains of P and C.  
%              
%     ----->0------>[ P ]----------+----> 
%           |-                     |
%           |                      |
%     <-----+-------[ C ]<---------0<---- 
% 
% LOOPTRANSFER is a PSTRUCT containing fields for all input/output 
% transfer functions of the feedback loop at each point in the domain of P 
% and C.  It also has fields Poles and Stable that contain the closed-loop 
% poles and a flag for the stability of the feedback loop. Additional 
% details on the fields of LOOPTRANSFER can be found in the documentation 
% for DynamicSystem/LOOPSENS.
%
% See also: loopsens.

% Define P and C on a common domain
% (domunion will lift to pss if required)
[Pext,Cext] = domunion(P,C);
szP = size(Pext);
niv = Pext.Domain.NumIV;
nad = numel(szP)-2;
if nad>0
    % TODO AH 9/25/14 - Implementation for array dimensions
    % LOOPSENS has the quirk that array dimensions are contained in the fields
    % of the resuls, instead of the structure itself, i.e. the results for
    % a 2x3xAD1xAD2 SS will be returned as a 1x1 struct with fields that 
    % are 2x3xAD1xAD2 in size.This behaviour is different from the behaviour 
    % seen in e.g. robustperf where the result is a AD1xAD2 struct.
    error('LOOPSENS does not allow PSS with array dimensions.');
end
LIVData = Pext.Domain.LIVData;
LIVData = LIVData';

% Reorder as [row col AD IV]
IVidx = 1:niv;
ADidx = (1:nad)+niv;
if nad>0
    Pdata = permute(Pext.Data,[ADidx IVidx]);
    Cdata = permute(Cext.Data,[ADidx IVidx]);
else
    Pdata = Pext.Data;
    Cdata = Cext.Data;
end

% LTI/LOOPSENS at each point in IV domain
npts = prod(LIVData);
idx = [cellstr(repmat(':',2+nad,1))' 1];
for i=1:npts
    idx{end} = i;
    if i==1
        LTF = loopsens( Pdata(idx{:}), Cdata(idx{:}) );
    else
        LTF(idx{:}) = loopsens( Pdata(idx{:}), Cdata(idx{:}) );
    end
end

% Reshape 1-by-npts struct array LFT into dimensions [row col IV]
LTF = reshape(LTF,[1 1 LIVData]);
LTF = pstruct(LTF, Pext.Domain);


%% Old Code
% 
% % Define P and C on a common domain
% % (domunion will lift to pss if required)
% [Pext,Cext] = domunion(P,C);
% szP = size(Pext.Data);
% if numel(szP) == 2
%     szP = [szP 1 1];
% end
% 
% % Reorder as [row col AD IV]
% niv = Pext.Domain.NumIV;
% nad = numel(szP)-niv-2;
% IVidx = 3:2+niv;
% ADidx = 3+niv:2+niv+nad;
% Pdata = permute(Pext.Data,[ADidx IVidx]-2);
% Cdata = permute(Cext.Data,[ADidx IVidx]-2);
% szPdata = size(Pdata);
% 
% % LTI/LOOPSENS at each point in IV domain
% npts = prod(szP(IVidx));
% idx = [cellstr(repmat(':',2+nad,1))' 1];
% % loopsens acts on LTI arrays to put them into a struct with array fields.
% % Hence LTF will be a 1xnpts struct array. Array dimensions are moved into
% % the fields.
% for i=1:npts
%     idx{end} = i;
%     LTF(i) = loopsens( Pdata(idx{:}), Cdata(idx{:}) );
% end
% 
% % Reshape 1-by-npts struct array LFT into dimensions [row col IV]
% LTF = reshape(LTF,[1 1 szP(IVidx)]);
% % Reorder dimensions as [row col IV AD]
% LTF = permute(LTF,[1 2 IVidx ADidx]);
% 
% LTF = pstruct(LTF, Pext.Domain);

