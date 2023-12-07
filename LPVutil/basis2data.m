function [Gdata,RBx,BFdatax,Pdatax,RBy,BFdatay,Pdatay] = basis2data(G,Xb,Yb)

% 1/14/2014 TODO PJS: Need to review this function to ensure that
% we are properly handling array dimnesions (public vs. private data)

% Get plant domain and size
Dom = G.Domain;
nmod = prod(Dom.LIVData);
nparg = G.Domain.NumIV;
IVNg = G.Domain.IVName;
Gdata = G.Data(:,:,:);

% Get data for Xb basis
BFx = Xb.BasisFunction;
Px  = Xb.Partials;
IVNx = Xb.IVName;
nbasisx = size(BFx,1);
nparx = numel(IVNx);
if ~all( ismember(IVNx,IVNg) )
    % Make sure bases functions depend on a subset of the IVs of G
    error('Xb has IV dependence that G does not have.');
end

% Construct Xb data:
% BFdatax: Basis functions used in storage function Nbasis-by-1-by-nmod
% Pdatax: Gradient of the basis functions wrt. the parameters:
%         Nbasis-by-Nparameters-by-nmod
% RBx: Upper and lower bounds on the parameter rates: Nparameter-by-2
BFdatax = zeros(nbasisx,1,nmod);
Pdatax = zeros(nbasisx,nparx,nmod);
RBx = zeros(nparx,2);

% Evaluate basis function and partials
for imod = 1:nmod
    % Get parameter value
    if nparg ==0
        pvaluec = {[]};
    else
        pvaluec = num2cell(Dom(imod));  % single index into RGRID gives value
    end
    
    % XXX 8/17/15 PJS Convert for-loop to single LPVSPLIT Call
    BFdatax(:,1,imod)=double(lpvsplit(BFx,IVNg,pvaluec,'value'));
    Pdatax(:,:,imod)=double(lpvsplit(Px,IVNg,pvaluec,'value'));        
%     for ibase = 1:nbasisx
%         % Call to lpvsplit is using a range call with a singleton pvaluec
%         BFdatax(ibase,1,imod)=double(lpvsplit(BFx(ibase),IVNg,pvaluec,'value'));
%         Pdatax(ibase,:,imod)=double(lpvsplit(Px(ibase,:),IVNg,pvaluec,'value'));
%     end
end

% Populate Rate bound data structure for Xb
for ipar = 1:nparx
    tf2 = strcmp( IVNx{ipar} , IVNg);
    RBx(ipar,:) = G.Domain.IVRateBounds(tf2,:);
end


% Handle Yb
BFdatay = [];
Pdatay= [];
RBy = [];
if nargin ==3
    
    % Get data for Yb basis
    BFy = Yb.BasisFunction;
    Py  = Yb.Partials;
    IVNy = Yb.IVName;
    nbasisy = size(BFy,1);
    npary = numel(IVNy);
    if ~all( ismember(IVNy,IVNg) )
        % Make sure basis functions depend on a subset of the IVs of G
        error('Yb has IV dependence that G does not have.');
    end
    
    % Construct Yb data:
    % BFdatay: Basis functions used in storage function Nbasis-by-1-by-nmod
    % Pdatay: Gradient of the basis functions wrt. the parameters:
    %         Nbasis-by-Nparameters-by-nmod
    % RBy: Upper and lower bounds on the parameter rates: Nparameter-by-2
    BFdatay = zeros(nbasisy,1,nmod);
    Pdatay = zeros(nbasisy,npary,nmod);
    RBy = zeros(npary,2);
    
    % Evaluate basis function and partials
    for imod = 1:nmod
        % Get parameter value
        if nparg ==0
            pvaluec = cell(0,1);
        else
            pvaluec = num2cell(Dom(imod));  % single index into RGRID gives value
        end
        
        % XXX 8/17/15 PJS Convert for-loop to single LPVSPLIT Call
        BFdatay(:,1,imod)=double(lpvsplit(BFy,IVNg,pvaluec,'value'));
        Pdatay(:,:,imod)=double(lpvsplit(Py,IVNg,pvaluec,'value'));
%         for ibase = 1:nbasisy
%             % Call to lpvsplit is using a range call with a singleton pvaluec
%             BFdatay(ibase,1,imod)=double(lpvsplit(BFy(ibase),IVNg,pvaluec,'value'));
%             Pdatay(ibase,:,imod)=double(lpvsplit(Py(ibase,:),IVNg,pvaluec,'value'));
%         end
    end
    
    % Populate Rate bound data structure for Yb
    for ipar = 1:npary
        tf2 = strcmp( IVNy{ipar} , IVNg);
        RBy(ipar,:) = G.Domain.IVRateBounds(tf2,:);
    end
end