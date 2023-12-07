function [omega,defaultpl] = polefreqauto(G,unames,polelist,omega)
% This function is a utility that automatically generates a
% frequency grid and a pole list for IQC functions.

if isequal(omega,'nogenfresp') || isempty(omega)
    % Autoselect frequency grid.
    zp = zpk(G);
    zp = [zp.z{:};zp.p{:}];
    zp = zp(abs(zp)>1e5*eps & imag(zp)>=0,:);
    fzp = abs(zp);
    
    if isempty(zp)
        fmax = 100;
        fmin = 0.01;
    else
        fmin = min(fzp)/20;
        fmax = max(fzp)*20;
    end
    if isempty(omega)
        %[~,~,omega] = genfresp(G,4,{fmin,fmax});
        [~,omega] = sigma(G,{fmin,fmax});
    end
elseif iscell(omega)
    if omega{1} == 0
        fmin = omega{2};
    else
        fmin = omega{1};
    end
    fmax = omega{2};
    %[~,~,omega] = genfresp(G,4,{fmin,fmax});
    [~,omega] = sigma(G,{fmin,fmax});
elseif isa(omega,'double')
    if omega(1) == 0
        fmin = omega(2);
    else
        fmin = omega(1);
    end
    fmax = omega(end);
end


% Set default values:
if isa(polelist,'double') && ~isempty(polelist)
    fg = polelist;
else
    % XXX what is a good heuristic?
        
    %     % Option #1 - min, max and geometric mean
    %     fg = [-fmin -sqrt(fmin*fmax) -fmax];
    
    % Option #2 - N equally spaced points on a log scale:
    Npts = 3;
    lfmin = log(fmin); 
    lfmax = log(fmax);
    
    lfg = linspace(lfmin,lfmax,Npts+2);
    fg = -exp(lfg(2:end-1));    
    fg = [fg inf];
end

% Create a default pole list for all the uncertainties:
defaultpl = unames;
defaultpl(:,2) = {fg};
defaultpl = defaultpl';
defaultpl = struct(defaultpl{:});


% Replace default values with user specified pole lists:
if isa(polelist,'struct')
    fn = fieldnames(polelist);    
    for r = 1:numel(fn)
        defaultpl.(fn{r}) = polelist.(fn{r});
    end
end