function out = connect(varargin)
%CONNECT  Block-diagram interconnections of PSS
%
% CONNECT overloads the standard MATLAB function to work with PSS objects.
% See DynamicSystem/connect for details.
%
% See also: connect, append, series, parallel, sysic, feedback, lft.

if isa(varargin{end},'char') || isa(varargin{end},'cell')
        
    % Cycle through inputs
    varcell = cell(0,2);
    for i = 1:numel(varargin)
        vi = varargin{i};
        
        if ~isa(vi,'char') && ~isa(vi,'cell')
            varcell = [varcell;{vi,i}];
        end
        
    end
    
    % Populate all items with same domain
    [varcell{:,1}] = domunion(varcell{:,1});
    Dom = varcell{1,1}.Domain;
    
    % Get underlying data
    for i = 1:numel(varcell(:,1))
        varargin{varcell{i,2}} = varcell{i,1}.Data;
    end
else 
    % when syntax connect(blksys,connections,in,out); is used
    % in,out are not characters

    % Get domain:
    Dom = varargin{1}.Domain;
    
    % Get underlying data:
    varargin{1} = varargin{1}.Data;
end
% Apply CONNECT
out = connect(varargin{:});

% Package output into a pss that uses the common domain
out = pss(out,Dom);

