function [flg,listidx] = propmatch(prop,varargin)
% prop is a character string
% varargin{i} is a cell array of strings. It is assumed that
%    the same string does not appear in more than one list of strings.
%
% flg = 0 (No match), 1 (ambiguous match), 2 (partial match), 
%       3 (exact match)
% listidx = If there is a match (flg=2 or flg=3) then listidx(1) and
%    listidx(2) are the list and list entry where the match occurs.
%    listidx is returned as empty if there is no match (flg=0) or
%    ambiguous match (flg=1)
%

Nlist = length(varargin);
lprop = length(prop);
flg = 0;
listidx = [];
for i=1:Nlist
    % Case-insensitive partial match
    list = varargin{i};
    idx = find( strncmpi(prop,list,lprop) );
    lidx = length(idx);
    if lidx==1 
        if length(prop)==length(list{idx})
            % Exact match
            flg = 3;
            listidx = [i idx];
            return
        elseif flg==0
            % This is the first partial match
            flg = 2;
            listidx = [i idx];
        elseif flg==2
            % Ambiguity with a previous partial match
            flg = 1;
            return
        end
    elseif lidx>=2
        for j=1:length(idx)
           if length(prop)==length(list{idx(j)});
               % Exact Match
               flg = 2;
               listidx = [i idx(j)];
               return       
           end
        end   
        
        % Ambiguous match
        flg=1;
        return
    end
    
end

