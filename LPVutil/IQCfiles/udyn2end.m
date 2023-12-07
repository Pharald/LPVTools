function [M,allPi,mublk,udyncnt,udynsize,udyncop] = udyn2end(G)
%
% MUSYNID
[M,b,Ulist] = lftdata(G);
% Code to reorder UDYNs to the end (which is what IQCSOLVE NEEDS)
allPi = cell(0,1);
mublk = zeros(0,2);  % build this for the non-IQC blocks, which are moved to top
udynsize = zeros(0,2);
udyncop = zeros(0,1);
udynridx = zeros(0,1); % this index resorts the rows of DELTA (cols of M)
udyncidx = zeros(0,1); % this index resorts the cols of DELTA (rows of M)
muridx = zeros(0,1); % this index resorts the rows of DELTA (cols of M)
mucidx = zeros(0,1); % this index resorts the cols of DELTA (rows of M)
rpt = 1;
cpt = 1;
udyncnt = 0;  % count the number of UDYN blocks
for i=1:numel(Ulist)
    sz = Ulist(i).Size;
    cop = Ulist(i).Occurrences;
    switch Ulist(i).Type
        case 'udyn'
            %if cop==1
                udyncnt = udyncnt + 1;
                udynridx = [udynridx rpt:rpt+sz(1)*cop-1];
                udyncidx = [udyncidx cpt:cpt+sz(2)*cop-1];
                rpt = rpt+sz(1)*cop;
                cpt = cpt+sz(2)*cop;
                %  UserData will have two fields: i) handle to IQC setup
                %  function and ii) parameters.
                allPi{udyncnt,1} = G.Uncertainty.(Ulist(i).Name).UserData;
                udynsize = [udynsize; sz];
                udyncop = [udyncop; cop];
%             else
%                 % i don't think we know how to deal with general, repeated
%                 % nonlinearities (via IQCs or otherwise)
%                 error('Currently cannot address repeated UDYN blocks.');
            %end
        case 'ureal'
            mublk = [mublk;-cop 0];
            muridx = [muridx rpt:rpt+cop-1];
            mucidx = [mucidx cpt:cpt+cop-1];
            rpt = rpt + cop;
            cpt = cpt + cop;
        case 'ucomplex'
            mublk = [mublk;cop 0];
            muridx = [muridx rpt:rpt+cop-1];
            mucidx = [mucidx cpt:cpt+cop-1];
            rpt = rpt + cop;
            cpt = cpt + cop;
        case {'ultidyn' 'ucomplexm'}
            if isequal(sz,[1 1])
                mublk = [mublk;cop 0];
                muridx = [muridx rpt:rpt+cop-1];
                mucidx = [mucidx cpt:cpt+cop-1];
                rpt = rpt + cop;
                cpt = cpt + cop;
            elseif cop==1
                mublk = [mublk;sz];
                muridx = [muridx rpt:rpt+sz(1)-1];
                mucidx = [mucidx cpt:cpt+sz(2)-1];
                rpt = rpt + sz(1);
                cpt = cpt + sz(2);
            else
                error('Currently cannot address repeated, non-scalar ULTIDYN or UCOMPLEXM blocks.');
            end
    end
end
szG = size(G);
yidx = cpt+(0:szG(1)-1);
uidx = rpt+(0:szG(2)-1);
M = M([mucidx udyncidx yidx],[muridx udynridx uidx]);
% At this point, M, allPi, mublk are almost ready to go into IQCsolve.
% Depending on STAB/PERF, we would need another IQC
% for the performance IQC (which doesn't have to be gain).