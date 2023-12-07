% UPFRD   Create an uncertain parameter-varying frequency response data model
%
%   S = UPFRD(Data,Domain) creates an uncertain parameter-varying
%   frequency response data model defined on an N-dimensional rectangular
%   grid. Domain is an RGRID object that  specifies the N independent
%   variables and the rectangular grid domain. Data is an N-dimensional
%   uncertain frequency response data (UFRD) array.
%   Data(:,:,i1,...,iN) is the uncertain frequency response data evaluated
%   at the point Domain(i1,....,iN).
% 
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 1-by-1 UFRD defined on a 1-dimensional grid
%   IVData = linspace(2,20,10);
%   Domain = rgrid('a',IVData);
%   omeg = logspace(-1,2,30);
%   unc = ureal('unc',10) 
%   usys = rss(1,1,2,10)*unc;
%   Data = ufrd(usys,omeg);
%   S = upfrd(Data,Domain)
% 
%   % Overlay Bode plots at each independent variable
%   bode(S);
%
% See also: ufrd, rgrid, pgrid, pmat, pss, pfrd, upmat, upss, pstruct.

% TODO PJS 5/1/2011: This implementation requires frequency to be constant
% across the domain. Allow it to be varying?

classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?pss,?pfrd,?uss,?ufrd,...
        ?umat,?ureal,?ucomplex,?ucomplexm,?ultidyn,?udyn}) upfrd
    % Class definition for uncertain parameter-varying FRD
    
    properties (Dependent)
        % Public data with array dimensions trailing after indep vars
        Data;
        Domain;
    end
    
    properties (Hidden = true)
        DomainPrivate = rgrid;
        DataPrivate = [];
    end
    
    properties (Dependent)
        InputName;
        OutputName;
        Frequency;
        NominalValue;
        Uncertainty;
        Parameter;
    end
    
    methods
        
        % Constructor
        function obj = upfrd(a1,a2,varargin)
            if nargin==1
                if isa(a1,'upfrd');
                    obj = a1;
                elseif ( isa(a1,'frd') || isa(a1,'ufrd') )
                    obj = upfrd(ufrd(a1),rgrid);
                elseif isa(a1,'pfrd')
                    obj = upfrd(ufrd(a1.DataPrivate),a1.DomainPrivate);
                else
                    error('Unsupported conversion to UPFRD.');
                end
            elseif nargin==2 && isa(a2,'rgrid')
                % Set properties of a UPFRD
                niv = a2.NumIV;
                szMD = size(a1);
                szad = szMD(2+niv+1:end);
                obj.DomainPrivate = rgrid(a2,szad);
                obj.DataPrivate =  a1;
                
                % Use isvalid to perform all error checking
                [pflag,errstr] = isvalid(obj);
                if pflag==0
                    error(errstr);
                end
            elseif nargin~=0
                % upfrd(Reponse,Freq,Ts) where Reponse is a PMAT
                Response = a1;
                Freq = a2;
                if nargin==3
                    Ts = varargin{3};
                else
                    Ts = 0;
                end
                obj.DataPrivate = ufrd(Response.DataPrivate,Freq,Ts);
                obj.DomainPrivate = Response.DomainPrivate;
            end
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if UPFRD object is valid.
            errstr = [];
            pflag = 1;
            if ~isa(obj.DomainPrivate,'rgrid')
                pflag = 0;
                errstr = ['Domain must be an rgrid object.'];
            elseif ~isa(obj.DataPrivate,'ufrd')
                pflag = 0;
                errstr = ['Data must be a UFRD array.'];
            else
                szDomain = size(obj.DomainPrivate);
                szData = size(obj.DataPrivate);
                lszD = length(szDomain);
                lszM = length(szData);
                szDomain = [szDomain ones(1,lszM-2-lszD)];
                szData = [szData ones(1,lszD-lszM+2)];
                if ~( all( szDomain == szData(3:end) ) || ...
                        ( isempty(obj.DomainPrivate) && ndims(obj.DataPrivate)==2)   )
                    pflag = 0;
                    errstr = ['Dimensions of Domain and Data are incompatible.'];
                end
            end
            
            % TODO PJS 4/30/2011: Anything else to check?
            % Check for varying frequency grids?
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % InputName
        function out = get.InputName(obj)
            out = obj.DataPrivate.InputName;
        end
        function obj = set.InputName(obj,Val)
            obj.DataPrivate.InputName = Val;
        end
        
        % OutputName
        function out = get.OutputName(obj)
            out = obj.DataPrivate.OutputName;
        end
        function obj = set.OutputName(obj,Val)
            obj.DataPrivate.OutputName = Val;
        end
        
        % Frequency
        function out = get.Frequency(obj)
            out = obj.DataPrivate.Frequency;
        end
        function obj = set.Frequency(obj,Val)
            obj.DataPrivate.Frequency = Val;
        end
        
        % get
        function SA = get.Data(obj)
            [IVidx,ADidx] = upidx(obj.DomainPrivate);
            % Ordering: Parameter vars, Array dims
            tmp = [IVidx' ADidx'];
            if numel(tmp)==0 ||  numel(tmp)==1
                SA = obj.DataPrivate;
            else
                SA = permute(obj.DataPrivate,tmp);
            end
        end
        
        % get
        function Domain = get.Domain(obj)
            DomainPrivate = obj.DomainPrivate;
            IVNamePrivate = DomainPrivate.IVName;
            IVDataPrivate = DomainPrivate.IVData;
            IVRateBoundsPrivate = DomainPrivate.IVRateBounds;
            
            ArrayName = DomainPrivate.ArrayName;
            idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
            Domain = rgrid( IVNamePrivate(~idx), IVDataPrivate(~idx),...
                IVRateBoundsPrivate(~idx,:));
        end
        
        % XXX should user be able to set domain parameters by hand.
        % Set Domain
        function out = set.Domain(obj,Domain)
            InData = obj.Data;
            out = upfrd(InData,Domain);
        end
        
        % get
        function out = get.NominalValue(obj)
            Domain = obj.Domain;
            Data = obj.Data.NominalValue;
            out = pfrd(Data,Domain);
        end
        
        % get
        function out = get.Uncertainty(obj)
            out = obj.Data.Uncertainty;
        end
        
        % Display
        function s=display(obj)
            % Make Header
            niv = obj.Domain.NumIV;
            S = obj.DataPrivate;
            szm = size(S);
            nfreq = length(S.Frequency);
            Ts = S.Ts;
            
            % Make header string
            [adcs,nad] = ad2char(obj.DomainPrivate);
            if nad == 0
                s1 = 'UPFRD with ';
            else
                s1 = [adcs ' array of UPFRDs with '];
            end
            % AH - 4/8/14 - FRDs don't keep track of the number of states.
%           s1 = [s1 int2str(ns) ' States, '];
            s1 = [s1 int2str(szm(1)) ' Outputs, '];
            s1 = [s1 int2str(szm(2)) ' Inputs, '];
            if Ts==0
                s1 = [s1 'Continuous System, '];
            else
                s1 = [s1 'Discrete System, Ts = ' num2str(Ts) ', '];
            end
            s1 = [s1 int2str(nfreq) ' Frequency points.'];
            
            % Make display string
            if niv>0
                dispStr = display(obj.Domain);
                tmp = 'The UPFRD consists of the following blocks:';
                dispStr = char(tmp,dispStr(2:end,:));
            else
                dispStr = '';
            end
            
            % Make display string for uncertainty
            [~,~,BlkStruct]=lftdata(obj.DataPrivate);
            Nblk = length(BlkStruct);
            for i=1:Nblk
                % Get short text description for each block
                % XXX PJS: Is getDescription a public function?
                Namei = BlkStruct(i).Name;
                Copiesi = BlkStruct(i).Occurrences;
                Blocki = obj.Data.Blocks.(Namei);
                D = getDescription(Blocki,Copiesi);
                dispStr = char(dispStr,['  ' D]);
            end
            
            s = char(s1,dispStr);
            if nargout==0
                disp(s)
                
%                 % TODO PJS 4/4/2011: Revisit. The code below displays the
%                 % matrix data if the the PMAT has only 1 IV point.
%                 if prod( szm(3:end) ) == 1 && ~isempty(obj)
%                     %             disp([inputname(1) ' = ']);
%                     %             disp(obj.Data)
%                     eval([inputname(1) ' = obj.Data'])
%                     
%                 end
            end
        end
        
        
        
%         function display(obj)
%             % Make Header
%             niv = obj.Domain.NumIV;
%             S = obj.DataPrivate;
%             szm = size(S);
%             nfreq = length(S.Frequency);
%             Ts = S.Ts;
%             if szm(1)==1
%                 ioStr = ['1 Output, '];
%             else
%                 ioStr = [int2str(szm(1)) ' Outputs, '];
%             end
%             if szm(2)==1
%                 ioStr = [ioStr '1 Input, '];
%             else
%                 ioStr = [ioStr int2str(szm(2)) ' Inputs, '];
%             end
%             if Ts==0
%                 ioStr = [ioStr 'Continuous System, '];
%                 %elseif Ts==-1
%                 %    ioStr = [ioStr 'Discrete System, unspecified Ts, '];
%             else
%                 ioStr = [ioStr 'Discrete System, Ts = ' num2str(Ts) ', '];
%             end
%             if nfreq==1
%                 ioStr = [ioStr '1 Frequency point, '];
%             else
%                 ioStr = [ioStr int2str(nfreq) ' Frequency points, '];
%             end
%             
%             if niv==1
%                 ioStr = [ioStr '1 IV'];
%             else
%                 ioStr = [ioStr int2str(niv) ' IVs'];
%             end
%             hdr = [ioStr];
%             
%             % Make header string
%             [adcs,nad] = ad2char(obj.DomainPrivate);
%             if nad==0
%                 disp(['UPFRD: ' hdr]);
%             else
%                 disp([adcs ' array of UPFRD: ' hdr]);
%             end
%             
%             % Make display string
%             if niv>0
%                 dispStr = display(obj.Domain);
%             else
%                 dispStr = '';
%             end
%             disp(dispStr)
%             
%             % TODO PJS 4/30/2011: Revisit. The code below displays the
%             % FRD if the the PFRD has only 1 IV point.
%             if prod( szm(3:end) ) == 1 && ~isempty(obj)
%                 disp([inputname(1) ' = ']);
%                 S
%                 %disp(S);
%             end
%             % GJB 26Jan12
%             % Handle Display of Uncertainty
%             zobjuss = obj.Data;
%             if isuncertain(zobjuss)
%                 warning off
%                 zobjussS = struct(zobjuss(:,:,1));
%                 warning on
%                 BD = getSummary(zobjussS.Data_.Blocks);
%                 nblk = numel(BD);
%                 if nblk==0
%                     disp('No uncertain block')
%                 else
%                     disp('The model uncertainty consists of the following blocks:')
%                     for ct=1:nblk
%                         disp(['  ' BD{ct}]);
%                     end
%                 end
%             end
%         end
        
    end % end of methods
    methods (Access = private)
        % privatesize
        function varargout = privatesize(m,arg2)
            out = size(m.DataPrivate);
            if length(out)==3
                out = [out 1];
            end
            
            if nargin==1
                arg2 = nan;
            end
            
            nout = nargout;
            if nout==2
                % Direct call to csize with nout = 2 put sthe product of the column
                % and IV dims into the second output. Instead SIZE with two
                % output args should return only the row/col dimensions.
                tmpout = csize(out,arg2,nargin,3);
                varargout = {tmpout{1}, tmpout{2}};
            else
                varargout = csize(out,arg2,nargin,nout);
            end
        end
    end
    
end % end of classdef





