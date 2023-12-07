function [a1,a2,a3,a4] = lpvblksim(t,xt,ut,Message,...
    nU,ABMat,CDMat,IVData,dIVData,niv,LIVData,x0)
% LPVBLKSIM   Utility function for Simulink simulations with LPVTools blocks
switch Message
    case 0  % init
        if isa(dIVData,'double') && isnan(dIVData)
            % PLFTSS Case
            nX = size(ABMat{4},1);
            nY = size(CDMat{4},1);
        else
            % PSS Case
            nX = size(ABMat,1);
            nY = size(CDMat,1);
        end
        
        tmp = simsizes;
        tmp.NumContStates = nX;
        tmp.NumDiscStates = 0;
        tmp.NumOutputs = nY;
        tmp.NumInputs = nU+niv;
        tmp.NumSampleTimes = 1;
        
        % XXX TODO: Check for direct feedthrough
        %dmat = CDMat(:,nX+1:end,:);
        tmp.DirFeedthrough = true; % any(dmat(:)~=0);
        
        a1 = simsizes(tmp);
        a2 = x0;        % x(0), initial condition, nx1 column
        a3 = [];        % always,
        a4 = [0 0];     % Sample time, offset
    case 1  % xdot
        Uinput = ut(1:nU);
        pValue = ut(nU+1:end);
        
        % Project pValue into the domain 
        % TODO PJS 5/27/2011: How should we handle parameters outside
        %   the specified domain?  See comments at bottom of file.        
        R = zeros(niv,2);
        for i=1:niv
            R(i,:) = [IVData{i}(1) IVData{i}(end)];
            pValue(i) = min( max(pValue(i),R(i,1)), R(i,2) );
        end        
                
        if isa(dIVData,'double') && isnan(dIVData)
            % PLFTSS Case
            M11 = ABMat{1};
            M12 = ABMat{2};
            M21 = ABMat{3};
            M22 = ABMat{4};
            r = size(M11,1);
            
            % Normalize parameter value to [-1,1].
            % TODO PJS 1/24/2014: The normalization assumes centered nominal!
            pValue = 2*(pValue - R(:,1))./(R(:,2) - R(:,1)) -1;            
            Del = pValue( LIVData{1} );
            
            M21Del = lrscale(M21,[],Del(:));
            M11Del = lrscale(M11,[],Del(:));
            AB = M22 + M21Del*( (eye(r)-M11Del)\M12 );
        else
            % PSS Case
            AB = ndLinInterp(ABMat,niv,LIVData,IVData,dIVData,pValue);
        end
        
        a1 = AB*[xt;Uinput];
        
    case 3  % output
        Uinput = ut(1:nU);
        pValue = ut(nU+1:end);
        
        % Project pValue into the domain 
        % TODO PJS 5/27/2011: How should we handle parameters outside
        %   the specified domain?  See comments at bottom of file.        
        R = zeros(niv,2);
        for i=1:niv
            R(i,:) = [IVData{i}(1) IVData{i}(end)];
            pValue(i) = min( max(pValue(i),R(i,1)), R(i,2) );
        end        
                        
        if isa(dIVData,'double') && isnan(dIVData)
            % PLFTSS Case
            M11 = CDMat{1};
            M12 = CDMat{2};
            M21 = CDMat{3};
            M22 = CDMat{4};
            r = size(M11,1);
            
            % Normalize parameter value to [-1,1].
            % TODO PJS 1/24/2014: The normalization assumes centered nominal!
            pValue = 2*(pValue - R(:,1))./(R(:,2) - R(:,1)) -1;            
            Del = pValue( LIVData{2} );
            
            M21Del = lrscale(M21,[],Del(:));
            M11Del = lrscale(M11,[],Del(:));            
            CD = M22 + M21Del*( (eye(r)-M11Del)\M12 );
        else
            % PSS Case
            CD = ndLinInterp(CDMat,niv,LIVData,IVData,dIVData,pValue);
        end
        
        a1 = CD*[xt;Uinput];
end



%    if t==0 && all(ut==0)
%        % TODO PJS 5/23/2011: If there are feedback loops, Simulink calls
%        % this function with ut=0 at t=0. Simulink is probably trying to
%        % determine signal dimensions.  This causes an error in ndLinInterp
%        % if pValue = 0 is not in the domain of sys.  Temporarily output
%        % zero in this case.  Need to revisit this.
%         a1 = zeros(size(CDMat,1),1);
%    else
%         CD = ndLinInterp(CDMat,niv,LIVData,IVData,dIVData,pValue);
%         a1 = CD*[xt;Uinput];
%    end

%         try
%             % TODO PJS 5/27/2011: The issue appears to be more complicated.
%             % On the BACT example this gets called at t=0 with Uinput non-zero
%             % but pValue = 0 (eventhough the input pValue is specified to be
%             % a non-zero constant).  In fact, even if the input parameter signal
%             % is a constant [0.5; 125], this function gets called at nonzero
%             % times with pValue not equal [0.5;125]. In the BACT example
%             % its getting called with pValue outside the Domain range
%             % and this is causing an error.
%             CD = ndLinInterp(CDMat,niv,LIVData,IVData,dIVData,pValue);
%             a1 = CD*[xt;Uinput];
%         catch
%             pValue
%             if pValue(1)>=0.5
%                 keyboard
%             end
%             a1 = zeros(size(CDMat,1),1);
%         end

