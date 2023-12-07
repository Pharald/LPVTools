function [y,t,x,u,ptrajout] = lpvimpulse(P,ptraj,T,opt)
% LPVIMPULSE  Parameter dependent impulse response for PLFTSS objects
% 
% [Y,T,X,U,PTRAJOUT] = LPVIMPULSE(SYS,PTRAJ) computes the parameter dependent
% impulse response of SYS. SYS is a PLFTSS with Ny outputs, Nx states,
% Nu inputs, and N independent variables IVName1,...,IVNameN. PTRAJ is a 
% struct which defines the time-variation of the parameters (independent 
% variables) in SYS. The field PTRAJ.time contains a sorted row vector of 
% time-values. PTRAJ must also have a field for each independend variable 
% in SYS, such that PTRAJ.IVName1, ... ,PTRAJ.IVNameN each contain a row 
% vector of parameter trajectories corresponding to PTRAJ.time. 
% The output Y is a length(T)-by-NY-by-Nu matrix such that Y(:,i,j)  
% corresponds to the i-th output of SYS due to a step command in the j-th 
% input channel. Similarly X is a length(T)-by-Nx-by-Nu matrix describing
% the state trajectories of SYS, U is a length(T)-by-Nu-by-Nu matrix 
% describing the trajectory of the inputs to SYS, and T is a column vector 
% of time values corresponding to Y, X and U. PTRAJOUT contains the 
% corresponding parameter trajectories.
%
% LPVIMPULSE(SYS,PTRAJ) generates plots of the parameter dependent impulse 
% response of SYS.
% 
% [Y,T,X,U,PTRAJOUT] = LPVIMPULSE(SYS,ptraj,TFINAL) simulates the impulse
% response up to the time TFINAL.
%
% [Y,T,X,U,PTRAJOUT] = LPVIMPULSE(SYS,ptraj,T) simulates the impulse response
% using a user supplied time vector T.
% 
% See also: lpvlsim, lpvinitial, lpvimpulse, lsim.

% TODO - AH - 11/14/13 - Need to add functionality that pads time vector 
% when it is too course to capture very fast dynamics. Choose time vector
% padding based on magnitude of fastest pole?

narginchk(2,4)
% Error checking
if isfield(ptraj,'time')
    time = ptraj.time;
else
    error('ptraj must contain a field "time"')
end
opt = [];
x0 = [];
if nargin ==2   
    T = time(end);
    if nargout == 0
        lpvimpulse(P,ptraj,T);
    else
        [y,t,x,u,ptrajout] = lpvimpulse(P,ptraj,T);
    end    
    return
elseif nargin>=3
    if isnumeric(T) && isscalar(T)
        tIn = [0;T];        
    elseif isnumeric(T) && isvector(T)
        tIn = T(:);
    elseif ~isnumeric(T)
        opt = T;
        tIn = ptraj.time;
    end   
end

ptrajout = [];
% Deal with array dimensions. Array dimensions only allowed if nargout = 0
szP = size(P);
if nargout == 0
    
    szP = size(P);
    ny = szP(1);
    nu = szP(2);
    L.type = '()';
    if numel(szP)>2
        % Array dimensions present. Need to loop through the array.
        nAD = prod(szP(3:end));
    else
        nAD = 1;
    end
    
    % Add a time point at t(1)+1e-6
    tIn = [tIn(1);tIn(1)+1e-6;tIn(1)+2e-6;tIn(2:end)];
    for i = 1:nAD
        % Grab model at i-th array index
        L.subs = {':',':',i};
        Pi = subsref(P,L);
        
        for k = 1:nu
            % Loop through input dimensions and do axis-by-axis step
            % commands.
            
            % Add impulse to u by setting u(t(1)) = 1 and u(t(1)+1e-6)= 0;       
            uIn = zeros(numel(tIn),size(P,2));
            uIn(2,k) = 2/1e-6;            
            [yik,tik,~,~,~] = lpvlsim(Pi,ptraj,uIn,tIn,x0,opt);
            
            for j = 1:ny
                subplot(ny,nu,(j-1)*ny+k)
                plot(tik,yik(:,j));
                if j ==1
                    title('Input');
                end
                if j == ny
                    xlabel('Time [sec]');
                end
                if k==1
                    ylabel('Output');
                end                                
            end
        end                     
        
        if nAD>1
            % Set plot holds to on/off if looping through array
            if i==1
                for j=1:ny*nu
                    % switch to "hold on" in all subplots
                    subplot(ny,nu,j);
                    hold on;
                end
            end
            if i==nAD
                for j=1:ny*nu
                    % switch to "hold off" in all subplots
                    subplot(ny,nu,j);
                    hold off;
                end
            end
        end        
        
    end
    
else
    if numel(szP)>2
        error(['The "lpvimpulse" command operates on a non-array dimensional'...
            ' model when used with output arguments.'])
    else                
        
        [ny,nu] =size(P);
        nx = size(P.Data.a,1);
        % Add a time point at t(1)+1e-6
        tIn = [tIn(1);tIn(1)+1e-6;tIn(1)+2e-6;tIn(2:end)];
        for k = 1:nu
            % Loop through input dimensions and do axis-by-axis step
            % commands.
            
            % Add impulse to u by setting u(t(1)) = 1 and u(t(1)+1e-6)= 0;           
            uIn = zeros(numel(tIn),size(P,2));
            uIn(2,k) = 2/1e-6;
            [yik,tik,xik,uik,ptrajoutik] = lpvlsim(P,ptraj,uIn,tIn,x0,opt);
            if k ==1
                y = zeros(numel(tik),ny,nu);
                x = zeros(numel(tik),nx,nu);
                u = zeros(numel(tik),nu,nu);                
            end
            y(:,:,k) = yik;
            x(:,:,k) = xik;
            u(:,:,k) = uik;
            ptrajout = [ptrajout ptrajoutik];
        end
        t = tik;
    end
end
