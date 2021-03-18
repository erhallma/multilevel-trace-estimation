function [mu,v,lvl,Nl] = mlmcTrace(A,f,lvl,npilot,vtol,nmax,d)
    % [mu,v,lvl,Nl] = mlmcTrace(A,f,lvl,npilot,vtol,ctol)
    %   Estimates tr(f(A)) using multilevel trace estimation
    %   INPUTS:
    %       A: a symmetric matrix or function handle, spectrum in [-1,1]
    %       f: a function
    %       lvl: a vector containing the levels to be used. If 
    %           lvl = [n], then n will be the maximum degree of the 
    %           interpolating Chebyshev polynomial and the other levels
    %           will be selected automatically using a pilot sample. 
    %       npilot: size of the pilot sample 
    %       vtol: target variance of the estimator
    %       nmax: determines computational budget to be nmax*n, 
    %           equivalent of nmax samples at highest level
    %       d: size(A,2). 
    %   OUTPUTS:
    %       mu : the trace estimate
    %       v  : the variance estimate
    %       lvl: the levels used
    %       Nl : the number of samples taken at each level
    
    if nargin < 5 || isempty(vtol), vtol = 0; end
    if nargin < 6 || isempty(nmax), nmax = 100; end
    
    if ~isa(A,'function_handle')
        d = size(A,2); 
        A = @(x) A*x; 
    end
    
    n = lvl(end); 
    lmax = length(lvl); 
    
    % find Chebyshev coefficients
    a = chebyFit(f,n); 
    
    % level selection if necessary
    if lmax == 1
        [lvl,Ns,sums] = levelSelection_pilot(A,d,a,npilot,vtol,nmax); 
        lmax = length(lvl); 
        Nl(1:lmax) = 0; Nl(end) = npilot; 
        suml(1:2,1:lmax) = 0; suml(:,end) = sums; 
        dNl = max(0,Ns-Nl); 
        dNl = ceil(dNl/3); % allow more time for variance estimates
        dNl = max(2,dNl); % use m >= 2 to allow for variance estimation
    else
        Nl(1:lmax) = 0; 
        suml(1:2,1:lmax) = 0;
        dNl(1:lmax) = npilot; 
    end
    
    Cl = cost_fun(lvl); 
    ctol = nmax*Cl(end); 

    % main loop
    while sum(dNl) > 0         
        % new samples
        for l = 1:lmax           
            if dNl(l) > 0
                al = a(1:lvl(l)+1); 
                if l > 1
                    al(1:lvl(l-1)+1) = 0; 
                end
                z = rademacher(d,dNl(l)); 
                s = chebyVal(A,z,al); 
                Nl(l)     = Nl(l) + dNl(l); 
                suml(1,l) = suml(1,l) + sum(s); 
                suml(2,l) = suml(2,l) + dot(s,s);           
            end
        end
        
        % update average and variance
        ml =        suml(1,:)./Nl; 
        Vl = max(0, suml(2,:)./Nl - ml.^2); 
        
        Vl = Vl.*(Nl)./(Nl-1); %unbiased variance estimator
        
        % update optimal number of samples
        rho    = sum(sqrt(Vl.*Cl)); 
        lambda = min(rho/vtol, ctol/rho);
        Ns     = lambda*sqrt(Vl./Cl);
        dNl    = max(0,Ns-Nl); 
        
        % don't go over budget 
        cRem = ctol - dot(Nl,Cl); 
        x    = min(1, cRem/dot(dNl,Cl));
        dNl  = max(0, ceil(dNl*x/2)); % more time for variance estimates
    end
    % end main loop
     
    
    % finally, estimate mean and standard error
    mu = sum(suml(1,:)./Nl);
    v  = sum(Vl./Nl);
end 


% -------------------
% Auxiliary functions
% -------------------

function a = chebyFit(f,n)
    % degree-n Chebyshev interpolation on interval [-1,1]
    x = cos(pi*(0:n)/n); 
    fx = feval(f,x)/(2*n); 
    g = real(fft(fx([1:n+1 n:-1:2]))); 
    a = [g(1), g(2:n)+g(2*n:-1:n+2), g(n+1)]; 
end

function z = rademacher(d,n)
    z = 2*(rand(d,n)>0.5)-1; 
end

function c = cost_fun(n)
    c = n; 
end

