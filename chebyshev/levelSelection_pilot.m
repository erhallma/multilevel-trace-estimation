function [L,Ns,sums] = levelSelection_pilot(A,d,a,npilot,vtol,nmax)
    % [L,Ns,sums] = levelSelection_pilot(A,d,a,npilot,vtol,nmax)
    % Automated level selection for MLMC trace estimation
    %  on the basis of a pilot sample
    % INPUTS 
    %  A: Symmetric matrix or function handle
    %  d: size(A,2)
    %  a: a vector [a0,...,an] of Chebyshev coefficients
    %  npilot: number of pilot samples
    %  vtol: target variance
    %  nmax: computational budget in terms of highest level
    % OUTPUTS
    %  L: a vector of the selected levels [l1,...,lk]
    %  Nl: recommended number of samples at each level
    %  sums: sum(x) and sum(x.^2) at pilot level lk. 
    
    n = length(a) - 1; 
    C = zeros(1,n+1); 
    T = zeros(1,n+1);
    V = zeros(1,n+1); 
    sums = zeros(1,2); 
    
    ctol = nmax*cost_fun(n); 
    
    % Pilot sample
    z = rademacher(d,npilot); 
    [~,X] = chebyVal(A,z,a); 
        
    for j = 0:n
        csx = cumsum(X(:,1:j+1),2,'reverse');
        varj = var(csx); 
        vx = sqrt(varj*cost_fun(j)); 
        cx = vx + [0,C(1:j)]; 
        if j == n % extra constraint at final step
            mus = min(cx/vtol,ctol./cx); 
            mx  = (vx.*mus)/cost_fun(n);
            cx(mx<npilot) = Inf; 
        end
        cx(2) = Inf;        % avoid level 0. 
        [m,ix] = min(cx); 
        C(j+1) = m; 
        T(j+1) = ix-2; 
        V(j+1) = vx(ix); 
    end
    pilotVals = csx(:,ix); 
    sums(1) = sum(pilotVals); 
    sums(2) = sum(pilotVals.^2); 

    
    % Backtrack to find optimal levels
    L = n; 
    ix = T(n+1); 
    while ix > -1
        L(end+1) = ix; %#ok<AGROW>
        ix = T(ix+1); 
    end
    L = fliplr(L); % vector of optimal levels
    
    % Recommended number of samples
    V = V(L+1); 
    mu = min(sum(V)/vtol,ctol/sum(V)); 
    Cl = cost_fun(L); 
    Ns = mu * (V./Cl); 
          
end


% -------------------
% Auxiliary functions
% -------------------

function z = rademacher(d,n)
    z = 2*(rand(d,n)>0.5)-1; 
end

function c = cost_fun(n)
    c = n; 
end
