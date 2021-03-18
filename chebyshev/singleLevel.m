function [muEst,vEst] = singleLevel(Afun,d,a,N)
    % [m,v] = singleLevel(Afun,d,a,N)
    % Stochastic trace estimation for a symmetric matrix
    % INPUTS: 
    %  Afun: a function handle for a matrix A, Afun(x) = A*x
    %  d: the number of columns of A
    %  a: a vector of Chebyshev coefficients ordered [0,...,n]
    %  N: the number of samples
    % OUTPUTS:
    %  m: the sample mean
    %  v: an estimate of the sample variance 
    %       (sqrt(v) estimates the standard error)
    
    z = rademacher(d,N); 
    ys = chebyVal(Afun,z,a); 
    
    muEst = mean(ys); 
    vEst = var(ys)/N; 
end


% Auxiliary functions
function z = rademacher(d,n)
    z = 2*(rand(d,n)>0.5)-1; 
end
