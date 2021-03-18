function s = chebyPol(xs,a)
    % s = chebyPol(xs,a)
    % evaluates Chebyshev polynomial at a given set of points
    % INPUTS
    %  xs: points to evaluate
    %  a: vector of Chebyshev coefficients ordered as [0,...,n]
    % OUTPUTS
    %  s: p(xs), where p() is the polynomial in the Chebyshev basis
    
    s = a(1)*ones(size(xs)); 
    if length(a) == 1, return, end
    
    x = xs; 
    s = s + a(2)*x; 
    z = 1; 
    
    for i = 3:length(a)
        tmp = 2*x.*xs-z; 
        z = x; 
        x = tmp; 
        
        s = s + a(i)*x; 
    end
end