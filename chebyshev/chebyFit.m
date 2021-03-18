function a = chebyFit(f,n)
    % a = chebyFit(f,n)
    % Chebyshev interpolation on interval [-1,1]
    % INPUTS
    %  f: a function
    %  n: degree of interpolating polynomial
    % OUTPUTS
    %  a: row vector of coefficients [0,...,n]
    
    x = cos(pi*(0:n)/n); 
    fx = feval(f,x)/(2*n); 
    g = real(fft(fx([1:n+1 n:-1:2]))); 
    a = [g(1), g(2:n)+g(2*n:-1:n+2), g(n+1)]; 
end