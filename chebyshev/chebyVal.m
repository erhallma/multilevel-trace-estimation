function [s,c] = chebyVal(A,z,a)
    % [s,c] = chebyVal(A,z,a)
    % evaluates z'*p(A)*z for symmetric matrix A
    % INPUTS
    %  A: matrix or function handle
    %  z: initial vector or matrix of column vectors
    %  a: coefficients of Chebyshev polynomial
    % OUTPUTS
    %  s: The output of the bilinear form
    %  c: The individual coefficients z'*Tj(A)*z
    %      Rows of c correspond to columns of z
    
    if ~isa(A,'function_handle')
        A = @(x) A*x; 
    end
    
    c(:,1) = a(1)*dot(z,z)';    
    if length(a) == 1, s = c; return, end
    
    zj = A(z);  % z1
    zjm1 = z;   % z0
    c(:,2) = a(2)*dot(z,zj)'; 
    
    for j = 2:(length(a)-1)
        %next Chebyshev vector
        tmp  = 2*A(zj) - zjm1; 
        zjm1 = zj;         % z_{j-1}
        zj   = tmp;        % z_j
        
        % increment
        c(:,j+1) = a(j+1)*dot(z,zj)'; 
    end
    s = sum(c,2)'; 
end