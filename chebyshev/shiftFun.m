function [Afun,fTilde] = shiftFun(A,f,ea,isSym,lambda)
    % Shifts A to implicitly have spectrum [-1,1]
    %  and shift the function f accordingly
    %  so that f(A) = fTilde(Afun)
    % INPUTS
    %  A: a matrix
    %  f: a function
    %  ea: bounds [lmin,lmax] or [smin,smax] on the spectrum of A
    %  isSym: true for A, false to use A'*A instead
    %  lambda: regularization, apply to (A+lambda*I) 
    %           or (A'*A + lambda*I) instead of A or A'*A
    %  NOTE: if this is done, then ea should be for the regularized 
    %           matrix rather than the original
    
    if nargin < 4 || isempty(isSym)
        isSym = issymmetric(A); 
    end
    if nargin < 5 || isempty(lambda)
        lambda = 0; 
    end

    % Create Afun, an implicit symmetric matrix
    %  with spectrum contained in [-1,1]
    %  fTilde(Afun) = f(A)
    if isSym
        b = 2/diff(ea); 
        c = sum(ea)/diff(ea); 
        Afun = @(z)  b*(A*z + lambda*z) - c*z;
    else
        b = 2/diff(ea.^2); 
        c = sum(ea.^2)/diff(ea.^2); 
        Afun = @(z) b*(A'*(A*z) + lambda*z) - c*z;
    end
    
    fTilde = @(x)f((x+c)/b); 
end