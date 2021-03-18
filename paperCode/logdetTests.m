% logdetTests.m
% Estimating the log determinant of matrices
% taken from the SuiteSparse database
% Test cases: 
%   thermotechTC, boneS01, ecology2, thermal2
% Used to produce Table 5.4 in paper

core

tests = {'thermomech_TC','boneS01','ecology2','thermal2'};
ns    = [75,150,60,60]; 
eas   = [ 4.5397e-04, 0.030549;
          0.00284726, 4.84756e4;
          1.2e-06, 80;
          1.6e-06, 7.8];

f = @log; 
m = 30;         % 30 samples
npilot = 5;     % 5 pilot samples

rng(1)
for i = 1:length(tests)
    test = tests{i}; 
    fprintf("Test: %s\n", test)
    
    load(test)
    A = Problem.A; 
    d = size(A,2); 
    
    ea = eas(i,:); 
    [Afun,ft] = shiftFun(A,f,ea); 
    
    n = ns(i);
    a = chebyFit(ft,n);
    
    % multilevel
    tic 
    [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,n,npilot,0,m,d);
    tE = toc; 
    fprintf("ML time: %.4f\n", tE)
    fprintf("ML mean: %.8e\n", muML)
    fprintf("ML serr: %.8e\n", sqrt(vML))
    fprintf("Budget: %d\n", dot(lvl,Nl))
    
    % single level
    tic
    [muSL,vSL] = singleLevel(Afun,d,a,m);
    tE = toc;
    fprintf("SL Time: %.4f\n", tE)
    fprintf("SL mean: %.8e\n", muSL)
    fprintf("SL serr: %.8e\n", sqrt(vSL))
    fprintf("Budget: %d\n\n", m*n)
    
end
