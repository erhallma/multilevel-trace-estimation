% estradaTests.m
% Estimating the Estrada Index of graphs
% taken from the SuiteSparse database
% Test cases: 
%   fe_4elt2, Erdos02, Roget
% Used to produce Table 5.5 in paper

core

tests = {'fe_4elt2','Erdos02','Roget'};
ns = [15,20,20]; 
eas = [ -3.33, 6.28;
        -21.4726, 25.8415;
        -6.4415, 12.0273];

f = @exp; 
m = 100;        % 100 samples
npilot = 10;    % 10 pilot samples

rng(1)
for i = 1:length(tests)
    test = tests{i}; 
    fprintf("Test: %s\n", test)
    
    load(test)
    A = Problem.A; 
    A = A + A'; 
    A(A~=0) = 1; 
    A = A - diag(diag(A)); 
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
    fprintf("Budget: %d\n", m*n)
    
end
