% nuclearTests.m
% Estimating the nuclear norm of matrices
% taken from the SuiteSparse database
% Test cases: 
%   California, FA, Erdos02, fe_4elt2, deter3, ukerbe1
% Used to produce Table 5.3 in paper


tests = {'California','FA','Erdos02','fe_4elt2','deter3','ukerbe1'}; 
%lmb   = 1e-10; % for arXiV version
lmb = 0; % for LAA version
f     = @sqrt; 

ns     = [100,300,100,70,70,20];
eas =  [0, 21.56;
        0, 3.424;
        0, 25.842; 
        2.87e-05, 6.28;
        0.033, 10.30; 
        0, 3.1314]; 
     
m      = 50; 
nPilot = 10; 

rng(1)
for ix = 1:length(tests)
    test = tests{ix}; 
    load(test)
    
    A = Problem.A; 
    if size(A,2) > size(A,1)
        A = A'; 
    end
    d = size(A,2); 
    
    % regularize if the matrix is low-rank
    ea = eas(ix,:); 
    if min(ea) == 0
        lambda = lmb*max(ea); 
    else
        lambda = 0; 
    end
    ea = sqrt(ea.^2 + lambda); 

    [Afun,ft] = shiftFun(A,f,ea,false,lambda); 
    n = ns(ix); 
    a = chebyFit(ft,n); 
    
    % multilevel
    tic
    [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,n,nPilot,0,m,d); 
    tE = toc; 
    
    disp([lvl;Nl])
    fprintf("Time Elapsed: %.2f\n", tE)
    fprintf("Nuclear Norm Estimate: %.4f\n", muML)
    fprintf("Standard Error: %.4f\n", sqrt(vML))
    
    % single level
    tic
    [muSL,vSL] = singleLevel(Afun,d,a,m); 
    tE = toc;
    fprintf("Time Elapsed: %.2f\n", tE)
    fprintf("Nuclear Norm Estimate: %.4f\n", muSL)
    fprintf("Standard Error: %.4f\n", sqrt(vSL))
    
end
