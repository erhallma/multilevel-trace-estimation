% degreeTestFA.m
% Testing the behavior of the multilevel method
%  and single level method using different 
%  degrees for the Chebyshev approximation
% Test problem: 
%   nuclear norm of FA
% Used to produce Figure 5.3 in paper
core

load FA.mat
load FA_SVD.mat

A = Problem.A; 
d = size(A,2);
sa = S.s; 

lmb = 1e-10; 
lambda = lmb*max(sa); 

sb = sqrt(sa.^2 + lambda); 
ea = [min(sb), max(sb)]; 
t = sum(sa); % the true nuclear norm

f = @sqrt; 
[Afun,ft] = shiftFun(A,f,ea,false,lambda);

ns = [25,50,75,100,150,200,250,300,350]; 
nk = length(ns); 
m = 50; 
npilot = 10; 
 
nTrials = 100; 

LVLS = cell(nk,nTrials); 
TS   = zeros(nk,nTrials); 
SDS  = zeros(nk,nTrials); 
ERR  = zeros(nk,nTrials);  
EST  = zeros(nk,nTrials); 

TSL   = zeros(nk,nTrials); 
SDSL  = zeros(nk,nTrials); 
ERRSL = zeros(nk,nTrials);  
ESTSL = zeros(nk,nTrials); 

PS = zeros(1,nk); 
ES = zeros(1,nk); 

% run the trials
rng(1) 
for j = 1:length(ns)
    n = ns(j); 
    fprintf("Degree %d\n", n)
    a = chebyFit(ft,n); 
    
    p = @(x) chebyPol(x,a); 
    b = 2/diff(ea.^2); 
    c = sum(ea.^2)/diff(ea.^2); 
    sc = b*sb.^2 - c; 
    pn = sum(p(sc)); 
    PS(j) = pn; 
    ES(j) = abs(pn-t)/t; %relative errors of Chebyshev approximation

    for i = 1:nTrials
        fprintf("Trial %d\n", i)
        % Version 1: automated level selection
        tic
        [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,n,npilot,0,m,d);
        TS(j,i) = toc; 
        LVLS{j,i} = [lvl;Nl]; 
        SDS(j,i) = sqrt(vML); 
        ERR(j,i) = abs(muML-t)/t;  
        EST(j,i) = muML; 

        % Version 2: single level method
        tic
        [muSL,vSL] = singleLevel(Afun,d,a,m);
        TSL(j,i) = toc; 
        SDSL(j,i) = sqrt(vSL); 
        ERRSL(j,i) = abs(muSL-t)/t; 
        ESTSL(j,i) = muSL; 

    end
end

% Figure 5.3, right: plot of the errors
PML = prctile(ERR,[25 50 75],2); 
PSL = prctile(ERRSL,[25 50 75],2); 
figure
g1 = semilogy(ns,PML(:,2),'k','linewidth',2); hold on
g2 = plot(ns,PSL(:,2),'k--','linewidth',2);
f1 = fill([ns,fliplr(ns)],[PML(:,1);flipud(PML(:,3))],'k');
f1.FaceAlpha = 0.2; 
f1.LineStyle = 'none'; 

f2 = fill([ns,fliplr(ns)],[PSL(:,1);flipud(PSL(:,3))],'k');
f2.FaceAlpha = 0.2; 
f2.LineStyle = 'none'; 

g3 = plot(ns,ES,'k:','linewidth',2);

ax = gca; 
ax.FontSize = 14; 
xlabel("Polynomial degree", 'FontSize',16)
ylabel("Relative error of trace estimate", 'FontSize',16)
lgd = legend([g1,g2,g3]); 
lgd.FontSize = 16; 
lgd.String = {'Multilevel','Single level','Exact error'}; 

print('FAdegree1','-dpng')

% Figure 5.3, left: plot of standard errors
figure
plot(ns,std(EST,[],2),'k','linewidth',2), hold on
plot(ns,std(ESTSL,[],2),'k--','linewidth',2)
ax = gca; 
ax.FontSize = 14; 
yl = ax.YLim; 
plot(ns,ES*t,'k:','linewidth',2)
ylim([0,max(yl)])
xlabel("Polynomial degree", 'FontSize',16)
ylabel("Approximate standard error", 'FontSize',16)
lgd = legend('Multilevel','Single level','Exact error'); 
lgd.FontSize = 16; 
lgd.Location = 'southwest'; 

print('FAdegree2','-dpng')
