% levelTestFA.m
% Testing the behavior of the multilevel method
%  using automated vs fixed level selection
% Test problem: 
%   nuclear norm of FA
% Used to produce Table 5.2 and Figures 5.1 and 5.2 in paper
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

n = 300; 
a = chebyFit(ft,n); 
m = 50; 
npilot = 10; 

L1 = [3 30 300]; 
L2 = [1:10,13,20,32,48,66,300]; 
L3 = [1:15,29,300]; 
 
nTrials = 100; 
LVLS = cell(3,nTrials); 
LVL = zeros(1,nTrials);

SDS = zeros(4,nTrials); 
ERR = zeros(4,nTrials); 
TS  = zeros(4,nTrials); 
EST = zeros(4,nTrials); 

rng(1)
for i = 1:nTrials
    fprintf("Trial %d\n", i)
    % Version 1: automated level selection
    tic
    [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,n,npilot,0,m,d);
    TS(1,i) = toc; 
    LVLS{1,i} = [lvl;Nl]; 
    LVL(i) = length(lvl); 
    
    SDS(1,i) = sqrt(vML); 
    ERR(1,i) = abs(muML-t)/t;  
    EST(1,i) = muML; 
    
    % Version 2: Straight single level
    tic
    [muSL,vSL] = singleLevel(Afun,d,a,m);
    TS(2,i) = toc; 
    SDS(2,i) = sqrt(vSL); 
    ERR(2,i) = abs(muSL-t)/t; 
    EST(2,i) = muSL; 
    
    % Version 3: Three levels
    tic
    [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,L1,npilot,0,m,d);
    LVLS{2,i} = [lvl;Nl];
    TS(3,i) = toc; 
    SDS(3,i) = sqrt(vML); 
    ERR(3,i) = abs(muML-t)/t;  
    EST(3,i) = muML; 
    
    % Version 4: Seventeen levels
    tic
    [muML,vML,lvl,Nl] = mlmcTrace(Afun,ft,L3,npilot,0,m,d);
    LVLS{3,i} = [lvl;Nl];
    TS(4,i) = toc; 
    SDS(4,i) = sqrt(vML); 
    ERR(4,i) = abs(muML-t)/t;  
    EST(4,i) = muML; 
end

% Fig 5.1, right: cumulative budget expenditure
figure, hold on
yline(m*n,'k--')
for i = 1:100
    L = LVLS{1,i}; 
    plot(L(1,:),cumsum(L(1,:).*L(2,:)),'k','linewidth',1); 
end
ax = gca; 
ax.FontSize = 14; 
xlabel('Polynomial degree','fontsize',16)
ylabel('Cumulative budget','fontsize',16)
print('FAbudget','-dpng')


% Fig 5.1, left: Histogram of the number of levels used
figure
h = histogram(LVL,'numBins',20,'FaceColor','k','FaceAlpha',0.1);
ax = gca; 
ax.FontSize = 14; 
xlabel('Number of levels','fontsize',16)
ylabel('Frequency','fontsize',16)
print('FAfreq','-dpng')

% Table 5.2: samples for one level in particular
L = LVLS{1,98}; 
disp(L)

% Look at second-largest levels
B = zeros(1,100); 
for i = 1:100
    L = LVLS{1,i}; 
    B(i) = L(1,end-1); 
end
fprintf("Median: %.2f, Max: %d\n", median(B), max(B)); 

% Fig 5.2: boxplot of the estimates over 100 trials
figure
boxplot(EST',{'Auto','L = 1','L = 3','L = 17'},'Color','k','Symbol','k+')
yline(t,'k--','linewidth',1)
ax = gca;
ax.FontSize = 14; 
xlabel('Levels','fontsize',16)
ylabel('Estimate','fontsize',16)
print('FAests','-dpng')

% Standard error estimates
std(EST,[],2) % taken over the 100 trials
median(SDS,2) % estimates of each trial
