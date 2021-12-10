% twoLevelFA.m
% Testing theoretical bounds for when 
%   2-level method outperforms 
%   single-level method
% Test problem: 
%   nuclear norm of FA
% Used to produce Figure 5.4 of revised paper
core

load FA.mat
load FA_SVD.mat

A = Problem.A; 
d = size(A,2);
sa = S.s; 
ea = [min(sa), max(sa)]; 
t = sum(sa); % the true nuclear norm

f = @sqrt; 
[Afun,ft] = shiftFun(A,f,ea,false);

n = 300; 
a = chebyFit(ft,n); 
p = @(x) chebyPol(x,a); 
b = 2/diff(ea.^2); 
c = sum(ea.^2)/diff(ea.^2); 
sx = b*sa.^2 - c; % Shifted eigenvalues of A'A, lying in [-1,1]


lmax = 150; %maximum degree of polynomial
LPAR = zeros(1,lmax); 
LPAR2 = zeros(1,lmax); 

normp2 = max(p(sx)); 
normpf = norm(p(sx)); 
r = normpf^2/normp2^2; % stable rank
ar = a; 

for l = 2:lmax
    % Split the Chebyshev polynomial in two
    ar(1:l+1) = 0; 
    qr = @(x) chebyPol(x,ar); 

    % Compute lambda
    l1 = max(qr(sx))/normp2; 
    l2 = norm(qr(sx))/normpf; 

    LPAR(l) = max(l1,l2); 

    % Compute lambda2
    LPAR2(l) = sqrt(d/r)*sum(abs(ar))/normp2; 

end

% Make a plot of theoretical values
alpha = 0.7; 
figure
semilogy(2:lmax,LPAR(2:end),'k','LineWidth',2), hold on
semilogy(2:lmax,LPAR2(2:end),'k--','LineWidth',2)
xlabel('Degree $$\ell$','FontSize',18,'Interpreter','latex')
ylabel('$$\lambda$$','fontsize',24,'Interpreter','latex')
yline(0.1,'color',[alpha,alpha,alpha],'LineWidth',4)

lgd = legend; 
lgd.String = {'$$\lambda$$', '$$\lambda_2^\prime$$'};
lgd.FontSize = 18; 
lgd.Interpreter = 'latex'; 
ax = gca; 
ax.FontSize = 14; 

n1 = find(LPAR<0.1,2,'first');
n2 = find(LPAR2<0.1,2,'first');

print('FAlambda','-dpng')

%%
% Now run some trials 
rng(1) 

m = 50; 
mpilot = 5; 
ls = [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 125, 150]; 
trials = 100; 

PML = zeros(3,length(ls)); 
mus = zeros(1,trials); 

for i = 1:length(ls)
    l = ls(i); 
    fprintf("l = %d\n", l)
    lvl = [l,n]; 
    for j = 1:trials
        mu = mlmcTrace(Afun,ft,lvl,mpilot,0,m,d);
        mus(j) = mu; 
    end
    errs = abs(mus - t)/t; 
    PML(:,i) = prctile(errs,[25,50,75]); 
end

fprintf("Single Level\n")
for j = 1:trials
    mu = singleLevel(Afun,d,a,m);
    mus(j) = mu; 
end
errs = abs(mus-t)/t; 
PSL = prctile(errs,[25,50,75]); 


%% Plot

figure
g1 = plot(ls,PML(2,:),'k','linewidth',2); hold on

%g2 = plot(ns,PSL(),'k--','linewidth',2);
f1 = fill([ls,fliplr(ls)],[PML(1,:),fliplr(PML(3,:))],'k');
f1.FaceAlpha = 0.2; 
f1.LineStyle = 'none'; 

g2 = yline(PSL(2),'k--','LineWidth',2);
%yline(PSL(1),'k:','LineWidth',2)
%yline(PSL(3),'k:','LineWidth',2)

% f2 = fill(ls([1,end,end,1]),PSL([1,1,3,3]),'k');
% f2.FaceAlpha = 0.2; 
% f2.LineStyle = 'none'; 

ax = gca; 
ax.FontSize = 14; 
xlabel('Degree $$\ell$','FontSize',18,'Interpreter','latex')
ylabel("Relative error of trace estimate", ...
        'FontSize',18,'Interpreter','latex')

lgd = legend([g1,g2]);
lgd.String = {'Two levels','Single level'}; 
lgd.Location = 'best'; 
lgd.FontSize = 16; 

print('FAtwoLevel','-dpng')