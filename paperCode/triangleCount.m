% triangleCount.m
% Variance reduction techniques for 
% estimating tr(A^3)
% Test cases: 
%   ca-GrQc, wiki-Vote
% Used to produce Figure 5.4 in paper

core
rng(1)

tests = {'ca-GrQc.txt','wiki-Vote.txt'}; 
names = {'ca-GrQc','wiki-Vote'};

for ix = 1:2
    test = tests{ix};
    name = names{ix}; 

    T = readtable(test,'ReadVariableNames',false); 
    G = graph(T{:,1},T{:,2}); 
    Alarge = adjacency(G); 

    s = sum(Alarge); 
    v = (s~=0); 
    A = Alarge(v,v); 

    for i = 1:size(A)
        A(i,i) = 0; 
    end


    % Compute the true answer
    fprintf("Computing trace...\n")
    tic
    A2 = A^2; 
    t = dot(A(:),A2(:)); 
    toc
    x = [trace(A); nnz(A)];

    % Generate random vectors
    m = 1000; 
    ntrials = 100; 
    N = size(A,2);

    est  = zeros(ntrials,m); 
    estx = zeros(ntrials,m); 

    tic
    for itn = 1:ntrials
        Z = rademacher(N,m); 

        Z1 = A*Z;
        C = [dot(Z,Z1); dot(Z1,Z1)]; 
        b = dot(Z1,A*Z1); 

        est(itn,:)  = cumsum(b)./(1:m); 

        for i = 1:m % control variates
            Ci = C(:,1:i);
            c3i = b(1:i); 
            alpha = c3i/Ci; 
            c3i = c3i - alpha*Ci; 
            estx(itn,i) = mean(c3i) + dot(alpha,x); 
        end
    end
    toc

    % Error analysis and plots
    err  = abs(est-t)/t; 
    errx = abs(estx-t)/t; 

    pct  = prctile(err,  [25,50,75]); 
    pctx = prctile(errx, [25,50,75]); 
    
    xpts = [2:10,20:10:100,200:100:1000];
    zpts = fliplr(xpts); 

    figure
    loglog(xpts,pct(2,xpts),'k--','linewidth',2), hold on
    plot(xpts,pctx(2,xpts),'k','linewidth',2)
    f = fill([xpts,zpts],[pct(1,xpts),pct(3,zpts)],'k'); 
    f.FaceAlpha = 0.1; 
    f.LineStyle = ':'; 
    xlim([2,m])

    fx = fill([xpts,zpts],[pctx(1,xpts),pctx(3,zpts)],'k');
    fx.FaceAlpha = 0.1; 
    fx.LineStyle = '-'; 

    xlabel('Number of samples $m$','fontsize',16','interpreter','latex')
    ylabel('Relative error of trace estimate','fontsize',16, ...
                'interpreter','latex')

    lgd = legend("Hutchinson's","Control Variates"); 
    lgd.FontSize = 16; 

    print(name,'-dpng');
    
    
    fprintf("Median variance ratio: %.4f\n", median(pctx(2,:)./pct(2,:)))
end


% -------------------
% Auxiliary functions
% -------------------

function y = rademacher(n,d)
    y = 2*(rand(n,d)>0.5)-1; 
end
