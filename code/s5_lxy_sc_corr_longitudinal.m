%% Longitudinal analysis, corr CT maturation rate with mean of sc-defined neighbors
load('./data/CT_info.mat','delta_CT_rate');
load('./data/SC_individual.mat');
[N_sub,~] = size(delta_CT_rate);

rall = [];
mean_delta_neighbor_all = [];
for i_sub = 1:N_sub
    G = A(:,:,i_sub);
    [~,n]=size(G);
    mean_delta_neighbor=[];
    CT_rate_subj = delta_CT_rate(i_sub,:)';
    
    for i = 1:n
        Index=find(G(:,i));
        mean_delta_neighbor(i,1) = mean(CT_rate_subj(Index));
    end
    
    stat2 = regstats(CT_rate_subj,mean_delta_neighbor);
    if stat2.tstat.beta(2) >= 0
        r = sqrt(stat2.adjrsquare);
    else
        r = -sqrt(stat2.adjrsquare);
    end
    p = stat2.tstat.pval(2);
    rall(i_sub) = r;
    mean_delta_neighbor_all(:,i_sub) = mean_delta_neighbor;
end
save('./data/rall.mat','rall');
save('./data/mean_delta_neighbor_all.mat','mean_delta_neighbor_all');
