%% Longitudinal analysis, corr CT maturation rate with mean of sc-defined neighbors
load('/Data/CT_info.mat','delta_CT_rate','longitudinal_sub');
N_sub = length(longitudinal_sub);
SC_matrix = '/Data/trk_matrix/';
Dir_output = '/Data/r_emp/';
for i_sub = 1:N_sub
    SC_pathall{i_sub} = strcat(SC_matrix,longitudinal_sub{i_sub,1},'.mat');
    Output{i_sub} = strcat(Dir_output,longitudinal_sub{i_sub,1},'_r.mat');
end
rall = [];
mean_delta_neighbor_all = [];
for i_sub = 1:N_sub
    load(SC_pathall{i_sub},'G');
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
    save(Output{i_sub},'r','p','stat2','mean_delta_neighbor')
    rall(i_sub) = r;
    mean_delta_neighbor_all(:,i_sub) = mean_delta_neighbor;
end
save('/Data/rall.mat','rall');
save('/Data/mean_delta_neighbor_all.mat','mean_delta_neighbor_all');
