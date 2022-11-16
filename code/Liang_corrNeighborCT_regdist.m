% =========================================================================
% This procedure is used to calculate the spatial correlation between the nodal CT maturation degree 
% and the mean of its directly connected neighbors.We regressed out the effect of nodal 
% mean Euclidean distance to its connected neighbors from the mean CT maturation degree.
%
% Written by Xinyuan Liang, xyliang@mail.bnu.edu.cn
% State Key Laboratory of Cognitive Neuroscience and Learning &
% IDG/McGovern Institute of Brain Research, 
% Beijing Normal University,
% Beijing, PR China.
% =========================================================================
ROI_num = {'125','250','500'};

for roi = 1:length(ROI_num)
    GroupSC_path = strcat('F:\data\CBDP\',ROI_num{roi},'\Group_sc.mat');
    Tvector_path = strcat('F:\data\CBDP\',ROI_num{roi} ,'\TVector.txt');
    dist_path = strcat('F:\data\CBDP\',ROI_num{roi},'\dist.mat');
    savepath = strcat('F:\data\CBDP\results\',ROI_num{roi});
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    
    load(GroupSC_path);
    load(Tvector_path);
    load(dist_path);
    
    [~,n] = size(G); % group sc size
    TVector = -TVector;% greater positive values indicate more significant cortical thinning

    mean_neighbor = []; % mean t-value of connected neighbors
    mean_delta_neighbor = []; % regressed
    mean_nei_dist = []; % mean Euclidean distance from a given node to its connected neighbors
    for i = 1:n
        Index = find(G(:,i));
        node_distance = dist(:,i);
        mean_neighbor(i,1) = mean(TVector(Index));
        mean_nei_dist(i,1) = mean(node_distance(Index));
    end
    stat = regstats(mean_neighbor,mean_nei_dist);
    mean_delta_neighbor = stat.r + stat.beta(1);
    
    stat2 = regstats(TVector,mean_delta_neighbor);
    if stat2.tstat.beta(2) >= 0
        r_adj = sqrt(stat2.adjrsquare);
    else
        r_adj = -sqrt(stat2.adjrsquare);
    end
    p = stat2.tstat.pval(2);
    save(fullfile(savepath,'r_wholebrain_regdist.mat'),'stat','stat2','r_adj','p','mean_delta_neighbor');
end

    
    
   