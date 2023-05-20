% ===================================================================================================
% This procedure is used to calculate the spatial correlation between the nodal CT maturation degree 
% and the mean of its directly connected neighbors.
% Written by Xinyuan Liang, xyliang@mail.bnu.edu.cn
% State Key Laboratory of Cognitive Neuroscience and Learning &
% IDG/McGovern Institute of Brain Research, 
% Beijing Normal University,
% Beijing, PR China.
% ===================================================================================================

ROI_num = {'125','250','500'};
for roi = 1:length(ROI_num)
    GroupSC_path = strcat('F:\data\CBDP\',ROI_num{roi},'\Group_sc.mat');
    Tvector_path = strcat('F:\data\CBDP\',ROI_num{roi} ,'\TVector.txt');
    savepath = strcat('F:\data\CBDP\results\',ROI_num{roi});
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    
    load(GroupSC_path);
    load(Tvector_path);
    [~,n] = size(G); % group sc size
    TVector = -TVector; % greater positive values indicate more significant cortical thinning
    mean_delta_neighbor = [];
    for i = 1:n
        Index = find(G(:,i));
        mean_delta_neighbor(i,1) = mean(TVector(Index));
    end
    stat1 = regstats(TVector,mean_delta_neighbor);
    if stat1.tstat.beta(2) >= 0
        r_adj = sqrt(stat1.adjrsquare);
    else
        r_adj = -sqrt(stat1.adjrsquare);
    end
    p = stat1.tstat.pval(2);
    save(fullfile(savepath,'r_wholebrain.mat'),'r_adj','p','mean_delta_neighbor');
end

    



    
    
   
