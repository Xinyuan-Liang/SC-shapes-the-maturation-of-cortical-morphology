% ================================================================================================
% Using a random walk procedure to describe the diffusion properties of the white matter network 
% Written by Xinyuan Liang, xyliang@mail.bnu.edu.cn
% State Key Laboratory of Cognitive Neuroscience and Learning &
% IDG/McGovern Institute of Brain Research, 
% Beijing Normal University,
% Beijing, PR China.
% ================================================================================================
%get transition_probabilities Discrete time
function [TP] = s6_Liang_random_walk(GroupSC_path,step_num)
    load(GroupSC_path);
    degree = sum(G); % node degree
    node = length(degree); % node number
    D = diag(degree);
    P = inv(D) * G; % transition matrix
    p0 = eye(node);
    ini = p0;% initial distribution of random walkers
    TP = [];% transition probability 

    for i = 1:step_num
        ini = ini * P;
        TP(i,:,:) = ini;
    end
end

