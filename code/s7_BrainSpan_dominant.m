load('data\gene\Brainspan_neocortex.mat');
load('data\gene\neurodev_process.mat');
load('data\gene\dominant_regions.mat');
 
% dominant region
Index_dom = find([dominant{:,2}]);
dom_structure = {dominant{Index_dom,1}};
notdom_structure = setdiff(dominant(:,1),dom_structure);

diff_PCperm = []; % Differences between the first principal components (permutation)
diff_PC = []; % Differences between the first principal components (observed)
field = fieldnames(neurodev_process); % four neurodevelopment processes
for i = 1:length(field)
    disp(i)
    proc = field{i};
    gene_proc = getfield(neurodev_process, proc);
    gene_remain_Index = ismember(brainspan_genename, gene_proc);
    domIndex = ismember(sample_info.structure_name, dom_structure);
    notdomIndex = ismember(sample_info.structure_name, notdom_structure);
    domGene = brainspan(gene_remain_Index, domIndex); % gene * dominant samples
    notdomGene = brainspan(gene_remain_Index, notdomIndex);
    Data = zscore( [ domGene'; notdomGene' ] ); 

    [ ~, PC ] = pca( Data, 'Centered', false );
    if corr( PC( :, 1 ), mean( Data, 2 ) ) < 0
        PC = -PC;
    end
    PC = PC(:,1);
    
%     PC = ( PC( :, 1 ) - min( PC( :, 1 ) ) )/( max( PC( :, 1 ) ) - min( PC( :, 1 ) ) ); % for plot
    domAge = log2(sample_info.Day(domIndex));
    notdomAge = log2(sample_info.Day(notdomIndex));
    domPC = PC(1:length(domAge));
    notdomPC = PC((length(domAge)+1):end); 
%      csvwrite( [ 'data\results\Brainspan', filesep, 'DevelopmentTrajectory-', proc, '-dom.csv' ], [ domAge'; domPC' ] );
%      csvwrite( [ 'data\results\Brainspan',filesep,'DevelopmentTrajectory-', proc, '-notdom.csv' ], [ notdomAge'; notdomPC' ] );

    % PC diff between 6y-14y
    Index_dom_6_14 = find(domAge > log2(2456) & domAge < log2(5376)); 
    Index_notdom_6_14 = find(notdomAge > log2(2456) & notdomAge < log2(5376)); 

    domPC_6_14 = mean(domPC(Index_dom_6_14));
    notdomPC_6_14 = mean(notdomPC(Index_notdom_6_14));
    diff_PC(i) = domPC_6_14 - notdomPC_6_14;

    % permutation
    Index_gene_per = find(gene_remain_Index == 0);
    generemain_num = length(find(gene_remain_Index == 1)); % overlay gene number
    nperm = 1000; % permutation number

    for j = 1 : nperm
        myresample = randsample(Index_gene_per,generemain_num);
        domGene_perm = brainspan(myresample,domIndex);
        notdomGene_perm = brainspan(myresample,notdomIndex);
        Dataperm = zscore( [ domGene_perm'; notdomGene_perm' ] );
        [ ~, PCperm ] = pca( Dataperm, 'Centered', false );
        if corr( PCperm( :, 1 ), mean( Dataperm, 2 ) ) < 0
            PCperm = -PCperm;
        end
        PCperm = PCperm(:,1);
        domPCperm = PCperm(1:length(domAge));
        notdomPCperm = PCperm((length(domAge) + 1):end); 
        domPCperm_6_14 = mean(domPCperm(Index_dom_6_14));
        notdomPCperm_6_14 = mean(notdomPCperm(Index_notdom_6_14));
        diff_PCperm(j,i) = domPCperm_6_14 - notdomPCperm_6_14;
    end
    if diff_PC(i) > 0
        p_perm(i) = length(find(diff_PCperm(:,i) > diff_PC(i))) / nperm;
    else
        p_perm(i) = length(find(diff_PCperm(:,i) < diff_PC(i))) / nperm;
    end
    
end


