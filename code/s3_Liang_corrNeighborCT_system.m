ROI_num = '500';
load(strcat('data\CBDP\',ROI_num,'\Group_sc.mat'));
load(strcat('data\CBDP\',ROI_num ,'\TVector.txt'));
TVector = -TVector;% greater positive values indicate more significant cortical thinning

% load Yeo7 and Mesulam label
load(strcat('data\CBDP\',ROI_num,'\Yeo7_label_',ROI_num,'.mat'));
load(strcat('data\CBDP\',ROI_num,'\Mesulam_',ROI_num,'_label.mat'));

Outpath = strcat('data\results\CBDP\',ROI_num);

[nodenum,~] = size(G);
nei_Yeo = zeros(nodenum,7);
nei_sys4 = zeros(nodenum,4);

%% yeo7
for i = 1:nodenum
    index = find(G(:,i));
    neighbor_T = TVector(index);
    neighbor_label = Yeo7_label(index);
    for j = 1:7
        Index_label = find(neighbor_label == j);
        nei_Yeo(i,j) = mean(neighbor_T(Index_label));
    end
end

for i = 1:7
    Index = find(Yeo7_label == i);
    T = TVector(Index); % T vectors for ith system
    for j = 1:7
        nei = nei_Yeo(Index,j);
        if length(find(~isnan(nei))) < 10
            r_yeo(i,j) = 0;
            p_yeo(i,j) = 1;
        else
           [r,p] = corr(T,nei,'rows','complete');
           data_num = length(T) - length(find(isnan(nei)));
           adjrs = 1-(1-r^2) * (data_num-1) / (data_num-2); %adjust rsquare
           if adjrs < 0
               adjrs = 0; 
           end
           
           if r > 0
               r_yeo(i,j) = sqrt(adjrs);
           else
               r_yeo(i,j) = -sqrt(adjrs);
           end  
           p_yeo(i,j) = p;
        end
    end
end
      
%% Mesulam
for i = 1:nodenum
    index = find(G(:,i));
    neighbor_T = TVector(index);
    neighbor_label = label(index);
    for j = 1:4
        Index_label = find(neighbor_label == j);
        nei_sys4(i,j) = mean(neighbor_T(Index_label));
    end
end

for i = 1:4
    Index = find(label == i);
    T = TVector(Index);
    for j = 1:4
        nei = nei_sys4(Index,j);
        if length(find(~isnan(nei))) < 10
            r_sys4(i,j) = 0;
            p_sys4(i,j) = 1;
        else
            [r,p] = corr(T,nei,'rows','complete');
            data_num = length(T) - length(find(isnan(nei)));
            adjrs = 1-(1-r^2)*(data_num-1) / (data_num-2);
            if adjrs < 0
               adjrs = 0; 
            end 
            
            if r > 0
                r_sys4(i,j) = sqrt(adjrs);
            else
                r_sys4(i,j) = -sqrt(adjrs);
            end  
           p_sys4(i,j) = p;
        end
    end
end

save(strcat(Outpath,'\','r_sys.mat'),'r_sys4','p_sys4','nei_sys4','r_yeo','p_yeo','nei_Yeo');   
