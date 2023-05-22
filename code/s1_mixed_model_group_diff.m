%===========================================================================
% This procedure is used to estimated the statistical differences in CT 
% between the child and adolescent groups by a mixed linear analysis.
%===========================================================================

    Roinum = '500';
    load('child_info.mat');
    load('cortical_thickness.mat');
    [~,ROInumber] = size(cortical_thickness);
    for i = 1:ROInumber
        table_model.depen_var = cortical_thickness(:,i);
        lme1 = fitlme(table_model,'depen_var ~ 1 + group + sex + (1|subname) + (-1 + group|subname)');
        group_pValue(i) = lme1.Coefficients.pValue(2);
        group_beta(i) = lme1.Coefficients.Estimate(2);
        group_tValue(i) = lme1.Coefficients.tStat(2);
    end
    group_pValue = group_pValue';
    group_beta = group_beta';
    group_tValue = group_tValue';
    Outpath = strcat('data\CBDP\',Roinum);
    Outmat = strcat(Outpath,'group_diff.mat');
    save (Outmat, 'group_tValue', 'group_pValue','group_beta');
    outgroup_tValue_txt = strcat(Outpath,'TVector.txt');
    save (outgroup_tValue_txt, 'group_tValue', '-ascii');
    outgroup_pValue_txt = strcat(Outpath,'PVector.txt');
    save (outgroup_pValue_txt, 'group_pValue', '-ascii');
    outgroup_betaValue_txt = strcat(Outpath,'betaVector.txt');
    save (outgroup_betaValue_txt, 'group_beta', '-ascii');

