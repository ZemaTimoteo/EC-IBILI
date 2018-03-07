%% Statistic analysis

% Initialization
struct = G_AllSubs_perConds_FDR;
GCI_1 = struct.baseline;  % Condition 1
GCI_2 = struct.Nback1;    % Condition 2
[datachan,~,ntrials] = size(GCI_1);

% Espected values
avgGCI_1 = mean(GCI_1,3);
avgGCI_2 = mean(GCI_2,3);

% Difference between conditions
diff_GCI = avgGCI_1 - avgGCI_2;
diff_GCI(find(diff_GCI==0)) = NaN;

% ----- FWE correction - maximum permutation statistics --------
permut   = 1000; % number of perumtation

for i=1:permut
    % ... Building permutation matrix ...
    perGCI_1 = permute_GCI(GCI_1,GCI_2);
    perGCI_2 = permute_GCI(GCI_1,GCI_2);
    
    % ... Expected values of permutation matrix ...
    GCsur_scores_1 = mean(perGCI_1,3);
    GCsur_scores_2 = mean(perGCI_2,3);
    
    % ... Max value of the differences of the expected values ...
    diff_shufGC_Net    = GCsur_scores_1_Net - GCsur_scores_2_Net;
    maxGC_sur_Net(i)   = max(max(abs(diff_shufGC_Net))); % max value of shuffle distribution
end
    
    % ----- Plotting null distribution -----------------------------
    % ... threshold ...
    ite           = 0;
    thresh        = 0;
    maxGC_sur_Net = sort(maxGC_sur_Net);
    while thresh <0.90
        ite    = ite+1;    % ----- Family-wise Errors correction ---------------------------

        thresh = sum(maxGC_sur_Net(1:ite))/sum(maxGC_sur_Net);
    end
    thresh = maxGC_sur_Net(ite-1);
    % ... plot ...
    figure(21)
    histfit(maxGC_sur_Net);
    line([thresh thresh],[0 20],'Color','g'); % percentil of 0.95%
    
    % ----- Applying FWE to data ----------------------------------
    diff_GCI(find(diff_GCI<=thresh & diff_GCI>=-thresh)) = 0;
    figure(22)
    plot_pw(diff_GCI);
    title('GCI values: Condition 1 - Condition 2 (SingSubj)');
    caxis([min(min(diff_GCI)) max(max(diff_GCI))])
    colormap(c_map)

