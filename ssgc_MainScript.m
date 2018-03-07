%%% ========= CONNECTIVITY ANALYSIS SCRIPT ========= %%%
%%% ============= SSGC ============= %%%
%%% ======================================== %%%

% Remember to have NeuroElf installed

%% Load Workspace
% % % clear;
% % % load('workspace.mat');
% % % i=1


%% Configuration
close all;

% addpath('mvgc_new')

%% Initalization and First Subject Selection

% subjectslist = importdata('danielalist.mat');
% n_sub = length(subjectslist);


%% Group VOI file
% txtPath = 'C:\Users\alexa\Desktop\DATA_Daniela\ANALYSIS\VOIClusterPeakTable_1back2backvsBaseline_RFX_FDR005_max300';
% peaks = importClusterPeakTable([ txtPath '.txt'] );
% voiFile = xff(fullfile(dataRootAnalysis,'GroupSphericalVOIs_4mm.voi'));
voiFile = xff('GroupVOIs_RFX_FDR005_17sub_plus_plus.voi');

n_rois = voiFile.NrOfVOIs;
voiNames = voiFile.VOINames;
n_con = (n_rois*(n_rois-1))/2; % #Connections



nBlocks = 1;
[nRois,blockDur,nSubj] = size(DataMVGC_AllCond);
try
    ResultsSSGC = ssgc_demo_function( blockDur , ... % number of observations per trial
        DataMVGC_AllCond , ... % data
        nRois , ... % number of rois
        0.05 , ... % significance level for significance test
        1 , ... % maximum model order for model order estimation
        'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
        nSubj , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
        nBlocks, ... % number of blocks per conditions
        cond_names{cond} , ... % condition names (or plot title in case of #conds = 1)
        0 , ... % Display model order plot
        0 , ... % Display results plot
        '', ... % 'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        0); % checkConsistency (0 without check and 1 with check -mandatory nnrois<blocksize )
    
    G_AllSubs_perConds_pw_FDR = ResultsSSGC.G;
    
    SSsig_AllSubs_perConds_pw_FDR = ResultsSSGC.sig;
    
    
catch ME
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCM failed')
end



%% State-space - SSGC toolbox all ROIS - each conditions - without blocks
failures = 0;  % blockDur(1) =15;
DataMVGC_II = DataMVGC;
for cond = 1:n_cond
    fprintf('Performing SSGC for condition %s... \n',cond_names{cond});
    [nnrois,~,nblocks]= size(DataMVGC_II.(cond_names{cond}));
    % normalization 1
    DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*nblocks);
    DataMVGC_II.(cond_names{cond}) = DataMVGC_II.(cond_names{cond}) - (ones(size(DataMVGC_II.(cond_names{cond})'))*diag(mean(DataMVGC_II.(cond_names{cond}),2)))';
    % normalization 2
%     DataMVGC_II.(cond_names{cond}) = demean(DataMVGC_II.(cond_names{cond})); % no constant term
%     DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*nblocks);
    
    [~,blockDur(cond),nblocks] = size(DataMVGC_II.(cond_names{cond}));
     sizeblocks = blockDur(cond); %1000 = sizeblocks;
%     try        
        ResultsSSGC = ssgc_demo_function( sizeblocks , ... % number of observations per trial
            DataMVGC_II.(cond_names{cond})(:,1:sizeblocks) , ... % data (alteraçao do nºpontos)
            nnrois , ... % number of rois
            0.05 , ... % significance level for significance test
            1 , ... % maximum model order for model order estimation
            'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
            [] , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
            nblocks, ... % number of blocks per conditions
            cond_names{cond} , ... % condition names (or plot title in case of #conds = 1)
            0 , ... % Display model order plot
            0 , ... % Display results plot
            '', ... % 'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
            0); % checkConsistency (0 without check and 1 with check -mandatory nnrois<blocksize )

        G_AllSubs_perConds_FDR.(cond_names{cond}) = ResultsSSGC.G;        
        SSsig_AllSubs_perConds_FDR.(cond_names{cond}) = ResultsSSGC.sig;

        
        
%     catch ME
%         disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCM failed')
%         failures = failures + 1;
%     end
end

% Plot baseline & C1back & C2back
figure(); clf;
subplot(1,3,1);
plot_pw(SSsig_AllSubs_perConds_FDR.baseline);
title('Pairwise-conditional SSGC - baseline');
subplot(1,3,2);
plot_pw(SSsig_AllSubs_perConds_FDR.C1back);
title('Pairwise-conditional SSGC - C1Back');
subplot(1,3,3);
plot_pw(SSsig_AllSubs_perConds_FDR.C2back);
title('Pairwise-conditional SSGC - C2Back');

% Plot baseline - C1back - C2back
SSsig_AllSubs_perConds_FDR.combination = ... % SSsig_AllSubs_perConds_FDR.baseline + ...
                                         SSsig_AllSubs_perConds_FDR.C1back + ...
                                         SSsig_AllSubs_perConds_FDR.C2back;
figure()
plot_pw(SSsig_AllSubs_perConds_FDR.combination);
title('Pairwise-conditional SSGC - Combination');


% % % figure(); clf;
% % % subplot(1,2,1);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.Baseline);
% % % title('Pairwise-conditional SSGC - Baseline');
% % % subplot(1,2,2);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.Imaginacao);
% % % title('Pairwise-conditional SSGC - Imaginacao');
% % % 
% % % 
% % % % Plot baseline - C1back - C2back
% % % SSsig_AllSubs_perConds_FDR.combination = SSsig_AllSubs_perConds_FDR.Baseline + ...
% % %                                          SSsig_AllSubs_perConds_FDR.Imaginacao;
% % %                                      
% % % figure(); clf;
% % % plot_pw(SSsig_AllSubs_perConds_FDR.combination);
% % % title('Pairwise-conditional SSGC - Combination');

%% State-space - SSGC toolbox all ROIS - each conditions - with blocks
failures = 0; % blockDur(1) =15; % n_cond = 1;
DataMVGC_II = DataMVGC;
clear G_AllSubs_perConds_FDR ResultsSSGC SSsig_AllSubs_perConds_FDR
for cond = 1:n_cond
    fprintf('Performing SSGC for condition %s... \n',cond_names{cond});
    [nnrois,blockDur(cond),nblocks]= size(DataMVGC_II.(cond_names{cond}));
    % normalization
    DataMVGC_II.(cond_names{cond}) = demean(DataMVGC_II.(cond_names{cond})); % no constant term
    spFactor = 40; aux = (blockDur(cond)*spFactor)*(floor(nblocks/spFactor));
    DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*nblocks);
    DataMVGC_II.(cond_names{cond})(:,aux+1:end) = [];
    DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*spFactor,floor(nblocks/spFactor));
    [~,blockDur(cond),nblocks] = size(DataMVGC_II.(cond_names{cond}));
    sizeblocks = blockDur(cond);
    
    for blo=1:nblocks
        fprintf('.. Block number: %d \n',blo);
        
        try
            ResultsSSGC = ssgc_demo_function_v2( blo, ...  % number of block
                DataMVGC_II.(cond_names{cond})(:,1:sizeblocks,blo) , ... % data (alteraçao do nºpontos)
                1 , ... % maximum model order for model mation
                []); % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
            
            G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = ResultsSSGC.G;
            
        catch ME
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed')
            failures = failures + 1;
        end
    end
    
    % Significance test
    G_AllSubs_perConds_FDR.(cond_names{cond}) = nanmean(G_AllSubs_perConds_FDR.(cond_names{cond}),3);
    ResultsSSGC.pval = mvgc_pval(G_AllSubs_perConds_FDR.(cond_names{cond}),1,sizeblocks,nblocks,1,1,nnrois-2,''); % take careful note of arguments!
    SSsig_AllSubs_perConds_FDR.(cond_names{cond})  = significance(ResultsSSGC.pval,0.05,'FDR');
    
end

% Plot baseline & C1back & C2back
% figure(); clf;
% subplot(1,1,1);
% plot_pw(mean(SSsig_AllSubs_perConds_FDR.baseline,3));
% % % title('Pairwise-conditional SSGC - baseline');
% % % subplot(1,3,2);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.C1back);
% % % title('Pairwise-conditional SSGC - C1Back');
% % % subplot(1,3,3);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.C2back);
% % % title('Pairwise-conditional SSGC - C2Back');
% % 
% % % Plot baseline - C1back - C2back
% % % SSsig_AllSubs_perConds_FDR.combination = ... % SSsig_AllSubs_perConds_FDR.baseline + ...
% % %                                          SSsig_AllSubs_perConds_FDR.C1back + ...
% % %                                          SSsig_AllSubs_perConds_FDR.C2back;
% % % 
% % % figure(); clf;
% % % subplot(1,2,1);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.Baseline);
% % % title('Pairwise-conditional SSGC - Baseline');
% % % subplot(1,2,2);
% % % plot_pw(SSsig_AllSubs_perConds_FDR.Imaginacao);
% % % title('Pairwise-conditional SSGC - Imaginacao');
% % % 
% % % 
% % % % Plot baseline - C1back - C2back
% % % SSsig_AllSubs_perConds_FDR.combination = SSsig_AllSubs_perConds_FDR.Baseline + ...
% % %                                          SSsig_AllSubs_perConds_FDR.Imaginacao;
% % %                                      
% % % figure(); clf;
% % % plot_pw(SSsig_AllSubs_perConds_FDR.combination);
% % % title('Pairwise-conditional SSGC - Combination');




%% State-space - SSGC toolbox - each condition - pairwise
failures = 0;

for j = 1: size(C,1)
    fprintf('\n\n Combination %i\n\n',j)
    
    for cond = 1:n_cond
        fprintf('Performing SSGC for condition %s... \n',cond_names{cond});
        nblocks = size(DataMVGC.(cond_names{cond})(C(j,:),:,:),3);
        
        try
            ResultsSSGC.(['C' num2str(C(j,:),'%i%i')]) = ssgc_demo_function( blockDur(cond) , ... % number of observations per trial
                DataMVGC.(cond_names{cond})(C(j,:),:,:) , ... % data
                2 , ... % number of rois
                0.05 , ... % significance level for significance test
                1 , ... % maximum model order for model order estimation
                'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
                size( DataMVGC_II.(cond_names{cond}),3 ) , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                nblocks, ... % number of blocks per conditions
                cond_names{cond} , ... % condition names (or plot title in case of #conds = 1)
                0 , ... % Display model order plot
                0 , ... % Display results plot
                '', ... % 'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
                0); % checkConsistency (0 without check and 1 with check -mandatory nnrois<blocksize )
            
                G_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,1),C(j,2)) = ResultsSSGC.(['C' num2str(C(j,:),'%i%i')]).G(1,2);
                G_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,2),C(j,1)) = ResultsSSGC.(['C' num2str(C(j,:),'%i%i')]).G(2,1);

                SSsig_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,1),C(j,2)) = ResultsSSGC.(['C' num2str(C(j,:),'%i%i')]).sig(1,2);
                SSsig_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,2),C(j,1)) = ResultsSSGC.(['C' num2str(C(j,:),'%i%i')]).sig(2,1);

           
        catch ME
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCM failed')
            failures = failures + 1;
        end
    end
end

