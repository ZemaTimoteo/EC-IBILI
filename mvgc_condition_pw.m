%%% ========= CONNECTIVITY ANALYSIS SCRIPT ========= %%%
%%% ============= PC and GCM ============= %%%
%%% ======================================== %%%

% Remember to have NeuroElf installed

%% Load Workspace
% % % clear;
% % % load('workspace.mat');
% % % i=1


%% Configuration
close all;

addpath('mvgc_new')

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


%% GCM

% for all conditions

cond = 1;
% sig_AllCond = zeros(n_rois,n_rois);
% F_AllCond = zeros(n_rois,n_rois,n_sub);
for c = 1:3
    DataMVGC_II.(cond_names{c}) = DataMVGC.(cond_names{c})(13:end, :,:);
end

voiIdx_II = 1:size(DataMVGC.(cond_names{cond}), 1);
voiNames_II = voiNames(13:end);
voiIdx_II(1:12) = [];

nnrois = size(DataMVGC_II.(cond_names{1}),1);

counter = 1;


for i = 1:nnrois
    
    %         [ksstat,cval] = mvgc_kpss(DataMVGC_II.(cond_names{cond})(:,:,i),0.05);
    
    if ( find(DataMVGC_II.(cond_names{cond})(i,:,:)==0) )
        idxsRemoval(counter) = i;
        counter = counter + 1;

    end
        
end


% remove idxs
for c = 1:n_cond
        DataMVGC_II.(cond_names{c})(idxsRemoval,:,:)=[];

end
         voiNames_II(idxsRemoval)=[];
        voiIdx_II(idxsRemoval) = [];

nnrois = size(DataMVGC_II.(cond_names{1}),1);

C = combnk(1:nnrois,2); % Possible Combinations


%% MVGC for each condition of all subjects

failures = 0;

for j = 1: size(C,1)
    fprintf('Combination %i\n',j)
    
    for cond = 1:n_cond
        fprintf('Performing MVGC for condition %s... \n',cond_names{cond});
        
        
        try
            ResultsMVGC.(['C' num2str(C(j,:),'%i%i')]) = mvgc_demo_function( blockDur(cond) , ... % number of observations per trial
                1/(TR/1000) , ... % sample rate (Hz)
                DataMVGC.(cond_names{cond})(C(j,:),:,:) , ... % data
                2 , ... % number of rois
                0.001 , ... % significance level for significance test
                trialIndex(cond)-1 , ... % number of trials
                1 , ... % maximum model order for model order estimation
                'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
                size( DataMVGC_II.(cond_names{cond}),3 ) , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                cond_names{cond} , ... % condition names (or plot title in case of #conds = 1)
                0 , ... % Display model order plot
                0 , ... % Display results plot
                0); % checkConsistency (0 without check and 1 with check -mandatory nnrois<blocksize )
            
            
                F_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,1),C(j,2)) = ResultsMVGC.(['C' num2str(C(j,:),'%i%i')]).F(1,2);
                F_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,2),C(j,1)) = ResultsMVGC.(['C' num2str(C(j,:),'%i%i')]).F(2,1);

                sig_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,1),C(j,2)) = ResultsMVGC.(['C' num2str(C(j,:),'%i%i')]).sig(1,2);
                sig_AllSubs_perConds_pw_FDR.(cond_names{cond})(C(j,2),C(j,1)) = ResultsMVGC.(['C' num2str(C(j,:),'%i%i')]).sig(2,1);

        
        catch ME
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCM failed')
            failures = failures + 1;
        end
    end
end

%%


   figure
        plot_pw_voiNames(sig_AllSubs_AllConds_pw_FDR, voiNames_II);
        colorbar
        title('All Subs All Conds pw FDR');
        

%%
for i = 1:n_cond
   figure
        plot_pw_voiNames(sig_AllSubs_perConds_pw_FDR.(cond_names{i}), voiNames_II);
        colorbar
        title(['All Subs per Conds pw FDR: ' cond_names{i}]);
        
end




%% Surrugate
fprintf('--- Performing Surrugate... --- \n');

% --- time-reversed strategy ---
fprintf('   -> Time-reversed strategy \n');
for i=1:size(X,3)
    X_sur(:,:,i) = fliplr(X(:,:,i)); % generate time-reversed data
    
    % morder selection - past/future horizons for canonical correlations (pf)
    momax = 10;
    icregmode = 'OLS';
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X(:,:,i),momax,icregmode);
    pf(i)   = 2*moAIC;
    morder  = 3;%moBIC;
    
    % building state-space
    [m_sur,A_sur(:,:,i),C_sur(:,:,i),K_sur(:,:,i),V_sur(:,:,i),z_sur,e_sur] ...
                                            = s4sid_CCA(X_sur(:,:,i),...   % time-series
                                                        pf(i),...   
                                                        []);               % model order                                                                                            
                                                    
    % buildng GC for surrugate
    G_sur(:,:,i) = iss_PWGC(A_sur(:,:,i),C_sur(:,:,i),K_sur(:,:,i),V_sur(:,:,i)); % valores de GC steady-space
    fprintf('\n\n');
end



%% Plots
fprintf('--- Plots... --- \n\n');

% Plot original
figure(); clf;
subplot(1,3,1);
plot_pw(mean(G,3,'omitnan'));
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])

% Plot original & surrogate (time.reversed) ((û - u) / ome)
figure(); clf;
subplot(1,3,1);
plot_pw(F_corr);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval_sur);
title('p-values');
subplot(1,3,3);
plot_pw(sig_sur);
title(['Significant at p = ' num2str(alpha)])
  