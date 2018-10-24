%% Load data
load('..\ResultsMatrix\5var_HighCorr_SSGC_timeRev_SNR200a_ChannelGap10-200_Morder1.mat')

%%

n_cond= 3; 
cond_names = {'simulate05','simulate250','simulateOriginal'};

% variable *data*
f1_results = data.results.Fone;

% final matrix
results05 = [];
results250 = [];
resultsOrig = [];

  
for p = 1 : 9 % num of points
    for b = 1 : 9 % num of blocks
        results05(p,b) = f1_results{p,b}.(cond_names{1});
        results250(p,b) = f1_results{p,b}.(cond_names{2});
        resultsOrig(p,b) = f1_results{p,b}.(cond_names{3});
    end
end






    