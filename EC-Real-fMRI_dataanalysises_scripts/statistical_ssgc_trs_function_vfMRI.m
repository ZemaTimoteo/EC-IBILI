%% =============================================================================================================

% ---- Statistical Validation over the SSGC data by TRS - using t-test ------------

%% ... Initialization of variables ...

tRev      = '2';                    % '1' - only the MVGC statistics, '2' - only the paired t-test, 'both' - both tests
graph1    = [];
FWE       = 'true';                 % Perform FWE - 'true' or 'false'
nhyp      = nodes*nodes - nodes;    % number of hypothesis
alpha     = 0.05;
mth       = 'none';                  % FWE method - 'bonf' Bonferroni, - 'FDR' FDR, - 'none' no FWE correction
plots     = 0;
%% ... Statistical testing ...

if strcmp(tRev,'1')
    % ... with surrugate (time.reversed) ((û - u) / ome) ...
    avg_sur = mean(Results.G_AllSubs_FDR_sur{aux,v},3,'omitnan');
    std_sur = std(Results.G_AllSubs_FDR_sur{aux,v},0,3,'omitnan');
    G_corr = (mean(Results.G_AllSubs_FDR{aux,v},3,'omitnan') - avg_sur) ./ std_sur;
    nobs = size(test.rshapData{aux,v},2); ntrials = size(test.rshapData{aux,v},3);
    pval_sur = mvgc_pval(G_corr,morder,nobs,ntrials,1,1,nvars-2,'chi2'); % take careful note of arguments!
    sig_sur  = significance(pval_sur,alpha,'FDR');
    if plots == 1
        figure(cond+10)
        title('Z-scores for GCreal vs GCsurrogate tRev. (subtraction)');
        plot_pw(sig_sur);
    end
end

if strcmp(tRev,'2')
    % ... Paired t-test for all channels ...    
    cat_G_AllSubs     = cat(3,Results.G_AllSubs_FDR{aux,v});
    cat_G_AllSubs_sur = cat(3,Results.G_AllSubs_FDR_sur{aux,v});    
    [h,p,ci,stats]    = ttest(cat_G_AllSubs,cat_G_AllSubs_sur,alpha,'both',3);
    
    % ... Family-wise Errors correction ...
    if strcmp(mth,'FDR') % FDR
        [h,alpha1] = fdr_bh_ECt(p,alpha,true);
        stats.tstat = stats.tstat .* h;
        
        % Convert z-scores
        mu       = 0;
        sigma    = 1;
        prob     = tcdf(stats.tstat,stats.df(1,2));     % CDF of the Students t-distribution with df degrees of freedom
        z_scor   = norminv(prob,mu,sigma);              % Inverse cdf of the univariate standard normal distribution
        [blo,j]  = find(z_scor == Inf);
        infinit  = find(z_scor == Inf);
    end
       
    if strcmp(mth,'bonf') % Bonferroni
        % Convert z-scores
        mu       = 0;
        sigma    = 1;
        prob     = tcdf(stats.tstat,stats.df(1,2));     % CDF of the Students t-distribution with df degrees of freedom
        z_scor   = norminv(prob,mu,sigma);              % Inverse cdf of the univariate standard normal distribution
        [blo,j]  = find(z_scor == Inf);
        infinit  = find(z_scor == Inf);
        
        alpha1 = alpha/nhyp;
        Auxiliar calculus
        z = @(p) -sqrt(2) * erfcinv(p*2);
        cval = z(1-alpha1); % critical value - z-value corrected for FWE
        z_scor(find(abs(z_scor)<cval)) = 0;
    end
    
    if strcmp(mth,'none') % no FWE correction
        FWE      = 'false';
        % Convert z-scores
        mu       = 0;
        sigma    = 1;
        prob     = tcdf(stats.tstat,stats.df(1,2));     % CDF of the Students t-distribution with df degrees of freedom
        z_scor   = norminv(prob,mu,sigma);              % Inverse cdf of the univariate standard normal distribution
%         [blo,j]  = find(z_scor == Inf);
%         infinit  = find(z_scor == Inf);
        
        z = @(p) -sqrt(2) * erfcinv(p*2);
        cval = z(1-alpha); % critical value - z-value corrected for FWE
%         z_scor(find(abs(z_scor)<cval)) = 0; % elimina os valores abaixo do limiar Critical Value
        HighCritVal = zeros(nodes);
        HighCritVal(find(abs(z_scor)>cval)) = 1;
        
        clear mu sigma prob cval
    end
    
    % ... Plot with FWE correction ...
    if plots == 1
        figure()
        plot_pw(z_scor);
        str = ['Z-scores - GCreal vs GCsurr - tRev - paired t-test'];
        title(str);
        caxis([min(min(z_scor)) max(max(z_scor))]);
        colormap(color_map);
    end
    
    % ... save z-scores ...
    Results.z_scores{aux,v}    = abs(z_scor);
    Results.HighCritVal{aux,v} = HighCritVal;
end

test.stat.tRev      = tRev;                 % '1' - only the MVGC statistics, '2' - only the paired t-test, 'both' - both tests
test.stat.graph1    = graph1;
test.stat.FWE       = FWE;                  % Perform FWE - 'true' or 'false'
test.stat.nhyp      = nhyp;                 % number of hypothesis
test.stat.alpha     = alpha;
test.stat.mth       = mth;                  % FWE method - 'bonf' Bonferroni, - 'FDR' FDR, - 'none'

clear tRev graph1 FWE nhyp alpha mth str stop plots
