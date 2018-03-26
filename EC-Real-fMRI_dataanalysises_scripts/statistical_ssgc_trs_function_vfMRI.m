%% =============================================================================================================

% ---- Statistical Validation over the SSGC data by TRS - using t-test ------------

%% ... Initialization of variables ...

tRev      = '2';                    % '1' - only the MVGC statistics, '2' - only the paired t-test, 'both' - both tests
graph1    = [];
FWE       = 'true';                 % Perform FWE - 'true' or 'false'
nhyp      = nodes*nodes - nodes;    % number of hypothesis
alpha     = 0.05;
mth       = 'FDR';                  % FWE method - 'bonf' Bonferroni, - 'FDR' FDR
plots     = 1;
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
    [h,p,ci,stats] = ttest(Results.G_AllSubs_FDR{aux,v},Results.G_AllSubs_FDR_sur{aux,v},0.05,'both',3);
    
    % ... Family-wise Errors correction ...
    if strcmp(mth,'FDR') % FDR
        [h,alpha1] = fdr_bh_ECt(p,0.05,true);
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
        alpha1 = alpha/nhyp;
        Auxiliar calculus
        z = @(p) -sqrt(2) * erfcinv(p*2);
        cval = z(1-alpha1); % critical value - z-value corrected for FWE
        z_scor(find(abs(z_scor)<cval)) = 0;
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
    Results.z_scores{aux,u} = z_scor;
    
end
