%% Surrogate script for SSGC and MVGC

fprintf('');
tRev      = '2';                    % '1' - only the MVGC statistics, '2' - only the paired t-test, 'both' - both tests
graph1    = [];
FWE       = 'true';                 % Perform FWE - 'true' or 'false'
nhyp      = nodes*nodes - nodes;    % number of hypothesis
alpha     = 0.05;
mth       = 'FDR';                  % FWE method - 'bonf' Bonferroni, - 'FDR' FDR




%% =============================================================================================================
% ---- time-reversed ------------

if strcmp(surrogate,'timeRev')
    clear G_AllSubs_perConds_FDR_sur
    stop = 'false';        % 'false' condition defining to perform statistical test
    
    for blo = 1:nblocks
        fprintf('        - Block number: %d \n',blo);
        
        % ... condicao de que deu erro ao realizar o SSGC ...
        if isnan(Results.G_AllSubs_FDR{aux,i}(2,1,blo)) | Results.G_AllSubs_FDR{aux,i}(2,1,blo) == 0  
            stop = 'true';
            if simul
                z_scores.(cond_names{cond}) = zeros(nvars);
            end
            continue
        end
        
        % ... State-space for Time-Reversed Time-series ...
        % -> generate time-reversed data
        X_sur = fliplr(test.rshapData{aux,u}(:,:,blo)); 
        
        % -> model order
        if mordermax == 1
            moAIC = 1;
            moBIC = 1;
        else
            momax = mordermax;icregmode = 'OLS';[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit_v2(X_sur,momax,icregmode,false,[]);
        end
        pf = 2*moAIC; morder = moAIC; % when using AIC
        % pf = moBIC; morder = moBIC; % when using BIC
        
        % -> applied with SSGC
        try
            Results.G_AllSubs_FDR_sur{aux,i}(:,:,blo) = diag(diagonal);
            
            % a) building state-space
            [m_sur,A_sur,C_sur,K_sur,V_sur,z_sur,e_sur] = s4sid_CCA(X_sur,...   % time-series
                                                          pf,...
                                                          [],...               % model order
                                                          []);
            
            % b) buildng GC for surrugate
            Results.G_AllSubs_FDR_sur{aux,i}(:,:,blo)   = iss_PWGC(A_sur,C_sur,K_sur,V_sur); % valores de GC state-space
            
        catch
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed - Surrogate')
            %                     G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = zeros(var,var);
            failures_sur = failures_sur + 1;
        end

        
    end
    
    
    % ... Statistical Validation over the SSGC data by TRS ...
    if strcmp(stop,'false')
        if strcmp(tRev,'1')
            % ... with surrugate (time.reversed) ((û - u) / ome) ...
            avg_sur = mean(G_AllSubs_perConds_FDR_sur.(cond_names{cond}),3,'omitnan');
            std_sur = std(G_AllSubs_perConds_FDR_sur.(cond_names{cond}),0,3,'omitnan');
            G_corr = (mean( G_AllSubs_perConds_FDR.(cond_names{cond}),3,'omitnan') - avg_sur) ./ std_sur;
            morder = 1; nobs = size(DataMVGC_II.(cond_names{cond}),2); ntrials = size(DataMVGC_II.(cond_names{cond}),3);
            pval_sur = mvgc_pval(G_corr,morder,nobs,ntrials,1,1,nvars-2,'chi2'); % take careful note of arguments!
            sig_sur  = significance(pval_sur,alpha,'FDR');
            if plots == 1
                figure(cond+10)
                title('Z-scores for GCreal vs GCsurrogate tRev. (subtraction)');
                plot_pw(sig_sur);
            end
        end
        
        if strcmp(tRev,'2')
            %... paired t-test for all channels ...
            [h,p,ci,stats] = ttest(G_AllSubs_perConds_FDR.(cond_names{cond}),G_AllSubs_perConds_FDR_sur.(cond_names{cond}),0.05,'both',3);
            
            if strcmp(mth,'FDR') % Family-wise Errors correction - FDR
                [h,alpha1] = fdr_bh_ECt(p,0.05,true);
                stats.tstat = stats.tstat .* h;
            end
            
            %... Convert z-scores ...
            mu     = 0;
            sigma  = 1;
            prob   = tcdf(stats.tstat,stats.df(1,2));     % CDF of the Students t-distribution with df degrees of freedom
            z_scor = norminv(prob,mu,sigma);              % Inverse cdf of the univariate standard normal distribution
            [blo,j]       = find(z_scor == Inf);
            infinit     = find(z_scor == Inf);
            
            
            if strcmp(mth,'bonf') % Family-wise Errors correction - Bonferroni
                alpha1 = alpha/nhyp;
                Auxiliar calculus
                z = @(p) -sqrt(2) * erfcinv(p*2);
                cval = z(1-alpha1); % critical value - z-value corrected for FWE
                z_scor(find(abs(z_scor)<cval)) = 0;
            end
            
            %... Plot with FWE correction ...
            if plots == 1
                figure()
                plot_pw(z_scor);
                str = ['Z-scores - GCreal vs GCsurr - tRev - paired t-test - Cond: ',(cond_names{cond})];
                title(str);
                caxis([min(min(z_scor)) max(max(z_scor))]);
                colormap(color_map);
            end
            if simul
                z_scores.(cond_names{cond}) = z_scor;
            end
        end
    end
end





