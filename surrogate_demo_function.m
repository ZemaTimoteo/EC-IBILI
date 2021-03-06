%% Surrogate script for SSGC and MVGC

fprintf('    - Performing Surrugate... \n');
tRev      = '2';                    % '1' - only the MVGC statistics, '2' - only the paired t-test, 'both' - both tests
graph1    = [];
FWE       = 'true';                 % Perform FWE - 'true' or 'false'
nhyp      = nvars*nvars - nvars;    % number of hypothesis
alpha     = 0.05;
mth       = 'FDR';                  % FWE method - 'bonf' Bonferroni, - 'FDR' FDR




%% =============================================================================================================
% ---- time-reversed ------------

if strcmp(surrogate,'timeRev')
    clear G_AllSubs_perConds_FDR_sur
    for cond = 1:n_cond
        stop      = 'false';        % 'false' to perform statistical test
        
        for i = 1:nblocks
            fprintf('         -Condition %d',cond);
            fprintf(' - bloco number: %d\n',i);
            if isnan(G_AllSubs_perConds_FDR.(cond_names{cond})(2,1,i)) | G_AllSubs_perConds_FDR.(cond_names{cond})(2,1,i) == 0  % condicao de que deu erro ao realizar o SSGC
                stop = 'true';
                if simul
                    z_scores.(cond_names{cond}) = zeros(nvars);
                end
                continue
            end
            % State-space for t-rev.
            X_sur = fliplr(DataMVGC_II.(cond_names{cond})(:,:,i)); % generate time-reversed data
            momax = mordermax;icregmode = 'OLS';[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit_v2(X_sur,momax,icregmode,false,[]);
            pf = 2*moAIC; morder = moAIC; % when using AIC
%             pf = moBIC; morder = moBIC; % when using BIC
            
            
            
            
            if test == 1 % ----- MVGC trials --------
                try
                    G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,nblocks) = nan(nnrois);
                    ResultsMVGC = mvgc_demo_function( blockDur(cond) , ... % number of observations per trial
                        fs(cond) , ...  % sample rate (Hz)
                        X_sur  , ... % data
                        nnrois , ... % number of rois
                        0.05 , ... % significance level for significance test
                        1 , ... % number of trials
                        mordermax , ... % maximum model order for model order estimation
                        'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
                        'AIC', ... % model order method
                        [] , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                        'All Conditions' , ... % condition names (or plot title in case of #conds = 1)
                        0 , ... % Display model order plot
                        0); % Display results plot
                    G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = ResultsMVGC.F;
                catch
                    G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = zeros(5);
                    failures = failures + 1;
                end
            end
            
            
            
            
            if test == 3 % ----- SSGC trials ------
                try
                    diagonal  = nan(var,1)';
                    G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = diag(diagonal);
                    % building state-space
                    [m_sur,A_sur,C_sur,K_sur,V_sur,z_sur,e_sur] = s4sid_CCA(X_sur,...   % time-series
                        pf,...
                        [],...               % model order
                        []);
                    % buildng GC for surrugate
                    G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = iss_PWGC(A_sur,C_sur,K_sur,V_sur); % valores de GC steady-space
                catch
                    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed - Surrogate')
%                     G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = zeros(var,var);
                    failures = failures + 1;                   
                end
            end
            
            
            
            if test == 4 % ----- Coherence trials ------
                G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,nblocks) = nan(nnrois);
                try
                    Coh = zeros(nnrois);
                    for comb = 1:size(c,1)
                        [tempCoh,tempF] = mscohere(X_sur(c(comb,1),:)',X_sur(c(comb,2),:)',[],[],[],fs(cond));
                        Coh(c(comb,2),c(comb,1)) = mean(tempCoh(tempF<0.15));
                    end
                catch ME
                    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Coherence failed')
                    failures = failures + 1;
                end
                
                G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = atanh(Coh);
                
                label_ec = [0 0 0 0 0; ...
                            2 0 0 0 0; ...
                            2 0 0 0 0; ...
                            2 0 0 0 0; ...
                            0 0 0 2 0];
            end
            
        end
        
        
        
        if strcmp(stop,'false')
            if strcmp(tRev,'1')
                % ... with surrugate (time.reversed) ((� - u) / ome) ...
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
                [i,j]       = find(z_scor == Inf);
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
end



%% =============================================================================================================
% ---- permutation --------------
if strcmp(surrogate,'perm')
    
    sizeblocks = 300; blockDur(cond) = sizeblocks;
    U = DataMVGC_II.(cond_names{cond})(:,1:sizeblocks);  %  multi-trial time series data
    [n,m,nblocks] = size(U);  %  multi-trial time series data
    
    p =  1;       %  model order (number of lags)
    bsize   = p;   % permutation block size (default: use model order)
    nsamps  = 500; %   number of permutations
    npermt  = floor(m/bsize); % number of blocks
    
    if npermt*bsize ~= m
        oldm = m;
        m = npermt*bsize;
        U = U(:,1:m,:);
    end
    
    for cond = 1:1
        G_AllSubs_perConds_FDR_perm.(cond_names{cond}) = nan(n,n,nsamps);
        SSsig_AllSubs_perConds_FDR_perm.(cond_names{cond}) = nan(n,n,nsamps);
        permutation = 1;
        fprintf('Condition %d',cond,' Number of permutation: %d\n',nsamps);
        for i=1:nsamps
            fprintf('%d / %d - ',permutation,nsamps);
            Uperm = U(:,randperm(npermt));
            try
                ResultsSSGC_perm = ssgc_demo_function( sizeblocks , ... % number of observations per trial
                    Uperm , ... % data (altera�ao do n�pontos)
                    nnrois , ... % number of rois
                    0.05 , ... % significance level for significance test
                    p , ... % maximum model order for model order estimation
                    'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
                    [] , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                    nblocks, ... % number of blocks per conditions
                    cond_names{cond} , ... % condition names (or plot title in case of #conds = 1)
                    0 , ... % Display model order plot
                    0 , ... % Display results plot
                    '', ... % 'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
                    0); % checkConsistency (0 without check and 1 with check -mandatory nnrois<blocksize )
                
                G_AllSubs_perConds_FDR_perm.(cond_names{cond})(:,:,i) = ResultsSSGC_perm.G;
                SSsig_AllSubs_perConds_FDR_perm.(cond_names{cond})(:,:,i) = ResultsSSGC_perm.sig;
                distrib(i) = max(max(G_AllSubs_perConds_FDR_perm.(cond_names{cond})(:,:,i)));
                permutation = 1 + permutation;
            catch ME
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GCM failed')
                failures = failures + 1;
            end
        end
    end
    
    % Plotting null distribution
    % ... threshold ...
    ite           = 0;
    thresh        = 0;
    distrib = sort(distrib);
    while thresh <0.95
        ite    = ite+1;    % Family-wise Errors correction
        
        thresh = sum(distrib(1:ite))/sum(distrib);
    end
    thresh = distrib(ite-1);
    % ... plot ...
    figure()
    histfit(distrib);
    line([thresh thresh],[0 20],'Color','g'); % percentil of 0.95%
    % ... applying FWE to data ...
    if exist('ResultsSSGC') == 0 
        ResultsSSGC.G = ResultsMVGC.F;
    end
    ResultsSSGC.G(find(ResultsSSGC.G<=thresh & ResultsSSGC.G>=-thresh)) = 0;
    figure()
    plot_pw(ResultsSSGC.G);
    title('GCI values: permutation statistics');
    caxis([min(min(ResultsSSGC.G)) max(max(ResultsSSGC.G))])
    %     colormap(c_map)
    
end




%% =============================================================================================================
% ---- permutation with blocks --------------
if strcmp(surrogate,'permBlock')
    for cond = 1:n_cond
        %  multi-trial time series data
        U = DataMVGC_II.(cond_names{cond});
        [n,blockDur(cond),nblocks] = size(U);
        sizeblocks = blockDur(cond);
        % permutation intilization
        p =  1;        % model order (number of lags)
        bsize   = p;   % permutation block size (default: use model order)
        nsamps  = floor(1000/nblocks); % number of permutations
        npermt  = floor(sizeblocks/bsize); % number of blocks
        if npermt*bsize ~= sizeblocks
            oldm = sizeblocks;
            blockDur(cond) = npermt*bsize;
            U = U(:,1:blockDur(cond),:);
        end
        
        
        fprintf('Condition %d - ',cond)
        fprintf(' Number of permutation: %d\n',nsamps);
        value = 1;
        for blo=1:nblocks
            fprintf('... %d BLOCK ...\n',blo);
            bloco{blo} = ['bloco' num2str(blo)];
            G_AllSubs_perConds_FDR_perm.(cond_names{cond}).(bloco{blo}) = nan(n,n,nsamps);
            SSsig_AllSubs_perConds_FDR_perm.(cond_names{cond}).(bloco{blo}) = nan(n,n,nsamps);
            for i=1:nsamps
                fprintf('Cond %d, Block number: %d -> Permt %d / %d - ',cond,blo,i,nsamps);
                Uperm = U(:,randperm(npermt),blo);
                try
                    ResultsSSGC_perm = ssgc_demo_function_v2( blo, ...  % number of block
                        Uperm(:,:) , ... % data (altera�ao do n�pontos)
                        p , ... % maximum model order for model mation
                        []); % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                    
                    G_AllSubs_perConds_FDR_perm.(cond_names{cond}).(bloco{blo})(:,:,i) = ResultsSSGC_perm.G;
                    distrib(value,cond) = max(max(ResultsSSGC_perm.G));
                    value = value+1;
                catch ME
                    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed\n\n\n\n\n\n')
                    failures = failures + 1;
                end
                
            end
        end
    end
    for cond = 1:n_cond
        
        % Plotting null distribution
        % ... threshold ...
        ite     = 0;
        thresh  = 0;
        alfa    = 0.05;
        
        if FWE % Family-wise error correction (Bonferroni)
            alfa = alfa/nhyp;
        end
        if  exist('distrib')
            dis(:,cond) = sort(distrib(:,cond));
            while thresh <1-alfa
                ite    = ite+1;    % ----- Family-wise Errors correction ---------------------------
                
                thresh = sum(dis(1:ite,1))/(sum(dis(:,1)));
            end
            try
                thresh = dis(ite-1,cond);
            catch
                Gcond = zeros(size(G_AllSubs_perConds_FDR.(cond_names{cond}),1),size(G_AllSubs_perConds_FDR.(cond_names{cond}),1));
                z_scores.(cond_names{cond}) = Gcond;         % save variable
            end
            % ... plot ...
            if graph1
                figure(55)
                histfit(dis(:,cond));
                line([thresh thresh],[0 20],'Color','g'); % percentil of 0.95%
            end
            % ... applying FWE to data ...
            Gcond = mean(G_AllSubs_perConds_FDR.(cond_names{cond}),3);
            Gcond(find(Gcond<=thresh & Gcond>=-thresh)) = 0;
            if graph1
                figure()
                plot_pw(Gcond);
                title(['GCI values: permutation statistics for ' (cond_names{cond})]);
                caxis([min(min(Gcond)) max(max(Gcond))])
                colormap(c_map)
            end
            % save variable
        else
            Gcond = zeros(size(G_AllSubs_perConds_FDR.(cond_names{cond}),1),size(G_AllSubs_perConds_FDR.(cond_names{cond}),1));
        end
        
        if simul
            z_scores.(cond_names{cond}) = Gcond;
        end
        
    end
end
