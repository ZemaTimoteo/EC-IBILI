% ===== Using Dual-Extended kalman filter ========

%% load
clear all, clc
fprintf('--- Load... --- \n\n');
addpath(genpath('C:\Users\admin\Documents\Tiago\IBILI\Code\Toolboxes'))
load('C:\Users\admin\Documents\Tiago\IBILI\Data\testes\test_ti_simul\testData_noVC.mat')

%% pre-processing 
fprintf('--- Pre-processing... --- \n');
% --- struct ---
[CH,m,nt]        = size(X);
inp_model.methd_ord = 'BIC'; % Method of selecting the model order
inp_model.momax  = 10;
inp_model.srate  = 1000;
acmaxlags        = 500;
timeStim         = 1;

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
graph     = 'true'; % do the plot ('true') it does not ('false')


% --- filtering ---
fprintf('   -> Filtering \n');
figure()
plot (X(2,:,1))
hold on,
X = reshape(X,CH,m*nt);
X = filtra(X,1,250,inp_model.srate);
% X = filtra_notch(X,50,inp_model.srate);
X = reshape(X,CH,m,nt);
plot (X(2,:,1))
hold off,
% build data
% data_plot = iddata(X(:,:,500)',[],1/1000);


% --- data ---
inp_model.data   = X(:,:,1)';
[m,CH,nt]        = size(inp_model.data);


% --- model order selection ---
inp_model.order = 2;%morder_ECt( X , inp_model );


% --- downsampling ---
fprintf('   -> Downsampling \n');
new_fs   = 250;
n        = inp_model.srate/new_fs;   % choose a point each n-th datapoints
% data_seg = inp_model.data;%reshape (inp_model.data,CH,m*nt);
inp_model.data  = downsample(inp_model.data,n);
inp_model.srate = new_fs;
% inp_model.data  = data_seg; %reshape (data_seg,inp_model.srate*timeStim,CH);
    

% --- GC pre-proc. ---
% fprintf('   -> Non-stati. removal \n');
% X = mvdetrend(X);      % detrend of data
% X = demean(X,'true');  % remove of ensemble mean and ensemble varia

fprintf('\n');
%% Dual Kalman Filter - parameter estimation
fprintf('--- Dual Kalman Filter... --- \n\n');
[A_full,SIG_full] = DEKF3(inp_model);
[CH,L,T] = size(A_full);
A_full = reshape(A_full,CH,CH,L/CH,T);

%% consistency test - Seth based
fprintf('--- Consistency... --- \n');
ite = size(X,1)+1;
sum = 0;
for ii = 1:ite:size(SIG_full,3)
    if ii+ite-1<=size(SIG_full,3)
        aux = SIG_full(:,:,ii:ii+ite-1);
        cons = consistency(X(:,ii:ii+ite-1,1),mean(aux,3));
        
        figure(5), plot(ii,cons,'-ro');  % plot
        hold on
        grid on
        title('Consistency statistics');
        ylabel('Consistency %');
        xlabel('Evaluation points (each ii = 50 datapoints)');
        
        if cons>= 0.8
            sum = 1+sum;
        end
    end
end
    fprintf('     -> %d points with consistency bigger than 0.8 \n\n',sum);

%% calculate GC time-domain
fprintf('--- GC (Time-domain)... --- \n\n');
lag_avg = 50;
for t = inp_model.order : T % Number of time points
    F(:,:,t) = nan(size(A_full,1));
    fprintf('Point number %d ... ',t);
    % Generate Autocov coef. 
    [G,info] = var_to_autocov(A_full(:,:,:,t),SIG_full(:,:,t),acmaxlags);
    
%     if graph
%         compare(data,sys)
%     end
    
    try
        var_info(info,true); % report results (and bail out on error)
    catch
        continue
    end
	% Generate GC values  
    n = size(G,1);
   
    % full regression
    owstate = warn_supp;
    [~,SIG] = autocov_to_var(G);
    warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
    
    if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
    LSIG = log(diag(SIG));
    
    for j = 1:n
        % reduced regression
        jo = [1:j-1 j+1:n]; % omit j
        owstate = warn_supp;
        [~,SIGj] = autocov_to_var(G(jo,jo,:));
        warn_test(owstate,     sprintf('in reduced regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
        if warn_if(isbad(SIGj),sprintf('in reduced regression for target node %d - regression failed',j)), continue; end % show-stopper!
        LSIGj = log(diag(SIGj));
        for ii=1:n-1;
            i = jo(ii);
            F(i,j,t) = LSIGj(ii)-LSIG(i);
        end
    end
        assert(~isbad(F(:,:,t),false),'GC calculation failed');
end

%% Surrugate
fprintf('--- Surrogates... --- \n');

% --- time-reversed strategy ---
% fprintf('   -> Time-reversed strategy \n');
% for i=1:size(inp_model.order,3)
%     nsamps = 1;
%     sur.data(nsamps,:,:,i) = fliplr(inp_model.data(:,:,i)); % generate time-reversed data
% end
% 

% --- permutation strategy ---
fprintf('   -> Permutation strategy \n');
bsize     = [];     % permutation test block size (uses model order if =[])
nsamps    = 100;    % number of permutations for permutation test
regmode   = 'OLS';  
acdectol  = [];
U = inp_model.data'; [n,m,N] = size(U);
nblocks   = floor(m/inp_model.order); % number of blocks

if nblocks*inp_model.order ~= m
    oldm = m;
    m = nblocks*inp_model.order;
    U = U(:,1:m,:);
    fprintf(2,'WARNING: truncating sequence length by %d observations\n',oldm-m);
end
for j = 1:n
    jo  = [1:j-1 j+1:n]; % omit j
%     sur.data = U;
    UUj = reshape(U(j,:,:),1,inp_model.order,nblocks,N); % stack blocks
    
    for s = 1:nsamps
        fprintf('PWCGC from node %d: permutation test sample %d of %d\n',j,s,nsamps);
        for r = 1:N
            sur.data(j,:,r,s) = reshape(UUj(:,:,randperm(nblocks),r),1,m); % permute blocks and unstack
        end
    end
end


% --- implementation ---
len = size(sur.data,2);
surr.order = inp_model.order;
surF       = nan(nsamps,n,n,len);
sur_ff     = nan(nsamps,n,n);
sur_std    = nan(nsamps,n,n);

for v= 1:nsamps
    fprintf('Perm. GC number %d ... \n',v);
    surr.data = sur.data(:,:,:,v)';
    [surA_full,surSIG_full] = DEKF3(surr);      % Dual Kalman Filter - parameter estimation
    
    % [CH,L,T] = size(surA_full);
    surA_full = reshape(surA_full,CH,CH,L/CH,len);
    
    for t = surr.order : len % Number of time points
        surF_aux(:,:,t) = nan(size(surA_full,1));
        % Generate Autocov coef.
        [surG,info] = var_to_autocov(surA_full(:,:,:,t),surSIG_full(:,:,t),acmaxlags);
        
        %     if graph
        %         compare(data,sys)
        %     end
        
        try
            var_info(info,true); % report results (and bail out on error)
        catch
            continue
        end
        % Generate GC values
        n = size(surG,1);
        
        % full regression
        owstate = warn_supp;
        [~,surSIG] = autocov_to_var(surG);
        warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
        
        if warn_if(isbad(surSIG),'in full regression - regression failed'), return; end % show-stopper!
        surLSIG = log(diag(surSIG));
        
        for j = 1:n
            % reduced regression
            jo = [1:j-1 j+1:n]; % omit j
            owstate = warn_supp;
            [~,surSIGj] = autocov_to_var(surG(jo,jo,:));
            warn_test(owstate,     sprintf('in reduced regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
            if warn_if(isbad(surSIGj),sprintf('in reduced regression for target node %d - regression failed',j)), continue; end % show-stopper!
            surLSIGj = log(diag(surSIGj));
            for ii=1:n-1;
                i = jo(ii);
                surF_aux(i,j,t) = surLSIGj(ii)-surLSIG(i);
            end
        end
        assert(~isbad(surF_aux(:,:,t),false),'GC calculation failed');
    end
    surF(v,:,:,:)  = surF_aux;
    sur_ff(v,:,:)  = mean(surF_aux(:,:,1:end),3,'omitnan');
    sur_std(v,:,:) = std(surF_aux(:,:,1:end),0,3,'omitnan');
end

fprintf('\n');
%% Significance test
fprintf('--- Significance test... --- \n');
ff = mean(F(:,:,1:end),3,'omitnan');
ff(find(ff==0))= nan;
sur_ff(find(sur_ff==0))= nan;
ntrials = 1;
morder  = inp_model.order;
nobs    = 1000;%lag_avg;
nvars   = CH;


% --- without surrugate ---
% pval = mvgc_pval(ff,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
% sig  = significance(pval,alpha,mhtc);


% --- with surrugate (time.reversed) ((û - u) / ome) ---
% fprintf('   -> Time-reversed strategy \n');

% fff = (ff - sur_ff) ./ sur_std;
% pval = mvgc_pval(fff,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
% sig  = significance(pval,alpha,mhtc);


% --- with surrugate (time.reversed) t-test  ---
% [datachan,~,ntrials] = size(F);
% 
% % Building antisymmetric matrix
% % GCI_dat_Net=zeros(datachan,datachan,ntrials);
% % GCI_sur_Net=zeros(datachan,datachan,ntrials);
% % 
% % for i=1:datachan
% %     for j=1:datachan
% %         if j~=i
% %             GCI_dat_Net(i,j,:) = F(i,j,:) - F(j,i,:);
% %             GCI_sur_Net(i,j,:) = surF(i,j,:) - surF(j,i,:);
% %         end
% %     end
% % end
% % 
% % figure(11)
% % imagesc(GCI_dat_Net(:,:,2));
% % colormap(gray);
% % figure(12),
% % imagesc(GCI_sur_Net(:,:,2));
% % colormap(gray);
% 
% % paired t-test for all channels
% [h,p,ci,stats] = ttest(F,surF,0.05,'both',3);
% % ... plot p-values ...
%     figure(13)
%     p_aux = p;
%     p_aux(find(p>0.05)) = 1;
%     plot_pw(p_aux);
%     figure(14)
%     plot_pw(p);
%     title('P-values for GCreal vs GCsurrogate');
% 
% % Convert z-scores
% mu     = 0;
% sigma  = 1;
% prob   = tcdf(stats.tstat,stats.df(1,1));  % CDF of the Students t-distribution with df degrees of freedom
% z_scor = norminv(prob,mu,sigma);         % Inverse cdf of the univariate standard normal distribution
% % z_scor = icdf('Normal',prob,mu,sigma);     % Inverse cdf of the univariate standard normal distribution
% [i,j]       = find(z_scor == Inf);
% infinit     = find(z_scor == Inf);
% if size(i,2)~=0 && size(j,2)~=0
%     diagonal    = diag(z_scor(j,i));
%     z_scor(infinit) = diagonal*(-1);
% end
% % Family-wise Errors correction
% nhyp   = (datachan*datachan)/2 - datachan;
% alpha  = 0.05;
% alpha1 = alpha/nhyp; % Bonferroni method - alpha value
% % ... Auxiliar calculus ...df
% z = @(p) -sqrt(2) * erfcinv(p*2);
% 
% cval = z(1-alpha1) % critical value - z-value corrected for FWE
% z_scor(find(abs(z_scor)<cval)) = 0;
% 
% % Plot with FWE correction
% figure(15)
% plot_pw(z_scor);
% title('Z-scores for GCreal vs GCsurrogate');
% caxis([min(min(z_scor)) max(max(z_scor))])
% colormap(color_map)


% --- with surrugate (permutation) ---
fprintf('   -> Permutation strategy \n');
% Permutation test significance test (adjusting for multiple hypotheses).

pval_p = empirical_pval(ff,sur_ff);
sig_p  = significance(pval_p,alpha,mhtc);

% Unbiased GC values
avg_sur_ff = mean(nanmean(nanmean(sur_ff)));
UNff = ff-avg_sur_ff;

% Theoretical significance test (adjusting for multiple hypotheses).
pval_t = mvgc_pval(UNff,morder,nobs,ntrials,1,1,nvars-2,tstat);
sig_t  = significance(pval_t,alpha,mhtc);
    
    
%% Plots
fprintf('--- Plots... --- \n\n');
% Plot original & surrogate (time.reversed)
% figure(5); clf;
% subplot(1,3,1);
% plot_pw(ff);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])

% Plot surrogate (permutation)
figure(); clf;
subplot(2,3,1);
plot_pw(ff);
title('Pairwise-conditional GC');
subplot(2,3,2);
plot_pw(pval_t);
title({'p-values';'(theoretical)'});
subplot(2,3,3);
plot_pw(sig_t);
title({['Significant at p = ' num2str(alpha)];'(theoretical)'});
subplot(2,3,5);
plot_pw(pval_p);
title({'p-values';'(perm test)'});
subplot(2,3,6);
plot_pw(sig_p);
title({['Significant at p = ' num2str(alpha)];'(perm test)'});

fprintf('--- ... Finish --- \n');
