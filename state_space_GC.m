% === Modelling state-space ===

%% load
clear all, clc
fprintf('--- Load... --- \n\n');
addpath(genpath('C:\Users\admin\Documents\Tiago\IBILI\Code\Toolboxes'))
load('C:\Users\admin\Documents\Tiago\IBILI\Data\testes\test_ti_simul\testData_noVC.mat')

%% Initialization
[nvars,nobs,ntrials] = size (X);
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
graph     = 'true'; % do the plot ('true') it does not ('false')

%% teste

% y = X; % Data
% K = qq_coisa;  % Kalman gain
% 
% % for i=2:fim
% %     x(i+1) = A*x(i) + u(i);    % state transition equation
% %     y(i)   = C*x(i) + v(i);    % observation equation
% % end
% 
% 
% uvRES = cat(1,u,v);
% noiseRES = cov (uvRES,uvRES'); % noise covariance matrix
% 
% z = qq_coisa; % state variable
% 
% % state space
% for i=2:fim
%     z(i)   = (((eye(qq_tamanho) - A*z)^-1)*Kz)*e(i);
%     z(i+1) = A*z(i) + K*e(i);                                      % state transition equation
%     B      = A - K*C;
%     z(i+1) = B*z(i) + K*y(i);
%     H(z)   = eye(qq_tamanho) + C*((eye(qq_tamanho) - A*z)^-1)*K*z; % transfer function
%     B(z)   = inv(H);
%     y(i)   = H(z)*e(i);
%     y(i)   = C*z(i) + e(i);                                        % observation equation
%     [S,H]  = var_to_cpsd(A,SIG,fres);                              % spectral factorisation
%     S      = tsdata_to_cpsd(X,fres,'WELCH');                       % spectral estimation
%     %  discrete-time algebraic Riccati equations - DARE 
%     [~,~,K,E] = dare_ECt(A,C,Q,R,S);
%     
% end
% 
% R = cov(K*e(i),e(i));
% Q = K*R*K';
% S = K*R;

%% Using Dual-Extended kalman filter
% load
% clear all, clc
% load('C:\Users\admin\Documents\Tiago\IBILI\Data\testes\test1_ti_simul\testData_noVC.mat')

% % struct
% inp_model.data   = X(:,:,1)';
% inp_model.methd_ord = 'BIC'; % Method of selecting the model order
% inp_model.momax  = 10;
% inp_model.srate  = 1000;
% acmaxlags        = 500;
% 
% tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
% alpha     = 0.05;   % significance level for significance test
% mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
% graph     = 'true'; % do the plot ('true') it does not ('false')
% 
% % build data
% % data_plot = iddata(X(:,:,500)',[],1/1000);
% 
% % model order selection
% inp_model.order = 3;%morder_ECt( X , inp_model );
% 
% % Dual Kalman Filter - parameter estimation
% [A_full,SIG_full] = DEKF3(inp_model);
% [CH,L,T] = size(A_full);
% A_full = reshape(A_full,CH,CH,L/CH,T);
% 
% % consistency test - Seth based
% for ii = 1:10:size(SIG_full,3)
%     aux = SIG_full(:,:,ii:ii+9);
%     cons = consistency(X(:,ii:ii+9,1),mean(aux,3));
%    
%     figure(5), plot(ii,cons,'-ro');  % plot
%     hold on
%     grid on
%     title('Consistency statistics');
%     ylabel('Consistency %');
%     xlabel('Evaluation points (each ii = 50 datapoints)');
% end
% 
% % calculate GC time-domain
% lag_avg = 50;
% for t = inp_model.order : T % Number of time points
%     F(:,:,t) = nan(5);
%     fprintf('Trial number %d ... \n\n',t);
%     % -------- Generate Autocov coef. ------------ 
%     [G,info] = var_to_autocov(A_full(:,:,:,t),SIG_full(:,:,t),acmaxlags);
%     
%     if graph
%         compare(data,sys)
%     end
%     
%     try
%         var_info(info,true); % report results (and bail out on error)
%     catch
%         continue
%     end
% 	% -------- Generate GC values ----------------  
%     n = size(G,1);
%    
%     % full regression
%     owstate = warn_supp;
%     [~,SIG] = autocov_to_var(G);
%     warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
%     
%     if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
%     LSIG = log(diag(SIG));
%     
%     for j = 1:n
%         % reduced regression
%         jo = [1:j-1 j+1:n]; % omit j
%         owstate = warn_supp;
%         [~,SIGj] = autocov_to_var(G(jo,jo,:));
%         warn_test(owstate,     sprintf('in reduced regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
%         if warn_if(isbad(SIGj),sprintf('in reduced regression for target node %d - regression failed',j)), continue; end % show-stopper!
%         LSIGj = log(diag(SIGj));
%         for ii=1:n-1;
%             i = jo(ii);
%             F(i,j,t) = LSIGj(ii)-LSIG(i);
%         end
%     end
%         assert(~isbad(F(:,:,t),false),'GC calculation failed');
% end
% 
% 
% 
%     % Significance test using theoretical null distribution, adjusting for multiple
%     % hypotheses.
%     ff = mean(F(:,:,1:end),3,'omitnan');
%     ntrials = 1;
%     morder  = inp_model.order;
%     nobs    = 1000;%lag_avg;
%     nvars   = 5; 
%     pval = mvgc_pval(ff,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
%     sig  = significance(pval,alpha,mhtc);
%     
%     figure(2); clf;
%     subplot(1,3,1);
%     plot_pw(ff);
%     title('Pairwise-conditional GC');
%     subplot(1,3,2);
%     plot_pw(pval);
%     title('p-values');
%     subplot(1,3,3);
%     plot_pw(sig);
%     title(['Significant at p = ' num2str(alpha)])



%% Building state-space - Matlab System Identification Toolbox
% for i = 1:10
%     data = iddata(X(:,:,i)',[],1/1000);%X0',[],1); % generate data
%     opt  = ssestOptions('Focus','stability','EstCovar',true,'Display','off'); % options
%     nx   = 3;%[1:10];
%     sys  = ssest(data,...
%                  nx, ...     % state space model order / dimension
%                  'Ts',1/1000, ...
%                  'Form','canonical', ...
%                  opt);       % build data state-space
%     % sys2 = ssregest(data,nx);
%     % opt.Regularization.Lambda = 10;
%     % opt = ssestOptions('SearchMethod','lm');
%     % mr = ssest(data,3,'form','modal','Ts',data.Ts,opt);
%     % bodeplot(sys)
%     figure()
%     compare(data,sys)                      % quality control
%     
%     A(:,:,i) = sys.A;                      % state transition matrix
%     C(:,:,i) = sys.C;                      % observation matrix
%     K(:,:,i) = sys.K;                      % inovation form
%     B(:,:,i) = sys.A-sys.K*sys.C;          % Inverse of trasnfer function
%     
%     % rhoaa = specnorm(AA);
%     rhob(i) = specnorm(B(:,:,i));
%     
%     % if rhoa > 1
%     %     warning('The state space does not have a minimum phase (rhoa = %d)',rhob);
%     % end
%     if rhob(i) > 1
%         warning('The state space does not have a minimum phase (rhob(%d) = %d.3f)',i,rhob);
%     end
%     
%     % Calcular a matriz de covariancia dos residuos
%     [Estr,Rstr] = resid(data,sys);
%     E(:,i)    = mean(Estr.OutputData);
%     V(:,:,i)  = abs(cov(E(:,i)*E(:,i)'));
%     
% end


%%  Building state-space - SSGC toolbox
fprintf('--- Building state-space... --- \n\n');

for i=1:10
    % morder selection - past/future horizons for canonical correlations (pf)
    fprintf('State-space for trial: %d \n\n',i);
    momax = 10;
    icregmode = 'OLS';
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X(:,:,i),momax,icregmode);
    pf(i) = 2*moAIC;
    morder    = 3;%moBIC;
    
    [m,A(:,:,i),C(:,:,i),K(:,:,i),V(:,:,i),z,e] = s4sid_CCA(X(:,:,i),...  % time-series
                                                            pf(i),...   
                                                            []);          % model order                                                                                                       
    fprintf('\n\n');
end

%% consistency test - Seth based
fprintf('--- Consistency... --- \n');
sum = 0;
for ii = 1:size(V,3)
    cons = consistency(X(:,:,ii),K(:,:,ii).*V(:,:,ii));
    figure(14), plot(ii,cons,'-ro');  % plot
    hold on
    grid on
    title('Consistency statistics');
    ylabel('Consistency %');
    xlabel('Evaluation points (each ii = 50 datapoints)');
    
    if cons>= 0.8
        sum = 1+sum;
    end
end
fprintf('     -> %d points with consistency bigger than 0.8 \n\n',sum);

%% GC for state-space
fprintf('--- GC for ss ... --- \n\n');

for i = 1:10
    G(:,:,i) = iss_PWGC(A(:,:,i),C(:,:,i),K(:,:,i),V(:,:,i)); % valores de GC steady-space
end

%% Significance test
fprintf('--- Significance test... --- \n');
pval = mvgc_pval(mean(G,3,'omitnan'),morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

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

%% Surrugate Significance test
% --- with surrugate (time.reversed) ((û - u) / ome) ---
avg_sur = mean(G_sur,3,'omitnan');
std_sur = std(G_sur,0,3,'omitnan');
    
F_corr = (mean(G,3,'omitnan') - avg_sur) ./ std_sur;
pval_sur = mvgc_pval(F_corr,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig_sur  = significance(pval_sur,alpha,mhtc);


% --- with surrugate (time.reversed) t-test  ---
% Building antisymmetric matrix
% GCI_dat_Net=zeros(datachan,datachan,ntrials);
% GCI_sur_Net=zeros(datachan,datachan,ntrials);
% 
% for i=1:datachan
%     for j=1:datachan
%         if j~=i
%             GCI_dat_Net(i,j,:) = F(i,j,:) - F(j,i,:);
%             GCI_sur_Net(i,j,:) = surF(i,j,:) - surF(j,i,:);
%         end
%     end
% end
% 
% figure(11)
% imagesc(GCI_dat_Net(:,:,2));
% colormap(gray);
% figure(12),
% imagesc(GCI_sur_Net(:,:,2));
% colormap(gray);

% paired t-test for all channels
[h,p,ci,stats] = ttest(G,G_sur,0.05,'both',3);

% ... plot p-values ...
figure()
p_aux = p;
p_aux(find(p>0.05)) = 1;
plot_pw(p_aux);
figure(14)
plot_pw(p);
title('P-values for GCreal vs GCsurrogate');

% Convert z-scores
mu     = 0;
sigma  = 1;
prob   = tcdf(stats.tstat,9);%stats.df(1,1));  % CDF of the Students t-distribution with df degrees of freedom
z_scor = norminv(prob,mu,sigma);           % Inverse cdf of the univariate standard normal distribution
% z_scor = icdf('Normal',prob,mu,sigma);   % Inverse cdf of the univariate standard normal distribution
[i,j]       = find(z_scor == Inf);
infinit     = find(z_scor == Inf);
% if size(i,2)~=0 && size(j,2)~=0
%     diagonal    = diag(z_scor(j,i));
%     z_scor(infinit) = diagonal*(-1);
% end
% Family-wise Errors correction
nhyp   = (nvars*nvars) - nvars;
alpha  = 0.05;
alpha1 = alpha/nhyp; % Bonferroni method - alpha value
% ... Auxiliar calculus ...df
z = @(p) -sqrt(2) * erfcinv(p*2);

cval = z(1-alpha1) % critical value - z-value corrected for FWE
z_scor(find(abs(z_scor)<cval)) = 0;

% Plot with FWE correction
figure()
plot_pw(z_scor);
title('Z-scores for GCreal vs GCsurrogate');
caxis([min(min(z_scor)) max(max(z_scor))])
colormap(color_map)

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