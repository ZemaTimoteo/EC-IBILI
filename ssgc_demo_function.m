function [ outputs ] = ssgc_demo_function(nobs,X,nvars,alpha,momax,mhtc,n_sub,nblocks,cond_name,graphs1,graphs2,tstat,checkConsistency)
% SSGC TOOLBOX based

%% Building state-space

% morder selection - past/future horizons for canonical correlations (pf)
icregmode = 'OLS';
if isempty(n_sub)
    [~,~,moAIC,~] = tsdata_to_infocrit_v2(X,momax,icregmode,[],[]);
else
    [~,~,moAIC,~] = tsdata_to_infocrit_v2(X,momax,icregmode,[],n_sub);
end
pf = 2*moAIC;
morder = moAIC;
%     morder    = 3;%moBIC;

failSSGC_build = 0;
i_failSS = [];
ii = [];

for i=1:nblocks
    try
    [outputs.m,outputs.A(:,:,i),outputs.C(:,:,i),outputs.K(:,:,i),outputs.V(:,:,i),~,~] = s4sid_CCA(X(:,:,i),...  % time-series
                                                            pf,...   
                                                            []);          % model order     
    catch ME
        failSSGC_build = failSSGC_build + 1;
        if i<nblocks
            i_failSS = [i_failSS i];
        end
        if i == nblocks
            [a,aa,~]=size(outputs.A);[c,cc,~]=size(outputs.C);[k,kk,~]=size(outputs.K);[v,vv,~]=size(outputs.V);
            outputs.A(:,:,i) = zeros(a,aa);
            outputs.C(:,:,i) = zeros(c,cc);
            outputs.K(:,:,i) = zeros(k,kk);
            outputs.V(:,:,i) = zeros(v,vv);
            i_failSS = [i_failSS i];            
        end
    end  
end
if size(i_failSS,2)>0
    fprintf('       - Failed to build SS in %d percent\n',round((failSSGC_build/nblocks)*100));
    outputs.A(:,:,i_failSS)=[];
    outputs.C(:,:,i_failSS)=[];
    outputs.K(:,:,i_failSS)=[];
    outputs.V(:,:,i_failSS)=[];
end


%% consistency test - Seth based
% if checkConsistency
%     fprintf('--- Consistency... --- \n');
%     sum = 0;
%     for ii = 1:size(outputs.V,3)
%         cons = consistency(X(:,:,ii),outputs.K(:,:,ii).*outputs.V(:,:,ii));
%         figure(14), plot(ii,cons,'-ro');  % plot
%         hold on
%         grid on
%         title('Consistency statistics');
%         ylabel('Consistency %');
%         xlabel('Evaluation points (each ii = 50 datapoints)');
%         
%         if cons>= 0.8
%             sum = 1+sum;
%         end
%     end
%     fprintf('     -> %d points with consistency bigger than 0.8 \n\n',sum);
% end

%% GC for state-space
% fprintf('--- GC for ss ... --- \n\n');

failSSGC_GC = 0;
i_failGC = [];

for i = 1:size(outputs.A,3)
%     try
%     fprintf('GC State-space for block: %d \n\n',i);
    outputs.G(:,:,i) = iss_PWGC(outputs.A(:,:,i),outputs.C(:,:,i),outputs.K(:,:,i),outputs.V(:,:,i)); % valores de GC steady-space
                                                      
%     catch ME
%         disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SSGC failed !')
%         failSSGC_GC = failSSGC_GC + 1;
%         if i<size(outputs.A,3)
%             i_failGC = [i_failGC i];
%             ii=i;
%         end
%       
%     end
end
if size(i_failGC,2)>0
    value = round((i_failGC/nblocks)*100);
    fprintf('       - Failed to perform GC-SS in %d percent\n',value);
end

%% Significance test
% fprintf('--- Significance test... --- \n');
outputs.pval = mvgc_pval(nanmean(outputs.G,3),morder,nobs,nblocks,1,1,nvars-2,tstat); % take careful note of arguments!
outputs.sig  = significance(outputs.pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.
if graphs2
 
    figure(); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(outputs.pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(outputs.sig);
    title(['Significant at p = ' num2str(alpha)])
    
    % Title for the figure
    suptitle(cond_name);
end

