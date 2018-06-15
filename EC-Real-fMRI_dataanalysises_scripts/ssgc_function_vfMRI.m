function [ outputs ] = ssgc_demo_function_v2(trial,X,momax,n_sub)
% SSGC TOOLBOX based

%% Model order
% morder selection - past/future horizons for canonical correlations (pf)
icregmode = 'OLS';
plotMorder = true;

if momax == 1
    moAIC = momax;
    moBIC = momax;
else
    if isempty(n_sub)
        [aic,bic,moAIC,moBIC] = tsdata_to_infocrit_v2(X,momax,icregmode,false,[]);
    else
        [aic,bic,moAIC,moBIC] = tsdata_to_infocrit_v2(X,momax,icregmode,false,n_sub);
    end
end


% Plot information criteria.
% % if plotMorder
%     figure(5); clf;
%     plot_tsdata([aic bic]',{'AIC','BIC'},1/fs);
%     title('Model order estimation');
% end

pf = 2*moAIC; morder = moAIC; % when using AIC
% pf = moBIC; morder = moBIC % when using BIC
outputs.morder = morder;
%     morder    = 3;%moBIC;

%% Build linear state-space

outputs.failSSGC_build = 0;
i_failSS = [];
ii = [];
% try
    [outputs.m,outputs.A(:,:),outputs.C(:,:),...
        outputs.K(:,:),outputs.V(:,:),~,~] = s4sid_CCA(X(:,:),...  % time-series
        pf,...
        [],...     % model order SVC
        []);          
% catch ME
%     fprintf(' - Failed to build SS in block %d', trial);
% end


%% Build nonlinear state-space


%%  Consistency statistics
outputs.cons = consistency(X,outputs.V);

%% GC for state-space
outputs.G(:,:) = iss_PWGC(outputs.A(:,:),outputs.C(:,:),outputs.K(:,:),outputs.V(:,:)); % valores de GC steady-space


