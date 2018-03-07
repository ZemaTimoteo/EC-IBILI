function [ outputs ] = mvgc_demo_function( nobs , fs , X , nvars , alpha, ntrials , momax , mhtc  , morder, n_sub , cond_name , graphs1 , graphs2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parameters

% ntrials   = 1;     % number of trials
% nobs      = 300;   % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

% morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
% momax     = 5;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
% alpha     = 0.01;   % significance level for significance test
% mhtc      = 'None';  % multiple hypothesis test correction (see routine 'significance')

% fs        = 0.5;    % sample rate (Hz)
% fres      = [];     % frequency resolution (empty for automatic calculation)

% seed      = 0;      % random seed (0 for unseeded)

outputs = cell(0);

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

% ptic('\n*** tsdata_to_infocrit\n');
if isempty(n_sub)
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit_v2(X,momax,icregmode,false,[]);
else
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit_v2(X,momax,icregmode,false,n_sub);
end
% ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.
if graphs1 && momax > 1
    figure(); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');
end

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

% Select model order.
if strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

outputs.morder = morder;

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

% ptic('\n*** tsdata_to_var... ');
if isempty(n_sub)
    [A,SIG] = tsdata_to_var(X,morder,regmode);
else
    [A,SIG] = tsdata_to_var(X,morder,regmode);
end
% ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.
%% Consistency
outputs.cons = consistency(X,SIG);

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

% ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
% ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

% ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
% ptoc;
outputs.F = F;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

outputs.pval = pval;
outputs.sig = sig;

% Plot time-domain causal graph, p-values and significance.
if graphs2
 
    figure(); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])
    
    % Title for the figure
    suptitle(cond_name);

end
% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

outputs.cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',outputs.cd);

%% Significant In / out flow
outputs.totalConnections = sum(sig(~isnan(sig)));
outputs.outConnections = nansum(sig,1);
outputs.inConnections = nansum(sig,2);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% % Calculate spectral pairwise-conditional causalities at given frequency
% % resolution - again, this only requires the autocovariance sequence.
% 
% % ptic('\n*** autocov_to_spwcgc... ');
% f = autocov_to_spwcgc(G,fres);
% % ptoc;
% 
% outputs.f = f;
% % Check for failed spectral GC calculation
% 
% assert(~isbad(f,false),'spectral GC calculation failed');
% 
% % Plot spectral causal graph.
% if graphs2
%     figure(); clf;
%     plot_spw(f,fs);
% end
% 
% %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)
% 
% % Check that spectral causalities average (integrate) to time-domain
% % causalities, as they should according to theory.
% 
% % fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
% Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
% mad = maxabs(F-Fint);
% madthreshold = 1e-5;
% if mad < madthreshold
%     fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
% else
%     fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
% end
%
% outputs.Fint = Fint;
% %%
% % <mvgc_demo.html back to top>


end

