% clear
% close all
% addpath('mvgc_new')

% %% Config (manual)
% fs      = 1000;         %In Hz
% t       = 200;          %In seconds
% l       = 20;           %In ms
% L       = t*fs;
% n       = 61 : L+60;
% nvars   = 5;
% j =1; blocksVect(j) = 25;
% test    = 3;            % 1 - 'MVGC', 2 - 'SSGC' reshape, 3  'SSGC' trials


%% Config (automatic - see spec_sens_testScript)
fs      = 1000;          %In Hz
t       = pointsVect(point); %In seconds
l       = 20;            %In ms
L       = t*fs;
n       = 61 : L+60;
nvars   = var;
% ntrials = 25;

%% Generate the simulated data
for trial = 1:blocksVect(block)
    fprintf('    - Performing trial %d... \n',trial);

    %% Model
    if var == 5
        y1_ti = zeros(L+60,1); y2_ti = zeros(L+60,1); y3_ti = zeros(L+60,1);
        y4_ti = zeros(L+60,1); y5_ti = zeros(L+60,1);
        w_ti  = randn(L+60,nvars); %  normal-distributed white noise
        
        
        y1_ti(n) = 0.95*sqrt(2)*y1_ti(n-l) - 0.9025*y1_ti(n-2*l) + w_ti(n,1);
        y2_ti(n) = 0.5*y1_ti(n-2*l) + w_ti(n,2);
        y3_ti(n) = -0.4*y1_ti(n-3*l) + w_ti(n,3);
        y4_ti(n) = -0.5*y1_ti(n-2*l) + 0.25*sqrt(2)*y4_ti(n-l) + ...
            0.25*sqrt(2)*y5_ti(n-l) + w_ti(n,4);
        y5_ti(n) = -0.25*sqrt(2)*y4_ti(n-l) + 0.25*sqrt(2)*y5_ti(n-l) + ...
            w_ti(n,5);
        
        %Discard first 60 points and concatenate
        X = [y1_ti(61:end) , y2_ti(61:end) , y3_ti(61:end) , y4_ti(61:end) , y5_ti(61:end) ];
        
        ord = [3.85 2.3 3 2.4 2.7]*2.5;  % for HRF

        clear y1_ti y2_ti y3_ti y4_ti y5_ti

        label_ec = [0 0 0 0 0; ...
                    2 0 0 0 0; ...
                    2 0 0 0 0; ...
                    2 0 0 0 2; ...
                    0 0 0 2 0];
        
    elseif var == 9 && varComplex == 0
        y1_ti = zeros(L+60,1); y2_ti = zeros(L+60,1); y3_ti = zeros(L+60,1);
        y4_ti = zeros(L+60,1); y5_ti = zeros(L+60,1); y6_ti = zeros(L+60,1);
        y7_ti = zeros(L+60,1); y8_ti = zeros(L+60,1); y9_ti = zeros(L+60,1);
        w_ti  = randn(L+60,nvars); %  normal-distributed white noise
        
        
        y1_ti(n) = 0.4088 * y3_ti(n-l) + 0.4236 * y4_ti(n-l) + 0.2114 * y1_ti(n-l) + ...
                    0.2590 * y3_ti(n-2*l) + 0.1965 * y4_ti(n-2*l) + 0.2148 * y1_ti(n-2*l) + ...
                    0.1434 * y3_ti(n-3*l) + 0.1787 * y4_ti(n-3*l) + 0.0976 * y1_ti(n-3*l) + ...
                    w_ti(n,1);
        y2_ti(n) = -0.3394 * y1_ti(n-l) + 0.2799 * y4_ti(n-l) + 0.3085 * y8_ti(n-l) - 0.3192 * y2_ti(n-l) + ...
                    0.2232 * y1_ti(n-2*l) + 0.1779 * y4_ti(n-2*l) - 0.2383 * y8_ti(n-2*l) - 0.1270 * y2_ti(n-2*l) + ...
                    0.1082 * y1_ti(n-3*l) + 0.1351 * y4_ti(n-3*l) + 0.2021 * y8_ti(n-3*l) - 0.1065 * y2_ti(n-3*l) + ...
                    w_ti(n,1);
        y3_ti(n) = 0.3437 * y4_ti(n-l) + 0.2944 * y3_ti(n-l) + ...
                    0.2008 * y4_ti(n-2*l) - 0.0958 * y3_ti(n-2*l) + ...
                    0.1826 * y4_ti(n-3*l) + 0.0950 * y3_ti(n-3*l) + ...
                    w_ti(n,1);
        y4_ti(n) = 0.4302 * y2_ti(n-l) + 0.4188 * y4_ti(n-l) + ...
                    0.2103 * y2_ti(n-2*l) + 0.1312 * y4_ti(n-2*l) + ...
                    0.2090 * y2_ti(n-3*l) - 0.1837 * y4_ti(n-3*l) + ...
                    w_ti(n,1);
        y5_ti(n) = 0.2851 * y3_ti(n-l) + 0.3160 * y7_ti(n-l) + 0.0827 * y5_ti(n-l) + ...
                    0.1597 * y3_ti(n-2*l) + 0.1989 * y7_ti(n-2*l) + 0.0965 * y5_ti(n-2*l) + ...
                    0.1943 * y3_ti(n-3*l) - 0.0698 * y7_ti(n-3*l) + 0.2197 * y5_ti(n-3*l) + ...
                    w_ti(n,1);
        y6_ti(n) = 0.3039 * y5_ti(n-l) - 0.0705 * y6_ti(n-l) + ...
                    0.2062 * y5_ti(n-2*l) + 0.3177 * y6_ti(n-2*l) + ...
                    0.1347 * y5_ti(n-3*l) - 0.1813 * y6_ti(n-3*l) + ...
                    w_ti(n,1);
        y7_ti(n) = 0.3030 * y6_ti(n-l) - 0.2186 * y7_ti(n-l) + ...
                    0.1895 * y6_ti(n-2*l) - 0.0805 * y7_ti(n-2*l) + ...
                    0.2270 * y6_ti(n-3*l) - 0.0404 * y7_ti(n-3*l) + ...
                    w_ti(n,1);
        y8_ti(n) = 0.3477 * y9_ti(n-l) - 0.1284 * y8_ti(n-l) + ...
                    0.1381 * y9_ti(n-2*l) - 0.2432 * y8_ti(n-2*l) + ...
                    0.1397 * y9_ti(n-3*l) - 0.1423 * y8_ti(n-3*l) + ...
                    w_ti(n,1);
        y9_ti(n) = 0.3037 * y7_ti(n-l) + 0.2810 * y8_ti(n-l) - 0.2083 * y9_ti(n-l) +  ...
                    0.1947 * y7_ti(n-1*l) - 0.1718 * y8_ti(n-2*l) + 0.6821 * y9_ti(n-2*l) + ...
                    0.2020 * y7_ti(n-3*l) + 0.1417 * y8_ti(n-3*l) - 0.2221 * y9_ti(n-3*l) + ...
                    w_ti(n,1);
                
        %Discard first 60 points and concatenate
        X = [y1_ti(61:end) , y2_ti(61:end) , y3_ti(61:end) , y4_ti(61:end) , ...
            y5_ti(61:end) , y6_ti(61:end) , y7_ti(61:end) , y8_ti(61:end) , ...
            y9_ti(61:end)];

        ord = [2.45 2.4 3.0 3.85 2.2 2.66 2.91 2.85 2.7]*2.5; % for the HRF

        clear y1_ti y2_ti y3_ti y4_ti y5_ti y6_ti y7_ti y8_ti y9_ti

        label_ec = [0 0 2 2 0 0 0 0 0; ...
                    2 0 0 2 0 0 0 2 0; ...
                    0 0 0 2 0 0 0 0 0; ...
                    0 2 0 0 0 0 0 0 0; ...
                    0 0 2 0 0 0 2 0 0; ...
                    0 0 0 0 2 0 0 0 0; ...
                    0 0 0 0 0 2 0 0 0; ...
                    0 0 0 0 0 0 0 0 2; ...
                    0 0 0 0 0 0 2 2 0];  
                
    elseif var == 9 && varComplex == 1
        y1_ti = zeros(L+60,1); y2_ti = zeros(L+60,1); y3_ti = zeros(L+60,1);
        y4_ti = zeros(L+60,1); y5_ti = zeros(L+60,1); y6_ti = zeros(L+60,1);
        y7_ti = zeros(L+60,1); y8_ti = zeros(L+60,1); y9_ti = zeros(L+60,1);
        w_ti  = randn(L+60,nvars); %  normal-distributed white noise
        
        y1_ti(n) = 0.95*sqrt(2)*y1_ti(n-l) - 0.9025*y1_ti(n-2*l) + w_ti(n,1);
        y2_ti(n) = 0.5*y1_ti(n-2*l) + w_ti(n,2);
        y3_ti(n) = -0.4*y1_ti(n-3*l) - 0.4381 * y8_ti(n-2*l) + ...
            0.4397 * y8_ti(n-3*l) + + w_ti(n,3);
        y4_ti(n) = -0.5*y1_ti(n-2*l) + 0.25*sqrt(2)*y4_ti(n-l) + ...
            0.25*sqrt(2)*y5_ti(n-l) + w_ti(n,4);
        y5_ti(n) = -0.25*sqrt(2)*y4_ti(n-l) + 0.25*sqrt(2)*y5_ti(n-l) + ...
            w_ti(n,5);
        y6_ti(n) = 0.3039 * y9_ti(n-l) - 0.9705 * y6_ti(n-l) + ...
            0.2062 * y9_ti(n-2*l) - 0.75*y4_ti(n-2*l) - 0.9177 * y6_ti(n-2*l) - ...
            0.1813 * y6_ti(n-3*l) + ...
            w_ti(n,1);
        y7_ti(n) = 0.3030 * y5_ti(n-l) - 0.2186 * y7_ti(n-l) - ...
            + 0.0805 * y7_ti(n-2*l) + ...
            0.2270 * y5_ti(n-3*l)  + ...
            w_ti(n,1);
        y8_ti(n) = 0.3477 * y3_ti(n-l) + ...
            0.4381 * y3_ti(n-2*l) + ...
            0.4397 * y3_ti(n-3*l) + ...
            w_ti(n,1);
        y9_ti(n) = - 0.2083 * y9_ti(n-l) +  0.25*sqrt(2)*y5_ti(n-l) + ...
            0.1947 * y3_ti(n-2*l) + 0.78*y5_ti(n-2*l)- 0.6821 * y9_ti(n-2*l) + ...
            0.2020 * y3_ti(n-3*l) - 0.2221 * y9_ti(n-3*l) + ...
            w_ti(n,1);
        
        %Discard first 60 points and concatenate
        X = [y1_ti(61:end) , y2_ti(61:end) , y3_ti(61:end) , y4_ti(61:end) , ...
            y5_ti(61:end) , y6_ti(61:end) , y7_ti(61:end) , y8_ti(61:end) , ...
            y9_ti(61:end)];
        
        ord = [2.45 2.4 3.0 3.85 2.2 2.66 2.91 2.85 2.7]*2.5; % for the HRF
        
        clear y1_ti y2_ti y3_ti y4_ti y5_ti y6_ti y7_ti y8_ti y9_ti
        
        label_ec = [0 0 0 0 0 0 0 0 0; ...
                    2 0 0 0 0 0 0 0 0; ...
                    2 0 0 0 0 0 0 0 0; ...
                    2 0 0 0 2 0 0 0 0; ...
                    0 0 0 2 0 0 0 0 0; ...
                    0 0 0 2 0 0 0 0 2; ...
                    0 0 0 0 2 0 0 0 0; ...
                    0 0 2 0 0 0 0 0 0; ...
                    0 0 2 0 2 0 0 0 0];  
        
    end
    
    
    %% Plot model
    if plots == 1
        figure (12)
        for i=1:var
            subplot(var,1,i)
            plot(0:1/1000:(size(X,1)-1)/1000,X(:,i))
            ylim([-5 5])
            if i==5, xlabel('Time(s)'), end
            title(sprintf('Source %i',i))
        end
    end
    %% Process HRF
    
    clear BOLD BOLD_true hrf
    if plots == 1; figure(); end
    for i = 1:nvars
        p = [ord(i) 16 1 1 6 0 32];
        hrf(:,i) = spm_hrf(1/1000,p);
        
        BOLD_true(:,i) = conv(X(:,i) , hrf(:,i));
        BOLD(:,i) = awgn(BOLD_true(:,i),ratioSNR);
        
        if plots == 1
            plot(BOLD_true(:,i))
            hold on
            if i==var, xlabel('Time(s)'), end
            title(sprintf('BOLD Simulated Signal',i))
        end
    end
    
    if plots == 1
        hold on
        plot(BOLD(:,2))
    end
    
    % Remove first 50s of data
    BOLD = BOLD((32+50)*fs+1:end,:);
    BOLD_Original(:,:,trial) = BOLD';

    clear X 
    %% HRF plot
    if plots == 1
        figure()
            plot(0:1/1000:(size(hrf,1)-1)/1000,hrf)
            title('HRF at 1000Hz')
            legend({'1','2','3','4','5','6','7','8','9'})
            xlabel('Time (s)')
            xlim([0 (size(hrf,1)-1)/1000])
            line([0 (size(hrf,1)-1)/1000],[0 0],'Color',[0.5 0.5 0.5])
    end

    %% Pos downsample
    BOLD_250(:,:,trial)      = downsample(BOLD,4)';
    BOLD_05(:,:,trial)       = downsample(BOLD,2000)';
    
end

    %% Correlation plot
figure()
if plots == 1
    subplot(1,3,1)
    X_dat_avg_original = mean(BOLD,3);
    R_original = corrcoef(X_dat_avg_original);
    imagesc(abs(R_original));colormap('gray');colorbar;
    xlabel('Channels'); ylabel('Channels')
    title('Original')
    caxis([0 1.0]);
    
    subplot (1,3,2)
    X_dat_avg_250 = mean(BOLD,3);
    R_250 = corrcoef(X_dat_avg_250);
    imagesc(abs(R_250));colormap('gray');colorbar;
    xlabel('Channels'); ylabel('Channels')
    title('250 Hz')
    caxis([0 1.0]);
    
    subplot (1,3,3)
    X_dat_avg_05 = mean(BOLD,3);
    R_05 = corrcoef(X_dat_avg_05);
    % Plot of matrix R
    imagesc(abs(R_05));colormap('gray');colorbar;
    xlabel('Channels'); ylabel('Channels')
    title('0.5 Hz')
    caxis([0 1.0]);
    
    titleP = ['Cross-Correlation - 9 Var'];
    p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',.012);
    set(p.th,'edgecolor',.5*[1 1 1]);
    
    R_origi = mean(mean(R_original));
    R_250 = mean(mean(R_250));
    R_05 = mean(mean(R_05));
end




%%
if test == 1 % MVGC
    %% Initialization
    clear G_AllSubs_perConds_FDR DataMVGC_II
    cond_names = {'simulate05','simulate250','simulateOriginal'};
    DataMVGC_II.simulate05 = BOLD_05;
    DataMVGC_II.simulate250 = BOLD_250;
    DataMVGC_II.simulateOriginal = BOLD_Original;
    n_cond = 3;
    fs     = [0.5 250 1000];
    
    %% test
    for cond = 1:n_cond
        failures = 0;
        [nnrois,blockDur(cond),nblocks]= size(DataMVGC_II.(cond_names{cond}));
        % ... normalization ...
        DataMVGC_II.(cond_names{cond}) = demean(DataMVGC_II.(cond_names{cond})); % no constant term
        fprintf('\n\n\n   => Performing MVGC for condition %s... \n',cond_names{cond});
        for blo=1:nblocks
            try
                G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = nan(nnrois);
                consist.(cond_names{cond})(block,blo,point) = nan;

                ResultsMVGC = mvgc_demo_function( blockDur(cond) , ... % number of observations per trial
                    fs(cond) , ...  % sample rate (Hz)
                    DataMVGC_II.(cond_names{cond})(:,:,blo)  , ... % data
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
                G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = ResultsMVGC.F;
                consist.(cond_names{cond})(block,blo,point) = ResultsMVGC.cons;

                
            catch
                fprintf('----> MVGC1 failed / Trial %i\n.',blo)
                G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = zeros(5);
                failures = failures + 1;
            end
            
        end
    end
    
end






if test == 2 % SSGC reshape in segments
    %% State-space - SSGC toolbox all ROIS - each conditions - with segments
    failures = 0;
    clear G_AllSubs_perConds_FDR DataMVGC_II
    cond_names = {'simulate05','simulate250','simulateOriginal'};
    n_cond = 3;
    DataMVGC_II.simulate05 = BOLD_05;
    DataMVGC_II.simulate250 = BOLD_250;
    DataMVGC_II.simulateOriginal = BOLD_Original;
    
    for cond = 1:n_cond
        fprintf('\n    - Performing SSGC for condition %s... \n',cond_names{cond});
        [nnrois,blockDur(cond),nblocks]= size(DataMVGC_II.(cond_names{cond}));
        % normalization
        DataMVGC_II.(cond_names{cond}) = demean(DataMVGC_II.(cond_names{cond})); % no constant term
        spFactor = 40; aux = floor((blockDur*spFactor)*(nblocks/spFactor));
        DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*nblocks);
        if aux < blockDur(cond)
            DataMVGC_II.(cond_names{cond})(:,aux+1:end) = [];
        end
        if floor(nblocks/spFactor) < 2
            Newblocks = 1;
        else
            Newblocks = floor(nblocks/spFactor);
        end
        %             DataMVGC_II.(cond_names{cond}) = reshape (DataMVGC_II.(cond_names{cond}),nnrois,blockDur(cond)*spFactor,floor(nblocks/spFactor));
        [~,blockDur(cond),nblocks] = size(DataMVGC_II.(cond_names{cond}));
        sizeblocks = blockDur(cond);
        
        for blo=1:nblocks
            fprintf('        - Block number: %d \n',blo);
            
            try
                ResultsSSGC = ssgc_demo_function_v2( blo, ...  % number of block
                    DataMVGC_II.(cond_names{cond})(:,1:sizeblocks,blo) , ... % data (alteraçao do nºpontos)
                    10 , ... % maximum model order for model mation
                    []); % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                
                G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = ResultsSSGC.G;
                morder.(cond_names{cond})(blo) = ResultsSSGC.morder;
                
            catch ME
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed')
                failures = failures + 1;
            end
        end
        
        % Significance test
        G_AllSubs_perConds_FDR.(cond_names{cond}) = nanmean(G_AllSubs_perConds_FDR.(cond_names{cond}),3);
        ResultsSSGC.pval = mvgc_pval(G_AllSubs_perConds_FDR.(cond_names{cond}),1,sizeblocks,nblocks,1,1,nnrois-2,''); % take careful note of arguments!
        SSsig_AllSubs.(cond_names{cond})  = significance(ResultsSSGC.pval,0.05,'FDR');
        
    end
end






if test == 3 % SSGC trials
    %% State-space - SSGC toolbox all ROIS - each conditions - trials
    failures = 0;
    clear G_AllSubs_perConds_FDR DataMVGC_II
    cond_names = {'simulate05','simulate250','simulateOriginal'};
    n_cond = 3;
    DataMVGC_II.simulate05 = BOLD_05;
    DataMVGC_II.simulate250 = BOLD_250;
    DataMVGC_II.simulateOriginal = BOLD_Original;
    
    for cond = 1:n_cond
        fprintf('    - Performing SSGC for condition %s... --- \n',cond_names{cond});
        [nnrois,blockDur(cond),nblocks]= size(DataMVGC_II.(cond_names{cond}));
        % ... normalization ...
        DataMVGC_II.(cond_names{cond}) = demean(DataMVGC_II.(cond_names{cond})); % no constant term
        
        [~,blockDur(cond),nblocks] = size(DataMVGC_II.(cond_names{cond}));
        sizeblocks = blockDur(cond);
        
        for blo=1:nblocks
            fprintf('        - Block number: %d \n',blo);
            diagonal  = nan(var,1)';
            consist.(cond_names{cond})(block,blo,point) = nan;
            G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = diag(diagonal);
            try
                ResultsSSGC = ssgc_demo_function_v2( blo, ...  % number of block
                    DataMVGC_II.(cond_names{cond})(:,1:sizeblocks,blo) , ... % data (alteraçao do nºpontos)
                    mordermax , ... % maximum model order for model mation
                    []); % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                
                G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = ResultsSSGC.G;
                consist.(cond_names{cond})(block,blo,point) = ResultsSSGC.cons;
                moOrder.(cond_names{cond})(block,blo,point)  = ResultsSSGC.morder;

            catch ME
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed')
                failures = failures + 1;
            end
        end
    end
    
end



%% Coherence
if test == 4
    failures = 0;
    clear DataMVGC_II G_AllSubs_perConds_FDR
    cond_names = {'simulate05','simulate250','simulateOriginal'};
    n_cond = 3;
    DataMVGC_II.simulate05       = BOLD_05;
    DataMVGC_II.simulate250      = BOLD_250;
    DataMVGC_II.simulateOriginal = BOLD_Original;
    fs = [0.5 250 1000];
    c  = combnk(1:5,2);
    
    for cond = 1:n_cond
        fprintf('    - Performing Coherence for condition %s... --- \n',cond_names{cond});
        [nnrois,~,nblocks]= size(DataMVGC_II.(cond_names{cond}));
        % ... normalization ...
        for blo=1:nblocks
            fprintf('        - Block number: %d \n',blo);
            G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = nan(nnrois);
            try
                Coh = zeros(nnrois);
                for comb = 1:size(c,1)
                    [tempCoh,tempF] = mscohere(DataMVGC_II.(cond_names{cond})(c(comb,1),:,blo)',DataMVGC_II.(cond_names{cond})(c(comb,2),:,blo)',[],[],[],fs(cond));
                    Coh(c(comb,2),c(comb,1)) = mean(tempCoh(tempF<0.15));
                end
            catch ME
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Coherence failed')
                failures = failures + 1;
            end
            
            G_AllSubs_perConds_FDR.(cond_names{cond})(:,:,blo) = atanh(Coh);
        end
    end
end




if test == 2 % SSGC
    surrugate_demo_function
end




% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================




if strcmp (test,'do')
    %% Pearson Correlation
    results_pc_1(:,:,trial) = corrcoef(BOLD_250);
    results_pc_2(:,:,trial) = corrcoef(BOLD_05);
    
    fprintf('--------> Ended %i trial\n',trial)
    
    
    %% Final Plots MVGC
    res1 = mean(results1,3);
    res2 = mean(results2,3);
    
    plot_GC = figure();
    subplot(1,2,1)
    plot_pw(res1)
    title('250Hz')
    subplot(1,2,2)
    plot_pw(res2)
    title('0.5Hz')
    
    saveas(plot_GC,'plot_GC.fig','fig'); % save figure for data
    
    %% Model Order
    res1_modelorder = [mean(results1_modelorder) std(results1_modelorder)];
    res2_modelorder = [mean(results2_modelorder) std(results2_modelorder)];
    
    %% Final Plots PC
    % load('colormap(-1,1).mat');
    %
    % res1_pc = mean(results_pc_1,3);
    % res2_pc = mean(results_pc_2,3);
    %
    % plot_PC = figure();
    %     subplot(1,2,1)
    %         plot_pw_pc(res1_pc,newColorMap)
    %         title('250Hz')
    %         colorbar
    %     subplot(1,2,2)
    %         plot_pw_pc(res2_pc,newColorMap)
    %         title('0.5Hz')
    %         colorbar
    %
    % saveas(plot_PC,'plot_PC.fig','fig'); % save figure for data
    
    %% Permutation Test
    % ResultsMVGC3 = mvgc_permtest_function( size(BOLD_250',2) , ... % number of observations per trial
    %                 250 , ...  % sample rate (Hz)
    %                 BOLD_250' , ... % data
    %                 nvars , ... % number of rois
    %                 0.05 , ... % significance level for significance test
    %                 1 , ... % number of trials
    %                 50 , ... % maximum model order for model order estimation
    %                 'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
    %                 [] , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
    %                 100, ... %number of permutations
    %                 'All Conditions' , ... % condition names (or plot title in case of #conds = 1)
    %                 0 , ... % Display model order plot
    %                 1); % Display results plot
    %
    % ResultsMVGC4 = mvgc_permtest_function( size(BOLD_05',2) , ... % number of observations per trial
    %                 0.5 , ...  % sample rate (Hz)
    %                 BOLD_05' , ... % data
    %                 nvars , ... % number of rois
    %                 0.05 , ... % significance level for significance test
    %                 1 , ... % number of trials
    %                 50 , ... % maximum model order for model order estimation
    %                 'FDR' , ... % multiple hypothesis test correction (see routine 'significance')
    %                 [] , ... % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
    %                 100, ... %number of permutations
    %                 'All Conditions' , ... % condition names (or plot title in case of #conds = 1)
    %                 0 , ... % Display model order plot
    %                 1); % Display results plot
    
    %% Save Workspace
    save('workspaceVARDataGenerator')
    
    %% lixo
    % %% Window Pearson Corr
    % window = 48;
    % n_rois = 5;
    % n_con = (n_rois*(n_rois-1))/2; % #Connections
    % C = combnk(1:n_rois,2); % Possible Combinations
    %
    % load('colormap(-1,1).mat');
    % mean_results = zeros(n_rois);
    % mean_resultsF = zeros(n_rois);
    % figure
    %
    % %Combinations Iteration
    % for l=1:n_con
    %
    %     temp1 = BOLD_pos_250(:,C(l,1));
    %     temp2 = BOLD_pos_250(:,C(l,2));
    %
    %     %Window Iteration
    %     for t = window : 45500
    %
    %         x = t-window+1 : t;
    %         temp3 = corrcoef(temp1(x),temp2(x));
    %
    %         Results.VOI_Corr{l}(t,1) = temp3(1,2);
    %         Results.VOI_CorrF{l}(t,1) = atanh(temp3(1,2)); %Fisher
    %
    %     end
    %
    %     a = Results.VOI_Corr{1,l};
    %     b = Results.VOI_CorrF{1,l};
    %
    %     subplot(n_con,1,l), plot(a);
    %
    %     mean_results(C(l,1),C(l,2)) = mean(a);
    %     mean_results(C(l,2),C(l,1)) = mean(a);
    %
    %     mean_resultsF(C(l,1),C(l,2)) = mean(b);
    %     mean_resultsF(C(l,2),C(l,1)) = mean(b);
    %
    % end
    %
    % %%
    % figure
    %     plot_pw_pc(mean_results,newColorMap)
    %     colorbar
    %
    % figure
    %     plot_pw_full(mean_results,newColorMap)
    %     colorbar
    %
end