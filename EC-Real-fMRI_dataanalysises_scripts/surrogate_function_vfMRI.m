%% =============================================================================================================


% ---- Perform time-reversed surrogate ------------

clear G_AllSubs_perConds_FDR_sur
stop = 'false';        % 'false' condition defining to perform statistical test

for blo = 1:nblocks
    fprintf('        - Block number: %d \n',blo);
    
    % ... condicao de que deu erro ao realizar o SSGC ...
    if isnan(Results.G_AllSubs_FDR{aux,u}(2,1,blo)) | Results.G_AllSubs_FDR{aux,u}(2,1,blo) == 0
        stop = 'true';
        z_scores.(cond_names{cond}) = zeros(nvars);
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
        Results.G_AllSubs_FDR_sur{aux,u}(:,:,blo) = diag(diagonal);
        
        % a) building state-space
        [m_sur,A_sur,C_sur,K_sur,V_sur,z_sur,e_sur] = s4sid_CCA(X_sur,...   % time-series
            pf,...
            [],...               % model order
            []);
        
        % b) buildng GC for surrugate
        Results.G_AllSubs_FDR_sur{aux,i}(:,:,blo)   = iss_PWGC(A_sur,C_sur,K_sur,V_sur); % valores de surrogate GC state-space
        clear m_sur A_sur C_sur K_sur V_sur z_sur e_sur pf X_sur
    catch
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed - Surrogate')
        %                     G_AllSubs_perConds_FDR_sur.(cond_names{cond})(:,:,i) = zeros(var,var);
        failures_sur = failures_sur + 1;
    end
    
    
end







