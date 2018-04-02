% =========================================================================
% Implementation of SSGC Script
% =========================================================================


%% State-space - SSGC toolbox all ROIS - each conditions - trials

[nodes,blockDur,nblocks]= size(BOLD_data);  % number of nodes, number of points, number of blocks

%         fprintf('    - Performing SSGC for condition %s... --- \n',cond_names{cond});
        
        for blo=1:nblocks
            fprintf('        - Block number: %d \n',blo);
            % ... normalization ...     
            BOLD_dataNorm = demean(BOLD_data(:,:,blo)); % no constant term
            
            % ... initialization ...
            diagonal  = nan(nodes,1)';
            Results.G_AllSubs_FDR{aux,i}(:,:,blo)   = diag(diagonal); 
            Results.consist{aux,i}(blo) = nan;            
            
            % ... trying to perform SSGC ...
            try
                ResultsSSGC = ssgc_function_vfMRI( blo, ... % number of block
                                           BOLD_dataNorm, ... % data (alteraçao do nºpontos)
                                           mordermax,     ... % maximum model order for model mation
                                           []);               % number of subjects (define only if n_sub > 1, otherwise n_sub = [])
                
                Results.G_AllSubs_FDR{aux,i}(:,:,blo)  = ResultsSSGC.G;
                Results.consist{aux,i}(blo)            = ResultsSSGC.cons;
                test.rshapData{aux,i}(:,:,blo)         = BOLD_dataNorm;  % save data normalized
                ResultsSSGC.Morder{aux,i}              = ResultsSSGC.morder;
                clear ResultsSSGC BOLD_dataNorm
            catch ME
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SSGC failed')
                failures = failures + 1;
            end
        end




