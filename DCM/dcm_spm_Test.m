%.................................................
% This batch script analyses the Attention to Visual Motion fMRI dataset
% available from the SPM website using DCM. 
% It will use 'spm_dcm_ui.m'
%.................................................


% % % % Directory containing the attention data
% % % %--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
% fprintf('%-40s:', 'Downloading Attention dataset...');
% urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/attention/attention.zip','attention.zip');
% unzip(fullfile(data_path,'attention.zip'));
data_path = fullfile(data_path,'attention');
fprintf(' %30s\n', '...done');
% ola
% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION & INFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factors = load(fullfile(data_path,'factors.mat'));

f = spm_select('FPList', fullfile(data_path,'functional'), '^snf.*\.img$'); % selecionar imagens

clear matlabbatch

% OUTPUT DIRECTORY
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% MODEL SPECIFICATION
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT    = 3.22;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans            = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name     = 'Photic';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset    = [factors.att factors.natt factors.stat];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name     = 'Motion';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset    = [factors.att factors.natt];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name     = 'Attention';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset    = [factors.att];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 10;

% MODEL ESTIMATION
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

% INFERENCE
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(3);
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Photic';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [1 0 0];
matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'Motion';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 1 0];
matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'Attention';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 0 1];

spm_jobman('run',matlabbatch);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLUMES OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

% EXTRACTING TIME SERIES: V5
%--------------------------------------------------------------------------
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{1}.spm.util.voi.session = 1; % session 1
matlabbatch{1}.spm.util.voi.name = 'V5';
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 3;  % "Motion" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.contrast = 4; % "Attention" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.mtype = 0; % inclusive
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [-36 -87 -3];
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: V1
%--------------------------------------------------------------------------
matlabbatch{2}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{2}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{2}.spm.util.voi.session = 1; % session 1
matlabbatch{2}.spm.util.voi.name = 'V1';
matlabbatch{2}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{2}.spm.util.voi.roi{1}.spm.contrast = 2;  % "Photic" T-contrast
matlabbatch{2}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{2}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{2}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.centre = [0 -93 18];
matlabbatch{2}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{2}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: SPC
%--------------------------------------------------------------------------
matlabbatch{3}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{3}.spm.util.voi.session = 1; % session 1
matlabbatch{3}.spm.util.voi.name = 'SPC';
matlabbatch{3}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{3}.spm.util.voi.roi{1}.spm.contrast = 4;  % "Attention" T-contrast
matlabbatch{3}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
matlabbatch{3}.spm.util.voi.roi{1}.spm.thresh = 0.001;
matlabbatch{3}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.centre = [-27 -84 36];
matlabbatch{3}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

spm_jobman('run',matlabbatch);


% DEFINE REGIONS FROM GC - in the xY structure type
%--------------------------------------------------------------------------
% ter de construir as regioes de interesse atraves da matriz de GC
% realizar essa struct. com base na identacao realizada por SPM
data = reshape(DataMVGC_AllCond,n_rois,n_vols*n_sub);
clear Sess xY DCM
for cond = 1:n_cond
    [n_rois,pnts,blocks] = size(DataMVGC.(cond_names{cond}));
    for roi =1:n_rois
        % --- Define rois ---
        %         if ~isempty(find(SSGCResults.(cond_names{cond})(roi,:)~=0) && isnan(SSGCResults.(cond_names{cond})(roi,:)))
        %             a=1;
        %         end
        xY.name  = voiNames(roi);  % nome da ROI
        xY.Ic    = 1;    %
        xY.Sess  = 1;  % sessao
        aux      = voiFile.L.DAT.VOI(roi).Voxels;
        xY.xyz   = mean(aux,1);   % centro da ROI, posicao xyz
        xY.def   = 'mask';   % string que define o tipo de análise
        %         xY.spec  = ;  %
        %         xY.str   = ;   % string com a origem da imagem
        %         xY.XYZmm = ; % posicoes relativas dos voxeis selecionados
        %         xY.X0    = ;    %
        
        y = data(roi,:)';
        [m n]   = size(y);
        if m > n
            [v s v] = svd(y'*y);
            s       = diag(s);
            v       = v(:,1);
            u       = y*v/sqrt(s(1));
        else
            [u s u] = svd(y*y');
            s       = diag(s);
            u       = u(:,1);
            v       = y'*u/sqrt(s(1));
        end
        d       = sign(sum(v));
        u       = u*d;
        v       = v*d;
        Y       = u*sqrt(s(1)/n);

        % --- set in structure ----
        xY.y    = y;                             % time series nos voxeis selecionados
        xY.yy   = transpose(mean(transpose(y))); % average (not in spm_regions)
        xY.u    = Y;                             % eigenvariate
        xY.v    = v;
        xY.s    = s;       
        
        % --- Define DCM struct for analysis ---
        DCM.xY(roi) = xY;
    end
    
    % build Inputs stimulus function in bins
    binsInt = intervalsPRT.(cond_names{cond})*16;  % inputs intervals of conditions in bins
    sfm = zeros(16*length(intervals),1);           % inputs or stimulus function matrix (sfm) in bins
    for i=1:size(binsInt,1)
        sfm(binsInt(i,1):binsInt(i,2)) = 1;
    end
    sfm = repmat(sfm,size(subjectslist,2),1);
    Sess.U(cond).u = sparse(sfm);
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC CAUSAL MODELLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear DCM

% SPECIFICATION DCMs "attentional modulation of backward/forward connection"
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

% Expects
%--------------------------------------------------------------------------
% DCM.a                              % switch on endogenous connections
% DCM.b                              % switch on bilinear modulations
% DCM.c                              % switch on exogenous connections
% DCM.d                              % switch on nonlinear modulations
% DCM.U                              % exogenous inputs
% DCM.Y.y                            % responses
% DCM.Y.X0                           % confounds
% DCM.Y.Q                            % array of precision components
% DCM.n                              % number of regions
% DCM.v                              % number of scans
% load(fullfile(data_path,'GLM','SPM.mat'));

% Options
%--------------------------------------------------------------------------
% DCM.options.two_state              % two regional populations (E and I)
% DCM.options.stochastic             % fluctuations on hidden states
% DCM.options.centre                 % mean-centre inputs
% DCM.options.nonlinear              % interactions among hidden states
% DCM.options.nograph                % graphical display
% DCM.options.induced                % switch for CSD data features
% DCM.options.P                      % starting estimates for parameters
% DCM.options.hidden                 % indices of hidden regions
% DCM.options.nmax                   % maximum number of (effective) nodes
% DCM.options.nN                     % maximum number of iterations


% Load regions of interest
%--------------------------------------------------------------------------
load(fullfile(data_path,'GLM','VOI_V1_1.mat'),'xY');
DCM.xY(1) = xY;
load(fullfile(data_path,'GLM','VOI_V5_1.mat'),'xY');
DCM.xY(2) = xY;
load(fullfile(data_path,'GLM','VOI_SPC_1.mat'),'xY');
DCM.xY(3) = xY;

DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

% Time series
%--------------------------------------------------------------------------
DCM.Y.dt  = datasetConfigs.TR;     % TR or repeat time
% DCM.Y.X0  = DCM.xY(1).X0;          % 
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

% Experimental inputs - escolher as condicoes a estudar
%--------------------------------------------------------------------------
DCM.U.dt   =  datasetConfigs.TR/16000;       % time bin length {seconds}
DCM.U.name = [cond_names'];                  % names of conditions
% DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...  % inputs or stimulus function matrix
%               SPM.Sess.U(2).u(33:end,1) ...
%               SPM.Sess.U(3).u(33:end,1)];
DCM.U.u = [];
for cond = 1:n_cond
    DCM.U.u    = [DCM.U.u Sess.U(cond).u];  % inputs or stimulus function matrix
end

% DCM parameters and options
%--------------------------------------------------------------------------
DCM.delays = repmat(DCM.Y.dt/2,DCM.n,1);
DCM.TE     = 0.03;          % echo time

DCM.options.nonlinear  = 0;  % Bilinear - '0', Nolinear - '1'
DCM.options.two_state  = 0;  % One State - '0', Two State - '1'
DCM.options.stochastic = 0;  % No - '0', Yes - '1'
DCM.options.nograph    = 1;

% Connectivity matrices for model with backward modulation
%--------------------------------------------------------------------------
% intrinsic (fixed) connection matrix (ANATOMIA)
DCM.a = eye(size(DCM.Y.y,2)) 
% for i = 1:size(DCM.Y.y,2)
DCM.a(find(SSGCResults.baseline~=0))=1;
% end
DCM.a = [1 1 0; 1 1 1; 0 1 1];                               
% input-dependent connection matrix (MODULACAO)
DCM.b = zeros(3,3,3);  DCM.b(2,1,2) = 1;  DCM.b(2,3,3) = 1;  
% input connection matrix (INFLUENCIA NA ANATOMIA)
DCM.c = [1 0 0; 0 0 0; 0 0 0];

DCM.d = zeros(3,3,0);                                        % 

save(fullfile(data_path,'GLM','DCM_mod_bwd.mat'),'DCM');

% Connectivity matrices for model with forward modulation
%--------------------------------------------------------------------------
DCM.b = zeros(3,3,3);  DCM.b(2,1,2) = 1;  DCM.b(2,1,3) = 1;

save(fullfile(data_path,'GLM','DCM_mod_fwd.mat'),'DCM');

% DCM Estimation
%--------------------------------------------------------------------------
clear matlabbatch

matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {...
    fullfile(data_path,'GLM','DCM_mod_bwd.mat'); ...
    fullfile(data_path,'GLM','DCM_mod_fwd.mat')};

spm_jobman('run',matlabbatch);

% Bayesian Model Comparison
%--------------------------------------------------------------------------
DCM_bwd = load('DCM_mod_bwd.mat','F');
DCM_fwd = load('DCM_mod_fwd.mat','F');
fprintf('Model evidence: %f (bwd) vs %f (fwd)\n',DCM_bwd.F,DCM_fwd.F);
