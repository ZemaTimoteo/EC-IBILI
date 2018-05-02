% =========================================================================
% Effective Connectiviy with SSGC and TRS for real fMRI data Script
% =========================================================================

close all
clc
clear

tic

set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"


%% ... import data ...

% data has to be structured as follows:
%   [field]                        [value]
%  - VOIS                       string in cell with the name of VOIS
%  - VOISGeoLoc                 Geolocalization of the coordintates of VOIS
%  - Name of Condition          Matrix (number of VOIS by Concatenated all subjects all runs all blocks) 
%  - [...] for all conditons     [...] equal to the upper line

% [data.filename,data.pathname,data.filterindex] = uigetfile({'*.*'},'Choose signal to analyse');
% path = fullfile(data.pathname , data.filename);
% data = [importdata(path)];
% meu computador
% data = load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\JP_EC_Data_&_Results\vois-data_vf.mat');
% computador IBILI
data = load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\JP_EC_Data_&_Results\vois-nosmoothing-data.mat');

%% ... Choose ROIS & GeoLocalization coordinates of ROIS...

ROIS_choice     = [1:4 9:10];
nodes           = size(ROIS_choice,2);
dat.VOISInterst = data.VOIS(ROIS_choice);
dat.VOISGeoLoc  = data.VOISGeoLoc(ROIS_choice,:);


%% ... data partitazion in conditions ...
fields ={'VOISGeoLoc','VOIS'};
data = rmfield(data,fields);
numberFields = length(fieldnames(data));
dat.nameCondit = fieldnames(data);

for i=1:numberFields
    cond =  ['cond',num2str(i)];
    dat.(cond) = data.(dat.nameCondit{i})(ROIS_choice,:);
end

clear data i cond numberFields ROIS_choice


%% ... Choose Conditions ...

prompt = {'Condition test A','Condition test B'};
title = 'Conditions to test';
dims = [1 35];
definput = {'1','2'};
answer   = inputdlg(prompt,title,dims,definput);

test.cond = {dat.nameCondit{str2num(answer{1})} , dat.nameCondit{str2num(answer{2})}};
test.Data = {dat.(['cond',num2str(answer{1})]) , dat.(['cond',num2str(answer{2})])};

clear definput dims title prompt

%% cicle for multiple shuffle
numberSHUF    = 200;                                 % number of Shuffles
initialiDiago = nan(nodes,1)';
aux = 1; signifResults{aux} = diag(initialiDiago);
aux = 2; signifResults{aux} = diag(initialiDiago);

clear initialiDiago aux

for pp=1:numberSHUF
    fprintf('\n\n\n   ===================-------------------================  ');
    fprintf('\n   ===================----shuffle %d-----================   ', pp);
    fprintf('\n   ===================-------------------================   \n\n');

    %% ... reorganization of the data ...
    DataTest{1} = test.Data{1};
    DataTest{2} = test.Data{2};
    
    if length(DataTest{1}) == 960
        % Data condition 1
        % - randomize subjects
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),96,10]);
        [~,~,s] = size(DataTest{1});
        test.shuffl{pp}  = randperm(s) ;
        DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),960,1]);
        % - rearrangement of data
        test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),48,20]);
        
        % Data condition 2
        % - randomize subjects
        DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),96,10]);
        DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
        DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),960,1]);
        % rearrangement of data
        test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),48,20]);
        
    elseif length(DataTest{1}) == 2880
        % Data condition 1
        % - randomize subjects
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),288,10]);
        [~,~,s] = size(DataTest{1});
        test.shuffl{pp}  = randperm(s) ;
        DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),2880,1]);
        
        % - rearrangement of data
        a = DataTest{1}(:,1:1152);
        b = DataTest{1}(:,1152+1:1152+864);
        c = DataTest{1}(:,1152+864+1:1152+864+864);
        
        test.rshapData{1,1} = reshape(a,[size(DataTest{1},1),48,24]);
        test.rshapData{1,2} = reshape(b,[size(DataTest{1},1),36,24]);
        test.rshapData{1,3} = reshape(c,[size(DataTest{1},1),36,24]);
        
        % Data condition 2
        % - randomize subjects
        DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),288,10]);
        DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
        DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),2880,1]);
        
        
        % - rearrangement of data
        aa = DataTest{2}(:,1:1152);
        bb = DataTest{2}(:,1152+1:1152+864);
        cc = DataTest{2}(:,1152+864+1:1152+864+864);
        
        test.rshapData{2,1} = reshape(aa,[size(DataTest{2},1),48,24]);
        test.rshapData{2,2} = reshape(bb,[size(DataTest{2},1),36,24]);
        test.rshapData{2,3} = reshape(cc,[size(DataTest{2},1),36,24]);
    end
    
    clear a aa b bb c cc s ss
%     test = rmfield(test,{'Data'});
    
    
    %% ... Performing SSGC ...
    
    fprintf('\n\n\n   ================   SSGC  TESTING   ================   \n\n');
    numberTests = size(test.rshapData,2);
    mordermax   = 1;
    failures = 0;
    
    for i=1:numberTests
        
        fprintf('\n\n   * Performing SSGC for Condition %s partition %d ... \n',test.cond{1},i);
        aux = 1;                              % defines test for condition 1
        BOLD_data = test.rshapData{aux,i};    % data identification
        SSGC_script                           % perform SSGC
        clear BOLD_data aux
        
        fprintf('\n\n   * Performing SSGC for Condition %s partition %d ... \n',test.cond{2},i);
        aux = 2;                             % defines test for condition 2
        BOLD_data = test.rshapData{aux,i};   % data identification
        SSGC_script                          % perform SSGC
        clear BOLD_data aux
        
    end
    
    
    %% ... Performing TRS for data ...
    
    fprintf('\n\n\n   ============ PERFORMING SURROGATE - TRS ===========   \n');
    failures_sur = 0;
    
    for u=1:numberTests
        
        fprintf('\n\n   * Surrogate TRS for Condition %s partition %d ... \n',test.cond{1},u);
        aux = 1;                              % defines test for condition 1
        surrogate_function_vfMRI              % perform TRS
        
        clear BOLD_data aux
        
        fprintf('\n\n   * Surrogate TRS for Condition %s partition %d ... \n',test.cond{2},u);
        aux = 2;                             % defines test for condition 2
        surrogate_function_vfMRI             % perform TRS
        clear BOLD_data aux
        
    end
    
    
    %% ... Statistical validation of SSGC vs TRS with t.test ...
    
    fprintf('\n\n\n   ========= SURROGATE STATISTICAL TESTING =========   \n\n');
    failures_sur = 0;
    
    for v=1:numberTests
        
        fprintf('\n\n   * Performing Statistic for Condition %s partition %d ... \n',test.cond{1},v);
        aux = 1;                              % defines test for condition 1
        statistical_ssgc_trs_function_vfMRI   % perform TRS

        clear aux 
        
        fprintf('\n\n   * Performing Statistic for Condition %s partition %d ... \n',test.cond{2},v);
        aux = 2;                             % defines test for condition 2
        statistical_ssgc_trs_function_vfMRI  % perform TRS
       
        clear aux
        
    end
    dat.mordermax        = mordermax;         % max model order
    Multiple_Zscores{pp} = Results.z_scores;  % Storage of Zscores
    
    clear moAIC moBIC morder mordermax i u v cat_G_AllSubs cat_G_AllSubs_sur ci alpha1 sigma stats z ...
        z_scor prob p mu j infinit h ci blo ans diagonal cval
    
    
%% ... Distribution of significancy ...  
    for ii=1:numberTests
        aux = 1;
        conectIndex                     = find(Multiple_Zscores{1,pp}{aux,ii}>0);
        signifResults{aux}(conectIndex) = signifResults{aux}(conectIndex)+Multiple_Zscores{1,pp}{aux,ii}(conectIndex); % add significative values to matrix of condition1
        
        aux = 2;
        conectIndex                     = find(Multiple_Zscores{1,pp}{aux,ii}>0);
        signifResults{aux}(conectIndex) = signifResults{aux}(conectIndex)+Multiple_Zscores{1,pp}{aux,ii}(conectIndex); % add significative values to matrix of condition2
    end
    
    clear aux conectIndex
end

%% ... Plots ...

% cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI')
set(0,'DefaultFigureVisible','on');  % all subsequent figures "off"
plot_testMultiShuffle


%% ... Head ROIS plot ...
% condition 1
aux = 1;
drawConnect_testMultiShuffle

% condition 2
aux = 2;
drawConnect_testMultiShuffle


%% ... Save ...

% mordermax   = 1;             % model order max
% 
% if test == 1; testType = 'MVGC_'; elseif test == 3; testType = 'SSGC_';end
% if varComplex == 1; complex = 'LowCorr_'; else complex = 'HighCorr_';end
% saveName = strcat(num2str(var),'var_',complex,testType,surrogate,'_SNR',...
%             num2str(ratioSNR),'a_ChannelGap',num2str(pointsVect(1)-50),...
%             '-',num2str(pointsVect(end)-50),'_Morder',...
%             num2str(mordermax),'.mat');
% cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix') % cd('C:\Users\admin\Documents\Tiago\IBILI\Results\Testes\MatrixResults_Test')
% save(saveName,'data')
% cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI') % cd('C:\Users\admin\Documents\Tiago\IBILI\Code\Test')

toc