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
% meu computador
% data = load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\JP_EC_Data_&_Results\vois-nosmoothing-data.mat');
data = load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\JP_EC_Data_&_Results\vois_shiftedtimecourse.mat');


%% ... Choose ROIS & GeoLocalization coordinates of ROIS...

ROIS_choice     = [1:10];
% ROIS_choice     = rois_aux;
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

%% ... Demean data ...
auxMean = [dat.cond1 dat.cond2 dat.cond3 dat.cond4 dat.cond5 dat.cond6];
meanAllCond = mean(auxMean,2);

dat.cond1 = dat.cond1 - meanAllCond;
dat.cond2 = dat.cond2 - meanAllCond;
dat.cond3 = dat.cond3 - meanAllCond;
dat.cond4 = dat.cond4 - meanAllCond;
dat.cond5 = dat.cond5 - meanAllCond;
dat.cond6 = dat.cond6 - meanAllCond;

clear meanAllCond auxMean

%% ... Choose Conditions ...

prompt = {'Condition test A','Condition test B'};
title = 'Conditions to test';
dims = [1 35];
definput = {'1','2'};
answer   = inputdlg(prompt,title,dims,definput);

test.cond = {dat.nameCondit{str2num(answer{1})} , dat.nameCondit{str2num(answer{2})}};
test.Data = {dat.(['cond',num2str(answer{1})]) , dat.(['cond',num2str(answer{2})])};

clear definput dims title prompt


%% preparation of cicle for multiple shuffle
numberSHUF = 100;   % number of Shuffles
typeOfTest = 'MT';  % ST - single trial analyses, MT - multiple trial analyses

initialiDiago = nan(nodes,1)';

if length(test.Data{1}) == 1600
    testnumber = 1;
elseif length(test.Data{1}) == 4800
    testnumber = 3;
end


for hh = 1:testnumber
    aux = 1;
    Z_scoreAllMat{aux,hh} = diag(initialiDiago);
    SignifBinary{aux,hh}  = diag(initialiDiago);
    % signifResults{aux} = diag(initialiDiago);
    aux = 2;
    Z_scoreAllMat{aux,hh} = diag(initialiDiago);
    SignifBinary{aux,hh}  = diag(initialiDiago);
    % signifResults{aux} = diag(initialiDiago);
end
clear initialiDiago aux

%% specification of binary approach

totalNodes = [1:nodes];
totalCombos = nchoosek(totalNodes,2);

AllMTResults.S{2,testnumber} = zeros(size(totalNodes,2),size(totalNodes,2));
AllMTResults.M{2,testnumber} = zeros(size(totalNodes,2),size(totalNodes,2));
AllMTResults.SignifBinary{2,testnumber} = zeros(size(totalNodes,2),size(totalNodes,2));
AllMTResults.Z_scoreAllMat{2,testnumber} = zeros(size(totalNodes,2),size(totalNodes,2));


%% cicle for permutation
for pp=1:numberSHUF
    
    fprintf('\n\n\n   ===================-------------------================  ');
    fprintf('\n   ===================--- shuffle %d -----================   ', pp);
    fprintf('\n   ===================-------------------================   \n\n');
    
    %% ... reorganization of the data ...
    DataTest{1} = test.Data{1};
    DataTest{2} = test.Data{2};
    
    % --- following paper - Fernandes, 2018
    switch typeOfTest
        case 'MT'
            
            % -------------- cond 1, 2, 3 e 4 ----------
            if length(DataTest{1}) == 960
                % Data condition 1
                % - randomize subjects
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),12,80]);
                [~,~,s] = size(DataTest{1});
                test.shuffl{pp}  = randperm(s) ;
                DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),960,1]);
                % - rearrangement of data
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),48,20]);
                
                % Data condition 2
                % - randomize subjects
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),12,80]);
                DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),960,1]);
                % rearrangement of data
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),48,20]);
                
            elseif length(DataTest{1}) == 1600
                % Data condition 1
                % - randomize subjects
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),20,80]);
                [~,~,s] = size(DataTest{1});
                test.shuffl{pp}  = randperm(s) ;
                DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),1600,1]);
                % - rearrangement of data
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),80,20]);
                
                % Data condition 2
                % - randomize subjects
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),20,80]);
                DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),1600,1]);
                % rearrangement of data
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),80,20]);
                
                % ---------------- cond 5-6 ---------------------------
            elseif length(DataTest{1}) == 2880
                % Data condition 1
                % - randomize subjects
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),12,240]);
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
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),12,240]);
                DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),2880,1]);
                
                
                % - rearrangement of data
                aa = DataTest{2}(:,1:1152);
                bb = DataTest{2}(:,1152+1:1152+864);
                cc = DataTest{2}(:,1152+864+1:1152+864+864);
                
                test.rshapData{2,1} = reshape(aa,[size(DataTest{2},1),48,24]);
                test.rshapData{2,2} = reshape(bb,[size(DataTest{2},1),36,24]);
                test.rshapData{2,3} = reshape(cc,[size(DataTest{2},1),36,24]);
                
            elseif length(DataTest{1}) == 4800
                % Data condition 1
                % - randomize subjects
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),20,240]);
                [~,~,s] = size(DataTest{1});
                test.shuffl{pp}  = randperm(s) ;
                DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
                DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),4800,1]);
                
                % - rearrangement of data
                a = DataTest{1}(:,1:1600*1);
                b = DataTest{1}(:,1600*1+1:1600*2);
                c = DataTest{1}(:,1600*2+1:1600*3);
                
                %                 a = DataTest{1}(:,1:1920);
                %                 b = DataTest{1}(:,1920+1:1920+1440);
                %                 c = DataTest{1}(:,1920+1440+1:1920+1440+1440);
                %
                test.rshapData{1,1} = reshape(a,[size(DataTest{1},1),80,20]);
                test.rshapData{1,2} = reshape(b,[size(DataTest{1},1),80,20]);
                test.rshapData{1,3} = reshape(c,[size(DataTest{1},1),80,20]);
                
                %                 test.rshapData{1,1} = reshape(a,[size(DataTest{1},1),80,24]);
                %                 test.rshapData{1,2} = reshape(b,[size(DataTest{1},1),60,24]);
                %                 test.rshapData{1,3} = reshape(c,[size(DataTest{1},1),60,24]);
                %
                % Data condition 2
                % - randomize subjects
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),20,240]);
                DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
                DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),4800,1]);
                
                
                % - rearrangement of data
                aa = DataTest{2}(:,1:1600*1);
                bb = DataTest{2}(:,1600*1+1:1600*2);
                cc = DataTest{2}(:,1600*2+1:1600*3);
                
                %                 aa = DataTest{2}(:,1:1920);
                %                 bb = DataTest{2}(:,1920+1:1920+1440);
                %                 cc = DataTest{2}(:,1920+1440+1:1920+1440+1440);
                
                test.rshapData{2,1} = reshape(aa,[size(DataTest{1},1),80,20]);
                test.rshapData{2,2} = reshape(bb,[size(DataTest{1},1),80,20]);
                test.rshapData{2,3} = reshape(cc,[size(DataTest{1},1),80,20]);
                
                %                 test.rshapData{2,1} = reshape(aa,[size(DataTest{2},1),80,24]);
                %                 test.rshapData{2,2} = reshape(bb,[size(DataTest{2},1),80,24]);
                %                 test.rshapData{2,3} = reshape(cc,[size(DataTest{2},1),60,24]);
            end
            
            clear a aa b bb c cc s ss
            %     test = rmfield(test,{'Data'});
            
            % --- follwing simple subject trial
        case'ST'
            
            % -------------- cond 1, 2, 3 e 4 ----------
            if length(DataTest{1}) == 960
                % Data condition 1 & 2
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),12,80]);
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),12,80]);
                
            elseif length(DataTest{1}) == 1600
                % Data condition 1 & 2
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),20,80]);
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),20,80]);
                
                % -------------- cond 5 e 6 ---------------
            elseif length(DataTest{1}) == 2880
                % Data condition 1 % 2
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),12,240]);
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),12,240]);
                
            elseif length(DataTest{1}) == 4800
                % Data condition 1 % 2
                test.rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),20,240]);
                test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),20,240]);
                
            end
    end
    
    auxrshapData = test.rshapData;
    test = rmfield(test,'rshapData');

    for tt=1:size(totalCombos)
        % selecting channels
        rois_aux = totalCombos(tt,:); % rois of interest
        test.rshapData{1,1} = auxrshapData{1,1}(rois_aux,:,:);
        test.rshapData{2,1} = auxrshapData{2,1}(rois_aux,:,:);
        
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
        
        dat.mordermax        = mordermax;           % max model order
        
        Multiple_Zscores{pp} = Results.z_scores;    % Storage of Zscores
        BinarySignif{pp}     = Results.HighCritVal; % Storage of Binary significant z_scores
        
        clear moAIC moBIC morder mordermax i u v cat_G_AllSubs cat_G_AllSubs_sur ci alpha1 sigma stats z ...
            z_scor prob p mu j infinit h ci blo ans diagonal cval
        
        
        %% ... Distribution of significancy ...
        for ii=1:numberTests
            aux = 1;
            Z_scoreAllMat{aux,ii}(rois_aux(1),rois_aux(2)) = Z_scoreAllMat{aux,ii}(rois_aux(1),rois_aux(2)) + Multiple_Zscores{1,pp}{aux,ii}(1,2); % add z_score values to matrix of condition1
            connectBinary = find(BinarySignif{1,pp}{aux,ii}(1,2)>0);
            if  isempty(connectBinary) == 1
                connectBinary = 0;
            end
            SignifBinary{aux,ii}(rois_aux(1),rois_aux(2))  = SignifBinary{aux,ii}(rois_aux(1),rois_aux(2)) + connectBinary; % add z_score values to matrix of condition1
            Z_scoreAllMat{aux,ii}(rois_aux(2),rois_aux(1)) = Z_scoreAllMat{aux,ii}(rois_aux(2),rois_aux(1)) + Multiple_Zscores{1,pp}{aux,ii}(2,1); % add z_score values to matrix of condition1
            connectBinary = find(BinarySignif{1,pp}{aux,ii}(2,1)>0);
            if  isempty(connectBinary) == 1
                connectBinary = 0;
            end
            SignifBinary{aux,ii}(rois_aux(2),rois_aux(1))  = SignifBinary{aux,ii}(rois_aux(2),rois_aux(1)) + connectBinary; % add z_score values to matrix of condition1
            
            aux = 2;
            Z_scoreAllMat{aux,ii}(rois_aux(1),rois_aux(2)) = Z_scoreAllMat{aux,ii}(rois_aux(1),rois_aux(2)) + Multiple_Zscores{1,pp}{aux,ii}(1,2); % add z_score values to matrix of condition1
            connectBinary = find(BinarySignif{1,pp}{aux,ii}(1,2)>0);
            if  isempty(connectBinary) == 1
                connectBinary = 0;
            end
            SignifBinary{aux,ii}(rois_aux(1),rois_aux(2))  = SignifBinary{aux,ii}(rois_aux(1),rois_aux(2)) + connectBinary; % add z_score values to matrix of condition1
            Z_scoreAllMat{aux,ii}(rois_aux(2),rois_aux(1)) = Z_scoreAllMat{aux,ii}(rois_aux(2),rois_aux(1)) + Multiple_Zscores{1,pp}{aux,ii}(2,1); % add z_score values to matrix of condition1
            connectBinary = find(BinarySignif{1,pp}{aux,ii}(2,1)>0);
            if  isempty(connectBinary) == 1
                connectBinary = 0;
            end
            SignifBinary{aux,ii}(rois_aux(2),rois_aux(1))  = SignifBinary{aux,ii}(rois_aux(2),rois_aux(1)) + connectBinary; % add z_score values to matrix of condition1            
        end
        
        clear aux conectIndex
    end
end

%% ... statistical measures of Permutations (Standart Deviantion & Mean) ...

% save sum variables
MTResults.Z_scoreAllMat = Z_scoreAllMat;
MTResults.SignifBinary = SignifBinary;
    
for ii=1:numberTests
    % condition 1
    aux = 1;
    auxMult_Zscores = zeros(nodes,nodes,pp);
    for j=1:pp
        auxMult_Zscores(:,:,j) = Multiple_Zscores{1,j}{aux,ii};
    end
    MTResults.S{aux,ii} = std(auxMult_Zscores,0,3);          % standart deviation
    MTResults.M{aux,ii} = MTResults.Z_scoreAllMat{aux,ii}/pp; % mean
    clear auxMultZscores

    % condition 2
    aux = 2;
    auxMult_Zscores = zeros(nodes,nodes,pp);
    for j=1:pp
        auxMult_Zscores(:,:,j) = Multiple_Zscores{1,j}{aux,ii};
    end
    MTResults.S{aux,ii} = std(auxMult_Zscores,0,3);          % standart deviation
    MTResults.M{aux,ii} = MTResults.Z_scoreAllMat{aux,ii}/pp; % mean

    clear auxMultZscores
end


% clear SignifBinary Z_scoreAllMat BinarySignif Multiple_Zscores

%% ... Plots ...

% cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI')
% for iii=1:numberTests    
%     set(0,'DefaultFigureVisible','on');  % all subsequent figures "off"
%     % plot_testMultiShuffle
%     signifResults = MTResults.Z_scoreAllMat;
%     plot_testMultiShuffle_v2
%     signifResults = MTResults.SignifBinary;
%     plot_testMultiShuffle_v2
%     signifResults = MTResults.S;
%     plot_testMultiShuffle_v2
%     signifResults = MTResults.M;
%     plot_testMultiShuffle_v2
%     
% end
%% ... Head ROIS plot ...
% % % % % condition 1
% % % % aux = 1;
% % % % drawConnect_testMultiShuffle
% % % % 
% % % % % condition 2
% % % % aux = 2;
% % % % drawConnect_testMultiShuffle
% % % % 

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