% =========================================================================
% Effective Connectiviy with SSGC and TRS for real fMRI data Script
% =========================================================================

close all
clc
clear
   
tic

set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

%% ... import data ...

[data.filename,data.pathname,data.filterindex] = uigetfile({'*.*'},'Choose signal to analyse');
path = fullfile(data.pathname , data.filename);
data = [importdata(path)];

%% ... data partitazion in conditions...

dat.VOIS = data.VOIS;
data = rmfield(data,'VOIS');
numberFields = length(fieldnames(data));
dat.nameCondit = fieldnames(data);

for i=1:numberFields
    cond =  ['cond',num2str(i)];
    dat.(cond) = data.(dat.nameCondit{i});
end

clear data i cond numberFields

%% ... Choose Conditions ...

prompt = {'Condition test A','Condition test B'};
title = 'Conditions to test';
dims = [1 35];
definput = {'1','2'};
answer   = inputdlg(prompt,title,dims,definput);

test.cond = {dat.nameCondit{str2num(answer{1})} , dat.nameCondit{str2num(answer{2})}};
test.Data = {dat.(['cond',num2str(answer{1})]) , dat.(['cond',num2str(answer{2})])};

clear answer definput dims title prompt

%% ... reorganization of the data ...

if length(test.Data{1}) == 960
    
   test111 = reshape(test.Data{1}, [10,32,30]);
    
    
elseif length(test.Data{1}) == 2880 
    % rearrangement for Data condition 1
    a = test.Data{1}(:,1:1152);
    b = test.Data{1}(:,1152+1:1152+864);
    c = test.Data{1}(:,1152+864+1:1152+864+864);

    test.rshapData{1,1} = reshape(a,[size(test.Data{1},1),48,24]);
    test.rshapData{1,2} = reshape(b,[size(test.Data{1},1),36,24]);
    test.rshapData{1,3} = reshape(c,[size(test.Data{1},1),36,24]);
    
    % rearrangement for Data condition 2
    aa = test.Data{2}(:,1:1152);
    bb = test.Data{2}(:,1152+1:1152+864);
    cc = test.Data{2}(:,1152+864+1:1152+864+864);

    test.rshapData{2,1} = reshape(aa,[size(test.Data{2},1),48,24]);
    test.rshapData{2,2} = reshape(bb,[size(test.Data{2},1),36,24]);
    test.rshapData{2,3} = reshape(cc,[size(test.Data{2},1),36,24]);
end

clear a aa b bb c cc
test = rmfield(test,{'Data'});

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
%% ... Performing TRS for SSGC obtained ...
fprintf('\n\n\n   ================ SURROGATE TESTING ================   \n\n');
failures_sur = 0;

for u=1:numberTests

    fprintf('\n\n   * Performing TRS for Condition %s partition %d ... \n',test.cond{1},i);    
    aux = 1;                              % defines test for condition 1
    BOLD_data = test.rshapData{aux,u};    % data identification
    surrogate_function_vfMRI              % perform TRS
                       
    clear BOLD_data aux

    fprintf('\n\n   * Performing SSGC for Condition %s partition %d ... \n',test.cond{2},i);
    aux = 2;                             % defines test for condition 2
    BOLD_data = test.rshapData{aux,u};   % data identification
    surrogate_function_vfMRI             % perform TRS
    clear BOLD_data aux
    
end


%% ... Testar Resultados ...
% save data in new variables
data.z_scores{point,block} = z_scores;
data.values{point,block} = [pointsVect(point)-50 blocksVect(block)];

% calculate sensi & speci
for cond = 1:n_cond
    [data.results.sensi{point,block}.(cond_names{cond}), ...
        data.results.speci{point,block}.(cond_names{cond}), ...
        data.results.Fone{point,block}.(cond_names{cond}),  ...
        data.results.bACC{point,block}.(cond_names{cond}),  ...
        data.results.ACC{point,block}.(cond_names{cond})]  ...
        = sensi_speci_ECt(data.z_scores{point,block}.(cond_names{cond}), ...
        label_ec, ...
        nvars, ...
        'SSGC');
end


clearvars -except pointsVect blocksVect point block test surrogate var ...
    simul data n_cond cond_names consistTest ratioSNR consist morder plots ...
    varComplex


data.cons = consist;  % data consistency


%% ... Save ...
mordermax   = 1;             % model order max

if test == 1; testType = 'MVGC_'; elseif test == 3; testType = 'SSGC_';end
if varComplex == 1; complex = 'LowCorr_'; else complex = 'HighCorr_';end
saveName = strcat(num2str(var),'var_',complex,testType,surrogate,'_SNR',...
            num2str(ratioSNR),'a_ChannelGap',num2str(pointsVect(1)-50),...
            '-',num2str(pointsVect(end)-50),'_Morder',...
            num2str(mordermax),'.mat');
cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix') % cd('C:\Users\admin\Documents\Tiago\IBILI\Results\Testes\MatrixResults_Test')
save(saveName,'data')
cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI') % cd('C:\Users\admin\Documents\Tiago\IBILI\Code\Test')

%% ... Plots ...
cd('C:\Users\tiago\Documents\MATLAB\EC-IBILI')
set(0,'DefaultFigureVisible','on');  % all subsequent figures "off"
F1scorePlot = 1;
% Plot sensi, speci, F1, ACC, bACC
plot_sensi_speci_demo_function

%%
toc