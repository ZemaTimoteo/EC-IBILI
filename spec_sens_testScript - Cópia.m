% =========================================================================
%  Specificity and Sensitivity for each fs
% =========================================================================

close all
clc
clear
   
tic

set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

%% ... Inicializacao ...
pointsVect  = [10 20 30 40 60 80 100 160 200] + 50; % 60 80 100 160 200]+50; % + 50 para excluir os primeiro 50 pontos afectados pela HRF [3000]; % [10 20 30 40 60 80 100 160 200]; % [40 80]; % [40:5:80]; % 
blocksVect  = [ 2  4 8 13 15 17  25  30  40]; %  [ 20]; % [ 2  4  8 13 15 17  25  30  40]; % [4 25]; % 40; %
test        = 1;             % 1 - MVGC, 3 - SSGC, 4 - Coherence
surrogate   = 'perm';   % 'timeRev' - for time-reversed, 'perm' - for permutation, 'permBlock' - for permutation by blocks
var         = 5;             % Number of Variables - '5' or '9'
varComplex  = 15;             % LowCorr '1' - true; HighCorr '0' - false
simul       = 'true';

consistTest = 'true';       % if consistentency test on: 'true'
ratioSNR    = 1000;          % Value of SNR
plots       = 0;            % 1 - 'true', 0 ou 2 - 'false'
cond_names  = {'simulate05','simulate250','simulateOriginal'};
n_cond      = 3;
% mordermax   = 1;            % model order max

%%
point=1;
while point<=size(pointsVect,2) % itterates on the number of points
    block=1;
    clear BOLD_Original BOLD_250 BOLD_05

    while block<=size(blocksVect,2)  % itterates on the number of blocks
        mordermax   = 1;             % model order max
        fprintf('\n\n\n***** -> Number of points (test %d out of %d) - ',point,size(pointsVect,2));
        fprintf(' N-blocks (test %d out of %d) ***** \n\n',block,size(blocksVect,2));

        %% ... Geracao de Resultados ...
        VARDataGenerator
        surrogate_demo_function
        
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
        
        block = block + 1;
        
        clearvars -except pointsVect blocksVect point block test surrogate var ...
        simul data n_cond cond_names consistTest ratioSNR consist morder plots ...
        varComplex
    end
    point = point + 1;
    
end

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