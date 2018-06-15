%% Analyse the SSGC + TRS metrics in a bivariate way, covering all nodes7

clear all
clc
close all

data = load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\JP_EC_Data_&_Results\vois_shiftedtimecourse.mat');


%% Initialization

totalNodes = [1:10];
totalCombos = nchoosek(totalNodes,2);
blocos = 1;

AllMTResults.S{2,blocos} = zeros(size(totalNodes,2),size(totalNodes,2));
AllMTResults.M{2,blocos} = zeros(size(totalNodes,2),size(totalNodes,2));
AllMTResults.SignifBinary{2,blocos} = zeros(size(totalNodes,2),size(totalNodes,2));
% AllMTResults.Z_scoreAllMat{2,blocos} = zeros(size(totalNodes,2),size(totalNodes,2));

for tt=1:size(totalCombos)
%     % choose ROISs
%     rois_aux = totalCombos(tt,:);
%     % run teste multishuffle
%     testeMultiShuffle
    
    % Save data in structure
    for jj=1:numberTests
        aux = 1;
        AllMTResults.S{aux,jj}(rois_aux(1),rois_aux(2)) = MTResults.S{aux,jj}(1,2);
        AllMTResults.M{aux,jj}(rois_aux(1),rois_aux(2)) = MTResults.M{aux,jj}(1,2);
        AllMTResults.SignifBinary{aux,jj}(rois_aux(1),rois_aux(2))  =  MTResults.SignifBinary{aux,jj}(1,2);
        AllMTResults.Z_scoreAllMat{aux,jj}(rois_aux(1),rois_aux(2)) =  MTResults.Z_scoreAllMat{aux,jj}(1,2);
        
        AllMTResults.S{aux,jj}(rois_aux(2),rois_aux(1)) = MTResults.S{aux,jj}(2,1);
        AllMTResults.M{aux,jj}(rois_aux(2),rois_aux(1)) = MTResults.M{aux,jj}(2,1);
        AllMTResults.SignifBinary{aux,jj}(rois_aux(2),rois_aux(1))  =  MTResults.SignifBinary{aux,jj}(2,1);
        AllMTResults.Z_scoreAllMat{aux,jj}(rois_aux(2),rois_aux(1)) =  MTResults.Z_scoreAllMat{aux,jj}(2,1);
        
        aux= 2;
        AllMTResults.S{aux,jj}(rois_aux(1),rois_aux(2)) = MTResults.S{aux,jj}(1,2);
        AllMTResults.M{aux,jj}(rois_aux(1),rois_aux(2)) = MTResults.M{aux,jj}(1,2);
        AllMTResults.SignifBinary{aux,jj}(rois_aux(1),rois_aux(2))  =  MTResults.SignifBinary{aux,jj}(1,2);
        AllMTResults.Z_scoreAllMat{aux,jj}(rois_aux(1),rois_aux(2)) =  MTResults.Z_scoreAllMat{aux,jj}(1,2);
        
        AllMTResults.S{aux,jj}(rois_aux(2),rois_aux(1)) = MTResults.S{aux,jj}(2,1);
        AllMTResults.M{aux,jj}(rois_aux(2),rois_aux(1)) = MTResults.M{aux,jj}(2,1);
        AllMTResults.SignifBinary{aux,jj}(rois_aux(2),rois_aux(1))  =  MTResults.SignifBinary{aux,jj}(2,1);
        AllMTResults.Z_scoreAllMat{aux,jj}(rois_aux(2),rois_aux(1)) =  MTResults.Z_scoreAllMat{aux,jj}(2,1);
    end
end

clear MTResults aux rois_aux

