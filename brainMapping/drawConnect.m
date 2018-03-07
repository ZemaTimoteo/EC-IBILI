% addpath('C:\Users\Admin\Documents\Projects\brain connectivity\connectivity\MVGC_Testing\BrainConnectionsDesign')


%% Import Data


% DATA

% dataConn = sig_AllSubs_AllConds_pw_FDR;
% dataConn = sig_AllSubs_perConds_pw_FDR.(cond_names{1});
% dataConn = sig_AllSubs_AllConds_pw_FDR;

dataConn = ResultsMVGC_AllCond_AllSub.sig;

nnrois = length(voiFile.L.DAT.VOI);

% COORDS and region tag

for i = 1:nnrois, 
    coords(i,:) = mean( voiFile.L.DAT.VOI(i).Voxels );
    
    tag_temp = voiFile.L.DAT.VOI(i).Name; 
    tag = strrep(tag_temp, '_', ' ');
    regionTag{i} = tag;
    
    
end

% roiSize (max - 200)
roiSize = ones(nnrois,1);


counter = 1;
% connStrength construcao da matrix que contem inform. dos canais e significancia
for i = 1: nnrois
    for j = 1:nnrois
        if (dataConn(i,j) == 1 )
            connStrength(counter,:) = [j i dataConn(i,j)] ;
            counter = counter +1 ;
        end
    end
end

showLabels = true;
exportBrain = false;

%% Number of connections in and out flow

from_count = hist(connStrength(:,1),1:nnrois);
to_count   = hist(connStrength(:,2),1:nnrois);

% norm count
from_count = round(from_count/norm(from_count) * 80);
to_count   = round(to_count/norm(to_count) * 80);

% from brain
brainConnectROIs(coords((from_count~=0),:), regionTag(from_count~=0), from_count(from_count~=0), connStrength, showLabels, 1, 'from' );


brainConnectROIs(coords((to_count~=0),:), regionTag(to_count~=0), to_count(to_count~=0), connStrength, showLabels, 1, 'to' );




%% Causal connections
brainConnect(coords, regionTag, roiSize, connStrength, showLabels,'connectivity all subs all conds', exportBrain );

%% Graph Theory

% % % % for j = 1: size(C,1)
% % % %         F_AllSubs_AllConds_pw_FDR(C(j,1),C(j,2)) = ResultsMVGC_pairs_AllCond.(['C' num2str(C(j,:),'%i%i')]).F(1,2);
% % % %         F_AllSubs_AllConds_pw_FDR(C(j,2),C(j,1)) = ResultsMVGC_pairs_AllCond.(['C' num2str(C(j,:),'%i%i')]).F(2,1);
% % % % end

for i=1:nnrois
    for j=1:nnrois
        % Connectivity Inflow.Outflow Relationships
        aux      = find(connStrength(:,1)==i & connStrength(:,3)==1);
        aux_II   = sum(F_AllSubs_AllConds_pw_FDR(connStrength(aux,1:2)));
        k_out(i) = aux_II(1,2)/nnrois;
        aux      = find(connStrength(:,1)==i & connStrength(:,3)==-1);
        aux_II   = sum(F_AllSubs_AllConds_pw_FDR(connStrength(aux,1:2)));
        k_in(i)  = aux_II(1,2)/nnrois;
        k_tot(i) = k_out(i) - k_in(i);
        
        % Clustering Coefficie
%         denom = (k_tot(i)*(k_tot(i)-1)-sum
%         Clust (i) = /nnrois;
        
        % Modulatory
        
        
        % Betweenness Centrality
    end
end