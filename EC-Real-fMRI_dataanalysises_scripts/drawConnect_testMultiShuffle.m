%% =============================================================================================================

% This scrips builds graph theory maps of the connectivity obtained in the
% metric.

%% Import Data


% ... DATA
dataConn = signifResults{1,aux};    % z_score matrix
% dataConn(find(dataConn>0))=1;     % uniformização
nnrois   = nodes;                   % number of ROIs

% ... Coordenates and region tag
regionTag = dat.VOISInterst;        % ROIs original names
coords    = dat.VOISGeoLoc;         % coordenates per region

% roiSize (max - 200)
roiSize = ones(nodes,1)*5;


% connStrength construcao da matrix que contem inform. dos canais e significancia
counter = 1;

for i = 1: nodes
    for j = 1:nodes
        if (dataConn(i,j) > 1 )
            connStrength(counter,:) = [j i dataConn(i,j)] ;
            counter = counter +1 ;
        end
    end
end

% define only the first nn connections - we can comment this if all are desired
nn = 6;
sortValue          = sortrows(connStrength,3,'descend');
sortValue(nn+1:end,:) = [];
connStrength = sortValue;

showLabels = true;
exportBrain = false;

%% Number of connections in and out flow

from_count = hist(connStrength(:,1),1:nodes);
to_count   = hist(connStrength(:,2),1:nodes);

% norm count
from_count = round(from_count/norm(from_count) * 40);
to_count   = round(to_count/norm(to_count) * 40);
both_count =  from_count + to_count;

% % % from brain
% % brainConnectROIs_vfMRI(coords((from_count~=0),:), regionTag(from_count~=0), from_count(from_count~=0), connStrength, showLabels, 1, 'from' );
% % 
% % % to brain
% % brainConnectROIs_vfMRI(coords((to_count~=0),:), regionTag(to_count~=0), to_count(to_count~=0), connStrength, showLabels, 1, 'to' );

% Causal connections
% brainConnectROIs_vfMRI(coords,regionTag, both_count(both_count~=0), connStrength, showLabels, 1, 'connectivity all subs all conds' );
brainConnectROIs_vfMRI(coords,regionTag, roiSize, connStrength, showLabels, 1, 'connectivity all subs all conds' );


%% Graph Theory

% % % % for j = 1: size(C,1)
% % % %         F_AllSubs_AllConds_pw_FDR(C(j,1),C(j,2)) = ResultsMVGC_pairs_AllCond.(['C' num2str(C(j,:),'%i%i')]).F(1,2);
% % % %         F_AllSubs_AllConds_pw_FDR(C(j,2),C(j,1)) = ResultsMVGC_pairs_AllCond.(['C' num2str(C(j,:),'%i%i')]).F(2,1);
% % % % end

% % % % for i=1:nnrois
% % % %     for j=1:nnrois
% % % %         % Connectivity Inflow.Outflow Relationships
% % % %         aux      = find(connStrength(:,1)==i & connStrength(:,3)==1);
% % % %         aux_II   = sum(F_AllSubs_AllConds_pw_FDR(connStrength(aux,1:2)));
% % % %         k_out(i) = aux_II(1,2)/nnrois;
% % % %         aux      = find(connStrength(:,1)==i & connStrength(:,3)==-1);
% % % %         aux_II   = sum(F_AllSubs_AllConds_pw_FDR(connStrength(aux,1:2)));
% % % %         k_in(i)  = aux_II(1,2)/nnrois;
% % % %         k_tot(i) = k_out(i) - k_in(i);
% % % %         
% % % %         % Clustering Coefficie
% % % % %         denom = (k_tot(i)*(k_tot(i)-1)-sum
% % % % %         Clust (i) = /nnrois;
% % % %         
% % % %         % Modulatory
% % % %         
% % % %         
% % % %         % Betweenness Centrality
% % % %     end
% % % % end