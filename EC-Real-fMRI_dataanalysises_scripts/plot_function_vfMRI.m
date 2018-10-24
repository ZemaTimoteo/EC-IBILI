  
%% =============================================================================================================

% ---- PLOT OF THE Statistical Validation over the SSGC data by TRS - using t-test ------------


figure()

% ... X & Y ticklabel
auxNames = char(dat.VOISInterst);
for j=1:nodes
    nodesNames(j*2+1,:) = auxNames(j,:);
end
% nodesNames(2,:) = '    ';

% ... plot the z_scores
for o=1:numberTests
    subplot(1,numberTests,o)
    imagesc (Results.z_scores{aux,o})
    colorbar
    
    title(strcat('Data Partition-',string(o)))
    axis('square');
    ylabel('VOIs'); xlabel('VOIs');
    set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
    set(gca,'XTickLabel',nodesNames)
    set(gca,'YTickLabel',nodesNames)
    ax = gca; ax.XTickLabelRotation = 45; %['   ';'';'   ';'';'   ';'';'   ';'';'   ']
    caxis([-.05 2.55]);
end


titleP = strcat('Statistical Validation of Connectivity for Condition: . ',dat.nameCondit{str2num(answer{aux})});
p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',-.12);
set(p.th,'edgecolor',.5*[1 1 1]);

clear titleP j o nodesNames auxNames ax p aux