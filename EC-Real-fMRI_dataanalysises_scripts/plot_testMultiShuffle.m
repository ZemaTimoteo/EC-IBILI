  
%% =============================================================================================================

% ---- PLOT OF THE multiShuffle Statistical Validation over the SSGC data by TRS - using t-test ------------


figure()

% ... X & Y ticklabel
auxNames = char(dat.VOISInterst);
for j=1:nodes
    nodesNames(j*2+1,:) = auxNames(j,:);
end
% nodesNames(2,:) = '    ';


% define max value for axis label
ConcatSignif = [signifResults{1,1}; signifResults{1,2}];
maxValue = max(max(ConcatSignif));

% ... plot the z_scores
for o=1:2
    subplot(1,2,o)
    imagesc (signifResults{1,o})
    colorbar
    
    title(strcat('Conditon: . ',dat.nameCondit{str2num(answer{o})}))
    axis('square');
    ylabel('VOIs'); xlabel('VOIs');
    set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
    set(gca,'XTickLabel',nodesNames)
    set(gca,'YTickLabel',nodesNames)
    ax = gca; ax.XTickLabelRotation = 45; %['   ';'';'   ';'';'   ';'';'   ';'';'   ']
    caxis([-.05 maxValue]);
end


% titleP = strcat('Statistical Validation of Connectivity for Condition: ',dat.nameCondit{str2num(answer{1})},' and Condition: .',dat.nameCondit{str2num(answer{2})});
% p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',-.12);
% set(p.th,'edgecolor',.5*[1 1 1]);

clear titleP j o nodesNames auxNames ax p aux