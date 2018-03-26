  
%% =============================================================================================================

% ---- Statistical Validation over the SSGC data by TRS - using t-test ------------


figure()
for cond=1:n_cond
    % Plot of the F1 score
    subplot(1,3,cond)
    imagesc (Fone(:,:,cond))
    colorbar
    title(cond_names{cond})
    axis('square');
    ylabel('Number of points'); xlabel('Number of blocks');
    set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
    set(gca,'XTickLabel',[' 0';'  ';' 2';'  ';' 4';'  ';' 8';'  ';'13';'  ';'15';'  ';'17';'  ';'25';'  ';'30';'  ';'40'])
    set(gca,'YTickLabel',['  0';'   ';' 10';'   ';' 20';'   ';' 30';'   ';' 40';'   ';' 60';'   ';' 80';'   ';'100';'   ';'160';'   ';'200'])
    ax = gca; ax.XTickLabelRotation = 45;
    caxis([-.05 1.05]);
end
titleP = ['F1 score test'];
p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',-.32);
set(p.th,'edgecolor',.5*[1 1 1]);