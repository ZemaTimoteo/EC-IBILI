%% ========================================================================
%  
%  This function allowed to plot the sensitivity and specificity as a
%  function of the parameters choosen by the user. This parameters are:
%       - Size of data;
%       - Number of channels;
%       - Number of trials;
%
%% ========================================================================

clc
option = 'b';
pointsVect  =  [10 20 30 40 60 80 100 160 200];

for cond = 1 : n_cond
    i = 0;
    clear condHigh
    %% ++++++++ Build data +++++++++++++++++++++++
    
    if consistTest
        cons = data.cons.(cond_names{cond});
        cons(find(data.cons.(cond_names{cond})==0)) = nan;
    end
    for po = 1 : size(pointsVect,2)
        for blo = 1 : size(blocksVect,2)
            sensi(po,blo,cond) = data.results.sensi{po,blo}.(cond_names{cond});
            speci(po,blo,cond) = data.results.speci{po,blo}.(cond_names{cond});
            Fone(po,blo,cond)  = data.results.Fone{po,blo}.(cond_names{cond});
            ACC(po,blo,cond)   = data.results.ACC{po,blo}.(cond_names{cond});            
%             bACC(po,blo,cond)  = data.results.bACC{po,blo}.(cond_names{cond});
            if sensi(po,blo,cond) > 0.55 && speci(po,blo,cond) > 0.55
                i=i+1;
                condHigh.value(i) = sensi(po,blo,cond)+speci(po,blo,cond);
                condHigh.po(i)    = po;
                condHigh.blo(i)   = blo;
            end
        end
        if consistTest
            consPLOT(cond,:,po) = reshape(cons(:,:,po),1,[]);
%             consPLOT(isnan(consPLOT(cond,:,po))) = [];
        end
    end
    
    if  exist('condHigh')
        pointsVect(condHigh.po);
        blocksVect(condHigh.blo);
        [M(cond),ind(cond)]= max(condHigh.value);
        fprintf('-> %d points & %d blocks are the Optimal parameters for condition %s,  \n\n',pointsVect(condHigh.po(ind(cond))),blocksVect(condHigh.blo(ind(cond))),(cond_names{cond}));
    end
    [pointMax,blockMax] = find(Fone(:,:,cond)==max(max(Fone(:,:,cond))));
    fprintf('-> %d points & %d blocks offer the highest F1 score %d for condition %s,  \n\n',pointsVect(pointMax(1)),blocksVect(blockMax(1)),max(max(Fone(:,:,cond))),(cond_names{cond}));
    
    %% ++++++++ Plot +++++++++++++++++++++++
    figure()
    if option == 'a'
        [X Y]     = meshgrid(pointsVect,blocksVect); 
        %     resol     = 10;
        %     ittPoints = (max(pointsVect) - min(pointsVect))/resol;
        %     ittBlocks = (max(blocksVect) - min(blocksVect))/resol;
        %     [Xq,Yq] = meshgrid(min(pointsVect):ittPoints:max(pointsVect),min(blocksVect):ittBlocks:max(blocksVect));
        
        % Plot of the sensitivity
        subplot(2,2,1)
        surf(Y,X,sensi(:,:,cond))
        colorbar
        title('Sensitivity test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');  
        zlim([-.05 1.05]);
        caxis([-.05 1.05]);
        
        % Plot of the specificity
        subplot(2,2,2)
        surf(Y,X,speci(:,:,cond))
        colorbar
        title('Specificity test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        zlim([-.05 1.05]);
        caxis([-.05 1.05]);
        
        % Plot of the F1 score
        subplot(2,2,3)
        %     Vq = interp2(X,Y,Fone(:,:,cond),Xq,Yq); % Interpolate at the query points.
        %     surf(Xq,Yq,Vq);
        surf(Y,X,Fone(:,:,cond))
        colorbar
        title('F1 score test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        zlim([-.05 1.05]);
        caxis([-.05 1.05]);
        
        % Plot of the Accuracy
        subplot(224)
        surf(Y,X,ACC(:,:,cond))
        colorbar
        title('Accuracy test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        zlim([-.05 1.05]);
        caxis([-.05 1.05]);

        
    elseif option == 'b'
        % Plot of the sensitivity
        subplot(2,2,1)
        imagesc (sensi(:,:,cond))
        colorbar
        title('Sensitivity test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])         
        set(gca,'XTickLabel',[' 0';'  ';' 2';'  ';' 4';'  ';' 8';'  ';'13';'  ';'15';'  ';'17';'  ';'25';'  ';'30';'  ';'40'])
        set(gca,'YTickLabel',['  0';'   ';' 10';'   ';' 20';'   ';' 30';'   ';' 40';'   ';' 60';'   ';' 80';'   ';'100';'   ';'160';'   ';'200'])        
        ax = gca; ax.XTickLabelRotation = 45;
        caxis([-.05 1.05]);
        
        % Plot of the specificity
        subplot(2,2,2)
        imagesc (speci(:,:,cond))
        colorbar
        title('Specificity test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
        set(gca,'XTickLabel',[' 0';'  ';' 2';'  ';' 4';'  ';' 8';'  ';'13';'  ';'15';'  ';'17';'  ';'25';'  ';'30';'  ';'40'])
        set(gca,'YTickLabel',['  0';'   ';' 10';'   ';' 20';'   ';' 30';'   ';' 40';'   ';' 60';'   ';' 80';'   ';'100';'   ';'160';'   ';'200'])        
        ax = gca; ax.XTickLabelRotation = 45;
        caxis([-.05 1.05]);
        
        % Plot of the F1 score
        subplot(2,2,3)
        imagesc (Fone(:,:,cond))
        colorbar
        title('F1 score test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
        set(gca,'XTickLabel',[' 0';'  ';' 2';'  ';' 4';'  ';' 8';'  ';'13';'  ';'15';'  ';'17';'  ';'25';'  ';'30';'  ';'40'])
        set(gca,'YTickLabel',['  0';'   ';' 10';'   ';' 20';'   ';' 30';'   ';' 40';'   ';' 60';'   ';' 80';'   ';'100';'   ';'160';'   ';'200'])        
        ax = gca; ax.XTickLabelRotation = 45;
        caxis([-.05 1.05]);
        
        % Plot of the Accuracy
        subplot(2,2,4)        
        imagesc (ACC(:,:,cond))
        colorbar
        title('Balanced Accuracy test')
        axis('square');
        ylabel('Number of points'); xlabel('Number of blocks');
        set(gca,'YTick',[0:0.5:9]);set(gca,'XTick',[0:0.5:9])
        set(gca,'XTickLabel',[' 0';'  ';' 2';'  ';' 4';'  ';' 8';'  ';'13';'  ';'15';'  ';'17';'  ';'25';'  ';'30';'  ';'40'])
        set(gca,'YTickLabel',['  0';'   ';' 10';'   ';' 20';'   ';' 30';'   ';' 40';'   ';' 60';'   ';' 80';'   ';'100';'   ';'160';'   ';'200'])        
        ax = gca; ax.XTickLabelRotation = 45;
        caxis([-.05 1.05]);
        
    end
    
    % Legend for the figure
    titleP = strjoin(cond_names(cond));
    titleP = [titleP ' - SNR = ' num2str(ratioSNR)];
    p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',.012);
    set(p.th,'edgecolor',.5*[1 1 1]);
    
    
end


if consistTest
    i = 1;
    figure()
    for po=1:size(pointsVect,2)
        if po == 2 || po == 5 || po == 8;
            subplot(1,3,i)
            thresh = 0.8;
            boxplot([consPLOT(1,:,po)',consPLOT(2,:,po)',consPLOT(3,:,po)'], ...
                'Notch','off','Labels',{'simul05','simul250','simulOrig'})
            ylim([-.05 1.05]);
            line([0 20],[thresh thresh],'Color','g');         % Legend for the figure
            %
            %         [consPLOT(1,:,po),consPLOT(2,:,po),consPLOT(3,:,po)],'Notch','on', ...
            %             'Labels',{'simul05','simul250','simulOrig'}
            i= i+1;
        end
    end
    
    titleP = ['Consistency: 20 s, 60 s & 160 s of data - SNR = ' num2str(ratioSNR)];
    p=mtit(titleP,'fontsize',30,'color',[0 0 0],'xoff',0.05,'yoff',.012);
    set(p.th,'edgecolor',.5*[1 1 1]);
end

if F1scorePlot
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
end
