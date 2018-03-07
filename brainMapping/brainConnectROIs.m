function [hf] = brainConnectROIs(coord, regionTag, roiSize, connStrength, showLabels, printBrain, direction )

%% VarIn and presets

if nargin < 6
    printBrain = 0;
end

% Number of connections
numbOfConn = size(connStrength, 1);

switch (direction) 
    case 'to'
        colMap = colormap(copper(max(roiSize)));
    case 'from'
        colMap = colormap(summer(max(roiSize)));
end

connStrenghColor (connStrength(:,3) == 1,1) = 'r';  % outflow

for t = 1: size(connStrength,1)
    i = connStrength(t,1);
    j = connStrength(t,2);
    if ~isempty(find((connStrength(:,1)==j & connStrength(:,2)==i))==1)
        connStrenghColor (t) = 'g';  % bi-flow
      
    end
end
% % % connStrenghColor (connStrength(:,3) == -1,1) = 'f'; % inflow




%% Mesh Load
% Read the GIfTI surface file
g = gifti('cortex_20484.surf.gii');

%% Display mesh

%% BOTTOM

% hf = figure('units','normalized','outerposition',[0 0 1 1]); Fullscreen
figure('position',[0,0, 800,800]);


    h = plot(g);
    set(h,'EdgeColor','none','FaceColor',[0.15 0.15 0.2])
    alpha(0.15);
    view(0,-90)
%     colormap(colMap)
%     colorbar;
    
% Anterior/Posterior, Right/Left
if showLabels
    text(0, 85, 30,...
        sprintf( ' %s', 'Ant.'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    text(0, -115, 30,...
        sprintf( ' %s', 'Post.'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    
    text(85, 0, 0,...
        sprintf( ' %s', 'Right'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    text(-85, 0, 0,...
        sprintf( ' %s', 'Left'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
end

hold on


for i = 1 : length(roiSize)
    plot3(coord(i,1),coord(i,2),coord(i,3), '.',...
        'Color',colMap((roiSize(i)),:),...
        'MarkerSize',roiSize(i)*10);
    
    if showLabels
        text(coord(i,1),coord(i,2),coord(i,3),...
            sprintf( ' %s', regionTag{i}),...
            'Color','black',...
            'BackgroundColor',[0.9 0.9 0.9] ,...
            'EdgeColor',[0.5 0.5 0.5] ,...
            'FontSize', 40,...
            'FontName','AvantGarde');
    end
    
end

hold on
%colMap = colormap(summer(numbOfConn));

if printBrain
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 200 100])
    print('-depsc2', 'BrainConn_Bottom_Exported.eps')
    close;
end

%% TOP
% 
figure('position',[0,0, 800,800]);

    h = plot(g);
    set(h,'EdgeColor','none','FaceColor',[0.15 0.15 0.2])
    alpha(0.15);
    view(0,90)
%     colormap(colMap)
%     colorbar;
    
% Anterior/Posterior, Right/Left
if showLabels
    text(0, 85, 30,...
        sprintf( ' %s', 'Ant.'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    text(0, -115, 30,...
        sprintf( ' %s', 'Post.'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    
    text(85, 0, 0,...
        sprintf( ' %s', 'Right'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
    text(-85, 0, 0,...
        sprintf( ' %s', 'Left'),...
        'Color','black',...
        'FontSize', 50,...
        'FontName','AvantGarde');
end

hold on


for i = 1 : length(roiSize)
    plot3(coord(i,1),coord(i,2),coord(i,3), '.',...
        'Color',colMap((roiSize(i)),:),...
        'MarkerSize',roiSize(i)*10);
    
    if showLabels
        text(coord(i,1),coord(i,2),coord(i,3),...
            sprintf( ' %s', regionTag{i}),...
            'Color','black',...
            'BackgroundColor',[0.9 0.9 0.9] ,...
            'EdgeColor',[0.5 0.5 0.5] ,...
            'FontSize', 40,...
            'FontName','AvantGarde');
    end
    
end

hold on
%colMap = colormap(summer(numbOfConn));



% % % % % % for i = 1:size(connStrength,1)
% % % % % %     %colorLine = connStrenghNrm(i);
% % % % % %     line(...
% % % % % %         [coord(connStrength(i,1),1); coord(connStrength(i,2),1)],...
% % % % % %         [coord(connStrength(i,1),2); coord(connStrength(i,2),2)],...
% % % % % %         [coord(connStrength(i,1),3); coord(connStrength(i,2),3)],...
% % % % % %         'LineWidth', 2,...
% % % % % %         'Color', connStrenghColor(i));
% % % % % % end
% % % % 
% % % % for i = 1:size(connStrength,1)
% % % %     mArrow3 (coord(connStrength(i,1),:),coord(connStrength(i,2),:),'color', 'red','stemWidth',0.2);
% % % % end
% % % % 
% % % % colormap(summer)
% % % % cb = colorbar;
% % % % set(cb, 'YTick',0:1/(numbOfConn-1):1,'YTickLabel', cellstr(num2str(sort(connStrength(:,3))))');

%% Save eps file

if printBrain
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 200 100])
    print('-depsc2', 'BrainConn_Top_Exported.eps')
    close;
end

end % end function