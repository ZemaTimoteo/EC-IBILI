%% FFT test - check if signal is change dramaticly after re-arrangment


%% ... reorganization of the data & Initialization ...
DataTest{1} = test.Data{1};
DataTest{2} = test.Data{2};
Fs = 0.5;                                    % Sampling frequency
T = 1/Fs;                                    % Sampling period 

%% Test data within trials subject specific

DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),96,10]);
[nodes,L,subjects] = size(DataTest{1});      % Length of signal
t = (0:L-1)*T;                               % Time vector
f = Fs*(0:(L/2))/L;                          % frequency vector
DataTest{1} = DataTest{1}-mean(DataTest{1},2);

for j=1:subjects
    figure()
    for i=1:nodes
        % fft process
        Y     = fft(DataTest{1}(i,:,1));
        P2    = abs(Y/L);
        P1{j} = P2(1:L/2+1);
        P1{j}(2:end-1) = 2*P1{j}(2:end-1);
        
        %downsampled for comparison
% %         P1{j} = downsample(P1{j},2);
%         L = size(P1{j},2);
%         f = Fs*(0:L-1)/L;
        
        subplot(3,2,i)
        plot(f,P1{j})
        title(strcat('Signal Amplitude Spectrum of variable .',num2str(i)))
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        clear P2 Y
    end
end



%% Test data within trial subject mix
% % numberSHUF = 10;

% % for pp=1:numberSHUF
    pp = 1;
    % ... reorganization of the data ...
    DataTest{1} = test.Data{1};
%     DataTest{2} = test.Data{2};
    
    if length(DataTest{1}) == 960
        % Data condition 1
        % - randomize subjects
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),96,10]);
        [~,~,s] = size(DataTest{1});
        shuffl{pp}  = randperm(s) ;
        DataTest{1}(:,:,shuffl{pp}) = DataTest{1}(:,:,:);
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),960,1]);
        % - rearrangement of data
        rshapData{1,1} = reshape(DataTest{1},[size(DataTest{1},1),48,20]);
        
        [nodes,Ll,subjects] = size(rshapData{1,1});     % Length of signal
        tt = (0:Ll-1)*T;                                % Time vector
        ff = Fs*(0:(Ll/2))/Ll;                          % frequency vector
        % fft process Data Condition 1
        for j=1:subjects
            figure()
            
            for i=1:nodes
                Yy     = fft(rshapData{1,1}(i,:,j));
                Pp2    = abs(Yy/Ll);
                Pp1{j} = Pp2(1:Ll/2+1);
                Pp1{j}(2:end-1) = 2*Pp1{j}(2:end-1);
                
                subplot(3,2,i)
                plot(ff,Pp1{j})
                title(strcat('Signal Amplitude Spectrum of variable .',num2str(i)))
                xlabel('f (Hz)')
                ylabel('|P1(f)|')
                clear Pp2 Yy
            end
        end

% % %         % Data condition 2
% % %         % - randomize subjects
% % %         DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),96,10]);
% % %         DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
% % %         DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),960,1]);
% % %         % rearrangement of data
% % %         test.rshapData{2,1} = reshape(DataTest{2},[size(DataTest{2},1),48,20]);
        
    elseif length(DataTest{1}) == 2880
        % Data condition 1
        % - randomize subjects
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),288,10]);
        [~,~,s] = size(DataTest{1});
        test.shuffl{pp}  = randperm(s) ;
        DataTest{1}(:,:,test.shuffl{pp}) = DataTest{1}(:,:,:);
        DataTest{1} = reshape(DataTest{1},[size(DataTest{1},1),2880,1]);
        
        % - rearrangement of data
        a = DataTest{1}(:,1:1152);
        b = DataTest{1}(:,1152+1:1152+864);
        c = DataTest{1}(:,1152+864+1:1152+864+864);
        
        test.rshapData{1,1} = reshape(a,[size(DataTest{1},1),48,24]);
        test.rshapData{1,2} = reshape(b,[size(DataTest{1},1),36,24]);
        test.rshapData{1,3} = reshape(c,[size(DataTest{1},1),36,24]);
        
% % %         % Data condition 2
% % %         % - randomize subjects
% % %         DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),288,10]);
% % %         DataTest{2}(:,:,test.shuffl{pp}) = DataTest{2}(:,:,:);
% % %         DataTest{2} = reshape(DataTest{2},[size(DataTest{2},1),2880,1]);
% % %         
% % %         
% % %         % - rearrangement of data
% % %         aa = DataTest{2}(:,1:1152);
% % %         bb = DataTest{2}(:,1152+1:1152+864);
% % %         cc = DataTest{2}(:,1152+864+1:1152+864+864);
% % %         
% % %         test.rshapData{2,1} = reshape(aa,[size(DataTest{2},1),48,24]);
% % %         test.rshapData{2,2} = reshape(bb,[size(DataTest{2},1),36,24]);
% % %         test.rshapData{2,3} = reshape(cc,[size(DataTest{2},1),36,24]);
    end
    
    clear a b c s
% % %     clear aa bb cc ss
    
% % % end




%% Statisitcal comparison between P1 and Pp1 coherence

% get data
DataTest{1} = test.Data{1};
DataTest{2} = test.Data{2};

% frequency and Period
Fs = 0.5;                                    % Sampling frequency
T = 1/Fs;                                    % Sampling period 

% re-arrange data per subject
DataPerSubjTest{1} = reshape(DataTest{1},[size(DataTest{1},1),96,10]);
[nodes,L,trialsPerSubj] = size(DataPerSubjTest{1});      % Length of signal

% re-arrange data per method
if length(DataTest{1}) == 960
    % - randomize subjects
    DataPerMethodTest{1} = reshape(DataTest{1},[size(DataTest{1},1),96,10]);
    [~,~,s] = size(DataPerMethodTest{1});
    shuffl  = randperm(s) ;
    DataPerMethodTest{1}(:,:,shuffl) = DataPerMethodTest{1}(:,:,:);
    DataPerMethodTest{1} = reshape(DataPerMethodTest{1},[size(DataPerMethodTest{1},1),960,1]);
    % - rearrangement of data
    DataPerMethodTest{1,1} = reshape(DataPerMethodTest{1},[size(DataPerMethodTest{1},1),48,20]);
    
    [~,Ll,trialsPerMethods] = size(DataPerMethodTest{1,1});     % Length of signal
    
elseif length(DataTest{1}) == 2880
    % Data condition 1
    % - randomize subjects
    DataPerMethodTest{1} = reshape(DataTest{1},[size(DataTest{1},1),288,10]);
    [~,~,s] = size(DataPerMethodTest{1});
    shuffl  = randperm(s) ;
    DataPerMethodTest{1}(:,:,shuffl) = DataPerMethodTest{1}(:,:,:);
    DataPerMethodTest{1} = reshape(DataPerMethodTest{1},[size(DataPerMethodTest{1},1),2880,1]);
    
    [~,Ll,trialsPerMethods] = size(DataPerMethodTest{1,1});     % Length of signal
end


%% test coherence
clear cxy coherencyAll

j=0;
auxJ = [1 0 2 0 3 0 4 0 5 0 6 0 7 0 8 0 9 0 10 0];

while j<trialsPerMethods-2
    for g=1:trialsPerSubj/2
        for i=1:nodes
%             dataA = downsample(DataPerSubjTest{1}(i,:,g),2);
%             dataA = DataPerSubjTest{1}(i,1:size(DataPerMethodTest{1}(i,:,j),2),g);
            dataA = DataPerSubjTest{1}(i,:,g);
%             dataAA = DataPerSubjTest{1}(i,:,g+5);
%             dataB = DataPerMethodTest{1}(i,:,j);
            dataB = [DataPerMethodTest{1}(i,:,j+1) DataPerMethodTest{1}(i,:,j+2)];
           
            %tests correlation & coherence
%             CoeffCorr{j,g} = corrcoef(dataA,dataB);
            [cxy{g,g},f] = mscohere (dataA,dataAA,[],[],[],0.5);
        end
    end
    j=j+2;
end

CoherencyAll = horzcat(cxy{:,:});
plot(f,CoherencyAll)
title('Magnitude-Squared Coherence')
xlabel('Frequency (Hz)')
grid

figure
h1 = histogram(CoherencyAll);





%%%%--- Exames a experimentar amanhã
% in time - fazer covariancia, cross-correlatio coef., 
% in frequency - statistical inference between 2 signals, coherence