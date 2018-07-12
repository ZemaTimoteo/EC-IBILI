%% Plot consistency with a boxplot
% Distribution of the consistency values across blocks with 20 s, 60 s and 160 s;
% for all three datasets of 5 var, 0.5 Hz sampled; 250 Hz sampled;
% 1000 Hz sampled, respectively.
%   - The top left plot depicts values of consistency for dataset with SNR = 1.
%  The top right plot depicts values of consistency for dataset with SNR = 200.
%  The bottom plot depicts values of consistency for dataset with SNR = 1000.
%  The horizontal green line, represents the threshold over which the data
%  is being consistent represented in the model.


%% loading
clear all
clc

load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1_20-25.mat')
consist_20(1,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR200_20-25.mat')
consist_20(2,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1000_20-25.mat')
consist_20(3,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1_60-25.mat')
consist_60(1,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR200_60-25.mat')
consist_60(2,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1000_60-25.mat')
consist_60(3,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1_160-25.mat')
consist_160(1,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR200_160-25.mat')
consist_160(2,:) = consist;
load('C:\Users\tiago\Documents\MATLAB\EC-IBILI\ResultsMatrix\Consistency\consist_SNR1000_160-25.mat')
consist_160(3,:) = consist;


cond_names  = {'simulate05','simulate250','simulateOriginal'};
snr_vector  = [1, 200, 1000];

clear consist

%%
% cons{number of points(20,60,160);sampled (0.5,250,Original);snr (1,200,1000)}

% consistency of 20 points
cons{1,1,1} = consist_20(1).simulate05 ; % snr = 1    fs=0.5
cons{1,1,2} = consist_20(2).simulate05 ; % snr = 200  fs=0.5
cons{1,1,3} = consist_20(3).simulate05 ; % snr = 1000 fs=0.5
cons{1,2,1} = consist_20(1).simulate250 ; % snr = 1    fs=250
cons{1,2,2} = consist_20(2).simulate250 ; % snr = 200  fs=250
cons{1,2,3} = consist_20(3).simulate250 ; % snr = 1000 fs=250
cons{1,3,1} = consist_20(1).simulateOriginal ; % snr = 1    fs=Original
cons{1,3,2} = consist_20(2).simulateOriginal ; % snr = 200  fs=Original
cons{1,3,3} = consist_20(3).simulateOriginal ; % snr = 1000 fs=Original

% consistency of 60 points
cons{2,1,1} = consist_60(1).simulate05 ; % snr = 1    fs=0.5
cons{2,1,2} = consist_60(2).simulate05 ; % snr = 200  fs=0.5
cons{2,1,3} = consist_60(3).simulate05 ; % snr = 1000 fs=0.5
cons{2,2,1} = consist_60(1).simulate250 ; % snr = 1    fs=250
cons{2,2,2} = consist_60(2).simulate250 ; % snr = 200  fs=250
cons{2,2,3} = consist_60(3).simulate250 ; % snr = 1000 fs=250
cons{2,3,1} = consist_60(1).simulateOriginal ; % snr = 1    fs=Original
cons{2,3,2} = consist_60(2).simulateOriginal ; % snr = 200  fs=Original
cons{2,3,3} = consist_60(3).simulateOriginal ; % snr = 1000 fs=Original

% consistency of 160 points
cons{3,1,1} = consist_160(1).simulate05 ; % snr = 1    fs=0.5
cons{3,1,2} = consist_160(2).simulate05 ; % snr = 200  fs=0.5
cons{3,1,3} = consist_160(3).simulate05 ; % snr = 1000 fs=0.5
cons{3,2,1} = consist_160(1).simulate250 ; % snr = 1    fs=250
cons{3,2,2} = consist_160(2).simulate250 ; % snr = 200  fs=250
cons{3,2,3} = consist_160(3).simulate250 ; % snr = 1000 fs=250
cons{3,3,1} = consist_160(1).simulateOriginal ; % snr = 1    fs=Original
cons{3,3,2} = consist_160(2).simulateOriginal ; % snr = 200  fs=Original
cons{3,3,3} = consist_160(3).simulateOriginal ; % snr = 1000 fs=Original

%% plotting
for snr_itt = 1:3
    for fs_itt=1:3 % itera pelo sampling
        for points_itt=1:3 % itera pelo numero de pontos no bloco
            cons{snr_itt,fs_itt,points_itt}(find(cons{snr_itt,fs_itt,points_itt})==0) = nan;
        end
    end
end

figure()

for snr_itt = 1:3
    subplot(1,3,snr_itt)
    
    thresh = 0.8;
    boxplot([cons{snr_itt,1,1}', ...
             cons{snr_itt,1,2}', ...
             cons{snr_itt,1,3}'], ...
            'Notch','off','Labels',{'20p','60p','120p'})
    ylim([0.75 1.05]);
    line([0 20],[thresh thresh],'Color','g');         % Legend for the figure
    title(['SNR = ' num2str(snr_vector(snr_itt))])
    axis('square');
    ylabel('Consistency'); 
end

figure()
for snr_itt = 1:3
    subplot(1,3,snr_itt)
    title('fs: 250 Hz')

    thresh = 0.8;
    boxplot([cons{snr_itt,2,1}', ...
             cons{snr_itt,2,2}', ...
             cons{snr_itt,2,3}'], ...
            'Notch','off','Labels',{'20p','60p','120p'})
    ylim([0.75 1.05]);
    line([0 20],[thresh thresh],'Color','g');         % Legend for the figure
    title(['SNR = ' num2str(snr_vector(snr_itt))])
    axis('square');
    ylabel('Consistency');     
end

figure()
for snr_itt = 1:3
    subplot(1,3,snr_itt)
    title('fs: Original')

    thresh = 0.8;
    boxplot([cons{snr_itt,3,1}', ...
             cons{snr_itt,3,2}', ...
             cons{snr_itt,3,3}'], ...
            'Notch','off','Labels',{'20p','60p','120p'})
    ylim([0.75 1.05]);
    line([0 20],[thresh thresh],'Color','g');         % Legend for the figure
    title(['SNR = ' num2str(snr_vector(snr_itt))])
    axis('square');
    ylabel('Consistency');     
end