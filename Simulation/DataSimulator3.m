clear
clc
close all
addpath('mvgc_new')

%% Config
fs = 1000; %In Hz
t = 200; %In seconds
l = 20; %In ms
L = t*fs;
n = 61 : L+60;

y1_ti = zeros(L+60,1);
y2_ti = zeros(L+60,1);
y3_ti = zeros(L+60,1);
y4_ti = zeros(L+60,1);
y5_ti = zeros(L+60,1);

w_ti  = randn(L+60,5); %  normal-distributed white noise

%% Model
y1_ti(n) = 0.95*sqrt(2)*y1_ti(n-l) - 0.9025*y1_ti(n-2*l) + w_ti(n,1);
y2_ti(n) = 0.5*y1_ti(n-2*l) + w_ti(n,2);
y3_ti(n) = -0.4*y1_ti(n-3*l) + w_ti(n,3);
y4_ti(n) = -0.5*y1_ti(n-2*l) + 0.25*sqrt(2)*y4_ti(n-l) + ...
    0.25*sqrt(2)*y5_ti(n-l) + w_ti(n,4);
y5_ti(n) = -0.25*sqrt(2)*y4_ti(n-l) + 0.25*sqrt(2)*y5_ti(n-l) + ...
    w_ti(n,5);

%Discard first 60 points and concatenate
X = [y1_ti(61:end) , y2_ti(61:end) , y3_ti(61:end) , y4_ti(61:end) , y5_ti(61:end) ];

X_250 = downsample(X,4);

X_05 = downsample(X,2000);

clear BOLD BOLD_250 hrf BOLD_05 hrf1 hrf2

for i = 1:5
    p = [6 16 1 1 6 0 32];
    
    hrf(:,i) = spm_hrf(1/1000,p);
    hrf1(:,i) = spm_hrf(1/250,p);
    hrf2(:,i) = spm_hrf(1/0.5,p);
    
    BOLD(:,i) = conv(X(:,i) , hrf(:,i));
    BOLD_250(:,i) = conv(X_250(:,i) , hrf1(:,i));
    BOLD_05(:,i) = conv(X_05(:,i) , hrf2(:,i));
end

% Remove first 50s of data
BOLD = BOLD(50*fs+1:end,:);
BOLD_250 = BOLD_250(50*250+1:end,:);
BOLD_05 = BOLD_05(50*0.5+1:end,:);

%% Pos downsample
BOLD_pos_250 = downsample(BOLD,4);
BOLD_pos_05 = downsample(BOLD,2000);

%%
window = 3;
n_rois = 5;
n_con = (n_rois*(n_rois-1))/2; % #Connections
C = combnk(1:n_rois,2); % Possible Combinations

%Window Iteration
for t = window : 91
    
    %Combinations Iteration
    for l=1:n_con
        
        temp1 = BOLD_05(:,C(l,1));
        temp2 = BOLD_05(:,C(l,2));
        
        x = t-window+1 : t;
        temp3 = corrcoef(temp1(x),temp2(x));
        
        Results.VOI_Corr{l}(t,1) = temp3(1,2);
        Results.VOI_CorrF{l}(t,1) = atanh(temp3(1,2)); %Fisher
        
    end

end
load('colormap(-1,1).mat');
mean_results = zeros(n_rois);

for j =1:n_con
    
    a = Results.VOI_Corr{1,j};
    
    
    mean_results(C(j,1),C(j,2)) = mean(a);
    mean_results(C(j,2),C(j,1)) = mean(a);  
end

figure
    plot_pw_full(mean_results,newColorMap)
    colorbar
%%
mean_results_full = corrcoef(BOLD_pos_05);

figure
    plot_pw_full(mean_results_full,newColorMap)
    colorbar
