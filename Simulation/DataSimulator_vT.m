% Simulation Model

% +++++++++++++++++++ MVAR model ++++++++++++++++++++++++++++++++++

%% Data without convulation
% Iniciallization
L       = 15;  % Number of samples for both models
fs      = 0.5;   % Sampling frequency
ntrials = 100;  % Number of trials
morder_ti = 3;
DataMVGC_II = struct();
y1_ti = zeros(L,1);
y2_ti = zeros(L,1);
y3_ti = zeros(L,1);
y4_ti = zeros(L,1);
y5_ti = zeros(L,1);
x     = zeros(5,L);

for ii=1:ntrials
    for n=morder_ti+1:L
        
        w_ti(n,:)  =[randn(1)  randn(1)  randn(1)  randn(1)  randn(1)]; %  normal-distributed white noise
        
        y1_ti(n) = 0.95*sqrt(2)*y1_ti(n-1) - 0.9025*y1_ti(n-2) + 10*w_ti(n,1);
        y2_ti(n) = 0.5*y1_ti(n-2) + 5*w_ti(n,2);
        y3_ti(n) = -0.4*y1_ti(n-3) + w_ti(n,3);
        y4_ti(n) = -0.5*y1_ti(n-2) + 0.25*sqrt(2)*y4_ti(n-1) + ...
            0.25*sqrt(2)*y5_ti(n-1) + 1.5*w_ti(n,4);
        y5_ti(n) = -0.25*sqrt(2)*y4_ti(n-1) + 0.25*sqrt(2)*y5_ti(n-1) + ...
            2*w_ti(n,5);
                     
        % Without volume conduction effect
        x(:,n,ii) = [y1_ti(n); y2_ti(n); y3_ti(n); y4_ti(n); y5_ti(n)];  % Without the volume conductions effect
        % With volume conduction effect
        %                 x(:,n,ii) = [y1_ti(n) + Vs(1); y2_ti(n) + Vs(2); ...
        %                              y3_ti(n) + Vs(3); y4_ti(n) + Vs(4); ...
        %                              y5_ti(n) + Vs(5)];       
    end
    
%     % +++++++++++++++++++ Generating Channels interactions ++++++++++++
%     V  = rand(5,6);           % Time-constant random mixing matrix of the sources
%     s  = 3*sprand(6,L,0.5);   % Intermittent interactions between channels
%     Vs = V*s;                 % Linear superposition of sparse uniformly distributed random sources
%     
%     % ----- Add channels interactions ------
%     x(:,:,ii) = x(:,:,ii) + Vs;
end

DataMVGC_II.baseline = x;


%% Data with convulution
% Initialization
fs = 0.5; %In Hz
t = 36; %In seconds
l = 30000; %In ms
L = t*fs;
ntrials = 100;
n = 61 : L+60;
morder_ti = 3;

y1_ti = zeros(L+60,1);
y2_ti = zeros(L+60,1);
y3_ti = zeros(L+60,1);
y4_ti = zeros(L+60,1);
y5_ti = zeros(L+60,1);

w_ti  = randn(L+60,5); %  normal-distributed white noise

for ii=1:ntrials
    for n=morder_ti+1:L
                
        y1_ti(n) = 0.95*sqrt(2)*y1_ti(n-1) - 0.9025*y1_ti(n-2) + 10*w_ti(n,1);
        y2_ti(n) = 0.5*y1_ti(n-2) + 5*w_ti(n,2);
        y3_ti(n) = -0.4*y1_ti(n-3) + w_ti(n,3);
        y4_ti(n) = -0.5*y1_ti(n-2) + 0.25*sqrt(2)*y4_ti(n-1) + ...
            0.25*sqrt(2)*y5_ti(n-1) + 1.5*w_ti(n,4);
        y5_ti(n) = -0.25*sqrt(2)*y4_ti(n-1) + 0.25*sqrt(2)*y5_ti(n-1) + ...
            2*w_ti(n,5);
        
        %Discard first 60 points and concatenate
        x(:,n,ii) = [y1_ti(n) ; y2_ti(n) ; y3_ti(n) ; y4_ti(n) ; y5_ti(n) ];
        
    end
end
BOLD = x;
% Convoluçao HRF
clear BOLD BOLD_250 hrf BOLD_05 hrf1 hrf2 DataMVGC_II
for ii=1:ntrials
    for i = 1:5
        p = [6 16 1 1 6 0 32];
        
        hrf(:,i) = spm_hrf(1/fs,p);
        BOLD(i,:,ii) = conv(x(i,:,ii)' , hrf(:,i));
    end
end

DataMVGC_II.baseline = BOLD(:,20:end,:);


