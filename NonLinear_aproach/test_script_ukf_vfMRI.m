%% Example - test script - ukf running

clear all
clc
close all

%%

% Variables definition -  [x,P]=ukf(fstate,x,P,hmeas,z,Q,R)

n=3;      %number of state - ?number of nodes?
q=0.1;    %std of process 
r=0.1;    %std of measurement

Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  

f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation

s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance

N=20;                                     % total dynamic steps - número total de pontos?
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual


zV = zeros(1,N);

%% aplicação de ukf
for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end

%% Plot Results
for k=1:3                                 
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end

%%
