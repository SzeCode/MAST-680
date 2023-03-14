%PCA and SINDy Implementation
%Will Sze
%40096561
clear all; close all; clc;

%Loading Data from extracted location of each camera
load cam1_1_xPos2.mat
load cam1_1_yPos2.mat
cam1_xPos = x_pos;
cam1_yPos = y_pos;

load cam2_1_xPos.mat
load cam2_1_yPos.mat
cam2_xPos = x_pos(1:length(cam1_xPos));
cam2_yPos = y_pos(1:length(cam1_yPos));

load cam3_1_xPos.mat
load cam3_1_yPos.mat
cam3_xPos = x_pos(1:length(cam1_xPos));
cam3_yPos = y_pos(1:length(cam1_yPos));

clear 'x_pos';
clear 'y_pos';

nfigure = 1;

%Plots all three camera y positions
figure(nfigure);
hold on
  plot(cam1_yPos);
  plot(cam2_yPos);
  plot(cam3_yPos);
hold off
nfigure++;

%Removes first few data
Start = 25;
X_og(1,:) = cam1_xPos(Start:end);
X_og(2,:) = cam1_yPos(Start:end);
X_og(3,:) = cam2_xPos(Start:end);
X_og(4,:) = cam2_yPos(Start:end);
X_og(5,:) = cam3_xPos(Start:end);
X_og(6,:) = cam3_yPos(Start:end);

for i = 1:6
  X_og_mean(i) = mean(X_og(i,:));
  X_PCA(i,:) = X_og(i,:) - X_og_mean(i);
end

%Plotting the mean subtracted data
figure(nfigure);
hold on
  plot(X_PCA(2,:),'r.','MarkerSize', 10);
  plot(X_PCA(4,:),'g.','MarkerSize', 10);
  plot(X_PCA(6,:),'b.','MarkerSize', 10);
hold off
nfigure++;


%--------- Applying PCA -----------%
A = 1/sqrt(length(X_og)-1)*X_PCA;

%A_AT = A*A'; %Alternative
%[eig_Vect, eig_Val] = eig(A_AT); %Alternative

[U,S,V] = svd(A,'econ');
Rank_S = 2; %First 2 principal components
Y_PCA = U(:,1:Rank_S)'*X_PCA;

%Ploting Original mean subtracted data and the data obtained from PCA
figure(nfigure);

  subplot (2, 1, 1)
  hold on
  plot(X_PCA(2,:),'r','MarkerSize', 10);
  plot(X_PCA(4,:),'g','MarkerSize', 10);
  plot(X_PCA(6,:),'b','MarkerSize', 10);
  hold off

  title("Original y Position");
  h = legend ("cam1","cam2","cam3");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

  subplot (2, 1, 2)
  hold on
  plot(Y_PCA(1,:),'r.','MarkerSize', 10);
  plot(Y_PCA(2,:),'b.','MarkerSize', 10);
  hold off

  title("PCA Coordinates");
  h = legend ("Unfiltered PCA Position","Unfiltered PCA Velocity");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

nfigure++;

%Comparing with a true sine curve
t = linspace(0,202,203);
y_sin = 90*sin(2*pi/40*t-pi/2-0.5);

figure(nfigure)

  hold on
  plot(Y_PCA(1,:),'r.','MarkerSize', 10);
  plot(y_sin,'b','MarkerSize', 10);

  %Filtering
  windowSize = 5;
  b = (1/windowSize)*ones(1,windowSize);
  a = 1;

  Y_PCA(1,:) = filter(b,a,Y_PCA(1,:));
  Y_PCA(2,:) = filter(b,a,Y_PCA(2,:));

  plot(Y_PCA(1,:),'g','MarkerSize', 10);
  hold off

  title("Comparison Between PCA Position, true sine curve and filtered PCA Position");
  h = legend ("PCA Position","True sine curve", "Filtered PCA Position");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

nfigure++;

%Plotting the filtered PCA Data
figure(nfigure);

  subplot (2, 1, 1)
  hold on
  plot(X_PCA(2,:),'r','MarkerSize', 10);
  plot(X_PCA(4,:),'g','MarkerSize', 10);
  plot(X_PCA(6,:),'b','MarkerSize', 10);
  hold off

  title("Original y Position");
  h = legend ("cam1","cam2","cam3");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

  subplot (2, 1, 2)
  hold on
  plot(Y_PCA(1,:),'r.','MarkerSize', 10);
  plot(Y_PCA(2,:),'b.','MarkerSize', 10);
  hold off

  title("PCA Coordinates");
  h = legend ("filtered PCA Position","filtered PCA Velocity");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");


nfigure++;

save("Y_PCA.mat",'Y_PCA');


%--------- Applying SINDy -----------%
Framerate = 1;
dt = 1/Framerate;
N_Start = 28;
N_End = 153;
X_SINDy(1,:) = Y_PCA(1,N_Start:(N_End-1));
X_SINDy(2,:) = Y_PCA(2,N_Start:(N_End-1));

Y_SINDy_dxdt = (Y_PCA(:,N_Start+1:(N_End))-Y_PCA(:,N_Start:(N_End-1)))/dt; %derivative
Y_SINDy = Y_PCA(:,N_Start+1:(N_End)) - Y_PCA(:,N_Start).*ones([size(Y_PCA(:,N_Start+1:(N_End)))]);

dxdt_Max = max(Y_SINDy_dxdt(1,:));
scaling = dxdt_Max/max(X_SINDy(2,:))
X_SINDy(2,:) = Y_PCA(2,N_Start:(N_End-1))*scaling;

%Plotting the principal components and their derivatives
figure(nfigure);

  hold on
  plot(X_SINDy(1,:),'r.','MarkerSize', 10);
  plot(X_SINDy(2,:),'b.','MarkerSize', 10);
  plot(Y_SINDy_dxdt(1,:),'g','MarkerSize', 10);
  hold off

  title("PCA Coordinates and derivative comparison");
  h = legend ("Filtered PCA Position","Filtered & Scaled PCA Velocity","Derivative PCA Position");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

nfigure++;

figure(nfigure);

  hold on
  plot(X_SINDy(1,:),'r.','MarkerSize', 10);
  plot(Y_SINDy_dxdt(1,:),'g.','MarkerSize', 10);
  hold off

  title("PCA Position and its derivative");
  h = legend ("PCA Position","Derivative PCA Position");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

nfigure++;

figure(nfigure);

  hold on
  plot(X_SINDy(2,:),'b.','MarkerSize', 10);
  plot(Y_SINDy_dxdt(2,:),'o.','MarkerSize', 10);
  hold off

  title("PCA Velocity and its derivative");
  h = legend ("PCA Velocity","Derivative PCA Velocity");
  xlabel ("Frame Number");
  ylabel ("Pixel Position");

nfigure++;

% Defining Dictionary (1, x1, x2, x1^2, x1*x2, x2^2);
Theta = ones([1,length(X_SINDy)]);
Theta(2:3,:) = X_SINDy;
Theta(4,:) = X_SINDy(1,:).^2;
Theta(5,:) = X_SINDy(1,:).*X_SINDy(2,:);
Theta(6,:) = X_SINDy(2,:).^2;

%Applying integration method
ThetaInt = zeros([size(Theta)]);
ThetaInt(:,1) = dt*Theta(:,1);

%reman sum
for i = 2:length(Theta(1,:))
  ThetaInt(:,i) = ThetaInt(:,i-1)+dt*Theta(:,i);
end

coeff = Y_SINDy*pinv(ThetaInt);

%SINDy method begins %adapted from Jason Bramburger's SINDy example code: https://github.com/jbramburger/DataDrivenDynSyst/tree/main/Identifying%20Nonlinear%20Dynamics
lam = 0.1;
k = 1;

Coeff_New = coeff;
Err = sum(sum((abs(coeff - Coeff_New)))); %Coefficient checking
while k == 1 || Err > 0;

  coeff = Coeff_New; %Save new coeff
  smallinds = (abs(coeff)<lam);
  Coeff_New(smallinds) = 0; %Set new coeff small indices to 0;

  for i = 1:2
    biginds = ~smallinds(i,:);
    Coeff_New(i,biginds) = Y_SINDy(i,:)*pinv(ThetaInt(biginds,:)); %Find new coeff for big indices;
  end

  k = k + 1;
  Err = sum(sum((abs(coeff - Coeff_New))));
end

%Display results
fprintf('Expected Output: \n')
fprintf('z1_dot = z2 \n')
fprintf('z2_dot = -Az2 - Bz1 \n')
fprintf('B is K/m (2pi/T)^2 should be around 0.025\n')

mons2 = {''; 'z1'; 'z2'; 'z1^2'; 'z1z2'; 'z2^2'};

fprintf('\nSINDy Output with lam = %d', lam)
fprintf(' and dt = %d:\n', dt)
fprintf('z1_dot = ')
for i = 1:length(Coeff_New)

  if Coeff_New(1,i) < -1e-5;
    mons2_get = mons2{i};
    fprintf(" - %d %s", abs(Coeff_New(1,i)),mons2_get);
  elseif Coeff_New(1,i) > 1e-5
    mons2_get = mons2{i};
    fprintf(" + %d %s", abs(Coeff_New(1,i)),mons2_get);
  end

end
fprintf('\n')


fprintf('z2_dot = ')
for i = 1:length(Coeff_New)

  if Coeff_New(2,i) < -1e-5
    mons2_get = mons2{i};
    fprintf(" - %d %s", abs(Coeff_New(2,i)),mons2_get);
  elseif Coeff_New(2,i) > 1e-5
    mons2_get = mons2{i};
    fprintf(" + %d %s", abs(Coeff_New(2,i)),mons2_get);
  end

end
fprintf('\n')

%Simulate results and plot
fprintf('\nSimulating SINDy discovered model...\n')
[t_sim,Y_PCA_Sim] = RK4_Sys_Inputs(@f, 0, [85,0], 126, 125, 0);
figure(nfigure);

hold on
plot(Y_PCA(1,N_Start:(N_End-1)),'r.','MarkerSize', 10);
plot(Y_PCA_Sim(:,1),'b','MarkerSize', 10);
hold off

title("PCA Position and Simulated SINDy model using RK4 for comparison");
h = legend ("Filtered PCA Position","RK4 Simulated SINDy Model");
xlabel ("Frame Number");
ylabel ("Pixel Position");

nfigure++;
