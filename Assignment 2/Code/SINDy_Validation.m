clear all;
close all;
clc;

load Y_PCA.mat

nfigure = 1;
%--------- SINDy Validation -----------%
Framerate = 60;
dt = 1/Framerate;


m = 1; K = 1; c = 0.1;
beta = sqrt(abs((c/(2*m))^2-(K/m)));
alpha = c/(2*m);

A = 10;
B = A*alpha/beta;

x = @(t) exp(-alpha*t)*(A*cos(beta*t)+B*sin(beta*t));
x_dot = @(t) exp(-alpha*t)*((-A*alpha+B*beta)*cos(beta*t)-(A*beta+B*alpha)*sin(beta*t));

t = 0;
for i = 1:1000
  X_Data(1,i) = x(t);
  X_Data(2,i) = x_dot(t);
  t = t + dt;
end

##A = 100;
##B = A*alpha/beta;
##t = 0;
##x = @(t) exp(-alpha*t)*(A*cos(beta*t)+B*sin(beta*t));
##x_dot = @(t) exp(-alpha*t)*((-A*alpha+B*beta)*cos(beta*t)-(A*beta+B*alpha)*sin(beta*t));
##for i = 501:1000
##  X_Data(1,i) = x(t);
##  X_Data(2,i) = x_dot(t);
##  t = t + dt;
##end


figure(nfigure);
  plot(X_Data(1,:));
  hold on
  plot(X_Data(2,:));
nfigure++;

X_SINDy(1,:) = X_Data(1,1:(end-1));
X_SINDy(2,:) = X_Data(2,1:(end-1));

Y_SINDy = (X_Data(:,2:(end))-X_Data(:,1:(end-1)))/dt; %derivative
##Y_SINDy(:,1:499) = (X_Data(:,2:(500))-X_Data(:,1:(500-1)))/dt; %derivative
##Y_SINDy(:,501:999) = (X_Data(:,502:(1000))-X_Data(:,501:(1000-1)))/dt; %derivative
%Y_SINDy = X_Data(:,2:end) - X_Data(:,1).*ones([size(X_Data(:,2:end))]);


figure(nfigure);

  hold on
  plot(X_SINDy(1,:),'r.','MarkerSize', 10);
  plot(Y_SINDy(1,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

figure(nfigure);

  hold on
  plot(X_SINDy(2,:),'b.','MarkerSize', 10);
  plot(Y_SINDy(2,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

% Defining Dictionary (1, x1, x2, x1^2, x1*x2, x2^2);
Theta = ones([1,length(X_SINDy(1,:))]);
Theta(2:3,:) = X_SINDy;
Theta(4,:) = X_SINDy(1,:).^2;
Theta(5,:) = X_SINDy(1,:).*X_SINDy(2,:);
Theta(6,:) = X_SINDy(2,:).^2;


##ThetaInt = zeros([size(Theta)]);
##ThetaInt(:,1) = dt*Theta(:,1);
##
##%reman sum
##for i = 2:length(Theta(:,1))
##  ThetaInt(:,i) = ThetaInt(:,i-1)+dt*Theta(:,i);
##end

coeff = Y_SINDy*pinv(Theta);

lam = 0.1;
k = 1;

Coeff_New = coeff;
Err = sum(sum((abs(coeff - Coeff_New))));
while k == 1 || Err > 0;

  coeff = Coeff_New;

  smallinds = (abs(coeff)<lam);

  Coeff_New(smallinds) = 0;

  for i = 1:2

    biginds = ~smallinds(i,:);

    Coeff_New(i,biginds) = Y_SINDy(i,:)*pinv(Theta(biginds,:));

  end

  k = k + 1;
  Err = sum(sum((abs(coeff - Coeff_New))));
end



##X_SINDy_dot = zeros([size(Y_SINDy)]);
##for i = 1:6
##
##  X_SINDy_dot = X_SINDy_dot + coeff(:,i)*Theta(i,:);
##
##end
fprintf('Using Derivative ------------------------- \n')
fprintf('Model Parameters: \n')
fprintf('Mass (m) = %d\n',m)
fprintf('Spring Constant (K) = %d\n',K)
fprintf('Damping Coefficient (c) = %d\n',c)

fprintf('\nExpected Output: \n')
fprintf('z1_dot = z2 \n')
fprintf('z2_dot = -%d z1',K/m)
fprintf(' - 0.1 z2 \n',c/m)


mons2 = {''; 'z1'; 'z2'; 'z1^2'; 'z1z2'; 'z2^2'};


fprintf('\nSINDy Output: \n')
fprintf('z1_dot = ')
for i = 1:length(Coeff_New)

  if Coeff_New(1,i) < -1e-5;
    mons2_get = mons2{i};
    fprintf(" - %d %s", abs(Coeff_New(1,i)),mons2_get);
  elseif Coeff_New(1,i) > 1e-5
    mons2_get = mons2{i};
    fprintf(" + %d %s", abs(Coeff_New(1,i)),mons2_get);
  else

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
  else


  end

end
fprintf('\n')




% ------- Integral implementation ---------%
noiseFactor = 0.1;
for i = 1:2
  for j = 1:(length(X_Data))
    noiseGen = 2*noiseFactor*rand(1)-(noiseFactor);
    X_Data_noise(i,j) = noiseGen+X_Data(i,j);
  end
end

X_SINDy = X_Data_noise(:,1:(end-1));
Y_SINDy = zeros(size(X_SINDy));
Y_SINDy = X_Data_noise(:,2:end) - X_Data_noise(:,1).*ones(size(X_Data_noise(:,2:end)));
%Y_SINDy = (X_Data(:,2:(end))-X_Data(:,1:(end-1)))/dt; %derivative

figure(nfigure);

  hold on
  plot(X_SINDy(1,:),'r.','MarkerSize', 10);
  plot(Y_SINDy(1,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

figure(nfigure);

  hold on
  plot(X_SINDy(2,:),'b.','MarkerSize', 10);
  plot(Y_SINDy(2,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

% Defining Dictionary (1, x1, x2, x1^2, x1*x2, x2^2);
Theta = ones([1,length(X_SINDy(1,:))]);
Theta(2:3,:) = X_SINDy;
Theta(4,:) = X_SINDy(1,:).^2;
Theta(5,:) = X_SINDy(1,:).*X_SINDy(2,:);
Theta(6,:) = X_SINDy(2,:).^2;


ThetaInt = zeros([size(Theta)]);
ThetaInt(:,1) = dt*Theta(:,1);

%reman sum
for i = 2:length(Theta(1,:))
  ThetaInt(:,i) = ThetaInt(:,i-1)+dt*Theta(:,i);
end

coeff = Y_SINDy*pinv(ThetaInt);

lam = 0.1;
k = 1;

Coeff_New = coeff;
Err = sum(sum((abs(coeff - Coeff_New))));
while k == 1 || Err > 0;

  coeff = Coeff_New;

  smallinds = (abs(coeff)<lam);

  Coeff_New(smallinds) = 0;

  for i = 1:2

    biginds = ~smallinds(i,:);

    Coeff_New(i,biginds) = Y_SINDy(i,:)*pinv(ThetaInt(biginds,:));

  end

  k = k + 1;
  Err = sum(sum((abs(coeff - Coeff_New))));
end



##X_SINDy_dot = zeros([size(Y_SINDy)]);
##for i = 1:6
##
##  X_SINDy_dot = X_SINDy_dot + coeff(:,i)*Theta(i,:);
##
##end

fprintf('\nIntegration SINDy begin -------------------------- \n')
fprintf('Model Parameters: \n')
fprintf('Mass (m) = %d\n',m)
fprintf('Spring Constant (K) = %d\n',K)
fprintf('Damping Coefficient (c) = %d\n',c)

fprintf('\nExpected Output: \n')
fprintf('z1_dot = z2 \n')
fprintf('z2_dot = -%d z1',K/m)
fprintf(' - 0.1 z2 \n',c/m)

mons2 = {''; 'z1'; 'z2'; 'z1^2'; 'z1z2'; 'z2^2'};


fprintf('\nSINDy Output: \n')
fprintf('z1_dot = ')
for i = 1:length(Coeff_New)

  if Coeff_New(1,i) < -1e-5;
    mons2_get = mons2{i};
    fprintf(" - %d %s", abs(Coeff_New(1,i)),mons2_get);
  elseif Coeff_New(1,i) > 1e-5
    mons2_get = mons2{i};
    fprintf(" + %d %s", abs(Coeff_New(1,i)),mons2_get);
  else

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
  else


  end

end
fprintf('\n')



windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

Y_PCA(1,:) = filter(b,a,Y_PCA(1,:));
Y_PCA(2,:) = filter(b,a,Y_PCA(2,:));


%----- Test
N_Start = 28;
N_End = 153;
X_SINDy = Y_PCA(:,N_Start:(N_End-1));
Y_SINDy = zeros(size(X_SINDy));
Y_SINDy = Y_PCA(:,(N_Start+1):N_End) - Y_PCA(:,N_Start).*ones(size(Y_PCA(:,(N_Start+1):N_End)));
%Y_SINDy = (X_Data(:,2:(end))-X_Data(:,1:(end-1)))/dt; %derivative

figure(nfigure);

  hold on
  plot(X_SINDy(1,:),'r.','MarkerSize', 10);
  plot(Y_SINDy(1,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

figure(nfigure);

  hold on
  plot(X_SINDy(2,:),'b.','MarkerSize', 10);
  plot(Y_SINDy(2,:),'g.','MarkerSize', 10);
  hold off

nfigure++;

% Defining Dictionary (1, x1, x2, x1^2, x1*x2, x2^2);
Theta = ones([1,length(X_SINDy(1,:))]);
Theta(2:3,:) = X_SINDy;
Theta(4,:) = X_SINDy(1,:).^2;
Theta(5,:) = X_SINDy(1,:).*X_SINDy(2,:);
Theta(6,:) = X_SINDy(2,:).^2;

dt = 1/60;
ThetaInt = zeros([size(Theta)]);
ThetaInt(:,1) = dt*Theta(:,1);

%reman sum
for i = 2:length(Theta(1,:))
  ThetaInt(:,i) = ThetaInt(:,i-1)+dt*Theta(:,i);
end

coeff = Y_SINDy*pinv(ThetaInt);

lam = 0.01;
k = 1;

Coeff_New = coeff;
Err = sum(sum((abs(coeff - Coeff_New))));
while k == 1 || Err > 0;

  coeff = Coeff_New;

  smallinds = (abs(coeff)<lam);

  Coeff_New(smallinds) = 0;

  for i = 1:2

    biginds = ~smallinds(i,:);

    Coeff_New(i,biginds) = Y_SINDy(i,:)*pinv(ThetaInt(biginds,:));

  end

  k = k + 1;
  Err = sum(sum((abs(coeff - Coeff_New))));
end



##X_SINDy_dot = zeros([size(Y_SINDy)]);
##for i = 1:6
##
##  X_SINDy_dot = X_SINDy_dot + coeff(:,i)*Theta(i,:);
##
##end

fprintf('\nActual Data -------------------------- \n')
fprintf('Model Parameters: \n')
fprintf('Mass (m) = %d\n',m)
fprintf('Spring Constant (K) = %d\n',K)
fprintf('Damping Coefficient (c) = %d\n',c)

fprintf('\nExpected Output: \n')
fprintf('z1_dot = z2 \n')
fprintf('z2_dot = -A z1 - B z2 \n')

mons2 = {''; 'z1'; 'z2'; 'z1^2'; 'z1z2'; 'z2^2'};

fprintf('\nSINDy Output: \n')
fprintf('z1_dot = ')
for i = 1:length(Coeff_New)

  if Coeff_New(1,i) < -1e-5;
    mons2_get = mons2{i};
    fprintf(" - %d %s", abs(Coeff_New(1,i)),mons2_get);
  elseif Coeff_New(1,i) > 1e-5
    mons2_get = mons2{i};
    fprintf(" + %d %s", abs(Coeff_New(1,i)),mons2_get);
  else

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
  else


  end

end
fprintf('\n')
