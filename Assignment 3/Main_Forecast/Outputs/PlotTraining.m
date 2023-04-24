clear all;
close all;
clc;

load('LorenzTRNDALL_step=3_rho=10LossN10000.mat');
load('Lorenz_rho=17_Model10.mat');
load('Lorenz_rho=35_Model10.mat');


nfigure = 1;
dt = Param(1);
N = Param(2);
Num_Hidden_Layers = Param(3);
Num_Neurons_per_Layer = Param(4);
steps = Param(5);
TrainingTime = Param(6)



rho = cast(rho, 'double');

if rho == 28
    DataStart = N;
elseif rho == 40
    DataStart = 2*N;
else
    DataStart = 0;
end

DataEndPred = 30;
DataEnd = DataStart+DataEndPred;
t = 0:dt:dt*(DataEndPred-1);



figure(nfigure)
  plot(xpred(1:DataEndPred,1),xpred(1:DataEndPred,3),'b--o');
  hold on
  plot(xtrue(1:DataEndPred,1),xtrue(1:DataEndPred,3),'r.','MarkerSize',10);
  hold off
nfigure++;

figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred(1:DataEndPred,1),'b--o','MarkerSize',12);
  hold on
  plot(xtrue(1:DataEndPred,1),'r.','MarkerSize',22);
  hold off
  title('rho = 10', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('x-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);

nfigure++;

figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred(1:DataEndPred,2),'b--o','MarkerSize',12);
  hold on
  plot(xtrue(1:DataEndPred,2),'r.','MarkerSize',22);
  hold off
  title('rho = 10', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('y-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred(1:DataEndPred,3),'b--o','MarkerSize',12);
  hold on
  plot(xtrue(1:DataEndPred,3),'r.','MarkerSize',22);
  hold off
  title('rho = 10', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('z-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

fprintf("\nRho 10 at n = 30 error is %d",norm(xpred(30,:)-xtrue(30,:)));


figure(nfigure)
  plot(xpred17(1:DataEndPred,1),xpred17(1:DataEndPred,3),'b--o');
  hold on
  plot(sol_17(1:DataEndPred,1),sol_17(1:DataEndPred,3),'r.','MarkerSize',10);
  hold off
nfigure++;

figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred17(1:DataEndPred,1),'b--o','MarkerSize',12);
  hold on
  plot(sol_17(1:DataEndPred,1),'r.','MarkerSize',22);
  hold off
  title('rho = 17', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('x-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

figure(nfigure)
axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred17(1:DataEndPred,2),'b--o','MarkerSize',12);
  hold on
  plot(sol_17(1:DataEndPred,2),'r.','MarkerSize',22);
  hold off
  title('rho = 17', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('y-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

figure(nfigure)
axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred17(1:DataEndPred,3),'b--o','MarkerSize',12);
  hold on
  plot(sol_17(1:DataEndPred,3),'r.','MarkerSize',22);
  hold off
  title('rho = 17', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('z-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

fprintf("\nRho 17 at n = 30 error is %d",norm(xpred17(30,:)-sol_17(30,:)));



figure(nfigure)

  plot(xpred35(1:DataEndPred,1),xpred35(1:DataEndPred,3),'b--o');
  hold on
  plot(sol_35(1:DataEndPred,1),sol_35(1:DataEndPred,3),'r.','MarkerSize',10);
  hold off
nfigure++;

figure(nfigure)
axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred35(1:DataEndPred,1),'b--o','MarkerSize',12);
  hold on
  plot(sol_35(1:DataEndPred,1),'r.','MarkerSize',22);
  hold off
  title('rho = 35', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('x-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:5:100);
nfigure++;

figure(nfigure)
axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred35(1:DataEndPred,2),'b--o','MarkerSize',12);
  hold on
  plot(sol_35(1:DataEndPred,2),'r.','MarkerSize',22);
  hold off
  title('rho = 35', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('y-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20);
  set(gca, 'xtick', 0:5:100);
nfigure++;

figure(nfigure)
axes('FontSize', 25, 'NextPlot', 'add');
  plot(xpred35(1:DataEndPred,3),'b--o','MarkerSize',12);
  hold on
  plot(sol_35(1:DataEndPred,3),'r.','MarkerSize',22);
  hold off
  title('rho = 35', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('z-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20);
  set(gca, 'xtick', 0:5:100);
nfigure++;

fprintf("\nRho 35 at n = 30 error is %d",norm(xpred35(30,:)-sol_35(30,:)));




figure(nfigure)
  plot(loss);
  %xlim([0 1])
  ylim([0 1]);
  title('Loss w.r.t Number of Epochs', 'fontsize', 20);
  xlabel('Iteration','fontsize',20);
  ylabel('Loss','fontsize',20);
nfigure++;
