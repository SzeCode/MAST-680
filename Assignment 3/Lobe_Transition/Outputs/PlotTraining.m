clear all;
close all;
clc;

load('LorenztJumpPred_rho=28Saw.mat');

nfigure = 1;
dt = Param(1);
N = Param(2);
Num_Hidden_Layers = Param(3);
Num_Neurons_per_Layer = Param(4);
TrainingTime = Param(5)



rho = cast(rho, 'double');

DataStart = 450;
DataEndPred = 100;
DataEnd = 1000;
t = 0:dt:dt*(DataEndPred-1);


figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(tpred(1:DataEnd,1),'b--o','MarkerSize',12);
  hold on
  plot(ttrue(1:DataEnd,1),'r.','MarkerSize',22);
  hold off
  title('Number of Iterates Remaining for the next Transition', 'fontsize', 20);
  xlabel('n','fontsize',40);
  ylabel('n remaining','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:100:1000);
  %xlim([470 530])
nfigure++;

figure(nfigure)
  axes('FontSize', 25, 'NextPlot', 'add');
  plot(1/10*tpred(1:DataEnd,1),'b--o','MarkerSize',2);
  hold on
  plot(FullSOl(1:DataEnd,1),'r.','MarkerSize',10);
  hold off
  title('Comparison Transition Predictions and X-position of Lorenz', 'fontsize', 20);
  xlabel('n','fontsize',40);
  %ylabel('x-position','fontsize',40);
  legend('Model Prediction','True Solution', 'fontsize',20)
  set(gca, 'xtick', 0:100:1000);
  %xlim([470 530])
nfigure++;




figure(nfigure)
  plot(loss);
  %xlim([0 1])
  ylim([0 1]);
  title('Loss w.r.t Number of Epochs', 'fontsize', 20);
  xlabel('Iteration','fontsize',20);
  ylabel('Loss','fontsize',20);
nfigure++;
