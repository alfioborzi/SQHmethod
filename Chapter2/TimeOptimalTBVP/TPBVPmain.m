%
% set up and solve TBVP 
%
close all; 
clear all;

ltype = {'b-','r--','m-.','b-*','r:','m-x'};

b=0.1;  % \nu
alpha = 0.01; % \mu

% solve 

 m=TPBVP(b,alpha);


 Tex = (1800*b./alpha).^0.2

% plot solution 

 figure(1)
 plot(m(1,:),m(6,:),ltype{1},'LineWidth',2);
 xlabel('t','FontSize',12);
 ylabel('u','FontSize',12)
 print('-depsc2', 'controlTBVP01.eps','-b0'); 
 
 figure(2)
 plot(m(1,:),m(2,:),ltype{1},'LineWidth',2); hold on;
 plot(m(1,:),m(3,:),ltype{2},'LineWidth',2)
 xlabel('t','FontSize',12)
 ylabel('y','FontSize',12)
 legend('y_1','y_2')
 print('-depsc2', 'trajectoryTBVP01.eps','-b0'); 

 