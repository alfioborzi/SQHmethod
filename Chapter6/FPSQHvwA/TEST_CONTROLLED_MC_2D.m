% test controlled MC 2 dim
clear all;
close all;

ltype = {'b-o','m-.','k:','r--'};


fname='data_traj.mat';
load(fname);   % all data needed for the simulation



% number of desired trajectories
Nsample=10;
Ntime = 300; % number of time points for the Milstein integration


T=T1-T0;

[ Nx1, Nx2, Nt] = size(u_1);

D1 = linspace(a,b,Nx1);
D2 = linspace(a,b,Nx2);
TT = linspace(0,T,Nt);

% interpolate controls: 

model

x0 = 1.0;  y0=1.0; %Starting point of the stochastic process


% Integrate the stochastic process
[y, t] = milstein_2D(m,x0,y0,u_1,u_2,Nsample,T,Ntime);


figure(1)

hold on


plot(y(:,:,1),y(:,:,2),'Linewidth',1);
%  xlabel('x_1')
%  ylabel('x_2')

 axis([a b a  b]);
 T=2;
 dt=0.01;

 xd=get_xd(T,dt);
 plot(xd(1,:),xd(2,:),'b--','Linewidth',2);
%  

print('-depsc2', 'trajs-vw-A.eps','-b0','-r300'); 
 print('-dpdf', 'trajs-vw-A.pdf','-b0','-r300');
 