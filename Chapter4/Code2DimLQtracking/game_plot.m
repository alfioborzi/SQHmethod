%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for differential Nash games - plots
%
close all;

ltype = {'b-','r--','m-.','b-*','r:','m-x'};

N=10001;
T=1.0;
time_vec = linspace( 0.0 , T , N );


load('u.mat');
u1S= u1;
u2S= u2;

load('y.mat');
yS=y;

load('j1.mat');
J1=J1k;
load('j2.mat');
J2=J2k;

load('psi.mat');
psiS=psik;

load('epsk.mat');
epsS=epsk;



figure(1)
plot(time_vec(2:end-1), u1S(1,2:end-1), ltype{1}, 'LineWidth', 2); hold on;
plot(time_vec(2:end-1), u1S(2,2:end-1), ltype{3}, 'LineWidth', 2);
xlabel('t')
%ylabel('u_1 ')
legend({'$u_{11}$','$u_{12}$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'control1SQHexp4.eps','-b0'); 
print('-dpdf', 'control1SQHexp4.pdf','-b0');

figure(2)
plot(time_vec(2:end-1), u2S(1,2:end-1), ltype{1}, 'LineWidth', 2); hold on;
plot(time_vec(2:end-1), u2S(2,2:end-1), ltype{3}, 'LineWidth', 2);
xlabel('t')
%ylabel('u_2')
legend({'$u_{21}$','$u_{22}$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'control2SQHexp4.eps','-b0'); 
print('-dpdf', 'control2SQHexp4.pdf','-b0');

figure(3)
plot(J1,ltype{1},'Linewidth',2); hold on; 
plot(J2,ltype{2},'Linewidth',2)
legend({'$J_{1}$','$J_{2}$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'J12SQHexp4.eps','-b0'); 
print('-dpdf', 'J12SQHexp4.pdf','-b0');

figure(4)
semilogy(psiS,ltype{1},'Linewidth',2);

print('-depsc2', 'psiSQHexp4.eps','-b0'); 
print('-dpdf', 'psiSQHexp4.pdf','-b0');

figure(5)
plot(time_vec ,y(1,:),ltype{1},'Linewidth',2); hold on; 
plot(time_vec ,y(2,:),ltype{2},'Linewidth',2); 
legend({'$y_{1}$','$y_{2}$'},'Interpreter','Latex','FontSize',12)

xlabel('t')
%ylabel('y')

print('-depsc2', 'ySQHexp4.eps','-b0'); 
print('-dpdf', 'ySQHexp4.pdf','-b0');

figure(6)
semilogy(epsS,ltype{1},'Linewidth',2); 

print('-depsc2', 'epsSQHexp4.eps','-b0'); 
print('-dpdf', 'epsSQHexp4.pdf','-b0');


