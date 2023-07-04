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

N=2501;
T=0.25;
time_vec = linspace( 0.0 , T , N );


load('u.mat');
u1S= u1;
u2S= u2;

load('y.mat');
yS=y;


figure(1)
plot(time_vec(2:end-1), u1S(1,2:end-1), ltype{1}, 'LineWidth', 2); hold on;
plot(time_vec(2:end-1), u1S(2,2:end-1), ltype{3}, 'LineWidth', 2);
xlabel('t')
%ylabel('u_1 ')
legend({'$u_{11}$','$u_{12}$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'control1SQHexp3.eps','-b0'); 
print('-dpdf', 'control1SQHexp3.pdf','-b0');

figure(2)
plot(time_vec(2:end-1), u2S(1,2:end-1), ltype{1}, 'LineWidth', 2); hold on;
plot(time_vec(2:end-1), u2S(2,2:end-1), ltype{3}, 'LineWidth', 2);
xlabel('t')
%ylabel('u_2')
legend({'$u_{21}$','$u_{22}$'},'Interpreter','Latex','FontSize',12)


print('-depsc2', 'control2SQHexp3.eps','-b0'); 
print('-dpdf', 'control2SQHexp3.pdf','-b0');


