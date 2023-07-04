%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% Glider - solution by SA method & outer loop for T
%
close all; 
clear all;

ltype = {'b-','r--','m-x','b-*','r:','m-.'};

OCP=struct('Nt',1000,'T',0.5);
DYN=struct('ND',4,'sigmaKC',0.5,'b',0.2,'alpha0',0.175,'K',0.1);

% parameters for augmented lagrangian method 
global A lambda
A = 1.0e+3; %B
lambda = 0.0;

% Initial conditions 
alpha0=DYN.alpha0;
% x,y,v,theta
y0=[0;0;1.0;0.25];

N=DYN.ND;
K=DYN.K;
y=zeros(N,OCP.Nt+1);
p=zeros(N,OCP.Nt+1);
u=zeros(2,OCP.Nt+1);

% Desired trajectory/terminal configuration
% yd=[0;0;0;0];

% Numerical parameters
Nt=OCP.Nt;
Nt1=Nt+1;

T=OCP.T;
dt=T/Nt;

TT = T; 

     
%Algorithm paramters
mmax=500;        %Maximum iteration number
kmax=10;        %Maximum iteration number
kappa=10^-12;        %Tolerance for convergence


% Two controls: alpha=u(1) and eta=u(2)
u=0.001*ones(2,Nt1);
u_old=u;

% 
y=forward(y0,u,TT,OCP,DYN);
p=backward(y,u,TT,OCP,DYN);
J_int=get_J(y,u,OCP,Nt1);

HXE=HPfunction(y(:,end),p(:,end),[u(1,end),u(2,end)],OCP,DYN);

fprintf('TT= %8.4f x= %8.4e y= %8.4e Jo=%12.6e \n',TT,y(1,Nt1),y(2,Nt1),J_int)


% step size for outer loop for T
st = 1.0e-4;
for m=1:mmax

HXEold = HXE;
fprintf('Outer loop m = %4i lambda= %8.4f HX=%8.4e TT= %8.4e \n',m,lambda,HXE,TT)

% SA loop for T fixed


for k=1:kmax
   
    for i=1:Nt
        
        u(1,i)=0.5*atan(tan(2*alpha0)*p(4,i)*K/(y(3,i)*(p(3,i))));
        
        v1 = 0.0;
        H1=HPfunction(y(:,i),p(:,i),[u(1,i),v1],OCP,DYN);
        v2 = 1.0;
        H2=HPfunction(y(:,i),p(:,i),[u(1,i),v2],OCP,DYN);
        
        u(2,i)=0.0;
        if (H2 > H1) 
            u(2,i)=1.0;
        end
    end 
    % Extend u to Nt1 for midpoint integration 
    u(:,Nt1)=u(:,Nt);   
    
    du=sum((u-u_old).^2,'all')*dt;
    u_old=u;
   
    y=forward(y0,u,TT,OCP,DYN);
    p=backward(y,u,TT,OCP,DYN);   
    J_int=get_J(y,u,OCP,Nt1);
    
  fprintf('k = %4i x= %8.4e y= %8.4e J= %12.6e \n',k,y(1,Nt1),y(2,Nt1),J_int)
  
     if(du<kappa)
        fprintf('converged at k = %4i \n',k)
        break;
   end 
    
end
% SA loop for T fixed - end

 Jk(m)=J_int;
    
 HXE=HPfunction(y(:,end),p(:,end),[u(1,end),u(2,end)],OCP,DYN);  
 
 lambda = lambda + A*y(2,end);
 
 lam(m)=lambda;
 
 % update TT
 TT = TT + st*HXE;
 dt = TT/Nt;
 
 if HXE*HXEold < 0
     disp('H changes sign');
 end
 HXEold=HXE;
 
 % print the actual p(T) and the one required for H(T)=0
 fprintf('pTy=%8.4e pTth=%8.4e \n ',(- A*y(2,Nt1) - lambda),(- cot(y(4,Nt1))))
 fprintf('\n')
 
 if abs(HXE) < 1.0e-2
     disp('H suff. small');
     break;
 end
   
end

    for ii = 1:Nt1    
    HX(ii)=HPfunction(y(:,ii),p(:,ii),[u(1,ii),u(2,ii)],OCP,DYN);
    end
 




%%%%%%%%%  % y = x,y,v,theta

figure(1)
plot(0:dt:TT,y(1,1:Nt1),ltype{1},'Linewidth',2)
xlabel('t')
ylabel('x')

print('-depsc2', 'xKC.eps','-b0'); 
%print('-dpdf', 'controlKC.pdf','-b0');

figure(2)
plot(0:dt:TT,y(2,1:Nt1),ltype{1},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('y')
print('-depsc2', 'yKC.eps','-b0'); 
%print('-dpdf', 'trajectoryKC.pdf','-b0');

figure(3)
plot(0:dt:TT,y(3,1:Nt1),ltype{1},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('v')

print('-depsc2', 'vKC.eps','-b0'); 
%print('-dpdf', 'adjtrajKC.pdf','-b0');

figure(4)
plot(0:dt:TT,y(4,1:Nt1),ltype{1},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('$\theta$','Interpreter','Latex','FontSize',12)

print('-depsc2', 'thKC.eps','-b0'); 
%print('-dpdf', 'adjtrajKC.pdf','-b0');


%%%%%%%%%


figure(5)
plot(0:dt:TT,u(1,1:Nt1),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:TT,u(2,1:Nt1),ltype{2},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('$\alpha , \, \eta$','Interpreter','Latex','FontSize',12)
legend({'$\alpha$','$\eta$'},'Interpreter','Latex','FontSize',12)
print('-depsc2', 'controlKC.eps','-b0'); 
% %print('-dpdf', 'controlKC.pdf','-b0');

% 
figure(6)
plot(0:dt:TT,p(1,1:Nt1),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:TT,p(2,1:Nt1),ltype{2},'Linewidth',2),
plot(0:dt:TT,p(3,1:Nt1),ltype{3},'Linewidth',2),
plot(0:dt:TT,p(4,1:Nt1),ltype{4},'Linewidth',2),
%axis([-inf inf -1 1])
xlabel('x/T')
ylabel('p')
legend({'$p_1$','$p_2$','$p_3$','$p_4$'},'Interpreter','Latex','FontSize',12)
% %print('-depsc2', 'adjtrajKC.eps','-b0'); 
% %print('-dpdf', 'adjtrajKC.pdf','-b0');
% 
figure(7)
plot(Jk,ltype{2},'Linewidth',2)
xlabel('iterations')
ylabel('J')
print('-depsc2', 'JKC.eps','-b0'); 

figure(8)
plot(0:dt:TT,HX,ltype{2},'Linewidth',2)
xlabel('t')
ylabel('H')
print('-depsc2', 'HKC.eps','-b0'); 

figure(8)
plot(lam,ltype{1},'Linewidth',2)
xlabel('iterations')
ylabel('\lambda')
print('-depsc2', 'lambdaKC.eps','-b0'); 


