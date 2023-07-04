%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for optimal control of the Bloch system
%
close all; 
clear all;

ltype = {'b-','b--','m-.','b-*','r:','m-x'};


%
%exp1
OCP=struct('Nt',1000,'T',1.0,'umin',-2,'umax',2,'nu',10^-7,'beta',10^-7);

% exp2
% OCP=struct('Nt',1000,'T',0.5,'umin',-2,'umax',2,'nu',10^-7,'beta',10^-3);


% initial condition 
y0=[0.0;0.0;1.0];
y0=y0/norm(y0);

% gg 
gg = 0.0;  % exp1 =0 
om = 1.0;
aa = 1.0;
dd = gg;

A=[-gg/2,-om,0;om,-gg/2,0; 0,0,-gg/2];

B1=[0,0,-aa;0,0,0;aa,0,0];
B2=[0,0,0;0,0,-aa;0,aa,0];

Dv=[0.0;0.0;-dd]; % D


% target and weight
nu=OCP.nu;
yd=[1.0;0.0;0.0];
yd=yd/norm(yd);

% admissible control values [umin,umax]
umin =OCP.umin;
umax =OCP.umax;
ux = [umin:0.1:umax]; % for search 

% time horizon and grid
T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

% SQH parameters
kmax=1000;        %Maximum iteration number
kappa=10^-10;              %Tolerance for convergence  
%Algorithm paramters
eta=10^-9;
zeta=0.8;
sigma=1.2;
epsilon=1.;      %Initial value for epsilon



% initial control 
u=0*ones(2,Nt+1);
u_old=u;

y=forward(A,B1,B2,Dv,y0,u,OCP);
p=backward(A,B1,B2,y,yd,u,OCP);
J_int=get_J(y,u,yd,OCP);

tic
for k=1:kmax
    
% save the value of the functional

    
y=forward(A,B1,B2,Dv,y0,u_old,OCP);    
p=backward(A,B1,B2,y,yd,u_old,OCP);
J_int=get_J(y,u_old,yd,OCP);
Jk(k)=J_int;


iFlag  = 1; 
epsMax = 1.0e+15;
%icount = 0;
while ( iFlag == 1 && epsilon <  epsMax )
    
%icount = icount + 1; 
% update control     
    for i=1:Nt+1
                

        uz1=max(min(OCP.umax,(2*epsilon*u_old(1,i)+(p(:,i)')*B1*y(:,i)-OCP.beta)/(2*epsilon+OCP.nu)),0);
        uz2=max(min(OCP.umax,(2*epsilon*u_old(2,i)+(p(:,i)')*B2*y(:,i)-OCP.beta)/(2*epsilon+OCP.nu)),0);
               
        us1=max(min(0,(2*epsilon*u_old(1,i)+(p(:,i)')*B1*y(:,i)+OCP.beta)/(2*epsilon+OCP.nu)),OCP.umin);
        us2=max(min(0,(2*epsilon*u_old(2,i)+(p(:,i)')*B2*y(:,i)+OCP.beta)/(2*epsilon+OCP.nu)),OCP.umin);
        
        uu1 = [uz1;uz2]' ; 
        H(1)=HPfunction(A,B1,B2,Dv,y(:,i),p(:,i),uu1,OCP) ... 
          - epsilon*(uz1-u_old(1,i))^2  - epsilon*(uz2-u_old(2,i))^2 ;  
      
        uu2 = [us1;us2]' ; 
        H(2)=HPfunction(A,B1,B2,Dv,y(:,i),p(:,i),uu2,OCP) ... 
          - epsilon*(us1-u_old(1,i))^2  - epsilon*(us2-u_old(2,i))^2 ;  

        uu3 = [uz1;us2]' ; 
        H(3)=HPfunction(A,B1,B2,Dv,y(:,i),p(:,i),uu3,OCP) ... 
          - epsilon*(uz1-u_old(1,i))^2  - epsilon*(us2-u_old(2,i))^2 ;  

      
        uu4 = [us1;uz2]' ; 
        H(4)=HPfunction(A,B1,B2,Dv,y(:,i),p(:,i),uu4,OCP) ... 
          - epsilon*(us1-u_old(1,i))^2  - epsilon*(uz2-u_old(2,i))^2 ;  

      
        [~,pos]=max(H);

     
        if pos == 1
        u(1,i)=uz1;
        u(2,i)=uz2;
        end

        if pos == 2
        u(1,i)=us1;
        u(2,i)=us2;
        end
        
        if pos == 3
        u(1,i)=uz1;
        u(2,i)=us2;
        end        

        if pos == 4
        u(1,i)=us1;
        u(2,i)=uz2;
        end        

    end 
  
   du=sum((u-u_old).^2,'all')*dt;
   
   y=forward(A,B1,B2,Dv,y0,u,OCP);
   Jt=get_J(y,u,yd,OCP);
   
   if(Jt-Jk(k)>-eta*du)
       epsilon=epsilon*sigma;    
       iFlag = 1; 
   else
       epsilon=epsilon*zeta;    
       u_old=u;
       iFlag = 0; 
   end
   
end


  fprintf('k = %4i du = %10.8e J = %10.8e \n',k,du,Jk(k))

  if(epsilon > epsMax) 
  fprintf('epsilon too large\n')
     break; 
  end
  if(du<kappa)
        fprintf('converged\n')
        break;
  end 
end 

% check norm of solution initial and final time
% ny0=y(:,1)'*y(:,1)
% nyf=y(:,end)'*y(:,end)

toc


% ------------------------ 

figure(1)
plot(0:dt:T,u(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:T,u(2,:),ltype{2},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('u')
legend({'$u_1$','$u_2$'},'Interpreter','Latex','FontSize',12)
print('-depsc2', 'BlochCu01.eps','-b0'); 
% print('-dpdf', 'BlochCu01.pdf','-b0');

figure(2)
plot(0:dt:T,y(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:T,y(2,:),ltype{2},'Linewidth',2),
plot(0:dt:T,y(3,:),ltype{3},'Linewidth',2),
axis([-inf inf -1 1])
xlabel('t')
ylabel('y')
legend({'$m_x$','$m_y$','$m_z$'},'Interpreter','Latex','FontSize',12)
print('-depsc2', 'trajectoryBloch01.eps','-b0'); 
% print('-dpdf', 'trajectoryBloch01.pdf','-b0');


figure(3)
plot(Jk,ltype{2},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J')
print('-depsc2', 'JBloch01.eps','-b0'); 

