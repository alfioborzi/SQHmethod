%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method displaying chattering control
%
close all; 
clear all;

ltype = {'b-','b--','m-.','b-*','r:','m-x'};

OCP=struct('Nt',2000,'T',1.0,'umin',-1,'umax',1,'nu',-1.0,'beta',0.0);

% initial condition 
y0=[0.0];

% parameters
aa = 0.0;
bb = 1.0;

A=[aa];
B=[bb];


% admissible control values [umin,umax]
umin =OCP.umin;
umax =OCP.umax;
ux = [umin:0.2:umax] ; % defines the discrete set in each component
Nu = length(ux) ;

% time horizon and grid
T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

% desired trajectory 
yd=0*ones(1,Nt+1);



% SQH parameters - orig
kmax=1000;        %Maximum iteration number
kappa=10^-12;              %Tolerance for convergence  
%Algorithm paramters
eta=10^-12;
zeta=0.8;
sigma=1.01;
epsilon=0.1;      %Initial value for epsilon


% initial control 
u=0.0*ones(1,Nt+1);
u_old=u;

y=forward(A,B,y0,u,OCP);
p=backward(A,B,y,yd,u,OCP);
J_int=get_J(y,u,yd,OCP);

tic
for k=1:kmax
    
% save the value of the functional

    
y=forward(A,B,y0,u_old,OCP);    
p=backward(A,B,y,yd,u_old,OCP);
J_int=get_J(y,u_old,yd,OCP);
Jk(k)=J_int;


iFlag  = 1; 
epsMax = 1.0e+15;

while ( iFlag == 1 && epsilon <  epsMax )
    
% update control     
    for i=1:Nt+1
                
        HPMax = -1.0e+10; 
        % search for max in the discrete set 
        for ii=1:Nu
                
        uus = ux(ii) ; 
        HP=HPfunction(A,B,y(:,i),p(:,i),uus,OCP) ... 
          - epsilon*(ux(ii)-u_old(1,i))^2  ;  
      
        if ( HP > HPMax ) 
           HPMax = HP;
           iim=ii;  
        end
      
        end
      
        u(1,i)=ux(iim);            

    end 
  
   du=sum((u-u_old).^2,'all')*dt;
   
   y=forward(A,B,y0,u,OCP);
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

fprintf('k = %4i tau = %10.8e eps = %10.8e J = %10.8e \n',k,du,epsilon,Jk(k))

  if(epsilon > epsMax) 
  fprintf('epsilon too large\n')
     break; 
  end
  if(du<kappa)
        fprintf('converged\n')
        break;
  end 
end 



toc


% ------------------------ 

figure(1)
plot(0:dt:T,u(1,:),ltype{1},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('u')

print('-depsc2', 'chatteringU01.eps','-b0'); 

figure(2)
plot(0:dt:T,y(1,:),ltype{1},'Linewidth',2)
axis([-inf inf -1 1])
xlabel('t')
ylabel('y')

print('-depsc2', 'chatteringTraj01.eps','-b0'); 


figure(3)
plot(Jk,ltype{1},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J')
print('-depsc2', 'chatteringJ01.eps','-b0'); 

