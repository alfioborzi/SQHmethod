%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for optimal control of Block system and path constraints
%

close all; 
clear all;

ltype = {'b-','b--','m-.','b-*','r:','m-x'};

%see sparsity
OCP=struct('Nt',500,'T',1.0,'umin',-2,'umax',2,'nu',10^-7, ... 
           'beta',10^-7,'alpha',10^-2);


% initial condition 
y0=[0.0;0.0;1.0];
y0=y0/norm(y0);

% gg = 0.001;
gg = 0.0;
om = 1.0;
aa = 1.0;

A=[-gg/2,-om,0;om,-gg/2,0; 0,0,-gg/2];

B1=[0,0,-aa;0,0,0;aa,0,0];
B2=[0,0,0;0,0,-aa;0,aa,0];

gv=[0.0;0.0;-gg];

% target and weight

yd=[1.0;0.0;0.0];
yd=yd/norm(yd);

% admissible control values [umin,umax] & grid 
umin =OCP.umin;
umax =OCP.umax;
ux = [umin:0.1:umax]; % for search 
Nu = length(ux) ;

alpha = OCP.alpha;
beta  = OCP.beta;
nu    = OCP.nu;

% time horizon and grid
T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

% SQH parameters
kmax=1000;                %Maximum iteration number
kappa=10^-8;              %Tolerance for convergence  
%Algorithm paramters
eta=10^-7;
zeta=0.9;
sigma=1.1;
epsilon=1.;      %Initial value for epsilon


% initial control 
u=0.0*ones(2,Nt+1);
u_old=u;

% for path constraint initialise ct = 0
ct=zeros(1,Nt+1);

kount = 0;
tic 

% outer loop 
for kk =1:10
    fprintf('Outer loop kk = %4i \n',kk)
    % plot
    
    
figure(4)
plot(0:dt:T,ct(:),ltype{1},'Linewidth',2); 
xlabel('t')
ylabel('c')
drawnow;
kp = 2;
if (kk == kp) % plot at kp iteration
    print('-depsc2', 'adjointMix01.eps','-b0'); 
end



y=forward(A,B1,B2,gv,y0,u,OCP);
p=backward(A,B1,B2,y,yd,u,ct,OCP);
J_int=get_J(y,u,ct,yd,OCP);


% check path constraint 
for i=1:Nt+1
pgc(i)= gc( y(:,i), u(:,i));
end
figure(5)
plot(0:dt:T, pgc(:),ltype{1},'Linewidth',2) ;
% title('path constraint g(y,u)');
xlabel('t')
ylabel('\phi')
drawnow;



% SQH loop 
for k=1:kmax
       
% save the value of the functional
    
y=forward(A,B1,B2,gv,y0,u_old,OCP);    
p=backward(A,B1,B2,y,yd,u_old,ct,OCP);
J_int=get_J(y,u_old,ct,yd,OCP);

kount = kount + 1;
Jk(kount)=J_int;

iFlag  = 1; 
epsMax = 1.0e+15;
%icount = 0;
while ( iFlag == 1 && epsilon <  epsMax )
    
% sqh update control by pre determined choice
    for i=1:Nt+1
        
    dldu = alpha*max(0,gc(y(:,i),u(:,i))+ct(1,i)/alpha )*dgdu(y,u)';
                
    uz1=max(min(OCP.umax,(2*epsilon*u_old(1,i)+(p(:,i)')*B1*y(:,i)-dldu(1)-beta)/(2*epsilon+nu)),0);
    uz2=max(min(OCP.umax,(2*epsilon*u_old(2,i)+(p(:,i)')*B2*y(:,i)-dldu(2)-beta)/(2*epsilon+nu)),0);
               
    us1=max(min(0,(2*epsilon*u_old(1,i)+(p(:,i)')*B1*y(:,i)-dldu(1)+beta)/(2*epsilon+nu)),OCP.umin);
    us2=max(min(0,(2*epsilon*u_old(2,i)+(p(:,i)')*B2*y(:,i)-dldu(2)+beta)/(2*epsilon+nu)),OCP.umin);
        
        uu1 = [uz1;uz2]' ; 
        H(1)=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uu1,ct(i),OCP) ... 
          - epsilon*(uz1-u_old(1,i))^2  - epsilon*(uz2-u_old(2,i))^2 ;  
      
        uu2 = [us1;us2]' ; 
        H(2)=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uu2,ct(i),OCP) ... 
          - epsilon*(us1-u_old(1,i))^2  - epsilon*(us2-u_old(2,i))^2 ;  

        uu3 = [uz1;us2]' ; 
        H(3)=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uu3,ct(i),OCP) ... 
          - epsilon*(uz1-u_old(1,i))^2  - epsilon*(us2-u_old(2,i))^2 ;  

      
        uu4 = [us1;uz2]' ; 
        H(4)=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uu4,ct(i),OCP) ... 
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
    % end update step 
    
%   % alternative : update control by inspection
%     for i=1:Nt+1
%                 
%         HPMax = -1.0e+10; 
%         % search for max in the discrete set 
%         for ii=1:Nu
%         for jj=1:Nu
%                 
%         uus = [ux(ii);ux(jj)]' ; 
%         HP=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uus,ct(i),OCP) ... 
%           - epsilon*(ux(ii)-u_old(1,i))^2  - epsilon*(ux(jj)-u_old(2,i))^2;  
%       
%         if ( HP > HPMax ) 
%            HPMax = HP;
%            iim=ii; jjm=jj; 
%         end
%       
%         end
%         end
%       
%         u(1,i)=ux(iim);
%         u(2,i)=ux(jjm);        
%     end   
%     
  
   du=sum((u-u_old).^2,'all')*dt;
   
   y=forward(A,B1,B2,gv,y0,u,OCP);
   Jt=get_J(y,u,ct,yd,OCP);
   
   if(Jt-Jk(k)>-eta*du)
       epsilon=epsilon*sigma;    
       iFlag = 1; 
   else
       epsilon=epsilon*zeta;    
       u_old=u;
       iFlag = 0; 
   end
   
end


  if(epsilon > epsMax) 
  fprintf('epsilon too large\n')
     break; 
  end
  if(du<kappa)
        fprintf('converged\n')
        break;
  end 
end 

fprintf('SQH: k = %4i du = %10.8e J = %10.8e \n',k,du,Jk(k))




% update ct & compute for stopping
ctx = ct;
ncstop = 0.;
for j=1:Nt+1
ctx(j) = max(0, alpha*gc(y(:,j),u(:,j)) + ct(1,j) );
ncstop = ncstop + dt*abs(min(-gc(y(:,j),u(:,j)),ct(1,j)));
end
ct = ctx;
  
  if(ncstop<kappa)
        fprintf('outer loop converged with tol = %10.8e\n',ncstop)
        break;
  end 

end

toc

% ------------------------ 

figure(1)
plot(0:dt:T,u(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:T,u(2,:),ltype{2},'Linewidth',2)
%axis([-inf inf -1 1])
xlabel('t')
ylabel('u')
legend({'$u_1$','$u_2$'},'Interpreter','Latex','FontSize',12)
print('-depsc2', 'MixCu01.eps','-b0'); 
% print('-dpdf', 'MixCu01.pdf','-b0');

figure(2)
plot(0:dt:T,y(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:dt:T,y(2,:),ltype{2},'Linewidth',2),
plot(0:dt:T,y(3,:),ltype{3},'Linewidth',2),
axis([-inf inf -1 1])
xlabel('t')
ylabel('y')
legend({'$m_x$','$m_y$','$m_z$'},'Interpreter','Latex','FontSize',12)
print('-depsc2', 'trajectoryMix01.eps','-b0'); 
% print('-dpdf', 'trajectoryMix01.pdf','-b0');


figure(3)
plot(Jk,ltype{2},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J')
print('-depsc2', 'mixedJ01.eps','-b0'); 

figure(4)
plot(0:dt:T,ct(:),ltype{1},'Linewidth',2)
xlabel('t')
ylabel('c')

% check path constraint 
for i=1:Nt+1
pgc(i)= gc( y(:,i), u(:,i));
end
figure(5)
plot(0:dt:T, pgc(:),ltype{1},'Linewidth',2) ;
% title('path constraint g(y,u)');
xlabel('t')
ylabel('\phi')
print('-depsc2', 'pathconstrMix01.eps','-b0'); 



