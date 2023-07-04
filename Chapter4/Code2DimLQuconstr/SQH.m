%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for differential Nash games
%

close all;
clear all;


ltype = {'b-','r--','m-.','b-*','r:','m-x'};
I = eye(2,2);


OCP=struct('dt',1e-4,'T',0.250,'nu1',0.01*I,'gamma1', 0.01*I,'nu2',0.01*I, ... 
    'gamma2', 0.01*I, 'alpha1', 0.1*I, 'alpha2', 1*I , ... 
    'A',[1, 0; 0,2],'B1',[1,0;0,-1],'B2',[2,-1;0,2], 'y0',[2;1], ...
    'umax',3,'umin',-3);



y0 =   OCP.y0;


A  =   OCP.A;
B1 =   OCP.B1;
B2 =   OCP.B2;

yd1 = [0;0];
yd2 = [0;0];


Nt=round(OCP.T/OCP.dt);

u1=0.0*ones(2,Nt+1);
u2=0.0*ones(2,Nt+1);


%Tolerance for convergence
kappa=10^-12;                
%Algorithm paramters
tol_psi = 10^-8;

zeta=0.95;          %original values: zeta=0.8; sigma=2;
sigma=1.05;


kmax=10000;        %Maximum iteration number

epsilon0 = 10.0;   %Initial guess for epsilon
epsilon = epsilon0;   %Initial guess for epsilon


u1_old=u1;
u2_old=u2;
du=1.0;

y=forward(A,B1,B2,y0,u1,u2,OCP);
y_old=y;

p1=backward1(A,B1,B2,y,yd1,u1,u2,OCP);
p2=backward2(A,B1,B2,y,yd2,u1,u2,OCP);


J1_plot(1,1)=get_J1(y,u1,yd1,OCP);
J1k(1)=J1_plot(1,1);

J2_plot(1,1)=get_J2(y,u2,yd2,OCP);
J2k(1)=J2_plot(1,1);

psi_old= 10;
count_updates=0;

tic
for k=1:kmax
    
    for i=1:Nt+1
        
        %Player 1
        uz1 = max(min( (2*epsilon*I + OCP.nu1)\(B1'*p1(:,i) + 2*epsilon*u1_old(:,i)), [OCP.umax; OCP.umax]), [OCP.umin;OCP.umin]);
        u1(:,i)=uz1;
      
        %Player 2
        uz2 = max(min((2*epsilon*I + OCP.nu2)\(B2'*p2(:,i) + 2*epsilon*u2_old(:,i)), [OCP.umax;OCP.umax]), [OCP.umin;OCP.umin]);
        u2(:,i)=uz2;
        
    end 

   du1=sum(sum((u1-u1_old).^2))*OCP.dt;
   du2=sum(sum((u2-u2_old).^2))*OCP.dt;
   du = sum (sum((u1-u1_old).^2) + sum((u2-u2_old).^2))*OCP.dt;
    
   y=forward(A,B1,B2,y0,u1,u2,OCP);
   
   J1_int=get_J1(y,u1,yd1,OCP); 
   J2_int=get_J2(y,u2,yd2,OCP); 
   
    
    %Nikaido-Isoda function
    y1 = forward(A,B1,B2,y0,u1_old,u2, OCP);
    y2 = forward(A,B1,B2,y0,u1,u2_old, OCP);
    J1_ni = get_J1(y1,u1_old,yd1,OCP);
    J2_ni = get_J2(y2,u2_old,yd2,OCP);
    
    
    psi = J1_int + J2_int - J1_ni - J2_ni;
   
       if (psi < -tol_psi*du) && (abs(psi) < abs(psi_old))

       count_updates=count_updates+1;
      
        epsilon=epsilon*zeta;   
       
       J1k(count_updates)=J1_int;
       J2k(count_updates)=J2_int;
          
       u1_old=u1;
       u2_old=u2;
       y_old=y;
       
         psi_old = psi;
         
       p1=backward1(A,B1,B2,y,yd1,u1,u2,OCP);
       p2=backward2(A,B1,B2,y,yd2,u1,u2,OCP);
       
   fprintf('| iter=%4i | J1= %10.8f | J2=%10.8f  | psi=%10.6e | du=%10.4e | eps=%10.4e |\n',k, J1_int, J2_int,psi,du,epsilon);
       
       else 
       u1=u1_old;
       u2=u2_old;
       y=y_old;
       epsilon=epsilon*sigma;      

   end


   if ( (du < kappa) && ( psi < 0 ) )
        u1=u1_old;
        u2=u2_old;
        fprintf('converged\n')
        break;
   end 
      
end 

fprintf('Number of successful updates = %4i \n', count_updates)
toc

time_vec=linspace(0,OCP.T, Nt+1);

figure(1)
plot(time_vec, u1, ltype{1}, 'Linewidth',2)
xlabel('x')
ylabel('u_1')

figure(2)
plot(time_vec, u2, ltype{1}, 'Linewidth',2)
xlabel('x')
ylabel('u_2')

figure(3)
plot(time_vec ,y,ltype{1},'Linewidth',2)
xlabel('x')
ylabel('y')

figure(10)
plot(time_vec ,p1,ltype{1},'Linewidth',2)
xlabel('x')
ylabel('p_1')


figure(11)
plot(time_vec ,p2,ltype{1},'Linewidth',2)
xlabel('x')
ylabel('p_2')


figure(5)
plot(J1k,ltype{2},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J_1')

figure(6)
plot(J2k,ltype{2},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J_2')


save('u.mat','u1','u2')
save('y.mat','y')


