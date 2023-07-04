%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% solving an elliptic optimal control problem with L^1 tracking term 
% and max(y,0) nonlinearity in the state equation on a square [0,a]x[0,b] 
%
% A robust SQH scheme to solve non-smooth PDE optimal control problems
% by Tim Breitenbach and Alfio Borzi'
%
% Problem specific parameters: 
% (a,b) interval length of the domain, 
% u_lo and u_up lower and upper bounds for the control, 
% alpha (weight quadratic cost term), 
% beta weight for the non-convex cost term, 
% N+1 number of intervals for the space discretisation, 
% kmax maximum number of iterations 

close all;
clear all; 

ltype = {'b-','r--','m-x','b-*','r:','m-.'};

OCP=struct('N',99,'a',0,'b',1,'u_lo',-40,'u_up',40,'alpha',0*10^-3,'beta',10^-5,'s',0);

%Tolerance for convergence
kappa=10^-6;                
%Algorithm paramters
eta=10^-7;
zeta=0.15;
sigma=10.;
kmax=5000;           %Maximum iteration number
epsilon=1;        %Initial guess for epsilon

h=1/(OCP.N+1); %Step size for the discretisation
[X1,X2]=meshgrid(OCP.a:h:OCP.b,OCP.a:h:OCP.b); %Meshgrid for the plots
[X1u,X2u]=meshgrid(OCP.a+h:h:OCP.b-h,OCP.a+h:h:OCP.b-h); %Meshgrid for the plots of the control (only interior points)

%Assemble the discretised Laplacian
n=OCP.N+1;              
N=(n+1)*(n+1)-4*n;
v=ones(N,1);
v_block=[ones(n-2,1);0];
v_block=kron(ones(N/(n-1),1),v_block); 
A_hilf=spdiags([(-1/h^2)*v_block (-1/h^2)*[v_block(N); v_block(1:N-1)]],[ -1  1 ],N,N);
A=spdiags([(-1/h^2)*v  (4/h^2)*v  (-1/h^2)*v],[-(n-1)  0  n-1],N,N);
A=A+A_hilf;

y=zeros(OCP.N+2,OCP.N+2);
yd=sin(2*pi*X1).*sin(2*pi*X2); %Desired state


u=zeros(OCP.N,OCP.N);  %Initial guess for the control 
y=forward_y(y,u,A,OCP);    %Calculate the corresponding state variable y
p=backward(y,yd,A,OCP);  %Calculate the corresponding adjoint state variable
u_old=u;
y_old=y;

J(1)=get_Jy(u,y,yd,OCP);%Calculate the corresponding cost functional value

count_updates=1;        %Counter for the updates; counter_updates -1 for the number of updates
time_calc=cputime;      %Start time measurement

for k=1:kmax
    %Loop over all grid points
    for i=1:OCP.N  
       for j=1:OCP.N
            %Pointwise minimisation of the corresponding augmented Hamiltonian
            u(i,j)=argMinH(u(i,j),p(i+1,j+1),epsilon,OCP); 
        end
    end
    y=forward_y(y,u,A,OCP);             %Calculation of the state equation with the updated control u
    J_int=get_Jy(u,y,yd,OCP);       %Corresponding cost functional value
    du=sum(sum((u_old-u).^2))*h*h;  %Evalute the norm square of the update
    if (J_int-J(count_updates)>-eta*du) 
        % If the improvemnt to the cost functional is not sufficiently large, 
        % take the old values of y and u, increase epsilon
        y=y_old;
        u=u_old;
        epsilon=epsilon*sigma;
    else   
        % If the improvment to the cost functional is sufficently large, 
        % take the new values of y, u, J, calculate p and decrease epsilon
        epsilon=epsilon*zeta;
        count_updates=count_updates+1;
        p=backward(y,yd,A,OCP);
        J(count_updates)=J_int;
        u_old=u;
        y_old=y; 
        % Print current values
          fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',k,... 
            count_updates-1,J_int,epsilon,du) 
    end
        
    if(du<kappa)   
        % If the norm square of the update is too small, then return the latest 
        % variable u that caused sufficient decrease of the cost functional 
        fprintf('Convergence achieved\n');
        u=u_old;
        break;
    end    
end
ct=cputime-time_calc;            %Stop time measurment for the calculation 
fprintf('CPU-time in seconds needed: %f\n',ct)

%Plot y and u
figure(1)
mesh(X1,X2,y);
  %title('The state variable');
  xlabel('$x_1$','interpreter','latex');
  ylabel('$x_2$','interpreter','latex');

  print('-depsc2', 'ellP4y.eps','-b0'); 
% print('-dpdf', 'ellP4y.pdf','-b0');

figure(2)
mesh(X1u,X2u,u);
  %title('The control variable');
  xlabel('$x_1$','interpreter','latex');
  ylabel('$x_2$','interpreter','latex');

  print('-depsc2', 'ellP4u.eps','-b0'); 
% print('-dpdf', 'ellP4u.pdf','-b0');

  
figure(3)
contour(X1u,X2u,u);
  %title('Contour of the control variable');
  xlabel('$x_1$','interpreter','latex');
  ylabel('$x_2$','interpreter','latex');
  
  print('-depsc2', 'ellP4ucontour.eps','-b0'); 
% print('-dpdf', 'ellP4ucontour.pdf','-b0');
  
%Plot convergence history
figure(4)
semilogx(J,ltype{3},'Linewidth',2)
% title('Convergence history of J');

 print('-depsc2', 'ellP4Jhist.eps','-b0'); 
% print('-dpdf', 'ellP4Jhist.pdf','-b0');


fprintf('\n')
fprintf('Check numerical optimality of the SQH solution\n')
% Check solution for optimality; 
numeric_optimality(u,y,yd,A,OCP)

