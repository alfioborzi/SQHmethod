%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% solving a parabolic bilinear optimal control problem on the
% space-time cylinder [a,b]x[0,T]
%
% A robust SQH scheme to solve non-smooth PDE optimal control problems
% by Tim Breitenbach and Alfio Borzi'
%
% Problem specific parameters: 
% N number of intervals for the discretisation of the space, 
% Nt for the discretation of the time, 
% (a,b) interval length of the domain, T time horizon, D diffusion constant, 
% u_lo and u_up lower and upper bound for the control, 
% alpha (weight quadratic cost term), 
% beta and s (weight discontinuous cost term) cost functional parameters

close all;
clear all; 

ltype = {'b-','r--','m-x','b-*','r:','m-.'};


OCP=struct('N',100,'Nt',200,'a',0,'b',1,'T',1,'D',1.0,'u_lo',-30,'u_up',50,'alpha',10^-6,'beta',10^-6,'s',20); 

%Tolerance for convergence
kappa=10^-8;                
%Algorithm paramters
eta=10^-5;
zeta=0.9;
sigma=1.1;
kmax=5000;         %Maximum iteration number
epsilon=1.0;         %Initial guess for epsilon

[X,T]=meshgrid(OCP.a:(OCP.b-OCP.a)/OCP.N:OCP.b,0:OCP.T/OCP.Nt:OCP.T); %Meshgrid for the plots
%Meshgrid for the plots of the control (only interior points)
[Xu,Tu]=meshgrid((OCP.a+(OCP.b-OCP.a)/OCP.N):(OCP.b-OCP.a)/OCP.N:(OCP.b-(OCP.b-OCP.a)/OCP.N),...
    0+OCP.T/OCP.Nt:OCP.T/OCP.Nt:OCP.T); 

f=1*ones(OCP.Nt*(OCP.N+1),1); %Values for the right hand side

hx=(OCP.b-OCP.a)/OCP.N; %Step length for the interval (a,b)
dt=OCP.T/OCP.Nt;        %Step length for the interval (0,T)

y0=zeros((OCP.N+1),1);  %Initial value for the state variable
yd=desiredState(OCP);   %Desired state

u=0*ones((OCP.N-1)*OCP.Nt,1);   %Initial guess for the control 
u_old=u;
y=forward_y(u,y0,f,OCP);        %Calculate the corresponding state variable y
y_old=y;
p=backward(y,yd,u,OCP);         %Calculate the corresponding adjoint state variable

Jk(1)=get_Jy(y,yd,u,OCP );      %Calculate the corresponding cost functional value

count_updates=1;                %Counter for the updates; counter_updates -1 for the number of updates
ctime=cputime;                  %Start time measurement

for k=1:kmax
    %Loop over all time points
    for j=1:OCP.Nt  
        %Loop over all interior space grid points
        for i=1:OCP.N-1 
            s=OCP.s;        
            pint=p((j-1)*(OCP.N+1)+1+i);
            yint=y((j-1)*(OCP.N+1)+1+i);
            uint=u((j-1)*(OCP.N-1)+i,1);
            %Candidates where the augmented Hamiltonian takes its minimum
            u1=min(max(OCP.u_lo,(2*epsilon*uint+pint*yint)/(2*(epsilon+0.5*OCP.alpha))),s);  
            u2=min(max(s,(2*epsilon*uint+(pint*yint-OCP.beta))/(2*(epsilon+0.5*OCP.alpha))),OCP.u_up);
            %Corresponding values of the augmented Hamiltonian
            H1=OCP.alpha*u1- pint*yint*u1+epsilon*(u1-uint)^2;                     
            H2=OCP.alpha*u2- pint*yint*u2+OCP.beta*u2+epsilon*(u2-uint)^2;
            %Choose the minimum value and where it is taken 
            [~,pos]=min([H1,H2]);  
           switch pos
               case 1
                   u((j-1)*(OCP.N-1)+i,1)=u1;
               case 2
                   u((j-1)*(OCP.N-1)+i,1)=u2;
           end
        end
    end
    du=sum((u-u_old).^2)*hx*dt;         %Evalute the norm square of the update
    y=forward_y(u,y0,f,OCP);            %Calculation of the state equation with the updated control u
    J_int=get_Jy(y,yd,u,OCP );          %Corresponding cost functional value
    if(J_int-Jk(count_updates)>-eta*du)  
        % If the improvemnt to the cost functional is not sufficiently large, 
        % take the old values of y and u, increase epsilon
        u=u_old;
        y=y_old;
        epsilon=epsilon*sigma;                            
    else                 
        % If the improvment to the cost functional is sufficently large, 
        % take the new values of y, u, J, calculate p and decrease epsilon
        count_updates=count_updates+1;
        p=backward(y,yd,u,OCP);
        epsilon=epsilon*zeta;                          
        Jk(count_updates)=J_int;
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
ct=cputime-ctime; %Stop time measurment for the calculation 
fprintf('CPU-time in seconds needed: %f\n',ct)

%Plot y and u
figure(1)
inty=reshape(y,[OCP.N+1,OCP.Nt+1])';
mesh(X,T,inty)
% title('The state variable');
xlabel('x');
ylabel('t');

  print('-depsc2', 'parabP6y.eps','-b0'); 
% print('-dpdf', 'parabP6y.pdf','-b0');

 
figure(2)
intu=reshape(u,[OCP.N-1,OCP.Nt])';
mesh(Xu,Tu,intu)
% title('The control variable');
xlabel('x');
ylabel('t');
  print('-depsc2', 'parabP6u.eps','-b0'); 
% print('-dpdf', 'parabP6u.pdf','-b0');

 
figure(3)
contour(Xu,Tu,intu)
% title('Contour of the control variable');
xlabel('x');
ylabel('t');

  print('-depsc2', 'parabP6ucontour.eps','-b0'); 
% print('-dpdf', 'parabP6contour.pdf','-b0');


%Plot convergence history 
figure(4)
semilogx(Jk,ltype{3},'Linewidth',2)
% title('Convergence history of J');

 print('-depsc2', 'parabP6Jhist.eps','-b0'); 
% print('-dpdf', 'parabP6Jhist.pdf','-b0');


% Check solution for optimality;
fprintf('\n')
fprintf('Check numerical optimality of the SQH solution\n')
numeric_optimality( u,y0,f,yd,OCP )


