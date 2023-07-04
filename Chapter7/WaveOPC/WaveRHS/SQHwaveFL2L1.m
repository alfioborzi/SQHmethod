%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
%
% solving a linear optimal control of 1D wave equation 
% d_tt y - v^2 d_xx y = u 
% initial conditions y(0)= y0 , y_t(0)=y1
% Dirichlet b.c. y=0
%
% minimise 
% J(y,f) = alpha/2*||y-y_d||_Q^2 + beta/2*||y-y_T||_O^2 
%                                + nu/2*||u||_Q^2  + gamma*||u||_1
%
% finite difference method

close all;
clear all;

ltype = {'b-','r--','m-x','b-*','r:','m-.'};
%-------------------------------------------------------------------------%
% Initialization
global L T v

L  = 10.0;           % domain size
Nx = 100;            % nr. subintervals
dx = L/Nx;          % mesh size
x(:,1) = (0:Nx)*dx; % grid

mpx = Nx/2;     % Mid point of x axis
                          
T  = 1.0;          % Total number of time steps
Nt = 1000;        % nr. time intervals 
dt = T/Nt;        % time step
t(:,1)= (0:Nt)*dt;

v = 10;             % Wave velocity
c = v*(dt/dx);      % for CFL condition < 1

fprintf('\n CFL %f  \n',c);

% cost functional weights
alpha = 1.0;
beta  = 1.0;
nu    = 1.0e-5;
gamma = 1.0e-3;

% control's constrains

F_lo = -30;
F_up =  30;

% initialise
U = zeros(Nt,Nx+1);  % U(t,x)  state
P = zeros(Nt,Nx+1);  % P(t,x)  adjoint
F = zeros(Nt,Nx+1);  % F(t,X)  control distributed
Y = zeros(Nt,Nx+1);  % Y(t,X)  target
J = [];

c2  = c*c;
dt2 = dt*dt; 


% 1. right-hand side = control & target
for j=1:Nt
    tt = t(j);
F(j,:) = frhs(x,tt);
Y(j,:) = target(x,tt); 
end

% SQH optimization
%-------------------------------------------------------------------------%
% tolerance for convergence
kappa=10^-8;                
% SQH algorithm paramters
eta=10^-5;
zeta=0.9;
sigma=1.1;
kmax=5000;         %Maximum iteration number
epsilon=1.0;         %Initial guess for epsilon
%-------------------------------------------------------------------------%

F_old=F;           % initial control
U = forWave(F,x,t,Nx,Nt);        %Calculate the corresponding state 
U_old=U;

P = adjWave(U,x,t,Nx,Nt,alpha,beta);   %Calculate the corresponding adjoint 

J(1) = JcostL1(U,Y,F,alpha,beta,nu,gamma,Nx,Nt);   %Calculate value cost functional 

count_updates=1;        % Counter for the SQH updates
ctime=cputime;          % Start clock

for k=1:kmax
    
    %Loop over all time points
    for j=1:Nt  
        %Loop over all interior space grid points
        for i=2:Nx
                  
            pp=P(j,i);
            uu=F(j,i);
            
 %  F(j,i)=min(max((2*epsilon*uu+pp)/(nu+2*epsilon),F_lo),F_up); % if gamma=0
            
            %Candidates where the augmented Hamiltonian takes its max 
            u1=min(max(F_lo,(2*epsilon*uu+pp+gamma)/(2*epsilon+nu)),0);  
            u2=min(max(0,(2*epsilon*uu+pp-gamma)/(2*epsilon+nu)),F_up);

            %Corresponding values of the augmented Hamiltonian
            H1=pp*u1 - 0.5*nu*u1^2 - gamma*abs(u1) - epsilon*(u1-uu)^2;                     
            H2=pp*u2 - 0.5*nu*u2^2 - gamma*abs(u2) - epsilon*(u2-uu)^2;
                       
            %Choose the max value and where it is taken 
            [~,pos]=max([H1,H2]);  
            switch pos
               case 1
                   F(j,i)=u1;
               case 2
                   F(j,i)=u2;
            end
        end
    end
    
    dF=sum(sum((F-F_old).^2))*dx*dt;         %Evalute the norm square of the update
    
     
    U = forWave(F,x,t,Nx,Nt);  %Calculation of the state equation with the updated control
  
    Jold = JcostL1(U,Y,F,alpha,beta,nu,gamma,Nx,Nt); %Corresponding cost functional value
    
    if(Jold-J(count_updates) > - eta * dF)  
        % If the improvemnt to the cost functional is not sufficiently large, 
        % take the old values of U and F, increase epsilon
        F=F_old;
        U=U_old;
        epsilon=epsilon*sigma;                            
    else                 
        % If the improvment to the cost functional is sufficently large, 
        % take the new values of U, F, J, calculate P and decrease epsilon
        count_updates=count_updates+1;
        
        P = adjWave(U,x,t,Nx,Nt,alpha,beta);  
         
        epsilon=epsilon*zeta;                          
        J(count_updates)=Jold;
        F_old=F;
        U_old=U; 
        
%   % Print current values
    fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',k,count_updates-1,Jold,epsilon,dF)

    end
    if(dF < kappa)      
        % If the norm square of the update is too small, then return the latest 
        % variable u that caused sufficient decrease of the cost functional 
        fprintf('Convergence achieved \n');
        F=F_old;
        break;
    end 
   
end


ct=cputime-ctime; %Stop time measurment for the calculation 
fprintf('CPU-time in seconds needed: %f\n',ct)


%-------------------------------------------------------------------------%
% Plot 
plot_times = [Nt/2 3*Nt/4];
lpt = length(plot_times);



% plot state
for i = 1:lpt
  figure(i);
  k = plot_times(i);
  plot(x,U(k,:),ltype{1},'linewidth',2); hold on;
  plot(x,Y(k,:),ltype{2},'linewidth',2);
  grid on;
%  axis([min(x) max(x) -3 3]);
  legend('state','target');
  xlabel('x');
%  ylabel('state');              
%   titlestring = ['Time step = ',num2str(k), ' Time = ',num2str(t(k))];
%   title(titlestring ,'fontsize',12);  
%
filename = sprintf('waveytk%d.eps',k);
print('-depsc2',filename,'-b0'); 
end

% plot control
for i = 1:lpt
  figure(i+lpt);
  k = plot_times(i);
  plot(x,F(k,:),'linewidth',2); 
  grid on;
 % axis([min(x) max(x) -3 3]);
 % legend('control');
  xlabel('x','fontSize',14);
%  ylabel('control','fontSize',14);              
%  titlestring = ['Time step = ',num2str(k), ' Time = ',num2str(t(k))];
%  title(titlestring ,'fontsize',12);                            
%
filename = sprintf('waveutk%d.eps',k);
print('-depsc2',filename,'-b0'); 
end

% plot adjoint
for i = 1:lpt
  figure(i+2*lpt);
  k = plot_times(i);
  plot(x,P(k,:),'linewidth',2); 
  grid on;
 % axis([min(x) max(x) -3 3]);
 % legend('adjoint');
  xlabel('x','fontSize',14);
  ylabel('adjoint','fontSize',14);              
  titlestring = ['Time step = ',num2str(k), ' Time = ',num2str(t(k)), 'second'];
  title(titlestring ,'fontsize',12);                            
end

tm=t(1:end-1);
%Plot u 
figure(18)
mesh(tm,x,F')
%title('The control variable');
ylabel('x');
xlabel('t');
view([135 30])
print('-depsc2', 'waveu01.eps','-b0'); 
% print('-dpdf', 'waveu01.pdf','-b0');


%Plot y 
figure(19)
mesh(tm,x,U')
%title('The state variable');
ylabel('x');
xlabel('t');
view([135 30])
print('-depsc2', 'wavey01.eps','-b0'); 
% print('-dpdf', 'wavey01.pdf','-b0');


figure(20)
% plot(J,ltype{2},'Linewidth',2)
semilogx(J,ltype{3},'Linewidth',2)

print('-depsc2', 'waveJrhs01.eps','-b0'); 
% print('-dpdf', 'waveJrhs01.pdf','-b0');

fprintf('\n CFL %f  \n',c);


% print('-depsc2', 'wave01.eps','-b0'); 
% print('-dpdf', 'wave01.pdf','-b0');
