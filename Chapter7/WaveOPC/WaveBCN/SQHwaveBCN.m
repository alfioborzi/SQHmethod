%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% solving a linear optimal control of 1D wave equation 
% d_tt y - v^2 d_xx y = f  
% initial conditions y(0)= y0 , y_t(0)=y1
%
% Neumann boundary control  d_n y=u 
%
% minimise 
% J(y,f) = alpha/2*||y-y_d||_Q^2 + beta/2*||y-y_T||_O^2 
%                                + nu/2*||u||_S^2  + gamma*||u||_1
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
                          
T  = 4;          % Total number of time steps
Nt = 1000;        % nr. time intervals 
dt = T/Nt;        % time step
t(:,1)= (0:Nt)*dt;

v = 10;             % Wave velocity
v2= v*v;
c = v*(dt/dx);      % for CFL condition < 1

fprintf('\n CFL %f  \n',c);

% cost functional weights
alpha = 0.0;
beta  = 1.0;
nu    = 1.0e-12;
gamma = 2.0e-1;

% control's constrains

F_lo = -0.04;
F_up =  0.04;

% initialise
U = zeros(Nt,Nx+1);  % U(t,x)  state
P = zeros(Nt,Nx+1);  % P(t,x)  adjoint
F = zeros(Nt,Nx+1);  % F(t,X)  forcing term
Y = zeros(Nt,Nx+1);  % Y(t,X)  target
CB = zeros(Nt,2);    % control left & right boundary


J = [];
E = [];

c2  = c*c;
dt2 = dt*dt; 


% 1. right-hand side, target, boundary controls
for j=1:Nt
    tt = t(j);
F(j,:) = frhs(x,tt);
Y(j,:) = target(x,tt); 
CB(j,1) = 0;
CB(j,2) = 0;
end


% SQH optimization
%-------------------------------------------------------------------------%
% tolerance for convergence
kappa=10^-10;                
% SQH algorithm paramters
eta=10^-5;
zeta=0.9;
sigma=1.1;
kmax=5000;         %Maximum iteration number
epsilon=1.0;         %Initial guess for epsilon
%-------------------------------------------------------------------------%

CB_old=CB;           % initial control
U = forWaveBCN(F,CB,x,t,Nx,Nt);        %Calculate the corresponding state 
U_old=U;

P = adjWaveBCN(U,x,t,Nx,Nt,alpha,beta);   %Calculate the corresponding adjoint 

J(1) = JcostBCN(U,Y,CB,alpha,beta,nu,gamma,Nx,Nt);   %Calculate value cost functional 

count_updates=1;        % Counter for the updates
ctime=cputime;          % Start clock

for k=1:kmax
    
    %Loop over all time points
    for j=2:Nt  
        
        % left boundary
        i=1;
                  
            pp=P(j,i);
            uu=CB(j,1);
            
            
            %Candidates where the augmented Hamiltonian takes its max 
            u1=min(max(F_lo,(2*epsilon*uu+v2*pp+gamma)/(2*epsilon+nu)),0);  
            u2=min(max(0,(2*epsilon*uu+v2*pp-gamma)/(2*epsilon+nu)),F_up);

            %Corresponding values of the augmented Hamiltonian
            H1=v2*pp*u1 - 0.5*nu*u1^2 - gamma*abs(u1) - epsilon*(u1-uu)^2;                     
            H2=v2*pp*u2 - 0.5*nu*u2^2 - gamma*abs(u2) - epsilon*(u2-uu)^2;
                       
            %Choose the max value and where it is taken 
            [~,pos]=max([H1,H2]);  
            switch pos
               case 1
                 CB(j,1)=u1;
               case 2
                 CB(j,1)=u2;
            end
        
            % right boundary
            i=Nx+1;
                  
            pp=P(j,i);
            uu=CB(j,2);
            
            
            %Candidates where the augmented Hamiltonian takes its max 
            u1=min(max(F_lo,(2*epsilon*uu+v2*pp+gamma)/(2*epsilon+nu)),0);  
            u2=min(max(0,(2*epsilon*uu+v2*pp-gamma)/(2*epsilon+nu)),F_up);

            %Corresponding values of the augmented Hamiltonian
            H1=v2*pp*u1 - 0.5*nu*u1^2 - gamma*abs(u1) - epsilon*(u1-uu)^2;                     
            H2=v2*pp*u2 - 0.5*nu*u2^2 - gamma*abs(u2) - epsilon*(u2-uu)^2;
                       
            %Choose the max value and where it is taken 
            [~,pos]=max([H1,H2]);  
            switch pos
               case 1
                 CB(j,2)=u1;
               case 2
                 CB(j,2)=u2;
            end
    end
    
    dF=sum(sum((CB-CB_old).^2))*dt;         %Evalute the norm square of the update
    
     
    U = forWaveBCN(F,CB,x,t,Nx,Nt);  %Calculation of the state equation with the updated control
  
    Jold = JcostBCN(U,Y,CB,alpha,beta,nu,gamma,Nx,Nt); %Corresponding cost functional value
    
    if(Jold-J(count_updates) > - eta * dF)  
        % If the improvemnt to the cost functional is not sufficiently large, 
        % take the old values of U and F, increase epsilon
        CB=CB_old;
        U=U_old;
        epsilon=epsilon*sigma;                            
    else                 
        % If the improvment to the cost functional is sufficently large, 
        % take the new values of U, F, J, calculate P and decrease epsilon
        count_updates=count_updates+1;
        
        P = adjWaveBCN(U,x,t,Nx,Nt,alpha,beta);  
         
        epsilon=epsilon*zeta;                          
        J(count_updates)=Jold;
        CB_old=CB;
        U_old=U; 
        
%   % Print current values
    fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',k,count_updates-1,Jold,epsilon,dF)

    end
    if(dF < kappa)      
        % If the norm square of the update is too small, then return the latest 
        % variable u that caused sufficient decrease of the cost functional 
        fprintf('Convergence achieved \n');
        CB=CB_old;
        break;
    end 
   
end


ct=cputime-ctime; %Stop time measurment for the calculation 
fprintf('CPU-time in seconds needed: %f\n',ct)

% Energy
for j=2:Nt
        E(j) = EnergyBCN(U,j,Nx);
end

%-------------------------------------------------------------------------%
% Plot 
plot_times = [Nt/2 Nt];
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
filename = sprintf('waveytkBCN%d.eps',k);
print('-depsc2',filename,'-b0'); 
end

% plot control
tm=t(1:end-1);


  figure(1+lpt);
  plot(tm,CB(:,1),'linewidth',2); 
  xlabel('t','fontSize',14);
%  title('control left');
  print('-depsc2','waveuBCNleft.eps','-b0'); 
  
  figure(2+lpt);
  plot(tm,CB(:,2),'linewidth',2);
  xlabel('t','fontSize',14);
%  title('control right');
  print('-depsc2','waveuBCNright.eps','-b0'); 


% % plot adjoint
% for i = 1:lpt
%   figure(i+2*lpt);
%   k = plot_times(i);
%   plot(x,P(k,:),'linewidth',2); 
%   grid on;
%  % axis([min(x) max(x) -3 3]);
%  % legend('adjoint');
%   xlabel('x','fontSize',14);
%   ylabel('adjoint','fontSize',14);              
%   titlestring = ['Time step = ',num2str(k), ' Time = ',num2str(t(k)), 'second'];
%   title(titlestring ,'fontsize',12);                            
% end


%Plot y 
figure(10)
mesh(tm,x,U')
%title('The state variable');
ylabel('x');
xlabel('t');
view([135 30])
print('-depsc2', 'waveyBCN01.eps','-b0'); 
% print('-dpdf', 'waveyBCN01.pdf','-b0');

figure(15)
% plot(J,ltype{2},'Linewidth',2)
semilogx(E,ltype{3},'Linewidth',2);

print('-depsc2', 'waveEbcn01.eps','-b0'); 
% print('-dpdf', 'waveEbcn01.pdf','-b0');



figure(20)
% plot(J,ltype{2},'Linewidth',2)
semilogx(J,ltype{3},'Linewidth',2);

print('-depsc2', 'waveJbcn01.eps','-b0'); 
% print('-dpdf', 'waveJbcn01.pdf','-b0');

fprintf('\n CFL %f  \n',c);
