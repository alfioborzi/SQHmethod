%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% solving a parabolic linear optimal control problem on the
% space-time cylinder [a,b]x[0,T] - discrete set of control values
% u in [u_lo: Du : u_up] where Du = (u_up - u_lo)/d_KU;
%
% Problem specific parameters: 
% N number of intervals for the discretisation of the space, 
% Nt for the discretation of the time, 
% (a,b) interval length of the domain, T time horizon, D diffusion constant, 
% u_lo and u_up lower and upper bound for the control, 
% alpha (weight quadratic cost term), 
% beta and s (weight discontinuous cost term) cost functional parameters

clear all;
close all; 

ltype = {'b-','r--','m-x','b-*','r:','m-.'};

OCP=struct('N',100,'Nt',200,'a',0,'b',1,'T',1,'D',1.0,'u_lo',-10,'u_up',10,'d_KU',10,'alpha',10^-3,'beta',0*10^-3,'P',0,'s',0);
d_KU=OCP.d_KU;

%Tolerance for convergence
kappa=10^-8;                
%Algorithm paramters
eta=10^-5;
zeta=0.9;
sigma=1.1;
kmax=200;         %Maximum iteration number
epsilon=1.0;         %Initial guess for epsilon

[X,T]=meshgrid(OCP.a:(OCP.b-OCP.a)/OCP.N:OCP.b,0:OCP.T/OCP.Nt:OCP.T);  %Meshgrid for the plots
%Meshgrid for the plots of the control (only interior points)
[Xu,Tu]=meshgrid((OCP.a+(OCP.b-OCP.a)/OCP.N):(OCP.b-OCP.a)/OCP.N:(OCP.b-(OCP.b-OCP.a)/OCP.N),0+OCP.T/OCP.Nt:OCP.T/OCP.Nt:OCP.T);

hx=(OCP.b-OCP.a)/OCP.N; 
dt=OCP.T/OCP.Nt; 

f=ones(OCP.Nt*(OCP.N+1),1); %Values for the right hand side
y0=zeros((OCP.N+1),1);
pT=zeros(OCP.N+1,1);
yd=desiredState3(OCP);



u=ones((OCP.N-1)*OCP.Nt,1);
u_old=u;
y=Forward(u,y0,f,OCP);
y_old=y;
p=Backward(y,yd,OCP);

Jk(1)=get_J(y,yd,u,OCP );

d=(OCP.b-OCP.a)*OCP.T*(1/(OCP.N*OCP.Nt));


count_updates=1;
ctime=cputime;
tic

for k=1:kmax
    for j=1:OCP.Nt
        for i=1:OCP.N-1
            u((j-1)*(OCP.N-1)+i,1)=argminH(i,j,d_KU,y,yd,p,u((j-1)*(OCP.N-1)+i,1),epsilon,f,OCP);      
        end
    end  
    du=sum((u-u_old).^2)*hx*dt;
   
   y=Forward(u,y0,f,OCP);
   J_int=get_J(y,yd,u,OCP );

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
        p=Backward(y,yd,OCP);
        epsilon=epsilon*zeta;                          
        Jk(count_updates)=J_int;
        u_old=u;
        y_old=y; 
        % Print current values
        fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',k,... 
            count_updates-1,J_int,epsilon,du) 
    end

   

    if(du<kappa && k > 200)      
        % If the norm square of the update is too small, then return the latest 
        % variable u that caused sufficient decrease of the cost functional 
        fprintf('Convergence achieved\n');
        u=u_old;
        break;
    end 
   
end
toc

ct=cputime-ctime; %Stop time measurment for the calculation 
fprintf('CPU-time in seconds needed: %f\n',ct)



%Plot y and u
figure(1)
inty=reshape(y,[OCP.N+1,OCP.Nt+1])';
mesh(X,T,inty)
% title('The state variable');
xlabel('x');
ylabel('t');

  print('-depsc2', 'parabP7y.eps','-b0'); 
% print('-dpdf', 'parabP7y.pdf','-b0');

 
figure(2)
intu=reshape(u,[OCP.N-1,OCP.Nt])';
mesh(Xu,Tu,intu)
% title('The control variable');
xlabel('x');
ylabel('t');
  print('-depsc2', 'parabP7u.eps','-b0'); 
% print('-dpdf', 'parabP7u.pdf','-b0');

 
figure(3)
contour(Xu,Tu,intu)
% title('Contour of the control variable');
xlabel('x');
ylabel('t');

  print('-depsc2', 'parabP7ucontour.eps','-b0'); 
% print('-dpdf', 'parabP7contour.pdf','-b0');


%Plot convergence history 
figure(4)
semilogx(Jk,ltype{3},'Linewidth',2)
% title('Convergence history of J');

 print('-depsc2', 'parabP7Jhist.eps','-b0'); 
% print('-dpdf', 'parabP7Jhist.pdf','-b0');

 figure(5)
 intyd=reshape(yd,[OCP.N+1,OCP.Nt+1])';
 mesh(X,T,intyd)
  xlabel('x');
  ylabel('t');
 view(-55,15);
