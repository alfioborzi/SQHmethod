%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for control of SDE via Fokker-Planck equation
%
% with u(x,t) control 
%
clear all; 
close all;

ltype = {'b-o','m-.','k:','r--'};

% Parameters

OCPA=struct('Nx',40,'Nt',80,'a',-2,'b',2,'T',2,'D',0.01,'umin',-4,'umax',4,...
            'alpha',10^(-4),'beta',0.0);

% Mesh sizes
dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);

T  = OCPA.T;
Nt = OCPA.Nt;
N  = OCPA.Nx;

x = linspace(OCPA.a,OCPA.b,N+1);
y = linspace(OCPA.a,OCPA.b,N+1);
t = linspace(0,T,Nt+1);

a = OCPA.a;
b = OCPA.b;

u_1=zeros(N+1,N+1,Nt+1);
u_2=zeros(N+1,N+1,Nt+1);

% Initialiaze controls

v(1,1:Nt+1) = 1.0*ones(1,Nt+1);
v(2,1:Nt+1) = pi*cos(pi*t(1:Nt+1)/2)/4;

w = 0.0*ones(2,Nt+1);
 v = 0.0*ones(2,Nt+1);


for i=1:Nt+1
    u_1(:,:,i)=1.0*(ones(N+1,N+1)*v(1,i)+((ones(N+1,1)*x)')*w(1,i));
    u_2(:,:,i)=1.0*(ones(N+1,N+1)*v(2,i)+(ones(N+1,1)*y)*w(2,i));
end

u_1_old=u_1;
u_2_old=u_2;



for k=1:N+1
for j=1:N+1
f0(k,j)=f00(x(k),y(j));
end
end

% Compute f

 f=fpsolve(u_1,u_2,@f00,OCPA);
 f_old=f;
 
 % showfp(f,OCPA,1);

p=adjfpsolve(u_1,u_2,@Vp,@xyd,OCPA);

% showadjfp(p,OCPA,8);


Jk(1)=Jfunc(f,u_1,u_2,@Vp,@xyd,OCPA);


kmax=200;
eta = 10^-8; 
sigma = 20;
zeta = 0.3;
epsi=100.0;

count_updates=1;
p_value(count_updates)=1;

Vs=(OCPA.b-OCPA.a)^2;



for k=1:kmax
    
    for i=1:Nt+1
   
        fint=f(2:N,2:N,i);
    
    %Aproximation of the derivative of the adjoint variable in the
    %augmented Hamiltonian
    dp1=(p(3:N+1,2:N,i)-p(1:N-1,2:N,i))/(2*dx);
    dp2=(p(2:N,3:N+1,i)-p(2:N,1:N-1,i))/(2*dx);
    

    u11=min(max(OCPA.umin,(2*epsi*u_1(2:N,2:N,i)-dp1)./(2*epsi+OCPA.alpha)),OCPA.umax);
    u22=min(max(OCPA.umin,(2*epsi*u_2(2:N,2:N,i)-dp2)./(2*epsi+OCPA.alpha)),OCPA.umax);
  
    
    
    for j=1:N-1
        for l=1:N-1
            
            u_1(j+1,l+1,i)=u11(j,l);
            u_2(j+1,l+1,i)=u22(j,l);

        end
    end
    
        
end
    
 
   
   intu=epsi*sum(sum(sum(((u_1-u_1_old).^2)+((u_2-u_2_old).^2))))*dt*(dx^2);

   
   f=fpsolve(u_1,u_2,@f00,OCPA);
   
   
   J_int = Jfunc(f,u_1,u_2,@Vp,@xyd,OCPA);
   
   dJ = J_int-Jk(count_updates);
   
   if(J_int-Jk(count_updates)>-eta*intu)
%       disp('no reduction');
       u_1=u_1_old;
       u_2=u_2_old;
       f=f_old;
       epsi=epsi*sigma;  
       if(epsi > 10^8) 
           break; 
       end
   else
 %      disp('ok reduction');
       count_updates=count_updates+1;
       
       p=adjfpsolve(u_1,u_2,@Vp,@xyd,OCPA);
       
       epsi=epsi*zeta;    
       
       Jk(count_updates)=J_int;
       
       u_1_old=u_1;
       u_2_old=u_2;
       f_old=f;
       
       p_value(count_updates)=sum(sum(p(:,:,1).*f0))*dx^2;
       
       % Print current values
    fprintf('k %i | k up %i | J %e | eps %e\n',k,... 
            count_updates-1,J_int,epsi) 
    
       
   end

   if(intu<10^-14)
      disp('convergence');
        u_1=u_1_old;
        u_2=u_2_old;
        break;
   end 


end

%Comparision of the Value function and the cost functional at the initial time
disp('Value function and the Cost functional at the initial time')
fprintf('int(p(0)f0),\t J\n\n')
for i=1:count_updates
    fprintf('%i \t %d \t %d\n',i,p_value(i),Jk(i))
end



% check total Prob
Intf0 = sum(sum(f(:,:,1)))*(dx^2);
IntfT = sum(sum(f(:,:,end)))*(dx^2);

[X,Y]=meshgrid(a:dx:b,a:dx:b);
 figure(1)
 i1=40;
 i2=60;
 subplot(2,2,1)
 mesh(X,Y,u_1(:,:,i1))
 axis([a b a b -4 4]);
 title('$u_1$','interpreter','latex');
 subplot(2,2,2)
 mesh(X,Y,u_2(:,:,i1))
 axis([a b a b -4 4]);
 title('$u_2$','interpreter','latex')
 subplot(2,2,3)
 mesh(X,Y,u_1(:,:,i2))
 axis([a b a b -4 4]);
 title('$u_1$','interpreter','latex');
 subplot(2,2,4)
 mesh(X,Y,u_2(:,:,i2))
 axis([a b a b -4 4]);
 title('$u_2$','interpreter','latex')
 
  print('-depsc2', 'controls-u.eps','-b0'); 
 print('-dpdf', 'controls-u.pdf','-b0');
 
 % ----- for print 
 figure(20)
 mesh(X,Y,u_1(:,:,i1))
 axis([a b a b -4 4]);
% title('$u_1$','interpreter','latex');
  print('-depsc2', 'controls-u1-t10.eps','-b0','-r300'); 
  
   figure(21)
 mesh(X,Y,u_2(:,:,i1))
 axis([a b a b -4 4]);
% title('$u_2$','interpreter','latex')
  print('-depsc2', 'controls-u2-t10.eps','-b0','-r300'); 
  
figure(22)
 mesh(X,Y,u_1(:,:,i2))
 axis([a b a b -4 4]);
% title('$u_1$','interpreter','latex');
  print('-depsc2', 'controls-u1-t15.eps','-b0','-r300'); 
  
  figure(23)
 mesh(X,Y,u_2(:,:,i2))
 axis([a b a b -4 4]);
% title('$u_2$','interpreter','latex') 
  print('-depsc2', 'controls-u2-t15.eps','-b0','-r300'); 

 
 % ----- for print 


%showadjfp(p,OCPA,2);

%showfp(f,OCPA,3);

figure(4)

plot(Jk,'Linewidth',2)

 print('-depsc2', 'J-u.eps','-b0','-r300'); 
 
 print('-dpdf', 'J-u.pdf','-b0','-r300');

% title('cost functional')


%Calculation of the mean value of f


for i=1:Nt+1
    
    fzwi=f(:,:,i);
    
    zwsumme1=0;
    zwsumme2=0;
    for j=1:N+1
         zwsumme1=zwsumme1+(a+dx*(j-1))*sum(fzwi(j,:))*dx^2;
         zwsumme2=zwsumme2+(a+dx*(j-1))*sum(fzwi(:,j))*dx^2;
    end
    mv(1,i)=zwsumme1;
    mv(2,i)=zwsumme2;
end

figure(5)
hold on
plot(mv(1,:),mv(2,:),ltype{1},'Linewidth',2);
 axis([a b a b]);
 
  [xt, yt] = xyd(t);
 plot(xt(:),yt(:),ltype{4},'Linewidth',2)
 
 print('-depsc2', 'mean-u.eps','-b0','-r300'); 
 print('-dpdf', 'mean-u.pdf','-b0','-r300');


 % title('f mean value')

%Save the control values
[a,b,m,h,x,y,Nt,t,dt,alpha,beta,nu] = parameters(1,OCPA);

ip = xyd(t(1)); 

T1 = T;
T0 = 0;
Nx = N;
fname = 'data_traj.mat';
save(fname,'a','b','Nx','T0','T1','Nt','u_1','u_2','ip');

MonteCarloSim(OCPA);

i=Nt;
uvel = u_1(:,:,i);
vvel = u_2(:,:,i);

figure(7)
quiver(uvel,vvel);

 

%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outf0] = f00(x,y)
    sigma=0.3;
    sigma2=2*sigma*sigma;
    x0 = -1.0;
    y0 = 0.0;
    outf0 = exp(-((x-x0).^2+(y-y0).^2)/(sigma2))/(pi*sigma2) ; 
end

function [outVp] = Vp(x,y,xt,yt)
    sigma=0.3;
    sigma2=2*sigma*sigma;
    Acoeff = 1.0d-3/(pi*sigma2);
outVp = - Acoeff * exp(-((x-xt)^2+(y-yt)^2)/(sigma2)) ; 
end


function [outX, outY] = xyd(t)
    x_t =  t - 1.0 ; 
    y_t =  sin(pi*t/2);
outX=x_t;
outY=y_t; 
end

function [ J ] = Jfunc(f,u_1,u_2,Vp,xyd,OCPA)

dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);
h = dx; 

T  = OCPA.T;
Nt = OCPA.Nt;
N  = OCPA.Nx;

x = linspace(OCPA.a,OCPA.b,N+1);
y = linspace(OCPA.a,OCPA.b,N+1);
t = linspace(0,T,Nt+1);


J = 0;

% space - time integration
for i=1:Nt
    [xt,yt] = xyd(t(i));

    for k=2:N
    for j=2:N
        
         dkji=0.5*OCPA.alpha*(u_1(k,j,i).^2+u_2(k,j,i).^2); 
            
    J = J + (Vp(x(k),y(j),xt,yt) + dkji)*f(k,j,i);
    end
   end
end

J = J * dx*dx*dt ; 

% space integration at t=T
i=Nt+1;
    [xt,yt] = xyd(t(i));
   for k=2:N
    for j=2:N
    J = J + Vp(x(k),y(j),xt,yt)*f(k,j,i)*dx*dx;
    end
   end

end








function [  ] = showfp( cd,OCPA, k )
a=OCPA.a;
b=OCPA.b;
T=OCPA.T;

dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);

N=OCPA.Nx;
Nt=OCPA.Nt;

[X,Y]=meshgrid(a:dx:b,a:dx:b);

    figure(k)
    
for i=1:Nt+1

    subplot(1,1,1);
    surf(X,Y,cd(:,:,i)')
    view(0,90);
    axis([a b a b -1 3]);
    pause(0.1)
end

end


function [  ] = showadjfp( cd,OCPA, k )
a=OCPA.a;
b=OCPA.b;
T=OCPA.T;

dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);

N=OCPA.Nx;
Nt=OCPA.Nt;

[X,Y]=meshgrid(a:dx:b,a:dx:b);

figure(k)


for i=Nt+1:-1:1
    subplot(1,1,1);
    mesh(X,Y,cd(:,:,i)')
    view(0,90);
    axis([a b a b -1 1]);
    pause(0.1)
end

end

function [] = MonteCarloSim(OCPA)

fname='data_traj.mat';
load(fname);   % all data needed for the simulation

D=OCPA.D;

% number of desired trajectories
Nsample=10;
Ntime = 300; % number of time points for the Milstein integration


T=T1-T0;

[ Nx1, Nx2, Nt1] = size(u_1);

D1 = linspace(a,b,Nx1);
D2 = linspace(a,b,Nx2);
TT = linspace(0,T,Nt1);

% interpolate controls: ???

model

%Starting points of the stochastic process
% as open loop
[x0, y0] = xyd(0);

% closed loop check
%x0=1;
%y0=1;

% Integrate the stochastic process
[y, t] = milstein_2D(m,x0,y0,u_1,u_2,Nsample,T,Ntime);



figure(6)

hold on

 plot(y(:,:,1),y(:,:,2),'Linewidth',1);
 hold on;
 axis([a b a b]);
 

 [xt, yt] = xyd(TT);
 
 plot(xt(:),yt(:),'--','Linewidth',1)
 
 print('-depsc2', 'trajs-u-B.eps','-b0','-r300'); 
 print('-dpdf', 'trajs-u-B.pdf','-b0','-r300');


end

