
function P = adjWave(U,x,t,Nx,Nt,alpha,beta)
%
global L T v

dx = L/Nx; 
dt = T/Nt; 
c = v*(dt/dx);    
c2  = c*c;
dt2 = dt*dt; 

% Adjoint problem
% d_tt p - v^2 d_xx p = -alpha*(y-y_d)
% terminal conditions p(T)= 0 , p_t(T)=beta*(y(T)-y_T); Dirichlet b.c. p=0
%-------------------------------------------------------------------------%
% 1. terminal condition p(T)=0
P(Nt,:) = 0*initp(x,T);
% 2. boundary conditions - homog. Dirichlet
tt = T; 
P(Nt,1) = 0*boundaryL(tt);
P(Nt,Nx+1) = 0*boundaryR(tt);
% 3. terminal condition p_t(T)
yt(1,:) = target(x,T); 
pt(1,:) = beta*(U(Nt,:) - yt(1,:));
yt(1,:) = - alpha*(U(Nt,:) - yt(1,:));
P(Nt-1,2:Nx)=(1-c2)*P(Nt,2:Nx)+0.5*c2*(P(Nt,1:Nx-1)+P(Nt,3:Nx+1)) ...
        - dt*pt(1,2:Nx) + 0.5*dt2*yt(1,2:Nx); 
% 4. boundary conditions - homog. Dirichlet
tt = T-dt; 
P(Nt-1,1) = 0*boundaryL(tt);
P(Nt-1,Nx+1) = 0*boundaryR(tt);
%-------------------------------------------------------------------------%
% Finite difference scheme
for j = Nt-1:-1:2  
% time step backward
     tt = t(j);
     yt(1,:) = target(x,tt); 
     yt(1,:) = -alpha*(U(j,:) - yt(1,:));
 P(j-1,2:Nx) = 2*(1-c2)*P(j,2:Nx)+c2*(P(j,1:Nx-1)+P(j,3:Nx+1)) ... 
               -P(j+1,2:Nx) + dt2*yt(1,2:Nx);    
           
% 2. boundary conditions - homog. Dirichlet
tt = t(j); 
P(j-1,1) = 0*boundaryL(tt);
P(j-1,Nx+1) = 0*boundaryR(tt);
end
