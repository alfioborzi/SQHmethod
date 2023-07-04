
function U = forWaveBCN(F,CB,x,t,Nx,Nt)
%
global L T v

dx = L/Nx; 
dt = T/Nt; 
c = v*(dt/dx);    
c2  = c*c;
dt2 = dt*dt; 

% Forward problem 
% d_tt y - v^2 d_xx y = F 
% initial conditions y(0)= y0 , y_t(0)=y1; Dirichlet 
% Neumann boundary control d_n y=u (CB)
%-------------------------------------------------------------------------%
% 1. initial condition y(0)=y0
U(1,:) = inity(x);

% 3. initial condition use y_t(0)=y1
yt(1,:) = inityt(x); 
U(2,2:Nx)=(1-c2)*U(1,2:Nx)+0.5*c2*(U(1,1:Nx-1)+U(1,3:Nx+1)) ...
         +dt*yt(1,2:Nx) + 0.5*dt2*F(1,2:Nx); 
% 4. boundary conditions - Neumann left
tt = dt; 
U(2,1)=(1-c2)*U(1,1)+0.5*c2*(2*dx*CB(1,1)+2*U(1,2)) ...
         +dt*yt(1,1) + 0.5*dt2*F(1,1); 
% U(2,1)=boundaryL(tt);
% 4. boundary conditions - Neumann right
U(2,Nx+1)=(1-c2)*U(1,Nx+1)+0.5*c2*(2*U(1,Nx)+2*dx*CB(1,2)) ...
         +dt*yt(1,Nx+1) + 0.5*dt2*F(1,Nx+1); 
% U(2,Nx+1) = boundaryR(tt);
%-------------------------------------------------------------------------%
% Finite difference scheme
for j = 2:Nt-1    
    tt = t(j);
% 1. time step forward
 U(j+1,2:Nx) = 2*(1-c2)*U(j,2:Nx)+c2*(U(j,1:Nx-1)+U(j,3:Nx+1)) ... 
               -U(j-1,2:Nx) + dt2*F(j,2:Nx);              
% 2. boundary conditions - Neumann left
U(j+1,1) = 2*(1-c2)*U(j,1)+c2*(2*dx*CB(j,1)+2*U(j,2)) ... 
               -U(j-1,1) + dt2*F(j,1);    

% 3. boundary conditions - Neumann right
 U(j+1,Nx+1) = 2*(1-c2)*U(j,Nx+1)+c2*(2*U(j,Nx)+2*dx*CB(j,2)) ... 
               -U(j-1,Nx+1) + dt2*F(j,Nx+1);     


end

end