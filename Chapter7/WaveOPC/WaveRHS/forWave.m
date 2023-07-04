
function U = forWave(F,x,t,Nx,Nt)
%
global L T v

dx = L/Nx; 
dt = T/Nt; 
c = v*(dt/dx);    
c2  = c*c;
dt2 = dt*dt; 

% Forward problem 
% d_tt y - v^2 d_xx y = u 
% initial conditions y(0)= y0 , y_t(0)=y1; Dirichlet b.c. y=0
%-------------------------------------------------------------------------%
% 1. initial condition y(0)=y0
U(1,:) = inity(x);
% 2. boundary conditions - Dirichlet
tt = 0; 
U(1,1)=boundaryL(tt);
U(1,Nx+1) = boundaryR(tt);

% 3. initial condition use y_t(0)=y1
yt(1,:) = inityt(x); 
U(2,2:Nx)=(1-c2)*U(1,2:Nx)+0.5*c2*(U(1,1:Nx-1)+U(1,3:Nx+1)) ...
         +dt*yt(1,2:Nx) + 0.5*dt2*F(1,2:Nx); 
% 4. boundary conditions - Dirichlet
tt = dt; 
U(2,1)=boundaryL(tt);
U(2,Nx+1) = boundaryR(tt);
%-------------------------------------------------------------------------%
% Finite difference scheme
for j = 2:Nt-1    
% 1. time step forward
 U(j+1,2:Nx) = 2*(1-c2)*U(j,2:Nx)+c2*(U(j,1:Nx-1)+U(j,3:Nx+1)) ... 
               -U(j-1,2:Nx) + dt2*F(j,2:Nx);              
% 2. boundary conditions - Dirichlet
tt = t(j); 
U(j+1,1)=boundaryL(tt);
U(j+1,Nx+1) = boundaryR(tt);
end

end