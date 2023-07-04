%
%The function uses an implicit Euler scheme to solve a bilinear equation with right hand side
%y'-D\Delta y +uy = f, see corresponding paper Subsection 4.2 for details

function [ y ] = forward_y(u,y0,f,OCP)

Nt=OCP.Nt;          %Number of intervals in time 
dt=OCP.T/Nt;
N=OCP.N;            %Number of intervals in space
h=(OCP.b-OCP.a)/N;
D=OCP.D;            %Diffusion constant

y=[y0;zeros((N+1)*Nt,1)];
%Vectors for assembling the one dimensional discretised Laplacian
e=ones(N-1,1);
ed=(D*dt/h^2)*e;

for j=2:Nt+1
    %Modify Laplacian due to the bilinear model
    A=spdiags([-ed e+2*ed+dt*u((j-2)*(N-1)+1:(j-1)*(N-1),1) -ed],-1:1, N-1,N-1); 
    %The update for each time step according to an implicit Euler scheme
    y((j-1)*(N+1)+2:j*(N+1)-1,1)=A\(y((j-2)*(N+1)+2:(j-1)*(N+1)-1,1)...
                                +dt*(f((j-2)*(N+1)+2:(j-1)*(N+1)-1,1)));  
end

end

