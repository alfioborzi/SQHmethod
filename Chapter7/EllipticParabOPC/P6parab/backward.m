%
%The function calculates the solution to the adjoint equation of a bilinear parabolic optimization problem
%-p'-D\Delta p +up =y-yd on (a,b)x(0,T) with an implicit Euler scheme, see
%corresponding paper Subsection 4.2 for details

function [ p] = backward(y,yd,u,OCP)

N=OCP.N;                %Number of intervals in space
h=(OCP.b-OCP.a)/N;
Nt=OCP.Nt;              %Number of intervals in time
dt=OCP.T/Nt;
D=OCP.D;                %Diffusion constant

p=zeros((N+1)*(Nt+1),1);
%Vectors for assembling the one dimensional discretised Laplacian
e=ones(N-1,1);
ed=(D*dt/h^2)*e;

for j=Nt:-1:1
    %Modify Laplacian due to the bilinear model
    A=spdiags([-ed e+2*ed+dt*u((j-1)*(N-1)+1:(j)*(N-1),1) -ed],-1:1, N-1,N-1); 
    %The update for each time step according to an implicit Euler scheme
    p((j-1)*(N+1)+2:j*(N+1)-1,1)=A\(p(j*(N+1)+2:(j+1)*(N+1)-1,1)...
                                +dt*(y((j)*(N+1)+2:(j+1)*(N+1)-1)-yd((j)*(N+1)+2:(j+1)*(N+1)-1))); 
end

end

