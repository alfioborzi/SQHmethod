function [ y ] = Forward(u,y0,f,OCP)
Nt=OCP.Nt;
dt=OCP.T/Nt;
N=OCP.N;
h=(OCP.b-OCP.a)/N;
D=OCP.D;

y=[y0;zeros((N+1)*Nt,1)];
e=ones(N-1,1);
ed=(D*dt/h^2)*e;
A=spdiags([-ed e+2*ed -ed],-1:1, N-1,N-1);

for j=2:Nt+1
   
   y((j-1)*(N+1)+2:j*(N+1)-1,1)=A\(y((j-2)*(N+1)+2:(j-1)*(N+1)-1,1)+dt*(u((j-2)*(N-1)+1:(j-1)*(N-1),1)+f((j-2)*(N+1)+2:(j-1)*(N+1)-1,1)));
end

end

