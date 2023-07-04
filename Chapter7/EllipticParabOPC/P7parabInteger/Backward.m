function [ p] = Backward(y,yd, OCP)
N=OCP.N;
h=(OCP.b-OCP.a)/N;
Nt=OCP.Nt;
dt=OCP.T/Nt;
D=OCP.D;

p=zeros((N+1)*(Nt+1),1);
e=ones(N-1,1);
ed=(D*dt/h^2)*e;
A=spdiags([-ed e+2*ed -ed],-1:1, N-1,N-1);

for j=Nt:-1:1
   p((j-1)*(N+1)+2:j*(N+1)-1,1)=A\(p(j*(N+1)+2:(j+1)*(N+1)-1,1)+dt*(y((j)*(N+1)+2:(j+1)*(N+1)-1)-yd((j)*(N+1)+2:(j+1)*(N+1)-1)));
end

end

