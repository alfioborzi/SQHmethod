function [p] = backward(A,B,y,yd,u,OCP)
dt=OCP.dt;
Nt=round(OCP.T/dt);
N=max(size(yd));
p(1:N,Nt+1)=y-yd;
I=eye(N,N);

for n=Nt:-1:1
    C=transpose(2*A+(u(n)+u(n+1))*B);
    p(1:N,n)=(I-0.25*dt*C)\(p(1:N,n+1)+0.25*dt*C*p(1:N,n+1));
end

end

