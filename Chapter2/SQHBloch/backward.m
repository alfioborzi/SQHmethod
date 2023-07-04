function [p] = backward(A,B1,B2,y,yd,u,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

N=max(size(yd));
p(1:N,Nt+1)=-(y(:,end)-yd);
I=eye(N,N);

% Midpoint in the bilinear case 
for n=Nt:-1:1
    C=transpose( 2*A + (u(1,n)+u(1,n+1))*B1 + (u(2,n)+u(2,n+1))*B2);
    p(1:N,n)=(I-0.25*dt*C)\(p(1:N,n+1)+0.25*dt*C*p(1:N,n+1));
end

end

